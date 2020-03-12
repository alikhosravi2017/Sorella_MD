//#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <time.h>
#include "gasdev.cpp"
#include "other/nr.h"
#define  LengthX 3.18198052
#define  LengthY 3.18198052
#define  LengthZ 4.5
//#define  displacement  1.5
#define  Number 54    //number of atoms should be   N=2*(n)^3 ==>> "n" should be integer
#define  Cutoff 5   //cut off
#define  dump_every 100
#define  log_every  100  
#define  thermostat_every  100000000
#define  h          0.0005
#define  desired_Temperature          20 //K
#define  kB_True         1.38064852e-23  //m2 kg s-2 K-1
#define  epsilon_True    1.666e-21 //J
#define  sigma_True      3.4e−10    //m
/* MD units 
  kB         1 
  epsilon    1
  sigma      1
 
 * */
using namespace std;
void Force(float [][3],float [][3]);
float KEnergy(float [][3]);
float PEnergy(float [][3]);
void velocity_verlet(float [][3],float [][3],float [][3], int);
void dump(float [][3],int, ofstream& );
void log(float[][3],float[][3],int, ofstream&);
void pbc(float[][3], int);
float Temperature(float [][3]);
void Thermostat_velocity_scaling(float [][3]);
int main(){
    
    float X[Number][3],V[Number][3],F[Number][3];
    float r,T,t,m,YMassVelocity=0,XMassVelocity=0 ,ZMassVelocity=0;
    float X_cm=0, Y_cm=0, Z_cm=0;
    int i,j,k,N,a,step=0,idum=-8; 
    ofstream M_output_file ("dump.xyz");
    ofstream E_output_file ("log.dat");
    cout<<"inter full time:\n ";
//     cin>>T;
    T = 1000000;

// Old positioning block!!!!!!!!!

//    k=sqrt(Number/2)+1;
////     m=(Length/2)/(sqrt(Number/2)+1);
//    m= displacement;
//    N=Number;
//        for(j=0;j<=Number/k;j++){
//            for(i=1;i<=k;i++){
//                if (N!=0){
//                X[Number-N][0]=i*m;
//                X[Number-N][1]=(j+1)*m;
//                X_cm += X[Number-N][0]/Number;
//                Y_cm += X[Number-N][1]/Number;
//                N=N-1;
//                }
//            }
//    }

//new positioning
    system( (" python create_argon.py -n '"+to_string(int(pow(Number/2,1./3)))+" '").c_str() ); // '"+name+"  '
    ifstream myfile ("initial_positions_"+to_string(Number)+"_atoms.xyz");
    float boxsize_X, boxsize_Y, boxsize_Z;
    myfile >> boxsize_X;
    myfile >> boxsize_Y;
    myfile >> boxsize_Z;
//    cout<<"yes!\t"<<boxsize_X<<"\t"<<boxsize_Y<<"\t"<<boxsize_Z<<"\n";
    for(j=0;j<Number;j++){
        myfile >> X[j][0];
        myfile >> X[j][1];
        myfile >> X[j][2];
        X_cm += X[j][0]/Number;
        Y_cm += X[j][1]/Number;
        Z_cm += X[j][2]/Number;
        }
    myfile.close();
//


//     cout<<X_cm<<"\t"<<Y_cm<<endl;
//     // move center of mass to the middle:
//    for(i=0;i<Number;i++){
//        X[i][0] += (LengthX/2-X_cm);
//        X[i][1] += (LengthY/2-Y_cm);
//        X[i][1] += (LengthZ/2-Z_cm);
//    }

    dump(X,step,M_output_file);
    //end of first place calculation...
    
    for(i=0;i<Number;i++){
        V[i][0]=NR::gasdev(idum);
        V[i][1]=NR::gasdev(idum);
        V[i][2]=NR::gasdev(idum);
    }
	for(i=0;i<Number;i++){                       //calculate Central Mass Velocity
		XMassVelocity+=V[i][0];
		YMassVelocity+=V[i][1];
		ZMassVelocity+=V[i][2];
	}
	for(i=0;i<Number;i++){                       //Central Mass Velocity=0
		V[i][2]=V[i][2]-ZMassVelocity/Number;
		V[i][1]=V[i][1]-YMassVelocity/Number;
		V[i][0]=V[i][0]-XMassVelocity/Number;
	}
	
	Thermostat_velocity_scaling(V);
	
	
	Force(X,F); // calculate a0
	
        /// MAIN ///
        clock_t tStart = clock();
        for(step=1;step<T;step+=1){
            velocity_verlet(V,X,F, step);
            dump(X,step,M_output_file);
            log(X,V,step,E_output_file);
            if(step%thermostat_every==0)        Thermostat_velocity_scaling(V);
        }
        
    
        printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
        
return 0;
}

void pbc(float X[][3], int i){
    X[i][0] -= LengthX*floor(float(X[i][0]/LengthX));
    X[i][1] -= LengthY*floor(float(X[i][1]/LengthY));
    X[i][2] -= LengthZ*floor(float(X[i][2]/LengthZ));

//     if(X[i][0]>Length)  X[i][0]+=(-1)*Length;    //boundary conditions...
//     if(X[i][0]<0)       X[i][0]+=Length;		 	//boundary conditions...	
//     if(X[i][1]>Length)  X[i][1]+=(-1)*Length;    //boundary conditions...
//     if(X[i][1]<0)       X[i][1]+=Length;         //boundary conditions...
}


void Force(float X[][3],float F[][3]){             //calculate acceleration
    int i,j;
    float r2,deltax,deltay, deltaz, solid_distanse_x,solid_distanse_y ,solid_distanse_z;
    float LengthX_2 = float(LengthX/2);
    float LengthY_2 = float(LengthY/2);
    float LengthZ_2 = float(LengthZ/2);

    for(i=0;i<Number;i++){
        F[i][0]=0;
        F[i][1]=0;
        F[i][2]=0;
//        cout<<X[i][0]<<X[i][1]<<"\n";
    }
//    cout<<"\n\n\n\n";
	for(i=0;i<Number;i++){
	    for(j=0;j<Number;j++){
                if(j!=i){ //this IF is unnecessary instead set i>j in the loop
                    
                    solid_distanse_x = abs(X[i][0]-X[j][0]);
                    deltax = -solid_distanse_x + LengthX*floor(float(solid_distanse_x/LengthX_2));
                                
                    solid_distanse_y = abs(X[i][1]-X[j][1]);
                    deltay=  -solid_distanse_y + LengthY*floor(float(solid_distanse_y/LengthX_2));

                    solid_distanse_z = abs(X[i][2]-X[j][2]);
                    deltaz=  -solid_distanse_z + LengthZ*floor(float(solid_distanse_z/LengthX_2));

                    r2=pow(deltax,2)+ pow(deltay,2) +pow(deltaz,2);
                    
                    
                    if(sqrt(r2)<Cutoff){
                        if (X[i][0]<X[j][0])     F[i][0]+=(deltax*48*(1/pow(r2,7)-1/(2*pow(r2,4))));   //  note that  r2=r*r
                        else                     F[i][0]+=(-1*deltax*48*(1/pow(r2,7)-1/(2*pow(r2,4))));

                        if (X[i][1]<X[j][1])     F[i][1]+=(deltay*48*(1/pow(r2,7)-1/(2*pow(r2,4))));
                        else                     F[i][1]+=(-1*deltay*48*(1/pow(r2,7)-1/(2*pow(r2,4))));

                        if (X[i][2]<X[j][2])     F[i][2]+=(deltaz*48*(1/pow(r2,7)-1/(2*pow(r2,4))));
                        else                     F[i][2]+=(-1*deltaz*48*(1/pow(r2,7)-1/(2*pow(r2,4))));
                    }
                }
            }
    }
}



// void Force(float X[][3],float F[][3]){             //calculate acceleration
//     int i,j,marzx,marzy;
//     float r2,deltax,deltay;
//     for(i=0;i<Number;i++){
// 		F[i][0]=0;
// 		F[i][1]=0;
//     }
// 	for(i=0;i<Number;i++){
// 	    for(j=0;j<Number;j++){
//     	    if(j!=i){
//         	    if(abs(X[i][0]-X[j][0])>Length/2){
//             	  deltax=Length-abs(X[i][0]-X[j][0]);
// 				  marzx=0;
// 				}
//             	else{
//                 	deltax=abs(X[i][0]-X[j][0]);
//                 	marzx=1;
// 				}
//             	if(abs(X[i][1]-X[j][1])>Length/2){
//             		deltay=Length-abs(X[i][1]-X[j][1]);
//               		marzy=0;
// 				}
//               	else{
//                 	deltay=abs(X[i][1]-X[j][1]);
//                 	marzy=1;
// 				}
//         		r2=pow(deltax,2)+pow(deltay,2);
//         		if(sqrt(r2)<Cutoff){
//                 	if (X[i][0]<X[j][0]) F[i][0]+=(pow(-1,marzx)*deltax*48*(1/pow(r2,7)-1/(2*pow(r2,4))));   //  note that  r2=r*r
//                 	else F[i][0]+=(pow(-1,marzx+1)*deltax*48*(1/pow(r2,7)-1/(2*pow(r2,4))));
//                 	if (X[i][1]<X[j][1])  F[i][1]+=(pow(-1,marzy)*deltay*48*(1/pow(r2,7)-1/(2*pow(r2,4))));
//                 	else F[i][1]+=(pow(-1,marzy+1)*deltay*48*(1/pow(r2,7)-1/(2*pow(r2,4))));
// 				}
// 			}
// 		}
//                
// 	}
// }

void Thermostat_velocity_scaling(float V[][3]){
    float lambda;
    int i;
    float temperature_True = float(epsilon_True)/kB_True;  // multiply temperature with this value to take in Kelvin

    float temperature_now =Temperature(V)*temperature_True;
//     cout<<"temperature_now="<<temperature_now<<endl;
    lambda = sqrt(desired_Temperature/temperature_now);
    for(i=0;i<Number;i++){
        V[i][0] *= lambda;
        V[i][1] *= lambda;
        V[i][2] *= lambda;
    }
    }

void velocity_verlet(float V[][3], float X[][3], float F[][3], int step){
    	                                       
    // verlet loop //
    int i;
    for(i=0;i<Number;i++){                             //calculate v half step
        V[i][0]+=(h/2)*(F[i][0]);
        V[i][1]+=(h/2)*(F[i][1]);
        V[i][2]+=(h/2)*(F[i][2]);
        X[i][0]+=V[i][0]*h;
        X[i][1]+=V[i][1]*h;
        X[i][2]+=V[i][2]*h;
        pbc(X,i);
    }
    Force(X,F);                                         //calculate a
    for(i=0;i<Number;i++){                               //calculate v half step
        V[i][0]+=(h/2)*(F[i][0]);
        V[i][1]+=(h/2)*(F[i][1]);
        V[i][2]+=(h/2)*(F[i][2]);
    }
        //cout<<step<<endl;
    }

void log(float X[][3],float V[][3],int step, ofstream& E_output_file){
    //cout<<step%log_every<<endl;
    if(step%log_every==0){
        float temrature_True = float(epsilon_True)/kB_True;  // multiply temperature with this value to take in Kelvin
        float kinetic_E   = KEnergy(V)*epsilon_True;       //calculate kinetic Energy
        float potential_E = PEnergy(X)*epsilon_True;    	//calculate potential Energy
        float temperature_now = Temperature(V)*temrature_True;
        E_output_file<<step<<"\t"<<kinetic_E<<"\t"<<potential_E<<"\t"<<kinetic_E+potential_E<<"\t"<<temperature_now<<endl;
        //cout<<"loged at step="<<step<<endl;
    }
    }

void dump(float X[][3], int step, ofstream& M_output_file){
    if(step%dump_every==0){
        M_output_file<<Number<<endl<<"atoms"<<endl; 
        for(int i=0;i<Number;i++){
            M_output_file<<1<<"\t"<<X[i][0]<<"\t"<<X[i][1]<<"\t"<<X[i][2]<<endl;
        }
//         cout<<"dumped at step="<<step<<endl;
    }
    }

float Temperature(float V[][3]){
    return KEnergy(V) * float(1)/(float(3)/2*Number);
    }

float KEnergy(float V[][3]){           //calculate kinetic Energy
    int i;
    float E=0;
    for(i=0;i<Number;i++){
    	E+=.5*(pow(V[i][0],2)+pow(V[i][1],2) +pow(V[i][2],2) );
    }
    return E;
    }
    
float PEnergy(float X[][3]){             //calculate potential Energy
    float r,deltax,deltay, deltaz,E=0;
    int j,i;
    for(i=0;i<Number;i++){
        for(j=0;j<Number;j++){
             if(j!=i){ // this IF is unnecessary use j=i instead
             	if(abs(X[i][0]-X[j][0])>LengthX/2){
            	  deltax=LengthX-abs(X[i][0]-X[j][0]);
				}
            	else{
                	deltax=abs(X[i][0]-X[j][0]);
				}


            	if(abs(X[i][1]-X[j][1])>LengthY/2){
            		deltay=LengthY-abs(X[i][1]-X[j][1]);
				}
              	else{
                	deltay=abs(X[i][1]-X[j][1]);
				}

				if(abs(X[i][2]-X[j][2])>LengthZ/2){
            		deltaz=LengthZ-abs(X[i][2]-X[j][2]);
				}
              	else{
                	deltaz=abs(X[i][2]-X[j][2]);
				}

        		r=sqrt(pow(deltax,2)+pow(deltay,2) +pow(deltaz,2)  );
                if(r<Cutoff){
                E+=4*(1/pow(r,12)-1/pow(r,6));
                }
		    }
		}
	}
    E=0.5*E;        // because of Uij=Uji
    
    return E;
    }

    
