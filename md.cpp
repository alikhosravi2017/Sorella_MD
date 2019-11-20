#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <time.h>
#include "gasdev.cpp"
#include "other/nr.h"
#define  Length 8
#define  displacement  1.5
#define  Number 20    //number of atoms
#define  Cutoff 2.5   //cut off
#define  dump_every 1000   
#define  log_every  100  
#define  thermostat_every  100
#define  h          0.00001
#define  desired_Temperature          160 //K
#define  kB_True         1.38064852e-23  //m2 kg s-2 K-1
#define  epsilon_True    1.65e-21 //J
#define  sigma_True      3.4e−10    //m
/* MD units 
  kB         1 
  epsilon    1
  sigma      1
 
 * */
using namespace std;
void Force(float [][2],float [][2]);
float KEnergy(float [][2]);
float PEnergy(float [][2]);
void time_integration_verlet(float [][2],float [][2],float [][2], int);
void dump(float [][2],int, ofstream& );
void log(float[][2],float[][2],int, ofstream&);
void pbc(float[][2], int);
float Temperature(float [][2]);
void Thermostat_velocity_scaling(float [][2]);
int main(){
    
    float X[Number][2],V[Number][2],F[Number][2];
    float r,T,t,m,YMassVelocity=0,XMassVelocity=0;
    float X_cm=0, Y_cm=0;
    int i,j,k,N,a,step=0,idum=-8; 
    ofstream M_output_file ("dump.xyz");
    ofstream E_output_file ("log.txt");
    cout<<"inter full time:\n ";
//     cin>>T;
    T = 100000;
    k=sqrt(Number/2)+1;                       	
//     m=(Length/2)/(sqrt(Number/2)+1);  
    m= displacement;
    N=Number;
        for(j=0;j<=Number/k;j++){                    
            for(i=1;i<=k;i++){                      
                if (N!=0){
                X[Number-N][0]=i*m;
                X[Number-N][1]=(j+1)*m;
                X_cm += X[Number-N][0]/Number;
                Y_cm += X[Number-N][1]/Number;
                N=N-1;
                }
            }
    }
//     cout<<X_cm<<"\t"<<Y_cm<<endl;
//     // move center of mass to the middle:
    for(i=0;i<Number;i++){
        X[i][0] += (Length/2-X_cm);
        X[i][1] += (Length/2-Y_cm);
    }

    dump(X,step,M_output_file);
    //end of first place calculation...
    
    for(i=0;i<Number;i++){
        V[i][0]=NR::gasdev(idum);
        V[i][1]=NR::gasdev(idum);
    }
	for(i=0;i<Number;i++){                       //calculate Central Mass Velocity
		YMassVelocity+=V[i][1];
		XMassVelocity+=V[i][0];
	}
	for(i=0;i<Number;i++){                       //Central Mass Velocity=0
		V[i][1]=V[i][1]-YMassVelocity/Number;
		V[i][0]=V[i][0]-XMassVelocity/Number;
	}
	
	Thermostat_velocity_scaling(V);
	
	
	Force(X,F); // calculate a0
	
        /// MAIN ///
        clock_t tStart = clock();
        for(step=1;step<T;step+=1){
            time_integration_verlet(V,X,F, step);
            dump(X,step,M_output_file);
            log(X,V,step,E_output_file);
            if(step%thermostat_every==0)        Thermostat_velocity_scaling(V);
        }
        
    
        printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
        
return 0;
}

void pbc(float X[][2], int i){
    X[i][0] -= Length*floor(float(X[i][0]/Length));
    X[i][1] -= Length*floor(float(X[i][1]/Length));
    
//     if(X[i][0]>Length)  X[i][0]+=(-1)*Length;    //boundary conditions...
//     if(X[i][0]<0)       X[i][0]+=Length;		 	//boundary conditions...	
//     if(X[i][1]>Length)  X[i][1]+=(-1)*Length;    //boundary conditions...
//     if(X[i][1]<0)       X[i][1]+=Length;         //boundary conditions...
}


void Force(float X[][2],float F[][2]){             //calculate acceleration
    int i,j;
    float r2,deltax,deltay,solid_distanse_x,solid_distanse_y;
    float Length_2 = float(Length/2);
    
    for(i=0;i<Number;i++){
        F[i][0]=0;
        F[i][1]=0;
    }
	for(i=0;i<Number;i++){
	    for(j=0;j<Number;j++){
                if(j!=i){
                    
                    solid_distanse_x = abs(X[i][0]-X[j][0]);
                    deltax = -solid_distanse_x + Length*floor(float(solid_distanse_x/Length_2));   
                                
                    solid_distanse_y = abs(X[i][1]-X[j][1]);
                    deltay=  -solid_distanse_y + Length*floor(float(solid_distanse_y/Length_2));  

                    r2=pow(deltax,2)+pow(deltay,2);
                    
                    
                    if(sqrt(r2)<Cutoff){
                        if (X[i][0]<X[j][0])     F[i][0]+=(deltax*48*(1/pow(r2,7)-1/(2*pow(r2,4)))  );   //  note that  r2=r*r
                        else                     F[i][0]+=(-1*deltax*48*(1/pow(r2,7)-1/(2*pow(r2,4))));
                        if (X[i][1]<X[j][1])     F[i][1]+=(deltay*48*(1/pow(r2,7)-1/(2*pow(r2,4))));
                        else                     F[i][1]+=(-1*deltay*48*(1/pow(r2,7)-1/(2*pow(r2,4))));
                    }
                }
            }
    }
}



// void Force(float X[][2],float F[][2]){             //calculate acceleration
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

void Thermostat_velocity_scaling(float V[][2]){
    float lambda;
    int i;
    float temrature_True = float(epsilon_True)/kB_True;  // multiply temperature with this value to take in Kelvin

    float temperature_now =Temperature(V)*temrature_True;
//     cout<<"temperature_now="<<temperature_now<<endl;
    lambda = sqrt(desired_Temperature/temperature_now);
    for(i=0;i<Number;i++){
        V[i][0] *= lambda;
        V[i][1] *= lambda;
    }
    }

void time_integration_verlet(float V[][2], float X[][2], float F[][2], int step){
    	                                       
    // verlet loop //
    int i;
    for(i=0;i<Number;i++){                             //calculate v half step
        V[i][0]+=(h/2)*(F[i][0]);
        V[i][1]+=(h/2)*(F[i][1]);
        X[i][0]+=V[i][0]*h;
        X[i][1]+=V[i][1]*h;
        pbc(X,i);
    }
    Force(X,F);                                         //calculate a
    for(i=0;i<Number;i++){                               //calculate v half step
        V[i][0]+=(h/2)*(F[i][0]);
        V[i][1]+=(h/2)*(F[i][1]);
    }
        //cout<<step<<endl;
    }

void log(float X[][2],float V[][2],int step, ofstream& E_output_file){
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

void dump(float X[][2], int step, ofstream& M_output_file){
    if(step%dump_every==0){
        M_output_file<<Number<<endl<<"atoms"<<endl; 
        for(int i=0;i<Number;i++){
            M_output_file<<1<<"\t"<<X[i][0]<<"\t"<<X[i][1]<<"\t"<<1<<endl;
        }
//         cout<<"dumped at step="<<step<<endl;
    }
    }

float Temperature(float V[][2]){
    return KEnergy(V) * float(1)/(float(3)/2*Number);
    }

float KEnergy(float V[][2]){           //calculate kinetic Energy
    int i;
    float E=0;
    for(i=0;i<Number;i++){
    	E+=.5*(pow(V[i][0],2)+pow(V[i][1],2));
    }
    return E;
    }
    
float PEnergy(float X[][2]){             //calculate potential Energy
    float r,deltax,deltay, E=0;
    int j,i;
    for(i=0;i<Number;i++){
        for(j=0;j<Number;j++){
             if(j!=i){
             	if(abs(X[i][0]-X[j][0])>Length/2){
            	  deltax=Length-abs(X[i][0]-X[j][0]);
				}
            	else{
                	deltax=abs(X[i][0]-X[j][0]);
				}
            	if(abs(X[i][1]-X[j][1])>Length/2){
            		deltay=Length-abs(X[i][1]-X[j][1]);
				}
              	else{
                	deltay=abs(X[i][1]-X[j][1]);
				}
        		r=sqrt(pow(deltax,2)+pow(deltay,2)); 	
                if(r<Cutoff){
                E+=4*(1/pow(r,12)-1/pow(r,6));
                }
		    }
		}
	}
    E=0.5*E;        // because of Uij=Uji
    
    return E;
    }

    
