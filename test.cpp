// stof example
//#include <iostream>   // cout
//#include <string>     // string, stof
//#include <sstream>
//#include <algorithm>
//#include <iterator>
//#include <vector>
#include <stdio.h>
#include <string.h>
//#include <iostream>


#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <cmath>

using namespace std;

int main ()
{
  string orbits ("686.97 4444.24   365.24 ");
  string::size_type sz;     // alias of size_t

  float mars = stof (orbits,&sz);
//  cout<<orbits.substr(sz);
  float earth = stof (orbits.substr(sz),&sz);
  float earth2 = stof (orbits.substr(sz).substr(sz),&sz);
  cout << "Fist= " << (mars) << "\t Second="<<earth<<"\t  Third= "<< earth2<<"\n";
//  cout << "Fist= " << stof("23.5") <<"\n";
//  cout<< strtok(orbits," ");
////  sss=
////  split2(orbits,sss);

//    string line="this ";
//    string line2="is a test";
//    int Number =34;
//    cout<<line+line2<<endl<<"initial_positions_"+to_string(Number)+"_atoms.xyz";
//    string::size_type sz;
//    char * pch;
//    ifstream myfile ("initial_positions_purenumbers.xyz");
////    if (myfile.is_open()){
////        while ( getline (myfile,line) ){
////
////            cout<<line<<endl;
////            pch = strtok (line," ,.-");
////            cout<<pch;
////        }
////        myfile.close();
////    }
////    else cout << "Unable to open file";
////
//
//  double a;
//  myfile >> a;
//  cout<<a<<endl;
//  myfile >> a;
//  cout<<a<<endl;
//  myfile >> a;
//  cout<<a<<endl;
//  myfile >> a;
//  cout<<a<<endl;
//  myfile >> a;
//  cout<<a<<endl;
//  myfile >> a;
//  cout<<a<<endl;
//  myfile >> a;
//  cout<<a<<endl;
//  myfile >> a;
//  cout<<a<<endl;

//  string aaa= "1213.4 34.0";
//  aaa >> a;
//  cout<<a<<endl;
//  aaa >> a;
//  cout<<a<<endl;
return 0;
}