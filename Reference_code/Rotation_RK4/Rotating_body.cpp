#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>

double func(double I1, double I2, double I3, double omega1,double omega2, double M1)
{
double omega_prim;
omega_prim=(M1-(I3-I2)*omega1*omega2)/I1;
return omega_prim;
}

int main()
{
std::ofstream printToFile;
printToFile.open ("output.dat");
double y0,x0,v0,w0,x,v,w,y;
int i;
double y1;
double t0,tf,N,h;
double K11,K12,K13;
double K21,K22,K23;
double K31,K32,K33;
double K41,K42,K43;
double I1=3500;
double I2=1000;
double I3=4200;
double M1=-1.2;
double M2=1.5;
double M3=13.5;
tf=100.;
t0=0.;
x0=0.0;
y0=0.1;
v0=-0.2;
w0=0.33;
N=10000;
h=(tf-t0)/N;
x=x0;
v=v0;
w=w0;
y=y0;
for (i=0;i<N;i++)
{
K11=h*func(I1,I2,I3,v,w,M1);//y
K12=h*func(I2,I3,I1,y,w,M2);//v
K13=h*func(I3,I1,I2,y,v,M3);//w

K21=h*func(I1,I2,I3,v+0.5*K12,w+0.5*K13,M1);
K22=h*func(I2,I3,I1,y+0.5*K11,w+0.5*K13,M2);
K23=h*func(I3,I1,I2,y+0.5*K11,v+0.5*K12,M3);

K31=h*func(I1,I2,I3,v+0.5*K22,w+0.5*K23,M1);
K32=h*func(I2,I3,I1,y+0.5*K21,w+0.5*K23,M2);
K33=h*func(I3,I1,I2,y+0.5*K21,v+0.5*K22,M3);

K41=h*func(I1,I2,I3,v+0.5*K32,w+0.5*K33,M1);
K42=h*func(I2,I3,I1,y+0.5*K31,w+0.5*K33,M2);
K43=h*func(I3,I1,I2,y+0.5*K31,v+0.5*K32,M3);

y=y+1./6.*(K11+2.*K21+2.*K31+K41);
v=v+1./6.*(K12+2.*K22+2.*K32+K42);
w=w+1./6.*(K13+2.*K23+2.*K33+K43);
//std::cout << x <<"\t"<< y <<"\t"<< v <<"\t"<< w <<"\t"<< h <<"\n";
printToFile << x <<"\t"<< y <<"\t"<< v <<"\t"<< w <<"\t"<< h <<"\n";
x=x+h;
}



//printToFile.close();
return 0;
}


