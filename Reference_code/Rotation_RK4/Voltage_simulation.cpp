#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>


void f(x)
{
double y=x;
return 0;
}

int main()
{
std::ofstream printToFile;
printToFile.open ("output.dat");
double y0,y1;
double x0=10.;
f(x0);
std::cout << "funkcja"<<f(x0) ;
//printToFile.close();
return 0;
}


