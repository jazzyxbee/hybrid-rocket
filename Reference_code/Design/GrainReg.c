typedef struct{
	double Mf;
	double Rf;
}Balistyka;

Balistyka Paliwo(double Massflowox,double Lcyl ,double Diam,double a,double n,double Furho)
{
double Asec,Gox,Rf,Mf;

Asec = pow(Diam,2.)/4.*pi;
//printf ("Asec = %lf \n ",Asec);
Gox = Massflowox/Asec;
//printf ("Gox = %lf \n",Gox);
Rf = a*pow(Gox,n);
//printf ("Rf = %lf \n ",Rf);

Mf = pi*Diam*Lcyl*Rf*Furho;

//printf ("Mf = %lf ",Mf);

Balistyka M={Mf,Rf};

return M;
}
