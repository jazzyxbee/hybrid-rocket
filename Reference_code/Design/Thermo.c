typedef struct{
	double Vtank;
	double Vpress;
	double Mpress;
}Press;

double Pe(double E, double gm1, double pr, double g, double gp1, double gm1g)
{
double ff1;
ff1=1.0/E-pow((gp1/2.0),(1.0/gm1))*pow(pr,(1.0/g))*sqrt((gp1/gm1)*(1.-pow(pr,gm1g)));
return ff1;
}

double dPe(double E, double gm1, double pr, double g, double gp1, double gm1g)
{
double dPe1;
double h=pow(10,-8.);
dPe1=(Pe(E,gm1,pr+h,g,gp1,gm1g)-Pe(E,gm1,pr,g,gp1,gm1g))/h;
return dPe1;
}

double pratio (double kappa,double E)
{
double gm1=kappa-1.0;
double gp1=kappa+1.0;
double gm1g=(kappa-1.0)/kappa;
double g=kappa;

//Metoda bisekcji
//double a1=1e-5;
//double a3=0.5;
//double a2=(a1+a3)/2.0;
double pr0=0.001;
double Err=0.1;
double pr;
//double delta,ff1,ff2;
//delta=100.;

/*printf("Solution of nonlinear equation for Pressure ratio using bisection:");
while( delta > 0.0001)
{

  pr=a1;
  ff1=1.0/E-pow((gp1/2.0),(1.0/gm1))*pow(pr,(1.0/g))*sqrt((gp1/gm1)*(1.-pow(pr,gm1g)));
  
  pr=a2;
  ff2=1.0/E-pow((gp1/2.0),(1.0/gm1))*pow(pr,(1.0/g))*sqrt((gp1/gm1)*(1.-pow(pr,gm1g)));

  	if ((ff1*ff2)>0.0)
	{
    	a1 = a2;
	}
 	 else
	{
    	a3 = a2;
	}
  
  a2=(a1+a3)/2.0;
  delta=(a3-a1)/a2;

}
pr=a2;
 */
int i=0;

		printf("Solution of nonlinear equation for Pressure ratio :");
		do
		{
		pr=pr0-Pe(E,gm1,pr0,g,gp1,gm1g)/dPe(E,gm1,pr0,g,gp1,gm1g);
		Err=fabs((pr-pr0)/pr0);
		pr0=pr;
		i=i+1;
		//printf("Pratio0 = %lf ,",pr0);
		}while(Err>0.0001);
		printf("Solution = %lf , after ",pr0);
		printf("%d iterations \n",i);

printf("%lf \n",pr);
return pr;
}

Press Design(double  Ppress,double Ptank,double Vliq ,double Filling )
{
	if ( Ppress <= Ptank )
	{    	
	printf("\n Error Ppress <= Ptank\n");
	system("pause");
    	}

	double gamma = 1.4;
 	double Mmol = 0.028;
 	double Tref = 300.;
 	double Vtank;
	double Texp,Rhop,Mpress,Vpress;
	Vtank=Vliq/Filling;
	//printf("V tank =%lf \n",Vtank);
 	Texp=Tref*pow((Ptank/Ppress),((gamma-1)/gamma));
	//printf("T exp =%lf \n",Texp);//Temperatura po rozprezeniu w zbiorniku
	Vpress=(Vtank*Ptank/Texp)/(Ppress/Tref-Ptank/Texp);  
	//Vpress=Vtank/(pow((Ptank/Ppress),(-1/gamma)) - 1);
	//printf("V press =%lf \n",Vpress);
	Rhop=Ppress/R2*Mmol/Tref;
	//printf("Rho press =%lf \n",Rhop);
 	Mpress=Vpress*Rhop;
	//printf("Mass press =%lf \n",Mpress);
	//printf("--------------------\n");
 	Press D = {Vtank, Vpress, Mpress};
	
	if ( Mpress <= 0 )
	{    	
	printf("\n Mpress <= 0\n");
	system("pause");
    	}

	return D;
}

