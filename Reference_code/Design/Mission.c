typedef struct{
	double Mprop;
	double Mox;
	double Mfu;
	double Ms;
	double Vox;
	double Vfu;
}Mass;

Mass Budget(double DeltaV, double Is, double X, double Mpay,double OF,double rho_ox, double rho_f)
{
	double K,Mprop,Mox,Mfu,Ms;
	double Vox,Vfu;
	K=exp(DeltaV/Is/g0);

	//Mprop=Mpay*(K-1.)/(1.+X-K*X);
	Mprop=Mpay*(1.-K)/(K*X-X-1.);
	//printf("MPROP--------- %lf \n",DeltaV/(Is*g0));
	//printf("MPROP--------- %lf \t %lf \t %lf \t %lf \n",Mprop,K,X,Is);
	Mox=Mprop*OF/(OF + 1.);
	//printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$Mox %lf",Mox);
	Mfu=Mprop-Mox;

	//Ms=Mprop/(1.-X)+Mpay; ORGINAL
	Ms=Mprop*X;//Moje!

	Vox=Mox/rho_ox;
	Vfu=Mfu/rho_f;
		if ( Mfu < 0.0 )
			{
		    	printf("\n Error : Masa paliwa ujemna\n");
			
			}
		if ( Mox < 0.0 )
			{
		    	printf("\n Error : Masa utleniacza ujemna \n");
			
			}
	Mass C = {Mprop, Mox, Mfu, Ms, Vox, Vfu};
	return C;
}

