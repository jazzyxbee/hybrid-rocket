typedef struct{
	double Ds;
	double Vs;
	double ts;
	double Ms;
}Size;

typedef struct{
	double L;
	double D;
	double Din;
	double Lcyl;
	double mcase;
	double EtaV;
	double tcase;
}Combustion;

double Lf(double X, double K, double a, double Mox, double n, double rho_f,double OFinit)
{
double Lf;

//FILE *ERR;
//ERR=fopen("Error.dat","w");
//Mox=1;

//for (X=0.01;X<=20;X=X+0.01)
//{
Lf=pi*X*X*a/K*pow((4.*(K*K)/(X*X)*Mox/pi),n)*rho_f-Mox/OFinit;
//fprintf(ERR,"%lf \t %lf \n",X,Lf);
//Lcyl = fsolve( @(X)( (pi .* X.^2 .* a ./ K .* ( 4 .* K.^2 ./ X.^2 .* Mox ./ pi ).^n .* Furho ) - Mox ./ OFinit ) , 0.001 );
//}

return Lf;
}


double dLf(double X, double K, double a, double Mox, double n, double rho_f,double OFinit)
{
double dLf1;
double h=pow(10,-8.);
dLf1=(Lf(X+h,K,a,Mox,n,rho_f,OFinit)-Lf(X,K,a,Mox,n,rho_f,OFinit))/h;
return dLf1;
}


Combustion Case(double Pmass,double LDinratio,double OFinit,double Massflowox,double Pburst,double rhocase,double sigmamax,double a, double n,double rho_f)
{
	double Lcyl,Din,D,Vcyl,L,tcase,mcase,EtaV;
	//double V;
	double K=LDinratio;
	double Mox = Massflowox;

	double Vprop = Pmass / rho_f;
	double X,X0,Err;
	Err=1.;
	X0=1.5;
	int i=0;
		//---rownanie nieliniowe---
		printf("Solution of nonlinear equation for Grain Length estimation :");
		do
		{
		X=X0-Lf(X0,K,a,Mox,n,rho_f,OFinit)/dLf(X0,K,a,Mox,n,rho_f,OFinit);
		Err=fabs((X-X0)/X0);

		X0=X;
			if(X<=0){X0=10.;}
		i=i+1;
		printf("Err = %lf ,",Err);
		printf("X0 = %lf ,",X0);

		}while(Err>0.000001);
		printf("number of iteration %d\n",i);

	//Lcyl = fsolve( @(X)( (pi .* X.^2 .* a ./ K .* ( 4 .* K.^2 ./ X.^2 .* Mox ./ pi ).^n .* Furho ) - Mox ./ OFinit ) , 0.001 );
	Lcyl=X0;//temp
//---rownanie nieliniowe---
	Din = Lcyl/K;

	D=sqrt(pow(Din,2.)+4.*Vprop/pi/Lcyl);
//printf("\n $$$$$$$$$$ Din =%lf \t K= %lf \t D= %lf \t L=%lf \n ",Din,K,D,Lcyl);
	Vcyl = pow(D,2.)/4.*pi*Lcyl;
	L=D+Lcyl;
	tcase=Pburst*D/2./sigmamax;

		if ( tcase < 2*pow(10.,-3.) ){tcase = 2*pow(10.,-3.);}	
	
		mcase=rhocase*tcase*pow(D,2.)*pi*(1+Lcyl/D);
	
	//V=Vcyl+pi*pow(D,3)/6; 

	EtaV = Vprop/ Vcyl; //???

	Combustion A={L, D, Din, Lcyl, mcase,EtaV,tcase};
	//Caseprop = struct ( "L" , L , "D" , D , "Din" , Din , "Lcyl" , Lcyl , "V" , V , "t" , t , "EtaV" , EtaV , "mcase" , mcase );

	/*if ( iscomplex( Din ) )
	    fprintf('\nBreak in MASSCASE2 per Din complesso\n');
	    keyboard()
	endif

	if ( iscomplex( mcase ) )
	    fprintf('\nBreak in MASSCASE2 per mcase complesso\n');
	    keyboard()
	endif*/
	return A;

}

Size Sphere(double Pmass,double volload,double rhop,double Pburst,double rhocase,double sigmamax)
{
	double Ds,Vs,ts,Ms;

	Vs=Pmass/rhop/volload;
	Ds=pow(6.*Vs/pi,0.33333);
	ts=Pburst*Ds/4./sigmamax;//thickness of the sphere.
	//printf("\n\n\n ts=%lf \t Pburst=%lf \t Ds=%lf \t sigmax=%lf \t Vs=%lf \n\n\n",ts,Pburst,Ds,sigmamax,Vs);
	
		if (ts<1.5*pow(10.,-3.))
		{
		ts=1.5*pow(10.,-3.);
		}

	Ms=rhocase*ts*pow(Ds,2.)*pi;

	Size F = {Ds,Vs,ts,Ms};
	return F;

}

double nozzle (double mprop, double E)
{
double mnoz;
//mnoz=1e-3;
//printf("\n\n\n!!!!!! %lf !!!!!!!\n\n\n",E);
mnoz =125.0*pow((mprop/5400),(2./3.))*pow(E*0.1,0.25);
//printf("\n\n\n!!!!!! %lf !!!!!!!\n\n\n",mnoz);
return mnoz;
}

double injector (double D, double t)
{
double minj;
double rhoal=2800.;
minj = pow(D,2.)/4.*pi*rhoal;
return minj;
}
