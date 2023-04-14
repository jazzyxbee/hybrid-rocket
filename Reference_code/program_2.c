#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>


#define NH 2
#define g0 9.81
#define R2 8.3144621
#define pi 3.141592

#include "Design/Mechanical.c"
#include "Design/Mission.c"
#include "Design/Thermo.c"
#include "Design/GrainReg.c"
#include "Aux/Interpolation.c"


//-------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
// -----------------------MAIN------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------
int main(void)
{

FILE *fw;
fw=fopen("wyniki.dat","w");
FILE *fp;
fp=fopen("graph3d.dat","w");
FILE *fr;
fr=fopen("DESIGN.dat","w");
FILE *ERR;
ERR=fopen("Error.dat","w");
FILE *CAD;
CAD=fopen("OpenCAD2","w");

double Tc;//K
double Cp,Cv,OF3;// specific heat for T 
double M;
double tcas,t1;//case thickness
//constants
double Mprop,L2,D2,Dox;//gas constants in different unit
double RM,Is;
double kappa,kappa1;
double DeltaOF;
double Mnoz,Minj,Minerts,Mtot,OFinit;
double Pt,Tt,ro_t,ut;
double ro_c;//kg/m^3
double Pe,Tal,cstar;
double Iv,Cf1;
//Mass Margin for inert and propeller
double margininerts;
double marginprop;
double Lcyl,Din,Dfin,Rin,Rfin;
double OFmean,D0,OF1,OF2,t,dt,OFmeanradius,Prm;

double	niu=1.0;//straty Is
double LDratio=5.;//Intial
double E=100;
double OF;//Initial
double rho_ox=1390;
//double rho_ox=1390;
double rho_f=920;
double Pc=1.0*pow(10.,6.); //atm
double Pa=0;
//double time=85;//sec
double propflowj,propflow;
double Moxx,Mcc;

for(OF3=4;OF3<=14;OF3=OF3+0.1)
{
OFinit=OF3;

for(LDratio=4;LDratio<=14;LDratio=LDratio+0.1)
{
OF=9;

printf("\n--------------------------------------------------- ITERACJA PETLI LD ----------------------- \n");
DeltaOF=100.;
	while (DeltaOF > 0.01)
	{
	printf("\n--------------------------------------------------- ITERACJA PETLI OF Delta OF = %lf ----------------------- \n",DeltaOF);

		//printf("OF przed gas= %lf \n", OF);
		Tc=Temp_c(OF);
		M=Mmol(OF)/1.e3;
		RM=R2/M;//(J/(K*mol)
		Cp=CP(OF);

		Cv=Cp-R2/(M*1e3);//meyers equation

		//printf ("\n M= %lf \t Cv = %lf \t Cp = %lf \t Tc = %lf \n",M,Cv,Cp,Tc);
		kappa1=Cp/Cv;
		kappa=Gamma(OF);
		printf("\t Kappa = %lf \t %lf \n", kappa,kappa1);
		//kappa=1.2962;

		fprintf (fw,"---------------Specific heat ratio---------------");
		fprintf (fw,"\n Kappa= %lf \t R/M = %lf \n",kappa,RM);
		//---------------rocket performance------------------//

		ro_c=rho_ox*rho_f*(1+OF)/(rho_f*OF+rho_ox);//g/cm3 - mieszanina paliwa z utleniaczem

		Pt=Pc*pow((2/(kappa+1)),(kappa/(kappa-1)));
		Tt=Tc*(2/(kappa+1));
		ro_t=ro_c*pow((2/(kappa+1)),(kappa/(kappa-1)));
		ut=sqrt(kappa*RM*Tt);

		fprintf (fw,"---------------Rocket Performance---------------");
		fprintf (fw,"\n Throat temperature %lf\t  Throat pressure %lf \n",Tt,Pt);
		fprintf (fw,"\n throat velocity %lf \t Throat density %lf\t density in CC %lf\n",ut,ro_t,ro_c);

	//-------------------------------------------------------SPECIMPULSE2------------------------------------------
			fprintf (fw,"---------------Specific Impulses---------------");

			Pe=Pc*pratio(kappa,E);
	
			//printf ("\n exit pressure [Pa] %lf \n",Pe);
			fprintf (fw,"\n exit pressure [Pa] %lf \n",Pe);
		
			Tal=sqrt(kappa*pow((2/(kappa+1)),((kappa+1)/(kappa-1))));// Vandenkerckhove function 
			cstar=1/Tal*sqrt(RM*Tc); // Characteristic velocity

			Cf1=Tal*sqrt(((2*kappa)/(kappa-1))*(1-pow(1/(Pc/Pe),(kappa-1)/kappa)))+E*(Pe-Pa)/Pc;
			Is=niu*Cf1*cstar/g0; //specific impulse
			Iv=Is*ro_c;
			fprintf (fw,"\n characteristic velocity %lf \n",cstar);
			fprintf (fw,"\n Specific Impulse %lf s \n Volumetric Imuplse %lf \n",Is,Iv);
			fprintf (fw,"\n thrust coefficient %lf \n",Cf1);
	//-------------------------------------------------------Mass Fun------------------------------------------
	//----------------------------------------------Obliczenia Wytrzymalosciowe---------------------
	double Ppress=15*pow(10.,6.);//Pressure of Pessurizers
	double Massflowox=1.00;//kg/s
	double X = 0.25;//STRUCT/PROPELANT
	double X2;//Nowy X - liczony
	double Ptank=Pc/0.4; 
	double Poper=Ptank;
	double Ullage=0.05; 
	double safety=2.0;
	double DeltaI=100.;//interations

	double DeltaV=300.00;//m/s
	double Mpay=1000;//[kg]
	double Filling=1.-Ullage;


	//Dane balistyczne
	double Marxa = 2.953e-5;//9.3684e-5;// 0.198e-5;// Marxman
	double Marxn = 0.52;//0.325;//0.325; // Marxman

	DeltaI=100.;
		while(DeltaI>0.01)
		{



		Mass C;
		C = Budget(DeltaV,Is,X,Mpay,OF,rho_ox,rho_f); 
		fprintf(fw," Mass propelant  %lf \t Mass Oxydizer %lf \n Mass Fuel %lf \t Mass Structure %lf \n Volume Oxydizer %lf \t Volume fuel %lf \n \n \n", C.Mprop, C.Mox, C.Mfu, C.Ms, C.Vox, C.Vfu);

		//parametry Pressurizers
		Press D;
		D = Design(Ppress,Ptank,C.Vox,Filling); 
		fprintf(fw," Presurizer tank diameter  %lf \t Amount of presurizer %lf \t Mass of presurizer %lf \n", D.Vtank, D.Vpress, D.Mpress);


		Size Ox;
		Ox = Sphere(C.Mox,Filling,rho_ox,Poper*safety,4430,724e6); 
		//Ox = Sphere(C.Mox,Filling,rho_ox,Poper*safety,1500,263e6); 
		fprintf(fw," Oxydizer tank diameter  %lf \t Oxydizer tank Volume %lf \t Oxydizer tank Mass %lf \n", Ox.Ds, Ox.Vs, Ox.Ms);
		Dox=Ox.Ds;

		Size Pr;
		Pr = Sphere(D.Mpress,Filling,Ppress/R2*0.028/300,Ppress*safety,4430,724e6); 
		fprintf(fw," Presurizer tank diameter  %lf \t Presurizer tank Volume %lf \t Presurizer tank Mass %lf \n", Pr.Ds, Pr.Vs, Pr.Ms);

		
		Combustion CC;
		CC = Case(C.Mfu,LDratio,OFinit,Massflowox,Pc*safety,2800,214e6,Marxa,Marxn,rho_f);//ALUMINUM
		//CC = Case(C.Mfu,LDratio,OFinit,Massflowox,Pc*safety,1600,223e6,Marxa,Marxn,rho_f);//carbon-epoxy
		fprintf(fw," Combustion Chamber Length  %lf \t Combstion Chamber diameter %lf \t Internal Diameter %lf \t CC Mass %lf \t Propelant volume ratio %lf \n", CC.L, CC.D, CC.Din, CC.mcase, CC.EtaV);

		tcas=CC.tcase;
		Lcyl=CC.Lcyl;
		Din=CC.Din;
		Dfin=CC.D;
		Mnoz=nozzle(C.Mprop,E);
		fprintf(fw," Nozzle mass  %lf \n", Mnoz);

		Minj=injector(CC.D/3,0.02);
		fprintf(fw," Injector mass  %lf \n", Minj);

		//Liczy mase calkowita
		  margininerts = 1.1;
		  marginprop = 1.05;
		
		Moxx=Ox.Ms;Mcc=CC.mcase;
		  Minerts=margininerts*(D.Mpress+Ox.Ms+Pr.Ms+CC.mcase+Mnoz+Minj);

		  //Minerts=margininerts*(Ox.Ms+CC.mcase+Mnoz+Minj);
		  Mprop=C.Mprop*marginprop;
		  Mtot=Mprop+Minerts;

		  fprintf(fw," Mass Inert %lf \n", Minerts);
		  fprintf(fw," Total Mass %lf \n", Mtot);

		  printf(" Inert/Propelant Old =  %lf ,", X);
		  X2=Minerts/Mprop;
		  printf(" Inert/Propelant New =  %lf ,", X2);
		  DeltaI=fabs((X2-X)/((X2+X)/2));
		  printf(" Inert/Propelant Difference =  %lf \n", DeltaI);
		  X=X2;

		printf("OFinitial=%lf \t LDratio=%lf\n",OFinit,LDratio);



		}


	int j;
	Rin=Din/2.;Rfin=Dfin/2.;
	D0=Din;
	Balistyka Pr1,Pr2;
	Pr1=Paliwo(Massflowox,Lcyl,Din,Marxa,Marxn,rho_f);
	Pr2=Paliwo(Massflowox,Lcyl,Dfin,Marxa,Marxn,rho_f);
	Prm=(Pr1.Rf+Pr2.Rf)/2.;
	dt=(Rfin-Rin)/100./Prm;
	t=0;
	OF1=0;
	j=0;
	propflowj=0;

			do
			{
			   Balistyka P;
			   P=Paliwo(Massflowox,Lcyl,D0,Marxa,Marxn,rho_f);
			   OF1=Massflowox/P.Mf+OF1;
			   D0=D0+2*P.Rf*dt;
			   t=t+dt;
			   j=j+1;
			   propflowj=P.Mf+Massflowox+propflowj;
			}
			while ( D0 < Dfin );
	
	propflow=propflowj/j;

	D0=Din;
	OFmean=OF1/j;	
	OFmeanradius=pow(Massflowox,(1-Marxn))*pow(pi,(Marxn-1))*(pow(Rfin,(2*Marxn))-pow(Rin,(2*Marxn)))/(4*Marxn*Marxn*Lcyl*rho_f*Marxa*(Rfin - Rin));
	OF2=(OF+OFmean)/2.;
	DeltaOF=fabs(OF-OFmean)/OF2;
	printf(" OFmean =  %lf , Old OF = %lf , Delta OF = %lf \n", OFmean,OF,DeltaOF);
	printf(" OFmeanradius =  %lf \n", OFmeanradius);
	OF=OF2;
	t1=tcas;
	L2=Lcyl;
	D2=Dfin;

	}
//printf("\n\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%lf  %lf!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n",OF3,OFinit);
//zapis do pliku OPENCAD
if(OF3<11.7 && OF3>11.6 && LDratio<4.1)
	{
	#include "Design/CAD.c"

	}

//D.Mpress+Ox.Ms+CC.mcase+Mnoz+Minj
fprintf(ERR,"%lf \t %lf \t %lf %lf \n",Moxx,Mcc,Mnoz,Minj);
//fprintf (fp,"%lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf  \n",OF3,L2/D2,Mtot,Minerts,Mprop,Minerts/Mprop,L2,D2,Is);
fprintf (fp,"%lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf  \n",OF3,LDratio,Mtot,Minerts,Mprop,Minerts/Mprop,L2,D2,Is);

}
fprintf (fp,"\n");
}

return 0;

}
