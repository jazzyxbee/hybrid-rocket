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

