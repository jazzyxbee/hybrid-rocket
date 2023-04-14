	
	double Thrust;
	double At,Ae;
	//Nozzle geometry
	Thrust=Is*propflow*ro_c;
	At=Thrust/(Cf1*Pc);
	Ae=At*E;	
fprintf (fr,"OF ratio =%lf \nLD ratio = %lf \nThrust= %lf [kN]\nA_t = %lf [m^2]\nA_e = %lf [m^2]\nM_inert/M_prop = %lf \nL = %lf [m]\nD_out = %lf [m]\nthick_case=%lf [mm] \nD_ox=%lf [m]\nIs = %lf [s] \n",OF3,LDratio,Thrust/1e3,At,Ae,Minerts/Mprop,L2,(D2+2*t1),t1*1e3,Dox,Is);

fprintf (CAD,"$fn = 150%c \n difference(){ cylinder(r=%lf, h=%lf, center=true)%c cylinder(r=%lf, h=%lf, center=true)%c }\n",59,(D2+2*t1)/2,Lcyl,59,D2/2,Lcyl*1.1,59);
fprintf (CAD,"difference(){ translate([0,0,%lf])sphere(r=%lf)%c cylinder(r=%lf, h=%lf, center=true)%c translate([0,0,%lf])sphere(r=%lf)%c }\n",Lcyl/2,(D2+2*t1)/2,59,(D2+2*t1)/2,Lcyl,59,Lcyl/2,D2/2,59);

fprintf (CAD,"translate([0,0,%lf])sphere(r=%lf)%c",Lcyl/2+(D2+2*t1)/2+(Dox/2),Dox/2,59);

