//
// Created by User on 09/11/2023.
//

#ifndef UNTITLED4_HYBRIDROCKETUPDATE9NOV_H
#define UNTITLED4_HYBRIDROCKETUPDATE9NOV_H

#include <iostream>
#include <stdio.h>
#include <cmath>
#include <iterator>
#include <vector>

double Pe(double E, double gm1, double pr, double g, double gp1, double gm1g)
{
    double ff1; ///why ff1?
    ff1=1.0/E-pow((gp1/2.0),(1.0/gm1))*pow(pr,(1.0/g))*sqrt((gp1/gm1)*(1.-pow(pr,gm1g)));
    return ff1;
}

double dPe(double E, double gm1, double pr, double g, double gp1, double gm1g)
{
    /// how does this work ? how does dividing by 10^-8 result in the de
    double dPe1;
    double h=pow(10,-8.);
    dPe1=(Pe(E,gm1,pr+h,g,gp1,gm1g)-Pe(E,gm1,pr,g,gp1,gm1g))/h;
    return dPe1;
}

double pratio (double kappa,double E) {
    //using newton raphson and (eq 3-25) from propulsion elements
    double gm1 = kappa - 1.0; /// why gm1 and gm1g?
    double gp1 = kappa + 1.0;
    double gm1g = (kappa - 1.0) / kappa;
    double g = kappa;
    double pr0=0.001;
    double Err=0.1;
    double pr;

    int i=0;

    //printf("Solution of nonlinear equation for Pressure ratio :");
    do
    {
        pr=pr0-Pe(E,gm1,pr0,g,gp1,gm1g)/dPe(E,gm1,pr0,g,gp1,gm1g);
        Err=fabs((pr-pr0)/pr0);
        pr0=pr;
        i=i+1;
        //printf("Pratio0 = %lf ,",pr0);
    }while(Err>0.0001);
    //printf("Solution = %lf , after ",pr0);
    //printf("%d iterations \n",i);

    //printf("%lf \n",pr);
    return pr;
}

double hybrid_rocket_thrust(double RI, double L, double T, double MOL, double kappa, double Pc, int E, double mDot_o) {
    // # kappa is specific heat ratios
    // # T is chamber temperature as a result of chemical reaction
    // # e is expansion ratio
    // # mDot_o is oxidizer flow rate
    // # equations from Rocket Propulsion Elements unless otherwise stated

    // -------------------------------Calculating Mdot propellant mass flow rate------------------------//
    // Mdot variables
    double a = 0.488;  // # regression rate coefficient of parrafin mm/s from Hybrid Rocket Fuel Regression Rate Data and Modeling
    double n = 0.62; // # regression rate pressure from Hybrid Rocket Fuel Regression Rate Data and Modeling
    int Pf = 900; // # fuel density of parrafin kg/m^3
    double A_b      = 0; // #  Combustion port surface area: L = port length  RI = internal radius of fuel grain (16-11)
    double G_o      = 0; // # Oxidizer mass velocity (16-2) gm/cm2 hence change RI to centimeters and mDot from kg to gm
    double rDot     = 0; // # Regression rate (16-10)
    double mDot_f   = 0; // # Mass Flow Rate of Fuel (16 - 11) kg/sec
    double mDot     = 0; // # Mass Flow Rate (6-2)

    A_b = 2 * M_PI * RI * L; // meters squared

    G_o = (mDot_o*1000) / (M_PI * pow((RI*100), 2));  //  gm/cm s
    rDot = a * pow(G_o, n); // millimeters per second
    mDot_f = (Pf * A_b * rDot/1000);/// kg/m^3 * m^2 * m/s
    mDot = mDot_o + mDot_f;

    // -------------------------------Calculating Ve Exit velocity------------------------------------//

    double R = 8314.41; // # Gas constant J / kmol * K
    double Pe = 0;  // # exit pressure from hybrid-rocket/Reference_code/program_2.c
    double v_e = 0; // # Calculate the exit velocity
    double Thrust = 0; // # Calculate the rocket thrust from Space-Propulsion-Analysis-and-Design-McGraw-Hill-(1995)

    Pe = Pc*pratio(kappa,E);
    v_e = sqrt((2 * kappa * R * T) / ((kappa - 1)*MOL) * (1.0 - pow((Pe / Pc), ((kappa - 1) / (kappa)))));
    Thrust = mDot * v_e;

    double Ip = v_e/9.81;

    /* variables can be printed here
    printf( "velocity %f\n",v_e) ;
    printf("Specific impulse V_e: %f \n", Ip);
    */

    return Thrust;
}
#endif //UNTITLED4_HYBRIDROCKETUPDATE9NOV_H
