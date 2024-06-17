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

struct dataPoint {
    double temperature; // oxygen fuel temperature, temperature inside chamber during reaction in Kelvin
    double gamma;       // kappa ratio of specific heats (constant pressure vs constant volume)
    double molecular_mass;  // mols used in reaction, molar mass (g/mol)
    double expansionRatio;
};

// Function to perform linear extrapolation
double linear_extrapolate(double x1, double y1, double x2, double y2, double x) {
    return y1 + ((y2 - y1) / (x2 - x1)) * (x - x1);
}

// Function to perform linear extrapolation based on two data points
dataPoint extrapolate_data_point(const dataPoint& dp1, const dataPoint& dp2, double expansionRatio) {
    double temperature_extrap = linear_extrapolate(dp1.expansionRatio, dp1.temperature,
                                                   dp2.expansionRatio, dp2.temperature, expansionRatio);
    double gamma_extrap = linear_extrapolate(dp1.expansionRatio, dp1.gamma,
                                             dp2.expansionRatio, dp2.gamma, expansionRatio);
    double molecular_mass_extrap = linear_extrapolate(dp1.expansionRatio, dp1.molecular_mass,
                                                      dp2.expansionRatio, dp2.molecular_mass, expansionRatio);
    return {temperature_extrap, gamma_extrap, molecular_mass_extrap, expansionRatio};
}


// Function for linear interpolation
double linear_interpolate(double x0, double y0, double x1, double y1, double x) {

    // from https://www.cuemath.com/linear-interpolation-formula/
    if (x0 == x1) {
        std::cerr << "Error: x0 and x1 cannot be the same value.\n";
        return -1;
    }
    return y0 + ((x - x0) * (y1 - y0) / (x1 - x0));
}

// Function to find the appropriate interval and perform interpolation to predict temperature, gamma, and molecular mass based on expansion ratio
dataPoint interpolate_data_point(const std::vector<dataPoint>& data, double expansionRatio) {
    // Extrapolate if values are out of range
    if (expansionRatio <= data.front().expansionRatio) {
        // Extrapolate based on the first two data points
        return extrapolate_data_point(data[0], data[1], expansionRatio);
    } else if (expansionRatio >= data.back().expansionRatio) {
        // Extrapolate based on the last two data points
        return extrapolate_data_point(data[data.size() - 2], data[data.size() - 1], expansionRatio);
    } else {
    }
    //Interpolate if values are within range
    for (size_t i = 0; i < data.size() - 1; ++i) {
        if ((data[i].expansionRatio <= expansionRatio && expansionRatio <= data[i + 1].expansionRatio) ||
            (data[i].expansionRatio >= expansionRatio && expansionRatio >= data[i + 1].expansionRatio)) {

            double temperature_interp = linear_interpolate(data[i].expansionRatio, data[i].temperature,
                                                           data[i + 1].expansionRatio, data[i + 1].temperature, expansionRatio);
            double gamma_interp = linear_interpolate(data[i].expansionRatio, data[i].gamma,
                                                     data[i + 1].expansionRatio, data[i + 1].gamma, expansionRatio);
            double molecular_mass_interp = linear_interpolate(data[i].expansionRatio, data[i].molecular_mass,
                                                              data[i + 1].expansionRatio, data[i + 1].molecular_mass, expansionRatio);
            return {temperature_interp, gamma_interp, molecular_mass_interp, expansionRatio};
        }
    }
    std::cerr << "Error: Expansion ratio value out of range for interpolation.\n";
    return {-1, -1, -1, -1}; // Return invalid dataPoint if out of range
}

double Pe(double E, double gm1, double pr, double g, double gp1, double gm1g)
{
    double ff1;
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
    double gm1 = kappa - 1.0;
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

double Mass_Flow_Rate(double RI, double L, double mDot_o){
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

    return mDot;
}

double exit_velocity(double E, double Pc,double kappa, double MOL, double T){
    // -------------------------------Calculating Ve Exit velocity------------------------------------//
    // # equations from Rocket Propulsion Elements unless otherwise stated
    double R = 8314.41; // # Gas constant J / kmol * K
    double Pe = 0;  // # exit pressure from hybrid-rocket/Reference_code/program_2.c
    double v_e = 0; // # Calculate the exit velocity

    Pe = Pc*pratio(kappa,E);
    v_e = sqrt((2 * kappa * R * T) / ((kappa - 1)*MOL) * (1.0 - pow((Pe / Pc), ((kappa - 1) / (kappa)))));

    return v_e;
}

double hybrid_rocket_thrust(const std::vector<dataPoint>& dataPoints, double gamma, double MOL, double RI, double L, double Pc, int E, double mDot_o) {

    double thrust = 0; // # Calculate the rocket thrust from Space-Propulsion-Analysis-and-Design-McGraw-Hill-(1995)
    double v_e = 0;

    double mDot = Mass_Flow_Rate(RI, L, mDot_o);

    // Interpolate a new data point based on expansion ratio
    dataPoint interpData = interpolate_data_point(dataPoints, E);
    if (interpData.temperature == -1) {
        std::cerr << "Error: Interpolation failed. Exiting thrust calculation.\n";
        return -1;
    }

    v_e = exit_velocity(E, Pc, interpData.gamma, interpData.molecular_mass,interpData.temperature);

    // Calculate thrust using the calculated exit velocity
    // Note: You may need to modify this calculation based on your requirements
    // Example calculation: Thrust = mDot * v_e
    thrust = mDot * v_e;

    printf("Interpolated Data - Temperature: %f, Gamma: %f, Molecular Mass: %f, Expansion Ratio: %f, Thrust: %f\n",
           interpData.temperature, interpData.gamma, interpData.molecular_mass, interpData.expansionRatio,thrust);


    double Ip = v_e/9.81;

    /* variables can be printed here
    printf( "velocity %f\n",v_e) ;
    printf("Specific impulse V_e: %f \n", Ip);
    */

    return thrust;
}
#endif //UNTITLED4_HYBRIDROCKETUPDATE9NOV_H
