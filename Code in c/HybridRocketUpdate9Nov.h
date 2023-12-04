//
// Created by User on 09/11/2023.
//

#ifndef UNTITLED4_HYBRIDROCKETUPDATE9NOV_H
#define UNTITLED4_HYBRIDROCKETUPDATE9NOV_H
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <iterator>
#include <vector>
#include <complex>

double pi = 3.14159265359;



double fuel_mass(double r, double h) {
    double paraffin_d = 900; // kg / m3\n
    return (pi * pow(r,2) * h) * paraffin_d; // pi * r ^ 2 * h * paraffin density     // should be changed to
}

double fuel_effective_grain_area(double r,double h){
    return 2*pi * r * h;//A = pi * r ^ 2     // isnt right 2Pi R1 H
}

//Solving for Pe using the Newton-Raphson method
// PROBLEM IS IN THIS FUNCTION RESULTING IN NAN VALUES AT e = 2 problem is pr goes negative and pow cannot go to a negative number
std::complex<double> Pe(double E, double gm1, std::complex <double> pr, double g, double gp1, double gm1g){
    // Calculate the exit pressure ratio from the given formula
   //return 1.0 / E - pow((gp1 / 2.0),(1.0 / gm1)) * pow(pr,(1.0 / g)) * sqrt(complex<double>(gp1 / gm1) * (1.0 - (pow(pr,gm1g))));
    // Check if pr is non-negative before performing calculations
   /* if (pr < 0) {
        // Handle the case when pr is negative to avoid complex results
        return std::nan("Invalid input");
    }*/

    // Calculate the exit pressure ratio from the given formula
    return 1.0 / E - pow((gp1 / 2.0), (1.0 / gm1)) * pow(pr, (1.0 / g)) * std::sqrt(std::complex<double>(gp1 / gm1) * (1.0 - pow(pr, gm1g)));
}

std::complex<double> dPe(double E, double gm1, std::complex <double> pr, double g, double gp1,double gm1g) {
// Calculate the numerical approximation of the derivative of the exit pressure ratio with respect to pressure ratio

// step size for numerical approximation
    double step = pow(10, -8);

// approximates the slope of the function at the point (pr, Pe(pr))
    std:: complex<double> numerator = Pe(E, gm1, pr + step, g, gp1, gm1g) - Pe(E, gm1, pr, g, gp1, gm1g); // bruh this is messed up hfuehgfouwehfoewhfoewhfouwehfuoewhfuoewhfuowehfwouefhouwehfwuoefhewuohfewuof
    std:: complex<double> denominator = step;

    // if the denominator = 0, set it to a small value
    if (abs(denominator) < 1e-10) {
        denominator = 1e-10;
    }

    return numerator / denominator;
}

std::complex <double> pratio(double kappa, double E){

    double gm1 = kappa - 1.0;
    double gp1 = kappa + 1.0;
    double gm1g = gm1 / kappa;
    double g = kappa;

// # Set initial conditions
    std::complex <double> pr0 = 0.001; // # Initial guess for the pressure ratio
    double Err = 0.1;  // # Initial error
    std:: complex <double> pr = 0.0; // # Initial pressure ratio

// # Iterate until the desired accuracy is achieved
    int i = 0;
    while (Err > 0.0001) {
// # Calculate the pressure ratio using Newton-Raphson method
        pr = pr0 - Pe(E, gm1, pr0, g, gp1, gm1g) / dPe(E, gm1, pr0, g, gp1, gm1g);// # Newton - Raphson method x(n + 1) = x(n) - f(x(n)) / f'(x(n))
        Err = abs((pr - pr0) / pr0);
        pr0 = pr;
        i += 1;
    }

// # print("Solution of nonlinear equation for Pressure ratio:", pr, "after", i, "iterations")
    return pr;
}

std::complex <double> hybrid_rocket_thrust(double r, double h, const std::vector<double>& OF, std::vector<double>& MOL , std::vector<double>& Gammas, int pressure,int expansionRatio) {

// # Fuel Mass
    double m_f = fuel_mass(r, h); // useless?????????????
    double Pc = pressure; //pressure ratio
    std::complex<double> Thrust = 0;

    double R = 8.3144598; // # J / mol * K

// # Molecular Weights
    double MW_fuel = 352;      // # g / mol, molecular weight of fuel(paraffin)
    double MW_oxidizer = 32;   // # g / mol, molecular weight of oxidizer(O2)

// # Oxidizer Flow Rate (O2)
    double m_o = 0.01429;

// # Mass flux oxidizer (total)
    double j_m = (m_o) / (pi * pow(r,2));

// # Regression rate
    double rdot = 0.304 * pow(j_m,0.527);

// # Area of the fuel grain
    double A_f = fuel_effective_grain_area(r,h);

// # Fuel Mass Flow Rate
    double mdot_f = 900 * A_f * rdot;

// # Mass Flow Rate
    std::complex <double> mdot_total = m_o + mdot_f;

// # Oxidizer to Fuel Ratio
    double OF_ratio = m_o / m_f; // useless???

    for (int i = 0; i<8; i++) {
        double t = OF[i];
        double k = Gammas[i];
        double m = MOL[i];

// # Convert Bar to Pa
        double Pc_Pa = Pc * 100000;

// # Area of the throat: Based on isentropic relations, assumes flow is choked
    std::complex <double> A_t = mdot_total / Pc_Pa * sqrt(8.314 * t / pow((k * 2 / (k + 1)), ((k + 1) / (k - 1)))); // # from rocket propulsion elements fig 3 - 5

// # Calculating the mass fraction of paraffin and o2
        double m_frac_paraffin = MW_fuel/(MW_fuel + MW_oxidizer);
        double m_frac_o2 = MW_oxidizer / (MW_oxidizer + MW_fuel);

// # Calculating average molecular weight
        double ave_mol_wt = (m_frac_paraffin * MW_fuel) + (m_frac_o2 * MW_oxidizer);

// # Using the ideal gas law to find the density of the mixture
        double rho_mixture = (Pc_Pa * ave_mol_wt) / (8.314 * t);
        double rho_oxidizer = (Pc_Pa * MW_oxidizer) / (8.314 * t); // possibly useless??

// # Flow velocity: V = m_dot / (rho * A)
    std::complex <double> V = mdot_total * (rho_mixture * A_t);

// # Speed of sound
        double a = sqrt(k * 8.314 * t);

// # Mach Number
    std::complex <double> M = V / a; // uselesss

//  # Exit Pressure Ratio
        std:: complex<double> Pe = pratio(k, expansionRatio);

// # Calculate the exit velocity
    std::complex <double> v_e = sqrt((t * 8.3144598) / m) * sqrt(((2 * k) / (k - 1))) * (1.0 - pow((Pe / Pc_Pa), ((k - 1) / k))); // double check correclty incoded

// # Calculate the rocket thrust
       Thrust = mdot_total * v_e;
    // # Return the calculated thrust

    return Thrust;

             }
        }



#endif //UNTITLED4_HYBRIDROCKETUPDATE9NOV_H
