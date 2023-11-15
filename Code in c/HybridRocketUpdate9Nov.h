//
// Created by User on 09/11/2023.
//

#ifndef UNTITLED4_HYBRIDROCKETUPDATE9NOV_H
#define UNTITLED4_HYBRIDROCKETUPDATE9NOV_H
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <iterator>


double pi = 3.14159265359;

double fuel_effective_grain_area(int r){
    return pi * pow(r,2);//A = pi * r ^ 2
}

double fuel_mass(double r, double h) {
    double paraffin_d = 900; // kg / m3\n
    return (pi * pow(r,2) * h) * paraffin_d; // pi * r ^ 2 * h * paraffin density
}

//Solving for Pe using the Newton-Raphson method
double Pe(double E, double gm1, double pr, double g, double gp1, double gm1g){
    // Calculate the exit pressure ratio from the given formula
    return 1.0 / E - (pow((gp1 / 2.0),(1.0 / gm1)) * (pow(pr,(1.0 / g)) * sqrt((gp1 / gm1) * (1.0 - (pow(pr,gm1g))))));
}

double dPe(double E, double gm1, double pr, double g, double gp1,double gm1g) {
// Calculate the numerical approximation of the derivative of the exit pressure ratio with respect to pressure ratio

// step size for numerical approximation
    int step = pow(10, -8);

// approximates the slope of the function at the point (pr, Pe(pr))
    double numerator = Pe(E, gm1, pr + step, g, gp1, gm1g) - Pe(E, gm1, pr, g, gp1, gm1g);
    double denominator = step;

    // if the denominator = 0, set it to a small value
    if (abs(denominator) < 1e-10) {
        denominator = 1e-10;
    }

    return numerator / denominator;
}

double pratio(double kappa, double E){

    double gm1 = kappa - 1.0;
    double gp1 = kappa + 1.0;
    double gm1g = gm1 / kappa;
    double g = kappa;

// # Set initial conditions
    double pr0 = 0.001; // # Initial guess for the pressure ratio
    double Err = 0.1;  // # Initial error
    double pr = 0.0; // # Initial pressure ratio

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

double hybrid_rocket_thrust(double r, double h, double T_array[256], double M_array[256], double G_array[256], double Pc_array[256]) {

// # Fuel Mass
    double m_f = fuel_mass(r, h);
    //printf("%f",m_f);

    double R = 8.3144598; // # J / mol * K

// # Molecular Weights
    int MW_fuel = 352;      // # g / mol, molecular weight of fuel(paraffin)
    int MW_oxidizer = 32;   // # g / mol, molecular weight of oxidizer(O2)

// # Oxidizer Flow Rate (O2)
    double m_o = 0.01429;

// # Mass flux oxidizer (total)
    double j_m = (m_o) / (pi * pow(r,2));

// # Regression rate
    double rdot = 0.304 * pow(j_m,0.527);

// # Area of the fuel grain
    double A_f = fuel_effective_grain_area(r);

// # Fuel Mass Flow Rate
    double mdot_f = 900 * A_f * rdot;

// # Mass Flow Rate
    double mdot_total = m_o + mdot_f;

// # Oxidizer to Fuel Ratio
    double OF_ratio = m_o / m_f;
    printf("%f",OF_ratio);

// make loops

    for (int tIndex = 0; T_array[tIndex] != -1; tIndex++) {
        double t = T_array[tIndex];

        for (int mIndex = 0; M_array[mIndex] != -1; mIndex++) {
            double m = M_array[mIndex];

            for (int kIndex = 0; G_array[kIndex] != -1; kIndex++) {
                double k = G_array[kIndex];

                for (int pcIndex = 0; Pc_array[pcIndex] != -1; pcIndex++) {
                    double Pc = Pc_array[pcIndex];

// # Convert Bar to Pa
                    double Pc_Pa = Pc * 100000;

// # Area of the throat
                    double A_t = mdot_total / Pc_Pa * sqrt(8.314 * t / pow((k * 2 / (k + 1)),(k + 1)) /(k - 1));//# from rocket propulsion elements fig 3 - 5

// # Calculating the mass fraction of paraffin and o2
                    double m_frac_paraffin = MW_fuel / (MW_fuel + MW_oxidizer);
                    double m_frac_o2 = MW_oxidizer / (MW_oxidizer + MW_fuel);

// # Calculating average molecular weight
                    double ave_mol_wt = (m_frac_paraffin * MW_fuel) + (m_frac_o2 * MW_oxidizer);

// # Using the ideal gas law to find the density of the mixture
                    double rho_mixture = (Pc_Pa * ave_mol_wt) / (8.314 * t);
                    double rho_oxidizer = (Pc_Pa * MW_oxidizer) / (8.314 * t);

// # Area of the throat: Based on isentropic relations, assumes flow is choked
                    A_t = mdot_total / Pc_Pa * sqrt(8.314 * t / pow((k * 2 / (k + 1)),((k + 1) / (k - 1)))); // # from rocket propulsion elements fig 3 - 5
                    printf("%d", A_t);

// # Flow velocity: V = m_dot / (rho * A)
                    double V = mdot_total * (rho_mixture * A_t);

// # Speed of sound
                    double a = sqrt(k * 8.314 * t);

// # Mach Number
                    double M = V / a;

//  # Expansion Ratio, from rocket propulsion elements example 3-1 (implies no losses)
                    double E = (1 / M) * pow(((1 + ((k - 1 / 2) * pow(M,2))) / (k + 1 / 2)),((k + 1) / (2 * (k - 1))));

//  # Exit Pressure Ratio
                    double Pe = pratio(k, E);

// # Calculate the exit velocity
                    double v_e = sqrt((t * 8.3144598) / m) * sqrt((2 * k / k - 1)) * (1 - pow((Pe / Pc_Pa),((k - 1) / k)));

// # Calculate the rocket thrust
                    double Thrust = (m_o + m_f) * v_e;

// # Return the calculated thrust

                    return Thrust;
                }
            }
        }
    }
}


#endif //UNTITLED4_HYBRIDROCKETUPDATE9NOV_H
