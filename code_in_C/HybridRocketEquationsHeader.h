//
// Created by User on 09/11/2023.
//

//#ifndef UNTITLED4_HYBRIDROCKETUPDATE9NOV_H
//#define UNTITLED4_HYBRIDROCKETUPDATE9NOV_H

#include <iostream>
//#include <stdio.h>
//#include <math.h>
//#include <iterator>
//#include <vector>
//#include <complex>



/*double function(double M, double y, double aRatio) {
    return (1 / pow(M, 2)) * pow(((2 / (y - 1)) * (1 + ((y - 1) / 2) * pow(M, 2))), ((y + 1) / (y - 1))) - pow(aRatio, 2);
}

double derivative(double M, double y) {
    double term1 = -2 * (y + 1) / ((y - 1) * M * M * M);
    double term2 = (y + 1) / (y - 1);
    double term3 = (1 + (y - 1) * M * M / 2);
    double term4 = (2 / (y - 1)) * (1 + (y - 1) * M * M / 2);

    return term1 * pow(term2 * term3, term4);
}

double newtonRaphson(double y, double aRatio, double initial_guess, double tolerance, int max_iterations) {
    double M = initial_guess;
    int iterations = 0;

    while (iterations < max_iterations) {
        double f = function(M, y, aRatio);
        double f_prime = derivative(M, y);
        double delta_M = -f / f_prime;

        M += delta_M;

        if (fabs(delta_M) < tolerance) {
            return M;
        }

        iterations++;
    }

    // If we reach here, the method did not converge
    return NAN; // Not a Number
}

double exitPressure(double y) {

    // # y = specific heat ratio gamma

    // # solve exit mach number
    double aRatio = 3;
    double initial_guess = 4; // Initial guess for M
    double tolerance = 1e-6;  // Tolerance for convergence
    int max_iterations = 100; // Maximum iterations

    double root1 = newtonRaphson(y, aRatio, initial_guess, tolerance, max_iterations);
    double root2 = newtonRaphson(y, aRatio, -initial_guess, tolerance, max_iterations);

    if (root1 != NAN) {
        printf("Root 1: %lf\n", root1);
    } else {
        printf("Root 1 not found.\n");
    }

    if (root2 != NAN) {
        printf("Root 2: %lf\n", root2);
    } else {
        printf("Root 2 not found.\n");
    }

    double Me = root1;

    // # return Pratio = exit pressure as inside pressure is assumed p0 = 1
    return 1 + pow((((y - 1) / 2) * pow(Me, 2)), (-y / (y - 1)));
}*/

/*double hybrid_rocket_thrust(double r, double L, double OF, double MOL, double y, int pressure) {

    // # equations from Rocket Propulsion Elements unless otherwise stated

    // -------------------------------Calculating Mdot propellant mass flow rate------------------------//

    // from Hybrid Rocket Fuel Regression Rate Data and Modeling
    // # regression rate coefficient of
    double a = 0.488;

    // from Hybrid Rocket Fuel Regression Rate Data and Modeling
    // # regression rate pressure
    double n = 0.62;

    // # Oxidizer Flow Rate (O2)
    double mDot_o = 0.01429;

    // # fuel density of parrafin kg/m^3
    int Pf = 900;

    // #  Combustion port surface area: L = port length  r = radius of fuel grain
    double A_f = 2 * M_PI * r * L;

    // # Oxidizer mass velocity (16-2)
    double G_o = mDot_o / M_PI * pow(r, 2);

    // # Regression rate (16-10)
    double rdot = a * pow(G_o, n);/// suspiciously small

    // # Mass Flow Rate of Fuel (16 - 11)
    double mDot_f = Pf * A_f * rdot;

    // # Mass Flow Rate (6-2)
    double mDot = mDot_o + mDot_f; /// too small?

    // -------------------------------Calculating Ve Exit velocity------------------------------------//

    // # Gas constant
    double R = 8.3144598; // J / mol * K

    // # absolute temperature at nozzle inlet
    double T1 = 1528.6;

    // # specific heat ratios
    double k = 1.2983;

    // # chamber pressure
    double Pc = pressure;

    /// # Convert Bar to Pa do I need this ????
   // Pc = Pc * 100000;

    // # Appendix 2 PROPERTIES OF THE EARTH’S STANDARD ATMOSPHERE with idealised flow p2 = p3
    // # exit pressure
    double Pe =  0.88700;    /// produces 31.5269 seems resonable

            /// exitPressure(y);  unsure if a method is required yet while we are still in the design process

    // # Calculate the exit velocity
    double v_e = sqrt(((2 * k) / (k - 1)) * (R * T1) * (1.0 - pow((Pe / Pc), ((k - 1) / k))));

    // # Calculate the rocket thrust from Space-Propulsion-Analysis-and-Design-McGraw-Hill-(1995
    double Thrust = mDot * v_e;

    // # Return the calculated thrust
    return Thrust;
}*/

int main() {
    std::cout << "Hello, World!" << std::endl;
    return 0;
}

//#endif //UNTITLED4_HYBRIDROCKETUPDATE9NOV_H
