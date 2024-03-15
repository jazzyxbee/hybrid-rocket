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
#include <complex>


double exitPressure(double y);

double exitMachNumber(double y);

double areaMachNumberRelation(double y, double m, double ratio);

double AM_EQN(double y, double m, double ratio);

double hybrid_rocket_thrust(double r, double L, double OF, double MOL , double y, int pressure) {

    // # equations from Rocket Propulsion Elements unless otherwise stated

    //initializing thrust to zero
    double Thrust = 0;

    // -------------------------------Calculating Mdot propellant mass flow rate------------------------//

    // from Hybrid Rocket Fuel Regression Rate Data and Modeling
    // # regression rate coefficient of
    double a = 0.488;

    // from Hybrid Rocket Fuel Regression Rate Data and Modeling
    // # regression rate pressure
    double n = 0.62;

    ///? forgot where from ask
    // # Oxidizer Flow Rate (O2)
    double mDot_o = 0.01429;

    ///? forgot where from ask
    // # fuel density of parrafin kg/m^3
    int Pf = 900;

    // #  Combustion port surface area: L = port length  R = inner radius
    double A_f = 2 * M_PI * r * L;

    // # Oxidizer mass velocity (16-2)
    double G_o = mDot_o / M_PI * pow(r,2);

    // # Regression rate (16-10)
    double rdot = a * pow(G_o,n);

    // # Mass Flow Rate of Fuel (16 - 11)
    double mDot_f = Pf * A_f * rdot;

    // # Mass Flow Rate (6-2)
    double mDot = mDot_o + mDot_f;

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

        // # exit pressure
        double Pe = exitPressure(y);  ///?

        // # Calculate the exit velocity
        double v_e =  sqrt(( ( 2*k )/( k-1 ) ) * (R*T1) * (1.0 - pow((Pe / Pc), ( (k-1)/k ))));

        // # Calculate the rocket thrust from Space-Propulsion-Analysis-and-Design-McGraw-Hill-(1995
        Thrust = mDot * v_e;

        // # Return the calculated thrust
        return Thrust;
    }

double exitPressure(double y) {

    // # y = specific heat ratio gamma

    // # solve exit mach number
    double Me = exitMachNumber(y);

    // # return Pratio
    return 1 + pow((((y - 1)/2) * pow(Me,2)), (-y/(y-1)));
}

// # calculate the exitMachNumber
double exitMachNumber(double y) {
   /// current issue: iterations dont reach below zero so no roots are ever found
    // # area ratio between any area of the nozzle over the throat: NASA area ratio - A/*A
    // # A/A* = {[(gam+1)/2]^-[(gam+1)/(gam-1)/2]} / M * [1 + M^2 * (gam-1)/2]^[(gam+1)/(gam-1)/2]
    // # assumed 1 for this test
    double aRatio = 1;
    double dM = 0.1;
    double M = 1e-6;
    double iConvSub = 0;
    double iterMax = 100;
    double stepMax = 100;

    for (int i = 0; i <= iterMax; i++) {
        double fj = AM_EQN(y, M, aRatio);
        double fjp1 = AM_EQN(y, M + dM, aRatio);

        for (int j = 0; j <= stepMax; j++) {
            if (fj * fjp1 > 0) {
                M = M + dM;
            } else if (fj * fjp1 < 0) {
                dM = dM * 0.1;
            }
        }

        if (fabs(fj - fjp1) <= iConvSub) { // Adjust tolerance as needed
            iConvSub = i;
        }
    }

    /// maybe needs to be float?
    double Msub = M;
    return Msub;
}

// # return estimated area mach number relation to find ideal exit Mach number
double AM_EQN(double y, double M, double aRatio) {
    return ( (1/pow(M,2)) * pow (( (2/(y-1)) * ( 1 + ((y-1)/2) * pow(M,2)))  ,  ((y+1) / (y-1)) ) - pow(aRatio,2) );
}


#endif //UNTITLED4_HYBRIDROCKETUPDATE9NOV_H
