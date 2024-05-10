
#include <cstdio>
#include "HybridRocketEquationsHeaderU75Check3.h"

int main() {
    ///U-75 values to aim for thrust = 11 KN, total Impulse = 1100 kN secs
    /// chamber pressure = 5 MPa, Expansion = 100,
    /// diameter = 380 mm, motor length = 1000 mm
    /// burn time = 30secs?, propellent weight = 104 (from calculation)

    // amroc validation numbers
    double Pc = 5000000; // pascals
    double l = 1;   // m
    double RI = 0.2;  // m
    double fuelMass = ( M_PI * pow(RI,2) * l ) * 920; // 104 kg
    double massFlowRate = fuelMass/30;  // 3.477 so at a 1.3 ratio 1.51173913 and 1.965260869
    double oD = 1.965260869;
    double E = 100; // expansion ratio
    double gammas = 1.09;
    double Mol = 0.25;  /// inversly porportional should mol be 0.25
    double OFtemps = 3200;


    double thrust = hybrid_rocket_thrust(RI, l, OFtemps, Mol, gammas, Pc, E,oD);
    printf("Thrust : %f \n", thrust);
}

