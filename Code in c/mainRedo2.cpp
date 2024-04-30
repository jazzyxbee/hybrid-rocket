
#include <cstdio>
#include "HybridRocketEquationsHeaderRedo2.h"

int main() {
    ///AMROC values to aim for thrust = 9000 lbs, vaccum Isp = 288 secs
    /// burn time = 85s, chamber pressure = 349 psia, Expansion = 75, propellent weight = 2670 lbs
    /// diameter = 22 inches, motor length = 69 inches

    /// current WIP investigate mass flow rate formulas and ensure the correct mass flow rate of 34.41 from calculations

    // amroc validation numbers
    double Pc = 349; // pascals N/m^2 average chamber pressure
    double l = 69;   // length of rocket in meters currently 1.756 m to validate with amroc rocket
    double RI = 22;  // amroc validation number
    double oD = 0; /// 0 while i test flow rate from pounds/secs = 2670 / 85 // flow rate needed for nominal 85 sec burn
    double E = 75; // expansion ratio
    double gammas = 1.197;
    double Mol = 24.39;
    double OFtemps = 2800;


    double thrust = hybrid_rocket_thrust(RI, l, OFtemps, Mol, gammas, Pc, E,oD);
    printf("Thrust : %f \n", thrust);
}

