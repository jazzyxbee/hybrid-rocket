
#include <cstdio>
#include "HybridRocketEquationsHeaderAquillaValidation.h"

int main() {
    ///AMROC values to aim for thrust = 9000 lbs, vaccum Isp = 288 secs
    /// burn time = 30s, chamber pressure = 5 MPa, Expansion = 75, propellent weight = 2670 lbs
    /// diameter = 22 inches, motor length = 69 inches, mixture ratio = 1.3

    // amroc validation numbers
    double Pc = 349; // pascals N/m^2 average chamber pressure (conversion 2.406e+6 pascals from 349 psia)
    double l = 69;   // length of rocket in meters currently 1.756 m to validate with amroc rocket
    double RI = 22;  // amroc validation number in inches
    double E = 75; // expansion ratio
    double gammas = 1.197;
    double Mol = 23;
    double OFtemps = 2800;
    double OxidizerFlowRate = 31.41; //14.25 kg as per belows calculations

    ///nitrous oxide used in amroc but its density is 1.22 so 1.7 of oxygen is a close approximation
    /// amroc rocket takes nominal 85 sec burn time
    /// volume of fuel = pi * r^2 * L * density of fuel = pi * 0.5588^2 * 1.756 * 920 (HTPB) = 1584.801 kg of fuel
    /// unit check = m^3 * kg/m^3
    ///an alterante way of validating this flow rate is the mass of the fuel is 2670 pounds = 1211 kg since the entire rocket isnt fuel
    ///this makes sense so the fuel rate is (1211kg/85s) 14.25 this way .
    /// I cannot confirm from this data what the individual flow rates are.


    double thrust = hybrid_rocket_thrust(RI, l, OFtemps, Mol, gammas, Pc, E,OxidizerFlowRate);
    printf("Thrust : %f \n", thrust);
}

