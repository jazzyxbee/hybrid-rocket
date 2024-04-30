
#include <stdio.h>
#include "HybridRocketEquationsHeader.h"
#include <iterator>
#include <sstream>
#include <algorithm>  // Include the algorithm header for sorting

int main() {

    //Our rocket numbers
    //int cP = 1;  // Chamber Pressures

    // size parameters
    //double RI = 0.5; // at time zero internal radius of fuel grain m
    //double RE = 1;  //external radius of rocket
    //double l = 1; // our rocket

    //double oD = 85.740; // maximum oxidizer flow rate kg/s
    /// previously we had 0.01429 as the flow rate which is the density of oxygen but also the max flow rate by chance
    //Double validated
    /// 1.429 kg/m^3 oxygen density
    /// convert from L to m^3 1000 L = 1 m^3
    /// 60000 L/s = 60 m^3/s
    /// 60 m^3/s * 1.429 kg/m^3 = 85.74 kg/s

    // amroc validation numbers
    double cP = 2.406e+6; // pascals N/m^2 average chamber pressure
    double l = 1.756;   // length of rocket in meters currently 1.756 m to validate with amroc rocket
    double RI = 0.5588;  // amroc validation number
    double oD = 16.344; // flow rate needed for nominal 85 sec burn

   /// 75 for amroc conformation
    int expansionRatio[] = {2,10,30,40,50,60,70,80,90,75};  // iterate through this to find different thrust
    // # Oxidizer Flow Rate (O2) (lbm/sec) max of testing chamber is 100 L/minute or 3.67436 lbs/s
    double oxidizerFlowRate[] = {oD/10, oD/5,oD/3.33333, oD/2.5, oD}; /// from Geraldine test chamber at robison capable of 60000 L/s
    ///nitrous oxide used in amroc but its density is 1.22 so 1.7 of oxygen is a close approximation
    /// amroc rocket takes nominal 85 sec burn time
    /// volume of fuel = pi * r^2 * L * density of fuel = pi * 0.5588^2 * 1.756 * 920 (HTPB) = 1584.801 kg of fuel
    /// unit check = m^3 * kg/m^3
    /// mdot f  is 2.3 kg/s so already we will burn through in s Volume/flowrate = T , so 1584.801/85 = Flow rate = 18.64 kg/s
    /// 18.64 kg/s - 2.3 kg/s = 16.344 kg/s from oxidizer flow rate

    //reading files
    // Initialize vectors to store data
    std::vector<double> OFtemps;  // oxygen fuel temperature, temperature inside chamber during reaction
    std::vector<double> gammas;   // kappa
    std::vector<double> Mol;      /// mols used in reaction, molar mass so is this equal to mdot total??
    std::vector<double> Thrust;

    // Sample data
    std::string data = "1528.6   1.2983  14.856\n"
                       "3017.45  1.1291  20.850\n"
                       "3118.76  1.1123  23.909\n"
                       "3074.36  1.1103  25.919\n"
                       "3012.54  1.1104  27.371\n"
                       "2947.19  1.1116  28.461\n"
                       "2880.86  1.1140  29.296";

    // Create a stringstream to read from the data string
    std::istringstream iss(data);

    // Process each line
    std::string line;
    while (std::getline(iss, line)) {
        std::istringstream lineStream(line);
        double temp, gamma, mol;

        // Read the three numbers from each line
        lineStream >> temp >> gamma >> mol;

        // Add the values to respective vectors
        OFtemps.push_back(temp);
        gammas.push_back(gamma);
        Mol.push_back(mol);
    }
    //for (int i =0; i < 10; i++) {
      //  printf("For Expansion ratio : %i\n", expansionRatio[i]);
        //for (int j = 0; j < OFtemps.size(); j++) {
          //  for (int z = 0; z < 4; z++) {
                //double thrust = hybrid_rocket_thrust(RI, l, OFtemps[j], Mol[j], gammas[j], cP, expansionRatio[i],oxidizerFlowRate[z]);
                ///amroc validation
                double thrust = hybrid_rocket_thrust(RI, l, OFtemps[3], Mol[3], gammas[3], cP, expansionRatio[9],oxidizerFlowRate[4]);
                Thrust.push_back(thrust);
            //}
        //}
    //}

    // returns largest value from the values
    auto max_element_iterator = std::max_element(Thrust.begin(), Thrust.end());
    if (max_element_iterator != Thrust.end()) {
        double largest_element = *max_element_iterator;
        //std::cout << "Largest Thrust: " << largest_element << std::endl;
    } else {
        std::cout << "Vector is empty." << std::endl;
    }
}

