
#include <stdio.h>
#include "HybridRocketEquationsHeader.h"
#include <iterator>
#include <sstream>
#include <fstream>
#include <algorithm>

int main() {

    //Our rocket numbers
    int Pc = 1;  // Chamber Pressures

    // size parameters
    double RI = 0.2; // at time zero internal radius of fuel grain m
    double l = 1; // our rocket
    double oD = 0.0238166666; // maximum oxidizer flow rate kg/s
    double volume = M_PI * RI * RI * l; // Volume of the cylindrical tank
    double mass = volume * 900; // around 113 kg

    //Validating mass flow rate based on maximum flow rate of Geraldine
    /*
    1.429 kg/m^3 oxygen density
    convert from L to m^3: 1000 L/min = 1 m^3 / min
    1 m^3/min * 1.429 kg/m^3 = 1.429 kg/min
    minutes to seconds = 1.429/60 =  0.0238166666 kg/s
    */

    int expansionRatio[] = {2,10,30,40,50,60,70,80,90,100};
    // # Oxidizer Flow Rate (O2) (lbm/sec) max of testing chamber is 100 L/minute or 3.67436 lbs/s
    double oxidizerFlowRate[] = {oD/10, oD/8,oD/6, oD/2, oD*2, oD*4, oD*8, oD*200}; /// from Geraldine test chamber at robison capable of 60000 L/s

    //reading files
    // Initialize vectors to store data
    std::vector<double> OFtemps;  // oxygen fuel temperature, temperature inside chamber during reaction
    std::vector<double> gammas;   // kappa
    std::vector<double> Mol;      // mols used in reaction, molar mass
    std::vector<double> Thrust;

    // Sample data
    std::string data = "1528.6   1.2983  14.856\n"
                       "3017.45  1.1291  20.850\n"
                       "3118.76  1.1123  23.909\n"
                       "3074.36  1.1103  25.919\n"
                       "3012.54  1.1104  27.371\n"
                       "2947.19  1.1116  28.461\n"
                       "2880.86  1.1140  29.296";

    // Create a stringStream to read from the data string
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

    std::ofstream outFile("thrust_data.csv");
    outFile << "ExpansionRatio,OxidizerFlowRate,Thrust\n"; //csv header

    for (int i =0; i < 10; i++) {
        printf("For Expansion ratio : %i\n", expansionRatio[i]);
        for (int j = 0; j < OFtemps.size(); j++) {
           for (int z = 0; z < 8; z++) {

               double thrust = hybrid_rocket_thrust(RI, l, OFtemps[j], Mol[j], gammas[j], Pc, expansionRatio[i],oxidizerFlowRate[z]);

                //printf("mol %f",Mol[j]);
                //printf("temp %f",OFtemps[j]);
                printf("Thrust : %f \n", thrust);
               Thrust.push_back(thrust);
               outFile << expansionRatio[i] << "," << oxidizerFlowRate[z] << "," << thrust << "\n"; // adds data to csv file
            }
        }
    }
    // returns largest value from the values
    auto max_element_iterator = std::max_element(Thrust.begin(), Thrust.end());
    if (max_element_iterator != Thrust.end()) {
        double largest_element = *max_element_iterator;
        std::cout << "Largest Thrust: " << largest_element << std::endl;
    } else {
        std::cout << "Vector is empty." << std::endl;
    }

}

