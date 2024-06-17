
#include <stdio.h>
#include "HybridRocketEquationsHeader.h"
#include <iterator>
#include <sstream>
#include <fstream>
#include <algorithm>


    // Create a vector of DataPoint structures to store the data
    std::vector<dataPoint> dataPoints;

void ReadInData() {
    //----------------------------------reading data---------------------------------------------------------------//
    // Sample data with expansion ratios
    std::string data = "1528.6   1.2983  14.856  1\n"
                       "3017.45  1.1291  20.850  2\n"
                       "3118.76  1.1123  23.909  3\n"
                       "3074.36  1.1103  25.919  4\n"
                       "3012.54  1.1104  27.371  5\n"
                       "2947.19  1.1116  28.461  6\n"
                       "2880.86  1.1140  29.296  7";

    // Create a stringStream to read data
    std::istringstream iss(data);

    // Add each value to the vector of DataPoint structures
    std::string line;
    while (std::getline(iss, line)) {
        std::istringstream lineStream(line);
        dataPoint dp;

        // Read the four numbers from each line
        lineStream >> dp.temperature >> dp.gamma >> dp.molecular_mass >> dp.expansionRatio;

        // Add the DataPoint to the vector
        dataPoints.push_back(dp);
    }
}

int main() {
    //Read the data into Struct
    ReadInData();

    //Rocket Parameters
    int Pc = 1;                                 // Chamber Pressures

    // size parameters
    double RI = 0.2;                            // at time zero internal radius of fuel grain m
    double l = 1;                               // length
    double oD = 0.0238166666l;                   // maximum oxidizer flow rate kg/s
    double volume = M_PI * RI * RI * l;         // Volume of the cylindrical tank
    double mass = volume * 900;                 // around 113 kg

    int expansionRatio[] = {2,10,30,40,50,60,70,80,90,100};
    // # Oxidizer Flow Rate (O2) (lbm/sec) max of testing chamber is 100 L/minute or 3.67436 lbs/s
    double oxidizerFlowRate[] = {oD/10, oD/8,oD/6, oD/2, oD*2, oD*4, oD*8, oD*10}; /// from Geraldine test chamber at robison capable of 60000 L/s

    // Initialize vectors to store data
    std::vector<double> Thrust;
    // create a thrust_data file to store results

    std::ofstream outFile("thrust_data.csv");
    outFile << "ExpansionRatio,OxidizerFlowRate,Thrust\n"; //csv header

    // Iterate over data points
    for (const auto& data : dataPoints) {
        // iterate through expansion ratio, vectored data, and oxidizer flow rates
        for (int i = 0; i < 10; i++) {
            //printf("For Expansion ratio : %i\n", expansionRatio[i]);
            for (int z = 0; z < 8; z++) {

                double thrust = hybrid_rocket_thrust(dataPoints, data.gamma,data.molecular_mass,RI, l, Pc, expansionRatio[i], oxidizerFlowRate[z]);
                //double thrust = hybrid_rocket_thrust(dataPoints, data.gamma,data.molecular_mass,RI, l, Pc, data.expansionRatio, oxidizerFlowRate[z]);

                //printf("mol %f",Mol[j]);
                //printf("temp %f",OFtemps[j]);
                //printf("Thrust : %f \n", thrust);
                Thrust.push_back(thrust);

                outFile << expansionRatio[i] << "," << oxidizerFlowRate[z] << "," << thrust
                        << "\n"; // adds data to csv file
            }
        }
    }

    // returns largest value from calculated Thrust for each expansion ratio
    auto max_element_iterator = std::max_element(Thrust.begin(), Thrust.end());
    if (max_element_iterator != Thrust.end()) {
        double largest_element = *max_element_iterator;
        std::cout << "Largest Thrust: " << largest_element << std::endl;
    } else {
        std::cout << "Vector is empty." << std::endl;
    }

    //Validating mass flow rate based on maximum flow rate of Geraldine
    /*
    1.429 kg/m^3 oxygen density
    convert from L to m^3: 1000 L/min = 1 m^3 / min
    1 m^3/min * 1.429 kg/m^3 = 1.429 kg/min
    minutes to seconds = 1.429/60 =  0.0238166666 kg/s
    */
}

