
#include <stdio.h>
#include "HybridRocketUpdate9Nov.h"
#include <iterator>
#include <sstream>
#include <algorithm>  // Include the algorithm header for sorting

int main() {

    // Chamber Pressures
    int cP = 1;
    // size parameters
    double r = 0.1;  // radii
    double l = 1.0;  // length

    //reading files
    // Initialize vectors to store data
    std::vector<double> OFtemps;
    std::vector<double> gammas;
    std::vector<double> Mol;

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

    double thrust = hybrid_rocket_thrust(r,l, OFtemps[0], Mol[0] , gammas[0],  cP);
    std::cout << thrust << std::endl;
    }

