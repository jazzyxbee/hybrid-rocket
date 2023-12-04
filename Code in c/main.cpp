
#include <stdio.h>
#include "HybridRocketUpdate9Nov.h"
#include <iterator>
#include <sstream>
#include <algorithm>  // Include the algorithm header for sorting

int main() {

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

    // code to check it is reading the correct values
    /*// Print the results
    std::cout << "OFtemps: ";
    for (const auto &value : OFtemps) {
        std::cout << value << " ";
    }
    std::cout << std::endl;

    std::cout << "gammas: ";
    for (const auto &value : gammas) {
        std::cout << value << " ";
    }
    std::cout << std::endl;

    std::cout << "Mol: ";
    for (const auto &value : Mol) {
        std::cout << value << " ";
    }
    std::cout << std::endl;
     */


    // Chamber Pressures

    int chamberPressures = 1;

    // calculates rocket thrust

    double radii[] = {0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1}; //# fuel grain radii  m
    double heights[] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};  //# fuel grain heights  m
    // make loops
    double expansion_ratio[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}; // should not be iterating here

    std::vector<std::complex<double>> thrustVals;

    for (int i = 0; i < sizeof(radii) / sizeof(radii[0]); i++) {

        double r = radii[i];

        for (int y = 0; y < sizeof(heights) / sizeof(heights[0]); y++) {

            for (double e: expansion_ratio) {
                double h = heights[y];
                std::complex<double> thrust = hybrid_rocket_thrust(r, h, OFtemps, gammas, Mol, chamberPressures, e);

                // Store thrust value in the vector
                thrustVals.push_back(thrust);
            }
        }

        // Sort the thrust values in descending order
       // std::sort(thrustVals.begin(), thrustVals.end(), std::greater<std::complex<double>>());

// Print the top 5 thrust values
    }
    size_t n = thrustVals.size();
    size_t k = std::min<size_t>(5, n);  // min ensures k is within the bounds of the vector size

    std::nth_element(thrustVals.begin(), thrustVals.begin() + k, thrustVals.end(),
                     [](const std::complex<double>& a, const std::complex<double>& b) {
                         return std::abs(a) > std::abs(b);
                     });

    // Output the top five values
    std::cout << "Top five values of thrustVals:" << std::endl;
    for (size_t i = 0; i < k; ++i) {
        std::cout << "Value " << i + 1 << ": " << thrustVals[i] << std::endl;
    }
}


