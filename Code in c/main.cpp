
#include <stdio.h>
#include "HybridRocketUpdate9Nov.h"
#include <iterator>
#include <vector>
#include <algorithm>  // Include the algorithm header for sorting

int main() {

    //reading files

    // OF temps
    printf("Select temps file: ");

    //char Tdata[256]; //previous method of obtaining file paths
    //scanf("%s", Tdata);

    //currently hard coding in file paths due to file finding errors
    char Tdata[256] = "C:/Users/User/OneDrive - Victoria University of Wellington - STUDENT/B.Thrust Vectoring project/CODE/HybridRocketUpdate9Nov/OFtemps.txt";

    printf("Selected file path: %s\n", Tdata);

    //opening file to read
    FILE *Tfile = fopen(Tdata, "r");

    //error checking file
    if (Tfile == NULL) {
        printf("Error: Unable to open file\n");
        exit(1);
    }

    //Creating an array to hold the values from the files
    double T_array[256];  // Assuming a maximum of 256 values
    int T_count = 0;
    while (fscanf(Tfile, "%lf", &T_array[T_count]) == 1) {
        T_count++;
    }
    fclose(Tfile);

    // Molar mass of exhaust gases. Same process as OF temps.

    printf("Select Molar Mass file: ");

    //char Mdata[256];
    //scanf("%s", Mdata);

    char Mdata[256] = "C:/Users/User/OneDrive - Victoria University of Wellington - STUDENT/B.Thrust Vectoring project/CODE/HybridRocketUpdate9Nov/MolarMasses.txt";

    printf("Selected file path: %s\n", Mdata);

    FILE *Mfile = fopen(Mdata, "r");

    if (Mfile == NULL) {
        printf("Error: Unable to open file\n");
        exit(1);
    }

    double M_array[256];  // Assuming a maximum of 256 values
    int M_count = 0;

    while (fscanf(Mfile, "%lf", &M_array[M_count]) == 1) {
        M_count++;
    }
    fclose(Mfile);

    // Gamma

    printf("Select Gamma file: ");

    //char Gdata[256];
    //scanf("%s", Gdata);

    char Gdata[256] = "C:/Users/User/OneDrive - Victoria University of Wellington - STUDENT/B.Thrust Vectoring project/CODE/HybridRocketUpdate9Nov/Gammas.txt";
    printf("Selected file path: %s\n", Gdata);

    FILE *Gfile = fopen(Gdata, "r");

    if (Gfile == NULL) {
        printf("Error: Unable to open file\n");
        exit(1);
    }

    double G_array[256];  // Assuming a maximum of 256 values
    int G_count = 0;

    while (fscanf(Gfile, "%lf", &G_array[G_count]) == 1) {
        G_count++;
    }

    fclose(Gfile);

    // Chamber Pressures

    printf("Select Chamber Pressures file: ");

    //char PcData[256];
    //scanf("%s", PcData);

    char PcData[256] = "C:/Users/User/OneDrive - Victoria University of Wellington - STUDENT/B.Thrust Vectoring project/CODE/HybridRocketUpdate9Nov/ChamberPressures.txt";

    printf("Selected file path: %s\n", PcData);

    FILE *Pcfile = fopen(PcData, "r");

    if (Pcfile == NULL) {
        printf("Error: Unable to open file\n");
        exit(1);
    }

    double Pc_array[256];  // Assuming a maximum of 256 values
    int Pc_count = 0;

    while (fscanf(Pcfile, "%lf", &Pc_array[Pc_count]) == 1) {
        Pc_count++;
    }
    fclose(Pcfile);

    // calculates rocket thrust

    double radii[] = {0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1} ; //# fuel grain radii  m
    double heights[] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};  //# fuel grain heights  m

    std::vector<double> thrustVals;

    for (int i = 0; i < sizeof(radii) / sizeof(radii[0]); i++) {

        double r = radii[i];

        for (int y = 0; y < sizeof(heights) / sizeof(heights[0]); y++) {

            double h = heights[y];
            double thrust = hybrid_rocket_thrust(r, h, T_array, M_array, G_array, Pc_array);
            printf("%f\n", thrust);

            // Store thrust value in the vector
            thrustVals.push_back(thrust);
        }
    }
    // Sort the thrust values in descending order
    std::sort(thrustVals.begin(), thrustVals.end(), std::greater<double>());

// Print the top 5 thrust values
    printf("Top 5 Thrust Values\n");
    for (int i = 0; i < std::min(5, static_cast<int>(thrustVals.size())); i++) {
        printf("Thrust = %f\n", thrustVals[i]);
    }
    }


