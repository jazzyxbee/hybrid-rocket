#include <iostream>
#include <cmath>
#include <fstream>

int main() {
    // This program can be validated with the solving from
    // https://www.youtube.com/watch?v=NFhiJWbNU1g&t=1856s

    // Initial parameters
    double t = 0; // seconds
    double v = 0; // meters per second
    //double Mr_t = 0.05 - 0.002*t; // starts at 50 grams
    double DMr_Dt = -0.002;  // kg/sec mass flow rate
    double Mr_atOneSecond = 0.048; //kg
    double Thrust = 15; // netwons
    double K = 0.02; // Friction

    // Open a file to write the data
    std::ofstream outfile("euler_data.csv");
    outfile << "Time,Velocity\n"; // Header

    // ODE

    // Substituting values

    //double Dv_t =  (Thrust - (Mr_t)*(9.8) - K*pow(v,2) - v*(DMr_Dt))/(Mr_t);

    //DE equals
    double V_t = (14.51 + 0.196*t - 0.02 * pow(v,2) + 0.002 * v)/(0.05 - 0.002*t) ;

    //solving using Eulers

    double dt = 0.01; // time step in seconds
    double t_end = 23.33; // end time for the simulation

    // Function to compute Mr(t)
    auto Mr_t = [](double t) { return 0.05 - 0.002 * t; };

    // Euler's method loop
    while (t <= t_end) {
        // Write current time and velocity to file
        outfile << t << "," << v << "\n";

        // Compute the derivative of v at time t
        double Mr = Mr_t(t);
        double Dv_t = (Thrust - Mr * 9.8 - K * pow(v, 2) - v * DMr_Dt) / Mr;

        // Update the velocity and time using Euler's method
        v = v + Dv_t * dt;
        t = t + dt;

        // Output the current time and velocity
        std::cout << "Time: " << t << " s, Velocity: " << v << " m/s" << std::endl;
    }

    //plotting
    // Plot the results using matplotlibcpp
    outfile.close();
    std::cout << "Data written to euler_data.csv" << std::endl;

    return 0;
}
