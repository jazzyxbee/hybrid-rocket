#include <iostream>
#include <cmath>
#include <fstream>

int main() {
    // Initial parameters
    // Using example values from https://www.youtube.com/watch?v=NFhiJWbNU1g&t=1856s
    double t = 0; // seconds
    double v = 0; // meters per second
    double y = 0; // meters
    double mDot = -0.002;  // kg/sec mass flow rate
    double fuelAt1Sec = 0.048; //kg
    double Thrust = 15; // Netwons
    double K = 0.02; // Friction

    // Open a file to write the data
    std::ofstream outfile("euler_data.csv");
    outfile << "Time,Velocity,Position\n"; // Header

    //solving using Eulers

    double dt = 0.01; // time step in seconds
    double t_end = 3; // end time for the simulation

    // Euler's method loop
    while (t <= t_end) {
        // Write current time, velocity, and position to file
        outfile << t << "," << v << "," << y << "\n";

        // Compute the derivative of v at time t
        double Mr = 0.05 - 0.002 * t; //Rocket mass

        // Check if Mr is very small or zero to avoid division by a small number
        if (Mr <= 0.0001) {
            std::cerr << "Mr is too small or zero at t = " << t << std::endl;
            break;
        }

        // Check if the fuel is depleted and set Thrust to 0
        if (Mr <= fuelAt1Sec) {
            Thrust = 0;
        }

        double Dv_dt;
        // Compute the derivative of v at time t
        if (Thrust > 0){
            Dv_dt = (Thrust - Mr * 9.8 - K * pow(v, 2) - v * mDot) / Mr;
        } else{
            Dv_dt = (-9.8 - (K / Mr) * pow(v, 2));
        }

        // Update the velocity and time using Euler's method
        v = v + Dv_dt * dt;
        y = y + v * dt;

        // Check if the position has reached or gone below 0
        if (y <= 0) {
            v = 0;
            y = 0;
        }

        t = t + dt;

        // Output the current time, velocity, and position
        std::cout << "Time: " << t << " s, Velocity: " << v << " m/s, Position: " << y << " m" << std::endl;
    }

    //plotting
    // Plot the results using matplotlibcpp
    outfile.close();
    std::cout << "Data written to euler_data.csv" << std::endl;

    return 0;
}
