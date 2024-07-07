#include <iostream>
#include <cmath>
#include <fstream>
#include <array>

// Function to convert degrees to radians
double deg2rad(double degrees) {
    return degrees * M_PI / 180.0;
}

// Function to convert radians to degrees
double rad2deg(double radians) {
    return radians * 180.0 / M_PI;
}

int main() {
    // Initial parameters
    // Using example values from https://www.youtube.com/watch?v=NFhiJWbNU1g&t=1856s
    double m0 = 0.05; //inital mass kg
    double mDot = 0.002;  // kg/sec mass flow rate
    double fuelAt1Sec = 0.048; //kg
    double thrust = 15; // Netwons
    double K = 0.02; // drag co=efficient
    double g = -9.81;
    double tburn = 2.5;

    double t = 0; // seconds
    double v = 0; // meters per second
    double y = 0; // vertical position in meters
    double x = 0; // horizontal position in meters
    double pitch = 90; // starts at 90 degrees, straight up.
    double pitchRate = -9.8; // gravity

    //solving using Eulers
    double dt = 0.01; // time step in seconds
    double t_end = 2.51; // end time for the simulation

    // Open a file to write the data
    std::ofstream outfile("euler_data.csv");
    outfile << "Time,Velocity,Position,PitchAngle\n"; // Header

    // Euler's method loop
    while (t <= t_end) {
        // Write current time, velocity, and position to file
        outfile << t << "," << v << "," << y << "\n";

        // Compute the derivative of v at time t
        double m = m0 - mDot * t; //Rocket mass

        // Check if Mr is very small or zero to avoid division by a small number
        if (m <= 0.0001) {
            std::cerr << "Mr is too small or zero at t = " << t << std::endl;
            break;
        }

        // Check if the fuel is depleted and set Thrust to 0
        if (m <= fuelAt1Sec) {
            thrust = 0;
        }

        // Calculate pitch angle in radians
        double pitchRad = deg2rad(pitch);

        // Calculate thrust components
        double thrustX = thrust * cos(pitchRad); // Horizontal component of thrust
        double thrustY = thrust * sin(pitchRad); // Vertical component of thrust

        // Compute drag force components
        double vx = v * cos(pitchRad); // Horizontal velocity component
        double vy = v * sin(pitchRad); // Vertical velocity component

        double dragX = K * pow(vx, 2); // Horizontal component of drag
        double dragY = K * pow(vy, 2); // Vertical component of drag

        double Dv_dt, Dv_dtX, Dv_dtY;
        // Compute the derivative of v at time t
        if (thrust > 0){ // before gravity turn
            Dv_dt = (thrust - m * 9.8 - K * pow(v, 2) - v * mDot) / m;
            Dv_dtX = (thrustX - dragX) / m;
            Dv_dtY = (thrustY - dragY - m * 9.8) / m;
        } else{ // after gravity turn
            Dv_dt = (-9.8 - (K / m) * pow(v, 2));
            Dv_dtX = -dragX / m;
            Dv_dtY = (-9.8 - dragY / m);
        }

        // Update velocities using Euler's method
        vx += Dv_dtX * dt;
        vy += Dv_dtY * dt;
        double vsquare = sqrt(vx * vx + vy * vy); // Total velocity

        // Update the velocity and time using Euler's method
        v = v + Dv_dt * dt; // velocity
        x = x + Dv_dtX * dt; // distance
        y = y + v * dt; // height

        // Check if the position has reached or gone below 0
        if (y <= 0) {
            v = 0;
            y = 0;
        }

        t = t + dt;

        // Update pitch angle based on the velocity and height
        if (y > 0) {
            pitch = 90; // Gravity turn approximation but it needs proper representation
        }

        // Output the current time, velocity, and position
        std::cout << "Time: " << t << " s, Velocity: " << v << " Vsquare: " << vsquare << " Acceleration via code " << Dv_dt <<" m/s, Height: " << y << " distance: " << x << " m" << " Pitch Angle: " << pitch <<std::endl;
    }

    //plotting
    // Plot the results using matplotlibcpp
    outfile.close();
    std::cout << "Data written to euler_data.csv" << std::endl;

    return 0;
}
