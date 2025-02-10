    #include <iostream>
    #include <cmath>
    #include <fstream>
    #include <array>
    #include <functional>
    #include "WhatItTakes.h"

    /// Equations and other methods are stored in WhatItTakes.h

    //Global Variables
    Vector3 thrust_local = {0, thrust, 0};   // Thrust in local frame

    // Atmospheric Conditions
    double g0 = 9.81; // gravity at sea level
    double rho0 = 12.93; // density of air kg/m^3 at sea level
    double hscale = 8500.0; // m, scale of rapid atmospheric change within Earths Atmosphere
    double Re = 6371000.0; // radius of earth from centre to surface
    double CL,CD;

    //function prototypes
    void rk4(Vector3& v, Vector3& h, double& t, double dt, double m, double wx, double wy, double wz);

    // Function to calculate the derivatives of the system
    std::tuple<Vector3, Vector3, Vector3>  derivatives(double t, Vector3 v, Vector3 h, double m, double wx,double wy, double wz) {

        /// figure out intergrating three velocitys and positions etc
        double v_magnitude = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);  // Rocket speed (magnitude of velocity)
        double g = g0 / pow((1 + h.y / Re), 2);    // Gravity at height
        double rho = rho0 * exp(-h.y / hscale);       // air density constant
        double drag = 0.5 * CD * rho * pow(v_magnitude, 2) * A;
        double lift = 0.5 * CL * rho * pow(v_magnitude, 2) * A;
        Vector3 acceleration = {0,0,0};
        Vector3 h_dot = {0,0,0};
        Vector3 dw_dot = {0,0,0};

            // Calculate the net force in each direction based on thrust, drag, lift, and gravity
           /*
            double drag_x = drag * (v.x / v_magnitude);  // Drag component in x
            double drag_y = drag * (v.y / v_magnitude);  // Drag component in y
            double drag_z = drag * (v.z / v_magnitude);  // Drag component in z

            double lift_x = lift * (v.x / v_magnitude);  // Lift component in x
            double lift_y = lift * (v.y / v_magnitude);  // Lift component in y
            double lift_z = lift * (v.z / v_magnitude);  // Lift component in z
*/
            double drag_x =  1;// Drag component in x
            double drag_y = 1;  // Drag component in y
            double drag_z = 1;  // Drag component in z

            double lift_x =1;
            double lift_y = 1;
            double lift_z = 1;

            acceleration.x = thrust_local.x / m - drag_x / m - g + lift_x/m;
            acceleration.y = thrust_local.y / m - drag_y / m - g + lift_y/m;
            acceleration.z = (lift_z - drag_z - g) / m;

            h_dot.x = v.x;
            h_dot.y = v.y;
            h_dot.z = v.z;

            // Get the torque based on thrust, lift, and drag
            auto torque = calculateTorque(lift, drag, radius);

            Inertia inertia = calculateInertia(diam, L, m);

            // Update angular velocities
            dw_dot = calculateAngularAccelerations(inertia.Ix,inertia.Iy,inertia.Iz, wx, wy, wz, torque[0], torque[1], torque[2]);


        return {acceleration, h_dot, dw_dot};
    }

    int main() {
        // Differential inputs
        double t = 0.0; // seconds
        double psi = deg2rad(10); // start at 10 degrees

        Quaternion orientation = {cos(psi / 2), 0, sin(psi / 2), 0}; // Rotation around y-axis

        Vector3 position = {0, 0, 0};   // Initial position
        Vector3 v = {0.1, 0.1, 0.1};         // Initial velocity


        // attempt to rectify code
        double wx = 0;
        double wy = 0;
        double wz = 0;
        Quaternion w = {0,wx,wy,wz};

        // Time parameters
        double dt = 1.0; // time per calculation
        double t_end = 1400.0; // time of simulation

        // Open a file to write data
        std::ofstream outfile("rocket_trajectory.csv");
        outfile << "Time,Velocity,Distance,Height,Angle,Vx,Vy,Psi2\n"; // Headers

        // Load data from the CSV file
        std::vector<CoeffData> data = loadCoeffData("xf-n0012-il-1000000.csv"); // Using NACA 0012 data for reynolds number 10^6
        if (data.empty()) {
            std::cerr << "Failed to load data." << std::endl;
            return 1;
        }

        // Simulation begins, RK4 method loop
        while (t <= t_end) {

            // calculate Lift co-efficient
            getCL_CD(rad2deg(psi), data, CL, CD);

            // Update rocket mass and thrust
            double m = (t <= tburn) ? (m0 - mDot * t) : mstruc;
            thrust_local.y = (t <= tburn) ? thrust: 0.0;

            // Update orientation with angular velocities
            Quaternion orientation_dot = rk4_quaternion_update(orientation,w,dt);
            //printf("Ori: %f ,  %f  , %f,   %f" , orientation.x,orientation.y,orientation.z,   orientation.w);

            // Rotate thrust vector to global frame
            Vector3 thrust_global = rotate_vector(orientation_dot, thrust_local.x,  thrust_local.y, thrust_local.z);

            //update_position using RK4 method
            // Perform RK4 integration step
            rk4(v, position, t, dt,  m,  wx,  wy, wz);

            // Write current time, velocity, and position to file
           // outfile << t << "," << v.y << "," << position.z << "," << rad2deg(psi) << "," << v.z << "," << v.x << "\n";

            // Output the current time, velocity, and position
            std::cout << "Time: " << t << " s, Velocity y : " << v.y << "   Vx :  " << v.x << ", Height: " << position.x << "  y: " << position.y << "  z:  " <<position.z << std::endl;

            // Stop simulation if the rocket hits the ground
            if (position.z < 0) {
                position.z = 0;
                break;
            }
        }

        outfile.close(); //Close file
        std::cout << "Data written to rocket_trajectory.csv" << std::endl;

        return 0;
    }

    Vector3 scaleVector(const Vector3& vec, double scalar) {
        return {vec.x * scalar, vec.y * scalar, vec.z * scalar};
    }


    // Helper function to sum two Vector3 objects
    Vector3 sumVectors(const Vector3& v1, const Vector3& v2) {
        return {v1.x + v2.x, v1.y + v2.y, v1.z+ v2.z};
    }

    // Function to perform the RK4 method Runge Kutta
    void rk4(Vector3& v, Vector3& h, double& t, double dt, double m, double wx, double wy, double wz) {
        // K1 beginning of the interval using Euler's method
        // K2 the midpoint of the interval
        // K3 again midpoint of the interval
        // K4 end of the interval
        // sum weighted v, h, psi, t
        //return {acceleration.x, acceleration.y, acceleration.z,
          //      h_dot.x, h_dot.y, h_dot.z, dw_dot[0],dw_dot[1],dw_dot[2]};


        auto k1 = derivatives(t, v, h, m, wx, wy, wz);

        // preforms h + k1[1] * dt / 2 but needs extra syntax due to Vector3 subtype
        auto vmid = sumVectors(v,scaleVector(std::get<0>(k1), dt / 2));
        auto hmid  = sumVectors(h,scaleVector(std::get<1>(k1), dt / 2));

        auto k2 = derivatives(t + dt / 2, vmid, hmid , m, wx, wy, wz);

        auto vmid2 = sumVectors(v,scaleVector(std::get<0>(k2), dt / 2));
        auto hmid2 = sumVectors(h,scaleVector(std::get<1>(k2), dt / 2));

        auto k3 = derivatives(t + dt / 2, vmid2, hmid2, m,  wx, wy, wz);

        auto vmid3 = sumVectors(v,scaleVector(std::get<0>(k3), dt));
        auto hmid3 = sumVectors(h,scaleVector(std::get<1>(k3), dt));

        auto k4 = derivatives(t + dt, vmid3, hmid3, m, wx, wy, wz);

        //middle points weighted by 1/3 and ends weighted by 1/6 to achieve 4th order accuracy
        // solves v += (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) * dt / 6.0;
        Vector3 Vpart1 = sumVectors(std::get<0>(k1),scaleVector(std::get<0>(k2),2));
        Vector3 Vpart2 = sumVectors(scaleVector(std::get<0>(k3), 2), scaleVector(std::get<0>(k4), dt / 6.0));
        v = sumVectors(v, sumVectors(Vpart1,Vpart2 )); // yn + 1

        // solves h += (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) * dt / 6.0;
        Vector3 hPart1 = sumVectors(std::get<1>(k1),scaleVector(std::get<1>(k2),2));
        Vector3 hPart2 = sumVectors(scaleVector(std::get<1>(k3), 2), scaleVector(std::get<1>(k4), dt / 6.0));
        h = sumVectors(h, sumVectors(hPart1,hPart2 )); // yn + 1

        t += dt;

    }


