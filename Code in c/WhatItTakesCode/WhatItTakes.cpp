#include <iostream>
#include <cmath>
#include <fstream>
#include <array>
#include <functional>
#include "WhatItTakes.h"

/// Equations and other methods are stored in WhatItTakes.h

//Global Variables
double T = thrust;

// Atmospheric Conditions
double g0 = 9.81; // gravity at sea level
double rho0 = 12.93; // density of air kg/m^3 at sea level
double hscale = 8500.0; // m, scale of rapid atmospheric change within Earths Atmosphere
double Re = 6371000.0; // radius of earth from centre to surface
double CL,CD;

//function prototypes
void rk4(double& v, double& h, double& psi, double& theta, double& phi, double& t, double dt, double m, double wx, double wy, double wz);

// Function to calculate the derivatives of the system
std::array<double, 10> derivatives(double t, double v, double h, double psi, double theta, double phi, double m, double wx,double wy, double wz) {
    double g = g0 / pow((1 + h / Re), 2);    // Gravity at height
    double rho = rho0 * exp(-h / hscale);       // air density constant
    double drag = 0.5 * CD * rho * pow(v, 2) * A;
    double lift = 0.5 * CL * rho * pow(v, 2) * A;
    double Dv_dt, h_dot, psi_dot, theta_dot;
    double dphi, dtheta, dpsi; // For Euler Angles

    if (h <= 15000 && T>0) { // before gravity turn currently at karman line possibly to high??
        Dv_dt = T / m - drag / m - g + lift/m;
        h_dot = v;

        // constant thrust results in no angle changes
        psi_dot = 0;
        theta_dot = 0;
        dphi = 0;
        dtheta = 0;
        dpsi = 0;

    } else { // after gravity turn
        // look at doing it differently here

        Dv_dt = T / m + drag / m + lift/m - g * cos(psi);

        // No thrust results in angle changes
        double phi_dot = g * sin(psi) / v; // defined here due to scope
        h_dot = v * cos(psi);
        theta_dot = (v * sin(psi)) / (Re + h);
        psi_dot = phi_dot - theta_dot;

        //double vx = v * cos(psi);
        //double vy = v * sin(psi);
        //psi_dot = psi + atan2(vx,vy);

        // Get the torque based on thrust, lift, and drag
        auto torque = calculateTorque(lift, drag, radius);

        Inertia inertia = calculateInertia(diam, L, m);

        // Update angular velocities
        std::array<double, 3> dw_dot = calculateAngularAccelerations(inertia.Ix,inertia.Iy,inertia.Iz, wx, wy, wz, torque[0], torque[1], torque[2]);

        // Integrate the angular velocity to angular displacement
        dphi = dw_dot[0]*t;
        dtheta = dw_dot[1]*t;
        dpsi = dw_dot[2]*t;

        // Compute the time derivatives of the Euler angles based on angular velocities
        //dphi = wx + sin(phi) * tan(theta) * wy + cos(phi) * tan(theta) * wz;
        //dtheta = cos(phi) * wy - sin(phi) * wz;
        //dpsi = (sin(phi) / cos(theta)) * wy + (cos(phi) / cos(theta)) * wz;

        // psi_dot = dpsi;
    }
    return {Dv_dt, h_dot, psi_dot, theta_dot,dphi,dtheta,dpsi};
}

int main() {

    // Differential inputs
    double t = 0.0; // seconds
    double v = 0.0; // meters per second
    double h = 0.0; // height position in meters
    double psi = deg2rad(10); // start at 10 degrees

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
        T = (t <= tburn) ? thrust : 0.0;

        // Components of velocity
        double vx = v *cos(psi);
        double vy = v *sin(psi);
        double psi2 = atan2(vy, vx);
        //printf("%f %n",   rad2deg(psi2));
        double vtot = sqrt(pow(vx,2) + pow(vy,2));

        // Calculate distance traveled along the Earth's surface
        double theta = (v * sin(psi)) / (Re + h) * dt;
        double dr = theta * Re / 1000;

        double phi; // idk how to acturally get this yet

        // Calculate pitching moment
        double pitchingMoment = calculatePitchingMoment(v, theta, rho0, diam, L);
        //double pitchingMoment = 0;

        // attempt to rectify code
        double wx = 1;
        double wy = 1;
        double wz = 1;
        // Perform RK4 integration step
        rk4( v,  h,  psi, theta, phi,  t, dt,  m,  wx,  wy, wz);

        // Write current time, velocity, and position to file
        outfile << t << "," << v << "," << dr << "," << h << "," << rad2deg(psi) << "," << vy << "," << vx << "," << rad2deg(psi2) << "\n";

        // Output the current time, velocity, and position
        std::cout << "Time: " << t << " s, Velocity: " << v << " m/s, Height: " << h << " m, Angle: " << rad2deg(psi) << " degrees" << "Pitching Moment: " << pitchingMoment << std::endl;

        // Stop simulation if the rocket hits the ground
        if (h < 0) {
            h = 0;
            break;
        }
    }

    outfile.close(); //Close file
    std::cout << "Data written to rocket_trajectory.csv" << std::endl;

    return 0;
}

// Function to perform the RK4 method Runge Kutta
void rk4(double& v, double& h, double& psi, double& theta, double& phi, double& t, double dt, double m, double wx, double wy, double wz) {
    // K1 beginning of the interval using Euler's method
    // K2 the midpoint of the interval
    // K3 again midpoint of the interval
    // K4 end of the interval
    // sum weighted v, h, psi, t

    auto k1 = derivatives(t, v, h, psi, theta, phi, m, wx, wy, wz);
    auto k2 = derivatives(t + dt / 2, v + k1[0] * dt / 2, h + k1[1] * dt / 2, psi + k1[2] * dt / 2,theta + k1[3] * dt / 2, phi + k1[4] * dt / 2, m, wx, wy, wz);
    auto k3 = derivatives(t + dt / 2, v + k2[0] * dt / 2, h + k2[1] * dt / 2, psi + k2[2] * dt / 2,theta + k2[3] * dt / 2, phi + k2[4] * dt / 2, m,  wx, wy, wz);
    auto k4 = derivatives(t + dt, v + k3[0] * dt, h + k3[1] * dt, psi + k3[2] * dt, theta + k3[3] * dt,phi + k3[4] * dt, m, wx, wy, wz);

    //middle points weighted by 1/3 and ends weighted by 1/6 to achieve 4th order accuracy
    v += (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) * dt / 6.0; // yn + 1
    h += (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) * dt / 6.0;
    psi += (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]) * dt / 6.0;
    theta += (k1[3] + 2 * k2[3] + 2 * k3[3] + k4[3]) * dt / 6.0;
    phi += (k1[4] + 2 * k2[4] + 2 * k3[4] + k4[4]) * dt / 6.0;
    t += dt;

    // Update angular position (Euler angles)
    double phi_increment = (k1[4] + 2 * k2[4] + 2 * k3[4] + k4[4]) * dt / 6.0;
    double theta_increment = (k1[5] + 2 * k2[5] + 2 * k3[5] + k4[5]) * dt / 6.0;
    double psi_increment = (k1[6] + 2 * k2[6] + 2 * k3[6] + k4[6]) * dt / 6.0;

    phi += phi_increment;
    theta += theta_increment;
    psi += psi_increment;
}

