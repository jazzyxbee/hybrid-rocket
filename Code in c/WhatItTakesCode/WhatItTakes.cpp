#include <iostream>
#include <cmath>
#include <fstream>
#include <array>
#include <functional>
#include "WhatItTakes.h"
#include <cassert>

// Atmospheric Conditions
double g0 = 9.81; // gravity at sea level
double rho = 12.93; // density of air kg/m^3 at sea level
double hscale = 8500.0; // m, scale of rapid atmospheric change within Earths Atmosphere
double Re = 6371000.0; // radius of earth from centre to surface
double CL,CD;

State derivatives(const State& state, double mass, double thrust, double CD, double A, double rho, const Vector3& gravity) {
    State d;

    // Thrust direction in body frame (+Z)
    Vector3 thrust_body = {0, 0, 1};

    Quaternion q = state.orientation;
    Quaternion thrust_quat = q * Quaternion{0, thrust_body.x, thrust_body.y, thrust_body.z} * Quaternion{q.w, -q.x, -q.y, -q.z};
    Vector3 thrust_direction = {thrust_quat.x, thrust_quat.y, thrust_quat.z};
    Vector3 thrust_force = thrust_direction * thrust;

    double v_mag = state.velocity.magnitude();
    Vector3 drag = state.velocity.normalize() * (-0.5 * CD * rho * v_mag * v_mag * A);

    Vector3 accel = (thrust_force + drag + gravity) / mass;

    d.position = state.velocity;
    d.velocity = accel;
    d.orientation = state.orientation.derivative(state.angular_velocity);
    d.angular_velocity = {0, 0, 0};  // For now, zero torque
    //std::cout << "accel: " << accel.x << ", " << accel.y << ", " << accel.z << std::endl;
    return d;
}

void rk4(State& state, double dt, double mass, double thrust, double CD, double A, double rho, const Vector3& gravity) {
    State k1 = derivatives(state, mass, thrust, CD, A, rho, gravity);
    State k2 = derivatives(state + k1 * (dt / 2), mass, thrust, CD, A, rho, gravity);
    State k3 = derivatives(state + k2 * (dt / 2), mass, thrust, CD, A, rho, gravity);
    State k4 = derivatives(state + k3 * dt, mass, thrust, CD, A, rho, gravity);

    state = state + (k1 + k2 * 2 + k3 * 2 + k4) * (dt / 6);
    state.orientation = state.orientation.normalize();
}

int main() {

    // Open a file to write data
    std::ofstream outfile("rocket_trajectory_Quaternions.csv");
    outfile << "Time,Vx,Vy,Vz,x,y,z\n"; // Headers

    // Load data from the CSV file
    std::vector<CoeffData> data = loadCoeffData("xf-n0012-il-1000000.csv"); // Using NACA 0012 data for reynolds number 10^6
    if (data.empty()) {
        std::cerr << "Failed to load data." << std::endl;
        return 1;
    }

    State state = {
            {0, 0, 0},    // Initial position
            {0, 0, 0},    // Initial velocity
            {1, 0, 0, 0}, // Orientation (identity quaternion)
            {0, 0, 0}     // Angular velocity
    };

    Vector3 gravity = {0, 0, -9.81};

    double t = 0.0, dt = 1, t_end = 1400.0;
    std::ofstream fout("trajectory.csv");
    fout << "time,x,y,z,vx,vy,vz\n";

    double psi = 10; /// maybe change this in the future just keeping it here to run code

    // calculate Lift co-efficient
    getCL_CD(rad2deg(psi), data, CL, CD);

    while (t < t_end) {
        // Update rocket mass and thrust
        double m = (t <= tburn) ? (m0 - mDot * t) : mstruc;
        thrust = (t <= tburn) ? thrust: 0;
        if(thrust==0){
            printf("ahhh");
        }

        std::cout << "t=" << t << " vel.z=" << state.velocity.z
                  << " accel.z=" << state.velocity.z
                  << " pos.z=" << state.position.z << std::endl;

        // Output rocket position to file
        fout << t << "," << state.position.x << "," << state.position.y << "," << state.position.z << " , " << state.velocity.x << " , "<< state.velocity.y << " , "<< state.velocity.z << "\n";
        //std::cout << t << "," << state.position.x << "," << state.position.y << "," << state.position.z << "\n";

        rk4(state, dt, m, thrust, CD, A, rho, gravity);
        t += dt;
    }

    fout.close();
    std::cout << "Simulation complete. Data saved to trajectory.csv\n";
    return 0;
}
