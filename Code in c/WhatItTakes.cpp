#include <iostream>
#include <cmath>
#include <fstream>
#include <array>
#include <functional>

// Function to convert degrees to radians
double deg2rad(double degrees) {
    return degrees * M_PI / 180.0;
}

// Function to convert radians to degrees
double rad2deg(double radians) {
    return radians * 180.0 / M_PI;
}

// Function to calculate the derivatives of the system
std::array<double, 4> derivatives(double t, double v, double h, double psi, double m, double T, double Re, double g0, double hscale, double rho0, double A, double CD, double CL, double AL) {
    double g = g0 / pow((1 + h / Re), 2);
    double rho = rho0 * exp(-h / hscale);
    double drag = (1 / 2) * rho * pow(v, 2) * A * CD;
    double lift = (1 / 2) * CL * rho * pow(v, 2) * AL;
    double Dv_dt, h_dot, psi_dot, theta_dot;

    if (h <= 1000000) { // before gravity turn currently at karman line possibly to high??
        Dv_dt = T / m - drag / m - g;
        h_dot = v;
        psi_dot = 0;
        theta_dot = 0;
    } else { // after gravity turn
        double phi_dot = g * sin(psi) / v; // defined here due to scope
        Dv_dt = T / m + drag / m - g * cos(psi);
        h_dot = v * cos(psi);
        theta_dot = (v * sin(psi)) / (Re + h);
        psi_dot = phi_dot - theta_dot;
    }

    return {Dv_dt, h_dot, psi_dot, theta_dot};
}

// Function to perform the RK4 method
void rk4(double& v, double& h, double& psi, double& t, double dt, double m, double T, double Re, double g0, double hscale, double rho0, double A, double CD, double CL,double AL) {
    // K1 beginning of the interval using Eulers method
    // K2 the midpoint of the interval
    // K3 again midpoint of the interval
    // K4 end of the interval
    // sum weighted v, h, psi, t

    auto k1 = derivatives(t, v, h, psi, m, T, Re, g0, hscale, rho0, A, CD, CL, AL);
    auto k2 = derivatives(t + dt / 2, v + k1[0] * dt / 2, h + k1[1] * dt / 2, psi + k1[2] * dt / 2, m, T, Re, g0, hscale, rho0, A, CD, CL, AL);
    auto k3 = derivatives(t + dt / 2, v + k2[0] * dt / 2, h + k2[1] * dt / 2, psi + k2[2] * dt / 2, m, T, Re, g0, hscale, rho0, A, CD, CL, AL);
    auto k4 = derivatives(t + dt, v + k3[0] * dt, h + k3[1] * dt, psi + k3[2] * dt, m, T, Re, g0, hscale, rho0, A, CD, CL, AL);

    //middle points weighted by 1/3 and ends weighted by 1/6 to achieve 4th order accuracy
    v += (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) * dt / 6.0; // yn + 1
    h += (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) * dt / 6.0;
    psi += (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]) * dt / 6.0;
    t += dt;
}

int main() {
    // Based on the Titan 2 and Gemini Spacecraft

    // Initial parameters
    double m = 0.0; // initial mass
    double mprop = 111130.0; // kg propellent mass 111130.0
    double mpl = 32000.0; // kg payload mass
    double mstruc = 6736.0; // kg structure mass
    double m0 = mprop + mpl + mstruc; // initialmass
    double tburn = 356.0; // burn time (s)
    double mDot = mprop / tburn; // kg/s propellent mass flow rate
    double thrust = 1900000.0; // Netwons
    double T = thrust;
    //double hturn = 1000.0; // experimenting with this but it still fails sometimes

    // Atmospheric Conditions
    double CD = 0.3; // drag coefficient
    double g0 = 9.81; // gravity at sea level
    double rho0 = 12.93; // density of air kg/m^3 at sea level
    double hscale = 8500.0; // m, scale of rapid atmospheric change within Earths Atmosphere

    // Size parameters
    double diam = 3.05; // m rocket diameter
    double A = M_PI * pow((diam / 2), 2); // area of the nose (frontal portion)
    double Re = 6371000.0; // radius of earth from centre to surface
    double AL = M_PI * pow((diam/2),2); // Lift area

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
    outfile << "Time,Velocity,Distance,Height,Angle\n"; // Headers

    // RK4 method loop
    while (t <= t_end) {
        // calculate Lift co-efficient
        double CL = psi * (2 * M_PI); // approximating it using thin air foil theory

        // Update rocket mass and thrust
        m = (t <= tburn) ? (m0 - mDot * t) : mstruc;
        T = (t <= tburn) ? thrust : 0.0;

        // Perform RK4 integration step
        rk4(v, h, psi, t, dt, m, T, Re, g0, hscale, rho0, A, CD,CL,AL);

        // Calculate distance traveled along the Earth's surface
        double theta = (v * sin(psi)) / (Re + h) * dt;
        double dr = theta * Re / 1000;

        // Write current time, velocity, and position to file
        outfile << t << "," << v << "," << dr << "," << h << "," << rad2deg(psi) << "\n";

        // Output the current time, velocity, and position
        std::cout << "Time: " << t << " s, Velocity: " << v << " m/s, Height: " << h << " m, Angle: " << rad2deg(psi) << " degrees" << std::endl;

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
