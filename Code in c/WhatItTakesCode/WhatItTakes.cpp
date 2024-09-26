#include <iostream>
#include <cmath>
#include <fstream>
#include <array>
#include <functional>
// Structure to hold angle, CL, and CD values
struct CoeffData {
    double angle;
    double CL;
    double CD;
};

struct Rocket {
    double m0;        // Initial mass
    double mprop;     // Propellant mass
    double mpl;       // Payload mass
    double mstruc;    // Structural mass
    double tburn;     // Burn time
    double thrust;    // Thrust force
    double mDot;      // Mass flow rate
    double diam;      // Diameter
    double A;         // Area (cross-sectional)
    double L;         // Length of rocket
};


// Function to load data from CSV file
std::vector<CoeffData> loadCoeffData(const std::string& filename) {
    std::vector<CoeffData> data;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return data;
    }

    std::string line;
    while (std::getline(file, line)) {
        CoeffData dp;
        if (sscanf(line.c_str(), "%lf,%lf,%lf", &dp.angle, &dp.CL, &dp.CD) == 3) {
            data.push_back(dp);
        }
    }

    file.close();
    return data;
}

void getCL_CD(double angle, std::vector<CoeffData>& data, double& CL, double& CD){
   // Interpolate the value of CL and CD based on a given angle CD and CL from data
   for(int i = 0; i < data.size()-1; i++){
       if(angle >= data[i].angle && angle < data[i+1].angle ) {
           double slopeCL = ((data[i+1].CL-data[i].CL)/(data[i+1].angle - data[i].angle));
           double slopeCD = ((data[i+1].CD-data[i].CD)/(data[i+1].angle - data[i].angle));
           CL = data[i].CL + (slopeCL) * (angle - data[i].angle);
           CD = data[i].CD + (slopeCD) * (angle - data[i].angle);
           return;
       }
   }
   //sets the values to zero if out of range
   CL = 0.0;
   CD = 0.0;
}

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
    double drag = (1 / 2) * CD * rho * pow(v, 2) * A;
    double lift = (1 / 2) * CL * rho * pow(v, 2) * AL;
    double Dv_dt, h_dot, psi_dot, theta_dot;

    if (h <= 15000 && T>0) { // before gravity turn currently at karman line possibly to high??
        Dv_dt = T / m - drag / m - g + lift/m;
        h_dot = v;
        psi_dot = 0;
        theta_dot = 0;

    } else { // after gravity turn
        // look at doing it differently here

        Dv_dt = T / m + drag / m + lift/m - g * cos(psi);
        double phi_dot = g * sin(psi) / v; // defined here due to scope
        h_dot = v * cos(psi);
        theta_dot = (v * sin(psi)) / (Re + h);
        psi_dot = phi_dot - theta_dot;

        //double vx = v * cos(psi);
        //double vy = v * sin(psi);
        //psi_dot = psi + atan2(vx,vy);
        //psuedo code
        //  use trig here to solve psi dot then let rk4 approximate the true value?
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

// Function to calculate pitching moment coefficient (C_M)
double calculateCM(double alpha){
    // Example formula replace this with empirical data or look-up tables
    return 0.1 * alpha; // A simple linear relationship for small angles
}

// Function to calculate the pitching moment
double calculatePitchingMoment(double velocity, double alpha, double rho, double diam, double L) {
    double S = M_PI*pow((diam/2),2); // reference area
    double CM = calculateCM(alpha); // moment co-effecient
    double dynamicPressure = 0.5 * rho * std::pow(velocity, 2);
    double M = CM * dynamicPressure * S * L;  // calculating the moment
    return M;
}

int main() {
    // Based on the Titan 2 and Gemini Spacecraft
    // Rocket parameters stored in a struct
    Rocket titan2 {
            .m0 = 111130 + 32000 + 6736,   // Initial mass = propellant + payload + structure
            .mprop = 111130.0,                   // kg propellent mass
            .mpl = 32000.0,                      // kg payload mass
            .mstruc = 6736.0,                    // kg structure mass
            .tburn = 356.0,                      // burn time in seconds
            .thrust = 1900000.0,                 // Newtons
            .mDot = titan2.mprop / titan2.tburn, // kg/s propellent mass flow rate
            .diam = 3.05,                        // rocket diameter meters
            .A = M_PI * pow((3.05 / 2), 2),      // m^2 (cross-sectional area or area of nose cone)
            .L = 33.0                            // meters (rocket length)
    };

    double m = 0.0; // resetting active mass to 0
    double T = titan2.thrust;

    // Atmospheric Conditions
    double g0 = 9.81; // gravity at sea level
    double rho0 = 12.93; // density of air kg/m^3 at sea level
    double hscale = 8500.0; // m, scale of rapid atmospheric change within Earths Atmosphere
    double Re = 6371000.0; // radius of earth from centre to surface

    // Differential inputs
    double t = 0.0; // seconds
    double v = 0.0; // meters per second
    double h = 0.0; // height position in meters
    double psi = deg2rad(10); // start at 10 degrees
    double CL,CD;

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

    // RK4 method loop
    while (t <= t_end) {
        // calculate Lift co-efficient
        getCL_CD(rad2deg(psi), data, CL, CD);
        //printf("%f",CL);

        // Update rocket mass and thrust
        m = (t <= titan2.tburn) ? (titan2.m0 - titan2.mDot * t) : titan2.mstruc;
        T = (t <= titan2.tburn) ? titan2.thrust : 0.0;

        // Perform RK4 integration step
        rk4(v, h, psi, t, dt, m, T, Re, g0, hscale, rho0, titan2.A, CD, CL, titan2.A);

        // Components of velocity
        double vx = v *cos(psi);
        double vy = v *sin(psi);
        double psi2 = atan2(vy, vx);
        //printf("%f %n",   rad2deg(psi2));

        double vtot = sqrt(pow(vx,2) + pow(vy,2));

        // Calculate distance traveled along the Earth's surface
        double theta = (v * sin(psi)) / (Re + h) * dt;
        double dr = theta * Re / 1000;

        // Calculate pitching moment
        //double pitchingMoment = calculatePitchingMoment(v, theta, rho0, titan2.diam, titan2.L);
        double pitchingMoment = 0;


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
