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
    double drag = (0.5) * CD * rho * pow(v, 2) * A;

    //printf("%f \n", drag);
    //double lift = (0.5) * CL * rho * pow(v, 2) * AL; //not accurate so removing for now
    double Dv_dt, h_dot, psi_dot, theta_dot;
    double phi_dot = 0; // defined here due to scope



  if(t <= 40) { // before gravity turn currently at karman line possibly to high??
          Dv_dt = T / m - drag / m - g ;
          h_dot = v ;
          theta_dot = 0;

    } else { // after gravity turn
        Dv_dt = T / m - drag / m - g * cos(psi);
        h_dot = v * cos(psi);
        theta_dot = (v * sin(psi)) / (Re + h);
        phi_dot = g * sin(psi) / v; // defined here due to scope
        psi_dot = phi_dot - theta_dot;
    }

    double targetPsi = psi; // default is "hold current heading"
    if (t > 144.611 && t < 167.67) {
        targetPsi = deg2rad(-90);
        // Simple proportional control law to steer to target psi
        double psi_rate_control = (targetPsi - psi) * 0.05;
        psi_dot = psi_rate_control;
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
    //SpaceX CRS11 flight Falcon 9 Flight Vehicle
    // Initial parameters
    // LOX/RP-1
    //Total values
    double m = 0.0; // initial mass
    double massTotal = 549054; //kg

    //first stage
    double mpropFirstStage = 395700.0; // kg propellent mass
    double tburnFirstStage = 162; // burn time (s)
    double emptyMassFirstStage = 25600; // kg
    double thrustAtSeaLevelFirstStage = 7607000.0; // Netwons
    double thrustInVaccumFirstStage = 8227000; // Newtons
    double Isp_FirstStage = 282; // s
    double mDotFirstStage = mpropFirstStage / tburnFirstStage; // kg/s propellent mass flow rate

    //second stage
    double mpropSecondStage = 92670; // kg
    double tburnSecondStage = 397.0; // burn time (s)
    double thrustInVaccumSecondStage = 981000; // Newtons
    double Isp_SecondStage = 348; // s
    double mDotSecondStage = mpropSecondStage / tburnFirstStage; // kg/s propellent mass flow rate

    //payload
    double mPayload = 22800.0; // kg payload mass
    double mstruc = massTotal -  mpropFirstStage - mpropSecondStage - mPayload;
    double m0 = mpropFirstStage + mpropSecondStage + mPayload + mstruc;  // initialmass

    //thrust levels
    double T = thrustAtSeaLevelFirstStage; // initial thrust

    // Atmospheric Conditions
    double g0 = 9.81; // gravity at sea level
    double rho0 = 12.93; // density of air kg/m^3 at sea level
    double hscale = 8500.0; // m, scale of rapid atmospheric change within Earths Atmosphere

    // Size parameters
    double diam = 3.7; // m rocket diameter
    double length = 70; // m length
    double A = M_PI/4 * pow(diam,2);
    double Re = 6371000.0; // radius of earth from centre to surface
    double AL = M_PI * pow((diam/2),2); // Lift area
    double payloadHeight = 13.1; // m
    double payloadDiameter = 5.2; // m

    // Differential inputs
    double t = 0.0; // seconds
    double v = 1.0; // meters per second
    double h = 0.0; // height position in meters
    double psi = deg2rad(8); // start at 90 degrees

    double CL,CD;

    // Time parameters
    double dt = 1.0; // time per calculation
    double t_end = 450.0; // time of simulation

    // Open a file to write data
    std::ofstream outfile("rocket_trajectory.csv");
    outfile << "Time,Velocity,Distance,Height,Angle,Vx,Vy,Psi2\n"; // Headers

    // Load data from the CSV file
    std::vector<CoeffData> data = loadCoeffData("xf-n0012-il-1000000.csv"); // Using NACA 0012 data for reynolds number 10^6
    if (data.empty()) {
        std::cerr << "Failed to load data." << std::endl;
        return 1;
    }

    bool rotatedForBoostback = false;

    // RK4 method loop
    while (t <= t_end) {
        // calculate Lift co-efficient
        //getCL_CD(rad2deg(psi), data, CL, CD); // always zero at the moment
        CD = 0.3; // validation with youtube
        CL = 0.2; //
        //printf("%f \n",CD);

        double fuelBurned = 0;

        // Update rocket mass and thrust
        if (t < 144.611) {
            //Rocket launch
            fuelBurned = mDotFirstStage * t;
            m = m0 - fuelBurned;
            T = (t <= tburnFirstStage) ? thrustAtSeaLevelFirstStage : 0.0;

        } else if (t > 144.611 && t < 162.168) {
            //Rocket second stage ejection, MECO, and rotation
            m = m0 - mPayload - mpropSecondStage - fuelBurned;
            T = 0;
        } else if (t > 162.168&& t < 210.77) {
            double deltaT = t - (162.168-144.611);
            //boostback burn ignition approximately three engines NASA Falcon9 Datasheet
            m = m0 - mPayload - mpropSecondStage - fuelBurned - (mDotFirstStage/3)*deltaT;
            T = thrustInVaccumFirstStage/10;
            //T = 0;
        } else if (t > 210.77 && t < 373) {
            //boostback burn cuts off
            T = 0;
        }else if (t>373){
            v = 0;
        }

        // Components of velocity     okay i dont know why but these two are diffenitley switched around and vx is vy etc
        // i need to figure out why this is happening
        //double vx = v * cos(psi);
        //double vy = v * sin(psi);
        double vy = v * cos(psi);
        double vx = v * sin(psi);
        double psi2 = atan2(vy, vx);

        // Perform RK4 integration step
        rk4(v, h, psi, t, dt, m, T, Re, g0, hscale, rho0, A, CD,CL,AL);

/*
        if (t>144.611 && t<167.67) {
            //  if (!rotatedForBoostback) {
            double targetPsi = deg2rad(130);  // This targets *opposite total velocity vector*
            psi += (targetPsi - psi) * 0.05;        // apply 5% of the rotation per step               // rotatedForBoostback = true;
        }
           // }

       } else if (t>167.67 && t<200) {
           //  if (!rotatedForBoostback) {
           double targetPsi = deg2rad(-90);  // This targets *opposite total velocity vector*
           psi += (targetPsi - psi) * 0.08;        // apply 5% of the rotation per step               // rotatedForBoostback = true;
           //}
       }else if (t>200 && t<250) {
           //  if (!rotatedForBoostback) {
           double targetPsi = deg2rad(-70);  // This targets *opposite total velocity vector*
           psi += (targetPsi - psi) * 0.08;        // apply 5% of the rotation per step               // rotatedForBoostback = true;
           //}
       }*/
            //printf("%f %n",   rad2deg(psi2));

        //double vtot = sqrt(pow(vx,2) + pow(vy,2));

        // Calculate distance traveled along the Earth's surface
        double theta = (v * sin(psi)) / (Re + h) * dt;
        double dr = theta * Re / 1000;

        // Write current time, velocity, and position to file
        outfile << t << "," << v << "," << dr << "," << h << "," << rad2deg(psi) << "," << vx << "," << vy << "," << rad2deg(psi2) << "\n";

        // Output the current time, velocity, and position
        //std::cout << "Time: " << t << " s, Velocity: " << v << "  Acceleration:   " << v*t << "  Height: " << h << " m, Angle: " << rad2deg(psi) << " degrees" << std::endl;

        // Stop simulation if the rocket hits the ground
        /*
        if (h < 0) {
          //  h = 0;
            break;
        }*/
    }

    outfile.close(); //Close file
    std::cout << "Data written to rocket_trajectory.csv" << std::endl;

    return 0;
}