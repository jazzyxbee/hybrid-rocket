//
// Created by User on 10/02/2025.
//

#ifndef WHATITTAKES_WHATITTAKES_H
#define WHATITTAKES_WHATITTAKES_H

typedef struct { double x, y, z; } Vector3;

// Based on the Titan 2 and Gemini Spacecraft
// -----------------------------------------------Rocket parameters stored in a struct--------------------------------
double mprop = 111130.0;                   // kg propellent mass
double mpl = 32000.0;                      // kg payload mass
double mstruc = 6736.0;                    // kg structure mass
double m0 = mprop + mpl + mstruc;   // Initial mass = propellant + payload + structure
double tburn = 356.0;                      // burn time in seconds
double thrust = 1900000.0;                 // Newtons
double mDot = 0.0; // kg/s propellent mass flow rate
double diam = 3.05;                        // rocket diameter meters
double radius = 3.05/2;
double A = M_PI * pow((3.05 / 2), 2);      // m^2 (cross-sectional area or area of nose cone)
double L = 33.0;                            // meters (rocket length)

// -----------------------------------------------------------Getting Co-Effecients-----------------------------------
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

// -----------------------------------------------------------Conversion formulas -----------------------------------------------------------
// Function to convert degrees to radians
double deg2rad(double degrees) {
    return degrees * M_PI / 180.0;
}

// Function to convert radians to degrees
double rad2deg(double radians) {
    return radians * 180.0 / M_PI;
}

// -----------------------------------------------------------ANGULAR EQUATIONS -----------------------------------------------------------
// Function to calculate angular accelerations based on current angular velocities and torques
Vector3 calculateAngularAccelerations(double Ix,double Iy, double Iz, double omega_x, double omega_y, double omega_z, double tau_x, double tau_y, double tau_z) {
    double w_x_dot, w_y_dot, w_z_dot;
    // Euler's equations of motion
    w_x_dot = (tau_x + (Iz - Iy) * omega_y * omega_z) / Ix;
    w_y_dot = (tau_y + (Ix - Iz) * omega_z * omega_x) / Iy;
    w_z_dot = (tau_z + (Iy - Ix) * omega_x * omega_y) / Iz;
    Vector3 w_dot = {w_x_dot,w_y_dot,w_z_dot};
    return {w_dot};
}

struct Inertia {
    double Ix;
    double Iy;
    double Iz;
};

Inertia calculateInertia(double diameter, double length, double mass) {
    // use a number as a placeholder
    // Calculate moments of inertia
    Inertia inertia;
    inertia.Ix = (1.0 / 12.0) * mass * (3 * pow(radius, 2) + pow(length, 2));
    inertia.Iy = (1.0 / 12.0) * mass * (3 * pow(radius, 2) + pow(length, 2));
    inertia.Iz = (1.0 / 2.0) * mass * pow(radius, 2);

    return inertia;
}

// Function to calculate torque
std::array<double, 3> calculateTorque(double lift, double drag, double radius) {
    double tau_x = lift * radius; // Assuming lift acts at a distance (e.g., half the diameter)
    double tau_y = -drag * radius; // Assuming drag acts in the opposite direction
    double tau_z = 0.0; // If there are no other moments about the z-axis

    return {tau_x, tau_y, tau_z};
}

// -----------------------------------------------------------PITCHING MOMENT EQUATIONS -----------------------------------------------------------

// Function to calculate pitching moment coefficient (C_M)
double calculateCM(double alpha){
    // Example formula replace this with empirical data or look-up tables
    return 0.1 * alpha; // A simple linear relationship for small angles
}

// Function to calculate the pitching moment
double calculatePitchingMoment(double velocity, double alpha, double rho, double diam, double L) {
    // maybe not needed
    double S = M_PI*pow((diam/2),2); // reference area
    double CM = calculateCM(alpha); // moment co-effecient
    double dynamicPressure = 0.5 * rho * std::pow(velocity, 2);
    double M = CM * dynamicPressure * S * L;  // calculating the moment
    return M;
}

//----------------------------------------------------------- Quaternion related code -----------------------------------------------------------PITCHING

typedef struct { double w, x, y, z; } Quaternion;

// Quaternion functions (multiply, conjugate, rotate vector, etc.)
Quaternion quaternion_multiply(Quaternion q1, Quaternion q2) {
    Quaternion result;
    result.w = q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z;
    result.x = q1.w * q2.x + q1.x * q2.w + q1.y * q2.z - q1.z * q2.y;
    result.y = q1.w * q2.y - q1.x * q2.z + q1.y * q2.w + q1.z * q2.x;
    result.z = q1.w * q2.z + q1.x * q2.y - q1.y * q2.x + q1.z * q2.w;
    return result;
}

Quaternion quaternion_conjugate(Quaternion q) {
    Quaternion result = {q.w, -q.x, -q.y, -q.z};
    return result;
}

// Update quaternion using angular velocities
Vector3 rotate_vector(Quaternion q, double vx, double vy, double vz) {
    Quaternion v = {0, vx, vy, vz};
    Quaternion q_conj = quaternion_conjugate(q);
    Quaternion rotated_q = quaternion_multiply(quaternion_multiply(q, v), q_conj);

    return {rotated_q.x, rotated_q.y, rotated_q.z};
}

Quaternion quaternion_scale(Quaternion q, double scale) {
    Quaternion result;
    result.w = q.w * scale;
    result.x = q.x * scale;
    result.y = q.y * scale;
    result.z = q.z * scale;
    return result;
}

Quaternion quaternion_add(Quaternion q1, Quaternion q2) {
    Quaternion result;
    result.w = q1.w + q2.w;
    result.x = q1.x + q2.x;
    result.y = q1.y + q2.y;
    result.z = q1.z + q2.z;
    return result;
}

// Normalize quaternion to avoid drift
Quaternion quaternion_normalize(Quaternion q) {
    double norm = sqrt(q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z);
    return {q.w / norm, q.x / norm, q.y / norm, q.z / norm};
}

// Runge-Kutta 4 method for quaternion update
Quaternion rk4_quaternion_update(Quaternion q, Quaternion omega_q, double dt) {
    Quaternion k1 = quaternion_scale(quaternion_multiply(q, omega_q), 0.5 * dt);
    Quaternion k2 = quaternion_scale(quaternion_multiply(quaternion_add(q, k1), omega_q), 0.5 * dt);
    Quaternion k3 = quaternion_scale(quaternion_multiply(quaternion_add(q, k2), omega_q), dt);
    Quaternion k4 = quaternion_scale(quaternion_multiply(quaternion_add(q, k3), omega_q), dt);

    Quaternion q_final = quaternion_add(q, quaternion_scale(k1, 1.0 / 6.0));
    q_final = quaternion_add(q_final, quaternion_scale(k2, 2.0 / 6.0));
    q_final = quaternion_add(q_final, quaternion_scale(k3, 2.0 / 6.0));
    q_final = quaternion_add(q_final, quaternion_scale(k4, 1.0 / 6.0));

    return quaternion_normalize(q_final);
}



#endif //WHATITTAKES_WHATITTAKES_H
