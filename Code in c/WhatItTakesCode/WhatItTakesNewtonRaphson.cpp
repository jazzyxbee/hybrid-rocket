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
    // Using example values from https://www.youtube.com/watch?v=22OCPbfY5SE&t=616s

    double m = 0.0; // initial mass
    double mprop = 111130.0; // kg propellent mass 111130.0
    double mpl = 32000.0; // kg payload mass
    double mstruc = 6736.0; // kg structure mass
    double m0 = mprop + mpl + mstruc; // initialmass
    double tburn = 356.0; // burn time (s)
    double mDot = mprop / tburn; // kg/s propellent mass flow rate
    double thrust = 1900000.0; // Netwons
    double T = thrust;

    double hturn = 1000 ;// pitchover height why is this known ???

    // Atmospheric Conditions
    double CD = 0.3; // drag coefficient
    double g0 = 9.81;
    double g = 0.0;
    double rho0 = 12.93; // density of air kg/m^3 at sea level
    long double rho = 0.0; // initialising
    double hscale = 8500.0; // m, scale of rapid atmospheric change within Earths Atmosphere

    //size parameters
    double diam = 3.05; // m rocket diameter
    double A = M_PI * pow((diam / 2), 2); // area of the nose (frontal portion)

    // differential inputs
    double t = 0.0; // seconds
    double v = 0.0; // meters per second
    double x = 0.0; // horizontal position in meters


    // angular parameters
    double psi_dot = 0.0;
    double theta_dot = 0.0;
    double h_dot = 0.0;
    double phi_dot = 0.0;
    double phi = 0.0;
    double psi = 10 *  M_PI / 180; //M_PI * 5/6; //10 * M_PI/180; //start at 10 degrees // cant be zero because then it never changes
    double theta = 0.0;
    double h = 0.0;
    double Re = 6371000.0; // radius of earth from centre to surface
    double drag = 0.0;

    //solving using Eulers
    double dt = 1; // time step in seconds
    double t_end = 1400.0; // end time for the simulation

    // Open a file to write the data
    std::ofstream outfile("rocket_trajectory.csv");
    outfile << "Time,Velocity,Distance,Height\n"; // Header

    // Euler's method loop
    while (t <= t_end) {

        //------------------------------Solving Vertical and Angular components-------------------------------------

        //calculating Drag and gravity
        rho = rho0 * exp(-h / hscale);
        g = g0 / pow((1 + h / Re), 2);
        drag = (1 / 2) * rho * pow(v, 2) * A * CD; // Vertical component of drag     DRAG IS ALWAYS ZERO

        //------------------------------Determine the thrust of rocket -------------------------------------
        // Check if the fuel is depleted and set Thrust to 0
       // T = (t <= tburn)? thrust : 0;

       if (t<=tburn){
           m = m0 - mDot * t; //Rocket mass
            T=thrust;
       }else{
           m = mstruc; //Rocket mass
           T = 0;
       }

        //------------------------------Solving De's-------------------------------------
        double Dv_dt;
        if (h<=100000) { // before gravity turn

            //m = m0 - mDot * t; //Rocket mass     moving this around did basically nothing so could return to regular one line if statement
            psi_dot = 0; // rate of change of psi
            theta_dot = 0; // rate of change of angle relative to Earths centre
            Dv_dt = T/m - drag/m - g ;
            h_dot = v; // change in height


        } else { // after gravity turn
           // m = mstruc; //Rocket mass
            phi_dot = g * sin(psi) / v;
            Dv_dt = T/m + (drag)/m  - (g * cos(psi)); // acceleration
            h_dot = v * cos(psi);
            theta_dot = (v * sin(psi)) / (Re + h);
            psi_dot = phi_dot - theta_dot;

        }

        // Update the velocity and time using Euler's method

        v = v + Dv_dt * dt;
        //h = h + v * dt;
        h = h + h_dot * dt; // height as well but calculated differently
        psi = psi + psi_dot * dt; // angle
        phi_dot = phi + phi_dot * dt;
        theta = theta + theta_dot * dt; //

        t = t + dt;

        ///When Rocket collides
        if (h<0){
            h = 0;
            break;
        }

        double dr = theta * Re/1000;
        //printf("%f",dr);
        // Write current time, velocity, and position to file
        outfile << t << "," << v << "," << dr << "," << h << "," << psi*180/M_PI <<"\n";

        // Output the current time, velocity, and position
        std::cout << "Time: " << t << " s, Velocity: " << v <<" Height: " << h << " m" << " Acceleration: " << Dv_dt << " angle: " << psi*180/M_PI <<std::endl;
    }

    //plotting
    // Plot the results using matplotlibcpp
    outfile.close();
    std::cout << "Data written to rocket_trajectory.csv" << std::endl;

    return 0;
}
