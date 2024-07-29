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
    double T = 0.0;
    double m = 0.0; // initial mass
    double mprop = 111130.0; // kg propellent mass
    double mpl = 32000.0; // kg payload mass
    double mstruc = 6736.0; // kg structure mass
    double m0 = mprop + mpl + mstruc; // initialmass
    double tburn = 356.0; // burn time (s)
    double mDot =  mprop/tburn; // kg/s propellent mass flow rate
    double thrust = 1900000.0; // Netwons

    //double hturn = 1000 // pitchover height why is this known ???

    // Atmospheric Conditions
    double CD = 0.3; // drag coefficient
    double g0 = -9.81;
    double g = 0.0;
    double rho0 = 12.93; // density of air kg/m^3 at sea level
    long double rho = 0.0; // initialising
    double hscale = 8500.0; // m, scale of rapid atmospheric change within Earths Atmosphere

    //size parameters
    double diam = 3.05; // m rocket diameter
    double A  = M_PI* pow((diam/2),2); // area of the nose (frontal portion)

    // differential inputs
    double t = 0.0; // seconds
    double v = 0.0; // meters per second
    double y = 0.0; // vertical position in meters
    double x = 0.0; // horizontal position in meters


    // angular parameters
    double psi_dot = 0.0;
    double theta_dot = 0.0;
    double h_dot = 0.0;
    double phi_dot = 0.0;
    double phi = 0.0;
    double deg = 0.0;
    double psi0 = 0.3*deg; // rad initial flight path angle
    double theta0 = 0.0; // rad inital downrange angle
    double h0 = 0.0; // initial height
    double psi = M_PI*5/6;
    double theta = 0.0;
    double h = 0.0;
    double Re = 6371000.0; // radius of earth from centre to surface

    //solving using Eulers
    double dt = 2; // time step in seconds
    double t_end = 1600.0; // end time for the simulation

    // Open a file to write the data
    std::ofstream outfile("rocket_trajectory.csv");
    outfile << "Time,Velocity,Distance,Height\n"; // Header


    // Euler's method loop
    while (t <= t_end) {
        //------------------------------Determine the thrust of rocket -------------------------------------
        // Check if the fuel is depleted and set Thrust to 0
        T = (t <= tburn) ? thrust : 0;

        //------------------------------Solving Vertical and Angular components-------------------------------------

        //why do this here why not just calculate individually
        // Calculate thrust components
        double thrustX = T * cos(psi); // Horizontal component of thrust
        double thrustY = T * sin(psi); // Vertical component of thrust

        // Compute drag force components
        double vx = v * cos(psi); // Horizontal velocity component
        double vy = v * sin(psi); // Vertical velocity component

        //calculating Drag and gravity
        rho = rho0 * exp(-h/hscale);
        printf("%Lf \n",rho);
        g = g0 /  pow((1 + h / Re),2);

        long double dragX = 0.5 * rho * pow(vx,2) * A * CD; // Horizontal component of drag
        long double dragY = 0.5 * rho * pow(vy,2) * A * CD; // Vertical component of drag

        //------------------------------Solving De's-------------------------------------
        double Dv_dt, Dv_dtX, Dv_dtY, Dv_dtdrag;
        if (T > 0){ // before gravity turn

            m = m0 - mDot * t; //Rocket mass
            psi_dot = 0; // rate of change of psi
            Dv_dtY = ((T - dragY) / m) + g; // acceleration
            //Dv_dtX = ((thrustX - dragX) / m); // acceleration
            theta_dot = 0; // rate of change of angle relative to Earths centre
            h_dot = v; // change in height

            //printf("no drag %f drag %f \n", Dv_dt, Dv_dtdrag);
            //Dv_dt = (thrust - m * 9.8 - CD * pow(v,2)- v * mDot) / m; // change to drag ... K * pow(v, 2   original DE

        } else{ // after gravity turn

            m = m0 - mDot * tburn; //Rocket mass

            phi_dot = g * sin(psi) / v;  //
            Dv_dtY = ((thrustY - dragY) / m) - g * cos(psi);
            //Dv_dtX = ((thrustX - dragX) / m); // acceleration
            h_dot = -v * cos(psi);
            theta_dot = v * sin(psi) / (Re + h);
            psi_dot = phi_dot - theta_dot;

            //Dv_dt = ((T-dragY) / m - 9.81); // original DE
        }

        // Update the velocity and time using Euler's method
        v = v + Dv_dtY * dt; // velocity //   v = v + (Dv_dtY + Dv_dtX) * dt; // velocity
        x = x + Dv_dtX * dt * dt; // distance
        y = y + Dv_dtY * dt; // height
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
        // Write current time, velocity, and position to file
        outfile << t << "," << v << "," << x << "," << h << "\n";

        // Output the current time, velocity, and position
        std::cout << "Time: " << t << " s, Velocity: " << v <<" Height: " << h << " m" << " Acceleration: " << Dv_dtY << " angle: " << psi*180/M_PI <<std::endl;
    }

    //plotting
    // Plot the results using matplotlibcpp
    outfile.close();
    std::cout << "Data written to rocket_trajectory.csv" << std::endl;

    return 0;
}
