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

std::array<double,4>derivitives(double t,double y0,double y1,double y2,double y3,double y4){
        double diam=0;//rocketdiameter
        double A=M_PI/(4*pow(diam,2));//m^2frontalarea
        double h=0;
        double g0=9.81;
        double rho= 1.293;
        double rho0=0;
        double Re=6371000;
        double D=0;
        double CD=0.5;
        double hscale=1;
        double hturn=30;
        double psi_dot=0;
        double v_dot=0;
        double theta_dot=0;
        double phi_dot=0;
        double h_dot=0;
        double tburn=2.5;
        double thrust=15;
        double m=0;

        double v = y0;
        double psi = y1;
        double theta = y2;
            h=y3;
            //determinegravityanddrag
            double g=g0/pow((1+h/Re),2);
            rho=rho0*exp(-h/hscale);
            D=0.5*rho*pow(v,2)*A*CD;

            if(tburn<t){
           //Computethederivativeofvattimet
                double m = 0.05-0.002*t;//Rocketmass

           //thrust=thrust
            }else{
                m = 0.05-0.002 * tburn;//Rocketmassoncefullhasallbeenburnt
                thrust=0;
            }

            if(h<=hturn){
                psi_dot=0;
                v_dot=thrust/m-D/m-g;
            theta_dot=0;
            h_dot=v;
        }else{
            phi_dot=g*sin(psi)/v;
            v_dot=thrust/m-D/m-g*cos(psi);
            h_dot=v*cos(psi);
            theta_dot=v*sin(psi)/(Re+h);
            psi_dot=phi_dot-theta_dot;
        }

        return{v_dot,psi_dot,theta_dot,h_dot};
}

int main() {
    // Initial parameters
    // Using example values from https://www.youtube.com/watch?v=NFhiJWbNU1g&t=1856s
    double m = 0;
    double m0 = 0.05; //inital mass kg
    double mDot = 0.002;  // kg/sec mass flow rate
    double fuelAt1Sec = 0.048; //kg
    double thrust = 15; // Netwons
    double K = 0.02; // drag co=efficient
    double g = -9.81;
    double tburn = 1;
    double rho= 1.293; // density of air kg/m^3
    double radius = 0.5; // meters
    double A  = M_PI * pow(radius,2); // area of the nose (frontal portion)

    double t = 0; // seconds
    double v = 0; // meters per second
    double y = 0; // vertical position in meters
    double pitch = 90; // starts at 90 degrees, straight up.
    double pitchRate = -9.8; // gravity

    //solving using Eulers
    double dt = 0.01; // time step in seconds
    double t_end = 2.51; // end time for the simulation

    // Open a file to write the data
    std::ofstream outfile("rocket_trajectory.csv");
    outfile << "Time,Velocity,Position,PitchAngle\n"; // Header

    // Euler's method loop
    while (t <= t_end) {

        if(t<=tburn){
            // Compute the derivative of v at time t
            m = m0 - mDot * t; //Rocket mass

            // Check if Mr is very small or zero to avoid division by a small number
            if (m <= 0.0001) {
                std::cerr << "Mr is too small or zero at t = " << t << std::endl;
                break;
            }

        } else {

            // Check if the fuel is depleted and set Thrust to 0
            m = m0 - mDot * tburn; //Rocket mass
            thrust = 0;

        }

        //calculating Drag
        double D = 0.5 * rho * pow(v,2) * A * K;

        double Dv_dt, Dv_dtdrag;
        // Compute the derivative of v at time t
        if (thrust > 0){ // before gravity turn
            Dv_dt = (thrust - m * 9.8 - D - v * mDot) / m; // change to drag ... K * pow(v, 2)
            Dv_dtdrag = (thrust - m * 9.8 - D ) / m;

            printf("no drag %f drag %f \n", Dv_dt, Dv_dtdrag);
        } else{ // after gravity turn
            Dv_dt = (-9.8 - (K / m) * pow(v, 2));
            Dv_dtdrag = (-9.8 - (D / m) );
        }

        // Update the velocity and time using Euler's method
        v = v + Dv_dt * dt; // velocity
        y = y + v * dt; // height

        // Check if the position has reached or gone below 0
        if (y <= 0) {
            v = 0;
            y = 0;
            Dv_dt = 0;
        }

        t = t + dt;

        // Update pitch angle based on the velocity and height
        if (y > 0) {
            pitch = 90; // Gravity turn approximation but it needs proper representation
        }

        // Write current time, velocity, and position to file
        outfile << t << "," << v << "," << y << "\n";

        // Output the current time, velocity, and position
        //std::cout << "Time: " << t << " s, Velocity: " << v <<" Height: " << y << " m" << " Acceleration: " << Dv_dt << " Pitch Angle: " << pitch <<std::endl;
    }

    //plotting
    // Plot the results using matplotlibcpp
    outfile.close();
    std::cout << "Data written to rocket_trajectory.csv" << std::endl;

    return 0;
}
