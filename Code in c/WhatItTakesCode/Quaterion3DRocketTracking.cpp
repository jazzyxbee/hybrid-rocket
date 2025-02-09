//
// Created by User on 09/01/2025.
//

#include <iostream>
#include <cmath>
#include <fstream>
#include <array>
#include <functional>

#include "Quaterion3DRocketTracking.h"

typedef struct { double w, x, y, z; } Quaternion;
typedef struct { double x, y, z; } Vector3;

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

// Update position using forces and acceleration
Vector3 update_position(Vector3 position, Vector3 velocity, Vector3 acceleration, double dt) {
    velocity.x += acceleration.x * dt;
    velocity.y += acceleration.y * dt;
    velocity.z += acceleration.z * dt;

    position.x += velocity.x * dt;
    position.y += velocity.y * dt;
    position.z += velocity.z * dt;

    return position;
}

Quaternion quaternion_derivative(Quaternion q, Quaternion omega_q) {
    return quaternion_scale(quaternion_multiply(q, omega_q), 0.5);
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

// eulers method update redundant because of RK4
Quaternion quaternion_update(Quaternion q, Quaternion omega_q, double dt) {
    Quaternion delta_q = quaternion_multiply(q, omega_q);
    delta_q = quaternion_scale(delta_q, 0.5 * dt);
    q = quaternion_add(q, delta_q);
    return quaternion_normalize(q);  // Normalize to maintain unit quaternion
}


int main() {
    Quaternion orientation = {1, 0, 0, 0}; // Initial orientation
    Vector3 position = {0, 0, 0};         // Initial position
    Vector3 velocity = {1, 1, 1};         // Initial velocity
    Vector3 thrust_local = {0, 0, 10};          // Thrust in local frame
    double mass = 1000.0;                 // Rocket mass
    double dt = 0.01;                     // Time step
    double t = 0, tf = 10.0;              // Simulation time
    double omegaX = 0;
    double omegaY = 0;
    double omegaZ = 0;
    Quaternion w = {0,omegaX,omegaY,omegaZ};

    while (t < tf) {
        // Update orientation with angular velocities
        Quaternion orientation_dot = rk4_quaternion_update(orientation,w,dt);
        //printf("Ori: %f ,  %f  , %f,   %f" , orientation.x,orientation.y,orientation.z,   orientation.w);

        // Rotate thrust vector to global frame
        Vector3 thrust_global = rotate_vector(orientation_dot, thrust_local.x,  thrust_local.y, thrust_local.z);

        // Compute acceleration
        Vector3 acceleration = {thrust_global.x / mass, thrust_global.y / mass, thrust_global.z / mass};

        //update_position();
        position = update_position(position, velocity, acceleration, dt);

        // Update position and velocity

        printf("Time: %.2f, Position: (%.2f, %.2f, %.2f)\n", t, position.x, position.y, position.z);
        t += dt;
    }

    return 0;
}