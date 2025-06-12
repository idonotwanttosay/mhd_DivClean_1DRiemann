#include "physics.hpp"
#include <cmath>
#include <iostream>
void initialize_riemann_problem(std::vector<Cell>& U, double x0, double dx) {
    // 1D Riemann problem - checking multiple sources for correct values
    // This appears to be a shock tube problem with magnetic field
    
    const double B0 = 5.0 / std::sqrt(4.0 * M_PI);  // This gives ~1.41
    
    for (size_t i = 0; i < U.size(); ++i) {
        double x = x0 + (static_cast<double>(i) - 1.0 + 0.5) * dx;
        Cell c;
        
        if (x < 0.0) {
            // Left state
            c.rho = 1.0;
            c.u = 10.0;   // High velocity to the right
            c.v = 0.0;    // No y-velocity initially
            c.w = 0.0;
            c.p = 20.0;   // High pressure
            c.bx = B0;    // Constant Bx = 5/√(4π)
            c.by = B0;    // Constant By = 5/√(4π)
            c.bz = 0.0;
        } else {
            // Right state  
            c.rho = 1.0;
            c.u = -10.0;  // High velocity to the left
            c.v = 0.0;    // No y-velocity initially
            c.w = 0.0;
            c.p = 1.0;    // Low pressure
            c.bx = B0;    // Same Bx
            c.by = B0;    // Same By
            c.bz = 0.0;
        }
        
        // Initialize psi to zero
        c.psi = 0.0;
        
        // Compute total energy
        double B2 = c.bx*c.bx + c.by*c.by + c.bz*c.bz;
        double v2 = c.u*c.u + c.v*c.v + c.w*c.w;
        c.e = c.p / (GAMMA - 1.0) + 0.5 * c.rho * v2 + 0.5 * B2;
        
        U[i] = c;
    }
    
    // Print initial conditions for verification
    static bool printed = false;
    if (!printed) {
        printed = true;
        std::cout << "\nInitial conditions:" << std::endl;
        std::cout << "Left state (x < 0):" << std::endl;
        std::cout << "  rho = 1.0, u = 10.0, v = 0.0, p = 20.0" << std::endl;
        std::cout << "  Bx = " << B0 << ", By = " << B0 << ", Bz = 0.0" << std::endl;
        std::cout << "Right state (x > 0):" << std::endl;
        std::cout << "  rho = 1.0, u = -10.0, v = 0.0, p = 1.0" << std::endl;
        std::cout << "  Bx = " << B0 << ", By = " << B0 << ", Bz = 0.0" << std::endl;
        std::cout << std::endl;
    }
}
