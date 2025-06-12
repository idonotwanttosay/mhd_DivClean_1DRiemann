#pragma once
#include <vector>

// Cell structure for 1D MHD with GLM divergence cleaning
struct Cell {
    // Primitive variables
    double rho;     // density
    double u, v, w; // velocity components
    double p;       // pressure
    double bx, by, bz; // magnetic field components
    double psi;     // GLM scalar field
    
    // Conserved variable (total energy)
    double e;       // total energy
};

static constexpr double GAMMA = 5.0/3.0;

void initialize_riemann_problem(std::vector<Cell>& U, double x0, double dx);
