#include "physics.hpp"
#include <cmath>

void initialize_riemann_problem(std::vector<Cell>& U, double x0, double dx) {
    for(size_t i = 0; i < U.size(); ++i) {
        double x = x0 + (static_cast<double>(i)-1.0+0.5)*dx; // ghost cells assumed
        Cell c;
        if (x < 0.0) {
            c.rho = 1.0;
            c.u   = 10.0;
            c.by  = 5.0 / std::sqrt(4.0 * M_PI);
            c.bx  = 5.0 / std::sqrt(4.0 * M_PI);
            c.p   = 20.0;
        } else {
            c.rho = 1.0;
            c.u   = -10.0;
            c.by  = 5.0 / std::sqrt(4.0 * M_PI);
            c.bx  = 5.0 / std::sqrt(4.0 * M_PI);
            c.p   = 0.0;
        }
        c.e   = c.p/(GAMMA-1.0) + 0.5*c.rho*c.u*c.u + 0.5*(c.bx*c.bx + c.by*c.by);
        c.psi = 0.0;
        U[i]  = c;
    }
}
