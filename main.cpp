// main.cpp - 1D GLM-MHD with proper By preservation
#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <iomanip>
#include <algorithm>
#include "physics.hpp"

static constexpr double CFL = 0.3;
static constexpr double CR = 0.18;

// State vector indices
enum { IRHO=0, IMOMX, IMOMY, IMOMZ, IENER, IBX, IBY, IBZ, IPSI, NVAR };

// Primitive variables from conservative
void cons2prim(const double U[NVAR], double& rho, double& u, double& v, double& w,
               double& p, double& Bx, double& By, double& Bz, double& psi) {
    rho = U[IRHO];
    u = U[IMOMX] / rho;
    v = U[IMOMY] / rho;
    w = U[IMOMZ] / rho;
    Bx = U[IBX];
    By = U[IBY];
    Bz = U[IBZ];
    psi = U[IPSI];
    
    double v2 = u*u + v*v + w*w;
    double B2 = Bx*Bx + By*By + Bz*Bz;
    double psi2 = psi*psi;
    p = (GAMMA - 1.0) * (U[IENER] - 0.5*rho*v2 - 0.5*B2 - 0.5*psi2);
    p = std::max(p, 1e-10);
}

// Conservative variables from primitive
void prim2cons(double rho, double u, double v, double w, double p,
               double Bx, double By, double Bz, double psi, double U[NVAR]) {
    U[IRHO] = rho;
    U[IMOMX] = rho * u;
    U[IMOMY] = rho * v;
    U[IMOMZ] = rho * w;
    U[IBX] = Bx;
    U[IBY] = By;
    U[IBZ] = Bz;
    U[IPSI] = psi;
    
    double v2 = u*u + v*v + w*w;
    double B2 = Bx*Bx + By*By + Bz*Bz;
    U[IENER] = p/(GAMMA-1.0) + 0.5*rho*v2 + 0.5*B2 + 0.5*psi*psi;
}

// Compute flux vector
void compute_flux(double rho, double u, double v, double w, double p,
                  double Bx, double By, double Bz, double psi, double ch,
                  double F[NVAR]) {
    double B2 = Bx*Bx + By*By + Bz*Bz;
    double ptot = p + 0.5*(B2 + psi*psi);
    
    F[IRHO] = rho * u;
    F[IMOMX] = rho*u*u + ptot - Bx*Bx;
    F[IMOMY] = rho*u*v - Bx*By;
    F[IMOMZ] = rho*u*w - Bx*Bz;
    
    double E = p/(GAMMA-1.0) + 0.5*rho*(u*u + v*v + w*w) + 0.5*B2
               + 0.5*psi*psi;
    F[IENER] = u*(E + ptot) - Bx*(u*Bx + v*By + w*Bz) + psi*Bx;
    
    F[IBX] = psi;  // GLM flux
    F[IBY] = u*By - v*Bx;
    F[IBZ] = u*Bz - w*Bx;
    F[IPSI] = ch*ch * Bx;  // GLM flux
}

// Rusanov (Lax-Friedrichs) flux - more diffusive but more stable
void rusanov_flux(const double UL[NVAR], const double UR[NVAR], double ch,
                  double F_interface[NVAR]) {
    // Left state
    double rhoL, uL, vL, wL, pL, BxL, ByL, BzL, psiL;
    cons2prim(UL, rhoL, uL, vL, wL, pL, BxL, ByL, BzL, psiL);
    
    // Right state
    double rhoR, uR, vR, wR, pR, BxR, ByR, BzR, psiR;
    cons2prim(UR, rhoR, uR, vR, wR, pR, BxR, ByR, BzR, psiR);
    
    // Wave speeds
    double a2L = GAMMA * pL / rhoL;
    double a2R = GAMMA * pR / rhoR;
    double CA2L = (BxL*BxL + ByL*ByL + BzL*BzL) / rhoL;
    double CA2R = (BxR*BxR + ByR*ByR + BzR*BzR) / rhoR;
    
    double cfL = std::sqrt(0.5*(a2L + CA2L + 
                 std::sqrt(std::max(0.0, (a2L+CA2L)*(a2L+CA2L) - 4*a2L*BxL*BxL/rhoL))));
    double cfR = std::sqrt(0.5*(a2R + CA2R + 
                 std::sqrt(std::max(0.0, (a2R+CA2R)*(a2R+CA2R) - 4*a2R*BxR*BxR/rhoR))));
    
    // Maximum wave speed (including GLM waves)
    double SL = std::min(uL - cfL, uR - cfR);
    double SR = std::max(uL + cfL, uR + cfR);
    double Smax = std::max({std::abs(SL), std::abs(SR), ch});
    
    // Left and right fluxes
    double FL[NVAR], FR[NVAR];
    compute_flux(rhoL, uL, vL, wL, pL, BxL, ByL, BzL, psiL, ch, FL);
    compute_flux(rhoR, uR, vR, wR, pR, BxR, ByR, BzR, psiR, ch, FR);
    
    // Rusanov flux
    for (int k = 0; k < NVAR; ++k) {
        F_interface[k] = 0.5 * (FL[k] + FR[k] - Smax * (UR[k] - UL[k]));
    }
}

// Write output
void write_output(const std::vector<double>& U, int nx, double x0, double dx,
                  const std::string& filename) {
    std::ofstream f(filename);
    f << std::scientific << std::setprecision(8);
    f << "x,rho,u,v,w,p,bx,by,bz,psi,divB\n";
    
    for (int i = 1; i < nx+1; ++i) {
        double x = x0 + (i-0.5)*dx;
        
        double rho, u, v, w, p, Bx, By, Bz, psi;
        cons2prim(&U[NVAR*i], rho, u, v, w, p, Bx, By, Bz, psi);
        
        // Compute div B
        double divB = 0.0;
        if (i > 1 && i < nx) {
            divB = (U[NVAR*(i+1) + IBX] - U[NVAR*(i-1) + IBX]) / (2.0*dx);
        }
        
        f << x << ',' << rho << ',' << u << ',' << v << ',' << w << ',' 
          << p << ',' << Bx << ',' << By << ',' << Bz << ',' << psi << ',' 
          << divB << '\n';
    }
}

int main() {
    // Problem parameters
    const int nx = 800;
    const double x0 = -0.5, x1 = 0.5;
    const double dx = (x1 - x0) / nx;
    const double t_end = 0.2;
    
    std::cout << "1D GLM-MHD solver (Dedner et al. 2002)" << std::endl;
    std::cout << "Method: Mixed GLM with Rusanov flux" << std::endl;
    std::cout << "Grid: nx = " << nx << ", dx = " << dx << std::endl;
    
    // Create output directory
    std::filesystem::create_directory("Result");
    
    // Allocate solution array (including ghost cells)
    std::vector<double> U((nx+2)*NVAR, 0.0);
    
    // Initialize
    std::vector<Cell> W(nx+2);
    initialize_riemann_problem(W, x0, dx);
    
    // Convert to conservative variables
    for (int i = 0; i < nx+2; ++i) {
        prim2cons(W[i].rho, W[i].u, W[i].v, W[i].w, W[i].p,
                  W[i].bx, W[i].by, W[i].bz, W[i].psi, &U[NVAR*i]);
    }
    
    // Time evolution
    double t = 0;
    int step = 0;
    int output_count = 0;
    
    // Initial output
    write_output(U, nx, x0, dx, "Result/out_state_0.csv");
    
    // Check initial By
    double By_initial = W[nx/2].by;
    std::cout << "Initial By = " << By_initial << std::endl;
    
    while (t < t_end) {
        // Boundary conditions (transmissive)
        for (int k = 0; k < NVAR; ++k) {
            U[k] = U[NVAR + k];
            U[NVAR*(nx+1) + k] = U[NVAR*nx + k];
        }
        
        // Maximum wave speed
        double max_speed = 0;
        for (int i = 1; i <= nx; ++i) {
            double rho, u, v, w, p, Bx, By, Bz, psi;
            cons2prim(&U[NVAR*i], rho, u, v, w, p, Bx, By, Bz, psi);
            
            double a2 = GAMMA * p / rho;
            double CA2 = (Bx*Bx + By*By + Bz*Bz) / rho;
            double cf2 = 0.5*(a2 + CA2 + std::sqrt(std::max(0.0, 
                         (a2+CA2)*(a2+CA2) - 4*a2*Bx*Bx/rho)));
            double cf = std::sqrt(cf2);
            max_speed = std::max(max_speed, std::abs(u) + cf);
        }
        
        // GLM wave speed
        double ch = std::max(max_speed, 1.0) * 1.1;
        double cp = ch / std::sqrt(CR);
        
        // Time step
        double dt = CFL * dx / std::max(max_speed, ch);
        if (t + dt > t_end) dt = t_end - t;
        
        // Copy current solution
        std::vector<double> U_old = U;
        
        // Update using Rusanov flux
        for (int i = 1; i <= nx; ++i) {
            double F_left[NVAR], F_right[NVAR];
            
            // Left interface flux
            rusanov_flux(&U_old[NVAR*(i-1)], &U_old[NVAR*i], ch, F_left);
            
            // Right interface flux
            rusanov_flux(&U_old[NVAR*i], &U_old[NVAR*(i+1)], ch, F_right);
            
            // Update
            for (int k = 0; k < NVAR; ++k) {
                U[NVAR*i + k] = U_old[NVAR*i + k] - dt/dx * (F_right[k] - F_left[k]);
            }
            
            // Apply mixed GLM correction to psi
            U[NVAR*i + IPSI] *= std::exp(-dt * ch*ch / (cp*cp));
        }
        
        // Update time
        t += dt;
        step++;
        
        // Output and diagnostics
        if (step % 50 == 0 || t >= t_end) {
            output_count++;
            write_output(U, nx, x0, dx, 
                        "Result/out_state_" + std::to_string(output_count*50) + ".csv");
            
            // Check By conservation
            double By_min = 1e10, By_max = -1e10;
            for (int i = 1; i <= nx; ++i) {
                double rho, u, v, w, p, Bx, By, Bz, psi;
                cons2prim(&U[NVAR*i], rho, u, v, w, p, Bx, By, Bz, psi);
                By_min = std::min(By_min, By);
                By_max = std::max(By_max, By);
            }
            
            std::cout << "Step " << step << ", t = " << std::setprecision(4) << t
                     << ", By: [" << By_min << ", " << By_max << "]"
                     << ", variation = " << std::scientific << (By_max - By_min) << std::endl;
        }
    }
    
    std::cout << "\nSimulation complete!" << std::endl;
    return 0;
}
