#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <iomanip>

static constexpr double GAMMA = 5.0/3.0;
static constexpr double CFL   = 0.3;
static constexpr double CR    = 0.18; // c_p^2 / c_h^2

struct Cell {
    double rho, u, p, e, bx, by, psi;
};

struct Flux {
    double rho, momx, E, Bx, By, psi;
};

// Fast magnetosonic speed
static double fast_speed(double rho,double p,double bx,double by){
    double a2 = GAMMA * p / rho;
    double b2 = (bx*bx + by*by)/rho;
    double bx2 = (bx*bx)/rho;
    double term = std::sqrt(std::max(0.0,(a2 + b2)*(a2 + b2) - 4*a2*bx2));
    return std::sqrt(0.5*(a2 + b2 + term));
}

static Flux hll_flux(const Cell& L,const Cell& R,double ch){
    double B2L = L.bx*L.bx + L.by*L.by;
    double B2R = R.bx*R.bx + R.by*R.by;
    double ptL = L.p + 0.5*B2L;
    double ptR = R.p + 0.5*B2R;
    double EL = L.p/(GAMMA-1.0) + 0.5*L.rho*L.u*L.u + 0.5*B2L;
    double ER = R.p/(GAMMA-1.0) + 0.5*R.rho*R.u*R.u + 0.5*B2R;

    double cfL = fast_speed(L.rho,L.p,L.bx,L.by);
    double cfR = fast_speed(R.rho,R.p,R.bx,R.by);
    double SL = std::min(L.u - cfL, R.u - cfR);
    double SR = std::max(L.u + cfL, R.u + cfR);

    Flux F;
    if(SL>0){
        F.rho = L.rho * L.u;
        F.momx = L.rho*L.u*L.u + ptL - L.bx*L.bx;
        F.E   = (EL + ptL)*L.u - L.bx*(L.u*L.bx);
        F.Bx  = L.psi;
        F.By  = L.u*L.by - 0.0*L.bx; // v=0
        F.psi = ch*ch*L.bx;
    }else if(SR<0){
        F.rho = R.rho * R.u;
        F.momx = R.rho*R.u*R.u + ptR - R.bx*R.bx;
        F.E   = (ER + ptR)*R.u - R.bx*(R.u*R.bx);
        F.Bx  = R.psi;
        F.By  = R.u*R.by - 0.0*R.bx;
        F.psi = ch*ch*R.bx;
    }else{
        double FL_rho = L.rho*L.u;
        double FR_rho = R.rho*R.u;
        double FL_momx = L.rho*L.u*L.u + ptL - L.bx*L.bx;
        double FR_momx = R.rho*R.u*R.u + ptR - R.bx*R.bx;
        double FL_E = (EL + ptL)*L.u - L.bx*(L.u*L.bx);
        double FR_E = (ER + ptR)*R.u - R.bx*(R.u*R.bx);
        double FL_By = L.u*L.by;
        double FR_By = R.u*R.by;
        F.rho = (SR*FL_rho - SL*FR_rho + SL*SR*(R.rho - L.rho))/(SR-SL);
        F.momx = (SR*FL_momx - SL*FR_momx + SL*SR*(R.rho*R.u - L.rho*L.u))/(SR-SL);
        F.E   = (SR*FL_E - SL*FR_E + SL*SR*(ER-EL))/(SR-SL);
        F.Bx  = (SR*L.psi - SL*R.psi + SL*SR*(R.bx - L.bx))/(SR-SL);
        F.By  = (SR*FL_By - SL*FR_By + SL*SR*(R.by - L.by))/(SR-SL);
        F.psi = ch*ch*(SR*L.bx - SL*R.bx + SL*SR*(R.psi - L.psi))/(SR-SL);
    }
    return F;
}

static void write_field(const std::vector<Cell>& U,double x0,double dx,const std::string& prefix,int step){
    std::ofstream f(prefix+std::to_string(step)+".csv");
    for(size_t i=1;i<U.size()-1;++i){
        double x = x0 + (i-1+0.5)*dx;
        f<<x<<','<<U[i].rho<<','<<U[i].u<<','<<U[i].p<<','<<U[i].bx<<','<<U[i].by<<','<<U[i].psi<<'\n';
    }
}

int main(){
    const int nx=1024;
    const double x0=-0.5,x1=0.5;
    const double dx=(x1-x0)/nx;
    const double t_end=0.25;
    const int output_every=20;
    std::filesystem::create_directory("Result");

    std::vector<Cell> U(nx+2); // ghost cells
    for(int i=0;i<nx+2;++i){
        double x = x0 + (i-1+0.5)*dx;
        Cell c;
        if(x<0){
            c.rho=1.0;
            c.u=10.0;
            c.by=5.0/std::sqrt(4*M_PI);
            c.bx=5.0/std::sqrt(4*M_PI);
            c.p=20.0;
        }else{
            c.rho=1.0;
            c.u=-10.0;
            c.by=5.0/std::sqrt(4*M_PI);
            c.bx=5.0/std::sqrt(4*M_PI);
            c.p=0.0;
        }
        c.e = c.p/(GAMMA-1.0) + 0.5*c.rho*c.u*c.u + 0.5*(c.bx*c.bx + c.by*c.by);
        c.psi=0.0;
        U[i]=c;
    }

    auto apply_bc=[&](){
        U[0]=U[1];
        U[nx+1]=U[nx];
    };
    apply_bc();

    double t=0; int step=0;
    while(t<t_end-1e-12){
        double max_speed=0.0;
        for(int i=1;i<=nx;++i){
            double cf = fast_speed(U[i].rho,U[i].p,U[i].bx,U[i].by);
            max_speed = std::max(max_speed, std::abs(U[i].u)+cf);
        }
        double ch = 1.5*max_speed;
        max_speed = std::max(max_speed, ch);
        double dt = CFL*dx/max_speed;
        if(t+dt>t_end) dt=t_end-t;

        std::vector<Flux> F(nx+1);
        for(int i=0;i<=nx;i++){
            F[i]=hll_flux(U[i],U[i+1],ch);
        }
        std::vector<Cell> Un(nx+2);
        for(int i=1;i<=nx;++i){
            Un[i].rho = U[i].rho - dt/dx*(F[i].rho - F[i-1].rho);
            double mom = U[i].rho*U[i].u - dt/dx*(F[i].momx - F[i-1].momx);
            Un[i].bx  = U[i].bx - dt/dx*(F[i].Bx - F[i-1].Bx);
            Un[i].by  = U[i].by - dt/dx*(F[i].By - F[i-1].By);
            Un[i].psi = U[i].psi - dt/dx*(F[i].psi - F[i-1].psi);
            double E  = U[i].e  - dt/dx*(F[i].E   - F[i-1].E);
            Un[i].u   = mom/Un[i].rho;
            Un[i].e   = E;
        }
        U.swap(Un);
        apply_bc();

        // GLM divergence cleaning
        std::vector<Cell> Uc = U;
        for(int i=1;i<=nx;++i){
            double divB = (U[i+1].bx - U[i-1].bx)/(2*dx);
            Uc[i].psi = U[i].psi - dt*(ch*ch*divB + (1.0/CR)*U[i].psi);
        }
        for(int i=1;i<=nx;++i){
            double dpsi_dx = (Uc[i+1].psi - Uc[i-1].psi)/(2*dx);
            Uc[i].bx = U[i].bx - dt*dpsi_dx;
        }
        U.swap(Uc);
        apply_bc();
        // update pressure from energy
        for(int i=1;i<=nx;++i){
            double B2 = U[i].bx*U[i].bx + U[i].by*U[i].by;
            double kin = 0.5*U[i].rho*U[i].u*U[i].u;
            U[i].p = (U[i].e - kin - 0.5*B2)*(GAMMA-1.0);
            if(U[i].p < 1e-8) U[i].p = 1e-8;
        }

        t += dt; step++;
        if(step%output_every==0 || t>=t_end-1e-12){
            write_field(U,x0,dx,"Result/out_state_",step);
            std::cout<<"step "<<step<<" t="<<t<<"\n";
        }
    }
    return 0;
}
