#pragma once
#include <vector>

struct Cell {
    double rho, u, p, e, bx, by, psi;
};

static constexpr double GAMMA = 5.0/3.0;

void initialize_riemann_problem(std::vector<Cell>& U, double x0, double dx);
