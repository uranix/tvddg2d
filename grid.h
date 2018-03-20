#pragma once

#include <vector>

template<class param>
struct grid {
    const int Nx, Ny;
    const double Lx, Ly;
    const double hx, hy;

    std::vector<param> data;

    grid(const int Nx, const int Ny, const double Lx, const double Ly)
        : Nx(Nx), Ny(Ny), Lx(Lx), Ly(Ly), hx(Lx/Nx), hy(Ly/Ny), data(Nx * Ny)
    { }

    const param operator()(int i, int j) const {
        return data[i * Ny + j];
    }

    param &operator()(int i, int j) {
        return data[i * Ny + j];
    }

    template<class PROBLEM>
    void fill(const PROBLEM &prob) {
        for (int i = 0; i < Nx; i++)
            for (int j = 0; j < Ny; j++) {
                double x = (i + 0.5) * hx;
                double y = (j + 0.5) * hy;
                (*this)(i, j) = prob.param(x, y);
            }
    }
};
