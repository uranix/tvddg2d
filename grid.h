#pragma once

#include <vector>
#include <cmath>

template<class PROBLEM>
struct grid {
    const int Nx, Ny;
    const double Lx, Ly;
    const double hx, hy;

    using param_t = typename PROBLEM::param_t;
    using cond_t = typename PROBLEM::cond_t;

    std::vector<param_t> data;
    std::vector<cond_t> xconddat, yconddat;

    grid(const int Nx, const int Ny, const double Lx, const double Ly)
        : Nx(Nx), Ny(Ny)
        , Lx(Lx), Ly(Ly)
        , hx(Lx/Nx), hy(Ly/Ny)
        , data(Nx * Ny)
        , xconddat((Nx - 1) * Ny), yconddat(Nx * (Ny - 1))
    { }

    const param_t &operator()(int i, int j) const {
        return data[i * Ny + j];
    }

    param_t &operator()(int i, int j) {
        return data[i * Ny + j];
    }

    // condition between (i, j) and (i+1, j) cell
    const cond_t &xcond(int i, int j) const {
        return xconddat[i * Ny + j];
    }

    cond_t &xcond(int i, int j) {
        return xconddat[i * Ny + j];
    }

    // condition between (i, j) and (i, j+1) cell
    const cond_t &ycond(int i, int j) const {
        return yconddat[i * (Ny - 1) + j];
    }

    cond_t &ycond(int i, int j) {
        return yconddat[i * (Ny - 1) + j];
    }

    void fill(const PROBLEM &prob) {
        for (int i = 0; i < Nx; i++)
            for (int j = 0; j < Ny; j++) {
                double x = (i + 0.5) * hx;
                double y = (j + 0.5) * hy;
                (*this)(i, j) = prob.param(x, y);
            }
    }

    double mark(const double x1, const double y1, const double x2, const double y2, const typename PROBLEM::cond_t &cond) {
        if (x1 == x2) {
            int i = x1 / hx - 0.5;
            int j1 = y1 / hy + 0.5;
            int j2 = y2 / hy + 0.5;
            double fuzz = std::abs((i+1) * hx - x1) + std::abs(j1 * hy - y1) + std::abs(j2 * hy - y2);
            if (i < 0 or i > Nx - 2)
                throw std::out_of_range("x1 = x2 is out of range");
            if (j1 < 0 or j1 > Ny)
                throw std::out_of_range("y1 is out of range");
            if (j2 < 0 or j2 > Ny)
                throw std::out_of_range("y2 is out of range");
            if (j1 > j2)
                std::swap(j1, j2);
            for (int j = j1; j < j2; j++)
                xcond(i, j) = cond;
            return fuzz;
        }
        if (y1 == y2) {
            int j = y1 / hy - 0.5;
            int i1 = x1 / hx + 0.5;
            int i2 = x2 / hx + 0.5;
            double fuzz = std::abs(i1 * hx - x1) + std::abs(i2 * hx - x2) + std::abs((j+1) * hy - y1);
            if (j < 0 or j > Ny - 2)
                throw std::out_of_range("y1 = y2 is out of range");
            if (i1 < 0 or i1 > Nx)
                throw std::out_of_range("x1 is out of range");
            if (i2 < 0 or i2 > Nx)
                throw std::out_of_range("x2 is out of range");
            if (i1 > i2)
                std::swap(i1, i2);
            for (int i = i1; i < i2; i++)
                ycond(i, j) = cond;
            return fuzz;
        }
        throw std::invalid_argument("Neither x1 == x2 nor y1 == y2");
    }
};
