#pragma once

#include <functional>
#include <fstream>
#include <string>
#include <vector>

#include "matrix.h"
#include "cell.h"
#include "grid.h"

enum struct dir { X, Y };
enum struct side { L, R, U, D };

template<typename PROBLEM, int p>
struct layer {
    using vars_t = typename PROBLEM::vars_t;
    using param_t = typename PROBLEM::param_t;

    const grid<param_t> &g;

    using cell_t = cell<vars_t, p>;
    using mat_t = matrix<p>;
    std::vector<cell_t> data;

    layer(const grid<param_t> &g) : g(g), data(g.Nx * g.Ny) {
    }
    const cell_t operator()(int i, int j) const {
        return data[i * g.Ny + j];
    }
    cell_t &operator()(int i, int j) {
        return data[i * g.Ny + j];
    }
    void fill(std::function<vars_t(double, double)> f) {
        for (int i = 0; i < g.Nx; i++)
            for (int j = 0; j < g.Ny; j++) {
                cell_t &c = (*this)(i, j);
                for (int ii = 0; ii <= p; ii++)
                    for (int jj = 0; jj <= p; jj++) {
                        double x = (i + mat_t::s[ii]) * g.hx;
                        double y = (j + mat_t::s[jj]) * g.hy;
                        c(ii, jj) = f(x, y);
                    }
                c.extrapolate();
            }
    }
    void extrapolate() {
        for (int i = 0; i < g.Nx; i++)
            for (int j = 0; j < g.Ny; j++)
                (*this)(i, j).extrapolate();
    }
    void put(std::ofstream &f, float v) const {
        union {
            float f;
            uint32_t u;
        } x;
        x.f = v;
        uint32_t data = __builtin_bswap32(x.u);
        f.write(reinterpret_cast<const char *>(&data), 4);
    }
    void save(const std::string &fn) const {
        std::ofstream f(fn, std::ios::binary);
        const double eps = 1e-10;

        f << "# vtk DataFile Version 3.0\nLayer\nBINARY\nDATASET RECTILINEAR_GRID\n";
        f << "DIMENSIONS " << (p+3)*g.Nx << " " << (p+3)*g.Ny << " 1\n";
        f << "X_COORDINATES " << (p+3)*g.Nx << " float\n";
        for (int i = 0; i < g.Nx; i++) {
            put(f, (i + eps)*g.hx);
            for (int j = 0; j <= p; j++)
                put(f, (i + mat_t::s[j])*g.hx);
            put(f, (i+1)*g.hx);
        }
        f << "Y_COORDINATES " << (p+3)*g.Ny << " float\n";
        for (int i = 0; i < g.Ny; i++) {
            put(f, (i + eps)*g.hy);
            for (int j = 0; j <= p; j++)
                put(f, (i + mat_t::s[j])*g.hy);
            put(f, (i+1)*g.hy);
        }
        f << "Z_COORDINATES 1 float\n";
        put(f, 0);
        f << "CELL_DATA " << ((p+3)*g.Nx - 1) * ((p+3)*g.Ny - 1) << "\n";
        f << "POINT_DATA " << ((p+3)*g.Nx) * ((p+3)*g.Ny) << "\n";
        for (int m = 0; m < vars_t::dim; m++) {
            f << "SCALARS " << vars_t::get_component_name(m) << " float\nLOOKUP_TABLE default\n";
            for (int j = 0; j < g.Ny; j++)
                for (int jj = -1; jj <= p+1; jj++)
                    for (int i = 0; i < g.Nx; i++) {
                        auto &c = (*this)(i, j);
                        for (int ii = -1; ii <= p+1; ii++)
                                put(f, c(ii, jj)[m]);
                    }
        }
    }
};
