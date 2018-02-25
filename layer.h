#pragma once

#include <functional>
#include <fstream>
#include <string>
#include <vector>

#include "matrix.h"
#include "cell.h"
#include "grid.h"

namespace dir {
    struct x     {};
    struct y     {};
    struct down  {};
    struct up    {};
    struct left  {};
    struct right {};
};

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
            }
    }
    void extrapolate() {
        for (int i = 0; i < g.Nx; i++)
            for (int j = 0; j < g.Ny; j++)
                (*this)(i, j).extrapolate();
    }
    void put(std::ofstream &f, float v) const {
        uint32_t data = __builtin_bswap32(*reinterpret_cast<uint32_t *>(&v));
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

template<typename PROBLEM, int p>
struct flux_layer {

    using param_t = typename PROBLEM::param_t;
    using vars_t = typename PROBLEM::vars_t;
    using cell_t = cell<vars_t, p>;

    const grid<param_t> &g;

    using flux_cell_t = flux_cell<vars_t, p>;
    using mat_t = matrix<p>;
    std::vector<flux_cell_t> data;

    flux_layer(const grid<param_t> &g) : g(g), data(g.Nx * g.Ny) {
    }
    const flux_cell_t operator()(int i, int j) const {
        return data[i * g.Ny + j];
    }
    flux_cell_t &operator()(int i, int j) {
        return data[i * g.Ny + j];
    }
    template<class problem>
    void compute_ho(const layer<problem, p> &lay, const problem &prob, double t) {
        auto &fl = *this;
        for (int i = 1; i < g.Nx; i++)
            for (int j = 0; j < g.Ny; j++)
                for (int k = 0; k <= p; k++) {
                    const std::pair<vars_t, vars_t> &ff = prob.template riemman<dir::x>(
                            lay(i-1, j)(p, k), lay(i, j)(0, k), g(i-1, j), g(i, j)
                        );
                    fl(i-1, j).F[p+1][k] = ff.first;
                    fl(i  , j).F[  0][k] = ff.second;
                }
        for (int i = 0; i < g.Nx; i++)
            for (int j = 1; j < g.Ny; j++)
                for (int k = 0; k <= p; k++) {
                    const std::pair<vars_t, vars_t> &ff = prob.template riemman<dir::y>(
                            lay(i, j-1)(k, p), lay(i, j)(k, 0), g(i, j-1), g(i, j)
                        );
                    fl(i, j-1).G[k][p+1] = ff.first;
                    fl(i, j  ).G[k][  0] = ff.second;
                }
        for (int j = 0; j < g.Ny; j++)
            for (int k = 0; k <= p; k++) {
                fl(0     , j).F[0  ][k] = prob.template bc<dir:: left>(lay(0     , j)( -1, k), g(0     , j), t);
                fl(g.Nx-1, j).F[p+1][k] = prob.template bc<dir::right>(lay(g.Nx-1, j)(p+1, k), g(g.Nx-1, j), t);
            }
        for (int i = 0; i < g.Nx; i++)
            for (int k = 0; k <= p; k++) {
                fl(i, 0     ).G[k][0  ] = prob.template bc<dir::down >(lay(i, 0     )(k,  -1), g(i, 0     ), t);
                fl(i, g.Ny-1).G[k][p+1] = prob.template bc<dir::up   >(lay(i, g.Ny-1)(k, p+1), g(i, g.Ny-1), t);
            }
        for (int i = 0; i < g.Nx; i++)
            for (int j = 0; j < g.Ny; j++) {
                const cell_t &c = lay(i, j);
                flux_cell_t &fc = fl(i, j);
                const param_t &pp = g(i, j);
                for (int k = 0; k <= p; k++) {
                    for (int m = 1; m <= p; m++) {
                        vars_t v = mat_t::F2F[m][0] * fc.F[0][k] + mat_t::F2F[m][p+2] * fc.F[p+1][k];
                        for (int z = 0; z <= p; z++)
                            v += mat_t::F2F[m][z+1] * prob.F(c(z, k), pp);
                        fc.F[m][k] = v;
                    }
                    for (int m = 1; m <= p; m++) {
                        vars_t v = mat_t::F2F[m][0] * fc.G[k][0] + mat_t::F2F[m][p+2] * fc.G[k][p+1];
                        for (int z = 0; z <= p; z++)
                            v += mat_t::F2F[m][z+1] * prob.G(c(k, z), pp);
                        fc.G[k][m] = v;
                    }
                }
            }
    }
};
