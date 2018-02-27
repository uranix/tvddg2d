#pragma once

#include "layer.h"

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
    void compute_lo(const layer<PROBLEM, p> &lay, const PROBLEM &prob, double t) {
        auto &fl = *this;
        for (int i = 1; i < g.Nx; i++)
            for (int j = 0; j < g.Ny; j++)
                for (int k = 0; k <= p; k++) {
                    const std::pair<vars_t, vars_t> &ff = prob.riemman(
                            dir::X, lay(i-1, j)(p, k), lay(i, j)(0, k), g(i-1, j), g(i, j)
                        );
                    fl(i-1, j).F[p+1][k] = ff.first;
                    fl(i  , j).F[  0][k] = ff.second;
                }
        for (int i = 0; i < g.Nx; i++)
            for (int j = 1; j < g.Ny; j++)
                for (int k = 0; k <= p; k++) {
                    const std::pair<vars_t, vars_t> &ff = prob.riemman(
                            dir::Y, lay(i, j-1)(k, p), lay(i, j)(k, 0), g(i, j-1), g(i, j)
                        );
                    fl(i, j-1).G[k][p+1] = ff.first;
                    fl(i, j  ).G[k][  0] = ff.second;
                }
        for (int j = 0; j < g.Ny; j++)
            for (int k = 0; k <= p; k++) {
                fl(0     , j).F[0  ][k] = prob.bc(side::L, lay(0     , j)( -1, k), g(0     , j), t);
                fl(g.Nx-1, j).F[p+1][k] = prob.bc(side::R, lay(g.Nx-1, j)(p+1, k), g(g.Nx-1, j), t);
            }
        for (int i = 0; i < g.Nx; i++)
            for (int k = 0; k <= p; k++) {
                fl(i, 0     ).G[k][0  ] = prob.bc(side::D, lay(i, 0     )(k,  -1), g(i, 0     ), t);
                fl(i, g.Ny-1).G[k][p+1] = prob.bc(side::U, lay(i, g.Ny-1)(k, p+1), g(i, g.Ny-1), t);
            }
        for (int i = 0; i < g.Nx; i++)
            for (int j = 0; j < g.Ny; j++) {
                const cell_t &c = lay(i, j);
                flux_cell_t &fc = fl(i, j);
                const param_t &pp = g(i, j);

                for (int k = 0; k <= p; k++)
                    for (int m = 1; m <= p; m++) {
                        const std::pair<vars_t, vars_t> &ff = prob.riemman(
                                dir::X, c(m-1, k), c(m, k), pp, pp
                            );

                        fc.F[m][k] = ff.first;

                        const std::pair<vars_t, vars_t> &gg = prob.riemman(
                                dir::Y, c(k, m-1), c(k, m), pp, pp
                            );

                        fc.G[k][m] = gg.first;
                    }
            }
    }
    void compute_ho(const layer<PROBLEM, p> &lay, const PROBLEM &prob, double t) {
        auto &fl = *this;
        for (int i = 1; i < g.Nx; i++)
            for (int j = 0; j < g.Ny; j++)
                for (int k = 0; k <= p; k++) {
                    const std::pair<vars_t, vars_t> &ff = prob.riemman(
                            dir::X, lay(i-1, j)(p+1, k), lay(i, j)(-1, k), g(i-1, j), g(i, j)
                        );
                    fl(i-1, j).F[p+1][k] = ff.first;
                    fl(i  , j).F[  0][k] = ff.second;
                }
        for (int i = 0; i < g.Nx; i++)
            for (int j = 1; j < g.Ny; j++)
                for (int k = 0; k <= p; k++) {
                    const std::pair<vars_t, vars_t> &ff = prob.riemman(
                            dir::Y, lay(i, j-1)(k, p+1), lay(i, j)(k, -1), g(i, j-1), g(i, j)
                        );
                    fl(i, j-1).G[k][p+1] = ff.first;
                    fl(i, j  ).G[k][  0] = ff.second;
                }
        for (int j = 0; j < g.Ny; j++)
            for (int k = 0; k <= p; k++) {
                fl(0     , j).F[0  ][k] = prob.bc(side::L, lay(0     , j)( -1, k), g(0     , j), t);
                fl(g.Nx-1, j).F[p+1][k] = prob.bc(side::R, lay(g.Nx-1, j)(p+1, k), g(g.Nx-1, j), t);
            }
        for (int i = 0; i < g.Nx; i++)
            for (int k = 0; k <= p; k++) {
                fl(i, 0     ).G[k][0  ] = prob.bc(side::D, lay(i, 0     )(k,  -1), g(i, 0     ), t);
                fl(i, g.Ny-1).G[k][p+1] = prob.bc(side::U, lay(i, g.Ny-1)(k, p+1), g(i, g.Ny-1), t);
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
    void compute(const layer<PROBLEM, p> &lay, const PROBLEM &prob, double t) {
        return compute_lo(lay, prob, t);
    }
};
