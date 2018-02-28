#pragma once

#include "layer.h"

#include "Profiler.h"

template<typename vec>
vec minmod(const vec &x, const vec &y) {
    const auto &sx = x.cwiseSign();
    const auto &sy = y.cwiseSign();
    return 0.5 * (sx + sy).cwiseProduct(x.cwiseAbs().cwiseMin(y.cwiseAbs()));
}

namespace scheme {
    enum { LO = 0, HO = 1, TVD = 2 };
}

template<typename PROBLEM, int p, int sch>
struct flux_layer {

    using param_t = typename PROBLEM::param_t;
    using vars_t = typename PROBLEM::vars_t;

    using vec_t = typename PROBLEM::vec_t;
    using mat_t = typename PROBLEM::mat_t;
    using cell_t = cell<vars_t, p>;

    using flux_cell_t = flux_cell<vars_t, p>;
    using flux_cell_lo_t = flux_cell_lo<vars_t, vec_t, mat_t, p>;
    using quad_t = quadrature<p>;

    const grid<param_t> &g;
    std::vector<flux_cell_lo_t> low_order_data;
    std::vector<flux_cell_t> high_order_data;

    flux_layer(const grid<param_t> &g)
        : g(g), low_order_data(g.Nx * g.Ny), high_order_data(g.Nx * g.Ny)
    {
    }
    const flux_cell_lo_t &low_order(int i, int j) const {
        return low_order_data[i * g.Ny + j];
    }
    flux_cell_lo_t &low_order(int i, int j) {
        return low_order_data[i * g.Ny + j];
    }
    const flux_cell_t &high_order(int i, int j) const {
        return high_order_data[i * g.Ny + j];
    }
    flux_cell_t &high_order(int i, int j) {
        return high_order_data[i * g.Ny + j];
    }
    void compute_eigen(dir d, const vars_t &uL, const vars_t &uR, const param_t &pp, const PROBLEM &prob, invars<vec_t, mat_t> &I) {
        vars_t uM = 0.5 * (uL + uR);
        I.iW = prob.iW(d, uM, pp);
        I.L  = prob.L(d, uM, pp);
        I.W  = prob.W(d, uM, pp);
        I.D  = I.W * (uR - uL);
    }
    void compute_lo(const layer<PROBLEM, p> &lay, const PROBLEM &prob, double t) {
        PROFILE_ME;
        #pragma omp parallel
        {
            #pragma omp for nowait collapse(3)
            for (int i = 1; i < g.Nx; i++)
                for (int j = 0; j < g.Ny; j++)
                    for (int k = 0; k <= p; k++) {
                        const vars_t &uL = lay(i-1, j)(p, k);
                        const vars_t &uR = lay(i,   j)(0, k);
                        const param_t &pL = g(i-1, j);
                        const param_t &pR = g(i  , j);

                        const std::pair<vars_t, vars_t> &ff = prob.riemman(dir::X, uL, uR, pL, pR);

                        flux_cell_lo_t &fL = low_order(i-1, j);
                        flux_cell_lo_t &fR = low_order(i  , j);

                        fL.F[p+1][k] = ff.first;
                        fR.F[  0][k] = ff.second;

                        compute_eigen(dir::X, uL, uR, pL, prob, fL.Ix[p+1][k]);
                        compute_eigen(dir::X, uL, uR, pR, prob, fR.Ix[  0][k]);
                    }
            #pragma omp for nowait collapse(3)
            for (int i = 0; i < g.Nx; i++)
                for (int j = 1; j < g.Ny; j++)
                    for (int k = 0; k <= p; k++) {
                        const vars_t &uL = lay(i, j-1)(k, p);
                        const vars_t &uR = lay(i, j  )(k, 0);
                        const param_t &pL = g(i, j-1);
                        const param_t &pR = g(i, j  );

                        const std::pair<vars_t, vars_t> &ff = prob.riemman(dir::Y, uL, uR, pL, pR);

                        flux_cell_lo_t &fL = low_order(i, j-1);
                        flux_cell_lo_t &fR = low_order(i, j  );

                        fL.G[k][p+1] = ff.first;
                        fR.G[k][  0] = ff.second;

                        compute_eigen(dir::Y, uL, uR, pL, prob, fL.Iy[k][p+1]);
                        compute_eigen(dir::Y, uL, uR, pR, prob, fR.Iy[k][  0]);
                    }
            #pragma omp for nowait collapse(2)
            for (int j = 0; j < g.Ny; j++)
                for (int k = 0; k <= p; k++) {
                    const vars_t &uL = lay(0,      j)(0, k);
                    const vars_t &uR = lay(g.Nx-1, j)(p, k);
                    const param_t &pL = g(0,      j);
                    const param_t &pR = g(g.Nx-1, j);

                    flux_cell_lo_t &fL = low_order(0,      j);
                    flux_cell_lo_t &fR = low_order(g.Nx-1, j);

                    fL.F[0  ][k] = prob.bc(side::L, uL, pL, t);
                    fR.F[p+1][k] = prob.bc(side::R, uR, pR, t);

                    compute_eigen(dir::X, uL, uL, pL, prob, fL.Ix[  0][k]);
                    compute_eigen(dir::X, uR, uR, pR, prob, fR.Ix[p+1][k]);
                }
            #pragma omp for nowait collapse(2)
            for (int i = 0; i < g.Nx; i++)
                for (int k = 0; k <= p; k++) {
                    const vars_t &uL = lay(i, 0     )(k, 0);
                    const vars_t &uR = lay(i, g.Ny-1)(k, p);
                    const param_t &pL = g(i, 0     );
                    const param_t &pR = g(i, g.Ny-1);

                    flux_cell_lo_t &fL = low_order(i, 0     );
                    flux_cell_lo_t &fR = low_order(i, g.Ny-1);

                    fL.G[k][0  ] = prob.bc(side::D, uL, pL, t);
                    fR.G[k][p+1] = prob.bc(side::U, uR, pR, t);

                    compute_eigen(dir::Y, uL, uL, pL, prob, fL.Iy[k][  0]);
                    compute_eigen(dir::Y, uR, uR, pR, prob, fR.Iy[k][p+1]);
                }
            #pragma omp for nowait collapse(2)
            for (int i = 0; i < g.Nx; i++)
                for (int j = 0; j < g.Ny; j++) {
                    const cell_t &c = lay(i, j);
                    flux_cell_lo_t &fc = low_order(i, j);
                    const param_t &pp = g(i, j);

                    for (int k = 0; k <= p; k++)
                        for (int m = 1; m <= p; m++) {
                            const vars_t &uL = c(m-1, k);
                            const vars_t &uR = c(m  , k);
                            auto &Ix = fc.Ix[m][k];

                            compute_eigen(dir::X, uL, uR, pp, prob, Ix);

                            // CIR flux, F_{j+1/2} = 1/2 (F(U_j) + F(U_{j+1}) - 1/2 iW |L| W (uR - uL))
                            fc.F[m][k] = 0.5 * (prob.F(dir::X, uL, pp) + prob.F(dir::X, uR, pp) - Ix.iW * Ix.L.cwiseAbs().cwiseProduct(Ix.D));
                        }

                    for (int k = 0; k <= p; k++)
                        for (int m = 1; m <= p; m++) {
                            const vars_t &uL = c(k, m-1);
                            const vars_t &uR = c(k, m  );
                            auto &Iy = fc.Iy[k][m];

                            compute_eigen(dir::Y, uL, uR, pp, prob, Iy);

                            fc.G[k][m] = 0.5 * (prob.F(dir::Y, uL, pp) + prob.F(dir::Y, uR, pp) - Iy.iW * Iy.L.cwiseAbs().cwiseProduct(Iy.D));
                        }
                }
        }
    }
    void compute_ho(const layer<PROBLEM, p> &lay, const PROBLEM &prob, double t) {
        PROFILE_ME;
        #pragma omp parallel
        {
            #pragma omp for nowait collapse(3)
            for (int i = 1; i < g.Nx; i++)
                for (int j = 0; j < g.Ny; j++)
                    for (int k = 0; k <= p; k++) {
                        const std::pair<vars_t, vars_t> &ff = prob.riemman(
                                dir::X, lay(i-1, j)(p+1, k), lay(i, j)(-1, k), g(i-1, j), g(i, j)
                            );
                        high_order(i-1, j).F[p+1][k] = ff.first;
                        high_order(i  , j).F[  0][k] = ff.second;
                    }
            #pragma omp for nowait collapse(3)
            for (int i = 0; i < g.Nx; i++)
                for (int j = 1; j < g.Ny; j++)
                    for (int k = 0; k <= p; k++) {
                        const std::pair<vars_t, vars_t> &ff = prob.riemman(
                                dir::Y, lay(i, j-1)(k, p+1), lay(i, j)(k, -1), g(i, j-1), g(i, j)
                            );
                        high_order(i, j-1).G[k][p+1] = ff.first;
                        high_order(i, j  ).G[k][  0] = ff.second;
                    }
            #pragma omp for nowait collapse(2)
            for (int j = 0; j < g.Ny; j++)
                for (int k = 0; k <= p; k++) {
                    high_order(0     , j).F[0  ][k] = prob.bc(side::L, lay(0     , j)( -1, k), g(0     , j), t);
                    high_order(g.Nx-1, j).F[p+1][k] = prob.bc(side::R, lay(g.Nx-1, j)(p+1, k), g(g.Nx-1, j), t);
                }
            #pragma omp for nowait collapse(2)
            for (int i = 0; i < g.Nx; i++)
                for (int k = 0; k <= p; k++) {
                    high_order(i, 0     ).G[k][0  ] = prob.bc(side::D, lay(i, 0     )(k,  -1), g(i, 0     ), t);
                    high_order(i, g.Ny-1).G[k][p+1] = prob.bc(side::U, lay(i, g.Ny-1)(k, p+1), g(i, g.Ny-1), t);
                }
        }
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < g.Nx; i++)
            for (int j = 0; j < g.Ny; j++) {
                const cell_t &c = lay(i, j);
                flux_cell_t &fc = high_order(i, j);
                const param_t &pp = g(i, j);
                for (int k = 0; k <= p; k++) {
                    for (int m = 1; m <= p; m++) {
                        vars_t v = quad_t::F2F[m][0] * fc.F[0][k] + quad_t::F2F[m][p+2] * fc.F[p+1][k];
                        for (int z = 0; z <= p; z++)
                            v += quad_t::F2F[m][z+1] * prob.F(dir::X, c(z, k), pp);
                        fc.F[m][k] = v;
                    }
                    for (int m = 1; m <= p; m++) {
                        vars_t v = quad_t::F2F[m][0] * fc.G[k][0] + quad_t::F2F[m][p+2] * fc.G[k][p+1];
                        for (int z = 0; z <= p; z++)
                            v += quad_t::F2F[m][z+1] * prob.F(dir::Y, c(k, z), pp);
                        fc.G[k][m] = v;
                    }
                }
            }
    }
    void limit_fluxes() {
        PROFILE_ME;
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < g.Nx; i++)
            for (int j = 0; j < g.Ny; j++) {
                for (int ii = 0; ii <= p+1; ii++)
                    for (int jj = 0; jj <= p; jj++) {
                        const auto &Ix = low_order(i, j).Ix[ii][jj];
                        const auto &g0 = Ix.D;
                        const auto &gL = ii > 0 ? low_order(i, j).Ix[ii-1][jj].D
                            : i > 0 ? low_order(i-1, j).Ix[p][jj].D : g0;
                        const auto &gR = ii <= p ? low_order(i, j).Ix[ii+1][jj].D
                            : i < g.Nx - 1 ? low_order(i+1, j).Ix[1][jj].D : g0;
                        const auto &mL = minmod(g0, gL);
                        const auto &mR = minmod(g0, gR);
                        const auto &phi0 = 0.5 * (Ix.L.cwiseAbs().cwiseProduct(mL + mR) + Ix.L.cwiseProduct(mL - mR));
                        const auto &phi = minmod<decltype(Ix.D)>(Ix.W * (high_order(i, j).F[ii][jj] - low_order(i, j).F[ii][jj]), phi0);
                        high_order(i, j).F[ii][jj] = low_order(i, j).F[ii][jj] + Ix.iW * phi;
                    }
                for (int ii = 0; ii <= p; ii++)
                    for (int jj = 0; jj <= p+1; jj++) {
                        const auto &Iy = low_order(i, j).Iy[ii][jj];
                        const auto &g0 = Iy.D;
                        const auto &gL = jj > 0 ? low_order(i, j).Iy[ii][jj-1].D
                            : j > 0 ? low_order(i, j-1).Iy[ii][p].D : g0;
                        const auto &gR = jj <= p ? low_order(i, j).Iy[ii][jj+1].D
                            : j < g.Ny - 1 ? low_order(i, j+1).Iy[ii][1].D : g0;
                        const auto &mL = minmod(g0, gL);
                        const auto &mR = minmod(g0, gR);
                        const auto &phi0 = 0.5 * (Iy.L.cwiseAbs().cwiseProduct(mL + mR) + Iy.L.cwiseProduct(mL - mR));
                        const auto &phi = minmod<decltype(Iy.D)>(Iy.W * (high_order(i, j).G[ii][jj] - low_order(i, j).G[ii][jj]), phi0);
                        high_order(i, j).G[ii][jj] = low_order(i, j).G[ii][jj] + Iy.iW * phi;
                    }
            }
    }
    void compute(const layer<PROBLEM, p> &lay, const PROBLEM &prob, double t) {
        if (sch != scheme::LO)
            compute_ho(lay, prob, t);
        if (sch != scheme::HO)
            compute_lo(lay, prob, t);
        if (sch == scheme::TVD)
            limit_fluxes();
    }
    const flux_cell_t &operator()(int i, int j) const {
        if (sch == scheme::LO)
            return low_order(i, j);
        return high_order(i, j);
    }
};
