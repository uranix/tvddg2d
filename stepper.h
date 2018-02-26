#pragma once

#include "layer.h"

struct base_stepper {
    template<typename PROBLEM, int p>
    void euler_step(const layer<PROBLEM, p> &from, layer<PROBLEM, p> &to,
            const flux_layer<PROBLEM, p> &flx, const double dt)
    {
        const auto &g = from.g;
        using mat_t = matrix<p>;

        for (int i = 0; i < g.Nx; i++)
            for (int j = 0; j < g.Ny; j++) {
                const auto &oc = from(i, j);
                auto &c = to(i, j);
                auto &fc = flx(i, j);

                for (int ii = 0; ii <= p; ii++)
                    for (int jj = 0; jj <= p; jj++) {
                        c(ii, jj) = oc(ii, jj)
                            - dt * (fc.F[ii+1][jj] - fc.F[ii][jj]) / (g.hx * mat_t::w[ii])
                            - dt * (fc.G[ii][jj+1] - fc.G[ii][jj]) / (g.hy * mat_t::w[jj]);
                    }
            }
    }
};

template<typename PROBLEM, int p, int ord> struct stepper;

template<typename PROBLEM, int p>
struct stepper<PROBLEM, p, 1> : public base_stepper {
    layer<PROBLEM, p> lay;
    flux_layer<PROBLEM, p> flx;

    stepper(const grid<typename PROBLEM::param_t> &g)
        : lay(g), flx(g)
    {
    }

    void advance(const PROBLEM &prob, const double dt, const double t) {
        flx.compute_ho(lay, prob, t);

        euler_step(lay, lay, flx, dt);

        lay.extrapolate();
    }
};

template<typename PROBLEM, int p>
struct stepper<PROBLEM, p, 2> : public base_stepper {
    layer<PROBLEM, p> lay;
    layer<PROBLEM, p> laymid;
    flux_layer<PROBLEM, p> flx;

    using mat_t = matrix<p>;

    stepper(const grid<typename PROBLEM::param_t> &g)
        : lay(g), laymid(g), flx(g)
    {
    }

    void advance(const PROBLEM &prob, const double dt, const double t) {
        flx.compute_ho(lay, prob, t);
        euler_step(lay, laymid, flx, dt);

        laymid.extrapolate();

        flx.compute_ho(laymid, prob, t + dt);
        euler_step(laymid, laymid, flx, dt);

        const auto &g = lay.g;

        for (int i = 0; i < g.Nx; i++)
            for (int j = 0; j < g.Ny; j++) {
                auto &u = lay(i, j);
                const auto &u2 = laymid(i, j);
                for (int ii = 0; ii <= p; ii++)
                    for (int jj = 0; jj <= p; jj++)
                        u(ii, jj) = 0.5 * u(ii, jj) + 0.5 * u2(ii, jj);
            }
        lay.extrapolate();
    }
};
