#pragma once

#include "layer.h"
#include "flux_layer.h"

template<typename PROBLEM, int p>
struct base_stepper {
    layer<PROBLEM, p> lay;
    base_stepper(const grid<typename PROBLEM::param_t> &g) : lay(g) { }
    // y := a y + b x
    void combine(const layer<PROBLEM, p> &x, layer<PROBLEM, p> &y, const double alpha, const double beta) {
        PROFILE_ME;
        const auto &g = x.g;

        if (alpha != 0) {
            for (int i = 0; i < g.Nx; i++)
                for (int j = 0; j < g.Ny; j++) {
                    const auto &xc = x(i, j);
                    auto &yc = y(i, j);

                    for (int ii = 0; ii <= p; ii++)
                        for (int jj = 0; jj <= p; jj++)
                            yc(ii, jj) = alpha * yc(ii, jj) + beta * xc(ii, jj);
                }
        } else {
            for (int i = 0; i < g.Nx; i++)
                for (int j = 0; j < g.Ny; j++) {
                    const auto &xc = x(i, j);
                    auto &yc = y(i, j);

                    for (int ii = 0; ii <= p; ii++)
                        for (int jj = 0; jj <= p; jj++)
                            yc(ii, jj) = beta * xc(ii, jj);
                }
        }
    }
    template<int sch>
    void euler_step(const layer<PROBLEM, p> &from, layer<PROBLEM, p> &to,
            const flux_layer<PROBLEM, p, sch> &flx, const double dt)
    {
        PROFILE_ME;
        const auto &g = from.g;
        using quad_t = quadrature<p>;

        for (int i = 0; i < g.Nx; i++)
            for (int j = 0; j < g.Ny; j++) {
                const auto &oc = from(i, j);
                auto &c = to(i, j);
                auto &fc = flx(i, j);

                for (int ii = 0; ii <= p; ii++)
                    for (int jj = 0; jj <= p; jj++) {
                        c(ii, jj) = oc(ii, jj)
                            - dt * (fc.F[ii+1][jj] - fc.F[ii][jj]) / (g.hx * quad_t::w[ii])
                            - dt * (fc.G[ii][jj+1] - fc.G[ii][jj]) / (g.hy * quad_t::w[jj]);
                    }
            }
    }
    double estimate_timestep(const PROBLEM &prob, const double C) {
        const auto &g = lay.g;

        double axMax = 0, ayMax = 0;

        for (int i = 0; i < g.Nx; i++)
            for (int j = 0; j < g.Ny; j++)
                for (int ii = 0; ii <= p; ii++)
                    for (int jj = 0; jj <= p; jj++) {
                        const auto &Lx = prob.L(dir::X, lay(i, j)(ii, jj), g(i, j));
                        const auto &Ly = prob.L(dir::Y, lay(i, j)(ii, jj), g(i, j));
                        double ax = Lx.cwiseAbs().maxCoeff();
                        double ay = Ly.cwiseAbs().maxCoeff();
                        if (ax > axMax)
                            axMax = ax;
                        if (ay > ayMax)
                            ayMax = ay;
                    }

        return C / (axMax / g.hx + ayMax / g.hy);
    }
};

template<typename PROBLEM, int p, int ord, int sch> struct stepper;

template<typename PROBLEM, int p, int sch>
struct stepper<PROBLEM, p, 1, sch> : public base_stepper<PROBLEM, p> {
    flux_layer<PROBLEM, p, sch> flx;

    stepper(const grid<typename PROBLEM::param_t> &g)
        : base_stepper<PROBLEM, p>(g), flx(g)
    {
    }

    void advance(const PROBLEM &prob, const double dt, const double t) {
        auto &lay = this->lay;
        flx.compute(lay, prob, t);

        this->euler_step(lay, lay, flx, dt);

        lay.extrapolate();
    }
};

template<typename PROBLEM, int p, int sch>
struct stepper<PROBLEM, p, 2, sch> : public base_stepper<PROBLEM, p> {
    layer<PROBLEM, p> laymid;
    flux_layer<PROBLEM, p, sch> flx;

    using quad_t = quadrature<p>;

    stepper(const grid<typename PROBLEM::param_t> &g)
        : base_stepper<PROBLEM, p>(g), laymid(g), flx(g)
    {
    }

    void advance(const PROBLEM &prob, const double dt, const double t) {
        auto &lay = this->lay;
        flx.compute(lay, prob, t);
        this->euler_step(lay, laymid, flx, dt);

        laymid.extrapolate();

        flx.compute(laymid, prob, t + dt);
        this->euler_step(laymid, laymid, flx, dt);

        this->combine(laymid, lay, 0.5, 0.5);
        lay.extrapolate();
    }
};

template<typename PROBLEM, int p, int sch>
struct stepper<PROBLEM, p, 3, sch> : public base_stepper<PROBLEM, p> {
    layer<PROBLEM, p> laymid;
    flux_layer<PROBLEM, p, sch> flx;

    using quad_t = quadrature<p>;

    stepper(const grid<typename PROBLEM::param_t> &g)
        : base_stepper<PROBLEM, p>(g), laymid(g), flx(g)
    {
    }

    void advance(const PROBLEM &prob, const double dt, const double t) {
        auto &lay = this->lay;
        flx.compute(lay, prob, t);
        this->euler_step(lay, laymid, flx, dt);
        laymid.extrapolate();

        flx.compute(laymid, prob, t + dt);
        this->euler_step(laymid, laymid, flx, dt);
        this->combine(lay, laymid, 0.25, 0.75);
        laymid.extrapolate();

        flx.compute(laymid, prob, t + 0.5*dt);
        this->euler_step(laymid, laymid, flx, dt);
        this->combine(laymid, lay, 1./3, 2./3);
        lay.extrapolate();
    }
};
