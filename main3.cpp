#include "stepper.h"

#include "Profiler.h"

#include <iostream>

#include "elast.h"

struct waves : public elast::equations {
    vars_t initial(double x, double y) const {
        double r2 = std::pow(x-0.5, 2) + std::pow(y-0.5, 2);
        double h = 0 * std::exp(-100 * r2);
        return vars_t(0, 0, h, h, 0);
    }
    param_t param(double x, double y) const {
        return param_t{1.0, 1.0, 1.0, x, y};
    }
    vars_t bc(side s, const vars_t &U, const param_t &p, double t) const {
        if (s == side::L) {
            dir d = dir::X;
            return iW(d, U, p) * L(d, U, p).cwiseMin(0).cwiseProduct(W(d, U, p) * U);
        }
        if (s == side::R) {
            dir d = dir::X;
            return iW(d, U, p) * L(d, U, p).cwiseMax(0).cwiseProduct(W(d, U, p) * U);
        }
        if (s == side::D) {
            dir d = dir::Y;
            return iW(d, U, p) * L(d, U, p).cwiseMin(0).cwiseProduct(W(d, U, p) * U);
        }
        if (s == side::U) {
            vars_t UR = U;

            double fx, fy;
            fx = 0;
            double s = (p.x - 0.5) / 0.05;
            fy = (std::abs(s) < 1) ? std::exp(-50*t - s*s) : 0;

            UR[4] = 2*fx - U[4];
            UR[3] = 2*fy - U[3];

            return riemman(dir::Y, U, UR, p, p).first;
        }
        throw;
    }
};

template<int sch>
void run() {
    PROFILE_ME;
    constexpr int p = 2;
    constexpr int order = 3;

    const char *sname[3] = {"LO", "HO", "TVD"};

    const std::string prefix("elast_" + std::string(sname[sch]) + ".");

    waves prob;
    grid<waves::param_t> g(51, 51, 1.0, 1.0);

    g.fill(prob);

    stepper<waves, p, order, sch> stp(g);

    stp.lay.fill(prob);
    stp.lay.save(prefix + "0.vtk");
    std::cout << "Now " << prefix << std::endl;

    double t = 0;
    const double tmax = 1;
    const double dtout = tmax / 300;
    double tout = dtout;
    const double C = 0.05;
    int step = 1;
    while (t < tmax) {
        double dt = stp.estimate_timestep(prob, C);
        stp.advance(prob, dt, t);
        if (t > tout) {
            std::cout << "t = " << t << std::endl;
            stp.lay.save(prefix + std::to_string(step) + ".vtk");
            tout += dtout;
        }
        t += dt;
        step++;
    }
}

int main() {
    PROFILE_ME;

//    run<scheme::LO>();
//    run<scheme::HO>();
    run<scheme::TVD>();

    PROFILE_END;
    std::cout << profiler::getInstance() << std::endl;

    return 0;
}

