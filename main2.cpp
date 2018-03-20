#include "stepper.h"

#include "Profiler.h"

#include <iostream>

#include "shallow.h"

struct drop : public shallow::equations {
    vars_t initial(double x, double y) const {
        double x0 = 0.2, y0 = 0.2;
        double r2 = std::pow(x - x0, 2) + std::pow(y - y0, 2);
        double s2 = std::pow(0.05, 2);
        double h = 1 + 3 * std::exp(-r2 / s2);
        return vars_t(h, 0, 0);
    }
    vars_t bc(side s, const vars_t &U, const param_t &p, double t) const {
        if (s == side::L) {
            vars_t UL = U;
            UL[1] = -U[1];
            return riemman(dir::X, UL, U, p, p).first;
        }
        if (s == side::R) {
            vars_t UR = U;
            UR[1] = -U[1];
            return riemman(dir::X, U, UR, p, p).first;
        }
        if (s == side::D) {
            vars_t UL = U;
            UL[2] = -U[2];
            return riemman(dir::Y, UL, U, p, p).first;
        }
        if (s == side::U) {
            vars_t UR = U;
            UR[2] = -U[2];
            return riemman(dir::Y, U, UR, p, p).first;
        }
        return vars_t();
    }
};

template<int sch>
void run() {
    PROFILE_ME;
    constexpr int p = 4;
    constexpr int order = 3;

    const char *sname[3] = {"LO", "HO", "TVD"};

    const std::string prefix("shallow_" + std::string(sname[sch]) + ".");

    drop prob;
    grid<drop::param_t> g(30, 30, 1.0, 1.0);

    stepper<drop, p, order, sch> stp(g);

    stp.lay.fill(prob);
    stp.lay.save(prefix + "0.vtk");
    std::cout << "Now " << prefix << std::endl;

    double t = 0;
    const double tmax = 1;
    const double dtout = tmax / 300;
    double tout = dtout;
    const double C = 0.02;
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

    run<scheme::LO>();
//    run<scheme::HO>();
    run<scheme::TVD>();

    PROFILE_END;
    std::cout << profiler::getInstance() << std::endl;

    return 0;
}

