#include "stepper.h"

#include "Profiler.h"

#include <iostream>

#include "elast.h"

struct waves : public elast::equations {
    vars_t initial(double x, double y) const {
        return vars_t(0, 0, 0, 0, 1);
    }
    param_t param(double x, double y) const {
        return param_t{1., 1., 1., x, y};
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
            dir d = dir::Y;
            return iW(d, U, p) * L(d, U, p).cwiseMax(0).cwiseProduct(W(d, U, p) * U);
        }
        throw;
    }
};

template<int sch>
void run() {
    PROFILE_ME;
    constexpr int p = 1;
    constexpr int order = 3;

    const char *sname[3] = {"LO", "HO", "TVD"};

    const std::string prefix("ivan_" + std::string(sname[sch]) + ".");

    waves prob;
    grid<waves::param_t> g(20, 20, 1.0, 1.0);

    g.fill(prob);

    stepper<waves, p, order, sch> stp(g);

    stp.lay.fill(prob);
    stp.lay.save(prefix + "0.vtk");
    std::cout << "Now " << prefix << std::endl;

    double t = 0;
    const double tmax = 0.1;
    const double dtout = tmax / 300;
    double tout = dtout;
    const double C = 0.14;
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
    run<scheme::HO>();
    run<scheme::TVD>();

    PROFILE_END;
    std::cout << profiler::getInstance() << std::endl;

    return 0;
}

