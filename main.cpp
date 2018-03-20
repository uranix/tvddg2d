#include "stepper.h"

#include "Profiler.h"

#include <iostream>

#include "elast.h"

struct waves : public elast::equations {
    vars_t initial(double x, double y) const {
        double x0 = 0.2;
        double sig = 0.01;
        return vars_t(-2, 0, 4, 2, 0) * std::exp(-0.5 * std::pow((x-x0)/sig, 2));
    }
    param_t param(double x, double y) const {
        if (std::abs(x - 1) < 0.5 and std::abs(y - 0.5) < 0.1)
            return param_t{1., 100., 200., x, y};
        return param_t{1., 1., 2., x, y};
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
            vars_t UL = U;

            UL[4] = -U[4];
            UL[3] = -U[3];

            return riemman(dir::Y, UL, U, p, p).first;
        }
        if (s == side::U) {
            vars_t UR = U;

            UR[4] = -U[4];
            UR[3] = -U[3];

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

    const std::string prefix("leveque_" + std::string(sname[sch]) + ".");

    waves prob;
    grid<waves::param_t> g(100, 50, 2.0, 1.0);

    g.fill(prob);

    stepper<waves, p, order, sch> stp(g);

    stp.lay.fill(prob);
    stp.lay.save(prefix + "0.vtk");
    std::cout << "Now " << prefix << std::endl;

    double t = 0;
    const double tmax = 1;
    const double dtout = tmax / 300;
    double tout = dtout;
    const double C = 0.04;
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
//    run<scheme::TVD>();

    PROFILE_END;
    std::cout << profiler::getInstance() << std::endl;

    return 0;
}

