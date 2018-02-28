#include "stepper.h"

#include "Profiler.h"

#include <iostream>
#include <Eigen/Core>

struct swvars : public Eigen::Vector3d {
    static constexpr int dim = 3;
    swvars() { *this << 0, 0, 0; }
    swvars(double h, double hu, double hv) {
        *this << h, hu, hv;
    }
    template<class Other>
    swvars(const Other &o) : Eigen::Vector3d(o) { }
    template<class Other>
    swvars &operator=(const Other &o) {
        Eigen::Vector3d::operator=(o);
        return *this;
    }
    static const char *get_component_name(int i) {
        if (i == 0) return "h";
        if (i == 1) return "hu";
        if (i == 2) return "hv";
        return "wtf";
    }
};

struct param {
};

struct shallow {
    using vars_t = swvars;
    using param_t = param;

    using vec_t = Eigen::Vector3d;
    using mat_t = Eigen::Matrix3d;

    vars_t initial(double x, double y) const {
        double x0 = 0.2, y0 = 0.2;
        double r2 = std::pow(x - x0, 2) + std::pow(y - y0, 2);
        double s2 = std::pow(0.05, 2);
        double h = 1 + 3 * std::exp(-r2 / s2);
        return vars_t(h, 0, 0);
    }
    std::pair<vars_t, vars_t> riemman(dir d, const vars_t &UL, const vars_t &UR, const param_t &pL, const param_t &pR) const {
        const vars_t UM = 0.5 * (UL + UR);

        const auto &Wm = W(d, UM, pL);
        const auto &Lm = L(d, UM, pL);
        const auto &iWm = iW(d, UM, pL);

        vars_t f = 0.5 * (F(d, UL, pL) + F(d, UR, pR) + iWm * Lm.cwiseAbs().cwiseProduct(Wm * (UL - UR)));
        return std::make_pair(f, f);
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

    vars_t F(dir d, const vars_t &U, const param_t &p) const {
        const double h = U[0];
        const double u = U[1] / h;
        const double v = U[2] / h;

        if (d == dir::X)
            return vars_t(U[1], U[1]*u + 0.5*h*h, U[2]*u);
        else
            return vars_t(U[2], U[1]*v, U[2]*v + 0.5*h*h);
    }
    mat_t W(dir d, const vars_t &U, const param_t &p) const {
        const double h = U[0];
        const double u = U[1] / h;
        const double v = U[2] / h;
        const double c = std::sqrt(h);

        mat_t ret;

        if (d == dir::X)
            ret << -v, 0, 1, -c-u, 1, 0, c-u, 1, 0;
        else
            ret << -u, 1, 0, -c-v, 0, 1, c-v, 0, 1;

        return ret;
    }
    mat_t iW(dir d, const vars_t &U, const param_t &p) const {
        const double h = U[0];
        const double u = U[1] / h;
        const double v = U[2] / h;
        const double c = std::sqrt(h);

        mat_t ret;

        if (d == dir::X)
            ret << 0, -1, 1, 0, c-u, c+u, 2*c, -v, v;
        else
            ret << 0, -1, 1, 2*c, -u, u, 0, c-v, c+v;

        return (0.5 / c) * ret;
    }
    vec_t L(dir d, const vars_t &U, const param_t &p) const {
        const double h = U[0];
        const double u = U[1] / h;
        const double v = U[2] / h;
        const double c = std::sqrt(h);
        if (d == dir::X)
            return vec_t(u, u-c, u+c);
        else
            return vec_t(v, v-c, v+c);
    }
};

template<int sch>
void run() {
    PROFILE_ME;
    constexpr int p = 3;
    constexpr int order = 3;

    const char *sname[3] = {"LO", "HO", "TVD"};

    const std::string prefix("shallow_" + std::string(sname[sch]) + ".");

    shallow prob;
    grid<param> g(50, 50, 1.0, 1.0);

    stepper<shallow, p, order, sch> stp(g);
    stp.lay.fill(prob);
    stp.lay.save(prefix + "0.vtk");
    std::cout << "Now " << prefix << std::endl;

    double t = 0;
    const double tmax = 0.3;
    const double C = 0.02;
    int step = 1;
    while (t < tmax) {
        double dt = stp.estimate_timestep(prob, C);
        stp.advance(prob, dt, t);
        if (step % 10 == 0) {
            std::cout << "t = " << t << std::endl;
            stp.lay.save(prefix + std::to_string(step) + ".vtk");
        }
        t += dt;
        step++;
    }
}

int main() {
    PROFILE_ME;

    run<scheme::LO>();
    run<scheme::HO>();
    run<scheme::TVD>();

    PROFILE_END;
    std::cout << profiler::getInstance() << std::endl;

    return 0;
}

