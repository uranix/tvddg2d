#include "stepper.h"

#include "Profiler.h"

#include <iostream>
#include <Eigen/Core>

struct gdvars : public Eigen::Vector3d {
    static constexpr int dim = 3;
    gdvars() { *this << 0, 0, 0; }
    gdvars(double rho, double u, double eps) {
//        *this << rho, rho * u, rho * (eps + 0.5*u*u);
        *this << rho, u, eps;
    }
    template<class Other>
    gdvars(const Other &o) : Eigen::Vector3d(o) { }
    template<class Other>
    gdvars &operator=(const Other &o) {
        Eigen::Vector3d::operator=(o);
        return *this;
    }
    static const char *get_component_name(int i) {
        if (i == 0) return "density";
        if (i == 1) return "momentum";
        if (i == 2) return "energy";
        return "wtf";
    }
};

struct param {
};

struct gasdyn {
    using vars_t = gdvars;
    using param_t = param;

    using vec_t = Eigen::Vector3d;
    using mat_t = Eigen::Matrix3d;

    vars_t initial(double x, double y) const {
        if (x > .2 && x < .4 && y > .2 && y < .4)
            return gdvars(1, 1, 1);
        else
            return gdvars(0, 0, 0);
    }
    std::pair<vars_t, vars_t> riemman(dir d, const vars_t &UL, const vars_t &UR, const param_t &pL, const param_t &pR) const {
        vars_t f = F(d, UL, pL);
        return std::make_pair(f, f);
    }
    vars_t bc(side s, const vars_t &U, const param_t &p, double t) const {
        switch (s) {
            case side::L: return F(dir::X, vars_t(), p);
            case side::R: return F(dir::X, U, p);
            case side::D: return F(dir::Y, vars_t(), p);
            case side::U: return F(dir::Y, U, p);
        }
        throw std::invalid_argument("What the side");
    }

    vars_t F(dir d, const vars_t &U, const param_t &p) const {
        return (d == dir::X ? 1 : 2) * U;
    }
    mat_t W(dir d, const vars_t &U, const param_t &p) const {
        return mat_t::Identity();
    }
    mat_t iW(dir d, const vars_t &U, const param_t &p) const {
        return mat_t::Identity();
    }
    vec_t L(dir d, const vars_t &U, const param_t &p) const {
        return (d == dir::X ? 1 : 2) * vec_t(1, 1, 1);
    }
};

template<int sch>
void run() {
    PROFILE_ME;
    constexpr int p = 3;
    constexpr int order = 3;

    const char *sname[3] = {"LO", "HO", "TVD"};

    const std::string prefix("res_" + std::string(sname[sch]) + ".");

    gasdyn prob;
    grid<param> g(100, 200, 1.0, 1.5);

    stepper<gasdyn, p, order, sch> stp(g);
    stp.lay.fill([&prob] (double x, double y) { return prob.initial(x, y); });
    stp.lay.save(prefix + "0.vtk");

    double t = 0;
    const double dt = 0.03 * std::min(g.hx, g.hy / 2);
    const double tmax = 10 * dt;
    int step = 1;
    while (t < tmax) {
        stp.advance(prob, dt, t);
        if (step % 50 == 0) {
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

