#include "stepper.h"

#include <iostream>
#include <Eigen/Core>

struct vars : public Eigen::Vector3d {
    static constexpr int dim = 3;
    vars() { *this << 0, 0, 0; }
    vars(double rho, double u, double eps) {
//        *this << rho, rho * u, rho * (eps + 0.5*u*u);
        *this << rho, u, eps;
    }
    template<class Other>
    vars(const Other &o) : Eigen::Vector3d(o) { }
    template<class Other>
    vars &operator=(const Other &o) {
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
    using vars_t = vars;
    using param_t = param;
    vars_t initial(double x, double y) const {
        if (x > .2 && x < .4 && y > .2 && y < .8)
            return vars(1, 1, 1);
        else
            return vars(0, 0, 0);
    }
    std::pair<vars_t, vars_t> riemman(dir d, const vars_t &UL, const vars_t &UR, const param_t &pL, const param_t &pR) const {
        vars_t f = (d == dir::X) ? F(UL, pL) : G(UL, pL);
        return std::make_pair(f, f);
    }
    vars_t bc(side s, const vars_t &U, const param_t &p, double t) const {
        switch (s) {
            case side::L: return F(vars_t(), p);
            case side::R: return F(U, p);
            case side::D: return G(vars_t(), p);
            case side::U: return G(U, p);
        }
        throw std::invalid_argument("What the side");
    }

    vars_t F(const vars_t &U, const param_t &p) const {
        return 1 * U;
    }
    vars_t G(const vars_t &U, const param_t &p) const {
        return 0 * U;
    }
};

int main() {
    constexpr int p = 1;
    constexpr int order = 2;
    gasdyn prob;
    grid<param> g(50, 1, 1.0, 1.0);
    stepper<gasdyn, p, order> stp(g);
    stp.lay.fill([&prob] (double x, double y) { return prob.initial(x, y); });
    stp.lay.save("res.0.vtk");

    double t = 0;
    const double dt = 0.0001;
    const double tmax = 0.1;
    int step = 1;
    while (t < tmax) {
        stp.advance(prob, dt, t);
        if (step % 10 == 0)
            stp.lay.save("res." + std::to_string(step) + ".vtk");
        t += dt;
        step++;
    }
    return 0;
}
