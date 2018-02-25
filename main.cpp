#include "layer.h"

#include <iostream>
#include <Eigen/Core>

struct vars : public Eigen::Vector3d {
    static constexpr int dim = 3;
    vars() { }
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

struct problem {
    using vars_t = vars;
    using param_t = param;
    vars_t initial(double x, double y) const {
        if (x < .3 && y < .3)
            return vars(x, y, sin(x) * cos(y));
        else
            return vars(0, 0, 0);
    }
    template<typename direction>
    std::pair<vars_t, vars_t> riemman(const vars_t &UL, const vars_t &UR, const param_t &pL, const param_t &pR) const {
        return std::make_pair(vars_t(), vars_t());
    }
    template<typename side>
    vars_t bc(const vars_t &U, const param_t &p, double t) const {
        return vars_t();
    }
    vars_t F(const vars_t &U, const param_t &p) const {
        return vars_t();
    }
    vars_t G(const vars_t &U, const param_t &p) const {
        return vars_t();
    }
};

int main() {
    problem prob;
    grid<param> g(10, 30, 1.0, 3.0);
    layer<problem, 3> lay(g);
    lay.fill([&prob] (double x, double y) { return prob.initial(x, y); });
    lay.extrapolate();
    lay.save("foo.vtk");
    flux_layer<problem, 3> flx(g);
    flx.compute_ho(lay, prob, 0);
    return 0;
}
