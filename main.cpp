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

int main() {
    grid g(10, 30, 1.0, 3.0);
    layer<vars, 3> lay(g);
    lay.fill([] (double x, double y) -> vars {
            if (x < .3 && y < .3)
                return vars(x, y, sin(x) * cos(y));
            else
                return vars(0, 0, 0);
        });
    lay.extrapolate();
    lay.save("foo.vtk");
    return 0;
}
