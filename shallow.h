#pragma once

#include <Eigen/Core>

namespace shallow {

struct vars : public Eigen::Vector3d {
    static constexpr int dim = 3;

    vars() { *this << 0, 0, 0; }
    vars(double h, double hu, double hv) {
        *this << h, hu, hv;
    }
    template<class Other>
    vars(const Other &o) : Eigen::Vector3d(o) { }
    template<class Other>
    vars &operator=(const Other &o) {
        Eigen::Vector3d::operator=(o);
        return *this;
    }
    static const char *get_component_name(int i) {
        if (i == 0) return "h";
        if (i == 1) return "hu";
        if (i == 2) return "hv";
        return nullptr;
    }
};

struct param { };

struct equations {
    using vars_t = vars;
    using param_t = param;

    using vec_t = Eigen::Vector3d;
    using mat_t = Eigen::Matrix3d;

    std::pair<vars_t, vars_t> riemman(dir d, const vars_t &UL, const vars_t &UR, const param_t &pL, const param_t &pR) const {
        const vars_t UM = 0.5 * (UL + UR);

        const auto &Wm = W(d, UM, pL);
        const auto &Lm = L(d, UM, pL);
        const auto &iWm = iW(d, UM, pL);

        vars_t f = 0.5 * (F(d, UL, pL) + F(d, UR, pR) + iWm * Lm.cwiseAbs().cwiseProduct(Wm * (UL - UR)));
        return std::make_pair(f, f);
    }
    vars_t F(dir d, const vars_t &U, const param_t &) const {
        const double h = U[0];
        const double u = U[1] / h;
        const double v = U[2] / h;

        if (d == dir::X)
            return vars_t(U[1], U[1]*u + 0.5*h*h, U[2]*u);
        else
            return vars_t(U[2], U[1]*v, U[2]*v + 0.5*h*h);
    }
    // The left eigenvectors of dF/dU by rows
    mat_t W(dir d, const vars_t &U, const param_t &) const {
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
    // The inverse of W
    mat_t iW(dir d, const vars_t &U, const param_t &) const {
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
    // Eigenvalues of dF/dU (same order as in W)
    vec_t L(dir d, const vars_t &U, const param_t &) const {
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

};
