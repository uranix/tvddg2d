#pragma once

#include <Eigen/Core>

namespace elast {

struct vars : public Eigen::Vector5d {
    static constexpr int dim = 5;

    vars() { *this << 0, 0, 0, 0, 0; }
    vars(double vx, double vy, double sxx, double syy, double sxy) {
        *this << vx << vy << sxx << syy << sxy;
    }
    template<class Other>
    vars(const Other &o) : Eigen::Vector5d(o) { }
    template<class Other>
    vars &operator=(const Other &o) {
        Eigen::Vector5d::operator=(o);
        return *this;
    }
    static const char *get_component_name(int i) {
        if (i == 0) return "vx";
        if (i == 1) return "vy";
        if (i == 2) return "sxx";
        if (i == 3) return "syy";
        if (i == 4) return "sxy";
        return nullptr;
    }
};

struct param {
    double rho, mu, lambda;
};

struct equations {
    using vars_t = vars;
    using param_t = param;

    using vec_t = Eigen::Vector5d;
    using mat_t = Eigen::Matrix5d;

    std::pair<vars_t, vars_t> riemman(dir d, const vars_t &UL, const vars_t &UR, const param_t &pL, const param_t &pR) const {
        if (d == dir::X) {
            throw
            vars_t f = 0.5 * (F(d, UL, pL) + F(d, UR, pR) + iWm * Lm.cwiseAbs().cwiseProduct(Wm * (UL - UR)));
            return std::make_pair(f, f);
        }

        throw
    }
    vars_t F(dir d, const vars_t &U, const param_t &p) const {
        const double vx = U[0];
        const double vy = U[1];
        const double sxx = U[2];
        const double syy = U[3];
        const double sxy = U[4];

        const double m = p.mu;
        const double l2m = p.lambda + 2 * p.mu;
        const double ir = 1. / p.rho;

        if (d == dir::X)
            return vars_t(-ir*sxx, -ir*sxy, -l2m*vx, -2*m*vx, -m*vy);
        else
            return vars_t(-ir*sxy, -ir*syy, -2*m*vy, -l2m*vy, -m*vx);
    }
    // The left eigenvectors of dF/dU by rows
    mat_t W(dir d, const vars_t &, const param_t &p) const {
        const double vx = U[0];
        const double vy = U[1];
        const double sxx = U[2];
        const double syy = U[3];
        const double sxy = U[4];

        const double r = p.rho;
        const double m = p.mu;
        const double l2m = p.lambda + 2 * p.mu;

        mat_t ret;

        const double as = std::sqrt(m * r);
        const double ap = std::sqrt(l2m * r);

        if (d == dir::X)
            ret <<
                0,  0, -2*m, l2m, 0,
                0,  as, 0, 0, 1,
                0, -as, 0, 0, 1,
                 ap, 0, 1, 0, 0,
                -ap, 0, 1, 0, 0
            ;
        else
            ret <<
                0,  0, l2m, -2*m, 0,
                 as, 0, 0, 0, 1,
                -as, 0, 0, 0, 1,
                0,  ap, 0, 1, 0,
                0, -ap, 0, 1, 0
            ;

        return ret;
    }
    // The inverse of W
    mat_t iW(dir d, const vars_t &, const param_t &p) const {
        const double vx = U[0];
        const double vy = U[1];
        const double sxx = U[2];
        const double syy = U[3];
        const double sxy = U[4];

        const double r = p.rho;
        const double m = p.mu;
        const double l2m = p.lambda + 2 * p.mu;

        mat_t ret;

        const double as = std::sqrt(m * r);
        const double ap = std::sqrt(l2m * r);

        if (d == dir::X)
            ret <<
                0, 0, 0, 0.5/ap, -0.5/ap,
                0, 0.5/as, -0.5/as, 0, 0,
                0, 0, 0, 0.5, 0.5,
                1/l2m, 0, 0, m/l2m, m/l2m,
                0, 0.5, 0.5, 0, 0
            ;
        else
            ret <<
                0, 0.5/as, -0.5/as, 0, 0,
                0, 0, 0, 0.5/ap, -0.5/ap,
                1/l2m, 0, 0, m/l2m, m/l2m,
                0, 0, 0, 0.5, 0.5,
                0, 0.5, 0.5, 0, 0
            ;

        return (0.5 / c) * ret;
    }
    // Eigenvalues of dF/dU (same order as in W)
    vec_t L(dir, const vars_t &, const param_t &p) const {
        const double r = p.rho;
        const double m = p.mu;
        const double l2m = p.lambda + 2 * p.mu;

        vec_t ret;

        const double cs = std::sqrt(m / r);
        const double cp = std::sqrt(l2m / r);

        ret << 0, -cs, cs, -cp, cp;

        return ret;
    }
};

};
