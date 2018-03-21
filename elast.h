#pragma once

#include <Eigen/Core>

namespace elast {

struct vars : public Eigen::Matrix<double, 5, 1> {
    static constexpr int dim = 5;

    vars() { *this << 0, 0, 0, 0, 0; }
    vars(double vx, double vy, double sxx, double syy, double sxy) {
        *this << vx, vy, sxx, syy, sxy;
    }
    template<class Other>
    vars(const Other &o) : Eigen::Matrix<double, 5, 1>(o) { }
    template<class Other>
    vars &operator=(const Other &o) {
        Eigen::Matrix<double, 5, 1>::operator=(o);
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

struct cond {
    // nu < 0 -> nu = inf
    double nu;
    explicit cond(const double nu = -1) : nu(nu) { }
};

struct equations {
    using vars_t = vars;
    using param_t = param;
    using cond_t = cond;

    using vec_t = Eigen::Matrix<double, 5, 1>;
    using mat_t = Eigen::Matrix<double, 5, 5>;

    std::pair<vars_t, vars_t> riemman(dir d, cond cnd, const vars_t &UL, const vars_t &UR, const param_t &pL, const param_t &pR) const {
        if (d == dir::X) {
            const double vxL  = UL[0];
            const double vyL  = UL[1];
            const double sxxL = UL[2];
            //const double syyL = UL[3];
            const double sxyL = UL[4];

            const double rL = pL.rho;
            const double mL = pL.mu;
            const double lL = pL.lambda;
            const double l2mL = lL + 2*mL;
            const double srL = lL / l2mL;

            const double asL = std::sqrt(rL * mL);
            const double apL = std::sqrt(rL * l2mL);

            const double vxR  = UR[0];
            const double vyR  = UR[1];
            const double sxxR = UR[2];
            //const double syyR = UR[3];
            const double sxyR = UR[4];

            const double rR = pR.rho;
            const double mR = pR.mu;
            const double lR = pR.lambda;
            const double l2mR = lR + 2*mR;
            const double srR = lR / l2mR;

            const double asR = std::sqrt(rR * mR);
            const double apR = std::sqrt(rR * l2mR);

            // nu = nunum / nuden
            // nu = inf  -> nunum = 1, nuden = 0
            // nu != inf -> nunum = nu, nuden = 1
            double nunum, nuden;

            if (cnd.nu < 0) {
                nunum = 1;
                nuden = 0;
            } else {
                nunum = cnd.nu;
                nuden = 1;
            }

            const double c1 = (apL*(-sxxL + sxxR + apR*(-vxL + vxR)))/(apL + apR);
            const double c5 = (apR*( sxxL - sxxR + apL*(-vxL + vxR)))/(apL + apR);

            const double c2 = -(asL*(nunum*(sxyL - sxyR) + asR*(nuden*sxyL + nunum*(vyL - vyR))))/(asL*asR*nuden + (asL + asR)*nunum);
            const double c4 =  (asR*(nunum*(sxyL - sxyR - asL*vyL + asL*vyR) - (asL*nuden*sxyR)))/(asL*asR*nuden + (asL + asR)*nunum);

            vars_t ULs = UL + vars_t( c1 / apL,  c2 / asL, c1, c1 * srL, c2);
            vars_t URs = UR + vars_t(-c5 / apR, -c4 / asR, c5, c5 * srR, c4);

            vars_t fL = F(dir::X, ULs, pL);
            vars_t fR = F(dir::X, URs, pR);

            return std::make_pair(fL, fR);
        }

        vars_t ULrot(UL[1], -UL[0], UL[3], UL[2], -UL[4]);
        vars_t URrot(UR[1], -UR[0], UR[3], UR[2], -UR[4]);
        std::pair<vars_t, vars_t> flxs = riemman(dir::X, cnd, ULrot, URrot, pL, pR);
        const auto &fLrot = flxs.first;
        const auto &fRrot = flxs.second;

        vars_t fL(-fLrot[1], fLrot[0], fLrot[3], fLrot[2], -fLrot[4]);
        vars_t fR(-fRrot[1], fRrot[0], fRrot[3], fRrot[2], -fRrot[4]);

        return std::make_pair(fL, fR);
    }
    vars_t F(dir d, const vars_t &U, const param_t &p) const {
        const double vx = U[0];
        const double vy = U[1];
        const double sxx = U[2];
        const double syy = U[3];
        const double sxy = U[4];

        const double m = p.mu;
        const double l = p.lambda;
        const double l2m = l + 2*m;
        const double ir = 1. / p.rho;

        if (d == dir::X)
            return vars_t(-ir*sxx, -ir*sxy, -l2m*vx, -l*vx, -m*vy);
        else
            return vars_t(-ir*sxy, -ir*syy, -l*vy, -l2m*vy, -m*vx);
    }
    // The left eigenvectors of dF/dU by rows
    mat_t W(dir d, const vars_t &, const param_t &p) const {
        const double r = p.rho;
        const double m = p.mu;
        const double l = p.lambda;
        const double l2m = l + 2*m;

        const double sr = l / l2m;

        mat_t ret;

        const double as = std::sqrt(m * r);
        const double ap = std::sqrt(l2m * r);

        if (d == dir::X)
            ret <<
                0,   0,  -sr, 1,  0,
                0,   as,  0,  0,  1,
                0,  -as,  0,  0,  1,
                 ap, 0,   1,  0,  0,
                -ap, 0,   1,  0,  0
            ;
        else
            ret <<
                 0,   0,  1, -sr,  0,
                 as,  0,  0,  0,   1,
                -as,  0,  0,  0,   1,
                 0,   ap, 0,  1,   0,
                 0,  -ap, 0,  1,   0
            ;

        return ret;
    }
    // The inverse of W
    mat_t iW(dir d, const vars_t &, const param_t &p) const {
        const double r = p.rho;
        const double m = p.mu;
        const double l = p.lambda;
        const double l2m = l + 2*m;

        const double sr = l / l2m;

        mat_t ret;

        const double as = std::sqrt(m * r);
        const double ap = std::sqrt(l2m * r);

        if (d == dir::X)
            ret <<
                0,     0,      0,       0.5/ap,  -0.5/ap,
                0,     0.5/as, -0.5/as, 0,        0,
                0,     0,      0,       0.5,      0.5,
                1,     0,      0,       0.5*sr,   0.5*sr,
                0,     0.5,    0.5,     0,        0
            ;
        else
            ret <<
                0,     0.5/as, -0.5/as, 0,        0,
                0,     0,      0,       0.5/ap,  -0.5/ap,
                1,     0,      0,       0.5*sr,   0.5*sr,
                0,     0,      0,       0.5,      0.5,
                0,     0.5,    0.5,     0,        0
            ;

        return ret;
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
