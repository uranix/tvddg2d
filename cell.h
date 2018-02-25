#pragma once

#ifndef SANITY_CHECKS
# define SANITY_CHECKS 0
#endif

#ifndef NDEBUG
# undef SANITY_CHECKS
# define SANITY_CHECKS 1
#endif

#if SANITY_CHECKS
# warning "Sanity checks are enabled"
#endif

#include "matrix.h"

template<typename vars, int _p>
struct cell {
	static constexpr int p = _p;
    using mat_t = matrix<_p>;
#if SANITY_CHECKS
    bool extrapolated;
#endif
    cell() {
#if SANITY_CHECKS
        extrapolated = false;
#endif
    }
	vars U[_p+3][_p+3];
    void extrapolate() {
#if SANITY_CHECKS
        if (extrapolated)
            throw std::logic_error("double call to extrapolate");
#endif

        for (int i = 0; i <= _p; i++) {
            vars v1 = mat_t::mu[0] * access(i,  0);
            vars v2 = mat_t::mu[0] * access(i, _p);
            for (int j = 1; j <= _p; j++) {
                v1 += mat_t::mu[j] * access(i,    j);
                v2 += mat_t::mu[j] * access(i, _p-j);
            }
            access(i, _p+1) = v1;
            access(i,   -1) = v2;
        }

        for (int j = -1; j <= _p+1; j++) {
            vars v3 = mat_t::mu[0] * access( 0, j);
            vars v4 = mat_t::mu[0] * access(_p, j);
            for (int i = 1; i <= _p; i++) {
                v3 += mat_t::mu[i] * access(   i, j);
                v4 += mat_t::mu[i] * access(_p-i, j);
            }
            access(_p+1, j) = v3;
            access(  -1, j) = v4;
        }
#if SANITY_CHECKS
        extrapolated = true;
#endif
    }
    void touch() {
#if SANITY_CHECKS
        if (!extrapolated)
            throw std::logic_error("unnessesary call to touch");
#endif
#if SANITY_CHECKS
        extrapolated = false;
#endif
    }
    const vars &access(int i, int j) const {
        return U[i+1][j+1];
    }
    vars &access(int i, int j) {
        return U[i+1][j+1];
    }
    const vars &operator()(int i, int j) const {
#if SANITY_CHECKS
        bool inside = (i >= 0 && i <= _p) && (j >= 0 && j <= _p);
        if (!inside && !extrapolated)
            throw std::out_of_range("Data is not extrapolated to the bounds");
        bool outside = (i < -1 || i > _p+1) || (j < -1 || j > _p+1);
        if (outside)
            throw std::out_of_range("Cell subindex out of range: (" +
                    std::to_string(i) + ", " + std::to_string(j) + ")");
#endif
        return access(i, j);
    }
    vars &operator()(int i, int j) {
        return const_cast<vars &>(
                static_cast<const cell<vars, _p> &>(*this)(i, j)
            );
    }
};
