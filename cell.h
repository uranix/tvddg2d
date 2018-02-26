#pragma once

#include "matrix.h"

template<typename vars, int p>
struct cell {
    using mat_t = matrix<p>;
	vars U[p+3][p+3];
    void extrapolate() {
        for (int i = 0; i <= p; i++) {
            vars v1 = mat_t::mu[0] * access(i, 0);
            vars v2 = mat_t::mu[0] * access(i, p);
            for (int j = 1; j <= p; j++) {
                v1 += mat_t::mu[j] * access(i,   j);
                v2 += mat_t::mu[j] * access(i, p-j);
            }
            access(i, p+1) = v1;
            access(i,  -1) = v2;
        }

        for (int j = -1; j <= p+1; j++) {
            vars v3 = mat_t::mu[0] * access(0, j);
            vars v4 = mat_t::mu[0] * access(p, j);
            for (int i = 1; i <= p; i++) {
                v3 += mat_t::mu[i] * access(  i, j);
                v4 += mat_t::mu[i] * access(p-i, j);
            }
            access(p+1, j) = v3;
            access( -1, j) = v4;
        }
    }
    const vars &access(int i, int j) const {
        return U[i+1][j+1];
    }
    vars &access(int i, int j) {
        return U[i+1][j+1];
    }
    const vars &operator()(int i, int j) const {
#ifdef NDEBUG
        bool outside = (i < -1 || i > p+1) || (j < -1 || j > p+1);
        if (outside)
            throw std::out_of_range("Cell subindex out of range: (" +
                    std::to_string(i) + ", " + std::to_string(j) + ")");
#endif
        return access(i, j);
    }
    vars &operator()(int i, int j) {
        return const_cast<vars &>(
                static_cast<const cell<vars, p> &>(*this)(i, j)
            );
    }
};

/*
  ^ j
  |
    +-G-+-G-+-G-+
  p F   F   F   F
    +-G-+-G-+-G-+
    F   F   F   F
    +-G-+-G-+-G-+
  0 F   F   F   F
    +-G-+-G-+-G-+
      0       p   -> i
*/
template<typename vars, int p>
struct flux_cell {
    vars F[p+2][p+1];
    vars G[p+1][p+2];
};
