#ifndef MATHX_LAMBERT_HPP
#define MATHX_LAMBERT_HPP

#include "export.hpp"
#include <complex>

namespace mathx {

//! Lambert W function on the real axis.
//! branch = 0  -> principal branch W0(x), x >= -1/e
//! branch = -1 -> lower branch W_{-1}(x), -1/e <= x < 0
MATHX_API double lambert_w(double x, int branch = 0);

//! Lambert W function for complex arguments (principal branch W0).
MATHX_API std::complex<double> lambert_w(const std::complex<double>& z);

//! Inverse: solves x * exp(x) = w  (equivalent to lambert_w).
MATHX_API double lambert_w_solve(double w, int branch = 0);

} // namespace mathx

#endif
