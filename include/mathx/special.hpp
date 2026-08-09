#ifndef MATHX_SPECIAL_HPP
#define MATHX_SPECIAL_HPP

#include "export.hpp"
#include <complex>

namespace mathx {

//! Gamma function (real argument).
MATHX_API double gamma(double x);

//! Log-Gamma function (real argument).
MATHX_API double lgamma(double x);

//! Digamma function (real argument).
MATHX_API double digamma(double x);

//! Error function.
MATHX_API double erf(double x);

//! Complementary error function.
MATHX_API double erfc(double x);

//! Gaussian (normal) cumulative distribution function.
MATHX_API double norm_cdf(double x);

//! Riemann zeta function (real argument).
MATHX_API double zeta(double s);

//! Beta function B(x, y).
MATHX_API double beta(double x, double y);

} // namespace mathx

#endif
