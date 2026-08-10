#ifndef MATHX_INTEGRATION_HPP
#define MATHX_INTEGRATION_HPP

#include "export.hpp"
#include <functional>

namespace mathx {

//! Adaptive Simpson quadrature.
//! Returns an approximation of int_a^b f(x) dx with estimated absolute
//! error below tol (clamped to a minimum based on machine precision).
MATHX_API double integrate_simpson(
    const std::function<double(double)>& f,
    double a, double b,
    double tol = 1e-10,
    int max_depth = 50);

//! Composite trapezoidal rule with n subintervals.
MATHX_API double integrate_trapezoid(
    const std::function<double(double)>& f,
    double a, double b, int n = 1000);

//! Composite Simpson's rule with n subintervals (n must be even).
MATHX_API double integrate_simpson_n(
    const std::function<double(double)>& f,
    double a, double b, int n = 1000);

//! Gauss-Legendre quadrature with n nodes (5, 10, 16 or 32 supported).
MATHX_API double integrate_gauss_legendre(
    const std::function<double(double)>& f,
    double a, double b, int n = 16);

//! Romberg integration to given tolerance.
MATHX_API double integrate_romberg(
    const std::function<double(double)>& f,
    double a, double b,
    double tol = 1e-10);

} // namespace mathx

#endif
