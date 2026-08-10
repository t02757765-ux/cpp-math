#ifndef MATHX_ROOTS_HPP
#define MATHX_ROOTS_HPP

#include "export.hpp"
#include <functional>

namespace mathx {

//! Bisection method. Requires f(a) and f(b) to have opposite signs.
//! Returns an approximate root in [a, b].
MATHX_API double bisect(const std::function<double(double)>& f,
                        double a, double b,
                        double tol = 1e-12, int max_iter = 100);

//! Regula falsi (false position) method.
MATHX_API double regula_falsi(const std::function<double(double)>& f,
                              double a, double b,
                              double tol = 1e-12, int max_iter = 100);

//! Secant method. Requires two starting values x0, x1.
MATHX_API double secant(const std::function<double(double)>& f,
                        double x0, double x1,
                        double tol = 1e-12, int max_iter = 100);

//! Newton-Raphson method using the provided derivative.
MATHX_API double newton(const std::function<double(double)>& f,
                        const std::function<double(double)>& df,
                        double x0,
                        double tol = 1e-12, int max_iter = 100);

//! Newton's method with a finite-difference derivative.
MATHX_API double newton_fd(const std::function<double(double)>& f,
                           double x0,
                           double tol = 1e-12, int max_iter = 100);

//! Brent's method (robust, combines bisection, secant and inverse
//! quadratic interpolation). Requires bracketing.
MATHX_API double brent(const std::function<double(double)>& f,
                       double a, double b,
                       double tol = 1e-12, int max_iter = 100);

} // namespace mathx

#endif
