#ifndef MATHX_ODE_HPP
#define MATHX_ODE_HPP

#include "export.hpp"
#include <functional>
#include <utility>
#include <vector>

namespace mathx {

//! First-order ODE: dy/dt = f(t, y), y(t0) = y0.
using OdeFunction = std::function<double(double t, double y)>;

//! Solve dy/dt = f(t, y) with the (explicit) Euler method.
//! Returns pairs (t, y). n is the number of steps.
MATHX_API std::vector<std::pair<double, double>> ode_euler(
    const OdeFunction& f, double y0,
    double t0, double t1, int n = 1000);

//! Solve with the classical 4th-order Runge-Kutta method (RK4).
MATHX_API std::vector<std::pair<double, double>> ode_rk4(
    const OdeFunction& f, double y0,
    double t0, double t1, int n = 1000);

//! Solve with an adaptive step-size RK45 (Dormand-Prince) to tolerance tol.
MATHX_API std::vector<std::pair<double, double>> ode_rk45(
    const OdeFunction& f, double y0,
    double t0, double t1, double tol = 1e-8);

//! Predictor-corrector (Heun's method).
MATHX_API std::vector<std::pair<double, double>> ode_heun(
    const OdeFunction& f, double y0,
    double t0, double t1, int n = 1000);

} // namespace mathx

#endif
