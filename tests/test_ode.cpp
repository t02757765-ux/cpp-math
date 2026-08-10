#include <cassert>
#include <cmath>
#include <iostream>

#include "mathx/ode.hpp"
#include "mathx/constants.hpp"

int main() {
    // dy/dt = y, y(0)=1 -> y(t) = e^t
    const auto exp_ode = [](double, double y) { return y; };

    const auto sol_rk4 = mathx::ode_rk4(exp_ode, 1.0, 0.0, 1.0, 1000);
    assert(std::fabs(sol_rk4.back().second - std::exp(1.0)) < 1e-6);

    const auto sol_rk45 = mathx::ode_rk45(exp_ode, 1.0, 0.0, 1.0, 1e-10);
    assert(std::fabs(sol_rk45.back().second - std::exp(1.0)) < 1e-8);

    const auto sol_heun = mathx::ode_heun(exp_ode, 1.0, 0.0, 1.0, 10000);
    assert(std::fabs(sol_heun.back().second - std::exp(1.0)) < 1e-6);

    const auto sol_euler = mathx::ode_euler(exp_ode, 1.0, 0.0, 1.0, 10000);
    assert(std::fabs(sol_euler.back().second - std::exp(1.0)) < 1e-3);

    // dy/dt = -y, y(0)=1 -> y(t) = e^{-t}; RK45 must handle negative too
    const auto dec_ode = [](double, double y) { return -y; };
    const auto sol_dec = mathx::ode_rk45(dec_ode, 1.0, 0.0, 2.0, 1e-10);
    assert(std::fabs(sol_dec.back().second - std::exp(-2.0)) < 1e-8);

    // sin(t): y' = cos(t), y(0)=0 -> y(t)=sin(t)
    const auto sin_ode = [](double t, double) { return std::cos(t); };
    const auto sol_sin = mathx::ode_rk45(sin_ode, 0.0, 0.0, mathx::constants::pi / 2.0, 1e-10);
    assert(std::fabs(sol_sin.back().second - 1.0) < 1e-8);

    std::cout << "ode tests passed\n";
    return 0;
}
