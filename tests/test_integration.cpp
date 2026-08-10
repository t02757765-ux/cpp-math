#include <cassert>
#include <cmath>
#include <iostream>

#include "mathx/integration.hpp"
#include "mathx/constants.hpp"

int main() {
    const auto f = [](double x) { return x * x; };            // int = b^3/3 - a^3/3
    const auto g = [](double x) { return std::sin(x); };      // int_0^pi = 2
    const auto h = [](double x) { return std::exp(-x * x); }; // ~sqrt(pi)

    assert(std::fabs(mathx::integrate_simpson(f, 0.0, 3.0) - 9.0) < 1e-9);
    assert(std::fabs(mathx::integrate_simpson(g, 0.0, mathx::constants::pi) - 2.0) < 1e-9);
    assert(std::fabs(mathx::integrate_trapezoid(f, 0.0, 3.0, 10000) - 9.0) < 1e-5);
    assert(std::fabs(mathx::integrate_simpson_n(f, 0.0, 3.0, 10000) - 9.0) < 1e-9);
    assert(std::fabs(mathx::integrate_gauss_legendre(g, 0.0, mathx::constants::pi, 16) - 2.0) < 1e-12);
    assert(std::fabs(mathx::integrate_romberg(f, 0.0, 3.0) - 9.0) < 1e-9);
    // Gaussian integral
    assert(std::fabs(mathx::integrate_gauss_legendre(h, -6.0, 6.0, 32)
                     - mathx::constants::sqrt_pi) < 1e-9);

    std::cout << "integration tests passed\n";
    return 0;
}
