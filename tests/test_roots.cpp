#include <cassert>
#include <cmath>
#include <iostream>

#include "mathx/roots.hpp"
#include "mathx/constants.hpp"

int main() {
    const auto f = [](double x) { return x * x - 2.0; };   // roots +/- sqrt(2)
    const auto df = [](double x) { return 2.0 * x; };
    const double sqrt2 = 1.4142135623730951;

    assert(std::fabs(mathx::bisect(f, 1.0, 2.0) - sqrt2) < 1e-10);
    assert(std::fabs(mathx::regula_falsi(f, 1.0, 2.0) - sqrt2) < 1e-10);
    assert(std::fabs(mathx::secant(f, 1.0, 2.0) - sqrt2) < 1e-10);
    assert(std::fabs(mathx::newton(f, df, 1.5) - sqrt2) < 1e-10);
    assert(std::fabs(mathx::newton_fd(f, 1.5) - sqrt2) < 1e-10);
    assert(std::fabs(mathx::brent(f, 1.0, 2.0) - sqrt2) < 1e-10);

    // sin(x) = 0 near 3.0 -> pi
    const auto s = [](double x) { return std::sin(x); };
    assert(std::fabs(mathx::brent(s, 3.0, 4.0) - mathx::constants::pi) < 1e-10);
    // Newton on cos: root at pi/2
    const auto c = [](double x) { return std::cos(x); };
    const auto dc = [](double x) { return -std::sin(x); };
    assert(std::fabs(mathx::newton(c, dc, 1.5) - mathx::constants::pi / 2.0) < 1e-10);

    // Unbracketed -> NaN
    assert(std::isnan(mathx::bisect(f, 1.0, 1.1)));

    std::cout << "roots tests passed\n";
    return 0;
}
