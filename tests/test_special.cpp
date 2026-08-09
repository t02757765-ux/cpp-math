#include <cassert>
#include <cmath>
#include <iostream>

#include "mathx/special.hpp"

namespace {
constexpr double kPi = 3.14159265358979323846;
constexpr double kSqrtPi = 1.77245385090551602729;
}

int main() {
    // gamma(n) = (n-1)!
    assert(std::fabs(mathx::gamma(1.0) - 1.0) < 1e-9);
    assert(std::fabs(mathx::gamma(2.0) - 1.0) < 1e-9);
    assert(std::fabs(mathx::gamma(5.0) - 24.0) < 1e-8);
    assert(std::fabs(mathx::gamma(0.5) - kSqrtPi) < 1e-9);

    // erf (A&S approx, max error ~1.5e-7)
    assert(std::fabs(mathx::erf(0.0)) < 1e-12);
    assert(std::fabs(mathx::erf(std::numeric_limits<double>::infinity()) - 1.0) < 1e-7);
    assert(std::fabs(mathx::erf(1.0) - 0.8427007929497149) < 2e-7);

    // erfc(0) = 1
    assert(std::fabs(mathx::erfc(0.0) - 1.0) < 1e-12);

    // norm_cdf(0) = 0.5
    assert(std::fabs(mathx::norm_cdf(0.0) - 0.5) < 1e-12);

    // zeta(2) = pi^2/6
    assert(std::fabs(mathx::zeta(2.0) - (kPi * kPi / 6.0)) < 1e-9);
    // zeta(4) = pi^4/90
    assert(std::fabs(mathx::zeta(4.0) - (kPi * kPi * kPi * kPi / 90.0)) < 1e-8);

    // beta(1,1) = 1
    assert(std::fabs(mathx::beta(1.0, 1.0) - 1.0) < 1e-9);
    // beta(0.5, 0.5) = pi
    assert(std::fabs(mathx::beta(0.5, 0.5) - kPi) < 1e-8);

    // digamma(1) = -gamma
    assert(std::fabs(mathx::digamma(1.0) + 0.5772156649015329) < 1e-9);

    std::cout << "special tests passed\n";
    return 0;
}
