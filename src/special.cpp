#include "mathx/special.hpp"

#include <cmath>
#include <limits>

namespace mathx {

namespace {
constexpr double kPi = 3.14159265358979323846;
constexpr double kEulerMascheroni = 0.57721566490153286060;
constexpr double kSqrtPi = 1.77245385090551602729;

//! Lanczos approximation coefficients (g = 7).
constexpr double kLanczos[9] = {
    0.99999999999980993,
    676.5203681218851,
    -1259.1392167224028,
    771.32342877765313,
    -176.61502916214059,
    12.507343278686905,
    -0.13857109526572012,
    9.9843695780195716e-6,
    1.5056327351493116e-7,
};

double lanczos_gamma(double x) {
    if (x < 0.5) {
        // Reflection formula.
        const double s = kPi / std::sin(kPi * x);
        return s / lanczos_gamma(1.0 - x);
    }
    x -= 1.0;
    double a = kLanczos[0];
    for (int i = 1; i < 9; ++i) {
        a += kLanczos[i] / (x + static_cast<double>(i));
    }
    const double t = x + 7.5;
    return std::sqrt(2.0 * kPi) * std::pow(t, x + 0.5) * std::exp(-t) * a;
}
} // namespace

double gamma(double x) {
    if (std::isnan(x)) return x;
    if (std::isinf(x)) return x > 0 ? x : std::numeric_limits<double>::quiet_NaN();
    // Poles at non-positive integers.
    if (x <= 0.0 && x == std::floor(x)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return lanczos_gamma(x);
}

double lgamma(double x) {
    if (x <= 0.0 && x == std::floor(x)) {
        return std::numeric_limits<double>::infinity();
    }
    return std::log(std::fabs(gamma(x)));
}

double digamma(double x) {
    if (x <= 0.0 && x == std::floor(x)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    double d = 0.0;
    // Shift for large x: psi(x) = psi(x + n) - sum 1/(x+k).
    while (x < 8.0) {
        d -= 1.0 / x;
        x += 1.0;
    }
    // Asymptotic series.
    const double inv = 1.0 / x;
    const double inv2 = inv * inv;
    d += std::log(x) - 0.5 * inv
       - inv2 * (1.0 / 12.0)
       + inv2 * inv2 * (1.0 / 120.0)
       - inv2 * inv2 * inv2 * (1.0 / 252.0);
    return d;
}

double erf(double x) {
    if (std::fabs(x) < 1e-8) {
        // erf(x) ~ 2x/sqrt(pi), exact at 0.
        return (2.0 / kSqrtPi) * x;
    }
    // Abramowitz & Stegun 7.1.26 (max error 1.5e-7).
    const double t = 1.0 / (1.0 + 0.3275911 * std::fabs(x));
    const double y = 1.0 - (((((1.061405429 * t - 1.453152027) * t)
        + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t
        * std::exp(-x * x);
    return x >= 0.0 ? y : -y;
}

double erfc(double x) {
    return 1.0 - erf(x);
}

double norm_cdf(double x) {
    return 0.5 * (1.0 + erf(x / std::sqrt(2.0)));
}

double beta(double x, double y) {
    return gamma(x) * gamma(y) / gamma(x + y);
}

double zeta(double s) {
    if (s == 1.0) {
        return std::numeric_limits<double>::infinity();
    }
    if (s <= 0.0) {
        // Functional equation: zeta(s) = 2^s pi^(s-1) sin(pi s / 2) Gamma(1 - s) zeta(1 - s)
        const double term = std::pow(2.0, s) * std::pow(kPi, s - 1.0)
                          * std::sin(0.5 * kPi * s) * gamma(1.0 - s);
        return term * zeta(1.0 - s);
    }
    // Euler-Maclaurin summation for the right half-plane.
    constexpr int N = 100;
    double sum = 0.0;
    for (int n = 1; n < N; ++n) {
        sum += std::pow(static_cast<double>(n), -s);
    }
    const double xN = static_cast<double>(N);
    const double term1 = std::pow(xN, 1.0 - s) / (s - 1.0);
    const double term2 = 0.5 * std::pow(xN, -s);
    const double term3 = std::pow(xN, -1.0 - s) * s / 12.0;
    return sum + term1 + term2 + term3;
}

} // namespace mathx
