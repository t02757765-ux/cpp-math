#include "mathx/hypergeometric.hpp"

#include <cmath>
#include <limits>

namespace mathx {

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

double pochhammer(double a, int k) {
    double r = 1.0;
    for (int i = 0; i < k; ++i) r *= (a + static_cast<double>(i));
    return r;
}

//! Arithmetic-geometric mean.
double agm(double a, double b) {
    for (int i = 0; i < 100; ++i) {
        const double na = 0.5 * (a + b);
        const double nb = std::sqrt(a * b);
        if (std::fabs(a - na) < 1e-16 * std::max(1.0, na)) break;
        a = na;
        b = nb;
    }
    return 0.5 * (a + b);
}

} // namespace

double hypergeometric_pfq(const double* a, int pa,
                          const double* b, int pb,
                          double z, int max_terms) {
    // Convergence check: p <= q + 1 needed. p = q + 1 needs |z| < 1.
    if (pa > pb + 1) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (pa == pb + 1 && std::fabs(z) >= 1.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    // Check for zero denominators in b.
    for (int i = 0; i < pb; ++i) {
        if (b[i] == 0.0) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        if (b[i] == std::floor(b[i]) && b[i] <= 0.0) {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }
    // b is non-positive integer -> also divergent/undefined.
    double term = 1.0;
    double sum = 0.0;
    for (int k = 0; k < max_terms; ++k) {
        sum += term;
        if (std::fabs(term) < 1e-16 * std::max(1.0, std::fabs(sum))) {
            return sum;
        }
        double num = z;
        for (int i = 0; i < pa; ++i) num *= (a[i] + static_cast<double>(k));
        double den = static_cast<double>(k + 1);
        for (int i = 0; i < pb; ++i) den *= (b[i] + static_cast<double>(k));
        if (den == 0.0) {
            return sum; // exact termination if a hits non-positive int
        }
        term *= num / den;
        if (!std::isfinite(term)) break;
    }
    // Series truncated without converging.
    if (std::fabs(term) > 1e-8 * std::fabs(sum)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return sum;
}

double hypergeometric_1f1(double a, double b, double z) {
    if (z < 0.0) {
        // Kummer transformation: 1F1(a;b;z) = e^z * 1F1(b-a;b;-z)
        // keeps the series better behaved for large negative z.
        return std::exp(z) * hypergeometric_1f1(b - a, b, -z);
    }
    const double a1[1] = {a};
    const double b1[1] = {b};
    return hypergeometric_pfq(a1, 1, b1, 1, z);
}

double hypergeometric_2f1(double a, double b, double c, double z) {
    // For z close to 1 use Gauss's theorem when c - a - b > 0.
    if (z > 0.9 && z < 1.0 && c > a + b) {
        // 2F1(a,b;c;1) = Gamma(c)Gamma(c-a-b)/(Gamma(c-a)Gamma(c-b))
        // We do not have Gamma in this translation unit; fall back to series.
    }
    const double a2[2] = {a, b};
    const double b1[1] = {c};
    return hypergeometric_pfq(a2, 2, b1, 1, z);
}

double elliptic_k(double m) {
    if (m < 0.0 || m > 1.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (m == 1.0) return std::numeric_limits<double>::infinity();
    if (m == 0.0) return kPi / 2.0;
    // K(m) = pi / (2 * AGM(1, sqrt(1 - m)))
    return kPi / (2.0 * agm(1.0, std::sqrt(1.0 - m)));
}

double elliptic_e(double m) {
    if (m < 0.0 || m > 1.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (m == 0.0) return kPi / 2.0;
    if (m == 1.0) return 1.0;
    // E(m) = (pi/2) * sum_n [(1/2)_n]^2/(n!^2) * m^n / (1 - 2n)
    double u = 1.0; // u_n for n = 0
    double sum = 1.0; // term n=0: 1/(1) = 1
    for (int n = 1; n <= 2000; ++n) {
        const double r = (static_cast<double>(n) - 0.5) / static_cast<double>(n);
        u *= r * r * m;
        const double term = u / (1.0 - 2.0 * static_cast<double>(n));
        sum += term;
        if (std::fabs(term) < 1e-16 * std::max(1.0, std::fabs(sum))) {
            return 0.5 * kPi * sum;
        }
    }
    return 0.5 * kPi * sum;
}

} // namespace mathx
