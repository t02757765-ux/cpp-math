#include "mathx/bessel.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace mathx {

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kEulerGamma = 0.577215664901532860606512090082402431;
constexpr double kSqrtPi = 1.772453850905516027298167483341145183;

//! Series for J_nu(x), valid for moderate |x|.
double j_series(double nu, double x) {
    const double x2 = x * x;
    double term = std::pow(0.5 * x, nu);
    for (int i = 0; i < nu; ++i) {
        term /= static_cast<double>(i + 1); // divide by Gamma(nu+1) when nu integer
    }
    double sum = 0.0;
    double t = term;
    int k = 0;
    while (std::fabs(t) > 1e-18 * std::max(1.0, std::fabs(sum))) {
        sum += t;
        ++k;
        t *= -x2 / (4.0 * static_cast<double>(k) * (nu + static_cast<double>(k)));
        if (k > 200) break;
    }
    return sum;
}

//! Asymptotic expansion of J_0 / J_1 for large x.
void j_asymptotic(double x, double& j0, double& j1) {
    const double t = 8.0 * x;
    const double p0 = 1.0 - 9.0 / (2.0 * t * t) + 3675.0 / (8.0 * t * t * t * t);
    const double q0 = -1.0 / t + 75.0 / (8.0 * t * t * t);
    const double p1 = 1.0 + 15.0 / (2.0 * t * t) - 14175.0 / (8.0 * t * t * t * t);
    const double q1 = 3.0 / t - 525.0 / (8.0 * t * t * t);
    const double a = x - kPi / 4.0;
    const double ca = std::cos(a);
    const double sa = std::sin(a);
    const double s = std::sqrt(2.0 / (kPi * x));
    j0 = s * (ca * p0 - sa * q0);
    j1 = s * (ca * p1 - sa * q1);
}

//! Series for I_nu(x).
double i_series(double nu, double x) {
    const double x2 = x * x;
    double term = std::pow(0.5 * x, nu);
    for (int i = 0; i < nu; ++i) term /= static_cast<double>(i + 1);
    double sum = 0.0;
    double t = term;
    int k = 0;
    while (std::fabs(t) > 1e-18 * std::max(1.0, std::fabs(sum))) {
        sum += t;
        ++k;
        t *= x2 / (4.0 * static_cast<double>(k) * (nu + static_cast<double>(k)));
        if (k > 300) break;
    }
    return sum;
}

//! Asymptotic expansion for I_0, I_1 (large x).
void i_asymptotic(double x, double& i0, double& i1) {
    const double t = 8.0 * x;
    const double p0 = 1.0 + 1.0 / (2.0 * t) + 9.0 / (8.0 * t * t) + 75.0 / (8.0 * t * t * t);
    const double p1 = 1.0 - 3.0 / (2.0 * t) - 15.0 / (8.0 * t * t) - 105.0 / (8.0 * t * t * t);
    const double e = std::exp(x) / std::sqrt(2.0 * kPi * x);
    i0 = e * p0;
    i1 = e * p1;
}

//! Asymptotic expansion for K_0, K_1 (large x).
void k_asymptotic(double x, double& k0, double& k1) {
    const double t = 8.0 * x;
    const double p0 = 1.0 - 1.0 / (2.0 * t) + 9.0 / (8.0 * t * t) - 75.0 / (8.0 * t * t * t);
    const double p1 = 1.0 + 3.0 / (2.0 * t) + 15.0 / (8.0 * t * t) + 105.0 / (8.0 * t * t * t);
    const double e = std::exp(-x) * std::sqrt(kPi / (2.0 * x));
    k0 = e * p0;
    k1 = e * p1;
}

double y0_from_j(double x, double j0, double j1) {
    // Y0(x) = 2/pi [ J0(x) (ln(x/2) + gamma) - S ], S = sum_k (-1)^k H_k (x/2)^{2k}/(k!)^2
    (void)j1;
    const double xh = x * 0.5;
    double sum = 0.0;
    double term = 1.0; // k=0 term (H_0 = 0, contributes nothing)
    double hk = 0.0;
    for (int k = 1; k <= 60; ++k) {
        term *= -(xh * xh) / (static_cast<double>(k * k)); // (-1)^k (xh)^(2k)/(k!)^2
        hk += 1.0 / static_cast<double>(k);
        sum += term * hk;
    }
    return (2.0 / kPi) * (j0 * (std::log(xh) + kEulerGamma) - sum);
}

} // namespace

double bessel_j0(double x) {
    x = std::fabs(x);
    if (x < 8.0) return j_series(0.0, x);
    double a, b;
    j_asymptotic(x, a, b);
    return a;
}

double bessel_j1(double x) {
    const double ax = std::fabs(x);
    double r;
    if (ax < 8.0) {
        r = j_series(1.0, ax);
    } else {
        double a, b;
        j_asymptotic(ax, a, b);
        r = b;
    }
    return x < 0.0 ? -r : r;
}

double bessel_jn(int n, double x) {
    if (n < 0) {
        const int m = -n;
        const double r = bessel_jn(m, x);
        return (m % 2 == 0) ? r : -r;
    }
    if (n == 0) return bessel_j0(x);
    if (n == 1) return bessel_j1(x);
    if (std::fabs(x) < 8.0) return j_series(static_cast<double>(n), x);
    // Backward recurrence.
    double j0 = bessel_j0(x);
    double j1 = bessel_j1(x);
    for (int k = 1; k < n; ++k) {
        const double jk = (2.0 * static_cast<double>(k) / x) * j1 - j0;
        j0 = j1;
        j1 = jk;
    }
    return j1;
}

double bessel_y0(double x) {
    if (x <= 0.0) return std::numeric_limits<double>::quiet_NaN();
    if (x < 8.0) return y0_from_j(x, bessel_j0(x), bessel_j1(x));
    // Asymptotic.
    const double a = x - kPi / 4.0;
    const double p0 = 1.0 - 9.0 / (2.0 * (8.0 * x) * (8.0 * x));
    const double q0 = -1.0 / (8.0 * x);
    const double s = std::sqrt(2.0 / (kPi * x));
    return s * (std::sin(a) * p0 + std::cos(a) * q0);
}

double bessel_y1(double x) {
    if (x <= 0.0) return std::numeric_limits<double>::quiet_NaN();
    if (x < 8.0) {
        // Y1(x) = 2/pi [ J1(x) ln(x/2) + ... ] - 2/(pi x) - (x/pi) sum ...
        const double xh = x * 0.5;
        const double j1 = bessel_j1(x);
        double sum = 0.0;
        double term = xh * 0.5;
        double hk = 0.0;
        for (int k = 1; k <= 60; ++k) {
            hk += 1.0 / static_cast<double>(k);
            sum += term * (hk + (hk - 1.0 / static_cast<double>(k)));
            term *= -(xh * xh) / (static_cast<double>(k) * static_cast<double>(k + 1));
        }
        return (2.0 / kPi) * (j1 * (std::log(xh) + kEulerGamma) - 1.0 / x - sum);
    }
    const double a = x - kPi / 4.0;
    const double s = std::sqrt(2.0 / (kPi * x));
    const double p1 = 1.0 + 15.0 / (2.0 * (8.0 * x) * (8.0 * x));
    const double q1 = 3.0 / (8.0 * x);
    return s * (std::sin(a) * p1 + std::cos(a) * q1);
}

double bessel_yn(int n, double x) {
    if (n < 0) {
        const int m = -n;
        const double r = bessel_yn(m, x);
        return (m % 2 == 0) ? r : -r;
    }
    if (n == 0) return bessel_y0(x);
    if (n == 1) return bessel_y1(x);
    double y0 = bessel_y0(x);
    double y1 = bessel_y1(x);
    for (int k = 1; k < n; ++k) {
        const double yk = (2.0 * static_cast<double>(k) / x) * y1 - y0;
        y0 = y1;
        y1 = yk;
    }
    return y1;
}

double bessel_i0(double x) {
    x = std::fabs(x);
    if (x < 16.0) return i_series(0.0, x);
    double a, b;
    i_asymptotic(x, a, b);
    return a;
}

double bessel_i1(double x) {
    const double ax = std::fabs(x);
    double r;
    if (ax < 16.0) {
        r = i_series(1.0, ax);
    } else {
        double a, b;
        i_asymptotic(ax, a, b);
        r = b;
    }
    return x < 0.0 ? -r : r;
}

double bessel_in(int n, double x) {
    if (n < 0) {
        const int m = -n;
        const double r = bessel_in(m, x);
        return (m % 2 == 0) ? r : -r;
    }
    if (n == 0) return bessel_i0(x);
    if (n == 1) return bessel_i1(x);
    if (std::fabs(x) < 16.0) return i_series(static_cast<double>(n), x);
    double i0 = bessel_i0(x);
    double i1 = bessel_i1(x);
    for (int k = 1; k < n; ++k) {
        const double ik = (2.0 * static_cast<double>(k) / x) * i1 + i0;
        i0 = i1;
        i1 = ik;
    }
    return i1;
}

double bessel_k0(double x) {
    if (x <= 0.0) return std::numeric_limits<double>::quiet_NaN();
    if (x < 2.0) {
        // K0(x) = -[ln(x/2) + gamma] I0(x) + sum_k H_k (x/2)^{2k}/(k!)^2
        const double xh = x * 0.5;
        const double i0 = bessel_i0(x);
        double sum = 0.0;
        double term = 1.0;
        double hk = 0.0;
        for (int k = 1; k <= 60; ++k) {
            term *= (xh * xh) / (static_cast<double>(k * k));
            hk += 1.0 / static_cast<double>(k);
            sum += term * hk;
        }
        return -(std::log(xh) + kEulerGamma) * i0 + sum;
    }
    double a, b;
    k_asymptotic(x, a, b);
    return a;
}

double bessel_k1(double x) {
    if (x <= 0.0) return std::numeric_limits<double>::quiet_NaN();
    if (x < 2.0) {
        // K1(x) = 1/x + ln(x/2) I1(x) - (x/4) sum_k [psi(k+1)+psi(k+2)] (x/2)^{2k}/(k!(k+1)!)
        const double xh = x * 0.5;
        const double i1 = bessel_i1(x);
        double sum = 0.0;
        double term = 1.0; // k=0 term
        double hk = 0.0;   // H_{k+1}
        for (int k = 0; k <= 60; ++k) {
            const double psi1 = hk - kEulerGamma;          // psi(k+1)
            const double psi2 = (hk + 1.0 / (k + 1.0)) - kEulerGamma; // psi(k+2)
            if (k == 0) {
                sum += term * (psi1 + psi2);
                term *= (xh * xh) / (1.0 * 2.0);
                hk = 1.0;
            } else {
                sum += term * (psi1 + psi2);
                term *= (xh * xh) / (static_cast<double>(k + 1) * static_cast<double>(k + 2));
                hk += 1.0 / static_cast<double>(k + 1);
            }
        }
        return 1.0 / x + std::log(xh) * i1 - (x * 0.25) * sum;
    }
    double a, b;
    k_asymptotic(x, a, b);
    return b;
}

double bessel_kn(int n, double x) {
    if (n < 0) return bessel_kn(-n, x); // K_n symmetric in n
    if (n == 0) return bessel_k0(x);
    if (n == 1) return bessel_k1(x);
    double k0 = bessel_k0(x);
    double k1 = bessel_k1(x);
    for (int k = 1; k < n; ++k) {
        const double kk = (2.0 * static_cast<double>(k) / x) * k1 + k0;
        k0 = k1;
        k1 = kk;
    }
    return k1;
}

double airy_ai(double x) {
    // Series expansion of Ai(x) around 0, valid for moderate |x|.
    const double c0 = 0.355028053887817;
    const double c1 = 0.258819403792807;
    const double x3 = x * x * x;
    double s0 = c0, s1 = c1, t0 = c0, t1 = c1;
    for (int k = 1; k <= 60; ++k) {
        t0 *= x3 / (3.0 * static_cast<double>(k) * (3.0 * static_cast<double>(k) - 1.0));
        t1 *= x3 / (3.0 * static_cast<double>(k) * (3.0 * static_cast<double>(k) + 1.0));
        s0 += t0;
        s1 += t1;
        if (std::fabs(t0) < 1e-18 && std::fabs(t1) < 1e-18) break;
    }
    return s0 - s1 * x;
}

double airy_bi(double x) {
    // Series expansion of Bi(x) around 0, valid for moderate |x|.
    const double c0 = 0.614926627446001;
    const double c1 = 0.448288357353826;
    const double x3 = x * x * x;
    double s0 = c0, s1 = c1, t0 = c0, t1 = c1;
    for (int k = 1; k <= 60; ++k) {
        t0 *= x3 / (3.0 * static_cast<double>(k) * (3.0 * static_cast<double>(k) - 1.0));
        t1 *= x3 / (3.0 * static_cast<double>(k) * (3.0 * static_cast<double>(k) + 1.0));
        s0 += t0;
        s1 += t1;
        if (std::fabs(t0) < 1e-18 && std::fabs(t1) < 1e-18) break;
    }
    return s0 + s1 * x;
}

} // namespace mathx
