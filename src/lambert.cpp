#include "mathx/lambert.hpp"

#include <cmath>
#include <complex>
#include <limits>

namespace mathx {

namespace {

constexpr double kInvE = 0.3678794411714423215955237701614608;

//! Newton iteration for real Lambert W.
//! Converges for both branches when started from a reasonable guess.
double newton_w(double x, double w) {
    for (int i = 0; i < 100; ++i) {
        const double ew = std::exp(w);
        const double f  = w * ew - x;
        const double df = ew * (w + 1.0);
        if (df == 0.0) break;
        const double dw = f / df;
        w -= dw;
        if (std::fabs(dw) < 1e-14 * std::max(1.0, std::fabs(w))) break;
    }
    return w;
}

} // namespace

double lambert_w(double x, int branch) {
    if (x < -kInvE) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (std::fabs(x) < 1e-16) {
        return branch == 0 ? x : -std::numeric_limits<double>::infinity();
    }

    double w;
    if (branch == 0) {
        if (x < 0.0) {
            // Expansion in p = sqrt(2(e x + 1)) around the branch point.
            const double p = std::sqrt(2.0 * (x * std::exp(1.0) + 1.0));
            w = -1.0 + p - p * p / 3.0 + (11.0 / 72.0) * p * p * p;
        } else if (x < 1.0) {
            w = x;
        } else if (x < 14.0) {
            w = std::log(x);
        } else {
            const double l = std::log(x);
            w = l - std::log(l);
        }
    } else if (branch == -1) {
        if (x >= 0.0) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        // Good across the whole (-1/e, 0) interval for the lower branch.
        w = std::log(-x) - std::log(-std::log(-x));
    } else {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return newton_w(x, w);
}

double lambert_w_solve(double w, int branch) {
    return lambert_w(w, branch);
}

std::complex<double> lambert_w(const std::complex<double>& z) {
    using cplx = std::complex<double>;
    // Initial guess for the principal branch.
    cplx w;
    if (std::abs(z) < 0.7) {
        w = z;
    } else {
        cplx lz = std::log(z);
        w = lz - std::log(lz);
    }
    // Newton iteration in the complex plane.
    for (int i = 0; i < 200; ++i) {
        cplx ew = std::exp(w);
        cplx f  = w * ew - z;
        cplx df = ew * (w + 1.0);
        if (std::abs(df) == 0.0) break;
        cplx dw = f / df;
        w -= dw;
        if (std::abs(dw) < 1e-14 * std::max(1.0, std::abs(w))) break;
    }
    return w;
}

} // namespace mathx
