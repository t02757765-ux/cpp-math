#include "mathx/lambert.hpp"

#include <cmath>
#include <complex>
#include <limits>

namespace mathx {

namespace {

constexpr double kInvE = 0.36787944117144233;

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

//! Series expansion around the branch point x = -1/e (Corless et al.).
//! With p = sqrt(2(e x + 1)) this gives W0(x) = -1 + p - p^2/3 + ...
//! Accurate to machine precision for small p; Newton is numerically
//! unstable at the double root (w = -1), so we skip it there.
double branch_point_w(double x, int branch) {
    // Compute p = sqrt(2(e x + 1)) stably.  Using t = e*(x + 1/e) avoids the
    // cancellation/fma artifact at x = -1/e, where x + kInvE is exact (Sterbenz).
    double t = std::exp(1.0) * (x + kInvE);
    if (t < 0.0) t = 0.0;  // guard against rounding at x = -1/e
    const double p = std::sqrt(2.0 * t);
    const double p2 = p * p;
    const double p3 = p2 * p;
    const double p4 = p2 * p2;
    const double p5 = p3 * p2;
    const double p6 = p3 * p3;
    const double p7 = p4 * p3;
    double s = 1.0;
    if (branch == -1) s = -1.0;
    // W = -1 + s*p - p^2/3 + s*(11/72)p^3 - (43/540)p^4 + ...
    double w = -1.0
        + s * p
        - p2 / 3.0
        + s * (11.0 / 72.0) * p3
        - (43.0 / 540.0) * p4
        + s * (769.0 / 17280.0) * p5
        - (221.0 / 8505.0) * p6
        + s * (680863.0 / 43545600.0) * p7;
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
    const double p = std::sqrt(2.0 * std::exp(1.0) * (x + kInvE)); // stable near branch point
    if (branch == 0) {
        if (x < 0.0) {
            // Newton is unstable near the double root (w = -1); use the
            // series directly when p is small.
            if (p < 1e-3) {
                return branch_point_w(x, 0);
            }
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
        // Also refine near the branch point where Newton breaks down.
        if (p < 1e-3) {
            return branch_point_w(x, -1);
        }
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
