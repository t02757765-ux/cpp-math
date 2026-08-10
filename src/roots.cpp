#include "mathx/roots.hpp"

#include <cmath>
#include <limits>

namespace mathx {

namespace {
double rel_tol(double tol) {
    return std::max(tol, std::numeric_limits<double>::epsilon() * 16.0);
}
} // namespace

double bisect(const std::function<double(double)>& f,
              double a, double b, double tol, int max_iter) {
    double fa = f(a);
    double fb = f(b);
    if (fa == 0.0) return a;
    if (fb == 0.0) return b;
    if (fa * fb > 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double eps = rel_tol(tol);
    for (int i = 0; i < max_iter; ++i) {
        const double m = 0.5 * (a + b);
        const double fm = f(m);
        if (fm == 0.0 || (b - a) < 2.0 * eps * std::max(1.0, std::fabs(m))) {
            return m;
        }
        if (fa * fm < 0.0) {
            b = m;
            fb = fm;
        } else {
            a = m;
            fa = fm;
        }
    }
    return 0.5 * (a + b);
}

double regula_falsi(const std::function<double(double)>& f,
                    double a, double b, double tol, int max_iter) {
    double fa = f(a);
    double fb = f(b);
    if (fa == 0.0) return a;
    if (fb == 0.0) return b;
    if (fa * fb > 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double eps = rel_tol(tol);
    for (int i = 0; i < max_iter; ++i) {
        const double x = (a * fb - b * fa) / (fb - fa);
        const double fx = f(x);
        if (fx == 0.0 || (b - a) < eps * std::max(1.0, std::fabs(x))) {
            return x;
        }
        if (fa * fx < 0.0) {
            b = x;
            fb = fx;
        } else {
            a = x;
            fa = fx;
        }
    }
    return (a * fb - b * fa) / (fb - fa);
}

double secant(const std::function<double(double)>& f,
              double x0, double x1, double tol, int max_iter) {
    const double eps = rel_tol(tol);
    double f0 = f(x0);
    if (f0 == 0.0) return x0;
    for (int i = 0; i < max_iter; ++i) {
        const double f1 = f(x1);
        if (f1 == 0.0) return x1;
        const double denom = f1 - f0;
        if (denom == 0.0) return x1;
        const double x2 = x1 - f1 * (x1 - x0) / denom;
        if (std::fabs(x2 - x1) <= eps * std::max(1.0, std::fabs(x2))) {
            return x2;
        }
        x0 = x1;
        f0 = f1;
        x1 = x2;
    }
    return x1;
}

double newton(const std::function<double(double)>& f,
              const std::function<double(double)>& df,
              double x0, double tol, int max_iter) {
    const double eps = rel_tol(tol);
    double x = x0;
    for (int i = 0; i < max_iter; ++i) {
        const double fx = f(x);
        if (fx == 0.0) return x;
        const double dfx = df(x);
        if (dfx == 0.0) return x;
        const double dx = fx / dfx;
        x -= dx;
        if (std::fabs(dx) <= eps * std::max(1.0, std::fabs(x))) {
            return x;
        }
    }
    return x;
}

double newton_fd(const std::function<double(double)>& f,
                 double x0, double tol, int max_iter) {
    const double eps = rel_tol(tol);
    double x = x0;
    for (int i = 0; i < max_iter; ++i) {
        const double fx = f(x);
        if (fx == 0.0) return x;
        const double h = std::sqrt(std::numeric_limits<double>::epsilon())
                       * std::max(1.0, std::fabs(x));
        const double dfx = (f(x + h) - f(x - h)) / (2.0 * h);
        if (dfx == 0.0) return x;
        const double dx = fx / dfx;
        x -= dx;
        if (std::fabs(dx) <= eps * std::max(1.0, std::fabs(x))) {
            return x;
        }
    }
    return x;
}

double brent(const std::function<double(double)>& f,
             double a, double b, double tol, int max_iter) {
    const double eps = rel_tol(tol);
    double fa = f(a);
    double fb = f(b);
    if (fa == 0.0) return a;
    if (fb == 0.0) return b;
    if (fa * fb > 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    double c = a, fc = fa, d = b - a, e = d;
    for (int i = 0; i < max_iter; ++i) {
        if (std::fabs(fc) < std::fabs(fb)) {
            a = b; b = c; c = a;
            fa = fb; fb = fc; fc = fa;
        }
        const double tol1 = 2.0 * eps * std::fabs(b) + 0.5 * tol;
        const double xm = 0.5 * (c - b);
        if (std::fabs(xm) <= tol1 || fb == 0.0) {
            return b;
        }
        if (std::fabs(e) >= tol1 && std::fabs(fa) > std::fabs(fb)) {
            double s = fb / fa;
            double p, q, r;
            if (a == c) {
                // Secant.
                p = 2.0 * xm * s;
                q = 1.0 - s;
            } else {
                // Inverse quadratic interpolation.
                q = fa / fc;
                r = fb / fc;
                p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (p > 0.0) q = -q;
            p = std::fabs(p);
            if (2.0 * p < std::min(3.0 * xm * q - std::fabs(tol1 * q), std::fabs(e * q))) {
                e = d;
                d = p / q;
            } else {
                d = xm;
                e = d;
            }
        } else {
            d = xm;
            e = d;
        }
        a = b;
        fa = fb;
        if (std::fabs(d) > tol1) {
            b += d;
        } else {
            b += (xm > 0.0 ? tol1 : -tol1);
        }
        fb = f(b);
        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
            c = a; fc = fa;
            d = b - a; e = d;
        }
    }
    return b;
}

} // namespace mathx
