#include "mathx/ode.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace mathx {

std::vector<std::pair<double, double>> ode_euler(
    const OdeFunction& f, double y0, double t0, double t1, int n) {
    if (n <= 0) return {};
    const double h = (t1 - t0) / static_cast<double>(n);
    std::vector<std::pair<double, double>> out;
    out.reserve(n + 1);
    double t = t0;
    double y = y0;
    out.emplace_back(t, y);
    for (int i = 0; i < n; ++i) {
        y += h * f(t, y);
        t += h;
        out.emplace_back(t, y);
    }
    return out;
}

std::vector<std::pair<double, double>> ode_heun(
    const OdeFunction& f, double y0, double t0, double t1, int n) {
    if (n <= 0) return {};
    const double h = (t1 - t0) / static_cast<double>(n);
    std::vector<std::pair<double, double>> out;
    out.reserve(n + 1);
    double t = t0;
    double y = y0;
    out.emplace_back(t, y);
    for (int i = 0; i < n; ++i) {
        const double k1 = f(t, y);
        const double k2 = f(t + h, y + h * k1);
        y += 0.5 * h * (k1 + k2);
        t += h;
        out.emplace_back(t, y);
    }
    return out;
}

std::vector<std::pair<double, double>> ode_rk4(
    const OdeFunction& f, double y0, double t0, double t1, int n) {
    if (n <= 0) return {};
    const double h = (t1 - t0) / static_cast<double>(n);
    std::vector<std::pair<double, double>> out;
    out.reserve(n + 1);
    double t = t0;
    double y = y0;
    out.emplace_back(t, y);
    for (int i = 0; i < n; ++i) {
        const double k1 = f(t, y);
        const double k2 = f(t + 0.5 * h, y + 0.5 * h * k1);
        const double k3 = f(t + 0.5 * h, y + 0.5 * h * k2);
        const double k4 = f(t + h, y + h * k3);
        y += h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        t += h;
        out.emplace_back(t, y);
    }
    return out;
}

std::vector<std::pair<double, double>> ode_rk45(
    const OdeFunction& f, double y0, double t0, double t1, double tol) {
    std::vector<std::pair<double, double>> out;
    out.reserve(1024);
    double t = t0;
    double y = y0;
    out.emplace_back(t, y);
    double h = (t1 - t0) / 16.0;
    if (h == 0.0) return out;

    while (t < t1) {
        if (t + h > t1) h = t1 - t;
        // Dormand-Prince coefficients (RK5(4)7FM).
        const double k1 = f(t, y);
        const double k2 = f(t + 0.2 * h, y + h * (0.2 * k1));
        const double k3 = f(t + 0.3 * h, y + h * (0.075 * k1 + 0.225 * k2));
        const double k4 = f(t + 0.8 * h, y + h * (0.9777777777777777 * k1
                            - 3.7333333333333334 * k2 + 3.5555555555555554 * k3));
        const double k5 = f(t + 8.0 / 9.0 * h, y + h * (2.9525986892242035 * k1
                            - 11.595793324188385 * k2 + 9.822892851699975 * k3
                            - 0.29080932784636487 * k4));
        const double k6 = f(t + h, y + h * (2.8462752525252526 * k1
                            - 10.757575757575758 * k2 + 8.906422717743473 * k3
                            + 0.2784090909090909 * k4 - 0.2735313036020583 * k5));
        const double k7 = f(t + h, y + h * (0.09114583333333333 * k1
                            + 0.4492362982920919 * k3 + 0.6510416666666666 * k4
                            - 0.3223761792452830 * k5 + 0.13095238095238096 * k6));

        const double y5 = y + h * (35.0 / 384.0 * k1
                             + 500.0 / 1113.0 * k3 + 125.0 / 192.0 * k4
                             - 2187.0 / 6784.0 * k5 + 11.0 / 84.0 * k6);
        const double y4 = y + h * (5179.0 / 57600.0 * k1
                             + 7571.0 / 16695.0 * k3 + 393.0 / 640.0 * k4
                             - 92097.0 / 339200.0 * k5 + 187.0 / 2100.0 * k6
                             + 1.0 / 40.0 * k7);

        const double err = std::fabs(y5 - y4);
        if (err < tol) {
            t += h;
            y = y5;
            out.emplace_back(t, y);
        }
        // Step-size control: safe factor 0.9, power 1/5.
        double factor = 0.9 * std::pow(tol / (err + 1e-300), 0.2);
        factor = std::min(5.0, std::max(0.2, factor));
        h *= factor;
    }
    return out;
}

} // namespace mathx
