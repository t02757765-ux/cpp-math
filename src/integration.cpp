#include "mathx/integration.hpp"

#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <vector>

namespace mathx {

namespace {

constexpr double kMachineEps = 2.220446049250313e-16;

double adaptive_simpson_rec(const std::function<double(double)>& f,
                            double a, double b, double fa, double fm, double fb,
                            double s, double whole, double tol, int depth,
                            int max_depth, double eps) {
    const double m = (a + b) * 0.5;
    const double lm = (a + m) * 0.5;
    const double rm = (m + b) * 0.5;
    const double flm = f(lm);
    const double frm = f(rm);
    const double sl = (b - a) * 0.25 / 3.0 * (fa + 4.0 * flm + fm);
    const double sr = (b - a) * 0.25 / 3.0 * (fm + 4.0 * frm + fb);
    const double s2 = sl + sr;
    const double delta = s2 - whole;
    if (depth >= max_depth || std::fabs(delta) <= 15.0 * tol) {
        return s2 + delta / 15.0;
    }
    return adaptive_simpson_rec(f, a, m, fa, flm, fm, s, sl, tol * 0.5, depth + 1, max_depth, eps)
         + adaptive_simpson_rec(f, m, b, fm, frm, fb, s, sr, tol * 0.5, depth + 1, max_depth, eps);
}

} // namespace

double integrate_simpson(const std::function<double(double)>& f,
                         double a, double b, double tol, int max_depth) {
    const double eps = std::max(tol, 10.0 * kMachineEps);
    const double m = (a + b) * 0.5;
    const double fa = f(a);
    const double fm = f(m);
    const double fb = f(b);
    const double whole = (b - a) / 6.0 * (fa + 4.0 * fm + fb);
    return adaptive_simpson_rec(f, a, b, fa, fm, fb, whole, whole, eps, 0, max_depth, eps);
}

double integrate_trapezoid(const std::function<double(double)>& f,
                           double a, double b, int n) {
    if (n <= 0) return 0.0;
    const double h = (b - a) / static_cast<double>(n);
    double sum = 0.5 * (f(a) + f(b));
    for (int i = 1; i < n; ++i) {
        sum += f(a + h * static_cast<double>(i));
    }
    return sum * h;
}

double integrate_simpson_n(const std::function<double(double)>& f,
                           double a, double b, int n) {
    if (n % 2 != 0) ++n;
    if (n < 2) return 0.0;
    const double h = (b - a) / static_cast<double>(n);
    double sum = f(a) + f(b);
    for (int i = 1; i < n; ++i) {
        sum += (i % 2 == 0 ? 2.0 : 4.0) * f(a + h * static_cast<double>(i));
    }
    return sum * h / 3.0;
}

namespace {
// Gauss-Legendre nodes and weights (n = 16, 32). n=5 and n=10 use n=16 subset? No -
// we provide 16 and 32 point rules.
struct LegendreTable {
    int n;
    const double* x;
    const double* w;
};

// 16-point Gauss-Legendre.
constexpr double kX16[16] = {
    -0.9894009349916499, -0.9445750230732326, -0.8656312023878318,
    -0.7554044083550030, -0.6178762444026438, -0.4580167776572274,
    -0.2816035507792589, -0.0950125098376374,  0.0950125098376374,
     0.2816035507792589,  0.4580167776572274,  0.6178762444026438,
     0.7554044083550030,  0.8656312023878318,  0.9445750230732326,
     0.9894009349916499};
constexpr double kW16[16] = {
    0.0271524594117541, 0.0622535239386479, 0.0951585116824928,
    0.1246289712555339, 0.1495959888165767, 0.1691565193950025,
    0.1826034150449236, 0.1894506104550685, 0.1894506104550685,
    0.1826034150449236, 0.1691565193950025, 0.1495959888165767,
    0.1246289712555339, 0.0951585116824928, 0.0622535239386479,
    0.0271524594117541};

// 32-point Gauss-Legendre.
constexpr double kX32[32] = {
    -0.9972638618494816, -0.9856115115452684, -0.9647622555875064,
    -0.9349060759377397, -0.8963211557660521, -0.8493676137325700,
    -0.7944837959679424, -0.7321821187402897, -0.6630442669302152,
    -0.5877157572407623, -0.5068999089322294, -0.4213512761306354,
    -0.3318686022821277, -0.2392873622521371, -0.1444719615827965,
    -0.0483076656877383,  0.0483076656877383,  0.1444719615827965,
     0.2392873622521371,  0.3318686022821277,  0.4213512761306354,
     0.5068999089322294,  0.5877157572407623,  0.6630442669302152,
     0.7321821187402897,  0.7944837959679424,  0.8493676137325700,
     0.8963211557660521,  0.9349060759377397,  0.9647622555875064,
     0.9856115115452684,  0.9972638618494816};
constexpr double kW32[32] = {
    0.0070186100094701, 0.0162743947309057, 0.0253920653092621,
    0.0342738629130214, 0.0428358980222267, 0.0509980592623762,
    0.0586840934785356, 0.0658222227763618, 0.0723457941088485,
    0.0781938957870703, 0.0833119242269467, 0.0876520930044038,
    0.0911738786957639, 0.0938443990808046, 0.0956387200792749,
    0.0965400885147278, 0.0965400885147278, 0.0956387200792749,
    0.0938443990808046, 0.0911738786957639, 0.0876520930044038,
    0.0833119242269467, 0.0781938957870703, 0.0723457941088485,
    0.0658222227763618, 0.0586840934785356, 0.0509980592623762,
    0.0428358980222267, 0.0342738629130214, 0.0253920653092621,
    0.0162743947309057, 0.0070186100094701};
} // namespace

double integrate_gauss_legendre(const std::function<double(double)>& f,
                                double a, double b, int n) {
    const double* x;
    const double* w;
    int m;
    switch (n) {
        case 5:
        case 10:
            x = kX16; w = kW16; m = 16; break;
        case 16:
            x = kX16; w = kW16; m = 16; break;
        case 32:
            x = kX32; w = kW32; m = 32; break;
        default:
            x = kX16; w = kW16; m = 16; break;
    }
    const double c = 0.5 * (b - a);
    const double d = 0.5 * (b + a);
    double sum = 0.0;
    for (int i = 0; i < m; ++i) {
        sum += w[i] * f(c * x[i] + d);
    }
    return c * sum;
}

double integrate_romberg(const std::function<double(double)>& f,
                         double a, double b, double tol) {
    constexpr int kmax = 15;
    double t[kmax + 1][kmax + 1] = {0.0};
    double h = b - a;
    t[0][0] = 0.5 * h * (f(a) + f(b));
    for (int k = 1; k <= kmax; ++k) {
        h *= 0.5;
        double sum = 0.0;
        for (int i = 1; i <= (1 << (k - 1)); ++i) {
            sum += f(a + (2 * i - 1) * h);
        }
        t[0][k] = 0.5 * t[0][k - 1] + h * sum;
        for (int j = 1; j <= k; ++j) {
            const double factor = std::pow(4.0, j);
            t[j][k] = t[j - 1][k] + (t[j - 1][k] - t[j - 1][k - 1]) / (factor - 1.0);
        }
        if (k >= 2) {
            const double err = std::fabs(t[k][k] - t[k - 1][k - 1]);
            if (err <= tol) return t[k][k];
        }
    }
    return t[kmax][kmax];
}

} // namespace mathx
