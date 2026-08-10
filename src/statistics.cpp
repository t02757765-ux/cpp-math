#include "mathx/statistics.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace mathx {

double mean(const std::vector<double>& x) {
    if (x.empty()) return 0.0;
    return std::accumulate(x.begin(), x.end(), 0.0) / static_cast<double>(x.size());
}

double variance(const std::vector<double>& x) {
    const std::size_t n = x.size();
    if (n < 2) return 0.0;
    const double m = mean(x);
    double s = 0.0;
    for (double v : x) {
        const double d = v - m;
        s += d * d;
    }
    return s / static_cast<double>(n - 1);
}

double variance_pop(const std::vector<double>& x) {
    if (x.empty()) return 0.0;
    const double m = mean(x);
    double s = 0.0;
    for (double v : x) {
        const double d = v - m;
        s += d * d;
    }
    return s / static_cast<double>(x.size());
}

double stddev(const std::vector<double>& x) {
    return std::sqrt(variance(x));
}

double stddev_pop(const std::vector<double>& x) {
    return std::sqrt(variance_pop(x));
}

double median(std::vector<double> x) {
    if (x.empty()) return 0.0;
    std::sort(x.begin(), x.end());
    const std::size_t n = x.size();
    if (n % 2 == 1) return x[n / 2];
    return 0.5 * (x[n / 2 - 1] + x[n / 2]);
}

double min(const std::vector<double>& x) {
    if (x.empty()) return 0.0;
    return *std::min_element(x.begin(), x.end());
}

double max(const std::vector<double>& x) {
    if (x.empty()) return 0.0;
    return *std::max_element(x.begin(), x.end());
}

double sum(const std::vector<double>& x) {
    return std::accumulate(x.begin(), x.end(), 0.0);
}

double product(const std::vector<double>& x) {
    return std::accumulate(x.begin(), x.end(), 1.0,
                           [](double a, double b) { return a * b; });
}

double geometric_mean(const std::vector<double>& x) {
    if (x.empty()) return 0.0;
    double s = 0.0;
    for (double v : x) {
        if (v <= 0.0) return 0.0;
        s += std::log(v);
    }
    return std::exp(s / static_cast<double>(x.size()));
}

double harmonic_mean(const std::vector<double>& x) {
    if (x.empty()) return 0.0;
    double s = 0.0;
    for (double v : x) {
        if (v == 0.0) return 0.0;
        s += 1.0 / v;
    }
    return static_cast<double>(x.size()) / s;
}

double percentile(std::vector<double> x, double p) {
    if (x.empty()) return 0.0;
    std::sort(x.begin(), x.end());
    const std::size_t n = x.size();
    if (p <= 0.0) return x[0];
    if (p >= 100.0) return x[n - 1];
    const double rank = (p / 100.0) * static_cast<double>(n - 1);
    const std::size_t lo = static_cast<std::size_t>(std::floor(rank));
    const std::size_t hi = static_cast<std::size_t>(std::ceil(rank));
    if (lo == hi) return x[lo];
    const double frac = rank - static_cast<double>(lo);
    return x[lo] + frac * (x[hi] - x[lo]);
}

double covariance(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size() || x.size() < 2) return 0.0;
    const double mx = mean(x);
    const double my = mean(y);
    double s = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        s += (x[i] - mx) * (y[i] - my);
    }
    return s / static_cast<double>(x.size() - 1);
}

double correlation(const std::vector<double>& x, const std::vector<double>& y) {
    const double sx = stddev(x);
    const double sy = stddev(y);
    if (sx == 0.0 || sy == 0.0) return 0.0;
    return covariance(x, y) / (sx * sy);
}

double skewness(const std::vector<double>& x) {
    const std::size_t n = x.size();
    if (n < 3) return 0.0;
    const double m = mean(x);
    const double s = stddev_pop(x);
    if (s == 0.0) return 0.0;
    double m3 = 0.0;
    for (double v : x) {
        const double d = (v - m) / s;
        m3 += d * d * d;
    }
    return m3 / static_cast<double>(n);
}

double kurtosis(const std::vector<double>& x) {
    const std::size_t n = x.size();
    if (n < 4) return 0.0;
    const double m = mean(x);
    const double s = stddev_pop(x);
    if (s == 0.0) return 0.0;
    double m4 = 0.0;
    for (double v : x) {
        const double d = (v - m) / s;
        m4 += d * d * d * d;
    }
    return m4 / static_cast<double>(n) - 3.0;
}

double range(const std::vector<double>& x) {
    if (x.empty()) return 0.0;
    return max(x) - min(x);
}

} // namespace mathx
