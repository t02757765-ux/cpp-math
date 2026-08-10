#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "mathx/statistics.hpp"

int main() {
    const std::vector<double> x{1.0, 2.0, 3.0, 4.0, 5.0};
    assert(std::fabs(mathx::mean(x) - 3.0) < 1e-12);
    assert(std::fabs(mathx::variance(x) - 2.5) < 1e-12);
    assert(std::fabs(mathx::variance_pop(x) - 2.0) < 1e-12);
    assert(std::fabs(mathx::stddev(x) - std::sqrt(2.5)) < 1e-12);
    assert(std::fabs(mathx::median(x) - 3.0) < 1e-12);
    assert(std::fabs(mathx::min(x) - 1.0) < 1e-12);
    assert(std::fabs(mathx::max(x) - 5.0) < 1e-12);
    assert(std::fabs(mathx::sum(x) - 15.0) < 1e-12);
    assert(std::fabs(mathx::product(x) - 120.0) < 1e-12);
    assert(std::fabs(mathx::geometric_mean(x) - std::pow(120.0, 0.2)) < 1e-12);
    assert(std::fabs(mathx::harmonic_mean(x) - (5.0 / (1.0 + 0.5 + 1.0/3.0 + 0.25 + 0.2))) < 1e-12);
    assert(std::fabs(mathx::range(x) - 4.0) < 1e-12);
    assert(std::fabs(mathx::percentile(x, 50.0) - 3.0) < 1e-12);

    // Perfect linear correlation
    const std::vector<double> y{2.0, 4.0, 6.0, 8.0, 10.0};
    assert(std::fabs(mathx::covariance(x, y) - 5.0) < 1e-12);
    assert(std::fabs(mathx::correlation(x, y) - 1.0) < 1e-12);

    // Skewness of symmetric data ~ 0
    const std::vector<double> sym{-1.0, -1.0, 0.0, 1.0, 1.0};
    assert(std::fabs(mathx::skewness(sym)) < 1e-9);

    // Median of even-length array
    const std::vector<double> even{1.0, 2.0, 3.0, 4.0};
    assert(std::fabs(mathx::median(even) - 2.5) < 1e-12);

    std::cout << "statistics tests passed\n";
    return 0;
}
