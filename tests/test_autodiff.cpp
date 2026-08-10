#include <cassert>
#include <cmath>
#include <iostream>

#include "mathx/autodiff.hpp"

int main() {
    using mathx::Dual;

    // f(x) = x^2, f'(2) = 4
    auto f1 = [](const Dual& x) { return x * x; };
    auto d1 = mathx::derivative(f1, 2.0);
    assert(std::fabs(d1.second - 4.0) < 1e-12);

    // f(x) = sin(x), f'(0) = 1
    auto f2 = [](const Dual& x) { return mathx::sin(x); };
    auto d2 = mathx::derivative(f2, 0.0);
    assert(std::fabs(d2.second - 1.0) < 1e-12);

    // f(x) = exp(x^2), f'(1) = 2e
    auto f3 = [](const Dual& x) { return mathx::exp(x * x); };
    auto d3 = mathx::derivative(f3, 1.0);
    assert(std::fabs(d3.second - 2.0 * std::exp(1.0)) < 1e-9);

    // f(x) = ln(x) * sqrt(x), f'(2)
    auto f4 = [](const Dual& x) { return mathx::log(x) * mathx::sqrt(x); };
    auto d4 = mathx::derivative(f4, 2.0);
    // analytic: (ln x)(1/(2 sqrt x)) + sqrt(x)/x  at x=2
    const double analytic = std::log(2.0) / (2.0 * std::sqrt(2.0)) + std::sqrt(2.0) / 2.0;
    assert(std::fabs(d4.second - analytic) < 1e-9);

    // f(x) = tanh(x), f'(0) = 1
    auto f5 = [](const Dual& x) { return mathx::tanh(x); };
    auto d5 = mathx::derivative(f5, 0.0);
    assert(std::fabs(d5.second - 1.0) < 1e-12);

    // f(x) = pow(x, 3) via Dual pow
    auto f6 = [](const Dual& x) { return mathx::pow(x, 3.0); };
    auto d6 = mathx::derivative(f6, 2.0);
    assert(std::fabs(d6.second - 12.0) < 1e-9);

    // value must be preserved too
    assert(std::fabs(d1.first - 4.0) < 1e-12);
    assert(std::fabs(d2.first - 0.0) < 1e-12);

    // mixed arithmetic: f(x) = 3x^2 + 2x + 1, f'(5) = 32
    auto f7 = [](const Dual& x) { return 3.0 * x * x + 2.0 * x + 1.0; };
    auto d7 = mathx::derivative(f7, 5.0);
    assert(std::fabs(d7.second - 32.0) < 1e-12);

    std::cout << "autodiff tests passed\n";
    return 0;
}
