#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

#include "mathx/polynomial.hpp"

int main() {
    // x^2 - 3x + 2 = (x-1)(x-2)
    mathx::Polynomial p({2.0, -3.0, 1.0});
    assert(p.degree() == 2);
    assert(std::fabs(p(0.0) - 2.0) < 1e-12);
    assert(std::fabs(p(1.0)) < 1e-12);
    assert(std::fabs(p(2.0)) < 1e-12);

    const auto roots = p.real_roots();
    assert(roots.size() == 2);
    assert((std::fabs(roots[0] - 1.0) < 1e-9 && std::fabs(roots[1] - 2.0) < 1e-9) ||
           (std::fabs(roots[0] - 2.0) < 1e-9 && std::fabs(roots[1] - 1.0) < 1e-9));

    // derivative: 2x - 3
    const auto dp = p.derivative();
    assert(std::fabs(dp(0.0) + 3.0) < 1e-12);
    assert(std::fabs(dp(1.0) + 1.0) < 1e-12);

    // integral of x^2 - 3x + 2 => x^3/3 - 3x^2/2 + 2x, at x=1 => 1/3 - 3/2 + 2 = 5/6
    const auto ip = p.integral();
    assert(std::fabs(ip(1.0) - (5.0 / 6.0)) < 1e-12);

    // addition / multiplication
    mathx::Polynomial q({1.0, 1.0});   // x + 1
    const auto pq = p * q;              // (x^2 - 3x + 2)(x + 1) = x^3 - 2x^2 - x + 2
    assert(pq.degree() == 3);
    assert(std::fabs(pq(2.0)) < 1e-12);
    assert(std::fabs(pq(0.0) - 2.0) < 1e-12);

    // division: p / (x-1) = x - 2
    mathx::Polynomial d({-1.0, 1.0});
    const auto quo = p / d;
    assert(std::fabs(quo(3.0) - 1.0) < 1e-12);

    // cubic with three real roots: (x-1)(x-2)(x-3) = x^3 - 6x^2 + 11x - 6
    mathx::Polynomial c({-6.0, 11.0, -6.0, 1.0});
    const auto croots = c.real_roots();
    assert(croots.size() == 3);
    for (double r : croots) {
        assert(std::fabs(c(r)) < 1e-9);
    }

    // complex roots of x^2 + 1
    mathx::Polynomial c2({1.0, 0.0, 1.0});
    const auto c2roots = c2.roots();
    assert(c2roots.size() == 2);
    for (const auto& z : c2roots) {
        assert(std::fabs(std::abs(c2(z))) < 1e-9);
    }

    std::cout << "polynomial tests passed\n";
    return 0;
}
