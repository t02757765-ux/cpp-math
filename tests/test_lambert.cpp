#include <cassert>
#include <cmath>
#include <iostream>

#include "mathx/lambert.hpp"

int main() {
    // W(0) = 0
    assert(std::fabs(mathx::lambert_w(0.0) - 0.0) < 1e-12);
    // W(e) = 1
    assert(std::fabs(mathx::lambert_w(std::exp(1.0)) - 1.0) < 1e-12);
    // W(1) = Omega constant ~0.56714329
    const double omega = 0.5671432904097838729999686622;
    assert(std::fabs(mathx::lambert_w(1.0) - omega) < 1e-9);
    // W(-1/e) = -1 (both branches)
    const double inv_e = -1.0 / std::exp(1.0);
    assert(std::fabs(mathx::lambert_w(inv_e, 0) + 1.0) < 1e-9);
    assert(std::fabs(mathx::lambert_w(inv_e, -1) + 1.0) < 1e-9);
    // W_0(x) and W_{-1}(x) at -0.2 differ.
    const double w0 = mathx::lambert_w(-0.2, 0);
    const double wm1 = mathx::lambert_w(-0.2, -1);
    assert(std::fabs(w0 - wm1) > 1e-3);
    // Round-trip: W(x) * exp(W(x)) == x
    const double x = 3.7;
    const double w = mathx::lambert_w(x);
    assert(std::fabs(w * std::exp(w) - x) < 1e-12);
    // Complex principal branch: W(1 + 2i).
    const std::complex<double> z(1.0, 2.0);
    const std::complex<double> wz = mathx::lambert_w(z);
    const std::complex<double> r = wz * std::exp(wz) - z;
    assert(std::abs(r) < 1e-12);

    std::cout << "lambert tests passed\n";
    return 0;
}
