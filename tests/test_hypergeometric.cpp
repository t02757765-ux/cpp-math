#include <cassert>
#include <cmath>
#include <iostream>

#include "mathx/hypergeometric.hpp"
#include "mathx/constants.hpp"

int main() {
    // 1F1(1; 1; z) = e^z
    assert(std::fabs(mathx::hypergeometric_1f1(1.0, 1.0, 1.0) - mathx::constants::e) < 1e-9);
    // 2F1(1, 1; 1; z) = 1/(1-z)
    assert(std::fabs(mathx::hypergeometric_2f1(1.0, 1.0, 1.0, 0.5) - 2.0) < 1e-9);
    // 2F1(-1, b; c; z) is a polynomial: 2F1(-1,b;c;z) = 1 - b/c z
    assert(std::fabs(mathx::hypergeometric_2f1(-1.0, 2.0, 3.0, 0.5)
                     - (1.0 - (2.0 / 3.0) * 0.5)) < 1e-9);

    // K(0) = pi/2, E(0) = pi/2
    assert(std::fabs(mathx::elliptic_k(0.0) - mathx::constants::pi / 2.0) < 1e-12);
    assert(std::fabs(mathx::elliptic_e(0.0) - mathx::constants::pi / 2.0) < 1e-12);
    // E(1) = 1
    assert(std::fabs(mathx::elliptic_e(1.0) - 1.0) < 1e-12);
    // K(0.5) ~ 1.854074677301372
    assert(std::fabs(mathx::elliptic_k(0.5) - 1.854074677301372) < 1e-9);
    // E(0.5) ~ 1.350643881047675
    assert(std::fabs(mathx::elliptic_e(0.5) - 1.350643881047675) < 1e-9);

    std::cout << "hypergeometric tests passed\n";
    return 0;
}
