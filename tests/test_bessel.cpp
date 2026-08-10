#include <cassert>
#include <cmath>
#include <iostream>

#include "mathx/bessel.hpp"
#include "mathx/constants.hpp"

int main() {
    // J_0 known values
    assert(std::fabs(mathx::bessel_j0(0.0) - 1.0) < 1e-9);
    assert(std::fabs(mathx::bessel_j0(1.0) - 0.7651976865579666) < 1e-6);
    // J_1
    assert(std::fabs(mathx::bessel_j1(0.0)) < 1e-12);
    assert(std::fabs(mathx::bessel_j1(1.0) - 0.4400505857449335) < 1e-6);
    // J_2 = 2/x J_1 - J_0
    assert(std::fabs(mathx::bessel_jn(2, 2.0)
                     - (2.0 / 2.0 * mathx::bessel_j1(2.0) - mathx::bessel_j0(2.0))) < 1e-9);
    // Y_0, Y_1
    assert(std::fabs(mathx::bessel_y0(1.0) - 0.0882569642156770) < 1e-6);
    assert(std::fabs(mathx::bessel_y1(1.0) - (-0.7812128213002887)) < 1e-6);
    // I_0(0)=1, I_1(0)=0
    assert(std::fabs(mathx::bessel_i0(0.0) - 1.0) < 1e-9);
    assert(std::fabs(mathx::bessel_i1(0.0)) < 1e-12);
    assert(std::fabs(mathx::bessel_i0(1.0) - 1.2660658777520084) < 1e-6);
    // K_0(1) ~ 0.421024, K_1(1) ~ 0.601907
    assert(std::fabs(mathx::bessel_k0(1.0) - 0.4210244382407083) < 1e-6);
    assert(std::fabs(mathx::bessel_k1(1.0) - 0.6019072301972346) < 1e-6);
    // Airy Ai(0) = 0.355028053887817, Bi(0) = 0.614926627446001
    assert(std::fabs(mathx::airy_ai(0.0) - 0.3550280538878172) < 1e-9);
    assert(std::fabs(mathx::airy_bi(0.0) - 0.614926627446001) < 1e-9);

    // constants sanity
    assert(std::fabs(mathx::constants::pi - 3.141592653589793) < 1e-15);
    assert(std::fabs(mathx::constants::e - 2.718281828459045) < 1e-15);

    std::cout << "bessel tests passed\n";
    return 0;
}
