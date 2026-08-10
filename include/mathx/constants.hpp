#ifndef MATHX_CONSTANTS_HPP
#define MATHX_CONSTANTS_HPP

#include <cstdint>

namespace mathx {

//! Fundamental mathematical constants with full double precision.
namespace constants {

//! Pi = 3.14159...
constexpr double pi = 3.141592653589793238462643383279502884;
//! Tau = 2*pi
constexpr double tau = 6.283185307179586476925286766559005768;
//! Euler's number e
constexpr double e = 2.718281828459045235360287471352662498;
//! sqrt(2)
constexpr double sqrt2 = 1.414213562373095048801688724209698079;
//! sqrt(3)
constexpr double sqrt3 = 1.732050807568877293527446341505872367;
//! sqrt(pi)
constexpr double sqrt_pi = 1.772453850905516027298167483341145183;
//! Golden ratio phi
constexpr double phi = 1.618033988749894848204586834365638118;
//! ln(2)
constexpr double ln2 = 0.693147180559945309417232121458176568;
//! ln(10)
constexpr double ln10 = 2.302585092994045684017991454684364208;
//! Euler-Mascheroni constant gamma
constexpr double euler_gamma = 0.577215664901532860606512090082402431;
//! Catalan's constant G
constexpr double catalan = 0.915965594177219015054603514932384111;
//! Apery's constant zeta(3)
constexpr double apery = 1.202056903159594285399738161511449991;
//! Feigenbaum's delta
constexpr double feigenbaum_delta = 4.669201609102990671853203820466201617;
//! Feigenbaum's alpha
constexpr double feigenbaum_alpha = 2.502907875095892822283902873218215786;
//! Khinchin's constant
constexpr double khinchin = 2.685452001065306445309714835481795694;
//! Glaisher-Kinkelin constant
constexpr double glaisher = 1.282427129100622636875342568869791728;
//! Gauss's constant
constexpr double gauss = 0.834626841674073186281429732799046809;
//! Conway's constant (look-and-say)
constexpr double conway = 1.303577269034296391257099112152551891;

} // namespace constants

} // namespace mathx

#endif
