#ifndef MATHX_BESSEL_HPP
#define MATHX_BESSEL_HPP

#include "export.hpp"
#include <complex>

namespace mathx {

//! Bessel function of the first kind, order 0: J_0(x).
MATHX_API double bessel_j0(double x);
//! Bessel function of the first kind, order 1: J_1(x).
MATHX_API double bessel_j1(double x);
//! Bessel function of the first kind, integer order n: J_n(x).
MATHX_API double bessel_jn(int n, double x);

//! Bessel function of the second kind, order 0: Y_0(x).
MATHX_API double bessel_y0(double x);
//! Bessel function of the second kind, order 1: Y_1(x).
MATHX_API double bessel_y1(double x);
//! Bessel function of the second kind, integer order n: Y_n(x).
MATHX_API double bessel_yn(int n, double x);

//! Modified Bessel function of the first kind, order 0: I_0(x).
MATHX_API double bessel_i0(double x);
//! Modified Bessel function of the first kind, order 1: I_1(x).
MATHX_API double bessel_i1(double x);
//! Modified Bessel function of the first kind, integer order n: I_n(x).
MATHX_API double bessel_in(int n, double x);

//! Modified Bessel function of the second kind, order 0: K_0(x).
MATHX_API double bessel_k0(double x);
//! Modified Bessel function of the second kind, order 1: K_1(x).
MATHX_API double bessel_k1(double x);
//! Modified Bessel function of the second kind, integer order n: K_n(x).
MATHX_API double bessel_kn(int n, double x);

//! Airy function Ai(x).
MATHX_API double airy_ai(double x);
//! Airy function Bi(x).
MATHX_API double airy_bi(double x);

} // namespace mathx

#endif
