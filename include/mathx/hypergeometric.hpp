#ifndef MATHX_HYPERGEOMETRIC_HPP
#define MATHX_HYPERGEOMETRIC_HPP

#include "export.hpp"

namespace mathx {

//! Generalized hypergeometric function pFq evaluated by series summation.
//! a, b are the parameter lists. Converges for |z| < 1 (with 2F1), and
//! for p <= q anywhere.  Returns NaN if the series does not converge.
MATHX_API double hypergeometric_pfq(const double* a, int pa,
                                    const double* b, int pb,
                                    double z,
                                    int max_terms = 10000);

//! Kummer's confluent hypergeometric function 1F1(a; b; z).
MATHX_API double hypergeometric_1f1(double a, double b, double z);

//! Gauss hypergeometric function 2F1(a, b; c; z).
MATHX_API double hypergeometric_2f1(double a, double b, double c, double z);

//! Complete elliptic integral of the first kind K(m).
//! Argument is the parameter m (0 <= m <= 1).
MATHX_API double elliptic_k(double m);

//! Complete elliptic integral of the second kind E(m).
//! Argument is the parameter m (0 <= m <= 1).
MATHX_API double elliptic_e(double m);

} // namespace mathx

#endif
