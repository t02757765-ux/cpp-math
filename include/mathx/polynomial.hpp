#ifndef MATHX_POLYNOMIAL_HPP
#define MATHX_POLYNOMIAL_HPP

#include "export.hpp"
#include <complex>
#include <cstddef>
#include <vector>

namespace mathx {

//! Polynomial with double coefficients, stored in ascending order
//! (index 0 = constant term).
class MATHX_API Polynomial {
public:
    Polynomial() = default;
    explicit Polynomial(std::vector<double> coeffs);

    //! Degree of the polynomial (-1 for the zero polynomial).
    int degree() const;
    //! Leading coefficient.
    double leading() const;
    //! Evaluate at x using Horner's scheme.
    double operator()(double x) const;
    //! Evaluate at a complex argument.
    std::complex<double> operator()(const std::complex<double>& x) const;
    //! Coefficient of x^k.
    double coeff(std::size_t k) const;
    //! Trim high-order zero coefficients.
    void trim();

    //! Symbolic derivative.
    Polynomial derivative() const;
    //! Symbolic integral (integration constant 0).
    Polynomial integral() const;

    //! All real roots (for degree <= 3 exact; higher degrees via companion
    //! matrix / Durand-Kerner). May be approximate.
    std::vector<double> real_roots() const;

    //! All complex roots (Durand-Kerner method).
    std::vector<std::complex<double>> roots() const;

    std::vector<double> coeffs_;

    //! Size (number of stored coefficients).
    std::size_t size() const { return coeffs_.size(); }
};

MATHX_API Polynomial operator+(const Polynomial& a, const Polynomial& b);
MATHX_API Polynomial operator-(const Polynomial& a, const Polynomial& b);
MATHX_API Polynomial operator*(const Polynomial& a, const Polynomial& b);
MATHX_API Polynomial operator*(double s, const Polynomial& p);
MATHX_API Polynomial operator/(const Polynomial& a, const Polynomial& b); // quotient
MATHX_API Polynomial operator%(const Polynomial& a, const Polynomial& b); // remainder
MATHX_API std::vector<double> polynomial_divmod(const Polynomial& a,
                                                const Polynomial& b,
                                                Polynomial& quotient,
                                                Polynomial& remainder);

} // namespace mathx

#endif
