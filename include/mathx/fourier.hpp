#ifndef MATHX_FOURIER_HPP
#define MATHX_FOURIER_HPP

#include "export.hpp"
#include <complex>
#include <cstddef>
#include <vector>

namespace mathx {

//! Discrete Fourier Transform (naive O(n^2)).
//! Forward transform. Output has the same size as the input.
MATHX_API std::vector<std::complex<double>> dft(
    const std::vector<std::complex<double>>& x);

//! Inverse Discrete Fourier Transform.
MATHX_API std::vector<std::complex<double>> idft(
    const std::vector<std::complex<double>>& X);

//! Fast Fourier Transform (iterative Cooley-Tukey, O(n log n)).
//! Input size must be a power of two, otherwise an exception is thrown.
MATHX_API std::vector<std::complex<double>> fft(
    const std::vector<std::complex<double>>& x);

//! Inverse FFT.
MATHX_API std::vector<std::complex<double>> ifft(
    const std::vector<std::complex<double>>& X);

//! Real-input convenience overloads.
MATHX_API std::vector<std::complex<double>> fft(
    const std::vector<double>& x);

MATHX_API std::vector<std::complex<double>> dft(
    const std::vector<double>& x);

} // namespace mathx

#endif
