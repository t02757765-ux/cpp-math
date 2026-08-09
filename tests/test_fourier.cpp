#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

#include "mathx/fourier.hpp"

int main() {
    // FFT must equal DFT for a power-of-two signal.
    const std::size_t n = 16;
    std::vector<std::complex<double>> x(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double t = static_cast<double>(i);
        x[i] = {std::sin(2.0 * t) + 0.5 * std::cos(t), t * 0.1};
    }

    const auto X_fft = mathx::fft(x);
    const auto X_dft = mathx::dft(x);
    for (std::size_t k = 0; k < n; ++k) {
        assert(std::abs(X_fft[k] - X_dft[k]) < 1e-9);
    }

    // Round-trip: ifft(fft(x)) == x.
    const auto y = mathx::ifft(X_fft);
    for (std::size_t k = 0; k < n; ++k) {
        assert(std::abs(y[k] - x[k]) < 1e-9);
    }

    // Real-input overload.
    const auto z = mathx::fft(std::vector<double>{1.0, 0.0, 0.0, 0.0});
    assert(std::fabs(z[0].real() - 1.0) < 1e-12);

    // Non-power-of-two must throw.
    bool threw = false;
    try {
        mathx::fft(std::vector<std::complex<double>>(3));
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    assert(threw);

    std::cout << "fourier tests passed\n";
    return 0;
}
