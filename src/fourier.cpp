#include "mathx/fourier.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>

namespace mathx {

namespace {
constexpr double kPi = 3.14159265358979323846;

template <typename T>
std::vector<std::complex<double>> to_complex(const std::vector<T>& x) {
    return std::vector<std::complex<double>>(x.begin(), x.end());
}

bool is_power_of_two(std::size_t n) {
    return n != 0 && (n & (n - 1)) == 0;
}
} // namespace

std::vector<std::complex<double>> dft(const std::vector<std::complex<double>>& x) {
    const std::size_t n = x.size();
    std::vector<std::complex<double>> X(n);
    for (std::size_t k = 0; k < n; ++k) {
        std::complex<double> sum(0.0, 0.0);
        for (std::size_t t = 0; t < n; ++t) {
            const double angle = -2.0 * kPi * static_cast<double>(k * t) / static_cast<double>(n);
            sum += x[t] * std::polar(1.0, angle);
        }
        X[k] = sum;
    }
    return X;
}

std::vector<std::complex<double>> idft(const std::vector<std::complex<double>>& X) {
    const std::size_t n = X.size();
    std::vector<std::complex<double>> x(n);
    for (std::size_t t = 0; t < n; ++t) {
        std::complex<double> sum(0.0, 0.0);
        for (std::size_t k = 0; k < n; ++k) {
            const double angle = 2.0 * kPi * static_cast<double>(k * t) / static_cast<double>(n);
            sum += X[k] * std::polar(1.0, angle);
        }
        x[t] = sum / static_cast<double>(n);
    }
    return x;
}

std::vector<std::complex<double>> fft(const std::vector<std::complex<double>>& x) {
    const std::size_t n = x.size();
    if (!is_power_of_two(n)) {
        throw std::invalid_argument("mathx::fft: input size must be a power of two");
    }
    if (n == 0) return {};

    // Bit-reversal permutation.
    std::vector<std::complex<double>> a = x;
    for (std::size_t i = 1, j = 0; i < n; ++i) {
        std::size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;
        if (i < j) {
            std::swap(a[i], a[j]);
        }
    }

    for (std::size_t len = 2; len <= n; len <<= 1) {
        const double angle = -2.0 * kPi / static_cast<double>(len);
        const std::complex<double> wlen(std::cos(angle), std::sin(angle));
        for (std::size_t i = 0; i < n; i += len) {
            std::complex<double> w(1.0, 0.0);
            for (std::size_t j = 0; j < len / 2; ++j) {
                std::complex<double> u = a[i + j];
                std::complex<double> v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
    return a;
}

std::vector<std::complex<double>> ifft(const std::vector<std::complex<double>>& X) {
    std::vector<std::complex<double>> x = X;
    for (auto& v : x) {
        v = std::conj(v);
    }
    x = fft(x);
    const std::size_t n = x.size();
    for (auto& v : x) {
        v = std::conj(v) / static_cast<double>(n);
    }
    return x;
}

std::vector<std::complex<double>> fft(const std::vector<double>& x) {
    return fft(to_complex(x));
}

std::vector<std::complex<double>> dft(const std::vector<double>& x) {
    return dft(to_complex(x));
}

} // namespace mathx
