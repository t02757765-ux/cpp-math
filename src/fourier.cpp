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
    if (n == 0) return X;

    // Precompute the n-th roots of unity to avoid std::polar per element.
    std::vector<double> cs(n);
    const double inv = 2.0 * kPi / static_cast<double>(n);
    for (std::size_t t = 0; t < n; ++t) {
        const double a = -inv * static_cast<double>(t);
        cs[t] = std::sin(a);  // sin(-a) = -sin(a), store sin only
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (n > 256)
#endif
    for (std::size_t k = 0; k < n; ++k) {
        std::complex<double> sum(0.0, 0.0);
        for (std::size_t t = 0; t < n; ++t) {
            // e^{-i k t inv} = cos(k*theta) - i sin(k*theta)
            const double theta = -inv * static_cast<double>(k * t % n);
            sum += x[t] * std::complex<double>(std::cos(theta), std::sin(theta));
        }
        X[k] = sum;
    }
    return X;
}

std::vector<std::complex<double>> idft(const std::vector<std::complex<double>>& X) {
    const std::size_t n = X.size();
    std::vector<std::complex<double>> x(n);
    if (n == 0) return x;

    const double inv = 2.0 * kPi / static_cast<double>(n);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (n > 256)
#endif
    for (std::size_t t = 0; t < n; ++t) {
        std::complex<double> sum(0.0, 0.0);
        for (std::size_t k = 0; k < n; ++k) {
            const double theta = inv * static_cast<double>(k * t % n);
            sum += X[k] * std::complex<double>(std::cos(theta), std::sin(theta));
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
    if (n == 1) return x;

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

    // Precompute twiddle factors W_n^j = exp(-2 pi i j / n), j in [0, n/2).
    // Avoiding per-butterfly cos/sin and the serial w *= wlen chain lets the
    // compiler vectorize the inner loop.
    std::vector<std::complex<double>> tw(n / 2);
    const double base = -2.0 * kPi / static_cast<double>(n);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (n > 65536)
#endif
    for (std::size_t j = 0; j < n / 2; ++j) {
        const double a0 = base * static_cast<double>(j);
        tw[j] = std::complex<double>(std::cos(a0), std::sin(a0));
    }

    for (std::size_t len = 2; len <= n; len <<= 1) {
        const std::size_t half = len / 2;
        const std::size_t step = n / len;  // stride into twiddle table
#ifdef _OPENMP
        // Independent butterfly blocks for large stages.
#pragma omp parallel for schedule(static) if (len >= n / 4)
#endif
        for (std::size_t i = 0; i < n; i += len) {
            for (std::size_t j = 0; j < half; ++j) {
                const std::complex<double>& w = tw[j * step];
                const std::complex<double> u = a[i + j];
                const std::complex<double> v = a[i + j + half] * w;
                a[i + j] = u + v;
                a[i + j + half] = u - v;
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
