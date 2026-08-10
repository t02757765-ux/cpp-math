// Benchmark for cpp-math core kernels.
// Compile:  g++ -O3 -fopenmp -Iinclude benchmarks/bench_cpp.cpp -L. -lmathx
#include <cstdio>
#include <chrono>
#include <cmath>
#include <vector>

#include <mathx/mathx.hpp>

namespace {
double now_ms() {
    return std::chrono::duration<double, std::milli>(
        std::chrono::high_resolution_clock::now().time_since_epoch()).count();
}

//! MFLOP/s of a matrix multiply A(MxK) * B(KxN).
void bench_matmul() {
    const std::size_t n = 512;
    mathx::Matrix A(n, n), B(n, n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            A(i, j) = std::sin(static_cast<double>(i + j));
            B(i, j) = std::cos(static_cast<double>(i - j));
        }
    // Warmup
    volatile auto tmp = A * B;
    const int reps = 5;
    const double t0 = now_ms();
    for (int r = 0; r < reps; ++r) {
        volatile auto C = A * B;
    }
    const double dt = now_ms() - t0;
    const double flops = 2.0 * n * n * n * reps;
    printf("matmul %zux%zu        : %8.2f ms   %10.1f MFLOP/s\n",
           n, n, dt / reps, flops / dt / 1e3);
}

void bench_fft() {
    const std::size_t n = 1 << 16; // 65536
    std::vector<double> x(n);
    for (std::size_t i = 0; i < n; ++i) x[i] = std::sin(0.01 * static_cast<double>(i));
    const int reps = 20;
    const double t0 = now_ms();
    for (int r = 0; r < reps; ++r) {
        volatile auto X = mathx::fft(x);
    }
    const double dt = now_ms() - t0;
    printf("fft  %zu pts         : %8.3f ms  (%.1f x 5n log n FLOP)\n",
           n, dt / reps, 5.0 * n * std::log2(static_cast<double>(n)) / (dt / reps) / 1e6);
}

void bench_special() {
    const int reps = 2000000;
    const double t0 = now_ms();
    double acc = 0.0;
    for (int i = 0; i < reps; ++i) {
        const double x = 0.5 + 0.000001 * static_cast<double>(i % 1000);
        acc += mathx::lambert_w(x) + mathx::gamma(x) + mathx::bessel_j0(x);
    }
    const double dt = now_ms() - t0;
    printf("special funcs x3     : %8.3f ms  (%6.2f ns / call)\n",
           dt, dt * 1e6 / (3.0 * reps));
    printf("(sink %.3f)\n", acc);
}

void bench_dft() {
    const std::size_t n = 1024;
    std::vector<double> x(n);
    for (std::size_t i = 0; i < n; ++i) x[i] = std::sin(static_cast<double>(i));
    const int reps = 3;
    const double t0 = now_ms();
    for (int r = 0; r < reps; ++r) {
        volatile auto X = mathx::dft(x);
    }
    const double dt = now_ms() - t0;
    printf("dft   %zu pts (O(n^2)): %8.2f ms\n", n, dt / reps);
}
} // namespace

int main() {
    printf("=== cpp-math benchmark (OpenMP %s) ===\n",
#ifdef _OPENMP
    "on"
#else
    "off"
#endif
    );
    bench_matmul();
    bench_fft();
    bench_special();
    bench_dft();
    return 0;
}
