# cpp-math

**Advanced C++ mathematics library** — a NumPy/SymPy-inspired toolkit of special functions, transforms, numerical analysis and linear algebra, compiled as a shared library (`.dll` / `.so`).

📚 **Kapsamlı dökümantasyon:** [`docs/`](docs/index.md) — başlangıç, API referansı, modüller, entegrasyon ve sürüm geçmişi.

## Features

- **Special functions**
  - Lambert W function (real: W₀ and W₋₁ branches, plus complex W₀)
  - Gamma / log-Gamma / digamma
  - Error function `erf`, complementary error function `erfc`
  - Normal CDF, Riemann zeta, Beta function
- **Bessel & Airy**
  - Bessel functions of the first/second kind (J, Y) and modified (I, K)
  - Integer orders n via recurrence relations
  - Airy functions Ai(x) and Bi(x)
- **Transforms**
  - Discrete Fourier Transform (DFT / IDFT)
  - Fast Fourier Transform (iterative Cooley–Tukey, O(n log n))
  - Real-input overloads
- **Numerical analysis**
  - Quadrature: adaptive Simpson, trapezoid, Gauss-Legendre, Romberg
  - Root finding: bisection, regula falsi, secant, Newton, Brent
  - ODE solvers: Euler, Heun, RK4, adaptive Dormand-Prince (RK45)
- **Linear algebra**
  - Dense `Matrix` with arithmetic and matrix–vector products
  - Gaussian elimination with partial pivoting (`solve`), determinant
- **Polynomials**
  - `Polynomial` class: evaluation, derivative, integral, roots (exact ≤ 3, Durand-Kerner ≥ 4), division
- **Statistics**
  - Mean, variance, stddev, median, percentiles, covariance, correlation, skewness, kurtosis
- **Constants**
  - π, e, golden ratio, Euler-Mascheroni, Catalan, Apery, Feigenbaum and more

## Building

Requires a C++17 compiler and CMake ≥ 3.16.

```sh
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .
```

This produces `libmathx.so` on Linux/macOS and `mathx.dll` (plus `.lib`) on Windows.

### Enable tests

```sh
cmake .. -DCMAKE_BUILD_TYPE=Release -DMATHX_BUILD_TESTS=ON
cmake --build .
ctest --output-on-failure
```

## Usage

```cpp
#include <mathx/mathx.hpp>

// Lambert W
double w = mathx::lambert_w(1.0);            // 0.56714...

// Special functions
double g = mathx::gamma(5.0);                // 24
double z = mathx::zeta(2.0);                 // pi^2/6
double j = mathx::bessel_j0(1.0);            // 0.76519...

// FFT
auto X = mathx::fft(std::vector<double>{1, 2, 3, 4});

// Numerical integration: int_0^1 x^2 dx
double area = mathx::integrate_gauss_legendre(
    [](double x) { return x * x; }, 0.0, 1.0);

// Root finding: x^2 - 2 = 0
double r = mathx::brent([](double x) { return x*x - 2.0; }, 1.0, 2.0);

// Polynomial roots: x^2 - 3x + 2
mathx::Polynomial p({2.0, -3.0, 1.0});
auto roots = p.real_roots();                 // {1, 2}

// Solve A x = b
mathx::Matrix A(2, 2, {4, 2, 1, 3});
auto x = mathx::solve(A, {6, 4});            // {1, 1}

// ODE: dy/dt = y, y(0)=1 -> y(t) = e^t
auto sol = mathx::ode_rk45([](double, double y) { return y; }, 1.0, 0.0, 1.0);

// Statistics
std::vector<double> data{1, 2, 3, 4, 5};
double m = mathx::mean(data);                // 3
```

## Installation / consumption from another project

```sh
cmake --install . --prefix /path/to/prefix
```

Then in your `CMakeLists.txt`:

```cmake
find_package(mathx CONFIG REQUIRED)
target_link_libraries(your_target PRIVATE mathx::mathx)
```

See [docs/integration.md](docs/integration.md) for all three integration methods.

## Project layout

```
include/mathx/   public headers
src/             implementation (.cpp)
tests/           unit tests (CTest)
docs/            documentation
```

## Documentation

| Document | Description |
|---|---|
| [docs/index.md](docs/index.md) | Overview and module map |
| [docs/getting-started.md](docs/getting-started.md) | Install, build, first program |
| [docs/api-reference.md](docs/api-reference.md) | Complete API reference |
| [docs/modules.md](docs/modules.md) | Mathematics and algorithms behind each module |
| [docs/integration.md](docs/integration.md) | Using the library in your own projects |
| [docs/changelog.md](docs/changelog.md) | Version history |

## Roadmap

- Bessel functions, elliptic integrals, Airy functions (WIP)
- Numerical root-finding and integration (WIP)
- Complex matrices and Eigen-style decomposition
- OpenMP / SIMD acceleration
- Python bindings (pybind11)

## License

MIT
