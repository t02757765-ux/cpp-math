# cpp-math

**Advanced C++ mathematics library** — a NumPy/SymPy-inspired toolkit of special functions, transforms and linear algebra, compiled as a shared library (`.dll` / `.so`).

## Features

- **Special functions**
  - Lambert W function (real: W₀ and W₋₁ branches, plus complex W₀)
  - Gamma / log-Gamma / digamma
  - Error function `erf`, complementary error function `erfc`
  - Normal CDF
  - Riemann zeta function
  - Beta function
- **Transforms**
  - Discrete Fourier Transform (DFT / IDFT)
  - Fast Fourier Transform (iterative Cooley–Tukey, O(n log n))
  - Real-input overloads
- **Linear algebra**
  - Dense `Matrix` with arithmetic and matrix–vector products
  - Gaussian elimination with partial pivoting (`solve`)
  - Determinant

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

// FFT
auto X = mathx::fft(std::vector<double>{1, 2, 3, 4});

// Solve A x = b
mathx::Matrix A(2, 2, {4, 2, 1, 3});
auto x = mathx::solve(A, {6, 4});            // {1, 1}
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

## Project layout

```
include/mathx/   public headers
src/             implementation (.cpp)
tests/           unit tests (C++ assert based, CTest)
```

## Roadmap

- Bessel functions, elliptic integrals, Airy functions
- Numerical root-finding and integration (quadrature)
- Complex matrices and Eigen-style decomposition
- OpenMP / SIMD acceleration
- Python bindings

## License

MIT
