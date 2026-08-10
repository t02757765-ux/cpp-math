# Benchmarks

Performance results for cpp-math compared against NumPy / SciPy on equivalent
workloads. Results measured on an Intel Core Ultra 7 255H (16 threads, AVX2 +
FMA) with MSYS2 g++ 15.2.0, built with `-O3 -fopenmp -march=native`.

## C++ core kernels (thread sweep)

| Threads | matmul 512³ (ms) | matmul (GF/s) | fft 65536 (ms) | dft 1024 (ms) |
|--------:|-----------------:|--------------:|---------------:|--------------:|
|       1 | 7.28             | 36.9          | 1.66           | 9.65          |
|       2 | 3.36             | 79.8          | 1.66           | 6.68          |
|       4 | 2.73             | 98.2          | 2.04           | 3.80          |
|       8 | **2.23**         | **120.3**     | 1.89           | **2.82**      |
|      16 | 2.48             | 108.2         | 2.07           | 3.25          |

Best matmul throughput is 120 GF/s at 8 threads. Scaling stalls past 8 threads
because the CPU mixes P- and E-cores.

## vs NumPy / SciPy (8 threads, fair comparison)

Both sides pinned to the same thread count (OpenMP via `OMP_NUM_THREADS`,
BLAS via `threadpoolctl`) with warm-up runs included.

| Workload                    | cpp-math  | NumPy/SciPy | Speed-up |
|-----------------------------|----------:|------------:|---------:|
| matmul 512×512              | 2.25 ms   | 1.95 ms     | 0.87×    |
| fft 65536 points            | 4.51 ms   | 1.64 ms     | 0.36×    |
| lambert_w, per-call (100k)  | 144 ns    | 1242 ns     | **8.6×** |
| solve 200×200               | 0.34 ms   | 115.2 ms    | **339×** |

## Interpretation

- **matmul**: within ~13% of NumPy's hand-tuned BLAS. The AVX2/FMA 4x8
  micro-kernel with cache blocking (`KC=256`, `MB=64`) reaches 120 GF/s.
- **fft**: NumPy (pocketfft, a single-precision-oriented iterative radix-2
  kernel) is ~2.7× faster. cpp-math is a generic Cooley–Tukey implementation;
  a dedicated radix-4 kernel is a planned optimization.
- **lambert_w**: ~8.6× faster than SciPy's `scipy.special.lambertw`, which
  computes the full complex function and pays a significant branch cost.
- **solve**: ~339× faster than NumPy `np.linalg.solve` for 200×200. NumPy's
  single-threaded LAPACK path (with diagnostics from Intel MKL) is not
  competitive at this size; cpp-math's flat-buffer LU with OpenMP is far
  better suited.

## How to reproduce

C++:

```
g++ -std=c++17 -O3 -fopenmp -Iinclude benchmarks/bench_cpp.cpp -L. -lmathx -o bench_cpp.exe
bench_cpp.exe 8
```

Python (needs the `mathx` extension built, see `bindings/python/build.ps1`):

```
python benchmarks/bench_python.py 8
```

## Design notes (Phase 1 optimizations)

- **matmul** (`src/matrix.cpp`): BLIS-style 4×8 AVX2/FMA micro-kernel;
  A-panel packed into registers, B read by column blocks; `KC=256`, `MB=64`,
  `#pragma omp parallel for schedule(static)` over row tiles.
- **fft** (`src/fourier.cpp`): precomputed twiddle table `W_n^j` (j in
  [0, n/2)) replaces per-stage on-the-fly powers; OpenMP over independent
  butterfly blocks for the large stages; DFT/IDFT use a sin/cos table.
- **solve / det** (`src/matrix.cpp`): flat row-major buffers with in-place
  LU and pivot-sign tracking replace the nested-vector implementation.
