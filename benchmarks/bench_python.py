"""Benchmark cpp-math against NumPy / scipy on equivalent workloads.

Requires the cpp-math Python binding to be importable as `mathx`.

Usage: python benchmarks/bench_python.py [threads]

The thread count is applied to both sides (OpenMP via OMP_NUM_THREADS,
BLAS via threadpoolctl) so the comparison is fair.
"""

import math
import os
import sys
import time

import numpy as np


def _add_mingw_dll_dir():
    # On Windows, CPython 3.8+ needs the MinGW runtime DLLs (libgomp,
    # libstdc++, ...) on the DLL search path for the mathx extension.
    # Only the toolchain used to build mathx (ucrt64) must be added; adding
    # more than one MinGW runtime at once corrupts the heap.
    if os.name == "nt":
        for d in (r"C:\msys64\ucrt64\bin",):
            if os.path.isdir(d):
                try:
                    os.add_dll_directory(d)
                except OSError:
                    pass


_add_mingw_dll_dir()


def set_threads(n):
    os.environ["OMP_NUM_THREADS"] = str(n)
    os.environ["OPENBLAS_NUM_THREADS"] = str(n)
    os.environ["MKL_NUM_THREADS"] = str(n)
    try:
        import threadpoolctl

        threadpoolctl.threadpool_limits(limits=n)
    except ImportError:
        pass


def timeit(fn, reps):
    fn()  # warmup
    t0 = time.perf_counter()
    for _ in range(reps):
        fn()
    return (time.perf_counter() - t0) / reps * 1e3  # ms


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else os.cpu_count() or 1
    set_threads(n)

    try:
        import mathx  # noqa
    except ImportError:
        print("mathx python binding not found; only numpy results shown")
        mathx = None

    print(f"=== threads = {n} ===")
    A = np.random.rand(512, 512)
    B = np.random.rand(512, 512)
    if mathx:
        A_c = mathx.Matrix(A.tolist())
        B_c = mathx.Matrix(B.tolist())
        t_c = timeit(lambda: A_c * B_c, 5)
        print(f"  mathx (C++/OpenMP): {t_c:8.2f} ms")
    t_n = timeit(lambda: A @ B, 5)
    print(f"  numpy (BLAS)       : {t_n:8.2f} ms")

    print("=== fft 65536 pts (ms) ===")
    x = np.sin(0.01 * np.arange(1 << 16))
    if mathx:
        t_c = timeit(lambda: mathx.fft(x.tolist()), 10)
        print(f"  mathx (C++): {t_c:8.3f} ms")
    t_n = timeit(lambda: np.fft.fft(x), 10)
    print(f"  numpy     : {t_n:8.3f} ms")

    print("=== special functions (ms) ===")
    xs = [0.5 + 1e-6 * (i % 1000) for i in range(1_000_000)]
    if mathx:
        t_c = timeit(lambda: [mathx.lambert_w(v) for v in xs[:100000]], 3)
        print(f"  mathx lambert_w : {t_c:8.2f} ms / 100k")
    from scipy.special import lambertw as scipy_lambertw

    t_n = timeit(lambda: [float(scipy_lambertw(v).real) for v in xs[:10000]], 3)
    print(f"  scipy lambert_w : {t_n:8.2f} ms / 10k")

    print("=== solve 200x200 (ms) ===")
    M = np.random.rand(200, 200) + 10 * np.eye(200)
    b = np.random.rand(200)
    if mathx:
        M_c = mathx.Matrix(M.tolist())
        t_c = timeit(lambda: mathx.solve(M_c, b.tolist()), 5)
        print(f"  mathx (C++): {t_c:8.2f} ms")
    t_n = timeit(lambda: np.linalg.solve(M, b), 5)
    print(f"  numpy     : {t_n:8.2f} ms")


if __name__ == "__main__":
    main()
