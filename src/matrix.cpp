#include "mathx/matrix.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>

#ifdef _OPENMP
#include <omp.h>
#endif

#if defined(__AVX2__)
#include <immintrin.h>
#endif

namespace mathx {

Matrix::Matrix(std::size_t rows, std::size_t cols, double value)
    : rows_(rows), cols_(cols), data_(rows * cols, value) {}

Matrix::Matrix(std::size_t rows, std::size_t cols,
               std::initializer_list<double> values)
    : rows_(rows), cols_(cols), data_(values) {
    if (data_.size() != rows * cols) {
        throw std::invalid_argument("mathx::Matrix: initializer size mismatch");
    }
}

double& Matrix::operator()(std::size_t i, std::size_t j) {
    return data_[i * cols_ + j];
}

double Matrix::operator()(std::size_t i, std::size_t j) const {
    return data_[i * cols_ + j];
}

Matrix Matrix::transpose() const {
    Matrix t(cols_, rows_);
    for (std::size_t i = 0; i < rows_; ++i) {
        for (std::size_t j = 0; j < cols_; ++j) {
            t(j, i) = (*this)(i, j);
        }
    }
    return t;
}

std::vector<double> Matrix::operator*(const std::vector<double>& v) const {
    if (v.size() != cols_) {
        throw std::invalid_argument("mathx::Matrix: vector size mismatch");
    }
    std::vector<double> out(rows_, 0.0);
    for (std::size_t i = 0; i < rows_; ++i) {
        double acc = 0.0;
        for (std::size_t j = 0; j < cols_; ++j) {
            acc += (*this)(i, j) * v[j];
        }
        out[i] = acc;
    }
    return out;
}

namespace {
void check_same_shape(const Matrix& a, const Matrix& b) {
    if (a.rows() != b.rows() || a.cols() != b.cols()) {
        throw std::invalid_argument("mathx::Matrix: shape mismatch");
    }
}

//! Scalar (portable) matmul: C(i,j) += A(i,k) * B(k,j), ikj order.
void matmul_scalar(const double* a, const double* b, double* c,
                   std::size_t M, std::size_t K, std::size_t N) {
    for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t k = 0; k < K; ++k) {
            const double v = a[i * K + k];
            if (v == 0.0) continue;
            const double* bk = b + k * N;
            double* ri = c + i * N;
            std::size_t j = 0;
            for (; j + 4 <= N; j += 4) {
                ri[j]     += v * bk[j];
                ri[j + 1] += v * bk[j + 1];
                ri[j + 2] += v * bk[j + 2];
                ri[j + 3] += v * bk[j + 3];
            }
            for (; j < N; ++j) {
                ri[j] += v * bk[j];
            }
        }
    }
}

#if defined(__AVX2__)
//! AVX2/FMA micro-kernel: computes a 4 x BN tile of C from A(4xKC) and
//! B(KC x N, row-major, stride ldb). Accumulators live in YMM registers.
template <std::size_t KC, std::size_t BN>
inline void avx2_microkernel(const double* a, std::size_t lda,
                             const double* b, std::size_t ldb,
                             double* c, std::size_t ldc) {
    constexpr std::size_t NV = BN / 4;
    __m256d acc[4][NV];
#pragma GCC unroll 1
    for (std::size_t i = 0; i < 4; ++i) {
#pragma GCC unroll 1
        for (std::size_t j = 0; j < NV; ++j) {
            acc[i][j] = _mm256_loadu_pd(c + i * ldc + 4 * j);
        }
    }
#pragma GCC unroll 1
    for (std::size_t k = 0; k < KC; ++k) {
        const double* bp = b + k * ldb;
        for (std::size_t i = 0; i < 4; ++i) {
            const __m256d av = _mm256_set1_pd(a[i * lda + k]);
#pragma GCC unroll 1
            for (std::size_t j = 0; j < NV; ++j) {
                acc[i][j] = _mm256_fmadd_pd(av, _mm256_loadu_pd(bp + 4 * j),
                                            acc[i][j]);
            }
        }
    }
#pragma GCC unroll 1
    for (std::size_t i = 0; i < 4; ++i) {
#pragma GCC unroll 1
        for (std::size_t j = 0; j < NV; ++j) {
            _mm256_storeu_pd(c + i * ldc + 4 * j, acc[i][j]);
        }
    }
}

//! Blocked AVX2 matmul with i/k/j cache blocking.
void matmul_avx2(const double* a, const double* b, double* c,
                 std::size_t M, std::size_t K, std::size_t N) {
    constexpr std::size_t KC = 256;  // B panel width (fits L2)
    constexpr std::size_t BN = 8;    // B micro-panel width
    constexpr std::size_t MB = 64;   // C rows per task

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i0 = 0; i0 < M; i0 += MB) {
        const std::size_t i1 = std::min(i0 + MB, M);
        // Process rows in groups of 4 (micro-kernel granularity).
        for (std::size_t i = i0; i < i1;) {
            const bool full_group = (i + 3 < i1);
            const std::size_t group_rows = full_group ? 4 : (i1 - i);

            for (std::size_t k0 = 0; k0 < K; k0 += KC) {
                const std::size_t k1 = std::min(k0 + KC, K);
                const bool full_k = (k1 - k0 == KC);

                if (full_group && full_k) {
                    // AVX path: 4x8 micro-kernel over full k-block.
                    std::size_t j0 = 0;
                    for (; j0 + BN <= N; j0 += BN) {
                        avx2_microkernel<256, 8>(a + i * K + k0, K,
                                                 b + k0 * N + j0, N,
                                                 c + i * N + j0, N);
                    }
                    // Trailing columns (< BN).
                    for (; j0 < N; ++j0) {
                        double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
                        for (std::size_t k = k0; k < k1; ++k) {
                            const double bk = b[k * N + j0];
                            s0 += a[(i + 0) * K + k] * bk;
                            s1 += a[(i + 1) * K + k] * bk;
                            s2 += a[(i + 2) * K + k] * bk;
                            s3 += a[(i + 3) * K + k] * bk;
                        }
                        c[(i + 0) * N + j0] += s0;
                        c[(i + 1) * N + j0] += s1;
                        c[(i + 2) * N + j0] += s2;
                        c[(i + 3) * N + j0] += s3;
                    }
                } else {
                    // Scalar path: group_rows rows over [k0, k1) and all N.
                    for (std::size_t r = 0; r < group_rows; ++r) {
                        const double* ar = a + (i + r) * K;
                        double* cr = c + (i + r) * N;
                        for (std::size_t j = 0; j < N; ++j) {
                            double sum = 0.0;
                            for (std::size_t k = k0; k < k1; ++k) {
                                sum += ar[k] * b[k * N + j];
                            }
                            cr[j] += sum;
                        }
                    }
                }
            }
            i += group_rows;
        }
    }
}
#endif  // __AVX2__
} // namespace

Matrix operator+(const Matrix& a, const Matrix& b) {
    check_same_shape(a, b);
    Matrix r(a.rows(), a.cols());
    for (std::size_t i = 0; i < a.rows(); ++i)
        for (std::size_t j = 0; j < a.cols(); ++j)
            r(i, j) = a(i, j) + b(i, j);
    return r;
}

Matrix operator-(const Matrix& a, const Matrix& b) {
    check_same_shape(a, b);
    Matrix r(a.rows(), a.cols());
    for (std::size_t i = 0; i < a.rows(); ++i)
        for (std::size_t j = 0; j < a.cols(); ++j)
            r(i, j) = a(i, j) - b(i, j);
    return r;
}

Matrix operator*(const Matrix& a, const Matrix& b) {
    if (a.cols() != b.rows()) {
        throw std::invalid_argument("mathx::Matrix: dimension mismatch for product");
    }
    const std::size_t M = a.rows();
    const std::size_t K = a.cols();
    const std::size_t N = b.cols();

    if (M == 0 || K == 0 || N == 0) {
        return Matrix(M, N);
    }

    Matrix r(M, N);
    const double* ad = a.data();
    const double* bd = b.data();
    double* rd = r.data();

    // Small matrices: avoid parallel/AVX overhead.
    if (M * K * N < 200000) {
        matmul_scalar(ad, bd, rd, M, K, N);
        return r;
    }

#if defined(__AVX2__)
    // Cache-blocked AVX2/FMA kernel (OpenMP inside).
    matmul_avx2(ad, bd, rd, M, K, N);
#else
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i = 0; i < M; ++i) {
        matmul_scalar(ad + i * K, bd, rd + i * N, 1, K, N);
    }
#endif
    return r;
}

Matrix operator*(double s, const Matrix& a) {
    Matrix r(a.rows(), a.cols());
    for (std::size_t i = 0; i < a.rows(); ++i)
        for (std::size_t j = 0; j < a.cols(); ++j)
            r(i, j) = s * a(i, j);
    return r;
}

std::vector<double> solve(const Matrix& A, const std::vector<double>& b) {
    if (A.rows() != A.cols() || A.cols() != b.size()) {
        throw std::invalid_argument("mathx::solve: dimension mismatch");
    }
    const std::size_t n = A.rows();
    if (n == 0) return {};

    // Flat (row-major) copy of A.
    std::vector<double> lu(A.data(), A.data() + n * n);
    std::vector<double> rhs = b;

    // LU decomposition with partial pivoting (in place).
    for (std::size_t col = 0; col < n; ++col) {
        std::size_t piv = col;
        double maxv = std::fabs(lu[col * n + col]);
        for (std::size_t r = col + 1; r < n; ++r) {
            if (std::fabs(lu[r * n + col]) > maxv) {
                maxv = std::fabs(lu[r * n + col]);
                piv = r;
            }
        }
        if (maxv < 1e-14) {
            throw std::runtime_error("mathx::solve: singular matrix");
        }
        if (piv != col) {
            for (std::size_t c = 0; c < n; ++c) {
                std::swap(lu[piv * n + c], lu[col * n + c]);
            }
            std::swap(rhs[piv], rhs[col]);
        }
        const double diag = lu[col * n + col];
        for (std::size_t r = col + 1; r < n; ++r) {
            const double factor = lu[r * n + col] / diag;
            if (factor == 0.0) continue;
            lu[r * n + col] = factor;  // store L part
            for (std::size_t c = col + 1; c < n; ++c) {
                lu[r * n + c] -= factor * lu[col * n + c];
            }
        }
    }

    // Forward substitution on rhs (L y = P b).
    for (std::size_t i = 0; i < n; ++i) {
        double sum = rhs[i];
        for (std::size_t j = 0; j < i; ++j) {
            sum -= lu[i * n + j] * rhs[j];
        }
        rhs[i] = sum;
    }
    // Back substitution (U x = y).
    for (std::size_t i = n; i-- > 0;) {
        double sum = rhs[i];
        for (std::size_t j = i + 1; j < n; ++j) {
            sum -= lu[i * n + j] * rhs[j];
        }
        rhs[i] = sum / lu[i * n + i];
    }
    return rhs;
}

double det(const Matrix& A) {
    if (A.rows() != A.cols()) {
        throw std::invalid_argument("mathx::det: matrix must be square");
    }
    const std::size_t n = A.rows();
    if (n == 0) return 1.0;

    std::vector<double> lu(A.data(), A.data() + n * n);

    double d = 1.0;
    for (std::size_t col = 0; col < n; ++col) {
        std::size_t piv = col;
        double maxv = std::fabs(lu[col * n + col]);
        for (std::size_t r = col + 1; r < n; ++r) {
            if (std::fabs(lu[r * n + col]) > maxv) {
                maxv = std::fabs(lu[r * n + col]);
                piv = r;
            }
        }
        if (maxv < 1e-14) return 0.0;
        if (piv != col) {
            for (std::size_t c = 0; c < n; ++c) {
                std::swap(lu[piv * n + c], lu[col * n + c]);
            }
            d = -d;
        }
        const double diag = lu[col * n + col];
        d *= diag;
        for (std::size_t r = col + 1; r < n; ++r) {
            const double factor = lu[r * n + col] / diag;
            if (factor == 0.0) continue;
            for (std::size_t c = col; c < n; ++c) {
                lu[r * n + c] -= factor * lu[col * n + c];
            }
        }
    }
    return d;
}

} // namespace mathx
