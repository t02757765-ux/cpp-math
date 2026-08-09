#include "mathx/matrix.hpp"

#include <algorithm>
#include <cmath>

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
    Matrix r(a.rows(), b.cols());
    for (std::size_t i = 0; i < a.rows(); ++i) {
        for (std::size_t k = 0; k < a.cols(); ++k) {
            const double v = a(i, k);
            if (v == 0.0) continue;
            for (std::size_t j = 0; j < b.cols(); ++j) {
                r(i, j) += v * b(k, j);
            }
        }
    }
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
    // Build augmented matrix.
    std::vector<std::vector<double>> a(n, std::vector<double>(n + 1));
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            a[i][j] = A(i, j);
        }
        a[i][n] = b[i];
    }

    for (std::size_t col = 0; col < n; ++col) {
        // Partial pivoting.
        std::size_t piv = col;
        double maxv = std::fabs(a[col][col]);
        for (std::size_t r = col + 1; r < n; ++r) {
            if (std::fabs(a[r][col]) > maxv) {
                maxv = std::fabs(a[r][col]);
                piv = r;
            }
        }
        if (maxv < 1e-14) {
            throw std::runtime_error("mathx::solve: singular matrix");
        }
        if (piv != col) {
            std::swap(a[piv], a[col]);
        }
        // Eliminate.
        const double diag = a[col][col];
        for (std::size_t r = col + 1; r < n; ++r) {
            const double factor = a[r][col] / diag;
            if (factor == 0.0) continue;
            for (std::size_t c = col; c <= n; ++c) {
                a[r][c] -= factor * a[col][c];
            }
        }
    }

    // Back substitution.
    std::vector<double> x(n, 0.0);
    for (std::size_t i = n; i-- > 0;) {
        double sum = a[i][n];
        for (std::size_t j = i + 1; j < n; ++j) {
            sum -= a[i][j] * x[j];
        }
        x[i] = sum / a[i][i];
    }
    return x;
}

double det(const Matrix& A) {
    if (A.rows() != A.cols()) {
        throw std::invalid_argument("mathx::det: matrix must be square");
    }
    const std::size_t n = A.rows();
    std::vector<std::vector<double>> a(n, std::vector<double>(n));
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            a[i][j] = A(i, j);

    double d = 1.0;
    for (std::size_t col = 0; col < n; ++col) {
        std::size_t piv = col;
        double maxv = std::fabs(a[col][col]);
        for (std::size_t r = col + 1; r < n; ++r) {
            if (std::fabs(a[r][col]) > maxv) {
                maxv = std::fabs(a[r][col]);
                piv = r;
            }
        }
        if (maxv < 1e-14) return 0.0;
        if (piv != col) {
            std::swap(a[piv], a[col]);
            d = -d;
        }
        const double diag = a[col][col];
        d *= diag;
        for (std::size_t r = col + 1; r < n; ++r) {
            const double factor = a[r][col] / diag;
            for (std::size_t c = col; c < n; ++c) {
                a[r][c] -= factor * a[col][c];
            }
        }
    }
    return d;
}

} // namespace mathx
