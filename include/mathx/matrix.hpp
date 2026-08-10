#ifndef MATHX_MATRIX_HPP
#define MATHX_MATRIX_HPP

#include "export.hpp"
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace mathx {

//! Simple dense matrix with row-major storage.
class MATHX_API Matrix {
public:
    Matrix() = default;
    Matrix(std::size_t rows, std::size_t cols,
           double value = 0.0);
    Matrix(std::size_t rows, std::size_t cols,
           std::initializer_list<double> values);

    std::size_t rows() const { return rows_; }
    std::size_t cols() const { return cols_; }
    bool empty() const { return rows_ == 0 || cols_ == 0; }

    double& operator()(std::size_t i, std::size_t j);
    double  operator()(std::size_t i, std::size_t j) const;

    //! Raw access to the underlying row-major data buffer.
    double* data() { return data_.data(); }
    const double* data() const { return data_.data(); }

    Matrix transpose() const;

    //! Matrix-vector product.
    std::vector<double> operator*(const std::vector<double>& v) const;

private:
    std::size_t rows_ = 0;
    std::size_t cols_ = 0;
    std::vector<double> data_;
};

MATHX_API Matrix operator+(const Matrix& a, const Matrix& b);
MATHX_API Matrix operator-(const Matrix& a, const Matrix& b);
MATHX_API Matrix operator*(const Matrix& a, const Matrix& b);
MATHX_API Matrix operator*(double s, const Matrix& a);

//! Solve the linear system A x = b via Gaussian elimination with partial
//! pivoting. Throws std::runtime_error if A is singular.
MATHX_API std::vector<double> solve(const Matrix& A,
                                    const std::vector<double>& b);

//! Determinant of a square matrix.
MATHX_API double det(const Matrix& A);

} // namespace mathx

#endif
