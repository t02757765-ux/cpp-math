#include <cassert>
#include <cmath>
#include <iostream>

#include "mathx/matrix.hpp"

int main() {
    mathx::Matrix A(2, 2, {4.0, 2.0, 1.0, 3.0});
    assert(std::fabs(A(0, 0) - 4.0) < 1e-12);

    // det(A) = 4*3 - 2*1 = 10
    assert(std::fabs(mathx::det(A) - 10.0) < 1e-9);

    // Solve A x = b,  b = {6, 4}  =>  x = {1, 1}
    const std::vector<double> b{6.0, 4.0};
    const auto x = mathx::solve(A, b);
    assert(std::fabs(x[0] - 1.0) < 1e-9);
    assert(std::fabs(x[1] - 1.0) < 1e-9);

    // A * A
    const auto A2 = A * A;
    assert(std::fabs(A2(0, 0) - 18.0) < 1e-9);
    assert(std::fabs(A2(0, 1) - 14.0) < 1e-9);
    assert(std::fabs(A2(1, 0) - 7.0) < 1e-9);
    assert(std::fabs(A2(1, 1) - 11.0) < 1e-9);

    // transpose
    const auto At = A.transpose();
    assert(std::fabs(At(0, 1) - 1.0) < 1e-9);
    assert(std::fabs(At(1, 0) - 2.0) < 1e-9);

    // matrix-vector product
    const std::vector<double> v{1.0, 2.0};
    const auto av = A * v;
    assert(std::fabs(av[0] - 8.0) < 1e-9);
    assert(std::fabs(av[1] - 7.0) < 1e-9);

    // singular matrix throws
    mathx::Matrix S(2, 2, {1.0, 2.0, 2.0, 4.0});
    bool threw = false;
    try {
        mathx::solve(S, std::vector<double>{1.0, 1.0});
    } catch (const std::runtime_error&) {
        threw = true;
    }
    assert(threw);

    std::cout << "matrix tests passed\n";
    return 0;
}
