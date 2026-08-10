// cpp-math Python bindings (pybind11).
// Build:
//   g++ -O3 -shared -fopenmp -std=c++17 -Iinclude \
//       -I$(pybind11)/include $(python3 -m pybind11 --includes) \
//       bindings/python/module.cpp src/*.cpp -o mathx$(python3-config --extension-suffix)
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>

#include <complex>
#include <vector>

#include <mathx/mathx.hpp>

namespace py = pybind11;

PYBIND11_MODULE(mathx, m) {
    m.doc() = "cpp-math: advanced C++ math library (NumPy/SymPy-inspired)";
    m.attr("__version__") = "0.3.0";

    // Lambert W
    m.def("lambert_w", [](double x, int branch) { return mathx::lambert_w(x, branch); },
          py::arg("x"), py::arg("branch") = 0, "Lambert W function (real, W0/W-1)");
    m.def("lambert_w_complex", [](std::complex<double> z) { return mathx::lambert_w(z); },
          "Lambert W function (complex, principal branch)");

    // Special functions
    m.def("gamma", &mathx::gamma);
    m.def("lgamma", &mathx::lgamma);
    m.def("digamma", &mathx::digamma);
    m.def("erf", &mathx::erf);
    m.def("erfc", &mathx::erfc);
    m.def("norm_cdf", &mathx::norm_cdf);
    m.def("zeta", &mathx::zeta);
    m.def("beta", &mathx::beta);

    // Bessel
    m.def("bessel_j0", &mathx::bessel_j0);
    m.def("bessel_j1", &mathx::bessel_j1);
    m.def("bessel_jn", &mathx::bessel_jn);
    m.def("bessel_y0", &mathx::bessel_y0);
    m.def("bessel_y1", &mathx::bessel_y1);
    m.def("bessel_yn", &mathx::bessel_yn);
    m.def("bessel_i0", &mathx::bessel_i0);
    m.def("bessel_i1", &mathx::bessel_i1);
    m.def("bessel_in", &mathx::bessel_in);
    m.def("bessel_k0", &mathx::bessel_k0);
    m.def("bessel_k1", &mathx::bessel_k1);
    m.def("bessel_kn", &mathx::bessel_kn);
    m.def("airy_ai", &mathx::airy_ai);
    m.def("airy_bi", &mathx::airy_bi);

    // Hypergeometric / elliptic
    m.def("hypergeometric_1f1", &mathx::hypergeometric_1f1);
    m.def("hypergeometric_2f1", &mathx::hypergeometric_2f1);
    m.def("elliptic_k", &mathx::elliptic_k);
    m.def("elliptic_e", &mathx::elliptic_e);

    // Fourier
    m.def("fft", [](std::vector<double> x) { return mathx::fft(x); });
    m.def("ifft", [](std::vector<std::complex<double>> x) { return mathx::ifft(x); });
    m.def("dft", [](std::vector<double> x) { return mathx::dft(x); });

    // Integration
    m.def("integrate_simpson", [](py::function f, double a, double b) {
        return mathx::integrate_simpson([&](double x) { return f(x).cast<double>(); }, a, b);
    });
    m.def("integrate_gauss_legendre", [](py::function f, double a, double b, int n) {
        return mathx::integrate_gauss_legendre([&](double x) { return f(x).cast<double>(); }, a, b, n);
    }, py::arg("f"), py::arg("a"), py::arg("b"), py::arg("n") = 16);

    // Roots
    m.def("bisect", [](py::function f, double a, double b) {
        return mathx::bisect([&](double x) { return f(x).cast<double>(); }, a, b);
    });
    m.def("newton", [](py::function f, double x0) {
        return mathx::newton_fd([&](double x) { return f(x).cast<double>(); }, x0);
    });
    m.def("brent", [](py::function f, double a, double b) {
        return mathx::brent([&](double x) { return f(x).cast<double>(); }, a, b);
    });

    // ODE
    m.def("ode_rk4", [](py::function f, double y0, double t0, double t1, int n) {
        auto sol = mathx::ode_rk4(
            [&](double t, double y) { return f(t, y).cast<double>(); },
            y0, t0, t1, n);
        py::list out;
        for (auto& p : sol) out.append(py::make_tuple(p.first, p.second));
        return out;
    }, py::arg("f"), py::arg("y0"), py::arg("t0"), py::arg("t1"), py::arg("n") = 1000);
    m.def("ode_rk45", [](py::function f, double y0, double t0, double t1, double tol) {
        auto sol = mathx::ode_rk45(
            [&](double t, double y) { return f(t, y).cast<double>(); },
            y0, t0, t1, tol);
        py::list out;
        for (auto& p : sol) out.append(py::make_tuple(p.first, p.second));
        return out;
    }, py::arg("f"), py::arg("y0"), py::arg("t0"), py::arg("t1"), py::arg("tol") = 1e-8);

    // Statistics
    m.def("mean", &mathx::mean);
    m.def("variance", &mathx::variance);
    m.def("stddev", &mathx::stddev);
    m.def("median", &mathx::median);
    m.def("percentile", &mathx::percentile);
    m.def("correlation", &mathx::correlation);
    m.def("skewness", &mathx::skewness);

    // Matrix
    py::class_<mathx::Matrix>(m, "Matrix")
        .def(py::init<>())
        .def(py::init([](py::list rows) {
            const std::size_t n = rows.size();
            if (n == 0) return mathx::Matrix();
            const std::size_t m = py::len(rows[0]);
            mathx::Matrix r(n, m);
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t j = 0; j < m; ++j)
                    r(i, j) = py::cast<double>(py::cast<py::list>(rows[i])[j]);
            return r;
        }))
        .def("shape", [](const mathx::Matrix& a) {
            return py::make_tuple(a.rows(), a.cols());
        })
        .def("__mul__", [](const mathx::Matrix& a, const mathx::Matrix& b) { return a * b; })
        .def("__rmul__", [](const mathx::Matrix& a, double s) { return s * a; })
        .def("__add__", [](const mathx::Matrix& a, const mathx::Matrix& b) { return a + b; })
        .def("__sub__", [](const mathx::Matrix& a, const mathx::Matrix& b) { return a - b; })
        .def("transpose", &mathx::Matrix::transpose)
        .def("det", &mathx::det)
        .def("__repr__", [](const mathx::Matrix& a) {
            std::string s = "Matrix(" + std::to_string(a.rows()) + "x"
                          + std::to_string(a.cols()) + ")";
            return s;
        });
    m.def("solve", [](const mathx::Matrix& A, const std::vector<double>& b) {
        return mathx::solve(A, b);
    });

    // Polynomial
    py::class_<mathx::Polynomial>(m, "Polynomial")
        .def(py::init([](std::vector<double> c) { return mathx::Polynomial(std::move(c)); }))
        .def("degree", &mathx::Polynomial::degree)
        .def("__call__", static_cast<double (mathx::Polynomial::*)(double) const>(&mathx::Polynomial::operator()))
        .def("derivative", &mathx::Polynomial::derivative)
        .def("integral", &mathx::Polynomial::integral)
        .def("real_roots", &mathx::Polynomial::real_roots)
        .def("roots", &mathx::Polynomial::roots);
}
