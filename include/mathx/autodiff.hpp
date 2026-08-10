#ifndef MATHX_AUTODIFF_HPP
#define MATHX_AUTODIFF_HPP

#include "export.hpp"
#include <cmath>
#include <functional>
#include <utility>

namespace mathx {

//! Dual number for forward-mode automatic differentiation.
//! Stores value v and derivative dv.
struct MATHX_API Dual {
    double v;
    double dv;

    Dual() : v(0.0), dv(0.0) {}
    Dual(double vv) : v(vv), dv(0.0) {}
    Dual(double vv, double dvv) : v(vv), dv(dvv) {}

    Dual operator+(const Dual& o) const { return Dual(v + o.v, dv + o.dv); }
    Dual operator-(const Dual& o) const { return Dual(v - o.v, dv - o.dv); }
    Dual operator*(const Dual& o) const { return Dual(v * o.v, dv * o.v + v * o.dv); }
    Dual operator/(const Dual& o) const { return Dual(v / o.v, (dv * o.v - v * o.dv) / (o.v * o.v)); }

    Dual operator+(double s) const { return Dual(v + s, dv); }
    Dual operator-(double s) const { return Dual(v - s, dv); }
    Dual operator*(double s) const { return Dual(v * s, dv * s); }
    Dual operator/(double s) const { return Dual(v / s, dv / s); }

    Dual& operator+=(const Dual& o) { v += o.v; dv += o.dv; return *this; }
    Dual& operator-=(const Dual& o) { v -= o.v; dv -= o.dv; return *this; }
    Dual& operator*=(const Dual& o) { dv = dv * o.v + v * o.dv; v *= o.v; return *this; }
    Dual& operator/=(const Dual& o) { dv = (dv * o.v - v * o.dv) / (o.v * o.v); v /= o.v; return *this; }
};

inline Dual operator+(double s, const Dual& d) { return d + s; }
inline Dual operator-(double s, const Dual& d) { return Dual(s - d.v, -d.dv); }
inline Dual operator*(double s, const Dual& d) { return d * s; }
inline Dual operator/(double s, const Dual& d) { return Dual(s / d.v, -s * d.dv / (d.v * d.v)); }

// Math functions on dual numbers.
inline Dual sin(const Dual& d) { return Dual(std::sin(d.v), std::cos(d.v) * d.dv); }
inline Dual cos(const Dual& d) { return Dual(std::cos(d.v), -std::sin(d.v) * d.dv); }
inline Dual tan(const Dual& d) { const double c = std::cos(d.v); return Dual(std::tan(d.v), d.dv / (c * c)); }
inline Dual exp(const Dual& d) { return Dual(std::exp(d.v), std::exp(d.v) * d.dv); }
inline Dual log(const Dual& d) { return Dual(std::log(d.v), d.dv / d.v); }
inline Dual log10(const Dual& d) { return Dual(std::log10(d.v), d.dv / (d.v * std::log(10.0))); }
inline Dual sqrt(const Dual& d) { const double s = std::sqrt(d.v); return Dual(s, d.dv / (2.0 * s)); }
inline Dual pow(const Dual& d, double e) { return Dual(std::pow(d.v, e), e * std::pow(d.v, e - 1.0) * d.dv); }
inline Dual pow(const Dual& d, const Dual& e) {
    const double v = std::pow(d.v, e.v);
    return Dual(v, v * (e.dv * std::log(d.v) + e.v * d.dv / d.v));
}
inline Dual atan(const Dual& d) { return Dual(std::atan(d.v), d.dv / (1.0 + d.v * d.v)); }
inline Dual asin(const Dual& d) { return Dual(std::asin(d.v), d.dv / std::sqrt(1.0 - d.v * d.v)); }
inline Dual acos(const Dual& d) { return Dual(std::acos(d.v), -d.dv / std::sqrt(1.0 - d.v * d.v)); }
inline Dual sinh(const Dual& d) { return Dual(std::sinh(d.v), std::cosh(d.v) * d.dv); }
inline Dual cosh(const Dual& d) { return Dual(std::cosh(d.v), std::sinh(d.v) * d.dv); }
inline Dual tanh(const Dual& d) { const double t = std::tanh(d.v); return Dual(t, d.dv * (1.0 - t * t)); }
inline Dual abs(const Dual& d) { return Dual(std::fabs(d.v), (d.v < 0.0 ? -1.0 : 1.0) * d.dv); }
inline Dual hypot(const Dual& a, const Dual& b) {
    const double h = std::hypot(a.v, b.v);
    return Dual(h, (a.v * a.dv + b.v * b.dv) / h);
}

//! Evaluate f at x and its first derivative using forward-mode AD.
//! The functor must accept and return Dual numbers.
template <typename F>
std::pair<double, double> derivative(const F& f, double x) {
    const Dual r = f(Dual(x, 1.0));
    return {r.v, r.dv};
}

//! Second derivative via two AD passes.
template <typename F>
std::pair<double, std::pair<double, double>> derivatives(const F& f, double x) {
    const double eps = 1e-8;
    const auto a = derivative(f, x - eps);
    const auto b = derivative(f, x + eps);
    return {a.first, {b.second, (b.second - a.second) / (2.0 * eps)}};
}

} // namespace mathx

#endif
