#include "mathx/polynomial.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace mathx {

Polynomial::Polynomial(std::vector<double> coeffs) : coeffs_(std::move(coeffs)) {
    trim();
}

int Polynomial::degree() const {
    if (coeffs_.empty()) return -1;
    const std::size_t n = coeffs_.size() - 1;
    for (std::size_t i = n + 1; i-- > 0;) {
        if (coeffs_[i] != 0.0) return static_cast<int>(i);
    }
    return -1;
}

double Polynomial::leading() const {
    const int d = degree();
    return d < 0 ? 0.0 : coeffs_[static_cast<std::size_t>(d)];
}

double Polynomial::coeff(std::size_t k) const {
    return k < coeffs_.size() ? coeffs_[k] : 0.0;
}

void Polynomial::trim() {
    while (!coeffs_.empty() && coeffs_.back() == 0.0) {
        coeffs_.pop_back();
    }
}

double Polynomial::operator()(double x) const {
    double r = 0.0;
    for (auto it = coeffs_.rbegin(); it != coeffs_.rend(); ++it) {
        r = r * x + *it;
    }
    return r;
}

std::complex<double> Polynomial::operator()(const std::complex<double>& x) const {
    std::complex<double> r(0.0, 0.0);
    for (auto it = coeffs_.rbegin(); it != coeffs_.rend(); ++it) {
        r = r * x + std::complex<double>(*it, 0.0);
    }
    return r;
}

Polynomial Polynomial::derivative() const {
    if (coeffs_.size() <= 1) return Polynomial();
    std::vector<double> d(coeffs_.size() - 1);
    for (std::size_t i = 1; i < coeffs_.size(); ++i) {
        d[i - 1] = coeffs_[i] * static_cast<double>(i);
    }
    return Polynomial(std::move(d));
}

Polynomial Polynomial::integral() const {
    std::vector<double> d(coeffs_.size() + 1, 0.0);
    for (std::size_t i = 0; i < coeffs_.size(); ++i) {
        d[i + 1] = coeffs_[i] / static_cast<double>(i + 1);
    }
    return Polynomial(std::move(d));
}

namespace {
constexpr double kPi = 3.141592653589793238462643383279502884;
} // namespace

std::vector<double> Polynomial::real_roots() const {
    const int d = degree();
    std::vector<double> roots;
    if (d < 0) return roots;
    if (d == 0) return roots;
    if (d == 1) {
        roots.push_back(-coeffs_[0] / coeffs_[1]);
        return roots;
    }
    if (d == 2) {
        const double a = coeffs_[2], b = coeffs_[1], c = coeffs_[0];
        const double disc = b * b - 4.0 * a * c;
        if (disc >= 0.0) {
            const double sq = std::sqrt(disc);
            roots.push_back((-b - sq) / (2.0 * a));
            if (disc > 0.0) roots.push_back((-b + sq) / (2.0 * a));
        }
        return roots;
    }
    if (d == 3) {
        // Depressed cubic via Cardano's formula.
        const double a = coeffs_[3], b = coeffs_[2], c = coeffs_[1], dd = coeffs_[0];
        const double p = (3.0 * a * c - b * b) / (3.0 * a * a);
        const double q = (2.0 * b * b * b - 9.0 * a * b * c + 27.0 * a * a * dd)
                       / (27.0 * a * a * a);
        const double disc = q * q / 4.0 + p * p * p / 27.0;
        const double shift = b / (3.0 * a);
        if (disc > 0.0) {
            const double u = std::cbrt(-q / 2.0 + std::sqrt(disc));
            const double v = std::cbrt(-q / 2.0 - std::sqrt(disc));
            roots.push_back(u + v - shift);
        } else if (disc < 0.0) {
            const double r = std::sqrt(-p * p * p / 27.0);
            const double theta = std::acos(-q / (2.0 * r));
            for (int k = 0; k < 3; ++k) {
                roots.push_back(2.0 * std::cbrt(r)
                                * std::cos((theta + 2.0 * kPi) / 3.0) - shift);
            }
        } else {
            const double u = std::cbrt(-q / 2.0);
            roots.push_back(2.0 * u - shift);
            roots.push_back(-u - shift);
        }
        return roots;
    }
    // Higher degree: use Durand-Kerner then filter real parts.
    const auto cs = this->roots();
    for (const auto& z : cs) {
        if (std::fabs(z.imag()) < 1e-8 * std::max(1.0, std::fabs(z.real()))) {
            roots.push_back(z.real());
        }
    }
    std::sort(roots.begin(), roots.end());
    return roots;
}

std::vector<std::complex<double>> Polynomial::roots() const {
    std::vector<std::complex<double>> result;
    const int d = degree();
    if (d < 1) return result;

    // Divide out zeros (constant == 0).
    Polynomial p = *this;
    int zero_count = 0;
    while (p.coeff(0) == 0.0 && p.degree() > 0) {
        ++zero_count;
        p.coeffs_.erase(p.coeffs_.begin());
    }
    for (int i = 0; i < zero_count; ++i) result.emplace_back(0.0, 0.0);

    const int nd = p.degree();
    if (nd == 0) return result;
    if (nd == 1) {
        result.emplace_back(-p.coeff(0) / p.coeff(1), 0.0);
        return result;
    }

    // Durand-Kerner simultaneous iteration.
    const std::size_t n = static_cast<std::size_t>(nd);
    const double scale = 1.0 + std::fabs(p.leading());
    std::vector<std::complex<double>> z(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double angle = 2.0 * kPi
                           * (static_cast<double>(i) + 0.5) / static_cast<double>(n);
        z[i] = std::complex<double>(0.4 * std::cos(angle), 0.4 * std::sin(angle)) * scale;
    }
    for (int iter = 0; iter < 300; ++iter) {
        bool converged = true;
        for (std::size_t i = 0; i < n; ++i) {
            std::complex<double> denom(1.0, 0.0);
            for (std::size_t j = 0; j < n; ++j) {
                if (i == j) continue;
                denom *= (z[i] - z[j]);
            }
            if (std::abs(denom) < 1e-300) {
                denom = std::complex<double>(1e-300, 0.0);
            }
            const std::complex<double> fz = p(z[i]);
            const std::complex<double> dz = fz / denom;
            z[i] -= dz;
            if (std::abs(dz) > 1e-13 * std::max(1.0, std::abs(z[i]))) {
                converged = false;
            }
        }
        if (converged) break;
    }
    for (const auto& zz : z) result.push_back(zz);
    return result;
}

Polynomial operator+(const Polynomial& a, const Polynomial& b) {
    std::vector<double> c(std::max(a.size(), b.size()), 0.0);
    for (std::size_t i = 0; i < c.size(); ++i) {
        c[i] = a.coeff(i) + b.coeff(i);
    }
    return Polynomial(std::move(c));
}

Polynomial operator-(const Polynomial& a, const Polynomial& b) {
    std::vector<double> c(std::max(a.size(), b.size()), 0.0);
    for (std::size_t i = 0; i < c.size(); ++i) {
        c[i] = a.coeff(i) - b.coeff(i);
    }
    return Polynomial(std::move(c));
}

Polynomial operator*(const Polynomial& a, const Polynomial& b) {
    const std::size_t n = a.size(), m = b.size();
    std::vector<double> c(n + m, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        if (a.coeff(i) == 0.0) continue;
        for (std::size_t j = 0; j < m; ++j) {
            c[i + j] += a.coeff(i) * b.coeff(j);
        }
    }
    return Polynomial(std::move(c));
}

Polynomial operator*(double s, const Polynomial& p) {
    std::vector<double> c(p.size());
    for (std::size_t i = 0; i < c.size(); ++i) c[i] = s * p.coeff(i);
    return Polynomial(std::move(c));
}

std::vector<double> polynomial_divmod(const Polynomial& a, const Polynomial& b,
                                      Polynomial& quotient, Polynomial& remainder) {
    const int da = a.degree();
    const int db = b.degree();
    if (db < 0) {
        throw std::domain_error("mathx::polynomial_divmod: division by zero polynomial");
    }
    if (da < db) {
        quotient = Polynomial();
        remainder = a;
        return {};
    }
    std::vector<double> q(static_cast<std::size_t>(da - db + 1), 0.0);
    std::vector<double> r = a.coeffs_;
    for (int k = da - db; k >= 0; --k) {
        const double lead = b.leading();
        q[static_cast<std::size_t>(k)] = r[static_cast<std::size_t>(k + db)] / lead;
        for (int j = 0; j <= db; ++j) {
            r[static_cast<std::size_t>(k + j)] -= q[static_cast<std::size_t>(k)] * b.coeff(static_cast<std::size_t>(j));
        }
    }
    quotient = Polynomial(std::move(q));
    remainder = Polynomial(std::move(r));
    return {};
}

Polynomial operator/(const Polynomial& a, const Polynomial& b) {
    Polynomial q, r;
    polynomial_divmod(a, b, q, r);
    return q;
}

Polynomial operator%(const Polynomial& a, const Polynomial& b) {
    Polynomial q, r;
    polynomial_divmod(a, b, q, r);
    return r;
}

} // namespace mathx
