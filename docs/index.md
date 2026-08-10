# cpp-math Dökümantasyonu

**cpp-math** — C++ ile yazılmış, NumPy ve SymPy'nin birleşimini hedefleyen ileri düzey matematik kütüphanesi. Özel fonksiyonlar, dönüşümler, sayısal analiz ve lineer cebir modülleri içerir; `.dll` / `.so` olarak derlenir.

## İçindekiler

1. [Başlangıç](getting-started.md) — kurulum, derleme, örnekler
2. [API Referansı](api-reference.md) — tüm fonksiyonların listesi
3. [Modüller](modules.md) — her modülün ayrıntılı açıklaması
4. [C++ Entegrasyonu](integration.md) — CMake ile başka projelerde kullanım
5. [Sürüm Geçmişi](changelog.md)

## Özellik Özeti

| Alan | Modül | İçerik |
|---|---|---|
| Özel fonksiyonlar | `special.hpp` | Gamma, digamma, erf, erfc, zeta, beta, norm CDF |
| Lambert W | `lambert.hpp` | W₀, W₋₁, kompleks W |
| Fourier | `fourier.hpp` | DFT, IDFT, FFT, IFFT |
| Bessel & Airy | `bessel.hpp` | J, Y, I, K türleri, Ai, Bi |
| Sayısal integral | `integration.hpp` | Simpson, trapez, Gauss-Legendre, Romberg |
| Kök bulma | `roots.hpp` | Bisection, regula-falsi, secant, Newton, Brent |
| Lineer cebir | `matrix.hpp` | Matris, çözücü, determinant |
| Polinomlar | `polynomial.hpp` | Dereceler, kökler, türev, integral, bölme |
| İstatistik | `statistics.hpp` | Ortalama, varyans, korelasyon, çeyreklikler |
| ODE çözücüler | `ode.hpp` | Euler, Heun, RK4, adaptif RK45 |
| Sabitler | `constants.hpp` | π, e, γ, altın oran ve daha fazlası |

## Hızlı Örnek

```cpp
#include <mathx/mathx.hpp>
#include <iostream>

int main() {
    // Lambert W
    std::cout << mathx::lambert_w(1.0) << "\n";          // 0.567143

    // FFT
    auto X = mathx::fft(std::vector<double>{1, 2, 3, 4});

    // Sayısal integral: int_0^1 e^x dx
    double area = mathx::integrate_simpson(
        [](double x) { return std::exp(x); }, 0.0, 1.0);
    std::cout << area << "\n";                           // 1.71828

    // Kök bulma: x^2 - 2 = 0
    double r = mathx::brent([](double x) { return x*x - 2.0; }, 1.0, 2.0);

    // Polinom kökleri
    mathx::Polynomial p({2.0, -3.0, 1.0});               // x^2 - 3x + 2
    auto roots = p.real_roots();                         // {1, 2}
}
```

## Modül Haritası

```
include/mathx/
├── export.hpp       -- DLL/SO ihracat makroları
├── mathx.hpp        -- ana başlık (her şeyi dahil eder)
├── constants.hpp    -- matematiksel sabitler
├── lambert.hpp      -- Lambert W
├── fourier.hpp      -- Fourier dönüşümleri
├── special.hpp      -- özel fonksiyonlar
├── bessel.hpp       -- Bessel ve Airy fonksiyonları
├── integration.hpp  -- sayısal integrasyon
├── roots.hpp        -- kök bulma algoritmaları
├── matrix.hpp       -- lineer cebir
├── polynomial.hpp   -- polinom cebri
├── statistics.hpp   -- istatistik
└── ode.hpp          -- diferansiyel denklem çözücüler

src/                 -- uygulama dosyaları
tests/               -- test programları
```
