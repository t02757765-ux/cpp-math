# Sürüm Geçmişi

## 0.3.0 — 2026-08-10

### Yeni modüller
- **`hypergeometric.hpp`** — ₁F₁, ₂F₁, genelleştirilmiş pFq + tam eliptik integraller K(m), E(m)
- **`autodiff.hpp`** — Header-only ileri mod otomatik türev (`Dual` sayıları, `derivative`, `derivatives`)
- **Python bağlamaları** (`bindings/python/`) — pybind11 tabanlı `mathx` modülü: lambert, özel fonksiyonlar, Bessel, hipergeometrik, FFT, integrasyon, kök bulma, ODE, istatistik, `Matrix`, `Polynomial`

### Performans
- `Matrix::operator*` — OpenMP paralelleştirme + 4x manuel açılım (~50 GFLOP/s, 512×512)
- `dft`/`idft` — büyük girdiler için OpenMP
- CMake: `MATHX_ENABLE_OPENMP` ve `MATHX_ENABLE_NATIVE` seçenekleri (`-march=native`, `-funroll-loops`)
- `benchmarks/` — C++ ve Python (NumPy/scipy karşılaştırmalı) benchmark'lar

### Düzeltmeler
- `elliptic_e` ikinci tür eliptik integral doğru kuvvet serisiyle yeniden yazıldı
- `Matrix` için `data()` erişimcisi eklendi (paralel çarpım için)

### Testler
- `test_hypergeometric`, `test_autodiff` eklendi — toplam 12 test programı

## 0.2.0 — 2026-08-10

### Yeni modüller
- **`constants.hpp`** — 18 matematiksel sabit (`pi`, `e`, `euler_gamma`, `catalan`, `apery`, ...)
- **`bessel.hpp`** — Bessel fonksiyonları J, Y, I, K (mertebe 0, 1, n) + Airy Ai/Bi
- **`integration.hpp`** — Adaptif Simpson, trapez, bileşik Simpson, Gauss-Legendre (16/32 nokta), Romberg
- **`roots.hpp`** — Bisection, regula-falsi, secant, Newton, Newton-FD, Brent
- **`polynomial.hpp`** — `Polynomial` sınıfı: değerlendirme, türev, integral, kökler (derece ≤ 3 kesin, ≥ 4 Durand-Kerner), bölme
- **`statistics.hpp`** — Ortalama, varyans, standart sapma, medyan, yüzdelik, kovaryans, korelasyon, çarpıklık, basıklık
- **`ode.hpp`** — Euler, Heun, RK4, adaptif Dormand-Prince (RK45)

### Düzeltmeler
- Kompleks Lambert W artık Newton yinelemesi kullanıyor (Halley bazen yanlış dala yakınsıyordu)
- K₀/K₁/Y₀ serileri doğru harmonik sayı formülleriyle yeniden yazıldı
- RK45 katsayıları düzeltildi (k5 zaman noktası 8/9) — metod artık 5. derece
- `erf(0)` artık tam 0 döndürüyor
- `Polynomial::real_roots` için Durand-Kerner kök filtreleme

### Testler
- 10 test programı (`tests/`), CMake `ctest` ile çalıştırılabilir
- Bessel, integral, kök, polinom, istatistik, ODE testleri eklendi

## 0.1.0 — 2026-08-10

İlk sürüm.

- Lambert W (W₀, W₋₁, kompleks)
- DFT / IDFT / FFT / IFFT
- Gamma, digamma, erf, erfc, norm_cdf, zeta, beta
- `Matrix` sınıfı, `solve`, `det`
- CMake yapılandırması (.dll/.so hedefi)
- 4 test programı
