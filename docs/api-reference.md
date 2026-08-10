# API Referansı

Bu belge cpp-math'in tüm public API'sini listeler. Tüm fonksiyonlar `mathx` namespace'i altındadır.

## Sabitler (`constants.hpp`)

`namespace mathx::constants`

| Sabit | Değer (kısaltılmış) | Açıklama |
|---|---|---|
| `pi` | 3.14159... | Pi |
| `tau` | 6.28318... | 2π |
| `e` | 2.71828... | Euler sayısı |
| `sqrt2` | 1.41421... | √2 |
| `sqrt3` | 1.73205... | √3 |
| `sqrt_pi` | 1.77245... | √π |
| `phi` | 1.61803... | Altın oran |
| `ln2` | 0.69314... | ln 2 |
| `ln10` | 2.30258... | ln 10 |
| `euler_gamma` | 0.57721... | Euler-Mascheroni sabiti |
| `catalan` | 0.91596... | Catalan sabiti |
| `apery` | 1.20205... | ζ(3) |
| `feigenbaum_delta` | 4.66920... | Feigenbaum δ |
| `feigenbaum_alpha` | 2.50290... | Feigenbaum α |
| `khinchin` | 2.68545... | Khinchin sabiti |
| `glaisher` | 1.28242... | Glaisher-Kinkelin sabiti |
| `gauss` | 0.83462... | Gauss sabiti |
| `conway` | 1.30357... | Conway sabiti |

## Lambert W (`lambert.hpp`)

| Fonksiyon | İmza | Açıklama |
|---|---|---|
| `lambert_w` | `double (double x, int branch=0)` | Gerçel Lambert W; `branch=0` W₀, `branch=-1` W₋₁ |
| `lambert_w` | `std::complex<double> (const std::complex<double>&)` | Kompleks W₀ |
| `lambert_w_solve` | `double (double w, int branch=0)` | `x·eˣ = w` denklemini çözer |

**Davranış:**
- `x < -1/e` için `NaN`
- `x → 0⁻` ve `branch=-1` için `-∞`
- `lambert_w(-1/e, branch)` her iki dalda da `-1` verir

## Fourier (`fourier.hpp`)

| Fonksiyon | Açıklama |
|---|---|
| `dft(x)` | Naive O(n²) ayrık Fourier dönüşümü |
| `idft(X)` | Ters DFT |
| `fft(x)` | İteratif Cooley-Tukey, O(n log n); n, 2'nin kuvveti olmalı |
| `ifft(X)` | Ters FFT |
| `fft(vector<double>)` | Gerçel girdi için FFT |
| `dft(vector<double>)` | Gerçel girdi için DFT |

**Notlar:**
- FFT için girdi boyutu 2'nin kuvveti değilse `std::invalid_argument` fırlatılır
- FFT/DFT çıktıları aynı boyuttadır
- `ifft(fft(x)) ≈ x` döner

## Özel Fonksiyonlar (`special.hpp`)

| Fonksiyon | Açıklama |
|---|---|
| `gamma(x)` | Gamma fonksiyonu (Lanczos yaklaşımı) |
| `lgamma(x)` | log|Γ(x)| |
| `digamma(x)` | Digamma ψ(x) |
| `erf(x)` | Hata fonksiyonu |
| `erfc(x)` | Tamamlayıcı hata fonksiyonu |
| `norm_cdf(x)` | Standart normal CDF |
| `zeta(s)` | Riemann zeta fonksiyonu |
| `beta(x, y)` | Beta fonksiyonu B(x,y) |

## Bessel ve Airy (`bessel.hpp`)

| Fonksiyon | Açıklama |
|---|---|
| `bessel_j0(x)` | J₀(x) |
| `bessel_j1(x)` | J₁(x) |
| `bessel_jn(n, x)` | Jₙ(x) |
| `bessel_y0(x)` | Y₀(x) |
| `bessel_y1(x)` | Y₁(x) |
| `bessel_yn(n, x)` | Yₙ(x) |
| `bessel_i0(x)` | I₀(x) |
| `bessel_i1(x)` | I₁(x) |
| `bessel_in(n, x)` | Iₙ(x) |
| `bessel_k0(x)` | K₀(x) |
| `bessel_k1(x)` | K₁(x) |
| `bessel_kn(n, x)` | Kₙ(x) |
| `airy_ai(x)` | Airy Ai(x) |
| `airy_bi(x)` | Airy Bi(x) |

**Notlar:**
- Y ve K fonksiyonları `x ≤ 0` için `NaN` döner
- J, Y, I, K için yüksek mertebeler yineleme ilişkileriyle hesaplanır
- Büyük |x| için asimptotik açılımlar kullanılır
- Ai ve Bi seri açılımı ile hesaplanır (moderate |x| için)

## Sayısal İntegral (`integration.hpp`)

| Fonksiyon | Açıklama |
|---|---|
| `integrate_simpson(f, a, b, tol=1e-10, depth=50)` | Adaptif Simpson |
| `integrate_trapezoid(f, a, b, n=1000)` | Bileşik trapez |
| `integrate_simpson_n(f, a, b, n=1000)` | Bileşik Simpson |
| `integrate_gauss_legendre(f, a, b, n=16)` | Gauss-Legendre (n=16 veya 32) |
| `integrate_romberg(f, a, b, tol=1e-10)` | Romberg |

Tüm fonksiyonlar `std::function<double(double)>` alır.

## Kök Bulma (`roots.hpp`)

| Fonksiyon | Açıklama |
|---|---|
| `bisect(f, a, b, tol, max_iter)` | İkiye bölme |
| `regula_falsi(f, a, b, tol, max_iter)` | Yanlış konum |
| `secant(f, x0, x1, tol, max_iter)` | Kesen yöntemi |
| `newton(f, df, x0, tol, max_iter)` | Newton-Raphson |
| `newton_fd(f, x0, tol, max_iter)` | Sonlu fark türevli Newton |
| `brent(f, a, b, tol, max_iter)` | Brent yöntemi |

**Notlar:**
- `bisect`, `regula_falsi`, `brent` için f(a) ve f(b) zıt işaretli olmalı; değilse `NaN`
- Varsayılanlar: `tol=1e-12`, `max_iter=100`

## Lineer Cebir (`matrix.hpp`)

### Sınıf `Matrix`

| Üye | Açıklama |
|---|---|
| `Matrix(rows, cols, value=0.0)` | Sıfır/değer matrisi |
| `Matrix(rows, cols, {values...})` | İlkelleştirme listesiyle |
| `rows()`, `cols()` | Boyutlar |
| `operator()(i, j)` | Eleman erişimi |
| `transpose()` | Devrik |
| `operator*(vector)` | Matris-vektör çarpımı |

### Serbest fonksiyonlar

| Fonksiyon | Açıklama |
|---|---|
| `operator+`, `operator-` | Toplama, çıkarma |
| `operator*` | Matris çarpımı, skaler çarpım |
| `solve(A, b)` | Gauss eleme ile A·x = b çözümü |
| `det(A)` | Determinant |

**Notlar:**
- Boyut uyuşmazlığında `std::invalid_argument`
- Tekil matriste `solve` → `std::runtime_error`, `det` → 0

## Polinomlar (`polynomial.hpp`)

### Sınıf `Polynomial`

Katsayılar artan derecede saklanır: `{c₀, c₁, c₂, ...}`.

| Üye | Açıklama |
|---|---|
| `Polynomial({c0, c1, ...})` | Kurucu |
| `degree()` | Derece |
| `leading()` | Baş katsayı |
| `operator()(x)` | Horner ile değerlendirme |
| `coeff(k)` | k. katsayı |
| `derivative()` | Sembolik türev |
| `integral()` | Sembolik integral |
| `real_roots()` | Gerçel kökler |
| `roots()` | Kompleks kökler (Durand-Kerner) |

### Serbest fonksiyonlar

| Fonksiyon | Açıklama |
|---|---|
| `operator+`, `operator-`, `operator*` | Toplama, çıkarma, çarpma |
| `operator* (double, Poly)` | Skaler çarpım |
| `operator/` | Polinom bölümü (bölüm) |
| `operator%` | Polinom bölümü (kalan) |
| `polynomial_divmod(a, b, q, r)` | Bölüm + kalan |

**Notlar:**
- Derece ≤ 3 için kesin formüller (karekök / Cardano)
- Derece ≥ 4 için Durand-Kerner iterasyonu
- Sıfır polinomla bölme `std::domain_error` fırlatır

## İstatistik (`statistics.hpp`)

| Fonksiyon | Açıklama |
|---|---|
| `mean(x)` | Aritmetik ortalama |
| `variance(x)` | Örneklem varyansı (n-1) |
| `variance_pop(x)` | Anakütle varyansı (n) |
| `stddev(x)` | Örneklem standart sapması |
| `stddev_pop(x)` | Anakütle standart sapması |
| `median(x)` | Medyan |
| `min(x)`, `max(x)` | Min / max |
| `sum(x)`, `product(x)` | Toplam / çarpım |
| `geometric_mean(x)` | Geometrik ortalama |
| `harmonic_mean(x)` | Harmonik ortalama |
| `percentile(x, p)` | p. yüzdelik |
| `covariance(x, y)` | Kovaryans |
| `correlation(x, y)` | Pearson korelasyonu |
| `skewness(x)` | Çarpıklık |
| `kurtosis(x)` | Fazlalık basıklık |
| `range(x)` | Aralık (max - min) |

## ODE Çözücüler (`ode.hpp`)

`OdeFunction = std::function<double(double t, double y)>` — dy/dt = f(t, y).

| Fonksiyon | Açıklama |
|---|---|
| `ode_euler(f, y0, t0, t1, n=1000)` | Euler (1. derece) |
| `ode_heun(f, y0, t0, t1, n=1000)` | Heun (2. derece) |
| `ode_rk4(f, y0, t0, t1, n=1000)` | Klasik RK4 |
| `ode_rk45(f, y0, t0, t1, tol=1e-8)` | Adaptif Dormand-Prince |

Tümü `std::vector<std::pair<double, double>>` (t, y) dizisi döner. RK45, adım boyutunu `tol`'a göre otomatik ayarlar.

## Hipergeometrik Fonksiyonlar (`hypergeometric.hpp`)

| Fonksiyon | Açıklama |
|---|---|
| `hypergeometric_1f1(a, b, z)` | Kummer'in birleşik hipergeometrik fonksiyonu ₁F₁(a; b; z) |
| `hypergeometric_2f1(a, b, c, z)` | Gauss hipergeometrik fonksiyonu ₂F₁(a, b; c; z) |
| `hypergeometric_pfq(...)` | Genelleştirilmiş pFq serisi |
| `elliptic_k(m)` | Birinci tür tam eliptik integral K(m) |
| `elliptic_e(m)` | İkinci tür tam eliptik integral E(m) |

**Notlar:**
- ₁F₁ ve ₂F₁ Taylor serisi toplamıyla hesaplanır; c ≤ 0 ise `NaN`
- |z| < 1 için ₂F₁ yakınsar (mutlak yakınsaklık çemberi)
- `elliptic_k` ve `elliptic_e` AGM / kuvvet serisi ile hesaplanır; m ∈ [0, 1]
- `K(0) = E(0) = π/2`, `E(1) = 1`

## Otomatik Türev (`autodiff.hpp`)

Header-only. İleri mod dual sayılar ile otomatik türev.

### Sınıf `Dual`

| Üye | Açıklama |
|---|---|
| `Dual(value, derivative=0.0)` | Kurucu |
| `value()` | Değer |
| `deriv()` | Türev |

### Serbest fonksiyonlar

| Fonksiyon | Açıklama |
|---|---|
| `derivative(f, x)` | `f`'in x noktasında türevi (`std::pair<double,double>` döner) |
| `derivatives(f, x, order)` | x'te ilk `order` türev vektörü |
| `Dual operator+,-,*,/` | Dual aritmetik |
| `mathx::sin/cos/exp/log/sqrt/pow/tanh(Dual)` | Dual üzerinde standart fonksiyonlar |

**Örnek:**

```cpp
auto f = [](const mathx::Dual& x) { return x * x; };
auto [v, d] = mathx::derivative(f, 3.0);  // v = 9, d = 6
```
