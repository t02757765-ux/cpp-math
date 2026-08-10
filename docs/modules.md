# Modüller

Bu bölüm her modülün arkasındaki matematiği ve algoritmaları açıklar.

## Lambert W Fonksiyonu

Lambert W, `W(x)·e^(W(x)) = x` denkleminin çözümüdür.

**Dallar:**
- **W₀** (principal): `x ≥ -1/e` için tanımlı, `W₀(x) ≥ -1`
- **W₋₁**: `-1/e ≤ x < 0` için tanımlı, `W₋₁(x) ≤ -1`

**Algoritma:**
- Küçük |x|: seri açılımıyla başlangıç tahmini
- Büyük x: `log(x) - log(log(x))` asimptotik tahmini
- Ardından Newton yinelemesi (gerçel) veya Halley yinelemesi (kompleks)
- 10⁻¹⁴ göreli hassasiyet

**Uygulamalar:** Diferansiyel denklemler, kombinatorik, fizik (Bose-Einstein yoğuşması), bilgisayar bilimi (veri yapıları analizi).

## Fourier Dönüşümleri

Ayrık Fourier dönüşümü:

```
X[k] = Σ_{t=0}^{n-1} x[t] · e^(-2πi·k·t/n)
```

**DFT:** Doğrudan O(n²) hesaplama. Küçük girdiler için basit ve doğru.

**FFT:** İteratif Cooley-Tukey algoritması, O(n log n).
1. Bit-ters çevirme permütasyonu
2. Butterfly aşamaları (len = 2, 4, 8, ..., n)
3. Girdi boyutu 2'nin kuvveti olmalı

**Ters dönüşümler:** `idft` ve `ifft` standart 1/n normalizasyonu kullanır.

**Uygulamalar:** Sinyal işleme, spektral analiz, polinom çarpımı, görüntü işleme.

## Özel Fonksiyonlar

### Gamma fonksiyonu

```
Γ(x) = ∫₀^∞ t^(x-1) e^(-t) dt
```

Lanczos yaklaşımı (g=7) kullanılır. Pozitif olmayan tam sayılarda kutuplar `NaN` döndürür. Yansıma formülü negatif argümanları işler.

### Riemann zeta

```
ζ(s) = Σ_{n=1}^∞ 1/n^s
```

- `s ≤ 0`: fonksiyonel denklem kullanılır
- `s > 1`: Euler-Maclaurin toplamı
- `s = 1`: kutup → `+∞`

### Digamma

```
ψ(x) = d/dx ln Γ(x)
```

Küçük x için asimptotik seri kullanılmadan önce kaydırma (recurrence) uygulanır.

## Bessel Fonksiyonları

Bessel diferansiyel denkleminin çözümleri:

```
x²y'' + xy' + (x² - ν²)y = 0     → J_ν, Y_ν  (birinci, ikinci tür)
x²y'' + xy' - (x² + ν²)y = 0     → I_ν, K_ν  (modifiye)
```

**Algoritmalar:**
- **Küçük |x|:** Kuvvet serisi açılımı
- **Büyük |x|:** Asimptotik açılımlar (Debye)
- **Yüksek mertebe:** Geriye (J) / ileriye (Y, I, K) yineleme ilişkileri
- **Y₀, Y₁, K₀, K₁:** Harmonik sayılar içeren seriler

**Uygulamalar:** Dalga yayılımı, ısı denklemi (silindirik), kuantum mekaniği, akustik.

### Airy fonksiyonları

`y'' - xy = 0` denkleminin çözümleri Ai(x) ve Bi(x). Seri açılımı ile hesaplanır.

## Sayısal İntegral

| Yöntem | Sıra | Açıklama |
|---|---|---|
| Trapez | O(h²) | Doğrusal interpolasyon |
| Simpson | O(h⁴) | Parabolik interpolasyon |
| Adaptif Simpson | O(h⁴) | Hata tahminine göre aralıkları böler |
| Gauss-Legendre | O(h^(2n)) | Düğümler ve ağırlıklar seçilir |
| Romberg | Yüksek | Richardson ekstrapolasyonu |

**Öneriler:**
- Düzgün fonksiyonlar için `integrate_gauss_legendre(..., 32)`
- Genel amaç için `integrate_simpson`
- Tepe noktaları olan fonksiyonlar için `integrate_romberg`

## Kök Bulma

| Yöntem | Yakınsama | Sınırlama |
|---|---|---|
| Bisection | Doğrusal | Aralık gerektirir |
| Regula falsi | Doğrusal | Tek taraflı yakınsayabilir |
| Secant | Altın oran (~1.618) | Başlangıç noktaları gerektirir |
| Newton | İkinci derece | Türev gerektirir |
| Brent | Süperdoğrusal | Aralık gerektirir, en sağlam |

**Brent yöntemi**, bisection, secant ve ters kuadratik interpolasyonu birleştirir; bu yüzden genel amaç için önerilir.

## Lineer Cebir

`Matrix` sınıfı satır-majör yoğun matris saklar. `solve` Gauss eleme + kısmi pivotlama kullanır:

1. Kısmi pivotlama (sayısal kararlılık)
2. İleri eleme (augmented matris)
3. Geri ikame

**Kararlılık:** Pivotlama olmadan birçok sistem yanlış sonuç verir. Bu uygulama her adımda en büyük mutlak değerli satırı seçer.

## Polinomlar

- Katsayılar artan derecede: `{c₀, c₁, c₂, ...}`
- Derece 1-3: kesin formüller (doğrusal, kuadratik, Cardano)
- Derece ≥ 4: Durand-Kerner eşzamanlı iterasyonu (kompleks kökler)
- `real_roots()` yalnızca imajiner kısmı ihmal edilebilir kökleri döndürür
- Horner şeması O(n) değerlendirme sağlar

## İstatistik

Tüm fonksiyonlar `std::vector<double>` kabul eder.

- **Örneklem vs anakütle:** `variance` n-1, `variance_pop` n böler
- **Yüzdelik:** En yakın sıralar arasında doğrusal interpolasyon (NumPy'nin `percentile` davranışı)
- **Korelasyon:** Pearson; sıfır standart sapmada 0 döner
- Boş girdi durumlarında 0 döner (güvenli varsayılan)

## ODE Çözücüler

Tek birinci mertebe ODE: `dy/dt = f(t, y)`, `y(t₀) = y₀`.

| Yöntem | Sıra | Açıklama |
|---|---|---|
| Euler | 1 | `y + h·f` |
| Heun | 2 | Orta nokta düzeltmeli |
| RK4 | 4 | 4 aşamalı klasik |
| Dormand-Prince (RK45) | 5 | 7 aşamalı, adaptif adım |

**RK45 adım kontrolü:** 4. ve 5. derece çözümlerin farkı yerel hatayı tahmin eder; adım boyutu `tol^(1/5)` oranında büyütülür/küçültülür (güvenlik faktörü 0.9, sınırlar [0.2, 5]).

Yüksek mertebeden sistemler için ODE'yi vektörel olarak kendiniz yazabilir veya birden fazla skaler denklemi ayrı çözebilirsiniz.
