# Başlangıç

Bu bölüm cpp-math'i kurmanızı ve ilk programınızı yazmanızı sağlar.

## Gereksinimler

- C++17 uyumlu derleyici (GCC ≥ 8, Clang ≥ 8, MSVC ≥ 2019)
- CMake ≥ 3.16 (opsiyonel, kütüphaneyi kurmak için)
- Platform: Windows (.dll), Linux (.so), macOS (.dylib)

## Derleme (CMake ile)

```sh
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .
```

Çıktılar:
- **Windows**: `mathx.dll` + `mathx.lib` (veya MinGW ile `libmathx.dll`)
- **Linux**: `libmathx.so`
- **macOS**: `libmathx.dylib`

### Yapılandırma seçenekleri

| Seçenek | Varsayılan | Açıklama |
|---|---|---|
| `MATHX_BUILD_SHARED` | `ON` | Paylaşımlı kütüphane (.dll/.so) derle |
| `MATHX_BUILD_TESTS` | `OFF` | Testleri derle |

Testleri çalıştırmak için:

```sh
cmake .. -DCMAKE_BUILD_TYPE=Release -DMATHX_BUILD_TESTS=ON
cmake --build .
ctest --output-on-failure
```

## Manuel Derleme (CMake olmadan)

Doğrudan `g++` ile tek başına derlemek için tüm kaynakları verin:

```sh
g++ -std=c++17 -Iinclude -shared -fPIC \
    src/*.cpp -o libmathx.so
```

Windows (MinGW) için:

```sh
g++ -std=c++17 -DMATHX_BUILDING_LIBRARY -Iinclude -shared \
    src/*.cpp -o mathx.dll -Wl,--out-implib,mathx.lib
```

## İlk Program

`main.cpp`:

```cpp
#include <mathx/mathx.hpp>
#include <iostream>

int main() {
    // Sabitler
    std::cout << "pi = " << mathx::constants::pi << "\n";

    // Lambert W
    std::cout << "W(1) = " << mathx::lambert_w(1.0) << "\n";

    // Özel fonksiyonlar
    std::cout << "gamma(5) = " << mathx::gamma(5.0) << "\n";   // 24
    std::cout << "zeta(2) = " << mathx::zeta(2.0) << "\n";     // pi^2/6

    // Matris
    mathx::Matrix A(2, 2, {4.0, 2.0, 1.0, 3.0});
    auto x = mathx::solve(A, {6.0, 4.0});                      // {1, 1}

    // Sayısal integral
    double a = mathx::integrate_gauss_legendre(
        [](double t) { return t * t; }, 0.0, 1.0);             // 1/3

    std::cout << "integral = " << a << "\n";
    return 0;
}
```

Derleyin ve çalıştırın:

```sh
g++ -std=c++17 -Iinclude main.cpp -Lbuild -lmathx -o main
./main            # Linux/macOS
set PATH=%PATH%;<dll dizini> && main   # Windows
```

## Kurulum (opsiyonel)

```sh
cmake --install . --prefix /usr/local
```

Kurulum sonrası `#include <mathx/mathx.hpp>` ile kullanabilirsiniz. Diğer projelerde `find_package` ile tüketme yöntemi için [C++ Entegrasyonu](integration.md) bölümüne bakın.
