# C++ Entegrasyonu

Bu bölüm cpp-math'i kendi C++ projelerinizde nasıl kullanacağınızı açıklar.

## Yöntem 1: CMake ile kurulum (önerilen)

Önce kütüphaneyi kurun:

```sh
cmake --install build --prefix /path/to/prefix
```

Sonra kendi `CMakeLists.txt` dosyanızda:

```cmake
cmake_minimum_required(VERSION 3.16)
project(my_app LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

# prefix'i CMake'e göster
list(APPEND CMAKE_PREFIX_PATH /path/to/prefix)

find_package(mathx CONFIG REQUIRED)

add_executable(my_app main.cpp)
target_link_libraries(my_app PRIVATE mathx::mathx)
```

`find_package` kurulan `mathxTargets.cmake` dosyasını bulur ve `mathx::mathx` hedefini içe aktarır.

## Yöntem 2: add_subdirectory ile

cpp-math'i projenize alt dizin olarak eklerseniz:

```cmake
add_subdirectory(path/to/cpp-math)

add_executable(my_app main.cpp)
target_link_libraries(my_app PRIVATE mathx)
```

## Yöntem 3: Doğrudan derleme

Kütüphane dosyasını ve başlıkları elle bağlayın:

```sh
# Linux/macOS
g++ -std=c++17 main.cpp -I/cpp-math/include -L/cpp-math/build -lmathx -o app

# Windows (MinGW)
g++ -std=c++17 main.cpp -Icpp-math\include -Lcpp-math\build -lmathx -o app.exe

# Windows (MSVC)
cl /std:c++17 main.cpp /Icpp-math\include mathx.lib
```

## Çalışma zamanında DLL bulma

**Windows:** `mathx.dll` dosyası, uygulamanın çalıştığı dizinde veya PATH'te olmalı:

```bat
set PATH=%PATH%;C:\cpp-math\build
my_app.exe
```

**Linux:** `LD_LIBRARY_PATH` ayarlayın veya `ldconfig` ile kaydedin:

```sh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cpp-math/build
./app
```

**Windows/DLL bağımlılıkları:** MinGW derlemesi, çalışma zamanında `libgcc_s_seh-1.dll` ve `libstdc++-6.dll` (MSYS2/ucrt64 bin dizininde) gerektirir. Uygulamanın yanına kopyalayın veya PATH'e ekleyin.

## Başlık organizasyonu

Tüm kütüphaneyi tek başlıkla dahil edebilirsiniz:

```cpp
#include <mathx/mathx.hpp>
```

veya yalnızca ihtiyacınız olan modülü:

```cpp
#include <mathx/lambert.hpp>
#include <mathx/fft.hpp>
#include <mathx/matrix.hpp>
```

## Başka dillerden çağırma

C++ kütüphanesi doğrudan başka dillerden çağrılamaz, ancak:

- **C**: `extern "C"` sarmalayıcı yazarak
- **Python**: `pybind11` veya ctypes ile
- **Rust**: FFI ile

Sarmalayıcı örneği (C):

```cpp
// wrapper.cpp
#include <mathx/mathx.hpp>
extern "C" {
    double mathx_lambert_w(double x) { return mathx::lambert_w(x); }
    double mathx_gamma(double x) { return mathx::gamma(x); }
}
```
