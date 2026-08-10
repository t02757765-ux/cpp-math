cmake_minimum_required(VERSION 3.16)
project(cpp_math VERSION 0.3.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

option(MATHX_BUILD_SHARED "Build as a shared library (.dll/.so)" ON)
option(MATHX_BUILD_TESTS "Build unit tests" OFF)
option(MATHX_ENABLE_OPENMP "Enable OpenMP parallel acceleration" ON)
option(MATHX_ENABLE_NATIVE "Enable -march=native / /arch:AVX2 optimizations" ON)

find_package(OpenMP QUIET)

# Export macros for Windows (__declspec) / Unix (visibility)
if(MATHX_BUILD_SHARED)
    add_library(mathx SHARED
        src/lambert.cpp
        src/fourier.cpp
        src/special.cpp
        src/matrix.cpp
        src/bessel.cpp
        src/integration.cpp
        src/roots.cpp
        src/statistics.cpp
        src/polynomial.cpp
        src/ode.cpp
        src/hypergeometric.cpp
    )
else()
    add_library(mathx STATIC
        src/lambert.cpp
        src/fourier.cpp
        src/special.cpp
        src/matrix.cpp
        src/bessel.cpp
        src/integration.cpp
        src/roots.cpp
        src/statistics.cpp
        src/polynomial.cpp
        src/ode.cpp
        src/hypergeometric.cpp
    )
endif()

target_include_directories(mathx PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_compile_definitions(mathx PRIVATE MATHX_BUILDING_LIBRARY)

if(MSVC)
    target_compile_options(mathx PRIVATE /W4)
else()
    target_compile_options(mathx PRIVATE -Wall -Wextra)
endif()

if(MATHX_ENABLE_OPENMP AND OpenMP_CXX_FOUND)
    target_link_libraries(mathx PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(mathx PUBLIC MATHX_HAS_OPENMP)
endif()

# Native/SIMD optimization flags (Release only).
if(MATHX_ENABLE_NATIVE AND NOT MSVC)
    target_compile_options(mathx PRIVATE -O3 -march=native -funroll-loops)
endif()

# Install rules so the library can be consumed by other projects
include(GNUInstallDirs)
install(TARGETS mathx
    EXPORT mathxTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT mathxTargets
    FILE mathxTargets.cmake
    NAMESPACE mathx::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/mathx
)

if(MATHX_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
