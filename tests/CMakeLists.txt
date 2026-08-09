cmake_minimum_required(VERSION 3.16)
project(cpp_math_tests LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(mathx_tests
    test_lambert.cpp
    test_fourier.cpp
    test_special.cpp
    test_matrix.cpp
)

target_link_libraries(mathx_tests PRIVATE mathx)

add_test(NAME mathx_tests COMMAND mathx_tests)
