cmake_minimum_required(VERSION 3.16)
project(cpp_math_tests LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

enable_testing()

function(mathx_add_test name)
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name} PRIVATE mathx)
    add_test(NAME ${name} COMMAND ${name})
endfunction()

mathx_add_test(test_lambert)
mathx_add_test(test_fourier)
mathx_add_test(test_special)
mathx_add_test(test_matrix)
mathx_add_test(test_bessel)
mathx_add_test(test_integration)
mathx_add_test(test_roots)
mathx_add_test(test_statistics)
mathx_add_test(test_polynomial)
mathx_add_test(test_ode)
