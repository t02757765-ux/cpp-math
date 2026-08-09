#ifndef MATHX_EXPORT_HPP
#define MATHX_EXPORT_HPP

#if defined(_WIN32) || defined(__CYGWIN__)
    #define MATHX_EXPORT __declspec(dllexport)
    #define MATHX_IMPORT __declspec(dllimport)
#elif defined(__GNUC__) && __GNUC__ >= 4
    #define MATHX_EXPORT __attribute__((visibility("default")))
    #define MATHX_IMPORT __attribute__((visibility("default")))
#else
    #define MATHX_EXPORT
    #define MATHX_IMPORT
#endif

#ifdef MATHX_BUILDING_LIBRARY
    #define MATHX_API MATHX_EXPORT
#else
    #define MATHX_API MATHX_IMPORT
#endif

#endif
