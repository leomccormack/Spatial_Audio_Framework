
cmake_minimum_required(VERSION 3.15)

project(saf_head)

# Project options/defaults
option(BUILD_SHARED_LIBS             "Build shared instead of static libraries."  OFF)
option(SAF_BUILD_TESTS               "Build SAF unit tests."                      ON)
option(SAF_BUILD_EXAMPLES            "Build SAF examples."                        ON)
option(SAF_BUILD_EXTRAS              "Build SAF extras."                          OFF)
option(SAF_ENABLE_SOFA_READER_MODULE "Enable the SAF SOFA READER module"          OFF)
option(SAF_ENABLE_TRACKER_MODULE     "Enable the SAF TRACKER module"              OFF)
option(SAF_ENABLE_HADES_MODULE       "Enable the SAF HADES module"                OFF)
option(SAF_USE_INTEL_IPP             "Use Intel IPP for the FFT, resampler, etc." OFF)
option(SAF_USE_FFTW                  "Use FFTW3 for the FFT."                     OFF)
option(SAF_ENABLE_SIMD               "Enable the use of SSE3, AVX2, AVX512"       OFF)
option(SAF_ENABLE_NETCDF             "Enable netcdf for the sofa reader module"   OFF)
option(SAF_USE_FAST_MATH_FLAG        "Enable -ffast-math compiler flag"           ON)
if(NOT DEFINED SAF_PERFORMANCE_LIB AND NOT APPLE)
    message(FATAL_ERROR "You must define the SAF_PERFORMANCE_LIB environment variable")
endif()
if(NOT DEFINED SAF_PERFORMANCE_LIB AND APPLE)
    set(SAF_PERFORMANCE_LIB "SAF_USE_APPLE_ACCELERATE_LP64")
endif()
#message(STATUS "PERFORMANCE LIB is ${SAF_PERFORMANCE_LIB}")
if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

# Configure CMake
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
if(MSVC)
    set(CMAKE_C_STANDARD 99)
else()
    set(CMAKE_C_STANDARD 11)
endif()
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS true)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Platform specific flags and definitions
if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /MP")
elseif(APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++") 
    message(STATUS "Targeting both arm64 and x86_64.")
    set(CMAKE_OSX_ARCHITECTURES "arm64;x86_64")
    message(STATUS "Targeting MacOS 12.0.")
    set(CMAKE_OSX_DEPLOYMENT_TARGET "12.0")
elseif(UNIX AND NOT APPLE AND NOT ANDROID)
    add_definitions(
        #-std=c++11
        -DLINUX=1   
    )
endif()

# add subdirectory for SAF
add_subdirectory(framework)

# Optionally, also add the SAF examples and/or unit testing program 
set(example_prefix "saf_example_")
if(SAF_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()
if(SAF_BUILD_EXTRAS)
    add_subdirectory(extras)
endif()
if(SAF_BUILD_TESTS)
    add_subdirectory(test)
endif()
