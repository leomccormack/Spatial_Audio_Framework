cmake_minimum_required(VERSION 3.15)
#include(CheckIncludeFiles) # For some reason this causes an error on windows...

# Check required variables are defined
if(NOT DEFINED SAF_PERFORMANCE_LIB AND NOT APPLE)
    message(FATAL_ERROR "You must define the SAF_PERFORMANCE_LIB environment variable")
endif()
if(NOT DEFINED SAF_PERFORMANCE_LIB AND APPLE)
    set(SAF_PERFORMANCE_LIB "SAF_USE_APPLE_ACCELERATE")
endif()

# Pull version number from the framework's main header file:
file(READ "include/saf.h" ver)
string(REGEX MATCH "SAF_VERSION_MAJOR ([0-9]*)" _ ${ver})
set(ver_major ${CMAKE_MATCH_1})
string(REGEX MATCH "SAF_VERSION_MINOR ([0-9]*)" _ ${ver})
set(ver_minor ${CMAKE_MATCH_1})
string(REGEX MATCH "SAF_VERSION_PATCH ([0-9]*)" _ ${ver})
set(ver_patch ${CMAKE_MATCH_1})
set(SAF_VERSION "${ver_major}.${ver_minor}.${ver_patch}")

# Display License
if(SAF_ENABLE_TRACKER_MODULE)
    set(SAF_LICENSE "GNU GPLv2")
else()
    set(SAF_LICENSE "ISC")
endif()

# SAF Banner:
message("
     _____     _____     _____
   (  ____ ) (  ___  ) (  ____ )
   | ( __ |/ | ( _ ) | | ( _  |/
   (_____  ) |  ___  | |  __)
   /|____) | | |   | | | |
   (_______) |_|   |_| |_|   (Version: ${SAF_VERSION}, License: ${SAF_LICENSE})
")

# Configure SAF
message(STATUS "Configuring the Spatial_Audio_Framework (SAF):")
project(saf VERSION ${SAF_VERSION} LANGUAGES C)
add_library(${PROJECT_NAME}) #STATIC
set_target_properties(${PROJECT_NAME}
PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION} 
)

# Display system and configuration info
message(STATUS "${CMAKE_SYSTEM_INFO_FILE}, name: ${CMAKE_SYSTEM_NAME}, arch: ${CMAKE_SYSTEM_PROCESSOR}, sys: ${CMAKE_SYSTEM}")
if(BUILD_SHARED_LIBS)
    message(STATUS "Configuring ${PROJECT_NAME} as a SHARED library...")
else()
    message(STATUS "Configuring ${PROJECT_NAME} as a STATIC library...")
endif()

# Framework include paths
target_include_directories(${PROJECT_NAME} 
PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
include(GNUInstallDirs)

# Framework source files
add_subdirectory(modules)
add_subdirectory(resources)


############################################################################
# Performance library dependencies
if( ${SAF_PERFORMANCE_LIB} MATCHES "SAF_USE_INTEL_MKL_LP64" OR ${SAF_PERFORMANCE_LIB} MATCHES "SAF_USE_INTEL_MKL_ILP64")
   
    # Default header path and libs, (which are platform dependent):
    if(NOT DEFINED INTEL_MKL_HEADER_PATH)
        if(MSVC)
            set(INTEL_MKL_HEADER_PATH "C:/Program Files (x86)/Intel/oneAPI/mkl/latest/include") 
        elseif(MSYS OR MINGW) 
            set(INTEL_MKL_HEADER_PATH "/c/Program Files (x86)/Intel/oneAPI/mkl/latest/include") 
        elseif(APPLE)
            set(INTEL_MKL_HEADER_PATH "/opt/intel/oneapi/mkl/latest/include")
        elseif(UNIX AND NOT APPLE)
            if(EXISTS "/opt/intel/oneapi/mkl/latest/include")
                set(INTEL_MKL_HEADER_PATH "/opt/intel/oneapi/mkl/latest/include")
            else()
                message(SEND_ERROR "Intel mkl.h default path not found, please specify!")
            endif()
        endif()
    endif()
    if ( ${SAF_PERFORMANCE_LIB} MATCHES "SAF_USE_INTEL_MKL_LP64" )
        target_compile_definitions(${PROJECT_NAME} PUBLIC SAF_USE_INTEL_MKL_LP64=1)

        if (NOT DEFINED INTEL_MKL_LIB)
            if(MSVC)
                set(INTEL_MKL_LIB "${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/Win64/lib/saf_mkl_custom_lp64.lib")
            elseif(MSYS OR MINGW)
                set(INTEL_MKL_LIB "${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/Win64/lib/saf_mkl_custom_lp64.lib")
            elseif(APPLE)
                set(INTEL_MKL_LIB "/usr/local/lib/libsaf_mkl_custom_lp64.dylib")
            elseif(UNIX AND NOT APPLE)
                set(INTEL_MKL_LIB "/usr/local/lib/libsaf_mkl_custom_lp64.so")
            endif()
        endif()

    else()
        target_compile_definitions(${PROJECT_NAME} PUBLIC SAF_USE_INTEL_MKL_ILP64=1)

        if (NOT DEFINED INTEL_MKL_LIB)
            if(MSVC)
                set(INTEL_MKL_LIB "${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/Win64/lib/saf_mkl_custom_ilp64.lib")
            elseif(MSYS OR MINGW)
                set(INTEL_MKL_LIB "${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/Win64/lib/saf_mkl_custom_ilp64.lib")
            elseif(APPLE)
                set(INTEL_MKL_LIB "/usr/local/lib/libsaf_mkl_custom_ilp64.dylib")
            elseif(UNIX AND NOT APPLE)
                set(INTEL_MKL_LIB "/usr/local/lib/libsaf_mkl_custom_ilp64.so")
            endif()
        endif()
    endif()

    if (NOT EXISTS "${INTEL_MKL_HEADER_PATH}/mkl.h")
        message(SEND_ERROR "Intel mkl.h not found")
    endif()
    message(STATUS "Using Intel MKL. Linking against: ${INTEL_MKL_LIB}")

    # Add target header paths and libs
    target_include_directories(${PROJECT_NAME} PUBLIC ${INTEL_MKL_HEADER_PATH} )
    target_link_libraries(${PROJECT_NAME} PUBLIC ${INTEL_MKL_LIB} )

elseif( ${SAF_PERFORMANCE_LIB} MATCHES "SAF_USE_OPEN_BLAS_AND_LAPACKE")
    target_compile_definitions(${PROJECT_NAME} PUBLIC SAF_USE_OPEN_BLAS_AND_LAPACKE=1)

    # find and link libraries
    if (NOT DEFINED OPENBLAS_LIBRARY)
        find_library(OPENBLAS_LIBRARY openblas HINTS /usr/lib/x86_64-linux-gnu /usr/lib/arm-linux-gnueabihf)
    endif()
    if (NOT DEFINED LAPACKE_LIBRARY)
        find_library(LAPACKE_LIBRARY lapacke HINTS /usr/lib/x86_64-linux-gnu /usr/lib/arm-linux-gnueabihf)
    endif()
    target_link_libraries(${PROJECT_NAME} PUBLIC ${OPENBLAS_LIBRARY} ${LAPACKE_LIBRARY})

    # Add header search paths if they are defined
    if(DEFINED OPENBLAS_HEADER_PATH)
        target_include_directories(${PROJECT_NAME} PUBLIC ${OPENBLAS_HEADER_PATH})
    endif()
    if(DEFINED LAPACKE_HEADER_PATH)
        target_include_directories(${PROJECT_NAME} PUBLIC ${LAPACKE_HEADER_PATH})
    endif()

    # Disable this particular warning... We know, but it's OK.
    if(UNIX)
        target_compile_options(${PROJECT_NAME} PRIVATE -Wno-incompatible-pointer-types)
    endif()
    
    message(STATUS "Using OpenBLAS. Linking against: ${OPENBLAS_LIBRARY}; ${LAPACKE_LIBRARY}")

elseif( ${SAF_PERFORMANCE_LIB} MATCHES "SAF_USE_ATLAS")
    message(STATUS "Using ATLAS...")
    message(SEND_ERROR "Not yet supported! Please contribute if you use this library.")

elseif(APPLE AND ${SAF_PERFORMANCE_LIB} MATCHES "SAF_USE_APPLE_ACCELERATE")
    message(STATUS "Using Apple Accelerate. Linking against the Accelerate framework...")
    target_compile_definitions(${PROJECT_NAME} PUBLIC SAF_USE_APPLE_ACCELERATE=1) 
    target_link_libraries(${PROJECT_NAME} PUBLIC "-framework Accelerate") 

else()
    message(SEND_ERROR "The specified SAF_PERFORMANCE_LIB is not supported")

endif()

############################################################################
# Enable/Disable Intel IPP
if(SAF_USE_INTEL_IPP)

    # Default header path and libs, (which are platform dependent):
    if(NOT DEFINED INTEL_IPP_HEADER_PATH)
        if(MSVC)
            set(INTEL_IPP_HEADER_PATH "C:/Program Files (x86)/Intel/oneAPI/ipp/latest/include") 
        elseif(MSYS OR MINGW) 
            set(INTEL_IPP_HEADER_PATH "/c/Program Files (x86)/Intel/oneAPI/ipp/latest/include") 
        elseif(APPLE)
            set(INTEL_IPP_HEADER_PATH "/opt/intel/oneapi/ipp/latest/include")
        elseif(UNIX AND NOT APPLE)
            if(EXISTS "/opt/intel/oneapi/ipp/latest/include")
                set(INTEL_IPP_HEADER_PATH "/opt/intel/oneapi/ipp/latest/include")
            else()
                message(SEND_ERROR "Intel ipp.h default path not found, please specify!")
            endif()
        endif()
    endif()

    if (NOT DEFINED INTEL_IPP_LIB)
        if(MSVC)
            set(INTEL_IPP_LIB "${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/Win64/lib/saf_ipp_custom.lib")
        elseif(MSYS OR MINGW)
            set(INTEL_IPP_LIB "${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/Win64/lib/saf_ipp_custom.lib")
        elseif(APPLE)
            set(INTEL_IPP_LIB "/usr/local/lib/libsaf_ipp_custom.dylib")
        elseif(UNIX AND NOT APPLE)
            set(INTEL_IPP_LIB "/usr/local/lib/libsaf_ipp_custom.so")
        endif()
    endif()


    if (NOT EXISTS "${INTEL_IPP_HEADER_PATH}/ipp.h")
        message(SEND_ERROR "Intel ipp.h not found")
    endif()
    message(STATUS "Using Intel IPP. Linking against: ${INTEL_IPP_LIB}")

    # Add target header paths and libs
    target_include_directories(${PROJECT_NAME} PUBLIC ${INTEL_IPP_HEADER_PATH} )
    target_compile_definitions(${PROJECT_NAME} PUBLIC SAF_USE_INTEL_IPP=1)
    target_link_libraries(${PROJECT_NAME} PUBLIC ${INTEL_IPP_LIB} )
endif()

############################################################################
# Enable/Disable FFTW
if(SAF_USE_FFTW)
    # Indicate to saf that the saf_sofa_reader module should be enabled
    target_compile_definitions(${PROJECT_NAME} PUBLIC SAF_USE_FFTW=${SAF_USE_FFTW})
        
    # Header path and libs are platform dependent
    if(MSVC)
        message(FATAL_ERROR "Incomplete list")
    elseif(MSYS OR MINGW)
        message(FATAL_ERROR "Incomplete list")
    elseif(APPLE) 
        message(FATAL_ERROR "Incomplete list")
    elseif(UNIX AND NOT APPLE) 
        find_library(FFTW_LIBRARY fftw3f HINTS /usr/local/lib)
        if (NOT FFTW_LIBRARY)
            message(FATAL_ERROR "FFTW_LIBRARY not found")
        endif()
        target_link_libraries(${PROJECT_NAME} PUBLIC ${FFTW_LIBRARY})
        message(STATUS "Linking FFTW: ${FFTW_LIBRARY}")
    endif() 
endif()

############################################################################
# Enable SIMD intrinsics
if(SAF_ENABLE_SIMD)
    # Note that you need to pass the appropriate compiler flags to actually
    # enable certain intrinsics, for example:
    # $ cmake -S . -B build -DCMAKE_CXX_FLAGS="-msse3 -mavx2" -DCMAKE_C_FLAGS="-msse3 -mavx2"
    message(STATUS "SIMD intrinsics support is enabled.")
    target_compile_definitions(${PROJECT_NAME} PUBLIC SAF_ENABLE_SIMD=1)
endif()

############################################################################
# Sofa reader module dependencies
if(SAF_ENABLE_SOFA_READER_MODULE)
    message(STATUS "saf_sofa_reader module enabled.")

    # Indicate to saf that the saf_sofa_reader module should be enabled
    target_compile_definitions(${PROJECT_NAME} PUBLIC SAF_ENABLE_SOFA_READER_MODULE=${SAF_ENABLE_SOFA_READER_MODULE})

    if(SAF_ENABLE_NETCDF)
        message(STATUS "netcdf enabled for saf_sofa_reader module. Linking against netcdf...")
        target_compile_definitions(${PROJECT_NAME} PUBLIC SAF_ENABLE_NETCDF=${SAF_ENABLE_NETCDF})
        
        # Header path and libs are platform dependent
        if(MSVC)
            # Link NetCDF dependencies (HDF5 and its dependencies: zlib,szip):
            if(SAF_HDF5_INCLUDE_DIR AND SAF_HDF5_LIB_DIR)
                target_include_directories(${PROJECT_NAME} PUBLIC ${SAF_HDF5_INCLUDE_DIR})
                set(SAF_HDF5_LIBS
                    ${SAF_HDF5_LIB_DIR}/libhdf5.lib
                    ${SAF_HDF5_LIB_DIR}/libzlib.lib
                    ${SAF_HDF5_LIB_DIR}/libszip.lib
                    ${SAF_HDF5_LIB_DIR}/libhdf5_hl.lib
                )
            else()
                target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/Win64/include)
                set(SAF_HDF5_LIBS
                    ${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/Win64/lib/libzlib.lib
                ) 
                set(SAF_HDF5_LIBS ${SAF_HDF5_LIBS}
                    ${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/Win64/lib/libhdf5.lib
                    ${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/Win64/lib/libszip.lib
                    ${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/Win64/lib/libhdf5_hl.lib
                )
            endif()
            target_link_libraries(${PROJECT_NAME} PRIVATE ${SAF_HDF5_LIBS})
            message(STATUS "Linking HDF5: ${SAF_HDF5_LIBS}")

            # Link NetCDF
            if(SAF_NETCDF_INCLUDE_DIR AND SAF_NETCDF_LIB_DIR)
                target_include_directories(${PROJECT_NAME} PUBLIC ${SAF_NETCDF_INCLUDE_DIR})
                set(SAF_NETCDF_LIBS ${SAF_NETCDF_LIB_DIR}/netcdf.lib)
            else()
                target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/Win64/include)
                set(SAF_NETCDF_LIBS ${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/Win64/lib/netcdf.lib)
            endif()
            target_link_libraries(${PROJECT_NAME} PRIVATE ${SAF_NETCDF_LIBS})
            message(STATUS "Linking NETCDF: ${SAF_NETCDF_LIBS}")
            
        elseif(MSYS OR MINGW)
            message(SEND_ERROR "Currently unsupported on MSYS")

        elseif(APPLE) 
            if (SAF_NETCDF_INCLUDE_DIR AND SAF_NETCDF_LIBRARIES)
                 message(STATUS "Using NetCDF include dir ${SAF_NETCDF_INCLUDE_DIR} and libraries ${SAF_NETCDF_LIBRARIES}")
                 target_include_directories(${PROJECT_NAME} PUBLIC ${SAF_NETCDF_INCLUDE_DIR})
                 target_link_libraries(${PROJECT_NAME} PRIVATE ${SAF_NETCDF_LIBRARIES})
             else()
                 target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/MacOSX/include)
                 target_link_libraries(${PROJECT_NAME}
                 PRIVATE
                     ${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/MacOSX/lib/libz.a
                 )
                 target_link_libraries(${PROJECT_NAME}
                 PRIVATE
                     ${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/MacOSX/lib/libnetcdf.a
                     ${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/MacOSX/lib/libhdf5.a
                     ${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/MacOSX/lib/libhdf5_hl.a
                 )
             endif()

        elseif(UNIX AND NOT APPLE)
            find_library(NETCDF_LIBRARY netcdf HINTS HINTS /usr/lib/x86_64-linux-gnu /usr/lib/arm-linux-gnueabihf)
            if (NOT NETCDF_LIBRARY)
                message(FATAL_ERROR "NETCDF_LIBRARY not found")
            endif()
            target_link_libraries(${PROJECT_NAME} PUBLIC ${NETCDF_LIBRARY})
            message(STATUS "Linking NETCDF: ${NETCDF_LIBRARY}")

        endif() 
    endif() 
else()
    message(STATUS "saf_sofa_reader module disabled.")
endif()


############################################################################
# Enable/Disable tracker module
if(SAF_ENABLE_TRACKER_MODULE)
    message(STATUS "saf_tracker module enabled. ")

    # Indicate to saf that the saf_tracker module should be enabled
    target_compile_definitions(${PROJECT_NAME} PUBLIC SAF_ENABLE_TRACKER_MODULE=${SAF_ENABLE_TRACKER_MODULE})

else()
    message(STATUS "saf_tracker module disabled.")
endif()

############################################################################
# Enable/Disable HADES module
if(SAF_ENABLE_HADES_MODULE)
    message(STATUS "saf_hades module enabled. ")

    # Indicate to saf that the saf_hades module should be enabled
    target_compile_definitions(${PROJECT_NAME} PUBLIC SAF_ENABLE_HADES_MODULE=${SAF_ENABLE_HADES_MODULE})

else()
    message(STATUS "saf_hades module disabled.")
endif()


############################################################################
# Extra compiler flags
if(UNIX)
    target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra)
endif()
if(SAF_ENABLE_FAST_MATH_FLAG)
    if((CMAKE_CXX_COMPILER_ID MATCHES "Clang") OR (CMAKE_CXX_COMPILER_ID MATCHES "GNU"))
        add_compile_options(-ffast-math)
    endif()
endif()


############################################################################
# set RPATHs
set(CMAKE_INSTALL_RPATH $ORIGIN)
