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
     _____     _____     ______
   (  ____ ) (  ___  ) (  ____  )
   | (    |/ | (   ) | | (    |/
   | (_____  | (___) | | (__
   (_____  ) |  ___  | |  __)
         ) | | (   ) | | (
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
if( ${SAF_PERFORMANCE_LIB} MATCHES "SAF_USE_INTEL_MKL" )
    target_compile_definitions(${PROJECT_NAME} PUBLIC SAF_USE_INTEL_MKL=1)

    # Default header path and libs, (which are platform dependent):
    if(NOT DEFINED INTEL_MKL_HEADER_PATH)
        if(MSVC)
            set(INTEL_MKL_HEADER_PATH "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/include") 
        elseif(MSYS OR MINGW) 
            set(INTEL_MKL_HEADER_PATH "/c/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/include") 
        elseif(APPLE)
            set(INTEL_MKL_HEADER_PATH "/opt/intel/compilers_and_libraries/mac/mkl/include")
        elseif(UNIX AND NOT APPLE)
            if(EXISTS "~/intel/compilers_and_libraries/linux/mkl/include")
                set(INTEL_MKL_HEADER_PATH "~/intel/compilers_and_libraries/linux/mkl/include")
            elseif(EXISTS "/opt/intel/compilers_and_libraries/linux/mkl/include")
                set(INTEL_MKL_HEADER_PATH "/opt/intel/compilers_and_libraries/linux/mkl/include")
            else()
                message(SEND_ERROR "Intel mkl.h default path not found, please specify!")
            endif()
        endif()
    endif()
    if (NOT DEFINED INTEL_MKL_LIB)
        if(MSVC) 
            set(INTEL_MKL_LIB "${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/Win64/lib/saf_mkl_custom.lib")  
        elseif(MSYS OR MINGW) 
            set(INTEL_MKL_LIB "${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/Win64/lib/saf_mkl_custom.lib") 
        elseif(APPLE)
            set(INTEL_MKL_LIB "/usr/local/lib/libsaf_mkl_custom.dylib") 
        elseif(UNIX AND NOT APPLE)
            set(INTEL_MKL_LIB "/usr/lib/libsaf_mkl_custom.so") 
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
    message(STATUS "Using openBLAS...")
    target_compile_definitions(${PROJECT_NAME} PUBLIC SAF_USE_OPEN_BLAS_AND_LAPACKE=1)

    # find and link libraries
    if (NOT DEFINED OPENBLAS_LIBRARY)
        find_library(OPENBLAS_LIBRARY openblas)
    endif()
    if (NOT DEFINED LAPACKE_LIBRARY)
        find_library(LAPACKE_LIBRARY lapacke)
    endif()
    target_link_libraries(${PROJECT_NAME} PUBLIC ${OPENBLAS_LIBRARY} ${LAPACKE_LIBRARY})

    # Add header search paths if they are defined
    if(DEFINED OPENBLAS_HEADER_PATH)
        target_include_directories(${PROJECT_NAME} PUBLIC ${OPENBLAS_HEADER_PATH})
    endif()
    if(DEFINED LAPACKE_HEADER_PATH)
        target_include_directories(${PROJECT_NAME} PUBLIC ${LAPACKE_HEADER_PATH})
    endif()

    # Disable this particular warning...
    if(UNIX)
        target_compile_options(${PROJECT_NAME} PRIVATE -Wno-incompatible-pointer-types)
    endif()

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
# Sofa reader module dependencies
if(SAF_ENABLE_SOFA_READER_MODULE)
    message(STATUS "saf_sofa_reader module enabled. Linking against netcdf...")

    # Indicate to saf that the saf_sofa_reader module should be enabled
    target_compile_definitions(${PROJECT_NAME} PUBLIC SAF_ENABLE_SOFA_READER_MODULE=${SAF_ENABLE_SOFA_READER_MODULE})
        
    # Header path and libs are platform dependent
    if(MSVC) 
        target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/Win64/include)
        target_link_libraries(${PROJECT_NAME} 
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/Win64/lib/netcdf.lib
            ${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/Win64/lib/libhdf5.lib
            ${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/Win64/lib/libzlib.lib
            ${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/Win64/lib/libszip.lib
            ${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/Win64/lib/libhdf5_hl.lib
        ) 
    elseif(MSYS OR MINGW)
        message(SEND_ERROR "Unsupported on MSYS because I can't get the damn thing to link no matter what I do..") 
        target_include_directories(${PROJECT_NAME} PUBLIC 
             /c/msys64/mingw64/include)
        #set (CMAKE_SHARED_LINKER_FLAGS "-lnetcdf")
        find_library(NETCDF_LIBRARY netcdf HINTS "/c/msys64/mingw64/lib/libnetcdf.dll.a")
        if (NOT NETCDF_LIBRARY)
            message(FATAL_ERROR "NETCDF_LIBRARY not found")
        endif()
        target_link_libraries(${PROJECT_NAME} INTERFACE ${NETCDF_LIBRARY})
         
    elseif(APPLE)
        target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/MacOSX/include)
        target_link_libraries(${PROJECT_NAME} 
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/MacOSX/lib/libnetcdf.a
            ${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/MacOSX/lib/libhdf5.a
            ${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/MacOSX/lib/libz.a 
            ${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/MacOSX/lib/libhdf5_hl.a
        )  
    elseif(UNIX AND NOT APPLE) 
        find_library(NETCDF_LIBRARY netcdf HINTS /usr/lib/x86_64-linux-gnu)
        if (NOT NETCDF_LIBRARY)
            message(FATAL_ERROR "NETCDF_LIBRARY not found")
        endif()
        target_link_libraries(${PROJECT_NAME} PUBLIC ${NETCDF_LIBRARY})
    endif() 

else()
    message(STATUS "saf_sofa_reader module disabled.")
endif()


############################################################################
# Enable/Disable tracker module
if(SAF_ENABLE_TRACKER_MODULE)
    message(STATUS "saf_tracker module enabled. ")

    # Indicate to saf that the saf_sofa_reader module should be enabled
    target_compile_definitions(${PROJECT_NAME} PUBLIC SAF_ENABLE_TRACKER_MODULE=${SAF_ENABLE_TRACKER_MODULE})

else()
    message(STATUS "saf_tracker module disabled.")
endif()


############################################################################
if(UNIX)
    target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra)
endif()


############################################################################
# set RPATHs
set(CMAKE_INSTALL_RPATH $ORIGIN)
