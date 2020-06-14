cmake_minimum_required(VERSION 3.15)

option(SAF_ENABLE_SOFA_READER "Enable/Disable the sofa reader module" OFF) 

# Check required variables are defined
if(NOT DEFINED SAF_PERF_LIB)
    message(FATAL_ERROR "You must set the SAF_PERF_LIB environment variable")
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

# SAF Banner:
message("
      _____     _____     ______
    (  ____ ) (  ___  ) (  ____  )
    | (    |/ | (   ) | | (    |/
    | (_____  | (___) | | (__
    (_____  ) |  ___  | |  __)
          ) | | (   ) | | (
    /|____) | | |   | | | |
    (_______) |_|   |_| |_|        (Version ${SAF_VERSION})
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
message(${STATUS_2} "${CMAKE_SYSTEM_INFO_FILE}, name: ${CMAKE_SYSTEM_NAME}, arch: ${CMAKE_SYSTEM_PROCESSOR}, sys: ${CMAKE_SYSTEM}")
if(BUILD_SHARED_LIBS)
    message(${STATUS_2} "Configuring ${PROJECT_NAME} as a SHARED library...")
else()
    message(${STATUS_2} "Configuring ${PROJECT_NAME} as a STATIC library...")
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

# Performance library dependencies
if( ${SAF_PERF_LIB} MATCHES "SAF_USE_INTEL_MKL" )
    message(${STATUS_2} "Using Intel MKL. Linking against saf_mkl_custom library...")
    
    # Indicate to SAF, that it should use Intel MKL as the performance library
    target_compile_definitions(${PROJECT_NAME} PUBLIC SAF_USE_INTEL_MKL=1)

    # Header path and libs are platform dependent
    if(MSVC)
        set(INTEL_MKL_HEADER_PATH "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/include")
        set(INTEL_MKL_LIB "C:/Windows/System32/saf_mkl_custom.dll") 
        set(INTEL_MKL_LIB_WIN_LOCAL "${CMAKE_CURRENT_SOURCE_DIR}/../_ext_libs/Win64/lib/saf_mkl_custom.lib") 
        target_link_libraries(${PROJECT_NAME} PRIVATE ${INTEL_MKL_LIB_WIN_LOCAL} )
    elseif(MSYS OR MINGW) 
        set(INTEL_MKL_HEADER_PATH "/c/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/include")
        set(INTEL_MKL_LIB "/c/Windows/System32/saf_mkl_custom.dll") 
        set(INTEL_MKL_LIB_WIN_LOCAL "${CMAKE_CURRENT_SOURCE_DIR}/../_ext_libs/Win64/lib/saf_mkl_custom.lib") 
        target_link_libraries(${PROJECT_NAME} PRIVATE ${INTEL_MKL_LIB_WIN_LOCAL} )
    elseif(APPLE)
        set(INTEL_MKL_HEADER_PATH "/opt/intel/compilers_and_libraries/mac/mkl/include")
        set(INTEL_MKL_LIB "/usr/local/lib/libsaf_mkl_custom.dylib")
    elseif(UNIX AND NOT APPLE)
        set(INTEL_MKL_HEADER_PATH "/home/leo/intel/compilers_and_libraries/linux/mkl/include")
        set(INTEL_MKL_LIB "/usr/lib/libsaf_mkl_custom.so")
    endif()

    # Add target header paths and libs
    target_include_directories(${PROJECT_NAME} PUBLIC ${INTEL_MKL_HEADER_PATH} )
    target_link_libraries(${PROJECT_NAME} PRIVATE ${INTEL_MKL_LIB} )
elseif( ${SAF_PERF_LIB} MATCHES "SAF_USE_OPEN_BLAS_AND_LAPACKE")
    message(${STATUS_2} "Using openBLAS...")
    message(SEND_ERROR "Not yet supported! Please contribute if you use this library.")
elseif( ${SAF_PERF_LIB} MATCHES "SAF_USE_ATLAS")
    message(${STATUS_2} "Using ATLAS...")
    message(SEND_ERROR "Not yet supported! Please contribute if you use this library.")
else()
    message(SEND_ERROR "Unsupported performance library") 
endif()

# Sofa reader module dependencies
if(SAF_ENABLE_SOFA_READER) 
    message(${STATUS_2} "saf_sofa_reader module enabled. Linking against netcdf...")

    # Indicate to saf that the saf_sofa_reader module should be enabled
    target_compile_definitions(${PROJECT_NAME} PUBLIC SAF_ENABLE_SOFA_READER=${SAF_ENABLE_SOFA_READER})
        
    # Header path and libs are platform dependent
    if(MSVC) 
        target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../_ext_libs/Win64/include)
        target_link_libraries(${PROJECT_NAME} 
        INTERFACE
            ${CMAKE_CURRENT_SOURCE_DIR}/../_ext_libs/Win64/lib/netcdf.lib
            ${CMAKE_CURRENT_SOURCE_DIR}/../_ext_libs/Win64/lib/libhdf5.lib
            ${CMAKE_CURRENT_SOURCE_DIR}/../_ext_libs/Win64/lib/libzlib.lib
            ${CMAKE_CURRENT_SOURCE_DIR}/../_ext_libs/Win64/lib/libszip.lib
            ${CMAKE_CURRENT_SOURCE_DIR}/../_ext_libs/Win64/lib/libhdf5_hl.lib
        ) 
    elseif(MSYS OR MINGW)
	    message(SEND_ERROR "Unsupported on MSYS because I can't get the damn thing to link no matter what I do..") 
        target_include_directories(${PROJECT_NAME} PUBLIC 
		    /c/msys64/mingw64/include)
		#set (CMAKE_SHARED_LINKER_FLAGS "-lnetcdf")
		find_library(NETCDF_LIBRARY protobuf HINTS "/c/msys64/mingw64/lib/libnetcdf.dll.a")
		target_link_libraries(${PROJECT_NAME} INTERFACE ${NETCDF_LIBRARY})
		 
    elseif(APPLE)
        target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../_ext_libs/MacOSX/include)
        target_link_libraries(${PROJECT_NAME} 
        INTERFACE
            ${CMAKE_CURRENT_SOURCE_DIR}/../_ext_libs/MacOSX/lib/libnetcdf.a
            ${CMAKE_CURRENT_SOURCE_DIR}/../_ext_libs/MacOSX/lib/libhdf5.a
            ${CMAKE_CURRENT_SOURCE_DIR}/../_ext_libs/MacOSX/lib/libz.a 
            ${CMAKE_CURRENT_SOURCE_DIR}/../_ext_libs/MacOSX/lib/libhdf5_hl.a
        )  
    elseif(UNIX AND NOT APPLE) 
        find_library(NETCDF_LIBRARY protobuf HINTS /usr/lib/x86_64-linux-gnu)
        target_link_libraries(${PROJECT_NAME} PUBLIC ${NETCDF_LIBRARY})
    endif() 

else()
    message(${STATUS_2} "saf_sofa_reader module disabled.")
endif()
    
# set RPATHs
if(UNIX AND NOT APPLE)
    set (CMAKE_SHARED_LINKER_FLAGS -rpath,/usr/lib)
elseif(APPLE)
    set (CMAKE_SHARED_LINKER_FLAGS -rpath,/usr/local/lib)
elseif(MSYS)
    set (CMAKE_SHARED_LINKER_FLAGS -rpath,/c/msys64/mingw64/lib)
elseif(MSVC) 
    set (CMAKE_SHARED_LINKER_FLAGS -rpath,C:/Windows/System32)
endif()
