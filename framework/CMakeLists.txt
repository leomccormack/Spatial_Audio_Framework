cmake_minimum_required(VERSION 3.15)

option(SAF_ENABLE_SOFA_READER "Enable/Disable the sofa reader module" OFF) 

# Check required variables are defined
if(NOT DEFINED SAF_PERF_LIB)
    message(FATAL_ERROR "You must set the SAF_PERF_LIB environment variable")
endif()

# Pull version number from main header file
file(READ "modules/include/saf.h" ver)
string(REGEX MATCH "SAF_VERSION_MAJOR ([0-9]*)" _ ${ver})
set(ver_major ${CMAKE_MATCH_1})
string(REGEX MATCH "SAF_VERSION_MINOR ([0-9]*)" _ ${ver})
set(ver_minor ${CMAKE_MATCH_1})
string(REGEX MATCH "SAF_VERSION_PATCH ([0-9]*)" _ ${ver})
set(ver_patch ${CMAKE_MATCH_1})
set(SAF_VERSION "${ver_major}.${ver_minor}.${ver_patch}")
message("-- SAF Version: ${SAF_VERSION}")

# Configure project
project(saf VERSION ${SAF_VERSION} LANGUAGES C)
add_library(${PROJECT_NAME} SHARED)
set_target_properties(${PROJECT_NAME}
PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION}
)

# Display system info
message("-- ${CMAKE_SYSTEM_INFO_FILE}, name: ${CMAKE_SYSTEM_NAME}, arch: ${CMAKE_SYSTEM_PROCESSOR}, sys: ${CMAKE_SYSTEM}")

# Framework source and include files are pulled from subdirectories
add_subdirectory(modules)
add_subdirectory(resources)
include(GNUInstallDirs)

# Performance library to employ
if( ${SAF_PERF_LIB} MATCHES "SAF_USE_INTEL_MKL" )
    message("-- Using Intel MKL. Linking against saf_mkl_custom library...")
    
    # Indicate to SAF, that it should use Intel MKL as the performance library
    #add_compile_definitions(SAF_USE_INTEL_MKL)
    target_compile_definitions(${PROJECT_NAME} PUBLIC SAF_USE_INTEL_MKL=1)

    # Header path and libs are platform dependent
    if(MSVC) 
        set(INTEL_MKL_HEADER_PATH "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/include")
        set(INTEL_MKL_LIB "C:/Windows/System32/saf_mkl_custom.dll") 
        set(INTEL_MKL_LIB_WIN_LOCAL "${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/Win64/lib/saf_mkl_custom.lib") 
        target_link_libraries(${PROJECT_NAME} PRIVATE ${INTEL_MKL_LIB_WIN_LOCAL} )
    elseif(MSYS OR MINGW) 
        set(INTEL_MKL_HEADER_PATH "/c/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/include")
        set(INTEL_MKL_LIB "/c/Windows/System32/saf_mkl_custom.dll") 
        set(INTEL_MKL_LIB_WIN_LOCAL "${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/Win64/lib/saf_mkl_custom.lib") 
        target_link_libraries(${PROJECT_NAME} PRIVATE ${INTEL_MKL_LIB_WIN_LOCAL} )
    elseif(APPLE)
        set(INTEL_MKL_HEADER_PATH "/opt/intel/compilers_and_libraries/mac/mkl/include")
        set(INTEL_MKL_LIB "/usr/local/lib/libsaf_mkl_custom.dylib")
    elseif(UNIX AND NOT APPLE)
        set(INTEL_MKL_HEADER_PATH "/opt/intel/compilers_and_libraries/mac/mkl/include")
        set(INTEL_MKL_LIB "/usr/lib/libsaf_mkl_custom.so")
    endif()

    # Add target header paths and libs
    target_include_directories(${PROJECT_NAME} PUBLIC ${INTEL_MKL_HEADER_PATH} )
    target_link_libraries(${PROJECT_NAME} PRIVATE ${INTEL_MKL_LIB} )
elseif( ${SAF_PERF_LIB} MATCHES "SAF_USE_OPEN_BLAS_AND_LAPACKE")
    message("-- Using openBLAS...")
    message(SEND_ERROR "Not yet supported! Please contribute if you use this library.")
elseif( ${SAF_PERF_LIB} MATCHES "SAF_USE_ATLAS")
    message("-- Using ATLAS...")
    message(SEND_ERROR "Not yet supported! Please contribute if you use this library.")
else()
    message(SEND_ERROR "Unsupported performance library") 
endif()

# Sofa reader module 
if(SAF_ENABLE_SOFA_READER) 
    message("-- saf_sofa_reader module enabled. Linking against netcdf...")

    # Indicate to saf that the saf_sofa_reader module should be enabled
    target_compile_definitions(${PROJECT_NAME} PRIVATE SAF_ENABLE_SOFA_READER=${SAF_ENABLE_SOFA_READER})
        
    # Header path and libs are platform dependent
    if(MSVC OR MSYS OR MINGW) 
        target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/Win64/include)
        target_link_libraries(${PROJECT_NAME} 
        PRIVATE 
            ${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/Win64/lib/netcdf.lib
            ${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/Win64/lib/libhdf5.lib
            ${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/Win64/lib/libzlib.lib
            ${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/Win64/lib/libszip.lib
            ${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/Win64/lib/libhdf5_hl.lib
        ) 
    elseif(APPLE)
        target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/MacOSX/include)
        target_link_libraries(${PROJECT_NAME} 
        PRIVATE 
            ${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/MacOSX/lib/libnetcdf.a
            ${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/MacOSX/lib/libhdf5.a
            ${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/MacOSX/lib/libz.a 
            ${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/MacOSX/lib/libhdf5_hl.a
        ) 
    elseif(UNIX AND NOT APPLE) 
        message(SEND_ERROR "Incomplete") 
    endif() 

else()
    message("-- saf_sofa_reader module disabled.")
endif()
    
 
# sets the _EXPORT #define for the target (library) to ensure the functions are exported (visible)
target_compile_definitions(${PROJECT_NAME}
PRIVATE 
    -D_EXPORT
)

# Specify the target to install 
# Set the export name <name>-config (does not need to match target name)
# Also specify where the .dylib/.so/.dll+.lib file should be installed
install(
    TARGETS ${PROJECT_NAME}
    EXPORT ${PROJECT_NAME}-config
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(
    EXPORT ${PROJECT_NAME}-config
    NAMESPACE ${PROJECT_NAME}::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
 
# Install header include folder
install(
DIRECTORY
    ${CMAKE_CURRENT_LIST_DIR}/modules/include/
    ${CMAKE_CURRENT_LIST_DIR}/resources/include/
DESTINATION
    ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
)

