cmake_minimum_required(VERSION 3.15)

# Pull version number from main header file
file(READ "modules/include/saf.h" ver)
string(REGEX MATCH "SAF_VERSION_MAJOR ([0-9]*)" _ ${ver})
set(ver_major ${CMAKE_MATCH_1})
string(REGEX MATCH "SAF_VERSION_MINOR ([0-9]*)" _ ${ver})
set(ver_minor ${CMAKE_MATCH_1})
string(REGEX MATCH "SAF_VERSION_PATCH ([0-9]*)" _ ${ver})
set(ver_patch ${CMAKE_MATCH_1})
set(SAF_VERSION "${ver_major}.${ver_minor}.${ver_patch}")
message("SAF Version: ${SAF_VERSION}")

# Configure project
project(saf VERSION ${SAF_VERSION} LANGUAGES C)
add_library(${PROJECT_NAME} SHARED)
set_target_properties(${PROJECT_NAME}
PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION}
)

# Framework source and include files pulled from subdirectories
add_subdirectory(modules)
add_subdirectory(resources)

include(GNUInstallDirs)

# User settings
add_compile_definitions(SAF_USE_INTEL_MKL)

set(INTEL_MKL_HEADER_PATH "/opt/intel/compilers_and_libraries/mac/mkl/include")
set(INTEL_MKL_LIB "/usr/local/lib/libsaf_mkl_custom.dylib")

target_include_directories(${PROJECT_NAME} 
PRIVATE 
    ${INTEL_MKL_HEADER_PATH}
)
target_link_libraries(${PROJECT_NAME}
PRIVATE 
    ${INTEL_MKL_LIB}
)

# sets the _EXPORT #define for the target (library)
# to ensure the functions are exported (visible)
target_compile_definitions(${PROJECT_NAME}
PRIVATE 
    -D_EXPORT
)

# specify the target to install (calculator library defined above)
# set the export name <name>-config (does not need to match target name)
# also specify where the .dylib/.so/.dll+.lib file should be installed
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
 
install(
DIRECTORY
    ${CMAKE_CURRENT_LIST_DIR}/modules/include/
    ${CMAKE_CURRENT_LIST_DIR}/resources/include/
DESTINATION
    ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
)
