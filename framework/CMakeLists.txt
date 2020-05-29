cmake_minimum_required(VERSION 3.15)

# Configure project
project(saf VERSION 1.0.0 LANGUAGES C)
add_library(${PROJECT_NAME} SHARED)
set_target_properties(${PROJECT_NAME}
PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION}
)

# Source
add_subdirectory(modules)
add_subdirectory(resources)

# Main include header
target_include_directories(${PROJECT_NAME}
PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>  
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

include(GNUInstallDirs)

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
    DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/include/ 
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})
