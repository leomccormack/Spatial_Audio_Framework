# Source
add_subdirectory(modules)
add_subdirectory(resources)

# Main include header
target_include_directories(${PROJECT_NAME}
PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>  
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)