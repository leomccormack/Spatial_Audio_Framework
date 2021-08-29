
# Source
target_sources(${PROJECT_NAME}
PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/saf_sofa_reader.c
    ${CMAKE_CURRENT_SOURCE_DIR}/libmysofa/mysofa.c
    ${CMAKE_CURRENT_SOURCE_DIR}/libmysofa/internal/hdf_dataobject.c
    ${CMAKE_CURRENT_SOURCE_DIR}/libmysofa/internal/hdf_fractalhead.c
    ${CMAKE_CURRENT_SOURCE_DIR}/libmysofa/internal/hdf_reader.c
    ${CMAKE_CURRENT_SOURCE_DIR}/libmysofa/internal/kdtree.c
    ${CMAKE_CURRENT_SOURCE_DIR}/libmysofa/internal/mysofa_internal.c
)
