# Source
target_sources(${PROJECT_NAME} 
PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/afSTFT_protoFilter.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/afSTFTlib.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/convhull_3d.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/_kiss_fft_guts.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/kiss_fft.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/kiss_fftr.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/md_malloc.c
)

target_include_directories(${PROJECT_NAME}
PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>  
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
