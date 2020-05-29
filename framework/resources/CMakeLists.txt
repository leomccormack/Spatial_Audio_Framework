# Source
target_sources(${PROJECT_NAME} 
PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}/afSTFT/afSTFT_protoFilter.h
    ${CMAKE_CURRENT_SOURCE_DIR}/afSTFT/afSTFTlib.c
    ${CMAKE_CURRENT_SOURCE_DIR}/convhull_3d/convhull_3d.c
    ${CMAKE_CURRENT_SOURCE_DIR}/kissFFT/_kiss_fft_guts.h
    ${CMAKE_CURRENT_SOURCE_DIR}/kissFFT/kiss_fft.c
    ${CMAKE_CURRENT_SOURCE_DIR}/kissFFT/kiss_fftr.c
    ${CMAKE_CURRENT_SOURCE_DIR}/md_malloc/md_malloc.c
)