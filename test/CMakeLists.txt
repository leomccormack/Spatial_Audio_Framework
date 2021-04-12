project(saf_test LANGUAGES CXX C)

message(STATUS "Configuring SAF unit testing program...")
add_executable(${PROJECT_NAME})
target_sources(${PROJECT_NAME} 
PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}/saf_test_wrapper.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/saf_test.c
    ${CMAKE_CURRENT_SOURCE_DIR}/timer/timer.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/unity/unity_internals.h
    ${CMAKE_CURRENT_SOURCE_DIR}/unity/unity.c 
)

target_include_directories(${PROJECT_NAME} 
PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/timer/>  
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/unity/>  
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# SAF
target_link_libraries(${PROJECT_NAME} PRIVATE saf)

# SAF examples tests
if(SAF_ENABLE_SOFA_READER_MODULE AND SAF_BUILD_EXAMPLES)
    target_compile_definitions(${PROJECT_NAME} PRIVATE SAF_ENABLE_EXAMPLES_TESTS=1) 
    target_link_libraries(${PROJECT_NAME} 
    PRIVATE 
        "${example_prefix}ambi_bin"
        "${example_prefix}ambi_dec"
        "${example_prefix}ambi_drc"
        "${example_prefix}ambi_enc"
        "${example_prefix}array2sh"
        "${example_prefix}beamformer"
        "${example_prefix}binauraliser"
        "${example_prefix}decorrelator"
        "${example_prefix}dirass"
        "${example_prefix}matrixconv"
        "${example_prefix}multiconv"
        "${example_prefix}panner"
        "${example_prefix}pitch_shifter"
        "${example_prefix}powermap"
        "${example_prefix}rotator"
        "${example_prefix}sldoa"
        "${example_prefix}spreader"
    ) 
else()
    target_compile_definitions(${PROJECT_NAME} PRIVATE SAF_ENABLE_EXAMPLES_TESTS=0) 
    message(STATUS "  Note: unit tests for the SAF examples have been disabled")
endif()

