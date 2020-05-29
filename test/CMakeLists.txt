project(test LANGUAGES CXX C)

find_package(saf REQUIRED CONFIG)

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

target_link_libraries(${PROJECT_NAME} PUBLIC saf)
