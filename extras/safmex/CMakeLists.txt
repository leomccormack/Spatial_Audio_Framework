message(STATUS "Configuring SAFMEX wrappers:")

cmake_minimum_required(VERSION 3.15)

if(SAF_MEX_CXX)
    add_subdirectory(CXX_API) # new experimental mex versions using c++ API
endif()

project(safmex LANGUAGES C)

if(NOT Matlab_ROOT_DIR)
    # Try to find MATLAB installations
    if(APPLE)
	file(GLOB matlab_candidates "/Applications/MATLAB_R20*.app")
    elseif(UNIX AND NOT APPLE) 
        file(GLOB matlab_candidates
            "/usr/local/MATLAB/R20*"
            "/usr/local/MATLAB/R20*/*"
            "/usr/local/MATLAB/*R20*"
        )
    elseif(MSVC)
        file(GLOB matlab_candidates
            "C:/Program Files/MATLAB/R20*"
            "C:/Program Files (x86)/MATLAB/R20*"
        )
    endif()

    # If MATLAB installation were found then we pick the latest
    if(matlab_candidates)
        list(SORT matlab_candidates) 
        list(REVERSE matlab_candidates)
        list(GET matlab_candidates 0 Matlab_ROOT_DIR)
        message(STATUS "Auto-detected latest MATLAB as: ${Matlab_ROOT_DIR}")
    else()
        message(WARNING "No MATLAB installations found automatically. Matlab_ROOT_DIR must be set with e.g.: -DMatlab_ROOT_DIR=/usr/local/MATLAB/R2025b")
    endif()
endif()

#find_package(Matlab REQUIRED COMPONENTS MEX_LIBRARY)
set(MATLAB_EXECUTABLE "${Matlab_ROOT_DIR}/bin/matlab")

# List of available safmex wrappers
set(safmex_list
    safmex_afSTFT
    safmex_faf_IIRFilterbank
    safmex_generateVBAPgainTable3D
    safmex_getSHcomplex
    safmex_getSHreal
    safmex_latticeDecorrelator
    safmex_qmf) 

if(SAF_ENABLE_TRACKER_MODULE)
    set(safmex_list ${safmex_list}
        safmex_tracker3d)
else()
    message(STATUS "  (Disabling safmex wrappers related to the saf_tracker module)")
endif()

# Function to build a MEX file using MATLAB's mex command
function(add_matlab_mex target src)
    set(MATLAB_EXECUTABLE "${Matlab_ROOT_DIR}/bin/matlab")
    set(SAF_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/../../framework/include")
    set(SAF_LIB_DIR "${PROJECT_BINARY_DIR}/../../../build/framework")

    # Mex extension is platform specific:
    if(APPLE)
        set(MEX_EXT "mexmaca64")
    elseif(UNIX AND NOT APPLE)
        set(MEX_EXT "mexa64")
    elseif(MSVC)
        set(MEX_EXT "mexw64")
    endif()

    # Create a string of mex arguments
    set(MEX_CMD_LIST mex)
    list(APPEND MEX_CMD_LIST -I${SAF_INCLUDE_DIR})
    list(APPEND MEX_CMD_LIST -L${SAF_LIB_DIR})
    list(APPEND MEX_CMD_LIST -lsaf)
    if( ${SAF_PERFORMANCE_LIB} MATCHES "SAF_USE_APPLE_ACCELERATE_LP64" OR ${SAF_PERFORMANCE_LIB} MATCHES "SAF_USE_APPLE_ACCELERATE_ILP64" )
        list(APPEND MEX_CMD_LIST LDFLAGS="$LDFLAGS -framework Accelerate")
    endif()
    if(SAF_ENABLE_TRACKER_MODULE)
        list(APPEND MEX_CMD_LIST -DSAF_ENABLE_TRACKER_MODULE=1)
    endif()
    list(APPEND MEX_CMD_LIST ${CMAKE_CURRENT_SOURCE_DIR}/${src})
    string(JOIN " " MEX_CMD ${MEX_CMD_LIST})
    #message(STATUS "MEX_CMD = ${MEX_CMD}")

    add_custom_command(
        OUTPUT ${target}.${MEX_EXT}
        COMMAND "${MATLAB_EXECUTABLE}" -batch "${MEX_CMD}"
        DEPENDS ${src} saf
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Building MEX: ${target}"
        VERBATIM
    )

    add_custom_target(${target} ALL DEPENDS ${target}.${MEX_EXT})
endfunction()

# Find number of wrappers
list(LENGTH safmex_list num_wrappers)
math(EXPR num_wrappers "${num_wrappers} - 1")
  
# Loop over the wrappers
foreach(wrapper ${safmex_list})
    message(STATUS "  Building MEX wrapper: ${wrapper}")
    add_matlab_mex(${wrapper} ${wrapper}.c)
endforeach()

# make install safmex
install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/*.mexmaca64
    DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/
    COMPONENT safmex
)
