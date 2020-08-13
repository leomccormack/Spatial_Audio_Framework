message(STATUS "Configuring SAFMEX wrappers:")

cmake_minimum_required(VERSION 3.15)

project(safmex LANGUAGES C)
find_package(Matlab)

# List of available safmex wrappers
set(safmex_list
    safmex_afSTFT
    safmex_faf_IIRFilterbank
    safmex_generateVBAPgainTable3D
    safmex_getSHcomplex
    safmex_getSHreal
    safmex_latticeDecorrelator
    safmex_qmf)

# Find number of wrappers
list(LENGTH safmex_list num_wrappers)
math(EXPR num_wrappers "${num_wrappers} - 1")
  
# Loop over the wrappers
foreach(current_wrapper RANGE ${num_wrappers})
    # Get current example string
    list(GET safmex_list ${current_wrapper} current_string)

    message(STATUS "  ${current_string}")
    matlab_add_mex(NAME ${current_string} SRC ${current_string}.c LINK_TO saf)

endforeach()

 