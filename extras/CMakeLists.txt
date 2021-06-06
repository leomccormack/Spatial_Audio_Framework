
if( ${SAF_PERFORMANCE_LIB} MATCHES "SAF_USE_INTEL_MKL_ILP64" )
    add_subdirectory(safmex)
else()
    message(STATUS "safmex disabled, since it requires \"-DSAF_PERFORMANCE_LIB=SAF_USE_INTEL_MKL_ILP64\". It is currently set to \"${SAF_PERFORMANCE_LIB}\".")
endif()