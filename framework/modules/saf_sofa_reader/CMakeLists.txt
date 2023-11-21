
# Source
target_sources(${PROJECT_NAME}
PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/saf_sofa_reader.c
	${CMAKE_CURRENT_SOURCE_DIR}/libmysofa/config.h
	${CMAKE_CURRENT_SOURCE_DIR}/libmysofa/hdf/btree.c
	${CMAKE_CURRENT_SOURCE_DIR}/libmysofa/hdf/dataobject.c
	${CMAKE_CURRENT_SOURCE_DIR}/libmysofa/hdf/fractalhead.c
	${CMAKE_CURRENT_SOURCE_DIR}/libmysofa/hdf/gcol.c
	${CMAKE_CURRENT_SOURCE_DIR}/libmysofa/hdf/gunzip.c
	${CMAKE_CURRENT_SOURCE_DIR}/libmysofa/hdf/reader.h
	${CMAKE_CURRENT_SOURCE_DIR}/libmysofa/hdf/superblock.c
	${CMAKE_CURRENT_SOURCE_DIR}/libmysofa/hrtf/kdtree.c
	${CMAKE_CURRENT_SOURCE_DIR}/libmysofa/hrtf/kdtree.h
	${CMAKE_CURRENT_SOURCE_DIR}/libmysofa/hrtf/mysofa.h
	${CMAKE_CURRENT_SOURCE_DIR}/libmysofa/hrtf/mysofa_export.h
	${CMAKE_CURRENT_SOURCE_DIR}/libmysofa/hrtf/reader.c
	${CMAKE_CURRENT_SOURCE_DIR}/libmysofa/hrtf/tools.c
	${CMAKE_CURRENT_SOURCE_DIR}/libmysofa/hrtf/tools.h
)
