/*******************************************************************************
* Copyright 1999-2017 Intel Corporation All Rights Reserved.
*
* The source code,  information  and material  ("Material") contained  herein is
* owned by Intel Corporation or its  suppliers or licensors,  and  title to such
* Material remains with Intel  Corporation or its  suppliers or  licensors.  The
* Material  contains  proprietary  information  of  Intel or  its suppliers  and
* licensors.  The Material is protected by  worldwide copyright  laws and treaty
* provisions.  No part  of  the  Material   may  be  used,  copied,  reproduced,
* modified, published,  uploaded, posted, transmitted,  distributed or disclosed
* in any way without Intel's prior express written permission.  No license under
* any patent,  copyright or other  intellectual property rights  in the Material
* is granted to  or  conferred  upon  you,  either   expressly,  by implication,
* inducement,  estoppel  or  otherwise.  Any  license   under such  intellectual
* property rights must be express and approved by Intel in writing.
*
* Unless otherwise agreed by Intel in writing,  you may not remove or alter this
* notice or  any  other  notice   embedded  in  Materials  by  Intel  or Intel's
* suppliers or licensors in any way.
*******************************************************************************/

/*
!  Content:
!     Intel(R) Math Kernel Library (Intel(R) MKL) interface for service routines
!******************************************************************************/

#ifndef _MKL_SERVICE_H_
#define _MKL_SERVICE_H_

#include <stdlib.h>
#include "mkl_types.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#if !defined(MKL_CALL_CONV)
#   if defined(__MIC__) || defined(__TARGET_ARCH_MIC)
#       define MKL_CALL_CONV
#   else
#       if defined(MKL_STDCALL)
#           define MKL_CALL_CONV __stdcall
#       else
#           define MKL_CALL_CONV __cdecl
#       endif
#   endif
#endif

#if !defined(_Mkl_Api)
#define _Mkl_Api(rtype,name,arg)    extern rtype MKL_CALL_CONV name arg;
#endif

_Mkl_Api(void,MKL_Get_Version,(MKLVersion *ver)) /* Returns information about the version of the Intel(R) MKL software */
#define mkl_get_version             MKL_Get_Version

_Mkl_Api(void,MKL_Get_Version_String,(char * buffer, int len)) /* Returns a string that contains Intel(R) MKL version information */
#define mkl_get_version_string      MKL_Get_Version_String

_Mkl_Api(void,MKL_Free_Buffers,(void)) /* Frees the memory allocated by the Intel(R) MKL Memory Manager */
#define mkl_free_buffers            MKL_Free_Buffers

_Mkl_Api(void,MKL_Thread_Free_Buffers,(void)) /* Frees the memory allocated by the Intel(R) MKL Memory Manager in the current thread only */
#define mkl_thread_free_buffers     MKL_Thread_Free_Buffers

_Mkl_Api(MKL_INT64,MKL_Mem_Stat,(int* nbuffers)) /* Intel(R) MKL Memory Manager statistical information. */
                                                 /* Returns an amount of memory, allocated by the Intel(R) MKL Memory Manager */
                                                 /* in <nbuffers> buffers. */
#define mkl_mem_stat                MKL_Mem_Stat

#define  MKL_PEAK_MEM_DISABLE       0
#define  MKL_PEAK_MEM_ENABLE        1
#define  MKL_PEAK_MEM_RESET        -1
#define  MKL_PEAK_MEM               2
_Mkl_Api(MKL_INT64,MKL_Peak_Mem_Usage,(int reset))    /* Returns the peak amount of memory, allocated by the Intel(R) MKL Memory Manager */
#define mkl_peak_mem_usage          MKL_Peak_Mem_Usage

_Mkl_Api(void*,MKL_malloc,(size_t size, int align)) /* Allocates the aligned buffer */
#define mkl_malloc                  MKL_malloc

_Mkl_Api(void*,MKL_calloc,(size_t num, size_t size, int align)) /* Allocates the aligned num*size - bytes memory buffer initialized by zeros */
#define mkl_calloc                  MKL_calloc

_Mkl_Api(void*,MKL_realloc,(void *ptr, size_t size)) /* Changes the size of memory buffer allocated by MKL_malloc/MKL_calloc */
#define mkl_realloc                  MKL_realloc

_Mkl_Api(void,MKL_free,(void *ptr))                 /* Frees the memory allocated by MKL_malloc() */
#define mkl_free                    MKL_free

_Mkl_Api(int,MKL_Disable_Fast_MM,(void))            /* Turns off the Intel(R) MKL Memory Manager */
#define  mkl_disable_fast_mm        MKL_Disable_Fast_MM

_Mkl_Api(void,MKL_Get_Cpu_Clocks,(unsigned MKL_INT64 *)) /* Gets CPU clocks */
#define mkl_get_cpu_clocks          MKL_Get_Cpu_Clocks

_Mkl_Api(double,MKL_Get_Cpu_Frequency,(void)) /* Gets CPU frequency in GHz */
#define mkl_get_cpu_frequency       MKL_Get_Cpu_Frequency

_Mkl_Api(double,MKL_Get_Max_Cpu_Frequency,(void)) /* Gets max CPU frequency in GHz */
#define mkl_get_max_cpu_frequency   MKL_Get_Max_Cpu_Frequency

_Mkl_Api(double,MKL_Get_Clocks_Frequency,(void)) /* Gets clocks frequency in GHz */
#define mkl_get_clocks_frequency    MKL_Get_Clocks_Frequency

_Mkl_Api(int,MKL_Set_Num_Threads_Local,(int nth))
#define mkl_set_num_threads_local   MKL_Set_Num_Threads_Local
_Mkl_Api(void,MKL_Set_Num_Threads,(int nth))
#define mkl_set_num_threads         MKL_Set_Num_Threads
_Mkl_Api(int,MKL_Get_Max_Threads,(void))
#define mkl_get_max_threads         MKL_Get_Max_Threads
_Mkl_Api(void,MKL_Set_Num_Stripes,(int nstripes))
#define mkl_set_num_stripes         MKL_Set_Num_Stripes
_Mkl_Api(int,MKL_Get_Num_Stripes,(void))
#define mkl_get_num_stripes         MKL_Get_Num_Stripes
_Mkl_Api(int,MKL_Domain_Set_Num_Threads,(int nth, int MKL_DOMAIN))
#define mkl_domain_set_num_threads  MKL_Domain_Set_Num_Threads
_Mkl_Api(int,MKL_Domain_Get_Max_Threads,(int MKL_DOMAIN))
#define mkl_domain_get_max_threads  MKL_Domain_Get_Max_Threads
_Mkl_Api(void,MKL_Set_Dynamic,(int bool_MKL_DYNAMIC))
#define mkl_set_dynamic             MKL_Set_Dynamic
_Mkl_Api(int,MKL_Get_Dynamic,(void))
#define mkl_get_dynamic             MKL_Get_Dynamic

/* Intel(R) MKL Progress routine */
#ifndef _MKL_PROGRESS_H_
#define _MKL_PROGRESS_H_
_Mkl_Api(int,MKL_PROGRESS, ( int* thread, int* step, char* stage, int lstage ))
_Mkl_Api(int,MKL_PROGRESS_,( int* thread, int* step, char* stage, int lstage ))
_Mkl_Api(int,mkl_progress, ( int* thread, int* step, char* stage, int lstage ))
_Mkl_Api(int,mkl_progress_,( int* thread, int* step, char* stage, int lstage ))
#endif /* _MKL_PROGRESS_H_ */

_Mkl_Api(int,MKL_Enable_Instructions,(int))
#define  mkl_enable_instructions    MKL_Enable_Instructions
#define  MKL_ENABLE_SSE4_2          0
#define  MKL_ENABLE_AVX             1
#define  MKL_ENABLE_AVX2            2
#define  MKL_ENABLE_AVX512_MIC      3
#define  MKL_ENABLE_AVX512          4
#define  MKL_ENABLE_AVX512_MIC_E1   5
#define  MKL_SINGLE_PATH_ENABLE     0x0600

/* Single Dynamic library interface */
#define MKL_INTERFACE_LP64          0x0
#define MKL_INTERFACE_ILP64         0x1
#define MKL_INTERFACE_GNU           0x2
_Mkl_Api(int,MKL_Set_Interface_Layer,(int code))
#define mkl_set_interface_layer     MKL_Set_Interface_Layer

/* Single Dynamic library threading */
#define MKL_THREADING_INTEL         0
#define MKL_THREADING_SEQUENTIAL    1
#define MKL_THREADING_PGI           2
#define MKL_THREADING_GNU           3
#define MKL_THREADING_TBB           4
_Mkl_Api(int,MKL_Set_Threading_Layer,(int code))
#define mkl_set_threading_layer     MKL_Set_Threading_Layer

typedef void (* XerblaEntry) (const char * Name, const int * Num, const int Len);
_Mkl_Api(XerblaEntry,mkl_set_xerbla,(XerblaEntry xerbla))

typedef int (* ProgressEntry) (int* thread, int* step, char* stage, int stage_len);
_Mkl_Api(ProgressEntry,mkl_set_progress,(ProgressEntry progress))

/* MIC service routines */
_Mkl_Api(int,MKL_MIC_Enable,(void))
#define mkl_mic_enable              MKL_MIC_Enable
_Mkl_Api(int,MKL_MIC_Disable,(void))
#define mkl_mic_disable             MKL_MIC_Disable

_Mkl_Api(int,MKL_MIC_Get_Device_Count,(void))
#define mkl_mic_get_device_count MKL_MIC_Get_Device_Count

typedef enum MKL_MIC_TARGET_TYPE {
    MKL_TARGET_NONE = 0, /* Undefine target */
    MKL_TARGET_HOST = 1, /* Host used as target */
    MKL_TARGET_MIC  = 2  /* MIC target */
} MKL_MIC_TARGET_TYPE;

#define MKL_MIC_DEFAULT_TARGET_TYPE MKL_TARGET_MIC
#define MKL_MIC_DEFAULT_TARGET_NUMBER 0
#define MKL_MIC_AUTO_WORKDIVISION   -1.0
#define MKL_MPI_PPN -1.0

_Mkl_Api(int,MKL_MIC_Get_Cpuinfo,(MKL_MIC_TARGET_TYPE target_type, int target_number, int *ncores, int *nthreads, double *freq))
#define mkl_mic_get_cpuinfo MKL_MIC_Get_Cpuinfo

_Mkl_Api(int,MKL_MIC_Get_Meminfo,(MKL_MIC_TARGET_TYPE target_type, int target_number, int *totalmem, int *freemem))
#define mkl_mic_get_meminfo MKL_MIC_Get_Meminfo

_Mkl_Api(int,MKL_MIC_Set_Workdivision,(MKL_MIC_TARGET_TYPE target_type,
                                       int target_number, double wd))
#define mkl_mic_set_workdivision    MKL_MIC_Set_Workdivision

_Mkl_Api(int,MKL_MIC_Get_Workdivision,(MKL_MIC_TARGET_TYPE target_type,
                                       int target_number, double *wd))
#define mkl_mic_get_workdivision    MKL_MIC_Get_Workdivision

_Mkl_Api(int,MKL_MIC_Set_Max_Memory,(MKL_MIC_TARGET_TYPE target_type,
                                     int target_number, size_t card_mem_mbytes))
#define mkl_mic_set_max_memory      MKL_MIC_Set_Max_Memory

_Mkl_Api(int,MKL_MIC_Free_Memory,(MKL_MIC_TARGET_TYPE target_type,
                                  int target_number))
#define mkl_mic_free_memory         MKL_MIC_Free_Memory

_Mkl_Api(int,MKL_MIC_Set_Offload_Report,(int enabled))
#define mkl_mic_set_offload_report  MKL_MIC_Set_Offload_Report

_Mkl_Api(int,MKL_MIC_Set_Device_Num_Threads,(MKL_MIC_TARGET_TYPE target_type,
                                             int target_number, int num_threads))
#define mkl_mic_set_device_num_threads MKL_MIC_Set_Device_Num_Threads

_Mkl_Api(int,MKL_MIC_Set_Resource_Limit,(double fraction))
#define mkl_mic_set_resource_limit    MKL_MIC_Set_Resource_Limit

_Mkl_Api(int,MKL_MIC_Get_Resource_Limit,(double *fraction))
#define mkl_mic_get_resource_limit    MKL_MIC_Get_Resource_Limit

#define MKL_MIC_DEFAULT_FLAGS       (0)
#define MKL_MIC_DISABLE_HOST_FALLBACK    (1 << 0)

_Mkl_Api(int,MKL_MIC_Get_Flags,(void))
#define mkl_mic_get_flags    MKL_MIC_Get_Flags

_Mkl_Api(int,MKL_MIC_Set_Flags,(int flag))
#define mkl_mic_set_flags    MKL_MIC_Set_Flags

#define MKL_MIC_SUCCESS            0
#define MKL_MIC_NOT_IMPL           1
#define MKL_MIC_HOST_FALLBACK      2
#define MKL_MIC_DISABLED           3
#define MKL_MIC_FAILED             -1
#define MKL_MIC_HOST_FALLBACK_DISABLED  -2

_Mkl_Api(int,MKL_MIC_Get_Status,(void))
#define mkl_mic_get_status    MKL_MIC_Get_Status

_Mkl_Api(void,MKL_MIC_Clear_Status,(void))
#define mkl_mic_clear_status    MKL_MIC_Clear_Status


/* Intel(R) MKL CBWR */
_Mkl_Api(int,MKL_CBWR_Get,(int))
#define mkl_cbwr_get                MKL_CBWR_Get
_Mkl_Api(int,MKL_CBWR_Set,(int))
#define mkl_cbwr_set                MKL_CBWR_Set
_Mkl_Api(int,MKL_CBWR_Get_Auto_Branch,(void))
#define mkl_cbwr_get_auto_branch    MKL_CBWR_Get_Auto_Branch

_Mkl_Api(int,MKL_Set_Env_Mode,(int))
#define mkl_set_env_mode            MKL_Set_Env_Mode

_Mkl_Api(int,MKL_Verbose,(int))
#define mkl_verbose                MKL_Verbose

#define MKL_EXIT_UNSUPPORTED_CPU    1
#define MKL_EXIT_CORRUPTED_INSTALL  2
#define MKL_EXIT_NO_MEMORY          3

typedef void (* MKLExitHandler)(int why);
_Mkl_Api(void,MKL_Set_Exit_Handler,(MKLExitHandler h));
#define mkl_set_exit_handler       MKL_Set_Exit_Handler

#define MKL_MIC_REGISTRATION_DISABLE 0
#define MKL_MIC_REGISTRATION_ENABLE  1
_Mkl_Api(void,MKL_MIC_register_memory,(int enable));
#define mkl_mic_register_memory    MKL_MIC_register_memory

#define MKL_MEM_MCDRAM 1

_Mkl_Api(int,MKL_Set_Memory_Limit,(int mem_type,size_t limit));
#define mkl_set_memory_limit MKL_Set_Memory_Limit

enum {
    MKL_BLACS_CUSTOM = 0,
    MKL_BLACS_MSMPI = 1,
    MKL_BLACS_INTELMPI = 2,
    MKL_BLACS_MPICH2 = 3,
    MKL_BLACS_LASTMPI = 4,
};
int MKL_Set_mpi(int vendor, const char *custom_library_name);
#define mkl_set_mpi MKL_Set_mpi

_Mkl_Api(void,MKL_Finalize,(void));
#define mkl_finalize MKL_Finalize

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _MKL_SERVICE_H_ */
