/*******************************************************************************
* Copyright 1999-2018 Intel Corporation.
*
* This software and the related documents are Intel copyrighted  materials,  and
* your use of  them is  governed by the  express license  under which  they were
* provided to you (License).  Unless the License provides otherwise, you may not
* use, modify, copy, publish, distribute,  disclose or transmit this software or
* the related documents without Intel's prior written permission.
*
* This software and the related documents  are provided as  is,  with no express
* or implied  warranties,  other  than those  that are  expressly stated  in the
* License.
*******************************************************************************/

/*
!  Content:
!      Intel(R) Math Kernel Library (Intel(R) MKL) C/C++ interface for
!	   Cluster Sparse Solver
!******************************************************************************/

#if !defined( __MKL_CLUSTER_SPARSE_SOLVER_H )

#include "mkl_types.h"

#define __MKL_CLUSTER_SPARSE_SOLVER_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

void cluster_sparse_solver(
     void *, const MKL_INT *, const MKL_INT *, const MKL_INT *, const MKL_INT *, const MKL_INT *,
     const void *, const MKL_INT *, const MKL_INT *, MKL_INT *, const MKL_INT *, MKL_INT *,
     const MKL_INT *,       void *,       void *, const int *, MKL_INT *);

void CLUSTER_SPARSE_SOLVER(
     void *, const MKL_INT *, const MKL_INT *, const MKL_INT *, const MKL_INT *, const MKL_INT *,
     const void *, const MKL_INT *, const MKL_INT *, MKL_INT *, const MKL_INT *, MKL_INT *,
     const MKL_INT *,       void *,       void *, const int *, MKL_INT *);

void cluster_sparse_solver_64(
     void *, const long long int *, const long long int *, const long long int *, const long long int *, const long long int *,
     const void *, const long long int *, const long long int *, long long int *, const long long int *, long long int *,
     const long long int *,       void *,       void *, const int *, long long int *);

void CLUSTER_SPARSE_SOLVER_64(
     void *, const long long int *, const long long int *, const long long int *, const long long int *, const long long int *,
     const void *, const long long int *, const long long int *, long long int *, const long long int *, long long int *,
     const long long int *,       void *,       void *, const int *, long long int *);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
