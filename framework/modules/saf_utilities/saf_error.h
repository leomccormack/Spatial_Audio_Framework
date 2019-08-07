/*
 * Copyright 2019 Leo McCormack
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH
 * REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY
 * AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT,
 * INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
 * LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
 * OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
 * PERFORMANCE OF THIS SOFTWARE.
 */

/*
 * Filename: saf_error.h
 * ---------------------
 * List of error and warning codes.
 *
 * Dependencies:
 *     none
 * Author, date created:
 *     Leo McCormack, 05.08.2019
 */

#ifndef __SAF_ERROR_H_INCLUDED__
#define __SAF_ERROR_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "stdio.h"

/* ========================================================================== */
/*                                    Enums                                   */
/* ========================================================================== */
    
/*
 * Enum: SAF_ERRORS
 * ----------------
 * Error and warning codes. Error codes are considered fatal. Whereas warnings
 * are given if alternative measures were conducted due to some kind of strange
 * behaviour, but the program may still continue.
 *
 * Error Codes:
 *     SAF_ERROR__NO_ERROR
 *         No error was encountered.
 *     SAF_ERROR__ILLEGAL_INPUT_VALUE
 *         One or more input variable is assigned an illegal value.
 *     SAF_ERROR__UNALLOCATED_FUNCTION_ARGUMENT
 *         One or more input/output variable is NULL.
 *     SAF_ERROR__FAILED_TO_BUILD_CONVEX_HULL
 *         findLsTriplets - Failed to build Convex Hull.
 * Warning Codes:
 *     SAF_WARNING__SOFA_FILE_NOT_FOUND
 *         loadSofaFile - sofa file was not found at the specified directory.
 *         rememeber to include the ".sofa" suffix. In this case, the HRIR data
 *         is returned as NULL. The default HRIR set could be loaded instead by
 *         your program, or prompt the user to enter the path again.
 *     SAF_WARNING__UNABLE_TO_COMPUTE_BESSEL_FUNCTION_AT_SPECIFIED_ORDER
 *         bessel_jn/bessel_in/bessel_yn/bessel_kn/hankel_hn1/hankel_hn2 -
 *         Unable to compute the spherical Bessel/Hankel function at the
 *         specified order and input value. In this case, the Bessel/Hankel
 *         functions are returned at the maximum order that was possible, and
 *         this maximum order is returned by the function.
 *     SAF_WARNING__FAILED_TO_COMPUTE_SVD
 *         utility_?svd/utility_?pinv - The SVD failed to converge, or the input
 *         matrix contained illegal values so no solution was attempted. In
 *         these cases the function will zero all output matrices and vectors.
 *     SAF_WARNING__FAILED_TO_COMPUTE_EVG
 *         utility_?seig/utility_?eigmp/utility_?eig - Failed to compute all of
 *         the eigenvalues, no eigenvectors have been computed, or the input
 *         matrix contained illegal values so no solution was attempted. In
 *         these cases the function will zero all output matrices and vectors.
 *     SAF_WARNING__FAILED_TO_SOLVE_LINEAR_EQUATION
 *         utility_?glslv/utility_?slslv - Input matrix was singular, solution
 *         not computed, or the input matrix contained illegal values so no
 *         solution was attempted. In these cases the function will zero the
 *         output matrix.
 *     SAF_WARNING__FAILED_TO_COMPUTE_CHOL
 *         utility_?chol - input matrix is not positive definite, and the
 *         Cholesky factorization could not be computed, or the input matrix
 *         contained illegal values so no solution was attempted. In these cases
 *         the function will zero the output matrix.
 */
typedef enum _SAF_ERRORS {
    
    /* errors */
    SAF_ERROR__NO_ERROR = 0,
    SAF_ERROR__ILLEGAL_INPUT_VALUE,
    SAF_ERROR__UNALLOCATED_FUNCTION_ARGUMENT,
    
    /* saf_vbap errors */
    SAF_ERROR__FAILED_TO_BUILD_CONVEX_HULL,
    
    /* saf_hrir warnings */
    SAF_WARNING__SOFA_FILE_NOT_FOUND,
    
    /* saf_sh warnings */
    SAF_WARNING__UNABLE_TO_COMPUTE_BESSEL_FUNCTION_AT_SPECIFIED_ORDER,
    
    /* saf_utilites warnings */
    SAF_WARNING__FAILED_TO_COMPUTE_SVD,
    SAF_WARNING__FAILED_TO_COMPUTE_EVG,
    SAF_WARNING__FAILED_TO_SOLVE_LINEAR_EQUATION,
    SAF_WARNING__FAILED_TO_COMPUTE_CHOL
    
} SAF_ERRORS;
  

/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/*
 * Function: saf_error_print
 * -------------------------
 * Checks current error/warning code. Used when in debug mode.
 *
 * If there is no error/warning (SAF_ERROR__NO_ERROR), then the function does
 * nothing.
 * If there is a warning code, then an appropriate warning message is printed,
 * and err is reset upon return (if needed).
 * If there is an error code, then an appropriate error message is printed, and
 * the program is terminated.
 *
 * Input Arguments:
 *     err - & address saf error code (see "SAF_ERRORS" enum)
 * returns
 *     SAF_ERROR__NO_ERROR
 */
#ifndef NDEBUG
SAF_ERRORS saf_error_print(SAF_ERRORS err);
#endif /* NDEBUG */
    
    
#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __SAF_ERROR_H_INCLUDED__ */
