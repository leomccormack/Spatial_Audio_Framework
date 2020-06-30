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

/**
 *@addtogroup Utilities
 *@{
 * @file saf_utility_error.h
 * @brief A list of error and warning codes
 *
 * @author Leo McCormack
 * @date 05.08.2019 
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
    
/**
 * Error and warning codes. Error codes are considered fatal. Whereas warnings
 * are given if alternative measures have taken place (due to some kind of
 * unexpected behaviour), but the program may still continue.
 */
typedef enum _SAF_ERRORS {
    /* ---------------------------------------------------------------------- */
    /**
     * No error was encountered.
     */
    SAF_ERROR__NO_ERROR = 0,
    /**
     * One or more input variable is assigned an illegal value.
     */
    SAF_ERROR__ILLEGAL_INPUT_VALUE,
    /**
     * One or more input/output variable is NULL.
     */
    SAF_ERROR__UNALLOCATED_FUNCTION_ARGUMENT, 
    /**
     * findLsTriplets - Failed to build Convex Hull.
     */
    SAF_ERROR__FAILED_TO_BUILD_CONVEX_HULL,
    
    /* ---------------------------------------------------------------------- */
    /**
     * loadSofaFile(): sofa file was not found at the specified directory.
     * Remember to include the ".sofa" suffix. In this case, the default HRIR
     * set is loaded instead.
     */
    SAF_WARNING__SOFA_FILE_NOT_FOUND,
    /**
     * bessel_jn(), bessel_in(), bessel_yn(), bessel_kn(), hankel_hn1(), or
     * hankel_hn2():
     * Unable to compute the spherical Bessel/Hankel function at the specified
     * order and input value. In this case, the Bessel/Hankel functions are
     * returned at the maximum order that was possible, and this maximum order
     * is returned by the function.
     */
    SAF_WARNING__UNABLE_TO_COMPUTE_BESSEL_FUNCTION_AT_SPECIFIED_ORDER,
    /**
     * utility_?svd/utility_?pinv - The SVD failed to converge, or the input
     * matrix contained illegal values so no solution was attempted. In these
     * cases the function will zero all output matrices and vectors.
     */
    SAF_WARNING__FAILED_TO_COMPUTE_SVD,
    /**
     * utility_?seig/utility_?eigmp/utility_?eig - Failed to compute all of the
     * eigenvalues, no eigenvectors have been computed, or the input matrix
     * contained illegal values so no solution was attempted. In these cases the
     * function will zero all output matrices and vectors.
     */
    SAF_WARNING__FAILED_TO_COMPUTE_EVG,
    /**
     * utility_?glslv/utility_?slslv - Input matrix was singular, solution not
     * computed, or the input matrix contained illegal values so no solution was
     * attempted. In these cases the function will zero the output matrix.
     */
    SAF_WARNING__FAILED_TO_SOLVE_LINEAR_EQUATION,
    /**
     * utility_?chol - input matrix is not positive definite, and the Cholesky
     * factorization could not be computed, or the input matrix contained
     * illegal values so no solution was attempted. In these cases the function
     * will zero the output matrix.
     */
    SAF_WARNING__FAILED_TO_COMPUTE_CHOL
    
} SAF_ERRORS;


/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/**
 * Checks current error/warning code, and prints out a message if needed (used
 * when in debug mode).
 *
 * If there is no error/warning (#SAF_ERROR__NO_ERROR), then the function does
 * nothing.
 * If there is a warning code, then an appropriate warning message is printed,
 * and err is reset upon return (if needed).
 * If there is an error code, then an appropriate error message is printed, and
 * the program is terminated.
 *
 * @param[in] err SAF error code (see #_SAF_ERRORS enum)
 * @returns #SAF_ERROR__NO_ERROR
 */
#ifndef NDEBUG
SAF_ERRORS saf_error_print(SAF_ERRORS err);
#endif /* NDEBUG */


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __SAF_ERROR_H_INCLUDED__ */

/**@} */ /* doxygen addtogroup Utilities */
