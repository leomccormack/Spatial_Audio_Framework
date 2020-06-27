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
 * @file saf_utility_error.c
 * @ingroup Utilities
 * @brief A list of error and warning codes
 *
 * @author Leo McCormack
 * @date 05.08.2019
 */

#include "saf_utility_error.h"

#ifndef NDEBUG
/*
 * NOTE: refer to "SAF_ERRORS" enum descriptions in "saf_error.h" for more
 *       detailed error/warning information.
 */
SAF_ERRORS saf_error_print(SAF_ERRORS err)
{
    switch(err){
        /* errors */
        case SAF_ERROR__NO_ERROR: break;
        case SAF_ERROR__ILLEGAL_INPUT_VALUE:
            fprintf(stderr, "%s", "SAF Error: One or more input variable was assigned an illegal value.\n");
            break;
        case SAF_ERROR__UNALLOCATED_FUNCTION_ARGUMENT:
            fprintf(stderr, "%s", "SAF Error: Memory for one or more input/output matrix/vector was not allocated.\n");
            break;
            
        /* saf_vbap errors */
        case SAF_ERROR__FAILED_TO_BUILD_CONVEX_HULL:
            fprintf(stderr, "%s", "SAF Error: Failed to build Convex Hull.\n");
            break;
            
        /* saf_hrir warnings */
        case SAF_WARNING__SOFA_FILE_NOT_FOUND:
            fprintf(stdout, "%s", "SAF Warning: Could not open SOFA file. Loading default HRIR data. \n");
            break;
            
        /* saf_sh warnings */
        case SAF_WARNING__UNABLE_TO_COMPUTE_BESSEL_FUNCTION_AT_SPECIFIED_ORDER:
            fprintf(stdout, "%s", "SAF Warning: Could not compute spherical Bessel/Hankel at specified order. \n");
            break;
        
        /* saf_utilites warnings */
        case SAF_WARNING__FAILED_TO_COMPUTE_SVD:
            fprintf(stdout, "%s", "SAF Warning: Could not compute SVD. Output matrices/vectors have been zeroed. \n");
            break;
        case SAF_WARNING__FAILED_TO_COMPUTE_EVG:
            fprintf(stdout, "%s", "SAF Warning: Could not compute EVD. Output matrices/vectors have been zeroed. \n");
            break;
        case SAF_WARNING__FAILED_TO_SOLVE_LINEAR_EQUATION:
            fprintf(stdout, "%s", "SAF Warning: Could not solve linear equation. Output matrix has been zeroed. \n");
            break;
        case SAF_WARNING__FAILED_TO_COMPUTE_CHOL:
            fprintf(stdout, "%s", "SAF Warning: Could not compute Choleksy Factorisation. Output matrix has been zeroed. \n");
            break;
    }
    
    return SAF_ERROR__NO_ERROR;
}
#endif /* NDEBUG */
