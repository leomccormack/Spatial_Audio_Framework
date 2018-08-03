/*
 Copyright 2016-2018 Leo McCormack
 
 Permission to use, copy, modify, and/or distribute this software for any purpose with or
 without fee is hereby granted, provided that the above copyright notice and this permission
 notice appear in all copies.
 
 THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO
 THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT
 SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR
 ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
 CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
 OR PERFORMANCE OF THIS SOFTWARE.
*/ 
/*
 * Filename:
 *     saf_utilities.h (include header)
 * Description:
 *     Contains a collection of useful memory allocation functions and cross-platform
 *     complex number wrappers. Optimised linear algebra routines utilising BLAS and LAPACK
 *     are also included.
 * Dependencies:
 *     Windows users only: Intel's MKL must be installed, which can be freely aquired via:
 *     https://software.intel.com/en-us/articles/free-ipsxe-tools-and-libraries
 *     Mac users only: saf_utilities will utilise Apple's Accelerate library.
 * Author, date created:
 *     Leo McCormack, 11.07.2016
 */

#ifndef __SAF_UTILITIES_H_INCLUDED__
#define __SAF_UTILITIES_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif

#ifndef MIN
  #define MIN(a,b) (( (a) < (b) ) ? (a) : (b))
#endif
#ifndef MAX
  #define MAX(a,b) (( (a) > (b) ) ? (a) : (b))
#endif
#ifndef PI
  #define PI ( 3.14159265358979323846264338327950288f )
#endif
#ifndef M_PI
  #define M_PI ( 3.14159265358979323846264338327950288f )
#endif

/* For multi-dimensional memory handling */
#include "../saf_utilities/saf_malloc.h"
#include "../saf_utilities/saf_calloc.h"
#include "../saf_utilities/saf_free.h"

/* For sorting vectors */
#include "../saf_utilities/saf_sort.h"

/* For BLAS/LAPACK functions, plus some other handy linear algebra functions */
#include "../saf_utilities/saf_veclib.h"

/* For cross-platform complex numbers wrapper */
#include "../saf_utilities/saf_complex.h"

/* For various presets for loudspeaker, microphone, and hydrophone arrays.  */
#include "../saf_utilities/saf_loudspeaker_presets.h"
#include "../saf_utilities/saf_sensorarray_presets.h"

#ifdef __cplusplus
}
#endif

#endif /* __SAF_UTILITIES_H_INCLUDED__ */
