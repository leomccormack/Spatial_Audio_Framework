/*
 * Copyright 2016-2018 Leo McCormack
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
 * @file saf_utility_veclib.c
 * @ingroup Utilities
 * @brief Wrappers for optimised linear algebra routines, utilising CBLAS and
 *        LAPACK, and/or SIMD intrinsics
 *
 * ## Dependencies
 *   A performance library comprising CBLAS and LAPACK routines is required by
 *   the module and, thus, also by SAF as a whole. Add one of the following
 *   FLAGS to your project's preprocessor definitions list in order to enable
 *   one of these suitable performance libraries, which must also be correctly
 *   linked to your project.
 *   - SAF_USE_INTEL_MKL_LP64:
 *       to enable Intel's Math Kernel Library with the Fortran LAPACK interface
 *   - SAF_USE_INTEL_MKL_ILP64
 *       same as SAF_USE_INTEL_MKL except using int64 and LAPACKE interface
 *   - SAF_USE_OPENBLAS_WITH_LAPACKE:
 *       to enable OpenBLAS with the LAPACKE interface
 *   - SAF_USE_APPLE_ACCELERATE:
 *       to enable the Accelerate framework with the Fortran LAPACK interface
 *   - SAF_USE_ATLAS:
 *       to enable ATLAS BLAS routines and ATLAS's CLAPACK interface
 *
 * @see More information can be found here:
 *      https://github.com/leomccormack/Spatial_Audio_Framework/docs
 *
 * @author Leo McCormack
 * @date 11.07.2016
 * @license ISC
 */

#include "saf_utilities.h"
#include "saf_externals.h"

/* Specify which LAPACK interface should be used: */
#if defined(SAF_USE_INTEL_MKL_LP64)
/* Note that Intel MKL LP64 supports Fortran LAPACK and LAPACKE interfaces: */
# define SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE /**< LAPACK interface */
#elif defined(SAF_USE_INTEL_MKL_ILP64)
/* Note that Intel MKL ILP64 will only work with the LAPACKE interface: */
# define SAF_VECLIB_USE_LAPACKE_INTERFACE        /**< LAPACK interface */
#elif defined(SAF_USE_OPEN_BLAS_AND_LAPACKE)
# define SAF_VECLIB_USE_LAPACKE_INTERFACE        /**< LAPACK interface */
#elif defined(SAF_USE_ATLAS)
# define SAF_VECLIB_USE_CLAPACK_INTERFACE        /**< LAPACK interface */
#elif defined(__APPLE__) && defined(SAF_USE_APPLE_ACCELERATE)
# define SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE /**< LAPACK interface */
#elif defined(SAF_USE_GSL)
# define SAF_VECLIB_USE_GSL_LINALG /**< No LAPACK interface, use alternatives */
#else
# error No LAPACK interface was specified!
#endif

/* Additional flags for Intel MKL */
#if defined(SAF_USE_INTEL_MKL_LP64) || defined(SAF_USE_INTEL_MKL_ILP64)
# ifndef SAF_INTEL_MKL_VML_MODE
/**
 * Taken from the Intel MKL VML developers reference documentation:
 * The flag(s) that can be passed to Intel VML functions:
 *  - Faster processing of denormalized inputs: VML_FTZDAZ_ON || VML_FTZDAZ_OFF
 *  - High Accuracy (VML_HA), the default mode
 *  - Low Accuracy (VML_LA), which improves performance by reducing accuracy of
 *    the two least significant bits
 *  - Enhanced Performance (VML_EP), which provides better performance at the
 *    cost of significantly reduced accuracy. Approximately half of the bits in
 *    the mantissa are correct
 *
 * Note that using the EP mode does not guarantee accurate processing of corner
 * cases and special values. Although the default accuracy is HA, LA is
 * sufficient in most cases. For applications that require less accuracy (for
 * example, media applications, some Monte Carlo simulations, etc.), the EP mode
 * may be sufficient.
 *
 * Note that this default SAF_INTEL_MKL_VML_MODE value can be overriden.
 */
#  define SAF_INTEL_MKL_VML_MODE ( VML_LA | VML_FTZDAZ_ON )
# endif
#endif

/* These are mainly just to remove compiler warnings: */
#if defined(__APPLE__) && defined(SAF_USE_APPLE_ACCELERATE)
  typedef __CLPK_integer        veclib_int;            /**< integer: 4-bytes */
  typedef __CLPK_real           veclib_float;          /**< real: 4-bytes */
  typedef __CLPK_doublereal     veclib_double;         /**< real: 8-bytes */
  typedef __CLPK_complex        veclib_float_complex;  /**< complex: 8-bytes */
  typedef __CLPK_doublecomplex  veclib_double_complex; /**< complex: 16-bytes */
#elif defined(SAF_USE_INTEL_MKL_LP64) || defined(SAF_USE_INTEL_MKL_ILP64)
# ifdef SAF_USE_INTEL_MKL_LP64
  typedef MKL_INT               veclib_int;            /**< integer: 4-bytes */
# else /* SAF_USE_INTEL_MKL_ILP64: */
  typedef MKL_INT               veclib_int;            /**< integer: 8-bytes */
# endif
  typedef float                 veclib_float;          /**< real: 4-bytes */
  typedef double                veclib_double;         /**< real: 8-bytes */
  typedef MKL_Complex8          veclib_float_complex;  /**< complex: 8-bytes */
  typedef MKL_Complex16         veclib_double_complex; /**< complex: 16-bytes */
#elif defined(SAF_USE_OPEN_BLAS_AND_LAPACKE)
  typedef lapack_int            veclib_int;            /**< integer: 4-bytes */
  typedef float                 veclib_float;          /**< real: 4-bytes */
  typedef double                veclib_double;         /**< real: 8-bytes */
  typedef lapack_complex_float  veclib_float_complex;  /**< complex: 8-bytes */
  typedef lapack_complex_double veclib_double_complex; /**< complex: 16-bytes */
#else
  typedef int                   veclib_int;            /**< integer: 4-bytes */
  typedef float                 veclib_float;          /**< real: 4-bytes */
  typedef double                veclib_double;         /**< real: 8-bytes */
  typedef float_complex         veclib_float_complex;  /**< complex: 8-bytes */
  typedef double_complex        veclib_double_complex; /**< complex: 16-bytes */
#endif


/* ========================================================================== */
/*                     Built-in CBLAS Functions (Level 0)                     */
/* ========================================================================== */

#ifdef SAF_USE_BUILT_IN_NAIVE_CBLAS
void cblas_scopy(const int N, const float *X, const int incX, float *Y, const int incY){
    int i;
#if defined(SAF_ENABLE_SIMD)
    if(incX==1 && incY==1){
        for(i=0; i<(N-3); i+=4)
            _mm_storeu_ps(Y+i, _mm_loadu_ps(X+i));
        for(;i<N; i++) /* The residual (if N was not divisable by 4): */
            Y[i] = X[i];
    }
    else
        for(i=0; i<N; i++)
            Y[i*incY] = X[i*incX];
#else
    for(i=0; i<N; i++)
        Y[i*incY] = X[i*incX];
#endif
}

void cblas_dcopy(const int N, const double *X, const int incX, double *Y, const int incY){
    int i,j;
    for(i=j=0; i<N; i+=incX, j+=incY)
        Y[i] = X[j];
}

void cblas_ccopy(const int N, const void *X, const int incX, void *Y, const int incY){
    int i,j;
    float_complex *cX, *cY;
    cX = (float_complex*)X;
    cY = (float_complex*)Y;
    for(i=j=0; i<N; i+=incX, j+=incY)
        cY[i] = cX[j];
}

void cblas_zcopy(const int N, const void *X, const int incX, void *Y, const int incY){
    int i,j;
    double_complex *cX, *cY;
    cX = (double_complex*)X;
    cY = (double_complex*)Y;
    for(i=j=0; i<N; i+=incX, j+=incY)
        cY[i] = cX[j];
}

void cblas_saxpy(const int N, const float alpha, const float* X, const int incX, float* Y, const int incY) {
    int i,j;
    for (i=j=0; i<N; i+=incX, j+=incY)
        Y[i] = alpha * X[i] + Y[i];
}

void cblas_daxpy(const int N, const double alpha, const double* X, const int incX, double* Y, const int incY) {
    int i,j;
    for (i=j=0; i<N; i+=incX, j+=incY)
        Y[i] = alpha * X[i] + Y[i];
}

float cblas_sdot(const int N, const float* X, const int incX, const float* Y, const int incY){
    int i;
    float ret;
#if defined(SAF_ENABLE_SIMD)
    float sum2;
    if(incX==1 && incY==1){
        sum2 = 0.0f;
        __m128 sum = _mm_setzero_ps();
        /* Sum over all elements in groups of 4 */
        for(i=0; i<(N-3); i+=4)
           sum = _mm_add_ps(sum, _mm_mul_ps(_mm_loadu_ps(X+i), _mm_loadu_ps(Y+i)));
        /* sum the first 2 elements with the second 2 elements and store in the first 2 elements */
        sum = _mm_add_ps(sum, _mm_movehl_ps(sum, sum)/* copy elements 3 and 4 to 1 and 2 */);
        /* sum the first 2 elements and store in the first element */
        sum = _mm_add_ss(sum, _mm_shuffle_ps(sum, sum, _MM_SHUFFLE(3, 2, 0, 1))/* shuffle so that element 2 goes to element 1 */);
        _mm_store_ss(&ret, sum); /* get the first element (the result) */
        for(;i<N; i++) /* The residual (if N was not divisable by 4): */
            sum2 += X[i] * Y[i];
        ret += sum2;
    }
    else
        for(i=0; i<N; i++)
            ret += X[i*incX] * Y[i*incY];
#else
    ret = 0.0f;
    for(i=0; i<N; i++)
        ret += X[i*incX] * Y[i*incY];
#endif
    return ret;
}

#endif


/* ========================================================================== */
/*                     Built-in CBLAS Functions (Level 3)                     */
/* ========================================================================== */

#ifdef SAF_USE_BUILT_IN_NAIVE_CBLAS
void cblas_sgemm(const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE TransA,
                     const CBLAS_TRANSPOSE TransB, const int M, const int N,
                     const int K, const float alpha, const float *A,
                     const int lda, const float *B, const int ldb,
                     const float beta, float *C, const int ldc)
{
    saf_assert(0, "INCOMPLETE and UNTESTED");
    int i,j,k;
    int _transA;
    int _transB;

    switch(Layout){
        case CblasRowMajor:
            _transA = (TransA == CblasNoTrans) ? SAF_FALSE : SAF_TRUE;
            _transB = (TransB == CblasNoTrans) ? SAF_FALSE : SAF_TRUE;
            break;
        case CblasColMajor:
            _transA = (TransA == CblasNoTrans) ? SAF_TRUE : SAF_FALSE;
            _transB = (TransB == CblasNoTrans) ? SAF_TRUE : SAF_FALSE;
            break;
    }

    if(!_transA && !_transB){
        for(i=0; i<M; i++){
            for(j=0; j<N; j++){
                C[i*ldc+j] = 0.0f;
                for(k=0; k<K; k++)
                    C[i*ldc+j] += A[i*lda+k] * B[k*ldb+j];
            }
        }
    }
    else if(_transA && !_transB){
    }
    else if(_transA && _transB){
    }
}
#endif


/* ========================================================================== */
/*                     Find Index of Min-Abs-Value (?iminv)                   */
/* ========================================================================== */

void utility_siminv
(
    const float* a,
    const int len,
    int* index
)
{
#if defined(SAF_USE_INTEL_IPP)
    float minVal;
    ippsMinAbsIndx_32f((Ipp32f*)a, len, &minVal, index);
#elif defined(SAF_USE_APPLE_ACCELERATE)
    float minVal;
    vDSP_Length ind_tmp;
    vDSP_minmgvi(a, 1, &minVal, &ind_tmp, (vDSP_Length)len);
    *index = (int)ind_tmp; 
#elif defined(SAF_USE_INTEL_MKL_LP64) || defined(SAF_USE_INTEL_MKL_ILP64)
    *index = (int)cblas_isamin(len, a, 1);
#else
    int i;
    *index = 0;
    float minVal=FLT_MAX;
    for(i=0; i<len; i++){
        if(fabsf(a[i])<minVal){
            minVal = fabsf(a[i]);
            *index = i;
        } 
    }
#endif
}

void utility_ciminv
(
    const float_complex* a,
    const int len,
    int* index
)
{
#if defined(SAF_USE_APPLE_ACCELERATE) /* Unfortunately requires a malloc call */
    float* abs_a;
    float minVal;
    abs_a = malloc1d(len*sizeof(float));
    vDSP_vdist((float*)a/*real*/, 2, (float*)a+1/*imag*/, 2, abs_a, 1, (vDSP_Length)len); /* cabsf */
    vDSP_Length ind_tmp;
    vDSP_minmgvi(abs_a, 1, &minVal, &ind_tmp, (vDSP_Length)len);
    *index = (int)ind_tmp;
    free(abs_a);
#elif defined(SAF_USE_INTEL_MKL_LP64) || defined(SAF_USE_INTEL_MKL_ILP64)
    *index = (int)cblas_icamin(len, a, 1);
#else
    int i;
    *index = 0;
    float minVal=FLT_MAX;
    for(i=0; i<len; i++){
        if(cabsf(a[i])<minVal){
            minVal = cabsf(a[i]);
            *index = i;
        } 
    }
#endif
}

void utility_diminv
(
    const double* a,
    const int len,
    int* index
)
{
#if defined(SAF_USE_INTEL_IPP)
    double minVal;
    ippsMinAbsIndx_64f((Ipp64f*)a, len, &minVal, index);
#elif defined(SAF_USE_APPLE_ACCELERATE)
    double minVal;
    vDSP_Length ind_tmp;
    vDSP_minmgviD(a, 1, &minVal, &ind_tmp, (vDSP_Length)len);
    *index = (int)ind_tmp;
#elif defined(SAF_USE_INTEL_MKL_LP64) || defined(SAF_USE_INTEL_MKL_ILP64)
    *index = (int)cblas_idamin(len, a, 1);
#else
    int i;
    *index = 0;
    double minVal=DBL_MAX;
    for(i=0; i<len; i++){
        if(fabs(a[i])<minVal){
            minVal = fabs(a[i]);
            *index = i;
        }
    }
#endif
}

void utility_ziminv
(
    const double_complex* a,
    const int len,
    int* index
)
{
#if defined(SAF_USE_APPLE_ACCELERATE) /* Unfortunately requires a malloc call */
    double* abs_a;
    double minVal;
    abs_a = malloc1d(len*sizeof(double));
    vDSP_vdistD((double*)a/*real*/, 2, (double*)a+1/*imag*/, 2, abs_a, 1, (vDSP_Length)len); /* cabs */
    vDSP_Length ind_tmp;
    vDSP_minmgviD(abs_a, 1, &minVal, &ind_tmp, (vDSP_Length)len);
    *index = (int)ind_tmp;
    free(abs_a);
#elif defined(SAF_USE_INTEL_MKL_LP64) || defined(SAF_USE_INTEL_MKL_ILP64)
    *index = (int)cblas_izamin(len, a, 1);
#else
    int i;
    *index = 0;
    double minVal=DBL_MAX;
    for(i=0; i<len; i++){
        if(cabs(a[i])<minVal){
            minVal = cabs(a[i]);
            *index = i;
        }
    }
#endif
}


/* ========================================================================== */
/*                     Find Index of Max-Abs-Value (?imaxv)                   */
/* ========================================================================== */

void utility_simaxv
(
    const float* a,
    const int len,
    int* index
)
{
#if defined(SAF_USE_INTEL_IPP)
    float maxVal;
    ippsMaxAbsIndx_32f((Ipp32f*)a, len, &maxVal, index);
#else
    *index = (int)cblas_isamax(len, a, 1);
#endif
}

void utility_cimaxv
(
    const float_complex* a,
    const int len,
    int* index
)
{
    *index = (int)cblas_icamax(len, a, 1);
}

void utility_dimaxv
(
    const double* a,
    const int len,
    int* index
)
{
#if defined(SAF_USE_INTEL_IPP)
    double maxVal;
    ippsMaxAbsIndx_64f((Ipp64f*)a, len, &maxVal, index);
#else
    *index = (int)cblas_idamax(len, a, 1);
#endif
}

void utility_zimaxv
(
    const double_complex* a,
    const int len,
    int* index
)
{
    *index = (int)cblas_izamax(len, a, 1);
}


/* ========================================================================== */
/*                              Vector-Abs (?vabs)                            */
/* ========================================================================== */

void utility_svabs
(
    const float* a,
    const int len,
    float* c
)
{
#if defined(SAF_USE_INTEL_IPP)
    ippsAbs_32f((Ipp32f*)a, (Ipp32f*)c, len);
#elif defined(SAF_USE_APPLE_ACCELERATE)
    vDSP_vabs(a, 1, c, 1, (vDSP_Length)len);
#elif defined(SAF_USE_INTEL_MKL_LP64) || defined(SAF_USE_INTEL_MKL_ILP64)
    vmsAbs(len, a, c, SAF_INTEL_MKL_VML_MODE);
#else
    int i;
    for(i=0; i<len; i++)
        c[i] = fabsf(a[i]);
#endif
}

void utility_cvabs
(
    const float_complex* a,
    const int len,
    float* c
)
{
#if defined(SAF_USE_INTEL_IPP)
    ippsMagnitude_32fc((Ipp32fc*)a, (Ipp32f*)c, len);
#elif defined(SAF_USE_APPLE_ACCELERATE)
    vDSP_vdist((float*)a/*real*/, 2, (float*)a+1/*imag*/, 2, c, 1, (vDSP_Length)len);
#elif defined(SAF_USE_INTEL_MKL_LP64) || defined(SAF_USE_INTEL_MKL_ILP64)
    vmcAbs(len, (MKL_Complex8*)a, c, SAF_INTEL_MKL_VML_MODE);
#else
    int i;
    for(i=0; i<len; i++)
        c[i] = cabsf(a[i]);
#endif
}


/* ========================================================================== */
/*                            Vector-Modulus (?vmod)                          */
/* ========================================================================== */

void utility_svmod
(
    const float* a,
    const float* b,
    const int len,
    float* c
)
{
#if defined(SAF_USE_APPLE_ACCELERATE)
    vvfmodf(c, a, b, &len);
#elif defined(SAF_USE_INTEL_MKL_LP64) || defined(SAF_USE_INTEL_MKL_ILP64)
    vmsFmod(len, a, b, c, SAF_INTEL_MKL_VML_MODE);
#else
    int i;
    for(i=0; i<len; i++)
        c[i] = fmodf(a[i], b[i]);
#endif
}


/* ========================================================================== */
/*                          Vector-Reciprocal (?vrecip)                       */
/* ========================================================================== */

void utility_svrecip
(
    const float* a,
    const int len,
    float* c
)
{
#if defined(SAF_USE_INTEL_IPP)
    ippsDivCRev_32f((Ipp32f*)a, 1.0f, (Ipp32f*)c, len);
#elif defined(SAF_USE_APPLE_ACCELERATE)
    float one;
    one = 1.0f;
    vDSP_svdiv(&one, a, 1, c, 1, (vDSP_Length)len);
#elif defined(SAF_USE_INTEL_MKL_LP64) || defined(SAF_USE_INTEL_MKL_ILP64)
    vmsInv(len, a, c, SAF_INTEL_MKL_VML_MODE);
#elif defined(SAF_ENABLE_SIMD)
    int i;
    i = 0;
# if defined(__AVX512F__)
    for(; i<(len-15); i+=16)
        _mm512_storeu_ps(c+i, _mm512_rcp14_ps(_mm512_loadu_ps(a+i)));
# endif
# if defined(__AVX__) && defined(__AVX2__)
    for(; i<(len-7); i+=8)
        _mm256_storeu_ps(c+i, _mm256_rcp_ps(_mm256_loadu_ps(a+i)));
# endif
# if defined(__SSE__) && defined(__SSE2__) && defined(__SSE3__)
    for(; i<(len-3); i+=4)
        _mm_storeu_ps(c+i, _mm_rcp_ps(_mm_loadu_ps(a+i)));
# endif
    for(;i<len; i++) /* The residual (if len was not divisable by the step size): */
        c[i] = 1.0f/a[i];
#else
    int i;
    for(i=0; i<len; i++)
        c[i] = 1.0f/a[i];
#endif
}


/* ========================================================================== */
/*                           Vector-Conjugate (?vconj)                        */
/* ========================================================================== */

void utility_cvconj
(
    const float_complex* a,
    const int len,
    float_complex* c
)
{
#if defined(SAF_USE_INTEL_IPP)
    ippsConj_32fc((Ipp32fc*)a, (Ipp32fc*)c, len);
#else
    cblas_ccopy(len, a, 1, c, 1);
    cblas_sscal(len, -1.0f, ((float*)c)+1, 2);
#endif
}

void utility_zvconj
(
    const double_complex* a,
    const int len,
    double_complex* c
)
{
#if defined(SAF_USE_INTEL_IPP)
    ippsConj_64fc((Ipp64fc*)a, (Ipp64fc*)c, len);
#else
    cblas_zcopy(len, a, 1, c, 1);
    cblas_dscal(len, -1.0, ((double*)c)+1, 2);
#endif
}


/* ========================================================================== */
/*                        Vector-Vector Copy (?vvcopy)                        */
/* ========================================================================== */

void utility_svvcopy
(
    const float* a,
    const int len,
    float* c
)
{
    cblas_scopy(len, a, 1, c, 1);
}

void utility_cvvcopy
(
    const float_complex* a,
    const int len,
    float_complex* c
)
{
    cblas_ccopy(len, a, 1, c, 1);
}

void utility_dvvcopy
(
    const double* a,
    const int len,
    double* c
)
{
    cblas_dcopy(len, a, 1, c, 1);
}

void utility_zvvcopy
(
    const double_complex* a,
    const int len,
    double_complex* c
)
{
    cblas_zcopy(len, a, 1, c, 1);
}


/* ========================================================================== */
/*                       Vector-Vector Addition (?vvadd)                      */
/* ========================================================================== */

void utility_svvadd
(
    const float* a,
    const float* b,
    const int len,
    float* c
)
{
    /* Checks: */
#ifndef NDEBUG
    saf_assert(c!=NULL, "'c' can no longer be NULL");
    saf_assert(a!=c && b!=c, "In-place operation is no longer supported, use e.g. cblas_saxby() instead.");
#endif

    /* The operation: */
#if defined(SAF_USE_INTEL_IPP)
    ippsAdd_32f((Ipp32f*)a, (Ipp32f*)b, (Ipp32f*)c, len);
#elif defined(SAF_USE_APPLE_ACCELERATE)
    vDSP_vadd(a, 1, b, 1, c, 1, (vDSP_Length)len);
#elif defined(SAF_USE_INTEL_MKL_LP64) || defined(SAF_USE_INTEL_MKL_ILP64)
    vmsAdd(len, a, b, c, SAF_INTEL_MKL_VML_MODE);
#elif defined(SAF_ENABLE_SIMD)
    int i;
    i = 0;
# if defined(__AVX512F__)
    for(; i<(len-15); i+=16)
        _mm512_storeu_ps(c+i, _mm512_add_ps(_mm512_loadu_ps(a+i), _mm512_loadu_ps(b+i)));
# endif
# if defined(__AVX__) && defined(__AVX2__)
    for(; i<(len-7); i+=8)
        _mm256_storeu_ps(c+i, _mm256_add_ps(_mm256_loadu_ps(a+i), _mm256_loadu_ps(b+i)));
# endif
# if defined(__SSE__) && defined(__SSE2__) && defined(__SSE3__)
    for(; i<(len-3); i+=4)
        _mm_storeu_ps(c+i, _mm_add_ps(_mm_loadu_ps(a+i), _mm_loadu_ps(b+i)));
# endif
    for(; i<len; i++) /* The residual (if len was not divisable by the step size): */
        c[i] = a[i] + b[i];
#elif defined(NDEBUG)
    int i;
    /* try to indirectly "trigger" some compiler optimisations */
    for(i=0; i<len-3; i+=4){
        c[i] = a[i] + b[i];
        c[i+1] = a[i+1] + b[i+1];
        c[i+2] = a[i+2] + b[i+2];
        c[i+3] = a[i+3] + b[i+3];
    }
    for(; i<len; i++) /* The residual (if len was not divisable by the step size): */
        c[i] = a[i] + b[i];
#else
    int j;
    for (j = 0; j < len; j++)
        c[j] = a[j] + b[j];
#endif
}

void utility_cvvadd
(
    const float_complex* a,
    const float_complex* b,
    const int len,
    float_complex* c
)
{
    /* Checks: */
#ifndef NDEBUG
    saf_assert(c!=NULL, "'c' can no longer be NULL");
    saf_assert(a!=c && b!=c, "In-place operation is no longer supported, use e.g. cblas_caxby() instead.");
#endif

    /* The operation: */
#if defined(SAF_USE_INTEL_IPP)
    ippsAdd_32fc((Ipp32fc*)a, (Ipp32fc*)b, (Ipp32fc*)c, len);
#elif defined(SAF_USE_APPLE_ACCELERATE)
    vDSP_vadd((float*)a, 1, (float*)b, 1, (float*)c, 1, /*re+im*/2*(vDSP_Length)len);
#elif defined(SAF_USE_INTEL_MKL_LP64) || defined(SAF_USE_INTEL_MKL_ILP64)
    vmcAdd(len, (MKL_Complex8*)a, (MKL_Complex8*)b, (MKL_Complex8*)c, SAF_INTEL_MKL_VML_MODE);
#elif defined(SAF_ENABLE_SIMD)
    int i, len2;
    float* sa, *sb, *sc;
    len2 = len*2;
    sa = (float*)a; sb = (float*)b; sc = (float*)c;
    i = 0;
# if defined(__AVX512F__)
    for(; i<(len2-15); i+=16)
        _mm512_storeu_ps(sc+i, _mm512_add_ps(_mm512_loadu_ps(sa+i), _mm512_loadu_ps(sb+i)));
# endif
# if defined(__AVX__) && defined(__AVX2__)
    for(; i<(len2-7); i+=8)
        _mm256_storeu_ps(sc+i, _mm256_add_ps(_mm256_loadu_ps(sa+i), _mm256_loadu_ps(sb+i)));
# endif
# if defined(__SSE__) && defined(__SSE2__) && defined(__SSE3__)
    for(; i<(len2-3); i+=4)
        _mm_storeu_ps(sc+i, _mm_add_ps(_mm_loadu_ps(sa+i), _mm_loadu_ps(sb+i)));
# endif
    for(; i<len2; i++) /* The residual (if len2 was not divisable by the step size): */
        sc[i] = sa[i] + sb[i];
#elif __STDC_VERSION__ >= 199901L && defined(NDEBUG)
    int i;
    /* try to indirectly "trigger" some compiler optimisations */
    for(i=0; i<len-3; i+=4){
        c[i] = a[i] + b[i];
        c[i+1] = a[i+1] + b[i+1];
        c[i+2] = a[i+2] + b[i+2];
        c[i+3] = a[i+3] + b[i+3];
    }
    for(; i<len; i++) /* The residual (if len was not divisable by the step size): */
        c[i] = a[i] + b[i];
#else
    int j;
    for (j = 0; j < len; j++)
        c[j] = ccaddf(a[j], b[j]);
#endif
}

void utility_dvvadd
(
    const double* a,
    const double* b,
    const int len,
    double* c
)
{
    /* Checks: */
#ifndef NDEBUG
    saf_assert(c!=NULL, "'c' can no longer be NULL");
    saf_assert(a!=c && b!=c, "In-place operation is no longer supported, use e.g. cblas_daxby() instead.");
#endif

    /* The operation: */
#if defined(SAF_USE_INTEL_IPP)
    ippsAdd_64f((Ipp64f*)a, (Ipp64f*)b, (Ipp64f*)c, len);
#elif defined(SAF_USE_APPLE_ACCELERATE)
    vDSP_vaddD(a, 1, b, 1, c, 1, (vDSP_Length)len);
#elif defined(SAF_USE_INTEL_MKL_LP64) || defined(SAF_USE_INTEL_MKL_ILP64)
    vmdAdd(len, a, b, c, SAF_INTEL_MKL_VML_MODE);
#elif defined(SAF_ENABLE_SIMD)
    int i;
    i = 0;
# if defined(__AVX512F__)
    for(; i<(len-7); i+=8)
        _mm512_storeu_pd(c+i, _mm512_add_pd(_mm512_loadu_pd(a+i), _mm512_loadu_pd(b+i)));
# endif
# if defined(__AVX__) && defined(__AVX2__)
    for(; i<(len-3); i+=4)
        _mm256_storeu_pd(c+i, _mm256_add_pd(_mm256_loadu_pd(a+i), _mm256_loadu_pd(b+i)));
# endif
# if defined(__SSE__) && defined(__SSE2__) && defined(__SSE3__)
    for(; i<(len-1); i+=2)
        _mm_storeu_pd(c+i, _mm_add_pd(_mm_loadu_pd(a+i), _mm_loadu_pd(b+i)));
# endif
    for(; i<len; i++) /* The residual (if len was not divisable by the step size): */
        c[i] = a[i] + b[i];
#else
    int j;
    for (j = 0; j < len; j++)
        c[j] = a[j] + b[j];
#endif
}

void utility_zvvadd
(
    const double_complex* a,
    const double_complex* b,
    const int len,
    double_complex* c
)
{
    /* Checks: */
#ifndef NDEBUG
    saf_assert(c!=NULL, "'c' can no longer be NULL");
    saf_assert(a!=c && b!=c, "In-place operation is no longer supported, use e.g. cblas_zaxby() instead.");
#endif

    /* The operation: */
#if defined(SAF_USE_INTEL_IPP)
    ippsAdd_64fc((Ipp64fc*)a, (Ipp64fc*)b, (Ipp64fc*)c, len);
#elif defined(SAF_USE_APPLE_ACCELERATE)
    vDSP_vaddD((double*)a, 1, (double*)b, 1, (double*)c, 1, /*re+im*/2*(vDSP_Length)len);
#elif defined(SAF_USE_INTEL_MKL_LP64) || defined(SAF_USE_INTEL_MKL_ILP64)
    vmzAdd(len, (MKL_Complex16*)a, (MKL_Complex16*)b, (MKL_Complex16*)c, SAF_INTEL_MKL_VML_MODE);
#elif defined(SAF_ENABLE_SIMD)
    int i, len2;
    double* sa, *sb, *sc;
    len2 = len*2;
    sa = (double*)a; sb = (double*)b; sc = (double*)c;
    i = 0;
# if defined(__AVX512F__)
    for(; i<(len2-7); i+=8)
        _mm512_storeu_pd(sc+i, _mm512_add_pd(_mm512_loadu_pd(sa+i), _mm512_loadu_pd(sb+i)));
# endif
# if defined(__AVX__) && defined(__AVX2__)
    for(; i<(len2-3); i+=4)
        _mm256_storeu_pd(sc+i, _mm256_add_pd(_mm256_loadu_pd(sa+i), _mm256_loadu_pd(sb+i)));
# endif
# if defined(__SSE__) && defined(__SSE2__) && defined(__SSE3__)
    for(; i<(len2-1); i+=2)
        _mm_storeu_pd(sc+i, _mm_add_pd(_mm_loadu_pd(sa+i), _mm_loadu_pd(sb+i)));
# endif
    for(; i<len2; i++) /* The residual (if len2 was not divisable by the step size): */
        sc[i] = sa[i] + sb[i];
#else
    int j;
    for (j = 0; j < len; j++)
        c[j] = ccadd(a[j], b[j]);
#endif
}


/* ========================================================================== */
/*                     Vector-Vector Subtraction (?vvsub)                     */
/* ========================================================================== */

void utility_svvsub
(
    const float* a,
    const float* b,
    const int len,
    float* c
)
{
    /* Checks: */
#ifndef NDEBUG
    saf_assert(c!=NULL, "'c' cannot be NULL");
    saf_assert(a!=c && b!=c, "In-place operation is not supported.");
#endif

    /* The operation: */
#if defined(SAF_USE_INTEL_IPP)
    ippsSub_32f((Ipp32f*)b, (Ipp32f*)a, (Ipp32f*)c, len); /* 'a' and 'b' are switched */
#elif defined(SAF_USE_APPLE_ACCELERATE)
    vDSP_vsub(b, 1, a, 1, c, 1, (vDSP_Length)len);        /* 'a' and 'b' are switched */
#elif defined(SAF_USE_INTEL_MKL_LP64) || defined(SAF_USE_INTEL_MKL_ILP64)
    vmsSub(len, a, b, c, SAF_INTEL_MKL_VML_MODE);
#elif defined(SAF_ENABLE_SIMD)
    int i;
    i = 0;
# if defined(__AVX512F__)
    for(; i<(len-15); i+=16)
        _mm512_storeu_ps(c+i, _mm512_sub_ps(_mm512_loadu_ps(a+i), _mm512_loadu_ps(b+i)));
# endif
# if defined(__AVX__) && defined(__AVX2__)
    for(; i<(len-7); i+=8)
        _mm256_storeu_ps(c+i, _mm256_sub_ps(_mm256_loadu_ps(a+i), _mm256_loadu_ps(b+i)));
# endif
# if defined(__SSE__) && defined(__SSE2__) && defined(__SSE3__)
    for(; i<(len-3); i+=4)
        _mm_storeu_ps(c+i, _mm_sub_ps(_mm_loadu_ps(a+i), _mm_loadu_ps(b+i)));
# endif
    for(; i<len; i++) /* The residual (if len was not divisable by the step size): */
        c[i] = a[i] - b[i];
#elif defined(NDEBUG)
    int i;
    /* try to indirectly "trigger" some compiler optimisations */
    for(i=0; i<len-3; i+=4){
        c[i] = a[i] - b[i];
        c[i+1] = a[i+1] - b[i+1];
        c[i+2] = a[i+2] - b[i+2];
        c[i+3] = a[i+3] - b[i+3];
    }
    for(; i<len; i++) /* The residual (if len was not divisable by the step size): */
        c[i] = a[i] - b[i];
#else
    int j;
    for (j = 0; j < len; j++)
        c[j] = a[j] - b[j];
#endif
}

void utility_cvvsub
(
    const float_complex* a,
    const float_complex* b,
    const int len,
    float_complex* c
)
{
    /* Checks: */
#ifndef NDEBUG
    saf_assert(c!=NULL, "'c' cannot be NULL");
    saf_assert(a!=c && b!=c, "In-place operation is not supported.");
#endif

    /* The operation: */
#if defined(SAF_USE_INTEL_IPP)
    ippsSub_32fc((Ipp32fc*)b, (Ipp32fc*)a, (Ipp32fc*)c, len);                         /* 'a' and 'b' are switched */
#elif defined(SAF_USE_APPLE_ACCELERATE)
    vDSP_vsub((float*)b, 1, (float*)a, 1, (float*)c, 1, /*re+im*/2*(vDSP_Length)len); /* 'a' and 'b' are switched */
#elif defined(SAF_USE_INTEL_MKL_LP64) || defined(SAF_USE_INTEL_MKL_ILP64)
    vmcSub(len, (MKL_Complex8*)a, (MKL_Complex8*)b, (MKL_Complex8*)c, SAF_INTEL_MKL_VML_MODE);
#elif defined(SAF_ENABLE_SIMD)
    int i, len2;
    float* sa, *sb, *sc;
    len2 = len*2;
    sa = (float*)a; sb = (float*)b; sc = (float*)c;
    i = 0;
# if defined(__AVX512F__)
    for(; i<(len2-15); i+=16)
        _mm512_storeu_ps(sc+i, _mm512_sub_ps(_mm512_loadu_ps(sa+i), _mm512_loadu_ps(sb+i)));
# endif
# if defined(__AVX__) && defined(__AVX2__)
    for(; i<(len2-7); i+=8)
        _mm256_storeu_ps(sc+i, _mm256_sub_ps(_mm256_loadu_ps(sa+i), _mm256_loadu_ps(sb+i)));
# endif
# if defined(__SSE__) && defined(__SSE2__) && defined(__SSE3__)
    for(; i<(len2-3); i+=4)
        _mm_storeu_ps(sc+i, _mm_sub_ps(_mm_loadu_ps(sa+i), _mm_loadu_ps(sb+i)));
# endif
    for(; i<len2; i++) /* The residual (if len2 was not divisable by the step size): */
        sc[i] = sa[i] - sb[i];
#elif __STDC_VERSION__ >= 199901L && defined(NDEBUG)
    int i;
    /* try to indirectly "trigger" some compiler optimisations */
    for(i=0; i<len-3; i+=4){
        c[i] = a[i] - b[i];
        c[i+1] = a[i+1] - b[i+1];
        c[i+2] = a[i+2] - b[i+2];
        c[i+3] = a[i+3] - b[i+3];
    }
    for(; i<len; i++) /* The residual (if len was not divisable by the step size): */
        c[i] = a[i] - b[i];
#else
    int j;
    for (j = 0; j < len; j++)
        c[j] = ccsubf(a[j], b[j]);
#endif
}

void utility_dvvsub
(
    const double* a,
    const double* b,
    const int len,
    double* c
)
{
    /* Checks: */
#ifndef NDEBUG
    saf_assert(c!=NULL, "'c' cannot be NULL");
    saf_assert(a!=c && b!=c, "In-place operation is not supported.");
#endif

    /* The operation: */
#if defined(SAF_USE_INTEL_IPP)
    ippsSub_64f((Ipp64f*)b, (Ipp64f*)a, (Ipp64f*)c, len);  /* 'a' and 'b' are switched */
#elif defined(SAF_USE_APPLE_ACCELERATE)
    vDSP_vsubD(b, 1, a, 1, c, 1, (vDSP_Length)len);        /* 'a' and 'b' are switched */
#elif defined(SAF_USE_INTEL_MKL_LP64) || defined(SAF_USE_INTEL_MKL_ILP64)
    vmdSub(len, a, b, c, SAF_INTEL_MKL_VML_MODE);
#elif defined(SAF_ENABLE_SIMD)
    int i;
    i = 0;
# if defined(__AVX512F__)
    for(; i<(len-7); i+=8)
        _mm512_storeu_pd(c+i, _mm512_sub_pd(_mm512_loadu_pd(a+i), _mm512_loadu_pd(b+i)));
# endif
# if defined(__AVX__) && defined(__AVX2__)
    for(; i<(len-3); i+=4)
        _mm256_storeu_pd(c+i, _mm256_sub_pd(_mm256_loadu_pd(a+i), _mm256_loadu_pd(b+i)));
# endif
# if defined(__SSE__) && defined(__SSE2__) && defined(__SSE3__)
    for(; i<(len-1); i+=2)
        _mm_storeu_pd(c+i, _mm_sub_pd(_mm_loadu_pd(a+i), _mm_loadu_pd(b+i)));
# endif
    for(;i<len; i++) /* The residual (if len was not divisable by the step size): */
        c[i] = a[i] - b[i];
#else
    int j;
    for (j = 0; j < len; j++)
        c[j] = a[j] - b[j];
#endif
}

void utility_zvvsub
(
    const double_complex* a,
    const double_complex* b,
    const int len,
    double_complex* c
)
{
    /* Checks: */
#ifndef NDEBUG
    saf_assert(c!=NULL, "'c' cannot be NULL");
    saf_assert(a!=c && b!=c, "In-place operation is not supported.");
#endif

    /* The operation: */
#if defined(SAF_USE_INTEL_IPP)
    ippsSub_64fc((Ipp64fc*)b, (Ipp64fc*)a, (Ipp64fc*)c, len);                             /* 'a' and 'b' are switched */
#elif defined(SAF_USE_APPLE_ACCELERATE)
    vDSP_vsubD((double*)b, 1, (double*)a, 1, (double*)c, 1, /*re+im*/2*(vDSP_Length)len); /* 'a' and 'b' are switched */
#elif defined(SAF_USE_INTEL_MKL_LP64) || defined(SAF_USE_INTEL_MKL_ILP64)
    vmzSub(len, (MKL_Complex16*)a, (MKL_Complex16*)b, (MKL_Complex16*)c, SAF_INTEL_MKL_VML_MODE);
#elif defined(SAF_ENABLE_SIMD)
    int i, len2;
    double* sa, *sb, *sc;
    len2 = len*2;
    sa = (double*)a; sb = (double*)b; sc = (double*)c;
    i = 0;
# if defined(__AVX512F__)
    for(; i<(len2-7); i+=8)
        _mm512_storeu_pd(sc+i, _mm512_sub_pd(_mm512_loadu_pd(sa+i), _mm512_loadu_pd(sb+i)));
# endif
# if defined(__AVX__) && defined(__AVX2__)
    for(; i<(len2-3); i+=4)
        _mm256_storeu_pd(sc+i, _mm256_sub_pd(_mm256_loadu_pd(sa+i), _mm256_loadu_pd(sb+i)));
# endif
# if defined(__SSE__) && defined(__SSE2__) && defined(__SSE3__)
    for(; i<(len2-1); i+=2)
        _mm_storeu_pd(sc+i, _mm_sub_pd(_mm_loadu_pd(sa+i), _mm_loadu_pd(sb+i)));
# endif
    for(;i<len2; i++) /* The residual (if len2 was not divisable by the step size): */
        sc[i] = sa[i] - sb[i];
#else
    int j;
    for (j = 0; j < len; j++)
        c[j] = ccsub(a[j], b[j]);
#endif
}


/* ========================================================================== */
/*                    Vector-Vector Multiplication (?vvmul)                   */
/* ========================================================================== */

void utility_svvmul
(
    const float* a,
    const float* b,
    const int len,
    float* c
)
{
    /* Checks: */
#ifndef NDEBUG
    saf_assert(c!=NULL, "'c' can no longer be NULL");
    saf_assert(a!=c && b!=c, "In-place operation is no longer supported.");
#endif

    /* The operation: */
#if defined(SAF_USE_INTEL_IPP)
    ippsMul_32f((Ipp32f*)a, (Ipp32f*)b, (Ipp32f*)c, len);
#elif defined(SAF_USE_APPLE_ACCELERATE)
    vDSP_vmul(a, 1, b, 1, c, 1, (vDSP_Length)len);
#elif defined(SAF_USE_INTEL_MKL_LP64) || defined(SAF_USE_INTEL_MKL_ILP64)
    vmsMul(len, a, b, c, SAF_INTEL_MKL_VML_MODE);
#elif defined(SAF_ENABLE_SIMD)
    int i;
    i = 0;
# if defined(__AVX512F__)
    for(; i<(len-15); i+=16)
        _mm512_storeu_ps(c+i, _mm512_mul_ps(_mm512_loadu_ps(a+i), _mm512_loadu_ps(b+i)));
# endif
# if defined(__AVX__) && defined(__AVX2__)
    for(; i<(len-7); i+=8)
        _mm256_storeu_ps(c+i, _mm256_mul_ps(_mm256_loadu_ps(a+i), _mm256_loadu_ps(b+i)));
# endif
# if defined(__SSE__) && defined(__SSE2__) && defined(__SSE3__)
    for(; i<(len-3); i+=4)
        _mm_storeu_ps(c+i, _mm_mul_ps(_mm_loadu_ps(a+i), _mm_loadu_ps(b+i)));
# endif
    for(;i<len; i++) /* The residual (if len was not divisable by the step size): */
        c[i] = a[i] * b[i];
#elif defined(NDEBUG)
    int i;
    /* try to indirectly "trigger" some compiler optimisations */
    for(i=0; i<len-3; i+=4){
        c[i] = a[i] * b[i];
        c[i+1] = a[i+1] * b[i+1];
        c[i+2] = a[i+2] * b[i+2];
        c[i+3] = a[i+3] * b[i+3];
    }
    for(; i<len; i++) /* The residual (if len was not divisable by the step size): */
        c[i] = a[i] * b[i];
#else
    int j;
    for (j = 0; j < len; j++)
        c[j] = a[j] * b[j];
#endif
}

void utility_cvvmul
(
    const float_complex* a,
    const float_complex* b,
    const int len,
    float_complex* c
)
{
    /* Checks: */
#ifndef NDEBUG
    saf_assert(c!=NULL, "'c' can no longer be NULL");
    saf_assert(a!=c && b!=c, "In-place operation is no longer supported.");
#endif

    /* The operation: */
#if defined(SAF_USE_INTEL_IPP)
    ippsMul_32fc((Ipp32fc*)a, (Ipp32fc*)b, (Ipp32fc*)c, len);
#elif defined(SAF_USE_APPLE_ACCELERATE)  /* Due to Apple "logic", this is unfortunately quite complicated, and probably slower than it should be... */
    /* Imaginary part of the output */
    vDSP_vmul((float*)a/*real*/, 2, (float*)b+1/*imag*/, 2, (float*)c+1/*imag*/, 2, (vDSP_Length)len);
    vDSP_vma((float*)a+1/*imag*/, 2, (float*)b/*real*/, 2, (float*)c+1/*imag*/, 2, (float*)c/*real*/, 2, (vDSP_Length)len); /* Use the real part of c as a temporary buffer */
    cblas_scopy(len, (float*)c/*real*/, 2, (float*)c+1/*imag*/, 2); /* Copy the imaginary part from the temporary buffer */
    /* Real part of the output */
    vDSP_vmmsb((float*)a/*real*/, 2, (float*)b/*real*/, 2, (float*)a+1/*imag*/, 2, (float*)b+1/*imag*/, 2, (float*)c/*real*/, 2, (vDSP_Length)len);
#elif defined(SAF_USE_INTEL_MKL_LP64) || defined(SAF_USE_INTEL_MKL_ILP64)
    vmcMul(len, (MKL_Complex8*)a, (MKL_Complex8*)b, (MKL_Complex8*)c, SAF_INTEL_MKL_VML_MODE);
#elif defined(SAF_ENABLE_SIMD)
    int i;
    float* sa, *sb, *sc;
    sa = (float*)a; sb = (float*)b; sc = (float*)c;
# if (defined(__AVX__) && defined(__AVX2__)) || defined(__AVX512F__) /* I couldn't figure out an alternative for addsub with AVX-512... */
    __m256i permute_ri = _mm256_set_epi32(6, 7, 4, 5, 2, 3, 0, 1);
    for(i=0; i<(len-3); i+=4){
        /* Load only the real parts of a */
        __m256 src1 = _mm256_moveldup_ps(_mm256_loadu_ps(sa+2*i)/*|a1|b1|a2|b2|a3|b3|a4|b4|*/); /*|a1|a1|a2|a2|a3|a3|a4|a4|*/
        /* Load real+imag parts of b */
        __m256 src2 = _mm256_loadu_ps(sb+2*i); /*|c1|d1|c2|d2|c3|d3|c4|d4|*/
        /* Multiply together */
        __m256 tmp1 = _mm256_mul_ps(src1, src2);
        /* Swap the real+imag parts of b to be imag+real instead: */
        __m256 b1 = _mm256_permutevar8x32_ps(src2, permute_ri);
        /* Load only the imag parts of a */
        src1 = _mm256_movehdup_ps(_mm256_loadu_ps(sa+2*i)/*|a1|b1|a2|b2|a3|b3|a4|b4|*/); /*|b1|b1|b2|b2|b3|b3|b4|b4|*/
        /* Multiply together */
        __m256 tmp2 = _mm256_mul_ps(src1, b1);
        /* Add even indices, subtract odd indices */
        _mm256_storeu_ps(sc+2*i, _mm256_addsub_ps(tmp1, tmp2));
    }
# elif defined(__SSE__) && defined(__SSE2__) && defined(__SSE3__)
    for(i=0; i<(len-1); i+=2){
        /* Load only the real parts of a */
        __m128 src1 = _mm_moveldup_ps(_mm_loadu_ps(sa+2*i)/*|a1|b1|a2|b2|*/); /*|a1|a1|a2|a2|*/
        /* Load real+imag parts of b */
        __m128 src2 = _mm_loadu_ps(sb+2*i); /*|c1|d1|c2|d2|*/
        /* Multiply together */
        __m128 tmp1 = _mm_mul_ps(src1, src2);
        /* Swap the real+imag parts of b to be imag+real instead: */
        __m128 b1 = _mm_shuffle_ps(src2, src2, _MM_SHUFFLE(2, 3, 0, 1));
        /* Load only the imag parts of a */
        src1 = _mm_movehdup_ps(_mm_loadu_ps(sa+2*i)/*|a1|b1|a2|b2|*/); /*|b1|b1|b2|b2|*/
        /* Multiply together */
        __m128 tmp2 = _mm_mul_ps(src1, b1);
        /* Add even indices, subtract odd indices */
        _mm_storeu_ps(sc+2*i, _mm_addsub_ps(tmp1, tmp2));
    }
# endif
    for(;i<len; i++){ /* The residual (if len was not divisable by the step size): */
        sc[2*i]   = sa[2*i] * sb[2*i]   - sa[2*i+1] * sb[2*i+1];
        sc[2*i+1] = sa[2*i] * sb[2*i+1] + sa[2*i+1] * sb[2*i];
    }
#elif __STDC_VERSION__ >= 199901L && defined(NDEBUG)
    int i;
    /* try to indirectly "trigger" some compiler optimisations */
    for(i=0; i<len-3; i+=4){
        c[i] = a[i] * b[i];
        c[i+1] = a[i+1] * b[i+1];
        c[i+2] = a[i+2] * b[i+2];
        c[i+3] = a[i+3] * b[i+3];
    }
    for(; i<len; i++) /* The residual (if len was not divisable by the step size): */
        c[i] = a[i] * b[i];
#elif __STDC_VERSION__ >= 199901L 
    int i;
    for (i = 0; i < len; i++)
        c[i] = a[i] * b[i];
#else
    int i;
    for (i = 0; i < len; i++)
        c[i] = ccmulf(a[i], b[i]);
#endif
}


/* ========================================================================== */
/*                     Vector-Vector Dot Product (?vvdot)                     */
/* ========================================================================== */

void utility_svvdot
(
    const float* a,
    const float* b,
    const int len,
    float* c
)
{
    c[0] = cblas_sdot (len, a, 1, b, 1);
}

void utility_cvvdot
(
    const float_complex* a,
    const float_complex* b,
    const int len,
    CONJ_FLAG flag,
    float_complex* c
)
{
    switch(flag){
        default:
        case NO_CONJ:
            cblas_cdotu_sub(len, a, 1, b, 1, c);
            break;
        case CONJ:
            cblas_cdotc_sub(len, a, 1, b, 1, c);
            break;
    }
}


/* ========================================================================== */
/*                       Vector-Scalar Product (?vsmul)                       */
/* ========================================================================== */

void utility_svsmul
(
    float* a,
    const float* s,
    const int len,
    float* c
)
{
#if defined(SAF_USE_INTEL_IPP)
    if(c==NULL)
        ippsMulC_32f_I((Ipp32f)s[0], (Ipp32f*)a, len);
    else
        ippsMulC_32f((Ipp32f*)a, (Ipp32f)s[0], (Ipp32f*)c, len);
#elif defined(SAF_USE_APPLE_ACCELERATE)
    if(c==NULL)
        cblas_sscal(len, s[0], a, 1);
    else
        vDSP_vsmul(a, 1, s, c, 1, (vDSP_Length)len);
#else
    if (c == NULL)
        cblas_sscal(len, s[0], a, 1);
    else {
        utility_svvcopy(a, len, c); 
        cblas_sscal(len, s[0], c, 1);
    }
#endif
}

void utility_cvsmul
(
    float_complex* a,
    const float_complex* s,
    const int len,
    float_complex* c
)
{
    if (c == NULL)
        cblas_cscal(len, s, a, 1);
    else {
        cblas_ccopy(len, a, 1, c, 1);
        cblas_cscal(len, s, c, 1);
    }
}

void utility_dvsmul
(
    double* a,
    const double* s,
    const int len,
    double* c
)
{
#if defined(SAF_USE_INTEL_IPP)
    if(c==NULL)
        ippsMulC_64f_I((Ipp64f)s[0], (Ipp64f*)a, len);
    else
        ippsMulC_64f((Ipp64f*)a, (Ipp64f)s[0], (Ipp64f*)c, len);
#elif defined(SAF_USE_APPLE_ACCELERATE)
    if(c==NULL)
        cblas_dscal(len, s[0], a, 1);
    else
        vDSP_vsmulD(a, 1, s, c, 1, (vDSP_Length)len);
#else
    if (c == NULL)
        cblas_dscal(len, s[0], a, 1);
    else {
        utility_dvvcopy(a, len, c);
        cblas_dscal(len, s[0], c, 1);
    }
#endif
}

void utility_zvsmul
(
    double_complex* a,
    const double_complex* s,
    const int len,
    double_complex* c
)
{
    if (c == NULL)
        cblas_zscal(len, s, a, 1);
    else {
        cblas_zcopy(len, a, 1, c, 1);
        cblas_zscal(len, s, c, 1);
    }
}


/* ========================================================================== */
/*                       Vector-Scalar Division (?vsdiv)                      */
/* ========================================================================== */

void utility_svsdiv
(
    const float* a,
    const float* s,
    const int len,
    float* c
)
{
    if(s[0] == 0.0f){
        memset(c, 0, len*sizeof(float));
        return;
    }
#if defined(SAF_USE_INTEL_IPP)
    ippsDivC_32f((Ipp32f*)a, (Ipp32f)s[0], (Ipp32f*)c, len);
#elif defined(SAF_USE_APPLE_ACCELERATE)
    vDSP_vsdiv(a, 1, s, c, 1, (vDSP_Length)len);
#else
    cblas_scopy(len, a, 1, c, 1);
    cblas_sscal(len, 1.0f/s[0], c, 1);
#endif
}


/* ========================================================================== */
/*                       Vector-Scalar Addition (?vsadd)                      */
/* ========================================================================== */

void utility_svsadd
(
    float* a,
    const float* s,
    const int len,
    float* c
)
{
#if defined(SAF_USE_INTEL_IPP)
    ippsAddC_32f((Ipp32f*)a, (Ipp32f)s[0], (Ipp32f*)c, len);
#elif defined(SAF_USE_APPLE_ACCELERATE)
    vDSP_vsadd(a, 1, s, c, 1, (vDSP_Length)len);
#elif defined(SAF_ENABLE_SIMD)
    int i;
# if defined(__AVX512F__)
    __m512 s16 = _mm512_set1_ps(s[0]);
    for(i=0; i<(len-15); i+=16)
        _mm512_storeu_ps(c+i, _mm512_add_ps(_mm512_loadu_ps(a+i), s16));
# elif defined(__AVX__) && defined(__AVX2__)
    __m256 s8 = _mm256_set1_ps(s[0]);
    for(i=0; i<(len-7); i+=8)
        _mm256_storeu_ps(c+i, _mm256_add_ps(_mm256_loadu_ps(a+i), s8));
# elif defined(__SSE__) && defined(__SSE2__) && defined(__SSE3__)
    __m128 s4 = _mm_load_ps1(s);
    for(i=0; i<(len-3); i+=4)
        _mm_storeu_ps(c+i, _mm_add_ps(_mm_loadu_ps(a+i), s4));
# endif
    for(;i<len; i++) /* The residual (if len was not divisable by the step size): */
        c[i] = a[i] + s[0];
#else
    int i;
    for(i=0; i<len; i++)
        c[i] = a[i] + s[0];
#endif
}


/* ========================================================================== */
/*                     Vector-Scalar Subtraction (?vssub)                     */
/* ========================================================================== */

void utility_svssub
(
    float* a,
    const float* s,
    const int len,
    float* c
)
{
#if defined(SAF_USE_INTEL_IPP)
    ippsSubC_32f((Ipp32f*)a, (Ipp32f)s[0], (Ipp32f*)c, len);
#elif defined(SAF_USE_APPLE_ACCELERATE)
    float inv_s;
    inv_s = -s[0];
    vDSP_vsadd(a, 1, &inv_s, c, 1, (vDSP_Length)len);
#elif defined(SAF_ENABLE_SIMD)
    int i;
# if defined(__AVX512F__)
    __m512 s16 = _mm512_set1_ps(s[0]);
    for(i=0; i<(len-15); i+=16)
        _mm512_storeu_ps(c+i, _mm512_sub_ps(_mm512_loadu_ps(a+i), s16));
# elif defined(__AVX__) && defined(__AVX2__)
    __m256 s8 = _mm256_set1_ps(s[0]);
    for(i=0; i<(len-7); i+=8)
        _mm256_storeu_ps(c+i, _mm256_sub_ps(_mm256_loadu_ps(a+i), s8));
# elif defined(__SSE__) && defined(__SSE2__) && defined(__SSE3__)
    __m128 s4 = _mm_load_ps1(s);
    for(i=0; i<(len-3); i+=4)
        _mm_storeu_ps(c+i, _mm_sub_ps(_mm_loadu_ps(a+i), s4));
# endif
    for(;i<len; i++) /* The residual (if len was not divisable by the step size): */
        c[i] = a[i] - s[0];
#else
    int i;
    for(i=0; i<len; i++)
        c[i] = a[i] - s[0];
#endif
}


/* ========================================================================== */
/*      Sparse-Vector to Compressed-Vector (Known Indices) (?sv2cv_inds)      */
/* ========================================================================== */

void utility_ssv2cv_inds
(
    const float* sv,
    const int* inds,
    const int len,
    float* cv
)
{
    int i;
#ifdef SAF_USE_APPLE_ACCELERATE /* Unfortunately requires a malloc call */
    /* Due to Apple "logic", we first need to add 1 to all of the indicies, since "vDSP_vgathr" is going to then subtract 1 from them all... */
    vDSP_Length* inds_vDSP;
    inds_vDSP = malloc1d(len*sizeof(vDSP_Length));
    for(i=0; i<len; i++)
        inds_vDSP[i] = (vDSP_Length)(inds[i] + 1);
    vDSP_vgathr(sv, inds_vDSP, 1, cv, 1, (vDSP_Length)len);
    free(inds_vDSP);
#elif defined(SAF_USE_INTEL_MKL_LP64) || defined(SAF_USE_INTEL_MKL_ILP64)
    veclib_int* inds_tmp;
    if(sizeof(veclib_int)==sizeof(int)) /* LP64 MKL */
        cblas_sgthr(len, sv, cv, (veclib_int*)inds);
    else{ /* ILP64 MKL */
        inds_tmp = malloc1d(len*sizeof(veclib_int)); /* Unfortunately requires a malloc call */
        for(i=0; i<len; i++)
            inds_tmp[i] = (veclib_int)inds[i];
        cblas_sgthr(len, sv, cv, (veclib_int*)inds_tmp);
        free(inds_tmp);
    }
#elif defined(NDEBUG)
    /* try to indirectly "trigger" some compiler optimisations */
    for(i=0; i<len-3; i+=4){
        cv[i] = sv[inds[i]];
        cv[i+1] = sv[inds[i+1]];
        cv[i+2] = sv[inds[i+2]];
        cv[i+3] = sv[inds[i+3]];
    }
    for(; i<len; i++) /* The residual (if len was not divisable by the step size): */
        cv[i] = sv[inds[i]];
#else
    for(i=0; i<len; i++)
        cv[i] = sv[inds[i]];
#endif
}


/* ========================================================================== */
/*                     Singular-Value Decomposition (?svd)                    */
/* ========================================================================== */

/** Data structure for utility_ssvd() */
typedef struct _utility_ssvd_data {
    int maxDim1, maxDim2;
    veclib_int currentWorkSize;
    float* a, *s, *u, *vt;
    float* work;
}utility_ssvd_data;

void utility_ssvd_create(void ** const phWork, int maxDim1, int maxDim2)
{
    *phWork = malloc1d(sizeof(utility_ssvd_data));
    utility_ssvd_data *h = (utility_ssvd_data*)(*phWork);

    h->maxDim1 = maxDim1;
    h->maxDim2 = maxDim2;
    h->currentWorkSize = 0;
    h->a = malloc1d(maxDim1*maxDim2*sizeof(float));
    h->s = malloc1d(SAF_MIN(maxDim2,maxDim1)*sizeof(float));
    h->u = malloc1d(maxDim1*maxDim1*sizeof(float));
    h->vt = malloc1d(maxDim2*maxDim2*sizeof(float));
    h->work = NULL;
}

void utility_ssvd_destroy(void ** const phWork)
{
    utility_ssvd_data *h = (utility_ssvd_data*)(*phWork);

    if(h!=NULL){
        free(h->a);
        free(h->s);
        free(h->u);
        free(h->vt);
        free(h->work);

        free(h);
        h=NULL;
        *phWork = NULL;
     }
}

void utility_ssvd
(
    void* const hWork,
    const float* A,
    const int dim1,
    const int dim2,
    float* U,
    float* S,
    float* V,
    float* sing
)
{
    utility_ssvd_data *h;
    veclib_int i, j, m, n, lda, ldu, ldvt, info;
    veclib_int lwork;
    float wkopt;

    m = dim1; n = dim2; lda = dim1; ldu = dim1; ldvt = dim2;

    /* Work struct */
    if(hWork==NULL)
        utility_ssvd_create((void**)&h, dim1, dim2);
    else{
        h = (utility_ssvd_data*)(hWork);
#ifndef NDEBUG
        saf_assert(dim1<=h->maxDim1, "dim1 exceeds the maximum length specified");
        saf_assert(dim2<=h->maxDim2, "dim2 exceeds the maximum length specified");
#endif
    }

    /* store in column major order */
    for(i=0; i<dim1; i++)
        for(j=0; j<dim2; j++)
            h->a[j*dim1+i] = A[i*dim2 +j];

    /* Query how much "work" memory is required */
    lwork = -1;
#if defined(SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE)
    sgesvd_( "A", "A", &m, &n, h->a, &lda, h->s, h->u, &ldu, h->vt, &ldvt, &wkopt, &lwork, &info );
#elif defined(SAF_VECLIB_USE_CLAPACK_INTERFACE)
    saf_print_error("No such implementation available in ATLAS CLAPACK");
#elif defined(SAF_VECLIB_USE_LAPACKE_INTERFACE)
    info = LAPACKE_sgesvd_work(CblasColMajor, 'A', 'A', m, n, h->a, lda, h->s, h->u, ldu, h->vt, ldvt, &wkopt, lwork);
#endif
    lwork = (veclib_int)wkopt;
    if(lwork>h->currentWorkSize){
        h->currentWorkSize = lwork;
        h->work = realloc1d(h->work, h->currentWorkSize*sizeof(float));
    }

    /* perform the singular value decomposition */
#if defined(SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE)
    sgesvd_( "A", "A", &m, &n, h->a, &lda, h->s, h->u, &ldu, h->vt, &ldvt, h->work, &lwork, &info );
#elif defined(SAF_VECLIB_USE_CLAPACK_INTERFACE)
    saf_print_error("No such implementation available in ATLAS CLAPACK");
#elif defined(SAF_VECLIB_USE_LAPACKE_INTERFACE)
    info = LAPACKE_sgesvd_work(CblasColMajor, 'A', 'A', m, n, h->a, lda, h->s, h->u, ldu, h->vt, ldvt, h->work, lwork);
#endif

    /* svd failed to converge */
    if( info != 0 ) {
        if (S != NULL)
            memset(S, 0, dim1*dim2*sizeof(float));
        if (U != NULL)
            memset(U, 0, dim1*dim1*sizeof(float));
        if (V != NULL)
            memset(V, 0, dim2*dim2*sizeof(float));
        if (sing != NULL)
            memset(sing, 0, SAF_MIN(dim1, dim2)*sizeof(float));
#ifndef NDEBUG
        /* The SVD failed to converge, or the input matrix contained illegal
         * values so no solution was attempted. In these cases this function
         * will zero all output matrices and/or vectors. */
        saf_print_warning("Could not compute SVD in utility_ssvd(). Output matrices/vectors have been zeroed.");
#endif
    }
    /* svd successful */
    else {
        if (S != NULL){
            memset(S, 0, dim1*dim2*sizeof(float));
            /* singular values on the diagonal MIN(dim1, dim2). The remaining elements are 0.  */
            for(i=0; i<SAF_MIN(dim1, dim2); i++)
                S[i*dim2+i] = h->s[i];
        }
        
        /*return as row-major*/
        if (U != NULL)
            for(i=0; i<dim1; i++)
                for(j=0; j<dim1; j++)
                    U[i*dim1+j] = h->u[j*dim1+i];
        
        /* lapack returns VT, i.e. row-major V already */
        if (V != NULL)
            for(i=0; i<dim2; i++)
                for(j=0; j<dim2; j++)
                    V[i*dim2+j] = h->vt[i*dim2+j];
        if (sing != NULL)
            for(i=0; i<SAF_MIN(dim1, dim2); i++)
                sing[i] = h->s[i];
    }

    if(hWork == NULL)
        utility_ssvd_destroy((void**)&h);
}

/** Data structure for utility_csvd() */
typedef struct _utility_csvd_data {
    int maxDim1, maxDim2;
    veclib_int currentWorkSize;
    float_complex* a, *u, *vt, *work;
    float *s, *rwork;
}utility_csvd_data;

void utility_csvd_create(void ** const phWork, int maxDim1, int maxDim2)
{
    *phWork = malloc1d(sizeof(utility_csvd_data));
    utility_csvd_data *h = (utility_csvd_data*)(*phWork);

    h->maxDim1 = maxDim1;
    h->maxDim2 = maxDim2;
    h->currentWorkSize = 0;
    h->a = malloc1d(maxDim1*maxDim2*sizeof(float_complex));
    h->s = malloc1d(SAF_MIN(maxDim2,maxDim1)*sizeof(float));
    h->u = malloc1d(maxDim1*maxDim1*sizeof(float_complex));
    h->vt = malloc1d(maxDim2*maxDim2*sizeof(float_complex));
    h->rwork = malloc1d(maxDim1*SAF_MAX(1, 5*SAF_MIN(maxDim2,maxDim1))*sizeof(float));
    h->work = NULL;
}

void utility_csvd_destroy(void ** const phWork)
{
    utility_csvd_data *h = (utility_csvd_data*)(*phWork);

    if(h!=NULL){
        free(h->a);
        free(h->s);
        free(h->u);
        free(h->vt);
        free(h->work);
        free(h->rwork);

        free(h);
        h=NULL;
        *phWork = NULL;
     }
}

void utility_csvd
(
    void* const hWork,
    const float_complex* A,
    const int dim1,
    const int dim2,
    float_complex* U,
    float_complex* S,
    float_complex* V,
    float* sing
)
{
    utility_csvd_data *h;
    veclib_int m, n, lda, ldu, ldvt, info;
    veclib_int lwork;
    float_complex wkopt;
#if defined(SAF_USE_INTEL_MKL_LP64) || defined(SAF_USE_INTEL_MKL_ILP64)
    const MKL_Complex8 calpha = {1.0f, 0.0f};
#else
    int i, j;
#endif

    m = dim1; n = dim2; lda = dim1; ldu = dim1; ldvt = dim2;

    /* Work struct */
    if(hWork==NULL)
        utility_csvd_create((void**)&h, dim1, dim2);
    else{
        h = (utility_csvd_data*)(hWork);
#ifndef NDEBUG
        saf_assert(dim1<=h->maxDim1, "dim1 exceeds the maximum length specified");
        saf_assert(dim2<=h->maxDim2, "dim2 exceeds the maximum length specified");
#endif
    }

    /* store in column major order (i.e. transpose) */
#if defined(SAF_USE_INTEL_MKL_LP64) || defined(SAF_USE_INTEL_MKL_ILP64)
    MKL_Comatcopy('R', 'T', dim1, dim2, calpha, (veclib_float_complex*)A, dim2, (veclib_float_complex*)h->a, dim1);
#else
    for(i=0; i<dim1; i++)
        for(j=0; j<dim2; j++)
            h->a[j*dim1+i] = A[i*dim2+j];
#endif

    /* Query how much "work" memory is required */
    lwork = -1;
#if defined(SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE)
    cgesvd_( "A", "A", (veclib_int*)&m, (veclib_int*)&n, (veclib_float_complex*)h->a, (veclib_int*)&lda, h->s, (veclib_float_complex*)h->u, (veclib_int*)&ldu,
            (veclib_float_complex*)h->vt, &ldvt, (veclib_float_complex*)&wkopt, &lwork, h->rwork, (veclib_int*)&info );
#elif defined(SAF_VECLIB_USE_CLAPACK_INTERFACE)
    saf_print_error("No such implementation available in ATLAS CLAPACK");
#elif defined(SAF_VECLIB_USE_LAPACKE_INTERFACE)
    info = LAPACKE_cgesvd_work(CblasColMajor, 'A', 'A', m, n, (veclib_float_complex*)h->a, lda, h->s, (veclib_float_complex*)h->u, ldu,
                               (veclib_float_complex*)h->vt, ldvt, (veclib_float_complex*)&wkopt, lwork, h->rwork);
#endif
    lwork = (veclib_int)(crealf(wkopt)+0.01f);
    if(lwork>h->currentWorkSize){
        h->currentWorkSize = lwork;
        h->work = realloc1d(h->work, h->currentWorkSize*sizeof(float_complex));
    }

    /* perform the singular value decomposition */
#if defined(SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE)
    cgesvd_( "A", "A", &m, &n, (veclib_float_complex*)h->a, &lda, h->s, (veclib_float_complex*)h->u, &ldu, (veclib_float_complex*)h->vt, &ldvt,
            (veclib_float_complex*)h->work, &lwork, h->rwork, &info);
#elif defined(SAF_VECLIB_USE_CLAPACK_INTERFACE)
    saf_print_error("No such implementation available in ATLAS CLAPACK");
#elif defined(SAF_VECLIB_USE_LAPACKE_INTERFACE)
    info = LAPACKE_cgesvd_work(CblasColMajor, 'A', 'A', m, n, (veclib_float_complex*)h->a, lda, h->s, (veclib_float_complex*)h->u, ldu, (veclib_float_complex*)h->vt, ldvt,
                              (veclib_float_complex*)h->work, lwork, h->rwork);
#endif

    /* svd failed to converge */
    if( info != 0 ) {
        if (S != NULL)
            memset(S, 0, dim1*dim2*sizeof(float_complex));
        if (U != NULL)
            memset(U, 0, dim1*dim1*sizeof(float_complex));
        if (V != NULL)
            memset(V, 0, dim2*dim2*sizeof(float_complex));
        if (sing != NULL)
            memset(sing, 0, SAF_MIN(dim1, dim2)*sizeof(float_complex));
#ifndef NDEBUG
        /* The SVD failed to converge, or the input matrix contained illegal
         * values so no solution was attempted. In these cases this function
         * will zero all output matrices and/or vectors. */
        saf_print_warning("Could not compute SVD in utility_csvd(). Output matrices/vectors have been zeroed.");
#endif
    }
    /* svd successful */
    else {
        if (S != NULL){
            /* singular values on the diagonal MIN(dim1, dim2). The remaining elements are 0.  */
            memset(S, 0, dim1*dim2*sizeof(float_complex));
            cblas_scopy(SAF_MIN(dim1, dim2), h->s, 1, (float*)S, 2*(dim2+1));
        }

        /*return as row-major*/
        if (U!=NULL){
#if defined(SAF_USE_INTEL_MKL_LP64) || defined(SAF_USE_INTEL_MKL_ILP64)
            MKL_Comatcopy('R', 'T', dim1, dim1, calpha, (veclib_float_complex*)h->u, dim1, (veclib_float_complex*)U, dim1);
#else
            for(i=0; i<dim1; i++)
                for(j=0; j<dim1; j++)
                    U[i*dim1+j] = h->u[j*dim1+i];
#endif
        }

        /* lapack returns V^T, i.e. row-major V already, but we need V^H! */
        if (V != NULL){
            cblas_ccopy(dim2*dim2, h->vt, 1, V, 1);
            cblas_sscal(dim2*dim2, -1.0f, ((float*)V)+1, 2); /* conj */
        }
        
        if (sing != NULL)
            cblas_scopy(SAF_MIN(dim1, dim2), h->s, 1, sing, 1);
    }

    if(hWork == NULL)
        utility_csvd_destroy((void**)&h);
}


/* ========================================================================== */
/*                 Symmetric Eigenvalue Decomposition (?seig)                 */
/* ========================================================================== */

/** Data structure for utility_sseig() */
typedef struct _utility_sseig_data {
    int maxDim;
    veclib_int currentWorkSize;
    float* w;
    float* a;
    float* work;
}utility_sseig_data;

void utility_sseig_create(void ** const phWork, int maxDim)
{
    *phWork = malloc1d(sizeof(utility_sseig_data));
    utility_sseig_data *h = (utility_sseig_data*)(*phWork);

    h->maxDim = maxDim;
    h->currentWorkSize = 0;
    h->w = malloc1d(maxDim*sizeof(float));
    h->a = malloc1d(maxDim*maxDim*sizeof(float));
    h->work = NULL;
}

void utility_sseig_destroy(void ** const phWork)
{
    utility_sseig_data *h = (utility_sseig_data*)(*phWork);

    if(h!=NULL){
        free(h->w);
        free(h->a);
        free(h->work);

        free(h);
        h=NULL;
        *phWork = NULL;
     }
}

void utility_sseig
(
    void* const hWork,
    const float* A,
    const int dim,
    int sortDecFLAG,
    float* V,
    float* D,
    float* eig
)
{
    utility_sseig_data *h;
    veclib_int i, j, n, lda, info;
    veclib_int lwork;
    float wkopt;

    n = dim;
    lda = dim;

    /* Work struct */
    if(hWork==NULL)
        utility_sseig_create((void**)&h, dim);
    else{
        h = (utility_sseig_data*)(hWork);
#ifndef NDEBUG
        saf_assert(dim<=h->maxDim, "dim exceeds the maximum length specified");
#endif
    }
    
    /* store in column major order (i.e. transpose) */
    for(i=0; i<dim; i++)
        for(j=0; j<dim; j++)
            h->a[i*dim+j] = A[j*dim+i];

    /* Query how much "work" memory is required */
    lwork = -1;
#if defined(SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE)
    ssyev_( "Vectors", "Upper", &n, h->a, &lda, h->w, &wkopt, &lwork, &info );
#elif defined(SAF_VECLIB_USE_CLAPACK_INTERFACE)
    saf_print_error("No such implementation available in ATLAS CLAPACK");
#elif defined(SAF_VECLIB_USE_LAPACKE_INTERFACE)
    info = LAPACKE_ssyev_work(CblasColMajor, 'V', 'U', n, h->a, lda, h->w, &wkopt, lwork);
#endif
    lwork = (veclib_int)wkopt;
    if(lwork>h->currentWorkSize){
        h->currentWorkSize = lwork;
        h->work = realloc1d(h->work, h->currentWorkSize*sizeof(float));
    }

    /* solve the eigenproblem */
#if defined(SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE)
    ssyev_( "Vectors", "Upper", &n, h->a, &lda, h->w, h->work, &lwork, &info );
#elif defined(SAF_VECLIB_USE_CLAPACK_INTERFACE)
    saf_print_error("No such implementation available in ATLAS CLAPACK");
#elif defined(SAF_VECLIB_USE_LAPACKE_INTERFACE)
    info = LAPACKE_ssyev_work(CblasColMajor, 'V', 'U', n, h->a, lda, h->w, h->work, lwork);
#endif
    
    /* output */
    if(D!=NULL)
        memset(D, 0, dim*dim*sizeof(float));
    if( info != 0 ) {
        /* failed to converge and find the eigenvalues */
        if(V!=NULL)
            memset(V, 0, dim*dim*sizeof(float));
#ifndef NDEBUG
        /* Failed to compute all of the eigenvalues, some eigenvectors have not
         * been computed, or the input matrix contained illegal values so no
         * solution was attempted. In these cases the function will zero all
         * output matrices and vectors. */
        saf_print_warning("Could not compute EVD in utility_sseig(). Output matrices/vectors have been zeroed.");
#endif
    }
    else{
        if(sortDecFLAG){
            for(i=0; i<dim; i++){
                if(V!=NULL)
                    for(j=0; j<dim; j++)
                        V[i*dim+j] = h->a[(dim-j-1)*dim+i]; /* transpose, back to row-major and reverse order */
                if(D!=NULL)
                    D[i*dim+i] = h->w[dim-i-1]; /* store along the diagonal, reversing the order */
                if(eig!=NULL)
                    eig[i] = h->w[dim-i-1];
            }
        }
        else{
            for(i=0; i<dim; i++){
                if(V!=NULL)
                    for(j=0; j<dim; j++)
                        V[i*dim+j] = h->a[j*dim+i]; /* transpose, back to row-major */
                if(D!=NULL)
                    D[i*dim+i] = h->w[i]; /* store along the diagonal */
                if(eig!=NULL)
                    eig[i] = h->w[i];
            }
        }
    }

    if(hWork == NULL)
        utility_sseig_destroy((void**)&h);
}

/** Data structure for utility_cseig() */
typedef struct _utility_cseig_data {
    int maxDim;
    veclib_int currentWorkSize;
    float* rwork;
    float* w;
    float_complex* a;
    float_complex* work;
}utility_cseig_data;

void utility_cseig_create(void ** const phWork, int maxDim)
{
    *phWork = malloc1d(sizeof(utility_cseig_data));
    utility_cseig_data *h = (utility_cseig_data*)(*phWork);

    h->maxDim = maxDim;
    h->currentWorkSize = SAF_MAX(1, 2*maxDim-1);
    h->rwork = malloc1d((3*maxDim-2)*sizeof(float));
    h->w = malloc1d(maxDim*sizeof(float));
    h->a = malloc1d(maxDim*maxDim*sizeof(float_complex));
    h->work = malloc1d(h->currentWorkSize*sizeof(float_complex));
}

void utility_cseig_destroy(void ** const phWork)
{
    utility_cseig_data *h = (utility_cseig_data*)(*phWork);

    if(h!=NULL){
        free(h->rwork);
        free(h->w);
        free(h->a);
        free(h->work);

        free(h);
        h=NULL;
        *phWork = NULL;
     }
}

void utility_cseig
(
    void* const hWork,
    const float_complex* A,
    const int dim,
    int sortDecFLAG,
    float_complex* V,
    float_complex* D,
    float* eig
)
{
    utility_cseig_data *h;
    veclib_int i, n, lda, info;
    veclib_int lwork;
    float_complex wkopt;
#if defined(SAF_USE_INTEL_MKL_LP64) || defined(SAF_USE_INTEL_MKL_ILP64)
    const MKL_Complex8 calpha = {1.0f, 0.0f};
#else
    int j;
#endif

    n = dim;
    lda = dim;
    
    /* Work struct */
    if(hWork==NULL)
        utility_cseig_create((void**)&h, dim);
    else{
        h = (utility_cseig_data*)(hWork);
#ifndef NDEBUG
        saf_assert(dim<=h->maxDim, "dim exceeds the maximum length specified");
#endif
    }
    
    /* store in column major order (i.e. transpose) */
#if defined(SAF_USE_INTEL_MKL_LP64) || defined(SAF_USE_INTEL_MKL_ILP64)
    MKL_Comatcopy('R', 'T', dim, dim, calpha, (veclib_float_complex*)A, dim, (veclib_float_complex*)h->a, dim);
#else
    for(i=0; i<dim; i++)
        for(j=0; j<dim; j++)
            h->a[i*dim+j] = A[j*dim+i];
#endif

    /* Query how much "work" memory is required */
    lwork = -1;
#if defined(SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE)
    cheev_( "Vectors", "Upper", &n, (veclib_float_complex*)h->a, &lda, h->w, (veclib_float_complex*)&wkopt, &lwork, h->rwork, &info );
#elif defined(SAF_VECLIB_USE_CLAPACK_INTERFACE)
    saf_print_error("No such implementation available in ATLAS CLAPACK");
#elif defined(SAF_VECLIB_USE_LAPACKE_INTERFACE)
    info = LAPACKE_cheev_work(CblasColMajor, 'V', 'U', n, (veclib_float_complex*)h->a, lda, h->w, (veclib_float_complex*)&wkopt, lwork, h->rwork);
#endif
    lwork = (veclib_int)crealf(wkopt);
    if(lwork>h->currentWorkSize){
        h->currentWorkSize = lwork;
        h->work = realloc1d(h->work, h->currentWorkSize*sizeof(float_complex));
    }

    /* solve the eigenproblem */
#if defined(SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE)
    cheev_( "Vectors", "Upper", &n, (veclib_float_complex*)h->a, &lda, h->w, (veclib_float_complex*)h->work, &lwork, h->rwork, &info );
#elif defined(SAF_VECLIB_USE_CLAPACK_INTERFACE)
    saf_print_error("No such implementation available in ATLAS CLAPACK");
#elif defined(SAF_VECLIB_USE_LAPACKE_INTERFACE)
    info = LAPACKE_cheev_work(CblasColMajor, 'V', 'U', n, (veclib_float_complex*)h->a, lda, h->w, (veclib_float_complex*)h->work, lwork, h->rwork);
#endif

    /* output */
    if(D!=NULL)
        memset(D, 0, dim*dim*sizeof(float_complex));
    if( info != 0 ) {
        /* failed to converge and find the eigenvalues */
        if(V!=NULL)
            memset(V, 0, dim*dim*sizeof(float_complex));
#ifndef NDEBUG
        /* Failed to compute all of the eigenvalues, some eigenvectors have not
         * been computed, or the input matrix contained illegal values so no
         * solution was attempted. In these cases the function will zero all
         * output matrices and vectors. */
        saf_print_warning("Could not compute EVD in utility_cseig(). Output matrices/vectors have been zeroed.");
#endif
    }
    
    /* transpose, back to row-major and reverse order if sortDecFlag==1 */
    else{
        if(sortDecFLAG && V!=NULL)
            for(i=0; i<(int)((float)dim/2.0f); i++)
                cblas_cswap(dim, &h->a[i*dim], 1, &h->a[(dim-i-1)*dim], 1);

        if(V!=NULL){
#if defined(SAF_USE_INTEL_MKL_LP64) || defined(SAF_USE_INTEL_MKL_ILP64)
            MKL_Comatcopy('R', 'T', dim, dim, calpha, (veclib_float_complex*)h->a, dim, (veclib_float_complex*)V, dim);
#else
            for(i=0; i<dim; i++)
                for(j=0; j<dim; j++)
                    V[i*dim+j] = h->a[j*dim+i];
#endif
        }

        if(sortDecFLAG){
            for(i=0; i<dim; i++) {
                if(D!=NULL)
                    D[i*dim+i] = cmplxf(h->w[dim-i-1], 0.0f); /* store along the diagonal, reversing the order */
                if(eig!=NULL)
                    eig[i] = h->w[dim-i-1];
            }
        }
        else {
            for(i=0; i<dim; i++){
                if(D!=NULL)
                    D[i*dim+i] = cmplxf(h->w[i], 0.0f); /* store along the diagonal */
                if(eig!=NULL)
                    eig[i] = h->w[i];
            }
        }
    }

    if(hWork == NULL)
        utility_cseig_destroy((void**)&h);
}


/* ========================================================================== */
/*                     Eigenvalues of Matrix Pair (?eigmp)                    */
/* ========================================================================== */

/** Data structure for utility_ceigmp() */
typedef struct _utility_ceigmp_data {
    int maxDim;
    veclib_int currentWorkSize;
    float_complex* a, *b, *vl, *vr, *alpha, *beta;
    float* rwork;
    float_complex* work;
}utility_ceigmp_data;

void utility_ceigmp_create(void** const phWork, int maxDim)
{
    *phWork = malloc1d(sizeof(utility_ceigmp_data));
    utility_ceigmp_data *h = (utility_ceigmp_data*)(*phWork);

    h->maxDim = maxDim;
    h->currentWorkSize = 4*maxDim;
    h->rwork = malloc1d(4*(h->currentWorkSize)*sizeof(float));
    h->a = malloc1d(maxDim*maxDim*sizeof(float_complex));
    h->b = malloc1d(maxDim*maxDim*sizeof(float_complex));
    h->vl = malloc1d(maxDim*maxDim*sizeof(float_complex));
    h->vr = malloc1d(maxDim*maxDim*sizeof(float_complex));
    h->alpha = malloc1d(maxDim*sizeof(float_complex));
    h->beta = malloc1d(maxDim*sizeof(float_complex));
    h->work = malloc1d(h->currentWorkSize*sizeof(float_complex));
}

void utility_ceigmp_destroy(void ** const phWork)
{
    utility_ceigmp_data *h = (utility_ceigmp_data*)(*phWork);

    if(h!=NULL){
        free(h->rwork);
        free(h->a);
        free(h->b);
        free(h->vl);
        free(h->vr);
        free(h->alpha);
        free(h->beta);
        free(h->work);

        free(h);
        h=NULL;
        *phWork = NULL;
     }
}

void utility_ceigmp
(
    void* const hWork,
    const float_complex* A,
    const float_complex* B,
    const int dim,
    float_complex* VL,
    float_complex* VR,
    float_complex* D
)
{
    utility_ceigmp_data *h;
    veclib_int i, j;
    veclib_int n, lda, ldb, ldvl, ldvr, info;
    veclib_int lwork;
    
    n = lda = ldb = ldvl = ldvr = dim;

    /* Work struct */
    if(hWork==NULL)
        utility_ceigmp_create((void**)&h, dim);
    else{
        h = (utility_ceigmp_data*)(hWork);
#ifndef NDEBUG
        saf_assert(dim<=h->maxDim, "dim exceeds the maximum length specified");
#endif
    }
    
    /* store in column major order */
    for(i=0; i<dim; i++)
        for(j=0; j<dim; j++)
            h->a[j*dim+i] = A[i*dim+j];
    for(i=0; i<dim; i++)
        for(j=0; j<dim; j++)
            h->b[j*dim+i] = B[i*dim+j];

    /* solve eigen problem */
    lwork = h->currentWorkSize;
#if defined(SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE)
    cggev_("V", "V", &n, (veclib_float_complex*)h->a, &lda, (veclib_float_complex*)h->b, &ldb, (veclib_float_complex*)h->alpha, (veclib_float_complex*)h->beta,
           (veclib_float_complex*)h->vl, &ldvl, (veclib_float_complex*)h->vr, &ldvr, (veclib_float_complex*)h->work, &lwork, h->rwork, &info);
#elif defined(SAF_VECLIB_USE_CLAPACK_INTERFACE)
    saf_print_error("No such implementation available in ATLAS CLAPACK");
#elif defined(SAF_VECLIB_USE_LAPACKE_INTERFACE)
    info = LAPACKE_cggev_work(CblasColMajor, 'V', 'V', n, (veclib_float_complex*)h->a, lda, (veclib_float_complex*)h->b, ldb, (veclib_float_complex*)h->alpha, (veclib_float_complex*)h->beta,
                              (veclib_float_complex*)h->vl, ldvl, (veclib_float_complex*)h->vr, ldvr, (veclib_float_complex*)h->work, lwork, h->rwork);
#endif
    
    if(D!=NULL)
        memset(D, 0, dim*dim*sizeof(float_complex));
    
    /* failed to converge and find the eigenvalues */
    if( info != 0 ) {
        if(VL!=NULL)
            memset(VL, 0, dim*dim*sizeof(float_complex));
        if(VR!=NULL)
            memset(VR, 0, dim*dim*sizeof(float_complex));
#ifndef NDEBUG
        /* Failed to compute all of the eigenvalues, some eigenvectors have not
         * been computed, or the input matrix contained illegal values so no
         * solution was attempted. In these cases the function will zero all
         * output matrices and vectors. */
        saf_print_warning("Could not compute EVD in utility_ceigmp(). Output matrices/vectors have been zeroed.");
#endif
    }
    /* transpose, back to row-major */
    else{
        if(D!=NULL)
            for(i=0; i<dim; i++)
                D[i*dim+i] = ccdivf(h->alpha[i],h->beta[i]);
        
        if(VL!=NULL)
            for(i=0; i<dim; i++)
                for(j=0; j<dim; j++)
                    VL[i*dim+j] = h->vl[j*dim+i];
        if(VR!=NULL)
            for(i=0; i<dim; i++)
                for(j=0; j<dim; j++)
                    VR[i*dim+j] = h->vr[j*dim+i];
    }

    if(hWork == NULL)
        utility_ceigmp_destroy((void**)&h);
}

/** Data structure for utility_zeigmp() */
typedef struct _utility_zeigmp_data {
    int maxDim;
    veclib_int currentWorkSize;
    double_complex* a, *b, *vl, *vr, *alpha, *beta;
    double* rwork;
    double_complex* work;
}utility_zeigmp_data;

void utility_zeigmp_create(void** const phWork, int maxDim)
{
    *phWork = malloc1d(sizeof(utility_zeigmp_data));
    utility_zeigmp_data *h = (utility_zeigmp_data*)(*phWork);

    h->maxDim = maxDim;
    h->currentWorkSize = 4*maxDim;
    h->rwork = malloc1d(4*(h->currentWorkSize)*sizeof(double));
    h->a = malloc1d(maxDim*maxDim*sizeof(double_complex));
    h->b = malloc1d(maxDim*maxDim*sizeof(double_complex));
    h->vl = malloc1d(maxDim*maxDim*sizeof(double_complex));
    h->vr = malloc1d(maxDim*maxDim*sizeof(double_complex));
    h->alpha = malloc1d(maxDim*sizeof(double_complex));
    h->beta = malloc1d(maxDim*sizeof(double_complex));
    h->work = malloc1d(h->currentWorkSize*sizeof(double_complex));
}

void utility_zeigmp_destroy(void ** const phWork)
{
    utility_zeigmp_data *h = (utility_zeigmp_data*)(*phWork);

    if(h!=NULL){
        free(h->rwork);
        free(h->a);
        free(h->b);
        free(h->vl);
        free(h->vr);
        free(h->alpha);
        free(h->beta);
        free(h->work);

        free(h);
        h=NULL;
        *phWork = NULL;
     }
}

void utility_zeigmp
(
    void* const hWork,
    const double_complex* A,
    const double_complex* B,
    const int dim,
    double_complex* VL,
    double_complex* VR,
    double_complex* D
)
{
    utility_zeigmp_data *h;
    veclib_int i, j;
    veclib_int n, lda, ldb, ldvl, ldvr, info;
    veclib_int lwork;
    
    n = lda = ldb = ldvl = ldvr = dim;

    /* Work struct */
    if(hWork==NULL)
        utility_zeigmp_create((void**)&h, dim);
    else{
        h = (utility_zeigmp_data*)(hWork);
#ifndef NDEBUG
        saf_assert(dim<=h->maxDim, "dim exceeds the maximum length specified");
#endif
    }
    
    /* store in column major order */
    for(i=0; i<dim; i++)
        for(j=0; j<dim; j++)
            h->a[j*dim+i] = A[i*dim+j];
    for(i=0; i<dim; i++)
        for(j=0; j<dim; j++)
            h->b[j*dim+i] = B[i*dim+j]; /* store in column major order */
    
    /* solve eigen problem */
    lwork = h->currentWorkSize;
#if defined(SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE)
    zggev_("V", "V", &n, (veclib_double_complex*)h->a, &lda, (veclib_double_complex*)h->b, &ldb, (veclib_double_complex*)h->alpha, (veclib_double_complex*)h->beta,
           (veclib_double_complex*)h->vl, &ldvl, (veclib_double_complex*)h->vr, &ldvr, (veclib_double_complex*)h->work, &lwork, h->rwork, &info);
#elif defined(SAF_VECLIB_USE_CLAPACK_INTERFACE)
    saf_print_error("No such implementation available in ATLAS CLAPACK");
#elif defined(SAF_VECLIB_USE_LAPACKE_INTERFACE)
    info = LAPACKE_zggev_work(CblasColMajor, 'V', 'V', n, (veclib_double_complex*)h->a, lda, (veclib_double_complex*)h->b, ldb, (veclib_double_complex*)h->alpha, (veclib_double_complex*)h->beta,
                              (veclib_double_complex*)h->vl, ldvl, (veclib_double_complex*)h->vr, ldvr, (veclib_double_complex*)h->work, lwork, h->rwork);
#endif
    
    if(D!=NULL)
        memset(D, 0, dim*dim*sizeof(double_complex));
    
    /* failed to converge and find the eigenvalues */
    if( info != 0 ) {
        if(VL!=NULL)
            memset(VL, 0, dim*dim*sizeof(double_complex));
        if(VR!=NULL)
            memset(VR, 0, dim*dim*sizeof(double_complex));
#ifndef NDEBUG
        /* Failed to compute all of the eigenvalues, some eigenvectors have not
         * been computed, or the input matrix contained illegal values so no
         * solution was attempted. In these cases the function will zero all
         * output matrices and vectors. */
        saf_print_warning("Could not compute EVD in utility_zeigmp(). Output matrices/vectors have been zeroed.");
#endif
    }
    /* transpose, back to row-major */
    else{
        if(D!=NULL)
            for(i=0; i<dim; i++)
                D[i*dim+i] = ccdiv(h->alpha[i],h->beta[i]);
        
        if(VL!=NULL)
            for(i=0; i<dim; i++)
                for(j=0; j<dim; j++)
                    VL[i*dim+j] = h->vl[j*dim+i];
        if(VR!=NULL)
            for(i=0; i<dim; i++)
                for(j=0; j<dim; j++)
                    VR[i*dim+j] = h->vr[j*dim+i];
    }

    if(hWork == NULL)
        utility_zeigmp_destroy((void**)&h);
}


/* ========================================================================== */
/*                       Eigenvalue Decomposition (?eig)                      */
/* ========================================================================== */

/** Data structure for utility_ceig() */
typedef struct _utility_ceig_data {
    int maxDim;
    veclib_int currentWorkSize;
    float_complex *w, *vl, *vr, *a;
    float* rwork;
    float_complex* work;
}utility_ceig_data;

void utility_ceig_create(void ** const phWork, int maxDim)
{
    *phWork = malloc1d(sizeof(utility_ceig_data));
    utility_ceig_data *h = (utility_ceig_data*)(*phWork);

    h->maxDim = maxDim;
    h->currentWorkSize = 0;
    h->rwork = malloc1d(4*maxDim*sizeof(float));
    h->w = malloc1d(maxDim*sizeof(float_complex));
    h->vl = malloc1d(maxDim*maxDim*sizeof(float_complex));
    h->vr = malloc1d(maxDim*maxDim*sizeof(float_complex));
    h->a = malloc1d(maxDim*maxDim*sizeof(float_complex));
    h->work = NULL;
}

void utility_ceig_destroy(void ** const phWork)
{
    utility_ceig_data *h = (utility_ceig_data*)(*phWork);

    if(h!=NULL){
        free(h->rwork);
        free(h->w);
        free(h->vl);
        free(h->vr);
        free(h->a);
        free(h->work);

        free(h);
        h=NULL;
        *phWork = NULL;
     }
}

void utility_ceig
(
    void* const hWork,
    const float_complex* A,
    const int dim,
    float_complex* VL,
    float_complex* VR,
    float_complex* D,
    float_complex* eig
)
{
    utility_ceig_data *h;
    veclib_int i, j, n, lda, ldvl, ldvr, info;
    veclib_int lwork;
    float_complex wkopt;

    n = lda = ldvl = ldvr = dim;

    /* Work struct */
    if(hWork==NULL)
        utility_ceig_create((void**)&h, dim);
    else{
        h = (utility_ceig_data*)(hWork);
#ifndef NDEBUG
        saf_assert(dim<=h->maxDim, "dim exceeds the maximum length specified");
#endif
    }

    /* store in column major order (i.e. transpose) */
    for(i=0; i<dim; i++)
        for(j=0; j<dim; j++)
            h->a[i*dim+j] = A[j*dim+i];

    /* Query how much "work" memory is required */
    lwork = -1;
#if defined(SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE)
    cgeev_( "Vectors", "Vectors", &n, (veclib_float_complex*)h->a, &lda, (veclib_float_complex*)h->w, (veclib_float_complex*)h->vl, &ldvl,
           (veclib_float_complex*)h->vr, &ldvr, (veclib_float_complex*)&wkopt, &lwork, h->rwork, &info );
#elif defined(SAF_VECLIB_USE_CLAPACK_INTERFACE)
    saf_print_error("No such implementation available in ATLAS CLAPACK");
#elif defined(SAF_VECLIB_USE_LAPACKE_INTERFACE)
    info = LAPACKE_cgeev_work(CblasColMajor, 'V', 'V', n, (veclib_float_complex*)h->a, lda, (veclib_float_complex*)h->w, (veclib_float_complex*)h->vl, ldvl,
                              (veclib_float_complex*)h->vr, ldvr, (veclib_float_complex*)&wkopt, lwork, h->rwork);
#endif
    lwork = (veclib_int)crealf(wkopt);
    if(lwork>h->currentWorkSize){
        h->currentWorkSize = lwork;
        h->work = realloc1d(h->work, h->currentWorkSize*sizeof(float_complex));
    }

    /* solve the eigenproblem */
#if defined(SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE)
    cgeev_( "Vectors", "Vectors", &n, (veclib_float_complex*)h->a, &lda, (veclib_float_complex*)h->w, (veclib_float_complex*)h->vl, &ldvl,
           (veclib_float_complex*)h->vr, &ldvr, (veclib_float_complex*)h->work, &lwork, h->rwork, &info );
#elif defined(SAF_VECLIB_USE_CLAPACK_INTERFACE)
    saf_print_error("No such implementation available in ATLAS CLAPACK");
#elif defined(SAF_VECLIB_USE_LAPACKE_INTERFACE)
    info = LAPACKE_cgeev_work(CblasColMajor, 'V', 'V', n, (veclib_float_complex*)h->a, lda, (veclib_float_complex*)h->w, (veclib_float_complex*)h->vl, ldvl,
                              (veclib_float_complex*)h->vr, ldvr, (veclib_float_complex*)h->work, lwork, h->rwork);
#endif

    /* output */
    if(D!=NULL)
        memset(D, 0, dim*dim*sizeof(float_complex));
    
    /* failed to converge and find the eigenvalues */
    if( info != 0 ) {
        if(VL!=NULL)
            memset(VL, 0, dim*dim*sizeof(float_complex));
        if(VR!=NULL)
            memset(VR, 0, dim*dim*sizeof(float_complex));
        if(eig!=NULL)
            memset(eig, 0, dim*sizeof(float_complex));
#ifndef NDEBUG
        /* Failed to compute all of the eigenvalues, some eigenvectors have not
         * been computed, or the input matrix contained illegal values so no
         * solution was attempted. In these cases the function will zero all
         * output matrices and vectors. */
        saf_print_warning("Could not compute EVD in utility_ceig(). Output matrices/vectors have been zeroed.");
#endif
    }
    /* transpose, back to row-major */
    else{
        for(i=0; i<dim; i++){
            if(VL!=NULL)
                for(j=0; j<dim; j++)
                    VL[i*dim+j] = h->vl[j*dim+i];
            if(VR!=NULL)
                for(j=0; j<dim; j++)
                    VR[i*dim+j] = h->vr[j*dim+i];
            if(D!=NULL)
                D[i*dim+i] = h->w[i]; /* store along the diagonal */
            if(eig!=NULL)
                eig[i] = h->w[i];
        }
    }

    if(hWork == NULL)
        utility_ceig_destroy((void**)&h);
}

/** Data structure for utility_zeig() */
typedef struct _utility_zeig_data {
    int maxDim;
    veclib_int currentWorkSize;
    double_complex *w, *vl, *vr, *a;
    double* rwork;
    double_complex* work;
}utility_zeig_data;

void utility_zeig_create(void ** const phWork, int maxDim)
{
    *phWork = malloc1d(sizeof(utility_zeig_data));
    utility_zeig_data *h = (utility_zeig_data*)(*phWork);

    h->maxDim = maxDim;
    h->currentWorkSize = 0;
    h->rwork = malloc1d(4*maxDim*sizeof(double));
    h->w = malloc1d(maxDim*sizeof(double_complex));
    h->vl = malloc1d(maxDim*maxDim*sizeof(double_complex));
    h->vr = malloc1d(maxDim*maxDim*sizeof(double_complex));
    h->a = malloc1d(maxDim*maxDim*sizeof(double_complex));
    h->work = NULL;
}

void utility_zeig_destroy(void ** const phWork)
{
    utility_zeig_data *h = (utility_zeig_data*)(*phWork);

    if(h!=NULL){
        free(h->rwork);
        free(h->w);
        free(h->vl);
        free(h->vr);
        free(h->a);
        free(h->work);

        free(h);
        h=NULL;
        *phWork = NULL;
     }
}

void utility_zeig
(
    void* const hWork,
    const double_complex* A,
    const int dim,
    double_complex* VL,
    double_complex* VR,
    double_complex* D,
    double_complex* eig
)
{
    utility_zeig_data *h;
    veclib_int i, j, n, lda, ldvl, ldvr, info;
    veclib_int lwork;
    double_complex wkopt;

    n = lda = ldvl = ldvr = dim;

    /* Work struct */
    if(hWork==NULL)
        utility_zeig_create((void**)&h, dim);
    else{
        h = (utility_zeig_data*)(hWork);
#ifndef NDEBUG
        saf_assert(dim<=h->maxDim, "dim exceeds the maximum length specified");
#endif
    }

    /* store in column major order (i.e. transpose) */
    for(i=0; i<dim; i++)
        for(j=0; j<dim; j++)
            h->a[i*dim+j] = A[j*dim+i];

    /* Query how much "work" memory is required */
    lwork = -1;
#if defined(SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE)
    zgeev_( "Vectors", "Vectors", &n, (veclib_double_complex*)h->a, &lda, (veclib_double_complex*)h->w, (veclib_double_complex*)h->vl, &ldvl,
           (veclib_double_complex*)h->vr, &ldvr, (veclib_double_complex*)&wkopt, &lwork, h->rwork, &info );
#elif defined(SAF_VECLIB_USE_CLAPACK_INTERFACE)
    saf_print_error("No such implementation available in ATLAS CLAPACK");
#elif defined(SAF_VECLIB_USE_LAPACKE_INTERFACE)
    info = LAPACKE_zgeev_work(CblasColMajor, 'V', 'V', n, (veclib_double_complex*)h->a, lda, (veclib_double_complex*)h->w, (veclib_double_complex*)h->vl, ldvl,
                              (veclib_double_complex*)h->vr, ldvr, (veclib_double_complex*)&wkopt, lwork, h->rwork);
#endif
    lwork = (veclib_int)creal(wkopt);
    if(lwork>h->currentWorkSize){
        h->currentWorkSize = lwork;
        h->work = realloc1d(h->work, h->currentWorkSize*sizeof(double_complex));
    }

    /* solve the eigenproblem */
#if defined(SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE)
    zgeev_( "Vectors", "Vectors", &n, (veclib_double_complex*)h->a, &lda, (veclib_double_complex*)h->w, (veclib_double_complex*)h->vl, &ldvl,
           (veclib_double_complex*)h->vr, &ldvr, (veclib_double_complex*)h->work, &lwork, h->rwork, &info );
#elif defined(SAF_VECLIB_USE_CLAPACK_INTERFACE)
    saf_print_error("No such implementation available in ATLAS CLAPACK");
#elif defined(SAF_VECLIB_USE_LAPACKE_INTERFACE)
    info = LAPACKE_zgeev_work(CblasColMajor, 'V', 'V', n, (veclib_double_complex*)h->a, lda, (veclib_double_complex*)h->w, (veclib_double_complex*)h->vl, ldvl,
                              (veclib_double_complex*)h->vr, ldvr, (veclib_double_complex*)h->work, lwork, h->rwork);
#endif

    /* output */
    if(D!=NULL)
        memset(D, 0, dim*dim*sizeof(double_complex));

    /* failed to converge and find the eigenvalues */
    if( info != 0 ) {
        if(VL!=NULL)
            memset(VL, 0, dim*dim*sizeof(double_complex));
        if(VR!=NULL)
            memset(VR, 0, dim*dim*sizeof(double_complex));
        if(eig!=NULL)
            memset(eig, 0, dim*sizeof(double_complex));
#ifndef NDEBUG
        /* Failed to compute all of the eigenvalues, some eigenvectors have not
         * been computed, or the input matrix contained illegal values so no
         * solution was attempted. In these cases the function will zero all
         * output matrices and vectors. */
        saf_print_warning("Could not compute EVD in utility_zeig(). Output matrices/vectors have been zeroed.");
#endif
    }
    /* transpose, back to row-major */
    else{
        for(i=0; i<dim; i++){
            if(VL!=NULL)
                for(j=0; j<dim; j++)
                    VL[i*dim+j] = h->vl[j*dim+i];
            if(VR!=NULL)
                for(j=0; j<dim; j++)
                    VR[i*dim+j] = h->vr[j*dim+i];
            if(D!=NULL)
                D[i*dim+i] = h->w[i]; /* store along the diagonal */
            if(eig!=NULL)
                eig[i] = h->w[i];
        }
    }

    if(hWork == NULL)
        utility_zeig_destroy((void**)&h);
}


/* ========================================================================== */
/*                       General Linear Solver (?glslv)                       */
/* ========================================================================== */

/** Data structure for utility_sglslv() */
typedef struct _utility_sglslv_data {
    int maxDim;
    int maxNCol;
    veclib_int* IPIV;
    float* a;
    float* b;
}utility_sglslv_data;

void utility_sglslv_create(void ** const phWork, int maxDim, int maxNCol)
{
    *phWork = malloc1d(sizeof(utility_sglslv_data));
    utility_sglslv_data *h = (utility_sglslv_data*)(*phWork);

    h->maxDim = maxDim;
    h->maxNCol = maxNCol;

    h->IPIV = malloc1d(maxDim*sizeof(veclib_int));
    h->a = malloc1d(maxDim*maxDim*sizeof(float));
    h->b = malloc1d(maxDim*maxNCol*sizeof(float));
}

void utility_sglslv_destroy(void ** const phWork)
{
    utility_sglslv_data *h = (utility_sglslv_data*)(*phWork);

    if(h!=NULL){
        free(h->IPIV);
        free(h->a);
        free(h->b);

        free(h);
        h=NULL;
        *phWork = NULL;
     }
}

void utility_sglslv
(
    void* const hWork,
    const float* A,
    const int dim,
    float* B,
    int nCol,
    float* X
)
{
    utility_sglslv_data *h;
    veclib_int n = dim, nrhs = nCol, lda = dim, ldb = dim, info;
#if !defined(SAF_VECLIB_USE_LAPACKE_INTERFACE) && !(defined(SAF_USE_INTEL_MKL_LP64) || defined(SAF_USE_INTEL_MKL_ILP64))
    veclib_int i, j;
#endif

    /* Work struct */
    if(hWork==NULL)
        utility_sglslv_create((void**)&h, dim, nCol);
    else {
        h = (utility_sglslv_data*)(hWork);
#ifndef NDEBUG
        saf_assert(dim<=h->maxDim, "dim exceeds the maximum length specified");
        saf_assert(nCol<=h->maxNCol, "nCol exceeds the maximum length specified");
#endif
    }

#if defined(SAF_VECLIB_USE_LAPACKE_INTERFACE)
    /* Copy locally */
    cblas_scopy(dim*dim, A, 1, h->a, 1);
    cblas_scopy(dim*nCol, B, 1, h->b, 1);
#else
    /* store in column major order */
# if defined(SAF_USE_INTEL_MKL_LP64) || defined(SAF_USE_INTEL_MKL_ILP64)
    MKL_Somatcopy('R', 'T', dim, dim, 1.0f, A, dim, h->a, dim);
    MKL_Somatcopy('R', 'T', dim, nCol, 1.0f, B, nCol, h->b, dim);
# else
    for(i=0; i<dim; i++)
        for(j=0; j<dim; j++)
            h->a[i*dim+j] = A[j*dim+i];
    for(i=0; i<dim; i++)
        for(j=0; j<nCol; j++)
            h->b[j*dim+i] = B[i*nCol+j];
# endif
#endif
    
    /* solve Ax = b for each column in b (b is replaced by the solution: x) */
#ifdef SAF_VECLIB_USE_CLAPACK_INTERFACE
    info = clapack_sgesv(CblasColMajor, n, nrhs, h->a, lda, h->IPIV, h->b, ldb);
#elif defined(SAF_VECLIB_USE_LAPACKE_INTERFACE)
    info = LAPACKE_sgesv_work(CblasRowMajor, n, nrhs, h->a, lda, h->IPIV, h->b, ldb);
#elif defined(SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE)
    sgesv_( &n, &nrhs, h->a, &lda, h->IPIV, h->b, &ldb, &info );
#endif
    
    if(info!=0){
        /* A is singular, solution not possible */
        memset(X, 0, dim*nCol*sizeof(float));
#ifndef NDEBUG
        /* Input matrix was singular, solution was unable to be computed, or the
         * input matrix contained illegal values so no solution was attempted.
         * In these cases the function will zero the output matrix. */
        saf_print_warning("Could not solve the linear equation in utility_sglslv(). Output matrices/vectors have been zeroed.");
#endif
    }
    else{
#if defined(SAF_VECLIB_USE_LAPACKE_INTERFACE)
        cblas_scopy(dim*nCol, h->b, 1, X, 1);
#else
        /* store solution in row-major order */
# if defined(SAF_USE_INTEL_MKL_LP64) || defined(SAF_USE_INTEL_MKL_ILP64)
        MKL_Somatcopy('R', 'T', nCol, dim, 1.0f, h->b, dim, X, nCol);
# else
        for(i=0; i<dim; i++)
            for(j=0; j<nCol; j++)
                X[i*nCol+j] = h->b[j*dim+i];
# endif
#endif
    }
    
    if(hWork == NULL)
        utility_sglslv_destroy((void**)&h);
}

/** Data structure for utility_cglslv() */
typedef struct _utility_cglslv_data {
    int maxDim;
    int maxNCol;
    veclib_int* IPIV;
    float_complex* a;
    float_complex* b;
}utility_cglslv_data;

void utility_cglslv_create(void ** const phWork, int maxDim, int maxNCol)
{
    *phWork = malloc1d(sizeof(utility_cglslv_data));
    utility_cglslv_data *h = (utility_cglslv_data*)(*phWork);

    h->maxDim = maxDim;
    h->maxNCol = maxNCol;

    h->IPIV = malloc1d(maxDim*sizeof(veclib_int));
    h->a = malloc1d(maxDim*maxDim*sizeof(float_complex));
    h->b = malloc1d(maxDim*maxNCol*sizeof(float_complex));
}

void utility_cglslv_destroy(void ** const phWork)
{
    utility_cglslv_data *h = (utility_cglslv_data*)(*phWork);

    if(h!=NULL){
        free(h->IPIV);
        free(h->a);
        free(h->b);

        free(h);
        h=NULL;
        *phWork = NULL;
     }
}

void utility_cglslv
(
    void* const hWork,
    const float_complex* A,
    const int dim,
    float_complex* B,
    int nCol,
    float_complex* X
)
{
    utility_cglslv_data *h;
    veclib_int i, j, n, nrhs, lda, ldb, info;

    n = dim;
    nrhs = nCol;
    lda = dim;
    ldb = dim;

    /* Work struct */
    if(hWork==NULL)
        utility_cglslv_create((void**)&h, dim, nCol);
    else {
        h = (utility_cglslv_data*)(hWork);
#ifndef NDEBUG
        saf_assert(dim<=h->maxDim, "dim exceeds the maximum length specified");
        saf_assert(nCol<=h->maxNCol, "nCol exceeds the maximum length specified");
#endif
    }

    /* store in column major order */
    for(i=0; i<dim; i++)
        for(j=0; j<dim; j++)
            h->a[j*dim+i] = A[i*dim+j];
    for(i=0; i<dim; i++)
        for(j=0; j<nCol; j++)
            h->b[j*dim+i] = B[i*nCol+j];
    
    /* solve Ax = b for each column in b (b is replaced by the solution: x) */
#if defined(SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE)
    cgesv_( &n, &nrhs, (veclib_float_complex*)h->a, &lda, h->IPIV, (veclib_float_complex*)h->b, &ldb, &info );
#elif defined(SAF_VECLIB_USE_CLAPACK_INTERFACE)
    info = clapack_cgesv(CblasColMajor, n, nrhs, (veclib_float_complex*)h->a, lda, h->IPIV, (veclib_float_complex*)h->b, ldb);
#elif defined(SAF_VECLIB_USE_LAPACKE_INTERFACE)
    info = LAPACKE_cgesv_work(CblasColMajor, n, nrhs, (veclib_float_complex*)h->a, lda, h->IPIV, (veclib_float_complex*)h->b, ldb);
#endif
    
    /* A is singular, solution not possible */
    if(info!=0){
        memset(X, 0, dim*nCol*sizeof(float_complex));
#ifndef NDEBUG
        /* Input matrix was singular, solution was unable to be computed, or the
         * input matrix contained illegal values so no solution was attempted.
         * In these cases the function will zero the output matrix. */
        saf_print_warning("Could not solve the linear equation in utility_cglslv(). Output matrices/vectors have been zeroed.");
#endif
    }
    /* store solution in row-major order */
    else{
        for(i=0; i<dim; i++)
            for(j=0; j<nCol; j++)
                X[i*nCol+j] = h->b[j*dim+i];
    }

    if(hWork == NULL)
        utility_cglslv_destroy((void**)&h);
}

/** Data structure for utility_dglslv() */
typedef struct _utility_dglslv_data {
    int maxDim;
    int maxNCol;
    veclib_int* IPIV;
    double* a;
    double* b;
}utility_dglslv_data;

void utility_dglslv_create(void ** const phWork, int maxDim, int maxNCol)
{
    *phWork = malloc1d(sizeof(utility_dglslv_data));
    utility_dglslv_data *h = (utility_dglslv_data*)(*phWork);

    h->maxDim = maxDim;
    h->maxNCol = maxNCol;
    h->IPIV = malloc1d(maxDim*sizeof(veclib_int));
    h->a = malloc1d(maxDim*maxDim*sizeof(double));
    h->b = malloc1d(maxDim*maxNCol*sizeof(double));
}

void utility_dglslv_destroy(void ** const phWork)
{
    utility_dglslv_data *h = (utility_dglslv_data*)(*phWork);

    if(h!=NULL){
        free(h->IPIV);
        free(h->a);
        free(h->b);

        free(h);
        h=NULL;
        *phWork = NULL;
     }
}

void utility_dglslv
(
    void* const hWork,
    const double* A,
    const int dim,
    double* B,
    int nCol,
    double* X
)
{
    utility_dglslv_data *h;
    veclib_int i, j, n, nrhs, lda, ldb, info;

    n = dim;
    nrhs = nCol;
    lda = dim;
    ldb = dim;

    /* Work struct */
    if(hWork==NULL)
        utility_dglslv_create((void**)&h, dim, nCol);
    else {
        h = (utility_dglslv_data*)(hWork);
#ifndef NDEBUG
        saf_assert(dim<=h->maxDim, "dim exceeds the maximum length specified");
        saf_assert(nCol<=h->maxNCol, "nCol exceeds the maximum length specified");
#endif
    }

    /* store in column major order */
    for(i=0; i<dim; i++)
        for(j=0; j<dim; j++)
            h->a[j*dim+i] = A[i*dim+j];
    for(i=0; i<dim; i++)
        for(j=0; j<nCol; j++)
            h->b[j*dim+i] = B[i*nCol+j];
    
    /* solve Ax = b for each column in b (b is replaced by the solution: x) */
#ifdef SAF_VECLIB_USE_CLAPACK_INTERFACE
    info = clapack_dgesv(CblasColMajor, n, nrhs, h->a, lda, h->IPIV, h->b, ldb);
#elif defined(SAF_VECLIB_USE_LAPACKE_INTERFACE)
    info = LAPACKE_dgesv_work(CblasColMajor, n, nrhs, h->a, lda, h->IPIV, h->b, ldb);
#elif defined(SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE)
    dgesv_( &n, &nrhs, h->a, &lda, h->IPIV, h->b, &ldb, &info );
#endif
    
    /* A is singular, solution not possible */
    if(info!=0){
        memset(X, 0, dim*nCol*sizeof(double));
#ifndef NDEBUG
        /* Input matrix was singular, solution was unable to be computed, or the
         * input matrix contained illegal values so no solution was attempted.
         * In these cases the function will zero the output matrix. */
        saf_print_warning("Could not solve the linear equation in utility_dglslv(). Output matrices/vectors have been zeroed.");
#endif
    }
    /* store solution in row-major order */
    else{
        for(i=0; i<dim; i++)
            for(j=0; j<nCol; j++)
                X[i*nCol+j] = h->b[j*dim+i];
    }

    if(hWork == NULL)
        utility_dglslv_destroy((void**)&h);
}

/** Data structure for utility_zglslv() */
typedef struct _utility_zglslv_data {
    int maxDim;
    int maxNCol;
    veclib_int* IPIV;
    double_complex* a;
    double_complex* b;
}utility_zglslv_data;

void utility_zglslv_create(void ** const phWork, int maxDim, int maxNCol)
{
    *phWork = malloc1d(sizeof(utility_zglslv_data));
    utility_zglslv_data *h = (utility_zglslv_data*)(*phWork);

    h->maxDim = maxDim;
    h->maxNCol = maxNCol; 
    h->IPIV = malloc1d(maxDim*sizeof(veclib_int));
    h->a = malloc1d(maxDim*maxDim*sizeof(double_complex));
    h->b = malloc1d(maxDim*maxNCol*sizeof(double_complex));
}

void utility_zglslv_destroy(void ** const phWork)
{
    utility_zglslv_data *h = (utility_zglslv_data*)(*phWork);

    if(h!=NULL){
        free(h->IPIV);
        free(h->a);
        free(h->b);

        free(h);
        h=NULL;
        *phWork = NULL;
     }
}

void utility_zglslv
(
    void* const hWork,
    const double_complex* A,
    const int dim,
    double_complex* B,
    int nCol,
    double_complex* X
)
{
    utility_zglslv_data *h;
    veclib_int i, j, n = dim, nrhs = nCol, lda = dim, ldb = dim, info;

    n = dim;
    nrhs = nCol;
    lda = dim;
    ldb = dim;

    /* Work struct */
    if(hWork==NULL)
        utility_zglslv_create((void**)&h, dim, nCol);
    else {
        h = (utility_zglslv_data*)(hWork);
#ifndef NDEBUG
        saf_assert(dim<=h->maxDim, "dim exceeds the maximum length specified");
        saf_assert(nCol<=h->maxNCol, "nCol exceeds the maximum length specified");
#endif
    }
    
    /* store in column major order */
    for(i=0; i<dim; i++)
        for(j=0; j<dim; j++)
            h->a[j*dim+i] = A[i*dim+j];
    for(i=0; i<dim; i++)
        for(j=0; j<nCol; j++)
            h->b[j*dim+i] = B[i*nCol+j];
    
    /* solve Ax = b for each column in b (b is replaced by the solution: x) */
#if defined(SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE)
    zgesv_( &n, &nrhs, (veclib_double_complex*)h->a, &lda, h->IPIV, (veclib_double_complex*)h->b, &ldb, &info );
#elif defined(SAF_VECLIB_USE_CLAPACK_INTERFACE)
    info = clapack_zgesv(CblasColMajor, n, nrhs, (veclib_double_complex*)h->a, lda, h->IPIV, (veclib_double_complex*)h->b, ldb);
#elif defined(SAF_VECLIB_USE_LAPACKE_INTERFACE)
    info = LAPACKE_zgesv_work(CblasColMajor, n, nrhs, (veclib_double_complex*)h->a, lda, h->IPIV, (veclib_double_complex*)h->b, ldb);
#endif
    
    /* A is singular, solution not possible */
    if(info!=0){
        memset(X, 0, dim*nCol*sizeof(double_complex));
#ifndef NDEBUG
        /* Input matrix was singular, solution was unable to be computed, or the
         * input matrix contained illegal values so no solution was attempted.
         * In these cases the function will zero the output matrix. */
        saf_print_warning("Could not solve the linear equation in utility_zglslv(). Output matrices/vectors have been zeroed.");
#endif
    }
    /* store solution in row-major order */
    else{
        for(i=0; i<dim; i++)
            for(j=0; j<nCol; j++)
                X[i*nCol+j] = h->b[j*dim+i];
    }

    if(hWork == NULL)
        utility_zglslv_destroy((void**)&h);
}


/* ========================================================================== */
/*                      General Linear Solver (?glslvt)                       */
/* ========================================================================== */

/** Data structure for utility_sglslvt() */
typedef struct _utility_sglslvt_data {
    int maxDim;
    int maxNCol;
    veclib_int* IPIV;
    float* a;
    float* b;
}utility_sglslvt_data;

void utility_sglslvt_create(void ** const phWork, int maxDim, int maxNCol)
{
    *phWork = malloc1d(sizeof(utility_sglslvt_data));
    utility_sglslvt_data *h = (utility_sglslvt_data*)(*phWork);

    h->maxDim = maxDim;
    h->maxNCol = maxNCol;
    h->IPIV = malloc1d(maxDim*sizeof(veclib_int));
    h->a = malloc1d(maxDim*maxDim*sizeof(float));
    h->b = malloc1d(maxDim*maxNCol*sizeof(float));
}

void utility_sglslvt_destroy(void ** const phWork)
{
    utility_sglslvt_data *h = (utility_sglslvt_data*)(*phWork);

    if(h!=NULL){
        free(h->IPIV);
        free(h->a);
        free(h->b);

        free(h);
        h=NULL;
        *phWork = NULL;
     }
}

void utility_sglslvt
(
    void* const hWork,
    const float* A,
    const int dim,
    float* B,
    int nCol,
    float* X
)
{
    utility_sglslvt_data *h;
    veclib_int n = nCol, nrhs = dim, lda = nCol, ldb = nCol, info;

    /* Work struct */
    if(hWork==NULL)
        utility_sglslvt_create((void**)&h, dim, nCol);
    else {
        h = (utility_sglslvt_data*)(hWork);
#ifndef NDEBUG
        saf_assert(dim<=h->maxDim, "dim exceeds the maximum length specified");
        saf_assert(nCol<=h->maxNCol, "nCol exceeds the maximum length specified");
#endif
    }

    /* store locally */
    cblas_scopy(dim*dim, A, 1, h->a, 1);
    cblas_scopy(dim*nCol, B, 1, h->b, 1);

    /* solve Ax = b for each column in b (b is replaced by the solution: x) */
#ifdef SAF_VECLIB_USE_CLAPACK_INTERFACE
    info = clapack_sgesv(CblasColMajor, n, nrhs, h->b, ldb, h->IPIV, h->a, lda);
#elif defined(SAF_VECLIB_USE_LAPACKE_INTERFACE)
    info = LAPACKE_sgesv_work(CblasColMajor, n, nrhs, h->b, ldb, h->IPIV, h->a, lda);
#elif defined(SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE)
    //sposv_("U", &n, &nrhs, h->b, &ldb, h->a, &lda, &info );
    sgesv_( &n, &nrhs, h->b, &ldb, h->IPIV, h->a, &lda, &info );
#endif

    if(info!=0){
        /* A is singular, solution not possible */
        memset(X, 0, dim*nCol*sizeof(float));
#ifndef NDEBUG
        /* Input matrix was singular, solution was unable to be computed, or the
         * input matrix contained illegal values so no solution was attempted.
         * In these cases the function will zero the output matrix. */
        saf_print_warning("Could not solve the linear equation in utility_sglslvt(). Output matrices/vectors have been zeroed.");
#endif
    }
    else
        cblas_scopy(dim*nCol, h->a, 1, X, 1);

    if(hWork == NULL)
        utility_sglslvt_destroy((void**)&h);
}


/* ========================================================================== */
/*                      Symmetric Linear Solver (?slslv)                      */
/* ========================================================================== */

/** Data structure for utility_sslslv() */
typedef struct _utility_sslslv_data {
    int maxDim;
    int maxNCol;
    float* a;
    float* b;
}utility_sslslv_data;

void utility_sslslv_create(void ** const phWork, int maxDim, int maxNCol)
{
    *phWork = malloc1d(sizeof(utility_sslslv_data));
    utility_sslslv_data *h = (utility_sslslv_data*)(*phWork);

    h->maxDim = maxDim;
    h->maxNCol = maxNCol;
    h->a = malloc1d(maxDim*maxDim*sizeof(float));
    h->b = malloc1d(maxDim*maxNCol*sizeof(float));
}

void utility_sslslv_destroy(void ** const phWork)
{
    utility_sslslv_data *h = (utility_sslslv_data*)(*phWork);

    if(h!=NULL){
        free(h->a);
        free(h->b);

        free(h);
        h=NULL;
        *phWork = NULL;
     }
}

void utility_sslslv
(
    void* const hWork,
    const float* A,
    const int dim,
    float* B,
    int nCol,
    float* X
)
{
    utility_sslslv_data *h;
    veclib_int i, j, n, nrhs, lda, ldb, info;

    n = dim;
    nrhs = nCol;
    lda = dim;
    ldb = dim;

    /* Work struct */
    if(hWork==NULL)
        utility_sslslv_create((void**)&h, dim, nCol);
    else {
        h = (utility_sslslv_data*)(hWork);
#ifndef NDEBUG
        saf_assert(dim<=h->maxDim, "dim exceeds the maximum length specified");
        saf_assert(nCol<=h->maxNCol, "nCol exceeds the maximum length specified");
#endif
    }
    
    /* store in column major order */
    for(i=0; i<dim; i++)
        for(j=0; j<dim; j++)
            h->a[j*dim+i] = A[i*dim+j];
    for(i=0; i<dim; i++)
        for(j=0; j<nCol; j++)
            h->b[j*dim+i] = B[i*nCol+j];
    
    /* solve Ax = b for each column in b (b is replaced by the solution: x) */
#ifdef SAF_VECLIB_USE_CLAPACK_INTERFACE
    info = clapack_sposv(CblasColMajor, CblasUpper, n, nrhs, h->a, lda, h->b, ldb);
#elif defined(SAF_VECLIB_USE_LAPACKE_INTERFACE)
    info = LAPACKE_sposv_work(CblasColMajor, CblasUpper, n, nrhs, h->a, lda, h->b, ldb);
#elif defined(SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE)
    sposv_( "U", &n, &nrhs, h->a, &lda, h->b, &ldb, &info );
#endif
    
    /* A is not symmetric positive definate, solution not possible */
    if(info!=0){
        memset(X, 0, dim*nCol*sizeof(float));
#ifndef NDEBUG
        /* Input matrix was singular, solution was unable to be computed, or the
         * input matrix contained illegal values so no solution was attempted.
         * In these cases the function will zero the output matrix. */
        saf_print_warning("Could not solve the linear equation in utility_sslslv(). Output matrices/vectors have been zeroed.");
#endif
    }
    /* store solution in row-major order */
    else{
        for(i=0; i<dim; i++)
            for(j=0; j<nCol; j++)
                X[i*nCol+j] = h->b[j*dim+i];
    }
    
    if(hWork == NULL)
        utility_sslslv_destroy((void**)&h);
}

/** Data structure for utility_cslslv() */
typedef struct _utility_cslslv_data {
    int maxDim;
    int maxNCol;
    float_complex* a;
    float_complex* b;
}utility_cslslv_data;

void utility_cslslv_create(void ** const phWork, int maxDim, int maxNCol)
{
    *phWork = malloc1d(sizeof(utility_cslslv_data));
    utility_cslslv_data *h = (utility_cslslv_data*)(*phWork);

    h->maxDim = maxDim;
    h->maxNCol = maxNCol;
    h->a = malloc1d(maxDim*maxDim*sizeof(float_complex));
    h->b = malloc1d(maxDim*maxNCol*sizeof(float_complex));
}

void utility_cslslv_destroy(void ** const phWork)
{
    utility_cslslv_data *h = (utility_cslslv_data*)(*phWork);

    if(h!=NULL){
        free(h->a);
        free(h->b);

        free(h);
        h=NULL;
        *phWork = NULL;
     }
}

void utility_cslslv
(
    void* const hWork,
    const float_complex* A,
    const int dim,
    float_complex* B,
    int nCol,
    float_complex* X
)
{
    utility_cslslv_data *h;
    veclib_int i, j, n, nrhs, lda, ldb, info;

    n = dim;
    nrhs = nCol;
    lda = dim;
    ldb = dim;

    /* Work struct */
    if(hWork==NULL)
        utility_cslslv_create((void**)&h, dim, nCol);
    else {
        h = (utility_cslslv_data*)(hWork);
#ifndef NDEBUG
        saf_assert(dim<=h->maxDim, "dim exceeds the maximum length specified");
        saf_assert(nCol<=h->maxNCol, "nCol exceeds the maximum length specified");
#endif
    }

    /* store in column major order */
    for(i=0; i<dim; i++)
        for(j=0; j<dim; j++)
            h->a[j*dim+i] = A[i*dim+j];
    for(i=0; i<dim; i++)
        for(j=0; j<nCol; j++)
    h->b[j*dim+i] = B[i*nCol+j];
    
    /* solve Ax = b for each column in b (b is replaced by the solution: x) */
#if defined(SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE)
    cposv_( "U", &n, &nrhs, (veclib_float_complex*)h->a, &lda, (veclib_float_complex*)h->b, &ldb, &info );
#elif defined(SAF_VECLIB_USE_CLAPACK_INTERFACE)
    info = clapack_cposv(CblasColMajor, CblasUpper, n, nrhs, (veclib_float_complex*)h->a, lda, (veclib_float_complex*)h->b, ldb);
#elif defined(SAF_VECLIB_USE_LAPACKE_INTERFACE)
    info = LAPACKE_cposv_work(CblasColMajor, CblasUpper, n, nrhs, (veclib_float_complex*)h->a, lda, (veclib_float_complex*)h->b, ldb);
#endif
    
    /* A is not symmetric positive definate, solution not possible */
    if(info!=0){
        memset(X, 0, dim*nCol*sizeof(float_complex));
#ifndef NDEBUG
        /* Input matrix was singular, solution was unable to be computed, or the
         * input matrix contained illegal values so no solution was attempted.
         * In these cases the function will zero the output matrix. */
        saf_print_warning("Could not solve the linear equation in utility_cslslv(). Output matrices/vectors have been zeroed.");
#endif
    }
    /* store solution in row-major order */
    else{
        for(i=0; i<dim; i++)
            for(j=0; j<nCol; j++)
                X[i*nCol+j] = h->b[j*dim+i];
    }
    
    if(hWork == NULL)
        utility_cslslv_destroy((void**)&h);
}


/* ========================================================================== */
/*                        Matrix Pseudo-Inverse (?pinv)                       */
/* ========================================================================== */

/** Data structure for utility_spinv() */
typedef struct _utility_spinv_data {
    int maxDim1, maxDim2;
    veclib_int currentWorkSize;
    float* a, *s, *u, *vt, *inva;
    float* work;
}utility_spinv_data;

void utility_spinv_create(void ** const phWork, int maxDim1, int maxDim2)
{
    *phWork = malloc1d(sizeof(utility_spinv_data));
    utility_spinv_data *h = (utility_spinv_data*)(*phWork);

    h->maxDim1 = maxDim1;
    h->maxDim2 = maxDim2;
    h->currentWorkSize = 0;
    h->a = malloc1d(maxDim1*maxDim2*sizeof(float));
    h->s = malloc1d(SAF_MIN(maxDim2,maxDim1)*sizeof(float));
    h->u = malloc1d(maxDim1*maxDim1*sizeof(float));
    h->vt = malloc1d(maxDim2*maxDim2*sizeof(float));
    h->inva = malloc1d(maxDim1*maxDim2*sizeof(float));
    h->work = NULL;
}

void utility_spinv_destroy(void ** const phWork)
{
    utility_spinv_data *h = (utility_spinv_data*)(*phWork);

    if(h!=NULL){
        free(h->a);
        free(h->s);
        free(h->u);
        free(h->vt);
        free(h->inva);
        free(h->work);

        free(h);
        h=NULL;
        *phWork = NULL;
     }
}

void utility_spinv
(
    void* const hWork,
    const float* inM,
    const int dim1,
    const int dim2,
    float* outM
)
{
    utility_spinv_data *h;
    veclib_int i, j, m, n, k, lda, ldu, ldvt, ld_inva, info;
    float ss;
    veclib_int lwork;
    float wkopt;
    
    m = lda = ldu = dim1;
    n = dim2;
    k = ldvt = m < n ? m : n;

    /* Work struct */
    if(hWork==NULL)
        utility_spinv_create((void**)&h, dim1, dim2);
    else{
        h = (utility_spinv_data*)(hWork);
#ifndef NDEBUG
        saf_assert(dim1<=h->maxDim1, "dim1 exceeds the maximum length specified");
        saf_assert(dim2<=h->maxDim2, "dim2 exceeds the maximum length specified");
#endif
    }
    
    /* store in column major order */
    for(i=0; i<m; i++)
        for(j=0; j<n; j++)
            h->a[j*m+i] = inM[i*n+j];

    /* Query how much "work" memory is required */
    lwork = -1;
#if defined(SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE)
    sgesvd_("S", "S", &m, &n, h->a, &lda, h->s, h->u, &ldu, h->vt, &ldvt, &wkopt, &lwork, &info);
#elif defined(SAF_VECLIB_USE_CLAPACK_INTERFACE)
    saf_print_error("No such implementation available in ATLAS CLAPACK");
#elif defined(SAF_VECLIB_USE_LAPACKE_INTERFACE)
    info = LAPACKE_sgesvd_work(CblasColMajor, 'S', 'S', m, n, h->a, lda, h->s, h->u, ldu, h->vt, ldvt, &wkopt, lwork);
#endif
    lwork = (veclib_int)wkopt;
    if(lwork>h->currentWorkSize){
        h->currentWorkSize = lwork;
        h->work = realloc1d(h->work, h->currentWorkSize*sizeof(float));
    }

    /* Perform the singular value decomposition */
#if defined(SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE)
    sgesvd_("S", "S", &m, &n, h->a, &lda, h->s, h->u, &ldu, h->vt, &ldvt, h->work, &lwork, &info );
#elif defined(SAF_VECLIB_USE_CLAPACK_INTERFACE)
    saf_print_error("No such implementation available in ATLAS CLAPACK");
#elif defined(SAF_VECLIB_USE_LAPACKE_INTERFACE)
    info = LAPACKE_sgesvd_work(CblasColMajor, 'S', 'S', m, n, h->a, lda, h->s, h->u, ldu, h->vt, ldvt, h->work, lwork);
#endif
    
    if( info != 0 ) {
        memset(outM, 0, dim1*dim2*sizeof(float));
#ifndef NDEBUG
        /* The SVD failed to converge, or the input matrix contained illegal
         * values so no solution was attempted. In these cases this function
         * will zero all output matrices and/or vectors. */
        saf_print_warning("Could not compute SVD in utility_spinv(). Output matrices/vectors have been zeroed.");
#endif
    }
    else{
        for(i=0; i<k; i++){
            if(h->s[i] > 1.0e-5f)
                ss=1.0f/h->s[i];
            else
                ss=h->s[i];
            cblas_sscal(m, ss, &h->u[i*m], 1);
        }
        ld_inva=n;
        cblas_sgemm(CblasColMajor, CblasTrans, CblasTrans, n, m, k, 1.0f,
                    h->vt, ldvt,
                    h->u, ldu, 0.0f,
                    h->inva, ld_inva);

        /* return in row-major order */
        for(i=0; i<m; i++)
            for(j=0; j<n; j++)
                outM[j*m+i] = h->inva[i*n+j];
    }
    
    /* clean-up */
    if(hWork == NULL)
        utility_spinv_destroy((void**)&h);
}

/** Data structure for utility_cpinv() */
typedef struct _utility_cpinv_data {
    int maxDim1, maxDim2;
    veclib_int currentWorkSize;
    float_complex* a, *u, *vt, *inva;
    float* s, *rwork;
    float_complex* work;
}utility_cpinv_data;

void utility_cpinv_create(void ** const phWork, int maxDim1, int maxDim2)
{
    *phWork = malloc1d(sizeof(utility_cpinv_data));
    utility_cpinv_data *h = (utility_cpinv_data*)(*phWork);

    h->maxDim1 = maxDim1;
    h->maxDim2 = maxDim2;
    h->currentWorkSize = 0;
    h->a = malloc1d(maxDim1*maxDim2*sizeof(float_complex));
    h->s = malloc1d(SAF_MIN(maxDim2,maxDim1)*sizeof(float));
    h->u = malloc1d(maxDim1*maxDim1*sizeof(float_complex));
    h->vt = malloc1d(maxDim2*maxDim2*sizeof(float_complex));
    h->inva = malloc1d(maxDim1*maxDim2*sizeof(float_complex));
    h->rwork = malloc1d(maxDim1*SAF_MAX(1, 5*SAF_MIN(maxDim2,maxDim1))*sizeof(float));
    h->work = NULL;
}

void utility_cpinv_destroy(void ** const phWork)
{
    utility_cpinv_data *h = (utility_cpinv_data*)(*phWork);

    if(h!=NULL){
        free(h->a);
        free(h->s);
        free(h->u);
        free(h->vt);
        free(h->inva);
        free(h->work);

        free(h);
        h=NULL;
        *phWork = NULL;
     }
}

void utility_cpinv
(
    void* const hWork,
    const float_complex* inM,
    const int dim1,
    const int dim2,
    float_complex* outM
)
{
    utility_cpinv_data *h;
    veclib_int i, j, m, n, k, lda, ldu, ldvt, ld_inva, info;
    float_complex  ss_cmplx;
    const float_complex calpha = cmplxf(1.0f, 0.0f); const float_complex cbeta = cmplxf(0.0f, 0.0f); /* blas */
    float ss;
    veclib_int lwork;
    float_complex wkopt;
    
    m = lda = ldu = dim1;
    n = dim2;
    k = ldvt = m < n ? m : n;

    /* Work struct */
    if(hWork==NULL)
        utility_cpinv_create((void**)&h, dim1, dim2);
    else{
        h = (utility_cpinv_data*)(hWork);
#ifndef NDEBUG
        saf_assert(dim1<=h->maxDim1, "dim1 exceeds the maximum length specified");
        saf_assert(dim2<=h->maxDim2, "dim2 exceeds the maximum length specified");
#endif
    }
    
    /* store in column major order */
    for(i=0; i<dim1; i++)
        for(j=0; j<dim2; j++)
            h->a[j*dim1+i] = inM[i*dim2 +j];

    /* Query how much "work" memory is required */
    lwork = -1;
#if defined(SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE)
    cgesvd_( "A", "A", &m, &n, (veclib_float_complex*)h->a, &lda, h->s, (veclib_float_complex*)h->u, &ldu, (veclib_float_complex*)h->vt, &ldvt,
            (veclib_float_complex*)&wkopt, &lwork, h->rwork, &info );
#elif defined(SAF_VECLIB_USE_CLAPACK_INTERFACE)
    saf_print_error("No such implementation available in ATLAS CLAPACK");
#elif defined(SAF_VECLIB_USE_LAPACKE_INTERFACE)
    info = LAPACKE_cgesvd_work(CblasColMajor, 'S', 'S', m, n, (veclib_float_complex*)h->a, lda, h->s, (veclib_float_complex*)h->u, ldu, (veclib_float_complex*)h->vt, ldvt,
                               (veclib_float_complex*)&wkopt, lwork, h->rwork);
#endif
    lwork = (veclib_int)(crealf(wkopt)+0.01f);
    if(lwork>h->currentWorkSize){
        h->currentWorkSize = lwork;
        h->work = realloc1d(h->work, h->currentWorkSize*sizeof(float_complex));
    }

    /* Perform the singular value decomposition */
#if defined(SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE)
    cgesvd_( "A", "A", &m, &n, (veclib_float_complex*)h->a, &lda, h->s, (veclib_float_complex*)h->u, &ldu, (veclib_float_complex*)h->vt, &ldvt,
            (veclib_float_complex*)h->work, &lwork, h->rwork, &info);
#elif defined(SAF_VECLIB_USE_CLAPACK_INTERFACE)
    saf_print_error("No such implementation available in ATLAS CLAPACK");
#elif defined(SAF_VECLIB_USE_LAPACKE_INTERFACE)
    info = LAPACKE_cgesvd_work(CblasColMajor, 'S', 'S', m, n, (veclib_float_complex*)h->a, lda, h->s, (veclib_float_complex*)h->u, ldu, (veclib_float_complex*)h->vt, ldvt,
                               (veclib_float_complex*)h->work, lwork, h->rwork);
#endif
    
    if( info != 0 ) {
        memset(outM, 0, dim1*dim2*sizeof(float_complex));
#ifndef NDEBUG
        /* The SVD failed to converge, or the input matrix contained illegal
         * values so no solution was attempted. In these cases this function
         * will zero all output matrices and/or vectors. */
        saf_print_warning("Could not compute SVD in utility_cpinv(). Output matrices/vectors have been zeroed.");
#endif
    }
    else{
        for(i=0; i<k; i++){
            if(h->s[i] > 1.0e-5f)
                ss=1.0f/h->s[i];
            else
                ss=h->s[i];
            ss_cmplx = cmplxf(ss, 0.0f);
            cblas_cscal(m, &ss_cmplx, &h->u[i*m], 1);
        }
        ld_inva=n;
        cblas_cgemm(CblasColMajor, CblasConjTrans, CblasConjTrans, n, m, k, &calpha,
                    h->vt, ldvt,
                    h->u, ldu, &cbeta,
                    h->inva, ld_inva);

        /* return in row-major order */
        for(i=0; i<m; i++)
            for(j=0; j<n; j++)
                outM[j*m+i] = h->inva[i*n+j];
    }

    /* clean-up */
    if(hWork == NULL)
        utility_cpinv_destroy((void**)&h);
}

/** Data structure for utility_dpinv() */
typedef struct _utility_dpinv_data {
    int maxDim1, maxDim2;
    veclib_int currentWorkSize;
    double* a, *s, *u, *vt, *inva;
    double* work;
}utility_dpinv_data;

void utility_dpinv_create(void ** const phWork, int maxDim1, int maxDim2)
{
    *phWork = malloc1d(sizeof(utility_dpinv_data));
    utility_dpinv_data *h = (utility_dpinv_data*)(*phWork);

    h->maxDim1 = maxDim1;
    h->maxDim2 = maxDim2;
    h->currentWorkSize = 0;
    h->a = malloc1d(maxDim1*maxDim2*sizeof(double));
    h->s = malloc1d(SAF_MIN(maxDim2,maxDim1)*sizeof(double));
    h->u = malloc1d(maxDim1*maxDim1*sizeof(double));
    h->vt = malloc1d(maxDim2*maxDim2*sizeof(double));
    h->inva = malloc1d(maxDim1*maxDim2*sizeof(double));
    h->work = NULL;
}

void utility_dpinv_destroy(void ** const phWork)
{
    utility_dpinv_data *h = (utility_dpinv_data*)(*phWork);

    if(h!=NULL){
        free(h->a);
        free(h->s);
        free(h->u);
        free(h->vt);
        free(h->inva);
        free(h->work);

        free(h);
        h=NULL;
        *phWork = NULL;
     }
}

void utility_dpinv
(
    void* const hWork,
    const double* inM,
    const int dim1,
    const int dim2,
    double* outM
)
{
    utility_dpinv_data *h;
    veclib_int i, j, m, n, k, lda, ldu, ldvt, ld_inva, info;
    double ss;
    veclib_int lwork;
    double wkopt;
    
    m = lda = ldu = dim1;
    n = dim2;
    k = ldvt = m < n ? m : n;

    /* Work struct */
    if(hWork==NULL)
        utility_dpinv_create((void**)&h, dim1, dim2);
    else{
        h = (utility_dpinv_data*)(hWork);
#ifndef NDEBUG
        saf_assert(dim1<=h->maxDim1, "dim1 exceeds the maximum length specified");
        saf_assert(dim2<=h->maxDim2, "dim2 exceeds the maximum length specified");
#endif
    }
    
    /* store in column major order */
    for(i=0; i<m; i++)
        for(j=0; j<n; j++)
            h->a[j*m+i] = inM[i*n+j];

    /* Query how much "work" memory is required */
    lwork = -1;
#if defined(SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE)
    dgesvd_("S", "S", &m, &n, h->a, &lda, h->s, h->u, &ldu, h->vt, &ldvt, &wkopt, &lwork, &info);
#elif defined(SAF_VECLIB_USE_CLAPACK_INTERFACE)
    saf_print_error("No such implementation available in ATLAS CLAPACK");
#elif defined(SAF_VECLIB_USE_LAPACKE_INTERFACE)
    info = LAPACKE_dgesvd_work(CblasColMajor, 'S', 'S', m, n, h->a, lda, h->s, h->u, ldu, h->vt, ldvt, &wkopt, lwork);
#endif
    lwork = (veclib_int)wkopt;
    if(lwork>h->currentWorkSize){
        h->currentWorkSize = lwork;
        h->work = realloc1d(h->work, h->currentWorkSize*sizeof(double));
    }

    /* Perform the singular value decomposition */
#if defined(SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE)
    dgesvd_("S", "S", &m, &n, h->a, &lda, h->s, h->u, &ldu, h->vt, &ldvt, h->work, &lwork, &info );
#elif defined(SAF_VECLIB_USE_CLAPACK_INTERFACE)
    saf_print_error("No such implementation available in ATLAS CLAPACK");
#elif defined(SAF_VECLIB_USE_LAPACKE_INTERFACE)
    info = LAPACKE_dgesvd_work(CblasColMajor, 'S', 'S', m, n, h->a, lda, h->s, h->u, ldu, h->vt, ldvt, h->work, lwork);
#endif
    
    if( info != 0 ) {
        memset(outM, 0, dim1*dim2*sizeof(double));
#ifndef NDEBUG
        /* The SVD failed to converge, or the input matrix contained illegal
         * values so no solution was attempted. In these cases this function
         * will zero all output matrices and/or vectors. */
        saf_print_warning("Could not compute SVD in utility_dpinv(). Output matrices/vectors have been zeroed.");
#endif
    }
    else{
        for(i=0; i<k; i++){
            if(h->s[i] > 1.0e-9)
                ss=1.0/h->s[i];
            else
                ss=h->s[i];
            cblas_dscal(m, ss, &h->u[i*m], 1);
        }
        ld_inva=n;
        cblas_dgemm(CblasColMajor, CblasTrans, CblasTrans, n, m, k, 1.0,
                    h->vt, ldvt,
                    h->u, ldu, 0.0,
                    h->inva, ld_inva);

        /* return in row-major order */
        for(i=0; i<m; i++)
            for(j=0; j<n; j++)
                outM[j*m+i] = h->inva[i*n+j];
    }

    /* clean-up */
    if(hWork == NULL)
        utility_dpinv_destroy((void**)&h);
}

/** Data structure for utility_zpinv() */
typedef struct _utility_zpinv_data {
    int maxDim1, maxDim2;
    veclib_int currentWorkSize;
    double_complex* a, *u, *vt, *inva;
    double* s, *rwork;
    double_complex* work;
}utility_zpinv_data;

void utility_zpinv_create(void ** const phWork, int maxDim1, int maxDim2)
{
    *phWork = malloc1d(sizeof(utility_zpinv_data));
    utility_zpinv_data *h = (utility_zpinv_data*)(*phWork);

    h->maxDim1 = maxDim1;
    h->maxDim2 = maxDim2;
    h->currentWorkSize = 0;
    h->a = malloc1d(maxDim1*maxDim2*sizeof(double_complex));
    h->s = malloc1d(SAF_MIN(maxDim2,maxDim1)*sizeof(double));
    h->u = malloc1d(maxDim1*maxDim1*sizeof(double_complex));
    h->vt = malloc1d(maxDim2*maxDim2*sizeof(double_complex));
    h->inva = malloc1d(maxDim1*maxDim2*sizeof(double_complex));
    h->rwork = malloc1d(maxDim1*SAF_MAX(1, 5*SAF_MIN(maxDim2,maxDim1))*sizeof(double));
    h->work = NULL;
}

void utility_zpinv_destroy(void ** const phWork)
{
    utility_zpinv_data *h = (utility_zpinv_data*)(*phWork);

    if(h!=NULL){
        free(h->a);
        free(h->s);
        free(h->u);
        free(h->vt);
        free(h->inva);
        free(h->work);

        free(h);
        h=NULL;
        *phWork = NULL;
     }
}

void utility_zpinv
(
    void* const hWork,
    const double_complex* inM,
    const int dim1,
    const int dim2,
    double_complex* outM
)
{
    utility_zpinv_data *h;
    veclib_int i, j, m, n, k, lda, ldu, ldvt, info, ld_inva;
    double_complex ss_cmplx;
    const double_complex calpha = cmplx(1.0, 0.0); const double_complex cbeta = cmplx(0.0, 0.0); /* blas */
    double ss;
    veclib_int lwork;
    double_complex wkopt;
    
    m = lda = ldu = dim1;
    n = dim2;
    k = ldvt = m < n ? m : n;

    /* Work struct */
    if(hWork==NULL)
        utility_zpinv_create((void**)&h, dim1, dim2);
    else{
        h = (utility_zpinv_data*)(hWork);
#ifndef NDEBUG
        saf_assert(dim1<=h->maxDim1, "dim1 exceeds the maximum length specified");
        saf_assert(dim2<=h->maxDim2, "dim2 exceeds the maximum length specified");
#endif
    }
    
    /* store in column major order */
    for(i=0; i<dim1; i++)
        for(j=0; j<dim2; j++)
            h->a[j*dim1+i] = inM[i*dim2 +j];

    /* Query how much "work" memory is required */
    lwork = -1;
#if defined(SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE)
    zgesvd_( "A", "A", &m, &n, (veclib_double_complex*)h->a, &lda, h->s, (veclib_double_complex*)h->u, &ldu, (veclib_double_complex*)h->vt, &ldvt,
            (veclib_double_complex*)&wkopt, &lwork, h->rwork, &info );
#elif defined(SAF_VECLIB_USE_CLAPACK_INTERFACE)
    saf_print_error("No such implementation available in ATLAS CLAPACK");
#elif defined(SAF_VECLIB_USE_LAPACKE_INTERFACE)
    info = LAPACKE_zgesvd_work(CblasColMajor, 'A', 'A', m, n, (veclib_double_complex*)h->a, lda, h->s, (veclib_double_complex*)h->u, ldu, (veclib_double_complex*)h->vt, ldvt,
                               (veclib_double_complex*)&wkopt, lwork, h->rwork);
#endif
    lwork = (veclib_int)(creal(wkopt)+0.01);
    if(lwork>h->currentWorkSize){
        h->currentWorkSize = lwork;
        h->work = realloc1d(h->work, h->currentWorkSize*sizeof(double_complex));
    }

    /* Perform the singular value decomposition */
#if defined(SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE)
    zgesvd_( "A", "A", &m, &n, (veclib_double_complex*)h->a, &lda, h->s, (veclib_double_complex*)h->u, &ldu, (veclib_double_complex*)h->vt, &ldvt,
            (veclib_double_complex*)h->work, &lwork, h->rwork, &info);
#elif defined(SAF_VECLIB_USE_CLAPACK_INTERFACE)
    saf_print_error("No such implementation available in ATLAS CLAPACK");
#elif defined(SAF_VECLIB_USE_LAPACKE_INTERFACE)
    info = LAPACKE_zgesvd_work(CblasColMajor, 'A', 'A', m, n, (veclib_double_complex*)h->a, lda, h->s, (veclib_double_complex*)h->u, ldu, (veclib_double_complex*)h->vt, ldvt,
                               (veclib_double_complex*)h->work, lwork, h->rwork);
#endif

    if( info != 0 ) {
        memset(outM, 0, dim1*dim2*sizeof(double_complex));
#ifndef NDEBUG
        /* The SVD failed to converge, or the input matrix contained illegal
         * values so no solution was attempted. In these cases this function
         * will zero all output matrices and/or vectors. */
        saf_print_warning("Could not compute SVD in utility_zpinv(). Output matrices/vectors have been zeroed.");
#endif
    }
    else{
        for(i=0; i<k; i++){
            if(h->s[i] > 1.0e-5)
                ss=1.0/h->s[i];
            else
                ss=h->s[i];
            ss_cmplx = cmplx(ss, 0.0);
            cblas_zscal(m, &ss_cmplx, &h->u[i*m], 1);
        }
        ld_inva=n;
        cblas_zgemm(CblasColMajor, CblasConjTrans, CblasConjTrans, n, m, k, &calpha,
                    h->vt, ldvt,
                    h->u, ldu, &cbeta,
                    h->inva, ld_inva);

        /* return in row-major order */
        for(i=0; i<m; i++)
            for(j=0; j<n; j++)
                outM[j*m+i] = h->inva[i*n+j];
    }
    
    /* clean-up */
    if(hWork == NULL)
        utility_zpinv_destroy((void**)&h);
}


/* ========================================================================== */
/*                       Cholesky Factorisation (?chol)                       */
/* ========================================================================== */

/** Data structure for utility_schol() */
typedef struct _utility_schol_data {
    int maxDim;
    float* a;
}utility_schol_data;

void utility_schol_create(void ** const phWork, int maxDim)
{
    *phWork = malloc1d(sizeof(utility_schol_data));
    utility_schol_data *h = (utility_schol_data*)(*phWork);

    h->maxDim = maxDim;
    h->a = malloc1d(maxDim*maxDim*sizeof(float));
}

void utility_schol_destroy(void ** const phWork)
{
    utility_schol_data *h = (utility_schol_data*)(*phWork);

    if(h!=NULL){
        free(h->a);

        free(h);
        h=NULL;
        *phWork = NULL;
     }
}

void utility_schol
(
    void* const hWork,
    const float* A,
    const int dim,
    float* X
)
{
    utility_schol_data *h;
    veclib_int i, j, info, n, lda;
 
    n = lda = dim;

    /* Work struct */
    if(hWork==NULL)
        utility_schol_create((void**)&h, dim);
    else{
        h = (utility_schol_data*)(hWork);
#ifndef NDEBUG
        saf_assert(dim<=h->maxDim, "dim exceeds the maximum length specified");
#endif
    }
    
    /* store in column major order */
    for(i=0; i<dim; i++)
        for(j=0; j<dim; j++)
            h->a[j*dim+i] = A[i*dim+j];
    
    /* a is replaced by solution */
#ifdef SAF_VECLIB_USE_CLAPACK_INTERFACE
    info = clapack_spotrf(CblasColMajor, CblasUpper, n, h->a, lda);
#elif defined(SAF_VECLIB_USE_LAPACKE_INTERFACE)
    info = LAPACKE_spotrf_work(CblasColMajor, CblasUpper, n, h->a, lda);
#elif defined(SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE)
    spotrf_( "U", &n, h->a, &lda, &info );
#endif
    
    /* A is not positive definate, solution not possible */
    if(info!=0){
        memset(X, 0, dim*dim*sizeof(float));
#ifndef NDEBUG
        /* Input matrix is not positive definite so the Cholesky factorisation
         * could not be computed, or the input matrix contained illegal values
         * so no solution was attempted. In these cases the function will zero
         * the output matrix. */
        saf_print_warning("Could not compute the Cholesky factorisation in utility_schol(). Output matrices/vectors have been zeroed.");
#endif
    }
    /* store solution in row-major order */
    else{
        for(i=0; i<dim; i++)
            for(j=0; j<dim; j++)
                X[i*dim+j] = j>=i ? h->a[j*dim+i] : 0.0f;
    }
    
    /* clean-up */
    if(hWork == NULL)
        utility_schol_destroy((void**)&h);
}

/** Data structure for utility_cchol() */
typedef struct _utility_cchol_data {
    int maxDim;
    float_complex* a;
}utility_cchol_data;

void utility_cchol_create(void ** const phWork, int maxDim)
{
    *phWork = malloc1d(sizeof(utility_cchol_data));
    utility_cchol_data *h = (utility_cchol_data*)(*phWork);

    h->maxDim = maxDim;
    h->a = malloc1d(maxDim*maxDim*sizeof(float_complex));
}

void utility_cchol_destroy(void ** const phWork)
{
    utility_cchol_data *h = (utility_cchol_data*)(*phWork);

    if(h!=NULL){
        free(h->a);

        free(h);
        h=NULL;
        *phWork = NULL;
     }
}

void utility_cchol
(
    void* const hWork,
    const float_complex* A,
    const int dim,
    float_complex* X
)
{
    utility_cchol_data *h;
    veclib_int i, j, info, n, lda;
    
    n = lda = dim;

    /* Work struct */
    if(hWork==NULL)
        utility_cchol_create((void**)&h, dim);
    else{
        h = (utility_cchol_data*)(hWork);
#ifndef NDEBUG
        saf_assert(dim<=h->maxDim, "dim exceeds the maximum length specified");
#endif
    }
    
    /* store in column major order */
    for(i=0; i<dim; i++)
        for(j=0; j<dim; j++)
            h->a[j*dim+i] = A[i*dim+j];
    
    /* a is replaced by solution */
#if defined(SAF_VECLIB_USE_CLAPACK_INTERFACE)
    info = clapack_cpotrf(CblasColMajor, CblasUpper, n, (veclib_float_complex*)h->a, lda);
#elif defined(SAF_VECLIB_USE_LAPACKE_INTERFACE)
    info = LAPACKE_cpotrf_work(CblasColMajor, CblasUpper, n, (veclib_float_complex*)h->a, lda);
#elif defined(SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE)
    cpotrf_( "U", &n, (veclib_float_complex*)h->a, &lda, &info );
#endif
    
    /* A is not positive definate, solution not possible */
    if(info!=0){
        memset(X, 0, dim*dim*sizeof(float_complex));
#ifndef NDEBUG
        /* Input matrix is not positive definite so the Cholesky factorisation
         * could not be computed, or the input matrix contained illegal values
         * so no solution was attempted. In these cases the function will zero
         * the output matrix. */
        saf_print_warning("Could not compute the Cholesky factorisation in utility_cchol(). Output matrices/vectors have been zeroed.");
#endif
    }
    /* store solution in row-major order */
    else{
        for(i=0; i<dim; i++)
            for(j=0; j<dim; j++)
                X[i*dim+j] = j>=i ? h->a[j*dim+i] : cmplxf(0.0f, 0.0f);
    }
    
    /* clean-up */
    if(hWork == NULL)
        utility_cchol_destroy((void**)&h);
}

/* ========================================================================== */
/*                        Determinant of a Matrix (?det)                      */
/* ========================================================================== */

/** Data structure for utility_sdet() */
typedef struct _utility_sdet_data {
    veclib_int maxN;
    veclib_int *IPIV;
    float *tmp;
}utility_sdet_data;

void utility_sdet_create(void ** const phWork, int maxN)
{
    *phWork = malloc1d(sizeof(utility_sdet_data));
    utility_sdet_data *h = (utility_sdet_data*)(*phWork);

    h->maxN = maxN;
    h->IPIV = malloc1d(maxN*sizeof(veclib_int));
    h->tmp = malloc1d(maxN*maxN*sizeof(float));
}

void utility_sdet_destroy(void ** const phWork)
{
    utility_sdet_data *h = (utility_sdet_data*)(*phWork);

    if(h!=NULL){
        free(h->IPIV);
        free(h->tmp);

        free(h);
        h=NULL;
        *phWork = NULL;
     }
}

float utility_sdet
(
    void* const hWork,
    float* A,
    int N
)
{
    utility_sdet_data *h;
    veclib_int i, j, INFO;
    float det;

    /* Simple cases: */
    if(N==2){
        return A[0]*A[3] - A[2]*A[1];
    }
    else if(N==3){
        return A[0] * ((A[4]*A[8]) - (A[7]*A[5])) -A[1] * (A[3]
           * A[8] - A[6] * A[5]) + A[2] * (A[3] * A[7] - A[6] * A[4]);
    }
    else if (N==4){
        return
        A[3] * A[6] * A[9] * A[12] - A[2] * A[7] * A[9] * A[12] -
        A[3] * A[5] * A[10] * A[12] + A[1] * A[7] * A[10] * A[12] +
        A[2] * A[5] * A[11] * A[12] - A[1] * A[6] * A[11] * A[12] -
        A[3] * A[6] * A[8] * A[13] + A[2] * A[7] * A[8] * A[13] +
        A[3] * A[4] * A[10] * A[13] - A[0] * A[7] * A[10] * A[13] -
        A[2] * A[4] * A[11] * A[13] + A[0] * A[6] * A[11] * A[13] +
        A[3] * A[5] * A[8] * A[14] - A[1] * A[7] * A[8] * A[14] -
        A[3] * A[4] * A[9] * A[14] + A[0] * A[7] * A[9] * A[14] +
        A[1] * A[4] * A[11] * A[14] - A[0] * A[5] * A[11] * A[14] -
        A[2] * A[5] * A[8] * A[15] + A[1] * A[6] * A[8] * A[15] +
        A[2] * A[4] * A[9] * A[15] - A[0] * A[6] * A[9] * A[15] -
        A[1] * A[4] * A[10] * A[15] + A[0] * A[5] * A[10] * A[15];
    }

    /* Otherwise, the arbitrary NxN solution: */
    /* Work struct */
    if(hWork==NULL)
        utility_sdet_create((void**)&h, N);
    else{
        h = (utility_sdet_data*)(hWork);
#ifndef NDEBUG
        saf_assert(N<=h->maxN, "N exceeds the maximum length specified");
#endif
    }

    /* Store in column major order */
    for(i=0; i<N; i++)
        for(j=0; j<N; j++)
            h->tmp[j*N+i] = A[i*N+j];

#if defined(SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE)
    sgetrf_((veclib_int*)&N, (veclib_int*)&N, h->tmp, (veclib_int*)&N, h->IPIV, &INFO);
#elif defined(SAF_VECLIB_USE_CLAPACK_INTERFACE)
    INFO = clapack_sgetrf(CblasColMajor, N, N, h->tmp, N, h->IPIV);
#elif defined(SAF_VECLIB_USE_LAPACKE_INTERFACE)
    INFO = LAPACKE_sgetrf_work(CblasColMajor, N, N, h->tmp, N, h->IPIV);
#endif

    if(INFO!=0) {
        det=0.0f;
#ifndef NDEBUG
        /* Input matrix was singular or contained illegal values so no solution
         * was attempted. In these cases the determinant is returned as 0. */
        saf_print_warning("Unable to compute determinant of input matrix. The function utility_sdet() returned 0. ");
#endif
    }
    else {
        det = 1.0f;
        for( i = 0; i < N; i++ ) {
            det*=h->tmp[i*N+i];
            if(h->IPIV[i] != i+1)
                det *= -1.0f;
        }
    }

    /* clean-up */
    if(hWork == NULL)
        utility_sdet_destroy((void**)&h);

    return det;
}

/** Data structure for utility_ddet() */
typedef struct _utility_ddet_data {
    int currentWorkSize;
    veclib_int maxN;
    veclib_int *IPIV;
    double *tmp, *TAU, *WORK;
}utility_ddet_data;

void utility_ddet_create(void ** const phWork, int maxN)
{
    *phWork = malloc1d(sizeof(utility_ddet_data));
    utility_ddet_data *h = (utility_ddet_data*)(*phWork);

    h->maxN = maxN;
    h->currentWorkSize = 0;
    h->IPIV = malloc1d(maxN*sizeof(veclib_int));
    h->tmp = malloc1d(maxN*maxN*sizeof(double));
    h->TAU = malloc1d(maxN*sizeof(double));
    h->WORK = NULL;
}

void utility_ddet_destroy(void ** const phWork)
{
    utility_ddet_data *h = (utility_ddet_data*)(*phWork);

    if(h!=NULL){
        free(h->IPIV);
        free(h->tmp);
        free(h->TAU);
        free(h->WORK);

        free(h);
        h=NULL;
        *phWork = NULL;
     }
}

double utility_ddet
(
    void* const hWork,
    double* A,
    int N
)
{
    utility_ddet_data *h;
    veclib_int i,j,INFO, LWORK, lwork3;
    double lwork2, det;

    /* Simple cases: */
    if(N==2){
        return A[0]*A[3] - A[2]*A[1];
    }
    else if(N==3){
        return A[0] * ((A[4]*A[8]) - (A[7]*A[5])) -A[1] * (A[3]
           * A[8] - A[6] * A[5]) + A[2] * (A[3] * A[7] - A[6] * A[4]);
    }
    else if (N==4){
        return
        A[3] * A[6] * A[9] * A[12] - A[2] * A[7] * A[9] * A[12] -
        A[3] * A[5] * A[10] * A[12] + A[1] * A[7] * A[10] * A[12] +
        A[2] * A[5] * A[11] * A[12] - A[1] * A[6] * A[11] * A[12] -
        A[3] * A[6] * A[8] * A[13] + A[2] * A[7] * A[8] * A[13] +
        A[3] * A[4] * A[10] * A[13] - A[0] * A[7] * A[10] * A[13] -
        A[2] * A[4] * A[11] * A[13] + A[0] * A[6] * A[11] * A[13] +
        A[3] * A[5] * A[8] * A[14] - A[1] * A[7] * A[8] * A[14] -
        A[3] * A[4] * A[9] * A[14] + A[0] * A[7] * A[9] * A[14] +
        A[1] * A[4] * A[11] * A[14] - A[0] * A[5] * A[11] * A[14] -
        A[2] * A[5] * A[8] * A[15] + A[1] * A[6] * A[8] * A[15] +
        A[2] * A[4] * A[9] * A[15] - A[0] * A[6] * A[9] * A[15] -
        A[1] * A[4] * A[10] * A[15] + A[0] * A[5] * A[10] * A[15];
    }

    /* Otherwise, the arbitrary NxN solution: */
    /* Work struct */
    if(hWork==NULL)
        utility_ddet_create((void**)&h, N);
    else{
        h = (utility_ddet_data*)(hWork);
#ifndef NDEBUG
        saf_assert(N<=h->maxN, "N exceeds the maximum length specified");
#endif
    }

    /* Store in column major order */
    for(i=0; i<N; i++)
        for(j=0; j<N; j++)
            h->tmp[j*N+i] = A[i*N+j];

    /* Query how much "work" memory is required */
    LWORK=-1;
#if defined(SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE)
    dgeqrf_(&N, &N, h->tmp, &N, h->TAU, &lwork2, &LWORK, &INFO);
#elif defined(SAF_VECLIB_USE_CLAPACK_INTERFACE)
    saf_print_error("No such implementation in ATLAS CLAPACK");
#elif defined(SAF_VECLIB_USE_LAPACKE_INTERFACE)
    INFO = LAPACKE_dgeqrf_work(CblasColMajor, N, N, h->tmp, N, h->TAU, &lwork2, LWORK);
#endif
    lwork3=(veclib_int)lwork2;
    if(lwork3>h->currentWorkSize){
        h->currentWorkSize = lwork3;
        h->WORK = realloc1d(h->WORK, h->currentWorkSize*sizeof(double));
    }

    /* Decompose matrix */
#if defined(SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE)
    dgeqrf_(&N, &N, h->tmp, &N, h->TAU, h->WORK, &lwork3, &INFO);
#elif defined(SAF_VECLIB_USE_CLAPACK_INTERFACE)
    saf_print_error("No such implementation in ATLAS CLAPACK");
#elif defined(SAF_VECLIB_USE_LAPACKE_INTERFACE)
    INFO = LAPACKE_dgeqrf_work(CblasColMajor, N, N, h->tmp, N, h->TAU, h->WORK, lwork3);
#endif

    if(INFO!=0) {
        det=0.0;
#ifndef NDEBUG
        /* Input matrix was singular or contained illegal values so no solution
         * was attempted. In these cases the determinant is returned as 0. */
        saf_print_warning("Unable to compute determinant of input matrix. The function utility_ddet() returned 0. ");
#endif
    }
    else {
        det=1;
        for (i=0;i<N;i++)
            det*=h->tmp[i*N+i];
        /* Householder algorithm */
        if(N%2==0)
            det*=-1;
    }

    /* clean-up */
    if(hWork == NULL)
        utility_ddet_destroy((void**)&h);

    return det;
}


/* ========================================================================== */
/*                           Matrix Inversion (?inv)                          */
/* ========================================================================== */

/** Data structure for utility_sinv() */
typedef struct _utility_sinv_data {
    int maxN;
    veclib_int *IPIV;
    float *WORK, *tmp;
}utility_sinv_data;

void utility_sinv_create(void ** const phWork, int maxN)
{
    *phWork = malloc1d(sizeof(utility_sinv_data));
    utility_sinv_data *h = (utility_sinv_data*)(*phWork);

    h->maxN = maxN;
    h->IPIV = malloc1d(maxN*sizeof(veclib_int));
    h->tmp = malloc1d(maxN*maxN*sizeof(float));
    h->WORK = malloc1d(maxN*maxN*sizeof(float));
}

void utility_sinv_destroy(void ** const phWork)
{
    utility_sinv_data *h = (utility_sinv_data*)(*phWork);

    if(h!=NULL){
        free(h->IPIV);
        free(h->tmp);
        free(h->WORK);

        free(h);
        h=NULL;
        *phWork = NULL;
     }
}

void utility_sinv
(
    void* const hWork,
    float* A,
    float* B,
    const int N
)
{
    utility_sinv_data *h;
    veclib_int i, j, N_;
    veclib_int LWORK;
    veclib_int INFO;

    LWORK = N*N;
    N_ = (veclib_int)N;

    /* Work struct */
    if(hWork==NULL)
        utility_sinv_create((void**)&h, N);
    else{
        h = (utility_sinv_data*)(hWork);
#ifndef NDEBUG
        saf_assert(N<=h->maxN, "N exceeds the maximum length specified");
#endif
    }

    /* Store in column major order */
    for(i=0; i<N; i++)
        for(j=0; j<N; j++)
            h->tmp[j*N+i] = A[i*N+j];

#if defined(SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE)
    sgetrf_((veclib_int*)&N_, (veclib_int*)&N_, h->tmp, (veclib_int*)&N_, h->IPIV, &INFO);
    sgetri_((veclib_int*)&N_, h->tmp, (veclib_int*)&N_, h->IPIV, h->WORK, &LWORK, &INFO);
#elif defined(SAF_VECLIB_USE_CLAPACK_INTERFACE)
    INFO = clapack_sgetrf(CblasColMajor, N_, N_, h->tmp, N_, h->IPIV);
    INFO = clapack_sgetri(CblasColMajor, N_, h->tmp, N_, h->IPIV);
#elif defined(SAF_VECLIB_USE_LAPACKE_INTERFACE)
    INFO = LAPACKE_sgetrf_work(CblasColMajor, N_, N_, h->tmp, N_, h->IPIV);
    INFO = LAPACKE_sgetri_work(CblasColMajor, N_, h->tmp, N_, h->IPIV, h->WORK, LWORK);
#endif

    if(INFO!=0) {
        memset(B, 0, N*N*sizeof(float));
#ifndef NDEBUG
        /* Input matrix was singular or contained illegal values so no inversion
         * was attempted. In these cases the function will zero the output
         * matrix. */
        saf_print_warning("Unable to compute the inverse of input matrix. The function utility_sinv() returned a matrix of zeros. ");
#endif
    }
    else {
        /* Output in row major order */
        for(i=0; i<N; i++)
            for(j=0; j<N; j++)
                B[j*N+i] = h->tmp[i*N+j];
    }
    
    /* clean-up */
    if(hWork == NULL)
        utility_sinv_destroy((void**)&h);
}

/** Data structure for utility_dinv() */
typedef struct _utility_dinv_data {
    int maxN;
    veclib_int *IPIV;
    double *WORK, *tmp;
}utility_dinv_data;

void utility_dinv_create(void ** const phWork, int maxN)
{
    *phWork = malloc1d(sizeof(utility_dinv_data));
    utility_dinv_data *h = (utility_dinv_data*)(*phWork);

    h->maxN = maxN;
    h->IPIV = malloc1d(maxN*sizeof(veclib_int));
    h->tmp = malloc1d(maxN*maxN*sizeof(double));
    h->WORK = malloc1d(maxN*maxN*sizeof(double));
}

void utility_dinv_destroy(void ** const phWork)
{
    utility_dinv_data *h = (utility_dinv_data*)(*phWork);

    if(h!=NULL){
        free(h->IPIV);
        free(h->tmp);
        free(h->WORK);

        free(h);
        h=NULL;
        *phWork = NULL;
     }
}

void utility_dinv
(
    void* const hWork,
    double* A,
    double* B,
    const int N
)
{
    utility_dinv_data *h;
    veclib_int i, j, N_;
    veclib_int LWORK;
    veclib_int INFO;

    LWORK = N*N;
    N_ = (veclib_int)N;

    /* Work struct */
    if(hWork==NULL)
        utility_dinv_create((void**)&h, N);
    else{
        h = (utility_dinv_data*)(hWork);
#ifndef NDEBUG
        saf_assert(N<=h->maxN, "N exceeds the maximum length specified");
#endif
    }

    /* Store in column major order */
    for(i=0; i<N; i++)
        for(j=0; j<N; j++)
            h->tmp[j*N+i] = A[i*N+j];

#if defined(SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE)
    dgetrf_((veclib_int*)&N_, (veclib_int*)&N_, h->tmp, (veclib_int*)&N_, h->IPIV, &INFO);
    dgetri_((veclib_int*)&N_, h->tmp, (veclib_int*)&N_, h->IPIV, h->WORK, &LWORK, &INFO);
#elif defined(SAF_VECLIB_USE_CLAPACK_INTERFACE)
    INFO = clapack_dgetrf(CblasColMajor, N_, N_, h->tmp, N_, h->IPIV);
    INFO = clapack_dgetri(CblasColMajor, N_, h->tmp, N_, h->IPIV);
#elif defined(SAF_VECLIB_USE_LAPACKE_INTERFACE)
    INFO = LAPACKE_dgetrf_work(CblasColMajor, N_, N_, h->tmp, N_, h->IPIV);
    INFO = LAPACKE_dgetri_work(CblasColMajor, N_, h->tmp, N_, h->IPIV, h->WORK, LWORK);
#endif

    if(INFO!=0) {
        memset(B, 0, N*N*sizeof(double));
#ifndef NDEBUG
        /* Input matrix was singular or contained illegal values so no inversion
         * was attempted. In these cases the function will zero the output
         * matrix. */
        saf_print_warning("Unable to compute the inverse of input matrix. The function utility_dinv() returned a matrix of zeros. ");
#endif
    }
    else {
        /* Output in row major order */
        for(i=0; i<N; i++)
            for(j=0; j<N; j++)
                B[j*N+i] = h->tmp[i*N+j];
    }
    
    /* clean-up */
    if(hWork == NULL)
        utility_dinv_destroy((void**)&h);
}

/** Data structure for utility_cinv() */
typedef struct _utility_cinv_data {
    int maxN;
    veclib_int *IPIV;
    float_complex *WORK, *tmp;
}utility_cinv_data;

void utility_cinv_create(void ** const phWork, int maxN)
{
    *phWork = malloc1d(sizeof(utility_cinv_data));
    utility_cinv_data *h = (utility_cinv_data*)(*phWork);

    h->maxN = maxN;
    h->IPIV = malloc1d(maxN*sizeof(veclib_int));
    h->tmp = malloc1d(maxN*maxN*sizeof(float_complex));
    h->WORK = malloc1d(maxN*maxN*sizeof(float_complex));
}

void utility_cinv_destroy(void ** const phWork)
{
    utility_cinv_data *h = (utility_cinv_data*)(*phWork);

    if(h!=NULL){
        free(h->IPIV);
        free(h->tmp);
        free(h->WORK);

        free(h);
        h=NULL;
        *phWork = NULL;
     }
}

void utility_cinv
(
    void* const hWork,
    float_complex* A,
    float_complex* B,
    const int N
)
{
    utility_cinv_data *h;
    veclib_int i, j, N_;
    veclib_int LWORK;
    veclib_int INFO;

    N_ = (veclib_int)N;
    LWORK = N*N;

    /* Work struct */
    if(hWork==NULL)
        utility_cinv_create((void**)&h, N);
    else{
        h = (utility_cinv_data*)(hWork);
#ifndef NDEBUG
        saf_assert(N<=h->maxN, "N exceeds the maximum length specified");
#endif
    }

    /* Store in column major order */
    for(i=0; i<N; i++)
        for(j=0; j<N; j++)
            h->tmp[j*N+i] = A[i*N+j];
    
#if defined(SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE)
    cgetrf_((veclib_int*)&N_, (veclib_int*)&N_, (veclib_float_complex*)h->tmp, (veclib_int*)&N_, h->IPIV, &INFO);
    cgetri_((veclib_int*)&N_, (veclib_float_complex*)h->tmp, (veclib_int*)&N_, h->IPIV, (veclib_float_complex*)h->WORK, &LWORK, &INFO);
#elif defined(SAF_VECLIB_USE_CLAPACK_INTERFACE)
    INFO = clapack_cgetrf(CblasColMajor, N_, N_, (veclib_float_complex*)h->tmp, N_, h->IPIV);
    INFO = clapack_cgetri(CblasColMajor, N_, (veclib_float_complex*)h->tmp, N_, h->IPIV);
#elif defined(SAF_VECLIB_USE_LAPACKE_INTERFACE)
    INFO = LAPACKE_cgetrf_work(CblasColMajor, N_, N_, (veclib_float_complex*)h->tmp, N_, h->IPIV);
    INFO = LAPACKE_cgetri_work(CblasColMajor, N_, (veclib_float_complex*)h->tmp, N_, h->IPIV, (veclib_float_complex*)h->WORK, LWORK);
#endif

    if(INFO!=0) {
        memset(B, 0, N*N*sizeof(float_complex));
#ifndef NDEBUG
        /* Input matrix was singular or contained illegal values so no inversion
         * was attempted. In these cases the function will zero the output
         * matrix. */
        saf_print_warning("Unable to compute the inverse of input matrix. The function utility_cinv() returned a matrix of zeros. ");
#endif
    }
    else {
        /* Output in row major order */
        for(i=0; i<N; i++)
            for(j=0; j<N; j++)
                B[j*N+i] = h->tmp[i*N+j];
    }

    /* clean-up */
    if(hWork == NULL)
        utility_cinv_destroy((void**)&h);
}
