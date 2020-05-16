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
 * @file saf_veclib.h
 * @brief Contains wrappers for optimised linear algebra routines, utilising
 *        CBLAS and LAPACK
 *
 * ## Dependencies
 *   A performance library comprising CBLAS and LAPACK routines is required by
 *   the module and, thus, also by the SAF framework as a whole. Add one of the
 *   following FLAGS to your project's preprocessor definitions list, in order
 *   to enable one of these suitable performance libraries, which must also be
 *   linked correctly to your project.
 *   - SAF_USE_INTEL_MKL:
 *       to enable Intel's Math Kernal Library with Fortran LAPACK interface
 *   - SAF_USE_ATLAS:
 *       to enable ATLAS BLAS routines and ATLAS's CLAPACK interface
 *   - SAF_USE_OPENBLAS_WITH_LAPACKE:
 *       to enable OpenBLAS with LAPACKE interface
 *
 * @see More information can be found here:
 *      https://github.com/leomccormack/Spatial_Audio_Framework
 * @note MacOSX users only: saf_utilities will employ Apple's Accelerate library
 *       by default, if none of the above FLAGS are defined.
 *
 * @author Leo McCormack
 * @date 11.07.2016
 */

#ifndef SAF_VECLIB_H_INCLUDED
#define SAF_VECLIB_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "saf_complex.h"
#include "saf_error.h"
#include <immintrin.h>
#ifdef __SSE__
# include <xmmintrin.h>
#endif
#ifdef __SSE2__
# include <emmintrin.h>
#endif
#ifdef __SSE3__
# include <pmmintrin.h>
#endif
#ifdef __SSSE3__
# include <tmmintrin.h>
#endif
#ifdef __SSE4_1__
# include <smmintrin.h>
#endif
#ifdef CBLAS_H
# define NO_TRANSPOSE (CblasNoTrans)
# define TRANSPOSE (CblasTrans)
# define CONJ_TRANSPOSE (CblasConjTrans)
  typedef enum CBLAS_TRANSPOSE TRANS_FLAG;
#else
  typedef enum _TRANS_FLAG{
    NO_TRANSPOSE = 1,   /**< Do not transpose */
    TRANSPOSE = 2,      /**< Transpose */
    CONJ_TRANSPOSE = 3  /**< Conjugate transpose / Hermition */
  }TRANS_FLAG;
#endif
typedef enum _CONJ_FLAG{
  NO_CONJ = 1,  /**< Do not take the conjugate */
  CONJ = 2      /**< Take the conjugate */
}CONJ_FLAG;

/* KEY:
 * ? -> s -> single floating-point precision
 * ? -> d -> double floating-point precision
 * ? -> c -> complex single floating-point precision
 * ? -> z -> complex double floating-point precision
 *
 * s -> scalar
 * v -> vector
 * m -> matrix */

/* ========================================================================== */
/*                     Find Index of Min-Abs-Value (?iminv)                   */
/* ========================================================================== */

/**
 * Single-precision, index of minimum absolute value in a vector, i.e.
 * \code{.m}
 *     [~,ind] = min(abs(a))
 * \endcode
 *
 * @param[in]  a     Input vector a; len x 1
 * @param[in]  len   Vector length
 * @param[out] index (&) Index of minimum value; 1 x 1
 */
void utility_siminv(/* Input Arguments */
                    const float* a,
                    const int len,
                    /* Output Arguments */
                    int* index);

/**
 * Single-precision, complex, index of maximum absolute value in a vector, i.e.
 * \code{.m}
 *     [~,ind] = min(abs(a))
 * \endcode
 *
 * @param[in]  a     Input vector a; len x 1
 * @param[in]  len   Vector length
 * @param[out] index (&) index of minimum value; 1 x 1
 */
void utility_ciminv(/* Input Arguments */
                    const float_complex* a,
                    const int len,
                    /* Output Arguments */
                    int* index);

    
/* ========================================================================== */
/*                     Find Index of Max-Abs-Value (?imaxv)                   */
/* ========================================================================== */

/**
 * Single-precision, index of maximum absolute value in a vector, i.e.
 * \code{.m}
 *     [~,ind] = max(abs(a))
 * \endcode
 *
 * @param[in]  a     Input vector a; len x 1
 * @param[in]  len   Vector length
 * @param[out] index (&) index of maximum value; 1 x 1
 */
void utility_simaxv(/* Input Arguments */
                    const float* a,
                    const int len,
                    /* Output Arguments */
                    int* index);

/**
 * Single-precision, complex, index of maximum absolute value in a vector, i.e.
 * \code{.m}
 *     [~,ind] = max(abs(a))
 * \endcode
 *
 * @param[in]  a     Input vector a; len x 1
 * @param[in]  len   Vector length
 * @param[out] index (&) index of maximum value; 1 x 1
 */
void utility_cimaxv(/* Input Arguments */
                    const float_complex* a,
                    const int len,
                    /* Output Arguments */
                    int* index);


/* ========================================================================== */
/*                              Vector-Abs (?vabs)                            */
/* ========================================================================== */

/**
 * Single-precision, absolute values of vector elements, i.e.
 * \code{.m}
 *     c = abs(a)
 * \endcode
 *
 * @param[in]  a   Input vector a; len x 1
 * @param[in]  len Vector length
 * @param[out] c   Output vector c; len x 1
 */
void utility_svabs(/* Input Arguments */
                   const float* a,
                   const int len,
                   /* Output Arguments */
                   float* c);

/**
 * Single-precision, complex, absolute values of vector elements, i.e.
 * \code{.m}
 *     c = abs(a)
 * \endcode
 *
 * @param[in]  a   Input vector a; len x 1
 * @param[in]  len Vector length
 * @param[out] c   Output vector c; len x 1
 */
void utility_cvabs(/* Input Arguments */
                   const float_complex* a,
                   const int len,
                   /* Output Arguments */
                   float* c);


/* ========================================================================== */
/*                        Vector-Vector Copy (?vvcopy)                        */
/* ========================================================================== */

/**
 * Single-precision, vector-vector copy, i.e.
 * \code{.m}
 *     c = a
 * \endcode
 *
 * @param[in]  a   Input vector a; len x 1
 * @param[in]  len Vector length
 * @param[out] c   Output vector c; len x 1
 */
void utility_svvcopy(/* Input Arguments */
                     const float* a,
                     const int len,
                     /* Output Arguments */
                     float* c);

/**
 * Single-precision, complex, vector-vector copy, i.e.
 * \code{.m}
 *     c = a
 * \endcode
 *
 * @param[in]  a   Input vector a; len x 1
 * @param[in]  len Vector length
 * @param[out] c   Output vector c; len x 1
 */
void utility_cvvcopy(/* Input Arguments */
                     const float_complex* a,
                     const int len,
                     /* Output Arguments */
                     float_complex* c);


/* ========================================================================== */
/*                       Vector-Vector Addition (?vvadd)                      */
/* ========================================================================== */

/**
 * Single-precision, vector-vector addition, i.e.
 * \code{.m}
 *     c = a+b, OR: a = a+b (if c==NULL)
 * \endcode
 *
 * @param[in]  a   Input vector a, and output if c==NULL; len x 1
 * @param[in]  b   Input vector b; len x 1
 * @param[in]  len Vector length
 * @param[out] c   Output vector c (set to NULL if you want 'a' as output);
 *                 len x 1
 */
void utility_svvadd(/* Input Arguments */
                    float* a,
                    const float* b,
                    const int len,
                    /* Output Arguments */
                    float* c);

/**
 * Single-precision, complex, vector-vector addition, i.e.
 * \code{.m}
 *     c = a+b, OR: a = a+b (if c==NULL)
 * \endcode
 *
 * @param[in]  a   Input vector a, and output if c==NULL; len x 1
 * @param[in]  b   Input vector b; len x 1
 * @param[in]  len Vector length
 * @param[out] c   Output vector c (set to NULL if you want 'a' as output);
 *                 len x 1
 */
void utility_cvvadd(/* Input Arguments */
                    float_complex* a,
                    const float_complex* b,
                    const int len,
                    /* Output Arguments */
                    float_complex* c);


/* ========================================================================== */
/*                     Vector-Vector Subtraction (?vvsub)                     */
/* ========================================================================== */

/**
 * Single-precision, vector-vector subtraction, i.e.
 * \code{.m}
 *     c = a-b, OR: a = a-b (if c==NULL)
 * \endcode
 *
 * @param[in]  a   Input vector a, and output if c==NULL; len x 1
 * @param[in]  b   Input vector b; len x 1
 * @param[in]  len Vector length
 * @param[out] c   Output vector c (set to NULL if you want 'a' as output);
 *                 len x 1
 */
void utility_svvsub(/* Input Arguments */
                    float* a,
                    const float* b,
                    const int len,
                    /* Output Arguments */
	                float* c);

/**
 * Single-precision, complex, vector-vector subtraction, i.e.
 * \code{.m}
 *     c = a-b, OR: a = a-b (if c==NULL)
 * \endcode
 *
 * @param[in]  a   Input vector a, and output if c==NULL; len x 1
 * @param[in]  b   Input vector b; len x 1
 * @param[in]  len Vector length
 * @param[out] c   Output vector c (set to NULL if you want 'a' as output);
 *                 len x 1
 */
void utility_cvvsub(/* Input Arguments */
                    float_complex* a,
                    const float_complex* b,
                    const int len,
                    /* Output Arguments */
                    float_complex* c);


/* ========================================================================== */
/*                    Vector-Vector Multiplication (?vvmul)                   */
/* ========================================================================== */

/**
 * Single-precision, element-wise vector-vector multiplication i.e.
 * \code{.m}
 *     c = a.*b
 * \endcode
 *
 * @param[in]  a   Input vector a; len x 1
 * @param[in]  b   Input vector b; len x 1
 * @param[in]  len Vector length
 * @param[out] c   Output vector c
 */
void utility_svvmul(/* Input Arguments */
                    float* a,
                    const float* b,
                    const int len,
                    /* Output Arguments */
                    float* c);

/**
 * Single-precision, complex, element-wise vector-vector multiplication i.e.
 * \code{.m}
 *     c = a.*b
 * \endcode
 *
 * @param[in]  a   Input vector a; len x 1
 * @param[in]  b   Input vector b; len x 1
 * @param[in]  len Vector length
 * @param[out] c   Output vector c; len x 1
 */
void utility_cvvmul(/* Input Arguments */
                    float_complex* a,
	                const float_complex* b,
	                const int len,
                    /* Output Arguments */
	                float_complex* c);


/* ========================================================================== */
/*                     Vector-Vector Dot Product (?vvdot)                     */
/* ========================================================================== */

/**
 * Single-precision, vector-vector dot product, i.e.
 * \code{.m}
 *     c = a*b^T, (where size(c) = [1  1])
 * \endcode
 *
 * @param[in]  a   Input vector a; len x 1
 * @param[in]  b   Input vector b; len x 1
 * @param[in]  len Vector length
 * @param[out] c   (&) output vector c; 1 x 1
 */
void utility_svvdot(/* Input Arguments */
                    const float* a,
                    const float* b,
                    const int len,
                    /* Output Arguments */
                    float* c);

/**
 * Single-precision, complex, vector-vector dot product, i.e.
 * \code{.m}
 *     c = a*b^T, (where size(c) = [1  1])
 * \endcode
 *
 * @param[in]  a    Input vector a; len x 1
 * @param[in]  b    Input vector b; len x 1
 * @param[in]  flag '0' do not take the conjugate of 'b', '1', take the
 *                  conjugate of 'b'.
 * @param[in]  len  Vector length
 * @param[out] c    (&) output vector c; 1 x 1
 */
void utility_cvvdot(/* Input Arguments */
                    const float_complex* a,
                    const float_complex* b,
                    const int len,
                    CONJ_FLAG flag,
                    /* Output Arguments */
                    float_complex* c);


/* ========================================================================== */
/*                       Vector-Scalar Product (?vsmul)                       */
/* ========================================================================== */

/**
 * Single-precision, multiplies each element in vector 'a' with a scalar 's',
 * i.e.
 * \code{.m}
 *     c = a.*s, OR: a = a.*s (if c==NULL)
 * \endcode
 *
 * @param[in]  a   Input vector a, and output if c==NULL; len x 1
 * @param[in]  s   (&) input scalar s; 1 x 1
 * @param[in]  len Vector length
 * @param[out] c   Output vector c (set to NULL if you want 'a' as output);
 *                 len x 1
 */
void utility_svsmul(/* Input Arguments */
                    float* a,
                    const float* s,
                    const int len,
                    /* Output Arguments */
                    float* c);

/**
 * Single-precision, complex, multiplies each element in vector 'a' with a
 * scalar 's', i.e.
 * \code{.m}
 *     c = a.*s, OR: a = a.*s (if c==NULL)
 * \endcode
 *
 * @param[in]  a   Input vector a, and output if c==NULL; len x 1
 * @param[in]  s   (&) input scalar s; 1 x 1
 * @param[in]  len Vector length
 * @param[out] c   Output vector c (set to NULL if you want 'a' as output);
 */
void utility_cvsmul(/* Input Arguments */
                    float_complex* a,
                    const float_complex* s,
                    const int len,
                    /* Output Arguments */
                    float_complex* c);


/* ========================================================================== */
/*                       Vector-Scalar Division (?vsdiv)                      */
/* ========================================================================== */

/**
 * Single-precision, divides each element in vector 'a' with a scalar 's',
 * i.e.
 * \code{.m}
 *     c = a./s, OR: a = a./s (if c==NULL)
 * \endcode
 *
 * @param[in]  a   Input vector a, and output if c==NULL; len x 1
 * @param[in]  s   (&) input scalar s; 1 x 1
 * @param[in]  len Vector length
 * @param[out] c   Output vector c (set to NULL if you want 'a' as output);
 *                 len x 1
 */
void utility_svsdiv(/* Input Arguments */
                    float* a,
                    const float* s,
                    const int len,
                    /* Output Arguments */
                    float* c);


/* ========================================================================== */
/*                       Vector-Scalar Addition (?vsadd)                      */
/* ========================================================================== */

/**
 * Single-precision, adds each element in vector 'a' with a scalar 's', i.e.
 * \code{.m}
 *     c = a+s, OR: a = a+s (if c==NULL)
 * \endcode
 *
 * @param[in]  a   Input vector a, and output if c==NULL; len x 1
 * @param[in]  s   (&) input scalar s; 1 x 1
 * @param[in]  len Vector length
 * @param[out] c   Output vector c (set to NULL if you want 'a' as output);
 *                 len x 1
 */
void utility_svsadd(/* Input Arguments */
                    float* a,
                    const float* s,
                    const int len,
                    /* Output Arguments */
                    float* c);


/* ========================================================================== */
/*                     Vector-Scalar Subtraction (?vssub)                     */
/* ========================================================================== */
    
/**
 * Single-precision, subtracts each element in vector 'a' with a scalar 's',
 * i.e.
 * \code{.m}
 *     c = a-s, OR: a = a-s (if c==NULL)
 * \endcode
 *
 * @param[in]  a   Input vector a, and output if c==NULL; len x 1
 * @param[in]  s   (&) input scalar s; 1 x 1
 * @param[in]  len Vector length
 * @param[out] c   Output vector c (set to NULL if you want 'a' as output);
 *                 len x 1
 */
void utility_svssub(/* Input Arguments */
                    float* a,
                    const float* s,
                    const int len,
                    /* Output Arguments */
                    float* c);


/* ========================================================================== */
/*                     Singular-Value Decomposition (?svd)                    */
/* ========================================================================== */
    
/**
 * Row-major, singular value decomposition: single precision, i.e.
 * \code{.m}
 *     [U,S,V] = svd(A); such that A = U*S*V.' = U*diag(sing)*V.'
 * \endcode
 *
 * @note 'S' contains the singular values along the diagonal, whereas 'sing' are
 *       the singular values as a vector. Also, V is returned untransposed!
 *       (like in Matlab)
 *
 * @param[in]  A    Input matrix; FLAT: dim1 x dim2
 * @param[in]  dim1 First dimension of matrix 'A'
 * @param[in]  dim2 Second dimension of matrix 'A'
 * @param[out] U    Left matrix (set to NULL if not needed); FLAT: dim1 x dim1
 * @param[out] S    Singular values along the diagonal min(dim1, dim2), (set to
 *                  NULL if not needed); FLAT: dim1 x dim2
 * @param[out] V    Right matrix (UNTRANSPOSED!) (set to NULL if not needed);
 *                  FLAT: dim2 x dim2
 * @param[out] sing Singular values as a vector, (set to NULL if not needed);
 *                  min(dim1, dim2) x 1
 */
void utility_ssvd(/* Input Arguments */
                  const float* A,
                  const int dim1,
                  const int dim2,
                  /* Output Arguments */
                  float* U,
                  float* S,
                  float* V,
                  float* sing);

/**
 * Row-major, singular value decomposition: single precision complex, i.e.
 * \code{.m}
 *     [U,S,V] = svd(A); such that A = U*S*V' = U*diag(sing)*V'
 * \endcode
 *
 * @note 'S' contains the singular values along the diagonal, whereas 'sing' are
 *       the singular values as a vector. Also, V is returned untransposed!
 *       (like in Matlab)
 *
 * @param[in]  A    Input matrix; FLAT: dim1 x dim2
 * @param[in]  dim1 First dimension of matrix 'A'
 * @param[in]  dim2 Second dimension of matrix 'A'
 * @param[out] U    Left matrix (set to NULL if not needed); FLAT: dim1 x dim1
 * @param[out] S    Singular values along the diagonal min(dim1, dim2), (set to
 *                  NULL if not needed); FLAT: dim1 x dim2
 * @param[out] V    Right matrix (UNTRANSPOSED!) (set to NULL if not needed);
 *                  FLAT: dim2 x dim2
 * @param[out] sing Singular values as a vector, (set to NULL if not needed);
 *                  min(dim1, dim2) x 1
 */
void utility_csvd(/* Input Arguments */
                  const float_complex* A,
                  const int dim1,
                  const int dim2,
                  /* Output Arguments */
                  float_complex* U,
                  float_complex* S,
                  float_complex* V,
                  float* sing);


/* ========================================================================== */
/*                 Symmetric Eigenvalue Decomposition (?seig)                 */
/* ========================================================================== */

/**
 * Row-major, eigenvalue decomposition of a SYMMETRIC matrix: single precision,
 * i.e.
 * \code{.m}
 *     [V,D] = eig(A); where A*V = V*D, and  A*V = V*diag(eig)
 * \endcode
 *
 * @note 'D' contains the eigen values along the diagonal, while 'eig' are the
 *       eigen values as a vector
 *
 * @param[in]  A           Input SYMMETRIC square matrix; FLAT: dim x dim
 * @param[in]  dim         Dimensions for square matrix 'A'
 * @param[in]  sortDecFLAG '1' sort eigen values and vectors in decending order.
 *                         '0' ascending
 * @param[out] V           Eigen vectors (set to NULL if not needed);
 *                         FLAT: dim x dim
 * @param[out] D           Eigen values along the diagonal (set to NULL if not
 *                         needed); FLAT: dim x dim
 * @param[out] eig         Eigen values not diagonalised (set to NULL if not
 *                         needed); dim x 1
 */
void utility_sseig(/* Input Arguments */
                   const float* A,
                   const int dim,
                   int sortDecFLAG,
                   /* Output Arguments */
                   float* V,
                   float* D,
                   float* eig);

/**
 * Row-major, eigenvalue decomposition of a SYMMETRIC/HERMITION matrix: single
 * precision complex, i.e.
 * \code{.m}
 *     [V,D] = eig(A); where A*V = V*D, and  A*V = V*diag(eig)
 * \endcode
 *
 * @note 'D' contains the eigen values along the diagonal, while 'eig' are the
 *       eigen values as a vector
 *
 * @param[in]  A           Input SYMMETRIC square matrix; FLAT: dim x dim
 * @param[in]  dim         Dimensions for square matrix 'A'
 * @param[in]  sortDecFLAG '1' sort eigen values and vectors in decending order.
 *                         '0' ascending
 * @param[out] V           Eigen vectors (set to NULL if not needed);
 *                         FLAT: dim x dim
 * @param[out] D           Eigen values along the diagonal (set to NULL if not
 *                         needed); FLAT: dim x dim
 * @param[out] eig         Eigen values not diagonalised (set to NULL if not
 *                         needed); dim x 1
 */
void utility_cseig(/* Input Arguments */
                   const float_complex* A,
                   const int dim,
                   int sortDecFLAG,
                   /* Output Arguments */
                   float_complex* V,
                   float_complex* D,
                   float* eig);


/* ========================================================================== */
/*                     Eigenvalues of Matrix Pair (?eigmp)                    */
/* ========================================================================== */

/**
 * Row-major, finds eigenvalues of a matrix pair using the QZ method, single
 * precision complex, i.e.
 * \code{.m}
 *     [VL,VR,D] = eig(A,B,'qz'); where A*VL = B*VL*VR
 * \endcode
 *
 * @param[in]  A   Input left square matrix; FLAT: dim x dim
 * @param[in]  B   Input right square matrix; FLAT: dim x dim
 * @param[in]  dim Dimensions for square matrices 'A' and 'B'
 * @param[out] VL  Left Eigen vectors (set to NULL if not needed);
 *                 FLAT: dim x dim
 * @param[out] VR  Right Eigen vectors (set to NULL if not needed);
 *                 FLAT: dim x dim
 * @param[out] D   Eigen values along the diagonal (set to NULL if not needed);
 *                 FLAT: dim x dim
 */
void utility_ceigmp(/* Input Arguments */
                    const float_complex* A,
                    const float_complex* B,
                    const int dim,
                    /* Output Arguments */
                    float_complex* VL,
                    float_complex* VR,
                    float_complex* D);

/**
 * Row-major, finds eigenvalues of a matrix pair using the QZ method, double
 * precision complex, i.e.
 * \code{.m}
 *     [VL,VR,D] = eig(A,B,'qz'); where A*VL = B*VL*VR
 * \endcode
 *
 * @param[in]  A   Input left square matrix; FLAT: dim x dim
 * @param[in]  B   Input right square matrix; FLAT: dim x dim
 * @param[in]  dim Dimensions for square matrices 'A' and 'B'
 * @param[out] VL  Left Eigen vectors (set to NULL if not needed);
 *                 FLAT: dim x dim
 * @param[out] VR  Right Eigen vectors (set to NULL if not needed);
 *                 FLAT: dim x dim
 * @param[out] D   Eigen values along the diagonal (set to NULL if not needed);
 *                 FLAT: dim x dim
 */
void utility_zeigmp(/* Input Arguments */
                    const double_complex* A,
                    const double_complex* B,
                    const int dim,
                    /* Output Arguments */
                    double_complex* VL,
                    double_complex* VR,
                    double_complex* D);


/* ========================================================================== */
/*                       Eigenvalue Decomposition (?eig)                      */
/* ========================================================================== */

/**
 * Row-major, eigenvalue decomposition of a NON-SYMMETRIC matrix: single
 * precision complex, i.e.
 * \code{.m}
 *     [VL,VR,D] = eig(A); where A*VR = VR*D, and  A*VR = VR*diag(eig)
 * \endcode
 *
 * @note 'D' contains the eigen values along the diagonal, while 'eig' are the
 *       eigen values as a vector
 *
 * @param[in]  A           Input NON-SYMMETRIC square matrix; FLAT: dim x dim
 * @param[in]  dim         Dimensions for square matrix 'A'
 * @param[out] VL          Left Eigen vectors (set to NULL if not needed);
 *                         FLAT: dim x dim
 * @param[out] VR          Right Eigen vectors (set to NULL if not needed);
 *                         FLAT: dim x dim
 * @param[out] D           Eigen values along the diagonal (set to NULL if not
 *                         needed); FLAT: dim x dim
 * @param[out] eig         Eigen values not diagonalised (set to NULL if not
 *                         needed); dim x 1
 */
void utility_ceig(/* Input Arguments */
                  const float_complex* A,
                  const int dim, 
                  /* Output Arguments */
                  float_complex* VL,
                  float_complex* VR,
                  float_complex* D,
                  float_complex* eig);

/**
 * Row-major, eigenvalue decomposition of a NON-SYMMETRIC matrix: double
 * precision complex, i.e.
 * \code{.m}
 *     [VL,VR,D] = eig(A); where A*VR = VR*D, and  A*VR = VR*diag(eig)
 * \endcode
 *
 * @note 'D' contains the eigen values along the diagonal, while 'eig' are the
 *       eigen values as a vector
 *
 * @param[in]  A           Input NON-SYMMETRIC square matrix; FLAT: dim x dim
 * @param[in]  dim         Dimensions for square matrix 'A'
 * @param[out] VL          Left Eigen vectors (set to NULL if not needed);
 *                         FLAT: dim x dim
 * @param[out] VR          Right Eigen vectors (set to NULL if not needed);
 *                         FLAT: dim x dim
 * @param[out] D           Eigen values along the diagonal (set to NULL if not
 *                         needed); FLAT: dim x dim
 * @param[out] eig         Eigen values not diagonalised (set to NULL if not
 *                         needed); dim x 1
 */
void utility_zeig(/* Input Arguments */
                  const double_complex* A,
                  const int dim,
                  /* Output Arguments */
                  double_complex* VL,
                  double_complex* VR,
                  double_complex* D,
                  double_complex* eig);


/* ========================================================================== */
/*                       General Linear Solver (?glslv)                       */
/* ========================================================================== */

/**
 * Row-major, general linear solver: single precision, i.e.
 * \code{.m}
 *     X = linsolve(A,B) = A\B; where, AX = B
 * \endcode
 *
 * @param[in]  A    Input square matrix; FLAT: dim x dim
 * @param[in]  dim  Dimensions for square matrix 'A'
 * @param[in]  B    Right hand side matrix; FLAT: dim x nCol
 * @param[in]  nCol Number of columns in right hand side matrix
 * @param[out] X    The solution; FLAT: dim x nCol
 */
void utility_sglslv(/* Input Arguments */
                    const float* A,
                    const int dim,
                    float* B,
                    int nCol,
                    /* Output Arguments */
                    float* X);

/**
 * Row-major, general linear solver: single precision complex, i.e.
 * \code{.m}
 *     X = linsolve(A,B) = A\B; where, AX = B
 * \endcode
 *
 * @param[in]  A    Input square matrix; FLAT: dim x dim
 * @param[in]  dim  Dimensions for square matrix 'A'
 * @param[in]  B    Right hand side matrix; FLAT: dim x nCol
 * @param[in]  nCol Number of columns in right hand side matrix
 * @param[out] X    The solution; FLAT: dim x nCol
 */
void utility_cglslv(/* Input Arguments */
                    const float_complex* A,
                    const int dim,
                    float_complex* B,
                    int nCol,
                    /* Output Arguments */
                    float_complex* X);

/**
 * Row-major, general linear solver: double precision, i.e.
 * \code{.m}
 *     X = linsolve(A,B) = A\B; where, AX = B
 * \endcode
 *
 * @param[in]  A    Input square matrix; FLAT: dim x dim
 * @param[in]  dim  Dimensions for square matrix 'A'
 * @param[in]  B    Right hand side matrix; FLAT: dim x nCol
 * @param[in]  nCol Number of columns in right hand side matrix
 * @param[out] X    The solution; FLAT: dim x nCol
 */
void utility_dglslv(/* Input Arguments */
                    const double* A,
                    const int dim,
                    double* B,
                    int nCol,
                    /* Output Arguments */
                    double* X);

/**
 * Row-major, general linear solver: double precision complex, i.e.
 * \code{.m}
 *     X = linsolve(A,B) = A\B; where, AX = B
 * \endcode
 *
 * @param[in]  A    Input square matrix; FLAT: dim x dim
 * @param[in]  dim  Dimensions for square matrix 'A'
 * @param[in]  B    Right hand side matrix; FLAT: dim x nCol
 * @param[in]  nCol Number of columns in right hand side matrix
 * @param[out] X    The solution; FLAT: dim x nCol
 */
void utility_zglslv(/* Input Arguments */
                    const double_complex* A,
                    const int dim,
                    double_complex* B,
                    int nCol,
                    /* Output Arguments */
                    double_complex* X);


/* ========================================================================== */
/*                      Symmetric Linear Solver (?slslv)                      */
/* ========================================================================== */

/**
 * Row-major, linear solver for SYMMETRIC positive-definate 'A': single
 * precision, i.e.
 * \code{.m}
 *     opts.LT=true
 *     X = linsolve(A,B, opts); where, AX = B
 * \endcode
 *
 * @param[in]  A    Input square SYMMETRIC positive-definate matrix; FLAT: dim x dim
 * @param[in]  dim  Dimensions for square matrix 'A'
 * @param[in]  B    Right hand side matrix; FLAT: dim x nCol
 * @param[in]  nCol Number of columns in right hand side matrix
 * @param[out] X    The solution; FLAT: dim x nCol
 */
void utility_sslslv(/* Input Arguments */
                    const float* A,
                    const int dim,
                    float* B,
                    int nCol,
                    /* Output Arguments */
                    float* X);

/**
 * Row-major, linear solver for HERMITIAN positive-definate 'A': single
 * precision complex, i.e.
 * \code{.m}
 *     opts.LT=true
 *     X = linsolve(A,B, opts); where, AX = B
 * \endcode
 *
 * @param[in]  A    Input square SYMMETRIC positive-definate matrix; FLAT: dim x dim
 * @param[in]  dim  Dimensions for square matrix 'A'
 * @param[in]  B    Right hand side matrix; FLAT: dim x nCol
 * @param[in]  nCol Number of columns in right hand side matrix
 * @param[out] X    The solution; FLAT: dim x nCol
 */
void utility_cslslv(/* Input Arguments */
                    const float_complex* A,
                    const int dim,
                    float_complex* B,
                    int nCol,
                    /* Output Arguments */
                    float_complex* X);


/* ========================================================================== */
/*                        Matrix Pseudo-Inverse (?pinv)                       */
/* ========================================================================== */

/**
 * Row-major, general matrix pseudo-inverse (the svd way): single precision,
 * i.e.
 * \code{.m}
 *     B = pinv(A)
 * \endcode
 *
 * @param[in]  A    Input matrix; FLAT: dim1 x dim2
 * @param[in]  dim1 Number of rows in 'A' / columns in 'B'
 * @param[in]  dim2 Number of columns in 'A' / rows in 'B'
 * @param[out] B    The solution; FLAT: dim2 x dim1
 */
void utility_spinv(/* Input Arguments */
                   const float* A,
                   const int dim1,
                   const int dim2,
                   /* Output Arguments */
                   float* B);

/**
 * Row-major, general matrix pseudo-inverse (the svd way): single precision
 * complex, i.e.
 * \code{.m}
 *     B = pinv(A)
 * \endcode
 *
 * @param[in]  A    Input matrix; FLAT: dim1 x dim2
 * @param[in]  dim1 Number of rows in 'A' / columns in 'B'
 * @param[in]  dim2 Number of columns in 'A' / rows in 'B'
 * @param[out] B    The solution; FLAT: dim2 x dim1
 */
void utility_cpinv(/* Input Arguments */
                   const float_complex* A,
                   const int dim1,
                   const int dim2,
                   /* Output Arguments */
                   float_complex* B);

/**
 * Row-major, general matrix pseudo-inverse (the svd way): double precision,
 * i.e.
 * \code{.m}
 *     B = pinv(A)
 * \endcode
 *
 * @param[in]  A    Input matrix; FLAT: dim1 x dim2
 * @param[in]  dim1 Number of rows in 'A' / columns in 'B'
 * @param[in]  dim2 Number of columns in 'A' / rows in 'B'
 * @param[out] B    The solution; FLAT: dim2 x dim1
 */
void utility_dpinv(/* Input Arguments */
                   const double* A,
                   const int dim1,
                   const int dim2,
                   /* Output Arguments */
                   double* B);

/**
 * Row-major, general matrix pseudo-inverse (the svd way): double precision
 * complex, i.e.
 * \code{.m}
 *     B = pinv(A)
 * \endcode
 *
 * @param[in]  A    Input matrix; FLAT: dim1 x dim2
 * @param[in]  dim1 Number of rows in 'A' / columns in 'B'
 * @param[in]  dim2 Number of columns in 'A' / rows in 'B'
 * @param[out] B    The solution; FLAT: dim2 x dim1
 */
void utility_zpinv(/* Input Arguments */
                   const double_complex* A,
                   const int dim1,
                   const int dim2,
                   /* Output Arguments */
                   double_complex* B);


/* ========================================================================== */
/*                       Cholesky Factorisation (?chol)                       */
/* ========================================================================== */

/**
 * Row-major, Cholesky factorisation of a symmetric matrix positive-definate
 * matrix: single precision, i.e.
 * \code{.m}
 *     X = chol(A); where A = X.'*X
 * \endcode
 *
 * @param[in]  A   Input square symmetric positive-definate matrix;
 *                 FLAT: dim x dim
 * @param[in]  dim Number of rows/colums in 'A'
 * @param[out] X   The solution; FLAT: dim x dim
 */
void utility_schol(/* Input Arguments */
                   const float* A,
                   const int dim,
                   /* Output Arguments */
                   float* X);

/**
 * Row-major, Cholesky factorisation of a hermitian matrix positive-definate
 * matrix: single precision complex, i.e.
 * \code{.m}
 *     X = chol(A); where A = X.'*X
 * \endcode
 *
 * @param[in]  A   Input square symmetric positive-definate matrix;
 *                 FLAT: dim x dim
 * @param[in]  dim Number of rows/colums in 'A'
 * @param[out] X   The solution; FLAT: dim x dim
 */
void utility_cchol(/* Input Arguments */
                   const float_complex* A,
                   const int dim,
                   /* Output Arguments */
                   float_complex* X);

    
/* ========================================================================== */
/*                           Matrix Inversion (?inv)                          */
/* ========================================================================== */

/**
 * Row-major, matrix inversion: single precision, i.e.
 * \code{.m}
 *     B = inv(A);
 * \endcode
 *
 * @param[in]  A   Input square matrix; FLAT: dim x dim
 * @param[out] B   Inverted square matrix; FLAT: dim x dim
 * @param[in]  dim size of matrix
 */
void utility_sinv(float* A,
                  float* B,
                  const int dim);

/**
 * Row-major, matrix inversion: double precision, i.e.
 * \code{.m}
 *     B = inv(A);
 * \endcode
 *
 * @param[in]  A   Input square matrix; FLAT: dim x dim
 * @param[out] B   Inverted square matrix; FLAT: dim x dim
 * @param[in]  dim size of matrix
 */
void utility_dinv(double* A,
                  double* B,
                  const int dim);

/**
 * Row-major, matrix inversion: double precision complex, i.e.
 * \code{.m}
 *     B = inv(A);
 * \endcode
 *
 * @param[in]  A   Input square matrix; FLAT: dim x dim
 * @param[out] B   Inverted square matrix; FLAT: dim x dim
 * @param[in]  dim size of matrix
 */
void utility_cinv(float_complex* A,
                  float_complex* B,
                  const int dim);


#ifdef __cplusplus
}/* extern "C" */
#endif  /* __cplusplus */

#endif /* SAF_VECLIB_H_INCLUDED */
