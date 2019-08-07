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

/*
 * Filename: saf_veclib.h
 * ----------------------
 * Contains wrappers for optimised linear algebra routines, utilising CBLAS and
 * LAPACK.
 *
 * Dependencies:
 *     A performance library comprising CBLAS and LAPACK functions is required
 *     by the framework.
 *     Add one of the following FLAGS to your project's preprocessor definitions
 *     list, in order to enable one of these suitable performance libraries,
 *     which must also be linked correctly to your project:
 *       SAF_USE_INTEL_MKL
 *         to enable Intel's Math Kernal Library
 *       SAF_USE_OPENBLAS_AND_REF_LAPACK
 *         to enable OpenBLAS and to use the reference implementation of LAPACK
 *     More information can be found here:
 *     https://github.com/leomccormack/Spatial_Audio_Framework
 *     Mac users only: saf_utilities will employ Apple's Accelerate library by
 *     default, if none of the above FLAGS are defined.
 * Author, date created:
 *     Leo McCormack, 11.07.2016
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
#ifdef CBLAS_H
# define NO_TRANSPOSE (CblasNoTrans)
# define TRANSPOSE (CblasTrans)
# define CONJ_TRANSPOSE (CblasConjTrans)
  typedef enum CBLAS_TRANSPOSE TRANS_FLAG;
#else
  typedef enum _TRANS_FLAG{
    NO_TRANSPOSE = 1,
    TRANSPOSE = 2,
    CONJ_TRANSPOSE = 3
  }TRANS_FLAG;
#endif
typedef enum _CONJ_FLAG{
  NO_CONJ = 1,
  CONJ = 2
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

/*
 * Function: utility_siminv
 * ------------------------
 * s, single-precision, index of minimum absolute value in a vector, i.e.
 *     [~,ind] = min(abs(a))
 *
 * Input Arguments:
 *     a     - input vector a; len x 1
 *     len   - vector length
 * Output Arguments:
 *     index - & index of minimum value; 1 x 1
 */
void utility_siminv(/* Input Arguments */
                    const float* a,
                    const int len,
                    /* Output Arguments */
                    int* index);

/*
 * Function: utility_ciminv
 * ------------------------
 * c, single-precision, complex, index of maximum absolute value in a vector,
 * i.e.
 *     [~,ind] = min(abs(a))
 *
 * Input Arguments:
 *     a     - input vector a; len x 1
 *     len   - vector length
 * Output Arguments:
 *     index - & index of minimum value; 1 x 1
 */
void utility_ciminv(/* Input Arguments */
                    const float_complex* a,
                    const int len,
                    /* Output Arguments */
                    int* index);

    
/* ========================================================================== */
/*                     Find Index of Max-Abs-Value (?imaxv)                   */
/* ========================================================================== */

/*
 * Function: utility_siminv
 * ------------------------
 * s, single-precision, index of maximum absolute value in a vector, i.e.
 *     [~,ind] = max(abs(a))
 *
 * Input Arguments:
 *     a     - input vector a; len x 1
 *     len   - vector length
 * Output Arguments:
 *     index - & index of maximum value; 1 x 1
 */
void utility_simaxv(/* Input Arguments */
                    const float* a,
                    const int len,
                    /* Output Arguments */
                    int* index);

/*
 * Function: utility_ciminv
 * ------------------------
 * c, single-precision, complex, index of maximum absolute value in a vector,
 * i.e.
 *     [~,ind] = max(abs(a))
 *
 * Input Arguments:
 *     a     - input vector a; len x 1
 *     len   - vector length
 * Output Arguments:
 *     index - & index of maximum value; 1 x 1
 */
void utility_cimaxv(/* Input Arguments */
                    const float_complex* a,
                    const int len,
                    /* Output Arguments */
                    int* index);


/* ========================================================================== */
/*                              Vector-Abs (?vabs)                            */
/* ========================================================================== */

/*
 * Function: utility_svabs
 * -----------------------
 * s, single-precision, absolute values of vector elements, i.e.
 *     c = abs(a)
 *
 * Input Arguments:
 *     a   - input vector a; len x 1
 *     len - vector length
 * Output Arguments:
 *     c   - output vector c; len x 1
 */
void utility_svabs(/* Input Arguments */
                   const float* a,
                   const int len,
                   /* Output Arguments */
                   float* c);

/*
 * Function: utility_cvabs
 * -----------------------
 * c, single-precision, complex, absolute values of vector elements, i.e.
 *     c = abs(a)
 *
 * Input Arguments:
 *     a   - input vector a; len x 1
 *     len - vector length
 * Output Arguments:
 *     c   - output vector c; len x 1
 */
void utility_cvabs(/* Input Arguments */
                   const float_complex* a,
                   const int len,
                   /* Output Arguments */
                   float* c);


/* ========================================================================== */
/*                        Vector-Vector Copy (?vvcopy)                        */
/* ========================================================================== */

/*
 * Function: utility_svvcopy
 * -------------------------
 * s, single-precision, vector-vector copy, i.e.
 *     c = a
 *
 * Input Arguments:
 *     a   - input vector a; len x 1
 *     len - vector length
 * Output Arguments:
 *     c   - output vector c; len x 1
 */
void utility_svvcopy(/* Input Arguments */
                     const float* a,
                     const int len,
                     /* Output Arguments */
                     float* c);

/*
 * Function: utility_cvvcopy
 * -------------------------
 * c, single-precision, complex, vector-vector copy, i.e.
 *     c = a
 *
 * Input Arguments:
 *     a   - input vector a; len x 1
 *     len - vector length
 * Output Arguments:
 *     c   - output vector c; len x 1
 */
void utility_cvvcopy(/* Input Arguments */
                     const float_complex* a,
                     const int len,
                     /* Output Arguments */
                     float_complex* c);


/* ========================================================================== */
/*                       Vector-Vector Addition (?vvadd)                      */
/* ========================================================================== */

/*
 * Function: utility_svvadd
 * ------------------------
 * s, single-precision, vector-vector addition, i.e.
 *     c = a+b, OR: a = a+b (if c=NULL)
 *
 * Input Arguments:
 *     a   - input vector a, and output if c==NULL; len x 1
 *     b   - input vector b; len x 1
 *     len - vector length
 * Output Arguments:
 *     c   - output vector c (set to NULL if you want 'a' as output); len x 1
 */
void utility_svvadd(/* Input Arguments */
                    float* a,
                    const float* b,
                    const int len,
                    /* Output Arguments */
                    float* c);

/*
 * Function: utility_cvvadd
 * ------------------------
 * c, single-precision, complex, vector-vector addition, i.e.
 *     c = a+b, OR: a = a+b (if c=NULL)
 *
 * Input Arguments:
 *     a   - input vector a, and output if c==NULL; len x 1
 *     b   - input vector b; len x 1
 *     len - vector length
 * Output Arguments:
 *     c   - output vector c (set to NULL if you want 'a' as output); len x 1
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

/*
 * Function: utility_svvsub
 * ------------------------
 * s, single-precision, vector-vector subtraction, i.e.
 *     c = a-b, OR: a = a-b (if c=NULL)
 *
 * Input Arguments:
 *     a   - input vector a, and output if c==NULL; len x 1
 *     b   - input vector b; len x 1
 *     len - vector length
 * Output Arguments:
 *     c   - output vector c (set to NULL if you want 'a' as output); len x 1
 */
void utility_svvsub(/* Input Arguments */
                    float* a,
                    const float* b,
                    const int len,
                    /* Output Arguments */
	                float* c);

/*
 * Function: utility_cvvsub
 * ------------------------
 * c, single-precision, complex, vector-vector subtraction, i.e.
 *     c = a-b, OR: a = a-b (if c=NULL)
 *
 * Input Arguments:
 *     a   - input vector a, and output if c==NULL; len x 1
 *     b   - input vector b; len x 1
 *     len - vector length
 * Output Arguments:
 *     c   - output vector c (set to NULL if you want 'a' as output); len x 1
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

/*
 * Function: utility_svvmul
 * ------------------------
 * s, single-precision, element-wise vector-vector multiplication i.e.
 *     c = a.*b, OR: a = a.*b (if c=NULL)
 *
 * Input Arguments:
 *     a   - input vector a, and output if c==NULL; len x 1
 *     b   - input vector b; len x 1
 *     len - vector length
 * Output Arguments:
 *     c   - output vector c (set to NULL if you want 'a' as output); len x 1
 */
void utility_svvmul(/* Input Arguments */
                    float* a,
                    const float* b,
                    const int len,
                    /* Output Arguments */
                    float* c);

/*
 * Function: utility_cvvmul
 * ------------------------
 * c, single-precision, complex, element-wise vector-vector multiplication i.e.
 *     c = a.*b, OR: a = a.*b (if c=NULL)
 *
 * Input Arguments:
 *     a   - input vector a, and output if c==NULL; len x 1
 *     b   - input vector b; len x 1
 *     len - vector length
 * Output Arguments:
 *     c   - output vector c (set to NULL if you want 'a' as output); len x 1
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

/*
 * Function: utility_svvdot
 * ------------------------
 * s, single-precision, vector-vector dot product, i.e.
 *     c = a*b^T, (where size(c) = 1 x 1)
 *
 * Input Arguments:
 *     a   - input vector a; len x 1
 *     b   - input vector b; len x 1
 *     len - vector length
 * Output Arguments:
 *     c   - & output vector c; 1 x 1
 */
void utility_svvdot(/* Input Arguments */
                    const float* a,
                    const float* b,
                    const int len,
                    /* Output Arguments */
                    float* c);

/*
 * Function: utility_cvvdot
 * ------------------------
 * c, single-precision, complex, vector-vector dot product, i.e.
 *     c = a*b^T, (where size(c) = 1 x 1)
 *
 * Input Arguments:
 *     a    - input vector a; len x 1
 *     b    - input vector b; len x 1
 *     len  - vector length
 *     flag - 0: do not conjugate a, 1: conjuate a
 * Output Arguments:
 *     c    - & output vector c; 1 x 1
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

/*
 * Function: utility_svsmul
 * ------------------------
 * s, single-precision, multiplies each element in vector 'a' with a scalar 's',
 * i.e.
 *     c = a.*s, OR: a = a.*s (if c=NULL)
 *
 * Input Arguments:
 *     a   - input vector a, and output if c==NULL; len x 1
 *     s   - & input scalar s; 1 x 1
 *     len - vector length
 * Output Arguments:
 *     c   - output vector c (set to NULL if you want 'a' as output); len x 1
 */
void utility_svsmul(/* Input Arguments */
                    float* a,
                    const float* s,
                    const int len,
                    /* Output Arguments */
                    float* c);

/*
 * Function: utility_cvsmul
 * ------------------------
 * c, single-precision, complex, multiplies each element in vector 'a' with a
 * scalar 's', i.e.
 *     c = a.*s, OR: a = a.*s (if c=NULL)
 *
 * Input Arguments:
 *     a   - input vector a, and output if c==NULL; len x 1
 *     s   - & input scalar s; 1 x 1
 *     len - vector length
 * Output Arguments:
 *     c   - output vector c (set to NULL if you want 'a' as output); len x 1
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

/*
 * Function: utility_svsdiv
 * ------------------------
 * s, single-precision, divides each element in vector 'a' with a scalar 's',
 * i.e.
 *     c = a./s, OR: a = a./s (if c=NULL)
 *
 * Input Arguments:
 *     a   - input vector a, and output if c==NULL; len x 1
 *     s   - & input scalar s; 1 x 1
 *     len - vector length
 * Output Arguments:
 *     c   - output vector c (set to NULL if you want 'a' as output); len x 1
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

/*
 * Function: utility_svsadd
 * ------------------------
 * s, single-precision, adds each element in vector 'a' with a scalar 's', i.e.
 *     c = a+s, OR: a = a+s (if c=NULL)
 *
 * Input Arguments:
 *     a   - input vector a, and output if c==NULL; len x 1
 *     s   - & input scalar s; 1 x 1
 *     len - vector length
 * Output Arguments:
 *     c   - output vector c (set to NULL if you want 'a' as output); len x 1
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
    
/*
 * Function: utility_svssub
 * ------------------------
 * s, single-precision, subtracts each element in vector 'a' with a scalar 's',
 * i.e.
 *     c = a-s, OR: a = a-s (if c=NULL)
 *
 * Input Arguments:
 *     a   - input vector a, and output if c==NULL; len x 1
 *     s   - & input scalar s; 1 x 1
 *     len - vector length
 * Output Arguments:
 *     c   - output vector c (set to NULL if you want 'a' as output); len x 1
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
    
/*
 * Function: utility_ssvd
 * ----------------------
 * s, row-major, singular value decomposition: single precision, i.e.
 *     [U,S,V] = svd(A); such that A = U*S*V.' = U*diag(sing)*V.'
 * Note 'S' contains the singular values along the diagonal, whereas 'sing' are
 * the singular values as a vector
 * Further Note: V is returned untransposed! (like in Matlab)
 *
 * Input Arguments:
 *     A    - input matrix; FLAT: dim1 x dim2
 *     dim1 - first dimension of matrix 'A'
 *     dim2 - second dimension of matrix 'A'
 * Output Arguments:
 *     U    - left matrix (set to NULL if not needed); FLAT: dim1 x dim1
 *     S    - singular values along the diagonal min(dim1, dim2), (set to NULL
 *            if not needed); FLAT: dim1 x dim2
 *     V    - right matrix (UNTRANSPOSED!) (set to NULL if not needed);
 *            FLAT: dim2 x dim2
 *     sing - singular values as a vector, (set to NULL if not needed);
 *            min(dim1, dim2) x 1
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

/*
 * Function: utility_csvd
 * ----------------------
 * c, row-major, singular value decomposition: single precision complex, i.e.
 *     [U,S,V] = svd(A); such that A = U*S*V' = U*diag(sing)*V'
 * Note 'S' contains the singular values along the diagonal, while 'sing' are
 * the singular values as a vector
 * Further Note: V is returned untransposed!
 *
 * Input Arguments:
 *     A    - input matrix; FLAT: dim1 x dim2
 *     dim1 - first dimension of matrix 'A'
 *     dim2 - second dimension of matrix 'A'
 * Output Arguments:
 *     U    - left matrix (set to NULL if not needed); FLAT: dim1 x dim1
 *     S    - singular values along the diagonal min(dim1, dim2), (set to NULL
 *            if not needed); FLAT: dim1 x dim2
 *     V    - right matrix (UNTRANSPOSED!) (set to NULL if not needed);
 *            FLAT: dim2 x dim2
 *     sing - singular values as a vector, (set to NULL if not needed);
 *            min(dim1, dim2) x 1
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

/*
 * Function: utility_sseig
 * -----------------------
 * s, row-major, eigenvalue decomposition of a SYMMETRIC matrix: single
 * precision, i.e.
 *     [V,D] = eig(A); where A*V = V*D, and  A*V = V*diag(eig)
 * Note: 'D' contains the eigen values along the diagonal, while 'eig' are the
 * eigen values as a vector
 *
 * Input Arguments:
 *     A           - input SYMMETRIC square matrix; FLAT: dim x dim
 *     dim         - dimensions for square matrix 'A'
 *     sortDecFLAG - 1: sort eigen values and vectors in decending order.
 *                   0: ascending
 * Output Arguments:
 *     V           - Eigen vectors (set to NULL if not needed); FLAT: dim x dim
 *     D           - Eigen values along the diagonal (set to NULL if not
 *                   needed); FLAT: dim x dim
 *     eig         - Eigen values not diagonalised (set to NULL if not needed);
 *                   dim x 1
 */
void utility_sseig(/* Input Arguments */
                   const float* A,
                   const int dim,
                   int sortDecFLAG,
                   /* Output Arguments */
                   float* V,
                   float* D,
                   float* eig);

/*
 * Function: utility_cseig
 * -----------------------
 * c, row-major, eigenvalue decomposition of a SYMMETRIC/HERMITION matrix:
 * single precision complex, i.e.
 *     [V,D] = eig(A); where A*V = V*D, and  A*V = V*diag(eig)
 * Note: 'D' contains the eigen values along the diagonal, while 'eig' are the
 * eigen values as a vector
 *
 * Input Arguments:
 *     A           - input SYMMETRIC/HERMITION square matrix; FLAT: dim x dim
 *     dim         - dimensions for square matrix 'A'
 *     sortDecFLAG - 1: sort eigen values and vectors in decending order.
 *                   0: ascending
 * Output Arguments:
 *     V           - Eigen vectors (set to NULL if not needed); FLAT: dim x dim
 *     D           - Eigen values along the diagonal (set to NULL if not
 *                   needed); FLAT: dim x dim
 *     eig         - Eigen values not diagonalised (set to NULL if not needed);
 *                   dim x 1
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
 
/*
 * Function: utility_ceigmp
 * ------------------------
 * c, row-major, finds eigenvalues of a matrix pair using the QZ method, single
 * precision complex, i.e.
 *     [VL,VR,D] = eig(A,B,'qz'); where A*VL = B*VL*VR
 *
 * Input Arguments:
 *     A   - input left square matrix; FLAT: dim x dim
 *     B   - input right square matrix; FLAT: dim x dim
 *     dim - dimensions for square matrices 'A' and 'B'
 * Output Arguments:
 *     VL  - Left Eigen vectors (set to NULL if not needed); FLAT: dim x dim
 *     VR  - Right Eigen vectors (set to NULL if not needed); FLAT: dim x dim
 *     D   - Eigen values along the diagonal (set to NULL if not needed);
 *           FLAT: dim x dim
 */
void utility_ceigmp(/* Input Arguments */
                    const float_complex* A,
                    const float_complex* B,
                    const int dim,
                    /* Output Arguments */
                    float_complex* VL,
                    float_complex* VR,
                    float_complex* D);
    
/*
 * Function: utility_zeigmp
 * ------------------------
 * z, row-major, finds eigenvalues of a matrix pair using the QZ method, double
 * precision complex, i.e.
 *     [VL,VR,D] = eig(A,B,'qz'); where A*VL = B*VL*VR
 *
 * Input Arguments:
 *     A   - input left square matrix; FLAT: dim x dim
 *     B   - input right square matrix; FLAT: dim x dim
 *     dim - dimensions for square matrices 'A' and 'B'
 * Output Arguments:
 *     VL  - Left Eigen vectors (set to NULL if not needed); FLAT: dim x dim
 *     VR  - Right Eigen vectors (set to NULL if not needed); FLAT: dim x dim
 *     D   - Eigen values along the diagonal (set to NULL if not needed);
 *           FLAT: dim x dim
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

/*
 * Function: utility_ceig
 * ----------------------
 * c, row-major, eigenvalue decomposition of a NON-SYMMETRIC matrix: single
 * precision complex, i.e.
 *     [VL,VR,D] = eig(A); where A*VR = VR*D, and  A*VR = VR*diag(eig)
 * Note: 'D' contains the eigen values along the diagonal, while 'eig' are the
 * eigen values as a vector
 *
 * Input Arguments:
 *     A           - input NON-SYMMETRIC square matrix; FLAT: dim x dim
 *     dim         - dimensions for square matrix 'A'
 *     sortDecFLAG - 1: sort eigen values and vectors in decending order.
 *                   0: ascending
 * Output Arguments:
 *     VL          - Left Eigen vectors (set to NULL if not needed);
 *                   FLAT: dim x dim
 *     VR          - Right Eigen vectors (set to NULL if not needed);
 *                   FLAT: dim x dim
 *     D           - Eigen values along the diagonal (set to NULL if not
 *                   needed); FLAT: dim x dim
 *     eig         - Eigen values not diagonalised (set to NULL if not needed);
 *                   dim x 1
 */
void utility_ceig(/* Input Arguments */
                  const float_complex* A,
                  const int dim,
                  int sortDecFLAG,
                  /* Output Arguments */
                  float_complex* VL,
                  float_complex* VR,
                  float_complex* D,
                  float* eig);


/* ========================================================================== */
/*                       General Linear Solver (?glslv)                       */
/* ========================================================================== */

/*
 * Function: utility_sglslv
 * ------------------------
 * s, row-major, general linear solver: single precision, i.e.
 *     X = linsolve(A,B) = A\B; where, AX = B
 *
 * Input Arguments:
 *     A    - input square matrix; FLAT: dim x dim
 *     dim  - dimensions for square matrix 'A'
 *     B    - right hand side matrix; FLAT: dim x nCol
 *     nCol - number of columns in right hand side matrix
 * Output Arguments:
 *     X    - the solution; FLAT: dim x nCol
 */
void utility_sglslv(/* Input Arguments */
                    const float* A,
                    const int dim,
                    float* B,
                    int nCol,
                    /* Output Arguments */
                    float* X);

/*
 * Function: utility_cglslv
 * ------------------------
 * c, row-major, general linear solver: single precision complex, i.e.
 *     X = linsolve(A,B) = A\B; where, AX = B
 *
 * Input Arguments:
 *     A    - input square matrix; FLAT: dim x dim
 *     dim  - dimensions for square matrix 'A'
 *     B    - right hand side matrix; FLAT: dim x nCol
 *     nCol - number of columns in right hand side matrix
 * Output Arguments:
 *     X    - the solution; FLAT: dim x nCol
 */
void utility_cglslv(/* Input Arguments */
                    const float_complex* A,
                    const int dim,
                    float_complex* B,
                    int nCol,
                    /* Output Arguments */
                    float_complex* X);

/*
 * Function: utility_dglslv
 * ------------------------
 * d, row-major, general linear solver: double precision, i.e.
 *     X = linsolve(A,B) = A\B; where, AX = B
 *
 * Input Arguments:
 *     A    - input square matrix; FLAT: dim x dim
 *     dim  - dimensions for square matrix 'A'
 *     B    - right hand side matrix; FLAT: dim x nCol
 *     nCol - number of columns in right hand side matrix
 * Output Arguments:
 *     X    - the solution; FLAT: dim x nCol
 */
void utility_dglslv(/* Input Arguments */
                    const double* A,
                    const int dim,
                    double* B,
                    int nCol,
                    /* Output Arguments */
                    double* X);

/*
 * Function: utility_zglslv
 * ------------------------
 * z, row-major, general linear solver: double precision complex, i.e.
 *     X = linsolve(A,B) = A\B; where, AX = B
 *
 * Input Arguments:
 *     A    - input square matrix; FLAT: dim x dim
 *     dim  - dimensions for square matrix 'A'
 *     B    - right hand side matrix; FLAT: dim x nCol
 *     nCol - number of columns in right hand side matrix
 * Output Arguments:
 *     X    - the solution; FLAT: dim x nCol
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

/*
 * Function: utility_sslslv
 * ------------------------
 * s, row-major, linear solver for SYMMETRIC positive-definate 'A': single
 * precision, i.e.
 *     opts.LT=true
 *     X = linsolve(A,B, opts); where, AX = B
 *
 * Input Arguments:
 *     A    - input square SYMMETRIC positive-definate matrix; FLAT: dim x dim
 *     dim  - dimensions for square matrix 'A'
 *     B    - right hand side matrix; FLAT: dim x nCol
 *     nCol - number of columns in right hand side matrix
 * Output Arguments:
 *     X    - the solution; FLAT: dim x nCol
 */
void utility_sslslv(/* Input Arguments */
                    const float* A,
                    const int dim,
                    float* B,
                    int nCol,
                    /* Output Arguments */
                    float* X);

/*
 * Function: utility_cslslv
 * ------------------------
 * c, row-major, linear solver for HERMITIAN positive-definate 'A': single
 * precision complex, i.e.
 *     opts.LT=true
 *     X = linsolve(A,B, opts); where, AX = B
 *
 * Input Arguments:
 *     A    - input square HERMITIAN positive-definate matrix; FLAT: dim x dim
 *     dim  - dimensions for square matrix 'A'
 *     B    - right hand side matrix; FLAT: dim x nCol
 *     nCol - number of columns in right hand side matrix
 * Output Arguments:
 *     X    - the solution; FLAT: dim x nCol
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

/*
 * Function: utility_spinv
 * -----------------------
 * s, row-major, general matrix pseudo-inverse (the svd way): single precision,
 * i.e.
 *     B = pinv(A)
 *
 * Input Arguments:
 *     A    - input matrix; FLAT: dim1 x dim2
 *     dim1 - number of rows in 'A' / columns in 'B'
 *     dim2 - number of columns in 'A' / rows in 'B'
 * Output Arguments:
 *     B    - the solution; FLAT: dim2 x dim1
 */
void utility_spinv(/* Input Arguments */
                   const float* A,
                   const int dim1,
                   const int dim2,
                   /* Output Arguments */
                   float* B);

/*
 * Function: utility_cpinv
 * -----------------------
 * c, row-major, general matrix pseudo-inverse (the svd way): single precision
 * complex, i.e.
 *     B = pinv(A)
 *
 * Input Arguments:
 *     A    - input matrix; FLAT: dim1 x dim2
 *     dim1 - number of rows in 'A' / columns in 'B'
 *     dim2 - number of columns in 'A' / rows in 'B'
 * Output Arguments:
 *     B    - the solution; FLAT: dim2 x dim1
 */
void utility_cpinv(/* Input Arguments */
                   const float_complex* A,
                   const int dim1,
                   const int dim2,
                   /* Output Arguments */
                   float_complex* B);

/*
 * Function: utility_dpinv
 * -----------------------
 * d, row-major, general matrix pseudo-inverse (the svd way): double precision,
 * i.e.
 *     B = pinv(A)
 *
 * Input Arguments:
 *     A    - input matrix; FLAT: dim1 x dim2
 *     dim1 - number of rows in 'A' / columns in 'B'
 *     dim2 - number of columns in 'A' / rows in 'B'
 * Output Arguments:
 *     B    - the solution; FLAT: dim2 x dim1
 */
void utility_dpinv(/* Input Arguments */
                   const double* A,
                   const int dim1,
                   const int dim2,
                   /* Output Arguments */
                   double* B);

/*
 * Function: utility_zpinv
 * -----------------------
 * z, row-major, general matrix pseudo-inverse (the svd way): double precision
 * complex, i.e.
 *     B = pinv(A)
 *
 * Input Arguments:
 *     A    - input matrix; FLAT: dim1 x dim2
 *     dim1 - number of rows in 'A' / columns in 'B'
 *     dim2 - number of columns in 'A' / rows in 'B'
 * Output Arguments:
 *     B    - the solution; FLAT: dim2 x dim1
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

/*
 * Function: utility_schol
 * -----------------------
 * s, row-major, Cholesky factorisation of a symmetric matrix positive-definate
 * matrix: single precision, i.e.
 *     X = chol(A); where A = X.'*X
 *
 * Input Arguments:
 *     A   - input square symmetric positive-definate matrix; FLAT: dim x dim
 *     dim - number of rows/colums in 'A'
 * Output Arguments:
 *     X   - the solution; FLAT: dim x dim
 */
void utility_schol(/* Input Arguments */
                   const float* A,
                   const int dim,
                   /* Output Arguments */
                   float* X);

/*
 * Function: utility_cchol
 * -----------------------
 * c, row-major, Cholesky factorisation of a hermitian matrix positive-definate
 * matrix: single precision complex, i.e.
 *     X = chol(A); where A = X.'*X
 *
 * Input Arguments:
 *     A   - input square hermitian positive-definate matrix; FLAT: dim x dim
 *     dim - number of rows/colums in 'A'
 * Output Arguments:
 *     X   - the solution; FLAT: dim x dim
 */
void utility_cchol(/* Input Arguments */
                   const float_complex* A,
                   const int dim,
                   /* Output Arguments */
                   float_complex* X);

    
/* ========================================================================== */
/*                           Matrix Inversion (?inv)                          */
/* ========================================================================== */
 
//DEPRECATED
//TODO: rewrite for row-major:

/* s, column-major, matrix inversion: single precision */
void utility_sinv(float * A,
                  const int  N);

/* d, column-major, matrix inversion: double precision */
void utility_dinv(double* A,
                  const int N);

/* c, column-major, matrix inversion: single precision complex */
void utility_cinv(float_complex * A,
                  const int N );

    
#ifdef __cplusplus
}/* extern "C" */
#endif  /* __cplusplus */

#endif /* SAF_VECLIB_H_INCLUDED */
