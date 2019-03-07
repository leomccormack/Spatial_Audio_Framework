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
 *     saf_veclib.h
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

#ifndef SAF_VECLIB_H_INCLUDED
#define SAF_VECLIB_H_INCLUDED

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
/* a performance library is required: */
#if defined(__APPLE__) && !defined(SAF_USE_INTEL_MKL)
  /* it is highly recommended to use (at least) the Accelerate performance library for Mac OSX */
  #include "Accelerate/Accelerate.h"
#else
  /* it is highly recommended to use Intel's MKL performance library for Windows & Mac OSX */
  #include "mkl.h"
#endif
#include "saf_complex.h"
#ifdef CBLAS_H
  #define NO_TRANSPOSE (CblasNoTrans)
  #define TRANSPOSE (CblasTrans)
  #define CONJ_TRANSPOSE (CblasConjTrans)
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

/*----------------------------- index of min-abs-value (?imaxv) -----------------------------*/

/* s, single-precision, index of minimum absolute value in a vector, [~,ind] = min(abs(a) */
void utility_siminv(const float* a,                /* input vector a; len x 1 */
                    const int len,                 /* vector length */
                    int* index);                   /* & index of value; 1 x 1 */

/* s, single-precision, complex, index of maximum absolute value in a vector, [~,ind] = min(abs(a) */
void utility_ciminv(const float_complex* a,        /* input vector a; len x 1 */
                    const int len,                 /* vector length */
                    int* index);                   /* & index of value; 1 x 1 */

/*----------------------------- index of max-abs-value (?imaxv) -----------------------------*/

/* s, single-precision, index of maximum absolute value in a vector, [~,ind] = max(abs(a) */
void utility_simaxv(const float* a,                /* input vector a; len x 1 */
                    const int len,                 /* vector length */
                    int* index);                   /* & index of value; 1 x 1 */

/* s, single-precision, complex, index of maximum absolute value in a vector, [~,ind] = max(abs(a) */
void utility_cimaxv(const float_complex* a,        /* input vector a; len x 1 */
                    const int len,                 /* vector length */
                    int* index);                   /* & index of value; 1 x 1 */

/*----------------------------------- vector-abs (?vabs) ------------------------------------*/

/* s, single-precision, absolute value of vector elements, c = fabsf(a) */
void utility_svabs(const float* a,                 /* input vector a; len x 1 */
                   const int len,                  /* vector length */
                   float* c);                      /* output vector c; len x 1 */

/* s, single-precision, complex, absolute value of vector elements, c = cabsf(a) */
void utility_cvabs(const float_complex* a,         /* input vector a; len x 1 */
                   const int len,                  /* vector length */
                   float* c);                      /* output vector c; len x 1 */

/*------------------------------ vector-vector copy (?vvcopy) -------------------------------*/

/* s, single-precision, vector-vector copy, c = a */
void utility_svvcopy(const float* a,               /* input vector a; len x 1 */
                     const int len,                /* vector length */
                     float* c);                    /* output vector c; len x 1 */

/* c, single-precision, complex, vector-vector copy, c = a */
void utility_cvvcopy(const float_complex* a,       /* input vector a; len x 1 */
                     const int len,                /* vector length */
                     float_complex* c);            /* output vector c; len x 1 */

/*------------------------- vector-vector multiplication (?vvmul) ---------------------------*/

/* s, single-precision, element-wise vector-vector multiplication, c = a.*b || a = a.*b */
void utility_svvmul(float* a,                      /* input vector a, and output if c==NULL; len x 1 */
                    const float* b,                /* input vector b; len x 1 */
                    const int len,                 /* length of vectors a,b,c */
                    float* c);                     /* output vector c; len x 1 */

/*--------------------------- vector-vector dot product (?vvdot) ----------------------------*/

/* s, single-precision, vector-vector dot product */
void utility_svvdot(const float* a,              /* input vector a; len x 1 */
                    const float* b,              /* input vector b; len x 1 */
                    const int len,               /* length of vectors a,b,c */
                    float* c);                   /* output vector c; len x 1 */

/* c, single-precision, complex, vector-vector dot product */
void utility_cvvdot(const float_complex* a,      /* input vector a; len x 1 */
                    const float_complex* b,      /* input vector b; len x 1 */
                    const int len,               /* length of vectors a,b,c */
                    CONJ_FLAG flag,              /* flag, should vector a should be conjugated */
                    float_complex* c);           /* output vector c; len x 1 */

/*------------------------------ vector-scalar product (?vsmul) -----------------------------*/

/* s, single-precision, multiplies each element in vector 'a' with a scalar 's' */
void utility_svsmul(float* a,               /* input vector, and output if c==NULL; len x 1 */
                    const float* s,         /* & scalar; 1 x 1 */
                    const int len,          /* vector length */
                    float* c);              /* output if c!=NULL; len x 1  */

/* c, single-precision, complex, multiplies each element in vector 'a' with a scalar 's' */
void utility_cvsmul(float_complex* a,       /* input vector, and output if c==NULL; len x 1 */
                    const float_complex* s, /* & scalar; 1 x 1 */
                    const int len,          /* vector length */
                    float_complex* c);      /* output if c!=NULL; len x 1  */

/*----------------------------- vector-scalar division (?vsdiv) -----------------------------*/

/* s, single-precision, divides each element in vector 'a' with a scalar 's' */
void utility_svsdiv(float* a,
                    const float* s,
                    const int len,
                    float* c);

/*----------------------------- vector-scalar addition (?vsadd) -----------------------------*/

/* s, single-precision, adds a scalar 's' to each element in vector 'a' */
void utility_svsadd(float* a,
                    const float* s,
                    const int len,
                    float* c);

/*---------------------------- vector-scalar subtraction (?vssub) ---------------------------*/

/* s, single-precision, subtracts a scalar 's' from each element in vector 'a' */
void utility_svssub(float* a,              /* vector; len x 1 */
                    const float* s,        /* scalar */
                    const int len,         /* length of vector */
                    float* c);             /* c = a-s, set to NULL for a = a-s; len x 1 */

/*---------------------------- singular-value decomposition (?svd) --------------------------*/

/* s, row-major, singular value decomposition: single precision */
void utility_ssvd(const float* A,          /* in matrix; flat: dim1 x dim2 */
                  const int dim1,          /* first dimension of A */
                  const int dim2,          /* second dimension of A */
                  float* U,                /* left matrix (set to NULL if not needed); flat: dim1 x dim1 */
                  float* S,                /* singular values along the diagonal min(dim1, dim2), (set to NULL if not needed); flat: dim1 x dim2 */
                  float* V,                /* right matrix (UNTRANSPOSED!) (set to NULL if not needed); flat: dim2 x dim2 */
                  float* sing);            /* singular values as a vector, (set to NULL if not needed); min(dim1, dim2) x 1 */

/* s, row-major, singular value decomposition: single precision complex */
void utility_csvd(const float_complex* A,  /* in matrix; flat: dim1 x dim2 */
                  const int dim1,          /* first dimension of A */
                  const int dim2,          /* second dimension of A */
                  float_complex* U,        /* left matrix (set to NULL if not needed); flat: dim1 x dim1 */
                  float_complex* S,        /* singular values along the diagonal min(dim1, dim2), (set to NULL if not needed); flat: dim1 x dim2 */
                  float_complex* V,        /* right matrix (UNTRANSPOSED!) (set to NULL if not needed); flat: dim2 x dim2 */
                  float* sing);            /* singular values as a vector, (set to NULL if not needed); min(dim1, dim2) x 1 */

/*------------------------ symmetric eigenvalue decomposition (?seig) -----------------------*/

/* s, row-major, eigenvalue decomposition of a symmetric matrix: single precision */
void utility_sseig(const float* A,         /* in symmetric square matrix; flat: dim x dim */
                   const int dim,          /* dimensions for the square matrix, A */
                   int sortDecFLAG,        /* 1: sort eigen values and vectors in decending order. 0: ascending */
                   float* V,               /* Eigen vectors (set to NULL if not needed); dim x dim */
                   float* D,               /* Eigen values along the diagonal (set to NULL if not needed); dim x dim */
                   float* eig);            /* Eigen values not diagonalised (set to NULL if not needed); dim x 1 */

/* c, row-major, eigenvalue decomposition of a symmetric/Hermitian matrix: single precision complex */
void utility_cseig(const float_complex* A, /* in symmetric square matrix; flat: dim x dim */
                   const int dim,          /* dimensions for the square matrix, A */
                   int sortDecFLAG,        /* 1: sort eigen values and vectors in decending order. 0: ascending */
                   float_complex* V,       /* Eigen vectors (set to NULL if not needed); dim x dim */
                   float_complex* D,       /* Eigen values along the diagonal (set to NULL if not needed); dim x dim */
                   float* eig);            /* Eigen values not diagonalised (set to NULL if not needed); dim x 1 */

/*-----------------------------  eigenvalue decomposition (?eig) ----------------------------*/

/* s, row-major, eigenvalue decomposition of a non-symmetric matrix: single precision complex */
void utility_ceig(const float_complex* A,  /* in nonsymmetric square matrix; flat: dim x dim */
                  const int dim,           /* dimensions for the square matrix, A */
                  int sortDecFLAG,         /* 1: sort eigen values and vectors in decending order. 0: ascending */
                  float_complex* VL,       /* Left Eigen vectors (set to NULL if not needed); dim x dim */
                  float_complex* VR,       /* Right Eigen vectors (set to NULL if not needed); dim x dim */
                  float_complex* D,        /* Eigen values along the diagonal (set to NULL if not needed); dim x dim */
                  float* eig);             /* Eigen values not diagonalised (set to NULL if not needed); dim x 1 */

/*------------------------------ general linear solver (?glslv) -----------------------------*/

/* s, row-major, general linear solver (AX=B): single precision */
void utility_sglslv(const float* A,          /* input square matrix; flat: dim x dim */
                    const int dim,           /* dimensions for the square matrix, A */
                    float* B,                /* right hand side matrix; flat: dim x nCol */
                    int nCol,                /* number of columns in right hand side matrix */
                    float* X);               /* the solution; dim x nCol */

/* c, row-major, general linear solver (AX=B): single precision complex */
void utility_cglslv(const float_complex* A,  /* input square matrix; flat: dim x dim */
                    const int dim,           /* dimensions for the square matrix, A */
                    float_complex* B,        /* right hand side matrix; flat: dim x nCol */
                    int nCol,                /* number of columns in right hand side matrix */
                    float_complex* X);       /* the solution; dim x nCol */

/*----------------------------- symmetric linear solver (?slslv) ----------------------------*/

/* s, row-major, linear solver (AX=B) for symmetric positive-definate 'A': single precision */
void utility_sslslv(const float* A,          /* square symmetric positive-definate matrix; flat: dim x dim */
                    const int dim,           /* dimensions for the square matrix, A */
                    float* B,                /* right hand side matrix; flat: dim x nCol */
                    int nCol,                /* number of columns in right hand side matrix */
                    float* X);               /* the solution; dim x nCol */

/* c, row-major, linear solver (AX=B) for symmetric positive-definate 'A': single precision complex */
void utility_cslslv(const float_complex* A,  /* square symmetric positive-definate matrix; flat: dim x dim */
                    const int dim,           /* dimensions for the square matrix, A */
                    float_complex* B,        /* right hand side matrix; flat: dim x nCol */
                    int nCol,                /* number of columns in right hand side matrix */
                    float_complex* X);       /* the solution; dim x nCol */

/*------------------------------- matrix pseudo-inverse (?pinv) -----------------------------*/

/* s, row-major, general matrix pseudo-inverse (the svd way): single precision */
void utility_spinv(const float* inM,          /* in matrix; flat:[dim1][dim2] */
                   const int dim1,
                   const int dim2,
                   float* outM);              /* out matrix; flat:[dim2][dim1] */

/* d, row-major, general matrix pseudo-inverse (the svd way): double precision */
void utility_dpinv(const double* inM,          /* in matrix; flat:[dim1][dim2] */
                   const int dim1,
                   const int dim2,
                   double* outM);              /* out matrix; flat:[dim2][dim1] */

/*-------------------------------- matrix inversion (?inv) ----------------------------------*/
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


#endif /* SAF_VECLIB_H_INCLUDED */
