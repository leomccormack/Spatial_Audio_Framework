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
 *     saf_veclib.c
 * Description:
 *     Contains a collection of useful memory allocation functions and cross-platform
 *     complex number wrappers. Optimised linear algebra routines utilising BLAS and LAPACK
 *     are also included.
 * Dependencies:
 *     Windows users only: custom Intel MKL '.lib/.dll' files are required.
 *     Mac users only: saf_utilities will utilise Apple's Accelerate library by default.
 *     However, Mac users may elect to use a custom Intel MKL '.dylib' instead.
 *     Further instructions for both Windows/Mac users can be found here:
 *     https://github.com/leomccormack/Spatial_Audio_Framework
 * Author, date created:
 *     Leo McCormack, 11.07.2016
 */

#include "saf_veclib.h"
#include "saf_complex.h"
#include "saf_sort.h"
#include <float.h>

#ifndef MIN
  #define MIN(a,b) (( (a) < (b) ) ? (a) : (b))
#endif
#ifndef MAX
  #define MAX(a,b) (( (a) > (b) ) ? (a) : (b))
#endif

/*----------------------------- index of min-abs-value (?iminv) -----------------------------*/

void utility_siminv(const float* a, const int len, int* index)
{
#if defined(__ACCELERATE__)
    float minVal;
    vDSP_Length ind_tmp;
    vDSP_minvi(a, 1, &minVal, &ind_tmp, len);
    *index = (int)ind_tmp;
#elif defined(INTEL_MKL_VERSION)
    *index = (int)cblas_isamin(len, a, 1);
#else
    UNTESTED
    int j;
    float minVal;
    minVal = 3.402823e+38f;
    for(j=0; j<len; j++){
        if(fabsf(a[j])<minVal){
            minVal = fabsf(a[j]);
            (*index) = j;
        }
    }
#endif
}

void utility_ciminv(const float_complex* a, const int len, int* index)
{
#if defined(__ACCELERATE__)
    int i;
    float* abs_a;
    float minVal;
    abs_a = malloc(len*sizeof(float));
    for(i=0; i<len; i++)
        abs_a[i] = cabsf(a[i]);
    vDSP_Length ind_tmp;
    vDSP_minvi(abs_a, 1, &minVal, &ind_tmp, len);
    *index = (int)ind_tmp;
    free(abs_a);
#elif defined(INTEL_MKL_VERSION)
    *index = (int)cblas_icamin(len, a, 1);
#else
    UNTESTED
    int j;
    float minVal;
    minVal = 3.402823e+38f;
    for(j=0; j<len; j++){
        if(cabsf(a[j])<maxVal){
            minVal = cabsf(a[j]);
            (*index) = j;
        }
    }
#endif
}

/*----------------------------- index of max-abs-value (?imaxv) -----------------------------*/

void utility_simaxv(const float* a, const int len, int* index)
{
#if defined(__ACCELERATE__) || defined(INTEL_MKL_VERSION)
    *index = (int)cblas_isamax(len, a, 1);
#else
    UNTESTED
    int j;
    float maxVal;
    maxVal = 1.175494e-38f;
    for(j=0; j<len; j++){
        if(fabsf(a[j])>maxVal){
            maxVal = fabsf(a[j]);
            (*index) = j;
        }
    }
#endif
}

void utility_cimaxv(const float_complex* a, const int len, int* index)
{
#if defined(__ACCELERATE__) || defined(INTEL_MKL_VERSION)
    *index = (int)cblas_icamax(len, a, 1);
#else
    UNTESTED
    int j;
    float maxVal;
    maxVal = 1.175494e-38f;
    for(j=0; j<len; j++){
        if(cabsf(a[j])>maxVal){
            maxVal = cabsf(a[j]);
            (*index) = j;
        }
    }
#endif
}

/*----------------------------------- vector-abs (?vabs) ------------------------------------*/

void utility_svabs(const float* a, const int len, float* c)
{
#ifdef INTEL_MKL_VERSION
    vsAbs(len, a, c);
#else
    int i;
    for(i=0; i<len; i++)
        c[i] = fabsf(a[i]);
#endif
}

void utility_cvabs(const float_complex* a, const int len, float* c)
{
#ifdef INTEL_MKL_VERSION
    vcAbs(len, (MKL_Complex8*)a, c);
#else
    int i;
    for(i=0; i<len; i++)
        c[i] = cabsf(a[i]);
#endif
}

/*------------------------------ vector-vector copy (?vvcopy) -------------------------------*/

void utility_svvcopy(const float* a, const int len, float* c)
{
#if defined(__ACCELERATE__) || defined(INTEL_MKL_VERSION)
    cblas_scopy(len, a, 1, c, 1);
#else
    memcpy(c, a, len*sizeof(float));
#endif
}

void utility_cvvcopy(const float_complex* a, const int len, float_complex* c)
{
#if defined(__ACCELERATE__) || defined(INTEL_MKL_VERSION)
    cblas_ccopy(len, a, 1, c, 1);
#else
    memcpy(c, a, len*sizeof(float_complex));
#endif
}

/*---------------------------- vector-vector addition (?vvadd) ------------------------------*/

void utility_svvadd(float* a, const float* b, const int len, float* c)
{
#ifdef __ACCELERATE__
	if (c == NULL) {
		float* tmp;
		tmp = malloc(len * sizeof(float));
		vDSP_vadd(a, 1, b, 1, tmp, 1, len);
		utility_svvcopy(tmp, len, a);
		free(tmp);
	}
	else
		vDSP_vadd(a, 1, b, 1, c, 1, len);
#elif INTEL_MKL_VERSION
	if (c == NULL) {
		float* tmp;
		tmp = malloc(len * sizeof(float));
		vsAdd(len, a, b, tmp);
		utility_svvcopy(tmp, len, a);
		free(tmp);
	}
	else
		vsAdd(len, a, b, c);
#else
	int i;
	for (i = 0; i < len; i++)
		c[i] = a[i] + b[i];
#endif
}

void utility_cvvadd(float_complex* a, const float_complex* b, const int len, float_complex* c)
{
#ifdef INTEL_MKL_VERSION
	if (c == NULL) {
		float_complex* tmp;
		tmp = malloc(len * sizeof(float_complex));
		vcAdd(len, (MKL_Complex8*)a, (MKL_Complex8*)b, (MKL_Complex8*)tmp);
		utility_cvvcopy(tmp, len, a);
		free(tmp);
	}
	else
		vcAdd(len, (MKL_Complex8*)a, (MKL_Complex8*)b, (MKL_Complex8*)c);
#else
	int i;
	for (i = 0; i < len; i++)
		c[i] = ccaddf(a[i], b[i]);
#endif
}

/*--------------------------- vector-vector subtraction (?vvsub) ----------------------------*/

void utility_svvsub(float* a, const float* b, const int len, float* c)
{
#ifdef __ACCELERATE__
	if (c == NULL) {
		float* tmp;
		tmp = malloc(len * sizeof(float));
		vDSP_vsub(a, 1, b, 1, tmp, 1, len);
		utility_svvcopy(tmp, len, a);
		free(tmp);
	}
	else
		vDSP_vsub(a, 1, b, 1, c, 1, len);
#elif INTEL_MKL_VERSION
	if (c == NULL) {
		float* tmp;
		tmp = malloc(len * sizeof(float));
		vsSub(len, a, b, tmp);
		utility_svvcopy(tmp, len, a);
		free(tmp);
	}
	else
		vsSub(len, a, b, c);
#else
	int i;
	for (i = 0; i < len; i++)
		c[i] = a[i] + b[i];
#endif
}

void utility_cvvsub(float_complex* a, const float_complex* b, const int len, float_complex* c)
{
#ifdef INTEL_MKL_VERSION
	if (c == NULL) {
		float_complex* tmp;
		tmp = malloc(len * sizeof(float_complex));
		vcSub(len, (MKL_Complex8*)a, (MKL_Complex8*)b, (MKL_Complex8*)tmp);
		utility_cvvcopy(tmp, len, a);
		free(tmp);
	}
	else
		vcSub(len, (MKL_Complex8*)a, (MKL_Complex8*)b, (MKL_Complex8*)c);
#else
	int i;
    if (c == NULL) {
        for (i = 0; i < len; i++)
            a[i] = ccaddf(a[i], b[i]);
    }
    else{
        for (i = 0; i < len; i++)
            c[i] = ccaddf(a[i], b[i]);
    }
#endif
}

/*------------------------- vector-vector multiplication (?vvmul) ---------------------------*/

void utility_svvmul(float* a, const float* b, const int len, float* c)
{
#ifdef __ACCELERATE__
    if(c==NULL){
        float* tmp;
        tmp=malloc(len*sizeof(float));
        vDSP_vmul(a, 1, b, 1, tmp, 1, len);
        utility_svvcopy(tmp, len, a);
        free(tmp);
    }
    else
        vDSP_vmul(a, 1, b, 1, c, 1, len);
#elif INTEL_MKL_VERSION
	if (c == NULL) {
		float* tmp;
		tmp = malloc(len * sizeof(float));
        vsMul(len, a, b, tmp);
		utility_svvcopy(tmp, len, a);
		free(tmp);
	}
	else
		vsMul(len, a, b, c);
#else
    int i;
    for (i = 0; i < len; i++)
        c[i] = a[i] * b[i];
#endif
}

void utility_cvvmul(float_complex* a, const float_complex* b, const int len, float_complex* c)
{
#ifdef INTEL_MKL_VERSION
	if (c == NULL) {
		float_complex* tmp;
		tmp = malloc(len * sizeof(float_complex));
		vcMul(len, (MKL_Complex8*)a, (MKL_Complex8*)b, (MKL_Complex8*)tmp);
		utility_cvvcopy(tmp, len, a);
		free(tmp);
	}
	else
		vcMul(len, (MKL_Complex8*)a, (MKL_Complex8*)b, (MKL_Complex8*)c);
#else
	int i;
    if (c == NULL) {
        for (i = 0; i < len; i++)
            a[i] = ccmulf(a[i], b[i]);
    }
    else{
        for (i = 0; i < len; i++)
            c[i] = ccmulf(a[i], b[i]);
    }
#endif
}

/*--------------------------- vector-vector dot product (?vvdot) ----------------------------*/

void utility_svvdot(const float* a, const float* b, const int len, float* c)
{
#if defined(__ACCELERATE__) || defined(INTEL_MKL_VERSION)
    c[0] = cblas_sdot (len, a, 1, b, 1);
#else
    int i;
    c[0] = 0.0f;
    for(i=0; i<len; i++)
        c[0] += a[i] * b[i];
#endif
}

void utility_cvvdot(const float_complex* a, const float_complex* b, const int len, CONJ_FLAG flag, float_complex* c)
{
#if defined(__ACCELERATE__) || defined(INTEL_MKL_VERSION)
    switch(flag){
        default:
        case NO_CONJ:
            cblas_cdotu_sub(len, a, 1, b, 1, c);
            break;
        case CONJ:
            cblas_cdotc_sub(len, a, 1, b, 1, c);
            break;
    }
#else
    int i;
    c[0] = cmplxf(0.0f,0.0f);
    switch(flag){
        default:
        case NO_CONJ:
            for(i=0; i<len; i++)
                c[0] = ccaddf(c[0], ccmulf(a[i], b[i]));
            break;
        case CONJ:
            for(i=0; i<len; i++)
                c[0] = ccaddf(c[0], ccmulf(conjf(a[i]), b[i]));
            break;
    }
#endif
}

/*------------------------------ vector-scalar product (?vsmul) -----------------------------*/

void utility_svsmul(float* a, const float* s, const int len, float* c)
{
#ifdef __ACCELERATE__
    if(c==NULL)
        cblas_sscal(len, s[0], a, 1);
    else
        vDSP_vsmul(a, 1, s, c, 1, len);
#elif INTEL_MKL_VERSION
    if (c == NULL)
        cblas_sscal(len, s[0], a, 1);
    else {
        memcpy(c, a, len*sizeof(float));
        cblas_sscal(len, s[0], c, 1);
    }
#else
    int i;
    for(i=0; i<len; i++)
        c[i] = a[i] * s[0];
#endif
}

void utility_cvsmul(float_complex* a, const float_complex* s, const int len, float_complex* c)
{
#if defined(__ACCELERATE__) || defined(INTEL_MKL_VERSION)
    if (c == NULL)
        cblas_cscal(len, s, a, 1);
    else {
        int i;
        for (i = 0; i<len; i++)
            c[i] = ccmulf(a[i], s[0]);
    }
#else
    int i;
    for(i=0; i<len; i++)
        c[i] = ccmulf(a[i], s[0]);
#endif
}

/*----------------------------- vector-scalar division (?vsdiv) -----------------------------*/

void utility_svsdiv(float* a, const float* s, const int len, float* c)
{
    if(s[0] == 0.0f){
        memset(c, 0, len*sizeof(float));
        return;
    }
#ifdef __ACCELERATE__
    vDSP_vsdiv(a, 1, s, c, 1, len);
#else
    int i;
    for(i=0; i<len; i++)
        c[i] = a[i] / s[0];
#endif
}

/*----------------------------- vector-scalar addition (?vsadd) -----------------------------*/

void utility_svsadd(float* a, const float* s, const int len, float* c)
{
#ifdef __ACCELERATE__
    vDSP_vsadd(a, 1, s, c, 1, len);
#else
    int i;
    for(i=0; i<len; i++)
        c[i] = a[i] + s[0];
#endif
}

/*---------------------------- vector-scalar subtraction (?vssub) ---------------------------*/

void utility_svssub(float* a, const float* s, const int len, float* c)
{
    int i;
    for(i=0; i<len; i++)
        c[i] = a[i] - s[0];
}

/*---------------------------- singular-value decomposition (?svd) --------------------------*/

void utility_ssvd(const float* A, const int dim1, const int dim2, float* U, float* S, float* V, float* sing)
{
    int i, j, m, n, lda, ldu, ldvt, info, lwork;
    m = dim1; n = dim2; lda = dim1; ldu = dim1; ldvt = dim2;
    float wkopt;
    float* a, *s, *u, *vt, *work;

    a = malloc(lda*n*sizeof(float));
    s = malloc(MIN(n,m)*sizeof(float));
    u = malloc(ldu*m*sizeof(float));
    vt = malloc(ldvt*n*sizeof(float));
    /* store in column major order */
    for(i=0; i<dim1; i++)
        for(j=0; j<dim2; j++)
            a[j*dim1+i] = A[i*dim2 +j];
    lwork = -1;
    sgesvd_( "A", "A", &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, &wkopt, &lwork,
           &info );
    lwork = (int)wkopt;
    work = (float*)malloc( lwork*sizeof(float) );
    sgesvd_( "A", "A", &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork,
           &info );
    if( info > 0 ) {
        /* svd failed to converge */
    }
    else {
        /* svd successful */
        if (S != NULL){
            memset(S, 0, dim1*dim2*sizeof(float));
            /* singular values on the diagonal MIN(dim1, dim2). The remaining elements are 0.  */
            for(i=0; i<MIN(dim1, dim2); i++)
                S[i*dim2+i] = s[i];
        }
        /*return as row-major*/
        if (U != NULL)
            for(i=0; i<dim1; i++)
                for(j=0; j<dim1; j++)
                    U[i*dim1+j] = u[j*dim1+i];
        
        /* lapack returns VT, i.e. row-major V already */
        if (V != NULL)
            for(i=0; i<dim2; i++)
                for(j=0; j<dim2; j++)
                    V[i*dim2+j] = vt[i*dim2+j];
        
        if (sing != NULL)
            for(i=0; i<MIN(dim1, dim2); i++)
                sing[i] = s[i];
    }
    
    free( (void*)a );
    free( (void*)s );
    free( (void*)u );
    free( (void*)vt );
    free( (void*)work );
}

void utility_csvd(const float_complex* A, const int dim1, const int dim2, float_complex* U, float_complex* S, float_complex* V, float* sing)
{
    int i, j, m, n, lda, ldu, ldvt, info, lwork;
    m = dim1; n = dim2; lda = dim1; ldu = dim1; ldvt = dim2;
    float_complex wkopt;
    float_complex* a, *u, *vt, *work;
    float* s, *rwork;
    
    a = malloc(lda*n*sizeof(float_complex));
    s = malloc(MIN(n,m)*sizeof(float));
    u = malloc(ldu*m*sizeof(float_complex));
    vt = malloc(ldvt*n*sizeof(float_complex));
    rwork = malloc(m*MAX(1, 5*MIN(n,m))*sizeof(float));
    /* store in column major order */
    for(i=0; i<dim1; i++)
        for(j=0; j<dim2; j++)
            a[j*dim1+i] = A[i*dim2 +j];
    lwork = -1;
#if defined(__APPLE__) && !defined(SAF_USE_INTEL_MKL)
    cgesvd_( "A", "A", (__CLPK_integer*)&m, (__CLPK_integer*)&n, (__CLPK_complex*)a, (__CLPK_integer*)&lda, s,
            (__CLPK_complex*)u, (__CLPK_integer*)&ldu, (__CLPK_complex*)vt, &ldvt, (__CLPK_complex*)&wkopt, &lwork, rwork, (__CLPK_integer*)&info );
#elif INTEL_MKL_VERSION
    cgesvd_( "A", "A", &m, &n, (MKL_Complex8*)a, &lda, s, (MKL_Complex8*)u, &ldu, (MKL_Complex8*)vt, &ldvt,
            (MKL_Complex8*)&wkopt, &lwork, rwork, &info );
#endif
    lwork = (int)(crealf(wkopt)+0.01f);
    work = malloc( lwork*sizeof(float_complex) );
#if defined(__APPLE__) && !defined(SAF_USE_INTEL_MKL)
    cgesvd_( "A", "A", &m, &n, (__CLPK_complex*)a, &lda, s, (__CLPK_complex*)u, &ldu, (__CLPK_complex*)vt, &ldvt,
            (__CLPK_complex*)work, &lwork, rwork, &info);
#elif INTEL_MKL_VERSION
    cgesvd_( "A", "A", &m, &n, (MKL_Complex8*)a, &lda, s, (MKL_Complex8*)u, &ldu, (MKL_Complex8*)vt, &ldvt,
            (MKL_Complex8*)work, &lwork, rwork, &info);
#endif
    if( info > 0 ) {
        /* svd failed to converge */
    }
    else {
        /* svd successful */
        if (S != NULL){
            memset(S, 0, dim1*dim2*sizeof(float_complex));
            /* singular values on the diagonal MIN(dim1, dim2). The remaining elements are 0.  */
            for(i=0; i<MIN(dim1, dim2); i++)
                S[i*dim2+i] = cmplxf(s[i], 0.0f);
        }
        /*return as row-major*/
        if (U != NULL)
            for(i=0; i<dim1; i++)
                for(j=0; j<dim1; j++)
                    U[i*dim1+j] = u[j*dim1+i];
        
        /* lapack returns VT, i.e. row-major V already */
        if (V != NULL)
            for(i=0; i<dim2; i++)
                for(j=0; j<dim2; j++)
                    V[i*dim2+j] = conjf(vt[i*dim2+j]); /* v^H */
        
        if (sing != NULL)
            for(i=0; i<MIN(dim1, dim2); i++)
                sing[i] = s[i];
    }
    
    free( (void*)a );
    free( (void*)s );
    free( (void*)u );
    free( (void*)vt );
    free( (void*)work );
    free(rwork);
}

/*------------------------ symmetric eigenvalue decomposition (?seig) -----------------------*/

void utility_sseig(const float* A, const int dim, int sortDecFLAG, float* V, float* D, float* eig)
{
    int i, j, n, lda, info, lwork;
    float wkopt;
    float* work, *w, *a;
    
    n = dim;
    lda = dim;
    w = malloc(dim*sizeof(float));
    a = malloc(dim*dim*sizeof(float));
    
    /* store in column major order (i.e. transpose) */
    for(i=0; i<dim; i++)
        for(j=0; j<dim; j++)
            a[i*dim+j] = A[j*dim+i];
    
    /* solve the eigenproblem */
    lwork = -1;
    ssyev_( "Vectors", "Upper", &n, a, &lda, w, &wkopt, &lwork, &info );
    lwork = (int)wkopt;
    work = (float*)malloc( lwork*sizeof(float) );
    ssyev_( "Vectors", "Upper", &n, a, &lda, w, work, &lwork, &info );

    /* output */
    if(D!=NULL)
        memset(D, 0, dim*dim*sizeof(float));
    if( info > 0 ) {
        /* failed to converge and find the eigenvalues */
        if(V!=NULL)
            memset(V, 0, dim*dim*sizeof(float));
    }
    else{
        if(sortDecFLAG){
            for(i=0; i<dim; i++){
                if(V!=NULL)
                    for(j=0; j<dim; j++)
                        V[i*dim+j] = a[(dim-j-1)*dim+i]; /* traspose, back to row-major and reverse order */
                if(D!=NULL)
                    D[i*dim+i] = w[dim-i-1]; /* store along the diagonal, reversing the order */
                if(eig!=NULL)
                    eig[i] = w[dim-i-1];
            }
        }
        else{
            for(i=0; i<dim; i++){
                if(V!=NULL)
                    for(j=0; j<dim; j++)
                        V[i*dim+j] = a[j*dim+i]; /* traspose, back to row-major */
                if(D!=NULL)
                    D[i*dim+i] = w[i]; /* store along the diagonal */
                if(eig!=NULL)
                    eig[i] = w[i];
            }
        }
    }
    
    free(w);
    free(a);
    free(work);
}

void utility_cseig(const float_complex* A, const int dim, int sortDecFLAG, float_complex* V, float_complex* D, float* eig)
{
    int i, j, n, lda, info, lwork;
    float_complex wkopt;
    float *w, *rwork;
    float_complex* a, *work;
    
    n = dim;
    lda = dim;
    w = malloc(dim*sizeof(float));
    a = malloc(dim*dim*sizeof(float_complex));
    rwork = malloc((3*n-2)*sizeof(float));
    
    /* store in column major order (i.e. transpose) */
    for(i=0; i<dim; i++)
        for(j=0; j<dim; j++)
            a[i*dim+j] = A[j*dim+i];
    
    /* solve the eigenproblem */
    lwork = -1;
#if defined(__APPLE__) && !defined(SAF_USE_INTEL_MKL)
    cheev_( "Vectors", "Upper", (__CLPK_integer*)&n, (__CLPK_complex*)a, (__CLPK_integer*)&lda,
           (__CLPK_real*)w, (__CLPK_complex*)&wkopt, (__CLPK_integer*)&lwork, (__CLPK_real*)rwork, (__CLPK_integer*)&info );
#elif INTEL_MKL_VERSION
    cheev_( "Vectors", "Upper", &n, (MKL_Complex8*)a, &lda, w, (MKL_Complex8*)&wkopt, &lwork, rwork, &info );
#endif
    lwork = (int)crealf(wkopt);
    work = (float_complex*)malloc( lwork*sizeof(float_complex) );
#if defined(__APPLE__) && !defined(SAF_USE_INTEL_MKL)
    cheev_( "Vectors", "Upper", (__CLPK_integer*)&n, (__CLPK_complex*)a, (__CLPK_integer*)&lda,
           (__CLPK_real*)w, (__CLPK_complex*)work, (__CLPK_integer*)&lwork, (__CLPK_real*)rwork, (__CLPK_integer*)&info );
#elif INTEL_MKL_VERSION
    cheev_( "Vectors", "Upper", &n, (MKL_Complex8*)a, &lda, w, (MKL_Complex8*)work, &lwork, rwork, &info );
#endif
    
    /* output */
    if(D!=NULL)
        memset(D, 0, dim*dim*sizeof(float_complex));
    if( info > 0 ) {
        /* failed to converge and find the eigenvalues */
        if(V!=NULL)
            memset(V, 0, dim*dim*sizeof(float_complex));
    }
    else{
        if(sortDecFLAG){
            for(i=0; i<dim; i++){
                if(V!=NULL)
                    for(j=0; j<dim; j++)
                        V[i*dim+j] = a[(dim-j-1)*dim+i]; /* traspose, back to row-major and reverse order */
                if(D!=NULL)
                    D[i*dim+i] = cmplxf(w[dim-i-1], 0.0f); /* store along the diagonal, reversing the order */
                if(eig!=NULL)
                    eig[i] = w[dim-i-1];
            }
        }
        else{
            for(i=0; i<dim; i++){
                if(V!=NULL)
                    for(j=0; j<dim; j++)
                        V[i*dim+j] = a[j*dim+i]; /* traspose, back to row-major */
                if(D!=NULL)
                    D[i*dim+i] = cmplxf(w[i], 0.0f); /* store along the diagonal */
                if(eig!=NULL)
                    eig[i] = w[i];
            }
        }
    }
    
    free(w);
    free(a);
    free(work);
}

/*-----------------------------  eigenvalue decomposition (?eig) ----------------------------*/

void utility_ceig(const float_complex* A, const int dim, int sortDecFLAG, float_complex* VL, float_complex* VR, float_complex* D, float* eig)
{
    int i, j, n, lda, ldvl, ldvr, info, lwork;
    float_complex wkopt;
    float_complex* work, *w, *vl, *vr, *a;
    float* rwork;
    
    n = lda = ldvl = ldvr = dim;
    rwork = malloc(2*dim*sizeof(float));
    w = malloc(dim*sizeof(float_complex));
    vl = malloc(dim*dim*sizeof(float_complex));
    vr = malloc(dim*dim*sizeof(float_complex));
    a = malloc(dim*dim*sizeof(float_complex));
    
    /* store in column major order (i.e. transpose) */
    for(i=0; i<dim; i++)
        for(j=0; j<dim; j++)
            a[i*dim+j] = A[j*dim+i];
    
    /* solve the eigenproblem */
    lwork = -1;
#if defined(__APPLE__) && !defined(SAF_USE_INTEL_MKL)
    cgeev_( "Vectors", "Vectors", &n, (__CLPK_complex*)a, &lda, (__CLPK_complex*)w, (__CLPK_complex*)vl,
           &ldvl, (__CLPK_complex*)vr, &ldvr, (__CLPK_complex*)&wkopt, &lwork, rwork, &info );
#elif INTEL_MKL_VERSION
    cgeev_( "Vectors", "Vectors", &n, (MKL_Complex8*)a, &lda, (MKL_Complex8*)w, (MKL_Complex8*)vl, &ldvl, (MKL_Complex8*)vr, &ldvr, (MKL_Complex8*)&wkopt, &lwork, rwork, &info );
#endif
    lwork = (int)crealf(wkopt);
    work = (float_complex*)malloc( lwork*sizeof(float_complex) );
#if defined(__APPLE__) && !defined(SAF_USE_INTEL_MKL)
    cgeev_( "Vectors", "Vectors", &n, (__CLPK_complex*)a, &lda, (__CLPK_complex*)w, (__CLPK_complex*)vl,
           &ldvl, (__CLPK_complex*)vr, &ldvr, (__CLPK_complex*)work, &lwork, rwork, &info );
#elif INTEL_MKL_VERSION
    cgeev_( "Vectors", "Vectors", &n, (MKL_Complex8*)a, &lda, (MKL_Complex8*)w, (MKL_Complex8*)vl, &ldvl, (MKL_Complex8*)vr, &ldvr, (MKL_Complex8*)work, &lwork, rwork, &info );
#endif
    
    /* sort the eigenvalues */
    float* wr;
    int* sort_idx;
    wr=malloc(dim*sizeof(float));
    sort_idx=malloc(dim*sizeof(int));
    for(i=0; i<dim; i++)
        wr[i] = crealf(w[i]);
    sortf(wr, NULL, sort_idx, dim, sortDecFLAG);
    
    /* output */
    if(D!=NULL)
        memset(D, 0, dim*dim*sizeof(float_complex));
    if( info > 0 ) {
        /* failed to converge and find the eigenvalues */
        if(VL!=NULL)
            memset(VL, 0, dim*dim*sizeof(float_complex));
        if(VR!=NULL)
            memset(VR, 0, dim*dim*sizeof(float_complex));
        if(eig!=NULL)
            memset(eig, 0, dim*sizeof(float));
    }
    else{
        for(i=0; i<dim; i++){
            if(VL!=NULL)
                for(j=0; j<dim; j++)
                    VL[i*dim+j] = vl[sort_idx[j]*dim+i]; /* transpose, back to row-major */
            if(VR!=NULL)
                for(j=0; j<dim; j++)
                    VR[i*dim+j] = vr[sort_idx[j]*dim+i]; /* transpose, back to row-major */
            if(D!=NULL)
                D[i*dim+i] = cmplxf(wr[i], 0.0f); /* store along the diagonal */
            if(eig!=NULL)
                eig[i] = wr[i];
        }
    }
    
    free(rwork);
    free(work);
    free(w);
    free(vl);
    free(vr);
    free(a);
    free(wr);
    free(sort_idx);
}

/*------------------------------ general linear solver (?glslv) -----------------------------*/

void utility_sglslv(const float* A, const int dim, float* B, int nCol, float* X)
{
    int i, j, n = dim, nrhs = nCol, lda = dim, ldb = dim, info;
    int* IPIV;
    IPIV = malloc(dim*sizeof(int));
    float* a, *b;
    a = malloc(dim*dim*sizeof(float));
    b = malloc(dim*nrhs*sizeof(float));
    
    /* store in column major order */
    for(i=0; i<dim; i++)
        for(j=0; j<dim; j++)
            a[j*dim+i] = A[i*dim+j];
    for(i=0; i<dim; i++)
        for(j=0; j<nCol; j++)
            b[j*dim+i] = B[i*nCol+j];
    
    /* solve Ax = b for each column in b (b is replaced by the solution: x) */
    sgesv_( &n, &nrhs, a, &lda, IPIV, b, &ldb, &info );
    
    if(info>0){
        /* A is singular, solution not possible */
        memset(X, 0, dim*nCol*sizeof(float));
    }
    else{
        /* store solution in row-major order */
        for(i=0; i<dim; i++)
            for(j=0; j<nCol; j++)
                X[i*nCol+j] = b[j*dim+i];
    }
    
    free(IPIV);
    free(a);
    free(b);
}

void utility_cglslv(const float_complex* A, const int dim, float_complex* B, int nCol, float_complex* X)
{
    int i, j, n = dim, nrhs = nCol, lda = dim, ldb = dim, info;
    int* IPIV;
    IPIV = malloc(dim*sizeof(int));
    float_complex* a, *b;
    a = malloc(dim*dim*sizeof(float_complex));
    b = malloc(dim*nrhs*sizeof(float_complex));
    
    /* store in column major order */
    for(i=0; i<dim; i++)
        for(j=0; j<dim; j++)
            a[j*dim+i] = A[i*dim+j];
    for(i=0; i<dim; i++)
        for(j=0; j<nCol; j++)
            b[j*dim+i] = B[i*nCol+j];
    
    /* solve Ax = b for each column in b (b is replaced by the solution: x) */
#if defined(__APPLE__) && !defined(SAF_USE_INTEL_MKL)
    cgesv_( &n, &nrhs, (__CLPK_complex*)a, &lda, IPIV, (__CLPK_complex*)b, &ldb, &info );
#elif INTEL_MKL_VERSION
    cgesv_( &n, &nrhs, (MKL_Complex8*)a, &lda, IPIV, (MKL_Complex8*)b, &ldb, &info );
#endif
    
    if(info>0){
        /* A is singular, solution not possible */
        memset(X, 0, dim*nCol*sizeof(float_complex));
    }
    else{
        /* store solution in row-major order */
        for(i=0; i<dim; i++)
            for(j=0; j<nCol; j++)
                X[i*nCol+j] = b[j*dim+i];
    }
    
    free(IPIV);
    free(a);
    free(b);
}

/*----------------------------- symmetric linear solver (?slslv) ----------------------------*/

void utility_sslslv(const float* A, const int dim, float* B, int nCol, float* X)
{
    int i, j, n = dim, nrhs = nCol, lda = dim, ldb = dim, info;
    float* a, *b;
    a = malloc(dim*dim*sizeof(float));
    b = malloc(dim*nrhs*sizeof(float));
    
    /* store in column major order */
    for(i=0; i<dim; i++)
        for(j=0; j<dim; j++)
            a[j*dim+i] = A[i*dim+j];
    for(i=0; i<dim; i++)
        for(j=0; j<nCol; j++)
            b[j*dim+i] = B[i*nCol+j];
    
    /* solve Ax = b for each column in b (b is replaced by the solution: x) */
    sposv_( "U", &n, &nrhs, a, &lda, b, &ldb, &info );
    
    if(info>0){
        /* A is not positive definate, solution not possible */
        memset(X, 0, dim*nCol*sizeof(float));
    }
    else{
        /* store solution in row-major order */
        for(i=0; i<dim; i++)
            for(j=0; j<nCol; j++)
                X[i*nCol+j] = b[j*dim+i];
    }
    
    free(a);
    free(b);
}

void utility_cslslv(const float_complex* A, const int dim, float_complex* B, int nCol, float_complex* X)
{
    int i, j, n = dim, nrhs = nCol, lda = dim, ldb = dim, info;
    float_complex* a, *b;
    a = malloc(dim*dim*sizeof(float_complex));
    b = malloc(dim*nrhs*sizeof(float_complex));
    
    /* store in column major order */
    for(i=0; i<dim; i++)
        for(j=0; j<dim; j++)
            a[j*dim+i] = A[i*dim+j];
    for(i=0; i<dim; i++)
        for(j=0; j<nCol; j++)
            b[j*dim+i] = B[i*nCol+j];
    
    /* solve Ax = b for each column in b (b is replaced by the solution: x) */
#if defined(__APPLE__) && !defined(SAF_USE_INTEL_MKL)
    cposv_( "U", &n, &nrhs, (__CLPK_complex*)a, &lda, (__CLPK_complex*)b, &ldb, &info );
#elif INTEL_MKL_VERSION
    cposv_( "U", &n, &nrhs, (MKL_Complex8*)a, &lda, (MKL_Complex8*)b, &ldb, &info );
#endif
    
    if(info>0){
        /* A is not positive definate, solution not possible */
        memset(X, 0, dim*nCol*sizeof(float_complex));
    }
    else{
        /* store solution in row-major order */
        for(i=0; i<dim; i++)
            for(j=0; j<nCol; j++)
                X[i*nCol+j] = b[j*dim+i];
    }
    
    free(a);
    free(b);
}



/*----------------------------- matrix pseudo-inverse (?pinv) -------------------------------*/

void utility_spinv(const float* inM, const int dim1, const int dim2, float* outM)
{
    
    int i, j, m, n, k, lda, ldu, ldvt, lwork, info;
    float* a, *s, *u, *vt, *inva, *work;
    float ss, wkopt;
    
    m = lda = ldu = dim1;
    n = dim2;
    k = ldvt = m < n ? m : n; 
    a = malloc(m*n*sizeof(float) );
    for(i=0; i<m; i++)
        for(j=0; j<n; j++)
            a[j*m+i] = inM[i*n+j]; /* store in column major order */
    s = (float*)malloc(k*sizeof(float));
    u = (float*)malloc(ldu*k*sizeof(float));
    vt = (float*)malloc(ldvt*n*sizeof(float));
    lwork = -1;
#if defined(__APPLE__) && !defined(SAF_USE_INTEL_MKL)
    sgesvd_("S", "S", (__CLPK_integer*)&m, (__CLPK_integer*)&n, a, (__CLPK_integer*)&lda,
            s, u, (__CLPK_integer*)&ldu, vt, (__CLPK_integer*)&ldvt, &wkopt, (__CLPK_integer*)&lwork, (__CLPK_integer*)&info);
#elif INTEL_MKL_VERSION
    sgesvd_("S", "S", &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, &wkopt, &lwork, &info);
#endif
    lwork = (int)wkopt;
    work = (float*)malloc(lwork*sizeof(float));
    
#if defined(__APPLE__) && !defined(SAF_USE_INTEL_MKL)
    sgesvd_("S", "S", (__CLPK_integer*)&m, (__CLPK_integer*)&n, a, (__CLPK_integer*)&lda,
            s, u, (__CLPK_integer*)&ldu, vt, (__CLPK_integer*)&ldvt, work, (__CLPK_integer*) &lwork, (__CLPK_integer*)&info); /* Compute SVD */
#elif INTEL_MKL_VERSION
    sgesvd_("S", "S", &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info ); /* Compute SVD */
#endif
    if( info > 0 ) {
        memset(outM, 0, dim1*dim2*sizeof(float));
        free((void*)a);
        free((void*)s);
        free((void*)u);
        free((void*)vt);
        free((void*)work);
        return; /*failed to converge, output 0s */
    }
    int incx=1;
    for(i=0; i<k; i++){
        if(s[i] > 1.0e-5f)
            ss=1.0f/s[i];
        else
            ss=s[i];
        cblas_sscal(m, ss, &u[i*m], incx);     
    }
    inva = (float*)malloc(n*m*sizeof(float));
    int ld_inva=n;
    cblas_sgemm( CblasColMajor, CblasTrans, CblasTrans, n, m, k, 1.0f,
                vt, ldvt,
                u, ldu, 0.0f,
                inva, ld_inva);
    for(i=0; i<m; i++)
        for(j=0; j<n; j++)
            outM[j*m+i] = inva[i*n+j]; /* return in row-major order */

    /* clean-up */
    free((void*)a);
    free((void*)s);
    free((void*)u);
    free((void*)vt);
    free((void*)inva);
    free((void*)work);
}

void utility_dpinv(const double* inM, const int dim1, const int dim2, double* outM)
{
    
    int i, j, m, n, k, lda, ldu, ldvt, lwork, info;
    double* a, *s, *u, *vt, *inva, *work;
    double ss, wkopt;
    
    m = lda = ldu = dim1;
    n = dim2;
    k = ldvt = m < n ? m : n;
    a = malloc(m*n*sizeof(double) );
    for(i=0; i<m; i++)
        for(j=0; j<n; j++)
            a[j*m+i] = inM[i*n+j]; /* store in column major order */
    s = (double*)malloc(k*sizeof(double));
    u = (double*)malloc(ldu*k*sizeof(double));
    vt = (double*)malloc(ldvt*n*sizeof(double));
    lwork = -1;
#if defined(__APPLE__) && !defined(SAF_USE_INTEL_MKL)
    dgesvd_("S", "S", (__CLPK_integer*)&m, (__CLPK_integer*)&n, a, (__CLPK_integer*)&lda,
            s, u, (__CLPK_integer*)&ldu, vt, (__CLPK_integer*)&ldvt, &wkopt, (__CLPK_integer*)&lwork, (__CLPK_integer*)&info);
#elif INTEL_MKL_VERSION
    dgesvd_("S", "S", &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, &wkopt, &lwork, &info);
#endif
    lwork = (int)wkopt;
    work = (double*)malloc(lwork*sizeof(double));
    
#if defined(__APPLE__) && !defined(SAF_USE_INTEL_MKL)
    dgesvd_("S", "S", (__CLPK_integer*)&m, (__CLPK_integer*)&n, a, (__CLPK_integer*)&lda,
            s, u, (__CLPK_integer*)&ldu, vt, (__CLPK_integer*)&ldvt, work, (__CLPK_integer*) &lwork, (__CLPK_integer*)&info); /* Compute SVD */
#elif INTEL_MKL_VERSION
    dgesvd_("S", "S", &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info ); /* Compute SVD */
#endif
    if( info > 0 ) {
        memset(outM, 0, dim1*dim2*sizeof(float));
        free((void*)a);
        free((void*)s);
        free((void*)u);
        free((void*)vt);
        free((void*)work);
        return; /*failed to converge, output 0s */
    }
    int incx=1;
    for(i=0; i<k; i++){
        if(s[i] > 1.0e-9f)
            ss=1.0f/s[i];
        else
            ss=s[i];
        cblas_dscal(m, ss, &u[i*m], incx);
    }
    inva = (double*)malloc(n*m*sizeof(double));
    int ld_inva=n;
    cblas_dgemm( CblasColMajor, CblasTrans, CblasTrans, n, m, k, 1.0f,
                vt, ldvt,
                u, ldu, 0.0f,
                inva, ld_inva);
    
    for(i=0; i<m; i++)
        for(j=0; j<n; j++)
            outM[j*m+i] = inva[i*n+j]; /* return in row-major order */
    
    /* clean-up */
    free((void*)a);
    free((void*)s);
    free((void*)u);
    free((void*)vt);
    free((void*)inva);
    free((void*)work);
}

/*------------------------------- Cholesky factorisation (?chol) -----------------------------*/

void utility_schol
(
    const float* A,
    const int dim,
    float* X
)
{
    int i, j, info, n, lda;
    n = lda = dim;
    float* a;
    a = malloc(dim*dim*sizeof(float));
    
    /* store in column major order */
    for(i=0; i<dim; i++)
        for(j=0; j<dim; j++)
            a[j*dim+i] = A[i*dim+j];
    
    /* a is replaced by solution */
    spotrf_( "U", &n, a, &lda, &info );
    
    if(info>0){
        /* A is not positive definate, solution not possible */
        memset(X, 0, dim*dim*sizeof(float));
    }
    else{
        /* store solution in row-major order */
        for(i=0; i<dim; i++)
            for(j=0; j<dim; j++)
                X[i*dim+j] = j>=i ? a[j*dim+i] : 0.0f;
    }
    
    free(a);
}

void utility_cchol
(
    const float_complex* A,
    const int dim,
    float_complex* X
)
{
    int i, j, info, n, lda;
    n = lda = dim;
    float_complex* a;
    a = malloc(dim*dim*sizeof(float_complex));
    
    /* store in column major order */
    for(i=0; i<dim; i++)
        for(j=0; j<dim; j++)
            a[j*dim+i] = A[i*dim+j];
    
    /* a is replaced by solution */
#if defined(__APPLE__) && !defined(SAF_USE_INTEL_MKL)
    cpotrf_( "U", &n, (__CLPK_complex*)a, &lda, &info );
#elif INTEL_MKL_VERSION
    cpotrf_( "U", &n, (MKL_Complex8*)a, &lda, &info );
#endif
    
    if(info>0){
        /* A is not positive definate, solution not possible */
        memset(X, 0, dim*dim*sizeof(float_complex));
    }
    else{
        /* store solution in row-major order */
        for(i=0; i<dim; i++)
            for(j=0; j<dim; j++)
                X[i*dim+j] = j>=i ? a[j*dim+i] : cmplxf(0.0f, 0.0f);
    }
    
    free(a);
}

/*-------------------------------- matrix inversion (?inv) ----------------------------------*/

// TODO: rewrite for row-major
void utility_sinv(float * A, const int N)
{
    int *IPIV;
    IPIV = malloc(N * sizeof(int));
    int LWORK = N*N;
    float *WORK;
    WORK = (float*)malloc(LWORK * sizeof(float));
    int INFO;
    
#if defined(__APPLE__) && !defined(SAF_USE_INTEL_MKL)
    sgetrf_((__CLPK_integer*)&N, (__CLPK_integer*)&N, A, (__CLPK_integer*)&N, (__CLPK_integer*)IPIV, (__CLPK_integer*)&INFO);
    sgetri_((__CLPK_integer*)&N, A, (__CLPK_integer*)&N, (__CLPK_integer*)IPIV, WORK, (__CLPK_integer*)&LWORK, (__CLPK_integer*)&INFO);
#elif INTEL_MKL_VERSION
    sgetrf_(&N, &N, A, &N, IPIV, &INFO);
    sgetri_(&N, A, &N, IPIV, WORK, &LWORK, &INFO);
#endif
    
    free(IPIV);
    free(WORK);
}

void utility_dinv(double* A, const int N)
{
    int *IPIV;
    IPIV = malloc( (N+1)*sizeof(int));
    int LWORK = N*N;
    double *WORK;
    WORK = malloc( LWORK*sizeof(double));
    int INFO;
    
#if defined(__APPLE__) && !defined(SAF_USE_INTEL_MKL)
    dgetrf_((__CLPK_integer*)&N, (__CLPK_integer*)&N, A, (__CLPK_integer*)&N, (__CLPK_integer*)IPIV, (__CLPK_integer*)&INFO);
    dgetri_((__CLPK_integer*)&N, A, (__CLPK_integer*)&N, (__CLPK_integer*)IPIV, WORK, (__CLPK_integer*)&LWORK, (__CLPK_integer*)&INFO);
#elif INTEL_MKL_VERSION
    dgetrf_(&N, &N, A, &N, IPIV, &INFO);
    dgetri_(&N, A, &N, IPIV, WORK, &LWORK, &INFO);
#endif
    
    free((void*)IPIV);
    free((void*)WORK);
}

void utility_cinv(float_complex * A, const int N)
{
    int *IPIV;
    IPIV = malloc(N * sizeof(int));
    int LWORK = N*N;
    float_complex *WORK;
    WORK = (float_complex*)malloc(LWORK * sizeof(float_complex));
    int INFO;
    
#if defined(__APPLE__) && !defined(SAF_USE_INTEL_MKL)
    cgetrf_((__CLPK_integer*)&N, (__CLPK_integer*)&N, (__CLPK_complex*)A, (__CLPK_integer*)&N, (__CLPK_integer*)IPIV, (__CLPK_integer*)&INFO);
    cgetri_((__CLPK_integer*)&N, (__CLPK_complex*)A, (__CLPK_integer*)&N, (__CLPK_integer*)IPIV, (__CLPK_complex*)WORK, (__CLPK_integer*)&LWORK, (__CLPK_integer*)&INFO);
#elif INTEL_MKL_VERSION
    cgetrf_(&N, &N, (MKL_Complex8*)A, &N, IPIV, &INFO);
    cgetri_(&N, (MKL_Complex8*)A, &N, IPIV, (MKL_Complex8*)WORK, &LWORK, &INFO);
#endif
    
    free(IPIV);
    free(WORK);
}



