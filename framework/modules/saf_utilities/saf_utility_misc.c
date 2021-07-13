/*
 * Copyright 2020 Leo McCormack
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
 * @file saf_utility_misc.c
 * @ingroup Utilities
 * @brief A collection of miscellaneous functions
 *
 * @author Leo McCormack
 * @date 29.01.2020
 * @license ISC
 */

#include "saf_utilities.h"
#include "saf_externals.h"

/**
 * Precomputed factorials for up to !15 (i.e. the "getSH" functions will employ
 * these up to 7th order) */
static const long double factorials_15[15] =
{1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0, 362880.0, 3628800.0, 39916800.0, 479001600.0, 6.2270208e9, 8.71782891e10};

/** Helper function for findCombinations() */
static void combinationUtil(int* arr, int* data, int start, int end, int index, int r, int** comb, int* nComb) {
    if (index == r) {
        (*nComb)++;
        (*comb) = realloc1d( (*comb), (*nComb)*r*sizeof(int));
        for (int j=0; j<r; j++)
            (*comb)[((*nComb)-1)*r+j] = data[j];
        return;
    }
    for (int i=start; i<=end && end-i+1 >= r-index; i++) {
        data[index] = arr[i];
        combinationUtil(arr, data, i+1, end, index+1, r, comb, nComb);
    }
}

void convert_0_360To_m180_180
(
    float* dirs_deg,
    int nDirs
)
{
    int i;
    for(i=0; i<nDirs; i++){
        if(dirs_deg[i*2]>180.0f)
            dirs_deg[i*2] = -360.0f + dirs_deg[i*2];
    }
}

int nextpow2
(
    int numsamp
)
{
    int npts_max;

    if (numsamp > INT_MAX)
        return 0;
    npts_max = 1;
    while( 1 ){
        npts_max *= 2;
        if (npts_max >= numsamp)
            return npts_max;
    }
}

void lagrangeWeights
(
    int N,
    float* x,
    int len_x,
    float* weights
)
{
    int k, i, l;

    for(l=0; l<len_x; l++){
        for (k=0; k<=N; k++)
            weights[k*len_x+l] = 1.0f;
        for (k=0; k<=N; k++){
            for(i=0; i<=N; i++)
                if(k!=i)
                    weights[i*len_x+l] *= ((x[l]-(float)k) / (float)(i-k));
        }
    }
}

void findERBpartitions
(
    float* centerFreq,
    int nBands,
    float maxFreqLim,
    int** erb_idx,
    float** erb_freqs,
    int* nERBBands
)
{
    int i, band, counter, next_erb_idx;
    float band_centreFreq, erb, erb_centre, tmp;

    band_centreFreq = (powf(2.0f, 1.0f/3.0f)+1.0f)/2.0f;
    free(*erb_idx);
    free(*erb_freqs);
    (*erb_idx) = malloc1d(sizeof(int));
    (*erb_freqs) = malloc1d(sizeof(float));
    (*erb_idx)[0] = 1;
    (*erb_freqs)[0] = centerFreq[0];
    counter = 0;
    next_erb_idx = 0;
    while((*erb_freqs)[counter]<maxFreqLim){
        erb = 24.7f + 0.108f * (*erb_freqs)[counter] * band_centreFreq;
        (*erb_idx) = realloc1d((*erb_idx), (counter+2)*sizeof(int));
        (*erb_freqs) = realloc1d((*erb_freqs), (counter+2)*sizeof(float));
        (*erb_freqs)[counter+1] = (*erb_freqs)[counter] + erb;
        erb_centre = FLT_MAX;
        /*  find closest band frequency as upper partition limit */
        for(band=0; band<nBands; band++){
            tmp =fabsf((*erb_freqs)[counter+1] - centerFreq[band]);
            if(tmp <erb_centre){
                erb_centre = tmp;
                next_erb_idx = band;
            }
        }
        (*erb_idx)[counter+1] = next_erb_idx + 1;
        if((*erb_idx)[counter+1] == (*erb_idx)[counter])
            (*erb_idx)[counter+1] = (*erb_idx)[counter+1]+1;
        (*erb_freqs)[counter+1] = centerFreq[(*erb_idx)[counter+1]-1];
        counter++;
    }
    /* last limit set at last band */
    (*erb_idx) = realloc1d((*erb_idx), (counter + 2) * sizeof(int));
    (*erb_freqs) = realloc1d((*erb_freqs), (counter + 2) * sizeof(float));
    (*erb_idx)[counter+1] = nBands;
    (*erb_freqs)[counter+1] = centerFreq[nBands-1];
    (*nERBBands) = counter+2;

    /* subtract 1 from the indices (the above is a direct port from Matlab...) */
    for(i=0; i<(*nERBBands); i++)
        (*erb_idx)[i] = (*erb_idx)[i] - 1;
}

void randperm
(
    int len,
    int* randperm
)
{
    int i, j, tmp;

    for (i = 0; i < len; i++)
        randperm[i] = i;
    for (i = 0; i < len; i++) {
        j = rand() % (len-i) + i;
        tmp = randperm[j];
        randperm[j] = randperm[i];
        randperm[i] = tmp;
    }
}

long double factorial(int n)
{
    int i;
    long double ff;
    if(n<15)
        return factorials_15[n];
    else{
        ff = 1.0;
        for(i = 1; i<=n; i++)
            ff *= (long double)i;
        return ff;
    }
}

float matlab_fmodf(float x, float y) {
    float tmp = fmodf(x, y);
    return tmp >= 0 ? tmp : tmp + y;
}

void cxcorr
(
    float* a,
    float* b,
    float* x_ab,
    size_t la,
    size_t lb
)
{
    int m, n, negFLAG, arg, len, lim;
    
    len = (int)(la + lb) - 1;
    memset(x_ab, 0, len*sizeof(float));
    for(m=1; m<=len; m++){
        arg = m-(int)la;
        if(arg<0){
            negFLAG = 1;
            lim = (int)la + arg;
        }
        else{
            negFLAG = 0;
            lim = (int)la - arg;
        }
        for(n=1; n<=lim; n++){
            if(negFLAG == 0)
                x_ab[m-1] += (a[arg+n-1] * b[n-1]);
            else
                x_ab[m-1] += (a[n-1] * b[n-arg-1]);
        }
    }
}

void rand_m1_1
(
    float* vector,
    int length
)
{
    int i;
    for(i=0; i<length; i++)
        vector[i] = (2.0f*rand()/(float)RAND_MAX)-1.0f;
}

void rand_cmplx_m1_1
(
    float_complex* vector,
    int length
)
{
    int i;
    for(i=0; i<length; i++)
        vector[i] = cmplxf((2.0f*rand()/(float)RAND_MAX)-1.0f, (2.0f*rand()/(float)RAND_MAX)-1.0f);
}

void rand_0_1
(
    float* vector,
    int length
)
{
    int i;
    for(i=0; i<length; i++)
        vector[i] = rand()/(float)RAND_MAX;
}

void convd
(
    double* x,
    double* h,
    int len_x,
    int len_h,
    double* y
)
{
    int i, j, h_start, x_start, x_end, len_y;

    len_y = len_h+len_x-1;
    memset(y, 0, len_y*sizeof(double));
    for (i=0; i<len_y; i++) {
        x_start = SAF_MAX(0,i-len_h+1);
        x_end   = SAF_MIN(i+1,len_x);
        h_start = SAF_MIN(i,len_h-1);
        for(j=x_start; j<x_end; j++)
            y[i] += h[h_start--]*x[j];
    }
}

void convz
(
    double_complex* x,
    double_complex* h,
    int len_x,
    int len_h,
    double_complex* y
)
{
    int i, j, h_start, x_start, x_end, len_y;

    len_y = len_h+len_x-1;
    memset(y, 0, len_y*sizeof(double_complex));
    for (i=0; i<len_y; i++) {
        x_start = SAF_MAX(0,i-len_h+1);
        x_end   = SAF_MIN(i+1,len_x);
        h_start = SAF_MIN(i,len_h-1);
        for(j=x_start; j<x_end; j++)
            y[i] = ccadd(y[i], ccmul(h[h_start--], x[j]));
    }
}

void polyd_v
(
    double* x,
    double* poly,
    int len_x
)
{
    int j,i;

    memset(poly, 0, (len_x+1)*sizeof(double));
    poly[0] = 1.0;
    for (j=0; j<len_x; j++){
        for(i=j+1; i>0; i--){
            poly[i] = poly[i] - x[j] * (poly[i-1]);
        }
    } 
}

void polyz_v
(
    double_complex* x,
    double_complex* poly,
    int len_x
)
{
    int j,i;

    memset(poly, 0, (len_x+1)*sizeof(double_complex));
    poly[0] = cmplx(1.0, 0.0);
    for (j=0; j<len_x; j++){
        for(i=j+1; i>0; i--){
            poly[i] = ccsub(poly[i], ccmul(x[j], poly[i-1]));
        }
    }
}

void polyd_m
(
    double* X,
    double_complex* poly,
    int size_x
)
{
    int j,i;
    double_complex* Xcmplx, *e;

    /* Characteristic polynomial */
    Xcmplx = malloc1d(size_x*size_x*sizeof(double_complex));
    e = malloc1d(size_x*sizeof(double_complex));
    for(i=0; i<size_x*size_x; i++)
        Xcmplx[i] = cmplx(X[i], 0.0);
    utility_zeig(NULL, Xcmplx, size_x, NULL, NULL, NULL, e);

    /* recursion formula */
    memset(poly, 0, (size_x+1)*sizeof(double_complex));
    poly[0] = cmplx(1.0, 0.0);
    for (j=0; j<size_x; j++){
        for(i=j+1; i>0; i--){
            poly[i] = ccsub(poly[i], ccmul(e[j], poly[i-1]));
        }
    }

    /* clean-up */
    free(Xcmplx);
    free(e);
}

float sumf
(
    float* values,
    int nValues
)
{
    float sum;
#if defined(SAF_USE_INTEL_IPP)
    ippsSum_32f((Ipp32f*)values, nValues, &sum, ippAlgHintNone);
#else
    int i;
    sum = 0.0f;
    for(i=0; i<nValues; i++)
        sum += values[i];
#endif
    return sum;
}

int anyLessThanf
(
    float* values,
    int nValues,
    float threshold
)
{
    int i;
    for(i=0; i<nValues; i++)
        if(values[i]<threshold)
            return 1;
    return 0;
}

void unique_i
(
    int* input,
    int nInputs,
    int** uniqueVals,
    int** uniqueInds,
    int* nUnique
)
{
    int i, j, k, nDups, foundDups, foundNewDup;
    int* dups, *nDuplicates_perInput;

    /* If only 1 input... */
    if(nInputs == 1){
        (*nUnique) = 1;
        if(uniqueVals!=NULL){
            (*uniqueVals) = malloc1d((*nUnique)*sizeof(int));
            (*uniqueVals)[0] = input[0];
        }
        if(uniqueInds!=NULL){
            (*uniqueInds) = malloc1d((*nUnique)*sizeof(int));
            (*uniqueInds)[0] = 0;
        }
    }

    /* prep */
    dups = malloc1d(nInputs*sizeof(int));
    nDuplicates_perInput = calloc1d(nInputs, sizeof(int));
    (*nUnique) = nInputs;

    /* Find duplicates */
    nDups = 0;
    for(i=0; i<nInputs; i++) {
        foundDups = 0;
        for(j=i+1; j<nInputs; j++) {
            if(input[i] == input[j]) {
                foundNewDup = 1;
                for(k=0; k<nDups; k++)
                    if(input[i]==dups[k])
                        foundNewDup = 0;

                /* input value is repeated: */
                nDuplicates_perInput[i]++;
                if(foundNewDup){
                    (*nUnique)--;  /* one less unique value... */
                    foundDups = 1;
                }
            }
        }
        /* Store repeated value, so "nUnique" is not decreased more than once for it */
        if(foundDups){
            dups[nDups] = input[i];
            nDups++;
        }
    }
    saf_assert((*nUnique)>-1, "Something very bad happened");
    free(dups);

    /* If no unique values were found */
    if((*nUnique)==0){
        (*uniqueVals) = NULL;
        (*uniqueInds) = NULL;
        (*nUnique) = 0;
        free(nDuplicates_perInput);
        return;
    }

    /* Save unique values and/or their indices */
    j=0;
    if(uniqueVals!=NULL)
        (*uniqueVals) = malloc1d((*nUnique)*sizeof(int));
    if(uniqueInds!=NULL)
        (*uniqueInds) = malloc1d((*nUnique)*sizeof(int));
    for(i=0; i<nInputs; i++) {
        if(nDuplicates_perInput[i] == 0) {
            /* If no duplicate exists, append... */
            if(uniqueVals!=NULL)
                (*uniqueVals)[j] = input[i];
            if(uniqueInds!=NULL)
                (*uniqueInds)[j] = i;
            j++;
        }
    }

    /* clean-up */
    free(nDuplicates_perInput);
}

void findCombinations
(
    int* arrValues,
    int nValues,
    int nElements,
    int** comb,
    int* nComb
)
{
    int* data;
    data = malloc1d(nElements*sizeof(int));
    saf_assert(*comb==NULL, "comb must be empty and NULL prior to calling findCombinations()");
    (*nComb) = 0;
    combinationUtil(arrValues, data, 0, nValues-1, 0, nElements, comb, nComb);
    free(data);
}
 
/* Based heavily on the Matlab script found here:
 * https://se.mathworks.com/matlabcentral/fileexchange/50413-generalized-matrix-exponential
 * Copyright (c) 2015, Kenneth Johnson (BSD-3-clause license) */
void gexpm
(
    float* D,
    int sizeD,
    int m1,
    float* Y
)
{
    int i, j, k;
    float tol, s, h2, h, hh, hhh, two;
    float** D_2, **D_3, **D_4, **D_5, **Dh, **Ym1, **Ym2;

    tol = FLT_EPSILON;

    /* Calculate Y = expm(D), Ym1 = Y-I.
     *
     * Scale and square: Y = expm(D/n)^n; n = 2^s (non-negative integer s)
     *
     * Pade approximation: expm(D/n) = R
     *   = (I-Dh+(2/5)*Dh^2-(1/15)*Dh^3)\(I+Dh+(2/5)*Dh^2+(1/15)*Dh^3)
     * where Dh = D*h; h = 1/(2*n) (h = integration half-step)
     *
     * Pade approximation error: R-expm(D/n) = (-2/1575)*(Dh)^7
     *
     * Set Y = R^n. Approximate absolute error:
     *   Y-expm(D) = R^n-expm(D/n)^n = n*R^(n-1)*(R-expm(D/n))
     *     = n*Y*R^(-1)*(-2/1575)*(Dh)^7
     * R^(-1) is close to I, so
     *   Y-expm(D) = n*Y*(-2/1575)*(Dh)^7 (approx)
     *
     * Large D (|D|>=1): Bound the relative error magnitude by tol,
     *   n*(2/1575)*|(Dh)^7| <= tol
     *
     * Small D (|D|<1): |Y-I| is of order 1 or less; the absolute error
     * Y-expm(D) is of order n*(-2/1575)*(Dh)^7. Bound the error magnitude
     * by tol*|D| to preserve relative accuracy of Ym1:
     *   n*(2/1575)*|(Dh)^7| <= tol*|D|
     *
     * Combine large/small Dh conditions conjunctively:
     *   n*(2/1575)*|(Dh)^7| <= tol*min(1,|D|)
     *
     * Substitute h = 1/(2*n):
     *   (2*n)^6 >= |D^7|/(1575*tol*min(1,|D|))
     * Substitute n = 2^s:
     *   s >= log2(|D^7|/(1575*tol*min(1,|D|)))/6-1
     * (Use the Frobenius norm for |...| to preserve symmetry of expm under
     * matrix transposition.)
     */
    D_2 = (float**)malloc2d(sizeD, sizeD, sizeof(float));
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, sizeD, sizeD, sizeD, 1.0f,
                D, sizeD,
                D, sizeD, 0.0f,
                FLATTEN2D(D_2), sizeD);
    D_3 = (float**)malloc2d(sizeD, sizeD, sizeof(float));
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, sizeD, sizeD, sizeD, 1.0f,
                FLATTEN2D(D_2), sizeD,
                D, sizeD, 0.0f,
                FLATTEN2D(D_3), sizeD);
    D_4 = (float**)malloc2d(sizeD, sizeD, sizeof(float));
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, sizeD, sizeD, sizeD, 1.0f,
                FLATTEN2D(D_3), sizeD,
                FLATTEN2D(D_3), sizeD, 0.0f,
                FLATTEN2D(D_4), sizeD);
    D_5 = (float**)malloc2d(sizeD, sizeD, sizeof(float));
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, sizeD, sizeD, sizeD, 1.0f,
                FLATTEN2D(D_4), sizeD,
                D, sizeD, 0.0f,
                FLATTEN2D(D_5), sizeD); 
    s = ceilf(log2f(Frob_norm(FLATTEN2D(D_5), sizeD, sizeD)/
                    (1575.0f*tol*SAF_MIN(1.0f,Frob_norm(D, sizeD, sizeD))))/6.0f-1.0f);
    s = SAF_MAX(s, 0.0f);

    /* Get Pade approximation for expm(D*h2) = Y =
     *   (I-Dh+(2/5)*Dh^2-(1/15)*Dh^3)\(I+Dh+(2/5)*Dh^2+(1/15)*Dh^3)
     * Ym1 = Y-I =
     *   (I-Dh+(2/5)*Dh^2-(1/15)*Dh^3)\(2*(Dh+(1/15)*Dh^3))
     * (Calculate Ym1, not Y, to avoid precision loss from dominant I
     * terms when Dh is small.) */
    h2 = powf(2.0f,-s);
    h = h2/2.0f;
    hh = h*h;
    hhh = hh*h;
    Dh = (float**)malloc2d(sizeD, sizeD, sizeof(float));
    memcpy(FLATTEN2D(Dh), D, sizeD*sizeD*sizeof(float));
    utility_svsmul(FLATTEN2D(Dh), &h, sizeD*sizeD, NULL);
    utility_svsmul(FLATTEN2D(D_2), &hh, sizeD*sizeD, NULL);
    utility_svsmul(FLATTEN2D(D_3), &hhh, sizeD*sizeD, NULL);
    Ym1 = (float**)malloc2d(sizeD, sizeD, sizeof(float));
    for(i=0; i<sizeD; i++)
        for(j=0; j<sizeD; j++)
            Ym1[i][j] = Dh[i][j] + (1.0f/15.0f)*D_3[i][j];
    Ym2 = (float**)malloc2d(sizeD, sizeD, sizeof(float));
    for(i=0; i<sizeD; i++){
        for(j=0; j<sizeD; j++){
            Ym2[i][j] = (2.0f/5.0f)*D_2[i][j]-Ym1[i][j];
            if(i==j)
                Ym2[i][j] += 1.0f;
        }
    }
    two = 2.0f;
    utility_svsmul(FLATTEN2D(Ym1), &two, sizeD*sizeD, NULL);
    utility_sglslv(NULL, FLATTEN2D(Ym2), sizeD, FLATTEN2D(Ym1), sizeD, FLATTEN2D(Ym1));

    /* Y = Ym1+I = expm(D)
     * Square Y (i.e., Y <-- Y*Y) s times to get Y = expm(D*2^s). */
    for (k = 0; k<(int)s; k++){
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, sizeD, sizeD, sizeD, 1.0f,
                    FLATTEN2D(Ym1), sizeD,
                    FLATTEN2D(Ym1), sizeD, 0.0f,
                    FLATTEN2D(Ym2), sizeD);
        for(i=0; i<sizeD; i++)
            for(j=0; j<sizeD; j++)
                Ym1[i][j] = Ym2[i][j] + 2.0f*Ym1[i][j]; /* (Ym1+I) <-- (Ym1+I)*(Ym1+I) */
    }
    memcpy(Y, FLATTEN2D(Ym1), sizeD*sizeD*sizeof(float));
    if (m1){}
    else{
        for(i=0; i<sizeD; i++)
           for(j=0; j<sizeD; j++)
               if(i==j)
                   Y[i*sizeD+j] += 1.0f;
    }

    /* clean-up */
    free(D_2);
    free(D_3);
    free(D_4);
    free(D_5);
    free(Dh);
    free(Ym1);
    free(Ym2);
}
