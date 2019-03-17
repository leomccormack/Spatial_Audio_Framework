/*
 Copyright 2018 Leo McCormack
 
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
 *     saf_hoa.c
 * Description:
 *     A Higher-order Ambisonics C library; largely derived from the MatLab library by
 *     Archontis Politis: https://github.com/polarch/Higher-Order-Ambisonics
 * Dependencies:
 *     saf_utilities, saf_sh, saf_vbap
 * Author, date created:
 *     Leo McCormack, 19.03.2018
 */

#include "saf_hoa.h"
#include "saf_hoa_internal.h"

/* Zotter, F., Frank, M. (2012). All-Round Ambisonic Panning and Decoding. Journal of the Audio Engineering Society, 60(10), 807?820. */
void getMaxREweights
(
    int order,
    float* a_n /* (order+1)^2 x (order+1)^2 */
)
{
    int n, i, idx, nSH;
    double x;
    double* ppm;
    
    x = cosf(137.9f*(M_PI/180.0f)/((float)order+1.51f));
    nSH = (order+1)*(order+1);
    memset(a_n, 0, nSH*nSH*sizeof(float));
    ppm = calloc((order+1),sizeof(double));
    idx = 0;
    for(n=0; n<=order; n++){
        unnorm_legendreP(n, &x, 1, ppm);
        /* store the first Legendre polynomial value for each order along the diagonal of a_n */
        for(i = 0; i<2*n+1; i++)
            a_n[(idx+i)*nSH + (idx+i)] = (float)ppm[0];
        idx += 2*n+1;
    }
    free(ppm);
}

void getAmbiDecoder
(
    float* ls_dirs_deg,
    int nLS,
    AMBI_DECODER_METHODS method,
    int order,
    float** decMtx
)
{
    int i, j, nSH;
    float* Y_ls;
    
    nSH = (order+1) * (order+1);
    Y_ls = NULL;
    free(*decMtx);
    (*decMtx) = malloc(nLS*nSH*sizeof(float));
    switch(method){
        default:
        case DECODER_DEFAULT:
        case DECODER_SAD:
            /* Sampling Ambisonic Decoder (SAD) is simply the loudspeaker spherical harmonic matrix scaled by the
             * number of loudspeakers. */
            getRSH(order, ls_dirs_deg, nLS, &Y_ls);
            for(i=0; i<nLS; i++)
                for(j=0; j<nSH; j++)
                    (*decMtx)[i*nSH+j] = Y_ls[j*nLS + i]/(float)nLS;
            free(Y_ls);
            break;
           
        case DECODER_MMD:
            /* Mode-Matching Decoder (MMD) is simply the psuedo inverse of the loudspeaker spherical harmonic
             * matrix. */
            getRSH(order, ls_dirs_deg, nLS, &Y_ls);
            utility_spinv(Y_ls, nSH, nLS, (*decMtx));
            free(Y_ls);
            break;
            
        case DECODER_EPAD:
            getEPAD(order, ls_dirs_deg, nLS, decMtx);
            break;
            
        case DECODER_ALLRAD:
            getAllRAD(order, ls_dirs_deg, nLS, decMtx);
            break;
    }
}

void getBinauralAmbiDecoder
(
    float_complex* hrtfs,
    float* hrtf_dirs_deg,
    int N_dirs,
    int N_bands,
    BINAURAL_AMBI_DECODER_METHODS method,
    int order,
    float* freqVector,
    float* itd_s,
    float* weights,
    float_complex* decMtx
)
{
    int nSH;
    
    nSH = (order+1) * (order+1);
    switch(method){
        default:
        case BINAURAL_DECODER_DEFAULT:
        case BINAURAL_DECODER_LS:
            getBinDecoder_LS(hrtfs, hrtf_dirs_deg, N_dirs, N_bands, order, weights, decMtx);
            break;
            
        case BINAURAL_DECODER_LSDIFFEQ:
            getBinDecoder_LSDIFFEQ(hrtfs, hrtf_dirs_deg, N_dirs, N_bands, order, weights, decMtx);
            break;
            
        case BINAURAL_DECODER_SPR:
            getBinDecoder_SPR(hrtfs, hrtf_dirs_deg, N_dirs, N_bands, order, weights, decMtx);
            break;
            
        case BINAURAL_DECODER_TAC:
            getBinDecoder_TAC(hrtfs, hrtf_dirs_deg, N_dirs, N_bands, order, freqVector, itd_s, weights, decMtx);
            //applyDiffCovMatching(hrtfs, hrtf_dirs_deg, N_dirs, N_bands, order, weights, decMtx);
            break;
            
        case BINAURAL_DECODER_MAGLS:
            getBinDecoder_MAGLS(hrtfs, hrtf_dirs_deg, N_dirs, N_bands, order, freqVector, weights, decMtx);
            //applyDiffCovMatching(hrtfs, hrtf_dirs_deg, N_dirs, N_bands, order, weights, decMtx);
            break;
    }
}

void applyDiffCovMatching
(
    float_complex* hrtfs,
    float* hrtf_dirs_deg,
    int N_dirs,
    int N_bands,
    int order,
    float* weights,
    float_complex* decMtx
)
{
    int i, nSH, band;
    float* Y_tmp;
    float_complex* W, *Y_na, *H_W, *H_ambi, *decMtx_diffMatched;
    float_complex C_ref[2][2], C_ambi[2][2], X[2][2], X_ambi[2][2], XH_Xambi[2][2], U[2][2], V[2][2], UX[2][2], VUX[2][2], M[2][2];
    const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);
    
    nSH = (order+1)*(order+1);
    
    /* integration weights */
    W = calloc(N_dirs*N_dirs, sizeof(float_complex));
    if(weights!=NULL)
        for(i=0; i<N_dirs; i++)
            W[i*N_dirs+i] = cmplxf(weights[i], 0.0f);
    else
        for(i=0; i<N_dirs; i++)
            W[i*N_dirs+i] = cmplxf(1.0f/(float)N_dirs, 0.0f);
    
    /* SH */
    Y_tmp = NULL;
    Y_na = malloc(nSH*N_dirs*sizeof(float_complex));
    getRSH(order, hrtf_dirs_deg, N_dirs, &Y_tmp);
    for(i=0; i<nSH*N_dirs; i++)
        Y_na[i] = cmplxf(Y_tmp[i], 0.0f);
    free(Y_tmp);
    
    /* apply diffuse-field coherence matching per band */
    H_W = malloc(2*N_dirs*sizeof(float_complex));
    H_ambi = malloc(2*N_dirs*sizeof(float_complex));
    decMtx_diffMatched = malloc(2*nSH*sizeof(float_complex));
    for(band=0; band<N_bands; band++){
        /* Diffuse-field responses */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, N_dirs, N_dirs, &calpha,
                    &hrtfs[band*2*N_dirs], N_dirs,
                    W, N_dirs, &cbeta,
                    H_W, N_dirs);
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, 2, 2, N_dirs, &calpha,
                    H_W, N_dirs,
                    &hrtfs[band*2*N_dirs], N_dirs, &cbeta,
                    (float_complex*)C_ref, 2);
        for(i=0; i<2; i++)
            C_ref[i][i] = cmplxf(crealf(C_ref[i][i]), 0.0f); /* force diagonal to be real */
        utility_cchol((float_complex*)C_ref, 2, (float_complex*)X);
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, N_dirs, nSH, &calpha,
                    &decMtx[band*2*nSH], nSH,
                    Y_na, N_dirs, &cbeta,
                    H_ambi, N_dirs);
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, N_dirs, N_dirs, &calpha,
                    H_ambi, N_dirs,
                    W, N_dirs, &cbeta,
                    H_W, N_dirs);
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, 2, 2, N_dirs, &calpha,
                    H_W, N_dirs,
                    H_ambi, N_dirs, &cbeta,
                    (float_complex*)C_ambi, 2);
        for(i=0; i<2; i++)
            C_ambi[i][i] = cmplxf(crealf(C_ambi[i][i]), 0.0f); /* force diagonal to be real */
        utility_cchol((float_complex*)C_ambi, 2, (float_complex*)X_ambi);
        
        /* SVD */
        cblas_cgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans, 2, 2, 2, &calpha,
                    (float_complex*)X_ambi, 2,
                    (float_complex*)X, 2, &cbeta,
                    (float_complex*)XH_Xambi, 2);
        utility_csvd((float_complex*)XH_Xambi, 2, 2, (float_complex*)U, NULL, (float_complex*)V, NULL);
        
        /* apply matching */
        cblas_cgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans, 2, 2, 2, &calpha,
                    (float_complex*)U, 2,
                    (float_complex*)X, 2, &cbeta,
                    (float_complex*)UX, 2);
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, 2, 2, &calpha,
                    (float_complex*)V, 2,
                    (float_complex*)UX, 2, &cbeta,
                    (float_complex*)VUX, 2);
        utility_cglslv((float_complex*)X_ambi, 2, (float_complex*)VUX, 2, (float_complex*)M);
        cblas_cgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans, 2, nSH, 2, &calpha,
                    (float_complex*)M, 2,
                    &decMtx[band*2*nSH], nSH, &cbeta,
                    decMtx_diffMatched, nSH);
        memcpy(&decMtx[band*2*nSH], decMtx_diffMatched, 2*nSH*sizeof(float_complex));
    }
    
    free(W);
    free(Y_na);
    free(H_W);
    free(H_ambi);
    free(decMtx_diffMatched);
}


















