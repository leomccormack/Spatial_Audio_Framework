/*
 * Copyright 2018 Leo McCormack
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
 * @file saf_hoa.c
 * @ingroup HOA
 * @brief Public source for the higher-order Ambisonics module (#SAF_HOA_MODULE)
 *
 * A collection of Ambisonics related functions. Many of which are derived from
 * the MATLAB library found in [1].
 *
 * @see [1] https://github.com/polarch/Higher-Order-Ambisonics
 *          Copyright (c) 2015, Archontis Politis, BSD-3-Clause License
 *
 * @author Leo McCormack
 * @date 19.03.2018
 * @license ISC
 */

#include "saf_hoa.h"
#include "saf_hoa_internal.h"

/* ========================================================================== */
/*                               Misc. Functions                              */
/* ========================================================================== */

void convertHOAChannelConvention
(
    float* insig,
    int order,
    int signalLength,
    HOA_CH_ORDER inConvention,
    HOA_CH_ORDER outConvention
)
{
    int i, nSH;

    nSH = ORDER2NSH(order);

    /* bypass, if 0th order, or no conversion required */
    if(order==0 || inConvention == outConvention)
        return;

    /* Convert FUMA to ACN */
    if(inConvention==HOA_CH_ORDER_FUMA && outConvention==HOA_CH_ORDER_ACN){
        cblas_sswap(signalLength, &insig[1*signalLength], 1, &insig[3*signalLength], 1); /* Swap X and Z */
        cblas_sswap(signalLength, &insig[1*signalLength], 1, &insig[2*signalLength], 1); /* Swap Z and Y */
    }
    /* Convert ACN to FUMA */
    else if(inConvention==HOA_CH_ORDER_ACN && outConvention==HOA_CH_ORDER_FUMA){
        cblas_sswap(signalLength, &insig[1*signalLength], 1, &insig[2*signalLength], 1); /* Swap Y and Z */
        cblas_sswap(signalLength, &insig[1*signalLength], 1, &insig[3*signalLength], 1); /* Swap Z and X */
    }
    /* Fill any remaining channels with zeros */
    for(i=4; i<nSH; i++)
        memset(&insig[i*signalLength], 0, signalLength * sizeof(float));
}

void convertHOANormConvention
(
    float* insig,
    int order,
    int signalLength,
    HOA_NORM inConvention,
    HOA_NORM outConvention
)
{
    int n, ch;

    if(order==0 || inConvention == outConvention)
        return; /* Nothing to do */

    if(inConvention==HOA_NORM_N3D){
        if(outConvention == HOA_NORM_SN3D){
            for (n = 0; n<order+1; n++)
                for (ch = ORDER2NSH(n-1); ch < ORDER2NSH(n); ch++)
                    cblas_sscal(signalLength, 1.0f/sqrtf(2.0f*(float)n+1.0f), &insig[ch*signalLength], 1);
        }
        else if(outConvention==HOA_NORM_FUMA){
            cblas_sscal(signalLength, 1.0f/sqrtf(2.0f), insig, 1);
            for (ch = 1; ch<4 /* 1st order only */; ch++)
                cblas_sscal(signalLength, 1.0f/sqrtf(3.0f), &insig[ch*signalLength], 1);
        }
    }
    else if(inConvention==HOA_NORM_SN3D){
        if(outConvention == HOA_NORM_N3D){
            for (n = 0; n<order+1; n++)
                for (ch = ORDER2NSH(n-1); ch<ORDER2NSH(n); ch++)
                    cblas_sscal(signalLength, sqrtf(2.0f*(float)n+1.0f), &insig[ch*signalLength], 1);
        }
        else if(outConvention==HOA_NORM_FUMA)
            cblas_sscal(signalLength, 1.0f/sqrtf(2.0f), insig, 1);
    }
    else if(inConvention==HOA_NORM_FUMA){
        if(outConvention == HOA_NORM_N3D){
            cblas_sscal(signalLength, sqrtf(2.0f), insig, 1);
            for (ch = 1; ch<4 /* 1st order only */; ch++)
                cblas_sscal(signalLength, sqrtf(3.0f), &insig[ch*signalLength], 1);
        }
        else if(outConvention == HOA_NORM_SN3D)
            cblas_sscal(signalLength, sqrtf(2.0f), insig, 1);
    }
}

void getRSH
(
    int N,
    float* dirs_deg,
    int nDirs,
    float* Y
)
{
    int i, nSH;
    float scale;
    float* dirs_rad;

    if(nDirs<1)
        return;
    
    nSH = (N+1)*(N+1);
    scale = sqrtf(4.0f*SAF_PI);
    
    /* convert [azi, elev] in degrees, to [azi, inclination] in radians */
    dirs_rad = malloc1d(nDirs*2*sizeof(float));
    for(i=0; i<nDirs; i++){
        dirs_rad[i*2+0] = dirs_deg[i*2+0] * SAF_PI/180.0f;
        dirs_rad[i*2+1] = SAF_PI/2.0f - (dirs_deg[i*2+1] * SAF_PI/180.0f);
    }
    
    /* get real-valued spherical harmonics */
    getSHreal(N, dirs_rad, nDirs, Y);
    
    /* remove sqrt(4*pi) term */
    utility_svsmul(Y, &scale, nSH*nDirs, NULL);
    
    free(dirs_rad);
}

void getRSH_recur
(
    int N,
    float* dirs_deg,
    int nDirs,
    float* Y
)
{
    int n, m, i, dir, index_n;
    float Nn0, Nnm;
    float sleg_n[8], sleg_n_1[8], sleg_n_2[8], ssin_el, sfactorials_n[15];
    float* leg_n, *leg_n_1, *leg_n_2, *sin_el, *factorials_n;

    if(nDirs<1)
        return;

    if(N<=7 && nDirs == 1){
        /* Single direction optimisation for up to 7th order */
        leg_n = sleg_n;
        leg_n_1 = sleg_n_1;
        leg_n_2 = sleg_n_2;
        sin_el = &ssin_el;
        factorials_n = sfactorials_n;
    }
    else{
        factorials_n = malloc1d((2*N+1)*sizeof(float));
        leg_n = malloc1d((N+1)*nDirs * sizeof(float));
        leg_n_1 = malloc1d((N+1)*nDirs * sizeof(float));
        leg_n_2 = malloc1d((N+1)*nDirs * sizeof(float));
        sin_el = malloc1d(nDirs * sizeof(float));
    } 
    index_n = 0;
    
    /* precompute factorials */
    for (i = 0; i < 2*N+1; i++)
        factorials_n[i] = (float)factorial(i);
    
    /* cos(inclination) = sin(elevation) */
    for (dir = 0; dir<nDirs; dir++)
        sin_el[dir] = sinf(dirs_deg[dir*2+1] * SAF_PI/180.0f);
    
    /* compute SHs with the recursive Legendre function */
    for (n = 0; n<N+1; n++) {
        if (n==0) {
            for (dir = 0; dir<nDirs; dir++)
                Y[n*nDirs+dir] = 1.0f;
            index_n = 1;
        }
        else {
            unnorm_legendreP_recur(n, sin_el, nDirs, leg_n_1, leg_n_2, leg_n); /* does NOT include Condon-Shortley phase term */
            
            Nn0 = sqrtf(2.0f*(float)n+1.0f);
            for (dir = 0; dir<nDirs; dir++){
                for (m = 0; m<n+1; m++) {
                    if (m==0)
                        Y[(index_n+n)*nDirs+dir] = Nn0  * leg_n[m*nDirs+dir];
                    else {
                        Nnm = Nn0* sqrtf( 2.0f * factorials_n[n-m]/factorials_n[n+m] );
                        Y[(index_n+n-m)*nDirs+dir] = Nnm * leg_n[m*nDirs+dir] * sinf((float)m * (dirs_deg[dir*2])*SAF_PI/180.0f);
                        Y[(index_n+n+m)*nDirs+dir] = Nnm * leg_n[m*nDirs+dir] * cosf((float)m * (dirs_deg[dir*2])*SAF_PI/180.0f);
                    }
                }
            }
            index_n += 2*n+1;
        }
        utility_svvcopy(leg_n_1, (N+1)*nDirs, leg_n_2);
        utility_svvcopy(leg_n,   (N+1)*nDirs, leg_n_1);
    }
    
    if(N>7 || nDirs > 1){
        free(factorials_n);
        free(leg_n);
        free(leg_n_1);
        free(leg_n_2);
        free(sin_el);
    }
}


/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

void getMaxREweights
(
    int order,
    int diagMtxFlag,
    float* a_n
)
{
    int n, i, idx, nSH;
    double x;
    double* ppm;
    
    x = cosf(137.9f*(SAF_PI/180.0f)/((float)order+1.51f));
    nSH = ORDER2NSH(order);
    if(diagMtxFlag)
        memset(a_n, 0, nSH*nSH*sizeof(float));
    else
        memset(a_n, 0, nSH*sizeof(float));
    ppm = calloc1d((order+1),sizeof(double));
    idx = 0;
    for(n=0; n<=order; n++){
        unnorm_legendreP(n, &x, 1, ppm);
        /* store the first Legendre polynomial value for each order along the
         * diagonal of a_n */
        for(i = 0; i<2*n+1; i++){
            if(diagMtxFlag)
                a_n[(idx+i)*nSH + (idx+i)] = (float)ppm[0];
            else
                a_n[(idx+i)] = (float)ppm[0];
        }
        idx += 2*n+1;
    }
    free(ppm);
}

void truncationEQ
(
    float* w_n,
    int order_truncated,
    int order_target,
    double* kr,
    int nBands,
    float softThreshold,
    float* gain
)
{
    /* For more information refer to:
     * Hold, C., Gamper, H., Pulkki, V., Raghuvanshi, N., & Tashev, I. J. (2019).
     * Improving Binaural Ambisonics Decoding by Spherical Harmonics Domain Tapering
     * and Coloration Compensation. ICASSP, IEEE International Conference on Acoustics,
     * Speech and Signal Processing - Proceedings.
     */
    double_complex* b_n_target = calloc1d(nBands*(order_target+1), sizeof(double_complex));
    double_complex* b_n_truncated = calloc1d(nBands*(order_truncated+1), sizeof(double_complex));
    double* p_target = calloc1d(nBands, sizeof(double));
    double* p_truncated = calloc1d(nBands, sizeof(double));

    sphModalCoeffs(order_target, kr, nBands, ARRAY_CONSTRUCTION_RIGID, 0., b_n_target);
    sphModalCoeffs(order_truncated, kr, nBands, ARRAY_CONSTRUCTION_RIGID, 0., b_n_truncated);

    /* p_full */
    for (int idxBand=0; idxBand<nBands; idxBand++)
        for (int idxOrder=0; idxOrder<=order_target; idxOrder++)
            p_target[idxBand] += (2.0*idxOrder + 1.0) * pow(cabs(b_n_target[idxBand*(order_target+1) + idxOrder]), 2.0);
    /* p_truncated from (tapered) modal weighting */
    for (int idxBand=0; idxBand<nBands; idxBand++)
        for (int idxOrder=0; idxOrder<=order_truncated; idxOrder++)
            p_truncated[idxBand] += w_n[idxOrder] * (2.0*idxOrder + 1.0) * pow(cabs(b_n_truncated[idxBand*(order_truncated+1) + idxOrder]), 2.0);

    /* inverse ratio is filter gain */
    for (int idxBand=0; idxBand < nBands; idxBand++) {
        p_target[idxBand] = 1.0/(4.0*SAF_PI) * sqrt(p_target[idxBand]);
        p_truncated[idxBand] = 1.0/(4.0*SAF_PI) * sqrt(p_truncated[idxBand]);
        gain[idxBand] = (float)(p_target[idxBand] / (p_truncated[idxBand]+2.23e-13));
    }

    /* soft clip to limit maximum gain */
    float clipFactor = powf(10.0f, softThreshold/20.0f);
    for (int idxBand=0; idxBand<nBands; idxBand++){
        gain[idxBand] = gain[idxBand]/clipFactor;  /* norm by threshold */
        if (gain[idxBand] > 1.0f)
            gain[idxBand] = 1.0f + tanhf(gain[idxBand] - 1.0f);  /* soft clip to 2 */
        gain[idxBand] = gain[idxBand] * clipFactor;  /* de-norm */
    }

    /* clean-up */
    free(b_n_target);
    free(b_n_truncated);
    free(p_target);
    free(p_truncated);
}

void getLoudspeakerDecoderMtx
(
    float* ls_dirs_deg,
    int nLS,
    LOUDSPEAKER_AMBI_DECODER_METHODS method,
    int order,
    int enableMaxReWeighting,
    float* decMtx
)
{
    int i, j, nSH;
    float scale;
    float* Y_ls, *a_n, *decMtx_maxrE;

    nSH = ORDER2NSH(order);
    scale = 1.0f/SQRT4PI;
 
    switch(method){
        default:
        case LOUDSPEAKER_DECODER_DEFAULT:
        case LOUDSPEAKER_DECODER_SAD:
            /* Sampling Ambisonic Decoder (SAD) is simply the loudspeaker
             * spherical harmonic matrix scaled by the number of loudspeakers.
             */
            Y_ls = malloc1d(nSH*nLS*sizeof(float));
            getRSH(order, ls_dirs_deg, nLS, Y_ls);
            cblas_sscal(nLS*nSH, scale, Y_ls, 1);
            for(i=0; i<nLS; i++)
                for(j=0; j<nSH; j++)
                    decMtx[i*nSH+j] = (4.0f*SAF_PI) * Y_ls[j*nLS + i]/(float)nLS;
            free(Y_ls);
            break;
           
        case LOUDSPEAKER_DECODER_MMD:
            /* Mode-Matching Decoder (MMD) is simply the psuedo inverse of the
             * loudspeaker spherical harmonic matrix. */
            Y_ls = malloc1d(nSH*nLS*sizeof(float));
            getRSH(order, ls_dirs_deg, nLS, Y_ls);
            cblas_sscal(nLS*nSH, scale, Y_ls, 1);
            utility_spinv(NULL, Y_ls, nSH, nLS, decMtx);
            free(Y_ls);
            break;
            
        case LOUDSPEAKER_DECODER_EPAD:
            getEPAD(order, ls_dirs_deg, nLS, decMtx);
            break;
            
        case LOUDSPEAKER_DECODER_ALLRAD:
            getAllRAD(order, ls_dirs_deg, nLS, decMtx);
            break;
    }
    
    /* Apply maxRE weighting */
    if(enableMaxReWeighting){
        a_n = malloc1d(nSH*nSH*sizeof(float));
        getMaxREweights(order, 1, a_n); /* 1: weights returned as diagonal matrix */
        decMtx_maxrE = malloc1d(nLS * nSH * sizeof(float));
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nLS, nSH, nSH, 1.0f,
                    decMtx, nSH,
                    a_n, nSH, 0.0f,
                    decMtx_maxrE, nSH);
        memcpy(decMtx, decMtx_maxrE, nLS*nSH*sizeof(float));
        
        free(a_n);
        free(decMtx_maxrE);
    }
}

void getBinauralAmbiDecoderMtx
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
    int enableDiffCovMatching,
    int enableMaxReWeighting,
    float_complex* decMtx
)
{
    int i, k, nSH;
    float *tmp;
    float_complex* a_n, *decMtx_rE;
    const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);
    
    nSH = ORDER2NSH(order);
    
    switch(method){
        default:
        case BINAURAL_DECODER_DEFAULT:
        case BINAURAL_DECODER_LS:       getBinDecoder_LS(hrtfs, hrtf_dirs_deg, N_dirs, N_bands, order, weights, decMtx); break;
        case BINAURAL_DECODER_LSDIFFEQ: getBinDecoder_LSDIFFEQ(hrtfs, hrtf_dirs_deg, N_dirs, N_bands, order, weights, decMtx); break;
        case BINAURAL_DECODER_SPR:      getBinDecoder_SPR(hrtfs, hrtf_dirs_deg, N_dirs, N_bands, order, weights, decMtx); break;
        case BINAURAL_DECODER_TA:       getBinDecoder_TA(hrtfs, hrtf_dirs_deg, N_dirs, N_bands, order, freqVector, itd_s, weights, decMtx); break;
        case BINAURAL_DECODER_MAGLS:    getBinDecoder_MAGLS(hrtfs, hrtf_dirs_deg, N_dirs, N_bands, order, freqVector, weights, decMtx); break;
    }
    
    /* apply Max RE weighting per bin */
    if(enableMaxReWeighting){
        tmp = malloc1d(nSH*nSH*sizeof(float));
        a_n = malloc1d(nSH*nSH*sizeof(float_complex));
        decMtx_rE = malloc1d(NUM_EARS*nSH*sizeof(float_complex));
        getMaxREweights(order, 1, tmp);
        for(i=0; i<nSH*nSH; i++)
            a_n[i] = cmplxf(tmp[i], 0.0f);
        for(k=0; k<N_bands; k++){
            cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, NUM_EARS, nSH, nSH, &calpha,
                        &decMtx[k*NUM_EARS*nSH], nSH,
                        a_n, nSH, &cbeta,
                        decMtx_rE, nSH);
            memcpy(&decMtx[k*NUM_EARS*nSH], decMtx_rE, NUM_EARS*nSH*sizeof(float_complex));
        }
        free(tmp);
        free(a_n);
        free(decMtx_rE);
    }
    
    /* apply diffuse-field coherence matching per bin */
    if(enableDiffCovMatching)
        applyDiffCovMatching(hrtfs, hrtf_dirs_deg, N_dirs, N_bands, order, weights, decMtx);
}

void getBinauralAmbiDecoderFilters
(
    float_complex* hrtfs,
    float* hrtf_dirs_deg,
    int N_dirs,
    int fftSize,
    float fs,
    BINAURAL_AMBI_DECODER_METHODS method,
    int order,
    float* itd_s,
    float* weights,
    int enableDiffCovMatching,
    int enableMaxReWeighting,
    float* decFilters
)
{
    int i, j, k, nBins, nSH;
    float* freqVector;
    float_complex* decMtx, *decMtx_bins;
    void* hSafFFT;
    
    /* frequency-vector */
    nBins = fftSize/2 + 1;
    freqVector = malloc1d(nBins*sizeof(float));
    getUniformFreqVector(fftSize, fs, freqVector);
    
    /* compute decoding matrix per bin */
    nSH = ORDER2NSH(order);
    decMtx = malloc1d(nBins*NUM_EARS*nSH*sizeof(float_complex));
    getBinauralAmbiDecoderMtx(hrtfs, hrtf_dirs_deg, N_dirs, nBins, method,
                              order, freqVector, itd_s, weights, enableDiffCovMatching,
                              enableMaxReWeighting, decMtx);
    
    /* ifft, to obtain time-domain filters */
    decMtx_bins = malloc1d(nBins*sizeof(float_complex));
    saf_rfft_create(&hSafFFT, fftSize);
    for(i=0; i<NUM_EARS; i++){
        for(j=0; j<nSH; j++){
            for(k=0; k<nBins; k++)
                decMtx_bins[k] = decMtx[k*NUM_EARS*nSH + i*nSH + j];
            saf_rfft_backward(hSafFFT, decMtx_bins, &decFilters[i*nSH*fftSize + j*fftSize]);
        }
    }
    
    saf_rfft_destroy(&hSafFFT);
    free(freqVector);
    free(decMtx);
    free(decMtx_bins);
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
    float_complex C_ref[NUM_EARS][NUM_EARS], C_ambi[NUM_EARS][NUM_EARS];
    float_complex X[NUM_EARS][NUM_EARS], X_ambi[NUM_EARS][NUM_EARS];
    float_complex XH_Xambi[NUM_EARS][NUM_EARS], U[NUM_EARS][NUM_EARS];
    float_complex V[NUM_EARS][NUM_EARS], UX[NUM_EARS][NUM_EARS];
    float_complex VUX[NUM_EARS][NUM_EARS], M[NUM_EARS][NUM_EARS];
    const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);
    
    nSH = ORDER2NSH(order);
    
    /* integration weights */
    W = calloc1d(N_dirs*N_dirs, sizeof(float_complex));
    if(weights!=NULL)
        for(i=0; i<N_dirs; i++)
            W[i*N_dirs+i] = cmplxf(weights[i], 0.0f);
    else
        for(i=0; i<N_dirs; i++)
            W[i*N_dirs+i] = cmplxf(1.0f/(float)N_dirs, 0.0f);
    
    /* SH */
    Y_tmp = malloc1d(nSH*N_dirs*sizeof(float));
    Y_na = malloc1d(nSH*N_dirs*sizeof(float_complex));
    getRSH(order, hrtf_dirs_deg, N_dirs, Y_tmp);
    for(i=0; i<nSH*N_dirs; i++)
        Y_na[i] = cmplxf(Y_tmp[i], 0.0f);
    free(Y_tmp);
    
    /* apply diffuse-field coherence matching per band */
    H_W = malloc1d(NUM_EARS*N_dirs*sizeof(float_complex));
    H_ambi = malloc1d(NUM_EARS*N_dirs*sizeof(float_complex));
    decMtx_diffMatched = malloc1d(NUM_EARS*nSH*sizeof(float_complex));
    for(band=0; band<N_bands-1 /* skip Nyquist */; band++){
        /* Diffuse-field responses */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, NUM_EARS, N_dirs, N_dirs, &calpha,
                    &hrtfs[band*NUM_EARS*N_dirs], N_dirs,
                    W, N_dirs, &cbeta,
                    H_W, N_dirs);
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, NUM_EARS, NUM_EARS, N_dirs, &calpha,
                    H_W, N_dirs,
                    &hrtfs[band*NUM_EARS*N_dirs], N_dirs, &cbeta,
                    (float_complex*)C_ref, NUM_EARS);
        for(i=0; i<NUM_EARS; i++)
            C_ref[i][i] = cmplxf(crealf(C_ref[i][i]), 0.0f); /* force diagonal to be real */
        utility_cchol(NULL, (float_complex*)C_ref, NUM_EARS, (float_complex*)X);
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, NUM_EARS, N_dirs, nSH, &calpha,
                    &decMtx[band*NUM_EARS*nSH], nSH,
                    Y_na, N_dirs, &cbeta,
                    H_ambi, N_dirs);
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, NUM_EARS, N_dirs, N_dirs, &calpha,
                    H_ambi, N_dirs,
                    W, N_dirs, &cbeta,
                    H_W, N_dirs);
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, NUM_EARS, NUM_EARS, N_dirs, &calpha,
                    H_W, N_dirs,
                    H_ambi, N_dirs, &cbeta,
                    (float_complex*)C_ambi, NUM_EARS);
        for(i=0; i<NUM_EARS; i++)
            C_ambi[i][i] = cmplxf(crealf(C_ambi[i][i]), 0.0f); /* force diagonal to be real */
        utility_cchol(NULL, (float_complex*)C_ambi, NUM_EARS, (float_complex*)X_ambi);
        
        /* SVD */
        cblas_cgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans, NUM_EARS, NUM_EARS, NUM_EARS, &calpha,
                    (float_complex*)X_ambi, NUM_EARS,
                    (float_complex*)X, NUM_EARS, &cbeta,
                    (float_complex*)XH_Xambi, NUM_EARS);
        utility_csvd(NULL, (float_complex*)XH_Xambi, NUM_EARS, NUM_EARS, (float_complex*)U, NULL, (float_complex*)V, NULL);
        
        /* apply matching */
        cblas_cgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans, NUM_EARS, NUM_EARS, NUM_EARS, &calpha,
                    (float_complex*)U, NUM_EARS,
                    (float_complex*)X, NUM_EARS, &cbeta,
                    (float_complex*)UX, NUM_EARS);
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, NUM_EARS, NUM_EARS, NUM_EARS, &calpha,
                    (float_complex*)V, NUM_EARS,
                    (float_complex*)UX, NUM_EARS, &cbeta,
                    (float_complex*)VUX, NUM_EARS);
        utility_cglslv(NULL, (float_complex*)X_ambi, NUM_EARS, (float_complex*)VUX, NUM_EARS, (float_complex*)M);
        cblas_cgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans, NUM_EARS, nSH, NUM_EARS, &calpha,
                    (float_complex*)M, NUM_EARS,
                    &decMtx[band*NUM_EARS*nSH], nSH, &cbeta,
                    decMtx_diffMatched, nSH);
        memcpy(&decMtx[band*NUM_EARS*nSH], decMtx_diffMatched, NUM_EARS*nSH*sizeof(float_complex));
    }
    
    free(W);
    free(Y_na);
    free(H_W);
    free(H_ambi);
    free(decMtx_diffMatched);
} 
