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
 * @file saf_hrir.c
 * @ingroup HRIR
 * @brief Public source for the HRIR/HRTF processing module (#SAF_HRIR_MODULE)
 *
 * A collection functions for processing head-related impulse-response (HRIR).
 * Including: estimation of the interaural time differences (ITDs), conversion
 * of HRIRs to HRTFs or filterbank coefficients; diffuse-field equalisation and
 * phase simplication; and HRTF interpolation.
 *
 * @author Leo McCormack
 * @date 12.12.2016
 * @license ISC
 */
 
#include "saf_hrir.h"
#include "../saf_utilities/saf_utilities.h"
#include "saf_externals.h"

/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

void estimateITDs
(
    float* hrirs /* N_dirs x NUM_EARS x hrir_len */,
    int N_dirs,
    int hrir_len,
    int fs,
    float* itds_s
)
{
    int i, n, j, k, maxIdx, xcorr_len;
    float maxVal, itd_bounds, fc, Q, K, KK, D, wn, Wz1[2], Wz2[2], b[3], a[3];
    float* xcorr_LR, *ir_L, *ir_R, *hrir_lpf;

    /* calculate LPF coefficients, 2nd order IIR design equations from DAFX (2nd ed) p50 */
    fc = 750.0f;
    Q = 0.7071f;
    K = tanf(SAF_PI * fc/(float)fs);
    KK = K * K; 
    D = KK * Q + K + Q;
    b[0] = (KK * Q) / D; b[1] = (2.0f * KK * Q) / D; b[2] = (KK * Q) / D;
    a[0] = 1.0f; a[1] = (2.0f * Q * (KK - 1.0f)) / D; a[2] = (KK * Q - K + Q) / D;
    
    /* determine the ITD via the cross-correlation between the LPF'd left and right HRIR signals */
    xcorr_len = 2*(hrir_len)-1;
    itd_bounds = sqrtf(2.0f)/2e3f;
    xcorr_LR = (float*)malloc1d(xcorr_len*sizeof(float));
    ir_L = (float*)malloc1d(hrir_len*sizeof(float));
    ir_R = (float*)malloc1d(hrir_len*sizeof(float));
    hrir_lpf = (float*)malloc1d(hrir_len*2*sizeof(float));
    for(i=0; i<N_dirs; i++){
        /* apply lpf */
        memset(Wz1, 0, 2*sizeof(float));
        memset(Wz2, 0, 2*sizeof(float));
        for (n=0; n<hrir_len; n++){
            for(j=0; j<NUM_EARS; j++){
                /* biquad difference equation (Direct form 2) */
                wn = hrirs[i*NUM_EARS*hrir_len + j*hrir_len + n] - a[1] * Wz1[j] - a[2] * Wz2[j];
                hrir_lpf[n*2+j]  = b[0] * wn + b[1]*Wz1[j] + b[2]*Wz2[j];
                
                /* shuffle delays */
                Wz2[j] = Wz1[j];
                Wz1[j] = wn;
            }
        }
        
        /* xcorr between L and R */
        for(k=0; k<hrir_len; k++){
            ir_L[k] = hrir_lpf[k*2+0];
            ir_R[k] = hrir_lpf[k*2+1];
        }
        cxcorr(ir_L, ir_R, xcorr_LR, hrir_len, hrir_len);
        maxIdx = 0;
        maxVal = 0.0f;
        for(j=0; j<xcorr_len; j++){
            if(xcorr_LR[j] > maxVal){
                maxIdx = j;
                maxVal = xcorr_LR[j];
            }
        }
        itds_s[i] = ((float)hrir_len-(float)maxIdx-1.0f)/(float)fs;
        itds_s[i] = itds_s[i] >  itd_bounds ?  itd_bounds : itds_s[i];
        itds_s[i] = itds_s[i] < -itd_bounds ? -itd_bounds : itds_s[i];
    }
    
    free(xcorr_LR);
    free(ir_L);
    free(ir_R);
    free(hrir_lpf);
}

void HRIRs2HRTFs_afSTFT
(
    float* hrirs, /* N_dirs x NUM_EARS x hrir_len */
    int N_dirs,
    int hrir_len,
    int hopsize,
    int LDmode,
    int hybridmode,
    float_complex* hrtf_fb /* nBands x NUM_EARS x N_dirs */
)
{
    /* convert the HRIRs to filterbank coefficients */
    afSTFT_FIRtoFilterbankCoeffs(hrirs, N_dirs, NUM_EARS, hrir_len, hopsize, LDmode, hybridmode, hrtf_fb);
}

void HRIRs2HRTFs_qmf
(
    float* hrirs, /* N_dirs x NUM_EARS x hrir_len */
    int N_dirs,
    int hrir_len,
    int hopsize,
    int hybridmode,
    float_complex* hrtf_fb /* nBands x NUM_EARS x N_dirs */
)
{
    /* convert the HRIRs to filterbank coefficients */
    qmf_FIRtoFilterbankCoeffs(hrirs, N_dirs, NUM_EARS, hrir_len, hopsize, hybridmode, hrtf_fb);
}

void HRIRs2HRTFs
(
    float* hrirs, /* N_dirs x NUM_EARS x hrir_len */
    int N_dirs,
    int hrir_len,
    int fftSize,
    float_complex* hrtfs /* (fftSize/2+1) x NUM_EARS x N_dirs  */
)
{
    int i, j, k, nBins;
    void* hSafFFT;
    float* hrir_pad;
    float_complex* hrtf;

    //TODO: if fftSize is shorter than hrir_len, maybe truncate based on the median peak?
    /* Perform FFT */
    nBins = fftSize/2 + 1;
    saf_rfft_create(&hSafFFT, fftSize);
    hrir_pad = calloc1d(fftSize, sizeof(float));
    hrtf = malloc1d(nBins*sizeof(float_complex));
    for(i=0; i<N_dirs; i++){
        for(j=0; j<NUM_EARS; j++){
            memcpy(hrir_pad, &hrirs[i*NUM_EARS*hrir_len+j*hrir_len], SAF_MIN(fftSize, hrir_len)*sizeof(float));
            saf_rfft_forward(hSafFFT, hrir_pad, hrtf);
            for(k=0; k<nBins; k++)
                hrtfs[k*NUM_EARS*N_dirs + j*N_dirs + i] = hrtf[k];
         }
    }
    
    saf_rfft_destroy(&hSafFFT);
    free(hrir_pad);
    free(hrtf);
}

void diffuseFieldEqualiseHRTFs
(
    int N_dirs,
    float* itds_s,
    float* centreFreq,
    int N_bands,
    float* weights,
    int applyEQ,
    int applyPhase,
    float_complex* hrtfs   /* N_bands x #NUM_EARS x N_dirs */
)
{
    /* Anything to do at all? */
    if(applyEQ + applyPhase)
    {
        int i, j, nd, band;
        float* ipd, *hrtf_diff, *_weights;
      
        /* diffuse-field equalise */
        if(applyEQ){
            hrtf_diff = calloc1d(N_bands*NUM_EARS, sizeof(float));
            if(weights == NULL){
                _weights = malloc1d(N_dirs*sizeof(float));
                for(int idx=0; idx < N_dirs; idx++)
                    _weights[idx] = 4.f*SAF_PI / (float)N_dirs;
            }
            else
                _weights = weights;
            for(band=0; band<N_bands; band++)
                for(i=0; i<NUM_EARS; i++)
                    for(j=0; j<N_dirs; j++)
                        hrtf_diff[band*NUM_EARS + i] += _weights[j]/(4.f*SAF_PI) * powf(cabsf(hrtfs[band*NUM_EARS*N_dirs + i*N_dirs + j]), 2.0f);
            for(band=0; band<N_bands; band++)
                for(i=0; i<NUM_EARS; i++)
                    hrtf_diff[band*NUM_EARS + i] = sqrtf(hrtf_diff[band*NUM_EARS + i]);
            for(band=0; band<N_bands; band++)
                for(i=0; i<NUM_EARS; i++)
                    for(nd=0; nd<N_dirs; nd++)
                        hrtfs[band*NUM_EARS*N_dirs + i*N_dirs + nd] = ccdivf(hrtfs[band*NUM_EARS*N_dirs + i*N_dirs + nd], cmplxf(hrtf_diff[band*NUM_EARS + i] + 2.23e-8f, 0.0f));
            free(hrtf_diff);
            if(weights==NULL)
                free(_weights);
        }
        
        /* Create complex HRTFs by introducing the interaural phase differences
         * (IPDs) to the HRTF magnitude responses */
        if(applyPhase){
            /* convert ITDs to phase differences -pi..pi */
            ipd = malloc1d(N_bands*N_dirs*sizeof(float));
            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, N_bands, N_dirs, 1, 1.0,
                        centreFreq, 1,
                        itds_s, 1, 0.0,
                        ipd, N_dirs);
            for(i=0; i<N_bands; i++)
                for(j=0; j<N_dirs; j++)
                    ipd[i*N_dirs+j] = (matlab_fmodf(2.0f*SAF_PI*ipd[i*N_dirs+j] + SAF_PI, 2.0f*SAF_PI) - SAF_PI)/2.0f; /* /2 here, not later */
  
            for(band=0; band<N_bands; band++){
                for(nd=0; nd<N_dirs; nd++){
                    hrtfs[band*NUM_EARS*N_dirs + 0*N_dirs + nd] = crmulf( cexpf(cmplxf(0.0f, ipd[band*N_dirs + nd])), cabsf(hrtfs[band*NUM_EARS*N_dirs + 0*N_dirs + nd]) );
                    hrtfs[band*NUM_EARS*N_dirs + 1*N_dirs + nd] = crmulf( cexpf(cmplxf(0.0f,-ipd[band*N_dirs + nd])), cabsf(hrtfs[band*NUM_EARS*N_dirs + 1*N_dirs + nd]) );
                }
            }
            free(ipd);
        }
    }
}
 
void interpHRTFs
(
    float_complex* hrtfs, /* N_bands x 2 x N_hrtf_dirs */
    float* itds,
    float* freqVector,
    float* interp_table,
    int N_hrtf_dirs,
    int N_bands,
    int N_interp_dirs,
    float_complex* hrtfs_interp /* pre-alloc, N_bands x 2 x N_interp_dirs */
)
{
    int i, band;
    float* itd_interp, *mags_interp, *ipd_interp;
    float** mags;
    float_complex* interp_table_cmplx;
    const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);

    if(itds==NULL || freqVector==NULL){
        /* prep */
        interp_table_cmplx = calloc1d(N_interp_dirs*N_hrtf_dirs, sizeof(float_complex));
        cblas_scopy(N_interp_dirs*N_hrtf_dirs, interp_table, 1, (float*)interp_table_cmplx, 2);

        /* interpolate HRTF spectra */
        for(band=0; band<N_bands; band++){
            cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasTrans, NUM_EARS, N_interp_dirs, N_hrtf_dirs, &calpha,
                        &hrtfs[band*NUM_EARS*N_hrtf_dirs], N_hrtf_dirs,
                        interp_table_cmplx, N_hrtf_dirs, &cbeta,
                        &hrtfs_interp[band*NUM_EARS*N_interp_dirs], N_interp_dirs);
        }

        /* clean-up */
        free(interp_table_cmplx);
    }
    else{
        /* prep */
        mags = (float**)malloc1d(N_bands*sizeof(float*));
        itd_interp = malloc1d(N_interp_dirs*sizeof(float));
        mags_interp = malloc1d(N_interp_dirs*NUM_EARS*sizeof(float));
        ipd_interp = malloc1d(N_interp_dirs*sizeof(float));

        /* calculate HRTF magnitudes */
        for(band=0; band<N_bands; band++){
            mags[band] = malloc1d(NUM_EARS * N_hrtf_dirs*sizeof(float));
            for(i=0; i< NUM_EARS * N_hrtf_dirs ; i++)
                mags[band][i] = cabsf(hrtfs[band*NUM_EARS * N_hrtf_dirs + i]);
        }

        /* interpolate ITDs */
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N_interp_dirs, 1, N_hrtf_dirs, 1.0f,
                    interp_table, N_hrtf_dirs,
                    itds, 1, 0.0f,
                    itd_interp, 1);
        for(band=0; band<N_bands; band++){
            /* interpolate HRTF magnitudes */
            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, N_interp_dirs, NUM_EARS, N_hrtf_dirs, 1.0f,
                        interp_table, N_hrtf_dirs,
                        mags[band], N_hrtf_dirs, 0.0f,
                        mags_interp, NUM_EARS);

            /* convert ITDs to phase differences -pi..pi */
            for(i=0; i<N_interp_dirs; i++)
                ipd_interp[i] = (matlab_fmodf(2.0f*SAF_PI*freqVector[band]*itd_interp[i] + SAF_PI, 2.0f*SAF_PI) - SAF_PI)/2.0f; /* /2 here, not later */

            /* reintroduce the interaural phase differences (IPD) */
            for(i=0; i<N_interp_dirs; i++){
                hrtfs_interp[band*NUM_EARS*N_interp_dirs + 0* N_interp_dirs +i] = ccmulf( cmplxf(mags_interp[i*NUM_EARS+0],0.0f), cexpf(cmplxf(0.0f, ipd_interp[i])) );
                hrtfs_interp[band*NUM_EARS*N_interp_dirs + 1* N_interp_dirs +i] = ccmulf( cmplxf(mags_interp[i*NUM_EARS+1],0.0f), cexpf(cmplxf(0.0f,-ipd_interp[i])) );
            }
        }

        /* clean-up */
        free(itd_interp);
        for(band=0; band<N_bands; band++)
            free(mags[band]);
        free(mags);
        free(mags_interp);
        free(ipd_interp);
    }
}

void binauralDiffuseCoherence
(
    float_complex* hrtfs, /* N_bands x 2 x N_hrtf_dirs */
    float* itds,
    float* freqVector,
    int N_hrtf_dirs,
    int N_bands,
    float* HRTFcoh
)
{
    int i, j;
    float* ipd;
    float_complex *hrtf_ipd_lr;
    
    /* convert ITDs to phase differences -pi..pi */
    ipd = malloc1d(N_bands*N_hrtf_dirs*sizeof(float));
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, N_bands, N_hrtf_dirs, 1, 1.0,
                freqVector, 1,
                itds, 1, 0.0,
                ipd, N_hrtf_dirs);
    for(i=0; i<N_bands; i++)
        for(j=0; j<N_hrtf_dirs; j++)
            ipd[i*N_hrtf_dirs+j] = (matlab_fmodf(2.0f*SAF_PI*ipd[i*N_hrtf_dirs+j] + SAF_PI, 2.0f*SAF_PI) - SAF_PI);
    
    /* compute complex coherence */
    hrtf_ipd_lr = calloc1d(N_bands, sizeof(float_complex));
    for(i=0; i<N_bands; i++){
        for(j=0; j<N_hrtf_dirs; j++)
            hrtf_ipd_lr[i] = ccaddf(hrtf_ipd_lr[i], crmulf(cexpf(crmulf(cmplxf(0.0f, 1.0f), ipd[i*N_hrtf_dirs+j])),
                                                    cabsf(hrtfs[i*NUM_EARS*N_hrtf_dirs + 0*N_hrtf_dirs + j])*
                                                    cabsf(hrtfs[i*NUM_EARS*N_hrtf_dirs + 1*N_hrtf_dirs + j])));
        hrtf_ipd_lr[i] = ccdivf(hrtf_ipd_lr[i], cmplxf((float)N_hrtf_dirs, 0.0f));
    }
    
    /* due to almost axisymmetry of ITD, the coherence is almost real */
    for(i=0; i<N_bands; i++)
        HRTFcoh[i] = crealf(hrtf_ipd_lr[i]) < 0.0f ? 0.0f : crealf(hrtf_ipd_lr[i]);
    HRTFcoh[0] = 1.0f; /* force 1 at DC */
    
    free(ipd);
    free(hrtf_ipd_lr);
}

void resampleHRIRs
(
    float* hrirs_in,
    int hrirs_N_dirs,
    int hrirs_in_len,
    int hrirs_in_fs,
    int hrirs_out_fs,
    int padToNextPow2,
    float** hrirs_out,
    int* hrirs_out_len
)
{
    int ch, hrirs_out_ld;
    float resample_factor;
#if defined(SAF_USE_INTEL_IPP) && 0 /* works fine on macOS, but not on MSVC. Not tried Linux. Use SPEEX for now... */
    Ipp64f pTime;
    const int history = 128;
    float *inBuffer, *outBuffer;
    int filterLength, pSize, numFilters, outL;
    IppStatus error;
#else
    unsigned int in_length, out_length;
    int ERROR_VAL, out_latency, nsample_proc;
    float *zeros;
    SpeexResamplerState *pRS;
#endif

    /* New HRIR length */
    resample_factor = (float)hrirs_out_fs / (float)hrirs_in_fs;
    (*hrirs_out_len) = (int)ceilf((float)hrirs_in_len * resample_factor);
    hrirs_out_ld = padToNextPow2 ? (int)pow(2.0, ceil(log((double)(*hrirs_out_len))/log(2.0))) : (*hrirs_out_len);

#if defined(SAF_USE_INTEL_IPP) && 0 /* works fine on macOS, but not on MSVC. Not tried Linux. Use SPEEX for now... */
    /* Initialise IPP resampler */
    error = ippsResamplePolyphaseFixedGetSize_32f(hrirs_in_fs, hrirs_out_fs, 2*(history-1), &pSize, &filterLength, &numFilters, ippAlgHintFast);
    saf_assert(!error, "IPP error");
    IppsResamplingPolyphaseFixed_32f* spec;
    spec = (IppsResamplingPolyphaseFixed_32f*)ippsMalloc_8u(pSize);
    error = ippsResamplePolyphaseFixedInit_32f(hrirs_in_fs, hrirs_out_fs, 2*(history-1), 0.98f, 12.0f, spec, ippAlgHintFast);
    saf_assert(!error, "IPP error");
    inBuffer = ippsMalloc_32f(hrirs_in_len + history * 2 + 2);
    outBuffer = ippsMalloc_32f(hrirs_out_ld + 2);
    ippsZero_32f(inBuffer, hrirs_in_len + history * 2 + 2);

    /* Apply IPP resampler */
    (*hrirs_out) = calloc1d(hrirs_N_dirs*NUM_EARS*(hrirs_out_ld), sizeof(float));
    for(ch=0; ch<hrirs_N_dirs*NUM_EARS; ch++){
        pTime = history;
        outL =  0;

        /* Apply resampling */
        ippsCopy_32f(hrirs_in + ch * hrirs_in_len, inBuffer + history, hrirs_in_len);
        saf_assert(!ippsResamplePolyphaseFixed_32f(inBuffer, hrirs_in_len, outBuffer, 1.0f, &pTime, &outL, spec), "IPP error");
        saf_assert(hrirs_out_ld==outL, "Not all samples were processed!");
        ippsCopy_32f(outBuffer, (*hrirs_out) + ch * (hrirs_out_ld), hrirs_out_ld);
    }

    (*hrirs_out_len) = hrirs_out_ld;

    /* Clean-up */
    ippsFree(spec);
    ippsFree(inBuffer);
    ippsFree(outBuffer);

#else
    /* Initialise SPEEX resampler */
    pRS = speex__resampler_init(1 /*one channel at a time*/, hrirs_in_fs, hrirs_out_fs, SPEEX_RESAMPLER_QUALITY_MAX, &ERROR_VAL);
    out_latency = speex__resampler_get_output_latency(pRS);
    zeros = calloc1d(out_latency, sizeof(float));

    /* Apply SPEEX resampler */
    (*hrirs_out) = calloc1d(hrirs_N_dirs*NUM_EARS*(hrirs_out_ld), sizeof(float));
    for(ch=0; ch<hrirs_N_dirs*NUM_EARS; ch++){
        speex__resampler_reset_mem(pRS);
        speex__resampler_skip_zeros(pRS);
        nsample_proc = 0;

        /* Pass the FIR through the resampler */
        in_length = hrirs_in_len;
        out_length = hrirs_out_ld;
        ERROR_VAL = speex__resampler_process_float((pRS), 0, hrirs_in + ch * hrirs_in_len, &in_length,
                                                   (*hrirs_out) + ch * (hrirs_out_ld), &out_length);
        nsample_proc += out_length; /* Current number of output samples processed */

        /* Pass through zeros to get the tail of the filter too */
        while(nsample_proc<(hrirs_out_ld)){
            in_length = out_latency;
            out_length = (hrirs_out_ld)-nsample_proc;
            ERROR_VAL = speex__resampler_process_float((pRS), 0, zeros, &in_length,
                                                       (*hrirs_out) + ch * (hrirs_out_ld) + nsample_proc, &out_length);
            nsample_proc += out_length;
        }
        saf_assert(nsample_proc==(hrirs_out_ld), "Not all samples were processed!");
    }

    (*hrirs_out_len) = hrirs_out_ld;

    /* Clean-up */
    speex__resampler_destroy(pRS);
    free(zeros);
#endif
}
