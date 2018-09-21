/*
 Copyright 2017-2018 Leo McCormack
 
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
 *     saf_hrir.c
 * Description:
 *     A collection of head-related impulse-response (HRIR)- related functions.
 * Dependencies:
 *     saf_utilities, afSTFTlib
 * Author, date created:
 *     Leo McCormack, 12.12.2016
 */

#include "saf_hrir.h"
#include "saf_hrir_internal.h"

static inline float matlab_fmodf(float x, float y) {
    float tmp = fmodf(x, y);
    return tmp >= 0 ? tmp : tmp + y;
}

void estimateITDs
(
    float* hrirs /*N_dirs x 2 x hrir_len*/,
    int N_dirs,
    int hrir_len,
    int fs,
    float** itds_s /* & */
)
{
    int i, n, j, k, maxIdx, xcorr_len;
    float maxVal, itd_bounds, fc, Q, K, KK, D, wn, Wz1[2], Wz2[2], b[3], a[3];
    float* xcorr_LR, *ir_L, *ir_R, *hrir_lpf;

    /* calculate LPF coefficients, 2nd order IIR design equations from DAFX (2nd ed) p50 */
    fc = 750.0f;
    Q = 0.7071f;
    K = tanf(M_PI * fc/(float)fs);
    KK = K * K; 
    D = KK * Q + K + Q;
	b[0] = (KK * Q) / D; b[1] = (2.0f * KK * Q) / D; b[2] = (KK * Q) / D;
	a[0] = 1.0f; a[1] = (2.0f * Q * (KK - 1.0f)) / D; a[2] = (KK * Q - K + Q) / D;
    
    /* determine the ITD via the cross-correlation between the LPF'd left and right HRIR signals */
    xcorr_len = 2*(hrir_len)-1;
    itd_bounds = sqrtf(2.0f)/2e3f;
    if((*itds_s)!=NULL)
        free((*itds_s));
    (*itds_s) = malloc(N_dirs*sizeof(float));
    xcorr_LR = (float*)malloc(xcorr_len*sizeof(float));
    ir_L = (float*)malloc(hrir_len*sizeof(float));
    ir_R = (float*)malloc(hrir_len*sizeof(float));
    hrir_lpf = (float*)malloc(hrir_len*2*sizeof(float));
    for(i=0; i<N_dirs; i++){
        /* apply lpf */
        memset(Wz1, 0, 2*sizeof(float));
        memset(Wz2, 0, 2*sizeof(float));
        for (n=0; n<hrir_len; n++){
            for(j=0; j<2; j++){
                /* biquad difference equation (Direct form 2) */
                wn = hrirs[i*2*hrir_len + j*hrir_len + n] - a[1] * Wz1[j] - a[2] * Wz2[j];
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
        (*itds_s)[i] = ((float)hrir_len-(float)maxIdx-1.0f)/(float)fs;
        (*itds_s)[i] = (*itds_s)[i]>itd_bounds  ? itd_bounds  : (*itds_s)[i];
        (*itds_s)[i] = (*itds_s)[i]<-itd_bounds ? -itd_bounds : (*itds_s)[i];
    }
    
    free(xcorr_LR);
    free(ir_L);
    free(ir_R);
    free(hrir_lpf);
}

void estimateIPDmanipCurve
(
    float* itds_s,
    int N_dirs,
    float* centreFreq,
    int N_bands,
    float c,
    float maxVal,
    float* phi_bands
)
{
    int i;
    float ITD_max, f1, steepness;
    
    ITD_max = FLT_MIN;
    for(i=0; i<N_dirs; i++)
        ITD_max = fabsf(itds_s[i]) > ITD_max ? fabsf(itds_s[i]) : ITD_max;
    f1 = 1.0f/ITD_max;
    steepness = 1.7f;
    for(i=0; i<N_bands; i++)
       phi_bands[i] = MIN(powf(f1,steepness)/(powf(centreFreq[i]+2.23e-9f,steepness)), maxVal);
}

/* A C implementation of a MatLab function by Archontis Politis; published with permission */
void HRIRs2FilterbankHRTFs
(
    float* hrirs, /*N_bands x 2 x hrir_len*/
    int N_dirs,
    int hrir_len,
    float* itds_s,
    float* centreFreq,
    int N_bands,
    int enablePhaseManipFLAG,
    float_complex** hrtf_fb /* &, N_bands x 2 x N_dirs */
)
{
    int i, j, nd, band;
    float* ipd, *hrtf_diff, *phi_bands;

    /* convert the HRIRs to filterbank coefficients */
    FIRtoFilterbankCoeffs(hrirs, N_dirs, NUM_EARS, hrir_len, N_bands, hrtf_fb);
#if 1
    /* estimate phase manipulation curve */
    if (enablePhaseManipFLAG) {
        phi_bands = malloc(N_bands*sizeof(float));
        estimateIPDmanipCurve(itds_s, N_dirs, centreFreq, N_bands, 343.0f, 1.15f, phi_bands);
    }
    
    /* convert ITDs to phase differences -pi..pi */
    ipd = malloc(N_bands*N_dirs*sizeof(float));
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, N_bands, N_dirs, 1, 1.0,
                centreFreq, 1,
                itds_s, 1, 0.0,
                ipd, N_dirs);
    for(i=0; i<N_bands; i++){
        for(j=0; j<N_dirs; j++){
            if (enablePhaseManipFLAG)
                ipd[i*N_dirs+j] = phi_bands[i]*(matlab_fmodf(2.0f*M_PI*ipd[i*N_dirs+j] + M_PI, 2.0f*M_PI) - M_PI)/2.0f; /* /2 here, not later */
            else
                ipd[i*N_dirs+j] = (matlab_fmodf(2.0f*M_PI*ipd[i*N_dirs+j] + M_PI, 2.0f*M_PI) - M_PI)/2.0f; /* /2 here, not later */
        }
    }
    
    /* diffuse-field equalise */
    hrtf_diff = calloc(N_bands*NUM_EARS, sizeof(float));
    for(band=0; band<N_bands; band++)
        for(i=0; i<NUM_EARS; i++)
            for(j=0; j<N_dirs; j++)
                hrtf_diff[band*NUM_EARS + i] += powf(cabsf((*hrtf_fb)[band*NUM_EARS*N_dirs + i*N_dirs + j]), 2.0f);
    for(band=0; band<N_bands; band++)
        for(i=0; i<NUM_EARS; i++)
            hrtf_diff[band*NUM_EARS + i] = sqrtf(hrtf_diff[band*NUM_EARS + i]/(float)N_dirs);
    for(band=0; band<N_bands; band++)
        for(i=0; i<NUM_EARS; i++)
            for(nd=0; nd<N_dirs; nd++)
                (*hrtf_fb)[band*NUM_EARS*N_dirs + i*N_dirs + nd] = ccdivf((*hrtf_fb)[band*NUM_EARS*N_dirs + i*N_dirs + nd], cmplxf(hrtf_diff[band*NUM_EARS + i], 0.0f));
    
    /* create complex HRTFs by introducing the interaural phase differences (IPDs) to the HRTF magnitude responses */
    for(band=0; band<N_bands; band++){
        for(nd=0; nd<N_dirs; nd++){
            (*hrtf_fb)[band*NUM_EARS*N_dirs + 0*N_dirs + nd] = crmulf( cexpf(cmplxf(0.0f, ipd[band*N_dirs + nd])), cabsf((*hrtf_fb)[band*NUM_EARS*N_dirs + 0*N_dirs + nd]) );
            (*hrtf_fb)[band*NUM_EARS*N_dirs + 1*N_dirs + nd] = crmulf( cexpf(cmplxf(0.0f,-ipd[band*N_dirs + nd])), cabsf((*hrtf_fb)[band*NUM_EARS*N_dirs + 1*N_dirs + nd]) );
        }
    }

    free(ipd);
    free(hrtf_diff);
    if (enablePhaseManipFLAG)
        free(phi_bands);
#endif
}

/* A C implementation of a MatLab function by Archontis Politis; published with permission */
void interpFilterbankHRTFs
(
    float_complex* hrtfs, /* N_bands x 2 x N_hrtf_dirs */
    float* itds,
    float* freqVector,
    float* vbap_gtable, 
    int N_hrtf_dirs,
    int N_bands,
    int N_interp_dirs,
    int enablePhaseManipFLAG,      /* 0: off, 1: on */
    float_complex* hrtfs_interp /* pre-alloc, N_bands x 2 x N_interp_dirs */
)
{
    int i, band;
    float* itd_interp, *mags_interp, *ipd_interp, *phi_bands;
    float** mags;
    
    mags = (float**)malloc(N_bands*sizeof(float*));
    itd_interp = malloc(N_interp_dirs*sizeof(float));
    mags_interp = malloc(N_interp_dirs*NUM_EARS*sizeof(float));
    ipd_interp = malloc(N_interp_dirs*sizeof(float));
    
    /* estimate phase manipulation curve */
    if (enablePhaseManipFLAG) {
        phi_bands = malloc(N_bands*sizeof(float));
        estimateIPDmanipCurve(itds, N_hrtf_dirs, freqVector, N_bands, 343.0f, 1.15f, phi_bands);
    }
    
    /* calculate HRTF magnitudes */
    for(band=0; band<N_bands; band++){
        mags[band] = malloc(NUM_EARS * N_hrtf_dirs*sizeof(float));
        for(i=0; i< NUM_EARS * N_hrtf_dirs ; i++)
            mags[band][i] = cabsf(hrtfs[band*NUM_EARS * N_hrtf_dirs + i]);
    }
    
    /* interpolate ITDs */
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N_interp_dirs, 1, N_hrtf_dirs, 1.0f,
                vbap_gtable, N_hrtf_dirs,
                itds, 1, 0.0f,
                itd_interp, 1);
    for(band=0; band<N_bands; band++){
        /* interpolate HRTF magnitudes */
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, N_interp_dirs, NUM_EARS, N_hrtf_dirs, 1.0f,
                    vbap_gtable, N_hrtf_dirs,
                    mags[band], N_hrtf_dirs, 0.0f,
                    mags_interp, NUM_EARS);
        
        /* convert ITDs to phase differences -pi..pi */
        for(i=0; i<N_interp_dirs; i++){
            if (enablePhaseManipFLAG)
                ipd_interp[i] = phi_bands[i]*(matlab_fmodf(2.0f*M_PI*freqVector[band]*itd_interp[i] + M_PI, 2.0f*M_PI) - M_PI)/2.0f; /* /2 here, not later */
            else
                ipd_interp[i] = (matlab_fmodf(2.0f*M_PI*freqVector[band]*itd_interp[i] + M_PI, 2.0f*M_PI) - M_PI)/2.0f; /* /2 here, not later */
        }
        
        /* reintroduce the interaural phase differences (IPD) */
        for(i=0; i<N_interp_dirs; i++){
            hrtfs_interp[band*NUM_EARS*N_interp_dirs + 0* N_interp_dirs +i] = ccmulf( cmplxf(mags_interp[i*NUM_EARS+0],0.0f), cexpf(cmplxf(0.0f, ipd_interp[i])) );
            hrtfs_interp[band*NUM_EARS*N_interp_dirs + 1* N_interp_dirs +i] = ccmulf( cmplxf(mags_interp[i*NUM_EARS+1],0.0f), cexpf(cmplxf(0.0f,-ipd_interp[i])) );
        }
    }

    free(itd_interp);
    for(band=0; band<N_bands; band++)
        free(mags[band]);
    free(mags);
    free(mags_interp);
    free(ipd_interp);
    if (enablePhaseManipFLAG)
        free(phi_bands);
}











