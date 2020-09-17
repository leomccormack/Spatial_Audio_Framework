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
 * @file saf_hrir_internal.c
 * @ingroup HRIR
 * @brief Internal source for the HRIR/HRTF processing module (#SAF_HRIR_MODULE)
 *
 * A collection of head-related impulse-response (HRIR) functions. Including
 * estimation of the interaural time differences (ITDs), conversion of HRIRs to
 * HRTF filterbank coefficients, and HRTF interpolation utilising amplitude-
 * normalised VBAP gains.
 *
 * @author Leo McCormack
 * @date 12.12.2016
 */

#include "saf_hrir.h"
#include "saf_hrir_internal.h"

/**
 * Passes input time-domain data through the afSTFT filterbank.
 *
 * Copyright (c) 2015 Juha Vilkamo, MIT license
 */
static void afAnalyse
(
    float* inTD/* nSamplesTD x nCH */,
    int nSamplesTD,
    int nCH,
    int hopSize,
    int hybridmode,
    float_complex* outTF /* out_nBands x nTimeslots x nCH */
)
{
    int t, ch, sample, band;
    void* hSTFT;
    float_complex*** FrameTF;
    float** tempFrameTD;
    int nTimeSlots, nBands;

    nBands = hopSize + (hybridmode ? 5 : 1);
    nTimeSlots = nSamplesTD/hopSize;

    /* allocate memory */
    afSTFT_create(&(hSTFT), nCH, 1, hopSize, 0, hybridmode, AFSTFT_TIME_CH_BANDS);
    FrameTF = (float_complex***)malloc3d(nTimeSlots, nCH, nBands, sizeof(float_complex));
    tempFrameTD = (float**)malloc2d(nCH, nSamplesTD, sizeof(float));

    /* perform TF transform */
    for(ch=0; ch<nCH; ch++)
        for(sample=0; sample<nSamplesTD; sample++)
            tempFrameTD[ch][sample] = inTD[sample* nCH + ch];
    afSTFT_forward(hSTFT, tempFrameTD, nSamplesTD, FrameTF);

    /* save result to output */
    for (band = 0; band < nBands; band++)
        for (t = 0; t < nTimeSlots; t++)
            for (ch = 0; ch < nCH; ch++)
                outTF[band * nTimeSlots * nCH + t * nCH + ch] = FrameTF[t][ch][band];

    /* clean-up */
    afSTFT_destroy(&hSTFT); 
    free(FrameTF);
    free(tempFrameTD);
}

/**
 * Passes input time-domain data through the QMF filterbank.
 */
static void qmfAnalyse
(
    float* inTD/* nSamplesTD x nCH */,
    int nSamplesTD,
    int nCH,
    int hopSize,
    int hybridmode,
    float_complex* outTF /* out_nBands x nTimeslots x nCH */
)
{
    int t, ch, sample, band;
    void* hQMF;
    float_complex*** FrameTF;
    float** tempFrameTD;
    int nTimeSlots, nBands;

    nBands = hopSize + (hybridmode ? 7 : 0);
    nTimeSlots = nSamplesTD/hopSize;

    /* allocate memory */
    qmf_create(&(hQMF), nCH, 1, hopSize, hybridmode, QMF_TIME_CH_BANDS);
    FrameTF = (float_complex***)malloc3d(nTimeSlots, nCH, nBands, sizeof(float_complex));
    tempFrameTD = (float**)malloc2d(nCH, nSamplesTD, sizeof(float));

    /* perform TF transform */
    for(ch=0; ch<nCH; ch++)
        for(sample=0; sample<nSamplesTD; sample++)
            tempFrameTD[ch][sample] = inTD[sample* nCH + ch];
    qmf_analysis(hQMF, tempFrameTD, nSamplesTD, FrameTF);

    /* save result to output */
    for (band = 0; band < nBands; band++)
        for (t = 0; t < nTimeSlots; t++)
            for (ch = 0; ch < nCH; ch++)
                outTF[band * nTimeSlots * nCH + t * nCH + ch] = FrameTF[t][ch][band];

    /* clean-up */
    qmf_destroy(&hQMF);
    free(FrameTF);
    free(tempFrameTD);
}

void FIRtoFilterbankCoeffs_afSTFT
(
    float* hIR /*N_dirs x nCH x ir_len*/,
    int N_dirs,
    int nCH,
    int ir_len,
    int hopSize,
    int hybridmode,
    float_complex* hFB /* nBands x nCH x N_dirs */
)
{
    int i, j, t, nd, nm, nTimeSlots, ir_pad, nBands;
    int* maxIdx;
    float maxVal, idxDel, irFB_energy, irFB_gain, phase;
    float* centerImpulse, *centerImpulseFB_energy, *ir;
    float_complex cross;
    float_complex* centerImpulseFB, *irFB;

    nBands = hopSize + (hybridmode ? 5 : 1);
    ir_pad = 1024;//+512;
    nTimeSlots = (ir_len+ir_pad)/hopSize;
    maxIdx = calloc1d(nCH,sizeof(int));
    centerImpulse = calloc1d(ir_len+ir_pad, sizeof(float));
    
    /* pick a direction to estimate the center of FIR delays */
    for(j=0; j<nCH; j++){
        maxVal = 2.23e-13f;
        for(i=0; i<ir_len; i++){
            if(hIR[j*ir_len + i] > maxVal){
                maxVal = hIR[j*ir_len + i];
                maxIdx[j] = i;
            }
        }
    }
    idxDel = 0.0f;
    for(j=0; j<nCH; j++)
        idxDel += (float)maxIdx[j];
    idxDel /= (float)nCH;
    idxDel = (idxDel + 1.5f);
    
    /* ideal impulse at mean delay */
    centerImpulse[(int)idxDel] = 1.0f;
    
    /* analyse impulse with the filterbank */
    centerImpulseFB = malloc1d(nBands*nTimeSlots*nCH*sizeof(float_complex));
    afAnalyse(centerImpulse, ir_len+ir_pad, 1, hopSize, hybridmode, centerImpulseFB);
    centerImpulseFB_energy = calloc1d(nBands, sizeof(float));
    for(i=0; i<nBands; i++)
        for(t=0; t<nTimeSlots; t++)
            centerImpulseFB_energy[i] += powf(cabsf(centerImpulseFB[i*nTimeSlots + t]), 2.0f);
    
    /* initialise FB coefficients */
    ir = calloc1d( (ir_len+ir_pad) * nCH, sizeof(float));
    irFB = malloc1d(nBands*nCH*nTimeSlots*sizeof(float_complex));
    for(nd=0; nd<N_dirs; nd++){
        for(j=0; j<ir_len; j++)
            for(i=0; i<nCH; i++)
                ir[j*nCH+i] = hIR[nd*nCH*ir_len + i*ir_len + j];
        afAnalyse(ir, ir_len+ir_pad, nCH, hopSize, hybridmode, irFB);
        for(nm=0; nm<nCH; nm++){
            for(i=0; i<nBands; i++){
                irFB_energy = 0;
                for(t=0; t<nTimeSlots; t++)
                    irFB_energy += powf(cabsf(irFB[i*nTimeSlots*nCH + t*nCH + nm]), 2.0f); /* out_nBands x nTimeslots x nCH */
                irFB_gain = sqrtf(irFB_energy/MAX(centerImpulseFB_energy[i], 2.23e-8f));
                cross = cmplxf(0.0f,0.0f);
                for(t=0; t<nTimeSlots; t++)
                    cross = ccaddf(cross, ccmulf(irFB[i*nTimeSlots*nCH + t*nCH + nm], conjf(centerImpulseFB[i*nTimeSlots + t])));
                phase = atan2f(cimagf(cross), crealf(cross));
                hFB[i*nCH*N_dirs + nm*N_dirs + nd] = crmulf( cexpf(cmplxf(0.0f, phase)), irFB_gain);
            }
        }
    }
    
    /* clean-up */
    free(maxIdx);
    free(centerImpulse);
    free(centerImpulseFB_energy);
    free(centerImpulseFB);
    free(ir);
    free(irFB);
}

void FIRtoFilterbankCoeffs_qmf
(
    float* hIR /*N_dirs x nCH x ir_len*/,
    int N_dirs,
    int nCH,
    int ir_len,
    int hopSize,
    int hybridmode,
    float_complex* hFB /* nBands x nCH x N_dirs */
)
{
    int i, j, t, nd, nm, nTimeSlots, ir_pad, nBands;
    int* maxIdx;
    float maxVal, idxDel, irFB_energy, irFB_gain, phase;
    float* centerImpulse, *centerImpulseFB_energy, *ir;
    float_complex cross;
    float_complex* centerImpulseFB, *irFB;

    nBands = hopSize + (hybridmode ? 7 : 0);
    ir_pad = 1024;//+512;
    nTimeSlots = (ir_len+ir_pad)/hopSize;
    maxIdx = calloc1d(nCH,sizeof(int));
    centerImpulse = calloc1d(ir_len+ir_pad, sizeof(float));

    /* pick a direction to estimate the center of FIR delays */
    for(j=0; j<nCH; j++){
        maxVal = 2.23e-13f;
        for(i=0; i<ir_len; i++){
            if(hIR[j*ir_len + i] > maxVal){
                maxVal = hIR[j*ir_len + i];
                maxIdx[j] = i;
            }
        }
    }
    idxDel = 0.0f;
    for(j=0; j<nCH; j++)
        idxDel += (float)maxIdx[j];
    idxDel /= (float)nCH;
    idxDel = (idxDel + 1.5f);

    /* ideal impulse at mean delay */
    centerImpulse[(int)idxDel] = 1.0f;

    /* analyse impulse with the filterbank */
    centerImpulseFB = malloc1d(nBands*nTimeSlots*nCH*sizeof(float_complex));
    qmfAnalyse(centerImpulse, ir_len+ir_pad, 1, hopSize, hybridmode, centerImpulseFB);
    centerImpulseFB_energy = calloc1d(nBands, sizeof(float));
    for(i=0; i<nBands; i++)
        for(t=0; t<nTimeSlots; t++)
            centerImpulseFB_energy[i] += powf(cabsf(centerImpulseFB[i*nTimeSlots + t]), 2.0f);

    /* initialise FB coefficients */
    ir = calloc1d( (ir_len+ir_pad) * nCH, sizeof(float));
    irFB = malloc1d(nBands*nCH*nTimeSlots*sizeof(float_complex));
    for(nd=0; nd<N_dirs; nd++){
        for(j=0; j<ir_len; j++)
            for(i=0; i<nCH; i++)
                ir[j*nCH+i] = hIR[nd*nCH*ir_len + i*ir_len + j];
        qmfAnalyse(ir, ir_len+ir_pad, nCH, hopSize, hybridmode, irFB);
        for(nm=0; nm<nCH; nm++){
            for(i=0; i<nBands; i++){
                irFB_energy = 0;
                for(t=0; t<nTimeSlots; t++)
                    irFB_energy += powf(cabsf(irFB[i*nTimeSlots*nCH + t*nCH + nm]), 2.0f); /* out_nBands x nTimeslots x nCH */
                irFB_gain = sqrtf(irFB_energy/MAX(centerImpulseFB_energy[i], 2.23e-8f));
                cross = cmplxf(0.0f,0.0f);
                for(t=0; t<nTimeSlots; t++)
                    cross = ccaddf(cross, ccmulf(irFB[i*nTimeSlots*nCH + t*nCH + nm], conjf(centerImpulseFB[i*nTimeSlots + t])));
                phase = atan2f(cimagf(cross), crealf(cross));
                hFB[i*nCH*N_dirs + nm*N_dirs + nd] = crmulf( cexpf(cmplxf(0.0f, phase)), irFB_gain);
            }
        }
    }

    /* clean-up */
    free(maxIdx);
    free(centerImpulse);
    free(centerImpulseFB_energy);
    free(centerImpulseFB);
    free(ir);
    free(irFB);
}
