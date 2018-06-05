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
 *     saf_hrir_internal.c
 * Description:
 *     A collection of head-related impulse-response (HRIR)- related functions.
 * Dependencies:
 *     saf_utilities, afSTFTlib
 * Author, date created:
 *     Leo McCormack, 12.12.2016
 */

#include "saf_hrir.h"
#include "saf_hrir_internal.h"

void cxcorr
(
    float* a,
    float* b,
    float* x_ab,
    size_t la,
    size_t lb
)
{
    int m, n, negFLAG, arg;
    size_t len, lim;
    
    len = la + lb - 1;
    memset(x_ab, 0, len*sizeof(float));
    for(m=1; m<=len; m++){
        arg = m-(int)la;
        if(arg<0){
            negFLAG = 1;
            lim = la + arg;
        }
        else{
            negFLAG = 0;
            lim = la - arg;
        }
        for(n=1; n<=lim; n++){
            if(negFLAG == 0)
                x_ab[m-1] += (a[arg+n-1] * b[n-1]);
            else
                x_ab[m-1] += (a[n-1] * b[n-arg-1]);
        }
    }
}

/* currently hard coded for a 128 hop size with hybrid mode enabled */
/* Copyright (c) 2015 Juha Vilkamo, MIT license */
static void afAnalyse
(
    float* inTD/* nSamplesTD x nCH */,
    int nSamplesTD,
    int nCH,
    float_complex* outTF /* out_nBands x nTimeslots x nCH */
)
{
    int t, ch, sample, band;
    void* hSTFT;
    complexVector** FrameTF;
    float** tempHopFrameTD;
    int nTimeSlots, hyrbidBands, hopSize;
    
    hopSize = 128;
    hyrbidBands =hopSize+5; 
    nTimeSlots = nSamplesTD/hopSize;
    
    /* allocate memory */
    afSTFTinit(&(hSTFT), hopSize, nCH, 1, 0, 1);
    FrameTF = (complexVector**)malloc2d(nTimeSlots, nCH, sizeof(complexVector));
    for(t=0; t<nTimeSlots; t++) {
        for(ch=0; ch< nCH; ch++) {
            FrameTF[t][ch].re = (float*)calloc(hyrbidBands, sizeof(float));
            FrameTF[t][ch].im = (float*)calloc(hyrbidBands, sizeof(float));
        }
    }
    tempHopFrameTD = (float**)malloc2d(nCH, hopSize, sizeof(float));
    
    /* perform TF transform */
    for ( t=0; t< nTimeSlots; t++) {
        for( ch=0; ch < nCH; ch++)
            for ( sample=0; sample < hopSize; sample++)
                tempHopFrameTD[ch][sample] = inTD[(sample + t*hopSize)*nCH + ch];
        afSTFTforward(hSTFT, (float**)tempHopFrameTD, (complexVector*)FrameTF[t]);
    }
    
    /* save result to output */
    for(band=0; band<hyrbidBands; band++)
        for ( t=0; t<nTimeSlots; t++)
            for( ch=0; ch < nCH; ch++)
                outTF[band*nTimeSlots*nCH + t*nCH + ch] = cmplxf(FrameTF[t][ch].re[band], FrameTF[t][ch].im[band]);
    
    /* clean-up */
    afSTFTfree(hSTFT);
    for (t = 0; t<nTimeSlots; t++) {
        for(ch=0; ch< nCH; ch++) {
            free(FrameTF[t][ch].re);
            free(FrameTF[t][ch].im);
        }
    }
    free2d((void**)FrameTF, nTimeSlots);
    free2d((void**)tempHopFrameTD, nCH);
    
}

void FIRtoFilterbankCoeffs
(
    float* hIR /*N_dirs x nCH x ir_len*/,
    int N_dirs,
    int nCH,
    int ir_len,
    int nBands,
    float_complex** hFB /* nBands x nCH x N_dirs */
)
{
    int i, j, t, nd, nm, nTimeSlots, ir_pad, hopSize;
    int* maxIdx;
    float maxVal, idxDel, irFB_energy, irFB_gain, phase;
    float* centerImpulse, *centerImpulseFB_energy, *ir;
    float_complex cross;
    float_complex* centerImpulseFB, *irFB;
    
    ir_pad = 1024+512;
    hopSize = 128;
    nTimeSlots = (ir_len+ir_pad)/hopSize;
    maxIdx = calloc(nCH,sizeof(int));
    centerImpulse = calloc(ir_len+ir_pad, sizeof(float));
    
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
    centerImpulseFB = malloc(nBands*nTimeSlots*nCH*sizeof(float_complex));
    afAnalyse(centerImpulse, ir_len+ir_pad, 1, centerImpulseFB);
    centerImpulseFB_energy = calloc(nBands, sizeof(float));
    for(i=0; i<nBands; i++)
        for(t=0; t<nTimeSlots; t++)
            centerImpulseFB_energy[i] += powf(cabsf(centerImpulseFB[i*nTimeSlots + t]), 2.0f);
    
    /* initialise FB coefficients */
    (*hFB) = malloc(nBands*nCH*N_dirs*sizeof(float_complex));
    ir = calloc( (ir_len+ir_pad) * nCH, sizeof(float));
    irFB = malloc(nBands*nCH*nTimeSlots*sizeof(float_complex));
    for(nd=0; nd<N_dirs; nd++){
        for(j=0; j<ir_len; j++)
            for(i=0; i<nCH; i++)
                ir[j*nCH+i] = hIR[nd*nCH*ir_len + i*ir_len + j];
        afAnalyse(ir, ir_len+ir_pad, nCH, irFB);
        for(nm=0; nm<nCH; nm++){
            for(i=0; i<nBands; i++){
                irFB_energy = 0;
                for(t=0; t<nTimeSlots; t++)
                    irFB_energy += powf(cabsf(irFB[i*nTimeSlots*nCH + t*nCH + nm]), 2.0f); /* out_nBands x nTimeslots x nCH */
                irFB_gain = sqrtf(irFB_energy/centerImpulseFB_energy[i]);
                cross = cmplxf(0.0f,0.0f);
                for(t=0; t<nTimeSlots; t++)
                    cross = ccaddf(cross, ccmulf(irFB[i*nTimeSlots*nCH + t*nCH + nm], conjf(centerImpulseFB[i*nTimeSlots + t])));
                phase = atan2f(cimagf(cross), crealf(cross));
                (*hFB)[i*nCH*N_dirs + nm*N_dirs  +nd] = crmulf( cexpf(cmplxf(0.0f, phase)), irFB_gain);
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
