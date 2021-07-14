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
 * @file saf_utility_decor.c
 * @ingroup Utilities
 * @brief A collection of signal decorrelators
 *
 * @author Leo McCormack
 * @date 30.07.2018
 * @license ISC
 */
 
#include "saf_utilities.h"
#include "saf_externals.h"

/**
 * Internal Lattice all-pass filter structure */
typedef struct _latticeAPF{
    int order;
    int filterLength;
    float** coeffs;
    float_complex* buffer;

}latticeAPF;

/**
 * Internal Lattice all-pass filter based decorrelator structure */
typedef struct _latticeDecor_data{
    int nCH;
    int nCutoffs;
    int nBands;
    int maxBufferLen;
    int* orders;
    int* TF_delays;
    latticeAPF** lttc_apf;
    float enComp_coeff;

    /* run-time */
    float** in_energy;
    float** decor_energy;
    float_complex*** delayBuffers;
    int* wIdx;
    int* rIdx;

}latticeDecor_data;

/**
 * Internal structure used by the transient Ducker */
typedef struct _transientDucker_data{
    int nCH;
    int nBands;
    float** transientDetector1;
    float** transientDetector2;

}transientDucker_data;

void getDecorrelationDelays
(
    int nChannels,  /* number of channels */
    float* freqs,   /* centre frequencies */
    int nFreqs,     /* number of elements in frequency vector */
    float fs,       /* host fs */
    int maxTFdelay, /* max number of time-slots to delay */
    int hopSize,    /* STFT hop size */
    int* delayTF    /* nFreq x nChannels */
)
{
    int band, ch;
    int* randperm_nCH;
    float maxMilliseconds, nChannelsf;
    float* delayRangeMax, *delayRangeMin, *tmp_delays, *delays;
    
    nChannelsf = (float)nChannels;
    randperm_nCH = malloc1d(nChannels*sizeof(int));
    delayRangeMax = malloc1d(nFreqs*sizeof(float));
    delayRangeMin = malloc1d(nFreqs*sizeof(float));
    tmp_delays = malloc1d(nChannels*sizeof(float));
    delays = malloc1d(nFreqs*nChannels*sizeof(float));
    maxMilliseconds = SAF_MIN(80.0f, (maxTFdelay-1.0f)*(float)hopSize/fs*1000.0f);
    for(band=0; band<nFreqs; band++){
        delayRangeMax[band] = SAF_MAX(7.0f, SAF_MIN(maxMilliseconds, 50.0f*1000.0f/(freqs[band]+2.23e-9f)));
        delayRangeMin[band] = SAF_MAX(3.0f, SAF_MIN(20.0f, 10.0f*1000.0f/(freqs[band]+2.23e-9f)));
    }
    for(band=0; band<nFreqs; band++)
        for(ch=0; ch<nChannels; ch++)
            delays[band*nChannels+ch] = (float)ch/nChannelsf + ((float)rand()/(float)RAND_MAX)/nChannelsf;
    for(band=0; band<nFreqs; band++){
        randperm(nChannels, randperm_nCH);
        memcpy(tmp_delays, &delays[band*nChannels], nChannels*sizeof(float));
        for(ch=0; ch<nChannels; ch++)
            delays[band*nChannels+ch] = tmp_delays[randperm_nCH[ch]];
    }
    for(band=0; band<nFreqs; band++){
        for(ch=0; ch<nChannels; ch++){
            delays[band*nChannels+ch] = delays[band*nChannels+ch]*(delayRangeMax[band]-delayRangeMin[band])+delayRangeMin[band];
            delayTF[band*nChannels+ch]  = SAF_MAX((int)(delays[band*nChannels+ch]/1000.0f*fs/(float)hopSize + 0.5f)-1,0);
        }
    }
    
    free(randperm_nCH);
    free(delayRangeMax);
    free(delayRangeMin);
    free(tmp_delays);
    free(delays);
}

void synthesiseNoiseReverb
(
    int nCH,
    float fs,
    float* t60,
    float* fcen_oct,
    int nBands,
    int flattenFLAG,
    float** rir_filt,
    int* rir_len
)
{
    int i, j, k, rir_filt_len, rir_filt_lout, filterOrder;
    float alpha, max_t60, t;
    float *rir, *fcut, *h_filt, *rir_filt_tmp;
    
    filterOrder = 800;
    
    /* find rir length */
    max_t60 = 0.0f;
    for(i=0; i<nBands; i++)
        max_t60 = max_t60 < t60[i] ? t60[i] : max_t60;
    rir_filt_len = (int)(max_t60*fs+0.5f); /* length of RIRs */
    rir_filt_lout = rir_filt_len + filterOrder/2; /* truncated output length */
    
    /* Generate noise and shape with exponentially decaying envelopes */
    rir = calloc1d(nCH*nBands*rir_filt_lout, sizeof(float));
    for(i=0; i<nCH; i++){
        for(j=0; j<nBands; j++){
            /* decay constants for t60 */
            alpha = 3.0f*logf(10.0f)/t60[j];
            for(k=0, t=0.0f; k<rir_filt_len; k++, t+=1.0f/fs)
                rir[i*nBands*rir_filt_lout + j*rir_filt_lout + k] = expf(-t*alpha) *     /* envelope */
                                                                    2.0f * ((float)rand()/(float)RAND_MAX-0.5f); /* whitenoise */
        }
    }

    /* get bank of FIRs filters - octave bands */
    fcut = malloc1d((nBands-1)*sizeof(float));
    h_filt = malloc1d(nBands*(filterOrder+1)*sizeof(float));
    getOctaveBandCutoffFreqs(fcen_oct, nBands, fcut);
    FIRFilterbank(filterOrder, fcut, (nBands-1), fs, WINDOWING_FUNCTION_HAMMING, 1, h_filt);
    
    /* filter RIRs with filterbank */
    (*rir_filt) = realloc1d((*rir_filt), nCH*rir_filt_lout*sizeof(float));
    memset((*rir_filt), 0, nCH*rir_filt_lout*sizeof(float));
    rir_filt_tmp = malloc1d(nBands*rir_filt_lout*sizeof(float));
    for(i=0; i<nCH; i++){
        fftfilt(&rir[i*nBands*rir_filt_lout], h_filt, rir_filt_lout, filterOrder+1, nBands, rir_filt_tmp);
        /* sum over bands */
        for(j=0; j<nBands; j++){
            for(k=0; k<rir_filt_lout; k++)
                (*rir_filt)[i*rir_filt_lout+k] += rir_filt_tmp[j*rir_filt_lout+k];
        }
    }
    
    /* equalise, to force flat magnitude response */
    if(flattenFLAG)
        for(i=0; i<nCH; i++)
            flattenMinphase(&((*rir_filt)[i*rir_filt_lout]), rir_filt_lout);
    
    /* remove filterbank delay */
    for(i=0; i<nCH; i++)
        memcpy(&((*rir_filt)[i*rir_filt_len]), &((*rir_filt)[i*rir_filt_lout + filterOrder/2]), rir_filt_len*sizeof(float));
    (*rir_len) = rir_filt_len;
    
    /* clean-up */
    free(rir);
    free(fcut);
    free(h_filt);
    free(rir_filt_tmp);
}

void latticeDecorrelator_create
(
    void** phDecor,
    float fs,
    int hopsize,
    float* freqVector,
    int nBands,
    int nCH,
    int* orders,
    float* freqCutoffs,
    int nCutoffs,
    int maxDelay,
    int lookupOffset,
    float enComp_coeff
)
{
    *phDecor = malloc1d(sizeof(latticeDecor_data));
    latticeDecor_data *h = (latticeDecor_data*)(*phDecor);
    int i, band, ch, o, filterIdx;

    h->nCH = nCH;
    h->nCutoffs = nCutoffs;
    h->nBands = nBands;

    /* alloc */
    h->orders = malloc1d(nCutoffs*sizeof(int));
    memcpy(h->orders, orders, nCutoffs*sizeof(int));
    h->TF_delays = malloc1d(nBands * nCH * sizeof(int));
    h->lttc_apf = (latticeAPF**)malloc2d(nBands, nCH, sizeof(latticeAPF));
    h->enComp_coeff = enComp_coeff;
    h->in_energy = (float**)calloc2d(nBands, nCH, sizeof(float));
    h->decor_energy = (float**)calloc2d(nBands, nCH, sizeof(float));

    /* Static delays */
    getDecorrelationDelays(h->nCH, freqVector, h->nBands, fs, maxDelay, hopsize, h->TF_delays);

    /* Find true max delay */
    maxDelay = 0;
    for(i=0; i<nBands*nCH; i++)
        maxDelay = h->TF_delays[i] > maxDelay ? h->TF_delays[i] : maxDelay;

    /* set-up allpass filters per band and channel */
    for(band=0; band<nBands; band++){
        filterIdx = -1;

        /* Find filter index for this band: */
        for(o=0; o<nCutoffs; o++){
            if(freqVector[band]<freqCutoffs[o]){
                filterIdx = o;
                break;
            }
        }

        /* Pull lattice allpass filter coefficients from database */
        for(ch=0; ch<nCH; ch++){
            if(filterIdx == -1){
                /* Not needed... */
                h->lttc_apf[band][ch].order = -1;
                h->lttc_apf[band][ch].filterLength = 0;
                h->lttc_apf[band][ch].coeffs = NULL;
                h->lttc_apf[band][ch].buffer = NULL;
            }
            else{
                h->lttc_apf[band][ch].order = h->orders[filterIdx];
                h->lttc_apf[band][ch].coeffs = (float**)malloc2d(2, h->orders[filterIdx], sizeof(float));
                h->lttc_apf[band][ch].buffer = calloc1d(h->lttc_apf[band][ch].order, sizeof(float_complex));

                /* numerator coefficients */
                switch(orders[filterIdx]){
                    case 20: memcpy(h->lttc_apf[band][ch].coeffs[0], __lattice_coeffs_o20[ch+lookupOffset], (h->orders[filterIdx])*sizeof(float)); break;
                    case 18: memcpy(h->lttc_apf[band][ch].coeffs[0], __lattice_coeffs_o18[ch+lookupOffset], (h->orders[filterIdx])*sizeof(float)); break;
                    case 16: memcpy(h->lttc_apf[band][ch].coeffs[0], __lattice_coeffs_o16[ch+lookupOffset], (h->orders[filterIdx])*sizeof(float)); break;
                    case 15: memcpy(h->lttc_apf[band][ch].coeffs[0], __lattice_coeffs_o15[ch+lookupOffset], (h->orders[filterIdx])*sizeof(float)); break;
                    case 14: memcpy(h->lttc_apf[band][ch].coeffs[0], __lattice_coeffs_o14[ch+lookupOffset], (h->orders[filterIdx])*sizeof(float)); break;
                    case 12: memcpy(h->lttc_apf[band][ch].coeffs[0], __lattice_coeffs_o12[ch+lookupOffset], (h->orders[filterIdx])*sizeof(float)); break;
                    case 10: memcpy(h->lttc_apf[band][ch].coeffs[0], __lattice_coeffs_o10[ch+lookupOffset], (h->orders[filterIdx])*sizeof(float)); break;
                    case 8:  memcpy(h->lttc_apf[band][ch].coeffs[0], __lattice_coeffs_o8[ch+lookupOffset],  (h->orders[filterIdx])*sizeof(float)); break;
                    case 6:  memcpy(h->lttc_apf[band][ch].coeffs[0], __lattice_coeffs_o6[ch+lookupOffset],  (h->orders[filterIdx])*sizeof(float)); break;
                    case 4:  memcpy(h->lttc_apf[band][ch].coeffs[0], __lattice_coeffs_o4[ch+lookupOffset],  (h->orders[filterIdx])*sizeof(float)); break;
                    case 3:  memcpy(h->lttc_apf[band][ch].coeffs[0], __lattice_coeffs_o3[ch+lookupOffset],  (h->orders[filterIdx])*sizeof(float)); break;
                    case 2:  memcpy(h->lttc_apf[band][ch].coeffs[0], __lattice_coeffs_o2[ch+lookupOffset],  (h->orders[filterIdx])*sizeof(float)); break;
                    default: saf_print_error("Unsupported filter order was specified"); break;
                }

                /* denominator coefficients */
                for(i=0; i<orders[filterIdx]; i++)
                    h->lttc_apf[band][ch].coeffs[1][i] = h->lttc_apf[band][ch].coeffs[0][orders[filterIdx]-i-1];
            }
        }
    }

    /* Run-time */
    h->maxBufferLen = maxDelay+1;
    h->delayBuffers = (float_complex***)calloc3d(nBands, nCH, h->maxBufferLen, sizeof(float_complex));
    h->wIdx = malloc1d(nBands*nCH*sizeof(int));
    for(band=0; band<nBands; band++)
        for(ch=0; ch<nCH; ch++)
            h->wIdx[band*nCH+ch] = h->TF_delays[band*nCH+ch];
    h->rIdx = calloc1d(nBands*nCH,sizeof(int));
}

void latticeDecorrelator_destroy
(
    void** phDecor
)
{
    latticeDecor_data *h = (latticeDecor_data*)(*phDecor);
    int band, ch;

    if(h!=NULL){
        free(h->orders);
        free(h->TF_delays);
        for(band=0; band <h->nBands; band++){
            for(ch=0; ch<h->nCH; ch++){
                free(h->lttc_apf[band][ch].buffer);
                free(h->lttc_apf[band][ch].coeffs);
            }
        }
        free(h->lttc_apf);

        free(h->delayBuffers);
        free(h->wIdx);
        free(h->rIdx);
        free(h);
        h=NULL;
        *phDecor = NULL;
    }
}

void latticeDecorrelator_reset
(
    void* hDecor
)
{
    latticeDecor_data *h = (latticeDecor_data*)(hDecor);
    int band, ch;

    memset(FLATTEN3D(h->delayBuffers), 0, h->nBands * h->nCH * h->maxBufferLen * sizeof(float_complex));
    for(band=0; band<h->nBands; band++)
        for(ch=0; ch<h->nCH; ch++)
            if(h->lttc_apf[band][ch].buffer!=NULL)
                memset(h->lttc_apf[band][ch].buffer, 0, h->lttc_apf[band][ch].order * sizeof(float_complex));
    memset(FLATTEN2D(h->in_energy), 0, h->nBands*h->nCH*sizeof(float));
    memset(FLATTEN2D(h->decor_energy), 0, h->nBands*h->nCH*sizeof(float));
}

void latticeDecorrelator_apply
(
    void* hDecor,
    float_complex*** inFrame,
    int nTimeSlots,
    float_complex*** decorFrame
)
{
    latticeDecor_data *h = (latticeDecor_data*)(hDecor);
    int band, ch, t, i;
    float_complex ytmp, xtmp;

    for(t=0; t<nTimeSlots; t++){
        for(band=0; band <h->nBands; band++){
            for(ch=0; ch<h->nCH; ch++){
                /* Apply fixed delay */
                h->delayBuffers[band][ch][h->wIdx[band*(h->nCH) + ch]] = inFrame[band][ch][t];
                decorFrame[band][ch][t] = h->delayBuffers[band][ch][h->rIdx[band*(h->nCH) + ch]];

                /* increment and wrap-around as needed */
                h->rIdx[band*(h->nCH) + ch]++;
                h->wIdx[band*(h->nCH) + ch]++;
                if( h->rIdx[band*(h->nCH) + ch] > h->TF_delays[band*(h->nCH) + ch] )
                    h->rIdx[band*(h->nCH) + ch] = 0;
                if( h->wIdx[band*(h->nCH) + ch] > h->TF_delays[band*(h->nCH) + ch] )
                    h->wIdx[band*(h->nCH) + ch] = 0;
            }
        }
    }

    /* Apply lattice allpass filters */
#if _MSC_VER >= 1900
    for(band=0; band <h->nBands; band++){
        for(ch=0; ch<h->nCH; ch++){
            if(h->lttc_apf[band][ch].buffer!=NULL){ /* only if filter is defined... */
                for(t=0; t<nTimeSlots; t++){
                    /* Compute energy of the input */
                    h->in_energy[band][ch] = (1.0f-h->enComp_coeff)*powf(cabsf(inFrame[band][ch][t]),2.0f) + (h->enComp_coeff)*h->in_energy[band][ch];

                    /* First tap in filter */
                    xtmp = decorFrame[band][ch][t];
                    ytmp = ccaddf(h->lttc_apf[band][ch].buffer[0], crmulf(xtmp, (h->lttc_apf[band][ch].coeffs[0][0])));
                    decorFrame[band][ch][t] = ytmp;

                    /* Energy compensation */
                    h->decor_energy[band][ch] = (1.0f-h->enComp_coeff)*powf(cabsf(decorFrame[band][ch][t]),2.0f) + (h->enComp_coeff)*h->decor_energy[band][ch];
                    decorFrame[band][ch][t] = crmulf(decorFrame[band][ch][t], SAF_MIN(sqrtf(h->in_energy[band][ch]/(h->decor_energy[band][ch]+2.23e-9f)), 1.0f));

                    /* propagate through the rest of the lattice filter structure */
                    for(i=0; i<h->lttc_apf[band][ch].order-1; i++){
                        h->lttc_apf[band][ch].buffer[i] = ccaddf(h->lttc_apf[band][ch].buffer[i+1],
                                                                 ccsubf(crmulf(xtmp, h->lttc_apf[band][ch].coeffs[0][i+1]),    /* numerator */
                                                                        crmulf(ytmp, h->lttc_apf[band][ch].coeffs[1][i+1])));  /* denominator */
                    }
                }
            }
        }
    }
#else
    for(band=0; band <h->nBands; band++){
        for(ch=0; ch<h->nCH; ch++){
            if(h->lttc_apf[band][ch].buffer!=NULL){ /* only if filter is defined */
                for(t=0; t<nTimeSlots; t++){
                    /* Compute energy of the input */
                    h->in_energy[band][ch] = (1.0f-h->enComp_coeff)*powf(cabsf(inFrame[band][ch][t]),2.0f) + (h->enComp_coeff)*h->in_energy[band][ch];

                    /* First tap in filter */
                    xtmp = decorFrame[band][ch][t];
                    ytmp = xtmp * (h->lttc_apf[band][ch].coeffs[0][0]) + h->lttc_apf[band][ch].buffer[0];
                    decorFrame[band][ch][t] = ytmp;

                    /* Energy compensation */
                    h->decor_energy[band][ch] = (1.0f-h->enComp_coeff)*powf(cabsf(decorFrame[band][ch][t]), 2.0f) + (h->enComp_coeff)*(h->decor_energy[band][ch]);
                    decorFrame[band][ch][t] *= SAF_MIN(sqrtf(h->in_energy[band][ch]/(h->decor_energy[band][ch]+2.23e-9f)), 1.0f);

                    /* propagate through the rest of the lattice filter structure */
                    for(i=0; i<h->lttc_apf[band][ch].order-1; i++){
                        h->lttc_apf[band][ch].buffer[i] = h->lttc_apf[band][ch].buffer[i+1] +
                                                          h->lttc_apf[band][ch].coeffs[0][i+1] * xtmp - /* numerator */
                                                          h->lttc_apf[band][ch].coeffs[1][i+1] * ytmp;  /* denominator */
                    }
                }
            }
        }
    }
#endif
}

void transientDucker_create
(
    void** phDucker,
    int nCH,
    int nBands
)
{
    *phDucker = malloc1d(sizeof(transientDucker_data));
    transientDucker_data *h = (transientDucker_data*)(*phDucker);

    h->nCH = nCH;
    h->nBands = nBands;
    h->transientDetector1 = (float**)calloc2d(nBands, nCH, sizeof(float));
    h->transientDetector2 = (float**)calloc2d(nBands, nCH, sizeof(float));
}

void transientDucker_destroy
(
    void** phDucker
)
{
    transientDucker_data *h = (transientDucker_data*)(*phDucker);

    if(h!=NULL){
        free(h->transientDetector1);
        free(h->transientDetector2);
        free(h);
        h=NULL;
        *phDucker = NULL;
    }
}

void transientDucker_apply
(
    void* hDucker,
    float_complex*** inFrame,
    int nTimeSlots,
    float alpha,
    float beta,
    float_complex*** residualFrame,
    float_complex*** transientFrame
)
{
    transientDucker_data *h = (transientDucker_data*)(hDucker);
    int band, i, t;
    float detectorEne, transientEQ;
    //const float alpha = 0.95f;
    //const float beta = 0.995f;

    for(band=0; band<h->nBands; band++){
        for(i=0; i<h->nCH; i++){
            for(t=0; t<nTimeSlots; t++){
                detectorEne = cabsf(inFrame[band][i][t]);
                detectorEne = detectorEne*detectorEne;
                h->transientDetector1[band][i] *= alpha;
                if(h->transientDetector1[band][i]<detectorEne)
                    h->transientDetector1[band][i] = detectorEne;
                h->transientDetector2[band][i] = h->transientDetector2[band][i]*beta + (1.0f-beta)*(h->transientDetector1[band][i]);
                if(h->transientDetector2[band][i] > h->transientDetector1[band][i])
                    h->transientDetector2[band][i] = h->transientDetector1[band][i];
                transientEQ = SAF_MIN(1.0f, 4.0f * (h->transientDetector2[band][i])/(h->transientDetector1[band][i]+2.23e-9f));
#ifdef _MSC_VER
                if(residualFrame!=NULL)
                    residualFrame[band][i][t] = crmulf(inFrame[band][i][t], transientEQ);
                if(transientFrame!=NULL)
                    transientFrame[band][i][t] = crmulf(inFrame[band][i][t], 1.0f-transientEQ);
#else
                if(residualFrame!=NULL)
                    residualFrame[band][i][t] = inFrame[band][i][t] * transientEQ;
                if(transientFrame!=NULL)
                    transientFrame[band][i][t] = inFrame[band][i][t] * (1.0f-transientEQ);
#endif
            }
        }
    }
}
