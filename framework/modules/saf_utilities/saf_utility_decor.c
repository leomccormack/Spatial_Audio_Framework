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
 * @author Leo McCormack
 * @date 30.07.2018
 */

#include "saf_utility_decor.h"
#include "saf_utilities.h"

/**
 * Random permutation of a vector of integers
 */
static void randperm(int len, int* randperm)
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
    maxMilliseconds = MIN(80.0f, (maxTFdelay-1.0f)*(float)hopSize/fs*1000.0f);
    for(band=0; band<nFreqs; band++){
        delayRangeMax[band] = MAX(7.0f, MIN(maxMilliseconds, 50.0f*1000.0f/(freqs[band]+2.23e-9f)));
        delayRangeMin[band] = MAX(3.0f, MIN(20.0f, 10.0f*1000.0f/(freqs[band]+2.23e-9f)));
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
            delayTF[band*nChannels+ch]  = MAX((int)(delays[band*nChannels+ch]/1000.0f*fs/(float)hopSize + 0.5f)-1,0); /* remove -1 for matlab version */
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

typedef struct _latticeAPF{
    int order;
    int filterLength;
    float** coeffs;
    float_complex* buffer;

}latticeAPF;

typedef struct _latticeDecor_data{
    int nCH;
    int nCutoffs;
    int nBands;
    int nXtimeslots;
    int maxBufferLen;
    int* orders;
    int* TF_delays;
    latticeAPF** lttc_apf;

    /* run-time */
    float_complex*** delayBuffers;
    int* wIdx;
    int* rIdx;

}latticeDecor_data;

void latticeDecorrelator_create
(
    void** phDecor,
    int nCH,
    int* orders,
    float* freqCutoffs,
    int* fixedDelays,
    int nCutoffs,
    float* freqVector,
    int lookupOffset,
    int nBands
)
{
    *phDecor = malloc1d(sizeof(latticeDecor_data));
    latticeDecor_data *h = (latticeDecor_data*)(*phDecor);
    int i, band, ch, o, filterIdx, maxDelay;

    h->nCH = nCH;
    h->nCutoffs = nCutoffs;
    h->nBands = nBands;

    /* alloc */
    h->orders = malloc1d(nCutoffs*sizeof(int));
    memcpy(h->orders, orders, nCutoffs*sizeof(int));
    h->TF_delays = malloc1d(nBands * sizeof(int));
    h->lttc_apf = (latticeAPF**)malloc2d(nBands, nCH, sizeof(latticeAPF));

    /* set-up delays and allpass filters per band and channel */
    maxDelay = 0;
    for(band=0; band<nBands; band++){
        filterIdx = -1;

        /* Find filter index for this band: */
        for(o=0; o<nCutoffs; o++){
            if(freqVector[band]<freqCutoffs[o]){
                filterIdx = o;
                break;
            }
        }

        /* Fixed frequency dependent delays */
        if(filterIdx==-1)
            h->TF_delays[band] = fixedDelays[nCutoffs]+1;
        else
            h->TF_delays[band] = fixedDelays[filterIdx]+1;

        /* keep track of maximum delay */
        if(h->TF_delays[band]>maxDelay)
            maxDelay = h->TF_delays[band];

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
                h->lttc_apf[band][ch].buffer = calloc1d(h->orders[filterIdx], sizeof(float_complex));

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
                    default: assert(0); /* Unsupported filter order specified */
                }

                /* denominator coefficients */
                for(i=0; i<orders[filterIdx]; i++)
                    h->lttc_apf[band][ch].coeffs[1][i] = h->lttc_apf[band][ch].coeffs[0][orders[filterIdx]-i-1];
            }
        }
    }

    /* Run-time */
    h->maxBufferLen = maxDelay;
    h->delayBuffers = (float_complex***)calloc3d(nBands, nCH, h->maxBufferLen, sizeof(float_complex));
    h->rIdx = malloc1d(nBands*sizeof(int));
    for(band=0; band<nBands; band++)
        h->rIdx[band] = h->TF_delays[band]-1;
    h->wIdx = calloc1d(nBands,sizeof(int));
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

    /* Apply fixed delay */
    for(t=0; t<nTimeSlots; t++){
        for(band=0; band <h->nBands; band++){
            for(ch=0; ch<h->nCH; ch++){
                h->delayBuffers[band][ch][h->wIdx[band]] = inFrame[band][ch][t];
                decorFrame[band][ch][t] = h->delayBuffers[band][ch][h->rIdx[band]];
            }

            /* increment and wrap-around as needed */
            h->rIdx[band]++;
            h->wIdx[band]++;
            if( h->rIdx[band] >= h->TF_delays[band] )
                h->rIdx[band] = 0;
            if( h->wIdx[band] >= h->TF_delays[band] )
                h->wIdx[band] = 0;
        }
    }

    /* Apply lattice allpass filters */
#if _MSC_VER >= 1900
    for(band=0; band <h->nBands; band++){
        for(ch=0; ch<h->nCH; ch++){
            if(h->lttc_apf[band][ch].buffer!=NULL){ /* only if filter is defined... */
                for(t=0; t<nTimeSlots; t++){
                    xtmp = decorFrame[band][ch][t];
                    ytmp = ccaddf(h->lttc_apf[band][ch].buffer[0], crmulf(xtmp, (h->lttc_apf[band][ch].coeffs[0][0])));
                    decorFrame[band][ch][t] = ytmp;

                    /* propagate through the rest of the lattice filter structure */
                    for(i=0; i<h->lttc_apf[band][ch].order-1; i++){
                        h->lttc_apf[band][ch].buffer[i] = ccaddf(h->lttc_apf[band][ch].buffer[i+1],
                                                                 ccsubf(crmulf(xtmp, h->lttc_apf[band][ch].coeffs[0][i+1]),   /* numerator */
                                                                        crmulf(ytmp, h->lttc_apf[band][ch].coeffs[1][i+1])));  /* denominator */
                    }
                }
            }
        }
    }
#else
    for(band=0; band <h->nBands; band++){
        for(ch=0; ch<h->nCH; ch++){
            if(h->lttc_apf[band][ch].buffer!=NULL){ /* only if filter is defined... */
                for(t=0; t<nTimeSlots; t++){
                    xtmp = decorFrame[band][ch][t];
                    ytmp = h->lttc_apf[band][ch].buffer[0] + xtmp * (h->lttc_apf[band][ch].coeffs[0][0]);
                    decorFrame[band][ch][t] = ytmp;
                    
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
