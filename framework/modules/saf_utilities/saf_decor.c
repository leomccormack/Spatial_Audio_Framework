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

/*
 * Filename: saf_decor.c
 * ---------------------
 * A collection of signal decorrelators.
 *
 * Dependencies:
 *     none
 * Author, date created:
 *     Leo McCormack, 30.07.2018
 */

#include "saf_decor.h"
#include "saf_utilities.h"

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

#if 0 // unfinished, untested
void synthesiseNoiseReverb
(
    int nCH,
    float fs,
    float* t60,
    float* fcen_oct,
    int nBands,
    int flattenFLAG,
    float* rir_filt,
    int* rir_len
)
{
    int i, j, k, rir_filt_len, rir_filt_lpad, rir_filt_lout, filterOrder, nBins;
    float alpha, max_t60, t;
    float *rir, *fcut, *h_filt, *rir_filt_tmp;
    
    filterOrder = 800;
    
    /* find rir length */
    max_t60 = 0.0f;
    for(i=0; i<nBands; i++)
        max_t60 = max_t60 < t60[i] ? t60[i] : max_t60;
    rir_filt_len = (int)(max_t60*fs+0.5f); /* length of RIRs */
    rir_filt_lpad = rir_filt_len + filterOrder; /* length after convolution with filterbank */
    rir_filt_lout = rir_filt_len + filterOrder/2; /* truncated output length */
    
    /* Generate Noise RIRs and shape with exponential envelopes */
    rir = calloc1d(nCH*nBands*rir_filt_lpad, sizeof(float));
    for(i=0; i<nCH; i++){
        for(j=0; j<nBands; j++){
            /* decay constants for t60 */
            alpha = 3.0f*logf(10.0f)/t60[j];
            for(k=0, t=0.0f; k<rir_filt_len; k++, t+=1.0f/fs)
                rir[i*nBands*rir_filt_lpad + j*rir_filt_lpad + k] = expf(-t*alpha) *     /* envelope */
                                                                  2.0f * ((float)rand()/(float)RAND_MAX-0.5f); /* whitenoise */
        }
    }
    
    /* get filterbank FIRs - octave bands */
    fcut = malloc1d((nBands-1)*sizeof(float));
    h_filt = malloc1d(nBands*(filterOrder+1)*sizeof(float));
    getOctaveBandCutoffFreqs(fcen_oct, nBands, fcut);
    FIRFilterbank(filterOrder, fcut, (nBands-1), fs, WINDOWING_FUNCTION_HAMMING, 1, h_filt);
    
    /* filter RIRs with filterbank */
    memset(rir_filt, 0, nCH*rir_filt_lout*sizeof(float));
    rir_filt_tmp = malloc1d(nBands*rir_filt_lpad*sizeof(float));
    for(i=0; i<nCH; i++){
        fftconv(&rir[i*nBands*rir_filt_lpad], h_filt, rir_filt_lpad, filterOrder+1, nBands, rir_filt_tmp);
        for(j=0; j<nBands; j++){
            for(k=0; k<rir_filt_lout; k++)
                rir_filt[i*rir_filt_lout] += rir_filt_tmp[j*rir_filt_lpad+k];
        }
    }
    (*rir_len) = rir_filt_lout;
    
    
    /*  */
    if(flattenFLAG){
        float_complex* tdi_f;
        void* hFFT;
        
        saf_rfft_create(&hFFT, (*rir_len));
        nBins = (*rir_len)/2+1;
        tdi_f = malloc1d(nBins*sizeof(float_complex));
        for(i=0; i<nCH; i++){
            saf_rfft_forward(hFFT, &rir_filt[i*rir_filt_lout], tdi_f);
            
            
            
            
        }
        
        
        
    }
    
    
    
}


#endif
