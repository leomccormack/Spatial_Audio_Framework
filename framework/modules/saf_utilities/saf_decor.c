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



