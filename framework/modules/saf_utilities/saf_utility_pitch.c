/*
 * Copyright 2019 Leo McCormack
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
 * @file saf_utility_pitch.c
 * @ingroup Utilities
 * @brief A collection of pitch shifting algorithms
 *
 * @author Leo McCormack
 * @date 04.05.2020
 * @license ISC
 */

#include "saf_utilities.h"

/**
 * Enables the use of the saf_utility_fft wrapper, instead of "smbFft" used in
 * the original implementation */
#define SMB_ENABLE_SAF_FFT

/* ========================================================================== */
/*                              SMB PitchShifter                              */
/* ========================================================================== */

#ifndef SMB_ENABLE_SAF_FFT
static void smbFft(float *fftBuffer, int fftFrameSize, long sign)
/*
    FFT routine, (C)1996 S.M.Bernsee. Sign = -1 is FFT, 1 is iFFT (inverse)
    Fills fftBuffer[0...2*fftFrameSize-1] with the Fourier transform of the
    time domain data in fftBuffer[0...2*fftFrameSize-1]. The FFT array takes
    and returns the cosine and sine parts in an interleaved manner, ie.
    fftBuffer[0] = cosPart[0], fftBuffer[1] = sinPart[0], asf. fftFrameSize
    must be a power of 2. It expects a complex input signal (see footnote 2),
    ie. when working with 'common' audio signals our input signal has to be
    passed as {in[0],0.,in[1],0.,in[2],0.,...} asf. In that case, the transform
    of the frequencies of interest is in fftBuffer[0...fftFrameSize].
*/
{
    float wr, wi, arg, *p1, *p2, temp;
    float tr, ti, ur, ui, *p1r, *p1i, *p2r, *p2i;
    long i, bitm, j, le, le2, k;

    for (i = 2; i < 2*fftFrameSize-2; i += 2) {
        for (bitm = 2, j = 0; bitm < 2*fftFrameSize; bitm <<= 1) {
            if (i & bitm) j++;
            j <<= 1;
        }
        if (i < j) {
            p1 = fftBuffer+i; p2 = fftBuffer+j;
            temp = *p1; *(p1++) = *p2;
            *(p2++) = temp; temp = *p1;
            *p1 = *p2; *p2 = temp;
        }
    }
    for (k = 0, le = 2; k < (int)(log(fftFrameSize)/log(2.)+.5); k++) {
        le <<= 1;
        le2 = le>>1;
        ur = 1.0;
        ui = 0.0;
        arg = SAF_PI / (le2>>1);
        wr = cos(arg);
        wi = sign*sin(arg);
        for (j = 0; j < le2; j += 2) {
            p1r = fftBuffer+j; p1i = p1r+1;
            p2r = p1r+le2; p2i = p2r+1;
            for (i = j; i < 2*fftFrameSize; i += le) {
                tr = *p2r * ur - *p2i * ui;
                ti = *p2r * ui + *p2i * ur;
                *p2r = *p1r - tr; *p2i = *p1i - ti;
                *p1r += tr; *p1i += ti;
                p1r += le; p1i += le;
                p2r += le; p2i += le;
            }
            tr = ur*wr - ui*wi;
            ui = ur*wi + ui*wr;
            ur = tr;
        }
    }
}
#endif

/**
 * Main structure for the SMB pitch shifter
 */
typedef struct _smb_pitchShift_data
{
    /* parameters */
    int fftFrameSize, osamp, nCH;
    float sampleRate, pitchShiftFactor;

    /* internals */
    void* hFFT;
    float* window;
    float** gInFIFO, **gOutFIFO;
#ifdef SMB_ENABLE_SAF_FFT
    float_complex** gFFTworksp_td;
    float_complex** gFFTworksp_fd;
#else
    float** gFFTworksp;
#endif
    float** gLastPhase, **gSumPhase;
    float** gOutputAccum;
    float** gAnaFreq, **gAnaMagn;
    float** gSynFreq, **gSynMagn;
    int* gRover;
    int stepSize, inFifoLatency;

}smb_pitchShift_data;

void smb_pitchShift_create
(
    void** hSmb,
    int nCH,
    int fftFrameSize,
    int osamp,
    float sampleRate
)
{
    *hSmb = malloc(sizeof(smb_pitchShift_data));
    smb_pitchShift_data *h = (smb_pitchShift_data*)(*hSmb);
    int i;

    h->fftFrameSize = fftFrameSize;
    h->osamp = osamp;
    h->sampleRate = sampleRate;
    h->nCH = nCH;
    h->pitchShiftFactor = 1.0f;

    /* internals */
    saf_fft_create(&(h->hFFT), fftFrameSize);
    h->stepSize = fftFrameSize/h->osamp;
    h->inFifoLatency = fftFrameSize - (h->stepSize);
    h->gRover = (int*)malloc1d(nCH * sizeof(int));
    for(i=0; i<nCH; i++)
        h->gRover[i] = h->inFifoLatency;
    h->window = (float*)malloc1d(fftFrameSize*sizeof(float));
    for (i = 0; i < fftFrameSize; i++)
        h->window[i] = -0.5f*cosf(2.0f*SAF_PI*(float)i/(float)fftFrameSize)+0.5f;
    h->gInFIFO = (float**)calloc2d(nCH,fftFrameSize,sizeof(float));
    h->gOutFIFO = (float**)calloc2d(nCH,fftFrameSize,sizeof(float));
#ifdef SMB_ENABLE_SAF_FFT
    h->gFFTworksp_td = (float_complex**)calloc2d(nCH,fftFrameSize,sizeof(float_complex));
    h->gFFTworksp_fd = (float_complex**)calloc2d(nCH,fftFrameSize,sizeof(float_complex));
#else
    h->gFFTworksp = (float**)calloc2d(nCH,2*fftFrameSize,sizeof(float));
#endif
    h->gLastPhase = (float**)calloc2d(nCH,fftFrameSize/2+1,sizeof(float));
    h->gSumPhase = (float**)calloc2d(nCH,fftFrameSize/2+1,sizeof(float));
    h->gOutputAccum = (float**)calloc2d(nCH,2*(h->fftFrameSize),sizeof(float));
    h->gAnaFreq = (float**)calloc2d(nCH,fftFrameSize,sizeof(float));
    h->gAnaMagn = (float**)calloc2d(nCH,fftFrameSize,sizeof(float));
    h->gSynFreq = (float**)malloc2d(nCH,fftFrameSize,sizeof(float));
    h->gSynMagn = (float**)malloc2d(nCH,fftFrameSize,sizeof(float));
}

void smb_pitchShift_destroy
(
    void ** const hSmb
)
{
    smb_pitchShift_data *h = (smb_pitchShift_data*)(*hSmb);
    if(h!=NULL){
        saf_fft_destroy(&(h->hFFT));
        free(h->window);
        free(h->gRover);
        free(h->gInFIFO);
        free(h->gOutFIFO);
#ifdef SMB_ENABLE_SAF_FFT
        free(h->gFFTworksp_td);
        free(h->gFFTworksp_fd);
#else
        free(h->gFFTworksp);
#endif
        free(h->gLastPhase);
        free(h->gSumPhase);
        free(h->gOutputAccum);
        free(h->gAnaFreq);
        free(h->gAnaMagn);
        free(h->gSynFreq);
        free(h->gSynMagn);
        free(h);
        h=NULL;
    }
}

/* Adapted from: http://blogs.zynaptiq.com/bernsee/pitch-shifting-using-the-ft/
 * Original Copyright notice:
 * COPYRIGHT 1999-2015 Stephan M. Bernsee <s.bernsee [AT] zynaptiq [DOT] com>
 *
 *                         The Wide Open License (WOL)
 *
 * Permission to use, copy, modify, distribute and sell this software and its
 * documentation for any purpose is hereby granted without fee, provided that
 * the above copyright notice and this license appear in all source copies.
 * THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY OF
 * ANY KIND. See http://www.dspguru.com/wol.htm for more information.
 */
void smb_pitchShift_apply
(
    void* hSmb,
    float pitchShift,
    int frameSize,
    float *indata,
    float *outdata
)
{
    smb_pitchShift_data *h = (smb_pitchShift_data*)(hSmb);
    float magn, phase, tmp, real, imag;
    float freqPerBin, expct;
    int ch, i, k, qpd, index, fftFrameSize2;

    /* set up some handy variables */
    fftFrameSize2 = h->fftFrameSize/2;
    freqPerBin = h->sampleRate/(float)h->fftFrameSize;
    expct = 2.0f*SAF_PI*(float)(h->stepSize)/(float)h->fftFrameSize;

    /* flush buffers */
    if(h->pitchShiftFactor!=pitchShift){
        h->pitchShiftFactor = pitchShift;
        for(ch=0; ch<h->nCH; ch++){
            memset(h->gOutputAccum[ch], 0, h->stepSize*sizeof(float));
            memset(h->gLastPhase[ch], 0, (h->fftFrameSize/2+1)*sizeof(float));
            memset(h->gSumPhase[ch], 0, (h->fftFrameSize/2+1)*sizeof(float));
        }
    }

    /* main processing loop */
    for(ch=0; ch<h->nCH; ch++){
        for (i = 0; i < frameSize; i++){
            /* As long as we have not yet collected enough data just read in */
            h->gInFIFO[ch][h->gRover[ch]] = indata[ch*frameSize+i];
            outdata[ch*frameSize+i] = h->gOutFIFO[ch][h->gRover[ch]-(h->inFifoLatency)];
            h->gRover[ch]++;

            /* now we have enough data for processing */
            if (h->gRover[ch] >= h->fftFrameSize) {
                h->gRover[ch] = h->inFifoLatency;

#ifdef SMB_ENABLE_SAF_FFT
                for (k = 0; k < h->fftFrameSize;k++)
                    h->gFFTworksp_td[ch][k] = cmplxf(h->gInFIFO[ch][k] * (float)h->window[k], 0.0f);
                saf_fft_forward(h->hFFT, h->gFFTworksp_td[ch], h->gFFTworksp_fd[ch]);
#else
                /* do windowing and re,im interleave */
                for (k = 0; k < h->fftFrameSize;k++) {
                    h->gFFTworksp[ch][2*k] = h->gInFIFO[ch][k] * h->window[k];
                    h->gFFTworksp[ch][2*k+1] = 0.;
                }

                /* do transform */
                smbFft(h->gFFTworksp[ch], h->fftFrameSize, -1);
#endif

                /* ***************** ANALYSIS ******************* */
                for (k = 0; k <= fftFrameSize2; k++) {
#ifdef SMB_ENABLE_SAF_FFT
                    real = crealf(h->gFFTworksp_fd[ch][k]);
                    imag = cimagf(h->gFFTworksp_fd[ch][k]);
#else
                    /* de-interlace FFT buffer */
                    real = h->gFFTworksp[ch][2*k];
                    imag = h->gFFTworksp[ch][2*k+1];
#endif
                    /* compute magnitude and phase */
                    magn = 2.0f*sqrtf(real*real + imag*imag);
                    phase = atan2f(imag,real);

                    /* compute phase difference */
                    tmp = phase - h->gLastPhase[ch][k];
                    h->gLastPhase[ch][k] = phase;

                    /* subtract expected phase difference */
                    tmp -= (float)k*expct;

                    /* map delta phase into +/- Pi interval */
                    qpd = (int)(tmp/SAF_PI);
                    if (qpd >= 0) qpd += qpd&1;
                    else qpd -= qpd&1;
                    tmp -= SAF_PI*(float)qpd;

                    /* get deviation from bin frequency from the +/- Pi interval */
                    tmp = (float)h->osamp*tmp/(2.0f*SAF_PI);

                    /* compute the k-th partials' true frequency */
                    tmp = (float)k*freqPerBin + tmp*freqPerBin;

                    /* store magnitude and true frequency in analysis arrays */
                    h->gAnaMagn[ch][k] = magn;
                    h->gAnaFreq[ch][k] = tmp;

                }

                /* ***************** PROCESSING ******************* */
                /* this does the actual pitch shifting */
                memset(h->gSynMagn[ch], 0, h->fftFrameSize*sizeof(float));
                memset(h->gSynFreq[ch], 0, h->fftFrameSize*sizeof(float));
                for (k = 0; k <= fftFrameSize2; k++) {
                    index = (int)((float)k*(h->pitchShiftFactor));
                    if (index <= fftFrameSize2) {
                        h->gSynMagn[ch][index] += (h->gAnaMagn[ch][k]);
                        h->gSynFreq[ch][index] = h->gAnaFreq[ch][k] * (h->pitchShiftFactor);
                    }
                }

                /* ***************** SYNTHESIS ******************* */
                /* this is the synthesis step */
                for (k = 0; k <= fftFrameSize2; k++) {
                    /* get magnitude and true frequency from synthesis arrays */
                    magn = h->gSynMagn[ch][k];
                    tmp = h->gSynFreq[ch][k];

                    /* subtract bin mid frequency */
                    tmp -= (float)k*freqPerBin;

                    /* get bin deviation from freq deviation */
                    tmp /= freqPerBin;

                    /* take osamp into account */
                    tmp = 2.0f*SAF_PI*tmp/(h->osamp);

                    /* add the overlap phase advance back in */
                    tmp += (float)k*expct;

                    /* accumulate delta phase to get bin phase */
                    h->gSumPhase[ch][k] += tmp;
                    phase = h->gSumPhase[ch][k];

#ifdef SMB_ENABLE_SAF_FFT
                    h->gFFTworksp_fd[ch][k] = cmplxf(magn*cosf(phase), magn*sinf(phase));
#else
                    /* get real and imag part and re-interleave */
                    h->gFFTworksp[ch][2*k] = magn*cos(phase);
                    h->gFFTworksp[ch][2*k+1] = magn*sin(phase);
#endif
                }

#ifdef SMB_ENABLE_SAF_FFT
                for (k = fftFrameSize2+1; k < h->fftFrameSize; k++)
                    h->gFFTworksp_fd[ch][k] = cmplxf(0.0f, 0.0f);

                saf_fft_backward(h->hFFT, h->gFFTworksp_fd[ch], h->gFFTworksp_td[ch]);

                for(k=0; k < h->fftFrameSize; k++)
                    h->gOutputAccum[ch][k] += 2.0f*(h->window[k])*crealf(h->gFFTworksp_td[ch][k])/((h->osamp)); 
#else
                /* zero negative frequencies */
                for (k = h->fftFrameSize+2; k < 2*(h->fftFrameSize); k++)
                    h->gFFTworksp[ch][k] = 0.;

                /* do inverse transform */
                smbFft(h->gFFTworksp[ch], h->fftFrameSize, 1);

                /* do windowing and add to output accumulator */
                for(k=0; k < h->fftFrameSize; k++)
                    h->gOutputAccum[ch][k] += 2.*(h->window[k])*(h->gFFTworksp[ch][2*k])/(fftFrameSize2*(h->osamp));
#endif
                for (k = 0; k < h->stepSize; k++)
                    h->gOutFIFO[ch][k] = h->gOutputAccum[ch][k];

                /* shift accumulator */
                memmove(h->gOutputAccum[ch], h->gOutputAccum[ch] + (h->stepSize), h->fftFrameSize*sizeof(float));

                /* move input FIFO */
                for (k = 0; k < h->inFifoLatency; k++)
                    h->gInFIFO[ch][k] = h->gInFIFO[ch][k+h->stepSize];
             }
         }
    }
}
