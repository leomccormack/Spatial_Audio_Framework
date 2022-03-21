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
 * @file: saf_utility_filters.c
 * @ingroup Utilities
 * @brief A collection of IIR/FIR filter and filterbank designs
 *
 * @author Leo McCormack
 * @date 01.03.2019
 * @license ISC
 */
 
#include "saf_utilities.h"
#include "saf_externals.h"

/**
 * Applies a windowing function (see #WINDOWING_FUNCTION_TYPES enum) of length
 * 'winlength', to vector 'x'.
 */
static void applyWindowingFunction
(
    WINDOWING_FUNCTION_TYPES type,
    int winlength,
    float* x
)
{
    int i, N;
    
    /* if winlength is odd -> symmetric window (mid index has value=1) */
    if ( !(winlength % 2 == 0) )
        N = winlength-1;
    /* otherwise, if winlength is even (index: winlength/2+1 = 1.0, but first
     * value != last value) */
    else
        N = winlength;
    
    switch(type){
        case WINDOWING_FUNCTION_RECTANGULAR:
            break;
            
        case WINDOWING_FUNCTION_HAMMING:
            for(i=0; i<winlength; i++)
                x[i] *= 0.54f - 0.46f * (cosf(2.0f*SAF_PI*(float)i/(float)N)); /* more wide-spread coefficient values */
            /* optimal equiripple coefficient values: */
            /*x[i] *= 0.53836f - 0.46164f * (cosf(2.0f*SAF_PI*(float)i/(float)N));*/
            break;
            
        case WINDOWING_FUNCTION_HANN:
            for(i=0; i<winlength; i++)
                x[i] *= 0.5f - 0.5f * (cosf(2.0f*SAF_PI*(float)i/(float)N));
            break;
            
        case WINDOWING_FUNCTION_BARTLETT:
            for(i=0; i<winlength; i++)
                x[i] *= 1.0f - 2.0f * fabsf((float)i-((float)N/2.0f))/(float)N;
            break;
            
        case WINDOWING_FUNCTION_BLACKMAN:
            for(i=0; i<winlength; i++){
                x[i] *= 0.42659f -
                        0.49656f *cosf(2.0f*SAF_PI*(float)i/(float)N) +
                        0.076849f*cosf(4.0f*SAF_PI*(float)i/(float)N);
            }
            break;
            
        case WINDOWING_FUNCTION_NUTTALL:
            for(i=0; i<winlength; i++){
                x[i] *= 0.355768f -
                        0.487396f*cosf(2.0f*SAF_PI*(float)i/(float)N) +
                        0.144232f*cosf(4.0f*SAF_PI*(float)i/(float)N) -
                        0.012604f*cosf(6.0f*SAF_PI*(float)i/(float)N);
            }
            break;
            
        case WINDOWING_FUNCTION_BLACKMAN_NUTTALL:
            for(i=0; i<winlength; i++){
                x[i] *= 0.3635819f -
                        0.4891775f*cosf(2.0f*SAF_PI*(float)i/(float)N) +
                        0.1365995f*cosf(4.0f*SAF_PI*(float)i/(float)N) +
                        0.0106411f*cosf(4.0f*SAF_PI*(float)i/(float)N);
            }
            break;
            
        case WINDOWING_FUNCTION_BLACKMAN_HARRIS:
            for(i=0; i<winlength; i++){
                x[i] *= 0.35875f -
                        0.48829f*cosf(2.0f*SAF_PI*(float)i/(float)N) +
                        0.14128f*cosf(4.0f*SAF_PI*(float)i/(float)N) +
                        0.01168f*cosf(4.0f*SAF_PI*(float)i/(float)N);
            }
            break;
    }
}

/**
 * Applies IIR filter of order 1 */
static void applyIIR_1
(
    float* in_signal,
    int nSamples,
    float* b,
    float* a,
    float* wz,
    float* out_signal
)
{
    int n;
    float wn;

    /*  difference equation (Direct form 2) */
    for (n=0; n<nSamples; n++){
        /* numerator */
        wn = in_signal[n] - (a[1] * wz[0]);

        /* denominator */
        out_signal[n] = b[0] * wn + b[1] * wz[0];

        /* shuffle delays */
        wz[0] = wn;
    }
}

/**
 * Applies IIR filter of order 2 */
static void applyIIR_2
(
    float* in_signal,
    int nSamples,
    float* b,
    float* a,
    float* wz,
    float* out_signal
)
{
    int n;
    float wn;

    /* Difference equation (Direct form 2) */
    for (n=0; n<nSamples; n++){
        /* numerator */
        wn = in_signal[n] - a[1] * wz[0] - a[2] * wz[1];

        /* denominator */
        out_signal[n] = b[0] * wn + b[1] * wz[0] + b[2] * wz[1];

        /* shuffle delays */
        wz[1] = wz[0];
        wz[0] = wn;
    }
}

/**
 * Applies IIR filter of order 3 */
static void applyIIR_3
(
    float* in_signal,
    int nSamples,
    float* b,
    float* a,
    float* wz,
    float* out_signal
)
{
    int n;
    float wn;

    /* Difference equation (Direct form 2) */
    for (n=0; n<nSamples; n++){
        /* numerator */
        wn = in_signal[n] - a[1] * wz[0] - a[2] * wz[1] - a[3] * wz[2];

        /* denominator */
        out_signal[n] = b[0] * wn + b[1] * wz[0] + b[2] * wz[1] + b[3] * wz[2];

        /* shuffle delays */
        wz[2] = wz[1];
        wz[1] = wz[0];
        wz[0] = wn;
    }
}


/* ========================================================================== */
/*                               Misc. Functions                              */
/* ========================================================================== */

void getWindowingFunction
(
    WINDOWING_FUNCTION_TYPES type,
    int winlength,
    float* win
)
{
    int i;
    for(i=0; i<winlength; i++)
        win[i] = 1.0f;
    applyWindowingFunction(type, winlength, win);
}

void getOctaveBandCutoffFreqs
(
    float* centreFreqs,
    int nCentreFreqs,
    float* cutoffFreqs
)
{
    int i;
    for(i=0; i<nCentreFreqs-1; i++)
        cutoffFreqs[i] = 2.0f*centreFreqs[i]/sqrtf(2.0f);
}

void flattenMinphase
(
    float* x,
    int len
)
{
    int i;
    float_complex* ctd_tmp, *tdi_f, *tdi_f_labs, *dt_min_f;
    void* hFFT;
    
    /* prep */
    ctd_tmp = malloc1d(len*sizeof(float_complex));
    tdi_f = malloc1d(len*sizeof(float_complex));
    tdi_f_labs = malloc1d(len*sizeof(float_complex));
    dt_min_f = malloc1d(len*sizeof(float_complex));
    saf_fft_create(&hFFT, len);
    
    /* fft */
    for(i=0; i<len; i++)
        ctd_tmp[i] = cmplxf(x[i], 0.0f);
    saf_fft_forward(hFFT, (float_complex*)ctd_tmp, (float_complex*)tdi_f);
    
    /* take log(cabs()) */
    for(i=0; i<len; i++)
        tdi_f_labs[i] = cmplxf(logf(cabsf(tdi_f[i])), 0.0f);
    
    /* Hilbert to acquire discrete-time analytic signal */
    hilbert(tdi_f_labs, len, dt_min_f);
    
    /* compute minimum-phase response, and apply to tdi_f to flatten it to unity magnitude */
    for(i=0; i<len; i++)
        dt_min_f[i] = ccdivf(tdi_f[i], cexpf(conjf(dt_min_f[i])));
    
    /* ifft */
    saf_fft_backward(hFFT, dt_min_f, ctd_tmp);
    
    /* overwrite input with EQ'd version */
    for(i=0; i<len; i++)
        x[i] = crealf(ctd_tmp[i]);
    
    /* tidy up */
    saf_fft_destroy(&hFFT);
    free(ctd_tmp);
    free(tdi_f);
    free(tdi_f_labs);
    free(dt_min_f);
}

void interpolateFiltersH
(
    int inFFTsize,
    int outFFTsize,
    int nFilters,
    float_complex* filters_in,
    float_complex* filters_out
)
{
    int i, j, nBins_in, nBins_out;
    float* M_ifft, *M_ifft_fl;
    float_complex* tmp;
    void* hFFT_in, *hFFT_out;

    nBins_in = inFFTsize/2 + 1;
    nBins_out = outFFTsize/2 + 1;
    saf_rfft_create(&hFFT_in, inFFTsize);
    saf_rfft_create(&hFFT_out, outFFTsize);
    M_ifft    = calloc1d(SAF_MAX(inFFTsize, outFFTsize), sizeof(float));
    M_ifft_fl = calloc1d(SAF_MAX(inFFTsize, outFFTsize), sizeof(float));
    tmp = malloc1d(SAF_MAX(nBins_in, nBins_out) * sizeof(float_complex));

    for(i=0; i<nFilters; i++){
        for(j=0; j<nBins_in; j++)
            tmp[j] = filters_in[j*nFilters+i];
        saf_rfft_backward(hFFT_in, tmp, M_ifft);

        /* flip */
        for(j=0; j<outFFTsize/2; j++){
            M_ifft_fl[j] = M_ifft[inFFTsize/2+j];
            M_ifft_fl[inFFTsize/2+j] = M_ifft[j];
        }
        saf_rfft_forward(hFFT_out, M_ifft_fl, tmp);
        for(j=0; j<nBins_out; j++)
            filters_out[j*nFilters+i] = tmp[j];
    }

    saf_rfft_destroy(&hFFT_in);
    saf_rfft_destroy(&hFFT_out);
    free(M_ifft);
    free(M_ifft_fl);
    free(tmp);
}

float convertBW2Q
(
    float BW
)
{
    return sqrtf(powf(2.0f, BW))/(powf(2.0f, BW)-1.0f);
}

float convertQ2BW
(
    float Q
)
{
    return logf(  (2.0f*Q*Q+1.0f)/(2.0f*Q*Q) + sqrtf( powf((2.0f*Q*Q+1.0f)/(Q*Q+2.23e-13f), 2.0f)/4.0f - 1.0f )) /logf(2.0f);
}


/* ========================================================================== */
/*                             IIR Filter Functions                           */
/* ========================================================================== */

void biQuadCoeffs
(
    BIQUAD_FILTER_TYPES filterType,  
    float fc,
    float fs,
    float Q,
    float gain_dB,
    float b[3],
    float a[3] 
)
{
    float K, KK, D, V0, A, w0, alpha, a0;
    
    a[0] = 1.0f;
    
    /* calculate the IIR filter coefficients */
    switch (filterType){
        case BIQUAD_FILTER_LPF:
            /* Filter design equations - DAFX (2nd ed) p50 */
            K = tanf(SAF_PI * fc/fs);
            KK = K * K;
            D = KK * Q + K + Q;
            b[0] = (KK * Q)/D;
            b[1] = (2.0f * KK * Q)/D;
            b[2] = b[0];
            a[1] = (2.0f * Q * (KK - 1.0f))/D;
            a[2] = (KK * Q - K + Q)/D;
            break;

        case BIQUAD_FILTER_LPF_EQCB:
            /* Filter design equations - https://webaudio.github.io/Audio-EQ-Cookbook/audio-eq-cookbook.html */
            w0 = 2.0f*SAF_PI*fc/fs;
            alpha = sinf(w0)/(2.0f*Q);
            b[0] = (1.0f - cosf(w0))/2.0f;
            b[1] = 1.0f - cosf(w0);
            b[2] = b[0];
            a0   = 1.0f + alpha;
            a[1] = -2.0f*cosf(w0);
            a[2] = 1.0f - alpha;

            /* Scale by a0, since applyBiQuadFilter() and applyIIR() assume a0 = 1.0 */
            b[0] /= a0;
            b[1] /= a0;
            b[2] /= a0;
            a[1] /= a0;
            a[2] /= a0;
            break;
            
        case BIQUAD_FILTER_HPF:
            /* Filter design equations - DAFX (2nd ed) p50 */
            K = tanf(SAF_PI * fc/fs);
            KK = K * K;
            D = KK * Q + K + Q;
            b[0] = (Q)/D;
            b[1] = -(2.0f * Q)/D;
            b[2] = b[0];
            a[1] = (2.0f * Q * (KK - 1.0f))/D;
            a[2] = (KK * Q - K + Q)/D;
            break;

        case BIQUAD_FILTER_HPF_EQCB:
            /* Filter design equations - https://webaudio.github.io/Audio-EQ-Cookbook/audio-eq-cookbook.html */
            w0 = 2.0f*SAF_PI*fc/fs;
            alpha = sinf(w0)/(2.0f*Q);
            b[0] = (1.0f + cosf(w0))/2.0f;
            b[1] = -(1.0f + cosf(w0));
            b[2] = b[0];
            a0   = 1.0f + alpha;
            a[1] = -2.0f*cosf(w0);
            a[2] = 1.0f - alpha;

            /* Scale by a0, since applyBiQuadFilter() and applyIIR() assume a0 = 1.0 */
            b[0] /= a0;
            b[1] /= a0;
            b[2] /= a0;
            a[1] /= a0;
            a[2] /= a0;
            break;
            
        case BIQUAD_FILTER_LOW_SHELF:
            /* Filter design equations - DAFX (2nd ed) p64 */
            K = tanf(SAF_PI * fc/fs);
            V0 = powf(10.0f, (gain_dB/20.0f));
            if (V0 < 1.0f)
                V0 = 1.0f/V0;
            KK = K * K;
            if (gain_dB > 0.0f){
                D = 1.0f + sqrtf(2.0f) * K + KK;
                b[0] = (1.0f + sqrtf(2.0f * V0) * K + V0 * KK)/D;
                b[1] = (2.0f*(V0*KK - 1.0f))/D;
                b[2] = (1.0f - sqrtf(2.0f * V0) * K + V0 * KK)/D;
                a[1] = (2.0f * (KK - 1.0f))/D;
                a[2] = (1.0f - sqrtf(2.0f) * K + KK)/D;
            }
            else{
                D = V0 + sqrtf(2.0f*V0)*K + KK;
                b[0] = (V0*(1.0f + sqrtf(2.0f)*K + KK))/D;
                b[1] = (2.0f*V0*(KK - 1.0f))/D;
                b[2] = (V0*(1.0f - sqrtf(2.0f)*K + KK))/D;
                a[1] = (2.0f * (KK - V0))/D;
                a[2] = (V0 - sqrtf(2.0f*V0)*K + KK)/D;
            }
            break;

        case BIQUAD_FILTER_LOW_SHELF_EQCB:
            /* Filter design equations - https://webaudio.github.io/Audio-EQ-Cookbook/audio-eq-cookbook.html */
            A = powf(10.0f, gain_dB/40.0f);
            w0 = 2.0f*SAF_PI*fc/fs;
            alpha = sinf(w0)/(2.0f*Q);
            b[0] = A*( (A+1.0f) - (A-1.0f) * cosf(w0) + 2.0f*sqrtf(A)*alpha );
            b[1] = 2.0f*A*( (A-1.0f) - (A+1.0f) * cosf(w0) );
            b[2] = A*( (A+1.0f) - (A-1.0f) * cosf(w0) - 2.0f*sqrtf(A)*alpha );
            a0   = (A+1.0f) + (A-1.0f) * cosf(w0) + 2.0f*sqrtf(A)*alpha;
            a[1] = -2.0f*( (A-1.0f) + (A+1.0f) * cosf(w0) );
            a[2] = (A+1.0f) + (A-1.0f) * cosf(w0) - 2.0f*sqrtf(A)*alpha;

            /* Scale by a0, since applyBiQuadFilter() and applyIIR() assume a0 = 1.0 */
            b[0] /= a0;
            b[1] /= a0;
            b[2] /= a0;
            a[1] /= a0;
            a[2] /= a0;
            break;
            
        case BIQUAD_FILTER_HI_SHELF:
            /* Filter design equations - DAFX (2nd ed) p64 */
            K = tanf(SAF_PI * fc/fs);
            V0 = powf(10.0f, (gain_dB/20.0f));
            if (V0 < 1.0f)
                V0 = 1.0f/V0;
            KK = K * K;
            if (gain_dB > 0.0f){
                D = 1.0f + sqrtf(2.0f) * K + KK;
                b[0] = (V0 + sqrtf(2.0f * V0) * K + KK)/D;
                b[1] = (2.0f*(KK - V0))/D;
                b[2] = (V0 - sqrtf(2.0f * V0) * K + KK)/D;
                a[1] = (2.0f*(KK - 1.0f))/D;
                a[2] = (1.0f - sqrtf(2.0f) * K + KK)/D;
            }
            else{
                D = 1.0f + sqrtf(2.0f*V0) * K + V0*KK;
                b[0] = (V0*(1.0f + sqrtf(2.0f)*K + KK))/D;
                b[1] = (2.0f*V0*(KK - 1.0f))/D;
                b[2] = (V0*(1.0f - sqrtf(2.0f)*K + KK))/D;
                a[1] = (2.0f * (V0*KK - 1.0f))/D;
                a[2] = (1.0f - sqrtf(2.0f*V0)*K + V0*KK)/D;
            }
            break;

        case BIQUAD_FILTER_HI_SHELF_EQCB:
            /* Filter design equations - https://webaudio.github.io/Audio-EQ-Cookbook/audio-eq-cookbook.html */
            A = powf(10.0f, gain_dB/40.0f);
            w0 = 2.0f*SAF_PI*fc/fs;
            alpha = sinf(w0)/(2.0f*Q);
            b[0] = A*( (A+1.0f) + (A-1.0f) * cosf(w0) + 2.0f*sqrtf(A)*alpha );
            b[1] = -2.0f*A*( (A-1.0f) + (A+1.0f) * cosf(w0) );
            b[2] = A*( (A+1.0f) + (A-1.0f) * cosf(w0) - 2.0f*sqrtf(A)*alpha );
            a0   = (A+1.0f) - (A-1.0f) * cosf(w0) + 2.0f*sqrtf(A)*alpha;
            a[1] = 2.0f*( (A-1.0f) - (A+1.0f) * cosf(w0) );
            a[2] = (A+1.0f) - (A-1.0f) * cosf(w0) - 2.0f*sqrtf(A)*alpha;

            /* Scale by a0, since applyBiQuadFilter() and applyIIR() assume a0 = 1.0 */
            b[0] /= a0;
            b[1] /= a0;
            b[2] /= a0;
            a[1] /= a0;
            a[2] /= a0;
            break;
            
        case BIQUAD_FILTER_PEAK:
            /* Filter design equations - DAFX (2nd ed) p66 */
            K = tanf(SAF_PI * fc/fs);
            V0 = powf(10.0f, (gain_dB/20.0f));
            KK = K * K;
            if (gain_dB > 0.0f){
                D = 1.0f + (K/Q) + KK;
                b[0] = (1.0f + (V0/Q) * K + KK)/D;
                b[1] = (2.0f*(KK - 1.0f))/D;
                b[2] = (1.0f - (V0/Q) * K + KK)/D;
                a[1] = b[1];
                a[2] = (1.0f - (K/Q) + KK)/D;
            }
            else {
                D = 1.0f + (K/(V0*Q)) + KK;
                b[0] = (1.0f + (K/Q) + KK)/D;
                b[1] = (2.0f*(KK - 1.0f))/D;
                b[2] = (1.0f - (K/Q) + KK)/D;
                a[1] = b[1];
                a[2] = (1.0f - (K/(V0*Q)) + KK)/D;
            }
            break;

        case BIQUAD_FILTER_PEAK_EQCB:
            /* Filter design equations - https://webaudio.github.io/Audio-EQ-Cookbook/audio-eq-cookbook.html */
            A = powf(10.0f, gain_dB/40.0f);
            w0 = 2.0f*SAF_PI*fc/fs;
            alpha = sinf(w0)/(2.0f*Q);
            b[0] = 1.0f + alpha*A;
            b[1] = -2.0f*cosf(w0);
            b[2] = 1.0f - alpha*A;
            a0   = 1.0f + alpha/A;
            a[1] = b[1];
            a[2] = 1.0f - alpha/A;

            /* Scale by a0, since applyBiQuadFilter() and applyIIR() assume a0 = 1.0 */
            b[0] /= a0;
            b[1] /= a0;
            b[2] /= a0;
            a[1] /= a0;
            a[2] /= a0;
            break;
    }
}

void applyBiQuadFilter
(
     float b[3],
     float a[3],
     float w_z_12[2],
     float* signal,
     int nSamples
)
{
    int n;
    float wn;
    
    /* biquad difference equation (Direct form 2) */
    for(n=0; n<nSamples; n++){
        wn = signal[n] - a[1] * w_z_12[0] - a[2] * w_z_12[1];
        signal[n] = b[0] * wn + b[1] * w_z_12[0] + b[2] * w_z_12[1];
        /* shuffle delays */
        w_z_12[1] = w_z_12[0];
        w_z_12[0] = wn;
    }
}

void evalBiQuadTransferFunction
(
    float b[3],
    float a[3],
    float* freqs,
    int nFreqs,
    float fs,
    int mag2dB,
    float* magnitude,
    float* phase_rad
)
{
    int ff;
    float w, denom_real, denom_imag, num_real, num_imag;
    
    for(ff=0; ff<nFreqs; ff++){
        w = tanf(SAF_PI * freqs[ff]/fs);
        /* substituting Euler, z = e^(-jwn) = cos(wn) + j*sin(wn), into:
         * H(z) = (b0 + b1*z^(-1) + b2*z^(-2)) / (1 + a1*z^(-1) + a2*z^(-2)): */
        denom_real = 1.0f + a[1]*cosf(w) + a[2]*cosf(2.0f*w);
        denom_imag = a[1]*sinf(w) + a[2]*sinf(2.0f*w);
        num_real = b[0] + b[1]*cosf(w) + b[2]*cosf(2.0f*w);
        num_imag = b[1]*sinf(w) + b[2]*sinf(2.0f*w);
        
        if(magnitude!=NULL){
            magnitude[ff] = sqrtf( (powf(num_real, 2.0f) + powf(num_imag, 2.0f)) / (powf(denom_real, 2.0f) + powf(denom_imag, 2.0f) + 2.23e-7f) ); 
            if(mag2dB)
                magnitude[ff] = 20.0f*log10f(magnitude[ff]);
        }
        if(phase_rad!=NULL)
            phase_rad[ff] = atan2f(num_imag,num_real) - atan2f(denom_imag, denom_real);
    }
}

void evalIIRTransferFunctionf
 (
    float* b_coeff, /* Note coeffs are *floats* */
    float* a_coeff,
    int nCoeffs,
    float* freqs,
    int nFreqs,
    float fs,
    int mag2dB,
    float* magnitude,
    float* phase_rad
)
{
    int ff;
    float w, x, a, b, c, d, h_re, h_im, cosx, sinx;
    float norm_frq = -2.0 * SAF_PI / fs; // -1 factored in from negated 'n' below

    /*
     * substitute Euler's  z = e^(jwn) = cos(wn) + j*sin(wn)
     * into
     *        [b0*z^(0) + b1*z^(-1) ... + bn*z^(-n)]
     * H(z) = --------------------------------------
     *        [a0*z^(0) + a1*z^(-1) ... + an*z^(-n)]
     */
    
    for(ff = 0; ff < nFreqs; ff++){
        w = freqs[ff] * norm_frq;
        
        // n = 0; real terms are b[0] and a[0], imag terms are 0
        a = b_coeff[0];    // num_re;   b_coeff[0] * cos(x); x = -1.f * n * w;
        b = 0.0;           // num_imag; b_coeff[0] * sin(x);
        c = a_coeff[0];    // den_re
        d = 0.0;           // den_imag
        
        /* sum over remaining numerator and denominator terms */
        for(int n = 1; n < nCoeffs; n++){
            x = (n * w);    // n * -1.f * w (-1.f applied in norm_frq)
            cosx = cosf(x); // alt: cosx = 1 - 2.f * powf(sinf(x/2.f), 2.f);
            sinx = sinf(x);
            a += b_coeff[n] * cosx;  // 'a'
            b += b_coeff[n] * sinx;  // 'b'
            c += a_coeff[n] * cosx;  // 'c'
            d += a_coeff[n] * sinx;  // 'd'
        }

        /* 1 / (c^2 + d^2 + eps) */
        double dvsr = 1.0 / (powf(c, 2.f) + powf(d, 2.f) + 2.23e-7f);
        
        if(magnitude!=NULL){
            /* sqrt((a^2 + b^2) / (c^2 + d^2)) */
            magnitude[ff] = (float)sqrt( (powf(a, 2.0f) + powf(b, 2.0f)) * dvsr );
            if(mag2dB)
                magnitude[ff] = 20.0f*log10f(magnitude[ff]);
        }
        
        if(phase_rad!=NULL) {
            /* complex division */
            h_re = (a*c + b*d) * dvsr;
            h_im = (b*c - a*d) * dvsr;
            phase_rad[ff] = (float)atan2(h_im, h_re);
        }
    }
}

void evalIIRTransferFunction
 (
    double* b_coeff, /* Note coeffs are *doubles* */
    double* a_coeff,
    int nCoeffs,
    float* freqs,
    int nFreqs,
    float fs,
    int mag2dB,
    float* magnitude,
    float* phase_rad
)
{
    int ff;
    float w;
    double x, a, b, c, d, h_re, h_im, cosx, sinx;
    float norm_frq = -2.0 * SAF_PI / fs; // -1 factored in from negated 'n' below
    
    /*
     * substitute Euler's  z = e^(jwn) = cos(wn) + j*sin(wn)
     * into
     *        [b0*z^(0) + b1*z^(-1) ... + bn*z^(-n)]
     * H(z) = --------------------------------------
     *        [a0*z^(0) + a1*z^(-1) ... + an*z^(-n)]
     */
    
    for(ff = 0; ff < nFreqs; ff++){
        w = freqs[ff] * norm_frq;
        
        // n = 0; real terms are b[0] and a[0], imag terms are 0
        a = b_coeff[0];    // num_re;   b_coeff[0] * cos(x); x = -1.f * n * w;
        b = 0.0;           // num_imag; b_coeff[0] * sin(x);
        c = a_coeff[0];    // den_re
        d = 0.0;           // den_imag
        
        /* sum over remaining numerator and denominator terms */
        for(int n = 1; n < nCoeffs; n++){
            x = (double)(n * w);     // n * -1.f * w (-1.f applied in norm_frq)
            cosx = 1 - 2.f * pow(sin(x/2.f), 2.f); // cos(x) = 1 - 2sin^2(x/2)
            sinx = sin(x);
            a += b_coeff[n] * cosx;  // 'a'
            b += b_coeff[n] * sinx;  // 'b'
            c += a_coeff[n] * cosx;  // 'c'
            d += a_coeff[n] * sinx;  // 'd'
        }

        /* 1 / (c^2 + d^2 + eps) */
        double dvsr = 1.0 / (pow(c, 2.f) + pow(d, 2.f) + 2.23e-17f);

        if(magnitude!=NULL){
            /* sqrt((a^2 + b^2) / (c^2 + d^2)) */
            magnitude[ff] = (float)sqrt( (pow(a, 2.0f) + pow(b, 2.0f)) * dvsr );
            if(mag2dB)
                magnitude[ff] = 20.0f*log10f(magnitude[ff]);
        }
        
        if(phase_rad!=NULL) {
            /* complex division */
            h_re = (a*c + b*d) * dvsr;
            h_im = (b*c - a*d) * dvsr;
            phase_rad[ff] = (float)atan2(h_im, h_re);
        }
    }
}


void applyIIR
(
    float* in_signal,
    int nSamples,
    int nCoeffs,
    float* b,
    float* a,
    float* wz,
    float* out_signal
)
{
    int n, i;
    float wn;
    
#if defined(SAF_USE_INTEL_IPP) && 0 /* Couldn't get this to give the same/correct result... */
    int pBufSize;
    Ipp32f* ba;
    IppsIIRState_32f* pIppIIR;
    Ipp8u * m_pBuf;
    ba = malloc1d(nCoeffs*2*sizeof(Ipp32f));
    for(i=0; i<nCoeffs; i++){
        ba[i] = b[i];
        ba[i+nCoeffs] = a[i];
    }
    ippsIIRGetStateSize_32f( nCoeffs-1, &pBufSize );
    m_pBuf = ippsMalloc_8u(pBufSize);

    pIppIIR = NULL;
    IppStatus error = ippsIIRInit_32f(&pIppIIR, ba, nCoeffs-1, NULL, m_pBuf);
    error = ippsIIR_32f( in_signal, out_signal, nSamples, pIppIIR );

    free(ba);
    ippsFree(m_pBuf);

    return;
#endif

    /* For compiler speed-ups */  
    switch(nCoeffs){
        case 1: saf_print_error("Just divide in_signal by b[0]...");
        case 2: applyIIR_1(in_signal, nSamples, b, a, wz, out_signal); return;
        case 3: applyIIR_2(in_signal, nSamples, b, a, wz, out_signal); return;
        case 4: applyIIR_3(in_signal, nSamples, b, a, wz, out_signal); return;
    }

    /*  difference equation (Direct form 2) */
    for (n=0; n<nSamples; n++){
        /* numerator */
        wn = in_signal[n];
        for (i=1; i<nCoeffs;i++)
            wn = wn - (a[i] * wz[i-1]);

        /* denominator */
        out_signal[n] = b[0] * wn;
        for (i=1; i<nCoeffs; i++)
            out_signal[n] = out_signal[n] + (b[i] * wz[i-1]);

        /* shuffle delays */
        switch(nCoeffs-1){
            case 10: wz[9] = wz[8]; /* fall through */
            case 9:  wz[8] = wz[7]; /* fall through */
            case 8:  wz[7] = wz[6]; /* fall through */
            case 7:  wz[6] = wz[5]; /* fall through */
            case 6:  wz[5] = wz[4]; /* fall through */
            case 5:  wz[4] = wz[3]; /* fall through */
            case 4:  wz[3] = wz[2]; /* fall through */
            case 3:  wz[2] = wz[1]; /* fall through */
            case 2:  wz[1] = wz[0]; /* fall through */
            case 1:  wz[0] = wn; break;
            default: saf_print_error("Unsupported number of IIR filter coefficients.");
        }
    }
}

void butterCoeffs
(
    BUTTER_FILTER_TYPES filterType,
    int order,
    float cutoff1,
    float cutoff2,
    float sampleRate,
    double* b_coeffs,
    double* a_coeffs
)
{
    int i, j, k, np, tmp_len, numStates, oddPolesFLAG, nCoeffs;
    double wlow, whi, w0, wl, w1, BW, Wn1, q;
    double den[3];
    double* c_state, *r, *b_coeffs_real;
    double** a_state, **bf_ss, **tmp1, **tmp2, **a_bili;
    double_complex kaT, kbT;
    double_complex den_cmplx[3];
    double_complex* proto, *proto_tmp, *a_coeffs_cmplx, *kern, *rcmplx, *b_coeffs_cmplx;

    wlow = (double)cutoff1/((double)sampleRate/2.0);
    whi = (double)cutoff2/((double)sampleRate/2.0);
    w0 = 4.0 * tan(SAF_PI*wlow/2.0);
    Wn1 = 0.0;

    /* Compute prototype for Nth order Butterworth analogue lowpass filter */
    if (order%2 != 0){/* ISODD */
        tmp_len = (int)((float)order/2.0f); /* floor */
        np = 2*tmp_len+1;
        proto = malloc1d(np*sizeof(double_complex));
        proto[np-1] = cmplx(-1.0,0.0);
    }
    else{ /* ISEVEN */
        tmp_len = order/2;
        np = 2*tmp_len;
        proto = malloc1d(np*sizeof(double_complex));
    }
    proto_tmp = malloc1d(np*sizeof(double_complex));
    for(i=1, j=0; i<=order-1; i+=2, j++)
        proto_tmp[j] = cexp(cmplx(0.0, SAF_PI*(double)i/(2.0*(double)order) + SAF_PI/2.0) );
    for (i=0; i<tmp_len; i++){
        proto[2*i] = proto_tmp[i];
        proto[2*i+1] = conj(proto_tmp[i]);
    }

    /* Transform prototype into state space  */
    numStates = np;
    cmplxPairUp(proto, proto_tmp, np);
    memcpy(proto, proto_tmp, np*sizeof(double_complex));
    free(proto_tmp);
    a_state = (double**)calloc2d(numStates,numStates,sizeof(double));
    c_state = malloc1d(numStates*sizeof(double));
    if (np%2 != 0){/* ISODD */
        a_state[0][0] = creal(proto[np-1]);
        c_state[0] = 1.0;
        np--;
        oddPolesFLAG = 1;
    }
    else
        oddPolesFLAG = 0;

    /* Adjust indices as needed */
    for(i=1; i<np; i+=2){
        polyz_v(&proto[i-1], den_cmplx, 2);
        for(j=0; j<3; j++)
            den[j] = creal(den_cmplx[j]);
        j = oddPolesFLAG ? i-1 : i-2;

        if(j==-1){
            a_state[0][0] = -den[1];
            a_state[0][1] = -den[2];
            a_state[1][0] = 1.0;
            a_state[1][1] = 0.0;
            c_state[0] = 0.0;
            c_state[1] = 1.0;
        }
        else{
            for(k=0; k<j+1; k++)
                a_state[j+1][k] = c_state[k];
            a_state[j+1][j+1] = -den[1];
            a_state[j+1][j+2] = -den[2];
            a_state[j+2][j+1] = 1.0;
            a_state[j+2][j+2] = 0.0;

            for(k=0; k<j+1; k++)
                c_state[k] = 0.0;
            c_state[j+1] = 0.0;
            c_state[j+2] = 1.0;
        }
    }

    /* Transform lowpass filter into the desired filter (while in state space) */
    bf_ss = NULL;
    switch(filterType){
        case BUTTER_FILTER_HPF:
            utility_dinv(NULL, FLATTEN2D(a_state), FLATTEN2D(a_state), numStates);
            /* fall through */
        case BUTTER_FILTER_LPF:
            bf_ss = (double**)malloc2d(numStates,numStates,sizeof(double));
            for(i=0; i<numStates; i++)
                for(j=0; j<numStates; j++)
                    bf_ss[i][j] = w0*(a_state[i][j]);
            break;
        case BUTTER_FILTER_BSF:
            utility_dinv(NULL, FLATTEN2D(a_state), FLATTEN2D(a_state), numStates);
            /* fall through */
        case BUTTER_FILTER_BPF:
            numStates = numStates*2;
            w1 = 4.0*tan(SAF_PI*whi/2.0);
            BW = w1 - w0;
            Wn1 = sqrt(w0*w1);
            q = Wn1/BW;
            bf_ss = (double**)calloc2d(numStates,numStates,sizeof(double));
            for(i=0; i<numStates/2; i++)
                for(j=0; j<numStates/2; j++)
                    bf_ss[i][j] = Wn1 * (a_state[i][j]) /q;
            for(i=numStates/2; i<numStates; i++)
                for(j=0; j<numStates/2; j++)
                    bf_ss[i][j] = (i-numStates/2) == j ? -Wn1 : 0.0;
            for(i=0; i<numStates/2; i++)
                for(j=numStates/2; j<numStates; j++)
                    bf_ss[i][j] = i == (j-numStates/2) ? Wn1 : 0.0;
            break;
    } 
    nCoeffs = numStates+1;

    /* Bilinear transformation to find the discrete equivalent of the filter */
    tmp1 = (double**)malloc2d(numStates,numStates,sizeof(double));
    tmp2 = (double**)malloc2d(numStates,numStates,sizeof(double));
    a_bili = (double**)malloc2d(numStates,numStates,sizeof(double));
    for(i=0; i<numStates; i++){
        for(j=0; j<numStates; j++){
            tmp1[i][j] = (i==j ? 1.0f : 0.0f) + bf_ss[i][j]*0.25;
            tmp2[i][j] = (i==j ? 1.0f : 0.0f) - bf_ss[i][j]*0.25;
        }
    }
    utility_dglslv(NULL, FLATTEN2D(tmp2), numStates, FLATTEN2D(tmp1), numStates, FLATTEN2D(a_bili));

    /* Compute the filter coefficients for the numerator and denominator */
    a_coeffs_cmplx = malloc1d(nCoeffs*sizeof(double_complex));
    polyd_m(FLATTEN2D(a_bili), a_coeffs_cmplx, numStates);
    rcmplx = NULL;
    r = NULL;
    switch(filterType){
        default:
            /* fall through */
        case BUTTER_FILTER_LPF:
            r = malloc1d(numStates*sizeof(double));
            for(i=0; i<numStates; i++)
                r[i] = -1.0;
            wl = 0.0;
            break;
        case BUTTER_FILTER_HPF:
            r = malloc1d(numStates*sizeof(double));
            for(i=0; i<numStates; i++)
                r[i] = 1.0;
            wl = SAF_PI;
            break;
        case BUTTER_FILTER_BPF:
            r = malloc1d(numStates*sizeof(double));
            wl = 2.0*atan2(Wn1, 4.0);
            for(i=0; i<order;i++)
                r[i] = 1.0;
            for(; i<2*order;i++)
                r[i] = -1.0;
            break;
        case BUTTER_FILTER_BSF:
            rcmplx = malloc1d(numStates*sizeof(double_complex));
            Wn1 = 2.0*atan2(Wn1,4.0);
            wl = 0.0;
            for(i=0; i<numStates;i++)
                rcmplx[i] = cexp(cmplx(0.0, Wn1*pow(-1.0,(double)i)));
            break;
    }
    b_coeffs_real = malloc1d(nCoeffs*sizeof(double));
    if(filterType == BUTTER_FILTER_BSF){
        b_coeffs_cmplx = malloc1d(nCoeffs*sizeof(double_complex));
        polyz_v(rcmplx, b_coeffs_cmplx, numStates);
        for(i=0; i<nCoeffs; i++)
            b_coeffs_real[i] = creal(b_coeffs_cmplx[i]);
        free(b_coeffs_cmplx);
    }
    else
        polyd_v(r, b_coeffs_real, numStates);
    kern = calloc1d(nCoeffs,sizeof(double_complex));
    kaT = cmplx(0.0,0.0);
    kbT = cmplx(0.0,0.0);
    for(i=0; i<nCoeffs; i++){
        kern[i] = cexp(cmplx(0.0,-wl*(double)i));
        kaT = ccadd(kaT, crmul(kern[i],creal(a_coeffs_cmplx[i])));
        kbT = ccadd(kbT, crmul(kern[i],b_coeffs_real[i]));
    }

    /* output */
    for(i=0; i<nCoeffs; i++){
        b_coeffs[i] = creal(crmul(ccdiv(kaT,kbT), b_coeffs_real[i]));
        a_coeffs[i] = creal(a_coeffs_cmplx[i]);
    }

    /* clean-up */
    free(proto);
    free(a_state);
    free(c_state);
    free(bf_ss);
    free(tmp1);
    free(tmp2);
    free(a_bili);
    free(a_coeffs_cmplx);
    free(b_coeffs_real);
    free(kern);
}

/** Main structure for the Favrot&Faller filterbank */
typedef struct _faf_IIRFB_data{
    int nBands;       /**< Number of bands in the filterbank */
    int nFilters;     /**< Number of filters used by the filterbank */
    int filtLen;      /**< Filter length */
    int filtOrder;    /**< Filter order (must be 1 or 3) */
    int maxNSamplesToExpect; /**< Maximum number of samples to expect to process
                              *   at a time */
    float** b_lpf;    /**< Numerator filter coeffs for low-pass filters */
    float** a_lpf;    /**< Denominator filter coeffs for low-pass filters */
    float** b_hpf;    /**< Numerator filter coeffs for high-pass filters */
    float** a_hpf;    /**< Denominator filter coeffs for high-pass filters */
    float*** wz_lpf;  /**< Delay buffers for low-pass filters */
    float*** wz_hpf;  /**< Delay buffers for high-pass filters */
    float*** wz_apf1; /**< Delay buffers for all-pass filter part 1 */
    float*** wz_apf2; /**< Delay buffers for all-pass filter part 2 */
    float* tmp;       /**< Temporary buffer; maxNSamplesToExpect x 1 */
    float* tmp2;      /**< Temporary buffer; maxNSamplesToExpect x 1 */

}faf_IIRFB_data;

void faf_IIRFilterbank_create
(
    void** phFaF,
    int order,
    float* fc,
    int nCutoffFreq,
    float sampleRate,
    int maxNumSamples
)
{
    *phFaF = malloc1d(sizeof(faf_IIRFB_data));
    faf_IIRFB_data *fb = (faf_IIRFB_data*)(*phFaF);
    double b_lpf[4], a_lpf[4], b_hpf[4], a_hpf[4], r[7], revb[4], reva[4], q[4];
    double tmp[7], tmp2[7];
    double_complex d1[3], d2[3], d1_num[3], d2_num[3];
    double_complex z[3], A[3][3], ztmp[7], ztmp2[7];
    int i, j, f, filtLen, d1_len, d2_len;

    saf_assert( (order==1) || (order==3), "Only odd number orders are supported, and 5th order+ is numerically unstable");
    saf_assert(nCutoffFreq>1, "Number of filterbank cut-off frequencies must be more than 1");
    filtLen = order + 1;
    fb->filtOrder = order;
    fb->filtLen = filtLen;

    /* Number of bands is always one more than the number of cut-off
     * frequencies */
    fb->nFilters = nCutoffFreq;
    fb->nBands = nCutoffFreq + 1;

    /* Allocate memory for filter coefficients and delay buffers */
    fb->b_hpf = (float**)malloc2d(nCutoffFreq, filtLen, sizeof(float));
    fb->a_hpf = (float**)malloc2d(nCutoffFreq, filtLen, sizeof(float));
    fb->b_lpf = (float**)malloc2d(nCutoffFreq, filtLen, sizeof(float));
    fb->a_lpf = (float**)malloc2d(nCutoffFreq, filtLen, sizeof(float));
    fb->wz_hpf = (float***)calloc3d(fb->nBands, nCutoffFreq, order, sizeof(float));
    fb->wz_lpf = (float***)calloc3d(fb->nBands, nCutoffFreq, order, sizeof(float));
    fb->wz_apf1 = (float***)calloc3d(fb->nBands, nCutoffFreq, order, sizeof(float));
    fb->wz_apf2 = (float***)calloc3d(fb->nBands, nCutoffFreq, order, sizeof(float));
    fb->maxNSamplesToExpect = maxNumSamples;
    fb->tmp = malloc1d(maxNumSamples*sizeof(float));
    fb->tmp2 = malloc1d(maxNumSamples*sizeof(float));

    /* Compute low-pass and complementary high-pass filter coefficients for each
     * cut-off frequency */
    for(f=0; f<nCutoffFreq; f++){
        /* Low-pass filter */
        butterCoeffs(BUTTER_FILTER_LPF, order, fc[f], 0.0f, sampleRate, (double*)b_lpf, (double*)a_lpf);

        /* IIR power complementary filter design (i.e. High-pass) */
        for(i=0; i<filtLen; i++){
            reva[i] = a_lpf[filtLen-i-1];
            revb[i] = b_lpf[filtLen-i-1];
        }
        convd(revb, b_lpf, filtLen, filtLen, tmp);
        convd(a_lpf, reva, filtLen, filtLen, tmp2);
        for(i=0; i<2*filtLen-1; i++)
            r[i] = tmp[i] - tmp2[i];
        q[0] = sqrt(-r[0]/-1.0);
        q[1] = -r[1]/(2.0*-1.0*q[0]);
        if(order==3){
            //q[3]=conj(-1.0*q[0]);
            //q[2]=conj(-1.0*q[1]);
            q[3] = -1.0*q[0];
            q[2] = -1.0*q[1];
        }
        for(i=0; i<filtLen; i++)
            q[i] =  b_lpf[i] - q[i];

        /* Find roots of polynomial  */
        if(order==1)
            z[0] = cmplx(-q[1]/q[0], 0.0);
        else if(order==3){
            memset(A, 0, 9*sizeof(double_complex));
            A[0][0] = cmplx(-q[1]/q[0], 0.0);
            A[0][1] = cmplx(-q[2]/q[0], 0.0);
            A[0][2] = cmplx(-q[3]/q[0], 0.0);
            A[1][0] = cmplx(1.0, 0.0);
            A[2][1] = cmplx(1.0, 0.0);
            utility_zeig(NULL, (double_complex*)A, 3, NULL, NULL, NULL, (double_complex*)z);
        }

        /* Separate the zeros inside the unit circle and the ones outside to
         * form the allpass functions */
        d1[0] = cmplx(1.0, 0.0);
        d2[0] = cmplx(1.0, 0.0);
        d1_len = d2_len = 1;
        for(i=0; i<order; i++){
            if (cabs(z[i]) < 1.0){
                ztmp[0] = cmplx(1.0, 0.0);
                ztmp[1] = crmul(z[i], -1.0);
                convz(d2,ztmp,d2_len,2,ztmp2);
                d2_len++;
                for(j=0; j<d2_len; j++)
                    d2[j] = ztmp2[j];
            }
            else{
                ztmp[0] = cmplx(1.0, 0.0);
                ztmp[1] = ccdiv(cmplx(-1.0, 0.0), conj(z[i]));
                convz(d1,ztmp,d1_len,2,ztmp2);
                d1_len++;
                for(j=0; j<d1_len; j++)
                    d1[j] = ztmp2[j];
            }
        }

        /* Convert coupled allpass filter to transfer function form (code from:
         * https://github.com/nsk1001/Scilab-functions written by Nagma Samreen
         * Khan) */
        for(i=0; i<d1_len; i++)
            d1_num[i] = conj(d1[d1_len-i-1]);
        for(i=0; i<d2_len; i++)
            d2_num[i] = conj(d2[d2_len-i-1]);
        convz(d1_num, d2, d1_len, d2_len, ztmp);
        convz(d2_num, d1, d2_len, d1_len, ztmp2);
        for(i=0; i<filtLen; i++){
            b_hpf[i] = -0.5 * creal(ccsub(ztmp[filtLen-i-1], ztmp2[filtLen-i-1]));
            a_hpf[i] = a_lpf[i];
        }

        /* Store in single precision for run-time */
        for(i=0; i<filtLen; i++){
            fb->b_hpf[f][i] = (float)b_hpf[i];
            fb->a_hpf[f][i] = (float)a_hpf[i];
            fb->b_lpf[f][i] = (float)b_lpf[i];
            fb->a_lpf[f][i] = (float)a_lpf[i];
        }
    }
}

void faf_IIRFilterbank_apply
(
    void* hFaF,
    float* inSig,
    float** outBands,
    int nSamples
)
{
    faf_IIRFB_data *fb = (faf_IIRFB_data*)(hFaF);
    int band,j;

    saf_assert(nSamples <= fb->maxNSamplesToExpect, "Number of samples exceeds the maximum number informed when calling faf_IIRFilterbank_create()");

    /* Copy input signal to all output bands/channels */
    for(band=0; band<fb->nBands; band++)
        memcpy(outBands[band], inSig, nSamples*sizeof(float));

    /* Band 0 */
    for (j = 0; j<fb->nFilters; j++)
        applyIIR(outBands[0], nSamples, fb->filtLen, fb->b_lpf[j], fb->a_lpf[j], fb->wz_lpf[0][j], outBands[0]);

    /* Band 1 */
    applyIIR(outBands[1], nSamples, fb->filtLen, fb->b_hpf[0], fb->a_hpf[0], fb->wz_hpf[1][0], outBands[1]);
    for (j = 1; j<fb->nFilters; j++)
        applyIIR(outBands[1], nSamples, fb->filtLen, fb->b_lpf[j], fb->a_lpf[j], fb->wz_lpf[1][j], outBands[1]);

    /* All-pass filters (bands 2..N-1) */
    for (band = 2; band < fb->nBands; band++){
        for (j=0; j<=band-2; j++){
            applyIIR(outBands[band], nSamples, fb->filtLen, fb->b_lpf[j], fb->a_lpf[j], fb->wz_apf1[band][j], fb->tmp);
            applyIIR(outBands[band], nSamples, fb->filtLen, fb->b_hpf[j], fb->a_hpf[j], fb->wz_apf2[band][j], fb->tmp2);
            utility_svvadd(fb->tmp, fb->tmp2, nSamples, outBands[band]);
        }
    }

    /* Bands 2..N-2 */
    for(band = 2; band< fb->nBands-1; band++){
        /* high-pass filter */
        applyIIR(outBands[band], nSamples, fb->filtLen, fb->b_hpf[band-1], fb->a_hpf[band-1], fb->wz_hpf[band][band-1], outBands[band]);

        /* low-pass filters */
        for(j=band; j<fb->nBands-1; j++)
            applyIIR(outBands[band], nSamples, fb->filtLen, fb->b_lpf[j], fb->a_lpf[j], fb->wz_lpf[band][j], outBands[band]);
    }

    /* Band N-1 */
    if (fb->nBands>2){
        applyIIR(outBands[fb->nBands-1], nSamples, fb->filtLen, fb->b_hpf[fb->nFilters-1], fb->a_hpf[fb->nFilters-1],
                 fb->wz_hpf[fb->nBands-1][fb->nFilters-1], outBands[fb->nBands-1]);
    }
}

void faf_IIRFilterbank_flushBuffers
(
    void* hFaF
)
{
    faf_IIRFB_data *fb = (faf_IIRFB_data*)(hFaF);

    memset(FLATTEN3D(fb->wz_hpf),  0, (fb->nBands) * (fb->nFilters) * (fb->filtOrder) * sizeof(float));
    memset(FLATTEN3D(fb->wz_lpf),  0, (fb->nBands) * (fb->nFilters) * (fb->filtOrder) * sizeof(float));
    memset(FLATTEN3D(fb->wz_apf1), 0, (fb->nBands) * (fb->nFilters) * (fb->filtOrder) * sizeof(float));
    memset(FLATTEN3D(fb->wz_apf2), 0, (fb->nBands) * (fb->nFilters) * (fb->filtOrder) * sizeof(float));
}

void faf_IIRFilterbank_destroy
(
    void** phFaF
)
{
    faf_IIRFB_data *fb = (faf_IIRFB_data*)(*phFaF);

    if(fb!=NULL){
        free(fb->b_hpf);
        free(fb->a_hpf);
        free(fb->b_lpf);
        free(fb->a_lpf);
        free(fb->wz_lpf);
        free(fb->wz_hpf);
        free(fb->wz_apf1);
        free(fb->wz_apf2);
        free(fb->tmp);
        free(fb->tmp2);
        fb=NULL;
        *phFaF = NULL;
    }
}


/* ========================================================================== */
/*                            FIR Filter Functions                            */
/* ========================================================================== */

void FIRCoeffs
(
    FIR_FILTER_TYPES filterType,
    int order,
    float fc1,
    float fc2, /* only needed for band-pass/stop */
    float fs,
    WINDOWING_FUNCTION_TYPES windowType,
    int scalingFLAG,
    float* h_filt
)
{
    int i, h_len;
    float ft1, ft2, h_sum, f0;
    float_complex h_z_sum;
    
    h_len = order + 1;
    ft1 = fc1/(fs*2.0f);
    
    /* compute filter weights */
    if(order % 2 == 0){
        /* if order is multiple of 2 */
        switch(filterType){
            case FIR_FILTER_LPF:
                for(i=0; i<h_len; i++)
                    h_filt[i] = i==order/2 ? 2.0f*ft1 : sinf(2.0f*SAF_PI*ft1*(float)(i-order/2)) / (SAF_PI*(float)(i-order/2));
                break;
                
            case FIR_FILTER_HPF:
                for(i=0; i<h_len; i++)
                    h_filt[i] = i==order/2 ? 1.0f - 2.0f*ft1 : -sinf(2.0f*ft1*SAF_PI*(float)(i-order/2)) / (SAF_PI*(float)(i-order/2));
                break;
                
            case FIR_FILTER_BPF:
                ft2 = fc2/(fs*2.0f);
                for(i=0; i<h_len; i++){
                    h_filt[i] = i==order/2 ? 2.0f*(ft2-ft1) :
                        sinf(2.0f*SAF_PI*ft2*(float)(i-order/2)) / (SAF_PI*(float)(i-order/2)) - sinf(2.0f*SAF_PI*ft1*(float)(i-order/2)) / (SAF_PI*(float)(i-order/2));
                }
                break;
                
            case FIR_FILTER_BSF:
                ft2 = fc2/(fs*2.0f);
                for(i=0; i<h_len; i++){
                    h_filt[i] = i==order/2 ? 1.0f - 2.0f*(ft2-ft1) :
                        sinf(2.0f*SAF_PI*ft1*(float)(i-order/2)) / (SAF_PI*(float)(i-order/2)) - sinf(2.0f*SAF_PI*ft2*(float)(i-order/2)) / (SAF_PI*(float)(i-order/2));
                }
                break;
        }
    }
    else
        saf_print_error("Please specify an even value for the filter 'order' argument"); 
    
    /* Apply windowing function */
    applyWindowingFunction(windowType, h_len, h_filt);
    
    /* Scaling, to ensure pass-band is truely at 1 (0dB).
     * [1] "Programs for Digital Signal Processing", IEEE Press John Wiley &
     *     Sons, 1979, pg. 5.2-1.
     */
    if(scalingFLAG){ 
        switch(filterType){
            case FIR_FILTER_LPF:
            case FIR_FILTER_BSF:
                h_sum = 0.0f;
                for(i=0; i<h_len; i++)
                    h_sum += h_filt[i];
                for(i=0; i<h_len; i++)
                    h_filt[i] /= h_sum;
                break;
                
            case FIR_FILTER_HPF:
                f0 = 1.0f;
                h_z_sum = cmplxf(0.0f, 0.0f);
                for(i=0; i<h_len; i++)
                    h_z_sum = ccaddf(h_z_sum, crmulf(cexpf(cmplxf(0.0f, -2.0f*SAF_PI*(float)i*f0/2.0f)), h_filt[i]));
                h_sum = cabsf(h_z_sum);
                for(i=0; i<h_len; i++)
                    h_filt[i] /= h_sum;
                break;
                
            case FIR_FILTER_BPF:
                f0 = (fc1/fs+fc2/fs)/2.0f;
                h_z_sum = cmplxf(0.0f, 0.0f);
                for(i=0; i<h_len; i++)
                    h_z_sum = ccaddf(h_z_sum, crmulf(cexpf(cmplxf(0.0f, -2.0f*SAF_PI*(float)i*f0/2.0f)), h_filt[i]));
                h_sum = cabsf(h_z_sum);
                for(i=0; i<h_len; i++)
                    h_filt[i] /= h_sum;
                break;
        }
    }
}

void FIRFilterbank
(
    int order,
    float* fc,  /* cut-off frequencies; nCutoffFreq x 1 */
    int nCutoffFreq,
    float sampleRate,
    WINDOWING_FUNCTION_TYPES windowType,
    int scalingFLAG,
    float* filterbank /* (nCutoffFreq+1) x (order+1) */
)
{
    int k, nFilt;
    
    /* Number of filters returned is always one more than the number of cut-off frequencies */
    nFilt = nCutoffFreq + 1;
    
    /* first and last bands are low-pass and high pass filters, using the first
     * and last cut-off frequencies in vector 'fc', respectively.  */
    FIRCoeffs(FIR_FILTER_LPF, order, fc[0], 0.0f, sampleRate, windowType, scalingFLAG, filterbank);
    FIRCoeffs(FIR_FILTER_HPF, order, fc[nCutoffFreq-1], 0.0f, sampleRate, windowType, scalingFLAG, &filterbank[(nFilt-1)*(order+1)]);
    
    /* the inbetween bands are then band-pass filters: */
    if(nCutoffFreq>1){
        for(k=1; k<nFilt-1; k++)
            FIRCoeffs(FIR_FILTER_BPF, order, fc[k-1], fc[k], sampleRate, windowType, scalingFLAG, &filterbank[k*(order+1)]);
    }
}

