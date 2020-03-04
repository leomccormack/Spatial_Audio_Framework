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
 * @file: saf_filters.c
 * @brief Contains a collection of filter design equations.
 *
 * @author Leo McCormack
 * @date 01.03.2019
 */

#include "saf_filters.h" 
#include "saf_utilities.h"

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
                x[i] *= 0.54f - 0.46f * (cosf(2.0f*M_PI*(float)i/(float)N)); /* more wide-spread coefficient values */
            /* optimal equiripple coefficient values: */
            /*x[i] *= 0.53836f - 0.46164f * (cosf(2.0f*M_PI*(float)i/(float)N));*/
            break;
            
        case WINDOWING_FUNCTION_HANN:
            for(i=0; i<winlength; i++)
                x[i] *= 0.5f - 0.5f * (cosf(2.0f*M_PI*(float)i/(float)N));
            break;
            
        case WINDOWING_FUNCTION_BARTLETT:
            for(i=0; i<winlength; i++)
                x[i] *= 1.0f - 2.0f * fabsf((float)i-((float)N/2.0f))/(float)N;
            break;
            
        case WINDOWING_FUNCTION_BLACKMAN:
            for(i=0; i<winlength; i++){
                x[i] *= 0.42659f -
                        0.49656f *cosf(2.0f*M_PI*(float)i/(float)N) +
                        0.076849f*cosf(4.0f*M_PI*(float)i/(float)N);
            }
            break;
            
        case WINDOWING_FUNCTION_NUTTALL:
            for(i=0; i<winlength; i++){
                x[i] *= 0.355768f -
                        0.487396f*cosf(2.0f*M_PI*(float)i/(float)N) +
                        0.144232f*cosf(4.0f*M_PI*(float)i/(float)N) -
                        0.012604f*cosf(6.0f*M_PI*(float)i/(float)N);
            }
            break;
            
        case WINDOWING_FUNCTION_BLACKMAN_NUTTALL:
            for(i=0; i<winlength; i++){
                x[i] *= 0.3635819f -
                        0.4891775f*cosf(2.0f*M_PI*(float)i/(float)N) +
                        0.1365995f*cosf(4.0f*M_PI*(float)i/(float)N) +
                        0.0106411f*cosf(4.0f*M_PI*(float)i/(float)N);
            }
            break;
            
        case WINDOWING_FUNCTION_BLACKMAN_HARRIS:
            for(i=0; i<winlength; i++){
                x[i] *= 0.35875f -
                        0.48829f*cosf(2.0f*M_PI*(float)i/(float)N) +
                        0.14128f*cosf(4.0f*M_PI*(float)i/(float)N) +
                        0.01168f*cosf(4.0f*M_PI*(float)i/(float)N);
            }
            break;
    }
}

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
    float* centreFreqs,  /* nCutoffFreq x 1 */
    int nCentreFreqs,
    float* cutoffFreqs   /* (nCutoffFreq+1) x 1 */
)
{
    int i;
    for(i=0; i<nCentreFreqs-1; i++)
        cutoffFreqs[i] = 2.0f*centreFreqs[i]/sqrtf(2.0f);
}

// TODO:
//    mid of log
//    exp(1)^((log(400)+log(2000))/2)

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
    float K, KK, D, V0;
    
    a[0] = 1.0f;
    
    /* calculate the IIR filter coefficients */
    switch (filterType){
        case BIQUAD_FILTER_LPF:
            /* Filter design equations - DAFX (2nd ed) p50 */
            K = tanf(M_PI * fc/fs);
            KK = K * K;
            D = KK * Q + K + Q;
            b[0] = (KK * Q)/D;
            b[1] = (2.0f * KK * Q)/D;
            b[2] = b[0];
            a[1] = (2.0f * Q * (KK - 1.0f))/D;
            a[2] = (KK * Q - K + Q)/D;
            break;
            
        case BIQUAD_FILTER_HPF:
            /* Filter design equations - DAFX (2nd ed) p50 */
            K = tanf(M_PI * fc/fs);
            KK = K * K;
            D = KK * Q + K + Q;
            b[0] = (Q)/D;
            b[1] = -(2.0f * Q)/D;
            b[2] = b[0];
            a[1] = (2.0f * Q * (KK - 1.0f))/D;
            a[2] = (KK * Q - K + Q)/D;
            break;
            
        case BIQUAD_FILTER_LOW_SHELF:
            /* Filter design equations - DAFX (2nd ed) p64 */
            K = tanf(M_PI * fc/fs);
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
            
        case BIQUAD_FILTER_HI_SHELF:
            /* Filter design equations - DAFX (2nd ed) p64 */
            K = tanf(M_PI * fc/fs);
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
            
        case BIQUAD_FILTER_PEAK:
            /* Filter design equations - DAFX (2nd ed) p66 */
            K = tanf(M_PI * fc/fs);
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
    float* magnitude_dB,
    float* phase_rad
)
{
    int ff;
    float w, denom_real, denom_imag, num_real, num_imag;
    
    for(ff=0; ff<nFreqs; ff++){
        w = tanf(M_PI * freqs[ff]/fs);
        /* substituting Euler, z = e^(-jwn) = cos(wn) + j*sin(wn), into:
         * H(z) = (b0 + b1*z^(-1) + b2*z^(-2)) / (1 + a1*z^(-1) + a2*z^(-2)): */
        denom_real = 1.0f + a[1]*cosf(w) + a[2]*cosf(2.0f*w);
        denom_imag = a[1]*sinf(w) + a[2]*sinf(2*w);
        num_real = b[0] + b[1]*cosf(w) + b[2]*cosf(2.0f*w);
        num_imag = b[1]*sinf(w) + b[2]*sinf(2.0f*w);
        
        if(magnitude_dB!=NULL){
            magnitude_dB[ff] = sqrtf( (powf(num_real, 2.0f) + powf(num_imag, 2.0f)) / (powf(denom_real, 2.0f) + powf(denom_imag, 2.0f)) );
            magnitude_dB[ff] = 20.0f*log10f(magnitude_dB[ff]);
        }
        if(phase_rad!=NULL)
            phase_rad[ff] = atan2f(num_imag,num_real) - atan2f(denom_imag, denom_real);
    }
}

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
                    h_filt[i] = i==order/2 ? 2.0f*ft1 : sinf(2.0f*M_PI*ft1*(float)(i-order/2)) / (M_PI*(float)(i-order/2));
                break;
                
            case FIR_FILTER_HPF:
                for(i=0; i<h_len; i++)
                    h_filt[i] = i==order/2 ? 1.0f - 2.0f*ft1 : -sinf(2.0f*ft1*M_PI*(float)(i-order/2)) / (M_PI*(float)(i-order/2));
                break;
                
            case FIR_FILTER_BPF:
                ft2 = fc2/(fs*2.0f);
                for(i=0; i<h_len; i++){
                    h_filt[i] = i==order/2 ? 2.0f*(ft2-ft1) :
                        sinf(2.0f*M_PI*ft2*(float)(i-order/2)) / (M_PI*(float)(i-order/2)) - sinf(2.0f*M_PI*ft1*(float)(i-order/2)) / (M_PI*(float)(i-order/2));
                }
                break;
                
            case FIR_FILTER_BSF:
                ft2 = fc2/(fs*2.0f);
                for(i=0; i<h_len; i++){
                    h_filt[i] = i==order/2 ? 1.0f - 2.0f*(ft2-ft1) :
                        sinf(2.0f*M_PI*ft1*(float)(i-order/2)) / (M_PI*(float)(i-order/2)) - sinf(2.0f*M_PI*ft2*(float)(i-order/2)) / (M_PI*(float)(i-order/2));
                }
                break;
        }
    }
    else
        assert(0); /* please specify an even value for the filter 'order' argument */
    
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
                    h_z_sum = ccaddf(h_z_sum, crmulf(cexpf(cmplxf(0.0f, -2.0f*M_PI*(float)i*f0/2.0f)), h_filt[i]));
                h_sum = cabsf(h_z_sum);
                for(i=0; i<h_len; i++)
                    h_filt[i] /= h_sum;
                break;
                
            case FIR_FILTER_BPF:
                f0 = (fc1/fs+fc2/fs)/2.0f;
                h_z_sum = cmplxf(0.0f, 0.0f);
                for(i=0; i<h_len; i++)
                    h_z_sum = ccaddf(h_z_sum, crmulf(cexpf(cmplxf(0.0f, -2.0f*M_PI*(float)i*f0/2.0f)), h_filt[i]));
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
