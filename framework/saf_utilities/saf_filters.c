/*
 Copyright 2019 Leo McCormack
 
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
 *     saf_filters.c
 * Description:
 *     Contains a collection of filter design equations.
 * Dependencies:
 *     Windows users only: Intel's MKL must be installed, which can be freely aquired via:
 *     https://software.intel.com/en-us/articles/free-ipsxe-tools-and-libraries
 *     Mac users only: saf_utilities will utilise Apple's Accelerate library.
 * Author, date created:
 *     Leo McCormack, 01.03.2019
 */

#include "saf_filters.h" 

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
                b[2] = (V0 - sqrtf(2.0 * V0) * K + KK)/D;
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
         * H(z) = (b0 + b1*z^(-1) + b2*z^(-2)) / (1 + a1*z^(-1) + a2*z^(-2)), gives: */
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

    












 
