/*
* Copyright ...
*
* Permission to use ...
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
* @file saf_utility_dvf.c
* @ingroup Utilities
* @brief Distance variation function filter coefficient data
*
*      Simone Spagnol, Erica Tavazzi, and Federico Avanzini. Distance rendering and perception of nearby
*      virtual sound sources with a near-field filter model. Applied Acoustics, 115:61–73, January 2017.
*
* @author Michael McCrea
* @date 20.02.2021
*/

#include "saf_utility_dvf.h"
#include "saf_utility_filters.h"

/**
 * Table 1: Coefficients for Eqs. (8), (13), and (14) for generating high-shelf coefficients
 */
const float p11[] = { 12.97f, 13.19f, 12.13f, 11.19f, 9.91f, 8.328f, 6.493f, 4.455f, 2.274f, 0.018f, -2.24f, -4.43f, -6.49f, -8.34f, -9.93f, -11.3f, -12.2f, -12.8f, -13.f };
const float p21[] = { -9.69f, 234.2f, -11.2f, -9.03f, -7.87f, -7.42f, -7.31f, -7.28f, -7.29f, -7.48f, -8.04f, -9.23f, -11.6f, -17.4f, -48.4f, 9.149f, 1.905f, -0.75f, -1.32f };
const float q11[] = { -1.14f, 18.48f, -1.25f, -1.02f, -0.83f, -0.67f, -0.5f, -0.32f, -0.11f, -0.13f, 0.395f, 0.699f, 1.084f, 1.757f, 4.764f, -0.64f, 0.109f, 0.386f, 0.45f };
const float q21[] = { 0.219f, -8.5f, 0.346f, 0.336f, 0.379f, 0.421f, 0.423f, 0.382f, 0.314f, 0.24f, 0.177f, 0.132f, 0.113f, 0.142f, 0.462f, -0.14f, -0.08f, -0.06f, -0.05f };
const float p12[] = { -4.39f, -4.31f, -4.18f, -4.01f, -3.87f, -4.1f, -3.87f, -5.02f, -6.72f, -8.69f, -11.2f, -12.1f, -11.1f, -11.1f, -9.72f, -8.42f, -7.44f, -6.78f, -6.58f };
const float p22[] = { 2.123f, -2.78f, 4.224f, 3.039f, -0.57f, -34.7f, 3.271f, 0.023f, -8.96f, -58.4f, 11.47f, 8.716f, 21.8f, 1.91f, -0.04f, -0.66f, 0.395f, 2.662f, 3.387f };
const float q12[] = { -0.55f, 0.59f, -1.01f, -0.56f, 0.665f, 11.39f, -1.57f, -0.87f, 0.37f, 5.446f, -1.13f, -0.63f, -2.01f, 0.15f, 0.243f, 0.147f, -0.18f, -0.67f, -0.84f };
const float q22[] = { -0.06f, -0.17f, -0.02f, -0.32f, -1.13f, -8.3f, 0.637f, 0.325f, -0.08f, -1.19f, 0.103f, -0.12f, 0.098f, -0.4f, -0.41f, -0.34f, -0.18f, 0.05f, 0.131f };
const float p13[] = { 0.457f, 0.455f, -0.87f, 0.465f, 0.494f, 0.549f, 0.663f, 0.691f, 3.507f, -27.4f, 6.371f, 7.032f, 7.092f, 7.463f, 7.453f, 8.101f, 8.702f, 8.925f, 9.317f };
const float p23[] = { -0.67f, 0.142f, 3404.f, -0.91f, -0.67f, -1.21f, -1.76f, 4.655f, 55.09f, 10336.f, 1.735f, 40.88f, 23.86f, 102.8f, -6.14f, -18.1f, -9.05f, -9.03f, -6.89f };
const float p33[] = { 0.174f, -0.11f, -1699.f, 0.437f, 0.658f, 2.02f, 6.815f, 0.614f, 589.3f, 16818.f, -9.39f, -44.1f, -23.6f, -92.3f, -1.81f, 10.54f, 0.532f, 0.285f, -2.08f };
const float q13[] = { -1.75f, -0.01f, 7354.f, -2.18f, -1.2f, -1.59f, -1.23f, -0.89f, 29.23f, 1945.f, -0.06f, 5.635f, 3.308f, 13.88f, -0.88f, -2.23f, -0.96f, -0.9f, -0.57f };
const float q23[] = { 0.699f, -0.35f, -5350.f, 1.188f, 0.256f, 0.816f, 1.166f, 0.76f, 59.51f, 1707.f, -1.12f, -6.18f, -3.39f, -12.7f, -0.19f, 1.295f, -0.02f, -0.08f, -0.4f };

const int numAz_table = sizeof(q23);

/* 8.75 centimeters, reference head size used to generate coeff lookup table */
const float a_0 = 0.0875;
/* this head size, TODO: make this a parameter */
const float a_head = 0.0875;
const float headDim = M_PI * (a_0 / a_head); // TODO: use saf_PI
const float sosDiv2PiA = 343 / (M_PI_2 * a_head);   // TODO: use saf_PI


/**
* Calculate high-shelf parameters from the lookup table coefficients (10 degree steps).
* Called twice per update as the returned values are subsequently interpolated to exact azimuth. */
static void calcHighShelfParams
(
    int i,          /* index into the coefficient table, dictated by azimuth */
    float rho,      /* normalized source distance */
    float* g0,      /* high shelf gain at DC */
    float* gInf,    /* high shelf gain at inf */
    float* fc       /* high shelf cutoff frequency */
)
{
    float rhoSq, fc_tmp;
    
    rhoSq = powf(rho, 2.0f);
    
    /*  Eq (8), (13) and (14) */
    *g0    = (p11[i] * rho   + p21[i]) / (rhoSq + q11[i] * rho + q21[i]);
    *gInf  = (p12[i] * rho   + p22[i]) / (rhoSq + q12[i] * rho + q22[i]);
    fc_tmp = (p13[i] * rhoSq + p23[i] * rho + p33[i]) / (rhoSq + q13[i] * rho + q23[i]);
    
    /* denormalize fc = fc * sos/(2pi*a) */
    *fc = fc_tmp * sosDiv2PiA;
}

static void calcIIRCoeffs
(
 float g0,      /* high shelf dc gain */
 float gInf,    /* high shelf high gain */
 float fc,      /* high shelf center freq */
 float fs,      /* sample rate */
 float* b0,     /* IIR coeffs */
 float* b1,
 float* a1
)
{
    float v0;
    float g0_mag;
    float tanF;
    float v0tanF;
    float a_c;
    float v;
    float va_c;
    
    v0     = db2mag(gInf);             /* Eq. (12), (10), and (11) */
    g0_mag = db2mag(g0);
    tanF   = tanf((headDim / fs) * fc);   // TODO: this /fs calc can be optimized out with precomputed head dimension
    v0tanF = v0 * tanF;
    a_c    = (v0tanF - 1.f) / (v0tanF + 1.f);
    
    v    = (v0 - 1) * 0.5;             /* Eq (10) */
    va_c = v * a_c;
    *b0  = g0_mag * (v - va_c + 1);    /* = V*(1 - a_c) + 1   */
    *b1  = g0_mag * (va_c - v + a_c);  /* = V*(a_c - 1) + a_c */
    *a1  = a_c;
}


/*
 *
 */
static void calcHighShelfCoeffs
(
 float theta,   /* ipsilateral azimuth, on the inter-aural axis [0, 180] (deg) */
 float rho,     /* distance, normalized to head radius, >= 1 */
 float fs,      /* sample rate */
 float* b0,     /* IIR coeffs */
 float* b1,
 float* a1
)
{
    int theta_idx_lower, theta_idx_upper;
    float ifac;
    float thetaDiv10;
    float g0;
    float gInf;
    float fc;
    // TODO: check pointer instantiation logic
    float* g0_1   = NULL; /* high shelf gain at DC */
    float* g0_2   = NULL;
    float* gInf_1 = NULL; /* high shelf gain at inf */
    float* gInf_2 = NULL;
    float* fc_1   = NULL; /* high shelf cutoff frequency */
    float* fc_2   = NULL;
    
    // TODO: range checking - clip theta and rho to valid range
    /* linearly interpolate DC gain, HF gain, center freq at theta */
    // TODO: rethink this indexing logic...
    thetaDiv10 = theta / 10.f;
    theta_idx_lower = (int)thetaDiv10;      /* because table is in 10 degree steps, floor(x/10) gets lower index */
    theta_idx_upper = theta_idx_lower + 1;
    if(theta_idx_upper == numAz_table) {    // TODO: if instead check theta_idx_upper => numAz_table, could clip the value > 180 here
        theta_idx_upper = theta_idx_lower;
        theta_idx_lower = theta_idx_lower - 1;
    }
    
    calcHighShelfParams(theta_idx_lower, rho, g0_1, gInf_1, fc_1);
    calcHighShelfParams(theta_idx_upper, rho, g0_2, gInf_2, fc_2);

    ifac = 1.f - thetaDiv10;                /* interpolation factor between table steps */
    g0   = interpolate_lin(*g0_1,   *g0_2,   ifac);
    gInf = interpolate_lin(*gInf_1, *gInf_2, ifac);
    fc   = interpolate_lin(*fc_1,   *fc_2,   ifac);
    
    calcIIRCoeffs(g0, gInf, fc, fs, b0, b1, a1);
}


static void applyDVF
(
    float theta,   /* ipsilateral azimuth, on the inter-aural axis [0, 180] (deg) */
    float rho,     /* distance, normalized to head radius, >= 1 */
    float* in_signal,
    int nSamples,
    float fs,
    float* wz,
    float* out_signal
)
{
    float b[2] = {0.f, 0.f};
    float a[2] = {1.f, 0.f};
    
    calcHighShelfCoeffs(theta, rho, fs, &b[0], &b[1], &a[1]);       // TODO: wacky pointer syntax
    applyIIR(in_signal, nSamples, 2, &b[0], &a[0], wz, out_signal); // TODO: wacky pointer syntax
}
