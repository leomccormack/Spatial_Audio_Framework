/*
* Copyright 2020-2021 Michael McCrea, Leo McCormack
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
 * @file saf_utility_dvf.c
 * @ingroup Utilities
 * @brief Distance variation function filter coefficient data [1].
 *
 * @see [1] S. Spagnol, E. Tavazzi, and F. Avanzini, “Distance rendering and
 *          perception of nearby virtual sound sources with a near-field filter
 *          model,” Applied Acoustics, vol. 115, pp. 61–73, Jan. 2017,
 *          doi: 10.1016/j.apacoust.2016.08.015.
 *
 * @author Michael McCrea
 * @date 20.02.2021
 * @license ISC
 */

#include "saf_utility_dvf.h"
#include "saf_utility_filters.h"

/**
 * Table 1: Coefficients for Eqs. (8), (13), and (14) for generating high-shelf coefficients
 */
static const double p11[19] = { 12.97, 13.19, 12.13, 11.19, 9.91, 8.328, 6.493, 4.455, 2.274, 0.018, -2.24, -4.43, -6.49, -8.34, -9.93, -11.3, -12.2, -12.8, -13.0 };
static const double p21[19] = { -9.69, 234.2, -11.2, -9.03, -7.87, -7.42, -7.31, -7.28, -7.29, -7.48, -8.04, -9.23, -11.6, -17.4, -48.4, 9.149, 1.905, -0.75, -1.32 };
static const double q11[19] = { -1.14, 18.48, -1.25, -1.02, -0.83, -0.67, -0.5, -0.32, -0.11, -0.13, 0.395, 0.699, 1.084, 1.757, 4.764, -0.64, 0.109, 0.386, 0.45 };
static const double q21[19] = { 0.219, -8.5, 0.346, 0.336, 0.379, 0.421, 0.423, 0.382, 0.314, 0.24, 0.177, 0.132, 0.113, 0.142, 0.462, -0.14, -0.08, -0.06, -0.05 };
static const double p12[19] = { -4.39, -4.31, -4.18, -4.01, -3.87, -4.1, -3.87, -5.02, -6.72, -8.69, -11.2, -12.1, -11.1, -11.1, -9.72, -8.42, -7.44, -6.78, -6.58 };
static const double p22[19] = { 2.123, -2.78, 4.224, 3.039, -0.57, -34.7, 3.271, 0.023, -8.96, -58.4, 11.47, 8.716, 21.8, 1.91, -0.04, -0.66, 0.395, 2.662, 3.387 };
static const double q12[19] = { -0.55, 0.59, -1.01, -0.56, 0.665, 11.39, -1.57, -0.87, 0.37, 5.446, -1.13, -0.63, -2.01, 0.15, 0.243, 0.147, -0.18, -0.67, -0.84 };
static const double q22[19] = { -0.06, -0.17, -0.02, -0.32, -1.13, -8.3, 0.637, 0.325, -0.08, -1.19, 0.103, -0.12, 0.098, -0.4, -0.41, -0.34, -0.18, 0.05, 0.131 };
static const double p13[19] = { 0.457, 0.455, -0.87, 0.465, 0.494, 0.549, 0.663, 0.691, 3.507, -27.4, 6.371, 7.032, 7.092, 7.463, 7.453, 8.101, 8.702, 8.925, 9.317 };
static const double p23[19] = { -0.67, 0.142, 3404., -0.91, -0.67, -1.21, -1.76, 4.655, 55.09, 10336., 1.735, 40.88, 23.86, 102.8, -6.14, -18.1, -9.05, -9.03, -6.89 };
static const double p33[19] = { 0.174, -0.11, -1699., 0.437, 0.658, 2.02, 6.815, 0.614, 589.3, 16818., -9.39, -44.1, -23.6, -92.3, -1.81, 10.54, 0.532, 0.285, -2.08 };
static const double q13[19] = { -1.75, -0.01, 7354., -2.18, -1.2, -1.59, -1.23, -0.89, 29.23, 1945., -0.06, 5.635, 3.308, 13.88, -0.88, -2.23, -0.96, -0.9, -0.57 };
static const double q23[19] = { 0.699, -0.35, -5350., 1.188, 0.256, 0.816, 1.166, 0.76, 59.51, 1707., -1.12, -6.18, -3.39, -12.7, -0.19, 1.295, -0.02, -0.08, -0.4 };
static const int numAz_table = sizeof(q23);

/* a_0 = 0.0875f; Reference head size, 8.75 centimeters, used in the generation of the coeff lookup table. */
/* a_head = 0.09096f; This head size, see note for head_radius in binauraliser_nf. */
static const float headDim = SAF_PI * (0.0875f / 0.09096f);             /* pi * (a_0 / a_head) */
static const float sosDiv2PiA = 343.0f / (2.0f * SAF_PI * 0.09096f);    /*  c / (2pi * a_head) */

/** Linear interpolation between two values */
static float interpolate_lin
(
    float a,
    float b,
    float ifac
)
{
    return a + (b-a) * ifac;
}

/** Covert decibels to a magnitude */
static float db2mag
(
    float dB
)
{
    return powf(10.f, dB / 20.f);
}

/*
 * Calculate high-shelf parameters, g0, gInf, fc, from the lookup table coefficients (10 degree steps).
 * Called twice per update as the returned values are subsequently interpolated to exact azimuth. */
void calcDVFShelfParams
(
    int i,          /* Index into the coefficient table, dictated by azimuth. */
    float rhoIn,    /* Normalized source distance. */
    float* g0,      /* High shelf gain at DC. */
    float* gInf,    /* High shelf gain at inf. */
    float* fc       /* High shelf cutoff frequency. */
)
{
    float fc_tmp;
    double rho = (double)rhoIn;
	double rhoSq = rho * rho;
    
    /*  Eq (8), (13) and (14) */
    *g0    = (float)((p11[i] * rho   + p21[i]) / (rhoSq + q11[i] * rho + q21[i]));
    *gInf  = (float)((p12[i] * rho   + p22[i]) / (rhoSq + q12[i] * rho + q22[i]));
    fc_tmp = (float)((p13[i] * rhoSq + p23[i] * rho + p33[i]) / (rhoSq + q13[i] * rho + q23[i]));
    
    /* denormalize fc = fc * sos/(2pi*a) */
    *fc = fc_tmp * sosDiv2PiA;
}

/*
 * Interpolate (linearly) the high shelf parameters generated by
 * calcDVFShelfParams() which is called twice to generate the high shelf
 * parameters for the nearest thetas in the lookup table. */
void interpDVFShelfParams
(
    float theta,   /* Ipsilateral azimuth, on the inter-aural axis [0, 180] (deg). */
    float rho,     /* Distance, normalized to head radius, >= 1. */
    /* output */
    float* iG0,
    float* iGInf,
    float* iFc
)
{
    int theta_idx_lower, theta_idx_upper;
    float ifac, thetaDiv10;
    float g0_1, g0_2;       /* high shelf gain at DC */
    float gInf_1, gInf_2;   /* high shelf gain at inf */
    float fc_1, fc_2;       /* high shelf cutoff frequency */
    
    // TBD: could add range checking, clipping theta and rho to valid range
    
    /* Linearly interpolate DC gain, HF gain, center freq at theta.
     * Table is in 10 degree steps, floor(x/10) gets lower index. */
    thetaDiv10 = theta / 10.f;
    theta_idx_lower = (int)thetaDiv10;
    theta_idx_upper = theta_idx_lower + 1;
    if(theta_idx_upper == numAz_table) {    // could alternatively check theta_idx_upper => numAz_table, clip the value > 180 here
        theta_idx_upper = theta_idx_lower;
        theta_idx_lower = theta_idx_lower - 1;
    }
    
    calcDVFShelfParams(theta_idx_lower, rho, &g0_1, &gInf_1, &fc_1);
    calcDVFShelfParams(theta_idx_upper, rho, &g0_2, &gInf_2, &fc_2);

    /* interpolation factor between table steps */
    ifac   = thetaDiv10 - theta_idx_lower;
    *iG0   = interpolate_lin(g0_1,   g0_2,   ifac);
    *iGInf = interpolate_lin(gInf_1, gInf_2, ifac);
    *iFc   = interpolate_lin(fc_1,   fc_2,   ifac);
}

/*
 * Calculate the DVF filter coefficients from shelving filter parameters.
 */
void dvfShelfCoeffs
(
    /* Input */
    float g0,      /* High shelf dc gain */
    float gInf,    /* High shelf high gain */
    float fc,      /* High shelf center freq */
    float fs,      /* Sample rate */
    /* Output */
    float* b0,     /* IIR coeffs */
    float* b1,
    float* a1
)
{
    float v0, g0_mag, tanF, v0tanF;
    float a_c, v, va_c;
    
    v0     = db2mag(gInf);              /* Eq. (12), (10), and (11) */
    g0_mag = db2mag(g0);                // optim TBD: revisit; does g0, gInf need to be in dB?
    tanF   = tanf((headDim / fs) * fc); // optim TBD: this /fs calc can be optimized out with precomputed head dimension
    v0tanF = v0 * tanF;
    a_c    = (v0tanF - 1.f) / (v0tanF + 1.f);
    
    v    = (v0 - 1.f) * 0.5f;           /* Eq (10) */
    va_c = v * a_c;
    *b0  = g0_mag * (v - va_c + 1.f);   /* = V*(1 - a_c) + 1   */
    *b1  = g0_mag * (va_c - v + a_c);   /* = V*(a_c - 1) + a_c */
    *a1  = a_c;
}

void calcDVFCoeffs
(
    float alpha,    /* Lateral angle, on the inter-aural axis [0, 180] (deg) */
    float rho,      /* Distance, normalized to head radius, >= 1 */
    float fs,       /* Sample rate */
    float* b,
    float* a
)
{
    float iG0, iGInf, iFc;
    interpDVFShelfParams(alpha, rho, &iG0, &iGInf, &iFc);
    dvfShelfCoeffs(iG0, iGInf, iFc, fs, &b[0], &b[1], &a[1]);
}

void doaToIpsiInteraural
(
    float azimuth,      /* azimuth wrt 0˚ forward-facing  (deg, [-360, 360]) */
    float elevation,    /* elevation from the horizontal plane (deg, [-180, 180]) */
    float* alphaLR,     /* (&) 2-element array of lateral angle alpha for left and right ear (deg, [0,180]) */
    float* betaLR       /* (&) 2-element array of vertical angle beta for left and right ear (deg, [0,90]) */
)
{
    float az_rad, el_rad, alpha_ipsi, beta_ipsi, alpha_ipsi_deg, beta_ipsi_deg;
    float sinaz, sinel, cosaz, cosel;
    
    az_rad = DEG2RAD(azimuth);
    el_rad = DEG2RAD(elevation);
    sinaz = sinf(az_rad);
    sinel = sinf(el_rad);
    cosaz = cosf(az_rad);
    cosel = cosf(el_rad);
    alpha_ipsi = SAF_PI/2.f - acosf(sinaz * cosel);
    beta_ipsi  = asinf(sinel / sqrtf(powf(sinel,2.f) + (powf(cosaz,2.f) * powf(cosel,2.f))));
    
    /* Convert interaural-polar coordinates back to spherical _ranges_. */
    if(beta_ipsi > SAF_PI/2.f) {
        alpha_ipsi = SAF_PI - alpha_ipsi;   /* [0,pi] */
        beta_ipsi = SAF_PI - beta_ipsi;     /* [0,pi/2] */
    }

    /* Convert to ipsilateral offset from the left ear: [-360, 360]->[0, 180] */
    alpha_ipsi = fabsf(SAF_PI/2.f - alpha_ipsi);
    if (alpha_ipsi > SAF_PI)
        alpha_ipsi = 2*SAF_PI - alpha_ipsi;
    
    alphaLR[0] = alpha_ipsi_deg = RAD2DEG(alpha_ipsi);
    alphaLR[1] = 180.f - alpha_ipsi_deg;
    
    if (betaLR != NULL) {
        betaLR[0] = beta_ipsi_deg = RAD2DEG(beta_ipsi);
        betaLR[1] = 180.f - beta_ipsi_deg;
    }
}
