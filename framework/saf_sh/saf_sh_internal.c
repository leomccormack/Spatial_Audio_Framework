/*
 Copyright 2016-2018 Leo McCormack
 
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
 *     saf_sh_internal.c
 * Description:
 *     A collection of spherical harmonic related functions. Some of which have been
 *     derived from the Matlab library by Archontis Politis; found here:
 *     https://github.com/polarch/Spherical-Harmonic-Transform
 *     and MATLAB code by Symeon Delikaris-Manias
 * Dependencies:
 *     saf_utilities
 * Author, date created:
 *     Leo McCormack, 22.05.2016
 */

#include "saf_sh.h"
#include "saf_sh_internal.h"

static unsigned long factorial(unsigned long f)
{
    if (f == 0)
        return 1;
    else
        return(f * factorial(f - 1));
}

/* A C implementation of a MatLab function by Symeon Delikaris-Manias; published with permission */
void ChebyshevPolyCoeff
(
    int     n,
    float * t_coeff
)
{
    int k, e, i;
    float* t_coeffm2, *t_coeffm1;

    if (n == 0)
        t_coeff[0] = 1.0f;
    else if (n == 1) {
        t_coeff[0] = 1.0f;
        t_coeff[1] = 0.0f;
    }
    else {
        t_coeffm2 = (float*)calloc((n + 1), sizeof(float));
        t_coeffm2[n] = 1.0f;
        t_coeffm1 = (float*)calloc((n + 1), sizeof(float));
        t_coeffm1[n - 1] = 1.0f;
        for (k = 2; k <= n; k++) {
            memset(t_coeff, 0, (n + 1) * sizeof(float));
            for (e = (n - k + 1); e <= n; e = e + 2)
                t_coeff[e - 1] = 2.0f*t_coeffm1[e] - t_coeffm2[e - 1];
            if ((k % 2) == 0)
                t_coeff[n] = powf(-1.0f, (float)k / 2.0f);
            if (k<n) {
                for (i = 0; i < n + 1; i++) { /* shuffle */
                    t_coeffm2[i] = t_coeffm1[i];
                    t_coeffm1[i] = t_coeff[i];
                }
            }
        }
        free(t_coeffm2);
        free(t_coeffm1);
    }
}

/* A C implementation of a MatLab function by Symeon Delikaris-Manias; published with permission */
void LegendrePolyCoeff
(
    int     n,
    float * p_coeff
)
{
    int k, e, i;
    float* p_coeffm2, *p_coeffm1;

    if (n == 0)
        p_coeff[0] = 1.0f;
    else if (n == 1) {
        p_coeff[0] = 1.0f;
        p_coeff[1] = 0.0f;
    }
    else {
        p_coeffm2 = (float*)calloc((n + 1), sizeof(float));
        p_coeffm2[n] = 1.0f;
        p_coeffm1 = (float*)calloc((n + 1), sizeof(float));
        p_coeffm1[n - 1] = 1.0f;
        for (k = 2; k <= n; k++) {
            memset(p_coeff, 0, (n + 1) * sizeof(float));
            for (e = (n - k + 1); e <= n; e = e + 2)
                p_coeff[e - 1] = (2.0f*(float)k - 1.0f)*p_coeffm1[e] + (1.0f - (float)k)*p_coeffm2[e - 1];
            p_coeff[n] += (1.0f - (float)k)*p_coeffm2[n];
            for (i = 0; i <= n; i++)
                p_coeff[i] = p_coeff[i] / (float)k;
            if (k<n) {
                for (i = 0; i < n + 1; i++) { /* shuffle */
                    p_coeffm2[i] = p_coeffm1[i];
                    p_coeffm1[i] = p_coeff[i];
                }
            }
        }
        free(p_coeffm2);
        free(p_coeffm1);
    }
}

/* A C implementation of a MatLab function by Symeon Delikaris-Manias; published with permission */
void dolph_chebyshev
(
    int                M,
    float            * d,
    int   type
)
{
    int n, i, j, k, q, s, m;
    float SNR, R, x0, theta0;
    float *x00, **P, **A, **C, **T, **PA, **CT, **PACT;
    float* p_coeff, *t_coeff;

    x00 = (float*)malloc((M + 1) * sizeof(float));
    P = (float**)malloc2d((M + 1), (M + 1), sizeof(float));
    A = (float**)calloc2d((M + 1), (M + 1), sizeof(float));
    C = (float**)calloc2d((M + 1), (M + 1), sizeof(float));
    T = (float**)calloc2d((M + 1), (M + 1), sizeof(float));
    PA = (float**)calloc2d((M + 1), (M + 1), sizeof(float));
    CT = (float**)calloc2d((M + 1), (M + 1), sizeof(float));
    PACT = (float**)calloc2d((M + 1), (M + 1), sizeof(float));
    p_coeff = (float*)calloc((M + 1), sizeof(float));
    t_coeff = (float*)calloc((2 * M + 1), sizeof(float));

    switch (type) {
        case 1: /* DESIRED_LOBE */
            SNR = 25.0f;
            R = powf(10.0f, SNR / 20.0f);
            x0 = coshf(1.0f / (2.0f*(float)M)*acoshf(R));
            //theta0 = 2.0f*acosf((1.f / x0)*cosf((float)PI / (4.0f*(float)M))) * (180.0f / (float)PI);
            break;
            
        case 0: /* MAIN_LOBE */
            theta0 = 60.0f* (float)PI / 180.0f;
            x0 = cosf((float)PI / (4.0f*(float)M)) / (cosf(theta0 / 2.0f));
            R = coshf(2.0f*(float)M*acoshf(x0));
            break;
            
        default:
            R = 0.0f;
            x0 = 0.0f;
            break;
    }
    for (n = 0; n <= M; n++)
        x00[n] = powf(x0, 2.0f*(float)n);
    for (n = 0; n <= M; n++) {
        LegendrePolyCoeff(n, p_coeff);
        memset(P[n], 0, (M + 1) * sizeof(float));
        for (i = 0; i <= n; i++)
            P[n][i] = p_coeff[n - i];
    }
    for (q = 0; q <= M; q++)
        for (s = 0; s <= M; s++)
            A[q][s] = ((1.0f - powf(-1.0f, (float)q + (float)s + 1.0f)) / ((float)q + (float)s + 1.0f));
    for (n = 0; n <= M; n++)
        for (m = 0; m <= n; m++)
            C[m][n] = (powf(2.0f, -(float)n) * (float)factorial(n)) / (float)(factorial(m) * (factorial(n - m)));
    ChebyshevPolyCoeff(M * 2, t_coeff);
    for (n = 0; n <= M; n++)
        T[n][n] = t_coeff[2 * (M - n)];
    for (i = 0; i <= M; i++) {
        for (j = 0; j <= M; j++) {
            for (k = 0; k <= M; k++) {
                PA[i][j] += P[i][k] * A[k][j];
                CT[i][j] += C[i][k] * T[k][j];
            }
        }
    }
    for (i = 0; i <= M; i++)
        for (j = 0; j <= M; j++)
            for (k = 0; k <= M; k++)
                PACT[i][j] += PA[i][k] * CT[k][j];
    memset(d, 0, (M + 1) * sizeof(float));
    for (i = 0; i <= M; i++)
        for (j = 0; j <= M; j++)
            d[i] += PACT[i][j] * x00[j] * (2.0f*(float)PI / R);

    free(x00);
    free2d((void**)P, (M + 1));
    free2d((void**)A, (M + 1));
    free2d((void**)C, (M + 1));
    free2d((void**)T, (M + 1));
    free2d((void**)PA, (M + 1));
    free2d((void**)CT, (M + 1));
    free2d((void**)PACT, (M + 1));
    free(p_coeff);
    free(t_coeff);
}

/* A C implementation of a MatLab function by Symeon Delikaris-Manias; published with permission */
void maxre3d
(
    int                M,
    float            * gm
)
{
    int i, nz, k, i0;
    float zmin, zmax, dz, zc, dz2, rE;
    float* z, **p;

    nz = 3 * M + 10;
    zmin = 0.5f;
    z = (float*)malloc(nz * sizeof(float));
    for (i = 0; i<nz; i++) 
        z[i] = (float)(i + 1) / (float)nz*(1.0f - zmin) + zmin;
    p = (float**)malloc(nz * sizeof(float*));
    for (i = 0; i<nz; i++) {
        p[i] = (float*)malloc((M + 2) * sizeof(float));
        legendreP(M + 1, z[i], p[i]);
    }
    k = 0;
    dz = 1.0f;

    /* Finding the smallest root of the Legendre polynomial gives rE, which is
        one of the first terms of the recurrence applied later */
    while ((dz>1e-7f) && (k<7)) {
        i0 = 0;
        for (i = 1; i<nz; i++)
            if ((p[i - 1][0] <= 0.0f) && (p[i][0] > 0.0f))
                i0 = MAX(i0, i);
        dz = z[i0 + 1] - z[i0];
        zc = (z[i0] * p[i0 + 1][0] - z[i0 + 1] * p[i0][0]) / (p[i0 + 1][0] - p[i0][0]);

        dz2 = MAX(zc - z[i0], z[i0 + 1] - zc);
        k++;
        free(z);
        for (i = 0; i<nz; i++)
			free(p[i]);
        free(p);
        nz = 14;
        zmin = zc - dz2;
        zmax = zc + dz2;
        z = (float*)malloc(nz * sizeof(float));
        for (i = 0; i<nz; i++)
            z[i] = (float)(i + 1) / (float)nz*(zmax - zmin) + zmin;
        p = (float**)malloc(nz * sizeof(float*));
        for (i = 0; i<nz; i++) {
            p[i] = (float*)malloc((M + 2) * sizeof(float));
            legendreP(M + 1, z[i], p[i]);
        }
    }
    rE = 0.0f;
    for (i = 0; i<nz; i++)
        rE += z[i];
    rE /= (float)nz;

    gm[0] = 1.0f;
    gm[1] = rE;
    for (i = 1; i<M; i++)
        gm[i + 1] = ((2.0f*(float)i + 1.0f)*rE*gm[i] - (float)i*gm[i - 1]) / ((float)i + 1.0f);
	
    free(z);
    for (i = 0; i<nz; i++)
        free(p[i]);
    free(p);
}







