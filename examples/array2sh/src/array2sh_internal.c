/*
 Copyright 2017-2018 Leo McCormack
 
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
 *     array2sh_internal.c
 * Description:
 *     Spatially encodes spherical or cylindrical sensor array signals into spherical harmonic
 *     signals utilising theoretical encoding filters.
 *     The algorithms within array2sh were pieced together and developed in collaboration
 *     with Symeon Delikaris-Manias.
 *     A more detailed explanation of the algorithms in array2sh can be found in:
 *     McCormack, L., Delikaris-Manias, S., Farina, A., Pinardi, D., and Pulkki, V.,
 *     “Real-time conversion of sensor array signals into spherical harmonic signals with
 *     applications to spatially localised sub-band sound-field analysis,” in Audio
 *     Engineering Society Convention 144, Audio Engineering Society, 2018.
 * Dependencies:
 *     saf_utilities, afSTFTlib, saf_sh
 * Author, date created:
 *     Leo McCormack, 13.09.2017
 */

#include "array2sh_internal.h"

static long factorial(int n)
{
    if (n == 0)
        return 1;
    else
        return(n * factorial(n-1));
}

static double_complex bessel_Hl2(int l, double z){
    double_complex Hl2;
#ifdef __APPLE__
    Hl2 = cmplx(jn(l,z), - yn(l,z));
#else
	Hl2 = cmplx(_jn(l, z), -_yn(l, z));
#endif
    return Hl2;
}

/* Spherical Bessel of the first kind j_l, l=-1..8 */
static double bessel_jl(int l, double z)
{
    double jl;
    z = MAX(z, 2.23e-13);

    /* if z << l, then an approximation of jl is more suitable, for higher values this approach breaks down */
    if (z < l/(z+2.23e-13)){
        if (l == -1){
            jl = sqrt(2.0/(M_PI*z))*cos(z);
            jl = jl*sqrt(M_PI/(2.0*z));
        }
        else {
            jl = pow(z,l) / ((double)factorial(2*l+1)/ (pow(2.0,l) * (double)factorial(l)));
        }
    }
    else {
        /* for higher 'z' values, the exact analytical solution is a better idea.
           (the insufficient numerical precision breaks this approach at lower 'z' values) */
        switch (l){
            case -1:
                /* reverts to negative half order cylindrical bessel J_l */
                jl = sqrt(2.0/(M_PI*z))*cos(z);
                jl = jl*sqrt(M_PI/(2.0*z));
                break;
            case 0:
                jl = sin(z)/z;
                break;
            case 1:
                jl = (sin(z) - z*cos(z))/pow(z,2.0);
                break;
            case 2:
                jl = -(pow(z,2.0)*sin(z) - 3.0*sin(z) + 3.0*z*cos(z))/pow(z,3.0);
                break;
            case 3:
                jl = (15.0*sin(z) + pow(z,3.0)*cos(z) - 6.0*pow(z,2.0)*sin(z) - 15.0*z*cos(z))/pow(z,4.0);
                break;
            case 4:
                jl = (105.0*sin(z) + 10.0*pow(z,3.0)*cos(z) - 45.0*pow(z,2.0)*sin(z) + pow(z,4.0)*sin(z) -
                      105.0*z*cos(z))/pow(z,5.0);
                break;
            case 5:
                jl = (945.0*sin(z) + 105.0*pow(z,3.0)*cos(z) - pow(z,5.0)*cos(z) - 420.0*pow(z,2.0)*sin(z) +
                      15.0*pow(z,4.0)*sin(z) - 945.0*z*cos(z))/pow(z,6.0);
                break;
            case 6:
                jl = -(21.0*pow(z,5.0)*cos(z) - 1260.0*pow(z,3.0)*cos(z) - 10395.0*sin(z) + 4725.0*pow(z,2.0)*sin(z)
                       - 210.0*pow(z,4.0)*sin(z) + pow(z,6.0)*sin(z) + 10395.0*z*cos(z))/pow(z,7.0);
                break;
            case 7:
                jl = (135135.0*sin(z) + 17325.0*pow(z,3.0)*cos(z) - 378.0*pow(z,5.0)*cos(z) + pow(z,7.0)*cos(z) -
                      62370.0*pow(z,2.0)*sin(z) + 3150.0*pow(z,4.0)*sin(z) - 28.0*pow(z,6.0)*sin(z) -
                      135135.0*z*cos(z))/pow(z,8.0);
                break;
            case 8:
                jl = (2027025.0*sin(z) + 270270.0*pow(z,3.0)*cos(z) - 6930.0*pow(z,5.0)*cos(z) + 36.0*pow(z,7.0)*cos(z) -
                      945945.0*pow(z,2.0)*sin(z) + 51975.0*pow(z,4.0)*sin(z) - 630.0*pow(z,6.0)*sin(z) +
                      pow(z,8.0)*sin(z) - 2027025.0*z*cos(z))/pow(z,9.0);
                break;
            default:
                jl = 2.23e-13;
                break;
        }
    }
    return jl;
}

/* Spherical Bessel of the third kind (Hankel of the second kind) h_l^(2), l=-1..8 */
static double_complex bessel_hl2(int l, double z)
{
    double_complex hl_2;
    z = MAX(z, 2.3e-13);
    switch (l){
        case -1:
            hl_2 = cmplx(sqrt(2.0/(M_PI*z))*cos(z), -sqrt(2.0/(M_PI*z))*sin(z));
            hl_2 = crmul(hl_2, sqrt(M_PI/(2.0*z)));
            break;
        case 0:
            hl_2 = ccmul(cexp(cmplx(0.0, -z)), cmplx(0.0, 1.0));
            hl_2 = ccdiv(hl_2, cmplx(z, 0.0));
            break;
        case 1:
            hl_2 = ccmul(cexp(cmplx(0.0, -z)), (cmplx(z, -1.0)));
            hl_2 = crmul(hl_2, -1.0);
            hl_2 = ccdiv(hl_2, cmplx(pow(z, 2.0), 0.0));
            break;
        case 2:
            hl_2 = ccmul(cexp(cmplx(0.0,-z)), (cmplx(3.0*z, pow(z,2.0) - 3.0)));
            hl_2 = crmul(hl_2, -1.0);
            hl_2 = ccdiv(hl_2, cmplx(pow(z, 3.0), 0.0));
            break;
        case 3:
            hl_2 = ccmul(cexp(cmplx(0.0,-z)), (cmplx( -pow(z,3.0) + 15.0*z, pow(z,2.0)*6.0 - 15.0)));
            hl_2 = crmul(hl_2, -1.0);
            hl_2 = ccdiv(hl_2, cmplx(pow(z, 4.0), 0.0));
            break;
        case 4:
            hl_2 = ccmul(cexp(cmplx(0.0,-z)), (cmplx(10.0*pow(z,3.0) - 105.0*z, pow(z,4.0)  - pow(z,2.0)*45.0 + 105.0)));
            hl_2 = ccdiv(hl_2, cmplx(pow(z, 5.0), 0.0));
            break;
        case 5:
            hl_2 = ccmul(cexp(cmplx(0.0,-z)), (cmplx(pow(z,5.0) - 105.0*pow(z,3.0) + 945.0*z, -pow(z,4.0)*15.0 +
                                                     pow(z,2.0)*420.0 - 945.0)));
            hl_2 = crmul(hl_2, -1.0);
            hl_2 = ccdiv(hl_2, cmplx(pow(z, 6.0), 0.0));
            break;
        case 6:
            hl_2 = ccmul(cexp(cmplx(0.0, -z)), (cmplx(21.0*pow(z, 5.0) - 1260.0*pow(z, 3.0) + 10395.0*z, pow(z, 6.0) -
                                                      pow(z, 4.0)*210.0 + pow(z, 2.0)*4725.0 - 10395.0)));
            hl_2 = crmul(hl_2, -1.0);
            hl_2 = ccdiv(hl_2, cmplx(pow(z, 7.0), 0.0));
            break;
        case 7:
            hl_2 = ccmul(cexp(cmplx(0.0, -z)), (cmplx(-pow(z, 7.0) + 378.0*pow(z, 5.0) - 17325.0*pow(z, 3.0) + 135135.0*z,
                                                      pow(z, 6.0)*28.0 - pow(z, 4.0)*3150.0 + pow(z, 2.0)*62370.0 - 135135.0)));
            hl_2 = crmul(hl_2, -1.0);
            hl_2 = ccdiv(hl_2, cmplx(pow(z, 8.0), 0.0));
            break;
        case 8:
            hl_2 = ccmul(cexp(cmplx(0.0,-z)), (cmplx(36.0*pow(z,7.0) - 6930.0*pow(z,5.0) + 270270.0*pow(z,3.0) - 2027025.0*z,
                                                     pow(z,8.0) - pow(z,6.0)*630.0  + pow(z,4.0)*51975.0 - pow(z,2.0)*945945.0
                                                     + 2027025.0)));
            hl_2 = ccdiv(hl_2, cmplx(pow(z, 9.0), 0.0));
            break;
        default:
            hl_2 = cmplx(2.23e-13, 0.0);
            break;
    }
    return hl_2;
}

static void array2sh_calculate_bN
(
    void* const hA2sh,
    float freqVector[HYBRID_BANDS]
)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    arrayPars* arraySpecs = (arrayPars*)(pData->arraySpecs);
    int band, n;
    double kr[HYBRID_BANDS], kR[HYBRID_BANDS];
    double sph2cyl_kR, admittance;
    admittance = arraySpecs->admittance;
    double_complex jl2, hl2, Jl2, Jl2_imag;
    
    for(band=0; band<HYBRID_BANDS; band++){
        kr[band] = 2.0*M_PI*(double)freqVector[band] * (double)arraySpecs->r / (double)pData->c;
        kR[band] = 2.0*M_PI*(double)freqVector[band] * (double)arraySpecs->R / (double)pData->c;
    }
    
    switch(arraySpecs->arrayType){
        case ARRAY_CYLINDRICAL:
            switch (arraySpecs->weightType){
                case WEIGHT_RIGID:
                    for(band=0; band<HYBRID_BANDS; band++){
                        for(n=0; n < SH_ORDER+1; n++){
                            if (n==0){
#ifdef __APPLE__
                                Jl2 = cmplx(-jn(1,kR[band])-jn(1,kR[band]),0.0);
                                Jl2_imag = cmplx(-jn(1,kR[band])- jn(1,kR[band]), -(-yn(1,kR[band]) - yn(1,kR[band])));
                                pData->bN[band][n] = crmul(cpow(cmplx(0.0,1.0),cmplx((double)n,0.0)),4.0*M_PI);
                                pData->bN[band][n] = ccmul( pData->bN[band][n], ccsub(cmplx(jn(0,kr[band]),0.0), ccmul(ccdiv(crmul(Jl2,0.5), (crmul(Jl2_imag,0.5))),  bessel_Hl2(0,kr[band]) )) );
#else
                                Jl2 = cmplx(-_jn(1,kR[band])-_jn(1,kR[band]),0.0);
                                Jl2_imag = cmplx(-_jn(1,kR[band])- _jn(1,kR[band]), -(-_yn(1,kR[band]) - _yn(1,kR[band])));
                                pData->bN[band][n] = crmul(cpow(cmplx(0.0,1.0),cmplx((double)n,0.0)),4.0*M_PI);
                                pData->bN[band][n] = ccmul( pData->bN[band][n], ccsub(cmplx(_jn(0,kr[band]),0.0), ccmul(ccdiv(crmul(Jl2,0.5), (crmul(Jl2_imag,0.5))),  bessel_Hl2(0,kr[band]) )) );
#endif
                            }
                            else {
#ifdef __APPLE__
                                Jl2 = cmplx(jn(n-1,kR[band])-jn(n+1,kR[band]),0.0);
                                Jl2_imag = cmplx(jn(n-1,kR[band])- jn(n+1,kR[band]), -(yn(n-1,kR[band]) - yn(n+1,kR[band])));
                                pData->bN[band][n] = crmul(cpow(cmplx(0.0,1.0),cmplx((double)n,0.0)),4.0*M_PI);
                                pData->bN[band][n] = ccmul( pData->bN[band][n], ccsub(cmplx(jn(n,kr[band]),0.0), ccmul(ccdiv(crmul(Jl2,0.5), (crmul(Jl2_imag,0.5))),  bessel_Hl2(n,kr[band]) )) );
#else
                                Jl2 = cmplx(_jn(n-1,kR[band])-_jn(n+1,kR[band]),0.0);
                                Jl2_imag = cmplx(_jn(n-1,kR[band])- _jn(n+1,kR[band]), -(_yn(n-1,kR[band]) - _yn(n+1,kR[band])));
                                pData->bN[band][n] = crmul(cpow(cmplx(0.0,1.0),cmplx((double)n,0.0)),4.0*M_PI);
                                pData->bN[band][n] = ccmul( pData->bN[band][n], ccsub(cmplx(_jn(n,kr[band]),0.0), ccmul(ccdiv(crmul(Jl2,0.5), (crmul(Jl2_imag,0.5))),  bessel_Hl2(n,kr[band]) )) );
#endif
                            }
                        }
                    }
                    break;
                    
                case WEIGHT_OPEN_OMNI:
                    for(band=0; band<HYBRID_BANDS; band++){
                        for(n=0; n < SH_ORDER+1; n++){
#ifdef __APPLE__
                            pData->bN[band][n] = crmul(cpow(cmplx(0.0,1.0),cmplx((double)n,0.0)), 4.0 * M_PI * jn(n, kr[band]));
#else
                            pData->bN[band][n] = crmul(cpow(cmplx(0.0,1.0),cmplx((double)n,0.0)), 4.0 * M_PI * _jn(n, kr[band]));
#endif
                        }
                    }
                    break;
                    
                 case WEIGHT_OPEN_CARD:
                    /* not supported */
                    break;
                    
                case WEIGHT_OPEN_DIPOLE:
                    /* not supported */
                    break;
            }
            break;
            
        case ARRAY_SPHERICAL:
            switch (arraySpecs->weightType){
                case WEIGHT_OPEN_OMNI:
                    for(band=0; band<HYBRID_BANDS; band++){
                        for(n=0; n < SH_ORDER+1; n++){
							pData->bN[band][n] = crmul(cpow(cmplx(0.0, 1.0), cmplx((double)n, 0.0)), 4.0*M_PI* bessel_jl(n, kr[band]));
                        }
                    }
                    break;
                    
                case WEIGHT_OPEN_CARD:
                    for(band=0; band<HYBRID_BANDS; band++){
                        for(n=0; n < SH_ORDER+1; n++){
							pData->bN[band][n] = crmul(cpow(cmplx(0.0, 1.0), cmplx((double)n, 0.0)), 4.0*M_PI);
							pData->bN[band][n] = ccmul(pData->bN[band][n], ccsub(cmplx(bessel_jl(n, kr[band]),0.0), crmul(cmplx(0.0f, 0.5), (bessel_jl(n - 1, kr[band]) - bessel_jl(n + 1, kr[band]) - 1.0 / kr[band] * bessel_jl(n, kr[band]))) ) );
                        }
                    }
                    break;
                    
                case WEIGHT_OPEN_DIPOLE:
                    for(band=0; band<HYBRID_BANDS; band++){
                        for(n=0; n < SH_ORDER+1; n++){
							pData->bN[band][n] = crmul(cpow(cmplx(0.0, 1.0), cmplx((double)n, 0.0)), 4.0*M_PI);
							pData->bN[band][n] = ccmul(pData->bN[band][n],  crmul( crmul(cmplx(0.0f, 0.5), -1.0), (bessel_jl(n - 1, kr[band]) - bessel_jl(n + 1, kr[band]) - 1.0 / kr[band] * bessel_jl(n, kr[band])) ) );
                        }
                    }
                    break;
                    
                case WEIGHT_RIGID:
                    for(band=0; band<HYBRID_BANDS; band++){
                        for(n=0; n < SH_ORDER+1; n++){
							sph2cyl_kR = MAX(sqrt(M_PI / (2.0*kR[band])), 2.23e-13);
							jl2 = ccadd(cmplx((bessel_jl(n - 1, kR[band]) - bessel_jl(n + 1, kR[band])) / sph2cyl_kR, 0.0), crmul((ccsub(cmplx(0.0, 2.0*admittance), cmplx(1.0 / kR[band],0.0))), bessel_jl(n, kR[band]) / sph2cyl_kR) );
							hl2 = ccdiv(ccsub(bessel_hl2(n - 1, kR[band]), bessel_hl2(n + 1, kR[band])), cmplx(sph2cyl_kR,0.0));
							hl2 = ccadd( hl2, ccmul((ccsub(cmplx(0.0, 2.0*admittance), cmplx( 1.0 / kR[band],0.0))), ccdiv(bessel_hl2(n, kR[band]), cmplx(sph2cyl_kR, 0.0))) );

							pData->bN[band][n] = crmul(cpow(cmplx(0.0, 1.0), cmplx((double)n, 0.0)), 4.0*M_PI);
							pData->bN[band][n] = ccmul(pData->bN[band][n], (ccsub(cmplx(bessel_jl(n, kr[band]),0.0), ccmul(ccdiv(jl2, hl2), bessel_hl2(n, kr[band])))));
                        }
                    }
                    break;
            }
            break;
    }
    
    for(band=0; band<HYBRID_BANDS; band++)
        for(n=0; n < SH_ORDER+1; n++)
            pData->bN[band][n] = ccdiv(pData->bN[band][n], cmplx(4.0*M_PI,0.0));
}


static void array2sh_reg_inv_bN
(
    void* const hA2sh
)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    arrayPars* arraySpecs = (arrayPars*)(pData->arraySpecs);
    int band, n;
    double regPar, g_lim, alpha, beta, f_n;
    
    regPar = (double)pData->regPar;
    
    for(band=0; band<HYBRID_BANDS; band++)
        for(n=0; n < SH_ORDER+1; n++)
            pData->bN_modal[band][n] = ccdiv(cmplx(1.0,0.0), (pData->bN[band][n]));
    
    switch(pData->regType){
        case REG_DAS:
            for(band=0; band<HYBRID_BANDS; band++){
                f_n = 0.0;
                for(n=0; n < SH_ORDER+1; n++)
                    f_n += (2.0*(double)n+1.0) * pow(cabs(pData->bN[band][n]), 2.0);
                beta = (1.0/ pow((double)(SH_ORDER)+1.0,2.0)) * f_n;
                for(n=0; n < SH_ORDER+1; n++)
					pData->bN_inv[band][n] = crmul(pData->bN_modal[band][n], pow(cabs(pData->bN[band][n]), 2.0) / beta);
            }
            break;
            
        case REG_SOFT_LIM:
            g_lim = sqrt(arraySpecs->Q)*pow(10.0,(regPar/20.0));
            for(band=0; band<HYBRID_BANDS; band++)
                for(n=0; n < SH_ORDER+1; n++)
                    pData->bN_inv[band][n] = crmul(pData->bN_modal[band][n], (2.0*g_lim*cabs(pData->bN[band][n]) / M_PI) * atan(M_PI / (2.0*g_lim*cabs(pData->bN[band][n]))) );
            break;
            
        case REG_TIKHONOV:
            alpha = sqrt(arraySpecs->Q)*pow(10.0,(regPar/20.0));
            for(band=0; band<HYBRID_BANDS; band++){
                for(n=0; n < SH_ORDER+1; n++){
                    beta = sqrt((1.0-sqrt(1.0-1.0/ pow(alpha,2.0)))/(1.0+sqrt(1.0-1.0/pow(alpha,2.0)))); /* Moreau & Daniel */
                    pData->bN_inv[band][n] = ccdiv(conj(pData->bN[band][n]), cmplx((pow(cabs(pData->bN[band][n]), 2.0) + pow(beta, 2.0)),0.0));
                }
            }
            break;
    }
}

static void array2sh_replicate_order
(
    void* const hA2sh
)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    int band, n, i;
    int o[SH_ORDER+2];
    
    for(n=0; n<SH_ORDER+2; n++)
        o[n] = n*n;
    for(band=0; band<HYBRID_BANDS; band++)
        for(n=0; n < SH_ORDER+1; n++)
            for(i=o[n]; i < o[n+1]; i++)
                pData->bN_inv_R[band][i] = pData->bN_inv[band][n];
}


void array2sh_initTFT
(
    void* const hA2sh
)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    arrayPars* arraySpecs = (arrayPars*)(pData->arraySpecs);
    int t, ch;
    
    if (pData->hSTFT != NULL){
        afSTFTfree(pData->hSTFT);
        pData->hSTFT = NULL;
        for (t = 0; t<TIME_SLOTS; t++) {
            for (ch = 0; ch< arraySpecs->Q; ch++) {
                free(pData->STFTInputFrameTF[t][ch].re);
                free(pData->STFTInputFrameTF[t][ch].im);
            }
        }
        free2d((void**)pData->STFTInputFrameTF, TIME_SLOTS);
        free2d((void**)pData->tempHopFrameTD, arraySpecs->Q);
    }
    if (pData->hSTFT == NULL){
        afSTFTinit(&(pData->hSTFT), HOP_SIZE, arraySpecs->newQ, NUM_SH_SIGNALS, 1, 1);
        pData->STFTInputFrameTF = (complexVector**)malloc2d(TIME_SLOTS, arraySpecs->newQ, sizeof(complexVector));
        for(t=0; t<TIME_SLOTS; t++) {
            for(ch=0; ch< arraySpecs->newQ; ch++) {
                pData->STFTInputFrameTF[t][ch].re = (float*)calloc(HYBRID_BANDS, sizeof(float));
                pData->STFTInputFrameTF[t][ch].im = (float*)calloc(HYBRID_BANDS, sizeof(float));
            } 
        }
        pData->tempHopFrameTD = (float**)malloc2d(arraySpecs->newQ, HOP_SIZE, sizeof(float));
        arraySpecs->Q = arraySpecs->newQ;
        pData->reinitSHTmatrixFLAG = 1; /* filters need to be updated too */
    }
}

void array2sh_calculate_sht_matrix
(
    void* const hA2sh
)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    arrayPars* arraySpecs = (arrayPars*)(pData->arraySpecs);
    int band, q, s;
    int n, nrhs, lda, ldb, info;
    int ipiv[NUM_SH_SIGNALS];
    float* steerTmp;
    float_complex diag_bN_inv_R[NUM_SH_SIGNALS][NUM_SH_SIGNALS];
    const float alpha = 1.0, beta = 0.0;
    const float_complex calpha = cmplxf(1.0f, 0.0f); const float_complex cbeta  = cmplxf(0.0f, 0.0f);
    double Y_col[NUM_SH_SIGNALS*MAX_NUM_SENSORS];
    double* Ytmp, *pinvY;
    
    steerTmp = (float*)malloc(NUM_SH_SIGNALS*sizeof(float));
    arraySpecs->R = MAX(arraySpecs->R, arraySpecs->r);

    /* calculate equalisation matrix */
    array2sh_calculate_bN(hA2sh, pData->freqVector);
    array2sh_reg_inv_bN(hA2sh);
    array2sh_replicate_order(hA2sh);
    
    /* calculate real-valued SH weights for each sensor direction */
    memset(Y_col, 0, NUM_SH_SIGNALS*MAX_NUM_SENSORS*sizeof(double));
    for (s = 0; s < NUM_SH_SIGNALS; s++) {
        memset(pData->Y[s], 0, MAX_NUM_SENSORS*sizeof(double));
        memset(pData->Y_cmplx[s], 0, MAX_NUM_SENSORS*sizeof(float_complex));
        memset(diag_bN_inv_R[s], 0, NUM_SH_SIGNALS*sizeof(float_complex));
    }
    Ytmp = malloc(NUM_SH_SIGNALS * (arraySpecs->Q) * sizeof(double));
    pinvY = malloc((arraySpecs->Q) * NUM_SH_SIGNALS * sizeof(double));
    for (q=0; q < arraySpecs->Q; q++){
        getSHreal(SH_ORDER, arraySpecs->sensorCoords_rad[q][0], M_PI / 2.0f - arraySpecs->sensorCoords_rad[q][1], steerTmp);
        for (s = 0; s < NUM_SH_SIGNALS; s++) {
            pData->Y[s][q] = (double)steerTmp[s] * sqrt(4.0*M_PI);
            Y_col[NUM_SH_SIGNALS*q + s] = pData->Y[s][q]; /* store col major order */
            //Ytmp[s * (arraySpecs->Q) + q] = pData->Y[s][q];
        }
    }
    
    //utility_dpinv(Ytmp, NUM_SH_SIGNALS, arraySpecs->Q, pinvY);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, NUM_SH_SIGNALS, NUM_SH_SIGNALS, MAX_NUM_SENSORS, alpha,
               (const double*)pData->Y, MAX_NUM_SENSORS,
               (const double*)pData->Y, MAX_NUM_SENSORS, beta,
               (double*)pData->YYT, NUM_SH_SIGNALS);
    n = NUM_SH_SIGNALS; nrhs = MAX_NUM_SENSORS; lda = NUM_SH_SIGNALS; ldb = NUM_SH_SIGNALS;
#ifdef __APPLE__
    dgesv_( (__CLPK_integer*)&n, (__CLPK_integer*)&nrhs, (__CLPK_doublereal*)pData->YYT, (__CLPK_integer*)&lda,
           (__CLPK_integer*)ipiv, (__CLPK_doublereal*)Y_col, (__CLPK_integer*)&ldb, (__CLPK_integer*)&info );
#else
    dgesv_( &n, &nrhs, (double*)pData->YYT, &lda, ipiv, (double*)Y_col, &ldb, &info );
#endif
    if( info > 0 ) {
        /* matrix pData->YYT is singular */
        for (q=0; q < arraySpecs->Q; q++){
            for (s = 0; s < NUM_SH_SIGNALS; s++) {
                /* Scale by number of sensors (assume a uniform sensor arrangment) */
                pData->Y_cmplx[s][q] = cmplxf((float)pData->Y[s][q]/(float)arraySpecs->Q, 0.0f);
            }
        }
    }
    else{
        for (q=0; q < arraySpecs->Q; q++){
            for (s = 0; s < NUM_SH_SIGNALS; s++) {
                pData->Y_cmplx[s][q] = cmplxf((float)Y_col[NUM_SH_SIGNALS*q + s], 0.0f); /* return to row major */
                //pData->Y_cmplx[s][q] = cmplxf((float)pinvY[NUM_SH_SIGNALS*q + s], 0.0f); /* transpose */
            }
        }
    }

    /* calculate spherical harmonic transform matrix */
    for (s = 0; s < NUM_SH_SIGNALS; s++)
        memset(pData->W[0][s], 0, MAX_NUM_SENSORS*sizeof(float_complex)); /* ignore DC */
    for(band = 1; band <HYBRID_BANDS; band++){
        for (s = 0; s < NUM_SH_SIGNALS; s++)
            diag_bN_inv_R[s][s] = cmplxf((float)creal(pData->bN_inv_R[band][s]), (float)cimag(pData->bN_inv_R[band][s]));
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, NUM_SH_SIGNALS, MAX_NUM_SENSORS, NUM_SH_SIGNALS, &calpha,
                    diag_bN_inv_R, NUM_SH_SIGNALS,
                    pData->Y_cmplx, MAX_NUM_SENSORS, &cbeta,
                    pData->W[band], MAX_NUM_SENSORS);
    }
    free((void*)steerTmp);
    free(Ytmp);
    free(pinvY);
}


void array2sh_calculate_mag_curves(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    int band, n;
    
    for(band = 0; band <HYBRID_BANDS-1; band++){
        for(n = 0; n <SH_ORDER+1; n++){
            pData->bN_inv_dB[band][n] = 20.0f * (float)log10(cabs(pData->bN_inv[band+1][n])); /* Ignore DC */
            pData->bN_modal_dB[band][n] = 20.0f * (float)log10(cabs(pData->bN_modal[band+1][n]));
        }
    }
}

void array2sh_evaluateSHTfilters(void* hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    arrayPars* arraySpecs = (arrayPars*)(pData->arraySpecs);
    int band, i, j, simOrder;
    double kr[HYBRID_BANDS-1];
    float* Y_grid_real;
    float_complex* Y_grid, *H_array, *Wshort;
    
    assert(pData->W != NULL);
    
    /* simulate the current array by firing 812 plane-waves around the surface of a theoretical sphere
     * and ascertaining the transfer function for each */
    simOrder = (int)(2.0f*M_PI*MAX_EVAL_FREQ_HZ*(arraySpecs->R)/pData->c)+1;
    for(i=0; i<HYBRID_BANDS-1; i++)
        kr[i] = 2.0*M_PI*(pData->freqVector[i+1/* ignore DC */])*(arraySpecs->R)/pData->c;
    H_array = malloc((HYBRID_BANDS-1) * (arraySpecs->Q) * 812*sizeof(float_complex));
    switch(arraySpecs->arrayType){
        case ARRAY_SPHERICAL:
            switch(arraySpecs->weightType){
                default:
                case WEIGHT_RIGID:
                    simulateSphArray(simOrder, kr, HYBRID_BANDS-1, (float*)arraySpecs->sensorCoords_rad, arraySpecs->Q, (float*)__geosphere_ico_9_0_dirs_deg, 812, ARRAY_CONSTRUCTION_RIGID, 0.0, H_array);
                    break;
                case WEIGHT_OPEN_OMNI:
                    simulateSphArray(simOrder, kr, HYBRID_BANDS-1, (float*)arraySpecs->sensorCoords_rad, arraySpecs->Q, (float*)__geosphere_ico_9_0_dirs_deg, 812, ARRAY_CONSTRUCTION_OPEN, 0.0, H_array);
                    break;
                case WEIGHT_OPEN_DIPOLE:
                    simulateSphArray(simOrder, kr, HYBRID_BANDS-1, (float*)arraySpecs->sensorCoords_rad, arraySpecs->Q, (float*)__geosphere_ico_9_0_dirs_deg, 812, ARRAY_CONSTRUCTION_DIRECTIONAL, 1.0, H_array);
                    break;
                case WEIGHT_OPEN_CARD:
                    simulateSphArray(simOrder, kr, HYBRID_BANDS-1, (float*)arraySpecs->sensorCoords_rad, arraySpecs->Q, (float*)__geosphere_ico_9_0_dirs_deg, 812, ARRAY_CONSTRUCTION_DIRECTIONAL, 0.5, H_array);
                    break;
            }
            break;
            
        case ARRAY_CYLINDRICAL:
            switch(arraySpecs->weightType){
                default:
                case WEIGHT_RIGID:
                    simulateCylArray(simOrder, kr, HYBRID_BANDS-1, (float*)arraySpecs->sensorCoords_rad, arraySpecs->Q, (float*)__geosphere_ico_9_0_dirs_deg, 812, ARRAY_CONSTRUCTION_RIGID, H_array);
                    break;
                case WEIGHT_OPEN_DIPOLE:
                case WEIGHT_OPEN_CARD:
                case WEIGHT_OPEN_OMNI:
                    simulateCylArray(simOrder, kr, HYBRID_BANDS-1, (float*)arraySpecs->sensorCoords_rad, arraySpecs->Q, (float*)__geosphere_ico_9_0_dirs_deg, 812, ARRAY_CONSTRUCTION_OPEN, H_array);
                    break;
            }
            break;
    }
    
    /* generate ideal (real) spherical harmonics to compare with */
    Y_grid_real = NULL;
    getRSH(SH_ORDER, (float*)__geosphere_ico_9_0_dirs_deg, 812, &Y_grid_real);
    Y_grid = malloc(NUM_SH_SIGNALS*812*sizeof(float_complex));
    for(i=0; i<NUM_SH_SIGNALS*812; i++)
        Y_grid[i] = cmplxf(Y_grid_real[i], 0.0f); /* "evaluateSHTfilters" function requires complex data type */
    
    /* compare the spherical harmonics obtained from encoding matrix 'W' with the ideal patterns */
    Wshort = malloc(HYBRID_BANDS*NUM_SH_SIGNALS*(arraySpecs->Q)*sizeof(float_complex));
    for(band=0; band<HYBRID_BANDS-1; band++)
        for(i=0; i<NUM_SH_SIGNALS; i++)
            for(j=0; j<(arraySpecs->Q); j++)
                Wshort[band*NUM_SH_SIGNALS*(arraySpecs->Q) + i*(arraySpecs->Q) + j] = pData->W[band+1/* skip DC */][i][j];
    evaluateSHTfilters(SH_ORDER, Wshort, arraySpecs->Q, HYBRID_BANDS-1, H_array, 812, Y_grid, pData->cSH, pData->lSH);

    free(Y_grid_real);
    free(Y_grid);
    free(H_array);
    free(Wshort);
}

void array2sh_createArray(void ** const hPars)
{
    arrayPars* pars = (arrayPars*)malloc(sizeof(arrayPars));
    if (pars == NULL) { return;/*error*/ }
    *hPars = (void*)pars;
}

void array2sh_destroyArray(void ** const hPars)
{
    arrayPars *pars = (arrayPars*)(*hPars);
    
    if(pars!=NULL) {
        free(pars);
        pars=NULL;
    }
}
 
void array2sh_initArray(void* const hPars, PRESETS preset, int firstInitFlag)
{
    arrayPars *pars = (arrayPars*)(hPars);
    int ch, i, Q;
    
    switch(preset){
	    default:
        case PRESET_DEFAULT:
            Q = (SH_ORDER+1)*(SH_ORDER+1); /* number of mics */
            pars->r = 0.042f; /* array radius */
            pars->R = 0.042f; /* radius of the sensors (incase they protrude from the surface of the array), (only for rigid arrays) */
            pars->admittance = 0.0f; /* acoustical admittance of the array material (only for rigid arrays) */
            pars->arrayType = ARRAY_SPHERICAL; /* spherical or cylindrical */
            pars->weightType = WEIGHT_RIGID; /* open or rigid, and directivity of the sensors (only for open arrays) */
            for(ch=0; ch<Q; ch++){
                for(i=0; i<2; i++){
                    pars->sensorCoords_rad[ch][i] = __default_coords_rad[ch][i]; /* spherical coordinates of the sensors, in radians */
                    pars->sensorCoords_deg[ch][i] = pars->sensorCoords_rad[ch][i] * (180.0f/M_PI);
                }
            }
            break;
#ifdef ENABLE_AALTO_HYDROPHONE_PRESET
        case PRESET_AALTO_HYDROPHONE:
            Q = 4;
            pars->r = 0.1555f;
            pars->R = 0.1555f;
            pars->arrayType = ARRAY_SPHERICAL;
            pars->weightType = WEIGHT_OPEN_OMNI;
            for(ch=0; ch<Q; ch++){
                for(i=0; i<2; i++){
                    pars->sensorCoords_rad[ch][i] = __Aalto_Hydrophone_coords_rad[ch][i];
                    pars->sensorCoords_deg[ch][i] = pars->sensorCoords_rad[ch][i] * (180.0f/M_PI);
                }
            }
            break;
#endif
#ifdef ENABLE_SENNHEISER_AMBEO_PRESET
        case PRESET_SENNHEISER_AMBEO:
            Q = 4;
            pars->r = 0.02f;
            pars->R = 0.02f;
            pars->arrayType = ARRAY_SPHERICAL;
            pars->weightType = WEIGHT_OPEN_CARD;
            for(ch=0; ch<Q; ch++){
                for(i=0; i<2; i++){
                    pars->sensorCoords_rad[ch][i] = __Sennheiser_Ambeo_coords_rad[ch][i];
                    pars->sensorCoords_deg[ch][i] = pars->sensorCoords_rad[ch][i] * (180.0f/M_PI);
                }
            }
            break;
#endif
#ifdef ENABLE_CORE_SOUND_TETRAMIC_PRESET
        case PRESET_CORE_SOUND_TETRAMIC:
            Q = 4;
            pars->r = 0.02f;
            pars->R = 0.02f;
            pars->arrayType = ARRAY_SPHERICAL;
            pars->weightType = WEIGHT_OPEN_CARD;
            for(ch=0; ch<Q; ch++){
                for(i=0; i<2; i++){
                    pars->sensorCoords_rad[ch][i] = __Core_Sound_TetraMic_coords_rad[ch][i];
                    pars->sensorCoords_deg[ch][i] = pars->sensorCoords_rad[ch][i] * (180.0f/M_PI);
                }
            }
            break;
#endif
#ifdef ENABLE_SOUND_FIELD_SPS200_PRESET
        case PRESET_SOUND_FIELD_SPS200:
            Q = 4;
            pars->r = 0.02f;
            pars->R = 0.02f;
            pars->arrayType = ARRAY_SPHERICAL;
            pars->weightType = WEIGHT_OPEN_CARD;
            for(ch=0; ch<Q; ch++){
                for(i=0; i<2; i++){
                    pars->sensorCoords_rad[ch][i] = __Sound_field_SPS200_coords_rad[ch][i];
                    pars->sensorCoords_deg[ch][i] = pars->sensorCoords_rad[ch][i] * (180.0f/M_PI);
                }
            }
            break;
#endif
#ifdef ENABLE_ZYLIA_1D_PRESET
        case PRESET_ZYLIA_1D:
            Q = 19;
            pars->r = 0.049f;
            pars->R = 0.049f;
            pars->arrayType = ARRAY_SPHERICAL;
            pars->weightType = WEIGHT_RIGID;
            for(ch=0; ch<Q; ch++){
                for(i=0; i<2; i++){
                    pars->sensorCoords_rad[ch][i] = __Zylia1D_coords_rad[ch][i];
                    pars->sensorCoords_deg[ch][i] = pars->sensorCoords_rad[ch][i] * (180.0f/M_PI);
                }
            }
            break;
#endif
#ifdef ENABLE_EIGENMIKE32_PRESET
        case PRESET_EIGENMIKE32:
            Q = 32;
            pars->r = 0.042f;
            pars->R = 0.042f;
            pars->admittance = 0.0f;
            pars->arrayType = ARRAY_SPHERICAL;
            pars->weightType = WEIGHT_RIGID;
            for(ch=0; ch<Q; ch++){
                for(i=0; i<2; i++){
                    pars->sensorCoords_rad[ch][i] = __Eigenmike32_coords_rad[ch][i];
                    pars->sensorCoords_deg[ch][i] = pars->sensorCoords_rad[ch][i] * (180.0f/M_PI);
                }
            }
            break;
#endif
#ifdef ENABLE_DTU_MIC_PRESET
        case PRESET_DTU_MIC:
            Q = 52;
            pars->r = 0.05f;
            pars->R = 0.05f;
            pars->admittance = 0.0f;
            pars->arrayType = ARRAY_SPHERICAL;
            pars->weightType = WEIGHT_RIGID;
            for(ch=0; ch<Q; ch++){
                for(i=0; i<2; i++){
                    pars->sensorCoords_rad[ch][i] = __DTU_mic_coords_rad[ch][i];
                    pars->sensorCoords_deg[ch][i] = pars->sensorCoords_rad[ch][i] * (180.0f/M_PI);
                }
            }
            break;
#endif 
    }
    
    /* Fill remaining slots with default coords */
    for(; ch<MAX_NUM_SENSORS; ch++){
        for(i=0; i<2; i++){
            pars->sensorCoords_rad[ch][i] = __default_SENSORcoords64_rad[ch][i];
            pars->sensorCoords_deg[ch][i] = pars->sensorCoords_rad[ch][i] * (180.0f/M_PI);
        }
    }
    
    /* For dynamically changing the number of TFT channels */
    if(firstInitFlag==1){
        pars->Q = Q;
        pars->newQ = pars->Q;
    }
    else
        pars->newQ = Q;
}




