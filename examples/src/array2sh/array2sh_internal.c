/*
 * Copyright 2017-2018 Leo McCormack
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
 * @file array2sh_internal.c
 * @brief Spatially encodes spherical microphone array signals into spherical
 *        harmonic signals (aka: Ambisonic signals) utilising theoretical
 *        encoding filters.
 *
 * The algorithms within array2sh were pieced together and developed in
 * collaboration with Symeon Delikaris-Manias and Angelo Farina.
 * A detailed explanation of the algorithms within array2sh can be found in [1].
 * Also included, is a diffuse-field equalisation option for frequencies past
 * aliasing, developed in collaboration with Archontis Politis, 8.02.2019
 *
 * @note Since the algorithms are based on theory, only array designs where
 *       there are analytical solutions available are supported. i.e. only
 *       spherical or cylindrical arrays, which have phase-matched sensors.
 *       For more information, the reader is referred to [2,3].
 * @test test__saf_example_array2sh()
 *
 * @see [1] McCormack, L., Delikaris-Manias, S., Farina, A., Pinardi, D., and
 *          Pulkki, V., "Real-time conversion of sensor array signals into
 *          spherical harmonic signals with applications to spatially localised
 *          sub-band sound-field analysis," in Audio Engineering Society
 *          Convention 144, Audio Engineering Society, 2018.
 * @see [2] Williams EG. Fourier acoustics: sound radiation and nearfield
 *          acoustical holography. Elsevier; 1999 Jun 10.
 * @see [3] Rafaely B. Fundamentals of spherical array processing. Berlin:
 *          Springer; 2015 Feb 18.
 *
 * @author Leo McCormack
 * @date 13.09.2017
 * @license ISC
 */

#include "array2sh_internal.h"

const int array2sh_defaultNumSensors = 4;

const float array2sh_defaultSensorsDirections[MAX_NUM_INPUTS][2] =
  { { 45.0f,    35.264f},
    { -45.0f,    -35.264f},
    { 135.0f,    -35.264f},
    { -135.0f,    35.264f},
    { -59.5674172539630f,    -69.7763953095265f},
    { -107.629493132161f,    10.0801907586439f},
    { 59.5651193014662f,    69.7750019669784f},
    { -72.3715136529922f,    -10.0909592824967f},
    { 120.418832045969f,    -69.7630732527928f},
    { 72.3705387467529f,    10.0908849630425f},
    { 10.5699823270029f,    17.3489915451677f},
    { -169.426851150349f,    17.3476277272366f},
    { -30.7728814706980f,    68.2546206344998f},
    { 101.529181734204f,    18.5653877408229f},
    { 149.250234431564f,    68.2619792676765f},
    { 18.9183978612362f,    -10.9197242075287f},
    { 78.4729799557400f,    -18.5558892226653f},
    { -18.9161109968480f,    10.9195275229078f},
    { -78.4723066816396f,    18.5565169424192f},
    { -161.085243560336f,    -10.9220975953767f},
    { -101.530063994217f,    -18.5648176173994f},
    { 161.085855454581f,    10.9237834974412f},
    { 30.7790837424989f,    -68.2562996173874f},
    { -149.255648697558f,    -68.2594892022576f},
    { 56.4638592172755f,    41.2645742763346f},
    { 46.4626904606488f,    24.5386896277818f},
    { -123.528839291756f,    41.2565889648213f},
    { 32.1932295461248f,    38.8027347956192f},
    { 133.534545076742f,    -24.5310453346374f},
    { -32.1935458377634f,    -38.8041545028931f},
    { -133.533477606995f,    24.5319420428588f},
    { -147.797881221742f,    38.7980990097138f},
    { -46.4626630742447f,    -24.5400749997299f},
    { 147.799544204027f,    -38.7965984357137f},
    { -56.4649576930965f,    -41.2658215291547f},
    { 123.530547454259f,    -41.2555316199835f},
    { 84.7385596473317f,    27.3132538556364f},
    { 27.4072860474693f,    4.67593317316088f},
    { -95.2612516264878f,    27.3036967398252f},
    { 10.0607591758466f,    62.2287426217158f},
    { 152.592710471291f,    -4.67007594116539f},
    { -10.0604693351594f,    -62.2317767049329f},
    { -152.592615409799f,    4.67188744845125f},
    { -169.918387889403f,    62.2272712110016f},
    { -27.4047652820400f,    -4.67684895484045f},
    { 169.918493477140f,    -62.2246416140826f},
    { -84.7390238025025f,    -27.3120311727342f},
    { 95.2622793653641f,    -27.3027426350066f},
    { 136.272846725200f,    -0.726139068263599f},
    { -1.05456977390436f,    -46.2698051615206f},
    { -43.7263256019408f,    -0.737666959763443f},
    { -91.0084089447853f,    43.7170295842697f},
    { -178.936314291048f,    46.2682278552731f},
    { 91.0100541861312f,    -43.7163112498585f},
    { 178.938490558618f,    -46.2666218304286f},
    { 88.9860902890987f,    43.7285680951606f},
    { 1.05095788619312f,    46.2682407577363f},
    { -88.9851906686018f,    -43.7278600822776f},
    { -136.273314433943f,    0.726484363994244f},
    { 43.7273191822901f,    0.733967099281741f},
    { 55.2321077221357f,    10.8187456245361f},
    { 13.0867814081219f,    34.0654299770448f},
    { -124.767898552651f,    10.8099467408727f},
    { 71.4751302534919f,    53.7984109774923f},
    { 166.905549455530f,    -34.0625862919868f},
    { -71.4765158530160f,    -53.8002805142714f},
    { -166.904587570338f,    34.0633858486074f},
    { -108.515564585463f,    53.7864186881779f},
    { -13.0869589355977f,    -34.0644403811701f},
    { 108.517274377512f,    -53.7853093207457f},
    { -55.2337689585088f,    -10.8200079671294f},
    { 124.767519417537f,    -10.8096092876561f},
    { -105.493359369733f,    -68.1345228940484f},
    { -111.151973933085f,    -5.71277564515985f},
    { 74.5244447329726f,    -68.1205741555505f},
    { -173.885132403874f,    -21.0401043953034f},
    { -68.8521018079939f,    5.70313885775370f},
    { 173.885050918764f,    21.0419786777516f},
    { 68.8513569996117f,    -5.70297394874866f},
    { 6.12041448982713f,    -21.0387333986971f},
    { 111.151515131908f,    5.71357243375474f},
    { -6.11913621739966f,    21.0397756114589f},
    { 105.487055344890f,    68.1331981897040f},
    { -74.5188619841755f,    68.1214716326340f},
    { 35.2822222163987f,    -15.1801841162885f},
    { -25.1732516880600f,    51.9826514887145f},
    { -144.722913967824f,    -15.1859965048827f},
    { 108.386507250934f,    33.8824333892584f},
    { -154.839566662519f,    -51.9856067440744f},
    { -108.385758587361f,    -33.8826426573830f},
    { 154.840828582605f,    51.9884126722963f},
    { -71.6110721306039f,    33.8714084774643f},
    { 25.1699762263473f,    -51.9843011897216f},
    { 71.6126173250777f,    -33.8705490198003f},
    { -35.2810594358470f,    15.1763090469469f},
    { 144.722474074749f,    15.1863112389402f},
    { -125.277803641269f,    -28.5561977838671f},
    { -146.317489805034f,    -30.4888082061806f},
    { 54.7220842889970f,    -28.5463137051834f},
    { -133.292600372594f,    -45.8199818171681f},
    { -33.6854220654234f,    30.4794807422719f},
    { 133.295169860873f,    45.8200619710773f},
    { 33.6830420572493f,    -30.4828689919959f},
    { 46.7098267598019f,    -45.8106991831118f},
    { 146.318736830907f,    30.4889588527258f},
    { -46.7107531920589f,    45.8077510168776f},
    { 125.278543193551f,    28.5567151663926f},
    { -54.7226679373581f,    28.5424675755949f},
    { -144.402199049652f,    54.7141577397485f},
    { 112.382216591782f,    -28.0100385212576f},
    { 35.5834085793474f,    54.7193853715013f},
    { -29.9155874729000f,    -19.6471282266303f},
    { 67.6163420587218f,    28.0211878299760f},
    { 29.9151078604395f,    19.6461577548916f},
    { -67.6164804090156f,    -28.0221553885221f},
    { 150.081679085309f,    -19.6396879622134f},
    { -112.381365317985f,    28.0105570315277f},
    { -150.081127534239f,    19.6412999777121f},
    { 144.402936805264f,    -54.7127597303097f},
    { -35.5837994140428f,    -54.7203576209200f},
    { 68.5348504928164f,    -52.8527273534571f},
    { -54.8218918027732f,    12.7570858971645f},
    { -111.464947340552f,    -52.8652964226689f},
    { 164.510735700738f,    34.1892226850382f},
    { -125.181243085238f,    -12.7710900404128f},
    { -164.510172672569f,    -34.1865270720665f},
    { 125.180742749445f,    12.7714932994183f},
    { -15.4967238792034f,    34.1855682977571f}};

/**
 * Takes the bNs computed up to N+1, and replicates them to be of length
 * (N+1)^2 (replicating the 1st order bNs 3 times, 2nd -> 5 times etc.)
 */
static void array2sh_replicate_order
(
    void* const hA2sh,
    int order
)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    int band, n, i;
    int o[MAX_SH_ORDER+2];
    
    for(n=0; n<order+2; n++)
        o[n] = n*n;
    for(band=0; band<HYBRID_BANDS; band++)
        for(n=0; n < order+1; n++)
            for(i=o[n]; i < o[n+1]; i++)
                pData->bN_inv_R[band][i] = pData->bN_inv[band][n];
}

void array2sh_initTFT
(
    void* const hA2sh
)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    array2sh_arrayPars* arraySpecs = (array2sh_arrayPars*)(pData->arraySpecs);
    int new_nSH, nSH;
    
    new_nSH = (pData->new_order+1)*(pData->new_order+1);
    nSH = (pData->order+1)*(pData->order+1);
    if(pData->hSTFT==NULL)
        afSTFT_create(&(pData->hSTFT), arraySpecs->newQ, new_nSH, HOP_SIZE, 0, 1, AFSTFT_BANDS_CH_TIME);
    else if(arraySpecs->newQ != arraySpecs->Q || nSH != new_nSH){
        afSTFT_channelChange(pData->hSTFT, arraySpecs->newQ, new_nSH);
        afSTFT_clearBuffers(pData->hSTFT); 
        pData->reinitSHTmatrixFLAG = 1; /* filters will need to be updated too */
    }
    arraySpecs->Q = arraySpecs->newQ;
}

void array2sh_calculate_sht_matrix
(
    void* const hA2sh
)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    array2sh_arrayPars* arraySpecs = (array2sh_arrayPars*)(pData->arraySpecs);
    int i, j, band, n, order, nSH;
    double alpha, beta, g_lim, regPar;
    double kr[HYBRID_BANDS], kR[HYBRID_BANDS];
    float sensorCoords_deg_local[MAX_NUM_SENSORS][2];
    float* Y_mic, *pinv_Y_mic;
    float_complex* pinv_Y_mic_cmplx, *diag_bN_inv_R;
    const float_complex calpha = cmplxf(1.0f, 0.0f); const float_complex cbeta  = cmplxf(0.0f, 0.0f);
    
    /* prep */
    order = pData->new_order;
    nSH = (order+1)*(order+1);
    arraySpecs->R = SAF_MIN(arraySpecs->R, arraySpecs->r);
    for(band=0; band<HYBRID_BANDS; band++){
        kr[band] = 2.0*SAF_PId*(pData->freqVector[band])*(arraySpecs->r)/pData->c;
        kR[band] = 2.0*SAF_PId*(pData->freqVector[band])*(arraySpecs->R)/pData->c;
    }
    for(i=0; i<arraySpecs->Q; i++){
        sensorCoords_deg_local[i][0] = arraySpecs->sensorCoords_deg[i][0];
        sensorCoords_deg_local[i][1] = arraySpecs->sensorCoords_deg[i][1];
    }
    
    /* Spherical harmponic weights for each sensor direction */
    Y_mic = malloc1d(nSH*(arraySpecs->Q)*sizeof(float));
    getRSH(order, (float*)sensorCoords_deg_local, arraySpecs->Q, Y_mic); /* nSH x Q */
    pinv_Y_mic = malloc1d( arraySpecs->Q * nSH *sizeof(float));
    utility_spinv(NULL, Y_mic, nSH, arraySpecs->Q, pinv_Y_mic);
    pinv_Y_mic_cmplx =  malloc1d((arraySpecs->Q) * nSH *sizeof(float_complex));
    for(i=0; i<(arraySpecs->Q)*nSH; i++)
        pinv_Y_mic_cmplx[i] = cmplxf(pinv_Y_mic[i], 0.0f);
    
    /* ------------------------------------------------------------------------------ */
    /* Encoding filters based on the regularised inversion of the modal coefficients: */
    /* ------------------------------------------------------------------------------ */
    if ( (pData->filterType==FILTER_SOFT_LIM) || (pData->filterType==FILTER_TIKHONOV) ){
        /* Compute modal responses */
        free(pData->bN);
        pData->bN = malloc1d((HYBRID_BANDS)*(order+1)*sizeof(double_complex));
        ARRAY2SH_ARRAY_TYPES arrayType = arraySpecs->arrayType;
        switch(arrayType){
            case ARRAY_CYLINDRICAL: {
                ARRAY2SH_WEIGHT_TYPES weightType = arraySpecs->weightType;
                switch (weightType){
                    case WEIGHT_RIGID_OMNI:   cylModalCoeffs(order, kr, HYBRID_BANDS, ARRAY_CONSTRUCTION_RIGID, pData->bN); break;
                    case WEIGHT_RIGID_CARD:   saf_print_error("weightType is not supported"); break;
                    case WEIGHT_RIGID_DIPOLE: saf_print_error("weightType is not supported"); break;
                    case WEIGHT_OPEN_OMNI:    cylModalCoeffs(order, kr, HYBRID_BANDS, ARRAY_CONSTRUCTION_OPEN, pData->bN);  break;
                    case WEIGHT_OPEN_CARD:    saf_print_error("weightType is not supported"); break;
                    case WEIGHT_OPEN_DIPOLE:  saf_print_error("weightType is not supported"); break;
                }
                break;
            }
            case ARRAY_SPHERICAL: {
                ARRAY2SH_WEIGHT_TYPES weightType = arraySpecs->weightType;
                switch (weightType){
                    case WEIGHT_OPEN_OMNI:   sphModalCoeffs(order, kr, HYBRID_BANDS, ARRAY_CONSTRUCTION_OPEN, 1.0, pData->bN); break;
                    case WEIGHT_OPEN_CARD:   sphModalCoeffs(order, kr, HYBRID_BANDS, ARRAY_CONSTRUCTION_OPEN_DIRECTIONAL, 0.5, pData->bN); break;
                    case WEIGHT_OPEN_DIPOLE: sphModalCoeffs(order, kr, HYBRID_BANDS, ARRAY_CONSTRUCTION_OPEN_DIRECTIONAL, 0.0, pData->bN); break;
                    case WEIGHT_RIGID_OMNI:
                    case WEIGHT_RIGID_CARD:
                    case WEIGHT_RIGID_DIPOLE:
                        /* if sensors are flushed with the rigid baffle: */
                        if(arraySpecs->R == arraySpecs->r )
                            sphModalCoeffs(order, kr, HYBRID_BANDS, ARRAY_CONSTRUCTION_RIGID, 1.0, pData->bN);

                        /* if sensors protrude from the rigid baffle: */
                        else{
                            if (arraySpecs->weightType == WEIGHT_RIGID_OMNI)
                                sphScattererModalCoeffs(order, kr, kR, HYBRID_BANDS, pData->bN);
                            else if (arraySpecs->weightType == WEIGHT_RIGID_CARD)
                                sphScattererDirModalCoeffs(order, kr, kR, HYBRID_BANDS, 0.5, pData->bN);
                            else if (arraySpecs->weightType == WEIGHT_RIGID_DIPOLE)
                                sphScattererDirModalCoeffs(order, kr, kR, HYBRID_BANDS, 0.0, pData->bN);
                        }
                        break;
                }
                break;
            }
        }
        
        for(band=0; band<HYBRID_BANDS; band++)
            for(n=0; n < order+1; n++)
                pData->bN[band*(order+1)+n] = ccdiv(pData->bN[band*(order+1)+n], cmplx(4.0*SAF_PId, 0.0)); /* 4pi term */

        /* direct inverse */
        regPar = pData->regPar;
        for(band=0; band<HYBRID_BANDS; band++)
            for(n=0; n < order+1; n++)
                pData->bN_modal[band][n] = ccdiv(cmplx(1.0,0.0), (pData->bN[band*(order+1)+n]));
        
        /* regularised inverse */
        if (pData->filterType == FILTER_SOFT_LIM){
            /* Bernschutz, B., Porschmann, C., Spors, S., Weinzierl, S., Versterkung, B., 2011. Soft-limiting der
             modalen amplitudenverst?rkung bei sph?rischen mikrofonarrays im plane wave decomposition verfahren.
             Proceedings of the 37. Deutsche Jahrestagung fur Akustik (DAGA 2011) */
            g_lim = sqrt(arraySpecs->Q)*pow(10.0,(regPar/20.0));
            for(band=0; band<HYBRID_BANDS; band++)
                for(n=0; n < order+1; n++)
                    pData->bN_inv[band][n] = crmul(pData->bN_modal[band][n], (2.0*g_lim*cabs(pData->bN[band*(order+1)+n]) / SAF_PId)
                                                     * atan(SAF_PId / (2.0*g_lim*cabs(pData->bN[band*(order+1)+n]))) );
        }
        else if(pData->filterType == FILTER_TIKHONOV){
            /* Moreau, S., Daniel, J., Bertet, S., 2006, 3D sound field recording with higher order ambisonics-objective
             measurements and validation of spherical microphone. In Audio Engineering Society Convention 120. */
            alpha = sqrt(arraySpecs->Q)*pow(10.0,(regPar/20.0));
            for(band=0; band<HYBRID_BANDS; band++){
                for(n=0; n < order+1; n++){
                    beta = sqrt((1.0-sqrt(1.0-1.0/ pow(alpha,2.0)))/(1.0+sqrt(1.0-1.0/pow(alpha,2.0))));
                    pData->bN_inv[band][n] = ccdiv(conj(pData->bN[band*(order+1)+n]), cmplx((pow(cabs(pData->bN[band*(order+1)+n]), 2.0) + pow(beta, 2.0)),0.0));
                }
            }
        }
        
        /* diag(filters) * Y */
        array2sh_replicate_order(hA2sh, order); /* replicate orders */
        
        diag_bN_inv_R = calloc1d(nSH*nSH, sizeof(float_complex));
        for(band=0; band<HYBRID_BANDS; band++){
            for(i=0; i<nSH; i++)
                diag_bN_inv_R[i*nSH+i] = cmplxf((float)creal(pData->bN_inv_R[band][i]), (float)cimag(pData->bN_inv_R[band][i]));
            cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nSH, (arraySpecs->Q), nSH, &calpha,
                        diag_bN_inv_R, nSH,
                        pinv_Y_mic_cmplx, nSH, &cbeta,
                        pData->W[band], MAX_NUM_SENSORS);
        }
        free(diag_bN_inv_R);
    }
    
    /* ------------------------------------------------------------- */
    /* Encoding filters based on a linear-phase filter-bank approach */
    /* ------------------------------------------------------------- */
    else if ( (pData->filterType==FILTER_Z_STYLE) || (pData->filterType==FILTER_Z_STYLE_MAXRE) ) {
        /* Zotter, F. A Linear-Phase Filter-Bank Approach to Process Rigid Spherical Microphone Array Recordings. */
        double normH;
        float f_lim[MAX_SH_ORDER+1];
        double H[HYBRID_BANDS][MAX_SH_ORDER+1];
        double_complex Hs[HYBRID_BANDS][MAX_SH_ORDER+1];
        
        /* find suitable cut-off frequencies */
        ARRAY2SH_WEIGHT_TYPES weightType = arraySpecs->weightType;
        switch (weightType){
            case WEIGHT_OPEN_OMNI:   sphArrayNoiseThreshold(order, arraySpecs->Q, arraySpecs->r, pData->c, ARRAY_CONSTRUCTION_OPEN, 1.0, pData->regPar, f_lim); break;
            case WEIGHT_OPEN_CARD:   sphArrayNoiseThreshold(order, arraySpecs->Q, arraySpecs->r, pData->c, ARRAY_CONSTRUCTION_OPEN_DIRECTIONAL, 0.5, pData->regPar, f_lim); break;
            case WEIGHT_OPEN_DIPOLE: sphArrayNoiseThreshold(order, arraySpecs->Q, arraySpecs->r, pData->c, ARRAY_CONSTRUCTION_OPEN_DIRECTIONAL, 0.0, pData->regPar, f_lim); break;
            case WEIGHT_RIGID_OMNI:
            case WEIGHT_RIGID_CARD:
            case WEIGHT_RIGID_DIPOLE:
                /* Currently no support for estimating the noise cut-off frequencies for rigid scatterers. */
                sphArrayNoiseThreshold(order, arraySpecs->Q, arraySpecs->r, pData->c, ARRAY_CONSTRUCTION_RIGID, 1.0, pData->regPar, f_lim); break;
        }
        
        /* design prototype filterbank */
        for(band=0; band<HYBRID_BANDS; band++){
            normH = 0.0;
            for (n=0; n<order+1; n++){
                if (n==0)
                    H[band][n] = 1.0/(1.0+ pow((double)(pData->freqVector[band]/f_lim[n]),2.0));
                else if (n==order)
                    H[band][n] = pow((double)(pData->freqVector[band]/f_lim[n-1]), (double)order+1.0 )  /
                                 (1.0 + pow((double)(pData->freqVector[band]/f_lim[n-1]), (double)order+1.0));
                else
                    H[band][n] = pow((double)(pData->freqVector[band]/f_lim[n-1]), (double)n+1.0 )  /
                                 (1.0 + pow((double)(pData->freqVector[band]/f_lim[n-1]), (double)n+1.0)) *
                                 (1.0 / (1.0 + pow((double)(pData->freqVector[band]/f_lim[n]), (double)n+2.0)));
                normH += H[band][n];
            }
            /* normalise */
            for (n=0; n<order+1; n++)
                H[band][n] = H[band][n]/normH;
        }
                
        /* compute inverse radial response */ 
        free(pData->bN);
        pData->bN = malloc1d((HYBRID_BANDS)*(order+1)*sizeof(double_complex));
        ARRAY2SH_ARRAY_TYPES arrayType = arraySpecs->arrayType;
        switch(arrayType){
            case ARRAY_CYLINDRICAL:{
                ARRAY2SH_WEIGHT_TYPES weightType = arraySpecs->arrayType;
                switch (weightType){
                    case WEIGHT_RIGID_OMNI:   cylModalCoeffs(order, kr, HYBRID_BANDS, ARRAY_CONSTRUCTION_RIGID, pData->bN); break;
                    case WEIGHT_RIGID_CARD:   /* not supported */ break;
                    case WEIGHT_RIGID_DIPOLE: /* not supported */ break;
                    case WEIGHT_OPEN_OMNI:    cylModalCoeffs(order, kr, HYBRID_BANDS, ARRAY_CONSTRUCTION_OPEN, pData->bN);  break;
                    case WEIGHT_OPEN_CARD:    /* not supported */ break;
                    case WEIGHT_OPEN_DIPOLE:  /* not supported */ break;
                }
                break;
            }
            case ARRAY_SPHERICAL:{
                ARRAY2SH_WEIGHT_TYPES weightType = arraySpecs->arrayType;
                switch (weightType){
                    case WEIGHT_OPEN_OMNI:   sphModalCoeffs(order, kr, HYBRID_BANDS, ARRAY_CONSTRUCTION_OPEN, 1.0, pData->bN); break;
                    case WEIGHT_OPEN_CARD:   sphModalCoeffs(order, kr, HYBRID_BANDS, ARRAY_CONSTRUCTION_OPEN_DIRECTIONAL, 0.5, pData->bN); break;
                    case WEIGHT_OPEN_DIPOLE: sphModalCoeffs(order, kr, HYBRID_BANDS, ARRAY_CONSTRUCTION_OPEN_DIRECTIONAL, 0.0, pData->bN); break;
                    case WEIGHT_RIGID_OMNI:
                    case WEIGHT_RIGID_CARD:
                    case WEIGHT_RIGID_DIPOLE:
                        /* if sensors are flushed with the rigid baffle: */
                        if(arraySpecs->R == arraySpecs->r )
                            sphModalCoeffs(order, kr, HYBRID_BANDS, ARRAY_CONSTRUCTION_RIGID, 1.0, pData->bN);
                        
                        /* if sensors protrude from the rigid baffle: */
                        else{
                            if (arraySpecs->weightType == WEIGHT_RIGID_OMNI)
                                sphScattererModalCoeffs(order, kr, kR, HYBRID_BANDS, pData->bN);
                            else if (arraySpecs->weightType == WEIGHT_RIGID_CARD)
                                sphScattererDirModalCoeffs(order, kr, kR, HYBRID_BANDS, 0.5, pData->bN);
                            else if (arraySpecs->weightType == WEIGHT_RIGID_DIPOLE)
                                sphScattererDirModalCoeffs(order, kr, kR, HYBRID_BANDS, 0.0, pData->bN);
                        }
                        break;
                }
                break;
            }
        }
        
        /* direct inverse (only required for GUI) */
        for(band=0; band<HYBRID_BANDS; band++)
            for(n=0; n < order+1; n++)
                pData->bN_modal[band][n] = ccdiv(cmplx(4.0*SAF_PId, 0.0), pData->bN[band*(order+1)+n]);

        /* phase shift */
        for(band=0; band<HYBRID_BANDS; band++)
            for (n=0; n<order+1; n++)
                Hs[band][n] = ccmul(cexp(cmplx(0.0, kr[band])), ccdiv(cmplx(4.0*SAF_PId, 0.0), pData->bN[band*(order+1)+n]));
        
        /* apply max-re order weighting and diffuse equalisation (not the same as "array2sh_apply_diff_EQ") */
        float* wn;
        double W[MAX_SH_ORDER+1][MAX_SH_ORDER+1];
        double EN, scale;
        int nSH_n;
        memset(W, 0, (MAX_SH_ORDER+1)*(MAX_SH_ORDER+1)*sizeof(double));
        for (n=0; n<order+1; n++){
            nSH_n = (n+1)*(n+1);
            wn = calloc1d(nSH_n*nSH_n, sizeof(float));
            if(pData->filterType==FILTER_Z_STYLE)
                for (i=0; i<n+1; i++)
                    wn[(i*i)*nSH_n+(i*i)] = 1.0f;
            else if(pData->filterType==FILTER_Z_STYLE_MAXRE)
                getMaxREweights(n, 1, wn);
            scale = 0.0;
            for (i=0; i<n+1; i++)
                scale += (double)(2*i+1)*pow((double)wn[(i*i)*nSH_n + (i*i)], 2.0);
            for (i=0; i<n+1; i++)
                W[i][n] = (double)wn[(i*i)*nSH_n + (i*i)]/ sqrt(scale);
            free(wn);
        }
        EN=W[0][n-1];
        for (n=0; n<order+1; n++)
            for (i=0; i<order+1; i++)
                W[i][n] /= EN;
        
        /* apply bandpass filterbank to the inverse array response to regularise it */
        double HW[HYBRID_BANDS];
        double H_np[HYBRID_BANDS][MAX_SH_ORDER+1];
        double W_np[MAX_SH_ORDER+1];
        for (n=0; n<order+1; n++){
            for(band=0; band< HYBRID_BANDS; band++)
                for (i=n, j=0; i<order+1; i++, j++)
                    H_np[band][j] = H[band][i];
            for (i=n, j=0; i<order+1; i++, j++)
                W_np[j] = W[n][i];
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, HYBRID_BANDS, 1, order+1-n, 1.0,
                        (const double*)H_np, MAX_SH_ORDER+1,
                        (const double*)W_np, MAX_SH_ORDER+1, 0.0,
                        (double*)HW, 1);
            for(band=0; band<HYBRID_BANDS; band++)
                pData->bN_inv[band][n] = crmul(Hs[band][n], HW[band]);
        }
        
        /* diag(filters) * Y */
        array2sh_replicate_order(hA2sh, order); /* replicate orders */
        diag_bN_inv_R = calloc1d(nSH*nSH, sizeof(float_complex));
        for(band=0; band<HYBRID_BANDS; band++){
            for(i=0; i<nSH; i++)
                diag_bN_inv_R[i*nSH+i] = cmplxf((float)creal(pData->bN_inv_R[band][i]), (float)cimag(pData->bN_inv_R[band][i])); /* double->single */
            cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nSH, (arraySpecs->Q), nSH, &calpha,
                        diag_bN_inv_R, nSH,
                        pinv_Y_mic_cmplx, nSH, &cbeta,
                        pData->W[band], MAX_NUM_SENSORS);
        }
        free(diag_bN_inv_R);
    }
     
    pData->order = order;
    
    if(pData->enableDiffEQpastAliasing)
        array2sh_apply_diff_EQ(hA2sh);
    
    free(Y_mic);
    free(pinv_Y_mic);
    free(pinv_Y_mic_cmplx);
}

/* Based on a MatLab script by Archontis Politis, 2019 */
void array2sh_apply_diff_EQ(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    array2sh_arrayPars* arraySpecs = (array2sh_arrayPars*)(pData->arraySpecs);
    int i, j, band, array_order, idxf_alias, nSH;
    float f_max, kR_max, f_alias, f_f_alias;
    double_complex* dM_diffcoh_s;
    const double_complex calpha = cmplx(1.0, 0.0); const double_complex cbeta  = cmplx(0.0, 0.0);
    double kr[HYBRID_BANDS];
    double* dM_diffcoh;
    float sensorCoords_rad_local[MAX_NUM_SENSORS][2];
    
    if(arraySpecs->arrayType==ARRAY_CYLINDRICAL)
        return; /* unsupported */
    
    /* prep */
    nSH = (pData->order+1)*(pData->order+1);
    dM_diffcoh = malloc1d((arraySpecs->Q)*(arraySpecs->Q)* (HYBRID_BANDS) * sizeof(double_complex));
    dM_diffcoh_s = malloc1d((arraySpecs->Q)*(arraySpecs->Q) * sizeof(double_complex));
    f_max = 20e3f;
    kR_max = 2.0f*SAF_PI*f_max*(arraySpecs->r)/pData->c;
    array_order = SAF_MIN((int)(ceilf(2.0f*kR_max)+0.01f), 28); /* Cap at around 28, as Bessels at 30+ can be numerically unstable */
    for(band=0; band<HYBRID_BANDS; band++)
        kr[band] = 2.0*SAF_PId*(pData->freqVector[band])*(arraySpecs->r)/pData->c;
    for(i=0; i<arraySpecs->Q; i++){
        sensorCoords_rad_local[i][0] = arraySpecs->sensorCoords_rad[i][0];
        sensorCoords_rad_local[i][1] = arraySpecs->sensorCoords_rad[i][1];
    }
    
    /* Get theoretical diffuse coherence matrix */
    ARRAY2SH_ARRAY_TYPES arrayType = arraySpecs->arrayType;
    switch(arrayType){
        case ARRAY_CYLINDRICAL:
            return; /* Unsupported */
            break;
        case ARRAY_SPHERICAL:{
            ARRAY2SH_WEIGHT_TYPES weightType = arraySpecs->weightType;
            switch (weightType){
                case WEIGHT_RIGID_OMNI: /* Does not handle the case where kr != kR ! */
                    sphDiffCohMtxTheory(array_order, (float*)sensorCoords_rad_local, arraySpecs->Q, ARRAY_CONSTRUCTION_RIGID, 1.0, kr, HYBRID_BANDS, dM_diffcoh);
                    break;
                case WEIGHT_RIGID_CARD:
                    sphDiffCohMtxTheory(array_order, (float*)sensorCoords_rad_local, arraySpecs->Q, ARRAY_CONSTRUCTION_RIGID_DIRECTIONAL, 0.5, kr, HYBRID_BANDS, dM_diffcoh);
                    break;
                case WEIGHT_RIGID_DIPOLE:
                    sphDiffCohMtxTheory(array_order, (float*)sensorCoords_rad_local, arraySpecs->Q, ARRAY_CONSTRUCTION_RIGID_DIRECTIONAL, 0.0, kr, HYBRID_BANDS, dM_diffcoh);
                    break;
                case WEIGHT_OPEN_OMNI:
                    sphDiffCohMtxTheory(array_order, (float*)sensorCoords_rad_local, arraySpecs->Q, ARRAY_CONSTRUCTION_OPEN, 1.0, kr, HYBRID_BANDS, dM_diffcoh);
                    break;
                case WEIGHT_OPEN_CARD:
                    sphDiffCohMtxTheory(array_order, (float*)sensorCoords_rad_local, arraySpecs->Q, ARRAY_CONSTRUCTION_OPEN_DIRECTIONAL, 0.5, kr, HYBRID_BANDS, dM_diffcoh);
                    break;
                case WEIGHT_OPEN_DIPOLE:
                    sphDiffCohMtxTheory(array_order, (float*)sensorCoords_rad_local, arraySpecs->Q, ARRAY_CONSTRUCTION_OPEN_DIRECTIONAL, 0.0, kr, HYBRID_BANDS, dM_diffcoh);
                    break;
            }
            break;
        }
    }

    /* determine band index for the spatial aliasing limit */
    f_alias = sphArrayAliasLim(arraySpecs->r, pData->c, pData->order);
    idxf_alias = 1;
    f_f_alias = 1e13f;
    for(band=0; band<HYBRID_BANDS; band++){
        if( fabsf(pData->freqVector[band]-f_alias) < f_f_alias){
            f_f_alias = fabsf(pData->freqVector[band]-f_alias);
            idxf_alias = band;
        }
    }

    /* baseline */
    for(i=0; i<arraySpecs->Q; i++)
        for(j=0; j<arraySpecs->Q; j++)
            dM_diffcoh_s[i*(arraySpecs->Q)+j] = cmplx(dM_diffcoh[i*(arraySpecs->Q)* (HYBRID_BANDS) + j*(HYBRID_BANDS) + (idxf_alias)], 0.0);
    for(i=0; i<nSH; i++)
        for(j=0; j<arraySpecs->Q; j++)
            pData->W_tmp[i*MAX_NUM_SENSORS+j]= cmplx((double)crealf(pData->W[idxf_alias][i][j]), (double)cimagf(pData->W[idxf_alias][i][j]));
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSH, (arraySpecs->Q), (arraySpecs->Q), &calpha,
                pData->W_tmp, MAX_NUM_SENSORS,
                dM_diffcoh_s, (arraySpecs->Q), &cbeta,
                pData->E_diff, MAX_NUM_SENSORS);
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, nSH, nSH, (arraySpecs->Q), &calpha,
                pData->E_diff, MAX_NUM_SENSORS,
                pData->W_tmp, MAX_NUM_SENSORS, &cbeta,
                pData->L_diff_fal, MAX_NUM_SH_SIGNALS);
    for(i=0; i<nSH; i++)
        pData->L_diff_fal[i*MAX_NUM_SH_SIGNALS+i] = crmul(pData->L_diff_fal[i*MAX_NUM_SH_SIGNALS+i], 1.0/(4.0*SAF_PId)); /* only care about the diagonal entries */

    /* diffuse-field equalise bands above aliasing. */
    for(band = SAF_MAX(idxf_alias,0)+1; band<HYBRID_BANDS; band++){
        for(i=0; i<arraySpecs->Q; i++)
            for(j=0; j<arraySpecs->Q; j++)
                dM_diffcoh_s[i*(arraySpecs->Q)+j] = cmplx(dM_diffcoh[i*(arraySpecs->Q)* (HYBRID_BANDS) + j*(HYBRID_BANDS) + (band)], 0.0);
        for(i=0; i<nSH; i++)
            for(j=0; j<arraySpecs->Q; j++)
                pData->W_tmp[i*MAX_NUM_SENSORS+j]= cmplx((double)crealf(pData->W[band][i][j]), (double)cimagf(pData->W[band][i][j]));
        cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSH, (arraySpecs->Q), (arraySpecs->Q), &calpha,
                    pData->W_tmp, MAX_NUM_SENSORS,
                    dM_diffcoh_s, (arraySpecs->Q), &cbeta,
                    pData->E_diff, MAX_NUM_SENSORS);
        cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, nSH, nSH, (arraySpecs->Q), &calpha,
                    pData->E_diff, MAX_NUM_SENSORS,
                    pData->W_tmp, MAX_NUM_SENSORS, &cbeta,
                    pData->L_diff, MAX_NUM_SH_SIGNALS);
        for(i=0; i<nSH; i++)
            for(j=0; j<nSH; j++)
                pData->L_diff[i*MAX_NUM_SH_SIGNALS+j] = i==j? csqrt(cradd(ccdiv(pData->L_diff_fal[i*MAX_NUM_SH_SIGNALS+j], crmul(pData->L_diff[i*MAX_NUM_SH_SIGNALS+j], 1.0/(4.0*SAF_PId))), 2.23e-10)): cmplx(0.0,0.0);
        cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSH, (arraySpecs->Q), nSH, &calpha,
                    pData->L_diff, MAX_NUM_SH_SIGNALS,
                    pData->W_tmp, MAX_NUM_SENSORS, &cbeta,
                    pData->W_diffEQ_tmp, MAX_NUM_SENSORS);
        for(i=0; i<nSH; i++)
            for(j=0; j<arraySpecs->Q; j++)
                pData->W[band][i][j] = cmplxf((float)creal(pData->W_diffEQ_tmp[i*MAX_NUM_SENSORS+j]), (float)cimag(pData->W_diffEQ_tmp[i*MAX_NUM_SENSORS+j]));
    }
    
    pData->evalStatus = EVAL_STATUS_NOT_EVALUATED;
    
    free(dM_diffcoh);
    free(dM_diffcoh_s);
}

void array2sh_calculate_mag_curves(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    int band, n;
    
    for(band = 0; band <HYBRID_BANDS; band++){
        for(n = 0; n <pData->order+1; n++){
            pData->bN_inv_dB[band][n] = 20.0f * (float)log10(cabs(pData->bN_inv[band][n]));
            pData->bN_modal_dB[band][n] = 20.0f * (float)log10(cabs(pData->bN_modal[band][n]));
        }
    }
}

void array2sh_evaluateSHTfilters(void* hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    array2sh_arrayPars* arraySpecs = (array2sh_arrayPars*)(pData->arraySpecs);
    int band, i, j, simOrder, order, nSH;
    double kr[HYBRID_BANDS];
    double kR[HYBRID_BANDS];
    float* Y_grid_real;
    float_complex* Y_grid, *H_array, *Wshort;
     
    saf_assert(pData->W != NULL, "The initCodec function must have been called prior to calling array2sh_evaluateSHTfilters()");
    
    strcpy(pData->progressBarText,"Simulating microphone array");
    pData->progressBar0_1 = 0.35f;
    
    /* simulate the current array by firing 812 plane-waves around the surface of a theoretical version of the array
     * and ascertaining the transfer function for each */
    simOrder = (int)(2.0f*SAF_PI*MAX_EVAL_FREQ_HZ*(arraySpecs->r)/pData->c)+1;
    for(band=0; band<HYBRID_BANDS; band++){
        kr[band] = 2.0*SAF_PId*(pData->freqVector[band])*(arraySpecs->r)/pData->c;
        kR[band] = 2.0*SAF_PId*(pData->freqVector[band])*(arraySpecs->R)/pData->c;
    }
    H_array = malloc1d((HYBRID_BANDS) * (arraySpecs->Q) * 812*sizeof(float_complex));
    ARRAY2SH_ARRAY_TYPES arrayType = arraySpecs->arrayType;
    switch(arrayType){
        case ARRAY_SPHERICAL:{
            ARRAY2SH_WEIGHT_TYPES weightType = arraySpecs->weightType;
            switch(weightType){
                default:
                case WEIGHT_RIGID_OMNI:
                    simulateSphArray(simOrder, kr, kR, HYBRID_BANDS, (float*)arraySpecs->sensorCoords_rad, arraySpecs->Q,
                                     (float*)__geosphere_ico_9_0_dirs_deg, 812, ARRAY_CONSTRUCTION_RIGID, 1.0, H_array);
                    break;
                case WEIGHT_RIGID_CARD:
                    simulateSphArray(simOrder, kr, kR, HYBRID_BANDS, (float*)arraySpecs->sensorCoords_rad, arraySpecs->Q,
                                     (float*)__geosphere_ico_9_0_dirs_deg, 812, ARRAY_CONSTRUCTION_RIGID_DIRECTIONAL, 0.5, H_array);
                    break;
                case WEIGHT_RIGID_DIPOLE:
                    simulateSphArray(simOrder, kr, kR, HYBRID_BANDS, (float*)arraySpecs->sensorCoords_rad, arraySpecs->Q,
                                     (float*)__geosphere_ico_9_0_dirs_deg, 812, ARRAY_CONSTRUCTION_RIGID_DIRECTIONAL, 0.0, H_array);
                    break;
                case WEIGHT_OPEN_OMNI:
                    simulateSphArray(simOrder, kr, NULL, HYBRID_BANDS, (float*)arraySpecs->sensorCoords_rad, arraySpecs->Q,
                                     (float*)__geosphere_ico_9_0_dirs_deg, 812, ARRAY_CONSTRUCTION_OPEN, 1.0, H_array);
                    break;
                case WEIGHT_OPEN_CARD:
                    simulateSphArray(simOrder, kr, NULL, HYBRID_BANDS, (float*)arraySpecs->sensorCoords_rad, arraySpecs->Q,
                                     (float*)__geosphere_ico_9_0_dirs_deg, 812, ARRAY_CONSTRUCTION_OPEN_DIRECTIONAL, 0.5, H_array);
                    break;
                case WEIGHT_OPEN_DIPOLE:
                    simulateSphArray(simOrder, kr, NULL, HYBRID_BANDS, (float*)arraySpecs->sensorCoords_rad, arraySpecs->Q,
                                     (float*)__geosphere_ico_9_0_dirs_deg, 812, ARRAY_CONSTRUCTION_OPEN_DIRECTIONAL, 0.0, H_array);
                    break;
            }
            break;
        }
        case ARRAY_CYLINDRICAL:{
            ARRAY2SH_WEIGHT_TYPES weightType = arraySpecs->weightType;
            switch(weightType){
                default:
                case WEIGHT_RIGID_OMNI:
                case WEIGHT_RIGID_CARD:
                case WEIGHT_RIGID_DIPOLE:
                    simulateCylArray(simOrder, kr, HYBRID_BANDS, (float*)arraySpecs->sensorCoords_rad, arraySpecs->Q, (float*)__geosphere_ico_9_0_dirs_deg, 812, ARRAY_CONSTRUCTION_RIGID, H_array);
                    break;
                case WEIGHT_OPEN_DIPOLE:
                case WEIGHT_OPEN_CARD:
                case WEIGHT_OPEN_OMNI:
                    simulateCylArray(simOrder, kr, HYBRID_BANDS, (float*)arraySpecs->sensorCoords_rad, arraySpecs->Q, (float*)__geosphere_ico_9_0_dirs_deg, 812, ARRAY_CONSTRUCTION_OPEN, H_array);
                    break;
            }
            break;
        }
    }
    
    strcpy(pData->progressBarText,"Evaluating encoding performance");
    pData->progressBar0_1 = 0.8f;
    
    /* generate ideal (real) spherical harmonics to compare with */
    order = pData->order;
    nSH = (order+1)*(order+1);
    Y_grid_real = malloc1d(nSH*812*sizeof(float));
    getRSH(order, (float*)__geosphere_ico_9_0_dirs_deg, 812, Y_grid_real);
    Y_grid = malloc1d(nSH*812*sizeof(float_complex));
    for(i=0; i<nSH*812; i++)
        Y_grid[i] = cmplxf(Y_grid_real[i], 0.0f); /* "evaluateSHTfilters" function requires complex data type */
    
    /* compare the spherical harmonics obtained from encoding matrix 'W' with the ideal patterns */
    Wshort = malloc1d(HYBRID_BANDS*nSH*(arraySpecs->Q)*sizeof(float_complex));
    for(band=0; band<HYBRID_BANDS; band++)
        for(i=0; i<nSH; i++)
            for(j=0; j<(arraySpecs->Q); j++)
                Wshort[band*nSH*(arraySpecs->Q) + i*(arraySpecs->Q) + j] = pData->W[band][i][j];
    evaluateSHTfilters(order, Wshort, arraySpecs->Q, HYBRID_BANDS, H_array, 812, Y_grid, pData->cSH, pData->lSH);

    free(Y_grid_real);
    free(Y_grid);
    free(H_array);
    free(Wshort);
}

void array2sh_createArray(void ** const hPars)
{
    array2sh_arrayPars* pars = (array2sh_arrayPars*)malloc1d(sizeof(array2sh_arrayPars));
    *hPars = (void*)pars;
}

void array2sh_destroyArray(void ** const hPars)
{
    array2sh_arrayPars *pars = (array2sh_arrayPars*)(*hPars);
    if(pars!=NULL) {
        free(pars);
        pars=NULL;
    }
}
 
void array2sh_initArray
(
    void* const hPars,
    ARRAY2SH_MICROPHONE_ARRAY_PRESETS preset,
    _Atomic_INT32* arrayOrder,
    int firstInitFlag
)
{
    array2sh_arrayPars *pars = (array2sh_arrayPars*)(hPars);
    int ch, i, Q;
    
    switch(preset){
        default:
        case MICROPHONE_ARRAY_PRESET_DEFAULT:
            (*arrayOrder) = 1;
            Q = 4;
            pars->r = 0.02f;
            pars->R = 0.02f;
            pars->arrayType = ARRAY_SPHERICAL;
            pars->weightType = WEIGHT_OPEN_CARD;
            for(ch=0; ch<Q; ch++){
                for(i=0; i<2; i++){
                    pars->sensorCoords_rad[ch][i] = __Sound_field_SPS200_coords_rad[ch][i];
                    pars->sensorCoords_deg[ch][i] = pars->sensorCoords_rad[ch][i] * (180.0f/SAF_PI);
                }
            }
            break;
        case MICROPHONE_ARRAY_PRESET_AALTO_HYDROPHONE:
            (*arrayOrder) = 1;
            Q = 4;
            pars->r = 0.173f;
            pars->R = 0.173f;
            pars->arrayType = ARRAY_SPHERICAL;
            pars->weightType = WEIGHT_OPEN_OMNI;
            for(ch=0; ch<Q; ch++){
                for(i=0; i<2; i++){
                    pars->sensorCoords_rad[ch][i] = __Aalto_Hydrophone_coords_rad[ch][i];
                    pars->sensorCoords_deg[ch][i] = pars->sensorCoords_rad[ch][i] * (180.0f/SAF_PI);
                }
            }
            break;
        case MICROPHONE_ARRAY_PRESET_SENNHEISER_AMBEO:
            (*arrayOrder) = 1;
            Q = 4;
            pars->r = 0.014f;
            pars->R = 0.014f;
            pars->arrayType = ARRAY_SPHERICAL;
            pars->weightType = WEIGHT_OPEN_CARD;
            for(ch=0; ch<Q; ch++){
                for(i=0; i<2; i++){
                    pars->sensorCoords_rad[ch][i] = __Sennheiser_Ambeo_coords_rad[ch][i];
                    pars->sensorCoords_deg[ch][i] = pars->sensorCoords_rad[ch][i] * (180.0f/SAF_PI);
                }
            }
            break;
        case MICROPHONE_ARRAY_PRESET_CORE_SOUND_TETRAMIC:
            (*arrayOrder) = 1;
            Q = 4;
            pars->r = 0.02f;
            pars->R = 0.02f;
            pars->arrayType = ARRAY_SPHERICAL;
            pars->weightType = WEIGHT_OPEN_CARD;
            for(ch=0; ch<Q; ch++){
                for(i=0; i<2; i++){
                    pars->sensorCoords_rad[ch][i] = __Core_Sound_TetraMic_coords_rad[ch][i];
                    pars->sensorCoords_deg[ch][i] = pars->sensorCoords_rad[ch][i] * (180.0f/SAF_PI);
                }
            }
            break;
        case MICROPHONE_ARRAY_PRESET_ZOOM_H3VR_PRESET:
            (*arrayOrder) = 1;
            Q = 4;
            pars->r = 0.012f;
            pars->R = 0.012f;
            pars->arrayType = ARRAY_SPHERICAL;
            pars->weightType = WEIGHT_OPEN_CARD;
            for(ch=0; ch<Q; ch++){
                for(i=0; i<2; i++){
                    pars->sensorCoords_rad[ch][i] = __Zoom_H3VR_coords_rad[ch][i];
                    pars->sensorCoords_deg[ch][i] = pars->sensorCoords_rad[ch][i] * (180.0f/SAF_PI);
                }
            }
            break;
        case MICROPHONE_ARRAY_PRESET_SOUND_FIELD_SPS200:
            (*arrayOrder) = 1;
            Q = 4;
            pars->r = 0.02f;
            pars->R = 0.02f;
            pars->arrayType = ARRAY_SPHERICAL;
            pars->weightType = WEIGHT_OPEN_CARD;
            for(ch=0; ch<Q; ch++){
                for(i=0; i<2; i++){
                    pars->sensorCoords_rad[ch][i] = __Sound_field_SPS200_coords_rad[ch][i];
                    pars->sensorCoords_deg[ch][i] = pars->sensorCoords_rad[ch][i] * (180.0f/SAF_PI);
                }
            }
            break;
        case MICROPHONE_ARRAY_PRESET_ZYLIA_1D:
            (*arrayOrder) = 3;
            Q = 19;
            pars->r = 0.049f;
            pars->R = 0.049f;
            pars->arrayType = ARRAY_SPHERICAL;
            pars->weightType = WEIGHT_RIGID_OMNI;
            for(ch=0; ch<Q; ch++){
                for(i=0; i<2; i++){
                    pars->sensorCoords_rad[ch][i] = __Zylia1D_coords_rad[ch][i];
                    pars->sensorCoords_deg[ch][i] = pars->sensorCoords_rad[ch][i] * (180.0f/SAF_PI);
                }
            }
            break;
        case MICROPHONE_ARRAY_PRESET_EIGENMIKE32:
            (*arrayOrder) = 4;
            Q = 32;
            pars->r = 0.042f;
            pars->R = 0.042f;
            pars->arrayType = ARRAY_SPHERICAL;
            pars->weightType = WEIGHT_RIGID_OMNI;
            for(ch=0; ch<Q; ch++){
                for(i=0; i<2; i++){
                    pars->sensorCoords_rad[ch][i] = __Eigenmike32_coords_rad[ch][i];
                    pars->sensorCoords_deg[ch][i] = pars->sensorCoords_rad[ch][i] * (180.0f/SAF_PI);
                }
            }
            break;
        case MICROPHONE_ARRAY_PRESET_EIGENMIKE64:
            (*arrayOrder) = 6;
            Q = 64;
            pars->r = 0.042f;
            pars->R = 0.042f;
            pars->arrayType = ARRAY_SPHERICAL;
            pars->weightType = WEIGHT_RIGID_OMNI;
            for(ch=0; ch<Q; ch++){
                for(i=0; i<2; i++){
                    pars->sensorCoords_rad[ch][i] = __Eigenmike64_coords_rad[ch][i];
                    pars->sensorCoords_deg[ch][i] = pars->sensorCoords_rad[ch][i] * (180.0f/SAF_PI);
                }
            }
            break;
        case MICROPHONE_ARRAY_PRESET_DTU_MIC:
            (*arrayOrder) = 6;
            Q = 52;
            pars->r = 0.05f;
            pars->R = 0.05f; 
            pars->arrayType = ARRAY_SPHERICAL;
            pars->weightType = WEIGHT_RIGID_OMNI;
            for(ch=0; ch<Q; ch++){
                for(i=0; i<2; i++){
                    pars->sensorCoords_rad[ch][i] = __DTU_mic_coords_rad[ch][i];
                    pars->sensorCoords_deg[ch][i] = pars->sensorCoords_rad[ch][i] * (180.0f/SAF_PI);
                }
            }
            break;
    }
    
    /* Fill remaining slots with default coords */
    for(; ch<MAX_NUM_SENSORS_IN_PRESET; ch++){
        for(i=0; i<2; i++){
            pars->sensorCoords_deg[ch][i] = __default_SENSORcoords128_deg[ch][i];
            pars->sensorCoords_rad[ch][i] = pars->sensorCoords_deg[ch][i] * (SAF_PI/180.0f);
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

