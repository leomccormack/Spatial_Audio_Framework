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
 *     ambi_dec_internal.c
 * Description:
 *     A frequency-dependent Ambisonic decoder for loudspeakers or headphones. Different
 *     decoder settings can be specified for the low and high frequencies. When utilising
 *     spherical harmonic signals derived from real microphone arrays, this implementation
 *     also allows the decoding order per frequency band to be specified. Optionally, a SOFA
 *     file may be loaded for personalised headphone listening.
 *     The algorithms utilised in this Ambisonic decoder were pieced together and developed
 *     in collaboration with Archontis Politis.
 * Dependencies:
 *     saf_utilities, afSTFTlib, saf_hoa, saf_vbap, saf_hrir, saf_sh
 * Author, date created:
 *     Leo McCormack, 07.12.2017
 */

#include "ambi_dec_internal.h"
#ifdef __APPLE__ /* lost my patience trying to statically build netcdf for Windows... */
  #define SAF_ENABLE_SOFA_READER
  #include "saf_sofa_reader.h"
#endif

static inline float matlab_fmodf(float x, float y) {
    float tmp = fmodf(x, y);
    return tmp >= 0 ? tmp : tmp + y;
}

void ambi_dec_initCodec
(
    void* const hAmbi
)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    codecPars* pars = pData->pars;
    int i, d, j, n, ng, nGrid_dirs, nSH_order;
    float* grid_dirs_deg, *Y, *M_dec_tmp, *g, *a, *e, *a_n;
    float a_avg[SH_ORDER], e_avg[SH_ORDER];
    
    M_dec_tmp = NULL;
    nGrid_dirs = 480; /* Minimum t-design of degree 30, has 480 points */
    g = malloc(pData->nLoudpkrs*sizeof(float));
    a = malloc(nGrid_dirs*sizeof(float));
    e = malloc(nGrid_dirs*sizeof(float));
    
    /* calculate loudspeaker decoding matrices */
    for( d=0; d<NUM_DECODERS; d++){
        getAmbiDecoder((float*)pData->loudpkrs_dirs_deg, pData->nLoudpkrs, pData->dec_method[d], SH_ORDER, &(M_dec_tmp));
        
        /* diffuse-field EQ for orders 1..SH_ORDER */
        for( n=1; n<=SH_ORDER; n++){
            /* truncate M_dec for each order */
            nSH_order = (n+1)*(n+1);
            free(pars->M_dec[d][n-1]); 
            pars->M_dec[d][n-1] = malloc(pData->nLoudpkrs * nSH_order * sizeof(float));
            free(pars->M_dec_cmplx[d][n-1]);
            pars->M_dec_cmplx[d][n-1] = malloc(pData->nLoudpkrs * nSH_order * sizeof(float_complex));
            for(i=0; i<pData->nLoudpkrs; i++){
                for(j=0; j<nSH_order; j++){
                    pars->M_dec[d][n-1][i*nSH_order+j] = M_dec_tmp[i*MAX_NUM_SH_SIGNALS +j]; /* for applying in the time domain, and... */
                    pars->M_dec_cmplx[d][n-1][i*nSH_order+j] = cmplxf(pars->M_dec[d][n-1][i*nSH_order+j], 0.0f); /* for the time-frequency domain */
                }
            }
            
            /* create dedicated maxrE weighted versions too (c'mon.. RAM is cheap nowadays) */
            a_n = malloc(nSH_order*nSH_order*sizeof(float));
            getMaxREweights(n, a_n); /* weights returned as diagonal matrix */
            free(pars->M_dec_maxrE[d][n-1]);
            pars->M_dec_maxrE[d][n-1] = malloc(pData->nLoudpkrs * nSH_order * sizeof(float));
            free(pars->M_dec_cmplx_maxrE[d][n-1]);
            pars->M_dec_cmplx_maxrE[d][n-1] = malloc(pData->nLoudpkrs * nSH_order * sizeof(float_complex));
            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, pData->nLoudpkrs, nSH_order, nSH_order, 1.0f,
                        pars->M_dec[d][n-1], nSH_order,
                        a_n, nSH_order, 0.0f,
                        pars->M_dec_maxrE[d][n-1], nSH_order); /* for applying in the time domain */
            for(i=0; i<pData->nLoudpkrs * nSH_order; i++)
                pars->M_dec_cmplx_maxrE[d][n-1][i] = cmplxf(pars->M_dec_maxrE[d][n-1][i], 0.0f); /* for the time-frequency domain */
            
            /* fire a plane-wave from each grid direction to find the total energy/amplitude (using non-maxrE weighted versions) */
            Y = malloc(nSH_order*sizeof(float));
            grid_dirs_deg = (float*)(&__Tdesign_degree_30_dirs_deg[0][0]);
            for(ng=0; ng<nGrid_dirs; ng++){
                getSHreal(n, grid_dirs_deg[ng*2]*M_PI/180.0f, M_PI/2.0f-grid_dirs_deg[ng*2+1]*M_PI/180.0f, Y);
                cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, pData->nLoudpkrs, 1, nSH_order, 1.0f,
                            pars->M_dec[d][n-1], nSH_order,
                            Y, nSH_order, 0.0f,
                            g, 1);
                a[ng] = e[ng] = 0.0f;
                for(i=0; i<pData->nLoudpkrs; i++){
                    a[ng] += g[i];
                    e[ng] += powf(g[i], 2.0f);
                }
            }
            
            /* determine the order+decoder dependent normalisation factor for energy+amplitude preserving decoding */
            a_avg[n-1] = e_avg[n-1] = 0.0f;
            for(ng=0; ng<nGrid_dirs; ng++){
                a_avg[n-1] += a[ng];
                e_avg[n-1] += e[ng];
            }
            a_avg[n-1] /= (float)nGrid_dirs;
            e_avg[n-1] /= (float)nGrid_dirs;
            pars->M_norm[d][n-1][0] = 1.0f/(a_avg[n-1]+2.23e-6f); /* use this to preserve omni amplitude */
            pars->M_norm[d][n-1][1] = sqrtf(1.0f/(e_avg[n-1]+2.23e-6f));  /* use this to preserve omni energy */
            free(a_n);
            free(Y);
        }
        free(M_dec_tmp);
        M_dec_tmp = NULL;
    }
    free(g);
    free(a);
    free(e);
}

void ambi_dec_initHRTFs
(
    void* const hAmbi
)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    codecPars* pars = pData->pars;
    int i, j, k;
    float* hrtf_vbap_gtable;
    
    /* load sofa file or load default hrir data */
    if(!pData->useDefaultHRIRsFLAG && pars->sofa_filepath!=NULL){
#ifdef SAF_ENABLE_SOFA_READER
        loadSofaFile(pars->sofa_filepath,
                     &pars->hrirs,
                     &pars->hrir_dirs_deg,
                     &pars->N_hrir_dirs,
                     &pars->hrir_len,
                     &pars->hrir_fs);
        if((pars->hrirs==NULL) || (pars->hrir_dirs_deg == NULL))
            pData->useDefaultHRIRsFLAG = 1;
#else
        pData->useDefaultHRIRsFLAG = 1;
#endif
    }
    if(pData->useDefaultHRIRsFLAG){
        /* load defaults */
        pars->N_hrir_dirs = __default_N_hrir_dirs;
        pars->hrir_len = __default_hrir_len;
        pars->hrir_fs = __default_hrir_fs;
        if(pars->hrirs != NULL)
            free(pars->hrirs);
        pars->hrirs = malloc(pars->N_hrir_dirs * 2 * (pars->hrir_len)*sizeof(float));
        for(i=0; i<pars->N_hrir_dirs; i++)
            for(j=0; j<2; j++)
                for(k=0; k< pars->hrir_len; k++)
                    pars->hrirs[i*2*(pars->hrir_len) + j*(pars->hrir_len) + k] = (float)__default_hrirs[i][j][k];
        if(pars->hrir_dirs_deg != NULL)
            free(pars->hrir_dirs_deg);
        pars->hrir_dirs_deg = malloc(pars->N_hrir_dirs * 2 * sizeof(float));
        for(i=0; i<pars->N_hrir_dirs; i++)
            for(j=0; j<2; j++)
                pars->hrir_dirs_deg[i*2+j] = (float)__default_hrir_dirs_deg[i][j];
    }
    
    /* estimate the ITDs for each HRIR */
    if(pars->itds_s!= NULL){
        free(pars->itds_s);
        pars->itds_s = NULL;
    }
    hrirlib_estimateITDs(pars->hrirs, pars->N_hrir_dirs, pars->hrir_len, pars->hrir_fs, &(pars->itds_s));
    
    /* generate VBAP gain table for the hrir_dirs */
    hrtf_vbap_gtable = NULL;
    pars->hrtf_vbapTableRes[0] = 2; /* azimuth resolution in degrees */
    pars->hrtf_vbapTableRes[1] = 5; /* elevation resolution in degrees */
    generateVBAPgainTable3D(pars->hrir_dirs_deg, pars->N_hrir_dirs, pars->hrtf_vbapTableRes[0], pars->hrtf_vbapTableRes[1], 1, 0,
                            &hrtf_vbap_gtable, &(pars->N_hrtf_vbap_gtable), &(pars->hrtf_nTriangles));
    if(hrtf_vbap_gtable==NULL){
        /* if generating vbap gain tabled failed, re-calculate with default HRIR set (which is known to triangulate correctly) */
        pData->useDefaultHRIRsFLAG = 1;
        ambi_dec_initHRTFs(hAmbi);
    }
    
    /* compress VBAP table (i.e. remove the zero elements) */
    if(pars->hrtf_vbap_gtableComp!= NULL){
        free(pars->hrtf_vbap_gtableComp);
        pars->hrtf_vbap_gtableComp = NULL;
    }
    if(pars->hrtf_vbap_gtableIdx!= NULL){
        free(pars->hrtf_vbap_gtableIdx);
        pars->hrtf_vbap_gtableIdx = NULL;
    }
    compressVBAPgainTable3D(hrtf_vbap_gtable, pars->N_hrtf_vbap_gtable, pars->N_hrir_dirs, &(pars->hrtf_vbap_gtableComp), &(pars->hrtf_vbap_gtableIdx));
    
    /* convert hrirs to filterbank coefficients */
    if(pars->hrtf_fb!= NULL){
        free(pars->hrtf_fb);
        pars->hrtf_fb = NULL;
    }
    hrirlib_HRIRs2FilterbankHRTFs(pars->hrirs, pars->N_hrir_dirs, pars->hrir_len, pars->itds_s, (float*)pData->freqVector, HYBRID_BANDS, &(pars->hrtf_fb));
    
    /* calculate magnitude responses */
    if(pars->hrtf_fb_mag!= NULL)
        free(pars->hrtf_fb_mag);
    pars->hrtf_fb_mag = malloc(HYBRID_BANDS*NUM_EARS* (pars->N_hrir_dirs)*sizeof(float));
    for(i=0; i<HYBRID_BANDS*NUM_EARS* (pars->N_hrir_dirs); i++)
        pars->hrtf_fb_mag[i] = cabsf(pars->hrtf_fb[i]);
    
    /* clean-up */
    if(hrtf_vbap_gtable!=NULL)
        free(hrtf_vbap_gtable);
}

void ambi_dec_initTFT
(
    void* const hAmbi
)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    int t, ch;

    /* free afSTFT + buffers, if already allocated */
    if (pData->hSTFT != NULL){
        afSTFTfree(pData->hSTFT);
        pData->hSTFT = NULL;
        if(pData->binauraliseLS){
            for (t = 0; t<TIME_SLOTS; t++) {
                for (ch = 0; ch< NUM_EARS; ch++) {
                    free(pData->STFTOutputFrameTF[t][ch].re);
                    free(pData->STFTOutputFrameTF[t][ch].im);
                }
            }
            free2d((void**)pData->tempHopFrameTD, MAX(MAX_NUM_SH_SIGNALS, NUM_EARS));
        }
        else{
            for (t = 0; t<TIME_SLOTS; t++) {
                for (ch = 0; ch< pData->nLoudpkrs; ch++) {
                    free(pData->STFTOutputFrameTF[t][ch].re);
                    free(pData->STFTOutputFrameTF[t][ch].im);
                }
            }
            free2d((void**)pData->tempHopFrameTD, MAX(MAX_NUM_SH_SIGNALS, pData->nLoudpkrs));
        }
        free2d((void**)pData->STFTOutputFrameTF, TIME_SLOTS);
    }
    
    /* reallocate afSTFT + buffers */
    if (pData->hSTFT == NULL){
        if(pData->new_binauraliseLS){
            afSTFTinit(&(pData->hSTFT), HOP_SIZE, MAX_NUM_SH_SIGNALS, NUM_EARS, 0, 1);
            pData->STFTOutputFrameTF = (complexVector**)malloc2d(TIME_SLOTS, NUM_EARS, sizeof(complexVector));
            for(t=0; t<TIME_SLOTS; t++) {
                for(ch=0; ch< NUM_EARS; ch++) {
                    pData->STFTOutputFrameTF[t][ch].re = (float*)calloc(HYBRID_BANDS, sizeof(float));
                    pData->STFTOutputFrameTF[t][ch].im = (float*)calloc(HYBRID_BANDS, sizeof(float));
                }
            }
            pData->tempHopFrameTD = (float**)malloc2d( MAX(MAX_NUM_SH_SIGNALS, NUM_EARS), HOP_SIZE, sizeof(float));
            pData->nLoudpkrs = pData->new_nLoudpkrs;
        }
        else{
            afSTFTinit(&(pData->hSTFT), HOP_SIZE, MAX_NUM_SH_SIGNALS, pData->new_nLoudpkrs, 0, 1);
            pData->STFTOutputFrameTF = (complexVector**)malloc2d(TIME_SLOTS, pData->new_nLoudpkrs, sizeof(complexVector));
            for(t=0; t<TIME_SLOTS; t++) {
                for(ch=0; ch< pData->new_nLoudpkrs; ch++) {
                    pData->STFTOutputFrameTF[t][ch].re = (float*)calloc(HYBRID_BANDS, sizeof(float));
                    pData->STFTOutputFrameTF[t][ch].im = (float*)calloc(HYBRID_BANDS, sizeof(float));
                }
            }
            pData->tempHopFrameTD = (float**)malloc2d( MAX(MAX_NUM_SH_SIGNALS, pData->new_nLoudpkrs), HOP_SIZE, sizeof(float));
            pData->nLoudpkrs = pData->new_nLoudpkrs;
        }
        pData->binauraliseLS = pData->new_binauraliseLS;
    }
}

void ambi_dec_interpHRTFs
(
    void* const hAmbi,
    float azimuth_deg,
    float elevation_deg,
    float_complex h_intrp[HYBRID_BANDS][NUM_EARS]
)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    codecPars* pars = pData->pars;
    int i, band;
    int aziIndex, elevIndex, N_azi, idx3d;
    float_complex ipd;
    float aziRes, elevRes, weights[1][3], itds3[3],  itdInterp[1];
    float magnitudes3[HYBRID_BANDS][3][NUM_EARS], magInterp[HYBRID_BANDS][NUM_EARS];

    /* find closest pre-computed VBAP direction */
    aziRes = (float)pars->hrtf_vbapTableRes[0];
    elevRes = (float)pars->hrtf_vbapTableRes[1];
    N_azi = (int)(360.0f / aziRes + 0.5f) + 1;
    aziIndex = (int)(matlab_fmodf(azimuth_deg + 180.0f, 360.0f) / aziRes + 0.5f);
    elevIndex = (int)((elevation_deg + 90.0f) / elevRes + 0.5f);
    idx3d = elevIndex * N_azi + aziIndex;
    for (i = 0; i < 3; i++)
        weights[0][i] = pars->hrtf_vbap_gtableComp[idx3d*3 + i];
    
    /* retrieve the 3 itds and hrtf magnitudes */
    for (i = 0; i < 3; i++) {
        itds3[i] = pars->itds_s[pars->hrtf_vbap_gtableIdx[idx3d*3+i]];
        for (band = 0; band < HYBRID_BANDS; band++) {
            magnitudes3[band][i][0] = pars->hrtf_fb_mag[band*NUM_EARS*(pars->N_hrir_dirs) + 0*(pars->N_hrir_dirs) + pars->hrtf_vbap_gtableIdx[idx3d*3+i]];
            magnitudes3[band][i][1] = pars->hrtf_fb_mag[band*NUM_EARS*(pars->N_hrir_dirs) + 1*(pars->N_hrir_dirs) + pars->hrtf_vbap_gtableIdx[idx3d*3+i]];
        }
    }
    
    /* interpolate hrtf magnitudes and itd seperately */
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 1, 1, 3, 1,
                (float*)weights, 3,
                (float*)itds3, 1, 0,
                (float*)itdInterp, 1);
    for (band = 0; band < HYBRID_BANDS; band++) {
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 1, 2, 3, 1,
                    (float*)weights, 3,
                    (float*)magnitudes3[band], 2, 0,
                    (float*)magInterp[band], 2);
    }
    
    /* reintroduce the interaural phase difference per band */
    for (band = 0; band < HYBRID_BANDS; band++) {
        ipd = cmplxf(0.0f, (matlab_fmodf(2.0f*PI* (pData->freqVector[band]) * itdInterp[0] + PI, 2.0f*PI) - PI) / 2.0f);
        h_intrp[band][0] = ccmulf(cmplxf(magInterp[band][0], 0.0f), cexpf(ipd));
        h_intrp[band][1] = ccmulf(cmplxf(magInterp[band][1], 0.0f), conjf(cexpf(ipd)));
    }
}

void ambi_dec_loadPreset(PRESETS preset, float dirs_deg[MAX_NUM_LOUDSPEAKERS][2], int* newNCH, int* nDims)
{
    float sum_elev;
    int ch, i, nCH;
    
    switch(preset){
        default:
        case PRESET_DEFAULT: 
#ifdef ENABLE_5PX_PRESET
        case PRESET_5PX:
            nCH = 5;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __5pX_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_7PX_PRESET
        case PRESET_7PX:
            nCH = 7;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __7pX_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_8PX_PRESET
        case PRESET_8PX:
            nCH = 8;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __8pX_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_9PX_PRESET
        case PRESET_9PX:
            nCH = 9;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __9pX_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_10PX_PRESET
        case PRESET_10PX:
            nCH = 10;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __10pX_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_11PX_PRESET
        case PRESET_11PX:
            nCH = 11;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __11pX_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_11PX_7_4_PRESET
        case PRESET_11PX_7_4:
            nCH = 11;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __11pX_7_4_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_13PX_PRESET
        case PRESET_13PX:
            nCH = 13;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __13pX_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_22PX_PRESET
        case PRESET_22PX:
            nCH = 22;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __22pX_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_AALTO_MCC_PRESET
        case PRESET_AALTO_MCC:
            nCH = 29;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Aalto_MCC_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_AALTO_APAJA_PRESET
        case PRESET_AALTO_APAJA:
            nCH = 29;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Aalto_Apaja_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_AALTO_APAJA2_PRESET
        case PRESET_AALTO_APAJA2:
            nCH = 39;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Aalto_Apaja2_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_AALTO_LR_PRESET
        case PRESET_AALTO_LR:
            nCH = 13;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Aalto_LR_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_DTU_AVIL_PRESET
        case PRESET_DTU_AVIL:
            nCH = 64;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __DTU_AVIL_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_T_DESIGN_4_PRESET
        case PRESET_T_DESIGN_4:
            nCH = 4;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_2_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_T_DESIGN_12_PRESET
        case PRESET_T_DESIGN_12:
            nCH = 12;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_4_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_T_DESIGN_24_PRESET
        case PRESET_T_DESIGN_24:
            nCH = 24;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_6_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_T_DESIGN_36_PRESET
        case PRESET_T_DESIGN_36:
            nCH = 36;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_8_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_T_DESIGN_48_PRESET
        case PRESET_T_DESIGN_48:
            nCH = 48;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_9_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_T_DESIGN_60_PRESET
        case PRESET_T_DESIGN_60:
            nCH = 60;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_10_dirs_deg[ch][i];
            break;
#endif
    }
    
    /* Fill remaining slots with default coords */
    for(; ch<MAX_NUM_LOUDSPEAKERS; ch++)
        for(i=0; i<2; i++)
            dirs_deg[ch][i] = default_LScoords64_rad[ch][i]* (180.0f/M_PI);
    
    /* specify new number of channels (for dynamically changing the number of TFT channels) */
    (*newNCH) = nCH;
    
    /* estimate number of dimensions. (Obviously fails if using 2D setups that are elevated.
     However, in these cases, triangulation should fail and revert to 2D anyway) */
    sum_elev = 0.0f;
    for(i=0; i<nCH; i++)
        sum_elev += fabsf(dirs_deg[i][1]);
    if(sum_elev < 0.01f)
        (*nDims) = 2;
    else
        (*nDims) = 3;
}
 







