/*
 Copyright 2018 Leo McCormack
 
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
 *     ambi_bin_internal.c
 * Description:
 *     A binaural Ambisonic decoder for reproducing ambisonic signals over headphones.
 *     Optionally, a SOFA file may be loaded for personalised headphone listening.
 * Dependencies:
 *     saf_utilities, afSTFTlib, saf_hrir, saf_sh
 * Author, date created:
 *     Leo McCormack, 14.04.2018
 */

#include "ambi_bin_internal.h"
 
#define SAF_ENABLE_SOFA_READER
#include "saf_sofa_reader.h" 

void ambi_bin_initCodec
(
    void* const hAmbi
)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    codecPars* pars = pData->pars;
    int i, j, k, nSH, band, t, nDirs_td;
    int* hrir_closest_idx;
    float scale;
    float* Y_td, *t_dirs;
    float_complex* M_dec_t, *hrtf_fb_short;
    const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);
    
    nSH = (pData->order+1)*(pData->order+1);
    
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
        pars->hrirs = malloc(pars->N_hrir_dirs * NUM_EARS * (pars->hrir_len)*sizeof(float));
        for(i=0; i<pars->N_hrir_dirs; i++)
            for(j=0; j<NUM_EARS; j++)
                for(k=0; k< pars->hrir_len; k++)
                    pars->hrirs[i*NUM_EARS*(pars->hrir_len) + j*(pars->hrir_len) + k] = (float)__default_hrirs[i][j][k];
        if(pars->hrir_dirs_deg != NULL)
            free(pars->hrir_dirs_deg);
        pars->hrir_dirs_deg = malloc(pars->N_hrir_dirs * NUM_EARS * sizeof(float));
        for(i=0; i<pars->N_hrir_dirs; i++)
            for(j=0; j<2; j++)
                pars->hrir_dirs_deg[i*2+j] = (float)__default_hrir_dirs_deg[i][j];
    }
    
    /* estimate the ITDs for each HRIR */
    if(pars->itds_s!= NULL){
        free(pars->itds_s);
        pars->itds_s = NULL;
    }
    estimateITDs(pars->hrirs, pars->N_hrir_dirs, pars->hrir_len, pars->hrir_fs, &(pars->itds_s));
    
    /* convert hrirs to filterbank coefficients */
    for (i=0; i<2; i++){ /* [0] WITHOUT, [1] WITH phase manipulation */
        if(pars->hrtf_fb[i]!= NULL){
            free(pars->hrtf_fb[i]);
            pars->hrtf_fb[i] = NULL; 
        }
#if USE_NEAREST_HRIRS
        HRIRs2FilterbankHRTFs(pars->hrirs, pars->N_hrir_dirs, pars->hrir_len, pars->itds_s, (float*)pData->freqVector, HYBRID_BANDS, i, &(pars->hrtf_fb[i]));
#else
        HRIRs2FilterbankHRTFs(pars->hrirs, pars->N_hrir_dirs, pars->hrir_len, pars->itds_s, (float*)pData->freqVector, HYBRID_BANDS, 0, &(pars->hrtf_fb[i]));
#endif
    }
    
    /* define t-design for this order */
    if(pData->order == 0)
        t = 2;
    else if(pData->order <= 1)
        t = 4*(pData->order);
    else
        t = 2*(pData->order);
    if (t<=20){
        nDirs_td = __Tdesign_nPoints_per_degree[t-1];
        t_dirs = (float*)__HANDLES_Tdesign_dirs_deg[t-1];
    }
    else {
        t_dirs = (float*)__Tdesign_degree_30_dirs_deg;
        nDirs_td = 480;
    }
    
    /* define M_dec_td (decoder to t-design) */
    Y_td = NULL;
    getRSH(pData->order, t_dirs, nDirs_td, &Y_td);
    float* Y_td_mrE,* a_n;
    M_dec_t = malloc(nSH*nDirs_td*sizeof(float_complex));
    if(pData->rE_WEIGHT){
        Y_td_mrE  = malloc(nSH*nDirs_td*sizeof(float));
        a_n = malloc(nSH*nSH*sizeof(float));
        getMaxREweights(pData->order, a_n);
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSH, nDirs_td, nSH, 1.0f,
                    a_n, nSH,
                    Y_td, nDirs_td, 0.0f,
                    Y_td_mrE, nDirs_td);
        memcpy(Y_td, Y_td_mrE, nSH*nDirs_td*sizeof(float));
        free(Y_td_mrE);
        free(a_n);
    }
    scale = 1.0f/(float)nDirs_td;
    for(i=0; i<nDirs_td*nSH; i++)
        M_dec_t[i] = cmplxf(Y_td[i] * scale, 0.0f);
    
#if USE_NEAREST_HRIRS
    /* M_dec(f) = H(f) * M_dec_td */
    hrir_closest_idx = malloc(nDirs_td*sizeof(int));
    hrtf_fb_short = malloc(NUM_EARS*nDirs_td*sizeof(float_complex));
    findClosestGridPoints(pars->hrir_dirs_deg, pars->N_hrir_dirs, t_dirs, nDirs_td, 1, hrir_closest_idx, NULL, NULL);
    memset(pars->M_dec, 0, 2*HYBRID_BANDS*NUM_EARS*MAX_NUM_SH_SIGNALS*sizeof(float_complex));
    for(j=0; j<2; j++){ /* [0] WITHOUT, [1] WITH phase manipulation */
        for(band=0; band<HYBRID_BANDS; band++){
            for(i=0; i<nDirs_td; i++){
                hrtf_fb_short[0*nDirs_td+i] = pars->hrtf_fb[j][band*NUM_EARS*(pars->N_hrir_dirs)+0*(pars->N_hrir_dirs) + hrir_closest_idx[i]];
                hrtf_fb_short[1*nDirs_td+i] = pars->hrtf_fb[j][band*NUM_EARS*(pars->N_hrir_dirs)+1*(pars->N_hrir_dirs) + hrir_closest_idx[i]];
            }
            cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasTrans, NUM_EARS, nSH, nDirs_td, &calpha,
                        hrtf_fb_short, nDirs_td,
                        M_dec_t, nDirs_td, &cbeta,
                        pars->M_dec[j][band], MAX_NUM_SH_SIGNALS);
        }
    }
    free(hrir_closest_idx);
    free(hrtf_fb_short);
    
    /* Alternatively, use triangular interpolation: */
#else
    /* calculate interpolation table */
    int N_gtable,nTriangles;
    float* gtable;
    gtable = NULL;
    generateVBAPgainTable3D_srcs(t_dirs, nDirs_td, pars->hrir_dirs_deg, pars->N_hrir_dirs, 0, 1, 0.0f, &gtable, &N_gtable, &nTriangles);
    VBAPgainTable2InterpTable(gtable, nDirs_td, pars->N_hrir_dirs);
    
    /* M_dec(f) = H(f) * M_dec_td */
    hrtf_fb_short = malloc(HYBRID_BANDS*NUM_EARS*nDirs_td*sizeof(float_complex));
    memset(pars->M_dec, 0, 2*HYBRID_BANDS*NUM_EARS*MAX_NUM_SH_SIGNALS*sizeof(float_complex));
    for(j=0; j<2; j++){ /* [0] WITHOUT, [1] WITH phase manipulation */
        interpFilterbankHRTFs(pars->hrtf_fb[j], pars->itds_s, pData->freqVector, gtable, pars->N_hrir_dirs, HYBRID_BANDS, nDirs_td, j, hrtf_fb_short);

        for(band=0; band<HYBRID_BANDS; band++){
            cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasTrans, NUM_EARS, nSH, nDirs_td, &calpha,
                        &hrtf_fb_short[band*NUM_EARS*nDirs_td], nDirs_td,
                        M_dec_t, nDirs_td, &cbeta,
                        pars->M_dec[j][band], MAX_NUM_SH_SIGNALS);
        }
    }
    free(gtable);
#endif
    
    free(M_dec_t);
    free(Y_td);
}

void ambi_bin_initTFT
(
    void* const hAmbi
)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);

	/* Initialise afSTFT */
    if(pData->hSTFT==NULL)
        afSTFTinit(&(pData->hSTFT), HOP_SIZE, pData->new_nSH, NUM_EARS, 0, 1);
    else /* Or change the number of channels */
        afSTFTchannelChange(pData->hSTFT, pData->new_nSH, NUM_EARS);
    pData->nSH = pData->new_nSH;
}






