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

//#if defined(__APPLE__) || defined(NDEBUG)
  #define SAF_ENABLE_SOFA_READER
  #include "saf_sofa_reader.h"
//#endif

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
    
    /* Convert HRIRs into Filterbank coefficients */ 
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
    if(pars->hrtf_fb!= NULL){
        free(pars->hrtf_fb);
        pars->hrtf_fb = NULL;
    }
    HRIRs2FilterbankHRTFs(pars->hrirs, pars->N_hrir_dirs, pars->hrir_len, pars->itds_s, (float*)pData->freqVector, HYBRID_BANDS, &(pars->hrtf_fb));
    
    /* calculate binaural ambisonic decoding matrix */
    t = 2*(pData->order+1);
    Y_td = NULL;
    nDirs_td = __Tdesign_nPoints_per_degree[t-1];
    t_dirs = (float*)__HANDLES_Tdesign_dirs_deg[t-1];
    getRSH(pData->order, t_dirs, nDirs_td, &Y_td);
    M_dec_t = malloc(nSH*nDirs_td*sizeof(float_complex));
    scale = 1.0f/(float)nDirs_td;
    for(i=0; i<nDirs_td*nSH; i++)
        M_dec_t[i] = cmplxf(Y_td[i] * scale, 0.0f);
    hrir_closest_idx = malloc(nDirs_td*sizeof(int));
    hrtf_fb_short = malloc(2*nDirs_td*sizeof(float_complex));
    memset(pars->M_dec, 0, HYBRID_BANDS*NUM_EARS*MAX_NUM_SH_SIGNALS*sizeof(float_complex));
    findClosestGridPoints(pars->hrir_dirs_deg, pars->N_hrir_dirs, t_dirs, nDirs_td, 1, hrir_closest_idx, NULL, NULL);
    for(band=0; band<HYBRID_BANDS; band++){
        for(i=0; i<nDirs_td; i++){
            hrtf_fb_short[0*nDirs_td+i] = pars->hrtf_fb[band*2*(pars->N_hrir_dirs)+0*(pars->N_hrir_dirs) + hrir_closest_idx[i]];
            hrtf_fb_short[1*nDirs_td+i] = pars->hrtf_fb[band*2*(pars->N_hrir_dirs)+1*(pars->N_hrir_dirs) + hrir_closest_idx[i]];
        }
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasTrans, NUM_EARS, nSH, nDirs_td, &calpha,
                    hrtf_fb_short, nDirs_td,
                    M_dec_t, nDirs_td, &cbeta,
                    pars->M_dec[band], MAX_NUM_SH_SIGNALS);
    }
    
    int dfdfdfdfdf[100] = {0};
    memcpy(dfdfdfdfdf,hrir_closest_idx, nDirs_td*sizeof(int));
    
    free(M_dec_t);
    free(hrir_closest_idx);
    free(hrtf_fb_short);
    free(Y_td);
}

void ambi_bin_initTFT
(
    void* const hAmbi
)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    int t, ch;

    /* free afSTFT + buffers, if already allocated */
    if (pData->hSTFT != NULL){
        afSTFTfree(pData->hSTFT);
        pData->hSTFT = NULL;
        for (t = 0; t<TIME_SLOTS; t++) {
            for (ch = 0; ch< pData->nSH; ch++) {
                free(pData->STFTInputFrameTF[t][ch].re);
                free(pData->STFTInputFrameTF[t][ch].im);
            }
        }
        free2d((void**)pData->tempHopFrameTD, MAX(pData->nSH, NUM_EARS));
        free2d((void**)pData->STFTInputFrameTF, TIME_SLOTS);
    }
    
    /* reallocate afSTFT + buffers */
    if (pData->hSTFT == NULL){
        afSTFTinit(&(pData->hSTFT), HOP_SIZE, pData->new_nSH, NUM_EARS, 0, 1);
        pData->STFTInputFrameTF = (complexVector**)malloc2d(TIME_SLOTS, pData->new_nSH, sizeof(complexVector));
        for(t=0; t<TIME_SLOTS; t++) {
            for(ch=0; ch< pData->new_nSH; ch++) {
                pData->STFTInputFrameTF[t][ch].re = (float*)calloc(HYBRID_BANDS, sizeof(float));
                pData->STFTInputFrameTF[t][ch].im = (float*)calloc(HYBRID_BANDS, sizeof(float));
            }
        }
        pData->tempHopFrameTD = (float**)malloc2d( MAX(pData->new_nSH, NUM_EARS), HOP_SIZE, sizeof(float));
        pData->nSH = pData->new_nSH;
    }
}






