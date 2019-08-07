/*
 * Copyright 2018 Leo McCormack
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

/*
 * Filename: ambi_bin_internal.c
 * -----------------------------
 * A binaural Ambisonic decoder for reproducing ambisonic signals over
 * headphones. The decoder supports sound-field rotation for head-tracking and
 * may also accomodate custom HRIR sets via the SOFA standard.
 *
 * Dependencies:
 *     saf_utilities, afSTFTlib, saf_hrir, saf_sh
 * Author, date created:
 *     Leo McCormack, 14.04.2018
 */

#include "ambi_bin_internal.h" 

void ambi_bin_initCodec
(
    void* const hAmbi
)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    codecPars* pars = pData->pars;
    int i, j, k, nSH, order, band;
    const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);
    
    order = pData->new_order;
    nSH = (order+1)*(order+1);
    
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
        free1d((void**)&(pars->hrirs));
        pars->hrirs = malloc1d(pars->N_hrir_dirs * NUM_EARS * (pars->hrir_len)*sizeof(float));
        for(i=0; i<pars->N_hrir_dirs; i++)
            for(j=0; j<NUM_EARS; j++)
                for(k=0; k< pars->hrir_len; k++)
                    pars->hrirs[i*NUM_EARS*(pars->hrir_len) + j*(pars->hrir_len) + k] = (float)__default_hrirs[i][j][k];
        free1d((void**)&(pars->hrir_dirs_deg));
        pars->hrir_dirs_deg = malloc1d(pars->N_hrir_dirs * NUM_EARS * sizeof(float));
        for(i=0; i<pars->N_hrir_dirs; i++)
            for(j=0; j<2; j++)
                pars->hrir_dirs_deg[i*2+j] = (float)__default_hrir_dirs_deg[i][j];
    }
    
    /* estimate the ITDs for each HRIR */
    free1d((void**)&(pars->itds_s));
    estimateITDs(pars->hrirs, pars->N_hrir_dirs, pars->hrir_len, pars->hrir_fs, &(pars->itds_s));
    
    /* convert hrirs to filterbank coefficients */
    free1d((void**)&(pars->hrtf_fb));
    HRIRs2FilterbankHRTFs(pars->hrirs, pars->N_hrir_dirs, pars->hrir_len, pars->itds_s, (float*)pData->freqVector, HYBRID_BANDS, 0, &(pars->hrtf_fb));

    /* get new decoder */
    float_complex* decMtx;
    decMtx = calloc(HYBRID_BANDS*NUM_EARS*nSH, sizeof(float_complex));
    switch(pData->method){
        default:
        case DECODING_METHOD_LS:
            getBinauralAmbiDecoder(pars->hrtf_fb, pars->hrir_dirs_deg, pars->N_hrir_dirs, HYBRID_BANDS, BINAURAL_DECODER_LS, order, pData->freqVector, pars->itds_s, NULL, decMtx);
            break;
        case DECODING_METHOD_LSDIFFEQ:
            getBinauralAmbiDecoder(pars->hrtf_fb, pars->hrir_dirs_deg, pars->N_hrir_dirs, HYBRID_BANDS, BINAURAL_DECODER_LSDIFFEQ, order, pData->freqVector, pars->itds_s, NULL, decMtx);
            break;
        case DECODING_METHOD_SPR:
            getBinauralAmbiDecoder(pars->hrtf_fb, pars->hrir_dirs_deg, pars->N_hrir_dirs, HYBRID_BANDS, BINAURAL_DECODER_SPR, order, pData->freqVector, pars->itds_s, NULL, decMtx);
            break;
        case DECODING_METHOD_TA:
            getBinauralAmbiDecoder(pars->hrtf_fb, pars->hrir_dirs_deg, pars->N_hrir_dirs, HYBRID_BANDS, BINAURAL_DECODER_TA, order, pData->freqVector, pars->itds_s, NULL, decMtx);
            break;
        case DECODING_METHOD_MAGLS:
            getBinauralAmbiDecoder(pars->hrtf_fb, pars->hrir_dirs_deg, pars->N_hrir_dirs, HYBRID_BANDS, BINAURAL_DECODER_MAGLS, order, pData->freqVector, pars->itds_s, NULL, decMtx);
            break;
    }
    
    /* Apply Phase Warping */
    if(pData->enablePhaseWarping){
        // COMING SOON
    }
    
    /* Apply Max RE */
    float* tmp;
    float_complex* a_n, *decMtx_rE;
    if (pData->enableMaxRE){
        tmp = malloc1d(nSH*nSH*sizeof(float));
        a_n = malloc1d(nSH*nSH*sizeof(float_complex));
        decMtx_rE = malloc1d(NUM_EARS*nSH*sizeof(float_complex));
        getMaxREweights(order, tmp);
        for(i=0; i<nSH*nSH; i++)
            a_n[i] = cmplxf(tmp[i], 0.0f);
        
        /* apply per band */
        for(band=0; band<HYBRID_BANDS; band++){
            cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, NUM_EARS, nSH, nSH, &calpha,
                        &decMtx[band*NUM_EARS*nSH], nSH,
                        a_n, nSH, &cbeta,
                        decMtx_rE, nSH);
            memcpy(&decMtx[band*NUM_EARS*nSH], decMtx_rE, NUM_EARS*nSH*sizeof(float_complex));
        }
        free(tmp);
        free(a_n);
        free(decMtx_rE);
    }
    
    /* Apply diffuse Matching/Correction */
    if (pData->enableDiffuseMatching)
        applyDiffCovMatching(pars->hrtf_fb, pars->hrir_dirs_deg, pars->N_hrir_dirs, HYBRID_BANDS, order, NULL, decMtx);
    
    /* replace current decoder */
    memset(pars->M_dec, 0, HYBRID_BANDS*NUM_EARS*MAX_NUM_SH_SIGNALS*sizeof(float_complex));
    for(band=0; band<HYBRID_BANDS; band++)
        for(i=0; i<NUM_EARS; i++)
            for(j=0; j<nSH; j++)
                pars->M_dec[band][i][j] = decMtx[band*2*nSH + i*nSH + j];
    free(decMtx);
     
    pData->order = order;
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

