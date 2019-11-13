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
 *     saf_utilities, afSTFTlib, saf_hrir, saf_vbap, saf_sh, saf_sofa_reader
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
    int i, j, nSH, order, band; 
    
    order = pData->new_order;
    nSH = (order+1)*(order+1);
    
    /* load sofa file or load default hrir data */
    if(!pData->useDefaultHRIRsFLAG && pars->sofa_filepath!=NULL){
        loadSofaFile(pars->sofa_filepath,
                     &(pars->hrirs),
                     &(pars->hrir_dirs_deg),
                     &(pars->N_hrir_dirs),
                     &(pars->hrir_len),
                     &(pars->hrir_fs));
    }
    else{
        loadSofaFile(NULL, /* setting path to NULL loads default HRIR data */
                     &(pars->hrirs),
                     &(pars->hrir_dirs_deg),
                     &(pars->N_hrir_dirs),
                     &(pars->hrir_len),
                     &(pars->hrir_fs));
    }
    
    /* estimate the ITDs for each HRIR */
    pars->itds_s = realloc1d(pars->itds_s, pars->N_hrir_dirs*sizeof(float));
    estimateITDs(pars->hrirs, pars->N_hrir_dirs, pars->hrir_len, pars->hrir_fs, pars->itds_s);
    
    /* convert hrirs to filterbank coefficients */ 
    pars->hrtf_fb = realloc1d(pars->hrtf_fb, HYBRID_BANDS * NUM_EARS * (pars->N_hrir_dirs)*sizeof(float_complex));
    HRIRs2FilterbankHRTFs(pars->hrirs, pars->N_hrir_dirs, pars->hrir_len, pars->hrtf_fb);
    diffuseFieldEqualiseHRTFs(pars->N_hrir_dirs, pars->itds_s, pData->freqVector, HYBRID_BANDS, pars->hrtf_fb);

    /* get new decoder */
    float_complex* decMtx;
    decMtx = calloc(HYBRID_BANDS*NUM_EARS*nSH, sizeof(float_complex));
    switch(pData->method){
        default:
        case DECODING_METHOD_LS:
            getBinauralAmbiDecoderMtx(pars->hrtf_fb, pars->hrir_dirs_deg, pars->N_hrir_dirs, HYBRID_BANDS,
                                      BINAURAL_DECODER_LS, order, pData->freqVector, pars->itds_s, NULL,
                                      pData->enableDiffuseMatching, pData->enableMaxRE, decMtx);
            break;
        case DECODING_METHOD_LSDIFFEQ:
            getBinauralAmbiDecoderMtx(pars->hrtf_fb, pars->hrir_dirs_deg, pars->N_hrir_dirs, HYBRID_BANDS,
                                      BINAURAL_DECODER_LSDIFFEQ, order, pData->freqVector, pars->itds_s, NULL,
                                      pData->enableDiffuseMatching, pData->enableMaxRE, decMtx);
            break;
        case DECODING_METHOD_SPR:
            getBinauralAmbiDecoderMtx(pars->hrtf_fb, pars->hrir_dirs_deg, pars->N_hrir_dirs, HYBRID_BANDS,
                                      BINAURAL_DECODER_SPR, order, pData->freqVector, pars->itds_s, NULL,
                                      pData->enableDiffuseMatching, pData->enableMaxRE, decMtx);
            break;
        case DECODING_METHOD_TA:
            getBinauralAmbiDecoderMtx(pars->hrtf_fb, pars->hrir_dirs_deg, pars->N_hrir_dirs, HYBRID_BANDS,
                                      BINAURAL_DECODER_TA, order, pData->freqVector, pars->itds_s, NULL,
                                      pData->enableDiffuseMatching, pData->enableMaxRE, decMtx);
            break;
        case DECODING_METHOD_MAGLS:
            getBinauralAmbiDecoderMtx(pars->hrtf_fb, pars->hrir_dirs_deg, pars->N_hrir_dirs, HYBRID_BANDS,
                                      BINAURAL_DECODER_MAGLS, order, pData->freqVector, pars->itds_s, NULL,
                                      pData->enableDiffuseMatching, pData->enableMaxRE, decMtx);
            break;
    }
    
    /* Apply Phase Warping */
    if(pData->enablePhaseWarping){
        // COMING SOON
    }
    
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

