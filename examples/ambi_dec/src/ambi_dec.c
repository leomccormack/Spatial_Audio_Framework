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
 * @file ambi_dec.c
 * @brief A frequency-dependent Ambisonic decoder for loudspeakers.
 *
 * Different decoder settings can be specified for the low and high frequencies.
 * When utilising spherical harmonic signals derived from real microphone
 * arrays, this implementation also allows the decoding order per frequency band
 * to be specified; of course, this may also be used creatively. Optionally, a
 * SOFA file may be loaded for personalised headphone listening.
 *
 * The algorithms utilised in this Ambisonic decoder were pieced together and
 * developed in collaboration with Archontis Politis.
 *
 * @author Leo McCormack
 * @date 07.12.2017
 */
 
#include "ambi_dec_internal.h"

void ambi_dec_create
(
    void ** const phAmbi
)
{
    ambi_dec_data* pData = (ambi_dec_data*)malloc1d(sizeof(ambi_dec_data));
    *phAmbi = (void*)pData;
    int i, j, ch, band;

    /* default user parameters */
    pData->masterOrder = pData->new_masterOrder = 1;
    for (band = 0; band<HYBRID_BANDS; band++)
        pData->orderPerBand[band] = 1;
    pData->useDefaultHRIRsFLAG = 1; /* pars->sofa_filepath must be valid to set this to 0 */
    loadLoudspeakerArrayPreset(LOUDSPEAKER_ARRAY_PRESET_T_DESIGN_24, pData->loudpkrs_dirs_deg, &(pData->new_nLoudpkrs), &(pData->loudpkrs_nDims));
    pData->nLoudpkrs = pData->new_nLoudpkrs;
    pData->chOrdering = CH_ACN;
    pData->norm = NORM_SN3D;
    pData->dec_method[0] = DECODING_METHOD_ALLRAD;
    pData->dec_method[1] = DECODING_METHOD_ALLRAD;
    pData->rE_WEIGHT[0] = 1;
    pData->rE_WEIGHT[1] = 1;
    pData->diffEQmode[0] = AMPLITUDE_PRESERVING;
    pData->diffEQmode[1] = ENERGY_PRESERVING;
    pData->transitionFreq = 800.0f;
    
    /* afSTFT stuff */
    pData->hSTFT = NULL;
    pData->STFTInputFrameTF = malloc1d(MAX_NUM_SH_SIGNALS * sizeof(complexVector));
    for(ch=0; ch< MAX_NUM_SH_SIGNALS; ch++) {
        pData->STFTInputFrameTF[ch].re = (float*)calloc1d(HYBRID_BANDS, sizeof(float));
        pData->STFTInputFrameTF[ch].im = (float*)calloc1d(HYBRID_BANDS, sizeof(float));
    }
    pData->tempHopFrameTD = (float**)malloc2d( MAX(MAX_NUM_SH_SIGNALS, MAX_NUM_LOUDSPEAKERS), HOP_SIZE, sizeof(float));
    pData->STFTOutputFrameTF = malloc1d(MAX_NUM_LOUDSPEAKERS * sizeof(complexVector));
    for(ch=0; ch< MAX_NUM_LOUDSPEAKERS; ch++) {
        pData->STFTOutputFrameTF[ch].re = (float*)calloc1d(HYBRID_BANDS, sizeof(float));
        pData->STFTOutputFrameTF[ch].im = (float*)calloc1d(HYBRID_BANDS, sizeof(float));
    }
    
    /* codec data */
    pData->progressBar0_1 = 0.0f;
    pData->progressBarText = malloc1d(AMBI_DEC_PROGRESSBARTEXT_CHAR_LENGTH*sizeof(char));
    strcpy(pData->progressBarText,"");
    pData->codecStatus = CODEC_STATUS_NOT_INITIALISED;
    pData->pars = (ambi_dec_codecPars*)malloc1d(sizeof(ambi_dec_codecPars));
    ambi_dec_codecPars* pars = pData->pars;
    for (i=0; i<NUM_DECODERS; i++){
        for(j=0; j<MAX_SH_ORDER; j++){
            pars->M_dec[i][j] = NULL;
            pars->M_dec_cmplx[i][j] = NULL;
            pars->M_dec_maxrE[i][j] = NULL;
            pars->M_dec_cmplx_maxrE[i][j] = NULL;
        }
    }
    pars->sofa_filepath = NULL;
    pars->hrirs = NULL;
    pars->hrir_dirs_deg = NULL;
    pars->hrtf_vbap_gtableIdx = NULL;
    pars->hrtf_vbap_gtableComp = NULL;
    pars->itds_s = NULL;
    pars->hrtf_fb = NULL;
    pars->hrtf_fb_mag = NULL;
    
    /* internal parameters */ 
    pData->binauraliseLS = pData->new_binauraliseLS = 0;
    
    /* flags */
    pData->procStatus = PROC_STATUS_NOT_ONGOING;
    pData->reinit_hrtfsFLAG = 1;
    for(ch=0; ch<MAX_NUM_LOUDSPEAKERS; ch++)
        pData->recalc_hrtf_interpFLAG[ch] = 1;

    /* set FIFO buffers */
    pData->FIFO_idx = 0;
    memset(pData->inFIFO, 0, MAX_NUM_SH_SIGNALS*FRAME_SIZE*sizeof(float));
    memset(pData->outFIFO, 0, MAX_NUM_LOUDSPEAKERS*FRAME_SIZE*sizeof(float));
}

void ambi_dec_destroy
(
    void ** const phAmbi
)
{
    ambi_dec_data *pData = (ambi_dec_data*)(*phAmbi);
    ambi_dec_codecPars *pars = pData->pars;
    int i, j, ch;
    
    if (pData != NULL) {
        /* not safe to free memory during intialisation/processing loop */
        while (pData->codecStatus == CODEC_STATUS_INITIALISING ||
               pData->procStatus == PROC_STATUS_ONGOING){
            SAF_SLEEP(10);
        }
        
        /* free afSTFT and buffers */
        if(pData->hSTFT!=NULL)
            afSTFTfree(pData->hSTFT);
        if(pData->STFTInputFrameTF!=NULL){
            for (ch = 0; ch< MAX_NUM_SH_SIGNALS; ch++) {
                free(pData->STFTInputFrameTF[ch].re);
                free(pData->STFTInputFrameTF[ch].im);
            }
        }
        if((pData->STFTOutputFrameTF!=NULL)){
            for (ch = 0; ch < MAX_NUM_LOUDSPEAKERS; ch++) {
                free(pData->STFTOutputFrameTF[ch].re);
                free(pData->STFTOutputFrameTF[ch].im);
            }
        }
        free(pData->STFTInputFrameTF);
        free(pData->STFTOutputFrameTF);
        free(pData->tempHopFrameTD);
        free(pars->hrtf_vbap_gtableComp);
        free(pars->hrtf_vbap_gtableIdx);
        free(pars->hrtf_fb);
        free(pars->hrtf_fb_mag);
        free(pars->itds_s);
        free(pars->hrirs);
        free(pars->hrir_dirs_deg);
        for (i=0; i<NUM_DECODERS; i++){
            for(j=0; j<MAX_SH_ORDER; j++){
                free(pars->M_dec[i][j]);
                free(pars->M_dec_cmplx[i][j]);
                free(pars->M_dec_maxrE[i][j]);
                free(pars->M_dec_cmplx_maxrE[i][j]);
            }
        }
        free(pData->progressBarText);
        free(pData);
        pData = NULL;
    }
}

void ambi_dec_init
(
    void * const hAmbi,
    int          sampleRate
)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    int band;
    
    /* define frequency vector */
    pData->fs = sampleRate;
    for(band=0; band <HYBRID_BANDS; band++){
        if(sampleRate == 44100)
            pData->freqVector[band] =  (float)__afCenterFreq44100[band];
        else /* Assume 48kHz */
            pData->freqVector[band] =  (float)__afCenterFreq48e3[band];
    }
}

void ambi_dec_initCodec
(
    void* const hAmbi
)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    ambi_dec_codecPars* pars = pData->pars;
    int i, ch, d, j, n, ng, nGrid_dirs, masterOrder, nSH_order, max_nSH, nLoudspeakers;
    float* grid_dirs_deg, *Y, *M_dec_tmp, *g, *a, *e, *a_n, *hrtf_vbap_gtable;;
    float a_avg[MAX_SH_ORDER], e_avg[MAX_SH_ORDER], azi_incl[2], sum_elev;
    
    if (pData->codecStatus != CODEC_STATUS_NOT_INITIALISED)
        return; /* re-init not required, or already happening */
    while (pData->procStatus == PROC_STATUS_ONGOING){
        /* re-init required, but we need to wait for the current processing loop to end */
        pData->codecStatus = CODEC_STATUS_INITIALISING; /* indicate that we want to init */
        SAF_SLEEP(10);
    }
    
    /* for progress bar */
    pData->codecStatus = CODEC_STATUS_INITIALISING;
    strcpy(pData->progressBarText,"Initialising");
    pData->progressBar0_1 = 0.0f;
    
    /* reinit afSTFT */
    masterOrder = pData->new_masterOrder;
    max_nSH = (masterOrder+1)*(masterOrder+1);
    nLoudspeakers = pData->new_nLoudpkrs;
    if(pData->hSTFT==NULL){
        if(pData->new_binauraliseLS)
            afSTFTinit(&(pData->hSTFT), HOP_SIZE, max_nSH, NUM_EARS, 0, 1);
        else
            afSTFTinit(&(pData->hSTFT), HOP_SIZE, max_nSH, nLoudspeakers, 0, 1);
        afSTFTclearBuffers(pData->hSTFT);
    }
    else{
        if(pData->new_binauraliseLS)
            afSTFTchannelChange(pData->hSTFT, max_nSH, NUM_EARS);
        else
            afSTFTchannelChange(pData->hSTFT, max_nSH, nLoudspeakers);
        afSTFTclearBuffers(pData->hSTFT);
    }
    pData->binauraliseLS = pData->new_binauraliseLS;
    pData->nLoudpkrs = nLoudspeakers;
    
    /* Quick and dirty check to find loudspeaker dimensionality */
    strcpy(pData->progressBarText,"Computing decoder");
    pData->progressBar0_1 = 0.2f;
    sum_elev = 0.0f;
    for(ch=0; ch < nLoudspeakers; ch++)
        sum_elev += fabsf(pData->loudpkrs_dirs_deg[ch][1]);
    if( (((sum_elev < 5.0f) && (sum_elev > -5.0f))) || (nLoudspeakers < 4) )
        pData->loudpkrs_nDims = 2;
    else
        pData->loudpkrs_nDims = 3;
    
    /* add virtual loudspeakers for 2D case */
    if (pData->loudpkrs_nDims == 2){
        assert(nLoudspeakers<=MAX_NUM_LOUDSPEAKERS-2);
        pData->loudpkrs_dirs_deg[nLoudspeakers][0] = 0.0f;
        pData->loudpkrs_dirs_deg[nLoudspeakers][1] = -90.0f;
        pData->loudpkrs_dirs_deg[nLoudspeakers+1][0] = 0.0f;
        pData->loudpkrs_dirs_deg[nLoudspeakers+1][1] = 90.0f;
        nLoudspeakers += 2;
    }
    
    /* prep */
    nGrid_dirs = 480; /* Minimum t-design of degree 30, has 480 points */
    g = malloc1d(nLoudspeakers*sizeof(float));
    a = malloc1d(nGrid_dirs*sizeof(float));
    e = malloc1d(nGrid_dirs*sizeof(float));
    
    /* calculate loudspeaker decoding matrices */
    for( d=0; d<NUM_DECODERS; d++){
        M_dec_tmp = malloc1d(nLoudspeakers * max_nSH * sizeof(float));
        switch(pData->dec_method[d]){
            case DECODING_METHOD_SAD:
                getLoudspeakerAmbiDecoderMtx((float*)pData->loudpkrs_dirs_deg, nLoudspeakers, LOUDSPEAKER_DECODER_SAD, masterOrder, 0, M_dec_tmp);
                break;
            case DECODING_METHOD_MMD:
                getLoudspeakerAmbiDecoderMtx((float*)pData->loudpkrs_dirs_deg, nLoudspeakers, LOUDSPEAKER_DECODER_MMD, masterOrder, 0, M_dec_tmp);
                break;
            case DECODING_METHOD_EPAD:
                getLoudspeakerAmbiDecoderMtx((float*)pData->loudpkrs_dirs_deg, nLoudspeakers, LOUDSPEAKER_DECODER_EPAD, masterOrder, 0, M_dec_tmp);
                break;
            case DECODING_METHOD_ALLRAD:
                getLoudspeakerAmbiDecoderMtx((float*)pData->loudpkrs_dirs_deg, nLoudspeakers, LOUDSPEAKER_DECODER_ALLRAD, masterOrder, 0, M_dec_tmp);
                break;
        }
        
        /* diffuse-field EQ for orders 1..masterOrder */
        for( n=1; n<=masterOrder; n++){
            /* truncate M_dec for each order */
            nSH_order = (n+1)*(n+1);
            free(pars->M_dec[d][n-1]);
            pars->M_dec[d][n-1] = malloc1d(nLoudspeakers* nSH_order * sizeof(float));
            free(pars->M_dec_cmplx[d][n-1]);
            pars->M_dec_cmplx[d][n-1] = malloc1d(nLoudspeakers * nSH_order * sizeof(float_complex));
            for(i=0; i<nLoudspeakers; i++){
                for(j=0; j<nSH_order; j++){
                    pars->M_dec[d][n-1][i*nSH_order+j] = M_dec_tmp[i*max_nSH +j]; /* for applying in the time domain, and... */
                    pars->M_dec_cmplx[d][n-1][i*nSH_order+j] = cmplxf(pars->M_dec[d][n-1][i*nSH_order+j], 0.0f); /* for the time-frequency domain */
                }
            }
            
            /* create dedicated maxrE weighted versions */
            a_n = malloc1d(nSH_order*nSH_order*sizeof(float));
            getMaxREweights(n, 1, a_n); /* weights returned as diagonal matrix */
            free(pars->M_dec_maxrE[d][n-1]);
            pars->M_dec_maxrE[d][n-1] = malloc1d(nLoudspeakers * nSH_order * sizeof(float));
            free(pars->M_dec_cmplx_maxrE[d][n-1]);
            pars->M_dec_cmplx_maxrE[d][n-1] = malloc1d(nLoudspeakers * nSH_order * sizeof(float_complex));
            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nLoudspeakers, nSH_order, nSH_order, 1.0f,
                        pars->M_dec[d][n-1], nSH_order,
                        a_n, nSH_order, 0.0f,
                        pars->M_dec_maxrE[d][n-1], nSH_order); /* for applying in the time domain */
            for(i=0; i<nLoudspeakers * nSH_order; i++)
                pars->M_dec_cmplx_maxrE[d][n-1][i] = cmplxf(pars->M_dec_maxrE[d][n-1][i], 0.0f); /* for the time-frequency domain */
            
            /* fire a plane-wave from each grid direction to find the total energy/amplitude (using non-maxrE weighted versions) */
            Y = malloc1d(nSH_order*sizeof(float));
            grid_dirs_deg = (float*)(&__Tdesign_degree_30_dirs_deg[0][0]);
            for(ng=0; ng<nGrid_dirs; ng++){
                azi_incl[0] = grid_dirs_deg[ng*2]*M_PI/180.0f;
                azi_incl[1] = M_PI/2.0f-grid_dirs_deg[ng*2+1]*M_PI/180.0f;
                getSHreal(n, (float*)azi_incl, 1,  Y);
                cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nLoudspeakers, 1, nSH_order, 1.0f,
                            pars->M_dec[d][n-1], nSH_order,
                            Y, nSH_order, 0.0f,
                            g, 1);
                a[ng] = e[ng] = 0.0f;
                for(i=0; i<nLoudspeakers; i++){
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
            
            /* remove virtual loudspeakers from the decoder */
            if (pData->loudpkrs_nDims == 2){
                pars->M_dec[d][n-1] = realloc1d(pars->M_dec[d][n-1], pData->nLoudpkrs * nSH_order * sizeof(float));
                pars->M_dec_cmplx[d][n-1] = realloc1d(pars->M_dec_cmplx[d][n-1], pData->nLoudpkrs * nSH_order * sizeof(float_complex));
                pars->M_dec_maxrE[d][n-1] = realloc1d(pars->M_dec_maxrE[d][n-1], pData->nLoudpkrs * nSH_order * sizeof(float));
                pars->M_dec_cmplx_maxrE[d][n-1] = realloc1d(pars->M_dec_cmplx_maxrE[d][n-1], pData->nLoudpkrs * nSH_order * sizeof(float_complex));
            }
        }
        free(M_dec_tmp);
    }
    
    /* update order */
    pData->masterOrder = pData->new_masterOrder;
    
    /* Binaural-related initialisations */
    if(pData->reinit_hrtfsFLAG){
        strcpy(pData->progressBarText,"Computing VBAP gain table");
        pData->progressBar0_1 = 0.4f;
        
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
        
        /* generate VBAP gain table for the hrir_dirs */
        hrtf_vbap_gtable = NULL;
        pars->hrtf_vbapTableRes[0] = 2; /* azimuth resolution in degrees */
        pars->hrtf_vbapTableRes[1] = 5; /* elevation resolution in degrees */
        generateVBAPgainTable3D(pars->hrir_dirs_deg, pars->N_hrir_dirs, pars->hrtf_vbapTableRes[0], pars->hrtf_vbapTableRes[1], 1, 0, 0.0f,
                                &hrtf_vbap_gtable, &(pars->N_hrtf_vbap_gtable), &(pars->hrtf_nTriangles));
        if(hrtf_vbap_gtable==NULL){
            /* if generating vbap gain tabled failed, re-calculate with default HRIR set (which is known to triangulate correctly) */
            pData->useDefaultHRIRsFLAG = 1;
            ambi_dec_initCodec(hAmbi);
        }
        
        /* compress VBAP table (i.e. remove the zero elements) */
        pars->hrtf_vbap_gtableComp = realloc1d(pars->hrtf_vbap_gtableComp, pars->N_hrtf_vbap_gtable * 3 * sizeof(float));
        pars->hrtf_vbap_gtableIdx  = realloc1d(pars->hrtf_vbap_gtableIdx,  pars->N_hrtf_vbap_gtable * 3 * sizeof(int));
        compressVBAPgainTable3D(hrtf_vbap_gtable, pars->N_hrtf_vbap_gtable, pars->N_hrir_dirs, pars->hrtf_vbap_gtableComp, pars->hrtf_vbap_gtableIdx);
        
        /* convert hrirs to filterbank coefficients */
        strcpy(pData->progressBarText,"Preparing HRIRs");
        pData->progressBar0_1 = 0.85f;
        pars->hrtf_fb = realloc1d(pars->hrtf_fb, HYBRID_BANDS * NUM_EARS * (pars->N_hrir_dirs)*sizeof(float_complex));
        HRIRs2FilterbankHRTFs(pars->hrirs, pars->N_hrir_dirs, pars->hrir_len, HOP_SIZE, 1, pars->hrtf_fb);
        diffuseFieldEqualiseHRTFs(pars->N_hrir_dirs, pars->itds_s, pData->freqVector, HYBRID_BANDS, pars->hrtf_fb);
        
        /* calculate magnitude responses */
        pars->hrtf_fb_mag = realloc1d(pars->hrtf_fb_mag, HYBRID_BANDS*NUM_EARS*(pars->N_hrir_dirs)*sizeof(float));
        for(i=0; i<HYBRID_BANDS*NUM_EARS* (pars->N_hrir_dirs); i++)
            pars->hrtf_fb_mag[i] = cabsf(pars->hrtf_fb[i]);
        
        /* clean-up */
        free(hrtf_vbap_gtable);
        pData->reinit_hrtfsFLAG = 0;
    }
    
    /* done! */
    strcpy(pData->progressBarText,"Done!");
    pData->progressBar0_1 = 1.0f;
    pData->codecStatus = CODEC_STATUS_INITIALISED;
    
    free(g);
    free(a);
    free(e);
}

void ambi_dec_process
(
    void  *  const hAmbi,
    float ** const inputs,
    float ** const outputs,
    int            nInputs,
    int            nOutputs,
    int            nSamples
)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    ambi_dec_codecPars* pars = pData->pars;
    int s, t, ch, ear, i, band, orderBand, nSH_band, decIdx, nSH;
    const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);

    /* local copies of user parameters */
    int nLoudspeakers, binauraliseLS, masterOrder;
    int orderPerBand[HYBRID_BANDS], rE_WEIGHT[NUM_DECODERS];
    float transitionFreq;
    AMBI_DEC_DIFFUSE_FIELD_EQ_APPROACH diffEQmode[NUM_DECODERS];
    AMBI_DEC_NORM_TYPES norm;
    AMBI_DEC_CH_ORDER chOrdering;
    masterOrder = pData->masterOrder;
    nSH = ORDER2NSH(masterOrder);
    nLoudspeakers = pData->nLoudpkrs;
    memcpy(orderPerBand, pData->orderPerBand, HYBRID_BANDS*sizeof(int));
    transitionFreq = pData->transitionFreq;
    memcpy(diffEQmode, pData->diffEQmode, NUM_DECODERS*sizeof(int));
    binauraliseLS = pData->binauraliseLS;
    norm = pData->norm;
    chOrdering = pData->chOrdering;
    memcpy(rE_WEIGHT, pData->rE_WEIGHT, NUM_DECODERS*sizeof(int));
    
    /* Loop over all samples */
    for(s=0; s<nSamples; s++){
        /* Load input signals into inFIFO buffer */
        for(ch=0; ch<MIN(nInputs,nSH); ch++)
            pData->inFIFO[ch][pData->FIFO_idx] = inputs[ch][s];
        for(; ch<nSH; ch++) /* Zero any channels that were not given */
            pData->inFIFO[ch][pData->FIFO_idx] = 0.0f;

        /* Pull output signals from outFIFO buffer */
        for(ch=0; ch<MIN(nOutputs, binauraliseLS ? NUM_EARS : MAX_NUM_LOUDSPEAKERS); ch++)
            outputs[ch][s] = pData->outFIFO[ch][pData->FIFO_idx];
        for(; ch<nOutputs; ch++) /* Zero any extra channels */
            outputs[ch][s] = 0.0f;

        /* Increment buffer index */
        pData->FIFO_idx++;

        /* Process frame if inFIFO is full and codec is ready for it */
        if (pData->FIFO_idx >= FRAME_SIZE && (pData->codecStatus == CODEC_STATUS_INITIALISED) ) {
            pData->FIFO_idx = 0;
            pData->procStatus = PROC_STATUS_ONGOING;

            /* Load time-domain data */
            switch(chOrdering){
                case CH_ACN:
                    convertHOAChannelConvention((float*)pData->inFIFO, masterOrder, FRAME_SIZE, HOA_CH_ORDER_ACN, HOA_CH_ORDER_ACN, (float*)pData->SHFrameTD);
                    break;
                case CH_FUMA:
                    convertHOAChannelConvention((float*)pData->inFIFO, masterOrder, FRAME_SIZE, HOA_CH_ORDER_FUMA, HOA_CH_ORDER_ACN, (float*)pData->SHFrameTD);
                    break;
            }

            /* account for input normalisation scheme */
            switch(norm){
                case NORM_N3D:  /* already in N3D, do nothing */
                    break;
                case NORM_SN3D: /* convert to N3D */
                    convertHOANormConvention((float*)pData->SHFrameTD, masterOrder, FRAME_SIZE, HOA_NORM_SN3D, HOA_NORM_N3D);
                    break;
                case NORM_FUMA: /* only for first-order, convert to N3D */
                    convertHOANormConvention((float*)pData->SHFrameTD, masterOrder, FRAME_SIZE, HOA_NORM_FUMA, HOA_NORM_N3D);
                    break;
            }

            /* Apply time-frequency transform (TFT) */
            for(t=0; t< TIME_SLOTS; t++) {
                for(ch = 0; ch < nSH; ch++)
                    utility_svvcopy(&(pData->SHFrameTD[ch][t*HOP_SIZE]), HOP_SIZE, pData->tempHopFrameTD[ch]);
                afSTFTforward(pData->hSTFT, (float**)pData->tempHopFrameTD, (complexVector*)pData->STFTInputFrameTF);
                for(band=0; band<HYBRID_BANDS; band++)
                    for(ch=0; ch < nSH; ch++)
                        pData->SHframeTF[band][ch][t] = cmplxf(pData->STFTInputFrameTF[ch].re[band], pData->STFTInputFrameTF[ch].im[band]);
            }

            /* Main processing: */
            /* Decode to loudspeaker set-up */
            memset(pData->outputframeTF, 0, HYBRID_BANDS*MAX_NUM_LOUDSPEAKERS*TIME_SLOTS*sizeof(float_complex));
            for(band=0; band<HYBRID_BANDS; band++){
                orderBand = MAX(MIN(orderPerBand[band], masterOrder),1);
                nSH_band = (orderBand+1)*(orderBand+1);
                decIdx = pData->freqVector[band] < transitionFreq ? 0 : 1; /* different decoder for low (0) and high (1) frequencies */
                if(rE_WEIGHT[decIdx]){
                    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nLoudspeakers, TIME_SLOTS, nSH_band, &calpha,
                                pars->M_dec_cmplx_maxrE[decIdx][orderBand-1], nSH_band,
                                pData->SHframeTF[band], TIME_SLOTS, &cbeta,
                                pData->outputframeTF[band], TIME_SLOTS);
                }
                else{
                    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nLoudspeakers, TIME_SLOTS, nSH_band, &calpha,
                                pars->M_dec_cmplx[decIdx][orderBand-1], nSH_band,
                                pData->SHframeTF[band], TIME_SLOTS, &cbeta,
                                pData->outputframeTF[band], TIME_SLOTS);
                }
                for(i=0; i<nLoudspeakers; i++){
                    for(t=0; t<TIME_SLOTS; t++){
                        if(diffEQmode[decIdx]==AMPLITUDE_PRESERVING)
                            pData->outputframeTF[band][i][t] = crmulf(pData->outputframeTF[band][i][t], pars->M_norm[decIdx][orderBand-1][0]);
                        else
                            pData->outputframeTF[band][i][t] = crmulf(pData->outputframeTF[band][i][t], pars->M_norm[decIdx][orderBand-1][1]);
                    }
                }
            }

            /* binauralise the loudspeaker signals */
            if(binauraliseLS){
                memset(pData->binframeTF, 0, HYBRID_BANDS*NUM_EARS*TIME_SLOTS * sizeof(float_complex));
                /* interpolate hrtfs and apply to each source */
                for (ch = 0; ch < nLoudspeakers; ch++) {
                    if(pData->recalc_hrtf_interpFLAG[ch]){
                        ambi_dec_interpHRTFs(hAmbi, pData->loudpkrs_dirs_deg[ch][0], pData->loudpkrs_dirs_deg[ch][1], pars->hrtf_interp[ch]);
                        pData->recalc_hrtf_interpFLAG[ch] = 0;
                    }
                    for (band = 0; band < HYBRID_BANDS; band++)
                        for (ear = 0; ear < NUM_EARS; ear++)
                            for (t = 0; t < TIME_SLOTS; t++)
                                pData->binframeTF[band][ear][t] = ccaddf(pData->binframeTF[band][ear][t], ccmulf(pData->outputframeTF[band][ch][t], pars->hrtf_interp[ch][band][ear]));
                }

                /* scale by sqrt(number of loudspeakers) */
                for (band = 0; band < HYBRID_BANDS; band++)
                    for (ear = 0; ear < NUM_EARS; ear++)
                        for (t = 0; t < TIME_SLOTS; t++)
                            pData->binframeTF[band][ear][t] = crmulf(pData->binframeTF[band][ear][t], 1.0f/sqrtf((float)nLoudspeakers));
            }


            /* inverse-TFT */
            for(t = 0; t < TIME_SLOTS; t++) {
                for(band = 0; band < HYBRID_BANDS; band++) {
                    if(binauraliseLS){
                        for (ch = 0; ch < NUM_EARS; ch++) {
                            pData->STFTOutputFrameTF[ch].re[band] = crealf(pData->binframeTF[band][ch][t]);
                            pData->STFTOutputFrameTF[ch].im[band] = cimagf(pData->binframeTF[band][ch][t]);
                        }
                    }
                    else {
                        for(ch = 0; ch < nLoudspeakers; ch++) {
                            pData->STFTOutputFrameTF[ch].re[band] = crealf(pData->outputframeTF[band][ch][t]);
                            pData->STFTOutputFrameTF[ch].im[band] = cimagf(pData->outputframeTF[band][ch][t]);
                        }
                    }
                }
                afSTFTinverse(pData->hSTFT, pData->STFTOutputFrameTF, pData->tempHopFrameTD);
                for(ch = 0; ch < MIN(binauraliseLS==1 ? NUM_EARS : nLoudspeakers, nOutputs); ch++)
                    utility_svvcopy(pData->tempHopFrameTD[ch], HOP_SIZE, &(pData->outFIFO[ch][t* HOP_SIZE]));
            }
        }
        else if(pData->FIFO_idx >= FRAME_SIZE){
            /* clear outFIFO if codec was not ready */
            pData->FIFO_idx = 0;
            memset(pData->outFIFO, 0, MAX_NUM_LOUDSPEAKERS*FRAME_SIZE*sizeof(float));
        }
    }

    pData->procStatus = PROC_STATUS_NOT_ONGOING;
}


/* Set Functions */

void ambi_dec_refreshSettings(void* const hAmbi)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    int ch;
    for(ch=0; ch<MAX_NUM_LOUDSPEAKERS; ch++)
        pData->recalc_hrtf_interpFLAG[ch] = 1;
    pData->reinit_hrtfsFLAG = 1;
    ambi_dec_setCodecStatus(hAmbi, CODEC_STATUS_NOT_INITIALISED);
}

void ambi_dec_setMasterDecOrder(void  * const hAmbi, int newValue)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    pData->new_masterOrder = MIN(MAX(newValue,1), MAX_SH_ORDER);
    ambi_dec_setCodecStatus(hAmbi, CODEC_STATUS_NOT_INITIALISED);
    /* FUMA only supports 1st order */
    if(pData->new_masterOrder!=MASTER_ORDER_FIRST && pData->chOrdering == CH_FUMA)
        pData->chOrdering = CH_ACN;
    if(pData->new_masterOrder!=MASTER_ORDER_FIRST && pData->norm == NORM_FUMA)
        pData->norm = NORM_SN3D;
}

void ambi_dec_setDecOrder(void  * const hAmbi, int newValue, int bandIdx)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    pData->orderPerBand[bandIdx] = MIN(MAX(newValue,1), pData->new_masterOrder);
}

void ambi_dec_setDecOrderAllBands(void  * const hAmbi, int newValue)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    int band;
    
    for(band=0; band<HYBRID_BANDS; band++)
        pData->orderPerBand[band] = MIN(MAX(newValue,1), pData->new_masterOrder);
}

void ambi_dec_setLoudspeakerAzi_deg(void* const hAmbi, int index, float newAzi_deg)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    if(newAzi_deg>180.0f)
        newAzi_deg = -360.0f + newAzi_deg;
    newAzi_deg = MAX(newAzi_deg, -180.0f);
    newAzi_deg = MIN(newAzi_deg, 180.0f);
    if(pData->loudpkrs_dirs_deg[index][0] != newAzi_deg){
        pData->loudpkrs_dirs_deg[index][0] = newAzi_deg;
        pData->recalc_hrtf_interpFLAG[index] = 1;
        ambi_dec_setCodecStatus(hAmbi, CODEC_STATUS_NOT_INITIALISED);
    }
}

void ambi_dec_setLoudspeakerElev_deg(void* const hAmbi, int index, float newElev_deg)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    newElev_deg = MAX(newElev_deg, -90.0f);
    newElev_deg = MIN(newElev_deg, 90.0f);
    if(pData->loudpkrs_dirs_deg[index][1] != newElev_deg){
        pData->loudpkrs_dirs_deg[index][1] = newElev_deg;
        pData->recalc_hrtf_interpFLAG[index] = 1;
        ambi_dec_setCodecStatus(hAmbi, CODEC_STATUS_NOT_INITIALISED);
    }
}

void ambi_dec_setNumLoudspeakers(void* const hAmbi, int new_nLoudspeakers)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    int ch; 
    pData->new_nLoudpkrs = new_nLoudspeakers > MAX_NUM_LOUDSPEAKERS ? MAX_NUM_LOUDSPEAKERS : new_nLoudspeakers;
    pData->new_nLoudpkrs = MAX(MIN_NUM_LOUDSPEAKERS, pData->new_nLoudpkrs);
    if(pData->nLoudpkrs != pData->new_nLoudpkrs){
        for(ch=0; ch<MAX_NUM_LOUDSPEAKERS; ch++)
            pData->recalc_hrtf_interpFLAG[ch] = 1;
        ambi_dec_setCodecStatus(hAmbi, CODEC_STATUS_NOT_INITIALISED);
    }
}

void ambi_dec_setBinauraliseLSflag(void* const hAmbi, int newState)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    
    pData->new_binauraliseLS = newState; 
    if(pData->new_binauraliseLS != pData->binauraliseLS)
        ambi_dec_setCodecStatus(hAmbi, CODEC_STATUS_NOT_INITIALISED);
}

void ambi_dec_setUseDefaultHRIRsflag(void* const hAmbi, int newState)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    
    if((!pData->useDefaultHRIRsFLAG) && (newState)){
        pData->useDefaultHRIRsFLAG = newState;
        pData->reinit_hrtfsFLAG = 1;
        ambi_dec_setCodecStatus(hAmbi, CODEC_STATUS_NOT_INITIALISED);
    }
}

void ambi_dec_setSofaFilePath(void* const hAmbi, const char* path)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    ambi_dec_codecPars* pars = pData->pars;
    
    pars->sofa_filepath = malloc1d(strlen(path) + 1);
    strcpy(pars->sofa_filepath, path);
    pData->useDefaultHRIRsFLAG = 0;
    pData->reinit_hrtfsFLAG = 1;
    ambi_dec_setCodecStatus(hAmbi, CODEC_STATUS_NOT_INITIALISED);
}

void ambi_dec_setOutputConfigPreset(void* const hAmbi, int newPresetID)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    int ch;
    
    loadLoudspeakerArrayPreset(newPresetID, pData->loudpkrs_dirs_deg, &(pData->new_nLoudpkrs), &(pData->loudpkrs_nDims));
    for(ch=0; ch<MAX_NUM_LOUDSPEAKERS; ch++)
        pData->recalc_hrtf_interpFLAG[ch] = 1;
    ambi_dec_setCodecStatus(hAmbi, CODEC_STATUS_NOT_INITIALISED);
}

void ambi_dec_setSourcePreset(void* const hAmbi, int newPresetID)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    int band, rangeIdx, curOrder, reverse;
    
    rangeIdx = 0;
    curOrder = 1;
    reverse = 0;
    switch(newPresetID){
        /* Ideal spherical harmonics will have SH_ORDER at all frequencies */
        case MIC_PRESET_IDEAL:
            for(band=0; band<HYBRID_BANDS; band++)
                pData->orderPerBand[band] = pData->masterOrder;
            break;
            
        /* For real microphone arrays, the maximum usable spherical harmonic order will depend on frequency  */
        case MIC_PRESET_ZYLIA:
            for(band=0; band<HYBRID_BANDS; band++){
                if(rangeIdx<2*(__Zylia_maxOrder-1)){
                    if(pData->freqVector[band]>__Zylia_freqRange[rangeIdx]){
                        if(!reverse)
                            curOrder++;
                        else
                            curOrder--;
                        reverse = (curOrder == __Zylia_maxOrder) || (reverse) ? 1 : 0;
                        rangeIdx++;
                    }
                }
                pData->orderPerBand[band] = MIN(pData->masterOrder,curOrder);
            }
            break;

        case MIC_PRESET_EIGENMIKE32:
            for(band=0; band<HYBRID_BANDS; band++){
                if(rangeIdx<2*(__Eigenmike32_maxOrder-1)){
                    if(pData->freqVector[band]>__Eigenmike32_freqRange[rangeIdx]){
                        if(!reverse)
                            curOrder++;
                        else
                            curOrder--;
                        reverse = (curOrder == __Eigenmike32_maxOrder) || (reverse) ? 1 : 0;
                        rangeIdx++;
                    }
                }
                pData->orderPerBand[band] = MIN(pData->masterOrder,curOrder);
            }
            break;

        case MIC_PRESET_DTU_MIC:
            for(band=0; band<HYBRID_BANDS; band++){
                if(rangeIdx<2*(__DTU_mic_maxOrder-1)){
                    if(pData->freqVector[band]>__DTU_mic_freqRange[rangeIdx]){
                        if(!reverse)
                            curOrder++;
                        else
                            curOrder--;
                        reverse = (curOrder == __DTU_mic_maxOrder) || (reverse) ? 1 : 0;
                        rangeIdx++;
                    }
                }
                pData->orderPerBand[band] = MIN(pData->masterOrder,curOrder);
            }
            break;
    }
}

void ambi_dec_setChOrder(void* const hAmbi, int newOrder)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    if((AMBI_DEC_CH_ORDER)newOrder != CH_FUMA || pData->new_masterOrder==MASTER_ORDER_FIRST) /* FUMA only supports 1st order */
        pData->chOrdering = (AMBI_DEC_CH_ORDER)newOrder;
}

void ambi_dec_setNormType(void* const hAmbi, int newType)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    if((AMBI_DEC_NORM_TYPES)newType != NORM_FUMA || pData->new_masterOrder==MASTER_ORDER_FIRST) /* FUMA only supports 1st order */
        pData->norm = (AMBI_DEC_NORM_TYPES)newType;
}

void ambi_dec_setDecMethod(void* const hAmbi, int index, int newID)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    pData->dec_method[index] = newID;
    ambi_dec_setCodecStatus(hAmbi, CODEC_STATUS_NOT_INITIALISED);
}

void ambi_dec_setDecEnableMaxrE(void* const hAmbi, int index, int newID)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    pData->rE_WEIGHT[index] = newID;
}

void ambi_dec_setDecNormType(void* const hAmbi, int index, int newID)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    pData->diffEQmode[index] = newID;
}

void ambi_dec_setTransitionFreq(void* const hAmbi, float newValue)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    pData->transitionFreq = CLAMP(newValue, AMBI_DEC_TRANSITION_MIN_VALUE, AMBI_DEC_TRANSITION_MAX_VALUE);
}


/* Get Functions */

AMBI_DEC_CODEC_STATUS ambi_dec_getCodecStatus(void* const hAmbi)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    return pData->codecStatus;
}

float ambi_dec_getProgressBar0_1(void* const hAmbi)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    return pData->progressBar0_1;
}

void ambi_dec_getProgressBarText(void* const hAmbi, char* text)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    memcpy(text, pData->progressBarText, AMBI_DEC_PROGRESSBARTEXT_CHAR_LENGTH*sizeof(char));
}

int ambi_dec_getMasterDecOrder(void  * const hAmbi)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    return pData->new_masterOrder;
}

int ambi_dec_getDecOrder(void  * const hAmbi, int bandIdx)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    return pData->orderPerBand[bandIdx];
}

int ambi_dec_getDecOrderAllBands(void  * const hAmbi)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    return pData->orderPerBand[0];
}

void ambi_dec_getDecOrderHandle
(
    void* const hAmbi,
    float** pX_vector,
    int** pY_values,
    int* pNpoints
)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    (*pX_vector) = &pData->freqVector[0];
    (*pY_values) = &pData->orderPerBand[0];
    (*pNpoints) = HYBRID_BANDS;
}

int ambi_dec_getNumberOfBands(void)
{
    return HYBRID_BANDS;
}

float ambi_dec_getLoudspeakerAzi_deg(void* const hAmbi, int index)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    return pData->loudpkrs_dirs_deg[index][0];
}

float ambi_dec_getLoudspeakerElev_deg(void* const hAmbi, int index)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    return pData->loudpkrs_dirs_deg[index][1];
}

int ambi_dec_getNumLoudspeakers(void* const hAmbi)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    return pData->new_nLoudpkrs;
}

int ambi_dec_getMaxNumLoudspeakers()
{
    return MAX_NUM_LOUDSPEAKERS;
}

int  ambi_dec_getNSHrequired(void* const hAmbi)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    return (pData->masterOrder+1)*(pData->masterOrder+1);
}

int ambi_dec_getBinauraliseLSflag(void* const hAmbi)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    return pData->new_binauraliseLS;
}

int ambi_dec_getUseDefaultHRIRsflag(void* const hAmbi)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    return pData->useDefaultHRIRsFLAG;
}

char* ambi_dec_getSofaFilePath(void* const hAmbi)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    ambi_dec_codecPars* pars = pData->pars;
    if(pars->sofa_filepath!=NULL)
        return pars->sofa_filepath;
    else
        return "no_file";
}

int ambi_dec_getChOrder(void* const hAmbi)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    return (int)pData->chOrdering;
}

int ambi_dec_getNormType(void* const hAmbi)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    return (int)pData->norm;
}

int ambi_dec_getDecMethod(void* const hAmbi, int index)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    return pData->dec_method[index];
}

int ambi_dec_getDecEnableMaxrE(void* const hAmbi, int index)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    return pData->rE_WEIGHT[index];
}

int ambi_dec_getDecNormType(void* const hAmbi, int index)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    return pData->diffEQmode[index];
}

float ambi_dec_getTransitionFreq(void* const hAmbi)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    return pData->transitionFreq;
}

int ambi_dec_getHRIRsamplerate(void* const hAmbi)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    ambi_dec_codecPars* pars = pData->pars;
    return pars->hrir_fs;
}

int ambi_dec_getDAWsamplerate(void* const hAmbi)
{
    ambi_dec_data *pData = (ambi_dec_data*)(hAmbi);
    return pData->fs;
}

int ambi_dec_getProcessingDelay()
{
    return FRAME_SIZE + 12*HOP_SIZE;
}



