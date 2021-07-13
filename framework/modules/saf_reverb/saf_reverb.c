/*
 * Copyright 2020 Leo McCormack
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
 * @file saf_reverb.c
 * @ingroup Reverb
 * @brief Public source for the reverb processing module (#SAF_REVERB_MODULE)
 *
 * A collection of reverb and room simulation algorithms.
 *
 * @author Leo McCormack
 * @date 06.05.2020
 * @license ISC
 */
 
#include "saf_reverb.h"
#include "saf_reverb_internal.h"

/* ========================================================================== */
/*                         IMS Shoebox Room Simulator                         */
/* ========================================================================== */

void ims_shoebox_create
(
    void** phIms,
    float roomDimensions[3],
    float* abs_wall,
    float lowestOctaveBand,
    int nOctBands,
    float c_ms,
    float fs
)
{
    *phIms = malloc1d(sizeof(ims_scene_data));
    ims_scene_data *sc = (ims_scene_data*)(*phIms);
    int i,j,band,wall;

    /* Shoebox dimensions */
    sc->room_dims[0] = roomDimensions[0];
    sc->room_dims[1] = roomDimensions[1];
    sc->room_dims[2] = roomDimensions[2];
    sc->c_ms = c_ms;

    /* Octave band centre frequencies */
    if(nOctBands>1){
        sc->nBands = nOctBands;
        sc->band_centerfreqs = malloc1d(nOctBands*sizeof(float));
        sc->band_centerfreqs[0] = lowestOctaveBand;
        for(band=1; band<nOctBands; band++)
            sc->band_centerfreqs[band] = sc->band_centerfreqs[band-1]*2.0f;
        sc->band_cutofffreqs = malloc1d((sc->nBands-1)*sizeof(float));
        getOctaveBandCutoffFreqs(sc->band_centerfreqs, sc->nBands, sc->band_cutofffreqs);
    }
    else { /* Broad-band operation */
        sc->nBands = 1;
        sc->band_centerfreqs = NULL;
        sc->band_cutofffreqs  = NULL;
    }

    /* Samplerate */
    sc->fs = fs;

    /* Absorption coeffients per wall and octave band */
    sc->abs_wall = (float**)malloc2d(sc->nBands, IMS_NUM_WALLS_SHOEBOX, sizeof(float));
    for(band=0; band<sc->nBands; band++)
        for(wall=0; wall<IMS_NUM_WALLS_SHOEBOX; wall++)
            sc->abs_wall[band][wall] = abs_wall[band*IMS_NUM_WALLS_SHOEBOX+wall];

    /* Default are no sources or receivers in the room */
    for(i=0; i<IMS_MAX_NUM_SOURCES; i++)
        sc->srcs[i].ID = IMS_UNASSIGNED; /* -1 indicates not in use */
    for(i=0; i<IMS_MAX_NUM_RECEIVERS; i++)
        sc->recs[i].ID = IMS_UNASSIGNED;
    sc->nSources = 0;
    sc->nReceivers = 0;

    /* ims_core_workspace per source / receiver combination */
    sc->hCoreWrkSpc = (voidPtr**)malloc2d(IMS_MAX_NUM_RECEIVERS, IMS_MAX_NUM_SOURCES, sizeof(voidPtr));
    for(i=0; i<IMS_MAX_NUM_RECEIVERS; i++)
        for(j=0; j<IMS_MAX_NUM_SOURCES; j++)
            sc->hCoreWrkSpc[i][j] = NULL;

    /* FIR Fiterbank */
    sc->H_filt = NULL;

    /* RIRs per source / receiver combination  */
    sc->rirs = (ims_rir**)malloc2d(IMS_MAX_NUM_RECEIVERS, IMS_MAX_NUM_SOURCES, sizeof(ims_rir));
    for(i=0; i<IMS_MAX_NUM_RECEIVERS; i++){
        for(j=0; j<IMS_MAX_NUM_SOURCES; j++){
            sc->rirs[i][j].data = NULL;
            sc->rirs[i][j].length = sc->rirs[i][j].nChannels = 0;
        }
    }

    /* Circular buffers (only used/allocated when applyEchogramTD() function is called for the first time) */
    memset(sc->wIdx, 0, IMS_MAX_NUM_RECEIVERS*IMS_MAX_NUM_SOURCES*2*sizeof(unsigned long));
    sc->circ_buffer[0] = NULL;
    sc->circ_buffer[1] = NULL;

    /* IIR Filterbank per source (only used/allocated when applyEchogramTD() function is called for the first time) */
    sc->hFaFbank = malloc1d(IMS_MAX_NUM_SOURCES*sizeof(voidPtr));
    sc->src_sigs_bands = malloc1d(IMS_MAX_NUM_SOURCES*sizeof(float**));
    for(j=0; j<IMS_MAX_NUM_SOURCES; j++){
        sc->hFaFbank[j] = NULL;
        sc->src_sigs_bands[j] = NULL;
    }

    /* Temp buffers for cross-fading (only used/allocated when applyEchogramTD() function is called for the first time) */
    sc->rec_sig_tmp[IMS_EG_CURRENT] = malloc1d(IMS_MAX_NUM_RECEIVERS*sizeof(float**));
    sc->rec_sig_tmp[IMS_EG_PREV] = malloc1d(IMS_MAX_NUM_RECEIVERS*sizeof(float**));
    for(j=0; j<IMS_MAX_NUM_RECEIVERS; j++){
        sc->rec_sig_tmp[IMS_EG_CURRENT][j] = NULL;
        sc->rec_sig_tmp[IMS_EG_PREV][j] = NULL;
    }
    memset(sc->applyCrossFadeFLAG, 0, IMS_MAX_NUM_RECEIVERS*IMS_MAX_NUM_SOURCES*sizeof(int));
    sc->interpolator_fIn = NULL;
    sc->interpolator_fOut = NULL;
    sc->tmp_frame = NULL;
    sc->framesize = -1;

    /* Lagrange interpolator look-up table */
    for(i=0; i<IMS_LAGRANGE_LOOKUP_TABLE_SIZE; i++)
        sc->lookup_fractions[i] = 1.0f/IMS_LAGRANGE_LOOKUP_TABLE_SIZE * (float)i;
    lagrangeWeights(IMS_LAGRANGE_ORDER, sc->lookup_fractions, IMS_LAGRANGE_LOOKUP_TABLE_SIZE, (float*)sc->lookup_H_frac);
}

void ims_shoebox_destroy
(
    void** phIms
)
{
    ims_scene_data *sc = (ims_scene_data*)(*phIms);
    int i,j;

    if(sc!=NULL){
        free(sc->band_centerfreqs);
        free(sc->band_cutofffreqs);
        free(sc->abs_wall);
        for(i=0; i<IMS_MAX_NUM_RECEIVERS; i++)
            for(j=0; j<IMS_MAX_NUM_SOURCES; j++)
                ims_shoebox_coreWorkspaceDestroy(&(sc->hCoreWrkSpc[i][j]));
        free(sc->hCoreWrkSpc);
        free(sc->H_filt);
        for(i=0; i<IMS_MAX_NUM_RECEIVERS; i++)
            for(j=0; j<IMS_MAX_NUM_SOURCES; j++)
                free(sc->rirs[i][j].data);
        free(sc->rirs);
        free(sc->circ_buffer[IMS_EG_CURRENT]);
        free(sc->circ_buffer[IMS_EG_PREV]);
        for(j=0; j<IMS_MAX_NUM_SOURCES; j++){
            faf_IIRFilterbank_destroy(&(sc->hFaFbank[j]));
            free(sc->src_sigs_bands[j]);
        }
        free(sc->hFaFbank);
        free(sc->src_sigs_bands);
        for(i=0; i<IMS_MAX_NUM_RECEIVERS; i++){
            free(sc->rec_sig_tmp[IMS_EG_CURRENT][i]);
            free(sc->rec_sig_tmp[IMS_EG_PREV][i]);
        }
        free(sc->rec_sig_tmp[IMS_EG_CURRENT]);
        free(sc->rec_sig_tmp[IMS_EG_PREV]);
        free(sc->interpolator_fIn);
        free(sc->interpolator_fOut);
        free(sc->tmp_frame);
        free(sc);
        sc=NULL;
        *phIms = NULL;
    }
}

void ims_shoebox_computeEchograms
(
    void* hIms,
    int maxN,
    float maxTime_ms
)
{
    ims_scene_data *sc = (ims_scene_data*)(hIms);
    ims_core_workspace* workspace;
    ims_pos_xyz src2, rec2;
    int src_idx, rec_idx, band;

    saf_assert(maxN<0 || maxTime_ms<0.0f, "one of these input arguments must be the same or greater than 0, and the other one must be less than 0.");
    saf_assert(maxN>=0 || maxTime_ms>0.0f, "one of these input arguments must be the same or greater than 0, and the other one must be less than 0.");

    /* Compute echograms for active source/receiver combinations */
    for(rec_idx = 0; rec_idx < IMS_MAX_NUM_RECEIVERS; rec_idx++){
        for(src_idx = 0; src_idx < IMS_MAX_NUM_SOURCES; src_idx++){
            if( (sc->srcs[src_idx].ID != IMS_UNASSIGNED) && (sc->recs[rec_idx].ID != IMS_UNASSIGNED) ){
                /* Change y coord for Receiver and Source to match convention
                 * used inside the coreInit function */
                rec2.x = sc->recs[rec_idx].pos.x;
                rec2.y = sc->room_dims[1] - sc->recs[rec_idx].pos.y;
                rec2.z = sc->recs[rec_idx].pos.z;
                src2.x = sc->srcs[src_idx].pos.x;
                src2.y = sc->room_dims[1] - sc->srcs[src_idx].pos.y;
                src2.z = sc->srcs[src_idx].pos.z;

                /* Workspace handle for this source/receiver combination */
                workspace = sc->hCoreWrkSpc[rec_idx][src_idx];

                /* Copy previous echograms */
                for(band=0; band<workspace->nBands; band++)
                    ims_shoebox_echogramCopy(workspace->hEchogram_abs[band], workspace->hPrevEchogram_abs[band]);

                /* Force refresh if target RIR length or max reflection order has changed */
                if(maxTime_ms>0.0f){
                    if(workspace->d_max != maxTime_ms)
                        workspace->refreshEchogramFLAG = 1;
                }
                else{
                    if(workspace->N_max != maxN)
                        workspace->refreshEchogramFLAG = 1;
                }

                /* Only update if it is required */
                if(workspace->refreshEchogramFLAG){
                    /* Compute echogram due to pure propagation (frequency-independent, omni-directional) */
                    if(maxTime_ms>0.0f)
                        ims_shoebox_coreInitT(workspace, sc->room_dims, src2, rec2, maxTime_ms, sc->c_ms);
                    else
                        ims_shoebox_coreInitN(workspace, sc->room_dims, src2, rec2, maxN, sc->c_ms);

                    /* Apply receiver directivities */
                    switch(sc->recs[rec_idx].type){
                        case RECEIVER_SH:
                            ims_shoebox_coreRecModuleSH(workspace, NSH2ORDER(sc->recs[rec_idx].nChannels));
                            break;
                    }

                    /* Apply boundary absorption per frequency band */
                    ims_shoebox_coreAbsorptionModule(workspace, sc->abs_wall);

                    /* Indicate that the echogram is now up to date, and that the RIR should be updated */
                    workspace->refreshEchogramFLAG = 0;
                    workspace->refreshRIRFLAG = 1;

                    /* Also indicate that applyTD() should cross-fade the next frame to void clicks */
                    sc->applyCrossFadeFLAG[rec_idx][src_idx] = 1;
                }
            }
        }
    }  
}

void ims_shoebox_renderRIRs
(
    void* hIms,
    int fractionalDelayFLAG
)
{
    ims_scene_data *sc = (ims_scene_data*)(hIms);
    ims_core_workspace* wrk;
    int src_idx, rec_idx;

    /* Compute FIR Filterbank coefficients (if this is the first time this
     * function is being called) */
    if(sc->H_filt==NULL){
        sc->H_filt = (float**)realloc2d((void**)sc->H_filt, sc->nBands, (IMS_FIR_FILTERBANK_ORDER+1), sizeof(float));
        FIRFilterbank(IMS_FIR_FILTERBANK_ORDER, sc->band_cutofffreqs, sc->nBands-1,
                      sc->fs, WINDOWING_FUNCTION_HAMMING, 1, FLATTEN2D(sc->H_filt));
    }

    /* Render RIRs for all active source/receiver combinations */
    for(rec_idx = 0; rec_idx < IMS_MAX_NUM_RECEIVERS; rec_idx++){
        for(src_idx = 0; src_idx < IMS_MAX_NUM_SOURCES; src_idx++){
            if( (sc->srcs[src_idx].ID!=IMS_UNASSIGNED) && (sc->recs[rec_idx].ID!=IMS_UNASSIGNED) ){

                /* Workspace handle for this source/receiver combination */
                wrk = sc->hCoreWrkSpc[rec_idx][src_idx];

                /* Only update if it is required */
                if(wrk->refreshRIRFLAG){
                    /* Render the RIRs for each band  */
                    ims_shoebox_renderRIR(wrk, fractionalDelayFLAG, sc->fs, sc->H_filt, &(sc->rirs[rec_idx][src_idx]));

                    wrk->refreshRIRFLAG = 0;
                }
            }
        }
    }
}

void ims_shoebox_applyEchogramTD
(
    void* hIms, 
    long receiverID,
    int nSamples,
    int fractionalDelaysFLAG
)
{
    ims_scene_data *sc = (ims_scene_data*)(hIms);
    ims_core_workspace* wrk;
    echogram_data *echogram_abs, *echogram_abs_0;
    int k, i, n, im, band, ch, rec_idx, src_idx, time_samples, wIdx_n;
    unsigned long rIdx;

    saf_assert(fractionalDelaysFLAG==0, "Untested!");
    saf_assert(nSamples <= IMS_MAX_NSAMPLES_PER_FRAME, "nSamples exceeds the maximum number that ims_shoebox_applyEchogramTD() can process at a time");

    /* Allocate circular buffers (if this is the first time this function is being called) */
    if(sc->circ_buffer[0] == NULL)
        sc->circ_buffer[0] = (float***)calloc3d(IMS_MAX_NUM_SOURCES, sc->nBands, IMS_CIRC_BUFFER_LENGTH, sizeof(float));
    if(sc->circ_buffer[1] == NULL)
        sc->circ_buffer[1] = (float***)calloc3d(IMS_MAX_NUM_SOURCES, sc->nBands, IMS_CIRC_BUFFER_LENGTH, sizeof(float));

    /* Also allocate signal buffers and filterbank handles (if this is the first time this function is being called) */
    for(src_idx = 0; src_idx < IMS_MAX_NUM_SOURCES; src_idx++){
        /* only allocate if source is active... */
        if( (sc->srcs[src_idx].ID != IMS_UNASSIGNED) && (sc->src_sigs_bands[src_idx] == NULL) ){
            if(sc->nBands>1){ /* ... and only create the filterbank if there is more than one band: */
                faf_IIRFilterbank_create(&(sc->hFaFbank[src_idx]), IMS_IIR_FILTERBANK_ORDER, sc->band_cutofffreqs,
                                         sc->nBands-1, sc->fs, IMS_MAX_NSAMPLES_PER_FRAME);
            }
            sc->src_sigs_bands[src_idx] = (float**)malloc2d(sc->nBands, IMS_MAX_NSAMPLES_PER_FRAME, sizeof(float)); 
        }
    }

    /* Find index corresponding to this receiver ID */
    rec_idx = IMS_UNASSIGNED;
    for(i=0; i<IMS_MAX_NUM_RECEIVERS; i++){
        if(sc->recs[i].ID == receiverID){
            rec_idx = i;
            break;
        }
    }
    saf_assert(rec_idx != IMS_UNASSIGNED, "Invalid receiverID");

    /* Allocate temporary buffer (if this is the first time this function is being called)  */
    if( (sc->recs[rec_idx].ID != IMS_UNASSIGNED) && (sc->rec_sig_tmp[IMS_EG_CURRENT][rec_idx] == NULL) ){
        sc->rec_sig_tmp[IMS_EG_CURRENT][rec_idx] = (float**)malloc2d(sc->recs[rec_idx].nChannels, IMS_MAX_NSAMPLES_PER_FRAME, sizeof(float));
        sc->rec_sig_tmp[IMS_EG_PREV][rec_idx] = (float**)malloc2d(sc->recs[rec_idx].nChannels, IMS_MAX_NSAMPLES_PER_FRAME, sizeof(float));
    }
    if(sc->framesize!=nSamples){
        sc->framesize = nSamples;
        sc->interpolator_fIn = realloc1d(sc->interpolator_fIn, nSamples*sizeof(float));
        sc->interpolator_fOut = realloc1d(sc->interpolator_fOut, nSamples*sizeof(float));
        sc->tmp_frame = realloc1d(sc->tmp_frame, nSamples*sizeof(float));
        for(i=0; i<nSamples; i++){
            sc->interpolator_fIn[i] = (i+1)*1.0f/(float)nSamples;
            sc->interpolator_fOut[i] = 1.0f-sc->interpolator_fIn[i];
        }
    }

    /* Initialise all receiver channels with zeros */
    for(ch=0; ch<sc->recs[rec_idx].nChannels; ch++) /* (looping over channels since this array is not guaranteed to be contiguous) */
        memset(sc->recs[rec_idx].sigs[ch], 0, nSamples * sizeof(float));
 
    /* Process all active sources for this specific receiver, directly in the time-domain */
    for(src_idx = 0; src_idx < IMS_MAX_NUM_SOURCES; src_idx++){
        if( (sc->srcs[src_idx].ID!=IMS_UNASSIGNED) && (sc->recs[rec_idx].ID!=IMS_UNASSIGNED) ){

            /* Broad-band operation */
            if(sc->nBands==1)
                memcpy(sc->src_sigs_bands[src_idx][0], sc->srcs[src_idx].sig, nSamples*sizeof(float));
            else /* OR: Pass source signal through the Favrot & Faller (power-complementary) IIR filterbank */
                faf_IIRFilterbank_apply(sc->hFaFbank[src_idx], sc->srcs[src_idx].sig, sc->src_sigs_bands[src_idx], nSamples);

            /* Workspace handle for this source/receiver combination */
            wrk = sc->hCoreWrkSpc[rec_idx][src_idx];

            /* k=0 is for the current echogram,
             * k=1 (when applyCrossFadeFLAG is enabled) is for the previous echogram */
            for(k=sc->applyCrossFadeFLAG[rec_idx][src_idx]; k>=0; k--){

                /* Loop over samples */
                for(n=0; n<nSamples; n++){

                    /* Determine write index */
                    wIdx_n = sc->wIdx[k][rec_idx][src_idx] & IMS_CIRC_BUFFER_LENGTH_MASK;

                    /* Since the time vector is the same across bands, it makes sense to determine the read-indices only once... */
                    if(k==1)
                        echogram_abs_0 = (echogram_data*)wrk->hPrevEchogram_abs[0];
                    else
                        echogram_abs_0 = (echogram_data*)wrk->hEchogram_abs[0];

                    /* Handle the special case of an empty echogram */
                    if(echogram_abs_0->numImageSources==0){
                        /* Set output to 0... */
                        for(ch=0; ch<sc->recs[rec_idx].nChannels; ch++)
                            sc->rec_sig_tmp[k][rec_idx][ch][n] = 0.0f;

                        /* Store current sample (per band) into the circular buffer */
                        for(band=0; band < sc->nBands; band++){
                            sc->circ_buffer[k][src_idx][band][wIdx_n] = sc->src_sigs_bands[src_idx][band][n];
                            if(sc->applyCrossFadeFLAG[rec_idx][src_idx]==0)
                                sc->circ_buffer[IMS_EG_PREV][src_idx][band][wIdx_n] = sc->src_sigs_bands[src_idx][band][n];
                        }

                        /* Increment write index */
                        sc->wIdx[k][rec_idx][src_idx]++;
                        if(sc->applyCrossFadeFLAG[rec_idx][src_idx]==0)
                            sc->wIdx[IMS_EG_PREV][rec_idx][src_idx]++;

                        continue; /* to next sample... */
                    }

                    /* Convert time from seconds to number of samples */
                    memset(echogram_abs_0->tmp1, 0, echogram_abs_0->numImageSources*sizeof(float));
                    cblas_saxpy(echogram_abs_0->numImageSources, (sc->fs), echogram_abs_0->time, 1, echogram_abs_0->tmp1, 1);

                    /* Determine read-indices and (optionally) also the interpolation weights */
                    if(fractionalDelaysFLAG){
                        /* Loop over all image sources, and determine the circular buffer read indices */
                        for(im=0; im <echogram_abs_0->numImageSources; im++){
                            /* Base read-index */
                            time_samples = (int)(echogram_abs_0->tmp1[im]);                /* FLOOR */
                            rIdx = IMS_CIRC_BUFFER_LENGTH-time_samples + sc->wIdx[k][rec_idx][src_idx] /* read index for this image source */
                                   + (IMS_LAGRANGE_ORDER/2);                               /* in order to correctly centre the filter */
                            echogram_abs_0->rIdx[im] = rIdx & IMS_CIRC_BUFFER_LENGTH_MASK; /* wrap-around if needed */
                        }

                        /* Find fractional parts */
                        utility_svmod(echogram_abs_0->tmp1, echogram_abs_0->ones_dummy, echogram_abs_0->numImageSources, echogram_abs_0->tmp2);

                        /* Find read-indices for interpolator */
                        for(im=0; im <echogram_abs_0->numImageSources; im++){
                            /* Centre the filter */
                            echogram_abs_0->tmp2[im] += (float)(IMS_LAGRANGE_ORDER/2);

                            /* Read-indices for lagrange interpolation */
                            for(i=1; i<IMS_LAGRANGE_ORDER; i++){
                                echogram_abs_0->rIdx_frac[i-1][im] = echogram_abs_0->rIdx[im] - i;
                                /* Wrap around if needed */
                                if(echogram_abs_0->rIdx_frac[i-1][im]<0)
                                    echogram_abs_0->rIdx_frac[i-1][im] += IMS_CIRC_BUFFER_LENGTH;
                            }
                        }

                        /* Compute interpolation weights */ // TODO: This line is around 50% CPU usage of the whole function, a look-up table would be faster...
                        lagrangeWeights(IMS_LAGRANGE_ORDER, echogram_abs_0->tmp2, echogram_abs_0->numImageSources, FLATTEN2D(echogram_abs_0->h_frac));
                    }
                    else{
                        /* Loop over all image sources, and determine the circular buffer read indices based on the nearest sample */
                        for(im=0; im <echogram_abs_0->numImageSources; im++){
                            time_samples = (int)(echogram_abs_0->tmp1[im] + 0.5f);         /* ROUND to nearest sample */
                            rIdx = IMS_CIRC_BUFFER_LENGTH-time_samples + sc->wIdx[k][rec_idx][src_idx]; /* read index for this image source */
                            echogram_abs_0->rIdx[im] = rIdx & IMS_CIRC_BUFFER_LENGTH_MASK; /* wrap-around if needed */
                        }
                    }

                    /* Loop over octave bands */
                    for(band=0; band < sc->nBands; band++){
                        /* Echogram for this source/receiver at this band */
                        if(k==1)
                            echogram_abs = (echogram_data*)wrk->hPrevEchogram_abs[band];
                        else
                            echogram_abs = (echogram_data*)wrk->hEchogram_abs[band];
                        saf_assert(echogram_abs_0->numImageSources == echogram_abs->numImageSources, "The below code is assuming that the number of image sources should be the same across octave bands!");

                        /* Pull values from the circular buffer corresponding to these read indices, and store this sparse vector as a "compressed" vector */
                        utility_ssv2cv_inds(sc->circ_buffer[k][src_idx][band], echogram_abs_0->rIdx, echogram_abs->numImageSources, echogram_abs->cb_vals[0]);

                        /* Apply interpolation if enabled */
                        if(fractionalDelaysFLAG){
                            /* Apply the weights corresponding to the first tap in the filter */
                            utility_svvmul(echogram_abs->cb_vals[0], echogram_abs_0->h_frac[0], echogram_abs->numImageSources, echogram_abs->tmp1);

                            /* Step through the rest of the filter */
                            for(i=1; i<IMS_LAGRANGE_ORDER; i++){
                                /* Pull values from the circular buffer corresponding to the read-indices for the current tap in the filter */
                                utility_ssv2cv_inds(sc->circ_buffer[k][src_idx][band], echogram_abs_0->rIdx_frac[i-1], echogram_abs->numImageSources, echogram_abs->cb_vals[0]);

                                /* Apply interpolation weights */
                                utility_svvmul(echogram_abs->cb_vals[0], echogram_abs_0->h_frac[i], echogram_abs->numImageSources, echogram_abs->tmp2);

                                /* Sum */
                                cblas_saxpy(echogram_abs->numImageSources, 1.0f, echogram_abs->tmp2, 1, echogram_abs->tmp1, 1);
                            }

                            /* Copy result */
                            cblas_scopy(echogram_abs->numImageSources, echogram_abs->tmp1, 1, echogram_abs->cb_vals[0], 1);
                        }

                        /* Replicate these circular buffer values for all output channels */
                        for(ch=1; ch<sc->recs[rec_idx].nChannels; ch++)
                            cblas_scopy(echogram_abs->numImageSources, echogram_abs->cb_vals[0], 1, echogram_abs->cb_vals[ch], 1);

                        /* Apply the echogram scalings to each image source - for all channels */
                        utility_svvmul(FLATTEN2D(echogram_abs->value), FLATTEN2D(echogram_abs->cb_vals),
                                       (sc->recs[rec_idx].nChannels)*echogram_abs->numImageSources, FLATTEN2D(echogram_abs->contrib));

                        /* Render frame */
                        for(ch=0; ch<sc->recs[rec_idx].nChannels; ch++)
                            sc->rec_sig_tmp[k][rec_idx][ch][n] = cblas_sdot(echogram_abs->numImageSources, echogram_abs->ones_dummy,
                                                                         1, echogram_abs->contrib[ch], 1);

                        /* Store current sample into the circular buffer */
                        sc->circ_buffer[k][src_idx][band][wIdx_n] = sc->src_sigs_bands[src_idx][band][n];
                        if(sc->applyCrossFadeFLAG[rec_idx][src_idx]==0)
                            sc->circ_buffer[IMS_EG_PREV][src_idx][band][wIdx_n] = sc->src_sigs_bands[src_idx][band][n];
                    }

                    /* Increment write index */
                    sc->wIdx[k][rec_idx][src_idx]++;
                    if(sc->applyCrossFadeFLAG[rec_idx][src_idx]==0)
                        sc->wIdx[IMS_EG_PREV][rec_idx][src_idx]++;

                } /* Loop over samples */
            } /* Loop over slots */

            /* Cross-fade between the buffers rendered using the previous and current echograms */
            if(sc->applyCrossFadeFLAG[rec_idx][src_idx]){
                for(ch=0; ch<sc->recs[rec_idx].nChannels; ch++){
                    /* Apply linear interpolator to fade in with the new echogram and fade out with the previous echogram */
                    utility_svvmul(sc->rec_sig_tmp[IMS_EG_CURRENT][rec_idx][ch], sc->interpolator_fIn,  nSamples, sc->tmp_frame);
                    cblas_saxpy(nSamples, 1.0f, sc->tmp_frame, 1, sc->recs[rec_idx].sigs[ch], 1);
                    utility_svvmul(sc->rec_sig_tmp[IMS_EG_PREV][rec_idx][ch],    sc->interpolator_fOut, nSamples, sc->tmp_frame);
                    cblas_saxpy(nSamples, 1.0f, sc->tmp_frame, 1, sc->recs[rec_idx].sigs[ch], 1); 
                }

                /* No longer need to cross-fade for future frames (unless the echograms change again that is...) */
                sc->applyCrossFadeFLAG[rec_idx][src_idx] = 0;
            }
            else
                for(ch=0; ch<sc->recs[rec_idx].nChannels; ch++)
                    cblas_saxpy(nSamples, 1.0f, sc->rec_sig_tmp[IMS_EG_CURRENT][rec_idx][ch], 1, sc->recs[rec_idx].sigs[ch], 1);

        } /* If source active */

    } /* Loop over sources */
}


/* set/get functions: */

void ims_shoebox_setRoomDimensions
(
    void* hIms,
    float new_roomDimensions[3]
)
{
    ims_scene_data *sc = (ims_scene_data*)(hIms);
    int rec_idx, src_idx;

    /* Only update if room dimensions are different */
    if( (sc->room_dims[0]!=new_roomDimensions[0]) ||
        (sc->room_dims[1]!=new_roomDimensions[1]) ||
        (sc->room_dims[2]!=new_roomDimensions[2]) )
    {
        sc->room_dims[0] = new_roomDimensions[0];
        sc->room_dims[1] = new_roomDimensions[1];
        sc->room_dims[2] = new_roomDimensions[2];

        /* Echograms must be re-initialised */
        for(rec_idx = 0; rec_idx < IMS_MAX_NUM_RECEIVERS; rec_idx++)
            for(src_idx = 0; src_idx < IMS_MAX_NUM_SOURCES; src_idx++)
                if( (sc->srcs[src_idx].ID != IMS_UNASSIGNED) && (sc->recs[rec_idx].ID != IMS_UNASSIGNED) )
                    ((ims_core_workspace*)(sc->hCoreWrkSpc[rec_idx][src_idx]))->refreshEchogramFLAG = 1;
    }
}

void ims_shoebox_setWallAbsCoeffs
(
    void* hIms,
    float* abs_wall
)
{
    ims_scene_data *sc = (ims_scene_data*)(hIms);
    int band, i, updateRequired, rec_idx, src_idx;
    updateRequired = 0;

    /* Only update if wall absorption coefficients are different */
    for(band=0; band<sc->nBands; band++){
        for(i=0; i<6; i++){
            if(sc->abs_wall[band][i] != abs_wall[band*6 + i]){
                sc->abs_wall[band][i] = abs_wall[band*6 + i];
                updateRequired = 1;
            }
        }
    }
    if(updateRequired){
        /* Echograms must be re-initialised */
        for(rec_idx = 0; rec_idx < IMS_MAX_NUM_RECEIVERS; rec_idx++)
            for(src_idx = 0; src_idx < IMS_MAX_NUM_SOURCES; src_idx++)
                if( (sc->srcs[src_idx].ID != IMS_UNASSIGNED) && (sc->recs[rec_idx].ID != IMS_UNASSIGNED) )
                    ((ims_core_workspace*)(sc->hCoreWrkSpc[rec_idx][src_idx]))->refreshEchogramFLAG = 1;
    }
}


/* add/remove/update functions: */

int ims_shoebox_addSource
(
    void* hIms,
    float src_xyz[3],
    float** pSrc_sig
)
{
    ims_scene_data *sc = (ims_scene_data*)(hIms);
    int i, rec, obj_idx;

    /* Increment number of sources */
    sc->nSources++;
    saf_assert(sc->nSources <= IMS_MAX_NUM_SOURCES, "Exceeded the maximum supported number of sources");

    /* Find an unoccupied object */
    obj_idx = IMS_UNASSIGNED;
    for(i=0; i<IMS_MAX_NUM_SOURCES; i++){
        /* an ID of '-1' indicates that it is free to use */
        if(sc->srcs[i].ID == IMS_UNASSIGNED){
            obj_idx = i;
            break;
        }
    }
    saf_assert(obj_idx != IMS_UNASSIGNED, "Ugly error");

    /* Assign unique ID */
    sc->srcs[obj_idx].ID = 0;
    for(i=0; i<IMS_MAX_NUM_SOURCES; i++)
        if(i!=obj_idx)
            if(sc->srcs[i].ID == sc->srcs[obj_idx].ID)
                sc->srcs[obj_idx].ID++; /* increment if ID is in use */

    //CHECK
    for(i=0; i<IMS_MAX_NUM_SOURCES; i++)
        if(i!=obj_idx)
            saf_assert(sc->srcs[obj_idx].ID != sc->srcs[i].ID, "Ugly error");

    /* Set source starting position and signal pointer */
    sc->srcs[obj_idx].pos.x = src_xyz[0];
    sc->srcs[obj_idx].pos.y = src_xyz[1];
    sc->srcs[obj_idx].pos.z = src_xyz[2]; 
    sc->srcs[obj_idx].sig = pSrc_sig == NULL ? NULL : *pSrc_sig;

    /* Create workspace for all receiver/source combinations, for this new source object */
    for(rec=0; rec<IMS_MAX_NUM_RECEIVERS; rec++)
        if(sc->recs[rec].ID!=IMS_UNASSIGNED)
            ims_shoebox_coreWorkspaceCreate(&(sc->hCoreWrkSpc[rec][obj_idx]), sc->nBands);

    return sc->srcs[obj_idx].ID;
}

int ims_shoebox_addReceiverSH
(
    void* hIms,
    int sh_order,
    float rec_xyz[3],
    float*** pSH_sigs
)
{
    ims_scene_data *sc = (ims_scene_data*)(hIms);
    int i, src, obj_idx;

    /* Increment number of receivers */
    sc->nReceivers++;
    saf_assert(sc->nReceivers <= IMS_MAX_NUM_RECEIVERS, "Exceeded the maximum supported number of receivers");

    /* Find an unoccupied object */
    obj_idx = IMS_UNASSIGNED;
    for(i=0; i<IMS_MAX_NUM_RECEIVERS; i++){
        /* an ID of '-1' indicates that it is free to use */
        if(sc->recs[i].ID == IMS_UNASSIGNED){
            obj_idx = i;
            break;
        }
    }
    saf_assert(obj_idx != IMS_UNASSIGNED, "Ugly error");

    /* Assign unique ID */
    sc->recs[obj_idx].ID = 0;
    for(i=0; i<IMS_MAX_NUM_RECEIVERS; i++)
        if(i!=obj_idx)
            if(sc->recs[i].ID == sc->recs[obj_idx].ID)
                sc->recs[obj_idx].ID++; /* increment if ID is in use */

    //CHECK
    for(i=0; i<IMS_MAX_NUM_RECEIVERS; i++)
        if(i!=obj_idx)
            saf_assert(sc->recs[obj_idx].ID != sc->recs[i].ID, "Ugly error");

    /* Set starting position, signal pointers, and indicate that this object is
     * a spherical harmonic receiver of order: "sh_order" */
    sc->recs[obj_idx].pos.x = rec_xyz[0];
    sc->recs[obj_idx].pos.y = rec_xyz[1];
    sc->recs[obj_idx].pos.z = rec_xyz[2];
    sc->recs[obj_idx].sigs = pSH_sigs == NULL ? NULL : *pSH_sigs;
    sc->recs[obj_idx].type = RECEIVER_SH;
    sc->recs[obj_idx].nChannels = ORDER2NSH(sh_order);

    /* Create workspace for all receiver/source combinations, for this new receiver object */
    for(src=0; src<IMS_MAX_NUM_SOURCES; src++)
        if(sc->srcs[src].ID!=IMS_UNASSIGNED)
            ims_shoebox_coreWorkspaceCreate(&(sc->hCoreWrkSpc[obj_idx][src]), sc->nBands);

    return sc->recs[obj_idx].ID;
}

void ims_shoebox_updateSource
(
    void* hIms,
    int sourceID,
    float new_position_xyz[3]
)
{
    ims_scene_data *sc = (ims_scene_data*)(hIms);
    ims_core_workspace* work;
    int i, rec, src_idx;

    saf_assert(sourceID >= 0, "Invalid sourceID");

    /* Find index corresponding to this source ID */
    src_idx = IMS_UNASSIGNED;
    for(i=0; i<IMS_MAX_NUM_SOURCES; i++){
        if(sc->srcs[i].ID == sourceID){
            src_idx = i;
            break;
        }
    }
    saf_assert(src_idx != IMS_UNASSIGNED, "Invalid sourceID");

    /* Check if source has actually moved */
    if( (new_position_xyz[0] != sc->srcs[src_idx].pos.x) ||
        (new_position_xyz[1] != sc->srcs[src_idx].pos.y) ||
        (new_position_xyz[2] != sc->srcs[src_idx].pos.z))
    {
       /* update source position */
        sc->srcs[src_idx].pos.x = new_position_xyz[0];
        sc->srcs[src_idx].pos.y = new_position_xyz[1];
        sc->srcs[src_idx].pos.z = new_position_xyz[2];

        /* All source/receiver combinations for this source index will need to
         * be refreshed */
        for(rec=0; rec<IMS_MAX_NUM_RECEIVERS; rec++){
            if(sc->recs[rec].ID!=IMS_UNASSIGNED){
                work = (ims_core_workspace*)(sc->hCoreWrkSpc[rec][src_idx]);
                work->refreshEchogramFLAG = 1;
            }
        }
    }
}

void ims_shoebox_updateReceiver
(
    void* hIms,
    int receiverID,
    float new_position_xyz[3]
)
{
    ims_scene_data *sc = (ims_scene_data*)(hIms);
    ims_core_workspace* work;
    int i, src, rec_idx;

    saf_assert(receiverID >= 0, "Invalid receiverID");

    /* Find index corresponding to this receiver ID */
    rec_idx = IMS_UNASSIGNED;
    for(i=0; i<IMS_MAX_NUM_RECEIVERS; i++){
        if(sc->recs[i].ID == receiverID){
            rec_idx = i;
            break;
        }
    }
    saf_assert(rec_idx != IMS_UNASSIGNED, "Invalid receiverID");

    /* Check if Receiver has actually moved */
    if( (new_position_xyz[0] != sc->recs[rec_idx].pos.x) ||
        (new_position_xyz[1] != sc->recs[rec_idx].pos.y) ||
        (new_position_xyz[2] != sc->recs[rec_idx].pos.z))
    {
        /* update receiver position */
        sc->recs[rec_idx].pos.x = new_position_xyz[0];
        sc->recs[rec_idx].pos.y = new_position_xyz[1];
        sc->recs[rec_idx].pos.z = new_position_xyz[2];

        /* All source/receiver combinations for this receiver index will need to
         * be refreshed */
        for(src=0; src<IMS_MAX_NUM_SOURCES; src++){
            if(sc->srcs[src].ID != IMS_UNASSIGNED){
                work = (ims_core_workspace*)(sc->hCoreWrkSpc[rec_idx][src]);
                work->refreshEchogramFLAG = 1;
            }
        }
    }
}

void ims_shoebox_removeSource
(
    void* hIms,
    int sourceID
)
{
    ims_scene_data *sc = (ims_scene_data*)(hIms);
    int i, obj_idx, rec;

    saf_assert(sourceID >= 0, "Invalid sourceID");

    /* Find index corresponding to this source ID */
    obj_idx = IMS_UNASSIGNED;
    for(i=0; i<IMS_MAX_NUM_SOURCES; i++){
        if(sc->srcs[i].ID == sourceID){
            obj_idx = i;
            break;
        }
    }
    saf_assert(obj_idx != IMS_UNASSIGNED, "Invalid sourceID");

    /* Set ID to -1 (invalid, so no longer rendered) */
    sc->srcs[obj_idx].ID = IMS_UNASSIGNED;

    /* Destroy workspace for all receiver/source combinations, for this dead source */
    for(rec=0; rec<IMS_MAX_NUM_RECEIVERS; rec++)
        if(sc->recs[rec].ID != IMS_UNASSIGNED)
            ims_shoebox_coreWorkspaceDestroy(&(sc->hCoreWrkSpc[rec][obj_idx]));

    /* De-increment number of sources */
    sc->nSources--;
}

void ims_shoebox_removeReceiver
(
    void* hIms,
    int receiverID
)
{
    ims_scene_data *sc = (ims_scene_data*)(hIms);
    int i, obj_idx, src;

    saf_assert(receiverID >= 0, "Invalid receiverID");

    /* Find index corresponding to this source ID */
    obj_idx = IMS_UNASSIGNED;
    for(i=0; i<IMS_MAX_NUM_RECEIVERS; i++){
        if(sc->recs[i].ID == receiverID){
            obj_idx = i;
            break;
        }
    }
    saf_assert(obj_idx != IMS_UNASSIGNED, "Invalid receiverID");

    /* Set ID to -1 (invalid, so no longer active) */
    sc->recs[obj_idx].ID = IMS_UNASSIGNED;

    /* Destroy workspace for all receiver/source combinations, for this dead receiver */
    for(src=0; src<IMS_MAX_NUM_SOURCES; src++)
        if(sc->srcs[src].ID != IMS_UNASSIGNED)
            ims_shoebox_coreWorkspaceDestroy(&(sc->hCoreWrkSpc[obj_idx][src]));

    /* De-increment number of receivers */
    sc->nReceivers--;
}
