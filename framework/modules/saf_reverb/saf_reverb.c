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
 * @brief Public source for the reverb processing module (#SAF_REVERB_MODULE)
 *
 * A collection of reverb and room simulation algorithms.
 *
 * @author Leo McCormack
 * @date 06.05.2020
 */
 
#include "saf_reverb.h"
#include "saf_reverb_internal.h"

/* ========================================================================== */
/*                         IMS Shoebox Room Simulator                         */
/* ========================================================================== */

void ims_shoebox_create
(
    void** phIms,
    int length,
    int width,
    int height,
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

    assert(nOctBands>1);

    /* Shoebox dimensions */
    sc->room_dimensions[0] = length;
    sc->room_dimensions[1] = width;
    sc->room_dimensions[2] = height;
    sc->c_ms = c_ms;

    /* Octave band centre frequencies */
    sc->nBands = nOctBands;
    sc->band_centerfreqs = malloc1d(nOctBands*sizeof(float));
    sc->band_centerfreqs[0] = lowestOctaveBand;
    for(band=1; band<nOctBands; band++)
        sc->band_centerfreqs[band] = sc->band_centerfreqs[band-1]*2.0f;
    sc->band_cutofffreqs = malloc1d((sc->nBands-1)*sizeof(float));
    getOctaveBandCutoffFreqs(sc->band_centerfreqs, sc->nBands, sc->band_cutofffreqs);
    sc->fs = fs;

    /* Absorption coeffients per wall and octave band */
    sc->abs_wall = (float**)malloc2d(nOctBands, IMS_NUM_WALLS_SHOEBOX, sizeof(float));
    for(band=0; band<nOctBands; band++)
        for(wall=0; wall<IMS_NUM_WALLS_SHOEBOX; wall++)
            sc->abs_wall[band][wall] = abs_wall[band*IMS_NUM_WALLS_SHOEBOX+wall];

    /* Default is no sources or receivers in the room */
    for(i=0; i<IMS_MAX_NUM_SOURCES; i++)
        sc->srcs[i].ID = -1; /* -1 indicates not in use */
    for(i=0; i<IMS_MAX_NUM_RECEIVERS; i++)
        sc->recs[i].ID = -1;
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

    /* Circular buffer */
    sc->wIdx = 0;
    sc->circ_buffer = NULL;

    /* IIR Filterbank per source */
    sc->hFaFbank = malloc1d(IMS_MAX_NUM_SOURCES*sizeof(voidPtr));
    sc->src_sigs_bands = malloc1d(IMS_MAX_NUM_SOURCES*sizeof(float**));
    for(j=0; j<IMS_MAX_NUM_SOURCES; j++){
        sc->hFaFbank[j] = NULL;
        sc->src_sigs_bands[j] = NULL;
    }
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
        free(sc->circ_buffer);
        for(j=0; j<IMS_MAX_NUM_SOURCES; j++){
            faf_IIRFilterbank_destroy(&(sc->hFaFbank[j]));
            free(sc->src_sigs_bands[j]);
        }
        free(sc->hFaFbank);
        free(sc->src_sigs_bands);
        free(sc);
        sc=NULL;
        *phIms = NULL;
    }
}

void ims_shoebox_computeEchograms
(
    void* hIms,
    float maxTime_ms
)
{
    ims_scene_data *sc = (ims_scene_data*)(hIms);
    ims_core_workspace* workspace;
    ims_pos_xyz src2, rec2;
    int src_idx, rec_idx;

    /* Compute echograms for active source/receiver combinations */
    for(rec_idx = 0; rec_idx < IMS_MAX_NUM_RECEIVERS; rec_idx++){
        for(src_idx = 0; src_idx < IMS_MAX_NUM_SOURCES; src_idx++){
            if( (sc->srcs[src_idx].ID != -1) && (sc->recs[rec_idx].ID != -1) ){
                /* Change y coord for Receiver and Source to match convention
                 * used inside the coreInit function */
                rec2.x = sc->recs[rec_idx].pos.x;
                rec2.y = (float)sc->room_dimensions[1] - sc->recs[rec_idx].pos.y;
                rec2.z = sc->recs[rec_idx].pos.z;
                src2.x = sc->srcs[src_idx].pos.x;
                src2.y = (float)sc->room_dimensions[1] - sc->srcs[src_idx].pos.y;
                src2.z = sc->srcs[src_idx].pos.z;

                /* Workspace handle for this source/receiver combination */
                workspace = sc->hCoreWrkSpc[rec_idx][src_idx];

                /* Only update if it is required */
                if(workspace->refreshEchogramFLAG){
                    /* Compute echogram due to pure propagation (frequency-independent) */
                    ims_shoebox_coreInit(workspace,
                                         sc->room_dimensions, src2, rec2, maxTime_ms, sc->c_ms);

                    /* Apply receiver directivities */
                    switch(sc->recs[rec_idx].type){
                        case RECEIVER_SH:
                            ims_shoebox_coreRecModuleSH(workspace, NSH2ORDER(sc->recs[rec_idx].nChannels));
                            break;
                    }

                    /* Apply boundary absoption per frequency band */
                    ims_shoebox_coreAbsorptionModule(workspace, sc->abs_wall);

                    /* Indicate that the echogram is now up to date, and that the RIR should now be updated */
                    workspace->refreshEchogramFLAG = 0;
                    workspace->refreshRIRFLAG = 1;
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
            if( (sc->srcs[src_idx].ID!=-1) && (sc->recs[rec_idx].ID!=-1) ){

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
    int fractionalDelaysFLAG // NOT IMPLEMENTED YET
)
{
    ims_scene_data *sc = (ims_scene_data*)(hIms);
    ims_core_workspace* wrk;
    echogram_data *echogram_abs;
    int i, n, im, band, ch, rec_idx, src_idx, time_samples, wIdx_n;
    float cb_val;
    unsigned int rIdx;

    assert(nSamples <= IMS_MAX_NSAMPLES_PER_FRAME);


    /* Allocate circular buffers and filterbank handles (if this is the first
     * time this function is being called) */
    if(sc->circ_buffer==NULL)
        sc->circ_buffer = (float***)calloc3d(IMS_MAX_NUM_SOURCES, sc->nBands, IMS_CIRC_BUFFER_LENGTH, sizeof(float));
    for(src_idx = 0; src_idx < IMS_MAX_NUM_SOURCES; src_idx++){
        /* Also, only allocate if source is active: */
        if( (sc->srcs[src_idx].ID != -1) && (sc->hFaFbank[src_idx] == NULL) ){
            faf_IIRFilterbank_create(&(sc->hFaFbank[src_idx]), IMS_IIR_FILTERBANK_ORDER, sc->band_cutofffreqs,
                                     sc->nBands-1, sc->fs, IMS_MAX_NSAMPLES_PER_FRAME);
            sc->src_sigs_bands[src_idx] = (float**)malloc2d(sc->nBands, IMS_MAX_NSAMPLES_PER_FRAME, sizeof(float));
        }
    }

    /* Find index corresponding to this receiver ID */
    rec_idx = -1;
    for(i=0; i<IMS_MAX_NUM_RECEIVERS; i++){
        if(sc->recs[i].ID == receiverID){
            rec_idx = i;
            break;
        }
    }
    assert(rec_idx != -1);

    /* Initialise buffer for receiver with zeros */
    memset(FLATTEN2D(sc->recs[rec_idx].sigs), sc->recs[rec_idx].nChannels * nSamples, sizeof(float));
 
    /* Process all active sources (for this specific receiver) directly in the
     * time-domain */
    for(src_idx = 0; src_idx < IMS_MAX_NUM_SOURCES; src_idx++){
        if( (sc->srcs[src_idx].ID!=-1) && (sc->recs[rec_idx].ID!=-1) ){
            /* Pass source signal through the Favrot & Faller IIR filterbank */
            faf_IIRFilterbank_apply(sc->hFaFbank[src_idx], sc->srcs[src_idx].sig, sc->src_sigs_bands[src_idx], nSamples);

            /* Workspace handle for this source/receiver combination */
            wrk = sc->hCoreWrkSpc[rec_idx][src_idx];

            /* Loop over samples */
            for(n=0; n<nSamples; n++){
                /* Number of image sources and time indices are the same across
                 * octave bands: */
                echogram_abs = (echogram_data*)wrk->hEchogram_abs[0];

                /* Determine write index */
                wIdx_n = sc->wIdx & IMS_CIRC_BUFFER_LENGTH_MASK;

                /* Loop over all image sources */
                for(im=0; im <echogram_abs->numImageSources; im++){
                    /* Determine read index */
                    time_samples = (int)(echogram_abs->time[im] * (sc->fs) + 0.5f); /* Round to nearest sample */
                    rIdx = IMS_CIRC_BUFFER_LENGTH-time_samples + sc->wIdx;          /* read index for this image source */
                    rIdx = rIdx & IMS_CIRC_BUFFER_LENGTH_MASK;                      /* wrap-around if needed */

                    /* Loop over octave bands */
                    for(band=0; band < sc->nBands; band++){
                        /* Echogram for this source/receiver at this band */
                        echogram_abs = (echogram_data*)wrk->hEchogram_abs[band];

                        /* Pull value from circular buffer at this read index */
                        cb_val = sc->circ_buffer[src_idx][band][rIdx]; // 5.3% of CPU TIME

                        /* For compiler optimisations up to 3rd order SH receiver*/
                        switch(sc->recs[rec_idx].type){
                            case RECEIVER_SH:
                                for(ch=16; ch<sc->recs[rec_idx].nChannels; ch++)  
                                    sc->recs[rec_idx].sigs[ch][n] += echogram_abs->value[im][ch] * cb_val;

                                switch(sc->recs[rec_idx].nChannels){
                                    case 64: case 49: case 36: case 25: case 16:
                                        sc->recs[rec_idx].sigs[15][n] += echogram_abs->value[im][15] * cb_val;
                                        sc->recs[rec_idx].sigs[14][n] += echogram_abs->value[im][14] * cb_val;
                                        sc->recs[rec_idx].sigs[13][n] += echogram_abs->value[im][13] * cb_val;
                                        sc->recs[rec_idx].sigs[12][n] += echogram_abs->value[im][12] * cb_val;
                                        sc->recs[rec_idx].sigs[11][n] += echogram_abs->value[im][11] * cb_val;
                                        sc->recs[rec_idx].sigs[10][n] += echogram_abs->value[im][10] * cb_val;
                                        sc->recs[rec_idx].sigs[9][n] += echogram_abs->value[im][9] * cb_val;
                                    case 9:
                                        sc->recs[rec_idx].sigs[8][n] += echogram_abs->value[im][8] * cb_val;
                                        sc->recs[rec_idx].sigs[7][n] += echogram_abs->value[im][7] * cb_val;
                                        sc->recs[rec_idx].sigs[6][n] += echogram_abs->value[im][6] * cb_val;
                                        sc->recs[rec_idx].sigs[5][n] += echogram_abs->value[im][5] * cb_val;
                                        sc->recs[rec_idx].sigs[4][n] += echogram_abs->value[im][4] * cb_val;
                                    case 4:
                                        sc->recs[rec_idx].sigs[3][n] += echogram_abs->value[im][3] * cb_val;
                                        sc->recs[rec_idx].sigs[2][n] += echogram_abs->value[im][2] * cb_val;
                                        sc->recs[rec_idx].sigs[1][n] += echogram_abs->value[im][1] * cb_val;
                                    case 1:
                                        sc->recs[rec_idx].sigs[0][n] += echogram_abs->value[im][0] * cb_val;
                                }
                                break;
                        }

                        /* Copy input to circular buffer */
                        sc->circ_buffer[src_idx][band][wIdx_n] = sc->src_sigs_bands[src_idx][band][n];
                    }
                }
 
                /* Increment write index */
                sc->wIdx++;
            }
        }
    }
}


/* add/remove/update functions: */

long ims_shoebox_addSource
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
    assert(sc->nSources <= IMS_MAX_NUM_SOURCES);

    /* Find an unoccupied object */
    obj_idx = -1;
    for(i=0; i<IMS_MAX_NUM_SOURCES; i++){
        /* an ID of '-1' indicates that it is free to use */
        if(sc->srcs[i].ID == -1){
            obj_idx = i;
            break;
        }
    }
    assert(obj_idx != -1);

    /* Assign unique ID */
    sc->srcs[obj_idx].ID = 0;
    for(i=0; i<IMS_MAX_NUM_SOURCES; i++)
        if(i!=obj_idx)
            if(sc->srcs[i].ID == sc->srcs[obj_idx].ID)
                sc->srcs[obj_idx].ID++; /* increment if ID is in use */

    //CHECK
    for(i=0; i<IMS_MAX_NUM_SOURCES; i++)
        if(i!=obj_idx)
            assert(sc->srcs[obj_idx].ID != sc->srcs[i].ID);

    /* Set source starting position and signal pointer */
    sc->srcs[obj_idx].pos.x = src_xyz[0];
    sc->srcs[obj_idx].pos.y = src_xyz[1];
    sc->srcs[obj_idx].pos.z = src_xyz[2]; 
    sc->srcs[obj_idx].sig = pSrc_sig == NULL ? NULL : *pSrc_sig;

    /* Create workspace for all receiver/source combinations, for this new source object */
    for(rec=0; rec<IMS_MAX_NUM_RECEIVERS; rec++)
        if(sc->recs[rec].ID!=-1)
            ims_shoebox_coreWorkspaceCreate(&(sc->hCoreWrkSpc[rec][obj_idx]), sc->nBands);

    return sc->srcs[obj_idx].ID;
}

long ims_shoebox_addReceiverSH
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
    assert(sc->nReceivers <= IMS_MAX_NUM_RECEIVERS);

    /* Find an unoccupied object */
    obj_idx = -1;
    for(i=0; i<IMS_MAX_NUM_RECEIVERS; i++){
        /* an ID of '-1' indicates that it is free to use */
        if(sc->recs[i].ID == -1){
            obj_idx = i;
            break;
        }
    }
    assert(obj_idx != -1);

    /* Assign unique ID */
    sc->recs[obj_idx].ID = 0;
    for(i=0; i<IMS_MAX_NUM_RECEIVERS; i++)
        if(i!=obj_idx)
            if(sc->recs[i].ID == sc->recs[obj_idx].ID)
                sc->recs[obj_idx].ID++; /* increment if ID is in use */

    //CHECK
    for(i=0; i<IMS_MAX_NUM_RECEIVERS; i++)
        if(i!=obj_idx)
            assert(sc->recs[obj_idx].ID != sc->recs[i].ID);

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
        if(sc->srcs[src].ID!=-1)
            ims_shoebox_coreWorkspaceCreate(&(sc->hCoreWrkSpc[obj_idx][src]), sc->nBands);

    return sc->recs[obj_idx].ID;
}

void ims_shoebox_updateSource
(
    void* hIms,
    long sourceID,
    float new_position_xyz[3]
)
{
    ims_scene_data *sc = (ims_scene_data*)(hIms);
    ims_core_workspace* work;
    int i, rec, src_idx;

    assert(sourceID >= 0);

    /* Find index corresponding to this source ID */
    src_idx = -1;
    for(i=0; i<IMS_MAX_NUM_SOURCES; i++){
        if(sc->srcs[i].ID == sourceID){
            src_idx = i;
            break;
        }
    }
    assert(src_idx != -1);

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
            if(sc->recs[rec].ID!=-1){
                work = (ims_core_workspace*)(sc->hCoreWrkSpc[rec][src_idx]);
                work->refreshEchogramFLAG = 1;
            }
        }
    }
}

void ims_shoebox_updateReceiver
(
    void* hIms,
    long receiverID,
    float new_position_xyz[3]
)
{
    ims_scene_data *sc = (ims_scene_data*)(hIms);
    ims_core_workspace* work;
    int i, src, rec_idx;

    assert(receiverID >= 0);

    /* Find index corresponding to this receiver ID */
    rec_idx = -1;
    for(i=0; i<IMS_MAX_NUM_RECEIVERS; i++){
        if(sc->recs[i].ID == receiverID){
            rec_idx = i;
            break;
        }
    }
    assert(rec_idx != -1);

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
            if(sc->srcs[src].ID != -1){
                work = (ims_core_workspace*)(sc->hCoreWrkSpc[rec_idx][src]);
                work->refreshEchogramFLAG = 1;
            }
        }
    }
}

void ims_shoebox_removeSource
(
    void* hIms,
    long sourceID
)
{
    ims_scene_data *sc = (ims_scene_data*)(hIms);
    int i, obj_idx, rec;

    assert(sourceID >= 0);

    /* Find index corresponding to this source ID */
    obj_idx = -1;
    for(i=0; i<IMS_MAX_NUM_SOURCES; i++){
        if(sc->srcs[i].ID == sourceID){
            obj_idx = i;
            break;
        }
    }
    assert(obj_idx != -1);

    /* Set ID to -1 (invalid, so no longer rendered) */
    sc->srcs[obj_idx].ID = -1;

    /* Destroy workspace for all receiver/source combinations, for this dead source */
    for(rec=0; rec<IMS_MAX_NUM_RECEIVERS; rec++)
        if(sc->recs[rec].ID != -1)
            ims_shoebox_coreWorkspaceDestroy(&(sc->hCoreWrkSpc[rec][obj_idx]));

    /* De-increment number of sources */
    sc->nSources--;
}

void ims_shoebox_removeReceiver
(
    void* hIms,
    long receiverID
)
{
    ims_scene_data *sc = (ims_scene_data*)(hIms);
    int i, obj_idx, src;

    assert(receiverID >= 0);

    /* Find index corresponding to this source ID */
    obj_idx = -1;
    for(i=0; i<IMS_MAX_NUM_RECEIVERS; i++){
        if(sc->recs[i].ID == receiverID){
            obj_idx = i;
            break;
        }
    }
    assert(obj_idx != -1);

    /* Set ID to -1 (invalid, so no longer active) */
    sc->recs[obj_idx].ID = -1;

    /* Destroy workspace for all receiver/source combinations, for this dead receiver */
    for(src=0; src<IMS_MAX_NUM_SOURCES; src++)
        if(sc->srcs[src].ID != -1)
            ims_shoebox_coreWorkspaceDestroy(&(sc->hCoreWrkSpc[obj_idx][src]));

    /* De-increment number of receivers */
    sc->nReceivers--;
}

