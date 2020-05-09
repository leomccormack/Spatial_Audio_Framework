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
 * @brief Public part of the reverb processing module (saf_reverb)
 *
 * ...
 *
 * @author Leo McCormack
 * @date 06.05.2020
 */
 
#include "saf_reverb.h"
#include "saf_reverb_internal.h"

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
    ims_scene_data *h = (ims_scene_data*)(*phIms);
    int i,j,band,wall;

    assert(nOctBands>1);

    /* Shoebox dimensions */
    h->room_dimensions[0] = length;
    h->room_dimensions[1] = width;
    h->room_dimensions[2] = height;
    h->c_ms = c_ms;

    /* Octave band centre frequencies */
    h->nBands = nOctBands;
    h->band_centerfreqs = malloc1d(nOctBands*sizeof(float));
    h->band_centerfreqs[0] = lowestOctaveBand;
    for(band=1; band<nOctBands; band++)
        h->band_centerfreqs[band] = h->band_centerfreqs[band-1];
    h->fs = fs;

    /* Absorption coeffients per wall and octave band */
    h->abs_wall = (float**)malloc2d(nOctBands, IMS_NUM_WALLS_SHOEBOX, sizeof(float));
    for(band=0; band<nOctBands; band++)
        for(wall=0; wall<IMS_NUM_WALLS_SHOEBOX; wall++)
            h->abs_wall[band][wall] = abs_wall[band*IMS_NUM_WALLS_SHOEBOX+wall];

    /* Default is no sources or receivers in the room */
    for(i=0; i<IMS_MAX_NUM_SOURCES; i++)
        h->src_IDs[i] = -1; /* -1 indicates not in use */
    for(i=0; i<IMS_MAX_NUM_RECEIVERS; i++)
        h->rec_IDs[i] = -1; 
    h->nSources = 0;
    h->nReceivers = 0;

    /* ims_core_workspace per source / receiver combination */
    h->hCoreWrkSpc = (voidPtr**)malloc2d(IMS_MAX_NUM_RECEIVERS, IMS_MAX_NUM_SOURCES, sizeof(voidPtr));
    for(i=0; i<IMS_MAX_NUM_RECEIVERS; i++)
        for(j=0; j<IMS_MAX_NUM_SOURCES; j++)
            h->hCoreWrkSpc[i][j] = NULL;
}

void ims_shoebox_destroy
(
    void** phIms
)
{
    ims_scene_data *h = (ims_scene_data*)(*phIms);
    int i,j;

    if(h!=NULL){
        free(h->band_centerfreqs);
        free(h->abs_wall);
        for(i=0; i<IMS_MAX_NUM_RECEIVERS; i++)
            for(j=0; j<IMS_MAX_NUM_SOURCES; j++)
                ims_shoebox_coreWorkspaceDestroy(&(h->hCoreWrkSpc[i][j]));
        free(h->hCoreWrkSpc);
        free(h);
        h=NULL;
        *phIms = NULL;
    }
}

void ims_shoebox_renderEchogramSH
(
    void* hIms,
    float maxTime_ms,
    int sh_order
)
{
    ims_scene_data *h = (ims_scene_data*)(hIms);
    ims_core_workspace* workspace;
    position_xyz src2, rec2;
    int src_idx, rec_idx;

    /* Compute echograms for active source/receiver combinations */
    for(rec_idx = 0; rec_idx < IMS_MAX_NUM_RECEIVERS; rec_idx++){
        for(src_idx = 0; src_idx < IMS_MAX_NUM_SOURCES; src_idx++){
            if( (h->src_IDs[src_idx]!=-1) && (h->rec_IDs[rec_idx]!=-1) ){
                /* Change y coord for Source to match convention used inside the coreInit function */
                rec2.x = h->rec_xyz[rec_idx].x;
                rec2.y = (float)h->room_dimensions[1] - h->rec_xyz[rec_idx].y;
                rec2.z = h->rec_xyz[rec_idx].z;
                src2.x = h->src_xyz[src_idx].x;
                src2.y = (float)h->room_dimensions[1] - h->src_xyz[src_idx].y;
                src2.z = h->src_xyz[src_idx].z;

                /* Workspace handle for this source/receiver combination */
                workspace = h->hCoreWrkSpc[rec_idx][src_idx];

                /* Only update if it is required */
                if(workspace->refreshEchogramFLAG){
                    /* Compute echogram due to pure propagation (frequency-independent) */
                    ims_shoebox_coreInit(workspace,
                                         h->room_dimensions, src2, rec2, maxTime_ms, h->c_ms);

                    /* Apply spherical harmonic directivities */
                    ims_shoebox_coreRecModuleSH(workspace, sh_order);

                    /* Apply boundary absoption per band */
                    ims_shoebox_coreAbsorptionModule(workspace, h->abs_wall);

                    /* Indicate that the echogram is now up to date */
                    workspace->refreshEchogramFLAG = 0;
                }
            }
        }
    }  
}

void ims_shoebox_renderSHRIRs
(
    void* hIms,
    int fractionalDelayFLAG
)
{
    ims_scene_data *h = (ims_scene_data*)(hIms);
    int src_idx, rec_idx;

    for(rec_idx = 0; rec_idx < IMS_MAX_NUM_RECEIVERS; rec_idx++){
        for(src_idx = 0; src_idx < IMS_MAX_NUM_SOURCES; src_idx++){
            if( (h->src_IDs[src_idx]!=-1) && (h->rec_IDs[rec_idx]!=-1) ){

            }
        }
    }

}



/* add/remove/update functions: */

long ims_shoebox_addSource
(
    void* hIms,
    float src_xyz[3]
)
{
    ims_scene_data *h = (ims_scene_data*)(hIms);
    int i, rec, slot_idx;

    /* Increment number of sources */
    h->nSources++;
    assert(h->nSources <= IMS_MAX_NUM_SOURCES);

    /* Find an unoccupied slot */
    slot_idx = -1;
    for(i=0; i<IMS_MAX_NUM_SOURCES; i++){
        /* an ID of '-1' indicates that it is free to use */
        if(h->src_IDs[i] == -1){
            slot_idx = i;
            break;
        }
    }
    assert(slot_idx != -1);

    /* Assign unique ID */
    h->src_IDs[slot_idx] = 0;
    for(i=0; i<IMS_MAX_NUM_SOURCES; i++)
        if(i!=slot_idx)
            if(h->src_IDs[i]==h->src_IDs[slot_idx])
                h->src_IDs[slot_idx]++; /* increment if ID is in use */

    //CHECK
    for(i=0; i<IMS_MAX_NUM_SOURCES; i++)
        if(i!=slot_idx)
            assert(h->src_IDs[slot_idx]!=h->src_IDs[i]);

    /* Set source starting position */
    for(i=0; i<3; i++)
        h->src_xyz[slot_idx].v[i] = src_xyz[i];

    /* Create workspace for all receiver/source combinations, for this new source slot */
    for(rec=0; rec<IMS_MAX_NUM_RECEIVERS; rec++)
        if(h->rec_IDs[rec]!=-1)
            ims_shoebox_coreWorkspaceCreate(&(h->hCoreWrkSpc[rec][slot_idx]), h->nBands);


    return h->src_IDs[slot_idx];
}

long ims_shoebox_addReceiver
(
    void* hIms,
    float rec_xyz[3]
)
{
    ims_scene_data *h = (ims_scene_data*)(hIms);
    int i, src, slot_idx; 

    /* Increment number of receivers */
    h->nReceivers++;
    assert(h->nReceivers <= IMS_MAX_NUM_RECEIVERS);

    /* Find an unoccupied slot */
    slot_idx = -1;
    for(i=0; i<IMS_MAX_NUM_RECEIVERS; i++){
        /* an ID of '-1' indicates that it is free to use */
        if(h->rec_IDs[i] == -1){
            slot_idx = i;
            break;
        }
    }
    assert(slot_idx != -1);

    /* Assign unique ID */
    h->rec_IDs[slot_idx] = 0;
    for(i=0; i<IMS_MAX_NUM_RECEIVERS; i++)
        if(i!=slot_idx)
            if(h->rec_IDs[i]==h->rec_IDs[slot_idx])
                h->rec_IDs[slot_idx]++; /* increment if ID is in use */

    //CHECK
    for(i=0; i<IMS_MAX_NUM_RECEIVERS; i++)
        if(i!=slot_idx)
            assert(h->rec_IDs[slot_idx]!=h->rec_IDs[i]);

    /* Set starting position */
    for(i=0; i<3; i++)
        h->rec_xyz[slot_idx].v[i] = rec_xyz[i];

    /* Create workspace for all receiver/source combinations, for this new receiver slot */
    for(src=0; src<IMS_MAX_NUM_SOURCES; src++)
        if(h->src_IDs[src]!=-1)
            ims_shoebox_coreWorkspaceCreate(&(h->hCoreWrkSpc[slot_idx][src]), h->nBands);

    return h->rec_IDs[slot_idx];
}

void ims_shoebox_updateSource
(
    void* hIms,
    long sourceID,
    float new_position_xyz[3]
)
{
    ims_scene_data *h = (ims_scene_data*)(hIms);
    ims_core_workspace* work;
    int i, rec, src_idx;

    assert(sourceID >= 0);

    /* Find index corresponding to this source ID */
    src_idx = -1;
    for(i=0; i<IMS_MAX_NUM_SOURCES; i++){
        if(h->src_IDs[i] == sourceID){
            src_idx = i;
            break;
        }
    }
    assert(src_idx != -1);

    /* update source position */
    for(i=0; i<3; i++)
        h->src_xyz[src_idx].v[i] = new_position_xyz[i];

    /* All source/receiver combinations for this source index will require refreshing */
    for(rec=0; rec<IMS_MAX_NUM_RECEIVERS; rec++){
        if(h->rec_IDs[rec]!=-1){
            work = (ims_core_workspace*)(h->hCoreWrkSpc[rec][src_idx]);
            work->refreshEchogramFLAG = 1;
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
    ims_scene_data *h = (ims_scene_data*)(hIms);
    ims_core_workspace* work;
    int i, src, rec_idx;

    assert(receiverID >= 0);

    /* Find index corresponding to this source ID */
    rec_idx = -1;
    for(i=0; i<IMS_MAX_NUM_RECEIVERS; i++){
        if(h->rec_IDs[i] == receiverID){
            rec_idx = i;
            break;
        }
    }
    assert(rec_idx != -1);

    /* update position */
    for(i=0; i<3; i++)
        h->rec_xyz[rec_idx].v[i] = new_position_xyz[i];

    /* All source/receiver combinations for this receiver index will require refreshing */
    for(src=0; src<IMS_MAX_NUM_SOURCES; src++){
        if(h->src_IDs[src]!=-1){
            work = (ims_core_workspace*)(h->hCoreWrkSpc[rec_idx][src]);
            work->refreshEchogramFLAG = 1;
        }
    }
}

void ims_shoebox_removeSource
(
    void* hIms,
    long sourceID
)
{
    ims_scene_data *h = (ims_scene_data*)(hIms);
    int i, slot_idx, rec;

    assert(sourceID >= 0);

    /* Find index corresponding to this source ID */
    slot_idx = -1;
    for(i=0; i<IMS_MAX_NUM_SOURCES; i++){
        if(h->src_IDs[i] == sourceID){
            slot_idx = i;
            break;
        }
    }
    assert(slot_idx != -1);

    /* Set ID to -1 (invalid, so no longer rendered) */
    h->src_IDs[slot_idx] = -1;

    /* Destroy workspace for all receiver/source combinations, for this dead source */
    for(rec=0; rec<IMS_MAX_NUM_RECEIVERS; rec++)
        if(h->rec_IDs[rec]!=-1)
            ims_shoebox_coreWorkspaceDestroy(&(h->hCoreWrkSpc[rec][slot_idx]));

    /* De-increment number of sources */
    h->nSources--;
}

void ims_shoebox_removeReceiver
(
    void* hIms,
    long receiverID
)
{
    ims_scene_data *h = (ims_scene_data*)(hIms);
    int i, slot_idx, src;

    assert(receiverID >= 0);

    /* Find index corresponding to this source ID */
    slot_idx = -1;
    for(i=0; i<IMS_MAX_NUM_RECEIVERS; i++){
        if(h->rec_IDs[i] == receiverID){
            slot_idx = i;
            break;
        }
    }
    assert(slot_idx != -1);

    /* Set ID to -1 (invalid, so no longer active) */
    h->rec_IDs[slot_idx] = -1;

    /* Destroy workspace for all receiver/source combinations, for this dead receiver */
    for(src=0; src<IMS_MAX_NUM_SOURCES; src++)
        if(h->src_IDs[src]!=-1)
            ims_shoebox_coreWorkspaceDestroy(&(h->hCoreWrkSpc[slot_idx][src]));

    /* De-increment number of receivers */
    h->nReceivers--;
}

