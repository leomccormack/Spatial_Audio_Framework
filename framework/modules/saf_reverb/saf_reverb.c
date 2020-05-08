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
        h->src_IDs[i] = -1;
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

    /* Compute echograms */
    for(rec_idx = 0; rec_idx < h->nReceivers; rec_idx++){
        /* Change y coord for receiver to match convention used inside the coreInit function */
        rec2.x = h->rec_xyz[rec_idx].x;
        rec2.y = (float)h->room_dimensions[1] - h->rec_xyz[rec_idx].y;
        rec2.z = h->rec_xyz[rec_idx].z;

        for(src_idx = 0; src_idx < h->nSources; src_idx++){
            /* Change y coord for Source to match convention used inside the coreInit function */
            src2.x = h->src_xyz[src_idx].x;
            src2.y = (float)h->room_dimensions[1] - h->src_xyz[src_idx].y;
            src2.z = h->src_xyz[src_idx].z;

            /* Workspace handle for this source/receiver combination */
            workspace = h->hCoreWrkSpc[rec_idx][src_idx];

            if(workspace->refreshEchogramFLAG){
                /* Compute echogram due to pure propagation (frequency-independent) */
                ims_shoebox_coreInit(workspace,
                                     h->room_dimensions, src2, rec2, maxTime_ms, h->c_ms);

                /* Apply spherical harmonic directivities */
                ims_shoebox_coreRecModuleSH(workspace, sh_order);

                /* Apply boundary absoption per band */
                ims_shoebox_coreAbsorptionModule(workspace, h->abs_wall);

                workspace->refreshEchogramFLAG = 0;
            }
        }
    } 
////TEST
//    int i;
//    float echograms_value[823];
//    float echograms_time[823];
//    reflOrder echograms_order[823];
//    position_xyz echograms_coords[823];
//
//    for(i=0; i<823; i++){
//        echograms_value[i] = h->echograms[0][0]->value[h->echograms[0][0].sortedIdx[i]][0];
//        echograms_time[i] = h->echograms[0][0]->time[h->echograms[0][0].sortedIdx[i]];
//        echograms_order[i] = h->echograms[0][0]->order[h->echograms[0][0].sortedIdx[i]];
//        echograms_coords[i] = h->echograms[0][0].coords[h->echograms[0][0].sortedIdx[i]];
//    }
}

void ims_shoebox_renderSHRIRs
(
    void* hIms,
    int fractionalDelayFLAG
)
{
    ims_scene_data *h = (ims_scene_data*)(hIms);
    int src_idx, rec_idx, band_idx;

    for(rec_idx = 0; rec_idx < h->nReceivers; rec_idx++){
        for(src_idx = 0; src_idx < h->nSources; src_idx++){
            for(band_idx = 0; band_idx < h->nBands; band_idx++){


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
    int i, rec;
    long newID;

    /* Append new source coordinates */
    h->nSources++;
    for(i=0; i<3; i++)
        h->src_xyz[h->nSources-1].v[i] = src_xyz[i];

    /* Assign unique ID */
    newID = 0;
    for(i=0; i<h->nSources-1; i++)
        if(h->src_IDs[i]==newID) /* check ID is not in use */
            newID++;
    h->src_IDs[h->nSources-1] = newID;

    /* Create workspace for all source/receiver combinations, for this new source */
    for(rec=0; rec<h->nReceivers; rec++)
        ims_shoebox_coreWorkspaceCreate(&(h->hCoreWrkSpc[rec][h->nSources-1]), h->nBands);

    return newID;
}

long ims_shoebox_addReceiver
(
    void* hIms,
    float rec_xyz[3]
)
{
    ims_scene_data *h = (ims_scene_data*)(hIms);
    int i, src;
    long newID;

    /* Append new source coordinates */
    h->nReceivers++;
    for(i=0; i<3; i++)
        h->rec_xyz[h->nReceivers-1].v[i] = rec_xyz[i];

    /* Assign unique ID */
    newID = 0;
    for(i=0; i<h->nReceivers-1; i++)
        if(h->rec_IDs[i]==newID) /* check ID is not in use */
            newID++;
    h->rec_IDs[h->nReceivers-1] = newID;

    /* Create workspace for all receiver/source combinations, for this new receiver */ 
    for(src=0; src<h->nSources; src++)
        ims_shoebox_coreWorkspaceCreate(&(h->hCoreWrkSpc[h->nReceivers-1][src]), h->nBands);

    return newID;
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
    int i, src_idx;

    /* Find index corresponding to this source ID */
    src_idx = -1;
    for(i=0; i<h->nSources; i++)
        if(h->src_IDs[i] == sourceID)
            src_idx = i;
    assert(src_idx != -1);

    /* update position */
    for(i=0; i<3; i++)
        h->src_xyz[src_idx].v[i] = new_position_xyz[i];

    /* All source/receiver combinations for this source index will require refreshing */
    for(i=0; i<h->nReceivers; i++){
        work = (ims_core_workspace*)(h->hCoreWrkSpc[i][src_idx]);
        work->refreshEchogramFLAG = 1;
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
    int i, rec_idx;

    /* Find index corresponding to this source ID */
    rec_idx = -1;
    for(i=0; i<h->nReceivers; i++)
        if(h->rec_IDs[i] == receiverID)
            rec_idx = i;
    assert(rec_idx != -1);

    /* update position */
    for(i=0; i<3; i++)
        h->rec_xyz[rec_idx].v[i] = new_position_xyz[i];

    /* All source/receiver combinations for this receiver index will require refreshing */
    for(i=0; i<h->nSources; i++){
        work = (ims_core_workspace*)(h->hCoreWrkSpc[rec_idx][i]);
        work->refreshEchogramFLAG = 1;
    }
}

void ims_shoebox_removeSource
(
    void* hIms,
    long sourceID
)
{
    ims_scene_data *h = (ims_scene_data*)(hIms);
    ims_core_workspace* work;
    voidPtr** new_hCoreWrkSpc;
    int i, src_idx;

    /* Find index corresponding to this source ID */
    src_idx = -1;
    for(i=0; i<h->nSources; i++)
        if(h->src_IDs[i] == sourceID)
            src_idx = i;
    assert(src_idx != -1);

//    /* Shift all the sources above this index down 1 slot, and truncate the src_xyz and src_IDs arrays */
//    for(i=src_idx; src_idx < h->nSources-1; src_idx++){
//        memcpy(&(h->src_xyz[i]), &(h->src_xyz[i+1]), sizeof(position_xyz));
//        memcpy(&(h->src_IDs[i]), &(h->src_IDs[i+1]), sizeof(long));
//    }
//    h->nSources--;
//    h->src_xyz = realloc1d(h->src_xyz, h->nSources*sizeof(position_xyz));
//    h->src_IDs = realloc1d(h->src_IDs, h->nSources*sizeof(long));

    /* Work space handles should also be shifted down */
//    new_hCoreWrkSpc =
//    for(i=src_idx; src_idx < h->nSources; src_idx++){
//
//    }


}
