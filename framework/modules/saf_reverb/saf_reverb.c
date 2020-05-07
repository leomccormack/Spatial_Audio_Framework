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

void ims_shoeboxroom_create
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
    int band, wall;

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
    h->abs_wall = (float**)malloc2d(nOctBands, NUM_WALLS_SHOEBOX, sizeof(float));
    for(band=0; band<nOctBands; band++)
        for(wall=0; wall<NUM_WALLS_SHOEBOX; wall++)
            h->abs_wall[band][wall] = abs_wall[band*NUM_WALLS_SHOEBOX+wall];

    /* Default is no sources or receivers in the room */
    h->src_xyz = NULL;
    h->rec_xyz = NULL;
    h->src_IDs = NULL;
    h->rec_IDs = NULL;
    h->nSources = 0;
    h->nReceivers = 0;

    /* ims_core_workspace per source / receiver */
    h->hCoreWrkSpc = NULL;
}

long ims_shoeboxroom_addSource
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
    h->src_xyz = realloc1d(h->src_xyz, h->nSources*sizeof(position_xyz));
    for(i=0; i<3; i++)
        h->src_xyz[h->nSources-1].v[i] = src_xyz[i];

    /* Assign unique ID */
    newID = 0;
    h->src_IDs = realloc1d(h->src_IDs, h->nSources*sizeof(long));
    for(i=0; i<h->nSources-1; i++)
        if(h->src_IDs[i]==newID) /* check ID is not in use */
            newID++;
    h->src_IDs[h->nSources-1] = newID;

    /* Create workspace for all source/receiver combinations, for this new source */
    h->hCoreWrkSpc = (voidPtr**)realloc2d((void**)h->hCoreWrkSpc, h->nReceivers, h->nSources, sizeof(voidPtr));
    for(rec=0; rec<h->nReceivers; rec++)
        ims_shoebox_coreWorkspaceCreate(&(h->hCoreWrkSpc[rec][h->nSources-1]), h->nBands);

    return newID;
}

long ims_shoeboxroom_addReceiver
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
    h->rec_xyz = realloc1d(h->rec_xyz, h->nReceivers*sizeof(position_xyz));
    for(i=0; i<3; i++)
        h->rec_xyz[h->nReceivers-1].v[i] = rec_xyz[i];

    /* Assign unique ID */
    newID = 0;
    h->rec_IDs = realloc1d(h->rec_IDs, h->nReceivers*sizeof(long));
    for(i=0; i<h->nReceivers-1; i++)
        if(h->rec_IDs[i]==newID) /* check ID is not in use */
            newID++;
    h->rec_IDs[h->nReceivers-1] = newID;

    /* Create workspace for all receiver/source combinations, for this new receiver */
    h->hCoreWrkSpc = (voidPtr**)realloc2d((void**)h->hCoreWrkSpc, h->nReceivers, h->nSources, sizeof(voidPtr));
    for(src=0; src<h->nSources; src++)
        ims_shoebox_coreWorkspaceCreate(&(h->hCoreWrkSpc[h->nReceivers-1][src]), h->nBands);

    return newID;
}

void ims_shoeboxroom_updateSource
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

void ims_shoeboxroom_updateReceiver
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

void ims_shoeboxroom_renderEchogramSH
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

void ims_shoeboxroom_renderSHRIRs
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
