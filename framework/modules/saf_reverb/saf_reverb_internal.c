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
 * @file saf_reverb_internal.c
 * @brief Internal part of the reverb processing module (saf_reverb)
 *
 * ...
 *
 * @author Leo McCormack
 * @date 06.05.2020
 */

#include "saf_reverb.h"
#include "saf_reverb_internal.h"


void ims_shoebox_coreWorkspaceCreate
(
    void** hWork
)
{
    *hWork = malloc1d(sizeof(ims_core_workspace));
    ims_core_workspace *h = (ims_core_workspace*)(*hWork);
    int i;

    h->d_max = 0.0f;
    h->lengthVec = 0;
    h->numImageSources = 0;
    memset(h->room, 0, 3*sizeof(int));
    for(i=0; i<3; i++){
        h->src.v[i] = -1; /* outside the room (forces reinit) */
        h->rec.v[i] = -1;
    }
    h->validIDs = NULL;
    h->II = h->JJ = h->KK = NULL;
    h->s_x = h->s_y = h->s_z = h->s_d = NULL;
    h->s_t = h->s_att = NULL;
}

void ims_shoebox_core
(
    void* hWork,
    int room[3],
    position_xyz src,
    position_xyz rec,
    float maxTime_s,
    float c_ms,
    echogram_data* echogram
)
{
    ims_core_workspace *h = (ims_core_workspace*)(hWork);
    position_xyz src_orig, rec_orig;
    int imsrc, vIdx;
    int ii, jj, kk;
    float d_max;

    d_max = maxTime_s*c_ms;

    /* move origin to the centre of the room */
    src_orig.x = src.x - (float)room[0]/2.0f;
    src_orig.y = (float)room[1]/2.0f - src.y;
    src_orig.z = src.z - (float)room[2]/2.0f;
    rec_orig.x = rec.x - (float)room[0]/2.0f;
    rec_orig.y = (float)room[1]/2.0f - rec.y;
    rec_orig.z = rec.z - (float)room[2]/2.0f;

    /* Update indices only if the maximum permitted delay or room dimensions have changed */
    if( (h->d_max != d_max) ||
        (h->room[0] != room[0]) || (h->room[1] != room[1]) || (h->room[2] != room[2]) )
    {
        h->d_max = d_max;
        memcpy(h->room, room, 3*sizeof(int));
        h->Nx = (int)(d_max/(float)room[0] + 1.0f); /* ceil */
        h->Ny = (int)(d_max/(float)room[1] + 1.0f); /* ceil */
        h->Nz = (int)(d_max/(float)room[2] + 1.0f); /* ceil */
        h->lengthVec = (2*(h->Nx)+1) * (2*(h->Ny)+1) * (2*(h->Nz)+1);

        /* i,j,k indices for calculation in x,y,z respectively */
        h->II = realloc1d(h->II, h->lengthVec*sizeof(float));
        h->JJ = realloc1d(h->JJ, h->lengthVec*sizeof(float));
        h->KK = realloc1d(h->KK, h->lengthVec*sizeof(float));
        ii = - h->Nx; jj = - h->Ny; kk = - h->Nz;
        for(imsrc = 0; imsrc<h->lengthVec; imsrc++){
            h->II[imsrc] = (float)ii;
            h->JJ[imsrc] = (float)jj;
            h->KK[imsrc] = (float)kk;
            ii++;
            if(ii>h->Nx){
                ii = -h->Nx;
                jj++;
            }
            if(jj>h->Ny){
                jj = -h->Ny;
                kk++;
            }
            if(kk>h->Nz){
                kk = -h->Nz;
            }
        }

        /* Re-allocate memory */
        h->validIDs = realloc1d(h->validIDs, h->lengthVec*sizeof(int));
        h->s_x = realloc1d(h->s_x, h->lengthVec*sizeof(float));
        h->s_y = realloc1d(h->s_y, h->lengthVec*sizeof(float));
        h->s_z = realloc1d(h->s_z, h->lengthVec*sizeof(float));
        h->s_d = realloc1d(h->s_d, h->lengthVec*sizeof(float));
        h->s_t = realloc1d(h->s_t, h->lengthVec*sizeof(float));
        h->s_att = realloc1d(h->s_att, h->lengthVec*sizeof(float));
    }

    /* Update echogram only if the source/receiver positions or room dimensions have changed */
    if( (h->rec.x != rec_orig.x) || (h->rec.y != rec_orig.y) || (h->rec.z != rec_orig.z) ||
        (h->src.x != src_orig.x) || (h->src.y != src_orig.y) || (h->src.z != src_orig.z) ||
        (h->room[0] != room[0]) || (h->room[1] != room[1]) || (h->room[2] != room[2]))
    {
        memcpy(h->room, room, 3*sizeof(int));
        memcpy(&(h->rec), &rec_orig, sizeof(position_xyz));
        memcpy(&(h->src), &src_orig, sizeof(position_xyz));

        /* image source coordinates with respect to receiver, and distance */
        for(imsrc = 0; imsrc<h->lengthVec; imsrc++){
            h->s_x[imsrc] = h->II[imsrc]*(float)room[0] + powf(-1.0f, h->II[imsrc])*src_orig.x - rec_orig.x;
            h->s_y[imsrc] = h->JJ[imsrc]*(float)room[1] + powf(-1.0f, h->JJ[imsrc])*src_orig.y - rec_orig.y;
            h->s_z[imsrc] = h->KK[imsrc]*(float)room[2] + powf(-1.0f, h->KK[imsrc])*src_orig.z - rec_orig.z;
            h->s_d[imsrc] = sqrtf(powf(h->s_x[imsrc], 2.0f) + powf(h->s_y[imsrc], 2.0f) + powf(h->s_z[imsrc], 2.0f));
        }

        /* Determine the indices where the distance is below the specified maximum */ 
        for(imsrc = 0, h->numImageSources = 0; imsrc<h->lengthVec; imsrc++){
            if(h->s_d[imsrc]<d_max){
                h->validIDs[imsrc] = 1;
                h->numImageSources++; /* (within maximum distance) */
            }
            else
                h->validIDs[imsrc] = 0;
        }

        /* Resize echogram struct if needed */
        if(echogram->numImageSources != h->numImageSources){
            echogram->numImageSources = h->numImageSources;
            echogram->value = (float**)realloc2d((void**)echogram->value, h->numImageSources, 1/*omni-pressure*/, sizeof(float));
            echogram->time = realloc1d(echogram->time, h->numImageSources*sizeof(float));
            echogram->order = realloc1d(echogram->order, h->numImageSources*sizeof(reflOrder));
            echogram->coords = realloc1d(echogram->coords, h->numImageSources*sizeof(position_xyz));
            echogram->sortedIdx = realloc1d(echogram->sortedIdx, h->numImageSources*sizeof(int));
        }

        /* Copy data into echogram struct */
        for(imsrc = 0, vIdx = 0; imsrc<h->lengthVec; imsrc++){
            if(h->validIDs[imsrc]){
                echogram->time[vIdx]     = h->s_d[imsrc]/c_ms;

                /* reflection propagation attenuation - if distance is <1m set
                 * at attenuation at 1 to avoid amplification */
                echogram->value[vIdx][0]   = h->s_d[imsrc]<=1 ? 1.0f : 1.0f / h->s_d[imsrc];

                /* Order */
                echogram->order[vIdx].v[0] = h->II[imsrc];
                echogram->order[vIdx].v[1] = h->JJ[imsrc];
                echogram->order[vIdx].v[2] = h->KK[imsrc];

                /* Coordinates */
                echogram->coords[vIdx].x = h->s_x[imsrc];
                echogram->coords[vIdx].y = h->s_y[imsrc];
                echogram->coords[vIdx].z = h->s_z[imsrc];
                vIdx++;
            }
        }

        /* Find indices to sort reflections according to propagation time (accending order) */
        sortf(echogram->time, NULL, echogram->sortedIdx, echogram->numImageSources, 0);
    }
}
