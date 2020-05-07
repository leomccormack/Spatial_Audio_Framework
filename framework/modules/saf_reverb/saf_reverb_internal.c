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


void ims_shoebox_echogramCreate
(
    void** hEcho
)
{
    *hEcho = malloc1d(sizeof(echogram_data));
    echogram_data *h = (echogram_data*)(*hEcho);

    h->numImageSources = 0;
    h->nChannels = 0;
    h->value = NULL;
    h->time = NULL;
    h->order = NULL;
    h->coords = NULL;
    h->sortedIdx = NULL;
}

void ims_shoebox_echogramResize
(
    void* hEcho,
    int numImageSources,
    int nChannels
)
{
    echogram_data *echogram = (echogram_data*)(hEcho);

    if(echogram->nChannels != nChannels ||  echogram->numImageSources != numImageSources){
        echogram->nChannels = nChannels;
        echogram->numImageSources = numImageSources;
        echogram->value = (float**)realloc2d((void**)echogram->value, numImageSources, nChannels, sizeof(float));
        echogram->time = realloc1d(echogram->time, numImageSources*sizeof(float));
        echogram->order = (int**)realloc2d((void**)echogram->order, numImageSources, 3, sizeof(int));
        echogram->coords = (float**)realloc2d((void**)echogram->coords, numImageSources, 3, sizeof(float));
        echogram->sortedIdx = realloc1d(echogram->sortedIdx, numImageSources*sizeof(int));
    }
}

void ims_shoebox_coreWorkspaceCreate
(
    void** hWork,
    int nBands
)
{
    *hWork = malloc1d(sizeof(ims_core_workspace));
    ims_core_workspace *h = (ims_core_workspace*)(*hWork);
    int i, band;

    /* locals */
    h->d_max = 0.0f;
    h->lengthVec = 0;
    h->numImageSources = 0;
//    for(i=0; i<NUM_WALLS_SHOEBOX; i++)
//        h->abs_wall[i] = -1.0f; /* illegal value (forces reinit) */
    memset(h->room, 0, 3*sizeof(int));
    for(i=0; i<3; i++){
        h->src.v[i] = -1; /* outside the room (forces reinit) */
        h->rec.v[i] = -1;
    }
    h->nBands = nBands;

    /* Internals */
    h->validIDs = NULL;
    h->II = h->JJ = h->KK = NULL;
    h->s_x = h->s_y = h->s_z = h->s_d = NULL;
    h->s_t = h->s_att = NULL;

    h->refreshEchogramFLAG = 1;

    /* Echograms */
    ims_shoebox_echogramCreate( &(h->hEchogram) );
    ims_shoebox_echogramCreate( &(h->hEchogram_rec) );
    h->hEchogram_abs = malloc1d(nBands*sizeof(voidPtr));
    for(band=0; band< nBands; band++)
        ims_shoebox_echogramCreate( &(h->hEchogram_abs[band]) );
}

void ims_shoebox_coreInit
(
    void* hWork,
    int room[3],
    position_xyz src,
    position_xyz rec,
    float maxTime_s,
    float c_ms
)
{
    ims_core_workspace *h = (ims_core_workspace*)(hWork);
    echogram_data *echogram = (echogram_data*)(h->hEchogram);
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

        /* Resize echogram container (only done if needed) */
        ims_shoebox_echogramResize(h->hEchogram, h->numImageSources, 1/*omni-pressure*/);

        /* Copy data into echogram struct */
        for(imsrc = 0, vIdx = 0; imsrc<h->lengthVec; imsrc++){
            if(h->validIDs[imsrc]){
                echogram->time[vIdx]     = h->s_d[imsrc]/c_ms;

                /* reflection propagation attenuation - if distance is <1m set
                 * at attenuation at 1 to avoid amplification */
                echogram->value[vIdx][0]   = h->s_d[imsrc]<=1 ? 1.0f : 1.0f / h->s_d[imsrc];

                /* Order */
                echogram->order[vIdx][0] = h->II[imsrc];
                echogram->order[vIdx][1] = h->JJ[imsrc];
                echogram->order[vIdx][2] = h->KK[imsrc];

                /* Coordinates */
                echogram->coords[vIdx][0] = h->s_x[imsrc];
                echogram->coords[vIdx][1] = h->s_y[imsrc];
                echogram->coords[vIdx][2] = h->s_z[imsrc];
                vIdx++;
            }
        }

        /* Find indices to sort reflections according to propagation time (accending order) */
        sortf(echogram->time, NULL, echogram->sortedIdx, echogram->numImageSources, 0);
    }
}

void ims_shoebox_coreRecModuleSH
(
    void* hWork,
    int sh_order
)
{
    ims_core_workspace *h = (ims_core_workspace*)(hWork);
    echogram_data *echogram = (echogram_data*)(h->hEchogram);
    echogram_data *echogram_rec = (echogram_data*)(h->hEchogram_rec);
    int i, j, nSH;
    float aziElev_rad[2];
    float* sh_gains;

    nSH = ORDER2NSH(sh_order);

    /* Resize container (only done if needed) */
    ims_shoebox_echogramResize(h->hEchogram_rec, echogram->numImageSources, nSH);

    /* Copy 'time', 'coord', 'order', except in accending order w.r.t propogation time  */
    for(i=0; i<echogram_rec->numImageSources; i++){
        echogram_rec->time[i] = echogram->time[echogram->sortedIdx[i]];
        for(j=0; j<3; j++){
            echogram_rec->order[i][j] = echogram->order[echogram->sortedIdx[i]][j];
            echogram_rec->coords[i][j] = echogram->coords[echogram->sortedIdx[i]][j];
        }
        echogram_rec->sortedIdx[i] = i;
    }

    /* Copy 'value' (the core omni-pressure), except also in accending order w.r.t propogation time */
    if(sh_order==0){
        for(i=0; i<echogram_rec->numImageSources; i++)
            echogram_rec->value[i][0] = echogram->value[echogram->sortedIdx[i]][0];
    }
    /* Impose spherical harmonic directivities onto 'value', and store in accending order w.r.t propogation time */
    else{
        sh_gains = malloc1d(nSH*sizeof(float));
        for(i=0; i<echogram_rec->numImageSources; i++){
            /* Cartesian coordinates to spherical coordinates */
            unitCart2Sph(echogram_rec->coords[i], (float*)aziElev_rad);
            aziElev_rad[1] = M_PI/2.0f-aziElev_rad[1]; /* AziElev to AziInclination convention */

            /* Apply spherical harmonic weights */
            getSHreal_recur(sh_order, (float*)aziElev_rad, 1, sh_gains);
            for(j=0; j<nSH; j++){
                echogram_rec->value[i][j] = sh_gains[j] * (echogram->value[echogram->sortedIdx[i]][0]);
            }
        }
        free(sh_gains);
    }
}

void ims_shoebox_coreAbsorptionModule
(
    void* hWork,
    float** abs_wall//,
    //void* hEchogram_abs
)
{
    ims_core_workspace *h = (ims_core_workspace*)(hWork);
    echogram_data *echogram_rec = (echogram_data*)(h->hEchogram_rec);
    echogram_data *echogram_abs;
    int i,j,band;
    float r_x[2], r_y[2], r_z[2];
    float abs_x, abs_y, abs_z, s_abs_tot;


    for(band=0; band < h->nBands; band++){
        echogram_abs = h->hEchogram_abs[band];

        /* Resize container (only done if needed) */
        ims_shoebox_echogramResize(h->hEchogram_abs[band], echogram_rec->numImageSources, echogram_rec->nChannels);

        /* Copy data */
        for(i=0; i<echogram_abs->numImageSources; i++){
            for(j=0; j<echogram_abs->nChannels; j++)
                echogram_abs->value[i][j] = echogram_rec->value[i][j];
            echogram_abs->time[i] = echogram_rec->time[i];
            for(j=0; j<3; j++){
                echogram_abs->order[i][j] = echogram_rec->order[i][j];
                echogram_abs->coords[i][j] = echogram_rec->coords[i][j];
            }
            echogram_abs->sortedIdx[i] = i;
        }

        /* Reflection coefficients given the absorption coefficients for x, y, z
         * walls per frequency */
        r_x[0] = sqrtf(1.0f - abs_wall[band][0]);
        r_x[1] = sqrtf(1.0f - abs_wall[band][1]);
        r_y[0] = sqrtf(1.0f - abs_wall[band][2]);
        r_y[1] = sqrtf(1.0f - abs_wall[band][3]);
        r_z[0] = sqrtf(1.0f - abs_wall[band][4]);
        r_z[1] = sqrtf(1.0f - abs_wall[band][5]);

        /* find total absorption coefficients by calculating the number of hits on
         * every surface, based on the order per dimension */
        for(i=0; i<echogram_abs->numImageSources; i++){
            /* Surfaces intersecting the x-axis */
            if(ISEVEN(echogram_abs->order[i][0]))
                abs_x = powf(r_x[0], (float)abs(echogram_abs->order[i][0])/2.0f) * powf(r_x[1], (float)abs(echogram_abs->order[i][0])/2.0f);
            else if (/* ISODD AND */echogram_abs->order[i][0]>0)
                abs_x = powf(r_x[0], ceilf((float)echogram_abs->order[i][0]/2.0f)) * powf(r_x[1], floorf((float)echogram_abs->order[i][0]/2.0f));
            else /* ISODD AND NEGATIVE */
                abs_x = powf(r_x[0], floorf((float)abs(echogram_abs->order[i][0])/2.0f)) * powf(r_x[1], ceilf((float)abs(echogram_abs->order[i][0])/2.0f));

            /* Surfaces intersecting the y-axis */
            if(ISEVEN(echogram_abs->order[i][1]))
                abs_y = powf(r_y[0], (float)abs(echogram_abs->order[i][1])/2.0f) * powf(r_y[1], (float)abs(echogram_abs->order[i][1])/2.0f);
            else if (/* ISODD AND */echogram_abs->order[i][1]>0)
                abs_y = powf(r_y[0], ceilf((float)echogram_abs->order[i][1]/2.0f)) * powf(r_y[1], floorf((float)echogram_abs->order[i][1]/2.0f));
            else /* ISODD AND NEGATIVE */
                abs_y = powf(r_y[0], floorf((float)abs(echogram_abs->order[i][1])/2.0f)) * powf(r_y[1], ceilf((float)abs(echogram_abs->order[i][1])/2.0f));

            /* Surfaces intersecting the y-axis */
            if(ISEVEN(echogram_abs->order[i][2]))
                abs_z = powf(r_z[0], (float)abs(echogram_abs->order[i][2])/2.0f) * powf(r_z[1], (float)abs(echogram_abs->order[i][2])/2.0f);
            else if (/* ISODD AND */echogram_abs->order[i][2]>0)
                abs_z = powf(r_z[0], ceilf((float)echogram_abs->order[i][2]/2.0f)) * powf(r_z[1], floorf((float)echogram_abs->order[i][2]/2.0f));
            else /* ISODD AND NEGATIVE */
                abs_z = powf(r_z[0], floorf((float)abs(echogram_abs->order[i][2])/2.0f)) * powf(r_z[1], ceilf((float)abs(echogram_abs->order[i][2])/2.0f));

            /* Apply Absorption */
            s_abs_tot = abs_x * abs_y * abs_z;
            for(j=0; j<echogram_abs->nChannels; j++)
                echogram_abs->value[i][j] *= s_abs_tot;
        }
    }
}


void ims_shoebox_renderRIR
(
    void* hWork,
    //void* hEchogram_abs,
    int fractionalDelayFLAG,
    float* rir // precompute, we know the maxlength and fs already, so should be able to do this for each workspace handle?
)
{
    ims_core_workspace *h = (ims_core_workspace*)(hWork);
    //echogram_data *echogram_abs = (echogram_data*)(hEchogram_abs);

    
}
