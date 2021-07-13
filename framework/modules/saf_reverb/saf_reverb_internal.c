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
 * @ingroup Reverb
 * @brief Internal source for the reverb processing module (#SAF_REVERB_MODULE)
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

void ims_shoebox_echogramCreate
(
    void** phEcho,
    int include_rt_vars
)
{
    *phEcho = malloc1d(sizeof(echogram_data));
    echogram_data *ec = (echogram_data*)(*phEcho);
    int i;

    saf_assert(include_rt_vars==0 || include_rt_vars==1, "include_rt_vars is a bool");

    /* Echogram data */
    ec->numImageSources = 0;
    ec->nChannels = 0;
    ec->value = NULL;
    ec->time = NULL;
    ec->order = NULL;
    ec->coords = NULL;
    ec->sortedIdx = NULL;

    /* Optional helper variables */
    ec->include_rt_vars = include_rt_vars;
    ec->tmp1 = NULL;
    ec->tmp2 = NULL;
    ec->rIdx = NULL;
    for(i=0; i<IMS_LAGRANGE_ORDER; i++)
        ec->rIdx_frac[i] = NULL;
    ec->h_frac = NULL;
    ec->cb_vals = NULL;
    ec->contrib = NULL;
    ec->ones_dummy = NULL;
}

void ims_shoebox_echogramResize
(
    void* hEcho,
    int numImageSources,
    int nChannels
)
{
    echogram_data *ec = (echogram_data*)(hEcho);
    int i;

    if(ec->nChannels != nChannels || ec->numImageSources != numImageSources){
        /* Resize echogram data */
        ec->nChannels = nChannels;
        ec->numImageSources = numImageSources;
        ec->value = (float**)realloc2d((void**)ec->value, nChannels, numImageSources, sizeof(float));
        ec->time = realloc1d(ec->time, numImageSources*sizeof(float));
        ec->order = (int**)realloc2d((void**)ec->order, numImageSources, 3, sizeof(int));
        ec->coords = realloc1d(ec->coords, numImageSources * sizeof(ims_pos_xyz));
        ec->sortedIdx = realloc1d(ec->sortedIdx, numImageSources*sizeof(int));

        /* Resize the optional helper variables (used for run-time speed-ups) */
        if(ec->include_rt_vars){
            ec->tmp1 = realloc1d(ec->tmp1, numImageSources*sizeof(float));
            ec->tmp2 = realloc1d(ec->tmp2, numImageSources*sizeof(float));
            ec->rIdx = realloc1d(ec->rIdx, numImageSources*sizeof(int));
            for(i=0; i<IMS_LAGRANGE_ORDER; i++)
                ec->rIdx_frac[i] = realloc1d(ec->rIdx_frac[i], numImageSources*sizeof(int));
            ec->h_frac = (float**)realloc2d((void**)ec->h_frac, IMS_LAGRANGE_ORDER+1, numImageSources, sizeof(float));
            ec->cb_vals = (float**)realloc2d((void**)ec->cb_vals, nChannels, numImageSources, sizeof(float));
            ec->contrib = (float**)realloc2d((void**)ec->contrib, nChannels, numImageSources, sizeof(float));
            ec->ones_dummy = realloc1d(ec->ones_dummy, numImageSources*sizeof(float));
            for(i=0; i<numImageSources; i++)
                ec->ones_dummy[i] = 1.0f;
        }
    }
}

void ims_shoebox_echogramCopy
(
    void* hEchoX,
    void* hEchoY
)
{
    echogram_data *ecX = (echogram_data*)(hEchoX);
    echogram_data *ecY = (echogram_data*)(hEchoY);
    int nChannels, numImageSources;

    saf_assert(hEchoX!=NULL && hEchoY!=NULL, "Echograms must be allocated first");
    if(hEchoX==hEchoY)
        return; /* no copying required... */

    if(ecX->nChannels==0 || ecX->numImageSources==0)
        return;

    /* Resize container 'Y' to be the same as container 'X' (if needed) */
    nChannels = ecX->nChannels;
    numImageSources = ecX->numImageSources;
    if(ecY->nChannels != nChannels || ecY->numImageSources != numImageSources)
        ims_shoebox_echogramResize(hEchoY, numImageSources, nChannels);

    /* Copy echogram data from X to Y */
    cblas_scopy(nChannels*numImageSources, FLATTEN2D(ecX->value), 1, FLATTEN2D(ecY->value), 1);
    cblas_scopy(numImageSources, ecX->time, 1, ecY->time, 1);
    memcpy(FLATTEN2D(ecY->order), FLATTEN2D(ecX->order), numImageSources*3*sizeof(int));
    memcpy(ecY->coords, ecX->coords, numImageSources*sizeof(ims_pos_xyz));
    memcpy(ecY->sortedIdx, ecX->sortedIdx, numImageSources*sizeof(int));
}

void ims_shoebox_echogramDestroy
(
    void** phEcho
)
{
    echogram_data *ec = (echogram_data*)(*phEcho);
    int i;

    if(ec!=NULL){
        /* Free echogram data */
        free(ec->value);
        free(ec->time);
        free(ec->order);
        free(ec->coords);
        free(ec->sortedIdx);

        /* Free the optional helper variables */
        if(ec->include_rt_vars){
            free(ec->tmp1);
            free(ec->tmp2);
            free(ec->rIdx);
            for(i=0; i<IMS_LAGRANGE_ORDER; i++)
                free(ec->rIdx_frac[i]);
            free(ec->h_frac);
            free(ec->cb_vals);
            free(ec->contrib);
            free(ec->ones_dummy);
        }

        free(ec);
        ec=NULL;
        *phEcho = NULL;
    }
}

void ims_shoebox_coreWorkspaceCreate
(
    void** phWork,
    int nBands
)
{
    ims_shoebox_coreWorkspaceDestroy(phWork);
    *phWork = malloc1d(sizeof(ims_core_workspace));
    ims_core_workspace *wrk = (ims_core_workspace*)(*phWork);
    int i, band;

    /* locals */
    wrk->d_max = -1.0f;
    wrk->N_max = -1;
    wrk->lengthVec = 0;
    wrk->numImageSources = 0;
    memset(wrk->room, 0, 3*sizeof(int));
    for(i=0; i<3; i++){
        wrk->src.v[i] = -1;
        wrk->rec.v[i] = -1;
    }
    wrk->nBands = nBands;

    /* Internals */
    wrk->validIDs = NULL;
    wrk->II  = wrk->JJ    = wrk->KK  = NULL;
    wrk->iII = wrk->iJJ   = wrk->iKK  = NULL;
    wrk->s_x = wrk->s_y   = wrk->s_z = wrk->s_d = NULL;
    wrk->s_t = wrk->s_att = NULL;
    wrk->s_ord = NULL;

    /* Echogram containers */
    wrk->refreshEchogramFLAG = 1;
    ims_shoebox_echogramCreate(&(wrk->hEchogram), 0);
    ims_shoebox_echogramCreate(&(wrk->hEchogram_rec), 0);
    wrk->hEchogram_abs = malloc1d(nBands*sizeof(voidPtr));
    wrk->hPrevEchogram_abs = malloc1d(nBands*sizeof(voidPtr));
    for(band=0; band<nBands; band++){
        ims_shoebox_echogramCreate(&(wrk->hEchogram_abs[band]), 1);
        ims_shoebox_echogramCreate(&(wrk->hPrevEchogram_abs[band]), 1);
    }

    /* Room impulse responses */
    wrk->refreshRIRFLAG = 1;
    wrk->rir_len_samples = 0;
    wrk->rir_len_seconds = 0.0f;
    wrk->rir_bands = (float***)malloc1d(nBands*sizeof(float**));
    for(band=0; band < nBands; band++)
        wrk->rir_bands[band] = NULL;
}

void ims_shoebox_coreWorkspaceDestroy
(
    void** phWork
)
{
    ims_core_workspace *wrk = (ims_core_workspace*)(*phWork);
    int band;

    if(wrk!=NULL){
        /* free internals */
        free(wrk->validIDs);
        free(wrk->II);
        free(wrk->JJ);
        free(wrk->KK);
        free(wrk->iII);
        free(wrk->iJJ);
        free(wrk->iKK);
        free(wrk->s_x);
        free(wrk->s_y);
        free(wrk->s_z);
        free(wrk->s_d);
        free(wrk->s_t);
        free(wrk->s_att);
        free(wrk->s_ord);

        /* Destroy echogram containers */
        ims_shoebox_echogramDestroy(&(wrk->hEchogram));
        ims_shoebox_echogramDestroy(&(wrk->hEchogram_rec));
        for(band=0; band< wrk->nBands; band++){
            ims_shoebox_echogramDestroy(&(wrk->hEchogram_abs[band]));
            ims_shoebox_echogramDestroy(&(wrk->hPrevEchogram_abs[band]));
        }
        free(wrk->hEchogram_abs);
        free(wrk->hPrevEchogram_abs);
 
        /* free rirs */
        for(band=0; band < wrk->nBands; band++)
            free(wrk->rir_bands[band]);

        free(wrk);
        wrk=NULL;
        *phWork = NULL;
    }
}

void ims_shoebox_coreInitT
(
    void* hWork,
    float room[3],
    ims_pos_xyz src,
    ims_pos_xyz rec,
    float maxTime_s,
    float c_ms
)
{
    ims_core_workspace *wrk = (ims_core_workspace*)(hWork);
    echogram_data *echogram = (echogram_data*)(wrk->hEchogram);
    ims_pos_xyz src_orig, rec_orig;
    int imsrc, vIdx;
    int ii, jj, kk;
    float d_max;

    d_max = maxTime_s*c_ms;

    /* move origin to the centre of the room */
    src_orig.x = src.x - room[0]/2.0f;
    src_orig.y = room[1]/2.0f - src.y;
    src_orig.z = src.z - room[2]/2.0f;
    rec_orig.x = rec.x - room[0]/2.0f;
    rec_orig.y = room[1]/2.0f - rec.y;
    rec_orig.z = rec.z - room[2]/2.0f;

    /* Update indices only if the maximum permitted delay or room dimensions have changed */
    if( (wrk->d_max != d_max) ||
        (wrk->room[0] != room[0]) || (wrk->room[1] != room[1]) || (wrk->room[2] != room[2]) )
    {
        wrk->Nx = (int)(d_max/room[0] + 1.0f); /* ceil */
        wrk->Ny = (int)(d_max/room[1] + 1.0f); /* ceil */
        wrk->Nz = (int)(d_max/room[2] + 1.0f); /* ceil */
        wrk->lengthVec = (2*(wrk->Nx)+1) * (2*(wrk->Ny)+1) * (2*(wrk->Nz)+1);

        /* i,j,k indices for calculation in x,y,z respectively */
        wrk->II = realloc1d(wrk->II, wrk->lengthVec*sizeof(float));
        wrk->JJ = realloc1d(wrk->JJ, wrk->lengthVec*sizeof(float));
        wrk->KK = realloc1d(wrk->KK, wrk->lengthVec*sizeof(float));
        ii = -(wrk->Nx); jj = -(wrk->Ny); kk = -(wrk->Nz);
        for(imsrc = 0; imsrc<wrk->lengthVec; imsrc++){
            wrk->II[imsrc] = (float)ii;
            wrk->JJ[imsrc] = (float)jj;
            wrk->KK[imsrc] = (float)kk;
            ii++;
            if(ii>wrk->Nx){
                ii = -(wrk->Nx);
                jj++;
            }
            if(jj>wrk->Ny){
                jj = -(wrk->Ny);
                kk++;
            }
            if(kk>wrk->Nz){
                kk = -(wrk->Nz);
            }
        }

        /* Re-allocate memory */
        wrk->validIDs = realloc1d(wrk->validIDs, wrk->lengthVec*sizeof(int));
        wrk->s_x = realloc1d(wrk->s_x, wrk->lengthVec*sizeof(float));
        wrk->s_y = realloc1d(wrk->s_y, wrk->lengthVec*sizeof(float));
        wrk->s_z = realloc1d(wrk->s_z, wrk->lengthVec*sizeof(float));
        wrk->s_d = realloc1d(wrk->s_d, wrk->lengthVec*sizeof(float));
        wrk->s_t = realloc1d(wrk->s_t, wrk->lengthVec*sizeof(float));
        wrk->s_att = realloc1d(wrk->s_att, wrk->lengthVec*sizeof(float));
    }

    /* Update echogram only if the source/receiver positions or room dimensions have changed */
    if( (wrk->d_max != d_max) ||
        (wrk->rec.x != rec_orig.x) || (wrk->rec.y != rec_orig.y) || (wrk->rec.z != rec_orig.z) ||
        (wrk->src.x != src_orig.x) || (wrk->src.y != src_orig.y) || (wrk->src.z != src_orig.z) ||
        (wrk->room[0] != room[0]) || (wrk->room[1] != room[1]) || (wrk->room[2] != room[2]))
    {
        wrk->d_max = d_max;
        wrk->room[0] = room[0];
        wrk->room[1] = room[1];
        wrk->room[2] = room[2];
        memcpy(&(wrk->rec), &rec_orig, sizeof(ims_pos_xyz));
        memcpy(&(wrk->src), &src_orig, sizeof(ims_pos_xyz));

        /* image source coordinates with respect to receiver, and distance */
        for(imsrc = 0; imsrc<wrk->lengthVec; imsrc++){
            wrk->s_x[imsrc] = wrk->II[imsrc]*room[0] + powf(-1.0f, wrk->II[imsrc])*src_orig.x - rec_orig.x;
            wrk->s_y[imsrc] = wrk->JJ[imsrc]*room[1] + powf(-1.0f, wrk->JJ[imsrc])*src_orig.y - rec_orig.y;
            wrk->s_z[imsrc] = wrk->KK[imsrc]*room[2] + powf(-1.0f, wrk->KK[imsrc])*src_orig.z - rec_orig.z;
            wrk->s_d[imsrc] = sqrtf(powf(wrk->s_x[imsrc], 2.0f) + powf(wrk->s_y[imsrc], 2.0f) + powf(wrk->s_z[imsrc], 2.0f));
        }

        /* Determine the indices where the distance is below the specified maximum */ 
        for(imsrc = 0, wrk->numImageSources = 0; imsrc<wrk->lengthVec; imsrc++){
            if(wrk->s_d[imsrc]<d_max){
                wrk->validIDs[imsrc] = 1;
                wrk->numImageSources++; /* (within maximum distance) */
            }
            else
                wrk->validIDs[imsrc] = 0;
        }

        /* Resize echogram container (only done if needed) */
        ims_shoebox_echogramResize(wrk->hEchogram, wrk->numImageSources, 1/*omni-pressure*/);

        /* Copy data into echogram struct */
        for(imsrc = 0, vIdx = 0; imsrc<wrk->lengthVec; imsrc++){
            if(wrk->validIDs[imsrc]){
                echogram->time[vIdx]     = wrk->s_d[imsrc]/c_ms;

                /* reflection propagation attenuation - if distance is <1m set
                 * attenuation to 1 to avoid amplification */
                echogram->value[0][vIdx]   = wrk->s_d[imsrc]<=1 ? 1.0f : 1.0f / wrk->s_d[imsrc];

                /* Order */
                echogram->order[vIdx][0] = (int)(wrk->II[imsrc] + (wrk->II[imsrc]>0 ? 0.1f : -0.1f)); /* round */
                echogram->order[vIdx][1] = (int)(wrk->JJ[imsrc] + (wrk->JJ[imsrc]>0 ? 0.1f : -0.1f));
                echogram->order[vIdx][2] = (int)(wrk->KK[imsrc] + (wrk->KK[imsrc]>0 ? 0.1f : -0.1f));

                /* Coordinates */
                echogram->coords[vIdx].x = wrk->s_x[imsrc];
                echogram->coords[vIdx].y = wrk->s_y[imsrc];
                echogram->coords[vIdx].z = wrk->s_z[imsrc];
                vIdx++;
            }
        }

        /* Find indices to sort reflections according to propagation time (accending order) */
        sortf(echogram->time, NULL, echogram->sortedIdx, echogram->numImageSources, 0);
    }
}

void ims_shoebox_coreInitN
(
    void* hWork,
    float room[3],
    ims_pos_xyz src,
    ims_pos_xyz rec,
    int maxN,
    float c_ms
)
{
    ims_core_workspace *wrk = (ims_core_workspace*)(hWork);
    echogram_data *echogram = (echogram_data*)(wrk->hEchogram);
    ims_pos_xyz src_orig, rec_orig;
    int imsrc;
    int ii, jj, kk;

    /* move origin to the centre of the room */
    src_orig.x = src.x - room[0]/2.0f;
    src_orig.y = room[1]/2.0f - src.y;
    src_orig.z = src.z - room[2]/2.0f;
    rec_orig.x = rec.x - room[0]/2.0f;
    rec_orig.y = room[1]/2.0f - rec.y;
    rec_orig.z = rec.z - room[2]/2.0f;

    /* Update indices only if the maximum reflection order has changed */
    if( (wrk->N_max != maxN) )
    {
        wrk->lengthVec = (2*(maxN)+1) * (2*(maxN)+1) * (2*(maxN)+1);

        /* i,j,k indices for calculation in x,y,z respectively */
        wrk->iII = realloc1d(wrk->iII, wrk->lengthVec*sizeof(int));
        wrk->iJJ = realloc1d(wrk->iJJ, wrk->lengthVec*sizeof(int));
        wrk->iKK = realloc1d(wrk->iKK, wrk->lengthVec*sizeof(int));
        wrk->s_ord = realloc1d(wrk->s_ord, wrk->lengthVec*sizeof(int));
        ii = -maxN; jj = -maxN; kk = -maxN;
        for(imsrc = 0; imsrc<wrk->lengthVec; imsrc++){
            wrk->iII[imsrc] = ii;
            wrk->iJJ[imsrc] = jj;
            wrk->iKK[imsrc] = kk;
            wrk->s_ord[imsrc] = abs(wrk->iII[imsrc]) + abs(wrk->iJJ[imsrc]) + abs(wrk->iKK[imsrc]);
            ii++;
            if(ii>maxN){
                ii = -(maxN);
                jj++;
            }
            if(jj>maxN){
                jj = -(maxN);
                kk++;
            }
            if(kk>maxN){
                kk = -(maxN);
            }
        }

        /* Cull the indices where the order is below the specified maximum */
        wrk->II = realloc1d(wrk->II, wrk->lengthVec*sizeof(float));
        wrk->JJ = realloc1d(wrk->JJ, wrk->lengthVec*sizeof(float));
        wrk->KK = realloc1d(wrk->KK, wrk->lengthVec*sizeof(float));
        for(imsrc = 0, wrk->numImageSources = 0; imsrc<wrk->lengthVec; imsrc++){
            if(wrk->s_ord[imsrc]<=maxN){
                wrk->II[wrk->numImageSources] = (float)wrk->iII[imsrc];
                wrk->JJ[wrk->numImageSources] = (float)wrk->iJJ[imsrc];
                wrk->KK[wrk->numImageSources] = (float)wrk->iKK[imsrc];
                wrk->numImageSources++; /* (within maximum order) */
            }
        }

        /* Re-allocate memory */
        wrk->s_x = realloc1d(wrk->s_x, wrk->numImageSources*sizeof(float));
        wrk->s_y = realloc1d(wrk->s_y, wrk->numImageSources*sizeof(float));
        wrk->s_z = realloc1d(wrk->s_z, wrk->numImageSources*sizeof(float));
        wrk->s_d = realloc1d(wrk->s_d, wrk->numImageSources*sizeof(float));
        wrk->s_t = realloc1d(wrk->s_t, wrk->numImageSources*sizeof(float));
        wrk->s_att = realloc1d(wrk->s_att, wrk->numImageSources*sizeof(float));
    }

    /* Update echogram only if the maximum reflection order, source/receiver positions or room dimensions have changed */
    if( (wrk->N_max != maxN) ||
        (wrk->rec.x != rec_orig.x) || (wrk->rec.y != rec_orig.y) || (wrk->rec.z != rec_orig.z) ||
        (wrk->src.x != src_orig.x) || (wrk->src.y != src_orig.y) || (wrk->src.z != src_orig.z) ||
        (wrk->room[0] != room[0]) || (wrk->room[1] != room[1]) || (wrk->room[2] != room[2]))
    {
        wrk->N_max = maxN;
        wrk->room[0] = room[0];
        wrk->room[1] = room[1];
        wrk->room[2] = room[2];
        memcpy(&(wrk->rec), &rec_orig, sizeof(ims_pos_xyz));
        memcpy(&(wrk->src), &src_orig, sizeof(ims_pos_xyz));

        /* image source coordinates with respect to receiver, and distance */
        for(imsrc = 0; imsrc<wrk->numImageSources; imsrc++){
            wrk->s_x[imsrc] = wrk->II[imsrc]*room[0] + powf(-1.0f, wrk->II[imsrc])*src_orig.x - rec_orig.x;
            wrk->s_y[imsrc] = wrk->JJ[imsrc]*room[1] + powf(-1.0f, wrk->JJ[imsrc])*src_orig.y - rec_orig.y;
            wrk->s_z[imsrc] = wrk->KK[imsrc]*room[2] + powf(-1.0f, wrk->KK[imsrc])*src_orig.z - rec_orig.z;
            wrk->s_d[imsrc] = sqrtf(powf(wrk->s_x[imsrc], 2.0f) + powf(wrk->s_y[imsrc], 2.0f) + powf(wrk->s_z[imsrc], 2.0f));
        }

        /* Resize echogram container (only done if needed) */
        ims_shoebox_echogramResize(wrk->hEchogram, wrk->numImageSources, 1/*omni-pressure*/);

        /* Copy data into echogram struct */
        for(imsrc = 0; imsrc<wrk->numImageSources; imsrc++){
            echogram->time[imsrc]     = wrk->s_d[imsrc]/c_ms;

            /* reflection propagation attenuation - if distance is <1m set
             * attenuation to 1 to avoid amplification */
            echogram->value[0][imsrc] = wrk->s_d[imsrc]<=1 ? 1.0f : 1.0f / wrk->s_d[imsrc];

            /* Order */
            echogram->order[imsrc][0] = (int)(wrk->II[imsrc] + (wrk->II[imsrc]>0 ? 0.1f : -0.1f)); /* round */
            echogram->order[imsrc][1] = (int)(wrk->JJ[imsrc] + (wrk->JJ[imsrc]>0 ? 0.1f : -0.1f));
            echogram->order[imsrc][2] = (int)(wrk->KK[imsrc] + (wrk->KK[imsrc]>0 ? 0.1f : -0.1f));

            /* Coordinates */
            echogram->coords[imsrc].x = wrk->s_x[imsrc];
            echogram->coords[imsrc].y = wrk->s_y[imsrc];
            echogram->coords[imsrc].z = wrk->s_z[imsrc];
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
    ims_core_workspace *wrk = (ims_core_workspace*)(hWork);
    echogram_data *echogram = (echogram_data*)(wrk->hEchogram);
    echogram_data *echogram_rec = (echogram_data*)(wrk->hEchogram_rec);
    int i, j, nSH;
    float aziElev_rad[2];
    float* sh_gains;

    nSH = ORDER2NSH(sh_order);

    /* Resize container (only done if needed) */
    ims_shoebox_echogramResize(wrk->hEchogram_rec, echogram->numImageSources, nSH);

    /* Copy 'time', 'coord', 'order', except in ascending order w.r.t propagation time  */
    for(i=0; i<echogram_rec->numImageSources; i++){
        echogram_rec->time[i] = echogram->time[echogram->sortedIdx[i]];
        echogram_rec->order[i][0] = echogram->order[echogram->sortedIdx[i]][0];
        echogram_rec->order[i][1] = echogram->order[echogram->sortedIdx[i]][1];
        echogram_rec->order[i][2] = echogram->order[echogram->sortedIdx[i]][2];
        echogram_rec->coords[i].v[0] = echogram->coords[echogram->sortedIdx[i]].v[0];
        echogram_rec->coords[i].v[1] = echogram->coords[echogram->sortedIdx[i]].v[1];
        echogram_rec->coords[i].v[2] = echogram->coords[echogram->sortedIdx[i]].v[2];

        echogram_rec->sortedIdx[i] = i; /* Should already sorted in ims_shoebox_coreInit() */
    }

    /* Copy 'value' (the core omni-pressure), except also in ascending order w.r.t propagation time */
    if(sh_order==0){
        for(i=0; i<echogram_rec->numImageSources; i++)
            echogram_rec->value[0][i] = echogram->value[0][echogram->sortedIdx[i]];
    }
    /* Impose spherical harmonic directivities onto 'value', and store in ascending order w.r.t propagation time */
    else{
        sh_gains = malloc1d(nSH*sizeof(float));
        for(i=0; i<echogram_rec->numImageSources; i++){
            /* Cartesian coordinates to spherical coordinates */
            unitCart2sph(echogram_rec->coords[i].v, 1, 0, (float*)aziElev_rad); 
            aziElev_rad[1] = SAF_PI/2.0f-aziElev_rad[1]; /* AziElev to AziInclination conversion */

            /* Apply spherical harmonic weights */
            getSHreal_recur(sh_order, (float*)aziElev_rad, 1, sh_gains);
            for(j=0; j<nSH; j++)
                echogram_rec->value[j][i] = sh_gains[j] * (echogram->value[0][echogram->sortedIdx[i]]);
        }
        free(sh_gains);
    }
}

void ims_shoebox_coreAbsorptionModule
(
    void* hWork,
    float** abs_wall
)
{
    ims_core_workspace *wrk = (ims_core_workspace*)(hWork);
    echogram_data *echogram_abs;
    int i,band,ch;
    float r_x[2], r_y[2], r_z[2];
    float abs_x, abs_y, abs_z, s_abs_tot;

    for(band=0; band < wrk->nBands; band++){
        echogram_abs = (echogram_data*)wrk->hEchogram_abs[band];

        /* Copy "receiver" echogram data into "absorption" echogram container */
        ims_shoebox_echogramCopy(wrk->hEchogram_rec, wrk->hEchogram_abs[band]);

        /* Then apply the wall absorption onto it, for this band... */

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
            if((echogram_abs->order[i][0]%2)==0) /* ISEVEN */
                abs_x = powf(r_x[0], (float)abs(echogram_abs->order[i][0])/2.0f) * powf(r_x[1], (float)abs(echogram_abs->order[i][0])/2.0f);
            else if (echogram_abs->order[i][0]>0) /* ISODD AND POSITIVE */
                abs_x = powf(r_x[0], ceilf((float)echogram_abs->order[i][0]/2.0f)) * powf(r_x[1], floorf((float)echogram_abs->order[i][0]/2.0f));
            else /* ISODD AND NEGATIVE */
                abs_x = powf(r_x[0], floorf((float)abs(echogram_abs->order[i][0])/2.0f)) * powf(r_x[1], ceilf((float)abs(echogram_abs->order[i][0])/2.0f));

            /* Surfaces intersecting the y-axis */
            if((echogram_abs->order[i][1]%2)==0) /* ISEVEN */
                abs_y = powf(r_y[0], (float)abs(echogram_abs->order[i][1])/2.0f) * powf(r_y[1], (float)abs(echogram_abs->order[i][1])/2.0f);
            else if (echogram_abs->order[i][1]>0) /* ISODD AND POSITIVE */
                abs_y = powf(r_y[0], ceilf((float)echogram_abs->order[i][1]/2.0f)) * powf(r_y[1], floorf((float)echogram_abs->order[i][1]/2.0f));
            else /* ISODD AND NEGATIVE */
                abs_y = powf(r_y[0], floorf((float)abs(echogram_abs->order[i][1])/2.0f)) * powf(r_y[1], ceilf((float)abs(echogram_abs->order[i][1])/2.0f));

            /* Surfaces intersecting the y-axis */
            if((echogram_abs->order[i][2]%2)==0) /* ISEVEN */
                abs_z = powf(r_z[0], (float)abs(echogram_abs->order[i][2])/2.0f) * powf(r_z[1], (float)abs(echogram_abs->order[i][2])/2.0f);
            else if (echogram_abs->order[i][2]>0) /* ISODD AND POSITIVE */
                abs_z = powf(r_z[0], ceilf((float)echogram_abs->order[i][2]/2.0f)) * powf(r_z[1], floorf((float)echogram_abs->order[i][2]/2.0f));
            else /* ISODD AND NEGATIVE */
                abs_z = powf(r_z[0], floorf((float)abs(echogram_abs->order[i][2])/2.0f)) * powf(r_z[1], ceilf((float)abs(echogram_abs->order[i][2])/2.0f));

            /* Apply total absorption */
            s_abs_tot = abs_x * abs_y * abs_z; 
            for(ch=0; ch<echogram_abs->nChannels; ch++)
                echogram_abs->value[ch][i] *= s_abs_tot;
        }
    }
}

void ims_shoebox_renderRIR
(
    void* hWork,
    int fractionalDelayFLAG,
    float fs,
    float** H_filt,
    ims_rir* rir
)
{
    ims_core_workspace *wrk = (ims_core_workspace*)(hWork);
    echogram_data *echogram_abs;
    float* temp;
    int i, j, refl_idx, band, rir_len_samples;
    float endtime, rir_len_seconds;

    /* Render RIR for each octave band */
    for(band=0; band<wrk->nBands; band++){
        echogram_abs = (echogram_data*)wrk->hEchogram_abs[band];

        /* Determine length of rir */
        endtime = echogram_abs->time[echogram_abs->numImageSources-1];
        rir_len_samples = (int)(endtime * fs + 1.0f) + 1; /* ceil + 1 */
        rir_len_seconds = (float)rir_len_samples/fs;

        /* Render rir */
        if(fractionalDelayFLAG){
            // TODO: implement
            saf_print_error("Not implemented yet!");
        }
        else{
            /* Resize RIR vector */
            wrk->rir_bands[band] = (float**)realloc2d((void**)wrk->rir_bands[band], echogram_abs->nChannels, rir_len_samples, sizeof(float));
            wrk->rir_len_samples = rir_len_samples;
            wrk->rir_len_seconds = rir_len_seconds;
            memset(FLATTEN2D(wrk->rir_bands[band]), 0, (echogram_abs->nChannels)*rir_len_samples*sizeof(float)); /* flush */

            /* Accumulate 'values' for each image source */
            for(i=0; i<echogram_abs->numImageSources; i++){
                refl_idx = (int)(echogram_abs->time[i]*fs+0.5f); /* round */
                for(j=0; j<echogram_abs->nChannels; j++)
                    wrk->rir_bands[band][j][refl_idx] += echogram_abs->value[j][i];
            }
        }
    }

    temp = malloc1d((wrk->rir_len_samples+IMS_FIR_FILTERBANK_ORDER)*sizeof(float));

    /* Resize rir->data if needed, then flush with 0s */
    echogram_abs = (echogram_data*)wrk->hEchogram_abs[0];
    if( (echogram_abs->nChannels!=rir->nChannels) || (wrk->rir_len_samples !=rir->length) ){
        rir->data = realloc1d(rir->data, echogram_abs->nChannels * (wrk->rir_len_samples) * sizeof(float));
        rir->length = wrk->rir_len_samples;
        rir->nChannels = echogram_abs->nChannels;
    }
    memset(rir->data, 0, echogram_abs->nChannels * (wrk->rir_len_samples) * sizeof(float));

    /* Apply filterbank to rir_bands and sum them up */
    for(band=0; band<wrk->nBands; band++){
        echogram_abs = (echogram_data*)wrk->hEchogram_abs[band];

        /* Apply the LPF (lowest band), HPF (highest band), and BPF (all other bands) */
        for(j=0; j<echogram_abs->nChannels; j++)
            fftconv(wrk->rir_bands[band][j], H_filt[band], wrk->rir_len_samples, IMS_FIR_FILTERBANK_ORDER+1, 1, temp);

        /* Sum */
        for(i=0; i<echogram_abs->nChannels; i++)
            cblas_saxpy(wrk->rir_len_samples, 1.0f, wrk->rir_bands[band][i], 1, &(rir->data[i*(wrk->rir_len_samples)]), 1);
    }

    free(temp);
}
