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


void ims_shoebox_coreWorkspaceCreate(void** hWork){
    *hWork = malloc1d(sizeof(ims_core_workspace));
    ims_core_workspace *h = (ims_core_workspace*)(*hWork);
    h->d_max = 0.0f;
    h->lengthVec = 0;
    memset(h->room, 0, 3*sizeof(int));
    memset(h->src.v, 0, 3*sizeof(int));
    memset(h->rec.v, 0, 3*sizeof(int));
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
    float c_ms
)
{
    ims_core_workspace *h = (ims_core_workspace*)(hWork);
    int i;
    float d_max;

    d_max = maxTime_s*c_ms;

    if( (h->d_max != d_max) || (h->room[0] != room[0]) || (h->room[1] != room[1]) || (h->room[2] != room[2]) ){
        h->d_max = d_max;
        memcpy(h->room, room, 3*sizeof(int));
        h->Nx = (int)(d_max/(float)room[0] + 1.0f); /* ceil */
        h->Ny = (int)(d_max/(float)room[1] + 1.0f); /* ceil */
        h->Nz = (int)(d_max/(float)room[2] + 1.0f); /* ceil */

        h->lengthVec = (2*(h->Nx)+1) * (2*(h->Ny)+1) * (2*(h->Nz)+1);

        

        i=5555;
    }


}
