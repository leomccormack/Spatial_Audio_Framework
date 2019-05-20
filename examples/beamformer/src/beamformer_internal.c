/*
 Copyright 2019 Leo McCormack
 
 Permission to use, copy, modify, and/or distribute this software for any purpose with or
 without fee is hereby granted, provided that the above copyright notice and this permission
 notice appear in all copies.
 
 THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO
 THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT
 SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR
 ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
 CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
 OR PERFORMANCE OF THIS SOFTWARE.
*/
/*
 * Filename:
 *     beamformer_internal.c
 * Description:
 *     Generates beamformers/virtual microphones in arbitrary directions. Several
 *     different beam pattern types are included.
 * Dependencies:
 *     saf_utilities, afSTFTlib, saf_sh
 * Author, date created:
 *     Leo McCormack, 17.05.2019
 */

#include "beamformer_internal.h"
#define SAF_ENABLE_SOFA_READER   
#include "saf_sofa_reader.h"

void beamformer_initTFT
(
    void* const hBeam
)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);

//    if(pData->hSTFT==NULL){
//        if(pData->new_binauraliseLS)
//            afSTFTinit(&(pData->hSTFT), HOP_SIZE, pData->new_nSH, NUM_EARS, 0, 1);
//        else
//            afSTFTinit(&(pData->hSTFT), HOP_SIZE, pData->new_nSH, pData->new_nLoudpkrs, 0, 1);
//    }
//    else{
//        if(pData->new_binauraliseLS)
//            afSTFTchannelChange(pData->hSTFT, pData->new_nSH, NUM_EARS);
//        else
//            afSTFTchannelChange(pData->hSTFT, pData->new_nSH, pData->new_nLoudpkrs);
//    }
    pData->nBeams = pData->new_nBeams;
    pData->nSH = pData->new_nSH;
}

