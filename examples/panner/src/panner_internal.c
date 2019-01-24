/*
 Copyright 2017-2018 Leo McCormack
 
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
 *     panner_internal.c
 * Description:
 *     A frequency-dependent 3D panner, based on the Vector-base Amplitude Panning (VBAP)
 *     method. Depending on the room, it may be beneficial to utilise amplitude-normalised
 *     gains for low frequencies, and energy-normalised gains for high frequencies; which
 *     this implemenation takes into account with one parameter "DTT". Set "DTT" to 0 for a
 *     normal room, 0.5 for listening room, and 1 for anechoic.
 * Dependencies:
 *     saf_utilities, saf_vbap, afSTFTlib
 * Author, date created:
 *     Leo McCormack, 25.09.2017
 */

#include "panner_internal.h" 

void panner_initGainTables(void* const hPan)
{
    panner_data *pData = (panner_data*)(hPan);
#ifndef FORCE_3D_LAYOUT
    int i;
    float sum_elev;
    
    /* determine dimensionality */
    sum_elev = 0.0f;
    for(i=0; i<pData->nLoudpkrs; i++)
        sum_elev += fabsf(pData->loudpkrs_dirs_deg[i][1]); 
    if(sum_elev < 0.01f)
        pData->output_nDims = 2;
    else
        pData->output_nDims = 3;
#endif
    
    /* generate VBAP gain table */
    if(pData->vbap_gtable!= NULL){
        free(pData->vbap_gtable);
        pData->vbap_gtable = NULL;
    } 
    pData->vbapTableRes[0] = 2;
    pData->vbapTableRes[1] = 5;
#ifdef FORCE_3D_LAYOUT
    pData->output_nDims = 3;
    generateVBAPgainTable3D((float*)pData->loudpkrs_dirs_deg, pData->nLoudpkrs, pData->vbapTableRes[0], pData->vbapTableRes[1], 1, 1, pData->spread_deg,
                            &(pData->vbap_gtable), &(pData->N_vbap_gtable), &(pData->nTriangles));
#else
    if(pData->output_nDims==2)
        generateVBAPgainTable2D((float*)pData->loudpkrs_dirs_deg, pData->nLoudpkrs, pData->vbapTableRes[0],
                                &(pData->vbap_gtable), &(pData->N_vbap_gtable), &(pData->nTriangles));
    else{
        generateVBAPgainTable3D((float*)pData->loudpkrs_dirs_deg, pData->nLoudpkrs, pData->vbapTableRes[0], pData->vbapTableRes[1], 1, 1, pData->spread_deg,
                                &(pData->vbap_gtable), &(pData->N_vbap_gtable), &(pData->nTriangles));
        if(pData->vbap_gtable==NULL){
            /* if generating vbap gain tabled failed, re-calculate with 2D VBAP */
            pData->output_nDims = 2;
            panner_initGainTables(hPan);
        }
    }
#endif
}

void panner_initTFT
(
    void* const hPan
)
{
    panner_data *pData = (panner_data*)(hPan);
    
    if(pData->hSTFT==NULL)
        afSTFTinit(&(pData->hSTFT), HOP_SIZE, pData->new_nSources, pData->new_nLoudpkrs, 0, 1);
    else
        afSTFTchannelChange(pData->hSTFT, pData->new_nSources, pData->new_nLoudpkrs);
    
    pData->nSources = pData->new_nSources;
    pData->nLoudpkrs = pData->new_nLoudpkrs;
}

void panner_loadPreset(PRESETS preset, float dirs_deg[MAX_NUM_INPUTS][2], int* newNCH, int* nDims)
{
    float sum_elev;
    int ch, i, nCH;
    
    switch(preset){
	    default:
        case PRESET_DEFAULT:
            nCH = 1;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = 0.0f;
            break;
#ifdef ENABLE_MONO_PRESET
        case PRESET_MONO:
            nCH = 1;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __mono_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_STEREO_PRESET
        case PRESET_STEREO:
            nCH = 2;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __stereo_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_5PX_PRESET
        case PRESET_5PX:
            nCH = 5;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __5pX_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_7PX_PRESET
        case PRESET_7PX:
            nCH = 7;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __7pX_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_8PX_PRESET
        case PRESET_8PX:
            nCH = 8;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __8pX_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_9PX_PRESET
        case PRESET_9PX:
            nCH = 9;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __9pX_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_10PX_PRESET
        case PRESET_10PX:
            nCH = 10;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __10pX_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_11PX_PRESET
        case PRESET_11PX:
            nCH = 11;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __11pX_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_11PX_7_4_PRESET
        case PRESET_11PX_7_4:
            nCH = 11;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __11pX_7_4_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_13PX_PRESET
        case PRESET_13PX:
            nCH = 13;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __13pX_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_22PX_PRESET
        case PRESET_22PX:
            nCH = 22;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __22pX_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_AALTO_MCC_PRESET
        case PRESET_AALTO_MCC:
            nCH = 29;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Aalto_MCC_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_AALTO_APAJA_PRESET
        case PRESET_AALTO_APAJA:
            nCH = 29;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Aalto_Apaja_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_AALTO_APAJA2_PRESET
        case PRESET_AALTO_APAJA2:
            nCH = 39;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Aalto_Apaja2_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_AALTO_LR_PRESET
        case PRESET_AALTO_LR:
            nCH = 13;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Aalto_LR_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_DTU_AVIL_PRESET
        case PRESET_DTU_AVIL:
            nCH = 64;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __DTU_AVIL_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_T_DESIGN_4_PRESET
        case PRESET_T_DESIGN_4:
            nCH = 4;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_2_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_T_DESIGN_12_PRESET
        case PRESET_T_DESIGN_12:
            nCH = 12;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_4_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_T_DESIGN_24_PRESET
        case PRESET_T_DESIGN_24:
            nCH = 24;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_6_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_T_DESIGN_36_PRESET
        case PRESET_T_DESIGN_36:
            nCH = 36;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_8_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_T_DESIGN_48_PRESET
        case PRESET_T_DESIGN_48:
            nCH = 48;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_9_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_T_DESIGN_60_PRESET
        case PRESET_T_DESIGN_60:
            nCH = 60;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_10_dirs_deg[ch][i];
            break;
#endif
    }
    
    /* Fill remaining slots with default coords */
    for(; ch<MAX_NUM_INPUTS; ch++){
        for(i=0; i<2; i++){
            dirs_deg[ch][i] = default_LScoords64_rad[ch][i]* (180.0f/M_PI);
        }
    }
    
    /* For dynamically changing the number of TFT channels */
    (*newNCH) = nCH;
    
    /* estimate number of dimensions. (Obviously fails if using 2D setups thare are on an angle.
       However, in these cases, triangulation should fail and revert to 2D anyway) */
    sum_elev = 0.0f;
    for(i=0; i<nCH; i++)
        sum_elev += fabsf(dirs_deg[i][1]);
    if(sum_elev < 0.01f)
        (*nDims) = 2;
    else
        (*nDims) = 3;
}

void panner_getPvalue
(
    float DTT,
    float f[HYBRID_BANDS],
    float p[HYBRID_BANDS]
)
{
    int band;
    float a1, a2, p0;
    
    a1 = 0.00045f;
    a2 = 0.000085f;
    
    for(band=0; band<HYBRID_BANDS; band++){
        p0 = 1.5f - 0.5f * cosf(4.7f*tanhf(a1*f[band]))*MAX(0.0f,1.0f-a2*f[band]);
        p[band] = (p0-2.0f)*sqrtf(DTT)+2.0f;
    }
}












 
