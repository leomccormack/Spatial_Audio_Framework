/*
 Copyright 2018 Leo McCormack
 
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
 *     upmix_internal.c
 * Description:
 *     A (soon to be) collection of upmixing algorithms. However, currently, only stereo to
 *     5.x is supported, utilising a modified version of the direct-ambient decomposition
 *     approach described in: Faller, C. (2006). Multiple-loudspeaker playback of stereo
 *     signals. Journal of the Audio Engineering Society, 54(11), 1051-1064.
 * Dependencies:
 *     saf_utilities, afSTFTlib, saf_vbap, saf_sh
 * Author, date created:
 *     Leo McCormack, 04.04.2018
 */

#include "upmix_internal.h"

/* a very lazy low-pass filter: */
//const float __diff_lpf[HYBRID_BANDS] = {
//    1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 0.996801706302619f, 0.987227283375627f, 0.971337974852030f, 0.949235418082441f, 0.921060994002885f, 0.886994922779284f, 0.847255111013416f, 0.802095757884293f, 0.751805729140895f, 0.696706709347165f, 0.637151144198580f, 0.573519986072457f, 0.506220257232778f, 0.435682446276712f, 0.362357754476674f, 0.286715209631955f, 0.209238665891419f, 0.130423708738145f, 0.0507744849335791f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0};

const float __diff_lpf[HYBRID_BANDS] = {
    1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0};

/* grp_idx start from 1 (matlab style), not 0 */
static void groupBands
(
    float freqVector[HYBRID_BANDS],
    float maxFreq,     /* past this frequency the bands are grouped into 1 */
    int USE_ERB_FLAG,  /* 0: group bands using Bark scale widths, 1: group bands using ERB  */
    int** grp_idx,     /* & grouped indices (start from 1); nGrpBands x 1 */
    float** grp_freqs, /* & erb frequencies; nGrpBands x 1 */
    int* nGrpBands     /* & number of grouped bands; 1 x 1 */
)
{
    int band, counter, next_grp_idx;
    float band_centreFreq, grp_f_width, grp_centre, tmp;
    
    if(!USE_ERB_FLAG)
        maxFreq = MIN(maxFreq, 18e3f);
    band_centreFreq = (powf(2.0f, 1.0f/3.0f)+1.0f)/2.0f;
    if((*grp_idx)!=NULL){
        free((void*)(*grp_idx));
        (*grp_idx) = NULL;
    }
    if((*grp_freqs)!=NULL){
        free((void*)(*grp_freqs));
        (*grp_freqs) = NULL;
    }
    (*grp_idx) = malloc(sizeof(int));
    (*grp_freqs) = malloc(sizeof(float));
    (*grp_idx)[0] = 1;
    (*grp_freqs)[0] = freqVector[0];
    counter = 0;
    next_grp_idx = 0;
    while((*grp_freqs)[counter]<maxFreq){
        if (USE_ERB_FLAG) /* ERB scale grouping */
            grp_f_width = 24.7f + 0.108f * (*grp_freqs)[counter] * band_centreFreq;
        else /* Bark scale grouping  */
            grp_f_width = 25.0f + 75.0f * powf(1.0f + 1.4f * powf((((*grp_freqs)[counter] * band_centreFreq)/1e3f),2.0f), 0.69f);
        (*grp_idx) = realloc((*grp_idx), (counter+2)*sizeof(int));
        (*grp_freqs) = realloc((*grp_freqs), (counter+2)*sizeof(float));
        (*grp_freqs)[counter+1] = (*grp_freqs)[counter] + grp_f_width;
        grp_centre = FLT_MAX;
        
        /*  find closest band frequency as upper partition limit */
        for(band=0; band<HYBRID_BANDS; band++){
            tmp =fabsf((*grp_freqs)[counter+1] - freqVector[band]);
            if(tmp <grp_centre){
                grp_centre = tmp;
                next_grp_idx = band;
            }
        }
        (*grp_idx)[counter+1] = next_grp_idx + 1;
        if((*grp_idx)[counter+1] == (*grp_idx)[counter])
            (*grp_idx)[counter+1] = (*grp_idx)[counter+1]+1;
        (*grp_freqs)[counter+1] = freqVector[(*grp_idx)[counter+1]-1];
        counter++;
    }
    
    /* set last limit as the last band */
    (*grp_idx) = realloc((*grp_idx), (counter + 2) * sizeof(int));
    (*grp_freqs) = realloc((*grp_freqs), (counter + 2) * sizeof(float));
    (*grp_idx)[counter+1] = HYBRID_BANDS;
    (*grp_freqs)[counter+1] = freqVector[HYBRID_BANDS-1];
    (*nGrpBands) = counter+2;
}


void upmix_initCodec
(
    void* const hUpmx
)
{
    upmix_data *pData = (upmix_data*)(hUpmx);
    codecPars* pars = pData->pars;
    int i, j;
    
    /* generate VBAP gain table for the grid */
    pars->vbap_azi_res = 1;
    pData->nLoudspeakers = MAX_NUM_OUTPUT_CHANNELS;
    for(i=0; i<pData->nLoudspeakers; i++)
        for(j=0; j<2; j++)
            pData->loudpkrs_dirs_deg[i][j] = __5pX_dirs_deg[i][j]; /* only stereo to 5.x is currently supported */
    free(pars->grid_vbap_gtable);
    generateVBAPgainTable2D((float*)pData->loudpkrs_dirs_deg, pData->nLoudspeakers, pars->vbap_azi_res , &(pars->grid_vbap_gtable), &(pars->grid_N_vbap_gtable), &(pars->grid_nPairs));
    
    /* define band grouping */
    pars->maxGrpFreq = MAX_GROUP_FREQ;
    groupBands(pData->freqVector, pars->maxGrpFreq, 1, &(pars->grp_idx), &(pars->grp_freqs), &(pars->nGrpBands));
    
    /* low-pass filter */
    memcpy(pars->diff_lpf, __diff_lpf, HYBRID_BANDS*sizeof(float));
    
    /* for averaging DoA estimate over time */
    free(pars->prev_est_dir);
    pars->prev_est_dir = calloc(pars->nGrpBands,sizeof(float));
}

















 
