/*
 * Copyright 2017-2018 Leo McCormack
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
 * @file saf_vbap.c
 * @brief Public part of the "saf_vbap" module
 *
 * VBAP functions largely derived from the MATLAB library by Archontis Politis,
 * found in [1].
 *
 * @see [1] https://github.com/polarch/Vector-Base-Amplitude-Panning
 *
 * @author Leo McCormack
 * @date 02.10.2017
 */

#include "saf_vbap.h"
#include "saf_vbap_internal.h" 

#define ENABLE_VBAP_DEBUGGING_CODE 0
#if ENABLE_VBAP_DEBUGGING_CODE
#define SAVE_PATH "../faces.txt"
#define SAVE_PATH2 "../vbapGains_compressed.txt"
#define SAVE_PATH3 "../vbapGains_table.txt"
#endif

void generateVBAPgainTable3D_srcs
(
    float* src_dirs_deg,
    int S,
    float* ls_dirs_deg,
    int L,
    int omitLargeTriangles,
    int enableDummies,
    float spread,
    float** gtable /* &: S x L */,
    int* N_gtable /* & S */,
    int* nTriangles
)
{
    int N_points, numOutVertices, numOutFaces;
    int* out_faces;
    float *out_vertices, *layoutInvMtx;
    int i, L_d;
    int needDummy[2] = {1, 1};
    float* ls_dirs_d_deg;
    
    /* find loudspeaker triangles */
    out_vertices = NULL;
    out_faces = NULL;
    if(enableDummies){
        /* scan the loudspeaker directions to see if dummies need to be added */
        for(i=0; i<L; i++){
            if(ls_dirs_deg[i*2+1] <= -ADD_DUMMY_LIMIT)
                needDummy[0] = 0;
            if(ls_dirs_deg[i*2+1] >=  ADD_DUMMY_LIMIT)
                needDummy[1] = 0;
        }
        if(needDummy[0] || needDummy[1]){
            /* add dummies to the extreme top/bottom as required */
            L_d = L+needDummy[0]+needDummy[1];
            ls_dirs_d_deg = malloc1d(L_d*2*sizeof(float));
            memcpy(ls_dirs_d_deg, ls_dirs_deg, L*2*sizeof(float));
            if (needDummy[0]){
                ls_dirs_d_deg[i*2+0] = 0.0f;
                ls_dirs_d_deg[i*2+1] = -90.0f;
                i++;
            }
            if (needDummy[1]){
                ls_dirs_d_deg[i*2+0] = 0.0f;
                ls_dirs_d_deg[i*2+1] = 90.0f;
            }
            
            /* triangulate while including the dummy loudspeaker directions */
            findLsTriplets(ls_dirs_d_deg, L_d, omitLargeTriangles, &out_vertices, &numOutVertices, &out_faces, &numOutFaces);
            free(ls_dirs_d_deg);
        }
        else /* triangulate as normal */
            findLsTriplets(ls_dirs_deg, L, omitLargeTriangles, &out_vertices, &numOutVertices, &out_faces, &numOutFaces);
    }
    else
        findLsTriplets(ls_dirs_deg, L, omitLargeTriangles, &out_vertices, &numOutVertices, &out_faces, &numOutFaces);
#if ENABLE_VBAP_DEBUGGING_CODE
    /* save faces and vertices for verification in matlab: */
    FILE* objfile = fopen(SAVE_PATH, "wt");
    fprintf(objfile, "faces = [\n");
    for (i = 0; i < numOutFaces; i++) {
        fprintf(objfile, " %u, %u, %u;\n",
                out_faces[3*i+0],
                out_faces[3*i+1],
                out_faces[3*i+2]);
    }
    fprintf(objfile, "];\n\n\n");
    fprintf(objfile, "vert = [\n");
    for (i = 0; i < numOutVertices; i++) {
        fprintf(objfile, " %f, %f, %f;\n",
                out_vertices[3*i+0],
                out_vertices[3*i+1],
                out_vertices[3*i+2]);
    }
    fprintf(objfile, "];\n\n\n");
    fclose(objfile);
#endif
    
    /* Invert matrix */
    layoutInvMtx = NULL;
    invertLsMtx3D(out_vertices, out_faces, numOutFaces, &layoutInvMtx);
    
    /* Calculate VBAP gains for each source position */
    N_points = S;
    vbap3D(src_dirs_deg, N_points, numOutVertices, out_faces, numOutFaces, spread, layoutInvMtx, gtable);
    if(enableDummies){
        if(needDummy[0] || needDummy[1]){
            /* remove the gains for the dummy loudspeakers, they have served their purpose and can now be laid to rest */
            for(i=0; i<N_points; i++)
                memcpy(&(*gtable)[i*L], &(*gtable)[i*numOutVertices], L*sizeof(float));
            (*gtable) = realloc((*gtable), N_points*L*sizeof(float));
        }
    }
    
    /* output */
    (*N_gtable) = N_points;
    (*nTriangles) = numOutFaces;
#if ENABLE_VBAP_DEBUGGING_CODE
    /* save gain table for verification in matlab: */
    FILE* objfile2 = fopen(SAVE_PATH3, "wt");
    fprintf(objfile2, "vbap_gtable = [\n");
    for (i = 0; i < N_points; i++) {
        for (j = 0; j < L; j++) {
            fprintf(objfile2, " %f", (*gtable)[3*L+j] );
            if(j<L-1)
                fprintf(objfile2,",");
        }
        fprintf(objfile2, ";\n");
    }
    fprintf(objfile2, "];\n\n\n");
    fclose(objfile2);
#endif
    
    /* clean up */
    free1d((void**)&(out_vertices));
    free1d((void**)&(out_faces));
    free1d((void**)&(layoutInvMtx));
}

void generateVBAPgainTable3D
(
    float* ls_dirs_deg,
    int L,
    int az_res_deg,
    int el_res_deg,
    int omitLargeTriangles,
    int enableDummies,
    float spread,
    float** gtable /* N_srcs x N_lspkrs  */,
    int* N_gtable,
    int* nTriangles
)
{
    int i, j, N_azi, N_ele, N_points, numOutVertices, numOutFaces;
    int* out_faces;
    float fi;
    float* azi, *ele, *src_dirs, *out_vertices, *layoutInvMtx;
    int L_d;
    int needDummy[2] = {1, 1};
    float* ls_dirs_d_deg;
    
    /* compute source directions for the grid */
    N_azi = (int)((360.0f/(float)az_res_deg) + 1.5f);
    N_ele = (int)((180.0f/(float)el_res_deg) + 1.5f);
    azi = malloc1d(N_azi * sizeof(float));
    ele = malloc1d(N_ele * sizeof(float));
    for(fi = -180.0f, i = 0; i<N_azi; fi+=(float)az_res_deg, i++)
        azi[i] = fi;
    for(fi = -90.0f,  i = 0; i<N_ele; fi+=(float)el_res_deg, i++)
        ele[i] = fi;
    src_dirs = malloc1d((N_azi*N_ele)*2*sizeof(float));
    for(i = 0; i<N_ele; i++){
        for(j=0; j<N_azi; j++){
            src_dirs[(i*N_azi + j)*2]   = azi[j];
            src_dirs[(i*N_azi + j)*2+1] = ele[i];
        }
    }
    
    /* find loudspeaker triangles */
    out_vertices = NULL;
    out_faces = NULL;
    if(enableDummies){
        /* scan the loudspeaker directions to see if dummies need to be added */
        for(i=0; i<L; i++){
            if(ls_dirs_deg[i*2+1] <= -ADD_DUMMY_LIMIT)
                needDummy[0] = 0;
            if(ls_dirs_deg[i*2+1] >=  ADD_DUMMY_LIMIT)
                needDummy[1] = 0;
        }
        if(needDummy[0] || needDummy[1]){
            /* add dummies to the extreme top/bottom as required */
            L_d = L+needDummy[0]+needDummy[1];
            ls_dirs_d_deg = malloc1d(L_d*2*sizeof(float));
            memcpy(ls_dirs_d_deg, ls_dirs_deg, L*2*sizeof(float));
            if (needDummy[0]){
                ls_dirs_d_deg[i*2+0] = 0.0f;
                ls_dirs_d_deg[i*2+1] = -90.0f;
                i++;
            }
            if (needDummy[1]){
                ls_dirs_d_deg[i*2+0] = 0.0f;
                ls_dirs_d_deg[i*2+1] = 90.0f;
            }
            
            /* triangulate while including the dummy loudspeakers */
            findLsTriplets(ls_dirs_d_deg, L_d, omitLargeTriangles, &out_vertices, &numOutVertices, &out_faces, &numOutFaces);
            free(ls_dirs_d_deg);
        }
        else /* triangulate as normal */
            findLsTriplets(ls_dirs_deg, L, omitLargeTriangles, &out_vertices, &numOutVertices, &out_faces, &numOutFaces);
    }
    else /* triangulate as normal */
        findLsTriplets(ls_dirs_deg, L, omitLargeTriangles, &out_vertices, &numOutVertices, &out_faces, &numOutFaces);
#if ENABLE_VBAP_DEBUGGING_CODE
    /* save faces and vertices for verification in matlab: */
    FILE* objfile = fopen(SAVE_PATH, "wt");
    fprintf(objfile, "faces = [\n");
    for (i = 0; i < numOutFaces; i++) {
        fprintf(objfile, " %u, %u, %u;\n",
                out_faces[3*i+0],
                out_faces[3*i+1],
                out_faces[3*i+2]);
    }
    fprintf(objfile, "];\n\n\n");
    fprintf(objfile, "vert = [\n");
    for (i = 0; i < numOutVertices; i++) {
        fprintf(objfile, " %f, %f, %f;\n",
                out_vertices[3*i+0],
                out_vertices[3*i+1],
                out_vertices[3*i+2]);
    }
    fprintf(objfile, "];\n\n\n");
    fclose(objfile);
#endif
    
    /* Invert matrix */
    layoutInvMtx = NULL;
    invertLsMtx3D(out_vertices, out_faces, numOutFaces, &layoutInvMtx);
    
    /* Calculate VBAP gains for each source position */
    N_points = N_azi*N_ele;
    vbap3D(src_dirs, N_points, numOutVertices, out_faces, numOutFaces, spread, layoutInvMtx,  gtable);
    
    /* remove the gains for the dummy loudspeakers, they have served their purpose and can now be laid to rest */
    if(enableDummies){
        if(needDummy[0] || needDummy[1]){
            for(i=0; i<N_points; i++)
                memcpy(&(*gtable)[i*L], &(*gtable)[i*numOutVertices], L*sizeof(float));
            (*gtable) = realloc((*gtable), N_points*L*sizeof(float));
        }
    }
    
    /* output */
    (*N_gtable) = N_points;
    (*nTriangles) = numOutFaces;
#if ENABLE_VBAP_DEBUGGING_CODE
    /* save gain table for verification in matlab: */
    FILE* objfile2 = fopen(SAVE_PATH3, "wt");
    fprintf(objfile2, "vbap_gtable = [\n");
    for (i = 0; i < N_points; i++) {
        for (j = 0; j < L; j++) {
            fprintf(objfile2, " %f", (*gtable)[i*L+j] );
            if(j<L-1)
                fprintf(objfile2,",");
        }
        fprintf(objfile2, ";\n");
    }
    fprintf(objfile2, "];\n\n\n");
    fclose(objfile2);
#endif
    
    /* clean up */
    free1d((void**)&(out_vertices));
    free1d((void**)&(out_faces));
    free1d((void**)&(layoutInvMtx));
    free(azi);
    free(ele);
}

void compressVBAPgainTable3D
(
    float* vbap_gtable,
    int nTable,
    int nDirs,
    float* vbap_gtableComp, /* nTable x 3  */
    int* vbap_gtableIdx     /* nTable x 3  */
)
{
    int i, j, nt;
    int idx_nt[3];
    float gains_nt[3];
    float gains_sum;
    
    memset(vbap_gtableComp, 0, nTable*3*sizeof(float));
    memset(vbap_gtableIdx, 0, nTable*3*sizeof(int));
    
    /* compress table by keeping only the non-zero gains and their indices, and also convert to AMPLITUDE NORMALISED */
    for(nt=0; nt<nTable; nt++){
        gains_sum = 0.0f;
        for(i=0, j=0; i<nDirs; i++){
            if(vbap_gtable[nt*nDirs+i]>0.0000001f){
                gains_nt[j] = vbap_gtable[nt*nDirs+i];
                gains_sum += gains_nt[j];
                idx_nt[j] = i;
                j++;
            }
        }
        //assert(j<4);
        for(i=0; i<j; i++){
            vbap_gtableComp[nt*3+i] = MAX(gains_nt[i]/gains_sum, 0.0f);
            vbap_gtableIdx[nt*3+i] = idx_nt[i];
        }
    }
#if ENABLE_VBAP_DEBUGGING_CODE
    /* save gain table for verification in matlab: */
    FILE* objfile = fopen(SAVE_PATH2, "wt");
    fprintf(objfile, "vbap_gtableComp = [\n");
    for (i = 0; i < nTable; i++) {
        fprintf(objfile, " %f, %f, %f;\n",
                vbap_gtableComp[3*i+0],
                vbap_gtableComp[3*i+1],
                vbap_gtableComp[3*i+2]);
    }
    fprintf(objfile, "];\n\n\n");
    fprintf(objfile, "vbap_gtableIdx = [\n");
    for (i = 0; i < nTable; i++) {
        fprintf(objfile, " %u, %u, %u;\n",
                vbap_gtableIdx[3*i+0],
                vbap_gtableIdx[3*i+1],
                vbap_gtableIdx[3*i+2]);
    }
    fprintf(objfile, "];\n\n\n");
    fclose(objfile);
#endif
}

void VBAPgainTable2InterpTable
(
    float* vbap_gtable,
    int nTable,
    int nDirs
)
{
    int i, j;
    float* gains_sum;
    
    gains_sum = calloc1d(nTable,sizeof(float));
    for(i=0; i<nTable; i++)
        for(j=0; j<nDirs; j++)
            gains_sum[i] += vbap_gtable[i*nDirs+j];
    for(i=0; i<nTable; i++)
        for(j=0; j<nDirs; j++)
            vbap_gtable[i*nDirs+j] /= gains_sum[i];
    
    free(gains_sum);
}

void generateVBAPgainTable2D_srcs
(
    float* src_dirs_deg,
    int S,
    float* ls_dirs_deg,
    int L,
    float** gtable /* S x L */,
    int* N_gtable,
    int* nPairs
)
{
    int i, N_points, numOutPairs;
    int* out_pairs;
    float *layoutInvMtx, *ls_vertices;
    
    out_pairs = NULL;
    findLsPairs(ls_dirs_deg, L, &out_pairs, &numOutPairs);
    ls_vertices = malloc1d(L*2*sizeof(float));
    for(i=0; i<L; i++){
        ls_vertices[i*2+0] = cosf(ls_dirs_deg[i*2]*M_PI/180.0f);
        ls_vertices[i*2+1] = sinf(ls_dirs_deg[i*2]*M_PI/180.0f);
    }
    
    /* Invert matrix */
    layoutInvMtx = NULL;
    invertLsMtx2D(ls_vertices, out_pairs, numOutPairs, &layoutInvMtx);
    
    /* Calculate VBAP gains for each source position */
    N_points = S;
    vbap2D(src_dirs_deg, N_points, L, out_pairs, numOutPairs, layoutInvMtx,  gtable);
    (*nPairs) = numOutPairs;
    (*N_gtable) = N_points;
    
    free(ls_vertices);
    free1d((void**)&(out_pairs));
    free1d((void**)&(layoutInvMtx));
}

void generateVBAPgainTable2D
(
    float* ls_dirs_deg,
    int L,
    int az_res_deg,
    float** gtable /* N_gtable x L  */,
    int* N_gtable,
    int* nPairs
)
{
    int i, N_azi, N_points, numOutPairs;
    int* out_pairs;
    float fi;
    float *src_dirs, *layoutInvMtx, *ls_vertices;
    
    /* compute directions for the grid */
    N_azi = (int)((360.0f/(float)az_res_deg) + 1.5f);
    src_dirs = malloc1d(N_azi * sizeof(float));
    for(fi = -180.0f, i = 0; i<N_azi; fi+=(float)az_res_deg, i++)
        src_dirs[i] = fi;
    out_pairs = NULL;
    findLsPairs(ls_dirs_deg, L, &out_pairs, &numOutPairs);
    ls_vertices = malloc1d(L*2*sizeof(float));
    for(i=0; i<L; i++){
        ls_vertices[i*2+0] = cosf(ls_dirs_deg[i*2]*M_PI/180.0f);
        ls_vertices[i*2+1] = sinf(ls_dirs_deg[i*2]*M_PI/180.0f);
    }
    
    /* Invert matrix */
    layoutInvMtx = NULL;
    invertLsMtx2D(ls_vertices, out_pairs, numOutPairs, &layoutInvMtx);
    
    /* Calculate VBAP gains for each source position */
    N_points = N_azi;
    vbap2D(src_dirs, N_points, L, out_pairs, numOutPairs, layoutInvMtx,  gtable);
    (*nPairs) = numOutPairs;
    (*N_gtable) = N_points;
    
    free(ls_vertices);
    free(src_dirs);
    free1d((void**)&(out_pairs));
    free1d((void**)&(layoutInvMtx));
}

/* Laitinen, M., Vilkamo, J., Jussila, K., Politis, A., Pulkki, V. (2014). Gain
 * normalisation in amplitude panning as a function of frequency and room
 * reverberance. 55th International Conference of the AES. Helsinki, Finland. */
void getPvalues
(
    float DTT,
    float* freq,
    int nFreq,
    float* pValues
)
{
    int i;
    float a1, a2, p0;
    
    a1 = 0.00045f;
    a2 = 0.000085f;
    for(i=0; i<nFreq; i++){
        p0 = 1.5f - 0.5f * cosf(4.7f*tanhf(a1*freq[i])) * MAX(0.0f, 1.0f-a2*freq[i]);
        pValues[i] = (p0-2.0f)*sqrtf(DTT)+2.0f;
    } 
}
