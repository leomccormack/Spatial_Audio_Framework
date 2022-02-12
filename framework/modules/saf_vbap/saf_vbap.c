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
 * @ingroup VBAP
 * @brief Public source for the VBAP/MDAP module (#SAF_VBAP_MODULE)
 *
 * VBAP functions largely derived from the MATLAB library found in [1].
 *
 * @see [1] https://github.com/polarch/Vector-Base-Amplitude-Panning
 *          Copyright (c) 2015, Archontis Politis, BSD-3-Clause License
 *
 * @author Leo McCormack
 * @date 02.10.2017
 * @license ISC
 */

#include "saf_vbap.h"
#include "saf_vbap_internal.h" 

/**
 * Enable to also export VBAP gain tables to a .txt file (for debugging
 * purposes) */
#define ENABLE_VBAP_DEBUGGING_CODE 0
#if ENABLE_VBAP_DEBUGGING_CODE
/** Save file path for exporting the face indices (to verify triangulation) */
# define SAVE_PATH "../faces.txt"
/** Save file path for exporting the compressed VBAP gain table */
# define SAVE_PATH2 "../vbapGains_compressed.txt"
/** Save file path for exporting the computed VBAP gain table */
# define SAVE_PATH3 "../vbapGains_table.txt"
#endif

/* ========================================================================== */
/*                               Misc. Functions                              */
/* ========================================================================== */

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
                memmove(&(*gtable)[i*L], &(*gtable)[i*numOutVertices], L*sizeof(float));
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
    free(out_vertices);
    free(out_faces);
    free(layoutInvMtx);
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
                memmove(&(*gtable)[i*L], &(*gtable)[i*numOutVertices], L*sizeof(float));
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
    free(out_vertices);
    free(out_faces);
    free(layoutInvMtx);
    free(src_dirs);
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
        //saf_assert(j<4);
        for(i=0; i<j; i++){
            vbap_gtableComp[nt*3+i] = SAF_MAX(gains_nt[i]/gains_sum, 0.0f);
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
        ls_vertices[i*2+0] = cosf(ls_dirs_deg[i*2]*SAF_PI/180.0f);
        ls_vertices[i*2+1] = sinf(ls_dirs_deg[i*2]*SAF_PI/180.0f);
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
    free(out_pairs);
    free(layoutInvMtx);
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
        ls_vertices[i*2+0] = cosf(ls_dirs_deg[i*2]*SAF_PI/180.0f);
        ls_vertices[i*2+1] = sinf(ls_dirs_deg[i*2]*SAF_PI/180.0f);
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
    free(out_pairs);
    free(layoutInvMtx);
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
        p0 = 1.5f - 0.5f * cosf(4.7f*tanhf(a1*freq[i])) * SAF_MAX(0.0f, 1.0f-a2*freq[i]);
        pValues[i] = (p0-2.0f)*sqrtf(DTT)+2.0f;
    } 
}


/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

void findLsTriplets
(
    float* ls_dirs_deg,
    int L,
    int omitLargeTriangles,
    float** out_vertices,
    int* numOutVertices,
    int** out_faces,
    int* numOutFaces
)
{
    int i, j, numValidFaces, nFaces;
    int* validFacesID, *valid_faces, *valid_faces2, *faces;
    float dotcc, aperture_lim;
    float vecs[3][3], cvec[3], centroid[3], a[3], b[3], abc[3];
    CH_FLOAT  rcoselev;
    ch_vertex* vertices;

    /* Build the convex hull for the points on the sphere - in this special case the
       result equals the Delaunay triangulation of the points */
    vertices = malloc1d(L*sizeof(ch_vertex));
    (*numOutVertices)  = L;
    (*out_vertices) = (float*)malloc1d(L*3*sizeof(float));
    for ( i = 0; i < L; i++) {
        (*out_vertices)[i*3+2] = (float)((CH_FLOAT)sin((double)ls_dirs_deg[i*2+1]* SAF_PId/180.0));
        rcoselev = (CH_FLOAT)cos((double)ls_dirs_deg[i*2+1]*SAF_PId/180.0);
        (*out_vertices)[i*3+0] = (float)(rcoselev * (CH_FLOAT)cos((double)ls_dirs_deg[i*2+0]*SAF_PId/180.0));
        (*out_vertices)[i*3+1] = (float)(rcoselev * (CH_FLOAT)sin((double)ls_dirs_deg[i*2+0]*SAF_PId/180.0));
        vertices[i].x = (*out_vertices)[i*3+0];
        vertices[i].y = (*out_vertices)[i*3+1];
        vertices[i].z = (*out_vertices)[i*3+2];
    }
    faces = NULL;
    convhull_3d_build(vertices, L, &faces, NULL, NULL, &nFaces);
#ifndef NDEBUG
    if(faces==NULL)
        saf_print_error("Failed to compute the Convex Hull of the specified vertices.");
#endif

#if 0
    int k, minIntVal, minIdx, nFaces;
    int tmp3[3], circface[3];

    /* circularily shift the indices to start from lowest value */
    for(i=0; i<nFaces; i++){
        minIntVal = L;
        minIdx = 0;
        for(j=0; j<3; j++){
            if (faces[i*3+j] < minIntVal){
                minIntVal = faces[i*3+j];
                minIdx=j;
            }
        }
        for(j=minIdx, k=0; j<minIdx+3; j++, k++)
            circface[k] = faces[i*3+(j % 3)];
        for(j=0; j<3; j++)
            faces[i*3+j] = circface[j];
    }

    /* sort indices in accending order for the first dimension */
    for (j=0; j<nFaces - 1; j++)  {
        for (i=0; i<nFaces - 1; i++) {
            if (faces[(i+1)*3+0] < faces[i*3+0])  {
                for(k=0; k<3; k++){
                    tmp3[k] = faces[i*3+k];
                    faces[i*3+k] = faces[(i+1)*3+k];
                    faces[(i+1)*3+k] = tmp3[k];
                }
            }
        }
    }

    /* sort indices in accending order for the second dimension */
    for (j=0; j<nFaces - 1; j++)  {
        for (i=0; i<nFaces - 1; i++) {
            if ( (faces[(i+1)*3+1] < faces[i*3+1]) && (faces[(i+1)*3+0] == faces[i*3+0]) ) {
                for(k=0; k<3; k++){
                    tmp3[k] = faces[i*3+k];
                    faces[i*3+k] = faces[(i+1)*3+k];
                    faces[(i+1)*3+k] = tmp3[k];
                }
            }
        }
    }
#endif

    /* Omit triplets if their normals and the centroid to the triplets have an angle larger than pi/2 */
    numValidFaces = 0;
    validFacesID = malloc1d(nFaces*sizeof(int));
    for(i=0; i<nFaces; i++){
        for(j=0; j<3; j++){
            vecs[0][j] = (*out_vertices)[faces[i*3+0]*3+j];
            vecs[1][j] = (*out_vertices)[faces[i*3+1]*3+j];
            vecs[2][j] = (*out_vertices)[faces[i*3+2]*3+j];
        }
        for(j=0; j<3; j++){
            a[j] = vecs[1][j]-vecs[0][j];
            b[j] = vecs[2][j]-vecs[1][j];
        }
        ccross(a, b, cvec);
        for(j=0; j<3; j++)
            centroid[j] = (vecs[0][j] + vecs[1][j] + vecs[2][j])/3.0f;
        dotcc = cvec[0] * centroid[0] + cvec[1] * centroid[1] + cvec[2] * centroid[2];
        if(acosf(SAF_MAX(SAF_MIN(dotcc,0.99999999f),-0.99999999f)/* avoids complex numbers */)<(SAF_PI/2.0f)){
            validFacesID[i] = 1;
            numValidFaces++;
        }
        else{
            validFacesID[i] = 0;
        }
    }
    valid_faces = malloc1d(numValidFaces*3*sizeof(int));
    for(i=0, j=0; i<nFaces; i++){
        if(validFacesID[i]==1){
            valid_faces[j*3+0] = faces[i*3+0];
            valid_faces[j*3+1] = faces[i*3+1];
            valid_faces[j*3+2] = faces[i*3+2];
            j++;
        }
    }
    free(validFacesID);

    /* Omit Trianges that have an aperture larger than APERTURE_LIMIT_DEG */
    valid_faces2 = NULL;
    if(omitLargeTriangles) {
        aperture_lim = APERTURE_LIMIT_DEG * SAF_PI/180.0f;
        nFaces = numValidFaces;
        numValidFaces = 0;
        validFacesID = malloc1d(nFaces*sizeof(int));
        for(i=0; i<nFaces; i++){
            for(j=0; j<3; j++){
                vecs[0][j] = (*out_vertices)[valid_faces[i*3+0]*3+j];
                vecs[1][j] = (*out_vertices)[valid_faces[i*3+1]*3+j];
                vecs[2][j] = (*out_vertices)[valid_faces[i*3+2]*3+j];
            }
            dotcc = vecs[0][0] * vecs[1][0] + vecs[0][1] * vecs[1][1] + vecs[0][2] * vecs[1][2];
            abc[0] = acosf(dotcc);
            dotcc = vecs[1][0] * vecs[2][0] + vecs[1][1] * vecs[2][1] + vecs[1][2] * vecs[2][2];
            abc[1] = acosf(dotcc);
            dotcc = vecs[2][0] * vecs[0][0] + vecs[2][1] * vecs[0][1] + vecs[2][2] * vecs[0][2];
            abc[2] = acosf(dotcc);
            if(abc[0]<aperture_lim && abc[1]<aperture_lim && abc[2]<aperture_lim){
                validFacesID[i] = 1;
                numValidFaces++;
            }
            else{
                validFacesID[i] = 0;
            }
        }
        valid_faces2 = malloc1d(numValidFaces*3*sizeof(int));
        for(i=0, j=0; i<nFaces; i++){
            if(validFacesID[i]==1){
                valid_faces2[j*3+0] = valid_faces[i*3+0];
                valid_faces2[j*3+1] = valid_faces[i*3+1];
                valid_faces2[j*3+2] = valid_faces[i*3+2];
                j++;
            }
        }
        free(validFacesID);
    }

    /* output valid faces */
    (*numOutFaces) = numValidFaces;
    (*out_faces) = (int*)malloc1d(numValidFaces*3*sizeof(int));
    if(omitLargeTriangles)
        memcpy((*out_faces), valid_faces2, numValidFaces*3*sizeof(int));
    else
        memcpy((*out_faces), valid_faces, numValidFaces*3*sizeof(int));

    /* clean-up */
    free(faces);
    free(vertices);
    free(valid_faces);
    if(omitLargeTriangles)
        free(valid_faces2);
}

void invertLsMtx3D
(
    float* U_spkr /* L x 3 */,
    int* ls_groups /* N_group x 3 */,
    int N_group,
    float** layoutInvMtx/* N_group x 9 */
)
{
    int i, j, n;
    float tempGroup[9];
    float tempInv[9];
    void* hSinv;

    /* pre-calculate inversions of the loudspeaker groups and store into matrix */
    (*layoutInvMtx) = malloc1d(N_group * 9 * sizeof(float));
    utility_sinv_create(&hSinv, 3);
    for(n=0; n<N_group; n++){
        /* get the unit vectors for the current group */
        for(i=0; i<3; i++)
            for(j=0; j<3; j++)
                tempGroup[j*3+i] = U_spkr[ls_groups[n*3+i]*3 + j]; /* ^T */

        /* get inverse of current group */
        utility_sinv(hSinv, tempGroup, tempInv, 3);

        /* store the vectorised inverse as a row in the output */
        cblas_scopy(9, tempInv, 1, (*layoutInvMtx) + n*9, 1);
    }
    utility_sinv_destroy(&hSinv);
}

void getSpreadSrcDirs3D
(
    float src_azi_rad,
    float src_elev_rad,
    float spread,
    int num_src,
    int num_rings_3d,
    float* U_spread
)
{
    int i, j, ns, nr;
    float theta, sin_theta, cos_theta, scale, spread_rad, ring_rad, U_spread_norm;
    float u[3], u_x_u[3][3], u_x[3][3], R_theta[3][3], uu2[3], spreadbase_ns[3];
    float* spreadbase;

    /* rotation matrix using the axis of rotation-angle definition (around source direction) */
    u[0] = cosf(src_elev_rad) * cosf(src_azi_rad);
    u[1] = cosf(src_elev_rad) * sinf(src_azi_rad);
    u[2] = sinf(src_elev_rad);
    u_x_u[0][0] = powf(u[0],2.0f); u_x_u[0][1] = u[0]*u[1]; u_x_u[0][2] = u[0]*u[2];
    u_x_u[1][0] = u[0]*u[1]; u_x_u[1][1] = powf(u[1],2.0f); u_x_u[1][2] = u[1]*u[2];
    u_x_u[2][0] = u[0]*u[2]; u_x_u[2][1] = u[1]*u[2]; u_x_u[2][2] = powf(u[2],2.0f);
    u_x[0][0] = 0.0f; u_x[0][1] = -u[2]; u_x[0][2] = u[1];
    u_x[1][0] = u[2]; u_x[1][1] = 0.0f; u_x[1][2] = -u[0];
    u_x[2][0] = -u[1]; u_x[2][1] = u[0]; u_x[2][2] = 0.0f;
    theta = 2.0f*SAF_PI/(float)num_src;
    sin_theta = sinf(theta);
    cos_theta = cosf(theta);
    for(i=0; i<3; i++)
        for(j=0; j<3; j++)
            R_theta[i][j] = sin_theta*u_x[i][j] + (1.0f-cos_theta)*u_x_u[i][j] + (i==j ? cos_theta : 0.0f);

    /*  create a ring of sources on the plane that is purpendicular to the source directions */
    spreadbase = calloc1d(num_src*3, sizeof(float));
    if ((src_elev_rad > SAF_PI/2.0f-0.01f ) || (src_elev_rad<-(SAF_PI/2.0f-0.01f)))
        spreadbase[0] = 1.0f;
    else{
        const float u2[3] = {0.0f, 0.0f, 1.0f};
        ccross(u, (float*)u2, uu2);
        scale = 0.0f;
        for(i=0; i<3; i++)
            scale += powf(uu2[i],2.0f);
        scale = sqrtf(scale);
        for(i=0; i<3; i++)
            spreadbase[i] = uu2[i]/scale;
    }

    /* get ring of directions by rotating the first vector around the source */
    for (ns = 1; ns<num_src; ns++){
        for(i=0; i<3; i++)
            spreadbase_ns[i] = spreadbase[(ns-1)*3+i];
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 1, 3, 1.0f,
                    (const float*)R_theta, 3,
                    spreadbase_ns, 1, 0.0f,
                    &spreadbase[ns*3], 1);
    }

    /* squeeze the perpendicular ring to the desired spread */
    spread_rad = (spread/2.0f)*SAF_PI/180.0f;
    ring_rad = spread_rad/(float)num_rings_3d;
    memset(U_spread, 0, num_rings_3d*num_src*3*sizeof(float));
    for(nr=0; nr<num_rings_3d; nr++)
        for (ns = 0; ns<num_src; ns++)
            for(i=0; i<3; i++)
                U_spread[(nr*num_src + ns)*3 + i] = u[i] + spreadbase[ns*3+i]*tanf(ring_rad*(float)(nr+1));

    /* normalise vectors to unity (based on first vector) */
    U_spread_norm = sqrtf(powf(U_spread[0],2.0f) + powf(U_spread[1],2.0f) + powf(U_spread[2],2.0f));
    for(i=0; i<num_rings_3d*num_src*3; i++)
        U_spread[i] /= U_spread_norm;

    /* append the original source direction at the end */
    for(i=0; i<3; i++)
        U_spread[(num_rings_3d*num_src)*3 + i] = u[i];

    free(spreadbase);
}


void vbap3D
(
    float* src_dirs,
    int src_num,
    int ls_num,
    int* ls_groups,
    int nFaces,
    float spread,
    float* layoutInvMtx,
    float** GainMtx
)
{
    int i, j, ns, nspr;
    float azi_rad, elev_rad, min_val, g_tmp_rms, gains_rms;
    float u[3], g_tmp[3], ls_invMtx_s[3];
    float* gains;

    (*GainMtx) = malloc1d(src_num*ls_num*sizeof(float));
    gains = malloc1d(ls_num*sizeof(float));

    /* MDAP (with spread) */
    if (spread > 0.1f) {
        const int nSpreadSrcs = 8;
        const int nRings = 1;
        float* U_spread;
        U_spread = malloc1d((nRings*nSpreadSrcs+1)*3*sizeof(float));
        for(ns=0; ns<src_num; ns++){
            azi_rad  = src_dirs[ns*2+0]*SAF_PI/180.0f;
            elev_rad = src_dirs[ns*2+1]*SAF_PI/180.0f;
            getSpreadSrcDirs3D(azi_rad, elev_rad, spread, nSpreadSrcs, nRings, U_spread);
            memset(gains, 0, ls_num*sizeof(float));
            for(nspr=0; nspr<(nRings*nSpreadSrcs+1); nspr++){
                u[0] = U_spread[nspr*3+0];
                u[1] = U_spread[nspr*3+1];
                u[2] = U_spread[nspr*3+2];
                for(i=0; i<nFaces; i++){
                    for(j=0; j<3; j++)
                        ls_invMtx_s[j] = layoutInvMtx[i*9+j];
                    utility_svvdot(ls_invMtx_s, u, 3, &g_tmp[0]);
                    for(j=0; j<3; j++)
                        ls_invMtx_s[j] = layoutInvMtx[i*9+j+3];
                    utility_svvdot(ls_invMtx_s, u, 3, &g_tmp[1]);
                    for(j=0; j<3; j++)
                        ls_invMtx_s[j] = layoutInvMtx[i*9+j+6];
                    utility_svvdot(ls_invMtx_s, u, 3, &g_tmp[2]);
                    min_val = 2.23e13f;
                    g_tmp_rms = 0.0;
                    for(j=0; j<3; j++){
                        min_val = SAF_MIN(min_val, g_tmp[j]);
                        g_tmp_rms +=  powf(g_tmp[j], 2.0f);
                    }
                    g_tmp_rms = sqrtf(g_tmp_rms);
                    if(min_val>-0.001){
                        for(j=0; j<3; j++)
                            gains[ls_groups[i*3+j]] += g_tmp[j]/g_tmp_rms;
                    }
                }
            }
            gains_rms = 0.0;
            for(i=0; i<ls_num; i++)
                gains_rms += powf(gains[i], 2.0f);
            gains_rms = sqrtf(gains_rms);
            for(i=0; i<ls_num; i++)
                (*GainMtx)[ns*ls_num+i] = SAF_MAX(gains[i]/gains_rms, 0.0f);
        }

        free(U_spread);
    }
    /* VBAP (no spread) */
    else{
        for(ns=0; ns<src_num; ns++){
            azi_rad  = src_dirs[ns*2+0]*SAF_PI/180.0f;
            elev_rad = src_dirs[ns*2+1]*SAF_PI/180.0f;
            u[0] = cosf(azi_rad)*cosf(elev_rad);
            u[1] = sinf(azi_rad)*cosf(elev_rad);
            u[2] = sinf(elev_rad);
            memset(gains, 0, ls_num*sizeof(float));
            for(i=0; i<nFaces; i++){
                for(j=0; j<3; j++)
                    ls_invMtx_s[j] = layoutInvMtx[i*9+j];
                utility_svvdot(ls_invMtx_s, u, 3, &g_tmp[0]);
                for(j=0; j<3; j++)
                    ls_invMtx_s[j] = layoutInvMtx[i*9+j+3];
                utility_svvdot(ls_invMtx_s, u, 3, &g_tmp[1]);
                for(j=0; j<3; j++)
                    ls_invMtx_s[j] = layoutInvMtx[i*9+j+6];
                utility_svvdot(ls_invMtx_s, u, 3, &g_tmp[2]);
                min_val = 2.23e13f;
                g_tmp_rms = 0.0;
                for(j=0; j<3; j++){
                    min_val = SAF_MIN(min_val, g_tmp[j]);
                    g_tmp_rms +=  powf(g_tmp[j], 2.0f);
                }
                g_tmp_rms = sqrtf(g_tmp_rms);
                if(min_val>-0.001){
                    for(j=0; j<3; j++)
                        gains[ls_groups[i*3+j]] = g_tmp[j]/g_tmp_rms;
                    break;
                }
            }
            gains_rms = 0.0;
            for(i=0; i<ls_num; i++)
                gains_rms += powf(gains[i], 2.0f);
            gains_rms = sqrtf(gains_rms);
            for(i=0; i<ls_num; i++)
                (*GainMtx)[ns*ls_num+i] = SAF_MAX(gains[i]/gains_rms, 0.0f);
        }
    }

    free(gains);
}

void findLsPairs
(
    float* ls_dirs_deg,
    int L,
    int** out_pairs,
    int* numOutPairs
)
{
    int n;
    float* ls_dirs_deg_tmp;
    int* idx_sorted;

    ls_dirs_deg_tmp = malloc1d(L*sizeof(float));
    idx_sorted = malloc1d(L*sizeof(int));
    for(n=0; n<L; n++)
        ls_dirs_deg_tmp[n] = ls_dirs_deg[n*2];

    /* find the loudspeaker pairs by sorting the angles */
    sortf(ls_dirs_deg_tmp, NULL, idx_sorted, L, 0);
    idx_sorted = realloc(idx_sorted, (L+1)*sizeof(int));
    idx_sorted[L] = idx_sorted[0];
    (*out_pairs) = malloc1d(L*2*sizeof(int));
    for(n=0; n<L; n++){
        (*out_pairs)[n*2] = idx_sorted[n];
        (*out_pairs)[n*2+1] = idx_sorted[n+1];
    }
    (*numOutPairs) = L;

    free(ls_dirs_deg_tmp);
    free(idx_sorted);
}

void invertLsMtx2D
(
    float* U_spkr /* L x 2 */,
    int* ls_pairs /* N_group x 2 */,
    int N_pairs,
    float** layoutInvMtx/* N_group x 4 */
)
{
    int i, j, n;
    float tempGroup[4];
    float tempInv[4];
    void* hSinv;

    /* pre-calculate inversions of the loudspeaker groups and store into matrix */
    (*layoutInvMtx) = malloc1d(N_pairs * 4 * sizeof(float));
    utility_sinv_create(&hSinv, 2);
    for(n=0; n<N_pairs; n++){
        /* get the unit vectors for the current group */
        for(i=0; i<2; i++)
            for(j=0; j<2; j++)
                tempGroup[j*2+i] = U_spkr[ls_pairs[n*2+i]*2 + j]; /* ^T */

        /* get inverse of current group */
        utility_sinv(hSinv, tempGroup, tempInv, 2);

        /* store the vectorised inverse as a row in the output */
        cblas_scopy(4, tempInv, 1, (*layoutInvMtx) + n*4, 1);
    }

    utility_sinv_destroy(&hSinv);
}

void vbap2D
(
    float* src_dirs,
    int src_num,
    int ls_num,
    int* ls_pairs,
    int N_pairs,
    float* layoutInvMtx,
    float** GainMtx
)
{
    int i, j, ns;
    float azi_rad, min_val, g_tmp_rms, gains_rms;
    float u[2], g_tmp[2], ls_invMtx_s[2];
    float* gains;

    (*GainMtx) = malloc1d(src_num*ls_num*sizeof(float));
    gains = malloc1d(ls_num*sizeof(float));
    for(ns=0; ns<src_num; ns++){
        azi_rad  = src_dirs[ns]*SAF_PI/180.0f;
        u[0] = cosf(azi_rad);
        u[1] = sinf(azi_rad);
        memset(gains, 0, ls_num*sizeof(float));
        for(i=0; i<N_pairs; i++){
            for(j=0; j<2; j++)
                ls_invMtx_s[j] = layoutInvMtx[i*4+j];
            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 1, 1, 2, 1.0,
                        ls_invMtx_s, 2,
                        u, 2, 0.0,
                        &g_tmp[0], 1);
            for(j=0; j<2; j++)
                ls_invMtx_s[j] = layoutInvMtx[i*4+j+2];
            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 1, 1, 2, 1.0,
                        ls_invMtx_s, 2,
                        u, 2, 0.0,
                        &g_tmp[1], 1);
            min_val = 2.23e13f;
            g_tmp_rms = 0.0;
            for(j=0; j<2; j++){
                min_val = SAF_MIN(min_val, g_tmp[j]);
                g_tmp_rms +=  powf(g_tmp[j], 2.0f);
            }
            g_tmp_rms = sqrtf(g_tmp_rms);
            if(min_val>-0.001){
                for(j=0; j<2; j++)
                    gains[ls_pairs[i*2+j]] = g_tmp[j]/g_tmp_rms;
            }
        }
        gains_rms = 0.0;
        for(i=0; i<ls_num; i++)
            gains_rms += powf(gains[i], 2.0f);
        gains_rms = sqrtf(gains_rms);
        for(i=0; i<ls_num; i++)
            (*GainMtx)[ns*ls_num+i] = SAF_MAX(gains[i]/gains_rms, 0.0f);
    }

    free(gains);
}
