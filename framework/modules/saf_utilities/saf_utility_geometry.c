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
 * @file saf_utility_geometry.c
 * @ingroup Utilities
 * @brief A collection of computational geometry related functions
 *
 * @author Leo McCormack
 * @date 03.07.2020
 */

#include "saf_utilities.h"
#include <stdio.h>
#include <stdlib.h>


float L2_norm(float v[3])
{
    return sqrtf(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

void crossProduct(float a[3], float b[3], float c[3]){
    c[0] = a[1]*b[2]-a[2]*b[1];
    c[1] = a[2]*b[0]-a[0]*b[2];
    c[2] = a[0]*b[1]-a[1]*b[0];
}

void convhull3d
(
    const float* vertices,
    const int nDirs,
    int** faces,
    int* nFaces
)
{
    int i;
    ch_vertex* ch_vertices;

    /* convert vertices to use "ch_vertex" format used by convhull_3d_build() */
    ch_vertices = malloc1d(nDirs*sizeof(ch_vertex));
    for(i = 0; i < nDirs; i++) {
        ch_vertices[i].z = (CH_FLOAT)vertices[i*3+2];
        ch_vertices[i].x = (CH_FLOAT)vertices[i*3];
        ch_vertices[i].y = (CH_FLOAT)vertices[i*3+1];
    }

    /* build convex hull */
    assert(*faces == NULL); /* nFaces not known yet, shouldn't be pre-allocated... */
    convhull_3d_build(ch_vertices, nDirs, faces, nFaces);

    /* clean-up */
    free(ch_vertices);
}

void sphDelaunay
(
    const float* dirs_deg,
    const int nDirs,
    int** faces,
    int* nFaces,
    float* vertices
)
{
    int i;
    float* vertices_tmp;
    float rcoselev;

    /* sph to unit Cart: */
    vertices_tmp = malloc1d(nDirs*3*sizeof(float));
    for(i = 0; i < nDirs; i++) {
        vertices_tmp[i*3+2] = sinf(dirs_deg[i*2+1]* SAF_PI/180.0f);
        rcoselev = cosf(dirs_deg[i*2+1]*SAF_PI/180.0f);
        vertices_tmp[i*3] = rcoselev * cosf(dirs_deg[i*2+0]*SAF_PI/180.0f);
        vertices_tmp[i*3+1] = rcoselev * sinf(dirs_deg[i*2+0]*SAF_PI/180.0f);
    }

    /* Delaunay triangulation of a spherical */
    convhull3d(vertices_tmp, nDirs, faces, nFaces);

    /* optionally, also output the vertices */
    if(vertices!=NULL)
        memcpy(vertices, vertices_tmp, nDirs*3*sizeof(float));

    /* clean up */
    free(vertices_tmp);
}

void sphVoronoi
(
    int* faces,
    int nFaces,
    float* vertices,
    int nDirs,
    voronoi_data* voronoi
)
{
    int n, m;
    int* duplicates;
    float r_12[3], r_13[3], r_normal[3];
    float norm;

    /* prep */
    voronoi->nVert = nFaces;
    voronoi->nFaces = nDirs;
    voronoi->vert = (float**)malloc2d(voronoi->nVert, 3, sizeof(float));
    voronoi->faces = (int**)malloc1d(voronoi->nFaces*sizeof(int*));
    voronoi->nPointsPerFace = malloc1d(voronoi->nFaces*sizeof(int));

    /* Calculate the voronoi vertices for each triangle, which for the unit
     * sphere are given by the unit normal vector of the triangle */
    for(n = 0; n<nFaces; n++){
        r_12[0] = vertices[faces[n*3+1]*3]   - vertices[faces[n*3]*3];
        r_12[1] = vertices[faces[n*3+1]*3+1] - vertices[faces[n*3]*3+1];
        r_12[2] = vertices[faces[n*3+1]*3+2] - vertices[faces[n*3]*3+2];
        r_13[0] = vertices[faces[n*3+2]*3]   - vertices[faces[n*3]*3];
        r_13[1] = vertices[faces[n*3+2]*3+1] - vertices[faces[n*3]*3+1];
        r_13[2] = vertices[faces[n*3+2]*3+2] - vertices[faces[n*3]*3+2];
        crossProduct(r_12, r_13, r_normal);

        norm = 1.0f/L2_norm(r_normal);
        utility_svsmul(r_normal, &norm, 3, r_normal);
        memcpy(voronoi->vert[n], r_normal, 3*sizeof(float));
    }

    /* Find duplicate vertices if any, due to two triangles sharing the same
     * circumscribed circle */
    duplicates = calloc1d(voronoi->nVert, sizeof(int));
    for(n = 0; n<voronoi->nVert; n++){
        if (duplicates[n] == 0 ){
            for (m = 0; m<voronoi->nVert; m++){
                if (n != m){
                    if (fabsf(voronoi->vert[n][0] - voronoi->vert[m][0]) < 1.0e-5f &&
                        fabsf(voronoi->vert[n][1] - voronoi->vert[m][1]) < 1.0e-5f &&
                        fabsf(voronoi->vert[n][2] - voronoi->vert[m][2]) < 1.0e-5f ){
                        duplicates[m] = n;
                    }
                }
            }
        }
    }

    int NOT_SORTED;
    int i, j, k, l, nFaceIdx, currentfaceIdx, currentvertIdx, currentvert, nSorted;
    int currentface[3];
    int* faceIdx, *sorted, *tempfacelist;
    faceIdx = NULL;
    sorted = NULL;
    tempfacelist = NULL;

    /* Calculate the voronoi polygons
     *
     * find the an ordered sequence of the triangles around each vertex and get
     * the proper sequence of the voronoi vertices that constitute a polygon */
    for (n = 0; n<voronoi->nFaces; n++){
        nFaceIdx=0;
        for(m=0; m<voronoi->nVert; m++)
            if(faces[m*3+0]==n || faces[m*3+1]==n || faces[m*3+2]==n)
                nFaceIdx++;
        faceIdx = realloc1d(faceIdx, nFaceIdx*sizeof(int));

        /* list of triangles that contain this specific vertex */
        i=0;
        for(m=0; m<voronoi->nVert; m++){
            if(faces[m*3+0]==n || faces[m*3+1]==n || faces[m*3+2]==n){
                faceIdx[i] = m;
                i++;
            }
        }

        /* Each triangle from the list contain the common vertex and two other
         * vertices - each triangle has two common vertices with each other. One
         * (brute) way of sorting the sequence is to pick one triangle, find the
         * neighbour triangle by finding their common
         * vertex, move to that triangle and iterate till all the number of
         * triangles have been checked. */
        k = 0;
        currentfaceIdx = faceIdx[k]; /* pick-up one of the triangles in the list */
        memcpy(currentface, &faces[currentfaceIdx*3], 3*sizeof(int)); /* the triangle's vertex indices */

        /* pick-up one of the vertices that is not the central one */
        currentvertIdx = -1;
        for(j=0; j<3; j++){
            if(currentface[j] != n){
                currentvertIdx = j;
                break;
            }
        }
        assert(currentvertIdx!=-1);
        currentvert = currentface[currentvertIdx];

        /* Prep */
        nSorted = 1;
        sorted = realloc1d(sorted, nSorted*sizeof(int));
        sorted[0] = faceIdx[k]; /* this is the list that keeps the the ordered triangles */
        NOT_SORTED = 1;
        while(NOT_SORTED){
            /* exclude the current triangle from the temporary list */
            tempfacelist = realloc1d(tempfacelist, (nFaceIdx-1)*sizeof(int));
            l=0;
            for(i=0; i<nFaceIdx; i++){
                if(faceIdx[i]!=currentfaceIdx){
                    tempfacelist[l] = faceIdx[i];
                    l++;
                }
            }
            assert(l==nFaceIdx-1);

            for (l = 0; l<nFaceIdx-1; l++){
                currentfaceIdx = tempfacelist[l];
                memcpy(currentface, &faces[currentfaceIdx*3], 3*sizeof(int));

                /* if the currentvert exists in the current triangles vertices... */
                if (currentface[0] == currentvert || currentface[1] == currentvert || currentface[2] == currentvert){
                    /* then it's the neighbour triangle, so store its index */
                    nSorted++;
                    sorted = realloc1d(sorted, nSorted*sizeof(int));
                    sorted[nSorted-1] = currentfaceIdx;

                    /* if the sorted list has the length of faceIdx, then we're done */
                    if (nSorted == nFaceIdx){
                        NOT_SORTED = 0;
                        break;
                    }

                    /* find the next vertex from current triangle that excludes the central one and the previous one */
                    for(j=0; j<3; j++)
                        if(currentface[j] != n && currentface[j] != currentvert)
                            currentvertIdx = j;
                    currentvert = currentface[currentvertIdx];
                    break;
                }
            }
        }

        /* remove the duplicate vertices from the list */
        for (i=0; i<nSorted; i++)
            if (duplicates[sorted[i]] != 0)
                sorted[i] = duplicates[sorted[i]];


        /* Identify unique IDs, and sort them in accending order */
        int nUnique;
        int* uniqueIDs;
        uniqueIDs = NULL;
        unique_i(sorted, nSorted, NULL, &uniqueIDs, &nUnique);
        sorti(uniqueIDs, uniqueIDs, NULL, nUnique, 0);

        /* Save */
        voronoi->faces[n] = malloc1d(nUnique*sizeof(int));
        for(i=0; i<nUnique; i++)
            voronoi->faces[n][i] = sorted[uniqueIDs[i]];
        voronoi->nPointsPerFace[n] = nUnique;

        /* clean-up */
        free(uniqueIDs);
    }

    /* clean-up */
    free(faceIdx);
    free(sorted);
    free(tempfacelist);
}

void sphVoronoiAreas
(
    voronoi_data* voronoi,
    float* areas
)
{
    int i, m, n, N_poly, tmp_i;
    int* face;
    float tmp, r_21_norm, r_23_norm;
    float *theta;
    float r_01[3], r_02[3], r_2x1[3], r_21[3], r_03[3], r_2x3[3], r_23[3];

    face = NULL;
    theta = NULL;
    for(m=0; m<voronoi->nFaces; m++){
        N_poly = voronoi->nPointsPerFace[m]; /* number of vertices in the polygon */
        face = realloc1d(face, N_poly*sizeof(int));
        theta = realloc1d(theta, N_poly*sizeof(float));
        memcpy(face, voronoi->faces[m], N_poly*sizeof(int)); /* current face */
 
        for(n=0; n<N_poly; n++){
            /* find vector between vertex origin and vertices 1 & 2 */
            memcpy(r_01, voronoi->vert[face[0]], 3*sizeof(float));
            memcpy(r_02, voronoi->vert[face[1]], 3*sizeof(float));

            /* find normal vector to the great circle of 1 & 2 */
            crossProduct(r_02, r_01, r_2x1);

            /* find tangent vector to great circle at 2 */
            crossProduct(r_2x1, r_02, r_21);

            /* find vector between vertex origin and vertex 3 */
            memcpy(r_03, voronoi->vert[face[2]], 3*sizeof(float));

            /* find normal vector to the great circle of 2 & 3 */
            crossProduct(r_02, r_03, r_2x3);

            /* find tangent vector to great circle at 2 */
            crossProduct(r_2x3, r_02, r_23);

            /* normalise tangent vectors */
            r_21_norm = 1.0f/L2_norm(r_21);
            utility_svsmul(r_21, &r_21_norm, 3, r_21);
            r_23_norm = 1.0f/L2_norm(r_23);
            utility_svsmul(r_23, &r_23_norm, 3, r_23);

            /* get angle between the normals */
            utility_svvdot(r_21, r_23, 3, &tmp);
            theta[n] = acosf(tmp);

            /* shift the vertex list by one position and repeat */
            tmp_i = face[0];
            for(i=1; i<N_poly; i++)
                face[i-1] = face[i];
            face[N_poly-1] = tmp_i;
        }

        /* compute area of spherical polygon by spherical excess */
        tmp = 0.0f;
        for(i=0; i<N_poly; i++)
            tmp += theta[i];
        areas[m] = tmp - ((float)N_poly-2.0f)*SAF_PI; 
    }
}

void getVoronoiWeights
(
    float* dirs_deg,
    int nDirs,
    int diagFLAG,
    float* weights
)
{
    int i, nFaces;
    int* faces;
    float* vertices, *areas;
    voronoi_data voronoi;

    /* Perform delaunay triangulation */
    faces = NULL;
    vertices = malloc1d(nDirs*3*sizeof(float));
    sphDelaunay(dirs_deg, nDirs, &faces, &nFaces, vertices);

    /* Get voronoi diagrams */
    sphVoronoi(faces, nFaces, vertices, nDirs, &voronoi);

    /* Compute areas of spherical voronoi polygons */
    areas = malloc1d(voronoi.nFaces * sizeof(float));
    sphVoronoiAreas(&voronoi, areas);

    /* Store weights */
    if(diagFLAG){
        memset(weights, 0, nDirs*nDirs*sizeof(float));
        for(i=0; i<nDirs; i++)
            weights[i*nDirs+i] = areas[i];
    }
    else
        memcpy(weights, areas, nDirs*sizeof(float));

    /* clean-up */
    free(faces);
    free(vertices);
    free(areas);
    for(i=0; i<voronoi.nFaces; i++)
        free(voronoi.faces[i]);
    free(voronoi.vert);
    free(voronoi.nPointsPerFace);
}
