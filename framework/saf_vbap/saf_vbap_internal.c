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
 *     saf_vbap_internal.c
 * Description:
 *    vbap functions largely derived from the MATLAB library by Archontis Politis,
 *    found here: https://github.com/polarch/Vector-Base-Amplitude-Panning
 * Dependencies:
 *     saf_utilities
 * Author, date created:
 *     Leo McCormack, 02.10.2017
 */

#include "saf_vbap.h"
#include "saf_vbap_internal.h"

/* highly recommended to use double floating point precision when triangulating large spherical meshes.*/
typedef double REAL;

typedef struct sort_REAL {
    REAL val;
    int idx;
}sort_REAL;

int cmp_asc_REAL(const void *a,const void *b) {
    struct sort_REAL *a1 = (struct sort_REAL*)a;
    struct sort_REAL *a2 = (struct sort_REAL*)b;
    if((*a1).val<(*a2).val)return -1;
    else if((*a1).val>(*a2).val)return 1;
    else return 0;
}

int cmp_desc_REAL(const void *a,const void *b) {
    struct sort_REAL *a1 = (struct sort_REAL*)a;
    struct sort_REAL *a2 = (struct sort_REAL*)b;
    if((*a1).val>(*a2).val)return -1;
    else if((*a1).val<(*a2).val)return 1;
    else return 0;
}

typedef struct sort_int {
    int val;
    int idx;
}sort_int;

int cmp_asc_int(const void *a,const void *b) {
    struct sort_int *a1 = (struct sort_int*)a;
    struct sort_int *a2 = (struct sort_int*)b;
    if((*a1).val<(*a2).val)return -1;
    else if((*a1).val>(*a2).val)return 1;
    else return 0;
}

int cmp_desc_int(const void *a,const void *b) {
    struct sort_int *a1 = (struct sort_int*)a;
    struct sort_int *a2 = (struct sort_int*)b;
    if((*a1).val>(*a2).val)return -1;
    else if((*a1).val<(*a2).val)return 1;
    else return 0;
}

static void please_sort_REAL
(
    REAL* in_vec,    /* vector[len] to be sorted */
    REAL* out_vec,   /* if NULL, then in_vec is sorted "in-place" */
    int* new_idices, /* set to NULL if you don't need them */
    int len,         /* number of elements in vectors, must be consistent with the input data */
    int descendFLAG  /* !1:ascending, 1:descending */
)
{
    int i;
    struct sort_REAL *data;
    
    data = malloc(len*sizeof(sort_REAL));
    for(i=0;i<len;i++) {
        data[i].val=in_vec[i];
        data[i].idx=i;
    }
    if(descendFLAG)
        qsort(data,len,sizeof(data[0]),cmp_desc_REAL);
    else
        qsort(data,len,sizeof(data[0]),cmp_asc_REAL);
    for(i=0;i<len;i++){
        if (out_vec!=NULL)
            out_vec[i] = data[i].val;
        else
            in_vec[i] = data[i].val; /* overwrite input vector */
        if(new_idices!=NULL)
            new_idices[i] = data[i].idx;
    }
    free(data);
}

static void please_sort_int
(
    int* in_vec,     /* vector[len] to be sorted */
    int* out_vec,    /* if NULL, then in_vec is sorted "in-place" */
    int* new_idices, /* set to NULL if you don't need them */
    int len,         /* number of elements in vectors, must be consistent with the input data */
    int descendFLAG  /* !1:ascending, 1:descending */
)
{
    int i;
    struct sort_int *data;
    
    data = malloc(len*sizeof(sort_int));
    for(i=0;i<len;i++) {
        data[i].val=in_vec[i];
        data[i].idx=i;
    }
    if(descendFLAG)
        qsort(data,len,sizeof(data[0]),cmp_desc_int);
    else
        qsort(data,len,sizeof(data[0]),cmp_asc_int);
    for(i=0;i<len;i++){
        if (out_vec!=NULL)
            out_vec[i] = data[i].val;
        else
            in_vec[i] = data[i].val; /* overwrite input vector */
        if(new_idices!=NULL)
            new_idices[i] = data[i].idx;
    }
    free(data);
}

/* calculates the determinent of a 4x4 matrix */
REAL convhull_det_4x4(REAL* m) {
    return
    m[3] * m[6] * m[9] * m[12] - m[2] * m[7] * m[9] * m[12] -
    m[3] * m[5] * m[10] * m[12] + m[1] * m[7] * m[10] * m[12] +
    m[2] * m[5] * m[11] * m[12] - m[1] * m[6] * m[11] * m[12] -
    m[3] * m[6] * m[8] * m[13] + m[2] * m[7] * m[8] * m[13] +
    m[3] * m[4] * m[10] * m[13] - m[0] * m[7] * m[10] * m[13] -
    m[2] * m[4] * m[11] * m[13] + m[0] * m[6] * m[11] * m[13] +
    m[3] * m[5] * m[8] * m[14] - m[1] * m[7] * m[8] * m[14] -
    m[3] * m[4] * m[9] * m[14] + m[0] * m[7] * m[9] * m[14] +
    m[1] * m[4] * m[11] * m[14] - m[0] * m[5] * m[11] * m[14] -
    m[2] * m[5] * m[8] * m[15] + m[1] * m[6] * m[8] * m[15] +
    m[2] * m[4] * m[9] * m[15] - m[0] * m[6] * m[9] * m[15] -
    m[1] * m[4] * m[10] * m[15] + m[0] * m[5] * m[10] * m[15];
}

/* Calculates the coefficients of the equation of a PLANE in 3D.
 * Copyright (c) 2014, George Papazafeiropoulos
 * All rights reserved. Distributed under the BSD (2-clause) license
 */
static void convhull_plane_3d
(
    REAL* p,
    REAL* c,
    REAL* d
)
{
    int i, j, k, l;
    int m;
    int* r;
    REAL sign, det, norm_c;
    REAL* pdiff, *pdiff_s;

    m=3;
    pdiff = malloc( (m-1)*m *sizeof(REAL));
    for(i=0; i<m-1; i++)
        for(j=0; j<m; j++)
            pdiff[i*m+j] = p[(i+1)*3+j] - p[i*3+j];
    memset(c, 0, m*sizeof(REAL));
    sign = 1.0;
    r=malloc(m*sizeof(int));
    for(i=0; i<m; i++)
        r[i] = i;
    pdiff_s = malloc((m-1)*(m-1)*sizeof(REAL));
    for(i=0; i<m; i++){
        for(j=0; j<m-1; j++){
            for(k=0, l=0; k<m; k++){
                if(r[k]!=i){
                    pdiff_s[j*(m-1)+l] = pdiff[j*m+k];
                    l++;
                }
            }
        }
        det = pdiff_s[0*(m-1)+0]*pdiff_s[1*(m-1)+1] -
              pdiff_s[1*(m-1)+0]*pdiff_s[0*(m-1)+1];
        c[i] = sign * det;
        sign *= -1.0;
    }
    norm_c = 0.0;
    for(i=0; i<m; i++)
        norm_c += (pow(c[i], 2.0));  
    norm_c = sqrt(norm_c);
    for(i=0; i<m; i++)
        c[i] /= norm_c;
    (*d) = 0.0;
    for(i=0; i<m; i++)
        (*d) += -p[i] * c[i];

    free(pdiff);
    free(r);
    free(pdiff_s);
}

static void convhull_ismember
(
    int* pLeft,
    int* pRight,
    int* pOut,     /* contains 0's and 1's; [nLeftElements]*/
    int nLeftElements,
    int nRightElements
)
{
    int i, j;
    memset(pOut, 0, nLeftElements*sizeof(int));
    for(i=0; i< nLeftElements; i++)
        for(j=0; j< nRightElements; j++)
            if(pLeft[i] == pRight[j] )
                pOut[i] = 1;
}

/* a stripped down C version of the 3D convex hull matlab implementation from here:
 * https://www.mathworks.com/matlabcentral/fileexchange/48509-computational-geometry-toolbox?focused=3851550&tab=example
 * (*out_faces) is returned as NULL, if triangulation fails *
 * Copyright (c) 2014, George Papazafeiropoulos
 * All rights reserved. Distributed under the BSD (2-clause) license
 */
static void convhull_3d
(
    REAL* vertices,
    int nVert,
    int** out_faces,
    int* nOut_faces
)
{
    int i, j, k, l, h;
    int d, nFaces, p;
    int* aVec, *faces;
    REAL dfi, v, max_p, min_p;
    REAL* points, *cf, *cfi, *df, *p_s, *span;
    
    /* if (*out_faces) returns as NULL, then the triangulation failed */
    d = 3; /* number of dimensions */ 
    span = malloc(d*sizeof(REAL));
    for(j=0; j<d; j++){
        max_p = 2.23e-13; min_p = 2.23e+13;
        for(i=0; i<nVert; i++){
            max_p = MAX(max_p, vertices[i*d+j]);
            min_p = MIN(min_p, vertices[i*d+j]);
        }
        span[j] = max_p - min_p;
    }
    points = malloc(nVert*(d+1)*sizeof(REAL));
    for(i=0; i<nVert; i++){
        for(j=0; j<3; j++)
            points[i*4+j] = vertices[i*3 + j];
        points[i*4+3] = 1.0f; /* add a last column of ones. Used only for determinant calculation */
    }
    
    /* The initial convex hull is a simplex with (d+1) facets */
    nFaces = d+1;
    faces = calloc(nFaces*d, sizeof(int));
    aVec = malloc(nFaces*sizeof(int));
    for(i=0; i<nFaces; i++)
        aVec[i] = i;
    
    /* Each column of cf contains the coefficients of a plane */
    cf = malloc(nFaces*d*sizeof(REAL));
    cfi = malloc(d*sizeof(REAL));
    df = malloc(nFaces*sizeof(REAL));
    p_s = malloc(d*d*sizeof(REAL));
    for(i=0; i<nFaces; i++){
        /* Set the indices of the points defining the face  */
        for(j=0, k=0; j<d+1; j++){
            if(aVec[j]!=i){
                faces[i*d+k] = aVec[j];
                k++;
            }
        }
        
        /* Calculate and store the plane coefficients of the face */
        for(j=0; j<d; j++)
            for(k=0; k<d; k++)
                p_s[j*d+k] = points[(faces[i*d+j])*(d+1) + k];
        
        /* Calculate and store the plane coefficients of the face */
        convhull_plane_3d(p_s, cfi, &dfi);
        for(j=0; j<d; j++)
            cf[i*d+j] = cfi[j];
        df[i] = dfi;
    }
    REAL *A;
    int *bVec, *fVec, *asfVec, *face_tmp;
    
    /* Check to make sure that faces are correctly oriented */
    bVec = malloc((d+1)*sizeof(int));
    for(i=0; i<d+1; i++)
        bVec[i] = i;
    
    /* A contains the coordinates of the points forming a simplex */
    A = calloc((d+1)*(d+1), sizeof(REAL));
    face_tmp = malloc(d*sizeof(int));
    fVec = malloc(d*sizeof(int));
    asfVec = malloc(d*sizeof(int));
    for(k=0; k<d+1; k++){
        /* Get the point that is not on the current face (point p) */
        for(i=0; i<d; i++)
            fVec[i] = faces[k*d+i];
        please_sort_int(fVec, NULL, NULL, d, 0); /* sort accending */
        p=k;
        for(i=0; i<d; i++)
            for(j=0; j<d+1; j++)
                A[i*(d+1)+j] = points[(faces[k*d+i])*(d+1) + j];
        for(; i<d+1; i++)
            for(j=0; j<d+1; j++)
                A[i*(d+1)+j] = points[p*(d+1)+j];
        
        /* det(A) determines the orientation of the face */
        v = convhull_det_4x4(A);
        
        /* Orient so that each point on the original simplex can't see the opposite face */
        if(v<0){
            /* Reverse the order of the last two vertices to change the volume */
            for(j=0; j<d; j++)
                face_tmp[j] = faces[k*d+j];
            for(j=0, l=d-2; j<d-1; j++, l++)
                faces[k*d+l] = face_tmp[d-j-1];
            
            /* Modify the plane coefficients of the properly oriented faces */
            for(j=0; j<d; j++)
                cf[k*d+j] = -cf[k*d+j];
            df[k] = -df[k];
            for(i=0; i<d; i++)
                for(j=0; j<d+1; j++)
                    A[i*(d+1)+j] = points[(faces[k*d+i])*(d+1) + j];
            for(; i<d+1; i++)
                for(j=0; j<d+1; j++)
                    A[i*(d+1)+j] = points[p*(d+1)+j];
        }
    }
    
    /* Coordinates of the center of the point set */
    REAL* meanp, *absdist, *reldist, *desReldist;
    meanp = calloc(d, sizeof(REAL));
    for(i=d+1; i<nVert; i++)
        for(j=0; j<d; j++)
            meanp[j] += points[i*(d+1)+j];
    for(j=0; j<d; j++)
        meanp[j] = meanp[j]/(REAL)(nVert-d-1);
    
    /* Absolute distance of points from the center */
    absdist = malloc((nVert-d-1)*d * sizeof(REAL));
    for(i=d+1, k=0; i<nVert; i++, k++)
        for(j=0; j<d; j++)
            absdist[k*d+j] = (points[i*(d+1)+j] -  meanp[j])/span[j];
    
    /* Relative distance of points from the center */
    reldist = calloc((nVert-d-1), sizeof(REAL));
    desReldist = malloc((nVert-d-1) * sizeof(REAL));
    for(i=0; i<(nVert-d-1); i++)
        for(j=0; j<d; j++)
            reldist[i] += pow(absdist[i*d+j], 2.0);
    
    /* Sort from maximum to minimum relative distance */
    int num_pleft, cnt;
    int* ind, *pleft;
    ind = malloc((nVert-d-1) * sizeof(int));
    pleft = malloc((nVert-d-1) * sizeof(int));
    please_sort_REAL(reldist, desReldist, ind, (nVert-d-1), 1);
    
    /* Initialize the vector of points left. The points with the larger relative
       distance from the center are scanned first. */
    num_pleft = (nVert-d-1);
    for(i=0; i<num_pleft; i++)
        pleft[i] = ind[i]+d+1;
    
    /* Loop over all remaining points that are not deleted. Deletion of points
       occurs every #iter2del# iterations of this while loop */
    memset(A, 0, (d+1)*(d+1) * sizeof(REAL));
    
    /* cnt is equal to the points having been selected without deletion of
       nonvisible points (i.e. points inside the current convex hull) */
    cnt=0;
    
    /* The main loop: */
    REAL detA;
    REAL* points_cf, *points_s;
    int* visible_ind, *visible, *nonvisible_faces, *f0, *face_s, *u, *gVec, *horizon, *hVec, *pp, *hVec_mem_face;
    int num_visible_ind, num_nonvisible_faces, n_newfaces, count, vis;
    int f0_sum, u_len, start, num_p, index, horizon_size1;
    int FUCKED;
    FUCKED = 0;
    u = horizon = NULL;
    nFaces = d+1;
    visible_ind = malloc(nFaces*sizeof(int));
    points_cf = malloc(nFaces*sizeof(REAL));
    points_s = malloc(d*sizeof(REAL));
    face_s = malloc(d*sizeof(int));
    gVec = malloc(d*sizeof(int));
    int HARD_STOP =0;
    while( (num_pleft>0) ){//}&& (HARD_STOP != 10 )){
        HARD_STOP++;
        /* i is the first point of the points left */
        i = pleft[0];
        
        /* Delete the point selected */
        for(j=0; j<num_pleft-1; j++)
            pleft[j] = pleft[j+1];
        num_pleft--;
        if(num_pleft == 0)
            free(pleft);
        else
            pleft = realloc(pleft, num_pleft*sizeof(int));
        
        /* Update point selection counter */
        cnt++;
        
        /* find visible faces */
        for(j=0; j<d; j++)
            points_s[j] = points[i*(d+1)+j];
        points_cf = realloc(points_cf, nFaces*sizeof(REAL));
        visible_ind = realloc(visible_ind, nFaces*sizeof(int));
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 1, nFaces, d, 1.0,
                     points_s, d,
                     cf, d, 0.0,
                     points_cf, nFaces);
        num_visible_ind = 0;
        for(j=0; j<nFaces; j++){
            if(points_cf[j] + df[j] > 0.0){
                num_visible_ind++; /* will sum to 0 if none are visible */
                visible_ind[j] = 1;
            }
            else
                visible_ind[j] = 0;
        }
        num_nonvisible_faces = nFaces - num_visible_ind;
        
        /* proceed if there are any visible faces */
        if(num_visible_ind!=0){
            /* Find visible face indices */
            visible = malloc(num_visible_ind*sizeof(int));
            for(j=0, k=0; j<nFaces; j++){
                if(visible_ind[j]==1){
                    visible[k]=j;
                    k++;
                }
            }
            
            /* Find nonvisible faces */
            nonvisible_faces = malloc(num_nonvisible_faces*d*sizeof(int));
            f0 = malloc(num_nonvisible_faces*d*sizeof(int));
            for(j=0, k=0; j<nFaces; j++){
                if(visible_ind[j]==0){
                    for(l=0; l<d; l++)
                        nonvisible_faces[k*d+l]= faces[j*d+l];
                    k++;
                }
            }
            
            /* Create horizon (count is the number of the edges of the horizon) */
            count=0;
            for(j=0; j<num_visible_ind; j++){
                /* visible face */
                vis = visible[j];
                for(k=0; k<d; k++)
                    face_s[k] = faces[vis*d+k];
                please_sort_int(face_s, NULL, NULL, d, 0);
                convhull_ismember(nonvisible_faces, face_s, f0, num_nonvisible_faces*d, d);
                u_len = 0;
                
                /* u are the nonvisible faces connected to the face v, if any */
                for(k=0; k<num_nonvisible_faces; k++){
                    f0_sum = 0;
                    for(l=0; l<d; l++)
                        f0_sum += f0[k*d + l];
                    if(f0_sum == d-1){
                        u_len++;
                        if(u_len==1)
                            u=malloc(u_len*sizeof(int));
                        else
                            u=realloc(u, u_len*sizeof(int));
                        u[u_len-1] = k;
                    }
                }
                for(k=0; k<u_len; k++){
                    /* The boundary between the visible face v and the k(th) nonvisible face connected to the face v forms part of the horizon */
                    count++;
                    if(count==1)
                        horizon = malloc(count*(d-1)*sizeof(int));
                    else
                        horizon = realloc(horizon, count*(d-1)*sizeof(int));
                    for(l=0; l<d; l++)
                        gVec[l] = nonvisible_faces[u[k]*d+l];
                    for(l=0, h=0; l<d; l++){
                        if(f0[u[k]*d+l]){
                            horizon[(count-1)*(d-1)+h] = gVec[l];
                            h++;
                        }
                    }
                }
                if(u_len!=0)
                    free(u);
            }
            horizon_size1 = count;
            for(j=0, l=0; j<nFaces; j++){
                if(!visible_ind[j]){
                    /* Delete visible faces */
                    for(k=0; k<d; k++)
                        faces[l*d+k] = faces[j*d+k];
                    
                    /* Delete the corresponding plane coefficients of the faces */
                    for(k=0; k<d; k++)
                        cf[l*d+k] = cf[j*d+k];
                    df[l] = df[j];
                    l++;
                }
            }
            
            /* Update the number of faces */
            nFaces = nFaces-num_visible_ind;
            faces = realloc(faces, nFaces*d*sizeof(int));
            cf = realloc(cf, nFaces*d*sizeof(REAL));
            df = realloc(df, nFaces*sizeof(REAL));
            
            /* start is the first row of the new faces */
            start=nFaces;
            
            /* Add faces connecting horizon to the new point */
            n_newfaces = horizon_size1;
            for(j=0; j<n_newfaces; j++){
                nFaces++;
                faces = realloc(faces, nFaces*d*sizeof(int));
                cf = realloc(cf, nFaces*d*sizeof(REAL));
                df = realloc(df, nFaces*sizeof(REAL));
                for(k=0; k<d-1; k++)
                    faces[(nFaces-1)*d+k] = horizon[j*(d-1)+k];
                faces[(nFaces-1)*d+(d-1)] = i;
                
                /* Calculate and store appropriately the plane coefficients of the faces */
                for(k=0; k<d; k++)
                    for(l=0; l<d; l++)
                        p_s[k*d+l] = points[(faces[(nFaces-1)*d+k])*(d+1) + l];
                convhull_plane_3d(p_s, cfi, &dfi);
                for(k=0; k<d; k++)
                    cf[(nFaces-1)*d+k] = cfi[k];
                df[(nFaces-1)] = dfi;
                if(nFaces > MAX_NUM_FACES){
                    FUCKED = 1;
                    nFaces = 0;
                    break; 
                }
            }
            
            /* Orient each new face properly */
            hVec = malloc( nFaces*sizeof(int));
            hVec_mem_face = malloc( nFaces*sizeof(int));
            for(j=0; j<nFaces; j++)
                hVec[j] = j;
            for(k=start; k<nFaces; k++){
                for(j=0; j<d; j++)
                    face_s[j] = faces[k*d+j];
                please_sort_int(face_s, NULL, NULL, d, 0);
                convhull_ismember(hVec, face_s, hVec_mem_face, nFaces, d);
                num_p = 0;
                for(j=0; j<nFaces; j++)
                    if(!hVec_mem_face[j])
                        num_p++;
                pp=malloc(num_p*sizeof(int));
                for(j=0, l=0; j<nFaces; j++){
                    if(!hVec_mem_face[j]){
                        pp[l] = hVec[j];
                        l++;
                    }
                }
                index = 0;
                detA = 0.0;
                
                /* While new point is coplanar, choose another point */
                while(detA==0.0){
                    for(j=0;j<d; j++)
                        for(l=0; l<d+1; l++)
                            A[j*(d+1)+l] = points[(faces[k*d+j])*(d+1) + l];
                    for(; j<d+1; j++)
                        for(l=0; l<d+1; l++)
                            A[j*(d+1)+l] = points[pp[index]*(d+1)+l];
                    index++;
                    detA = convhull_det_4x4(A);
                }
                
                /* Orient faces so that each point on the original simplex can't see the opposite face */
                if (detA<0.0){
                    /* If orientation is improper, reverse the order to change the volume sign */
                    for(j=0; j<d; j++)
                        face_tmp[j] = faces[k*d+j];
                    for(j=0, l=d-2; j<d-1; j++, l++)
                        faces[k*d+l] = face_tmp[d-j-1];
                    
                    /* Modify the plane coefficients of the properly oriented faces */
                    for(j=0; j<d; j++)
                        cf[k*d+j] = -cf[k*d+j];
                    df[k] = -df[k];
                    for(l=0; l<d; l++)
                        for(j=0; j<d+1; j++)
                            A[l*(d+1)+j] = points[(faces[k*d+l])*(d+1) + j];
                    for(; l<d+1; l++)
                        for(j=0; j<d+1; j++)
                            A[l*(d+1)+j] = points[pp[index]*(d+1)+j];
                }
                free(pp);
            }
            free(horizon);
            free(f0);
            free(nonvisible_faces);
            free(visible);
            free(hVec);
            free(hVec_mem_face);
        }
        if(FUCKED){
            break;
        }
    }
    
    /* output */
    if(FUCKED){
        (*out_faces) = NULL;
        (*nOut_faces) = nFaces;
    }
    else{
        (*out_faces) = (int*)malloc(nFaces*d*sizeof(int));
        memcpy((*out_faces),faces, nFaces*d*sizeof(int));
        (*nOut_faces) = nFaces;
    }
    
    /* clean-up */
    free(visible_ind);
    free(points_cf);
    free(points_s);
    free(face_s);
    free(gVec); 
    free(meanp);
    free(absdist);
    free(reldist);
    free(desReldist);
    free(ind);
    free(span);
    free(points);
    free(faces);
    free(aVec);
    free(cf);
    free(cfi);
    free(df);
    free(p_s);
    free(face_tmp);
    free(fVec);
    free(asfVec);
    free(bVec);
    free(A);
}

static void ccross(float a[3], float b[3], float c[3]){
    c[0] = a[1]*b[2]-a[2]*b[1];
    c[1] = a[2]*b[0]-a[0]*b[2];
    c[2] = a[0]*b[1]-a[1]*b[0];
}

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
    int i, j, k, numValidFaces, minIntVal, minIdx, nFaces;
    int tmp3[3], circface[3];
    int* validFacesID, *valid_faces, *valid_faces2, *faces;
    float dotcc, aperture_lim;
    float vecs[3][3], cvec[3], centroid[3], a[3], b[3], abc[3];
    REAL r, rcoselev;
    REAL* vertices;
    
    /* Find the convex hull of the points on the sphere - in this special case the
       result equals the Delaunay triangulation of the points */
    vertices = malloc(L*3*sizeof(REAL));  
    (*numOutVertices)  = L;
    (*out_vertices) = (float*)malloc(L*3*sizeof(float));
    for ( i = 0; i < L; i++) {
        if(L>1000)
            r = 1.0 + (REAL)rand()/(REAL)RAND_MAX*0.00001; /* add some noise to the data (it helps) */
        else
            r = 1.0 + (REAL)rand()/(REAL)RAND_MAX*0.001; /* add some noise to the data (it helps) */
        (*out_vertices)[i*3+2] = (float)(r* (REAL)sin((double)ls_dirs_deg[i*2+1]*M_PI/180.0));
        rcoselev = r*(REAL)cos((double)ls_dirs_deg[i*2+1]*M_PI/180.0);
        (*out_vertices)[i*3+0] = (float)(rcoselev * (REAL)cos((double)ls_dirs_deg[i*2+0]*M_PI/180.0));
        (*out_vertices)[i*3+1] = (float)(rcoselev * (REAL)sin((double)ls_dirs_deg[i*2+0]*M_PI/180.0));
        if(L>1000) {/* starts becoming unstable past this point */
            for(j=0; j<3; j++)
                vertices[i*3+j] = (REAL)(*out_vertices)[i*3+j] + (REAL)rand()/(REAL)RAND_MAX*0.000001; /* add some noise to the data (it helps) */
        }
        else if(L>100){
            for(j=0; j<3; j++)
                vertices[i*3+j] = (REAL)(*out_vertices)[i*3+j] + (REAL)rand()/(REAL)RAND_MAX*0.001; /* add some noise to the data (it helps) */
        }
        else{
            for(j=0; j<3; j++)
                vertices[i*3+j] = (REAL)(*out_vertices)[i*3+j] + (REAL)rand()/(REAL)RAND_MAX*0.01; /* add some noise to the data (it helps) */
        }
    }
    faces = NULL;
    convhull_3d(vertices, L, &faces, &nFaces);
    if(faces==NULL){
        /* out_faces returned as NULL, if triangulation failed */
        (*out_vertices) = NULL;
        (*numOutVertices) = 0;
        (*out_faces) = NULL;
        (*numOutFaces) = 0;
        return;
    }
    
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
#if 0
    /* sort indices in accending order for the third dimension */
    for (j=0; j<nFaces - 1; j++)  {
        for (i=0; i<nFaces - 1; i++) {
            if ( (faces[(i+1)*3+2] < faces[i*3+2]) && (faces[(i+1)*3+0] == faces[i*3+0]) && (faces[(i+1)*3+1] == faces[i*3+1])) {
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
    validFacesID = malloc(nFaces*sizeof(int));
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
        if(acosf(MAX(MIN(dotcc,0.99999999f),-0.99999999f)/* avoids complex numbers */)<(M_PI/2.0f)){
            validFacesID[i] = 1;
            numValidFaces++;
        }
        else{
            validFacesID[i] = 0;
        }
    }
    valid_faces = malloc(numValidFaces*3*sizeof(int));
    for(i=0, j=0; i<nFaces; i++){
        if(validFacesID[i]==1){
            valid_faces[j*3+0] = faces[i*3+0];
            valid_faces[j*3+1] = faces[i*3+1];
            valid_faces[j*3+2] = faces[i*3+2];
            j++;
        }
    }
    free(validFacesID);
    
    /* Omit Triangules that have an aperture larger than APERTURE_LIMIT_DEG */
    valid_faces2 = NULL;
    if(omitLargeTriangles) {
        aperture_lim = APERTURE_LIMIT_DEG * M_PI/180.0f;
        nFaces = numValidFaces;
        numValidFaces = 0;
        validFacesID = malloc(nFaces*sizeof(int));
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
        valid_faces2 = malloc(numValidFaces*3*sizeof(int));
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
    (*out_faces) = (int*)malloc(numValidFaces*3*sizeof(int));
    if(omitLargeTriangles)
        memcpy((*out_faces), valid_faces2, numValidFaces*3*sizeof(int));
    else
        memcpy((*out_faces), valid_faces, numValidFaces*3*sizeof(int));
    
    /* clean-up */
    if(faces!=NULL)
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
    float* tempGroup;
    float tempInv[9], eye3[9];
    
    for(i=0; i<3; i++)
        for(j=0; j<3; j++)
            eye3[i*3+j] = i==j ? 1.0f : 0.0f;
    tempGroup = malloc(9* sizeof(float));
    
    /* pre-calculate inversions of the loudspeaker groups and store into matrix */
    (*layoutInvMtx) = malloc(N_group * 9 * sizeof(float));
    for(n=0; n<N_group; n++){
        /* get the unit vectors for the current group */
        for(i=0; i<3; i++)
            for(j=0; j<3; j++)
                tempGroup[i*3+j] = U_spkr[ls_groups[n*3+i]*3 + j];
        
        /* get inverse of current group */
        utility_sinv(tempGroup,3);
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0,
                    eye3, 3,
                    tempGroup, 3, 0.0,
                    tempInv, 3);
        
        /* store the vectorized inverse as a row the output */
        for(i=0; i<3; i++)
            for(j=0; j<3; j++)
                (*layoutInvMtx)[n*9+(i*3+j)] = tempInv[j*3+i];
    }
    free(tempGroup);
}

void vbap3D
(
    float* src_dirs,
    int src_num,
    int ls_num,
    int* ls_groups,
    int nFaces,
    float* layoutInvMtx,
    float** GainMtx
)
{
    int i, j, ns;
    float azi_rad, elev_rad, min_val, g_tmp_rms, gains_rms;
    float u[3], g_tmp[3], ls_invMtx_s[3];
    float* gains;
    
    (*GainMtx) = malloc(src_num*ls_num*sizeof(float));
    gains = malloc(ls_num*sizeof(float)); 
    for(ns=0; ns<src_num; ns++){
        azi_rad  = src_dirs[ns*2+0]*M_PI/180.0f;
        elev_rad = src_dirs[ns*2+1]*M_PI/180.0f;
        u[0] = cosf(azi_rad)*cosf(elev_rad);
        u[1] = sinf(azi_rad)*cosf(elev_rad);
        u[2] = sinf(elev_rad);
        memset(gains, 0, ls_num*sizeof(float));
        for(i=0; i<nFaces; i++){
            for(j=0; j<3; j++)
                ls_invMtx_s[j] = layoutInvMtx[i*9+j];
            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 1, 1, 3, 1.0,
                        ls_invMtx_s, 3,
                        u, 3, 0.0,
                        &g_tmp[0], 1);
            for(j=0; j<3; j++)
                ls_invMtx_s[j] = layoutInvMtx[i*9+j+3];
            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 1, 1, 3, 1.0,
                        ls_invMtx_s, 3,
                        u, 3, 0.0,
                        &g_tmp[1], 1);
            for(j=0; j<3; j++)
                ls_invMtx_s[j] = layoutInvMtx[i*9+j+6];
            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 1, 1, 3, 1.0,
                        ls_invMtx_s, 3,
                        u, 3, 0.0,
                        &g_tmp[2], 1);
            min_val = 2.23e13f;
            g_tmp_rms = 0.0;
            for(j=0; j<3; j++){
                min_val = MIN(min_val, g_tmp[j]);
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
            (*GainMtx)[ns*ls_num+i] = MAX(gains[i]/gains_rms, 0.0f); 
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
    REAL* ls_dirs_deg_tmp;
    int* idx_sorted;
    
    ls_dirs_deg_tmp = malloc(L*sizeof(REAL));
    idx_sorted = malloc(L*sizeof(int));
    for(n=0; n<L; n++)
        ls_dirs_deg_tmp[n] = (REAL)ls_dirs_deg[n*2];
    
    /* find the loudspeaker pairs by sorting the angles */
    please_sort_REAL(ls_dirs_deg_tmp, NULL, idx_sorted, L, 0);
    idx_sorted = realloc(idx_sorted, (L+1)*sizeof(int));
    idx_sorted[L] = idx_sorted[0];
    (*out_pairs) = malloc(L*2*sizeof(int));
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
    float* tempGroup;
    float tempInv[4], eye2[4];
    
    for(i=0; i<2; i++)
        for(j=0; j<2; j++)
            eye2[i*2+j] = i==j ? 1.0f : 0.0f;
    tempGroup = malloc(4* sizeof(float));
    
    /* pre-calculate inversions of the loudspeaker groups and store into matrix */
    (*layoutInvMtx) = malloc(N_pairs * 4 * sizeof(float));
    for(n=0; n<N_pairs; n++){
        /* get the unit vectors for the current group */
        for(i=0; i<2; i++)
            for(j=0; j<2; j++)
                tempGroup[i*2+j] = U_spkr[ls_pairs[n*2+i]*2 + j];
        
        /* get inverse of current group */
        utility_sinv(tempGroup,2);
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, 2, 2, 1.0,
                    eye2, 2,
                    tempGroup, 2, 0.0,
                    tempInv, 2);
        
        /* store the vectorized inverse as a row the output */
        for(i=0; i<2; i++)
            for(j=0; j<2; j++)
                (*layoutInvMtx)[n*4+(i*2+j)] = tempInv[j*2+i];
    }
    free(tempGroup);
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
    
    (*GainMtx) = malloc(src_num*ls_num*sizeof(float));
    gains = malloc(ls_num*sizeof(float));
    for(ns=0; ns<src_num; ns++){
        azi_rad  = src_dirs[ns]*M_PI/180.0f;
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
                min_val = MIN(min_val, g_tmp[j]);
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
            (*GainMtx)[ns*ls_num+i] = MAX(gains[i]/gains_rms, 0.0f);
    }
    
    free(gains);
}








