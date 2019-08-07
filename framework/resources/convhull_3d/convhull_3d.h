/*
 Copyright (c) 2017-2018 Leo McCormack
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
*/
/*
 * Filename:
 *     convhull_3d.h
 * Description:
 *     An implementation of the 3-D quickhull algorithm.
 *     The code is largely derived from the "computational-geometry-toolbox"
 *     by George Papazafeiropoulos (c) 2014, originally distributed under
 *     the BSD (2-clause) license.
 *     Taken from: https://github.com/leomccormack/convhull_3d
 *     Reference: "The Quickhull Algorithm for Convex Hull, C. Bradford
 *                 Barber, David P. Dobkin and Hannu Huhdanpaa, Geometry
 *                 Center Technical Report GCG53, July 30, 1993"
 * Dependencies:
 *     cblas (optional for speed ups, especially for very large meshes)
 * Author, date created:
 *     Leo McCormack, 02.10.2017
 */

#ifndef CONVHULL_3D_INCLUDED
#define CONVHULL_3D_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

#ifdef CONVHULL_3D_USE_FLOAT_PRECISION
typedef float CH_FLOAT;
#else
typedef double CH_FLOAT;
#endif
typedef struct _ch_vertex {
    union {
        CH_FLOAT v[3];
        struct{
             CH_FLOAT x, y, z;
        };
    };
} ch_vertex;
typedef ch_vertex ch_vec3;
   
/* builds the convexhull, returning the face indices corresponding to "in_vertices" */
void convhull_3d_build(/* input arguments */
                       ch_vertex* const in_vertices,            /* vector of input vertices; nVert x 1 */
                       const int nVert,                         /* number of vertices */
                       /* output arguments */
                       int** out_faces,                         /* & of empty int*, output face indices; flat: nOut_faces x 3 */
                       int* nOut_faces);                        /* & of int, number of output face indices */
    
/* exports the vertices, face indices, and face normals, as an 'obj' file, ready for GPU */
void convhull_3d_export_obj(/* input arguments */
                            ch_vertex* const vertices,          /* vector of input vertices; nVert x 1 */
                            const int nVert,                    /* number of vertices */
                            int* const faces,                   /* face indices; flat: nFaces x 3 */
                            const int nFaces,                   /* number of faces in hull */
                            const int keepOnlyUsedVerticesFLAG, /* 0: exports in_vertices, 1: exports only used vertices  */
                            char* const obj_filename);          /* obj filename, WITHOUT extension */
    
/* exports the vertices, face indices, and face normals, as an 'm' file, for MatLab verification */
void convhull_3d_export_m(/* input arguments */
                          ch_vertex* const vertices,            /* vector of input vertices; nVert x 1 */
                          const int nVert,                      /* number of vertices */
                          int* const faces,                     /* face indices; flat: nFaces x 3 */
                          const int nFaces,                     /* number of faces in hull */
                          char* const m_filename);              /* m filename, WITHOUT extension */
    
/* reads an 'obj' file and extracts only the vertices */
void extractVerticesFromObjFile(/* input arguments */
                                char* const obj_filename,       /* obj filename, WITHOUT extension */
                                /* output arguments */
                                ch_vertex** out_vertices,       /* & of empty ch_vertex*, output vertices; out_nVert x 1 */
                                int* out_nVert);                /* & of int, number of vertices */

#ifdef __cplusplus
} /*extern "C"*/
#endif

#endif /* CONVHULL_3D_INCLUDED */

