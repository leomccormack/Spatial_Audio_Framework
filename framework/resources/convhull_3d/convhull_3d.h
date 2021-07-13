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

/**
 * @file: convhull_3d.h
 * @brief An implementation of the 3-D quickhull algorithm [1]
 *
 * The code is largely derived from the "computational-geometry-toolbox" by
 * George Papazafeiropoulos (c) 2014, originally distributed under the BSD
 * (2-clause) license. Taken from: https://github.com/leomccormack/convhull_3d
 *
 * ## Dependencies
 *   CBLAS (optional) for speed ups, especially for very large meshes
 *
 * @see [1] C. Bradford, Barber, David P. Dobkin and Hannu Huhdanpaa, "The
 *          Quickhull Algorithm for Convex Hull". Geometry Center Technical
 *          Report GCG53, July 30, 1993
 *
 * @author Leo McCormack
 * @date 02.10.2017
 * @license MIT
 */

#ifndef CONVHULL_3D_INCLUDED
#define CONVHULL_3D_INCLUDED

#if !defined(__cplusplus) && defined(_MSC_VER)
# pragma warning(disable : 4201)
#endif

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifdef CONVHULL_3D_USE_FLOAT_PRECISION
typedef float CH_FLOAT;
#else
typedef double CH_FLOAT;
#endif
/**
 * vertex structure, used by convhull_3d
 */
typedef struct _ch_vertex {
    union{
        struct{
            CH_FLOAT x, y, z;
        };
        CH_FLOAT v[3];
    };
} ch_vertex;
typedef ch_vertex ch_vec3;

/**
 * Builds the 3-D convexhull using the quickhull algorithm [1]
 *
 * @param[in]  in_vertices Vector of input vertices; nVert x 1
 * @param[in]  nVert       Number of vertices
 * @param[out] out_faces   (&) output face indices; FLAT: nOut_faces x 3
 * @param[out] out_cf      (&) contains the coefficients of the planes (set to
 *                         NULL if not wanted); FLAT: nOut_faces x 3
 * @param[out] out_df      (&) contains the constant terms of the planes (set to
 *                         NULL if not wanted); nOut_faces x 1
 * @param[out] nOut_faces  (&) number of output face indices
 *
 * @see [1] C. Bradford, Barber, David P. Dobkin and Hannu Huhdanpaa, "The
 *          Quickhull Algorithm for Convex Hull". Geometry Center Technical
 *          Report GCG53, July 30, 1993
 */
void convhull_3d_build(/* Input arguments */
                       ch_vertex* const in_vertices,
                       const int nVert,
                       /* Output arguments */
                       int** out_faces,
                       CH_FLOAT** out_cf,
                       CH_FLOAT** out_df,
                       int* nOut_faces);

/**
 * Builds the N-D convexhull using the quickhull algorithm [1]
 *
 * @param[in]  in_vertices Matrix of nVertices in 'd' dimensions; FLAT:nVert x d
 * @param[in]  nVert       Number of vertices
 * @param[in]  d           Number of dimensions
 * @param[out] out_faces   (&) output face indices; FLAT: nOut_faces x d
 * @param[out] out_cf      (&) contains the coefficients of the planes (set to
 *                         NULL if not wanted); FLAT: nOut_faces x d
 * @param[out] out_df      (&) contains the constant terms of the planes (set to
 *                         NULL if not wanted); nOut_faces x 1
 * @param[out] nOut_faces  (&) number of output face indices
 *
 * @see [1] C. Bradford, Barber, David P. Dobkin and Hannu Huhdanpaa, "The
 *          Quickhull Algorithm for Convex Hull". Geometry Center Technical
 *          Report GCG53, July 30, 1993
 */
void convhull_nd_build(/* Input arguments */
                       CH_FLOAT* const in_vertices,
                       const int nVert,
                       const int d,
                       /* Output arguments */
                       int** out_faces,
                       CH_FLOAT** out_cf,
                       CH_FLOAT** out_df,
                       int* nOut_faces);

/**
 * Exports the vertices, face indices, and face normals, as an '.obj' file, ready
 * for the GPU.
 *
 * @param[in] vertices                 Vector of input vertices; nVert x 1
 * @param[in] nVert                    Number of vertices
 * @param[in] faces                    Face indices; flat: nFaces x 3
 * @param[in] nFaces                   Number of faces in hull
 * @param[in] keepOnlyUsedVerticesFLAG '0' exports in_vertices, '1': exports
 *                                     Only used vertices
 * @param[in] obj_filename             *.obj filename, WITHOUT extension
 */
void convhull_3d_export_obj(/* Input arguments */
                            ch_vertex* const vertices,
                            const int nVert,
                            int* const faces,
                            const int nFaces,
                            const int keepOnlyUsedVerticesFLAG,
                            char* const obj_filename);

/**
 * Exports the vertices, face indices, and face normals, as an '.m' file, for
 * Matlab verification
 *
 * @param[in] vertices    Vector of input vertices; nVert x 1
 * @param[in] nVert       Number of vertices
 * @param[in] faces       Face indices; flat: nFaces x 3
 * @param[in] nFaces      Number of faces in hull
 * @param[in] m_filename  *.m filename, WITHOUT extension
 */
void convhull_3d_export_m(/* Input arguments */
                          ch_vertex* const vertices,
                          const int nVert,
                          int* const faces,
                          const int nFaces,
                          char* const m_filename);

/**
 * Reads an '.obj' file and extracts only the vertices.
 *
 * @param[in]  obj_filename *.obj filename, WITHOUT extension
 * @param[out] out_vertices (&) output vertices; out_nVert x 1
 * @param[out] out_nVert    (&) number of vertices
 */
void extractVerticesFromObjFile(/* Input arguments */
                                char* const obj_filename,
                                /* Output arguments */
                                ch_vertex** out_vertices,
                                int* out_nVert);        

#ifdef __cplusplus
} /*extern "C"*/
#endif /* __cplusplus */

#endif /* CONVHULL_3D_INCLUDED */

