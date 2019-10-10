(https://github.com/leomccormack/convhull_3d)

# convhull_3d

A header only C implementation of the 3-D Quickhull algorithm for building Convex Hulls. The code is also MSVC-C89 and C++ compiler compliant.

![](images/teapot_example.png)

## Getting Started

To use this 3-D Convex Hull implementation in a '.c' or '.cpp' file, first add the following:

```c
#define CONVHULL_3D_ENABLE
#include "convhull_3d.h"
```

Then specify the vertices, which can be optionally extracted from an '.obj' file using the following code:

```c
ch_vertex* vertices = NULL;
int nVertices;
extractVerticesFromObjFile(OBJ_FILE_NAME, &vertices, &nVertices);
/* Where 'vertices' is a vector of vertices [nVertices].x, .y, .z  */
```

where 'OBJ_FILE_NAME' is the '.obj' file path (without the extension).

Or they may be defined manually, for example: a random distribution of points on the unit sphere:

```c
int n = 936;
ch_vertex* vertices;
vertices = (ch_vertex*)malloc(n*sizeof(ch_vertex));
for (i = 0; i < n; i++) {
    float elev = rand()/(float)RAND_MAX * M_PI * 2.0;
    float azi = rand()/(float)RAND_MAX * M_PI * 2.0;
    vertices[i].x = cos(azi) * cos(elev) * rand()/(float)RAND_MAX;
    vertices[i].y = sin(azi) * cos(elev) * rand()/(float)RAND_MAX;
    vertices[i].z = sin(elev);
}
```

The Convex Hull may then be built and subsequently exported (including face normals) as an '.obj' file, using this code:

```c
int* faceIndices = NULL;
int nFaces;
convhull_3d_build(vertices, nVertices, &faceIndices, &nFaces);
/* Where 'faceIndices' is a flat 2D matrix [nFaces x 3] */
convhull_3d_export_obj(vertices, nVertices, faceIndices, nFaces, 1, OUTPUT_OBJ_FILE_NAME);
free(vertices);
free(faceIndices);
```

where 'OUTPUT_OBJ_FILE_NAME' is the output '.obj' file path (without the extension).

### Additional options

By default, the implementation uses double floating point precision to build the hull, while still exporting the results in single floating point precision. However, one may configure convhull_3d to use single precision to build the hull (which is less accurate and reliable, but quicker) by adding the following:
```c
#define CONVHULL_3D_USE_FLOAT_PRECISION /* (optional) */
#define CONVHULL_3D_ENABLE
#include "convhull_3d.h"
```

If your project has CBLAS linked, then you may also speed up the matrix multiplications by adding:
```c
#define CONVHULL_3D_USE_CBLAS /* (optional) */
#define CONVHULL_3D_USE_FLOAT_PRECISION /* (optional) */
#define CONVHULL_3D_ENABLE
#include "convhull_3d.h"
```

## Test

This repository contains files: 'test/test_convhull_3d.c' and 'test/test_script.m'. The former can be used to generate Convex Hulls of the '.obj' files located in the 'test/obj_files' folder, which can be subsequently verified in MatLab using the latter file; where the 'convhull_3d.h' implementation is compared with MatLab's built-in 'convhull' function, side-by-side. Furthermore, Visual Studio 2017 and Xcode project files have been included in the 'test' folder for convenience.

![](images/tdesign_5100_sph.png)
![](images/teapot_matlab.png)

## Examples

The 'test/test_convhull_3d.c' file may also serve as example usage of the convhull_3d implementation. The following images are of the original 'obj' files (left) and the corresponding Convex Hulls (right), depicted using Tim Maxwell's OBJ Viewer:

![](images/violincase_example.png)
![](images/sandal_example.png)
![](images/trumpet_example.png)

Convex Hulls of uniformly distributed points on a sphere (180, 840, 5100 points left-to-right):

![](images/sph_tdesigns.png)

## Future work

Possibly add the Delaunay triangulation and Voronoi diagram algorithms, as the additional effort to implement them is relatively low.

## License

The code is distributed under the MIT license, but contains code that was originally written for MatLab by George Papazafeiropoulos (c) 2014; which was distributed under the BSD (2-clause) license and can be found [here](https://se.mathworks.com/matlabcentral/fileexchange/48509-computational-geometry-toolbox?focused=3851286&tab=function).

## Contact

If you have any questions, or encounter any bugs, please email: leo.mccormack@aalto.fi
