# Performance Libraries Supported by SAF

The Spatial_Audio_Framework (SAF) requires any library (or combination of libraries) which supports the [CBLAS](https://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms#Implementations) and [LAPACK](https://en.wikipedia.org/wiki/LAPACK) standards

Currently, SAF supports [Intel MKL](https://software.intel.com/en-us/articles/free-ipsxe-tools-and-libraries), [OpenBLAS](https://github.com/xianyi/OpenBLAS), [Apple Accelerate](https://developer.apple.com/documentation/accelerate), and [ATLAS](http://math-atlas.sourceforge.net/). You may select which option you wish to use by adding one of the following pre-processor definitions:

```
SAF_USE_INTEL_MKL_LP64        # great option, but only for x86 architectures (using the LP64 config [int32])
SAF_USE_INTEL_MKL_ILP64       # great option, but only for x86 architectures (using the ILP64 config [int64])  
SAF_USE_OPEN_BLAS_AND_LAPACKE # good option, works on everything
SAF_USE_APPLE_ACCELERATE      # good option (x86 and ARM), faster than OpenBLAS, but MacOSX only & slower than MKL
SAF_USE_ATLAS                 # bad option (x86 and ARM), many LAPACK functions are missing
SAF_USE...                    # please get in touch if you use something else! :-)
```

## SAF_USE_INTEL_MKL_LP64 or SAF_USE_INTEL_MKL_ILP64

Intel MKL is perhaps the fastest library for **x86** platforms. It also includes an optimised FFT implementation and a number of additional vectorised utility functions which SAF is also able to make use of. 

**Intel MKL can be downloaded (royalty-free) [from here](https://software.intel.com/en-us/articles/free-ipsxe-tools-and-libraries)**.

Unfortunately, while you can indeed link directly to the Intel MKL static/shared libraries... these files are rather large. Therefore, you may prefer to use the MKL builder to build a custom MKL library, which contains only the routines that SAF actually needs.

Instructions regarding how to build a custom MKL library for Windows, MacOSX, and Linux users, can be found below.

### Linux/MacOSX users 

Run the following bash script (**sudo** privileges required):

```
cd scripts
sudo ./install-safmkl.sh [threaded|sequential] [lp64|ilp64] [path/to/mkl/tools/builder]
```

Add the following header search path to your project:

```
[/opt|~]/intel/oneapi/mkl/latest/include/                   # The new oneapi path
[/opt|~]/intel/compilers_and_libraries/linux/mkl/include    # Old path for Linux users
[/opt|~]/intel/compilers_and_libraries/mac/mkl/include      # Old path for MacOSX users
```

Then add the following linker flag to your project:

```
-L/usr/local/lib -lsaf_mkl_custom_lp64   # if you built the lp64 (32-bit integer) version
-L/usr/local/lib -lsaf_mkl_custom_ilp64  # if you built the ilp64 (64-bit integer) version
``` 

Finally, add the appropriate global definition:
```
SAF_USE_INTEL_MKL_LP64      # if you are linking with the lp64 (32-bit integer) version
SAF_USE_INTEL_MKL_ILP64     # if you are linking with the ilp64 (64-bit integer) version
```

### Windows users
 
Run the following batch script using **x64 Developer Command Prompt for VS.exe** (open as **administrator**):

```
cd scripts
install-safmkl.bat
```

Add the following header search path to your project:

```
C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/include
```

Add the following library search path to your project:

```
C:/YOUR/WORKING/DIRECTORY/Spatial_Audio_Framework/dependencies/Win64/lib
```

Link your project against the following library:
```
saf_mkl_custom_lp64.lib    # if you built the lp64 (32-bit integer) version
saf_mkl_custom_ilp64.lib   # if you built the ilp64 (64-bit integer) version
```

Finally, add the appropriate global definition:
```
SAF_USE_INTEL_MKL_LP64      # if you are linking with the lp64 (32-bit integer) version
SAF_USE_INTEL_MKL_ILP64     # if you are linking with the ilp64 (64-bit integer) version
```

### Other options

This custom MKL library (lp64) is also installed alongside the [SPARTA](http://research.spa.aalto.fi/projects/sparta_vsts/) VST plug-in suite. Athough, it will be noted that this may not be the most up-to-date version of the library.

You may of course also elect to not build this custom library and link directly, e.g. with:
```
/opt/intel/oneapi/mkl/latest/lib/libmkl_rt.[dylib|so]
```

This will select the interface and threading at runtime, based on environment variables.
You can check the options for your system with intel's [link-line-advisor](https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onemkl/link-line-advisor.html).
The options can be passed directly to CMake, e.g. with
```
-DSAF_PERFORMANCE_LIB=SAF_USE_INTEL_MKL_ILP64 -DINTEL_MKL_HEADER_PATH="/opt/intel/oneapi/mkl/latest/include" -DINTEL_MKL_LIB="/opt/intel/oneapi/mkl/latest/lib/libmkl_intel_ilp64.dylib;/opt/intel/oneapi/mkl/latest/lib/libmkl_sequential.dylib;/opt/intel/oneapi/mkl/latest/lib/libmkl_core.dylib"
```

Intel MKL also ships with recent [Anaconda](https://anaconda.org) distributions. The easiest way to install with conda is:
```
conda install mkl mkl-include
```
Please refer to the documentation to find out where conda installs packages on your platform.


## SAF_USE_OPEN_BLAS_AND_LAPACKE

The framework also supports [OpenBLAS](https://github.com/xianyi/OpenBLAS). However, unlike Intel MKL and Apple Accelerate, the OpenBLAS library does not offer an optimised DFT/FFT, and some vector-vector operations; such as element-wise multiplications and additions. Therefore, consider pairing OpenBLAS with **SAF_ENABLE_SIMD** (enabling SSE3 and AVX2) and **SAF_USE_INTEL_IPP** or **SAF_USE_FFTW**, in which case, the performance of SAF will become roughly on par with e.g. Apple Accelerate; but not quite as fast as Intel MKL.

For Debian/Ubuntu based Linux distributions the required libraries may be installed via:
```
sudo apt-get install liblapack3 liblapack-dev libopenblas-base libopenblas-dev liblapacke-dev
```
While Apple Accelerate is the more recommended option for MacOSX users, OpenBLAS may be installed and used with the following:
```
brew install openblas
export LDFLAGS="-L/usr/local/opt/openblas/lib"
export CPPFLAGS="-I/usr/local/opt/openblas/include"
```

Windows users can follow the instuctions found [here](https://github.com/xianyi/OpenBLAS/wiki/How-to-use-OpenBLAS-in-Microsoft-Visual-Studio).

Lastly, you may also instead build the required libraries yourself via CMake. For example:

```
git clone https://github.com/xianyi/OpenBLAS.git
# install a fortran compiler if you haven't already (for lapack support)
sudo apt-get install gfortran
# build static libs
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=0 -DBUILD_WITHOUT_LAPACK=0
make 
```

## SAF_USE_APPLE_ACCELERATE

Only for MacOSX users. Simply link against the [**Accelerate**](https://developer.apple.com/documentation/accelerate) framework and you're good to go. 

Note that Accelerate also includes an optimised DFT/FFT implementation and a number of additional vectorised utility functions which SAF is also able to make use of. 


## SAF_USE_ATLAS

Using [ATLAS](http://math-atlas.sourceforge.net/) is not recommended, since some LAPACK functions are missing. However, if you don't mind losing some SAF functionality, then this may still be a good choice for your particular project.
