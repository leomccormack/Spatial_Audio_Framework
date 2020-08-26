# Performance libraries supported by SAF

The Spatial_Audio_Framework (SAF) requires any library (or combination of libraries) which supports the [CBLAS](https://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms#Implementations) and [LAPACK](https://en.wikipedia.org/wiki/LAPACK) standards

Currently, SAF supports [Intel MKL](https://software.intel.com/en-us/articles/free-ipsxe-tools-and-libraries), [OpenBLAS](https://github.com/xianyi/OpenBLAS), [Apple Accelerate](https://developer.apple.com/documentation/accelerate), and [ATLAS](http://math-atlas.sourceforge.net/). You may select which option you wish to use by adding one of the following pre-processor definitions:

```
SAF_USE_INTEL_MKL             # great option, but only for x86 architectures    
SAF_USE_OPEN_BLAS_AND_LAPACKE # good option, works on everything
SAF_USE_APPLE_ACCELERATE      # good option (x86 and ARM), faster than OpenBLAS, but MacOSX only & slower than MKL
SAF_USE_ATLAS                 # bad option (x86 and ARM), many LAPACK functions are missing
SAF_USE...                    # we're also open to adding more alternatives; please get in touch if you have one! :-)
```

## SAF_USE_INTEL_MKL

Intel MKL is perhaps the fastest library for **x86** platforms. It also includes an optimised FFT implementation and a number of additional vectorised utility functions which SAF is also able to make use of. 

**Intel MKL can be downloaded (royalty-free) [from here](https://software.intel.com/en-us/articles/free-ipsxe-tools-and-libraries)**.

Unfortunately, while you can indeed link directly to the Intel MKL static/shared libraries... these files are rather large. Therefore, you may prefer to use the MKL builder to build a custom MKL library, which contains only the routines that SAF actually needs.

Instructions regarding how to build a custom MKL library for Windows, MacOSX, and Linux users, can be found below.

### Linux/MacOSX users 

Run the following bash script (**sudo** privileges required):

```
sudo ./install-safmkl.sh [threaded|sequential]
```

Add the following header search path to your project:

```
[/opt|~]/intel/compilers_and_libraries/linux/mkl/include    # Linux users
[/opt|~]/intel/compilers_and_libraries/mac/mkl/include      # MacOSX users
```

Then add the following linker flag to your project:

```
-L/usr/lib -lsaf_mkl_custom        # Linux users
-L/usr/local/lib -lsaf_mkl_custom  # MacOSX users
``` 

### Windows users
 
Run the following batch script using **x64 Developer Command Prompt for VS.exe** (open as **administrator**):

```
install-safmkl.bat
```

Add the following header search path to your project:

```
C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/include
```

Add the following library search path to your project:

```
C:/Users/[YOUR WORKING DIRECTORY]/SDKs/Spatial_Audio_Framework/dependencies/Win64/lib
```

link your project against the following library:
```
saf_mkl_custom.lib 
```

### Other options

This custom MKL library is also installed alongside the [SPARTA](http://research.spa.aalto.fi/projects/sparta_vsts/) VST plug-in suite. Athough, it will be noted that this may not be the most up-to-date version of the library.

You may of course also elect to not build this custom library and link directly with:
```
[/opt|~]/intel/compilers_and_libraries/linux/mkl/lib/mkl_rt.so                              # Linux users
[/opt|~]/intel/compilers_and_libraries/mac/mkl/lib/mkl_rt.dylib                             # MacOSX users
C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/lib/libmkl_rt.lib   # Windows users
```

Intel MKL also ships with recent [Anaconda](https://anaconda.org) distributions. The easiest way to install with conda is:
```
conda install mkl mkl-include
```
Please refer to the documentation to find out where conda installs packages on your platform.


## SAF_USE_OPEN_BLAS_AND_LAPACKE

For those who would like to use the framework with [OpenBLAS](https://github.com/xianyi/OpenBLAS), note that for Debian/Ubuntu based Linux distributions, the libraries may be installed via:

```
sudo apt-get install liblapack3 liblapack-dev libopenblas-base libopenblas-dev liblapacke-dev
```

However, if you are not on Linux, or would prefer to have more control over things, then the required libraries may also be built via CMake:

```
git clone https://github.com/xianyi/OpenBLAS.git
# install fortran compiler if you haven't already (for lapack support)
sudo apt-get install gfortran
# build static libs
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=0 -DBUILD_WITHOUT_LAPACK=0
make 
```

The static openblas library can then be found in:
```
OpenBLAS/build/lib/libopenblas.a
```

The required include files are then found in:
```
OpenBLAS/lapack-netlib/CBLAS/include
OpenBLAS/lapack-netlib/LAPACKE/include
```

## SAF_USE_APPLE_ACCELERATE

Only for MacOSX users. Simply add the [**Accelerate**](https://developer.apple.com/documentation/accelerate) framework to your XCode project, and you're good to go. 

Note that Accelerate also includes an optimised FFT implementation and a number of additional vectorised utility functions which SAF is also able to make use of. 


## SAF_USE_ATLAS

Using [ATLAS](http://math-atlas.sourceforge.net/) is not recommended, since some LAPACK functions are missing. However, if you don't mind losing some SAF functionality, then this may still be a good choice for your particular project.
