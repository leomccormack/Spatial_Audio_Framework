# Performance libraries supported by SAF

The Spatial_Audio_Framework (SAF) requires any library (or combination of libraries) which supports the CBLAS and LAPACK standards. 

Currently, SAF supports [Intel MKL](https://software.intel.com/en-us/articles/free-ipsxe-tools-and-libraries), [OpenBLAS](https://github.com/xianyi/OpenBLAS), [Apple Accelerate](https://developer.apple.com/documentation/accelerate), and [ALTAS](http://math-atlas.sourceforge.net/). You may select which option you wish to use by adding one of the following pre-processor definitions:

```
SAF_USE_INTEL_MKL
SAF_USE_OPEN_BLAS_AND_LAPACKE
SAF_USE_APPLE_ACCELERATE
SAF_USE_ATLAS
```

**MacOSX users only**: if you do not define one of the above flags, then SAF will automatically use [Apple Accelerate](https://developer.apple.com/documentation/accelerate) for CBLAS/LAPACK and also vDSP for the FFT. Although, it should be noted that Intel MKL is still the more recommended option, as it is generally faster than Apple Accelerate.


## OPTION: SAF_USE_INTEL_MKL

Intel MKL is perhaps the fastest library for **x86/amd64** platforms. It also has an optimised FFT implementation which SAF is able to use too. 

Note Intel MKL can be freely downloaded (royalty-free) [from here](https://software.intel.com/en-us/articles/free-ipsxe-tools-and-libraries).

It also ships with recent [Anaconda](https://anaconda.org) distributions. The easiest way to install with conda is:
```
conda install mkl mkl-include
```
Please refer to the documentation about where conda installs packages on your platform.

Unfortunately, while you can indeed link directly to the Intel MKL static/shared libraries... these files are rather large. Therefore, you may prefer to use the MKL builder to build a custom MKL library, which contains only the routines that SAF uses.

Instructions regarding how to build a custom MKL library for Windows, MacOSX, and Linux users, can be found below.

### MacOSX users 

Optionally, you may want to add the MKL global environment variables using this command:

```
source /opt/intel/compilers_and_libraries/mac/mkl/bin/mklvars.sh intel64
```

Copy the included **saf_mkl_list** file into the MKL **builder** folder:

```
sudo cp saf_mkl_list /opt/intel/compilers_and_libraries/mac/mkl/tools/builder
```

**EITHER** (The blue pill): generate and copy the **saf_mkl_custom.dylib** library to a system path folder:

```
cd /opt/intel/compilers_and_libraries/mac/mkl/tools/builder
sudo make intel64 interface=lp64 threading=sequential name=libsaf_mkl_custom export=saf_mkl_list
sudo cp libsaf_mkl_custom.dylib /usr/local/lib
```

**OR** (The red pill): you may instead build a threaded version of the library with:

```
cd /opt/intel/compilers_and_libraries/mac/mkl/tools/builder
sudo make intel64 interface=lp64 threading=parallel name=libsaf_mkl_custom export=saf_mkl_list

sudo cp libsaf_mkl_custom.dylib /usr/local/lib
sudo cp /opt/intel/compilers_and_libraries/mac/compiler/lib/libiomp5.dylib /usr/local/lib
```

Add the following header search path to your project:

```
/opt/intel/compilers_and_libraries/mac/mkl/include
```

Then add the following linker flag to your project:

```
-L/usr/local/lib -lsaf_mkl_custom 
```

### Linux users

Depending on the MKL version and how you configure the installer, the **mkl/tools/builder** folder can be in one of two places. Find out which, and keep this path consistent for further commands:

```
~/intel/compilers_and_libraries/linux/mkl/tools/builder
# Or it can be here:
/opt/intel/compilers_and_libraries/linux/mkl/tools/builder
# note, that if it is located in /opt/, then you will need to run the remaining commands with 'sudo'
```

Copy the included **saf_mkl_list** file into the MKL **builder** folder:
```
cp saf_mkl_list ~/intel/compilers_and_libraries/linux/mkl/tools/builder
```

**EITHER** (The blue pill): generate and copy the **saf_mkl_custom.so** library to a system path folder:

```
cd ~/intel/compilers_and_libraries/linux/mkl/tools/builder
make intel64 interface=lp64 threading=sequential name=libsaf_mkl_custom export=saf_mkl_list
sudo cp libsaf_mkl_custom.so /usr/lib
```

**OR** (The red pill): you may instead build a threaded version of the library with:

```
cd ~/intel/compilers_and_libraries/linux/mkl/tools/builder
make intel64 interface=lp64 threading=parallel name=libsaf_mkl_custom export=saf_mkl_list

sudo cp libsaf_mkl_custom.so /usr/lib
sudo cp ~/intel/compilers_and_libraries/linux/lib/intel64/libiomp5.so /usr/lib
```

Add the following header search path to your project:

```
~/intel/compilers_and_libraries/linux/mkl/include
```

Then add the following linker flag to your project:

```
-L/usr/lib -lsaf_mkl_custom
```

### Windows (64-bit) users

Open **x64 Developer Command Prompt for VS.exe** as **administrator**.

Run the following batch script

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

link your project against the following libraries (note that the second library is only needed if you built the threaded version):
```
saf_mkl_custom.lib
libiomp5md.lib
```

## OPTION: SAF_USE_OPEN_BLAS_AND_LAPACKE

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

## OPTION: SAF_USE_APPLE_ACCELERATE

Only for MacOSX. This is also the default performance library for Apple users if no performance library is specified. Simply add the [**Accelerate**](https://developer.apple.com/documentation/accelerate) framework to your XCode project, and you're good to go.


## OPTION: SAF_USE_ATLAS

Using [ALTAS](http://math-atlas.sourceforge.net/) is not recommended, since some LAPACK functions are missing. However, if you don't mind losing some SAF functionality, then this may still be a good choice for your particular project.
