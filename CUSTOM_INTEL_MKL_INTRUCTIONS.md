# Building and linking a custom Intel MKL library tailored for SAF 

The Spatial_Audio_Framework (SAF) requires any library (or combination of libraries) that supports the CBLAS and LAPACK standards. Personally, we like Intel MKL for this task, as it works well and also has an optimised FFT implementation that SAF can use too. 

Unfortunetly, while you can indeed link directly to the Intel MKL static/shared libraries, these are files stupidly huge. This is why we use the MKL builder to build a custom MKL library, which contains only the routines that SAF will actually use.

Instructions for Windows, MacOSX, and Linux users rocking x86_64/amd64 CPUs, can be found below.

## Windows (64-bit) users

1. Install [Intel MKL](https://software.intel.com/en-us/articles/free-ipsxe-tools-and-libraries). 

2. Open "x64 Developer Command Prompt for VS.exe" as **administrator**.

3. Run the following batch script

```
install-safmkl.bat
```

4. Add the following header search path to your project:

```
C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/include
```

5. Add the following library search path to your project:

```
C:/Users/[YOUR WORKING DIRECTORY]/SDKs/Spatial_Audio_Framework/dependencies/Win64/lib
```

6. link your project against the following libraries (note that the second library is only needed if you built the threaded version):
```
saf_mkl_custom.lib
libiomp5md.lib
```

7. Add "SAF_USE_INTEL_MKL" to your pre-processor definitions.

Note: If you built the threaded version of the library, then there are some more pre-processors defintions you can use. See below.


## MacOSX users 

1. Install [Intel MKL](https://software.intel.com/en-us/articles/free-ipsxe-tools-and-libraries). Optionally, you may want to also add the MKL global environment variables using this command:

```
source /opt/intel/compilers_and_libraries/mac/mkl/bin/mklvars.sh intel64
```

2. The required custom library may be obtained by first copying the included "**dependencies/saf_mkl_list**" file into the MKL "builder" folder:

```
sudo cp saf_mkl_list /opt/intel/compilers_and_libraries/mac/mkl/tools/builder
```

3. EITHER (The blue pill): generate and copy the "**saf_mkl_custom.dylib**" library to a system path folder:

```
cd /opt/intel/compilers_and_libraries/mac/mkl/tools/builder
sudo make intel64 interface=lp64 threading=sequential name=libsaf_mkl_custom export=saf_mkl_list
sudo cp libsaf_mkl_custom.dylib /usr/local/lib
```

3. OR (The red pill): you may instead build a threaded version of the library with:

```
cd /opt/intel/compilers_and_libraries/mac/mkl/tools/builder
sudo make intel64 interface=lp64 threading=parallel name=libsaf_mkl_custom export=saf_mkl_list

sudo cp libsaf_mkl_custom.dylib /usr/local/lib
sudo cp /opt/intel/compilers_and_libraries/mac/compiler/lib/libiomp5.dylib /usr/local/lib
```

4. Add the following header search path to your project:

```
/opt/intel/compilers_and_libraries/mac/mkl/include
```

5. Then add the following linker flags to your project (note that the second library is only needed if you built the threaded version):

```
-L/usr/local/lib -lsaf_mkl_custom
-L/usr/local/lib -liomp5    
```

6. Finally, add "SAF_USE_INTEL_MKL" to your project's pre-processor definitions.

Note: If you built the threaded version of the library, then there are some more pre-processors defintions you can use. See below.


## Linux users

1. Install [Intel MKL](https://software.intel.com/en-us/articles/free-ipsxe-tools-and-libraries). 

2. Depending on the MKL version and how you configure the installer, the "mkl/tools/builder" folder can be in one of two places. Find out which, and keep this path consistent for further commands:

```
~/intel/compilers_and_libraries/linux/mkl/tools/builder
# Or it can be here:
/opt/intel/compilers_and_libraries/linux/mkl/tools/builder
# note, that if it is located in '/opt/', you will need to run the remaining commands with 'sudo'
```

3. Copy the included "**dependencies/saf_mkl_list**" file into the MKL "builder" folder:
```
cp saf_mkl_list ~/intel/compilers_and_libraries/linux/mkl/tools/builder
```

4. EITHER (The blue pill): generate and copy the "**saf_mkl_custom.so**" library to a system path folder:

```
cd ~/intel/compilers_and_libraries/linux/mkl/tools/builder
make intel64 interface=lp64 threading=sequential name=libsaf_mkl_custom export=saf_mkl_list
sudo cp libsaf_mkl_custom.so /usr/lib
```

4. OR (The red pill): you may instead build a threaded version of the library with:

```
cd ~/intel/compilers_and_libraries/linux/mkl/tools/builder
make intel64 interface=lp64 threading=parallel name=libsaf_mkl_custom export=saf_mkl_list

sudo cp libsaf_mkl_custom.so /usr/lib
sudo cp ~/intel/compilers_and_libraries/linux/lib/intel64/libiomp5.so /usr/lib
```

5. Add the following header search path to your project:

```
~/intel/compilers_and_libraries/linux/include
```

6. Then add the following linker flags to your project (note that the second library is only needed if you built the threaded version):

```
-L/usr/lib -lsaf_mkl_custom
-L/usr/lib -liomp5
```

7. Add "SAF_USE_INTEL_MKL" to your pre-processor definitions.


## Intel MKL threading options (Optional)

If you built a threaded version of the custom library, then there are some additional options you may specify. More information can be found [here](https://software.intel.com/en-us/articles/recommended-settings-for-calling-intel-mkl-routines-from-multi-threaded-applications).

Note by default: MKL_NUM_THREADS = [number of CPU cores], MKL_DYNAMIC = 1. 

You may also change these at run time using the following functions:

```
/* for example: */
MKL_Set_Num_Threads(2);
MKL_Set_Dynamic(1);
```
