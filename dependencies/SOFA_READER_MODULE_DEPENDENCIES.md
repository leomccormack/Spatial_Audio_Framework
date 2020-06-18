
# Dependencies for the saf_sofa_reader module 

In order to use the built-in [SOFA](https://www.sofaconventions.org/mediawiki/index.php/SOFA_(Spatially_Oriented_Format_for_Acoustics)) reader module (**saf_sofa_reader**), your project must also link against the [netCDF](https://www.unidata.ucar.edu/software/netcdf/) library (including its dependencies). 

Note that the following preprocessor definition is also required to enable this optional module:

```
SAF_ENABLE_SOFA_READER_MODULE
```

Platform specific instructions on how to obtain and link these dependent libraries is provided below.


## Windows MSVC (64-bit) and MacOSX users

For convenience, the following statically built libraries are included in the **dependencies** folder; simply link your project against them:

```
libszip.lib; libzlib.lib; libhdf5.lib; libhdf5_hl.lib; netcdf.lib; # Win64
netcdf; hdf5; hdf5_hl; z; # MacOSX
```

Make sure to also add the appropriate 'include' and 'lib' directories to your project's header and library search paths, respectively.


## Linux (amd64) and Raspberry Pi (ARM) users

For Ubuntu/Debian based distros, you may install netCDF and its dependencies with:

```
sudo apt-get install libhdf5-dev libnetcdf-dev libnetcdff-dev
```

Add the directory of the **netcdf.h** file to your project's header search paths, and then add this linker flag:
```
-L/lib/x86_64-linux-gnu -lnetcdf  # (or wherever it was installed) 
```

## Windows MSYS2 (64-bit) users

The required libraries may be installed with:
```
pacman -S mingw-w64-x86_64-hdf5 mingw-w64-x86_64-netcdf
```