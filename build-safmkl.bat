:: "Batch script to build and install a custom MKL library for SAF. Run using: x64 Command Prompt for VS 2017." 
::
@echo off
cd /d "%~dp0"

CLS

ECHO ******************************************************************************************
ECHO ************** Custom MKL library installer for the Spatial_Audio_Framework **************
ECHO ******************************************************************************************
ECHO This batch script will build the required saf_mkl_custom.lib and saf_mkl_custom.dll files
ECHO and copy them into:
ECHO "Spatial_Audio_Framework/dependencies/Win64/lib/saf_mkl_custom.lib" 
ECHO "C:/Windows/System32/saf_mkl_custom.dll".  
ECHO You may choose between sequential and threaded versions of the library. The latter option
ECHO will also place libiomp5md.lib and libiomp5md.dll, in these same folders.
ECHO (Note, you will likely need to run this script in Administrator mode)

:: ========================================================================= ::
:: Check that MKL is installed
IF NOT EXIST "C:\Program Files (x86)\IntelSWTools\compilers_and_libraries\windows\mkl\tools\builder" (
    echo Intel MKL not installed on this machine. Note that Intel MKL can be freely acquired from here:
    echo "https://software.intel.com/en-us/articles/free-ipsxe-tools-and-libraries"
    EXIT /B
)
:: Copy saf_mkl_list to MKL builder folder
echo Copying saf_mkl_list into MKL builder folder
xcopy "dependencies\saf_mkl_list" "C:\Program Files (x86)\IntelSWTools\compilers_and_libraries\windows\mkl\tools\builder"


:: ========================================================================= ::
::CLS
ECHO.
ECHO 1.Sequential (recommended)
ECHO 2.Threaded/parallel
ECHO.

CHOICE /C 12 /M "Enter your choice:"

:: Note - list ERRORLEVELS in decreasing order
IF ERRORLEVEL 2 GOTO Parallel  
IF ERRORLEVEL 1 GOTO Sequential

:Parallel
ECHO.
ECHO Building threaded versions of the custom MKL library for SAF...
cd C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/tools/builder
nmake intel64 interface=lp64 threading=parallel name=saf_mkl_custom export=saf_mkl_list
ECHO.
ECHO Copying files to correct folders...
cd /d "%~dp0"
xcopy "C:\Program Files (x86)\IntelSWTools\compilers_and_libraries\windows\mkl\tools\builder\saf_mkl_custom.lib" "dependencies\Win64\lib" 
xcopy "C:\Program Files (x86)\IntelSWTools\compilers_and_libraries\windows\mkl\tools\builder\saf_mkl_custom.dll" "C:\Windows\System32\"
xcopy "C:\Program Files (x86)\IntelSWTools\compilers_and_libraries\windows\compiler\lib\intel64\libiomp5md.lib" "dependencies\Win64\lib"
xcopy "C:\Program Files (x86)\IntelSWTools\compilers_and_libraries\windows\redist\intel64\compiler\libiomp5md.dll" "C:\Windows\System32\"
GOTO End

:Sequential
ECHO.
ECHO Building sequential versions of the custom MKL library for SAF...
cd C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/tools/builder
nmake intel64 interface=lp64 threading=sequential name=saf_mkl_custom export=saf_mkl_list
ECHO.
ECHO Copying files to correct folders...
cd /d "%~dp0"
xcopy "C:\Program Files (x86)\IntelSWTools\compilers_and_libraries\windows\mkl\tools\builder\saf_mkl_custom.lib" "dependencies\Win64\lib" 
xcopy "C:\Program Files (x86)\IntelSWTools\compilers_and_libraries\windows\mkl\tools\builder\saf_mkl_custom.dll" "C:\Windows\System32\"
GOTO End

:End
