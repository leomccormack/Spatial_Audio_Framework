:: "Batch script to build and install a custom MKL library for SAF. Run using: x64 Command Prompt for VS 2017." 
::
@echo off
cd /d "%~dp0"

CLS

ECHO ******************************************************************************************
ECHO ************** Custom MKL library installer for the Spatial_Audio_Framework **************
ECHO ******************************************************************************************
ECHO.
ECHO This batch script will build the required saf_mkl_custom files
ECHO and copy them into:
ECHO   - "Spatial_Audio_Framework/dependencies/Win64/lib/saf_mkl_custom_[lp64|ilp64].lib" 
ECHO   - "C:/Windows/System32/saf_mkl_custom[lp64|ilp64].dll".  
ECHO You may choose between sequential and threaded versions of the library. Sequential is the
ECHO recommended option.
ECHO.
ECHO NOTE: you will need to run this script using the "x64 Command Prompt for VS.exe", which 
ECHO will require Administrator privileges!
ECHO.

:: ========================================================================= ::
:: Check that MKL is installed
IF NOT EXIST "C:\Program Files (x86)\Intel\oneAPI\mkl\latest\tools\builder" (
    echo Intel MKL not installed on this machine. Note that Intel MKL can be freely downloaded from here:
    echo "https://software.intel.com/content/www/us/en/develop/tools/oneapi/base-toolkit/download.html"
    EXIT /B
)
:: Copy saf_mkl_list to MKL builder folder
echo Copying saf_mkl_list into MKL builder folder
xcopy "saf_mkl_list" "C:\Program Files (x86)\Intel\oneAPI\mkl\latest\tools\builder"

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
ECHO 1.lp64 configuration (32-bit integers; recommended)
ECHO 2.ilp64 configuration (64-bit integers; required for safmex)
ECHO.

CHOICE /C 12 /M "Enter your choice:"

IF ERRORLEVEL 2 GOTO Parallel_ILP64   
IF ERRORLEVEL 1 GOTO Parallel_LP64


:Parallel_LP64
ECHO.
ECHO Building threaded versions of the custom MKL library (LP64 interface) for SAF...
cd C:/Program Files (x86)/Intel/oneAPI/mkl/latest/tools/builder
nmake intel64 interface=lp64 threading=parallel name=saf_mkl_custom_lp64 export=saf_mkl_list
ECHO.
ECHO Copying files to correct folders...
cd /d "%~dp0"
xcopy "C:\Program Files (x86)\Intel\oneAPI\mkl\latest\tools\builder\saf_mkl_custom_lp64.lib" "..\dependencies\Win64\lib" 
xcopy "C:\Program Files (x86)\Intel\oneAPI\mkl\latest\tools\builder\saf_mkl_custom_lp64.dll" "C:\Windows\System32\"
GOTO End

:Parallel_ILP64
ECHO.
ECHO Building threaded versions of the custom MKL library (ILP64 interface) for SAF...
cd C:/Program Files (x86)/Intel/oneAPI/mkl/latest/tools/builder
nmake intel64 interface=ilp64 threading=parallel name=saf_mkl_custom_ilp64 export=saf_mkl_list
ECHO.
ECHO Copying files to correct folders...
cd /d "%~dp0"
xcopy "C:\Program Files (x86)\Intel\oneAPI\mkl\latest\tools\builder\saf_mkl_custom_ilp64.lib" "..\dependencies\Win64\lib" 
xcopy "C:\Program Files (x86)\Intel\oneAPI\mkl\latest\tools\builder\saf_mkl_custom_ilp64.dll" "C:\Windows\System32\"
GOTO End

:Sequential
ECHO.
ECHO 1.lp64 configuration (32-bit integers; recommended)
ECHO 2.ilp64 configuration (64-bit integers; required for safmex)
ECHO.

CHOICE /C 12 /M "Enter your choice:"

IF ERRORLEVEL 2 GOTO Sequential_ILP64  
IF ERRORLEVEL 1 GOTO Sequential_LP64 

:Sequential_LP64
ECHO.
ECHO Building sequential versions of the custom MKL library (LP64 interface) for SAF...
cd C:/Program Files (x86)/Intel/oneAPI/mkl/latest/tools/builder
nmake intel64 interface=lp64 threading=sequential name=saf_mkl_custom_lp64 export=saf_mkl_list
ECHO.
ECHO Copying files to correct folders...
cd /d "%~dp0"
xcopy "C:\Program Files (x86)\Intel\oneAPI\mkl\latest\tools\builder\saf_mkl_custom_lp64.lib" "..\dependencies\Win64\lib" 
xcopy "C:\Program Files (x86)\Intel\oneAPI\mkl\latest\tools\builder\saf_mkl_custom_lp64.dll" "C:\Windows\System32\"
GOTO End

:Sequential_ILP64
ECHO.
ECHO Building sequential versions of the custom MKL library (ILP64 interface) for SAF...
cd C:/Program Files (x86)/Intel/oneAPI/mkl/latest/tools/builder
nmake intel64 interface=ilp64 threading=sequential name=saf_mkl_custom_ilp64 export=saf_mkl_list
ECHO.
ECHO Copying files to correct folders...
cd /d "%~dp0"
xcopy "C:\Program Files (x86)\Intel\oneAPI\mkl\latest\tools\builder\saf_mkl_custom_ilp64.lib" "..\dependencies\Win64\lib" 
xcopy "C:\Program Files (x86)\Intel\oneAPI\mkl\latest\tools\builder\saf_mkl_custom_ilp64.dll" "C:\Windows\System32\"
GOTO End

:End
