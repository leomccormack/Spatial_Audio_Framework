:: "Batch script to build and install a custom Intel IPP library for SAF" 
::
@echo off
cd /d "%~dp0"

CLS

ECHO ******************************************************************************************
ECHO ************** Custom IPP library installer for the Spatial_Audio_Framework **************
ECHO ******************************************************************************************
ECHO.
ECHO This batch script will build the required saf_ipp_custom library files
ECHO and copy them into:
ECHO   - "Spatial_Audio_Framework/dependencies/Win64/lib/saf_mkl_custom.lib" 
ECHO   - "C:/Windows/System32/saf_ipp_custom.dll".  
ECHO.
ECHO NOTE: you will need to run this script with Administrator privileges!
ECHO NOTE: python3 must also be installed on your system
ECHO.

:: ========================================================================= ::
:: Check that MKL is installed
IF NOT EXIST "C:\Program Files (x86)\Intel\oneAPI\ipp\latest\tools" (
    echo Intel IPP not installed on this machine. Note that Intel IPP can be freely downloaded from here:
    echo "https://software.intel.com/content/www/us/en/develop/tools/oneapi/base-toolkit/download.html"
    EXIT /B
)

:: Python 
echo Running build scrips
python "C:/Program Files (x86)/Intel/oneAPI/ipp/latest/tools/custom_library_tool_python/main.py" -c -g -n saf_ipp_custom -p "" -ff "saf_ipp_list.txt" -arch=intel64
CALL build_saf_ipp_custom_intel64.bat

:: Copying files
:: copy the generated .dll to a system PATH folder, and the .lib to somewhere local for you to link your application with, e.g.:
echo Copying files
xcopy "saf_ipp_custom.lib" "..\dependencies\Win64\lib" 
xcopy "saf_ipp_custom.dll" "C:\Windows\System32\"


:: clean everything up:
echo Cleaning up
del build_saf_ipp_custom_intel64.bat export.def main.c main.obj saf_ipp_custom.dll saf_ipp_custom.exp saf_ipp_custom.lib

GOTO End

:End
