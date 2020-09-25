@echo off
cd %~dp0

SET REV=R2018a& SET VER=v94& SET DLLVER=9_4
SET REV=R2018b& SET VER=v95& SET DLLVER=9_5
SET REV=R2019a& SET VER=v96& SET DLLVER=9_6

SET CYGWIN=C:\cygwin64\bin
SET RUNTIME=C:\Program Files\MATLAB\MATLAB Runtime\%VER%
SET LIBNAME=MMSE
SET PATH=%RUNTIME%\runtime\win64;%CYGWIN%

REM =====================================

echo ----------
echo [%~nx0]
echo CYGWIN  : %CYGWIN%
echo RUNTIME : %RUNTIME%
echo PATH    : %PATH%
echo ----------

REM =====================================

SET BUILD=x86_64-w64-mingw32-gcc.exe  ^
 -Wall  -O2 ^
 "%RUNTIME%\runtime\win64\libMatlabCppSharedLib%DLLVER%.dll" ^
 "%RUNTIME%\runtime\win64\mclmcrrt%DLLVER%.dll" ^
 -I "%RUNTIME%\extern\include" ^
    "%RUNTIME%\bin\win64\libmwblas.dll" ^
    "%RUNTIME%\bin\win64\libmwlapack.dll"

REM =====================================

%BUILD% -o mmse.exe ^
           mmse.c
