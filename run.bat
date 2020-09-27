@echo off
cd %~dp0

SET REV=R2018a& SET VER=v94& SET DLLVER=9_4
SET REV=R2018b& SET VER=v95& SET DLLVER=9_5
SET REV=R2019a& SET VER=v96& SET DLLVER=9_6
SET REV=R2020a& SET VER=v98& SET DLLVER=9_8

SET CYGWIN=C:\cygwin64\bin
SET RUNTIME=C:\Program Files\MATLAB\MATLAB Runtime\%VER%
SET LIBNAME=MMSE
:SET PATH=%RUNTIME%\runtime\win64;%CYGWIN%
SET PATH=%RUNTIME%\runtime\win64;%RUNTIME%\bin\win64;%CYGWIN%

SET EXE=mmse.exe

REM =====================================

echo ----------
echo [%~nx0]
echo CYGWIN  : %CYGWIN%
echo RUNTIME : %RUNTIME%
echo PATH    : %PATH%
echo EXE     : %EXE%
echo ----------

REM =====================================

%EXE%

REM =====================================

