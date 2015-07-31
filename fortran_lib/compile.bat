@echo off
setlocal
set /p fortranfile="Name of Fortan file:"%=%

for %%f in (%fortranfile%) do (
	f2py %fortranfile% -m %%~nf -h %%~nf.pyf --overwrite-signature
	f2py -c --fcompiler=g95 --compiler=mingw32 %%~nf.pyf %fortranfile%
)
