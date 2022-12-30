::delete previous files
del *.exe
del *.tmp
del *.obj
del *.asm
del *.pdb
del *.ptx
del *.cudafe1.cpp
del *.c
del *.res
del *.cubin
del *.fatbin
del *.gpu
del *.ii
del *.module_id
del *.optf

::Replace the path of vcvars64.bat on your computer
call "C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Auxiliary\Build\vcvars64.bat"

::Replace the path of cuda 64bit library path on your computer
set LIBDIR="/LIBPATH:C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.8\lib\x64"

set "FLAG=/c /D_HAS_EXCEPTIONS=0 /MP /MP8 /TP /w /W0 /WX- /nologo /I"." /Zc:wchar_t /Zc:inline /Ox /Ob2 /Oi /Ot /Oy /GL /GF /Gm- /MT /GS- /Gy /GR- /Qpar /arch:SSE2 /fp:precise /fp:except- /DNDEBUG /DCUDAx /DCUDA_NDEBUG /DENABLE_RELATEDNESS /D_CONSOLE /D_UNICODE /DUNICODE /std:c++20 /openmp -msse3 -msse4.1 -msse4.2 -mlzcnt -mpopcnt -Wwritable-strings"

::clean screen
cls

::compile CPU code

start /min /b clang-cl.exe %FLAG% mom2.cpp
start /min /b clang-cl.exe %FLAG% ml.cpp
start /min /b clang-cl.exe %FLAG% mom.cpp 
start /min /b clang-cl.exe %FLAG% dre.cpp mlbin.cpp matrix.cpp
start /min /b clang-cl.exe %FLAG% vcfpop.cpp global.cpp hash.cpp math2.cpp mathNEO.cpp mathSSE.cpp parameters.cpp misc.cpp file.cpp string2.cpp statistics.cpp spa.cpp ploidyinfer.cpp function.cpp load.cpp filter.cpp diversity.cpp haplotype.cpp conversion.cpp indstat.cpp dist.cpp pcoa.cpp clustering.cpp diff.cpp kinship.cpp relatedness.cpp amova.cpp popas.cpp structure.cpp structureNEO.cpp structureSSE.cpp ad.cpp menu.cpp structureCUDA.cu
start /min /b clang-cl.exe %FLAG% -mavx -mavx2 -mfma /arch:AVX2 mathAVX.cpp structureAVX.cpp
start /min /b clang-cl.exe %FLAG% -mavx -mavx2 -mfma -mavx512f -mavx512bw -mavx512 /arch:AVX512 math512.cpp structure512.cpp

clang-cl.exe %FLAG% ml2.cpp

::link
clang-cl.exe /Fe:vcfpop.msvc.clang.exe -MT -DSFML_STATIC vcfpop.obj dre.obj ml.obj ml2.obj mlbin.obj mom.obj mom2.obj global.obj hash.obj math2.obj mathNEO.obj mathSSE.obj mathAVX.obj math512.obj parameters.obj misc.obj file.obj string2.obj statistics.obj matrix.obj spa.obj ploidyinfer.obj function.obj load.obj filter.obj diversity.obj haplotype.obj conversion.obj indstat.obj dist.obj pcoa.obj clustering.obj diff.obj kinship.obj relatedness.obj amova.obj popas.obj structure.obj structureNEO.obj structureSSE.obj structureAVX.obj structure512.obj structureCUDA.obj ad.obj menu.obj cudart.lib zlibstat.lib kernel32.lib libomp.lib /link %LIBDIR%

::copy and delete
copy vcfpop.msvc.clang.exe ..\bin\vcfpop.msvc.clang.exe
echo "Compile complete"
del *.obj
del *.tmp
del *.asm
del *.pdb
del *.ptx
del *.cudafe1.cpp
del *.c
del *.res
del *.cubin
del *.fatbin
del *.gpu
del *.ii
del *.module_id
pause