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

set "FLAG=/c /D_HAS_EXCEPTIONS=0 /MP /MP8 /TP /w /W0 /WX- /nologo /I. /Zc:wchar_t /Zc:inline /Ox /Ob2 /Oi /Ot /Oy /GL /GF /Gm- /MT /GS- /Gy /GR- /Qpar /fp:precise /fp:except- /DNDEBUG /DCUDA /DCUDA_NDEBUG /DENABLE_RELATEDNESS /D_CONSOLE /D_UNICODE /DUNICODE /std:c++20 /openmp"

::clean screen
cls

::compile CPU code
cl.exe %FLAG% vcfpop.cpp ad.cpp amova.cpp clustering.cpp conversion.cpp diff.cpp dist.cpp diversity.cpp dre.cpp file.cpp filter.cpp function.cpp global.cpp slide.cpp haplotype.cpp hash.cpp indstat.cpp kinship.cpp load.cpp math2.cpp math512.cpp mathAVX.cpp mathNEO.cpp mathSSE.cpp matrix.cpp menu.cpp misc.cpp ml.cpp ml2.cpp mlbin.cpp mom.cpp mom2.cpp parameters.cpp pcoa.cpp ploidyinfer.cpp popas.cpp relatedness.cpp spa.cpp statistics.cpp string2.cpp structure.cpp structure512.cpp structureAVX.cpp structureNEO.cpp structureSSE.cpp 

::compile CUDA code
nvcc.exe -dlink -keep -I="." -DCUDA -DCUDA_NDEBUG structureCUDA.cu -w -O3 -m64 --gpu-architecture=sm_60 -Xcompiler "/EHsc /nologo /Zi"

::link
cl.exe /Fe:vcfpop.msvc.exe /MANIFEST /LTCG:incremental /NXCOMPAT /DYNAMICBASE /NDEBUG /MACHINE:X64 /OPT:REF /INCREMENTAL:NO /SUBSYSTEM:CONSOLE /OPT:ICF /ERRORREPORT:PROMPT /NOLOGO /TLBID:1 vcfpop.obj dre.obj ml.obj ml2.obj mlbin.obj mom.obj mom2.obj global.obj hash.obj math2.obj mathNEO.obj mathSSE.obj mathAVX.obj math512.obj parameters.obj misc.obj file.obj string2.obj statistics.obj matrix.obj spa.obj ploidyinfer.obj function.obj load.obj filter.obj diversity.obj slide.obj haplotype.obj conversion.obj indstat.obj dist.obj pcoa.obj clustering.obj diff.obj kinship.obj relatedness.obj amova.obj popas.obj structure.obj structureNEO.obj structureSSE.obj structureAVX.obj structure512.obj structureCUDA.obj a_dlink.obj ad.obj menu.obj cudart.lib zlibstat.lib kernel32.lib /link %LIBDIR%

::copy and delete
copy vcfpop.msvc.exe ..\bin\vcfpop.msvc.exe
echo "Compile complete"
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
pause