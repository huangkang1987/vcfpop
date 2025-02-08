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
del a_dlink.*

::Replace the path of vcvars64.bat on your computer
call "C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Auxiliary\Build\vcvars64.bat"

::Replace the path of cuda 64bit library path on your computer
set LIBDIR="/LIBPATH:C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.5\lib\x64"

set "FLAG=/c /D_HAS_EXCEPTIONS=0 -mno-fma -mno-fma4 /MT /TP /w /W0 /WX- /Zc:wchar_t /Ob2 /Oi /O2 /Ot /Oy /GL /GF /Gm- /GS- /Gy /GR- /Qpar /arch:SSE2 /fp:precise /fp:except- /DNDEBUG /DCUDA /DCUDA_NDEBUG /DENABLE_RELATEDNESS /D_CONSOLE /D_UNICODE /DUNICODE  /std:c17 /std:c++20 /openmp -msse3 -msse4.1 -msse4.2 -mlzcnt -mpopcnt -Wwritable-strings /Iarmadillo-14.2.2\include /Igcem-master\include /Istats-master\include /Izlib

::clean screen
cls

::compile CPU code

start /min /b clang-cl.exe %FLAG% mom2.cpp
start /min /b clang-cl.exe %FLAG% ml.cpp
start /min /b clang-cl.exe %FLAG% mom.cpp 
start /min /b clang-cl.exe %FLAG% dre.cpp mlbin.cpp matrix.cpp gwas.cpp block.cpp decay.cpp
start /min /b clang-cl.exe %FLAG% vcfpop.cpp global.cpp hash.cpp math2.cpp mathNEO.cpp mathSSE.cpp parameters.cpp misc.cpp file.cpp string2.cpp statistics.cpp spa.cpp ploidyinfer.cpp function.cpp load.cpp filter.cpp diversity.cpp slide.cpp haplotype.cpp conversion.cpp indstat.cpp dist.cpp pcoa.cpp clustering.cpp diff.cpp kinship.cpp relatedness.cpp amova.cpp popas.cpp structure.cpp structureNEO.cpp structureSSE.cpp ad.cpp menu.cpp
start /min /b clang-cl.exe %FLAG% -mavx -mavx2 /arch:AVX2 mathAVX.cpp structureAVX.cpp
start /min /b clang-cl.exe %FLAG% -mavx -mavx2 -mavx512 -mavx512f -mavx512bw -mavx512dq -mavx512vl /arch:AVX512 math512.cpp structure512.cpp

::compile CUDA code
start /min /b nvcc.exe -dlink -keep -Iarmadillo-14.2.2\include -Igcem-master\include -Istats-master\include -Izlib -DCUDA -DCUDA_NDEBUG structureCUDA.cpp -w -O3 -m64 -x cu -std c++20 -gencode=arch=compute_60,code=sm_60 -gencode=arch=compute_70,code=sm_70 -gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_80,code=sm_80 -gencode=arch=compute_86,code=sm_86 -gencode=arch=compute_89,code=sm_89 -gencode=arch=compute_90,code=sm_90

clang-cl.exe %FLAG% ml2.cpp

::link
clang-cl.exe /Fe:vcfpop.exe -DSFML_STATIC vcfpop.obj dre.obj ml.obj ml2.obj mlbin.obj mom.obj mom2.obj global.obj hash.obj math2.obj mathNEO.obj mathSSE.obj mathAVX.obj math512.obj parameters.obj misc.obj file.obj string2.obj statistics.obj matrix.obj spa.obj ploidyinfer.obj function.obj load.obj filter.obj diversity.obj slide.obj haplotype.obj conversion.obj indstat.obj dist.obj pcoa.obj clustering.obj diff.obj kinship.obj relatedness.obj amova.obj popas.obj structure.obj structureNEO.obj structureSSE.obj structureAVX.obj structure512.obj structureCUDA.obj a_dlink.obj ad.obj menu.obj gwas.obj block.obj decay.obj zlibstat.lib openblas.lib cudart_static.lib kernel32.lib Advapi32.lib libomp.lib /link %LIBDIR% 

::copy and delete
copy vcfpop.exe ..\bin\vcfpop.exe
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
del a_dlink.*
pause