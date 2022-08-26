call "C:\Program Files (x86)\Microsoft Visual Studio\2019\Professional\VC\Auxiliary\Build\vcvars64.bat"

cls

set "FLAG=/c /EHsc /MP /MP8 /TP /w /W0 /WX- /nologo /I"." /Zc:wchar_t /Zc:inline /Ox /Ob2 /Oi /Ot /Oy /GL /GF /Gm- /MT /GS- /Gy /GR- /Gv /Qpar /arch:AVX2 /fp:fast /fp:except- /D "NDEBUG" /D "ENABLE_RELATEDNESS" /D "_CONSOLE" /D "_UNICODE" /D "UNICODE" /std:c++17 /openmp -msse3 -msse4.1 -msse4.2 -mavx -mavx2 -mfma -mlzcnt -mpopcnt -Wwritable-strings"

clang-cl.exe %FLAG% vcfpop.cpp dre.cpp ml.cpp ml2.cpp mlbin.cpp mom.cpp mom2.cpp global.cpp hash.cpp math2.cpp mathNEO.cpp mathSSE.cpp mathAVX.cpp parameters.cpp misc.cpp file.cpp string2.cpp statistics.cpp matrix.cpp spa.cpp ploidyinfer.cpp function.cpp load.cpp filter.cpp diversity.cpp haplotype.cpp conversion.cpp indstat.cpp dist.cpp pcoa.cpp clustering.cpp diff.cpp kinship.cpp relatedness.cpp amova.cpp popas.cpp structure.cpp structureNEO.cpp structureSSE.cpp structureAVX.cpp ad.cpp menu.cpp

clang-cl.exe %FLAG% -mavx512f -mavx512bw -mavx512 /arch:AVX512 math512.cpp structure512.cpp

clang-cl.exe -MT -DSFML_STATIC vcfpop.obj dre.obj ml.obj ml2.obj mlbin.obj mom.obj mom2.obj global.obj hash.obj math2.obj mathNEO.obj mathSSE.obj mathAVX.obj math512.obj parameters.obj misc.obj file.obj string2.obj statistics.obj matrix.obj spa.obj ploidyinfer.obj function.obj load.obj filter.obj diversity.obj haplotype.obj conversion.obj indstat.obj dist.obj pcoa.obj clustering.obj diff.obj kinship.obj relatedness.obj amova.obj popas.obj structure.obj structureNEO.obj structureSSE.obj structureAVX.obj structure512.obj ad.obj menu.obj zlibstat.lib kernel32.lib libomp.lib

copy vcfpop.exe ..\bin\vcfpop.vs2019.clang.exe
del *.obj
del *.exe
pause