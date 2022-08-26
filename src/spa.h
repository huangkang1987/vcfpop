/* Analysis of Spatial Structure, TEST, DO NOT USE */

#pragma once
#include "vcfpop.h"

#pragma pack(push, 1)

#pragma pack(pop)

extern _thread bool* spa_valid;
extern _thread int spa_tn;
extern double* spa_x;								//Spatial pattern, TEST, don't use, n * (dim + 1)
extern int spa_n;									//Spatial pattern, TEST, don't use, n
extern int spa_np;									//Spatial pattern, TEST, don't use, n

/* Calculate spatical pattern */
TARGET void CalcSPA();

/* Test. SPA */
TARGET double SPA_Fij(double* a, int i);

TARGET double SPA_Likelihood(CPOINT& xx, void** Param);

TARGET double SPA_Likelihood(uint64 l, double* a);

TARGET double SPA_Hessian(uint64 l, double* G, double* H, double* a, double* a2);

TARGET void SPA_Locus(uint64 l, double* x, double* f, double* a, double* a2, double* am, double* g, double* h);

/* Calculate spatial pattern using multiple threads */
THREADH(SPAThread);