/* Analysis of Spatial Structure, TEST, DO NOT USE */

#pragma once
#include "vcfpop.h"

#pragma pack(push, 1)

#pragma pack(pop)

extern thread_local bool* spa_valid;
extern thread_local int spa_tn;
extern double* spa_x;								//Spatial pattern, TEST, don't use, n * (dim + 1)
extern int spa_n;									//Spatial pattern, TEST, don't use, n
extern int spa_np;									//Spatial pattern, TEST, don't use, n

/* Calculate spatical pattern */
template<typename REAL>
TARGET void CalcSPA();

/* Test. SPA */
template<typename REAL>
TARGET double SPA_Fij(double* a, int i);

template<typename REAL>
TARGET double SPA_Likelihood(void* Param, CPOINT& xx, rmat& G, rmat& H);

template<typename REAL>
TARGET double SPA_Likelihood(uint64 l, double* a);

template<typename REAL>
TARGET double SPA_Hessian(uint64 l, double* G, double* H, double* a, double* a2);

template<typename REAL>
TARGET void SPA_Locus(uint64 l, double* x, double* f, double* a, double* a2, double* am, double* g, double* h);

template<typename REAL>
TARGET void MatrixSPA(double* f, double* x, int spa_tn, int spa_np, double* a,
	int64 l, double* am, double* h, double* g, double* a2);

/* Calculate spatial pattern using multiple threads */
THREAD2H(SPAThread);
