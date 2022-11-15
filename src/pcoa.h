/* Principal Coordinate Analysis Functions */

#pragma once
#include "vcfpop.h"

template<typename REAL> struct PCOA;

#pragma pack(push, 1)

template<typename REAL>
struct PCOA
{
	REAL* D;								//Euclidean distance matrix
	int N;									//Number of objects
	REAL* U;								//Eigen-vector matrix
	REAL* V;								//Eigen-value matrix
	REAL Vt;								//Total variance
	int estimator;							//PCoA or PCA
	int type;								//Type of objects: 1 individual, 2 pop, 3 reg
	int maxp;								//Number of extracted coordinates

	/* Do nothing */
	TARGET PCOA();

	/* Destructor */
	TARGET ~PCOA();

	/* Perform PCoA */
	TARGET int CalcPCoA(int _maxp);

	/* Print PCOA */
	TARGET void PrintPCoA(FILE* fout, REAL* d, int n, int _est, int _type);
};

#pragma pack(pop)

extern MEMORY* pcoa_memory;					//Genetic distance memory class
extern void* pcoa_matrix_;					//Genetic distance array of genetic distance to perform PCoA
#define pcoa_matrix (*(REAL**)&pcoa_matrix_)

/* Calculate principal coordinate analysis */
template<typename REAL>
TARGET void CalcPCOA();

/* Calculate genetic distance for PCoA using multiple threads */
THREAD2H(PCoAThread);