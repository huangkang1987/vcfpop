/* Principal Coordinate Analysis Functions */

#pragma once
#include "vcfpop.h"

#pragma pack(push, 1)

struct PCOA
{
	double* D;								//Euclidean distance matrix
	int N;									//Number of objects
	double* U;								//Eigen-vector matrix
	double* V;								//Eigen-value matrix
	double Vt;								//Total variance
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
	TARGET void PrintPCoA(FILE* fout, double* d, int n, int _est, int _type);
};

#pragma pack(pop)

extern MEMORY* pcoa_memory;					//Genetic distance memory class
extern double* pcoa_matrix;					//Genetic distance array of genetic distance to perform PCoA

/* Calculate principal coordinate analysis */
TARGET void CalcPCOA();

/* Calculate genetic distance for PCoA using multiple threads */
THREADH(PCoAThread);