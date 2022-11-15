/* Kinship Functions */

#pragma once
#include "vcfpop.h"

template<typename REAL> struct KINSHIP;

#pragma pack(push, 1)

template<typename REAL>
struct KINSHIP
{
	REAL Ritland1996;						//Estimates for various kinship estimators
	REAL Loiselle1995;
	REAL Weir1996;
	int ABtype;								//Number of loci genotyped in both individuals
	int Atype;								//Number of loci genotyped in individual A
	int Btype;								//Number of loci genotyped in individual B

	TARGET static void ColumnPrintHeader();

	TARGET void ColumnPrintLine(int i, int j);

	TARGET static void MatrixPrintMatrixHeader(int k, int n);

	TARGET static void MatrixPrintRowHeader(int k, int i);

	TARGET void MatrixPrintCell(int k);

	TARGET void CalcKinship(IND<REAL>* a, IND<REAL>* b);
};

#pragma pack(pop)

extern void* kinship_buf_;						//Circle buffer for kinship estimation, NBUF
#define kinship_buf (*(KINSHIP<REAL>**)&kinship_buf_)

/* Calculate kinship coefficient */
template<typename REAL>
TARGET void CalcKinship();

/* Write column format kinship coefficient results in a guard thread */
THREAD2H(KinshipGuard1);

/* Write matrix format kinship coefficient results in a guard thread */
THREAD2H(KinshipGuard2);

/* Calculate kinship coefficient using multiple threads */
THREAD2H(KinshipThread);