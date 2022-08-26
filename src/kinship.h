/* Kinship Functions */

#pragma once
#include "vcfpop.h"

#pragma pack(push, 1)

struct KINSHIP
{
	double Ritland1996;						//Estimates for various kinship estimators
	double Loiselle1995;
	double Weir1996;
	int ABtype;								//Number of loci genotyped in both individuals
	int Atype;								//Number of loci genotyped in individual A
	int Btype;								//Number of loci genotyped in individual B

	TARGET static void ColumnPrintHeader();

	TARGET void ColumnPrintLine(int i, int j);

	TARGET static void MatrixPrintMatrixHeader(int k, int n);

	TARGET static void MatrixPrintRowHeader(int k, int i);

	TARGET void MatrixPrintCell(int k);

	TARGET void CalcKinship(IND* a, IND* b);
};

#pragma pack(pop)

extern KINSHIP* kinship_buf;						//Circle buffer for kinship estimation, NBUF

/* Calculate kinship coefficient */
TARGET void CalcKinship();

/* Calculate relatedness coefficient using multiple threads */
THREADH(RelatednessThread);

/* Write column format kinship coefficient results in a guard thread */
THREADH(KinshipGuard1);

/* Write matrix format kinship coefficient results in a guard thread */
THREADH(KinshipGuard2);

/* Calculate kinship coefficient using multiple threads */
THREADH(KinshipThread);