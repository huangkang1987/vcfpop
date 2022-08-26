/* Relatedness Functions */

#pragma once
#include "vcfpop.h"

#pragma pack(push, 1)

struct Huang2015ENTRY
{
	int ibs;								//IBS model index
	uint64 pattern;							//Genotype pair pattern
};

struct RELATEDNESS
{
	double Lynch1999;						//Estimates for various relatedness estimators
	double Wang2002;
	double Thomas2010;
	double Li1993;
	double Queller1989;
	double Huang2016A;
	double Huang2016B;
	double Milligan2003;
	double Anderson2007;
	double Huang2014;
	double Huang2015;
	double Ritland1996m;
	double Loiselle1995m;
	double Ritland1996;
	double Loiselle1995;
	double Weir1996;
	int ABtype;								//Number of loci genotyped in both individuals
	int Atype;								//Number of loci genotyped in individual A
	int Btype;								//Number of loci genotyped in individual B

	/* Write header row for relatedness estimation */
	TARGET static void ColumnPrintHeader();

	/* Write result row for relatedness estimation */
	TARGET void ColumnPrintLine(int i, int j);

	/* Write matrix format header for relatedness estimation */
	TARGET static void MatrixPrintMatrixHeader(int k, int n);

	/* Write matrix format row header for relatedness estimation */
	TARGET static void MatrixPrintRowHeader(int k, int i);

	/* Write matrix format grid for relatedness estimation */
	TARGET void MatrixPrintCell(int k);

	/* Calculate relatedness coefficient */
	TARGET void CalcRelatedness(IND* x, IND* y);

	/* Relatedness estimator Warpper */
	TARGET static double RelatednessEstimator(int k, IND* x, IND* y);

	/* Lynch 1999 relatedness estimator */
	TARGET static double R_Lynch1999(IND* x, IND* y);

	/* Wang 2002 relatedness estimator */
	TARGET static double R_Wang2002(IND* x, IND* y);

	/* Thomas 2010 relatedness estimator */
	TARGET static double R_Thomas2010(IND* x, IND* y);

	/* Li 1993 relatedness estimator */
	TARGET static double R_Li1993(IND* x, IND* y);

	/* Queller 1989 relatedness estimator */
	TARGET static double R_Queller1989(IND* x, IND* y);

	/* Huang 2016 relatedness estimator A */
	TARGET static double R_Huang2016A(IND* x, IND* y);

	/* Huang 2016 relatedness estimator B */
	TARGET static double R_Huang2016B(IND* x, IND* y);

	/* Initialize Anderson 2007 relatedness estimator */
	TARGET static void R_AndersonInitialize(IND* x, IND* y);

	/* Milligan 2003 relatedness estimator */
	TARGET static double R_Milligan2003(IND* x, IND* y);

	/* Anderson 2007 relatedness estimator */
	TARGET static double R_Anderson2007(IND* x, IND* y, bool confine);

	/* Calculate Anderson 2007 likelihood */
	TARGET static double L_Anderson(CPOINT& x, void** unusued);

	/* Ritland 1996 kinship estimator, convert into relatedness */
	TARGET static double R_Ritland1996(IND* x, IND* y, bool iscorrect, bool mulv);

	/* Loiselle 1995 kinship estimator, convert into relatedness */
	TARGET static double R_Loiselle1995(IND* x, IND* y, bool iscorrect, bool mulv);

	/* Weir 1996 kinship estimator, convert into relatedness */
	TARGET static double R_Weir1996(IND* x, IND* y, bool mulv);

	/* Huang 2014 relatedness estimator */
	TARGET static double R_Huang2014(IND* x, IND* y);

	/* Huang 2014 relatedness estimator : similarity index */
	TARGET static double S_Index(int* c, int* d, int ploidyx, int ploidyy);

	/* Huang 2014 relatedness estimator : get genotype pattern for reference individual */
	TARGET static int GetRefMode(int* a, int ploidy);

	/* Huang 2014 relatedness estimator : calculate relatedness */
	TARGET static double HuangMoment(GENOTYPE& gx, GENOTYPE& gy, int64 l, double& weight);

	/* Initialize Huang 2015 relatedness estimator */
	TARGET static void Huang2015_Initialize();

	/* Uninitialize Huang 2015 relatedness estimator */
	TARGET static void Huang2015_Uninitialize();

	/* Huang 2015 relatedness estimator */
	TARGET static double R_Huang2015(IND* x, IND* y);

	/* Calculate Huang 2015 likelihood */
	TARGET static double L_Huang2015(CPOINT& xx, void** unusued);

	/* Huang 2015 likelihood estimator: Match genotype-pair pattern and assign alleles */
	TARGET static void Huang2015_MatchAllele(int64 pattern, int* gx, int* gy, int* alleles, int p);
};

#pragma pack(pop)

extern _thread double* Anderson2007_Coef;			//Anderson 2007 relatedness estimator coefficients
extern _thread double* Huang2015_Coef;				//Huang 2015 relatedness estimator coefficients
extern TABLE<int, Huang2015ENTRY>* Huang2015_maps;	//Huang2015_maps[ploidylevel][hash] is a entry saves the ibs modex index and genotype pair pattern
extern RELATEDNESS* relatedness_buf;				//Circle buffer for relatedness estimation, NBUF
extern LOCSTAT2* relatedness_loc_stat2;				//Locus information temporatorily used for cpop

/* Calculate relatedness coefficient */
TARGET void CalcRelatedness();

/* Write column format relatedness coefficient results in a guard thread */
THREADH(RelatednessGuard1);

/* Write matrix format relatedness coefficient results in a guard thread */
THREADH(RelatednessGuard2);

/* Calculate relatedness coefficient using multiple threads */
THREADH(RelatednessThread);