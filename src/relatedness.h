/* Relatedness Functions */

#pragma once
#include "vcfpop.h"

struct Huang2015ENTRY;
template<typename REAL> struct RELATEDNESS;

#pragma pack(push, 1)

struct Huang2015ENTRY
{
	int ibs;								//IBS model index
	uint64 pattern;							//Genotype pair pattern
};

template<typename REAL>
struct RELATEDNESS
{
	REAL Lynch1999;						//Estimates for various relatedness estimators
	REAL Wang2002;
	REAL Thomas2010;
	REAL Li1993;
	REAL Queller1989;
	REAL Huang2016A;
	REAL Huang2016B;
	REAL Milligan2003;
	REAL Anderson2007;
	REAL Huang2014;
	REAL Huang2015;
	REAL Ritland1996m;
	REAL Loiselle1995m;
	REAL Ritland1996;
	REAL Loiselle1995;
	REAL Weir1996;
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
	TARGET void CalcRelatedness(IND<REAL>* x, IND<REAL>* y);

	/* Relatedness estimator Warpper */
	TARGET static double RelatednessEstimator(int k, IND<REAL>* x, IND<REAL>* y);

	/* Lynch 1999 relatedness estimator */
	TARGET static double R_Lynch1999(IND<REAL>* x, IND<REAL>* y);

	/* Wang 2002 relatedness estimator */
	TARGET static double R_Wang2002(IND<REAL>* x, IND<REAL>* y);

	/* Thomas 2010 relatedness estimator */
	TARGET static double R_Thomas2010(IND<REAL>* x, IND<REAL>* y);

	/* Li 1993 relatedness estimator */
	TARGET static double R_Li1993(IND<REAL>* x, IND<REAL>* y);

	/* Queller 1989 relatedness estimator */
	TARGET static double R_Queller1989(IND<REAL>* x, IND<REAL>* y);

	/* Huang 2016 relatedness estimator A */
	TARGET static double R_Huang2016A(IND<REAL>* x, IND<REAL>* y);

	/* Huang 2016 relatedness estimator B */
	TARGET static double R_Huang2016B(IND<REAL>* x, IND<REAL>* y);

	/* Initialize Anderson 2007 relatedness estimator */
	TARGET static void R_AndersonInitialize(IND<REAL>* x, IND<REAL>* y);

	/* Milligan 2003 relatedness estimator */
	TARGET static double R_Milligan2003(IND<REAL>* x, IND<REAL>* y);

	/* Anderson 2007 relatedness estimator */
	TARGET static double R_Anderson2007(IND<REAL>* x, IND<REAL>* y, bool confine);

	/* Calculate Anderson 2007 likelihood */
	TARGET static double L_Anderson(CPOINT& x, void** unusued);

	/* Ritland 1996 kinship estimator, convert into relatedness */
	TARGET static double R_Ritland1996(IND<REAL>* x, IND<REAL>* y, bool iscorrect, bool mulv);

	/* Loiselle 1995 kinship estimator, convert into relatedness */
	TARGET static double R_Loiselle1995(IND<REAL>* x, IND<REAL>* y, bool iscorrect, bool mulv);

	/* Weir 1996 kinship estimator, convert into relatedness */
	TARGET static double R_Weir1996(IND<REAL>* x, IND<REAL>* y, bool mulv);

	/* Huang 2014 relatedness estimator */
	TARGET static double R_Huang2014(IND<REAL>* x, IND<REAL>* y);

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
	TARGET static double R_Huang2015(IND<REAL>* x, IND<REAL>* y);

	/* Calculate Huang 2015 likelihood */
	TARGET static double L_Huang2015(CPOINT& xx, void** unusued);

	/* Huang 2015 likelihood estimator: Match genotype-pair pattern and assign alleles */
	TARGET static void Huang2015_MatchAllele(int64 pattern, int* gx, int* gy, int* alleles, int p);
};

#pragma pack(pop)

extern TABLE<int, Huang2015ENTRY>* Huang2015_maps;				//Huang2015_maps[ploidylevel][hash] is a entry saves the ibs modex index and genotype pair pattern
extern _thread double* Anderson2007_Coef;						//Anderson 2007 relatedness estimator coefficients
extern _thread double* Huang2015_Coef;							//Huang 2015 relatedness estimator coefficients
extern void* relatedness_buf_;									//Circle buffer for relatedness estimation, NBUF
#define relatedness_buf (*(RELATEDNESS<REAL>**)&relatedness_buf_)
extern void* relatedness_loc_stat_;							//Locus information temporatorily used for cpop
#define relatedness_loc_stat (*(LOCSTAT2<REAL>**)&relatedness_loc_stat_)

/* Calculate relatedness coefficient */
template<typename REAL>
TARGET void CalcRelatedness();

/* Write column format relatedness coefficient results in a guard thread */
THREAD2H(RelatednessGuard1);

/* Write matrix format relatedness coefficient results in a guard thread */
THREAD2H(RelatednessGuard2);

/* Calculate relatedness coefficient using multiple threads */
THREAD2H(RelatednessThread);