/* Genetic Distance Functions */

#pragma once
#include "vcfpop.h"

#pragma pack(push, 1)

template<typename REAL> struct INDGD;
template<typename REAL> struct GDIST;

template<typename REAL>
struct INDGD
{
	REAL ABtype;							//Number of loci both genotyped
	REAL Jx1, Jx2, Jxy;					//Nei1972; Nei1974;
	REAL Cavalli1967;						//Cavalli1967
	REAL t1, t2;							//Reynolds1983;
	REAL Nei1983;							//Nei1983
	REAL Euclidean;						//Euclidean
	REAL Goldstein1995;					//Goldstein1995
	REAL Roger1972;						//Roger1972
};

template<typename REAL>
struct GDIST
{
	double Nei1972;							//Genetic distance estimates
	double Cavalli1967;
	double Reynolds1983;
	double Nei1983;
	double Euclidean;
	double Goldstein1995;
	double Nei1974;
	double Roger1972;

	/* Converted from Fst by Slatkin's transform d = Fst/(1-Fst) */
	double Slatkin_Nei1973;					//Gst
	double Slatkin_Weir1984;				//ANOVA
	double Slatkin_Hudson1992;				//Mean allele difference
	double Slatkin_Slatkin1995;				//Allele size difference
	double Slatkin_Hedrick2005;				//G'st
	double Slatkin_Jost2008;				//D
	double Slatkin_Huang2021_homo;			//AMOVA homoploid
	double Slatkin_Huang2021_aneu;			//AMOVA aneuploid

	/* Converted from Fst by Reynolds's transform d = -ln(1 - Fst)*/
	double Reynolds_Nei1973;				//Gst
	double Reynolds_Weir1984;				//ANOVA
	double Reynolds_Hudson1992;				//Mean allele difference
	double Reynolds_Slatkin1995;			//Allele size difference
	double Reynolds_Hedrick2005;			//G'st
	double Reynolds_Jost2008;				//D
	double Reynolds_Huang2021_homo;			//AMOVA homoploid
	double Reynolds_Huang2021_aneu;		//AMOVA aneuploid

	/* Write column format header row for genetic distance estimation */
	TARGET static void ColumnFormatHeader();

	/* Write column format result row for genetic distance estimation */
	TARGET void ColumnFormatLine(int i, int j);

	/* Write matrix format header for genetic distance estimation */
	TARGET static void MatrixFormatHeader(int k, int n);

	/* Write matrix format row header for genetic distance estimation */
	TARGET static void MatrixFormatRowHeader(int k, int i);

	/* Write matrix format grid for genetic distance estimation */
	TARGET void MatrixFormatCell(int k);

	/* Use population/region allele frequency as the missing data */
	TARGET static void GetMissingFreq(GENOTYPE& gt, int64 l, REAL* p, int k);

	/* Calculate genetic distance between genotypes and save in gdtab */
	TARGET static void CacheIndGD();

	/* Calculate genetic distance between two individuals */
	TARGET void CalcGD(IND<REAL>* a, IND<REAL>* b, REAL* p1, REAL* p2);

	/* Calculate genetic distance between two populations/regions */
	TARGET void CalcGD(POP<REAL>* a, POP<REAL>* b, double* buf);
};

#pragma pack(pop)
extern int gdist_type;												//1 between inds, 2 between pops, 3 + between regions
extern int gdindex[N_GD_ESTIMATOR + 1];								//Index of ith used estimator
template<typename REAL>
extern GDIST<REAL>* gdist_buf;
template<typename REAL>
extern INDGD<REAL>** gd_tab;

/* Calculate genetic distance */
template<typename REAL>
TARGET void CalcDist();

/* Write column format genetic distance results in a guard thread */
THREAD2H(GeneticDistanceGuard1);

/* Write matrix format genetic distance results in a guard thread */
THREAD2H(GeneticDistanceGuard2);

/* Calculate genetic distance using multiple threads */
THREAD2H(GeneticDistanceThread);