/* Genetic Differentiation Functions */

#pragma once
#include "vcfpop.h"

template<typename REAL> struct FST;

#pragma pack(push, 1)

template<typename REAL>
struct FST
{
	/* Total value for all loci */
	double Nei1973T;							//Gst
	double Weir1984T;							//ANOVA
	double Hudson1992T;						//Mean allele difference
	double Slatkin1995T;						//Allele size difference
	double Hedrick2005T;						//G'st
	double Jost2008T;							//D
	double Huang2021homoT;					//AMOVA homoploid
	double Huang2021aneuT;					//AMOVA aneuploid

	/* Differentiation test by genotype distribution */
	double Genotype_GT;						//G-static
	int Genotype_DFT;						//degrees-of-freedom
	double Genotype_PT;						//P-val

	/* Differentiation test by allele distribution */
	double Allele_GT;							//G-static
	int Allele_DFT;							//Degrees-of-freedom
	double Allele_PT;							//P-val

	/* Value for each locus */
	double* Nei1973;							//Gst
	double* Weir1984;							//ANOVA
	double* Hudson1992;						//Mean allele difference
	double* Slatkin1995;						//Allele size difference
	double* Hedrick2005;						//G'st
	double* Jost2008;							//D
	double* Huang2021_homo;					//ANOVA homoploid
	double* Huang2021_aneu;					//ANOVA aneuploid

	/* Differentiation test by genotype distribution */
	double* Genotype_G;						//G-static
	int* Genotype_DF;						//degrees-of-freedom
	double* Genotype_P;						//P-val

	/* Differentiation test by allele distribution */
	double* Allele_G;							//G-static
	int* Allele_DF;							//Degrees-of-freedom
	double* Allele_P;							//P-val

	/* Fst estimator Warpper */
	TARGET static double FstEstimator(POP<REAL>** grps, int n, int e, double* each, double* buf);

	/* Estimate Fst and test differentiation for two pops */
	TARGET void CalcFst(POP<REAL>* a, POP<REAL>* b);

	/* Estimate Fst and test differentiation for multiple pops */
	TARGET void CalcFst(POP<REAL>** grps, int n);

	/* Uninitialize */
	TARGET void Uninitialize();

	/* Nei 1973 Fst estimator based on heterozgysotiy */
	TARGET static double Fst_Nei1973(POP<REAL>** grps, int n, double* each, double* buf, int64 _l = -1, double* frac1 = NULL, double* frac2 = NULL);

	/* Weir 1984 Fst estimator based on anova */
	TARGET static double Fst_Weir1984(POP<REAL>** grps, int n, double* each, int64 _l = -1, double* ffrac1 = NULL, double* frac2 = NULL);

	/* Hudson 1992 Fst estimator based on mean allele difference */
	TARGET static double Fst_Hudson1992(POP<REAL>** grps, int n, double* each, double* buf, int64 _l = -1, double* ffrac1 = NULL, double* frac2 = NULL);

	/* Slatkin 1995 Fst estimator based on allele size */
	TARGET static double Fst_Slatkin1995(POP<REAL>** grps, int n, double* each, double* buf, int64 _l = -1, double* ffrac1 = NULL, double* frac2 = NULL);

	/* Hedrick 2005 G'st */
	TARGET static double Fst_Hedrick2005(POP<REAL>** grps, int n, double* each, double* buf, int64 _l = -1, double* ffrac1 = NULL, double* frac2 = NULL);

	/* Jost 2008 D */
	TARGET static double Fst_Jost2008(POP<REAL>** grps, int n, double* each, double* buf, int64 _l = -1, double* ffrac1 = NULL, double* frac2 = NULL);

	/* Huang 2021 Fst estimator based on multi-level amova */
	TARGET static double Fst_Huang2021_homo(POP<REAL>** grps, int n, int layer, bool isiam, double* each);

	/* Huang 2021 Fst estimator based on multi-level amova */
	TARGET static double Fst_Huang2021_aneu(POP<REAL>** grps, int n, int layer, bool isiam, bool sumss, double* each, double* buf, int64 _l = -1, double* ffrac1 = NULL, double* frac2 = NULL);

	/* Hudson 1992 Fst estimator based on anova */
	TARGET static double dxy(POP<REAL>** grps, int n, double* buf, int64 l);

	/* Write results file in column format */
	TARGET static void ColumnFormat(FILE* fout);

	/* Write results file in matrix format */
	TARGET static void MatrixFormat(FILE* fout, FST* Fst, int n, int type);
};

#pragma pack(pop)

template<typename REAL>
extern FST<REAL>* fst_buf[6];							//Read/Write buffer for each type of fst
extern int fst_type;									//1 Among regions, 2 among pops, 3 among pops/regs in region, 4 between regions, 5 between pops 

/* Calculate genetic differentiation */
template<typename REAL>
TARGET void CalcDiff();

/* Calculate genetic differentiation using multiple threads */
THREAD2H(GeneticDifferentiationThread);