/* Genetic Diversity Functions */

#pragma once
#include "vcfpop.h"

#pragma pack(push, 1)

struct DIVERSITY
{
	int64 l;								//Locus id

	/* Add */
	double bmaf;							//Minor allele freq of biallelic locus
	double ptype;							//Genotype rate
	double pval;							//P-val of genotype distribution test
	double he;								//Expected heterozygosity
	double ho;								//Observed heterozygosity
	double pic;								//Polymorphic information contents
	double ae;								//Effective number of alleles
	double I;								//Shannon¡¯s Information Index

	double fis;								//Inbreeding coefficient
	double g;								//G-statistic in HWE test
	double df;								//Degrees of freedom in HWE test

	int k;									//Number of alleles, int
	int n;									//Number of individuals, int
	int nhaplo;								//Number of allele copies, int
	int v2i;								//Number of within-allele pairs to weight Ho, int

	byte minploidy;							//Min ploidy
	byte maxploidy;							//Max ploidy
	bool varploidy;							//Does ploidy level varying among individuals
	bool unusued;							//Unusued byte for alignment

	/* Multiply */
	double NE1P;							//Exclusion rate without known parent
	double NE2P;							//Exclusion rate with known parent
	double NEPP;							//Exclusion rate for parent pair
	double NEID;							//Exclusion nonrelatives in identity test
	double NESI;							//Exclusion full-sibs in identity test

	/* Initialize */
	TARGET DIVERSITY();

	/* Calculate diveristy indices */
	TARGET void CalcDiversity(int64 _l);

	/* Write header to the result file */
	TARGET static void WriteHeader(FILE* f);

	/* Write diversity of a locus to the result file */
	TARGET void WriteLocus(FILE* f, const char* name);
};

struct DIVSUM
{
	/* Add */
	double k;								//Number of alleles
	double n;								//Number of individuals
	double nhaplo;							//Number of allele copies
	double bmaf;							//Minor allele freq of biallelic locus
	double ptype;							//Genotype rate
	double pval;							//P-val of genotype distribution test
	double he;								//Expected heterozygosity
	double ho;								//Observed heterozygosity
	double pic;								//Polymorphic information contents
	double ae;								//Effective number of alleles
	double I;								//Shannon¡¯s Information Index
	double fis;								//Inbreeding coefficient
	int minploidy;							//Min ploidy
	int maxploidy;							//Max ploidy

	/* Count */
	int kc;									//Number of alleles
	int nc;									//Number of individuals
	int nhaploc;							//Number of allele copies
	int bmafc;								//Minor allele freq of biallelic locus
	int ptypec;								//Genotype rate
	int pvalc;								//P-val of genotype distribution test
	int hec;								//Expected heterozygosity
	int hoc;								//Observed heterozygosity
	int picc;								//Polymorphic information contents
	int aec;								//Effective number of alleles
	int Ic;									//Shannon¡¯s Information Index
	int fisc;								//Inbreeding coefficient

	/* Multiply */
	double NE1P;							//Exclusion rate without known parent
	double NE2P;							//Exclusion rate with known parent
	double NEPP;							//Exclusion rate for parent pair
	double NEID;							//Exclusion nonrelatives in identity test
	double NESI;							//Exclusion full-sibs in identity test

	LOCK lock;

	/* Initialize sum */
	TARGET void Init();

	/* Add loc to the sum */
	TARGET void Add(DIVERSITY& loc);

	/* Write sum */
	TARGET void Write(FILE* f, const char* name);
};

#pragma pack(pop)

extern DIVERSITY* diversity_buf;					//Circle buffer for diversity estimation, NBUF
extern DIVSUM diversity_sum;						//Diversity sum
extern int diversity_stage;							//Diversity level, 3 total, 2 pop, 1 reg

/* Add and write genetic diversity */
THREADH(DiversityGuard);

/* Calculate genetic diversity using multiple threads */
THREADH(DiversityThread);