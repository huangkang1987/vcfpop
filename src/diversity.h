/* Genetic Diversity Functions */

#pragma once
#include "vcfpop.h"

template<typename REAL> struct DIVERSITY;
template<typename REAL> struct DIVSUM;

#pragma pack(push, 1)

template<typename REAL>
struct DIVERSITY
{
	int64 l;								//Locus id

	/* Add */
	REAL bmaf;							//Minor allele freq of biallelic locus
	REAL ptype;							//Genotype rate
	REAL pval;							//P-val of genotype distribution test
	REAL he;								//Expected heterozygosity
	REAL ho;								//Observed heterozygosity
	REAL pic;								//Polymorphic information contents
	REAL ae;								//Effective number of alleles
	REAL I;								//Shannon¡¯s Information Index

	REAL fis;								//Inbreeding coefficient
	REAL g;								//G-statistic in HWE test
	REAL df;								//Degrees of freedom in HWE test

	int k;									//Number of alleles, int
	int n;									//Number of individuals, int
	int nhaplo;								//Number of allele copies, int
	int v2i;								//Number of within-allele pairs to weight Ho, int

	byte minploidy;							//Min ploidy
	byte maxploidy;							//Max ploidy
	bool varploidy;							//Does ploidy level varying among individuals
	bool unusued;							//Unusued byte for alignment

	/* Multiply */
	REAL NE1P;							//Exclusion rate without known parent
	REAL NE2P;							//Exclusion rate with known parent
	REAL NEPP;							//Exclusion rate for parent pair
	REAL NEID;							//Exclusion nonrelatives in identity test
	REAL NESI;							//Exclusion full-sibs in identity test

	/* Initialize */
	TARGET DIVERSITY();

	/* Calculate diveristy indices */
	TARGET void CalcDiversity(int64 _l);

	/* Write header to the result file */
	TARGET static void WriteHeader(FILE* f);

	/* Write diversity of a locus to the result file */
	TARGET void WriteLocus(FILE* f, const char* name);
};

template<typename REAL>
struct DIVSUM
{
	/* Add */
	REAL k;								//Number of alleles
	REAL n;								//Number of individuals
	REAL nhaplo;							//Number of allele copies
	REAL bmaf;							//Minor allele freq of biallelic locus
	REAL ptype;							//Genotype rate
	REAL pval;							//P-val of genotype distribution test
	REAL he;								//Expected heterozygosity
	REAL ho;								//Observed heterozygosity
	REAL pic;								//Polymorphic information contents
	REAL ae;								//Effective number of alleles
	REAL I;								//Shannon¡¯s Information Index
	REAL fis;								//Inbreeding coefficient
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
	REAL NE1P;							//Exclusion rate without known parent
	REAL NE2P;							//Exclusion rate with known parent
	REAL NEPP;							//Exclusion rate for parent pair
	REAL NEID;							//Exclusion nonrelatives in identity test
	REAL NESI;							//Exclusion full-sibs in identity test

	LOCK lock;
    
    /* Initialize sum */
    TARGET void Init();
    
    /* Uninitialize sum */
    TARGET void UnInit();

	/* Add loc to the sum */
	TARGET void Add(DIVERSITY<REAL>& loc);

	/* Write sum */
	TARGET void Write(FILE* f, const char* name);
};

#pragma pack(pop)

extern void* diversity_buf_;						//Circle buffer for diversity estimation, NBUF
#define diversity_buf (*(DIVERSITY<REAL>**)&diversity_buf_)

template<typename REAL>
extern DIVSUM<REAL> diversity_sum; 					//Diversity sum

extern int diversity_stage;							//Diversity level, 3 total, 2 pop, 1 reg

/* Add and write genetic diversity */
THREAD2H(DiversityGuard);

/* Calculate genetic diversity using multiple threads */
THREAD2H(DiversityThread);
