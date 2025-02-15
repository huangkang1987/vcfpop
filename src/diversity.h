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
	REAL bmaf;								//Minor allele freq of biallelic locus
	REAL ptype;								//Genotype rate
	REAL pval;								//P-val of genotype distribution test

	REAL ho;								//Total number of non-IBS of allele pairs (without replacement) within genotypes
	REAL he;								//Total number of non-IBS of allele pairs (with replacement) within populations
	REAL pi;								//Total number of non-IBS of allele pairs (without replacement) within populations
	REAL pic;								//Polymorphic information contents
	REAL ae;								//Effective number of alleles
	REAL I;									//Shannon's Information Index

	REAL how;								//Total number of allele pairs (without replacement) within genotypes
	REAL hew;								//Total number of allele pairs (with replacement) within populations
	REAL piw;								//Total number of allele pairs (without replacement) within populations
	REAL picw;								//Total number of allele pairs (with replacement) within populations
	REAL aew;								//Total number of allele pairs (with replacement) within populations
	REAL Iw;								//Total number of allele pairs (with replacement) within populations

	REAL g;									//G-statistic in HWE test
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
	REAL NE1P;								//Exclusion rate without known parent
	REAL NE2P;								//Exclusion rate with known parent
	REAL NEPP;								//Exclusion rate for parent pair
	REAL NEID;								//Exclusion nonrelatives in identity test
	REAL NESI;								//Exclusion full-sibs in identity test

	/* Initialize */
	TARGET DIVERSITY();

	/* Calculate diveristy indices */
	TARGET void CalcDiversity(int64 _l);

	/* Write header to the result file */
	TARGET static void WriteHeader(FILE* fout);

	/* Write diversity of a locus to the result file */
	TARGET void WriteLocus(FILE* fout, const char* name);
};

template<typename REAL>
struct DIVSUM
{
	/* Numerator */
	REAL k;									//Number of alleles
	REAL n;									//Number of individuals
	REAL nhaplo;							//Number of allele copies
	REAL bmaf;								//Minor allele freq of biallelic locus
	REAL ptype;								//Genotype rate
	REAL pval;								//P-val of genotype distribution test

	REAL ho;								//Expected heterozygosity
	REAL he;								//Observed heterozygosity
	REAL pi;								//Nucleotide diversity
	REAL pic;								//Polymorphic information contents
	REAL ae;								//Effective number of alleles
	REAL I;									//Shannon's Information Index

	int minploidy;							//Min ploidy
	int maxploidy;							//Max ploidy

	/* Denominator / Weight */
	int kc;									//Number of alleles
	int nc;									//Number of individuals
	int nhaploc;							//Number of allele copies
	int bmafc;								//Minor allele freq of biallelic locus
	int ptypec;								//Genotype rate
	int pvalc;								//P-val of genotype distribution test
	REAL how;								//Observed heterozygosity
	REAL hew;								//Expected heterozygosity
	REAL piw;								//Nucleotide diversity
	REAL picw;								//Polymorphic information contents
	REAL aew;								//Effective number of alleles
	REAL Iw;								//Shannon's Information Index

	/* Multiply */
	REAL NE1P;								//Exclusion rate without known parent
	REAL NE2P;								//Exclusion rate with known parent
	REAL NEPP;								//Exclusion rate for parent pair
	REAL NEID;								//Exclusion nonrelatives in identity test
	REAL NESI;								//Exclusion full-sibs in identity test

	LOCK lock;
    
    /* Initialize sum */
    TARGET void Init();
    
    /* Uninitialize sum */
    TARGET void UnInit();

	/* Add loc to the sum */
	TARGET void Add(DIVERSITY<REAL>& loc);

	/* Write sum */
	TARGET void Write(FILE* fout, const char* name);
};

#pragma pack(pop)

template<typename REAL>
extern DIVERSITY<REAL>* diversity_buf;				//Circle buffer for diversity estimation, NBUF
template<typename REAL>
extern DIVSUM<REAL> diversity_sum; 					//Diversity sum
extern int diversity_stage;							//Diversity level, 3 total, 2 pop, 1 reg

/* Add and write genetic diversity */
THREAD2H(DiversityGuard);

/* Calculate genetic diversity using multiple threads */
THREAD2H(DiversityThread);
