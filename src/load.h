/* Load File Functions */

#pragma once
#include "vcfpop.h"

#pragma pack(push, 1)

#pragma pack(pop)

extern int genotype_digit;
extern int genotype_extracol;
extern int genotype_missing;
extern int genotype_ambiguous;
extern bool load_complete;
extern bool usephase;
extern bool uselocpos;

extern INCBUFFER* load_buf;							//Circle buffer for loading vcf/bcf files, NBUF
extern char* vcf_header;							//Vcf header row

/* Check Genotype index is correct */
template<typename REAL>
TARGET void CheckGenotypeId();

/* Assign individual indid to population popid */
template<typename REAL>
TARGET void AssignIndSub(int indid, int popid, int& namelen, int& nind2);

/* Assign individuals to their populations */
template<typename REAL>
TARGET void AssignInds();

/* Load data from input files */
template<typename REAL>
TARGET void LoadFile();

/* Create genotype table for non-vcf/bcf input files */
TARGET void CreateGenoIndexTable(GENO_WRITER* wt = NULL);

/* load from Genepop input format */
THREAD2H(LoadGenepop);

/* load from Spagedi input format */
THREAD2H(LoadSpagedi);

/* load from Cervus input format */
THREAD2H(LoadCervus);

/* load from Arlequin input format */
THREAD2H(LoadArlequin);

/* load from Structure input format */
THREAD2H(LoadStructure);

/* load from PolyGene input format */
THREAD2H(LoadPolyGene);

/* load from PolyRelatedness input format */
THREAD2H(LoadPolyRelatedness);

/* load from PolyRelatedness input format */
THREAD2H(LoadGenoDive);

/* load from PolyRelatedness input format */
THREAD2H(LoadPlink);

/* Indexing alleles for non-vcf input, with allele identifier being the size */
TARGET void IndexAlleleLength();

/* Process lines from memory */
THREAD2H(LoadBCF);

/* Read lines from BCF file */
THREAD2H(LoadBCFGuard);

/* Process lines from memory */
THREAD2H(LoadVCF);

/* Read lines from VCF file */
THREAD2H(LoadVCFGuard);

/* Replace missing genotype with vmin */
TARGET void ReplaceMissingGenotypes();

/* Count ploidy level from VCF genotype string */
static forceinline TARGET int CountPloidy(char* str)
{
	if (str[0] == '\0') return 0;
	int count = 1;
	for (int i = 1; str[i]; ++i)
		count += (str[i] == '|') | (str[i] == '/');
	return count;
}

/* Count ploidy level from VCF genotype string */
static forceinline TARGET int CountPloidy(char* str, int len)
{
	int count = 1;
	switch (len)
	{
	case  0: return 0;
	default :
	{
		for (int i = 1; i < len - 1; ++i)
			count += (str[i] == '|') | (str[i] == '/');
		return count;
	}
	case 20: count += (str[18] == '|') | (str[18] == '/');
	case 19: count += (str[17] == '|') | (str[17] == '/');
	case 18: count += (str[16] == '|') | (str[16] == '/');
	case 17: count += (str[15] == '|') | (str[15] == '/');
	case 16: count += (str[14] == '|') | (str[14] == '/');
	case 15: count += (str[13] == '|') | (str[13] == '/');
	case 14: count += (str[12] == '|') | (str[12] == '/');
	case 13: count += (str[11] == '|') | (str[11] == '/');
	case 12: count += (str[10] == '|') | (str[10] == '/');
	case 11: count += (str[ 9] == '|') | (str[ 9] == '/');
	case 10: count += (str[ 8] == '|') | (str[ 8] == '/');
	case  9: count += (str[ 7] == '|') | (str[ 7] == '/');
	case  8: count += (str[ 6] == '|') | (str[ 6] == '/');
	case  7: count += (str[ 5] == '|') | (str[ 5] == '/');
	case  6: count += (str[ 4] == '|') | (str[ 4] == '/');
	case  5: count += (str[ 3] == '|') | (str[ 3] == '/');
	case  4: count += (str[ 2] == '|') | (str[ 2] == '/');
	case  3: count += (str[ 1] == '|') | (str[ 1] == '/');
	case  2: ;
	case  1: ;
	}
	return count;
}
