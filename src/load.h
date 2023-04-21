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
__forceinline TARGET int CountPloidy(char* str);

/* Count ploidy level from VCF genotype string */
__forceinline TARGET int CountPloidy(char* str, int len);