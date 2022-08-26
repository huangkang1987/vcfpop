/* Hi/* Load File Functions */

#pragma once
#include "vcfpop.h"

#pragma pack(push, 1)

#pragma pack(pop)

extern int genotype_digit;
extern int genotype_extracol;
extern int genotype_missing;
extern int genotype_ambiguous;

extern INCBUFFER* load_buf;							//Circle buffer for loading vcf/bcf files, NBUF
extern char* vcf_header;							//Vcf header row

/* Check Genotype index is correct */
TARGET void CheckGenotypeId();

/* Assign individual indid to population popid */
TARGET void AssignIndSub(int indid, int popid, int& namelen, int& nind2);

/* Assign individuals to their populations */
TARGET void AssignInds();

/* Load data from input files */
TARGET void LoadFile();

/* Create genotype table for non-vcf/bcf input files */
TARGET void CreateGenoIndexTable(GENO_WRITER* wt = NULL);

/* load from Genepop input format */
THREADH(LoadGenepop);

/* load from Spagedi input format */
THREADH(LoadSpagedi);

/* load from Cervus input format */
THREADH(LoadCervus);

/* load from Arlequin input format */
THREADH(LoadArlequin);

/* load from Structure input format */
THREADH(LoadStructure);

/* load from PolyGene input format */
THREADH(LoadPolyGene);

/* load from PolyRelatedness input format */
THREADH(LoadPolyRelatedness);

/* load from PolyRelatedness input format */
THREADH(LoadGenoDive);

/* Indexing alleles for non-vcf input, with allele identifier being the size */
TARGET void IndexAlleleLength();

/* Process lines from memory */
THREADH(LoadBCF);

/* Read lines from BCF file */
THREADH(LoadBCFGuard);

/* Process lines from memory */
THREADH(LoadVCF);

/* Read lines from VCF file */
THREADH(LoadVCFGuard);

/* Replace missing genotype with vmin */
TARGET void ReplaceMissingGenotypes();