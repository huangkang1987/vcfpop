/* File Conversion Functions */

#pragma once
#include "vcfpop.h"

#pragma pack(push, 1)

#pragma pack(pop)

extern MEMORY* conversion_memory;					//Memory class for conversion_string
extern MEMORY* conversion_memory2;					//Memory class for genotype_string and convert_buf
extern LIST<char*>* conversion_string;				//Genotype string for each genotype for file conversion
extern FILE* convert_file;							//Convert file pointer
extern char** convert_buf;							//Circle buffer for file conversion, NBUF
extern int64 convert_linesize;						//Max length of each line in converted file

/* Convert into other genotype formats */
TARGET void ConvertFile();

/* Write convert genepop genotypes in a guard thread */
THREADH(ConvertGenepopGuard);

/* Write convert arlequin genotypes in a guard thread */
THREADH(ConvertArlequinGuard);

/* Write convert genotypes in a guard thread */
THREADH(ConvertGuard);

/* Convert genotype string */
TARGET void PrepareGenotypeString(int format);

/* Convert into genepop format */
TARGET void ConvertGenepop(int ntot, bool& isfirst);

/* Convert individual genotypes into genepop format in multiple threads */
THREADH(ConvertGenepopInd);

/* Convert into spagedi format */
TARGET void ConvertSpagedi(int ntot, bool& isfirst);

/* Convert individual genotypes into spagedi format in multiple threads */
THREADH(ConvertSpagediInd);

/* Convert into cervus format */
TARGET void ConvertCervus(int ntot, bool& isfirst);

/* Convert individual genotypes into cervus format in multiple threads */
THREADH(ConvertCervusInd);

/* Convert into arlequin format */
TARGET void ConvertArlequin(int ntot, bool& isfirst);

/* Convert individual genotypes into arlequin format in multiple threads */
THREADH(ConvertArlequinInd);

/* Convert into structure format */
TARGET void ConvertStructure(int ntot, bool& isfirst);

/* Convert individual genotypes into structure format in multiple threads */
THREADH(ConvertStructureInd);

/* Convert into polygene format */
TARGET void ConvertPolygene(int ntot, bool& isfirst);

/* Convert individual genotypes into polygene format in multiple threads */
THREADH(ConvertPolygeneInd);

/* Convert into polyrelatedness format */
TARGET void ConvertPolyRelatedness(int ntot, bool& isfirst);

/* Convert individual genotypes into polyrelatedness format in multiple threads */
THREADH(ConvertPolyRelatednessInd);

/* Convert into genodive format */
TARGET void ConvertGenoDive(int ntot, bool& isfirst);

/* Convert individual genotypes into genodive format in multiple threads */
THREADH(ConvertGenoDiveInd);