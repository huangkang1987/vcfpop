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
template<typename REAL>
TARGET void ConvertFile();

/* Write convert genepop genotypes in a guard thread */
THREAD2H(ConvertGenepopGuard);

/* Write convert arlequin genotypes in a guard thread */
THREAD2H(ConvertArlequinGuard);

/* Write convert genotypes in a guard thread */
THREAD2H(ConvertGuard);

/* Convert genotype string */
TARGET void PrepareGenotypeString(int format);

/* Convert into genepop format */
template<typename REAL>
TARGET void ConvertGenepop(int ntot, bool& isfirst);

/* Convert individual genotypes into genepop format in multiple threads */
THREAD2H(ConvertGenepopInd);

/* Convert into spagedi format */
template<typename REAL>
TARGET void ConvertSpagedi(int ntot, bool& isfirst);

/* Convert individual genotypes into spagedi format in multiple threads */
THREAD2H(ConvertSpagediInd);

/* Convert into cervus format */
template<typename REAL>
TARGET void ConvertCervus(int ntot, bool& isfirst);

/* Convert individual genotypes into cervus format in multiple threads */
THREAD2H(ConvertCervusInd);

/* Convert into arlequin format */
template<typename REAL>
TARGET void ConvertArlequin(int ntot, bool& isfirst);

/* Convert individual genotypes into arlequin format in multiple threads */
THREAD2H(ConvertArlequinInd);

/* Convert into structure format */
template<typename REAL>
TARGET void ConvertStructure(int ntot, bool& isfirst);

/* Convert individual genotypes into structure format in multiple threads */
THREAD2H(ConvertStructureInd);

/* Convert into polygene format */
template<typename REAL>
TARGET void ConvertPolygene(int ntot, bool& isfirst);

/* Convert individual genotypes into polygene format in multiple threads */
THREAD2H(ConvertPolygeneInd);

/* Convert into polyrelatedness format */
template<typename REAL>
TARGET void ConvertPolyRelatedness(int ntot, bool& isfirst);

/* Convert individual genotypes into polyrelatedness format in multiple threads */
THREAD2H(ConvertPolyRelatednessInd);

/* Convert into genodive format */
template<typename REAL>
TARGET void ConvertGenoDive(int ntot, bool& isfirst);

/* Convert individual genotypes into genodive format in multiple threads */
THREAD2H(ConvertGenoDiveInd);

/* Convert into plink format */
template<typename REAL>
TARGET void ConvertPlink(int ntot, bool& isfirst);

/* Convert individual genotypes into plink format in multiple threads */
THREAD2H(ConvertPlinkInd);