/* Hash Functions */

#include "vcfpop.h"

/* Get crc32 for 64 bits data*/
uint crc32_u64(uint prev, uint64 val);

/* Get crc32 value for a genotype */
TARGET HASH HashGenotype(ushort* allele, int ploidy);

/* Get crc32 value for a string */
TARGET HASH HashString(char* str, int len = -1);

/* Get crc32 value for unsigned long integer */
TARGET HASH HashULong(uint64 val);

/* Get crc32 value for unsigned long integer */
TARGET uint64 Hash64ULong(uint64 val);

/* Get crc32 value for two unsigned integer */
TARGET HASH HashUInt(uint val, uint val2 = 0);

/* Initialize crypt table for Huang 2015 maximum-likelihood polyploid relatedness estimator */
TARGET void InitCryptTable();

/* Get IBS model of a genotype */
TARGET int GetSingleIBS(int* x, int ploidy);

/* Get hash of IBS mode of a pair of genotypes */
TARGET uint HashString32(char* s1, char* s2, char* s3);

/* Get IBS mode of a pair of genotypes */
TARGET uint GetHuang2015Hash(int* x, int* y, int p);

/* Get hash of a pair of genotype id */
TARGET HASH HashDyadGenotypeIndex(HASH ha);

/* Get hash of a haplotype */
template<typename REAL>
TARGET void HashHaplotype(IND<REAL>* ti, int64 st, int64 ed, HASH* hash, int& ploidy);

/* 32 bit Integer Hashing */
TARGET uint MurmurHash2(uint data, uint seed);

/* 64 bit Integer Hashing */
TARGET uint64 MurmurHash64(uint64 data, uint64 seed);

/* 64 bit Integer Hashing */
TARGET uint MurmurHash32(uint64 data, uint64 seed);

/* Mix high and low 32 bits */
TARGET uint Mix(uint64 x);