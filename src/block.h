/* Linkage Disequilibrium Decay Functions */

#pragma once
#include "vcfpop.h"

#pragma pack(push, 1)
template<typename REAL>
struct BLOCK_SINGLE
{
	char* chr;
	int id;
	int pos;
	int loc_st;
	int loc_ed;
};

template<typename REAL>
struct BLOCK_PAIR : DECAY_INTERVAL<REAL>
{
	BLOCK_SINGLE<REAL>* block1;
	BLOCK_SINGLE<REAL>* block2;

	/* Write table header */
	TARGET static void WriteHeader(FILE* fout);

	/* Write table entry */
	TARGET void Write(FILE* fout);
};

#pragma pack(pop)

template<typename REAL> extern BLOCK_SINGLE<REAL>*		block_singles;
template<typename REAL> extern BLOCK_PAIR<REAL>*		block_pairs;
extern vector<char*>									block_chromas;
extern umap<HASH, CHROM_PROP>							block_chrom_sted;

/* Calculate LD block */
template<typename REAL>
TARGET void CalcBlock();

/* Calculate LD block using multiple threads */
THREAD2H(BlockThread);