/* Linkage Disequilibrium Decay Functions */

#pragma once
#include "vcfpop.h"

#pragma pack(push, 1)

template<typename REAL> struct DECAY_INTERVAL;
template<typename REAL> struct DECAY;

template<typename REAL>
struct DECAY_INTERVAL
{
	int npair;
	int unusued;

	//weighted average
	double r2A;  //sum r2 * B
	double r2A2; //sum r2^2 * B
	double r2B;  //sum B			V1
	double r2B2; //sum B^2			V2

	double r2DeltaA;
	double r2DeltaA2;
	double r2DeltaB;
	double r2DeltaB2;

	double DpA;
	double DpA2;
	double DpB;
	double DpB2;

	double DeltapA;
	double DeltapA2;
	double DeltapB;
	double DeltapB2;

	/* Set zero */
	TARGET DECAY_INTERVAL();

	/* Add thread specific interval to global specific interval */
	TARGET void AddInterval(DECAY_INTERVAL<REAL>& interval);

	/* Add DECAY measures to thread specific interval */
	TARGET void AddDecay(DECAY<REAL>& ld);

	/* Write table header */
	TARGET static void WriteHeader(FILE* fout);

	/* Write table entry */
	TARGET void Write(FILE* fout, char* chrom, int id);
};

template<typename REAL>
struct DECAY
{
	int loc1;
	int loc2;
	int nhaplo;

	double r2A; //numerator
	double r2B; //denominaotr

	double r2DeltaA;
	double r2DeltaB;

	double DpA;
	double DpB;

	double DeltapA;
	double DeltapB;

	/* Calculate LD measures between loc1 and several loc2 */
	TARGET void CalcLD(int l1, int* l2, int* buf, int n, TABLE<uint64, int>* genopair);

	/* Write LD measures header */
	TARGET static void WriteHeader(FILE* fout);

	/* Write LD measures to a temp file */
	TARGET void Write(int id);
};

#pragma pack(pop)

template<typename REAL> extern DECAY_INTERVAL<REAL>*	decay_global_intervals;
template<typename REAL> extern DECAY_INTERVAL<REAL>**	decay_chrom_intervals;
extern umap<HASH, int>									decay_chrom_id;
extern vector<char*>									decay_chomas;
extern double											decay_interval_width;

/* Calculate LD decay */
template<typename REAL>
TARGET void CalcDecay();

/* Calculate LD decay using multiple threads */
THREAD2H(DecayThread);

/* Calculate allele frequencies for cpop<REAL> */
THREAD2H(DecayFreqThread);