/* Haplotype Extraction Functions */

#pragma once
#include "vcfpop.h"

struct QUICKSORT_PARAMETER;
struct HAPLO_DUMMY_HAPLOTYPE;
struct HAPLO_DUMMY_LOCUS;

#pragma pack(push, 1)

/* Quick sort locus */
struct QUICKSORT_PARAMETER
{
	int64 left;								//Left bound of quick sort range
	int64 right;							//Right bound of quick sort range
}; 

struct HAPLO_DUMMY_HAPLOTYPE
{
	ushort alleleid;						//Dummy allele index
	ushort* alleles;						//True allele array in [st-ed] variant, alloc in memory class

	/* Do nothing */
	TARGET HAPLO_DUMMY_HAPLOTYPE();

	/* Extract the ith haplotype from an individual */
	template<typename REAL>
	TARGET void ExtractHaplotype(int vi, IND<REAL>* ti, int64 st, int64 ed, int nvar, ushort aid, MEMORY& haplo_memory);

	/* Print information for an extracted locus */
	TARGET void PrintHaplotype(FILE* f1, int64 st, int64 ed);
};

struct HAPLO_DUMMY_LOCUS
{
	int64 chrid;							//Chrom id
	int64 stpos;							//Pos of start locus
	int64 st;								//Start locus id
	int64 ed;								//End locus id
	int nvar;								//Number of variants
	int hsize;								//Number of haplotypes
	int gsize;								//Number of genotypes
};

#pragma pack(pop)

extern atomic<int64> haplotype_contig;				//1st locus in the current contig 
extern LIST<QUICKSORT_PARAMETER> qslstack;			//Sort loci in haplotype extraction
extern LIST<HAPLO_DUMMY_LOCUS> haplotype_locus;		//Dummy locus information in haplotype extraction
extern SLOCUS* haplotype_nslocus;					//Extracted locus
extern uint64* locus_pos;							//Locus pos for small locus used in haplotype extraction
extern LOCN* locus_id;								//Locus id for small locus used in haplotype extraction

/* Check aneuploid in haplotype extraction */
THREAD2H(CheckAneuploid);

/* Perform haplotype extraction */
template<typename REAL>
TARGET void CalcHaplotype();

/* Quick sort locus by contig and position */
TARGET void QSLocus(int64 left, int64 right);

/* Quick sort locus in a contig */
THREADH(QSWorker);

/* Quick sort extracted locus by contig and position */
TARGET void QSHapLocus(int64 left, int64 right);

/* Quick sort extracted locus in a contig */
THREADH(QSHapWorker);

/* Get number of alleles and genotypes at a dummy locus */
template<typename REAL>
TARGET double GetDummyK(int64 st, int64 ed, TABLE<HASH, ushort>& hfidx, TABLE<HASH, ushort>& gfidx);

/* Create locus for haplotype extraction */
THREAD2H(CreateHaplotypeLocus);

/* Output locus for haplotype extraction */
THREAD2H(WriteHaplotypeLocus);