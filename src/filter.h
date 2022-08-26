/* Filter Functions */

#pragma once
#include "vcfpop.h"

#pragma pack(push, 1)

#pragma pack(pop)

extern bool diversity_filter;
extern bool genotype_filter;
extern bool individual_filter;
extern bool info_filter;
extern TABLE<HASH, TABLE<HASH, pair<LOCN, double>>> contig_diversity; //Max diversity indices in a shift-window
extern LOCN nfilter;								//Number of filtered loci
extern atomic<int> nfilterind;						//Number of filtered individuals

/* Applying individual and diversity filters */
TARGET void ApplyFilter();

/* Recursive set id, pops, inds for each region */
TARGET void SetVReg(int rl, int i);

/* Sort individuals by population index to rinds array */
TARGET void AssignVInds();

/* Marker individual filtered or not */
THREADH(MarkerIndividual);

/* Marker locus filtered or not */
THREADH(MarkerDiversity);

/* Remove individual fail to pass filter */
THREADH(RemoveIndividual);

/* Remove locus fail to pass filter */
THREADH(RemoveLocus);