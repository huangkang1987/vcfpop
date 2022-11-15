/* Ploidy Inference Functions */

#pragma once
#include "vcfpop.h"

template<typename REAL> struct SPF;

#pragma pack(push, 1)

template<typename REAL>
struct SPF
{
	//Used in ploidy infer
	int count;
	REAL val1[N_MAX_PLOIDY + 1];	//
	REAL val2[N_MAX_PLOIDY + 1];	//
};

#pragma pack(pop)

/* Calculate individual ploidy inference */
template<typename REAL>
TARGET void CalcPloidyInference();

/* Calculate ploidy inference using multiple threads */
THREAD2H(PloidyInferenceThread);
