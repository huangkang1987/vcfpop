/* Ploidy Inference Functions */

#pragma once
#include "vcfpop.h"

#pragma pack(push, 1)

#pragma pack(pop)

struct SPF
{
	//Used in ploidy infer
	int count;
	double val1[N_MAX_PLOIDY + 1];	//
	double val2[N_MAX_PLOIDY + 1];	//
};

/* Calculate individual ploidy inference */
TARGET void CalcPloidyInference();

/* Calculate ploidy inference using multiple threads */
THREADH(PloidyInferenceThread);
