/* Maximum-likelihood estimation of relatedness coefficient for polyploids */

#pragma once
#include "vcfpop.h"

/* Assign genotype pair pattern for polyploid maximum-likelihood relatedness estimation */
template<typename REAL>
TARGET void MLRelatednessAssign2(int sumploidy, REAL* f, int* a, double* c, int IBS);