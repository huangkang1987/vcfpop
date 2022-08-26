/* Maximum-likelihood estimation of relatedness coefficient for polyploids */

#pragma once
#include "vcfpop.h"

/* Assign genotype pair pattern for polyploid maximum-likelihood relatedness estimation */
TARGET void MLRelatednessAssign(int sumploidy, double* f, int* a, double* c, int IBS);