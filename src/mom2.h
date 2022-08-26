/* Method-of-moment estimation of relatedness coefficient for polyploids */

#pragma once
#include "vcfpop.h"

/* Assign allele and frequency of polyploid method-of-moment estimator */
TARGET void MOMRelatednessAssign8(int refmode, double* e, double* f, int* xx);