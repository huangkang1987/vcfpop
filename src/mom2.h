/* Method-of-moment estimation of relatedness coefficient for polyploids */

#pragma once
#include "vcfpop.h"

/* Assign allele and frequency of polyploid method-of-moment estimator */

template<typename REAL>
TARGET void MOMRelatednessAssign8(int refmode, double* e, REAL* f, int* xx);