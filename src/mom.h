/* Method-of-moment estimation of relatedness coefficient for polyploids */

#pragma once
#include "vcfpop.h"

/* Power of a real number */
TARGET double mp(double base, int index);

/* Assign allele and frequency of polyploid method-of-moment estimator */
TARGET void MOMRelatednessAssign(int p, int refmode, double* e, double* f, int* xx);

/* Assign allele and frequency of polyploid method-of-moment estimator */
TARGET void MOMRelatednessAssign1(int refmode, double* e, double* f, int* xx);

/* Assign allele and frequency of polyploid method-of-moment estimator */
TARGET void MOMRelatednessAssign2(int refmode, double* e, double* f, int* xx);

/* Assign allele and frequency of polyploid method-of-moment estimator */
TARGET void MOMRelatednessAssign3(int refmode, double* e, double* f, int* xx);

/* Assign allele and frequency of polyploid method-of-moment estimator */
TARGET void MOMRelatednessAssign4(int refmode, double* e, double* f, int* xx);

/* Assign allele and frequency of polyploid method-of-moment estimator */
TARGET void MOMRelatednessAssign5(int refmode, double* e, double* f, int* xx);

/* Assign allele and frequency of polyploid method-of-moment estimator */
TARGET void MOMRelatednessAssign6(int refmode, double* e, double* f, int* xx);

/* Assign allele and frequency of polyploid method-of-moment estimator */
TARGET void MOMRelatednessAssign7(int refmode, double* e, double* f, int* xx);
