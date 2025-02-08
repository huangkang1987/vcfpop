/* Method-of-moment estimation of relatedness coefficient for polyploids */

#pragma once
#include "vcfpop.h"

/* Power of a real number */
static forceinline TARGET double mp(double base, int index)
{
	double re = base;
	for(int i = 1; i < index; ++i)
		re *= base;
	return re;
}

/* Assign allele and frequency of polyploid method-of-moment estimator */
template<typename REAL>
TARGET void MOMRelatednessAssign(int p, int refmode, double* e, REAL* f, int* xx);

/* Assign allele and frequency of polyploid method-of-moment estimator */
template<typename REAL>
TARGET void MOMRelatednessAssign1(int refmode, double* e, REAL* f, int* xx);

/* Assign allele and frequency of polyploid method-of-moment estimator */
template<typename REAL>
TARGET void MOMRelatednessAssign2(int refmode, double* e, REAL* f, int* xx);

/* Assign allele and frequency of polyploid method-of-moment estimator */
template<typename REAL>
TARGET void MOMRelatednessAssign3(int refmode, double* e, REAL* f, int* xx);

/* Assign allele and frequency of polyploid method-of-moment estimator */
template<typename REAL>
TARGET void MOMRelatednessAssign4(int refmode, double* e, REAL* f, int* xx);

/* Assign allele and frequency of polyploid method-of-moment estimator */
template<typename REAL>
TARGET void MOMRelatednessAssign5(int refmode, double* e, REAL* f, int* xx);

/* Assign allele and frequency of polyploid method-of-moment estimator */
template<typename REAL>
TARGET void MOMRelatednessAssign6(int refmode, double* e, REAL* f, int* xx);

/* Assign allele and frequency of polyploid method-of-moment estimator */
template<typename REAL>
TARGET void MOMRelatednessAssign7(int refmode, double* e, REAL* f, int* xx);
