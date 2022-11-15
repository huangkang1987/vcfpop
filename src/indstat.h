/* Individual Statistics Functions */

#pragma once
#include "vcfpop.h"

#pragma pack(push, 1)

#pragma pack(pop)

/* Calculate individual statistics */
template<typename REAL>
TARGET void CalcIndstat();

/* Calculate individual statistics using multiple threads */
THREAD2H(IndividualStatisticsThread);