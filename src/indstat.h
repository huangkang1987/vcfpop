/* Individual Statistics Functions */

#pragma once
#include "vcfpop.h"

#pragma pack(push, 1)

#pragma pack(pop)

/* Calculate individual statistics */
TARGET void CalcIndstat();

/* Calculate individual statistics using multiple threads */
THREADH(IndividualStatisticsThread);