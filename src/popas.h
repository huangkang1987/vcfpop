/* Population Assignment Functions */

#pragma once
#include "vcfpop.h"

#pragma pack(push, 1)

#pragma pack(pop)

/* Calculate population assignment */
template<typename REAL>
TARGET void CalcAssignment();

/* Calculate population assignment using multiple threads */
THREAD2H(PopulationAssignmentThread);