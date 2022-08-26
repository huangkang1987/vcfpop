/* Population Assignment Functions */

#pragma once
#include "vcfpop.h"

#pragma pack(push, 1)

#pragma pack(pop)

/* Calculate population assignment */
TARGET void CalcAssignment();

/* Calculate population assignment using multiple threads */
THREADH(PopulationAssignmentThread);