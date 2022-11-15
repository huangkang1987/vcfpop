/* Menu */

#pragma once
#include "vcfpop.h"

/* Print help information */
TARGET bool PrintHelp();

/* Run benchmark for SIMD instruction set */
template<typename REAL>
TARGET void SimdBenchmark();
