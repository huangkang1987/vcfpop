/* AVX Instruction Set Functions */

#pragma once
#include "vcfpop.h"

#ifndef __aarch64__

/* Read allele id from decompressed bucket */
//TARGETAVX void ReadAidAVX(BAYESIAN_READER& rt, int size, __m256i* aid, __m256i* type);

/* Read allele frequency from decompressed bucket */
//TARGETAVX void ReadFreqAVX(BAYESIAN_READER& rt, int size, __m256i* slog, __m256d* prod, double* p, int K);

/* Update individual or allele origin when ancetral proportion is binary */
TARGETAVX void UpdateZBinaryAVX(int* Ni, ushort* Z);

/* Update a priori ancetral proportion for non-admix model */
TARGETAVX void UpdateQNoAdmixAVX(double* Q, int K, double* bufNK1, double* bufNK2, double* Gamma, RNGAVX& rng, bool locpriori, double* Base, ushort* Z);

/* Update a priori ancetral proportion by Metropolis-Hastings for admix model */
TARGETAVX void UpdateQMetroAVX(double* Q, int K, double* bufNK1, double* bufN1, double* bufN2, double* Base);

/* Update individual or allele origin for admix model */
TARGETAVX void UpdateZAdmixAVX(double* Q, int K, int64* Mi, int* Ni, double* bufNK1, double* bufNK2, double* Base, RNGAVX& rng);

#endif