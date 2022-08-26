/* SSE Instruction Set Functions */

#pragma once
#include "vcfpop.h"

#ifndef __aarch64__

/* Read allele id from decompressed bucket */
//TARGETSSE void ReadAidSSE(BAYESIAN_READER& rt, int size, __m128i* aid, __m128i* type);

/* Read allele frequency from decompressed bucket */
//TARGETSSE void ReadFreqSSE(BAYESIAN_READER& rt, int size, __m128i* slog, __m128d* prod, double* p, int K);

/* Update individual or allele origin when ancetral proportion is binary */
TARGETSSE void UpdateZBinarySSE(int* Ni, ushort* Z);

/* Update a priori ancetral proportion for non-admix model */
TARGETSSE void UpdateQNoAdmixSSE(double* Q, int K, double* bufNK1, double* bufNK2, double* Gamma, RNGSSE& rng, bool locpriori, double* Base, ushort* Z);

/* Update a priori ancetral proportion by Metropolis-Hastings for admix model */
TARGETSSE void UpdateQMetroSSE(double* Q, int K, double* bufNK1, double* bufN1, double* bufN2, double* Base);

/* Update individual or allele origin for admix model */
TARGETSSE void UpdateZAdmixSSE(double* Q, int K, int64* Mi, int* Ni, double* bufNK1, double* bufNK2, double* Base, RNGSSE& rng);

#endif