/* AVX512 Instruction Set Functions */

#pragma once
#include "vcfpop.h"

#ifndef __aarch64__

/* Read allele id from decompressed bucket */
//TARGET512 void ReadAid512(BAYESIAN_READER& rt, int size, __m512i& aid, byte& typeflag);

/* Read allele frequency from decompressed bucket */
//TARGET512 void ReadFreq512(BAYESIAN_READER& rt, int size, __m512i* slog, __m512d* prod, double* p, int K);

/* Update individual or allele origin when ancetral proportion is binary */
TARGET512 void UpdateZBinary512(int* Ni, ushort* Z);

/* Update a priori ancetral proportion for non-admix model */
TARGET512 void UpdateQNoAdmix512(double* Q, int K, double* bufNK1, double* bufNK2, double* Gamma, RNG512& rng, bool locpriori, double* Base, ushort* Z);

/* Update a priori ancetral proportion by Metropolis-Hastings for admix model */
TARGET512 void UpdateQMetro512(double* Q, int K, double* bufNK1, double* bufN1, double* bufN2, double* Base);

/* Update individual or allele origin for admix model */
TARGET512 void UpdateZAdmix512(double* Q, int K, int64* Mi, int* Ni, double* bufNK1, double* bufNK2, double* Base, RNG512& rng);

#endif