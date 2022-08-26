/* NEO Instruction Set Functions */

#pragma once
#include "vcfpop.h"

#ifdef __aarch64__

/* Read allele id from decompreNEOd bucket */
//TARGETNEO void ReadAidNEO(BAYESIAN_READER& rt, int size, uint64x2_t* aid, uint64x2_t* type);

/* Read allele frequency from decompreNEOd bucket */
//TARGETNEO void ReadFreqNEO(BAYESIAN_READER& rt, int size, uint64x2_t* slog, float64x2_t* prod, double* p, int K);

/* Update individual or allele origin when ancetral proportion is binary */
TARGETNEO void UpdateZBinaryNEO(int* Ni, ushort* Z);

/* Update a priori ancetral proportion for non-admix model */
TARGETNEO void UpdateQNoAdmixNEO(double* Q, int K, double* bufNK1, double* bufNK2, double* Gamma, RNGNEO& rng, bool locpriori, double* Base, ushort* Z);

/* Update a priori ancetral proportion by Metropolis-Hastings for admix model */
TARGETNEO void UpdateQMetroNEO(double* Q, int K, double* bufNK1, double* bufN1, double* bufN2, double* Base);

/* Update individual or allele origin for admix model */
TARGETNEO void UpdateZAdmixNEO(double* Q, int K, int64* Mi, int* Ni, double* bufNK1, double* bufNK2, double* Base, RNGNEO& rng);

#endif