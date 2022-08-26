/* AVX Instruction Set Functions */

#include "vcfpop.h"

#ifndef __aarch64__

/* Read allele id from decompressed bucket */
//TARGETAVX void ReadAidAVX(BAYESIAN_READER& rt, int size, __m256i* aid, __m256i* type)
#define ReadAidAVX \
{ \
	/* if data is empty */ \
	if (rt.nbits < size) [[unlikely]] \
	{ \
		/* move nbits data to gid */ \
		aid1 = rt.data256[0]; \
		aid2 = rt.data256[1]; \
 \
		/* remain number of bits to read */ \
		int rbits = size - rt.nbits; \
 \
		/* read 64 bits to data */ \
		rt.data256[0] = _mm256_i64gather_epi64((int64*)rt.pos, vindex, sizeof(int64)); \
		rt.data256[1] = _mm256_i64gather_epi64((int64*)rt.pos + structure_indnbytes * 4, vindex, sizeof(int64)); \
 \
		rt.pos++; \
 \
		/* read rbits from data and concate to higher bits in gid */ \
		aid1 = _mm256_or_si256(aid1, _mm256_slli_epi64(_mm256_and_si256(rt.data256[0], _mm256_set1_epi64x((1u << rbits) - 1u)), rt.nbits)); \
		aid2 = _mm256_or_si256(aid2, _mm256_slli_epi64(_mm256_and_si256(rt.data256[1], _mm256_set1_epi64x((1u << rbits) - 1u)), rt.nbits)); \
 \
		/*shift right */ \
		rt.data256[0] = _mm256_srli_epi64(rt.data256[0], rbits); \
		rt.data256[1] = _mm256_srli_epi64(rt.data256[1], rbits); \
 \
		rt.nbits = 64 - rbits; \
	} \
	else [[likely]] \
	{ \
		/*read size bits */ \
		aid1 = _mm256_and_si256(rt.data256[0], mask); \
		aid2 = _mm256_and_si256(rt.data256[1], mask); \
 \
		/*shift right */ \
		rt.data256[0] = _mm256_srli_epi64(rt.data256[0], size); \
		rt.data256[1] = _mm256_srli_epi64(rt.data256[1], size); \
 \
		rt.nbits -= size; \
	} \
 \
	/*typed is 0xFFFFFFFF */ \
	type1 = _mm256_cmpgt_epi64(mask, aid1); \
	type2 = _mm256_cmpgt_epi64(mask, aid2); \
 \
	/* set missing to 0 */ \
	aid1 = _mm256_and_si256(aid1, type1); \
	aid2 = _mm256_and_si256(aid2, type2); \
}

/* Read allele frequency from decompressed bucket */
/* TARGETAVX void ReadFreqAVX(BAYESIAN_READER& rt, int size, __m256i* slog, __m256d* prod, double* p, int K)
{
	__m256i mask1 = _mm256_set1_epi64x(0x7FF0000000000000);
	__m256i mask2 = _mm256_set1_epi64x(0x800FFFFFFFFFFFFF);
	__m256i mask3 = _mm256_set1_epi64x(0x3FF0000000000000);
	__m256i subv = _mm256_set1_epi64x(1023);
	__m256i maskff = _mm256_set1_epi64x(0xFFFFFFFFFFFFFFFF);
	__m256i maskunder = _mm256_set1_epi64x(0x1FF0000000000000); 
	__m256d maskone = _mm256_set1_pd(1.0);
	__m256i vindex = _mm256_set_epi64x(3 * structure_indnbytes, 2 * structure_indnbytes, structure_indnbytes, 0); */
#define ReadFreqAVX \
{ \
	__m256i aid1, aid2; \
	__m256i type1, type2; \
	__m256i a1, a2; \
	__m256i b1, b2; \
	__m256d freq1, freq2; \
	__m256i* prodi = (__m256i*)prod; \
 \
	ReadAidAVX \
 \
	for (int k = 0; k < K; ++k) \
	{ \
		/* get freq of 8 allele copies for 8 individuals in cluster k */ \
		freq1 = _mm256_mask_i64gather_pd(maskone, p + k * KT, aid1, _mm256_castsi256_pd(type1), sizeof(double)); \
		freq2 = _mm256_mask_i64gather_pd(maskone, p + k * KT, aid2, _mm256_castsi256_pd(type2), sizeof(double)); \
 \
		/* mul to prod */ \
		prod[k * 2 + 0] = _mm256_mul_pd(prod[k * 2 + 0], freq1); \
		prod[k * 2 + 1] = _mm256_mul_pd(prod[k * 2 + 1], freq2); \
\
		/* check underflow */ \
		a1 = _mm256_cmpgt_epi64(maskunder, _mm256_castpd_si256(prod[k * 2 + 0])); \
		a2 = _mm256_cmpgt_epi64(maskunder, _mm256_castpd_si256(prod[k * 2 + 1])); \
 \
		a1 = _mm256_or_si256(a1, a2); \
 \
		if (_mm256_testz_si256(a1, maskff)) [[likely]] continue; \
 \
		/* add exponent */ \
 \
		a1 = _mm256_and_si256(prodi[k * 4 + 0], mask1); \
		a2 = _mm256_and_si256(prodi[k * 4 + 1], mask1); \
 \
		b1 = _mm256_and_si256(prodi[k * 4 + 0], mask2); \
		b2 = _mm256_and_si256(prodi[k * 4 + 1], mask2); \
 \
		a1 = _mm256_srli_epi64(a1, 52); \
		a2 = _mm256_srli_epi64(a2, 52); \
 \
		prodi[k * 4 + 0] = _mm256_or_si256(b1, mask3); \
		prodi[k * 4 + 1] = _mm256_or_si256(b2, mask3); \
 \
		a1 = _mm256_sub_epi64(a1, subv); \
		a2 = _mm256_sub_epi64(a2, subv); \
 \
		slog[k * 4 + 0] = _mm256_add_epi64(slog[k * 4 + 0], a1); \
		slog[k * 4 + 1] = _mm256_add_epi64(slog[k * 4 + 1], a2); \
	} \
}

/* Update individual or allele origin when ancetral proportion is binary */
TARGETAVX void UpdateZBinaryAVX(int* Ni, ushort* Z)
{
	__m256i Kt = _mm256_set1_epi64x(KT);
	__m256i mask01 = _mm256_set1_epi64x(1);
	__m256i vindex = _mm256_set_epi64x(3 * structure_indnbytes, 2 * structure_indnbytes, structure_indnbytes, 0); 

	for (int i = 0, pad = 0; i < nind; i += STRUCTURE_NPACK)
	{
		if (i + STRUCTURE_NPACK >= nind) //fix i for the last several
		{
			pad = STRUCTURE_NPACK - (nind - i);
			i = nind - STRUCTURE_NPACK;
		}

		__m256i ZxKT1, ZxKT2;
		__m128i t1, t2;

		t1 = _mm_loadu_si128((__m128i*) & Z[i + 0]);
		t2 = _mm_loadu_si128((__m128i*) & Z[i + 4]);

		ZxKT1 = _mm256_cvtepi16_epi64(t1);
		ZxKT2 = _mm256_cvtepi16_epi64(t2);

		ZxKT1 = _mm256_mul_epu32(Kt, ZxKT1);
		ZxKT2 = _mm256_mul_epu32(Kt, ZxKT2);

		ZxKT1 = _mm256_slli_epi64(ZxKT1, 2);
		ZxKT2 = _mm256_slli_epi64(ZxKT2, 2);

		BAYESIAN_READER rt(i);
		int* ni = Ni;

		for (int64 l = 0; l < nloc; ++l, ni += GetLoc(l).k)
		{
			int size = structure_size[l];
			__m256i mask = _mm256_set1_epi64x((1u << size) - 1u); 
			__m256i aid1, aid2;
			__m256i type1, type2;
			__m256i npio2 = _mm256_set1_epi64x((int64)ni);
			__m256i a1, a2;
			__m256i b1, b2;

			a1 = _mm256_add_epi64(ZxKT1, npio2);
			a2 = _mm256_add_epi64(ZxKT2, npio2);

#define REP_MID \
			ReadAidAVX /*(rt, size, aid, type);*/ \
			\
			b1 = _mm256_slli_epi64(aid1, 2); \
			b2 = _mm256_slli_epi64(aid2, 2); \
			\
			b1 = _mm256_add_epi64(a1, b1); \
			b2 = _mm256_add_epi64(a2, b2); \
			\
			type1 = _mm256_and_si256(type1, mask01);\
			type2 = _mm256_and_si256(type2, mask01);\
			\
			switch (pad) \
			{ /* Z and Mi is updated before call, update Ni here */ \
			case 0: (*(int*)simd_u64(b1, 0)) += simd_u64(type1, 0); \
			case 1: (*(int*)simd_u64(b1, 1)) += simd_u64(type1, 1); \
			case 2: (*(int*)simd_u64(b1, 2)) += simd_u64(type1, 2); \
			case 3: (*(int*)simd_u64(b1, 3)) += simd_u64(type1, 3); \
			case 4: (*(int*)simd_u64(b2, 0)) += simd_u64(type2, 0); \
			case 5: (*(int*)simd_u64(b2, 1)) += simd_u64(type2, 1); \
			case 6: (*(int*)simd_u64(b2, 2)) += simd_u64(type2, 2); \
			case 7: (*(int*)simd_u64(b2, 3)) += simd_u64(type2, 3); \
			}

			switch (maxploidy)
			{
			case 10: REP_MID;
			case  9: REP_MID;
			case  8: REP_MID;
			case  7: REP_MID;
			case  6: REP_MID;
			case  5: REP_MID;
			case  4: REP_MID;
			case  3: REP_MID;
			case  2: REP_MID;
			case  1: REP_MID;
			}
#undef REP_MID
		}
	}
}

/* Update a priori ancetral proportion for non-admix model */
TARGETAVX void UpdateQNoAdmixAVX(double* Q, int K, double* bufNK1, double* bufNK2, double* Gamma, RNGAVX& rng, bool locpriori, double* Base, ushort* Z)
{
	__m256i mask1 = _mm256_set1_epi64x(0x7FF0000000000000);
	__m256i mask2 = _mm256_set1_epi64x(0x800FFFFFFFFFFFFF);
	__m256i mask3 = _mm256_set1_epi64x(0x3FF0000000000000);
	__m256i subv = _mm256_set1_epi64x(1023);
	__m256i maskff = _mm256_set1_epi64x(0xFFFFFFFFFFFFFFFF); 
	__m256i maskunder = _mm256_set1_epi64x(0x1FF0000000000000); 
	__m256d maskone = _mm256_set1_pd(1.0);
	__m256i vindex = _mm256_set_epi64x(3 * structure_indnbytes, 2 * structure_indnbytes, structure_indnbytes, 0);
	double* q = Q;

	for (int i = 0, pad = 0; i < nind; i += STRUCTURE_NPACK)
	{
		if (i + STRUCTURE_NPACK >= nind) //fix i for the last several
		{
			pad = STRUCTURE_NPACK - (nind - i);
			i = nind - STRUCTURE_NPACK;
		}

		//K elements, each save NPACK log or likelihood
		__m256i* slog = AlignSIMD((__m256i*)bufNK1);
		__m256d* prod = AlignSIMD((__m256d*)bufNK2);

		//slog and prod at K clusters
		OpenLog((int64*)slog, (double*)prod, K * STRUCTURE_NPACK);

		//add priori probability
		if (locpriori) for (int k = 0; k < K; ++k)
		{
			ChargeLog(*((int64*)slog + k * STRUCTURE_NPACK + 0), *((double*)prod + k * STRUCTURE_NPACK + 0), Gamma[ainds[i + 0]->popid * K + k]);
			ChargeLog(*((int64*)slog + k * STRUCTURE_NPACK + 1), *((double*)prod + k * STRUCTURE_NPACK + 1), Gamma[ainds[i + 1]->popid * K + k]);
			ChargeLog(*((int64*)slog + k * STRUCTURE_NPACK + 2), *((double*)prod + k * STRUCTURE_NPACK + 2), Gamma[ainds[i + 2]->popid * K + k]);
			ChargeLog(*((int64*)slog + k * STRUCTURE_NPACK + 3), *((double*)prod + k * STRUCTURE_NPACK + 3), Gamma[ainds[i + 3]->popid * K + k]);
			ChargeLog(*((int64*)slog + k * STRUCTURE_NPACK + 4), *((double*)prod + k * STRUCTURE_NPACK + 4), Gamma[ainds[i + 4]->popid * K + k]);
			ChargeLog(*((int64*)slog + k * STRUCTURE_NPACK + 5), *((double*)prod + k * STRUCTURE_NPACK + 5), Gamma[ainds[i + 5]->popid * K + k]);
			ChargeLog(*((int64*)slog + k * STRUCTURE_NPACK + 6), *((double*)prod + k * STRUCTURE_NPACK + 6), Gamma[ainds[i + 6]->popid * K + k]);
			ChargeLog(*((int64*)slog + k * STRUCTURE_NPACK + 7), *((double*)prod + k * STRUCTURE_NPACK + 7), Gamma[ainds[i + 7]->popid * K + k]);
		}

		double* p = Base;
		BAYESIAN_READER rt(i);
		for (int64 l = 0; l < nloc; ++l, p += GetLoc(l).k)
		{
			int size = structure_size[l];
			__m256i mask = _mm256_set1_epi64x((1u << size) - 1u);

			switch (maxploidy)
			{
			case 10: ReadFreqAVX/* (rt, size, slog, prod, p, K) */;
			case  9: ReadFreqAVX/* (rt, size, slog, prod, p, K) */;
			case  8: ReadFreqAVX/* (rt, size, slog, prod, p, K) */;
			case  7: ReadFreqAVX/* (rt, size, slog, prod, p, K) */;
			case  6: ReadFreqAVX/* (rt, size, slog, prod, p, K) */;
			case  5: ReadFreqAVX/* (rt, size, slog, prod, p, K) */;
			case  4: ReadFreqAVX/* (rt, size, slog, prod, p, K) */;
			case  3: ReadFreqAVX/* (rt, size, slog, prod, p, K) */;
			case  2: ReadFreqAVX/* (rt, size, slog, prod, p, K) */;
			case  1: ReadFreqAVX/* (rt, size, slog, prod, p, K) */;
			}
		}

		CloseLog((int64*)slog, (double*)prod, K * STRUCTURE_NPACK);

		__m256i ki[2];
		rng.PolyLog(prod, K, ki);

		for (int j = pad, ii = i + pad; j < STRUCTURE_NPACK; ++j, ++ii, q += K)
		{
			ushort k2 = simp_u64(ki, j);
			q[k2] = 1;
			Z[ii] = k2;
		}
	}
}

/* Update a priori ancetral proportion by Metropolis-Hastings for admix model */
TARGETAVX void UpdateQMetroAVX(double* Q, int K, double* bufNK1, double* bufN1, double* bufN2, double* Base)
{
	//add priori probability
	__m256i vindex = _mm256_set_epi64x(3 * structure_indnbytes, 2 * structure_indnbytes, structure_indnbytes, 0); 
	double* p = NULL, * q0 = Q, * b0 = bufNK1;
	OpenLog((int64*)bufN1, bufN2, nind);

	for (int i = 0, pad = 0; i < nind; i += STRUCTURE_NPACK, b0 += STRUCTURE_NPACK * K, q0 += STRUCTURE_NPACK * K)
	{
		if (i + STRUCTURE_NPACK >= nind) //fix i for the last several
		{
			pad = STRUCTURE_NPACK - (nind - i);
			i = nind - STRUCTURE_NPACK;
			q0 = Q + i * K;
			b0 = bufNK1 + i * K;
		}

		BAYESIAN_READER rt(i);
		p = Base;
		double* q1 = q0 + K, * q2 = q1 + K, * q3 = q2 + K, * q4 = q3 + K, * q5 = q4 + K, * q6 = q5 + K, * q7 = q6 + K;
		double* b1 = b0 + K, * b2 = b1 + K, * b3 = b2 + K, * b4 = b3 + K, * b5 = b4 + K, * b6 = b5 + K, * b7 = b6 + K;

		for (int64 l = 0; l < nloc; ++l, p += GetLoc(l).k)
		{
			int size = structure_size[l];
			__m256i mask = _mm256_set1_epi64x((1u << size) - 1u);
			__m256i aid1, aid2;
			__m256i type1, type2;
			__m256i paddr = _mm256_set1_epi64x((uint64)p);

#define REP_MID \
			ReadAidAVX/*(rt, size, aid, type)*/; \
			\
			aid1 = _mm256_slli_epi64(aid1, 3); \
			aid2 = _mm256_slli_epi64(aid2, 3); \
			\
			aid1 = _mm256_add_epi64(paddr, aid1); \
			aid2 = _mm256_add_epi64(paddr, aid2); \
			\
			switch (pad) \
			{ \
			case 0: if (simd_u64(type1, 0)) ChargeLog(*(int64*)&bufN1[i + 0], bufN2[i + 0], SumProdDiv(b0, q0, (double*)simd_u64(aid1, 0), KT, K)); \
			case 1: if (simd_u64(type1, 1)) ChargeLog(*(int64*)&bufN1[i + 1], bufN2[i + 1], SumProdDiv(b1, q1, (double*)simd_u64(aid1, 1), KT, K)); \
			case 2: if (simd_u64(type1, 2)) ChargeLog(*(int64*)&bufN1[i + 2], bufN2[i + 2], SumProdDiv(b2, q2, (double*)simd_u64(aid1, 2), KT, K)); \
			case 3: if (simd_u64(type1, 3)) ChargeLog(*(int64*)&bufN1[i + 3], bufN2[i + 3], SumProdDiv(b3, q3, (double*)simd_u64(aid1, 3), KT, K)); \
			case 4: if (simd_u64(type2, 0)) ChargeLog(*(int64*)&bufN1[i + 4], bufN2[i + 4], SumProdDiv(b4, q4, (double*)simd_u64(aid2, 0), KT, K)); \
			case 5: if (simd_u64(type2, 1)) ChargeLog(*(int64*)&bufN1[i + 5], bufN2[i + 5], SumProdDiv(b5, q5, (double*)simd_u64(aid2, 1), KT, K)); \
			case 6: if (simd_u64(type2, 2)) ChargeLog(*(int64*)&bufN1[i + 6], bufN2[i + 6], SumProdDiv(b6, q6, (double*)simd_u64(aid2, 2), KT, K)); \
			case 7: if (simd_u64(type2, 3)) ChargeLog(*(int64*)&bufN1[i + 7], bufN2[i + 7], SumProdDiv(b7, q7, (double*)simd_u64(aid2, 3), KT, K)); \
			}

			switch (maxploidy)
			{
			case 10: REP_MID;
			case  9: REP_MID;
			case  8: REP_MID;
			case  7: REP_MID;
			case  6: REP_MID;
			case  5: REP_MID;
			case  4: REP_MID;
			case  3: REP_MID;
			case  2: REP_MID;
			case  1: REP_MID;
			}
#undef REP_MID 
		}
	}
	CloseLog((int64*)bufN1, bufN2, nind);
}

/* Update individual or allele origin for admix model */
TARGETAVX void UpdateZAdmixAVX(double* Q, int K, int64* Mi, int* Ni, double* bufNK1, double* bufNK2, double* Base, RNGAVX& rng)
{
	__m256i mask01 = _mm256_set1_epi64x(1);
	__m256d maskone = _mm256_set1_pd(1.0);
	__m256d* freq = AlignSIMD((__m256d*)bufNK1);
	__m256d* qq = AlignSIMD((__m256d*)bufNK2);
	__m256d dmin = _mm256_set1_pd(MIN_FREQ);
	__m256i vindex = _mm256_set_epi64x(3 * structure_indnbytes, 2 * structure_indnbytes, structure_indnbytes, 0); 

	for (int i = 0, pad = 0; i < nind; i += STRUCTURE_NPACK)
	{
		if (i + STRUCTURE_NPACK >= nind) //fix i for the last several
		{
			pad = STRUCTURE_NPACK - (nind - i);
			i = nind - STRUCTURE_NPACK;
		}

		double* q0 = Q + i * K;
		int64* m0 = Mi + i * K, * m1 = m0 + K, * m2 = m1 + K, * m3 = m2 + K, * m4 = m3 + K, * m5 = m4 + K, * m6 = m5 + K, * m7 = m6 + K;
		__m256i vindex1 = _mm256_set_epi64x(3 * K, 2 * K, 1 * K, 0 * K);
		__m256i vindex2 = _mm256_set_epi64x(7 * K, 6 * K, 5 * K, 4 * K);
		
		for (int k = 0; k < K; ++k)
		{
			qq[k * 2 + 0] = _mm256_i64gather_pd(&q0[k], vindex1, sizeof(double));
			qq[k * 2 + 1] = _mm256_i64gather_pd(&q0[k], vindex2, sizeof(double));
		}

		int* ni = Ni;
		double* p = Base;
		BAYESIAN_READER rt(i);

		for (int64 l = 0; l < nloc; ++l, ni += GetLoc(l).k, p += GetLoc(l).k)
		{
			int size = structure_size[l];
			__m256i mask = _mm256_set1_epi64x((1u << size) - 1u);
			__m256i type[2];
			__m256i& type1 = type[0];
			__m256i& type2 = type[1];
			__m256i aid1, aid2;
			__m256i addr = _mm256_set1_epi64x((uint64)p);
			__m256d sumfreq[2];
			__m256d a1, a2;
			__m256i zid[2];
			__m256i& zid1 = zid[0];
			__m256i& zid2 = zid[1];
			ushort k2;

#define REP_MID \
			ReadAidAVX/*(rt, size, aid, type)*/; \
			\
			sumfreq[0] = _mm256_setzero_pd(); \
			sumfreq[1] = _mm256_setzero_pd(); \
			\
			for (int k = 0; k < K; ++k) \
			{ \
				a1 = _mm256_mask_i64gather_pd(maskone, p + k * KT, aid1, _mm256_castsi256_pd(type1), sizeof(double)); \
				a2 = _mm256_mask_i64gather_pd(maskone, p + k * KT, aid2, _mm256_castsi256_pd(type2), sizeof(double)); \
				\
				a1 = _mm256_mul_pd(qq[k * 2 + 0], a1); \
				a2 = _mm256_mul_pd(qq[k * 2 + 1], a2); \
				\
				freq[k * 2 + 0] = _mm256_add_pd(a1, dmin); \
				freq[k * 2 + 1] = _mm256_add_pd(a2, dmin); \
				\
				sumfreq[0] = _mm256_add_pd(sumfreq[0], freq[k * 2 + 0]); \
				sumfreq[1] = _mm256_add_pd(sumfreq[1], freq[k * 2 + 1]); \
			} \
			/* draw z */ \
			rng.Poly(freq, sumfreq, K, zid); \
			\
			type1 = _mm256_and_si256(type1, mask01); \
			type2 = _mm256_and_si256(type2, mask01); \
			\
			switch (pad) \
			{ /* update mi, ni */ \
			case 0: k2 = simd_u64(zid1, 0); m0[k2] += simd_u64(type1, 0); ni[k2 * KT + simd_u64(aid1, 0)] += simd_u64(type1, 0);  \
			case 1: k2 = simd_u64(zid1, 1); m1[k2] += simd_u64(type1, 1); ni[k2 * KT + simd_u64(aid1, 1)] += simd_u64(type1, 1);  \
			case 2: k2 = simd_u64(zid1, 2); m2[k2] += simd_u64(type1, 2); ni[k2 * KT + simd_u64(aid1, 2)] += simd_u64(type1, 2);  \
			case 3: k2 = simd_u64(zid1, 3); m3[k2] += simd_u64(type1, 3); ni[k2 * KT + simd_u64(aid1, 3)] += simd_u64(type1, 3);  \
			case 4: k2 = simd_u64(zid2, 0); m4[k2] += simd_u64(type2, 0); ni[k2 * KT + simd_u64(aid2, 0)] += simd_u64(type2, 0);  \
			case 5: k2 = simd_u64(zid2, 1); m5[k2] += simd_u64(type2, 1); ni[k2 * KT + simd_u64(aid2, 1)] += simd_u64(type2, 1);  \
			case 6: k2 = simd_u64(zid2, 2); m6[k2] += simd_u64(type2, 2); ni[k2 * KT + simd_u64(aid2, 2)] += simd_u64(type2, 2);  \
			case 7: k2 = simd_u64(zid2, 3); m7[k2] += simd_u64(type2, 3); ni[k2 * KT + simd_u64(aid2, 3)] += simd_u64(type2, 3);  \
			}
			
			switch (maxploidy)
			{
			case 10: REP_MID;
			case  9: REP_MID;
			case  8: REP_MID;
			case  7: REP_MID;
			case  6: REP_MID;
			case  5: REP_MID;
			case  4: REP_MID;
			case  3: REP_MID;
			case  2: REP_MID;
			case  1: REP_MID;
			}
#undef REP_MID 
		}
	}
}

#endif