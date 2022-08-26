/* SSE Instruction Set Functions */

#include "vcfpop.h"

#ifndef __aarch64__

/* Read allele id from decompressed bucket */
//ReadAidSSE(BAYESIAN_READER& rt, int size, __m128i* aid, __m128i* type)
#define ReadAidSSE \
{ \
	/* if data is empty */ \
	if (rt.nbits < size) [[unlikely]] \
	{ \
		/* move nbits data to gid */ \
		aid1 = rt.data128[0]; \
		aid2 = rt.data128[1]; \
		aid3 = rt.data128[2]; \
		aid4 = rt.data128[3]; \
 \
		/* remain number of bits to read */ \
		int rbits = size - rt.nbits; \
 \
		/* read 64 bits to data */ \
		uint64* tpos = rt.pos; \
 \
		rt.data64[0] = *tpos; tpos += structure_indnbytes; \
		rt.data64[1] = *tpos; tpos += structure_indnbytes; \
		rt.data64[2] = *tpos; tpos += structure_indnbytes; \
		rt.data64[3] = *tpos; tpos += structure_indnbytes; \
		rt.data64[4] = *tpos; tpos += structure_indnbytes; \
		rt.data64[5] = *tpos; tpos += structure_indnbytes; \
		rt.data64[6] = *tpos; tpos += structure_indnbytes; \
		rt.data64[7] = *tpos; tpos += structure_indnbytes; \
 \
		rt.pos++; \
 \
		/* read rbits from data and concate to higher bits in gid */ \
		__m128i maskrbits = _mm_set1_epi64x((1u << rbits) - 1u); \
		__m128i c1, c2, c3, c4; \
 \
		c1 = _mm_and_si128(rt.data128[0], maskrbits); \
		c2 = _mm_and_si128(rt.data128[1], maskrbits); \
		c3 = _mm_and_si128(rt.data128[2], maskrbits); \
		c4 = _mm_and_si128(rt.data128[3], maskrbits); \
 \
		c1 = _mm_slli_epi64(c1, rt.nbits); \
		c2 = _mm_slli_epi64(c2, rt.nbits); \
		c3 = _mm_slli_epi64(c3, rt.nbits); \
		c4 = _mm_slli_epi64(c4, rt.nbits); \
 \
		aid1 = _mm_or_si128(aid1, c1); \
		aid2 = _mm_or_si128(aid2, c2); \
		aid3 = _mm_or_si128(aid3, c3); \
		aid4 = _mm_or_si128(aid4, c4); \
 \
		/*shift right */ \
		rt.data128[0] = _mm_srli_epi64(rt.data128[0], rbits); \
		rt.data128[1] = _mm_srli_epi64(rt.data128[1], rbits); \
		rt.data128[2] = _mm_srli_epi64(rt.data128[2], rbits); \
		rt.data128[3] = _mm_srli_epi64(rt.data128[3], rbits); \
 \
		rt.nbits = 64 - rbits; \
	} \
	else [[likely]] \
	{ \
		/*read size bits */ \
		aid1 = _mm_and_si128(rt.data128[0], mask); \
		aid2 = _mm_and_si128(rt.data128[1], mask); \
		aid3 = _mm_and_si128(rt.data128[2], mask); \
		aid4 = _mm_and_si128(rt.data128[3], mask); \
 \
		/*shift right */ \
		rt.data128[0] = _mm_srli_epi64(rt.data128[0], size); \
		rt.data128[1] = _mm_srli_epi64(rt.data128[1], size); \
		rt.data128[2] = _mm_srli_epi64(rt.data128[2], size); \
		rt.data128[3] = _mm_srli_epi64(rt.data128[3], size); \
 \
		rt.nbits -= size; \
	} \
 \
	/*typed is 0xFFFFFFFF */ \
	type1 = _mm_cmpgt_epi64(mask, aid1); \
	type2 = _mm_cmpgt_epi64(mask, aid2); \
	type3 = _mm_cmpgt_epi64(mask, aid3); \
	type4 = _mm_cmpgt_epi64(mask, aid4); \
 \
	/* set missing to 0 */ \
	aid1 = _mm_and_si128(aid1, type1); \
	aid2 = _mm_and_si128(aid2, type2); \
	aid3 = _mm_and_si128(aid3, type3); \
	aid4 = _mm_and_si128(aid4, type4); \
}

/* Read allele frequency from decompressed bucket */
/*
TARGETSSE void ReadFreqSSE(BAYESIAN_READER& rt, int size, __m128i* slog, __m128d* prod, double* p, int K)
__m128i maskff = _mm_set1_epi64x(0xFFFFFFFFFFFFFFFF);
__m128i maskunder = _mm_set1_epi64x(0x1FF0000000000000);
__m128i mask1 = _mm_set1_epi64x(0x7FF0000000000000);
__m128i mask2 = _mm_set1_epi64x(0x800FFFFFFFFFFFFF);
__m128i mask3 = _mm_set1_epi64x(0x3FF0000000000000);
__m128i subv = _mm_set1_epi64x(1023);
__m128i addr_inc = _mm_set1_epi64x(KT * sizeof(double));
__m128d maskone = _mm_set1_pd(1.0);
*/
#define ReadFreqSSE \
{ \
	__m128i aid1, aid2, aid3, aid4; \
	__m128i a1, a2, a3, a4; \
	__m128i b1, b2, b3, b4; \
	__m128i type1, type2, type3, type4; \
	__m128d freq1, freq2, freq3, freq4; \
	__m128i* prodi = (__m128i*)prod; \
 \
	ReadAidSSE \
 \
	/* x sizeof(double) to load frequency */ \
	aid1 = _mm_slli_epi64(aid1, 3); \
	aid2 = _mm_slli_epi64(aid2, 3); \
	aid3 = _mm_slli_epi64(aid3, 3); \
	aid4 = _mm_slli_epi64(aid4, 3); \
 \
	__m128i addr = _mm_set1_epi64x((uint64)p); \
 \
	aid1 = _mm_add_epi64(aid1, addr); \
	aid2 = _mm_add_epi64(aid2, addr); \
	aid3 = _mm_add_epi64(aid3, addr); \
	aid4 = _mm_add_epi64(aid4, addr); \
 \
	for (int k = 0; k < K; ++k) \
	{ \
		/* get freq of 8 allele copies for 8 individuals in cluster k */ \
		freq1 = _mm_set_pd(*(double*)simd_u64(aid1, 1), *(double*)simd_u64(aid1, 0)); \
		freq2 = _mm_set_pd(*(double*)simd_u64(aid2, 1), *(double*)simd_u64(aid2, 0)); \
		freq3 = _mm_set_pd(*(double*)simd_u64(aid3, 1), *(double*)simd_u64(aid3, 0)); \
		freq4 = _mm_set_pd(*(double*)simd_u64(aid4, 1), *(double*)simd_u64(aid4, 0)); \
 \
		/* set missing freq to one */ \
		freq1 = _mm_blendv_pd(maskone, freq1, _mm_castsi128_pd(type1)); \
		freq2 = _mm_blendv_pd(maskone, freq2, _mm_castsi128_pd(type2)); \
		freq3 = _mm_blendv_pd(maskone, freq3, _mm_castsi128_pd(type3)); \
		freq4 = _mm_blendv_pd(maskone, freq4, _mm_castsi128_pd(type4)); \
 \
		/* mul to prod */ \
		prod[k * 4 + 0] = _mm_mul_pd(prod[k * 4 + 0], freq1); \
		prod[k * 4 + 1] = _mm_mul_pd(prod[k * 4 + 1], freq2); \
		prod[k * 4 + 2] = _mm_mul_pd(prod[k * 4 + 2], freq3); \
		prod[k * 4 + 3] = _mm_mul_pd(prod[k * 4 + 3], freq4); \
 \
		aid1 = _mm_add_epi64(aid1, addr_inc); \
		aid2 = _mm_add_epi64(aid2, addr_inc); \
		aid3 = _mm_add_epi64(aid3, addr_inc); \
		aid4 = _mm_add_epi64(aid4, addr_inc); \
 \
		a1 = _mm_cmpgt_epi64(maskunder, _mm_castpd_si128(prod[k * 4 + 0])); \
		a2 = _mm_cmpgt_epi64(maskunder, _mm_castpd_si128(prod[k * 4 + 1])); \
		a3 = _mm_cmpgt_epi64(maskunder, _mm_castpd_si128(prod[k * 4 + 2])); \
		a4 = _mm_cmpgt_epi64(maskunder, _mm_castpd_si128(prod[k * 4 + 3])); \
 \
		a1 = _mm_or_si128(a1, a2); \
		a3 = _mm_or_si128(a3, a4); \
		a1 = _mm_or_si128(a1, a3); \
 \
		/* check underflow */ \
		if (!(simd_u64(a1, 0) || simd_u64(a1, 1))) [[likely]] continue; \
 \
		/* add exponent */ \
		a1 = _mm_and_si128(prodi[k * 4 + 0], mask1); \
		a2 = _mm_and_si128(prodi[k * 4 + 1], mask1); \
		a3 = _mm_and_si128(prodi[k * 4 + 2], mask1); \
		a4 = _mm_and_si128(prodi[k * 4 + 3], mask1); \
 \
		b1 = _mm_and_si128(prodi[k * 4 + 0], mask2); \
		b2 = _mm_and_si128(prodi[k * 4 + 1], mask2); \
		b3 = _mm_and_si128(prodi[k * 4 + 2], mask2); \
		b4 = _mm_and_si128(prodi[k * 4 + 3], mask2); \
 \
		a1 = _mm_srli_epi64(a1, 52); \
		a2 = _mm_srli_epi64(a2, 52); \
		a3 = _mm_srli_epi64(a3, 52); \
		a4 = _mm_srli_epi64(a4, 52); \
 \
		prodi[k * 4 + 0] = _mm_or_si128(b1, mask3); \
		prodi[k * 4 + 1] = _mm_or_si128(b2, mask3); \
		prodi[k * 4 + 2] = _mm_or_si128(b3, mask3); \
		prodi[k * 4 + 3] = _mm_or_si128(b4, mask3); \
 \
		a1 = _mm_sub_epi64(a1, subv); \
		a2 = _mm_sub_epi64(a2, subv); \
		a3 = _mm_sub_epi64(a3, subv); \
		a4 = _mm_sub_epi64(a4, subv); \
 \
		slog[k * 4 + 0] = _mm_add_epi64(slog[k * 4 + 0], a1); \
		slog[k * 4 + 1] = _mm_add_epi64(slog[k * 4 + 1], a2); \
		slog[k * 4 + 2] = _mm_add_epi64(slog[k * 4 + 2], a3); \
		slog[k * 4 + 3] = _mm_add_epi64(slog[k * 4 + 3], a4); \
 \
	} \
}

/* Update individual or allele origin when ancetral proportion is binary */
TARGETSSE void UpdateZBinarySSE(int* Ni, ushort* Z)
{
	__m128i Kt = _mm_set1_epi64x(KT);
	__m128i mask01 = _mm_set1_epi64x(1);

	for (int i = 0, pad = 0; i < nind; i += STRUCTURE_NPACK)
	{
		if (i + STRUCTURE_NPACK >= nind) //fix i for the last several
		{
			pad = STRUCTURE_NPACK - (nind - i);
			i = nind - STRUCTURE_NPACK;
		}

		__m128i ZxKT1, ZxKT2, ZxKT3, ZxKT4;

		ZxKT1 = _mm_loadu_si128((__m128i*)&Z[i + 0]);
		ZxKT2 = _mm_loadu_si128((__m128i*)&Z[i + 2]);
		ZxKT3 = _mm_loadu_si128((__m128i*)&Z[i + 4]);
		ZxKT4 = _mm_loadu_si128((__m128i*)&Z[i + 6]);

		ZxKT1 = _mm_cvtepi16_epi64(ZxKT1);
		ZxKT2 = _mm_cvtepi16_epi64(ZxKT2);
		ZxKT3 = _mm_cvtepi16_epi64(ZxKT3);
		ZxKT4 = _mm_cvtepi16_epi64(ZxKT4);

		ZxKT1 = _mm_mul_epu32(Kt, ZxKT1);
		ZxKT2 = _mm_mul_epu32(Kt, ZxKT2);
		ZxKT3 = _mm_mul_epu32(Kt, ZxKT3);
		ZxKT4 = _mm_mul_epu32(Kt, ZxKT4);

		ZxKT1 = _mm_slli_epi64(ZxKT1, 2);
		ZxKT2 = _mm_slli_epi64(ZxKT2, 2);
		ZxKT3 = _mm_slli_epi64(ZxKT3, 2);
		ZxKT4 = _mm_slli_epi64(ZxKT4, 2);

		BAYESIAN_READER rt(i);
		int* ni = Ni;

		for (int64 l = 0; l < nloc; ++l, ni += GetLoc(l).k)
		{
			int size = structure_size[l];
			__m128i mask = _mm_set1_epi64x((1u << size) - 1u); 
			__m128i aid1, aid2, aid3, aid4;
			__m128i type1, type2, type3, type4;
			__m128i npio2 = _mm_set1_epi64x((int64)ni);
			__m128i a1, a2, a3, a4;
			__m128i b1, b2, b3, b4;

			a1 = _mm_add_epi64(ZxKT1, npio2);
			a2 = _mm_add_epi64(ZxKT2, npio2);
			a3 = _mm_add_epi64(ZxKT3, npio2);
			a4 = _mm_add_epi64(ZxKT4, npio2);

#define REP_MID \
			ReadAidSSE/*(rt, size, aid, type);*/ \
			\
			b1 = _mm_slli_epi64(aid1, 2); \
			b2 = _mm_slli_epi64(aid2, 2); \
			b3 = _mm_slli_epi64(aid3, 2); \
			b4 = _mm_slli_epi64(aid4, 2); \
			\
			b1 = _mm_add_epi64(a1, b1); \
			b2 = _mm_add_epi64(a2, b2); \
			b3 = _mm_add_epi64(a3, b3); \
			b4 = _mm_add_epi64(a4, b4); \
			\
			type1 = _mm_and_si128(type1, mask01);\
			type2 = _mm_and_si128(type2, mask01);\
			type3 = _mm_and_si128(type3, mask01);\
			type4 = _mm_and_si128(type4, mask01);\
			\
			switch (pad) \
			{ /* Z and Mi is updated before call, update Ni here */ \
			case 0: (*(int*)simd_u64(b1, 0)) += simd_u64(type1, 0); \
			case 1: (*(int*)simd_u64(b1, 1)) += simd_u64(type1, 1); \
			case 2: (*(int*)simd_u64(b2, 0)) += simd_u64(type2, 0); \
			case 3: (*(int*)simd_u64(b2, 1)) += simd_u64(type2, 1); \
			case 4: (*(int*)simd_u64(b3, 0)) += simd_u64(type3, 0); \
			case 5: (*(int*)simd_u64(b3, 1)) += simd_u64(type3, 1); \
			case 6: (*(int*)simd_u64(b4, 0)) += simd_u64(type4, 0); \
			case 7: (*(int*)simd_u64(b4, 1)) += simd_u64(type4, 1); \
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
TARGETSSE void UpdateQNoAdmixSSE(double* Q, int K, double* bufNK1, double* bufNK2, double* Gamma, RNGSSE& rng, bool locpriori, double* Base, ushort* Z)
{
	__m128i maskff = _mm_set1_epi64x(0xFFFFFFFFFFFFFFFF); 
	__m128i maskunder = _mm_set1_epi64x(0x1FF0000000000000);
	__m128i mask1 = _mm_set1_epi64x(0x7FF0000000000000);
	__m128i mask2 = _mm_set1_epi64x(0x800FFFFFFFFFFFFF);
	__m128i mask3 = _mm_set1_epi64x(0x3FF0000000000000);
	__m128i subv = _mm_set1_epi64x(1023);
	__m128i addr_inc = _mm_set1_epi64x(KT * sizeof(double));
	__m128d maskone = _mm_set1_pd(1.0);

	double* q = Q;

	for (int i = 0, pad = 0; i < nind; i += STRUCTURE_NPACK)
	{
		if (i + STRUCTURE_NPACK >= nind) //fix i for the last several
		{
			pad = STRUCTURE_NPACK - (nind - i);
			i = nind - STRUCTURE_NPACK;
		}

		//K elements, each save NPACK log or likelihood
		__m128i* slog = AlignSIMD((__m128i*)bufNK1);
		__m128d* prod = AlignSIMD((__m128d*)bufNK2);

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
			__m128i mask = _mm_set1_epi64x((1u << size) - 1u);

			switch (maxploidy)
			{
			case 10: ReadFreqSSE/* (rt, size, slog, prod, p, K) */;
			case  9: ReadFreqSSE/* (rt, size, slog, prod, p, K) */;
			case  8: ReadFreqSSE/* (rt, size, slog, prod, p, K) */;
			case  7: ReadFreqSSE/* (rt, size, slog, prod, p, K) */;
			case  6: ReadFreqSSE/* (rt, size, slog, prod, p, K) */;
			case  5: ReadFreqSSE/* (rt, size, slog, prod, p, K) */;
			case  4: ReadFreqSSE/* (rt, size, slog, prod, p, K) */;
			case  3: ReadFreqSSE/* (rt, size, slog, prod, p, K) */;
			case  2: ReadFreqSSE/* (rt, size, slog, prod, p, K) */;
			case  1: ReadFreqSSE/* (rt, size, slog, prod, p, K) */;
			}
		}

		CloseLog((int64*)slog, (double*)prod, K * STRUCTURE_NPACK);

		__m128i ki[4];
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
TARGETSSE void UpdateQMetroSSE(double* Q, int K, double* bufNK1, double* bufN1, double* bufN2, double* Base)
{
	//add priori probability
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
			__m128i mask = _mm_set1_epi64x((1u << size) - 1u);
			__m128i aid1, aid2, aid3, aid4;
			__m128i type1, type2, type3, type4;
			__m128i paddr = _mm_set1_epi64x((uint64)p);

#define REP_MID \
			ReadAidSSE/*(rt, size, aid, type)*/; \
			\
			aid1 = _mm_slli_epi64(aid1, 3); \
			aid2 = _mm_slli_epi64(aid2, 3); \
			aid3 = _mm_slli_epi64(aid3, 3); \
			aid4 = _mm_slli_epi64(aid4, 3); \
			\
			aid1 = _mm_add_epi64(paddr, aid1); \
			aid2 = _mm_add_epi64(paddr, aid2); \
			aid3 = _mm_add_epi64(paddr, aid3); \
			aid4 = _mm_add_epi64(paddr, aid4); \
			\
			switch (pad) \
			{ \
			case 0: if (simd_u64(type1, 0)) ChargeLog(*(int64*)&bufN1[i + 0], bufN2[i + 0], SumProdDiv(b0, q0, (double*)simd_u64(aid1, 0), KT, K)); \
			case 1: if (simd_u64(type1, 1)) ChargeLog(*(int64*)&bufN1[i + 1], bufN2[i + 1], SumProdDiv(b1, q1, (double*)simd_u64(aid1, 1), KT, K)); \
			case 2: if (simd_u64(type2, 0)) ChargeLog(*(int64*)&bufN1[i + 2], bufN2[i + 2], SumProdDiv(b2, q2, (double*)simd_u64(aid2, 0), KT, K)); \
			case 3: if (simd_u64(type2, 1)) ChargeLog(*(int64*)&bufN1[i + 3], bufN2[i + 3], SumProdDiv(b3, q3, (double*)simd_u64(aid2, 1), KT, K)); \
			case 4: if (simd_u64(type3, 0)) ChargeLog(*(int64*)&bufN1[i + 4], bufN2[i + 4], SumProdDiv(b4, q4, (double*)simd_u64(aid3, 0), KT, K)); \
			case 5: if (simd_u64(type3, 1)) ChargeLog(*(int64*)&bufN1[i + 5], bufN2[i + 5], SumProdDiv(b5, q5, (double*)simd_u64(aid3, 1), KT, K)); \
			case 6: if (simd_u64(type4, 0)) ChargeLog(*(int64*)&bufN1[i + 6], bufN2[i + 6], SumProdDiv(b6, q6, (double*)simd_u64(aid4, 0), KT, K)); \
			case 7: if (simd_u64(type4, 1)) ChargeLog(*(int64*)&bufN1[i + 7], bufN2[i + 7], SumProdDiv(b7, q7, (double*)simd_u64(aid4, 1), KT, K)); \
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
TARGETSSE void UpdateZAdmixSSE(double* Q, int K, int64* Mi, int* Ni, double* bufNK1, double* bufNK2, double* Base, RNGSSE& rng)
{
	__m128i mask01 = _mm_set1_epi64x(1);
	__m128d maskone = _mm_set1_pd(1.0);
	__m128d* freq = AlignSIMD((__m128d*)bufNK1);
	__m128d* qq = AlignSIMD((__m128d*)bufNK2);
	__m128d dmin = _mm_set1_pd(MIN_FREQ);
	__m128i addr_inc = _mm_set1_epi64x(KT * sizeof(double));

	for (int i = 0, pad = 0; i < nind; i += STRUCTURE_NPACK)
	{
		if (i + STRUCTURE_NPACK >= nind) //fix i for the last several
		{
			pad = STRUCTURE_NPACK - (nind - i);
			i = nind - STRUCTURE_NPACK;
		}

		double* q0 = Q + i * K, * q1 = q0 + K, * q2 = q1 + K, * q3 = q2 + K, * q4 = q3 + K, * q5 = q4 + K, * q6 = q5 + K, * q7 = q6 + K;
		int64* m0 = Mi + i * K, * m1 = m0 + K, * m2 = m1 + K, * m3 = m2 + K, * m4 = m3 + K, * m5 = m4 + K, * m6 = m5 + K, * m7 = m6 + K;

		for (int k = 0; k < K; ++k)
		{
			qq[k * 4 + 0] = _mm_set_pd(q1[k], q0[k]);
			qq[k * 4 + 1] = _mm_set_pd(q3[k], q2[k]);
			qq[k * 4 + 2] = _mm_set_pd(q5[k], q4[k]);
			qq[k * 4 + 3] = _mm_set_pd(q7[k], q6[k]);
		}

		int* ni = Ni;
		double* p = Base;
		BAYESIAN_READER rt(i);

		for (int64 l = 0; l < nloc; ++l, ni += GetLoc(l).k, p += GetLoc(l).k)
		{
			int size = structure_size[l];
			__m128i mask = _mm_set1_epi64x((1u << size) - 1u);
			__m128i type[4];
			__m128i aid1, aid2, aid3, aid4;
			__m128i& type1 = type[0];
			__m128i& type2 = type[1];
			__m128i& type3 = type[2];
			__m128i& type4 = type[3];
			__m128i t1, t2, t3, t4;
			__m128i addr = _mm_set1_epi64x((uint64)p);
			__m128d sumfreq[4];
			__m128d a1, a2, a3, a4;
			__m128i zid[4];
			__m128i& zid1 = zid[0];
			__m128i& zid2 = zid[1];
			__m128i& zid3 = zid[2];
			__m128i& zid4 = zid[3];
			ushort k2;

#define REP_MID \
			ReadAidSSE/*(rt, size, aid, type)*/; \
			\
			sumfreq[0] = _mm_setzero_pd(); \
			sumfreq[1] = _mm_setzero_pd(); \
			sumfreq[2] = _mm_setzero_pd(); \
			sumfreq[3] = _mm_setzero_pd(); \
			\
			t1 = _mm_slli_epi64(aid1, 3); \
			t2 = _mm_slli_epi64(aid2, 3); \
			t3 = _mm_slli_epi64(aid3, 3); \
			t4 = _mm_slli_epi64(aid4, 3); \
			\
			t1 = _mm_add_epi64(t1, addr); \
			t2 = _mm_add_epi64(t2, addr); \
			t3 = _mm_add_epi64(t3, addr); \
			t4 = _mm_add_epi64(t4, addr); \
			\
			for (int k = 0; k < K; ++k) \
			{ \
				a1 = _mm_set_pd(*(double*)simd_u64(t1, 1), *(double*)simd_u64(t1, 0)); \
				a2 = _mm_set_pd(*(double*)simd_u64(t2, 1), *(double*)simd_u64(t2, 0)); \
				a3 = _mm_set_pd(*(double*)simd_u64(t3, 1), *(double*)simd_u64(t3, 0)); \
				a4 = _mm_set_pd(*(double*)simd_u64(t4, 1), *(double*)simd_u64(t4, 0)); \
				\
				a1 = _mm_mul_pd(qq[k * 4 + 0], a1); \
				a2 = _mm_mul_pd(qq[k * 4 + 1], a2); \
				a3 = _mm_mul_pd(qq[k * 4 + 2], a3); \
				a4 = _mm_mul_pd(qq[k * 4 + 3], a4); \
				\
				freq[k * 4 + 0] = _mm_add_pd(a1, dmin); \
				freq[k * 4 + 1] = _mm_add_pd(a2, dmin); \
				freq[k * 4 + 2] = _mm_add_pd(a3, dmin); \
				freq[k * 4 + 3] = _mm_add_pd(a4, dmin); \
				\
				sumfreq[0] = _mm_add_pd(sumfreq[0], freq[k * 4 + 0]); \
				sumfreq[1] = _mm_add_pd(sumfreq[1], freq[k * 4 + 1]); \
				sumfreq[2] = _mm_add_pd(sumfreq[2], freq[k * 4 + 2]); \
				sumfreq[3] = _mm_add_pd(sumfreq[3], freq[k * 4 + 3]); \
				\
				t1 = _mm_add_epi64(t1, addr_inc); \
				t2 = _mm_add_epi64(t2, addr_inc); \
				t3 = _mm_add_epi64(t3, addr_inc); \
				t4 = _mm_add_epi64(t4, addr_inc); \
			} \
			/* draw z */ \
			rng.Poly(freq, sumfreq, K, zid); \
			\
			type1 = _mm_and_si128(type1, mask01);\
			type2 = _mm_and_si128(type2, mask01);\
			type3 = _mm_and_si128(type3, mask01);\
			type4 = _mm_and_si128(type4, mask01);\
			\
			switch (pad) \
			{ /* update mi, ni */ \
			case 0: k2 = simd_u64(zid1, 0); m0[k2] += simd_u64(type1, 0); ni[k2 * KT + simd_u64(aid1, 0)] += simd_u64(type1, 0); \
			case 1: k2 = simd_u64(zid1, 1); m1[k2] += simd_u64(type1, 1); ni[k2 * KT + simd_u64(aid1, 1)] += simd_u64(type1, 1); \
			case 2: k2 = simd_u64(zid2, 0); m2[k2] += simd_u64(type2, 0); ni[k2 * KT + simd_u64(aid2, 0)] += simd_u64(type2, 0); \
			case 3: k2 = simd_u64(zid2, 1); m3[k2] += simd_u64(type2, 1); ni[k2 * KT + simd_u64(aid2, 1)] += simd_u64(type2, 1); \
			case 4: k2 = simd_u64(zid3, 0); m4[k2] += simd_u64(type3, 0); ni[k2 * KT + simd_u64(aid3, 0)] += simd_u64(type3, 0); \
			case 5: k2 = simd_u64(zid3, 1); m5[k2] += simd_u64(type3, 1); ni[k2 * KT + simd_u64(aid3, 1)] += simd_u64(type3, 1); \
			case 6: k2 = simd_u64(zid4, 0); m6[k2] += simd_u64(type4, 0); ni[k2 * KT + simd_u64(aid4, 0)] += simd_u64(type4, 0); \
			case 7: k2 = simd_u64(zid4, 1); m7[k2] += simd_u64(type4, 1); ni[k2 * KT + simd_u64(aid4, 1)] += simd_u64(type4, 1); \
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
