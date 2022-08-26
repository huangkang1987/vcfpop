/* AVX512 Instruction Functions */

#include "vcfpop.h"

#ifndef __aarch64__

/* Read allele id from decompressed bucket */
//TARGET512 void ReadAid512(BAYESIAN_READER& rt, int size, __m512i& aid, byte& type)
#define ReadAid512 \
{ \
	/* if data is empty */ \
	if (rt.nbits < size) [[unlikely]] \
	{ \
		/* move nbits data to gid */ \
		aid = rt.data512; \
 \
		/* remain number of bits to read */ \
		int rbits = size - rt.nbits; \
 \
		/* read 64 bits to data */ \
		rt.data512 = _mm512_i64gather_epi64(vindex, (int64*)rt.pos, sizeof(int64)); \
 \
		rt.pos++; \
 \
		/* read rbits from data and concate to higher bits in gid */ \
		aid = _mm512_or_si512(aid,_mm512_slli_epi64(_mm512_and_si512(rt.data512, _mm512_set1_epi64((1u << rbits) - 1u)), rt.nbits)); \
 \
		/*shift right */ \
		rt.data512 = _mm512_srli_epi64(rt.data512, rbits); \
\
		rt.nbits = 64 - rbits; \
	} \
	else [[likely]] \
	{ \
		/*read size bits */ \
		aid = _mm512_and_si512(rt.data512, mask); \
 \
		/*shift right */ \
		rt.data512 = _mm512_srli_epi64(rt.data512, size); \
 \
		rt.nbits -= size; \
	} \
 \
	type = _mm512_cmpneq_epi64_mask(mask, aid); \
 \
	/*set missing to -1 */ \
	aid = _mm512_mask_mov_epi64(maskff, type, aid); \
}

/* Read allele frequency from decompressed bucket */
/*TARGET512 void ReadFreq512(BAYESIAN_READER& rt, int size, __m512i* slog, __m512d* prod, double* p, int K)
{
	__m512i vindex = _mm512_set_epi64(
		7 * structure_indnbytes, 6 * structure_indnbytes, 5 * structure_indnbytes, 4 * structure_indnbytes,
		3 * structure_indnbytes, 2 * structure_indnbytes, structure_indnbytes, 0);
	__m512i mask1 = _mm512_set1_epi64(0x7FF0000000000000);
	__m512i mask2 = _mm512_set1_epi64(0x800FFFFFFFFFFFFF);
	__m512i mask3 = _mm512_set1_epi64(0x3FF0000000000000);
	__m512i subv = _mm512_set1_epi64(1023);
	__m512i maskff = _mm512_set1_epi64(0xFFFFFFFFFFFFFFFF);
	__m512i maskunder = _mm512_set1_epi64(0x1FF0000000000000);
	__m512d maskone = _mm512_set1_pd(1.0);
*/
#define ReadFreq512 \
{ \
	__m512i aid, a, b; \
	__m512d freq; \
	__m512i* prodi = (__m512i*)prod; \
	byte type; \
 \
	ReadAid512 \
\
	for (int k = 0; k < K; ++k) \
	{ \
		/* get freq of 8 allele copies in NPACK individuals and cluster k, and set missing freq to one */ \
		freq = _mm512_mask_i64gather_pd(maskone, type, aid, p + k * KT, sizeof(double)); \
 \
		/* mul to prod */ \
		prod[k] = _mm512_mul_pd(prod[k], freq); \
 \
		/* check underflow */ \
		if (_mm512_cmplt_epi64_mask(maskunder, _mm512_castpd_si512(prod[k]))) [[likely]] continue; \
 \
		/* add component to slog */ \
		a = _mm512_and_si512(prodi[k], mask1); \
 \
		b = _mm512_and_si512(prodi[k], mask2); \
 \
		a = _mm512_srli_epi64(a, 52); \
 \
		prodi[k] = _mm512_or_si512(b, mask3); \
 \
		a = _mm512_sub_epi64(a, subv); \
 \
		slog[k] = _mm512_add_epi64(slog[k], a); \
	} \
}

/* Update individual or allele origin when ancetral proportion is binary */
TARGET512 void UpdateZBinary512(int* Ni, ushort* Z)
{
	__m512i Kt = _mm512_set1_epi64(KT);
	__m512i maskff = _mm512_set1_epi64(0xFFFFFFFFFFFFFFFF);
	__m512i vindex = _mm512_set_epi64(
		7 * structure_indnbytes, 6 * structure_indnbytes, 5 * structure_indnbytes, 4 * structure_indnbytes, 
		3 * structure_indnbytes, 2 * structure_indnbytes, structure_indnbytes, 0); 

	for (int i = 0, pad = 0; i < nind; i += STRUCTURE_NPACK)
	{
		if (i + STRUCTURE_NPACK >= nind) //fix i for the last several
		{
			pad = STRUCTURE_NPACK - (nind - i);
			i = nind - STRUCTURE_NPACK;
		}

		__m128i t1 = _mm_loadu_si128((__m128i*)&Z[i]);
		__m512i ZxKT = _mm512_cvtepi16_epi64(t1);
		ZxKT = _mm512_mul_epu32(Kt, ZxKT);
		ZxKT = _mm512_slli_epi64(ZxKT, 2);

		BAYESIAN_READER rt(i);
		int* ni = Ni;

		for (int64 l = 0; l < nloc; ++l, ni += GetLoc(l).k)
		{
			int size = structure_size[l];
			__m512i mask = _mm512_set1_epi64((1u << size) - 1u);
			__m512i aid;
			__m512i b;
			byte type;
			__m512i a = _mm512_add_epi64(ZxKT, _mm512_set1_epi64((int64)ni));

#define REP_MID \
			ReadAid512 /*(rt, size, aid, type);*/ \
			b = _mm512_slli_epi64(aid, 2); \
			\
			b = _mm512_add_epi64(a, b); \
			\
			switch (pad) \
			{ /* Z and Mi is updated before call, update Ni here */ \
			case 0: (*(int*)simd_u64(b, 0)) += (type &   1) != 0; \
			case 1: (*(int*)simd_u64(b, 1)) += (type &   2) != 0; \
			case 2: (*(int*)simd_u64(b, 2)) += (type &   4) != 0; \
			case 3: (*(int*)simd_u64(b, 3)) += (type &   8) != 0; \
			case 4: (*(int*)simd_u64(b, 4)) += (type &  16) != 0; \
			case 5: (*(int*)simd_u64(b, 5)) += (type &  32) != 0; \
			case 6: (*(int*)simd_u64(b, 6)) += (type &  64) != 0; \
			case 7: (*(int*)simd_u64(b, 7)) += (type & 128) != 0; \
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
TARGET512 void UpdateQNoAdmix512(double* Q, int K, double* bufNK1, double* bufNK2, double* Gamma, RNG512& rng, bool locpriori, double* Base, ushort* Z)
{
	__m512i mask1 = _mm512_set1_epi64(0x7FF0000000000000);
	__m512i mask2 = _mm512_set1_epi64(0x800FFFFFFFFFFFFF);
	__m512i mask3 = _mm512_set1_epi64(0x3FF0000000000000);
	__m512i subv = _mm512_set1_epi64(1023);
	__m512i vindex = _mm512_set_epi64(
		7 * structure_indnbytes, 6 * structure_indnbytes, 5 * structure_indnbytes, 4 * structure_indnbytes,
		3 * structure_indnbytes, 2 * structure_indnbytes, structure_indnbytes, 0);
	__m512i maskff = _mm512_set1_epi64(0xFFFFFFFFFFFFFFFF); 
	__m512i maskunder = _mm512_set1_epi64(0x1FF0000000000000); 
	__m512d maskone = _mm512_set1_pd(1.0); 

	double* q = Q;

	for (int i = 0, pad = 0; i < nind; i += STRUCTURE_NPACK)
	{
		if (i + STRUCTURE_NPACK >= nind) //fix i for the last several
		{
			pad = STRUCTURE_NPACK - (nind - i);
			i = nind - STRUCTURE_NPACK;
		}

		//K elements, each save NPACK log or likelihood
		__m512i* slog = AlignSIMD((__m512i*)bufNK1);
		__m512d* prod = AlignSIMD((__m512d*)bufNK2);

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
			__m512i mask = _mm512_set1_epi64((1u << size) - 1u); 

			switch (maxploidy)
			{
			case 10: ReadFreq512/*(rt, size, slog, prod, p, K)*/;
			case  9: ReadFreq512/*(rt, size, slog, prod, p, K)*/;
			case  8: ReadFreq512/*(rt, size, slog, prod, p, K)*/;
			case  7: ReadFreq512/*(rt, size, slog, prod, p, K)*/;
			case  6: ReadFreq512/*(rt, size, slog, prod, p, K)*/;
			case  5: ReadFreq512/*(rt, size, slog, prod, p, K)*/;
			case  4: ReadFreq512/*(rt, size, slog, prod, p, K)*/;
			case  3: ReadFreq512/*(rt, size, slog, prod, p, K)*/;
			case  2: ReadFreq512/*(rt, size, slog, prod, p, K)*/;
			case  1: ReadFreq512/*(rt, size, slog, prod, p, K)*/;
			}
		}

		CloseLog((int64*)slog, (double*)prod, K * STRUCTURE_NPACK);

		__m512i ki;
		rng.PolyLog(prod, K, ki);

		for (int j = pad, ii = i + pad; j < STRUCTURE_NPACK; ++j, ++ii, q += K)
		{
			ushort k2 = simd_u64(ki, j);
			q[k2] = 1;
			Z[ii] = k2;
		}
	}
}

/* Update a priori ancetral proportion by Metropolis-Hastings for admix model */
TARGET512 void UpdateQMetro512(double* Q, int K, double* bufNK1, double* bufN1, double* bufN2, double* Base)
{
	//add priori probability
	double* p = NULL, * q0 = Q, * b0 = bufNK1;
	__m512i maskff = _mm512_set1_epi64(0xFFFFFFFFFFFFFFFF);
	__m512i vindex = _mm512_set_epi64(
		7 * structure_indnbytes, 6 * structure_indnbytes, 5 * structure_indnbytes, 4 * structure_indnbytes, 
		3 * structure_indnbytes, 2 * structure_indnbytes, structure_indnbytes, 0); 

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
			__m512i mask = _mm512_set1_epi64((1u << size) - 1u); 
			__m512i aid;
			byte type;
			__m512i paddr = _mm512_set1_epi64((uint64)p);

#define REP_MID \
			ReadAid512 /*(rt, size, aid, type);*/ \
			\
			aid = _mm512_add_epi64(paddr, _mm512_slli_epi64(aid, 3)); \
			\
			switch(pad) \
			{ \
			case 0: if (type &   1) ChargeLog(*(int64*)&bufN1[i + 0], bufN2[i + 0], SumProdDiv(b0, q0, (double*)simd_u64(aid, 0), KT, K)); \
			case 1: if (type &   2) ChargeLog(*(int64*)&bufN1[i + 1], bufN2[i + 1], SumProdDiv(b1, q1, (double*)simd_u64(aid, 1), KT, K)); \
			case 2: if (type &   4) ChargeLog(*(int64*)&bufN1[i + 2], bufN2[i + 2], SumProdDiv(b2, q2, (double*)simd_u64(aid, 2), KT, K)); \
			case 3: if (type &   8) ChargeLog(*(int64*)&bufN1[i + 3], bufN2[i + 3], SumProdDiv(b3, q3, (double*)simd_u64(aid, 3), KT, K)); \
			case 4: if (type &  16) ChargeLog(*(int64*)&bufN1[i + 4], bufN2[i + 4], SumProdDiv(b4, q4, (double*)simd_u64(aid, 4), KT, K)); \
			case 5: if (type &  32) ChargeLog(*(int64*)&bufN1[i + 5], bufN2[i + 5], SumProdDiv(b5, q5, (double*)simd_u64(aid, 5), KT, K)); \
			case 6: if (type &  64) ChargeLog(*(int64*)&bufN1[i + 6], bufN2[i + 6], SumProdDiv(b6, q6, (double*)simd_u64(aid, 6), KT, K)); \
			case 7: if (type & 128) ChargeLog(*(int64*)&bufN1[i + 7], bufN2[i + 7], SumProdDiv(b7, q7, (double*)simd_u64(aid, 7), KT, K)); \
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
TARGET512 void UpdateZAdmix512(double* Q, int K, int64* Mi, int* Ni, double* bufNK1, double* bufNK2, double* Base, RNG512& rng)
{
	__m512d maskone = _mm512_set1_pd(1.0);
	__m512i maskff = _mm512_set1_epi64(0xFFFFFFFFFFFFFFFF);
	__m512d* freq = AlignSIMD((__m512d*)bufNK1);
	__m512d* qq = AlignSIMD((__m512d*)bufNK2);
	__m512d dmin = _mm512_set1_pd(MIN_FREQ);
	__m512i vindex = _mm512_set_epi64(
		7 * structure_indnbytes, 6 * structure_indnbytes, 5 * structure_indnbytes, 4 * structure_indnbytes,
		3 * structure_indnbytes, 2 * structure_indnbytes, structure_indnbytes, 0);

	for (int i = 0, pad = 0; i < nind; i += STRUCTURE_NPACK)
	{
		if (i + STRUCTURE_NPACK >= nind) //fix i for the last several
		{
			pad = STRUCTURE_NPACK - (nind - i);
			i = nind - STRUCTURE_NPACK;
		}

		double* q0 = Q + i * K;
		int64* m0 = Mi + i * K, * m1 = m0 + K, * m2 = m1 + K, * m3 = m2 + K, * m4 = m3 + K, * m5 = m4 + K, * m6 = m5 + K, * m7 = m6 + K;
		__m512i vindex1 = _mm512_set_epi64(7 * K, 6 * K, 5 * K, 4 * K, 3 * K, 2 * K, 1 * K, 0 * K);
		
		for (int k = 0; k < K; ++k)
			_mm512_storeu_pd((double*)&qq[k],
				_mm512_i64gather_pd(vindex1, &q0[k], sizeof(double)));

		int* ni = Ni;
		double* p = Base;
		BAYESIAN_READER rt(i);

		for (int64 l = 0; l < nloc; ++l, ni += GetLoc(l).k, p += GetLoc(l).k)
		{
			int size = structure_size[l];
			__m512i mask = _mm512_set1_epi64((1u << size) - 1u);
			__m512i aid, zid;
			__m512d sumfreq;
			byte type;
			ushort k2;

#define REP_MID \
			ReadAid512/*(rt, size, aid, type);*/ \
			\
			sumfreq = _mm512_setzero_pd(); \
			\
			for (int k = 0; k < K; ++k) \
			{ \
				freq[k] = _mm512_add_pd( \
					_mm512_mul_pd(_mm512_loadu_pd((double*)&qq[k]), _mm512_mask_i64gather_pd( \
						maskone, type, aid, p + k * KT, sizeof(double))), dmin); \
			  \
				sumfreq = _mm512_add_pd(sumfreq, freq[k]); \
			} \
			/* draw z */ \
			rng.Poly(freq, sumfreq, K, zid); \
			\
			switch (pad) \
			{ /* update mi, ni */ \
			case 0: k2 = simd_u64(zid, 0); m0[k2] += (type &   1) != 0; ni[k2 * KT + simd_u64(aid, 0)] += (type &   1) != 0;  \
			case 1: k2 = simd_u64(zid, 1); m1[k2] += (type &   2) != 0; ni[k2 * KT + simd_u64(aid, 1)] += (type &   2) != 0;  \
			case 2: k2 = simd_u64(zid, 2); m2[k2] += (type &   4) != 0; ni[k2 * KT + simd_u64(aid, 2)] += (type &   4) != 0;  \
			case 3: k2 = simd_u64(zid, 3); m3[k2] += (type &   8) != 0; ni[k2 * KT + simd_u64(aid, 3)] += (type &   8) != 0;  \
			case 4: k2 = simd_u64(zid, 4); m4[k2] += (type &  16) != 0; ni[k2 * KT + simd_u64(aid, 4)] += (type &  16) != 0;  \
			case 5: k2 = simd_u64(zid, 5); m5[k2] += (type &  32) != 0; ni[k2 * KT + simd_u64(aid, 5)] += (type &  32) != 0;  \
			case 6: k2 = simd_u64(zid, 6); m6[k2] += (type &  64) != 0; ni[k2 * KT + simd_u64(aid, 6)] += (type &  64) != 0;  \
			case 7: k2 = simd_u64(zid, 7); m7[k2] += (type & 128) != 0; ni[k2 * KT + simd_u64(aid, 7)] += (type & 128) != 0;  \
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