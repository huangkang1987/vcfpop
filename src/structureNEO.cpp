/* NEO Instruction Set Functions */

#include "vcfpop.h"

#ifdef __aarch64__

/* Read allele id from decompreNEOd bucket */
/*ReadAidNEO(BAYESIAN_READER& rt, int size, uint64x2_t* aid, uint64x2_t* type)
{ \
	uint64x2_t mask = vdupq_n_u64((1u << size) - 1u);
	uint64x2_t aid1, aid2, aid3, aid4; 
	uint64x2_t type1, type2, type3, type4; */
#define ReadAidNEO \
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
		uint64x2_t maskrbits = vdupq_n_u64((1u << rbits) - 1u); \
		uint64x2_t c1, c2, c3, c4; \
 \
		c1 = vandq_u64(rt.data128[0], maskrbits); \
		c2 = vandq_u64(rt.data128[1], maskrbits); \
		c3 = vandq_u64(rt.data128[2], maskrbits); \
		c4 = vandq_u64(rt.data128[3], maskrbits); \
 \
		int64x2_t nlbits = vdupq_n_s64(rt.nbits); \
 \
		c1 = vshlq_u64(c1, nlbits); /*left move rt.nbits (positive) */ \
		c2 = vshlq_u64(c2, nlbits); \
		c3 = vshlq_u64(c3, nlbits); \
		c4 = vshlq_u64(c4, nlbits); \
 \
		aid1 = vorrq_u64(aid1, c1); \
		aid2 = vorrq_u64(aid2, c2); \
		aid3 = vorrq_u64(aid3, c3); \
		aid4 = vorrq_u64(aid4, c4); \
 \
		/*shift right */ \
		nlbits = vdupq_n_s64(-rbits); \
		rt.data128[0] = vshlq_u64(rt.data128[0], nlbits);/* right move rbits (negative) */ \
		rt.data128[1] = vshlq_u64(rt.data128[1], nlbits); \
		rt.data128[2] = vshlq_u64(rt.data128[2], nlbits); \
		rt.data128[3] = vshlq_u64(rt.data128[3], nlbits); \
 \
		rt.nbits = 64 - rbits; \
	} \
	else [[likely]] \
	{ \
		/*read size bits */ \
		aid1 = vandq_u64(rt.data128[0], mask); \
		aid2 = vandq_u64(rt.data128[1], mask); \
		aid3 = vandq_u64(rt.data128[2], mask); \
		aid4 = vandq_u64(rt.data128[3], mask); \
 \
		/*shift right */ \
		rt.data128[0] = vshlq_u64(rt.data128[0], nsize); /* right move size (negative) */ \
		rt.data128[1] = vshlq_u64(rt.data128[1], nsize); \
		rt.data128[2] = vshlq_u64(rt.data128[2], nsize); \
		rt.data128[3] = vshlq_u64(rt.data128[3], nsize); \
 \
		rt.nbits -= size; \
	} \
 \
	/*typed is 0xFFFFFFFF */ \
	type1 = vcgtq_u64(mask, aid1); \
	type2 = vcgtq_u64(mask, aid2); \
	type3 = vcgtq_u64(mask, aid3); \
	type4 = vcgtq_u64(mask, aid4); \
 \
	/* set missing to 0 */ \
	aid1 = vandq_u64(aid1, type1); \
	aid2 = vandq_u64(aid2, type2); \
	aid3 = vandq_u64(aid3, type3); \
	aid4 = vandq_u64(aid4, type4); \
}

/* Read allele frequency from decompreNEOd bucket */
/*
TARGETNEO void ReadFreqNEO(BAYESIAN_READER& rt, int size, uint64x2_t* slog, float64x2_t* prod, double* p, int K)
{
	int64x2_t nsize = vdupq_n_s64(-size);
	uint64x2_t mask = vdupq_n_u64((1u << size) - 1u);
	uint64x2_t maskunder = vdupq_n_u64(0x1FF0000000000000);
	uint64x2_t mask1 = vdupq_n_u64(0x7FF0000000000000);
	uint64x2_t mask2 = vdupq_n_u64(0x800FFFFFFFFFFFFF);
	uint64x2_t mask3 = vdupq_n_u64(0x3FF0000000000000);
	uint64x2_t subv = vdupq_n_u64(1023);
	uint64x2_t addr_inc = vdupq_n_u64(KT * sizeof(double));
	float64x2_t maskone = vdupq_n_f64(1.0);
*/
#define ReadFreqNEO \
{ \
	uint64x2_t aid1, aid2, aid3, aid4; \
	uint64x2_t a1, a2, a3, a4; \
	uint64x2_t b1, b2, b3, b4; \
	uint64x2_t type1, type2, type3, type4; \
	float64x2_t freq1, freq2, freq3, freq4; \
	uint64x2_t* prodi = (uint64x2_t*)prod; \
 \
	ReadAidNEO \
 \
	/* x sizeof(double) to load frequency */ \
	aid1 = vshlq_n_u64(aid1, 3); \
	aid2 = vshlq_n_u64(aid2, 3); \
	aid3 = vshlq_n_u64(aid3, 3); \
	aid4 = vshlq_n_u64(aid4, 3); \
 \
	uint64x2_t addr = vdupq_n_u64((uint64)p); \
 \
	aid1 = vaddq_f64(aid1, addr); \
	aid2 = vaddq_f64(aid2, addr); \
	aid3 = vaddq_f64(aid3, addr); \
	aid4 = vaddq_f64(aid4, addr); \
 \
	for (int k = 0; k < K; ++k) \
	{ \
		double f[8] = { \
			*(double*)vgetq_lane_u64(aid1, 0), *(double*)vgetq_lane_u64(aid1, 1), \
			*(double*)vgetq_lane_u64(aid2, 0), *(double*)vgetq_lane_u64(aid2, 1), \
			*(double*)vgetq_lane_u64(aid3, 0), *(double*)vgetq_lane_u64(aid3, 1), \
			*(double*)vgetq_lane_u64(aid4, 0), *(double*)vgetq_lane_u64(aid4, 1) \
		}; \
        \
		/* get freq of 8 allele copies for 8 individuals in cluster k */ \
		freq1 = vld1q_f64(&f[0]); \
		freq2 = vld1q_f64(&f[2]); \
		freq3 = vld1q_f64(&f[4]); \
		freq4 = vld1q_f64(&f[6]); \
 \
		/* set missing freq to one */ \
		freq1 = vbslq_u64(type1, freq1, maskone); \
		freq2 = vbslq_u64(type2, freq2, maskone); \
		freq3 = vbslq_u64(type3, freq3, maskone); \
		freq4 = vbslq_u64(type4, freq4, maskone); \
 \
		/* mul to prod */ \
		prod[k * 4 + 0] = vmulq_f64(prod[k * 4 + 0], freq1); \
		prod[k * 4 + 1] = vmulq_f64(prod[k * 4 + 1], freq2); \
		prod[k * 4 + 2] = vmulq_f64(prod[k * 4 + 2], freq3); \
		prod[k * 4 + 3] = vmulq_f64(prod[k * 4 + 3], freq4); \
 \
		aid1 = vaddq_f64(aid1, addr_inc); \
		aid2 = vaddq_f64(aid2, addr_inc); \
		aid3 = vaddq_f64(aid3, addr_inc); \
		aid4 = vaddq_f64(aid4, addr_inc); \
 \
		a1 = vcgtq_u64(maskunder, (uint64x2_t)prod[k * 4 + 0]); \
		a2 = vcgtq_u64(maskunder, (uint64x2_t)prod[k * 4 + 1]); \
		a3 = vcgtq_u64(maskunder, (uint64x2_t)prod[k * 4 + 2]); \
		a4 = vcgtq_u64(maskunder, (uint64x2_t)prod[k * 4 + 3]); \
 \
		a1 = vorrq_u64(a1, a2); \
		a3 = vorrq_u64(a3, a4); \
		a1 = vorrq_u64(a1, a3); \
 \
		/* check underflow */ \
		if (!(vgetq_lane_u64(a1, 0) | vgetq_lane_u64(a1, 1))) [[likely]] continue; \
 \
		/* add exponent */ \
		a1 = vandq_u64(prodi[k * 4 + 0], mask1); \
		a2 = vandq_u64(prodi[k * 4 + 1], mask1); \
		a3 = vandq_u64(prodi[k * 4 + 2], mask1); \
		a4 = vandq_u64(prodi[k * 4 + 3], mask1); \
 \
		b1 = vandq_u64(prodi[k * 4 + 0], mask2); \
		b2 = vandq_u64(prodi[k * 4 + 1], mask2); \
		b3 = vandq_u64(prodi[k * 4 + 2], mask2); \
		b4 = vandq_u64(prodi[k * 4 + 3], mask2); \
 \
		a1 = vshrq_n_u64(a1, 52); \
		a2 = vshrq_n_u64(a2, 52); \
		a3 = vshrq_n_u64(a3, 52); \
		a4 = vshrq_n_u64(a4, 52); \
 \
		prodi[k * 4 + 0] = vorrq_u64(b1, mask3); \
		prodi[k * 4 + 1] = vorrq_u64(b2, mask3); \
		prodi[k * 4 + 2] = vorrq_u64(b3, mask3); \
		prodi[k * 4 + 3] = vorrq_u64(b4, mask3); \
 \
		a1 = vsubq_u64(a1, subv); \
		a2 = vsubq_u64(a2, subv); \
		a3 = vsubq_u64(a3, subv); \
		a4 = vsubq_u64(a4, subv); \
 \
		slog[k * 4 + 0] = vaddq_f64(slog[k * 4 + 0], a1); \
		slog[k * 4 + 1] = vaddq_f64(slog[k * 4 + 1], a2); \
		slog[k * 4 + 2] = vaddq_f64(slog[k * 4 + 2], a3); \
		slog[k * 4 + 3] = vaddq_f64(slog[k * 4 + 3], a4); \
 \
	} \
}

/* Update individual or allele origin when ancetral proportion is binary */
TARGETNEO void UpdateZBinaryNEO(int* Ni, ushort* Z)
{
	uint64x2_t Kt = vdupq_n_u64(KT);
	uint64x2_t mask01 = vdupq_n_u64(1);

	for (int i = 0, pad = 0; i < nind; i += STRUCTURE_NPACK)
	{
		if (i + STRUCTURE_NPACK >= nind) //fix i for the last several
		{
			pad = STRUCTURE_NPACK - (nind - i);
			i = nind - STRUCTURE_NPACK;
		}

		uint64x2_t ZxKT1, ZxKT2, ZxKT3, ZxKT4;

		uint64 z[8] = {
			Z[i + 0], Z[i + 1], Z[i + 2], Z[i + 3],
			Z[i + 4], Z[i + 5], Z[i + 6], Z[i + 7]
		};

		ZxKT1 = vld1q_u64(&z[0]);
		ZxKT2 = vld1q_u64(&z[2]);
		ZxKT3 = vld1q_u64(&z[4]);
		ZxKT4 = vld1q_u64(&z[6]);

		ZxKT1 = (uint64x2_t)vmulq_u32((uint32x4_t)Kt, (uint32x4_t)ZxKT1);
		ZxKT2 = (uint64x2_t)vmulq_u32((uint32x4_t)Kt, (uint32x4_t)ZxKT2);
		ZxKT3 = (uint64x2_t)vmulq_u32((uint32x4_t)Kt, (uint32x4_t)ZxKT3);
		ZxKT4 = (uint64x2_t)vmulq_u32((uint32x4_t)Kt, (uint32x4_t)ZxKT4);

		ZxKT1 = vshlq_n_u64(ZxKT1, 2);
		ZxKT2 = vshlq_n_u64(ZxKT2, 2);
		ZxKT3 = vshlq_n_u64(ZxKT3, 2);
		ZxKT4 = vshlq_n_u64(ZxKT4, 2);

		BAYESIAN_READER rt(i);
		int* ni = Ni;

		for (int64 l = 0; l < nloc; ++l, ni += GetLoc(l).k)
		{
			int size = structure_size[l];
			int64x2_t nsize = vdupq_n_s64(-size);
			uint64x2_t mask = vdupq_n_u64((1u << size) - 1u); 
			uint64x2_t aid1, aid2, aid3, aid4;
			uint64x2_t type1, type2, type3, type4;
			uint64x2_t npio2 = vdupq_n_u64((int64)ni);
			uint64x2_t a1, a2, a3, a4;
			uint64x2_t b1, b2, b3, b4;

			a1 = vaddq_f64(ZxKT1, npio2);
			a2 = vaddq_f64(ZxKT2, npio2);
			a3 = vaddq_f64(ZxKT3, npio2);
			a4 = vaddq_f64(ZxKT4, npio2);

#define REP_MID \
			ReadAidNEO/*(rt, size, aid, type);*/ \
			\
			b1 = vshlq_n_u64(aid1, 2); \
			b2 = vshlq_n_u64(aid2, 2); \
			b3 = vshlq_n_u64(aid3, 2); \
			b4 = vshlq_n_u64(aid4, 2); \
			\
			b1 = vaddq_f64(a1, b1); \
			b2 = vaddq_f64(a2, b2); \
			b3 = vaddq_f64(a3, b3); \
			b4 = vaddq_f64(a4, b4); \
			\
			type1 = vandq_u64(type1, mask01);\
			type2 = vandq_u64(type2, mask01);\
			type3 = vandq_u64(type3, mask01);\
			type4 = vandq_u64(type4, mask01);\
			\
			switch (pad) \
			{ /* Z and Mi is updated before call, update Ni here */ \
			case 0: (*(int*)vgetq_lane_u64(b1, 0)) += vgetq_lane_u64(type1, 0); \
			case 1: (*(int*)vgetq_lane_u64(b1, 1)) += vgetq_lane_u64(type1, 1); \
			case 2: (*(int*)vgetq_lane_u64(b2, 0)) += vgetq_lane_u64(type2, 0); \
			case 3: (*(int*)vgetq_lane_u64(b2, 1)) += vgetq_lane_u64(type2, 1); \
			case 4: (*(int*)vgetq_lane_u64(b3, 0)) += vgetq_lane_u64(type3, 0); \
			case 5: (*(int*)vgetq_lane_u64(b3, 1)) += vgetq_lane_u64(type3, 1); \
			case 6: (*(int*)vgetq_lane_u64(b4, 0)) += vgetq_lane_u64(type4, 0); \
			case 7: (*(int*)vgetq_lane_u64(b4, 1)) += vgetq_lane_u64(type4, 1); \
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
TARGETNEO void UpdateQNoAdmixNEO(double* Q, int K, double* bufNK1, double* bufNK2, double* Gamma, RNGNEO& rng, bool locpriori, double* Base, ushort* Z)
{
	uint64x2_t maskunder = vdupq_n_u64(0x1FF0000000000000);
	uint64x2_t mask1 = vdupq_n_u64(0x7FF0000000000000);
	uint64x2_t mask2 = vdupq_n_u64(0x800FFFFFFFFFFFFF);
	uint64x2_t mask3 = vdupq_n_u64(0x3FF0000000000000);
	uint64x2_t subv = vdupq_n_u64(1023);
	uint64x2_t addr_inc = vdupq_n_u64(KT * sizeof(double));
	float64x2_t maskone = vdupq_n_f64(1.0);

	double* q = Q;

	for (int i = 0, pad = 0; i < nind; i += STRUCTURE_NPACK)
	{
		if (i + STRUCTURE_NPACK >= nind) //fix i for the last several
		{
			pad = STRUCTURE_NPACK - (nind - i);
			i = nind - STRUCTURE_NPACK;
		}

		//K elements, each save NPACK log or likelihood
		uint64x2_t* slog = AlignSIMD((uint64x2_t*)bufNK1);
		float64x2_t* prod = AlignSIMD((float64x2_t*)bufNK2);

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
			int64x2_t nsize = vdupq_n_s64(-size);
			uint64x2_t mask = vdupq_n_u64((1u << size) - 1u);

			switch (maxploidy)
			{
			case 10: ReadFreqNEO/*(rt, size, slog, prod, p, K)*/;
			case  9: ReadFreqNEO/*(rt, size, slog, prod, p, K)*/;
			case  8: ReadFreqNEO/*(rt, size, slog, prod, p, K)*/;
			case  7: ReadFreqNEO/*(rt, size, slog, prod, p, K)*/;
			case  6: ReadFreqNEO/*(rt, size, slog, prod, p, K)*/;
			case  5: ReadFreqNEO/*(rt, size, slog, prod, p, K)*/;
			case  4: ReadFreqNEO/*(rt, size, slog, prod, p, K)*/;
			case  3: ReadFreqNEO/*(rt, size, slog, prod, p, K)*/;
			case  2: ReadFreqNEO/*(rt, size, slog, prod, p, K)*/;
			case  1: ReadFreqNEO/*(rt, size, slog, prod, p, K)*/;
			}
		}

		CloseLog((int64*)slog, (double*)prod, K * STRUCTURE_NPACK);

		uint64x2_t ki[4];
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
TARGETNEO void UpdateQMetroNEO(double* Q, int K, double* bufNK1, double* bufN1, double* bufN2, double* Base)
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
			int64x2_t nsize = vdupq_n_s64(-size);
			uint64x2_t mask = vdupq_n_u64((1u << size) - 1u);
			uint64x2_t aid1, aid2, aid3, aid4;
			uint64x2_t type1, type2, type3, type4;
			uint64x2_t paddr = vdupq_n_u64((uint64)p);

#define REP_MID \
			ReadAidNEO/*(rt, size, aid, type)*/; \
			\
			aid1 = vshlq_n_u64(aid1, 3); \
			aid2 = vshlq_n_u64(aid2, 3); \
			aid3 = vshlq_n_u64(aid3, 3); \
			aid4 = vshlq_n_u64(aid4, 3); \
			\
			aid1 = vaddq_f64(paddr, aid1); \
			aid2 = vaddq_f64(paddr, aid2); \
			aid3 = vaddq_f64(paddr, aid3); \
			aid4 = vaddq_f64(paddr, aid4); \
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
TARGETNEO void UpdateZAdmixNEO(double* Q, int K, int64* Mi, int* Ni, double* bufNK1, double* bufNK2, double* Base, RNGNEO& rng)
{
	uint64x2_t mask01 = vdupq_n_u64(1);
	float64x2_t* freq = AlignSIMD((float64x2_t*)bufNK1);
	float64x2_t* qq = AlignSIMD((float64x2_t*)bufNK2);
	float64x2_t dmin = vdupq_n_f64(MIN_FREQ);
	uint64x2_t addr_inc = vdupq_n_u64(KT * sizeof(double));

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
			double q[8] = { q0[k], q1[k], q2[k], q3[k], 
							q4[k], q5[k], q6[k], q7[k] };

			qq[k * 4 + 0] = vld1q_f64(&q[0]);
			qq[k * 4 + 1] = vld1q_f64(&q[2]);
			qq[k * 4 + 2] = vld1q_f64(&q[4]);
			qq[k * 4 + 3] = vld1q_f64(&q[6]);
		}

		int* ni = Ni;
		double* p = Base;
		BAYESIAN_READER rt(i);

		for (int64 l = 0; l < nloc; ++l, ni += GetLoc(l).k, p += GetLoc(l).k)
		{
			int size = structure_size[l];
			int64x2_t nsize = vdupq_n_s64(-size);
			uint64x2_t mask = vdupq_n_u64((1u << size) - 1u);
			uint64x2_t type[4];
			uint64x2_t aid1, aid2, aid3, aid4;
			uint64x2_t& type1 = type[0];
			uint64x2_t& type2 = type[1];
			uint64x2_t& type3 = type[2];
			uint64x2_t& type4 = type[3];
			uint64x2_t t1, t2, t3, t4;
			uint64x2_t addr = vdupq_n_u64((uint64)p);
			float64x2_t sumfreq[4];
			float64x2_t a1, a2, a3, a4;
			uint64x2_t zid[4];
			uint64x2_t& zid1 = zid[0];
			uint64x2_t& zid2 = zid[1];
			uint64x2_t& zid3 = zid[2];
			uint64x2_t& zid4 = zid[3];
			ushort k2;

#define REP_MID \
			ReadAidNEO/*(rt, size, aid, type)*/; \
			\
			sumfreq[0] = vdupq_n_f64(0); \
			sumfreq[1] = vdupq_n_f64(0); \
			sumfreq[2] = vdupq_n_f64(0); \
			sumfreq[3] = vdupq_n_f64(0); \
			\
			t1 = vshlq_n_u64(aid1, 3); \
			t2 = vshlq_n_u64(aid2, 3); \
			t3 = vshlq_n_u64(aid3, 3); \
			t4 = vshlq_n_u64(aid4, 3); \
			\
			t1 = vaddq_f64(t1, addr); \
			t2 = vaddq_f64(t2, addr); \
			t3 = vaddq_f64(t3, addr); \
			t4 = vaddq_f64(t4, addr); \
			\
			for (int k = 0; k < K; ++k) \
			{ \
				double t[8] = { \
					*(double*)vgetq_lane_u64(t1, 0), *(double*)vgetq_lane_u64(t1, 1), \
					*(double*)vgetq_lane_u64(t2, 0), *(double*)vgetq_lane_u64(t2, 1), \
					*(double*)vgetq_lane_u64(t3, 0), *(double*)vgetq_lane_u64(t3, 1), \
					*(double*)vgetq_lane_u64(t4, 0), *(double*)vgetq_lane_u64(t4, 1)  \
				}; \
				\
				a1 = vld1q_f64(&t[0]); \
				a2 = vld1q_f64(&t[2]); \
				a3 = vld1q_f64(&t[4]); \
				a4 = vld1q_f64(&t[6]); \
				\
				a1 = vmulq_f64(qq[k * 4 + 0], a1); \
				a2 = vmulq_f64(qq[k * 4 + 1], a2); \
				a3 = vmulq_f64(qq[k * 4 + 2], a3); \
				a4 = vmulq_f64(qq[k * 4 + 3], a4); \
				\
				freq[k * 4 + 0] = vaddq_f64(a1, dmin); \
				freq[k * 4 + 1] = vaddq_f64(a2, dmin); \
				freq[k * 4 + 2] = vaddq_f64(a3, dmin); \
				freq[k * 4 + 3] = vaddq_f64(a4, dmin); \
				\
				sumfreq[0] = vaddq_f64(sumfreq[0], freq[k * 4 + 0]); \
				sumfreq[1] = vaddq_f64(sumfreq[1], freq[k * 4 + 1]); \
				sumfreq[2] = vaddq_f64(sumfreq[2], freq[k * 4 + 2]); \
				sumfreq[3] = vaddq_f64(sumfreq[3], freq[k * 4 + 3]); \
				\
				t1 = vaddq_f64(t1, addr_inc); \
				t2 = vaddq_f64(t2, addr_inc); \
				t3 = vaddq_f64(t3, addr_inc); \
				t4 = vaddq_f64(t4, addr_inc); \
			} \
			/* draw z */ \
			rng.Poly(freq, sumfreq, K, zid); \
			\
			type1 = vandq_u64(type1, mask01);\
			type2 = vandq_u64(type2, mask01);\
			type3 = vandq_u64(type3, mask01);\
			type4 = vandq_u64(type4, mask01);\
			\
			switch (pad) \
			{ /* update mi, ni */ \
			case 0: k2 = vgetq_lane_u64(zid1, 0); m0[k2] += vgetq_lane_u64(type1, 0); ni[k2 * KT + vgetq_lane_u64(aid1, 0)] += vgetq_lane_u64(type1, 0); \
			case 1: k2 = vgetq_lane_u64(zid1, 1); m1[k2] += vgetq_lane_u64(type1, 1); ni[k2 * KT + vgetq_lane_u64(aid1, 1)] += vgetq_lane_u64(type1, 1); \
			case 2: k2 = vgetq_lane_u64(zid2, 0); m2[k2] += vgetq_lane_u64(type2, 0); ni[k2 * KT + vgetq_lane_u64(aid2, 0)] += vgetq_lane_u64(type2, 0); \
			case 3: k2 = vgetq_lane_u64(zid2, 1); m3[k2] += vgetq_lane_u64(type2, 1); ni[k2 * KT + vgetq_lane_u64(aid2, 1)] += vgetq_lane_u64(type2, 1); \
			case 4: k2 = vgetq_lane_u64(zid3, 0); m4[k2] += vgetq_lane_u64(type3, 0); ni[k2 * KT + vgetq_lane_u64(aid3, 0)] += vgetq_lane_u64(type3, 0); \
			case 5: k2 = vgetq_lane_u64(zid3, 1); m5[k2] += vgetq_lane_u64(type3, 1); ni[k2 * KT + vgetq_lane_u64(aid3, 1)] += vgetq_lane_u64(type3, 1); \
			case 6: k2 = vgetq_lane_u64(zid4, 0); m6[k2] += vgetq_lane_u64(type4, 0); ni[k2 * KT + vgetq_lane_u64(aid4, 0)] += vgetq_lane_u64(type4, 0); \
			case 7: k2 = vgetq_lane_u64(zid4, 1); m7[k2] += vgetq_lane_u64(type4, 1); ni[k2 * KT + vgetq_lane_u64(aid4, 1)] += vgetq_lane_u64(type4, 1); \
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
