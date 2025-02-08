/* NEO Bayesian clustering functions */

#include "vcfpop.h"

template struct BAYESIAN<double>;
template struct BAYESIAN<float >;

template TARGETNEO void BAYESIAN<double>::UpdateQNoAdmixNEO(int tid);
template TARGETNEO void BAYESIAN<float >::UpdateQNoAdmixNEO(int tid);

template TARGETNEO void BAYESIAN<double>::UpdateZNoAdmixNEO(int tid);
template TARGETNEO void BAYESIAN<float >::UpdateZNoAdmixNEO(int tid);

template TARGETNEO void BAYESIAN<double>::UpdateQMetroNEO<true >(int tid);
template TARGETNEO void BAYESIAN<double>::UpdateQMetroNEO<false>(int tid);
template TARGETNEO void BAYESIAN<float >::UpdateQMetroNEO<true >(int tid);
template TARGETNEO void BAYESIAN<float >::UpdateQMetroNEO<false>(int tid);

template TARGETNEO void BAYESIAN<double>::UpdateZAdmixNEO<true >(int tid);
template TARGETNEO void BAYESIAN<double>::UpdateZAdmixNEO<false>(int tid);
template TARGETNEO void BAYESIAN<float >::UpdateZAdmixNEO<true >(int tid);
template TARGETNEO void BAYESIAN<float >::UpdateZAdmixNEO<false>(int tid);

template TARGETNEO void BAYESIAN<double>::RecordNEO<true , true >(int tid);
template TARGETNEO void BAYESIAN<double>::RecordNEO<true , false>(int tid);
template TARGETNEO void BAYESIAN<double>::RecordNEO<false, true >(int tid);
template TARGETNEO void BAYESIAN<double>::RecordNEO<false, false>(int tid);
template TARGETNEO void BAYESIAN<float >::RecordNEO<true , true >(int tid);
template TARGETNEO void BAYESIAN<float >::RecordNEO<true , false>(int tid);
template TARGETNEO void BAYESIAN<float >::RecordNEO<false, true >(int tid);
template TARGETNEO void BAYESIAN<float >::RecordNEO<false, false>(int tid);

#ifdef __aarch64__

#ifndef _GENO_READERNEO
template<typename REAL>
struct GENO_READERNEO
{
	uint32x4_t data[16];					//Readed bits
	uint32x4_t msize;
	byte* pos;								//Current read pointer
	uint vindex[64];						//Offset of 16 loci
	byte size;								//Number of bits a genotype id used
	byte nbits;								//Number of bits remaining in data

	TARGETNEO GENO_READERNEO()
	{

	}

	TARGETNEO GENO_READERNEO(int indid, int64 l, int64 num, BUCKET* bucket = NULL)
	{
		//set pos and size
		SetZero(this, 1);
		num = std::min(num, 64LL);
            
		//set bucket from default bucket or assigned bucket
		if (bucket == NULL) bucket = &geno_bucket;

		OFFSET* offset = &bucket->offset.bucket[l];
		uint64 offset0 = offset[0].offset;
		pos = (byte*)(bucket->base_addr + offset0);

		for (int i = 0; i < num; ++i)
			vindex[i] = offset[i].offset - offset0;

		size = offset[0].size;

		msize = vdupq_n_u32((1u << size) - 1u);

		if (indid)
		{
			//skip indid * size bits
			int nreadbits = indid * size;

			//pos move nreadbits / 32
			pos += nreadbits >> 5;

			//remain bits to read
			nreadbits &= 31;

			//read 32 bits to data and read remain bits from data
			uint* data2 = (uint*)data;
			REP(64) data2[kk] = *(uint*)(pos + vindex[kk]);

			switch (nreadbits)
			{
			case  0:                                               break;
			case  1: REP(16) data[kk] = vshrq_n_u32(data[kk],  1); break;
			case  2: REP(16) data[kk] = vshrq_n_u32(data[kk],  2); break;
			case  3: REP(16) data[kk] = vshrq_n_u32(data[kk],  3); break;
			case  4: REP(16) data[kk] = vshrq_n_u32(data[kk],  4); break;
			case  5: REP(16) data[kk] = vshrq_n_u32(data[kk],  5); break;
			case  6: REP(16) data[kk] = vshrq_n_u32(data[kk],  6); break;
			case  7: REP(16) data[kk] = vshrq_n_u32(data[kk],  7); break;
			case  8: REP(16) data[kk] = vshrq_n_u32(data[kk],  8); break;
			case  9: REP(16) data[kk] = vshrq_n_u32(data[kk],  9); break;
			case 10: REP(16) data[kk] = vshrq_n_u32(data[kk], 10); break;
			case 11: REP(16) data[kk] = vshrq_n_u32(data[kk], 11); break;
			case 12: REP(16) data[kk] = vshrq_n_u32(data[kk], 12); break;
			case 13: REP(16) data[kk] = vshrq_n_u32(data[kk], 13); break;
			case 14: REP(16) data[kk] = vshrq_n_u32(data[kk], 14); break;
			case 15: REP(16) data[kk] = vshrq_n_u32(data[kk], 15); break;
			case 16: REP(16) data[kk] = vshrq_n_u32(data[kk], 16); break;
			case 17: REP(16) data[kk] = vshrq_n_u32(data[kk], 17); break;
			case 18: REP(16) data[kk] = vshrq_n_u32(data[kk], 18); break;
			case 19: REP(16) data[kk] = vshrq_n_u32(data[kk], 19); break;
			case 20: REP(16) data[kk] = vshrq_n_u32(data[kk], 20); break;
			case 21: REP(16) data[kk] = vshrq_n_u32(data[kk], 21); break;
			case 22: REP(16) data[kk] = vshrq_n_u32(data[kk], 22); break;
			case 23: REP(16) data[kk] = vshrq_n_u32(data[kk], 23); break;
			case 24: REP(16) data[kk] = vshrq_n_u32(data[kk], 24); break;
			case 25: REP(16) data[kk] = vshrq_n_u32(data[kk], 25); break;
			case 26: REP(16) data[kk] = vshrq_n_u32(data[kk], 26); break;
			case 27: REP(16) data[kk] = vshrq_n_u32(data[kk], 27); break;
			case 28: REP(16) data[kk] = vshrq_n_u32(data[kk], 28); break;
			case 29: REP(16) data[kk] = vshrq_n_u32(data[kk], 29); break;
			case 30: REP(16) data[kk] = vshrq_n_u32(data[kk], 30); break;
			case 31: REP(16) data[kk] = vshrq_n_u32(data[kk], 31); break;
			}

			pos += 4;

			//set nbits
			nbits = 32 - nreadbits;
		}
		else
		{
			//read 32 bits
			uint* data2 = (uint*)data;
			REP(64) data2[kk] = *(uint*)(pos + vindex[kk]);

			pos += 4;

			//set nbits
			nbits = 32;
		}
	}

	forceinline TARGETNEO void Read(uint32x4_t* gid)
	{
		// if data is empty
		if (nbits < size) [[unlikely]]
		{
			// move nbits data to gid
			memcpy(gid, data, 16 * sizeof(uint32x4_t));

			// remain number of bits to read
			int rbits = size - nbits;
			uint32x4_t tmask = vdupq_n_u32((1u << rbits) - 1u);

			// read 32 bits to data
			uint* data2 = (uint*)data;
			REP(64) data2[kk] = *(uint*)(pos + vindex[kk]);

			// read rbits from data and concate to higher bits in gid
			switch (nbits)
			{
			case  0: REP(16) gid[kk] = veorq_u32(gid[kk],             vandq_u32(data[kk], tmask)     ); break;
			case  1: REP(16) gid[kk] = veorq_u32(gid[kk], vshlq_n_u32(vandq_u32(data[kk], tmask),  1)); break;
			case  2: REP(16) gid[kk] = veorq_u32(gid[kk], vshlq_n_u32(vandq_u32(data[kk], tmask),  2)); break;
			case  3: REP(16) gid[kk] = veorq_u32(gid[kk], vshlq_n_u32(vandq_u32(data[kk], tmask),  3)); break;
			case  4: REP(16) gid[kk] = veorq_u32(gid[kk], vshlq_n_u32(vandq_u32(data[kk], tmask),  4)); break;
			case  5: REP(16) gid[kk] = veorq_u32(gid[kk], vshlq_n_u32(vandq_u32(data[kk], tmask),  5)); break;
			case  6: REP(16) gid[kk] = veorq_u32(gid[kk], vshlq_n_u32(vandq_u32(data[kk], tmask),  6)); break;
			case  7: REP(16) gid[kk] = veorq_u32(gid[kk], vshlq_n_u32(vandq_u32(data[kk], tmask),  7)); break;
			case  8: REP(16) gid[kk] = veorq_u32(gid[kk], vshlq_n_u32(vandq_u32(data[kk], tmask),  8)); break;
			case  9: REP(16) gid[kk] = veorq_u32(gid[kk], vshlq_n_u32(vandq_u32(data[kk], tmask),  9)); break;
			case 10: REP(16) gid[kk] = veorq_u32(gid[kk], vshlq_n_u32(vandq_u32(data[kk], tmask), 10)); break;
			case 11: REP(16) gid[kk] = veorq_u32(gid[kk], vshlq_n_u32(vandq_u32(data[kk], tmask), 11)); break;
			case 12: REP(16) gid[kk] = veorq_u32(gid[kk], vshlq_n_u32(vandq_u32(data[kk], tmask), 12)); break;
			case 13: REP(16) gid[kk] = veorq_u32(gid[kk], vshlq_n_u32(vandq_u32(data[kk], tmask), 13)); break;
			case 14: REP(16) gid[kk] = veorq_u32(gid[kk], vshlq_n_u32(vandq_u32(data[kk], tmask), 14)); break;
			case 15: REP(16) gid[kk] = veorq_u32(gid[kk], vshlq_n_u32(vandq_u32(data[kk], tmask), 15)); break;
			case 16: REP(16) gid[kk] = veorq_u32(gid[kk], vshlq_n_u32(vandq_u32(data[kk], tmask), 16)); break;
			case 17: REP(16) gid[kk] = veorq_u32(gid[kk], vshlq_n_u32(vandq_u32(data[kk], tmask), 17)); break;
			case 18: REP(16) gid[kk] = veorq_u32(gid[kk], vshlq_n_u32(vandq_u32(data[kk], tmask), 18)); break;
			case 19: REP(16) gid[kk] = veorq_u32(gid[kk], vshlq_n_u32(vandq_u32(data[kk], tmask), 19)); break;
			case 20: REP(16) gid[kk] = veorq_u32(gid[kk], vshlq_n_u32(vandq_u32(data[kk], tmask), 20)); break;
			case 21: REP(16) gid[kk] = veorq_u32(gid[kk], vshlq_n_u32(vandq_u32(data[kk], tmask), 21)); break;
			case 22: REP(16) gid[kk] = veorq_u32(gid[kk], vshlq_n_u32(vandq_u32(data[kk], tmask), 22)); break;
			case 23: REP(16) gid[kk] = veorq_u32(gid[kk], vshlq_n_u32(vandq_u32(data[kk], tmask), 23)); break;
			case 24: REP(16) gid[kk] = veorq_u32(gid[kk], vshlq_n_u32(vandq_u32(data[kk], tmask), 24)); break;
			case 25: REP(16) gid[kk] = veorq_u32(gid[kk], vshlq_n_u32(vandq_u32(data[kk], tmask), 25)); break;
			case 26: REP(16) gid[kk] = veorq_u32(gid[kk], vshlq_n_u32(vandq_u32(data[kk], tmask), 26)); break;
			case 27: REP(16) gid[kk] = veorq_u32(gid[kk], vshlq_n_u32(vandq_u32(data[kk], tmask), 27)); break;
			case 28: REP(16) gid[kk] = veorq_u32(gid[kk], vshlq_n_u32(vandq_u32(data[kk], tmask), 28)); break;
			case 29: REP(16) gid[kk] = veorq_u32(gid[kk], vshlq_n_u32(vandq_u32(data[kk], tmask), 29)); break;
			case 30: REP(16) gid[kk] = veorq_u32(gid[kk], vshlq_n_u32(vandq_u32(data[kk], tmask), 30)); break;
			case 31: REP(16) gid[kk] = veorq_u32(gid[kk], vshlq_n_u32(vandq_u32(data[kk], tmask), 31)); break;
			}




			//shift right
			switch (rbits)
			{
			case  0:                                               break;
			case  1: REP(16) data[kk] = vshrq_n_u32(data[kk],  1); break;
			case  2: REP(16) data[kk] = vshrq_n_u32(data[kk],  2); break;
			case  3: REP(16) data[kk] = vshrq_n_u32(data[kk],  3); break;
			case  4: REP(16) data[kk] = vshrq_n_u32(data[kk],  4); break;
			case  5: REP(16) data[kk] = vshrq_n_u32(data[kk],  5); break;
			case  6: REP(16) data[kk] = vshrq_n_u32(data[kk],  6); break;
			case  7: REP(16) data[kk] = vshrq_n_u32(data[kk],  7); break;
			case  8: REP(16) data[kk] = vshrq_n_u32(data[kk],  8); break;
			case  9: REP(16) data[kk] = vshrq_n_u32(data[kk],  9); break;
			case 10: REP(16) data[kk] = vshrq_n_u32(data[kk], 10); break;
			case 11: REP(16) data[kk] = vshrq_n_u32(data[kk], 11); break;
			case 12: REP(16) data[kk] = vshrq_n_u32(data[kk], 12); break;
			case 13: REP(16) data[kk] = vshrq_n_u32(data[kk], 13); break;
			case 14: REP(16) data[kk] = vshrq_n_u32(data[kk], 14); break;
			case 15: REP(16) data[kk] = vshrq_n_u32(data[kk], 15); break;
			case 16: REP(16) data[kk] = vshrq_n_u32(data[kk], 16); break;
			case 17: REP(16) data[kk] = vshrq_n_u32(data[kk], 17); break;
			case 18: REP(16) data[kk] = vshrq_n_u32(data[kk], 18); break;
			case 19: REP(16) data[kk] = vshrq_n_u32(data[kk], 19); break;
			case 20: REP(16) data[kk] = vshrq_n_u32(data[kk], 20); break;
			case 21: REP(16) data[kk] = vshrq_n_u32(data[kk], 21); break;
			case 22: REP(16) data[kk] = vshrq_n_u32(data[kk], 22); break;
			case 23: REP(16) data[kk] = vshrq_n_u32(data[kk], 23); break;
			case 24: REP(16) data[kk] = vshrq_n_u32(data[kk], 24); break;
			case 25: REP(16) data[kk] = vshrq_n_u32(data[kk], 25); break;
			case 26: REP(16) data[kk] = vshrq_n_u32(data[kk], 26); break;
			case 27: REP(16) data[kk] = vshrq_n_u32(data[kk], 27); break;
			case 28: REP(16) data[kk] = vshrq_n_u32(data[kk], 28); break;
			case 29: REP(16) data[kk] = vshrq_n_u32(data[kk], 29); break;
			case 30: REP(16) data[kk] = vshrq_n_u32(data[kk], 30); break;
			case 31: REP(16) data[kk] = vshrq_n_u32(data[kk], 31); break;
			}

			pos += 4;

			nbits = 32 - rbits;
		}
		else [[likely]]
		{
			//read size bits
			REP(16) gid[kk] = vandq_u32(data[kk], msize);

			//shift right
			switch (size)
			{
			case  0:                                              break;
			case  1: REP(16) data[kk] = vshrq_n_u32(data[kk], 1); break;
			case  2: REP(16) data[kk] = vshrq_n_u32(data[kk], 2); break;
			case  3: REP(16) data[kk] = vshrq_n_u32(data[kk], 3); break;
			case  4: REP(16) data[kk] = vshrq_n_u32(data[kk], 4); break;
			case  5: REP(16) data[kk] = vshrq_n_u32(data[kk], 5); break;
			case  6: REP(16) data[kk] = vshrq_n_u32(data[kk], 6); break;
			case  7: REP(16) data[kk] = vshrq_n_u32(data[kk], 7); break;
			case  8: REP(16) data[kk] = vshrq_n_u32(data[kk], 8); break;
			case  9: REP(16) data[kk] = vshrq_n_u32(data[kk], 9); break;
			case 10: REP(16) data[kk] = vshrq_n_u32(data[kk], 10); break;
			case 11: REP(16) data[kk] = vshrq_n_u32(data[kk], 11); break;
			case 12: REP(16) data[kk] = vshrq_n_u32(data[kk], 12); break;
			case 13: REP(16) data[kk] = vshrq_n_u32(data[kk], 13); break;
			case 14: REP(16) data[kk] = vshrq_n_u32(data[kk], 14); break;
			case 15: REP(16) data[kk] = vshrq_n_u32(data[kk], 15); break;
			case 16: REP(16) data[kk] = vshrq_n_u32(data[kk], 16); break;
			case 17: REP(16) data[kk] = vshrq_n_u32(data[kk], 17); break;
			case 18: REP(16) data[kk] = vshrq_n_u32(data[kk], 18); break;
			case 19: REP(16) data[kk] = vshrq_n_u32(data[kk], 19); break;
			case 20: REP(16) data[kk] = vshrq_n_u32(data[kk], 20); break;
			case 21: REP(16) data[kk] = vshrq_n_u32(data[kk], 21); break;
			case 22: REP(16) data[kk] = vshrq_n_u32(data[kk], 22); break;
			case 23: REP(16) data[kk] = vshrq_n_u32(data[kk], 23); break;
			case 24: REP(16) data[kk] = vshrq_n_u32(data[kk], 24); break;
			case 25: REP(16) data[kk] = vshrq_n_u32(data[kk], 25); break;
			case 26: REP(16) data[kk] = vshrq_n_u32(data[kk], 26); break;
			case 27: REP(16) data[kk] = vshrq_n_u32(data[kk], 27); break;
			case 28: REP(16) data[kk] = vshrq_n_u32(data[kk], 28); break;
			case 29: REP(16) data[kk] = vshrq_n_u32(data[kk], 29); break;
			case 30: REP(16) data[kk] = vshrq_n_u32(data[kk], 30); break;
			case 31: REP(16) data[kk] = vshrq_n_u32(data[kk], 31); break;
			}

			nbits -= size;
		}
	}
};

#endif

/* Update a priori ancetral proportion for non-admix model */
template<typename REAL>
TARGETNEO void BAYESIAN<REAL>::UpdateQNoAdmixNEO(int tid)
{
	if (tid == -1)
	{
		SetZero(Q, N * K);
		OpenLog((int64*)bufNK1, bufNK2, N * K * structure_nsubthread);

		//add priori probability
		double* buf1 = bufNK1, * buf2 = bufNK2;
		if (locpriori) for (int i = 0; i < N; ++i, buf1 += K, buf2 += K)
		{
			if (ainds<REAL>[i]->vt == 0) continue;
			ChargeLog((int64*)buf1, buf2, Gamma + ainds<REAL>[i]->popid * K, K);
		}

		//////////////////////////////////////////////////////////

		SetZero((int64*)l_atomic, 32);

		UpdateQNoAdmixNEO(0);

		//avoid thread-conflict
		for (int i = 1; i < structure_nsubthread; ++i)
			Add(bufNK1, bufNK1 + N * K * i, N * K);

		buf1 = bufNK1;
		REAL* q = Q;
		RNG<double> rng(seed + m, RNG_SALT_UPDATEQ);//double
		for (int i = 0; i < N; ++i, buf1 += K, q += K)
		{
			if (ainds<REAL>[i]->vt == 0) continue;
			ushort k2 = (ushort)rng.PolyLog(buf1, K);
			q[k2] = 1;
			Z[i] = k2;
		}
		return;
	}

	static float64x2_t maskoned = vdupq_n_f64(1.0);
	static float32x4_t maskones = vdupq_n_f32(1.0);
	static uint32x4_t mask00 = vdupq_n_u32(0);
	static uint32x4_t maskff = vdupq_n_u32(0xFFFFFFFF);
	static uint32x4_t mask24 = vdupq_n_u32(0xFFFFFF);
	static uint32x4_t mask01 = vdupq_n_u32(1);
	static uint32x4_t mask02 = vdupq_n_u32(2);
	alignas(16) static uint maskidx1[64] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63 };
	static uint32x4_t* maskidx = (uint32x4_t*)maskidx1;
	static uint PT_PLOIDYxNALLELES[150] = 									//Pattern index to ploidy level
	{ 0, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	atomic<int> thread_counter = 0;
#pragma omp parallel num_threads(structure_nsubthread)
	{
		int tid2 = thread_counter.fetch_add(1);

		for (int lsize = structure_loc_size_min; lsize <= structure_loc_size_max; ++lsize)
		{
			int64 lstart = structure_loc_lend[lsize - 1], lend = structure_loc_lend[lsize];
			int64 lend0 = (lend - lstart + 63) / 64;

			for (int64 ia = l_atomic[lsize].fetch_add(1); ia < lend0; ia = l_atomic[lsize].fetch_add(1))
			{
				int64 l = lstart + (ia * structure_loc_coprime64[lsize] % lend0 << 6);
				REAL* p = Freq + allele_freq_offset[l];
				int64 gtab_base = (int64)GetLoc(l).GetGtab();

				uint32x4_t gtab[16]; uint* gtab1 = (uint*)gtab;
				REP(64) gtab1[kk] = GetLocTabDiff(l + kk);

				uint32x4_t lmask[16], lmaskidx[16];
				uint64 lmask0 = lend - l >= 64 ? 0xFFFFFFFFFFFFFFFF : (1ull << (lend - l)) - 1ull;
				uint32x4_t lmaskt = vdupq_n_u32(lmask0);
				uint lmaskv1[4] = { 1, 2, 4, 8 }; uint32x4_t lmaskv2 = vld1q_u32(lmaskv1);
				REP(16)
				{
                    if (kk == 8) lmaskt = vdupq_n_u32(lmask0 >> 32);
					lmask[kk] = vcgtq_u32(vandq_u32(lmaskt, lmaskv2), mask00);
					lmaskt = vshrq_n_u32(lmaskt, 4);
					lmaskidx[kk] = vandq_u32(lmask[kk], maskidx[kk]);
				}

				GENO_READERNEO<REAL> rt(0, l, lend - l);

				uint32x4_t oindex[16];
				uint64* toffset = allele_freq_offset + l;
				uint64 oindex01 = toffset[0];
				int* lmaskidx2 = (int*)lmaskidx;

				REP(16) oindex[kk] = vld1q_u32(((uint[]) { (uint)(toffset[lmaskidx2[0 + (kk << 2)]] - oindex01), (uint)(toffset[lmaskidx2[1 + (kk << 2)]] - oindex01), (uint)(toffset[lmaskidx2[2 + (kk << 2)]] - oindex01), (uint)(toffset[lmaskidx2[3 + (kk << 2)]] - oindex01) }));

                uint32x4_t gtaddr[16], gtlow[16], gtploidy[16], gtals[16], allele[16], typed[16];
                uint64x2_t typed64[32];
				int* allele2 = (int*)allele, * gtploidy2 = (int*)gtploidy, * gtaddr2 = (int*)gtaddr;
				float64x2_t freq[32];

				int64* buf1 = (int64*)bufNK1 + N * K * tid2; double* buf2 = bufNK2 + N * K * tid2;

				for (int i = 0; i < N; ++i, buf1 += K, buf2 += K)
				{
					rt.Read(gtaddr);

					REP(16) gtaddr[kk] = vaddq_u32(gtab[kk], vshlq_n_u32(gtaddr[kk], 2));

					REP(16) gtlow[kk] = vld1q_u32(((uint[]) { *(uint*)(gtab_base + gtaddr2[0 + (kk << 2)]), *(uint*)(gtab_base + gtaddr2[1 + (kk << 2)]), *(uint*)(gtab_base + gtaddr2[2 + (kk << 2)]), *(uint*)(gtab_base + gtaddr2[3 + (kk << 2)]) }));

					REP(16) gtploidy[kk] = vandq_u32(lmask[kk], vshrq_n_u32(gtlow[kk], 24));

					REP(16) gtlow[kk] = vandq_u32(gtlow[kk], mask24);

					REP(16) gtploidy[kk] = vld1q_u32(((uint[]) { PT_PLOIDYxNALLELES[gtploidy2[0 + (kk << 2)]], PT_PLOIDYxNALLELES[gtploidy2[1 + (kk << 2)]], PT_PLOIDYxNALLELES[gtploidy2[2 + (kk << 2)]], PT_PLOIDYxNALLELES[gtploidy2[3 + (kk << 2)]] }));

					REP(16) gtals[kk] = vaddq_u32(gtaddr[kk], gtlow[kk]);

					int maxv = maxploidy;
					if (maxploidy != minploidy)
					{
						uint32x4_t maxv1 =
							vmaxq_u32(
								vmaxq_u32(
									vmaxq_u32(
										vmaxq_u32(gtploidy[0], gtploidy[1]),
										vmaxq_u32(gtploidy[2], gtploidy[3])),
									vmaxq_u32(
										vmaxq_u32(gtploidy[4], gtploidy[5]),
										vmaxq_u32(gtploidy[6], gtploidy[7]))),
								vmaxq_u32(
									vmaxq_u32(
										vmaxq_u32(gtploidy[8], gtploidy[9]),
										vmaxq_u32(gtploidy[10], gtploidy[11])),
									vmaxq_u32(
										vmaxq_u32(gtploidy[12], gtploidy[13]),
										vmaxq_u32(gtploidy[14], gtploidy[15]))));

						maxv = _neo_reduce_max_epi32(maxv1);
					}

					uint32x4_t ai = vdupq_n_u32(0);
					for (int a = 0; a < maxv; ++a, ai = vaddq_u32(ai, mask01))
					{
						REP(16)
						{
							typed[kk] = vcgtq_u32(gtploidy[kk], ai);

							if constexpr (std::is_same_v<REAL, double>)
							{
                                typed64[0 + (kk << 1)] = vld1q_u64(((uint64[]) { (uint64)vgetq_lane_s32(typed[kk], 0), (uint64)vgetq_lane_s32(typed[kk], 1) }));
                                typed64[1 + (kk << 1)] = vld1q_u64(((uint64[]) { (uint64)vgetq_lane_s32(typed[kk], 2), (uint64)vgetq_lane_s32(typed[kk], 3) }));
							}

							uint32x4_t als3 = vandq_u32(gtals[kk], typed[kk]);

                            allele[kk] = vld1q_u32(((uint []) { *(ushort*)(gtab_base + vgetq_lane_u32(als3, 0)), *(ushort*)(gtab_base + vgetq_lane_u32(als3, 1)), *(ushort*)(gtab_base + vgetq_lane_u32(als3, 2)), *(ushort*)(gtab_base + vgetq_lane_u32(als3, 3)) }));

							gtals[kk] = vaddq_u32(gtals[kk], mask02);

							allele[kk] = vaddq_u32(allele[kk], oindex[kk]);

							allele[kk] = vandq_u32(allele[kk], typed[kk]);
						}

						REAL* p2 = p;
						for (int k = 0; k < K; ++k, p2 += KT)
						{
							if constexpr (std::is_same_v<REAL, double>)
							{
								REP(32) freq[kk] = vbslq_f64(typed64[kk], vld1q_f64(((double []) { p2[allele2[0 + (kk << 1)]] , p2[allele2[1 + (kk << 1)]] })), maskoned);
							}
							else
							{
								REP(16)
								{
									float32x4_t v2 = vbslq_f64(typed[kk], vld1q_f32(((float []) { p2[allele2[0 + (kk << 2)]], p2[allele2[1 + (kk << 2)]], p2[allele2[2 + (kk << 2)]], p2[allele2[3 + (kk << 2)]] })), maskones);
									freq[0 + (kk << 1)] = vcvt_f64_f32(vget_low_f32 (v2));
									freq[1 + (kk << 1)] = vcvt_f64_f32(vget_high_f32(v2));
								}
							}

							for (int KK = sizeof(freq) / sizeof(freq[0]) / 2; KK >= 1; KK >>= 1)
								REP(KK) freq[kk] = vmulq_f64(freq[kk], freq[kk + KK]);

							ChargeLog(buf1[k], buf2[k], vgetq_lane_f64(freq[0], 0) * vgetq_lane_f64(freq[0], 1));
						}
					}
				}
			}
		}

		//avoid thread-conflict
		CloseLog((int64*)bufNK1 + N * K * tid2, bufNK2 + N * K * tid2, N * K);
	}
}

/* Update individual or allele origin when ancetral proportion is binary */
template<typename REAL>
TARGETNEO void BAYESIAN<REAL>::UpdateZNoAdmixNEO(int tid)
{
	if (tid == -1)
	{
		SetZero(Mi, N * K);
		SetZero(Ni, K * KT);

		for (int i = 0; i < N; ++i)
			Mi[i * K + Z[i]] = ainds<REAL>[i]->vt;

		//////////////////////////////////////////////////////////

		SetZero((int64*)l_atomic, 32);

		UpdateZNoAdmixNEO(0);

		return;
	}

	static uint32x4_t mask00 = vdupq_n_u32(0);
	static uint32x4_t maskff = vdupq_n_u32(0xFFFFFFFF);
	static uint32x4_t mask24 = vdupq_n_u32(0xFFFFFF);
	static uint32x4_t mask01 = vdupq_n_u32(1);
	static uint32x4_t mask02 = vdupq_n_u32(2);
	alignas(16) static uint maskidx1[64] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63 };
	static uint32x4_t* maskidx = (uint32x4_t*)maskidx1;
	static uint PT_PLOIDYxNALLELES[150] = 									//Pattern index to ploidy level
	{ 0, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

#pragma omp parallel num_threads(structure_nsubthread)
	{
		for (int lsize = structure_loc_size_min; lsize <= structure_loc_size_max; ++lsize)
		{
			int64 lstart = structure_loc_lend[lsize - 1], lend = structure_loc_lend[lsize];
			int64 lend0 = (lend - lstart + 63) / 64;

			for (int64 ia = l_atomic[lsize].fetch_add(1); ia < lend0; ia = l_atomic[lsize].fetch_add(1))
			{
				int64 l = lstart + (ia * structure_loc_coprime64[lsize] % lend0 << 6);
				int64 o = allele_freq_offset[l];
				int64 gtab_base = (int64)GetLoc(l).GetGtab();

				uint32x4_t gtab[16]; uint* gtab1 = (uint*)gtab;
				REP(64) gtab1[kk] = GetLocTabDiff(l + kk);

				uint32x4_t lmask[16], lmaskidx[16];
				uint64 lmask0 = lend - l >= 64 ? 0xFFFFFFFFFFFFFFFF : (1ull << (lend - l)) - 1ull;
				uint32x4_t lmaskt = vdupq_n_u32(lmask0);
				uint lmaskv1[4] = { 1, 2, 4, 8 }; uint32x4_t lmaskv2 = vld1q_u32(lmaskv1);
				REP(16)
				{
                    if (kk == 8) lmaskt = vdupq_n_u32(lmask0 >> 32);
					lmask[kk] = vcgtq_u32(vandq_u32(lmaskt, lmaskv2), mask00);
					lmaskt = vshrq_n_u32(lmaskt, 4);
					lmaskidx[kk] = vandq_u32(lmask[kk], maskidx[kk]);
				}

				GENO_READERNEO<REAL> rt(0, l, lend - l);

				uint32x4_t oindex[16];
				uint64* toffset = allele_freq_offset + l;
				uint64 oindex01 = toffset[0];
				int* lmaskidx2 = (int*)lmaskidx;

                REP(16) oindex[kk] = vld1q_u32(((uint[]) { (uint)(toffset[lmaskidx2[0 + (kk << 2)]] - oindex01), (uint)(toffset[lmaskidx2[1 + (kk << 2)]] - oindex01), (uint)(toffset[lmaskidx2[2 + (kk << 2)]] - oindex01), (uint)(toffset[lmaskidx2[3 + (kk << 2)]] - oindex01) }));

				uint32x4_t gtaddr[16], gtlow[16], gtploidy[16], gtals[16], allele[16], typed[16];
				int* typed2 = (int*)typed, * allele2 = (int*)allele, * gtploidy2 = (int*)gtploidy, * gtaddr2 = (int*)gtaddr;

				for (int i = 0; i < N; i++)
				{
					rt.Read(gtaddr);
                    
                    REP(16) gtaddr[kk] = vaddq_u32(gtab[kk], vshlq_n_u32(gtaddr[kk], 2));

                    REP(16) gtlow[kk] = vld1q_u32(((uint[]) { *(uint*)(gtab_base + gtaddr2[0 + (kk << 2)]), *(uint*)(gtab_base + gtaddr2[1 + (kk << 2)]), *(uint*)(gtab_base + gtaddr2[2 + (kk << 2)]), *(uint*)(gtab_base + gtaddr2[3 + (kk << 2)]) }));

                    REP(16) gtploidy[kk] = vandq_u32(lmask[kk], vshrq_n_u32(gtlow[kk], 24));

                    REP(16) gtlow[kk] = vandq_u32(gtlow[kk], mask24);

                    REP(16) gtploidy[kk] = vld1q_u32(((uint[]) { PT_PLOIDYxNALLELES[gtploidy2[0 + (kk << 2)]], PT_PLOIDYxNALLELES[gtploidy2[1 + (kk << 2)]], PT_PLOIDYxNALLELES[gtploidy2[2 + (kk << 2)]], PT_PLOIDYxNALLELES[gtploidy2[3 + (kk << 2)]] }));

                    REP(16) gtals[kk] = vaddq_u32(gtaddr[kk], gtlow[kk]);

					int maxv = maxploidy;
					if (maxploidy != minploidy)
					{
						uint32x4_t maxv1 =
							vmaxq_u32(
								vmaxq_u32(
									vmaxq_u32(
										vmaxq_u32(gtploidy[0], gtploidy[1]),
										vmaxq_u32(gtploidy[2], gtploidy[3])),
									vmaxq_u32(
										vmaxq_u32(gtploidy[4], gtploidy[5]),
										vmaxq_u32(gtploidy[6], gtploidy[7]))),
								vmaxq_u32(
									vmaxq_u32(
										vmaxq_u32(gtploidy[8], gtploidy[9]),
										vmaxq_u32(gtploidy[10], gtploidy[11])),
									vmaxq_u32(
										vmaxq_u32(gtploidy[12], gtploidy[13]),
										vmaxq_u32(gtploidy[14], gtploidy[15]))));
						
						maxv = _neo_reduce_max_epi32(maxv1);
					}

					int* ni = Ni + Z[i] * KT + o;
					uint32x4_t ai = vdupq_n_u32(0);
					for (int a = 0; a < maxv; ++a, ai = vaddq_u32(ai, mask01))
					{
						REP(16)
						{
							typed[kk] = vcgtq_u32(gtploidy[kk], ai);

							uint32x4_t als3 = vandq_u32(gtals[kk], typed[kk]);

							allele[kk] = vld1q_u32(((uint []) { *(ushort*)(gtab_base + vgetq_lane_u32(als3, 0)), *(ushort*)(gtab_base + vgetq_lane_u32(als3, 1)), *(ushort*)(gtab_base + vgetq_lane_u32(als3, 2)), *(ushort*)(gtab_base + vgetq_lane_u32(als3, 3)) }));

							gtals[kk] = vaddq_u32(gtals[kk], mask02);

							allele[kk] = vaddq_u32(allele[kk], oindex[kk]);

							allele[kk] = vandq_u32(allele[kk], typed[kk]);
						}

						REP(64) ni[allele2[kk]] -= typed2[kk];
					}
				}
			}
		}
	}
}

/* Update a priori ancetral proportion for non-admix model */
template<typename REAL>
template<bool fast_fp32>
TARGETNEO void BAYESIAN<REAL>::UpdateQMetroNEO(int tid)
{
	if (tid == -1)
	{
		RNG<double> rng(seed + m, RNG_SALT_UPDATEQ);//REAL
		REAL* bufi = (REAL*)bufNK1;
		REAL* q = NULL;

		for (int i = 0; i < N; ++i, bufi += K)
		{
			if (ainds<REAL>[i]->vt == 0) continue;
			if (locpriori) rng.Dirichlet(bufi, AlphaLocal + ainds<REAL>[i]->popid * K, K);
			else           rng.Dirichlet(bufi, Alpha, K);
		}

		OpenLog((int64*)bufN1, bufN2, N * structure_nsubthread);

		//////////////////////////////////////////////////////////

		SetZero((int64*)l_atomic, 32);

		g_fastsingle_val == 1 ? UpdateQMetroNEO<true >(0) : UpdateQMetroNEO<false>(0);

		//avoid thread-conflict
		for (int i = 1; i < structure_nsubthread; ++i)
			Add(bufN1, bufN1 + N * i, N);

		bufi = (REAL*)bufNK1; q = Q;
		for (int i = 0; i < N; ++i, q += K, bufi += K)
		{
			if (ainds<REAL>[i]->vt == 0) continue;
			if (bufN1[i] >= NZERO || rng.Uniform() < exp(bufN1[i]))
				SetVal(q, bufi, K);
		}
		return;
	}

	static float64x2_t maskoned = vdupq_n_f64(1.0);
	static float32x4_t maskones = vdupq_n_f32(1.0);
	static uint32x4_t mask00 = vdupq_n_u32(0);
	static uint32x4_t maskff = vdupq_n_u32(0xFFFFFFFF);
	static uint32x4_t mask24 = vdupq_n_u32(0xFFFFFF);
	static uint32x4_t mask01 = vdupq_n_u32(1);
	static uint32x4_t mask02 = vdupq_n_u32(2);
	alignas(16) static uint maskidx1[64] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63 };
	static uint32x4_t* maskidx = (uint32x4_t*)maskidx1;
	static uint PT_PLOIDYxNALLELES[150] = 									//Pattern index to ploidy level
	{ 0, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	atomic<int> thread_counter = 0;
#pragma omp parallel num_threads(structure_nsubthread)
	{
		int tid = thread_counter.fetch_add(1);

		for (int lsize = structure_loc_size_min; lsize <= structure_loc_size_max; ++lsize)
		{
			int64 lstart = structure_loc_lend[lsize - 1], lend = structure_loc_lend[lsize];
			int64 lend0 = (lend - lstart + 63) / 64;

			for (int64 ia = l_atomic[lsize].fetch_add(1); ia < lend0; ia = l_atomic[lsize].fetch_add(1))
			{
				int64 l = lstart + (ia * structure_loc_coprime64[lsize] % lend0 << 6);
				REAL* p = Freq + allele_freq_offset[l];
				int64 gtab_base = (int64)GetLoc(l).GetGtab();

				uint32x4_t gtab[16]; uint* gtab1 = (uint*)gtab;
				REP(64) gtab1[kk] = GetLocTabDiff(l + kk);

				uint32x4_t lmask[16], lmaskidx[16];
				uint64 lmask0 = lend - l >= 64 ? 0xFFFFFFFFFFFFFFFF : (1ull << (lend - l)) - 1ull;
				uint32x4_t lmaskt = vdupq_n_u32(lmask0);
				uint lmaskv1[4] = { 1, 2, 4, 8 }; uint32x4_t lmaskv2 = vld1q_u32(lmaskv1);
				REP(16)
				{
                    if (kk == 8) lmaskt = vdupq_n_u32(lmask0 >> 32);
					lmask[kk] = vcgtq_u32(vandq_u32(lmaskt, lmaskv2), mask00);
					lmaskt = vshrq_n_u32(lmaskt, 4);
					lmaskidx[kk] = vandq_u32(lmask[kk], maskidx[kk]);
				}

				GENO_READERNEO<REAL> rt(0, l, lend - l);

				uint32x4_t oindex[16];
				uint64* toffset = allele_freq_offset + l;
				uint64 oindex01 = toffset[0];
				int* lmaskidx2 = (int*)lmaskidx;

                REP(16) oindex[kk] = vld1q_u32(((uint[]) { (uint)(toffset[lmaskidx2[0 + (kk << 2)]] - oindex01), (uint)(toffset[lmaskidx2[1 + (kk << 2)]] - oindex01), (uint)(toffset[lmaskidx2[2 + (kk << 2)]] - oindex01), (uint)(toffset[lmaskidx2[3 + (kk << 2)]] - oindex01) }));
                
                uint32x4_t gtaddr[16], gtlow[16], gtploidy[16], gtals[16], allele[16], typed[16];
                uint64x2_t typed64[32];
				int* allele2 = (int*)allele, * gtploidy2 = (int*)gtploidy, * gtaddr2 = (int*)gtaddr;
				float64x2_t f1[32], f2[32];
				float32x4_t f1s[16], f2s[16];

				REAL* bufi = (REAL*)bufNK1, * q = Q;
				//avoid thread-conflict
				double* buf1 = bufN1 + N * tid, * buf2 = bufN2 + N * tid;

				for (int i = 0; i < N; ++i, q += K, bufi += K, buf1++, buf2++)
				{
                    rt.Read(gtaddr);
                    
                    REP(16) gtaddr[kk] = vaddq_u32(gtab[kk], vshlq_n_u32(gtaddr[kk], 2));

                    REP(16) gtlow[kk] = vld1q_u32(((uint[]) { *(uint*)(gtab_base + gtaddr2[0 + (kk << 2)]), *(uint*)(gtab_base + gtaddr2[1 + (kk << 2)]), *(uint*)(gtab_base + gtaddr2[2 + (kk << 2)]), *(uint*)(gtab_base + gtaddr2[3 + (kk << 2)]) }));

                    REP(16) gtploidy[kk] = vandq_u32(lmask[kk], vshrq_n_u32(gtlow[kk], 24));

                    REP(16) gtlow[kk] = vandq_u32(gtlow[kk], mask24);

                    REP(16) gtploidy[kk] = vld1q_u32(((uint[]) { PT_PLOIDYxNALLELES[gtploidy2[0 + (kk << 2)]], PT_PLOIDYxNALLELES[gtploidy2[1 + (kk << 2)]], PT_PLOIDYxNALLELES[gtploidy2[2 + (kk << 2)]], PT_PLOIDYxNALLELES[gtploidy2[3 + (kk << 2)]] }));

                    REP(16) gtals[kk] = vaddq_u32(gtaddr[kk], gtlow[kk]);
                    
					int maxv = maxploidy;
					if (maxploidy != minploidy)
					{
						uint32x4_t maxv1 =
							vmaxq_u32(
								vmaxq_u32(
									vmaxq_u32(
										vmaxq_u32(gtploidy[0], gtploidy[1]),
										vmaxq_u32(gtploidy[2], gtploidy[3])),
									vmaxq_u32(
										vmaxq_u32(gtploidy[4], gtploidy[5]),
										vmaxq_u32(gtploidy[6], gtploidy[7]))),
								vmaxq_u32(
									vmaxq_u32(
										vmaxq_u32(gtploidy[8], gtploidy[9]),
										vmaxq_u32(gtploidy[10], gtploidy[11])),
									vmaxq_u32(
										vmaxq_u32(gtploidy[12], gtploidy[13]),
										vmaxq_u32(gtploidy[14], gtploidy[15]))));
						
						maxv = _neo_reduce_max_epi32(maxv1);
					}

					uint32x4_t ai = vdupq_n_u32(0);
					for (int a = 0; a < maxv; ++a, ai = vaddq_u32(ai, mask01))
					{
						REP(16)
						{
							typed[kk] = vcgtq_u32(gtploidy[kk], ai);

							if constexpr (std::is_same_v<REAL, double> || !fast_fp32)
							{
                                typed64[0 + (kk << 1)] = vld1q_u64(((uint64[]) { (uint64)vgetq_lane_s32(typed[kk], 0), (uint64)vgetq_lane_s32(typed[kk], 1) }));
                                typed64[1 + (kk << 1)] = vld1q_u64(((uint64[]) { (uint64)vgetq_lane_s32(typed[kk], 2), (uint64)vgetq_lane_s32(typed[kk], 3) }));
							}

							uint32x4_t als3 = vandq_u32(gtals[kk], typed[kk]);

                            allele[kk] = vld1q_u32(((uint []) { *(ushort*)(gtab_base + vgetq_lane_u32(als3, 0)), *(ushort*)(gtab_base + vgetq_lane_u32(als3, 1)), *(ushort*)(gtab_base + vgetq_lane_u32(als3, 2)), *(ushort*)(gtab_base + vgetq_lane_u32(als3, 3)) }));

							gtals[kk] = vaddq_u32(gtals[kk], mask02);

							allele[kk] = vaddq_u32(allele[kk], oindex[kk]);

							allele[kk] = vandq_u32(allele[kk], typed[kk]);
						}

						if constexpr (std::is_same_v<REAL, double> || !fast_fp32)
						{
							REP(32) f1[kk] = vdupq_n_f64(0);
							REP(32) f2[kk] = vdupq_n_f64(0);
						}
						else
						{
							REP(16) f1s[kk] = vdupq_n_f32(0);
							REP(16) f2s[kk] = vdupq_n_f32(0);
						}

						REAL* p2 = p;
						for (int k = 0; k < K; ++k, p2 += KT)
						{
							if constexpr (std::is_same_v<REAL, double> || !fast_fp32)
							{
								float64x2_t pp[32], qq = vdupq_n_f64(q[k]), ii = vdupq_n_f64(bufi[k]);

								REP(32)
								{
									pp[kk] = vld1q_f64(((double []) { p2[allele2[0 + (kk << 1)]], p2[allele2[1 + (kk << 1)]] }));
									f1[kk] = vaddq_f64(f1[kk], vmulq_f64(pp[kk], ii));
									f2[kk] = vaddq_f64(f2[kk], vmulq_f64(pp[kk], qq));
								}
							}
							else
							{
								float32x4_t pp[16], v1[16], v2[16], qq = vdupq_n_f32(q[k]), ii = vdupq_n_f32(bufi[k]);

								REP(16)
								{
									pp[kk] = vld1q_f32(((float[]) { p2[allele2[0 + (kk << 2)]], p2[allele2[1 + (kk << 2)]], p2[allele2[2 + (kk << 2)]], p2[allele2[3 + (kk << 2)]] }));
									v1[kk] = vmulq_f32(pp[kk], ii);
									v2[kk] = vmulq_f32(pp[kk], qq);
									f1s[kk] = vaddq_f32(f1s[kk], v1[kk]);
									f2s[kk] = vaddq_f32(f2s[kk], v2[kk]);
								}
							}
						}

						if constexpr (std::is_same_v<REAL, float> && fast_fp32)
						{
							REP(16) f1s[kk] = vdivq_f32(f1s[kk], f2s[kk]);
							REP(16) f1s[kk] = vbslq_f32(typed[kk], f1s[kk], maskones);
							REP(16)
							{
								f1[0 + (kk << 1)] = vcvt_f64_f32(vget_low_f32(f1s[kk]));
								f1[1 + (kk << 1)] = vcvt_f64_f32(vget_high_f32(f1s[kk]));
							}
						}
						else
						{
							REP(32)
							{
								f1[kk] = vdivq_f64(f1[kk], f2[kk]);
								f1[kk] = vbslq_f64(typed64[kk], f1[kk], maskoned);
							}
						}

						for (int KK = sizeof(f1) / sizeof(f1[0]) / 2; KK >= 1; KK >>= 1)
							REP(KK) f1[kk] = vmulq_f64(f1[kk], f1[kk + KK]);

						ChargeLog(*(int64*)buf1, *buf2, vgetq_lane_f64(f1[0], 0) * vgetq_lane_f64(f1[0], 1));
					}
				}
			}
		}

		//avoid thread-conflict
		CloseLog((int64*)bufN1 + N * tid, bufN2 + N * tid, N);
	}
}

/* Update individual or allele origin when ancetral proportion */
template<typename REAL>
template<bool fast_fp32>
TARGETNEO void BAYESIAN<REAL>::UpdateZAdmixNEO(int tid)
{
	if (tid == -1)
	{
		SetZero(Mi, N * K * structure_nsubthread);
		SetZero(Ni, K * KT);

		//////////////////////////////////////////////////////////

		SetZero((int64*)l_atomic, 32);

		g_fastsingle_val == 1 ? UpdateZAdmixNEO<true >(0) : UpdateZAdmixNEO<false>(0);

		//avoid thread-conflict
		for (int i = 1; i < structure_nsubthread; ++i)
			Add(Mi, Mi + N * K * i, N * K);

		return;
	}

	static uint32x4_t mask00 = vdupq_n_u32(0);
	static uint32x4_t maskff = vdupq_n_u32(0xFFFFFFFF);
	static uint32x4_t mask24 = vdupq_n_u32(0xFFFFFF);
	static uint32x4_t mask01 = vdupq_n_u32(1);
	static uint32x4_t mask02 = vdupq_n_u32(2);
	alignas(16) static uint maskidx1[64] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63 };
	static uint32x4_t* maskidx = (uint32x4_t*)maskidx1;
	static float64x2_t minfreqd = vdupq_n_f64(MIN_FREQ);
	static float32x4_t minfreqs = vdupq_n_f32(MIN_FREQ);
	static uint PT_PLOIDYxNALLELES[150] = 									//Pattern index to ploidy level
	{ 0, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	atomic<int> thread_counter = 0;
#pragma omp parallel num_threads(structure_nsubthread)
	{
		int tid = thread_counter.fetch_add(1);

		//64 * K  * sizeof(double)
		float64x2_t* bufkd = (float64x2_t*)Align64((byte*)bufNK2 + (std::max(N, 64) * K * sizeof(double) + 63) * tid);
		float32x4_t* bufks = (float32x4_t*)Align64((byte*)bufNK2 + (std::max(N, 64) * K * sizeof(double) + 63) * tid);

		for (int lsize = structure_loc_size_min; lsize <= structure_loc_size_max; ++lsize)
		{
			int64 lstart = structure_loc_lend[lsize - 1], lend = structure_loc_lend[lsize];
			int64 lend0 = (lend - lstart + 63) / 64;

			for (int64 ia = l_atomic[lsize].fetch_add(1); ia < lend0; ia = l_atomic[lsize].fetch_add(1))
			{
				int64 l = lstart + (ia * structure_loc_coprime64[lsize] % lend0 << 6);
				REAL* p = Freq + allele_freq_offset[l];
				int* ni = Ni + allele_freq_offset[l];
				int64 gtab_base = (int64)GetLoc(l).GetGtab();

				uint32x4_t gtab[16]; uint* gtab1 = (uint*)gtab;
				REP(64) gtab1[kk] = GetLocTabDiff(l + kk);

				uint32x4_t lmask[16], lmaskidx[16];
				uint64 lmask0 = lend - l >= 64 ? 0xFFFFFFFFFFFFFFFF : (1ull << (lend - l)) - 1ull;
				uint32x4_t lmaskt = vdupq_n_u32(lmask0);
				uint lmaskv1[4] = { 1, 2, 4, 8 }; uint32x4_t lmaskv2 = vld1q_u32(lmaskv1);
				REP(16)
				{
                    if (kk == 8) lmaskt = vdupq_n_u32(lmask0 >> 32);
					lmask[kk] = vcgtq_u32(vandq_u32(lmaskt, lmaskv2), mask00);
					lmaskt = vshrq_n_u32(lmaskt, 4);
					lmaskidx[kk] = vandq_u32(lmask[kk], maskidx[kk]);
				}

				GENO_READERNEO<REAL> rt(0, l, lend - l);

				uint32x4_t oindex[16];
				uint64* toffset = allele_freq_offset + l;
				uint64 oindex01 = toffset[0];
				int* lmaskidx2 = (int*)lmaskidx;

                REP(16) oindex[kk] = vld1q_u32(((uint[]) { (uint)(toffset[lmaskidx2[0 + (kk << 2)]] - oindex01), (uint)(toffset[lmaskidx2[1 + (kk << 2)]] - oindex01), (uint)(toffset[lmaskidx2[2 + (kk << 2)]] - oindex01), (uint)(toffset[lmaskidx2[3 + (kk << 2)]] - oindex01) }));

				uint32x4_t gtaddr[16], gtlow[16], gtploidy[16], gtals[16], allele[16], typed[16];
				uint64x2_t allele64[32], k2[32], typed64[32];
				int* typed2 = (int*)typed, * allele2 = (int*)allele, * gtploidy2 = (int*)gtploidy, * gtaddr2 = (int*)gtaddr;
				uint64* k22 = (uint64*)&k2;

				RNGNEO<double> rngd; RNGNEO<float > rngs;

				if constexpr (std::is_same_v<REAL, double> || !fast_fp32)
					new (&rngd) RNGNEO<double>(seed + m * L + l, RNG_SALT_UPDATEZ);
				else
					new (&rngs) RNGNEO<float >(seed + m * L + l, RNG_SALT_UPDATEZ);

				REAL* q = Q;
				//avoid thread-conflict
				int64* mi = Mi + N * K * tid;

				for (int i = 0; i < N; ++i, q += K, mi += K)
				{
                    rt.Read(gtaddr);
                    
                    REP(16) gtaddr[kk] = vaddq_u32(gtab[kk], vshlq_n_u32(gtaddr[kk], 2));

                    REP(16) gtlow[kk] = vld1q_u32(((uint[]) { *(uint*)(gtab_base + gtaddr2[0 + (kk << 2)]), *(uint*)(gtab_base + gtaddr2[1 + (kk << 2)]), *(uint*)(gtab_base + gtaddr2[2 + (kk << 2)]), *(uint*)(gtab_base + gtaddr2[3 + (kk << 2)]) }));

                    REP(16) gtploidy[kk] = vandq_u32(lmask[kk], vshrq_n_u32(gtlow[kk], 24));

                    REP(16) gtlow[kk] = vandq_u32(gtlow[kk], mask24);

                    REP(16) gtploidy[kk] = vld1q_u32(((uint[]) { PT_PLOIDYxNALLELES[gtploidy2[0 + (kk << 2)]], PT_PLOIDYxNALLELES[gtploidy2[1 + (kk << 2)]], PT_PLOIDYxNALLELES[gtploidy2[2 + (kk << 2)]], PT_PLOIDYxNALLELES[gtploidy2[3 + (kk << 2)]] }));

                    REP(16) gtals[kk] = vaddq_u32(gtaddr[kk], gtlow[kk]);

					int maxv = maxploidy;

					uint32x4_t ai = vdupq_n_u32(0);
					for (int a = 0; a < maxv; ++a, ai = vaddq_u32(ai, mask01))
					{
						REP(16)
						{
							typed[kk] = vcgtq_u32(gtploidy[kk], ai);

							uint32x4_t als3 = vandq_u32(gtals[kk], typed[kk]);

                            allele[kk] = vld1q_u32(((uint []) { *(ushort*)(gtab_base + vgetq_lane_u32(als3, 0)), *(ushort*)(gtab_base + vgetq_lane_u32(als3, 1)), *(ushort*)(gtab_base + vgetq_lane_u32(als3, 2)), *(ushort*)(gtab_base + vgetq_lane_u32(als3, 3)) }));

							gtals[kk] = vaddq_u32(gtals[kk], mask02);

							allele[kk] = vaddq_u32(allele[kk], oindex[kk]);

							allele[kk] = vandq_u32(allele[kk], typed[kk]);
						}

						REP(16)
						{
							{
								allele64[0 + (kk << 1)] = vld1q_u64(((uint64[]) { vgetq_lane_u32(allele[kk], 0), vgetq_lane_u32(allele[kk], 1)}));
								allele64[1 + (kk << 1)] = vld1q_u64(((uint64[]) { vgetq_lane_u32(allele[kk], 2), vgetq_lane_u32(allele[kk], 3)}));
							}

							if constexpr (std::is_same_v<REAL, double> || !fast_fp32)
							{
                                typed64[0 + (kk << 1)] = vld1q_u64(((uint64[]) { (uint64)vgetq_lane_s32(typed[kk], 0), (uint64)vgetq_lane_s32(typed[kk], 1) }));
                                typed64[1 + (kk << 1)] = vld1q_u64(((uint64[]) { (uint64)vgetq_lane_s32(typed[kk], 2), (uint64)vgetq_lane_s32(typed[kk], 3) }));
							}
						}

						REAL* p2 = p;
						for (int k = 0; k < K; ++k, p2 += KT)
						{
							float64x2_t qd; float32x4_t qs;

							if constexpr (std::is_same_v<REAL, double>)
							{
								qd = vdupq_n_f64(q[k]);
								REP(32) bufkd[kk + (k << 5)] = vld1q_f64(((double []) { p2[allele2[0 + (kk << 1)]], p2[allele2[1 + (kk << 1)]] }));
								REP(32) bufkd[kk + (k << 5)] = vmulq_f64(bufkd[kk + (k << 5)], qd);
								REP(32) bufkd[kk + (k << 5)] = vaddq_f64(bufkd[kk + (k << 5)], minfreqd);
								REP(32) bufkd[kk + (k << 5)] = vandq_u64(bufkd[kk + (k << 5)], typed64[kk]);
							}
							else if constexpr (fast_fp32)
							{
								qs = vdupq_n_f32(q[k]);

								REP(16) bufks[kk + (k << 4)] = vld1q_f32(((float []) { p2[allele2[0 + (kk << 2)]], p2[allele2[1 + (kk << 2)]], p2[allele2[2 + (kk << 2)]], p2[allele2[3 + (kk << 2)]] }));
								REP(16) bufks[kk + (k << 4)] = vmulq_f32(bufks[kk + (k << 4)], qs);
								REP(16) bufks[kk + (k << 4)] = vaddq_f32(bufks[kk + (k << 4)], minfreqs);
								REP(16) bufks[kk + (k << 4)] = vandq_u32(bufks[kk + (k << 4)], typed[kk]);
							}
							else
							{
								qd = vdupq_n_f64(q[k]);
								REP(16)
								{
									float32x4_t v2 = vld1q_f32(((float []) { p2[allele2[0 + (kk << 2)]], p2[allele2[1 + (kk << 2)]], p2[allele2[2 + (kk << 2)]], p2[allele2[3 + (kk << 2)]] }));
									bufkd[0 + (kk << 1) + (k << 5)] = vcvt_f64_f32(vget_low_f32 (v2));
									bufkd[1 + (kk << 1) + (k << 5)] = vcvt_f64_f32(vget_high_f32(v2));
								}
								REP(32) bufkd[kk + (k << 5)] = vmulq_f64(bufkd[kk + (k << 5)], qd);
								REP(32) bufkd[kk + (k << 5)] = vaddq_f64(bufkd[kk + (k << 5)], minfreqd);
								REP(32) bufkd[kk + (k << 5)] = vandq_u64(bufkd[kk + (k << 5)], typed64[kk]);
							}
						}

						//draw cluster for each allele copy
						if constexpr (std::is_same_v<REAL, double> || !fast_fp32)
							rngd.Poly(bufkd, K, k2);
						else
							rngs.Poly(bufks, K, k2);

						//Update Mi
						REP(64) mi[k22[kk]] -= typed2[kk];

						REP(32) k2[kk] = vld1q_u64(((uint64[]) { k22[0 + (kk << 1)] * KT, k22[1 + (kk << 1)] * KT }));

						REP(32) k2[kk] = vaddq_u64(k2[kk], allele64[kk]);

						//ni[k2 * KT + als[a]]++;
						REP(64) ni[k22[kk]] -= typed2[kk];
					}
				}
			}
		}
	}
}

/* Record updated MCMC parameters */
template<typename REAL>
template<bool isadmix, bool fast_fp32>
TARGETNEO void BAYESIAN<REAL>::RecordNEO(int tid)
{
	if (tid == -1)
	{
		Add(MiSum, Mi, N * K);

		//////////////////////////////////////////////////////////

		SetZero((int64*)l_atomic, 32);

		switch ((int)binaryq * 10 + g_fastsingle_val)
		{
		case 01: BAYESIAN<REAL>::RecordNEO<true, true >(0); break;
		case 02: BAYESIAN<REAL>::RecordNEO<true, false>(0); break;
		case 11: BAYESIAN<REAL>::RecordNEO<false, true >(0); break;
		case 12: BAYESIAN<REAL>::RecordNEO<false, false>(0); break;
		}

		bufNK1[0] = Sum(bufNK1, structure_nsubthread);
		return;
	}

	static float64x2_t maskoned = vdupq_n_f64(1.0);
	static uint32x4_t mask00 = vdupq_n_u32(0);
	static uint32x4_t maskff = vdupq_n_u32(0xFFFFFFFF);
	static uint32x4_t mask24 = vdupq_n_u32(0xFFFFFF);
	static uint32x4_t mask01 = vdupq_n_u32(1);
	static uint32x4_t mask02 = vdupq_n_u32(2);
	alignas(16) static uint maskidx1[64] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63 };
	static uint32x4_t* maskidx = (uint32x4_t*)maskidx1;
	static uint PT_PLOIDYxNALLELES[150] = 									//Pattern index to ploidy level
	{ 0, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	atomic<int> thread_counter = 0;
#pragma omp parallel num_threads(structure_nsubthread)
	{
		int tid = thread_counter.fetch_add(1);

		int64 slog = 0; double prod = 1;
		OpenLog(slog, prod);

		for (int lsize = structure_loc_size_min; lsize <= structure_loc_size_max; ++lsize)
		{
			int64 lstart = structure_loc_lend[lsize - 1], lend = structure_loc_lend[lsize];
			int64 lend0 = (lend - lstart + 63) / 64;

			for (int64 ia = l_atomic[lsize].fetch_add(1); ia < lend0; ia = l_atomic[lsize].fetch_add(1))
			{
				int64 l = lstart + (ia * structure_loc_coprime64[lsize] % lend0 << 6);
				REAL* p = Freq + allele_freq_offset[l];
				int64 gtab_base = (int64)GetLoc(l).GetGtab();

				uint32x4_t gtab[16]; uint* gtab1 = (uint*)gtab;
				REP(64) gtab1[kk] = GetLocTabDiff(l + kk);

				uint32x4_t lmask[16], lmaskidx[16];
				uint64 lmask0 = lend - l >= 64 ? 0xFFFFFFFFFFFFFFFF : (1ull << (lend - l)) - 1ull;
				uint32x4_t lmaskt = vdupq_n_u32(lmask0);
				uint lmaskv1[4] = { 1, 2, 4, 8 }; uint32x4_t lmaskv2 = vld1q_u32(lmaskv1);
				REP(16)
				{
                    if (kk == 8) lmaskt = vdupq_n_u32(lmask0 >> 32);
					lmask[kk] = vcgtq_u32(vandq_u32(lmaskt, lmaskv2), mask00);
					lmaskt = vshrq_n_u32(lmaskt, 4);
					lmaskidx[kk] = vandq_u32(lmask[kk], maskidx[kk]);
				}

				GENO_READERNEO<REAL> rt(0, l, lend - l);

				uint32x4_t oindex[16];
				uint64* toffset = allele_freq_offset + l;
				uint64 oindex01 = toffset[0];
				int* lmaskidx2 = (int*)lmaskidx;

                REP(16) oindex[kk] = vld1q_u32(((uint[]) { (uint)(toffset[lmaskidx2[0 + (kk << 2)]] - oindex01), (uint)(toffset[lmaskidx2[1 + (kk << 2)]] - oindex01), (uint)(toffset[lmaskidx2[2 + (kk << 2)]] - oindex01), (uint)(toffset[lmaskidx2[3 + (kk << 2)]] - oindex01) }));

                uint32x4_t gtaddr[16], gtlow[16], gtploidy[16], gtals[16], allele[16], typed[16];
                uint64x2_t allele64[32], typed64[32];
				int* allele2 = (int*)allele, * gtploidy2 = (int*)gtploidy, * gtaddr2 = (int*)gtaddr;
				float64x2_t f2[32];
				float32x4_t f2s[16];

				REAL* q = Q;

				for (int i = 0; i < N; ++i, q += K)
				{
                    rt.Read(gtaddr);
                    
                    REP(16) gtaddr[kk] = vaddq_u32(gtab[kk], vshlq_n_u32(gtaddr[kk], 2));

                    REP(16) gtlow[kk] = vld1q_u32(((uint[]) { *(uint*)(gtab_base + gtaddr2[0 + (kk << 2)]), *(uint*)(gtab_base + gtaddr2[1 + (kk << 2)]), *(uint*)(gtab_base + gtaddr2[2 + (kk << 2)]), *(uint*)(gtab_base + gtaddr2[3 + (kk << 2)]) }));

                    REP(16) gtploidy[kk] = vandq_u32(lmask[kk], vshrq_n_u32(gtlow[kk], 24));

                    REP(16) gtlow[kk] = vandq_u32(gtlow[kk], mask24);

                    REP(16) gtploidy[kk] = vld1q_u32(((uint[]) { PT_PLOIDYxNALLELES[gtploidy2[0 + (kk << 2)]], PT_PLOIDYxNALLELES[gtploidy2[1 + (kk << 2)]], PT_PLOIDYxNALLELES[gtploidy2[2 + (kk << 2)]], PT_PLOIDYxNALLELES[gtploidy2[3 + (kk << 2)]] }));

                    REP(16) gtals[kk] = vaddq_u32(gtaddr[kk], gtlow[kk]);

					int maxv = maxploidy;
					if (maxploidy != minploidy)
					{
						uint32x4_t maxv1 =
							vmaxq_u32(
								vmaxq_u32(
									vmaxq_u32(
										vmaxq_u32(gtploidy[0], gtploidy[1]),
										vmaxq_u32(gtploidy[2], gtploidy[3])),
									vmaxq_u32(
										vmaxq_u32(gtploidy[4], gtploidy[5]),
										vmaxq_u32(gtploidy[6], gtploidy[7]))),
								vmaxq_u32(
									vmaxq_u32(
										vmaxq_u32(gtploidy[8], gtploidy[9]),
										vmaxq_u32(gtploidy[10], gtploidy[11])),
									vmaxq_u32(
										vmaxq_u32(gtploidy[12], gtploidy[13]),
										vmaxq_u32(gtploidy[14], gtploidy[15]))));
						
						maxv = _neo_reduce_max_epi32(maxv1);
					}

					uint32x4_t ai = vdupq_n_u32(0);
					for (int a = 0; a < maxv; ++a, ai = vaddq_u32(ai, mask01))
					{
						REP(16)
						{
							typed[kk] = vcgtq_u32(gtploidy[kk], ai);

							uint32x4_t als3 = vandq_u32(gtals[kk], typed[kk]);

                            allele[kk] = vld1q_u32(((uint []) { *(ushort*)(gtab_base + vgetq_lane_u32(als3, 0)), *(ushort*)(gtab_base + vgetq_lane_u32(als3, 1)), *(ushort*)(gtab_base + vgetq_lane_u32(als3, 2)), *(ushort*)(gtab_base + vgetq_lane_u32(als3, 3)) }));

							gtals[kk] = vaddq_u32(gtals[kk], mask02);

							allele[kk] = vaddq_u32(allele[kk], oindex[kk]);

							allele[kk] = vandq_u32(allele[kk], typed[kk]);
						}

						REP(16)
						{
                            {
                                allele64[0 + (kk << 1)] = vld1q_u64(((uint64[]) { vgetq_lane_u32(allele[kk], 0), vgetq_lane_u32(allele[kk], 1)}));
                                allele64[1 + (kk << 1)] = vld1q_u64(((uint64[]) { vgetq_lane_u32(allele[kk], 2), vgetq_lane_u32(allele[kk], 3)}));
                            }

							{
                                typed64[0 + (kk << 1)] = vld1q_u64(((uint64[]) { (uint64)vgetq_lane_s32(typed[kk], 0), (uint64)vgetq_lane_s32(typed[kk], 1) }));
                                typed64[1 + (kk << 1)] = vld1q_u64(((uint64[]) { (uint64)vgetq_lane_s32(typed[kk], 2), (uint64)vgetq_lane_s32(typed[kk], 3) }));
							}
						}

						if constexpr (std::is_same_v<REAL, double> || !fast_fp32)
							REP(32) f2[kk] = vdupq_n_f64(0);
						else
							REP(16) f2s[kk] = vdupq_n_f32(0);

						if constexpr (!isadmix)
						{
							REAL* p2 = p + Z[i] * KT;

							if constexpr (std::is_same_v<REAL, double> || !fast_fp32)
							{
								REP(32) f2[kk] = vld1q_f64(((double []) { p2[allele2[0 + (kk << 1)]], p2[allele2[1 + (kk << 1)]] }));
							}
							else
							{
								REP(16) f2s[kk] = vld1q_f32(((float []) { p2[allele2[0 + (kk << 2)]], p2[allele2[1 + (kk << 2)]], p2[allele2[2 + (kk << 2)]], p2[allele2[3 + (kk << 2)]] }));
							}
						}
						else
						{
							REAL* p2 = p;
							for (int k = 0; k < K; ++k, p2 += KT)
							{
								//disable fmadd
								if constexpr (std::is_same_v<REAL, double> || !fast_fp32)
								{
									float64x2_t pp[32], qq = vdupq_n_f64(q[k]);

									REP(32) pp[kk] = vld1q_f64(((double[]) { p2[allele2[0 + (kk << 1)]] , p2[allele2[1 + (kk << 1)]] }));
                                    
									REP(32) f2[kk] = vaddq_f64(f2[kk], vmulq_f64(pp[kk], qq));
								}
								else
								{
									float32x4_t pp[16], qq = vdupq_n_f32(q[k]);

									REP(16) pp[kk] = vld1q_f32(((float []) { p2[allele2[0 + (kk << 2)]], p2[allele2[1 + (kk << 2)]], p2[allele2[2 + (kk << 2)]], p2[allele2[3 + (kk << 2)]] }));
                                    
									REP(16) f2s[kk] = vaddq_f32(f2s[kk], vmulq_f32(pp[kk], qq));
								}
							}
						}

						if constexpr (std::is_same_v<REAL, float> && fast_fp32)
						{
							REP(16)
							{
								f2[0 + (kk << 1)] = vcvt_f64_f32(vget_low_f32 (f2s[kk]));
								f2[1 + (kk << 1)] = vcvt_f64_f32(vget_high_f32(f2s[kk]));
							}
						}

						REP(32) f2[kk] = vbslq_f64(typed64[kk], f2[kk], maskoned);

						for (int KK = sizeof(f2) / sizeof(f2[0]) / 2; KK >= 1; KK >>= 1)
							REP(KK) f2[kk] = vmulq_f64(f2[kk], f2[kk + KK]);

						ChargeLog(slog, prod, vgetq_lane_f64(f2[0], 0) * vgetq_lane_f64(f2[0], 1));
					}
				}
			}
		}

		CloseLog(slog, prod);
		bufNK1[tid] = prod;
	}
}

#else

template<typename REAL>
TARGETNEO void BAYESIAN<REAL>::UpdateQNoAdmixNEO(int tid) { }

template<typename REAL>
TARGETNEO void BAYESIAN<REAL>::UpdateZNoAdmixNEO(int tid) { }

template<typename REAL>
template<bool fast_fp32>
TARGETNEO void BAYESIAN<REAL>::UpdateQMetroNEO(int tid) { }

template<typename REAL>
template<bool fast_fp32>
TARGETNEO void BAYESIAN<REAL>::UpdateZAdmixNEO(int tid) { }

template<typename REAL>
template<bool isadmix, bool fast_fp32>
TARGETNEO void BAYESIAN<REAL>::RecordNEO(int tid) { }

#endif
