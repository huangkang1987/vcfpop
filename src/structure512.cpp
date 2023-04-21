/* AVX512 Bayesian clustering functions */

#include "vcfpop.h"

template struct BAYESIAN<double>;
template struct BAYESIAN<float >;

template TARGET512 void BAYESIAN<double>::UpdateQNoAdmix512(int tid);
template TARGET512 void BAYESIAN<float >::UpdateQNoAdmix512(int tid);

template TARGET512 void BAYESIAN<double>::UpdateZNoAdmix512(int tid);
template TARGET512 void BAYESIAN<float >::UpdateZNoAdmix512(int tid);

template TARGET512 void BAYESIAN<double>::UpdateQMetro512<true >(int tid);
template TARGET512 void BAYESIAN<double>::UpdateQMetro512<false>(int tid);
template TARGET512 void BAYESIAN<float >::UpdateQMetro512<true >(int tid);
template TARGET512 void BAYESIAN<float >::UpdateQMetro512<false>(int tid);

template TARGET512 void BAYESIAN<double>::UpdateZAdmix512<true >(int tid);
template TARGET512 void BAYESIAN<double>::UpdateZAdmix512<false>(int tid);
template TARGET512 void BAYESIAN<float >::UpdateZAdmix512<true >(int tid);
template TARGET512 void BAYESIAN<float >::UpdateZAdmix512<false>(int tid);

template TARGET512 void BAYESIAN<double>::Record512<true , true >(int tid);
template TARGET512 void BAYESIAN<double>::Record512<true , false>(int tid);
template TARGET512 void BAYESIAN<double>::Record512<false, true >(int tid);
template TARGET512 void BAYESIAN<double>::Record512<false, false>(int tid);
template TARGET512 void BAYESIAN<float >::Record512<true , true >(int tid);
template TARGET512 void BAYESIAN<float >::Record512<true , false>(int tid);
template TARGET512 void BAYESIAN<float >::Record512<false, true >(int tid);
template TARGET512 void BAYESIAN<float >::Record512<false, false>(int tid);

#ifndef __aarch64__

#ifndef _GENO_READER512
template<typename REAL>
struct GENO_READER512
{
	__m512i vindex[4];						//Offset of 16 loci
	__m512i data[4];						//Readed bits
	__m512i msize;
	byte* pos;							//Current read pointer
	byte size;								//Number of bits a genotype id used
	byte nbits;								//Number of bits remaining in data

	TARGET512 GENO_READER512()
	{

	}

	TARGET512 GENO_READER512(int indid, int64 l, int64 num, BUCKET* bucket = NULL)
	{
		//set pos and size
		SetZero(this, 1);
		num = Min(num, 64);

		//set bucket from default bucket or assigned bucket
		if (bucket == NULL) bucket = &geno_bucket;

		LIST<OFFSET>& offset = bucket->offset;
		uint64 offset0 = offset[l].offset;
		pos = bucket->base_addr + offset0;

		int* vindex2 = (int*)vindex;
		OFFSET* offset1 = &offset.bucket[l];
		for (int i = 0; i < num; ++i)
			vindex2[i] = offset1[i].offset - offset0;

		size = bucket->offset[l].size;
		msize = _mm512_set1_epi32((1u << size) - 1u);

		if (indid)
		{
			//skip indid * size bits
			int nreadbits = indid * size;

			//pos move nreadbits / 8
			pos += nreadbits >> 3;

			//remain bits to read
			nreadbits &= 31;

			//read 32 bits to data and read remain bits from data
			REP(4) data[kk] = _mm512_srli_epi32(_mm512_i32gather_epi32(vindex[kk], (void*)pos, 1), nreadbits);
		
			pos += 4;

			//set nbits
			nbits = 32 - nreadbits;
		}
		else
		{
			//read 32 bits
			//GCC 11 bug here
			REP(4) data[kk] = _mm512_i32gather_epi32(vindex[kk], (void*)pos, 1);
			
			pos += 4;

			//set nbits
			nbits = 32;
		}
	}

	__forceinline
	TARGET512 void Read(__m512i* gid)
	{
		// if data is empty
		if (nbits < size) [[unlikely]]
		{
			// remain number of bits to read
			int rbits = size - nbits;
			__m512i tmask = _mm512_set1_epi32((1u << rbits) - 1u);

			// move nbits data to gid
			REP(4) gid[kk] = data[kk];

			// read 32 bits to data
			REP(4) data[kk] = _mm512_i32gather_epi32(vindex[kk], (void*)pos, 1);

			// read rbits from data and concate to higher bits in gid
			REP(4) gid[kk] = _mm512_or_si512(gid[kk], _mm512_slli_epi32(_mm512_and_si512(data[kk], tmask), nbits));

			//shift right
			REP(4) data[kk] = _mm512_srli_epi32(data[kk], rbits);

			nbits = 32 - rbits;

			pos += 4;
		}
		else [[likely]]
		{
			//read size bits
			REP(4) gid[kk] = _mm512_and_epi32(data[kk], msize);

			//shift right
			REP(4) data[kk] = _mm512_srli_epi32(data[kk], size);

			nbits -= size;
		}
	}
};
#endif

#ifndef _GENO_READERSSE
template<typename REAL>
struct GENO_READERSSE
{
	__m128i data[16];						//Readed bits
	__m128i msize;
	byte* pos;								//Current read pointer
	uint vindex[64];						//Offset of 16 loci
	byte size;								//Number of bits a genotype id used
	byte nbits;								//Number of bits remaining in data

	TARGETSSE GENO_READERSSE()
	{

	}

	TARGETSSE GENO_READERSSE(int indid, int64 l, int64 num, BUCKET* bucket = NULL)
	{
		//set pos and size
		SetZero(this, 1);
		num = Min(num, 64);

		//set bucket from default bucket or assigned bucket
		if (bucket == NULL) bucket = &geno_bucket;

		OFFSET* offset = &bucket->offset.bucket[l];
		uint64 offset0 = offset[0].offset;
		pos = (byte*)(bucket->base_addr + offset0);

		for (int i = 0; i < num; ++i)
			vindex[i] = offset[i].offset - offset0;

		size = offset[0].size;

		msize = _mm_set1_epi32((1u << size) - 1u);

		if (indid)
		{
			//skip indid * size bits
			int nreadbits = indid * size;

			//pos move nreadbits / 32
			pos += nreadbits >> 5;

			//remain bits to read
			nreadbits &= 31;

			//read 32 bits to data and read remain bits from data
			REP(16) data[kk] = _mm_set_epi32(*(uint*)(pos + vindex[3 + (kk << 2)]), *(uint*)(pos + vindex[2 + (kk << 2)]), *(uint*)(pos + vindex[1 + (kk << 2)]), *(uint*)(pos + vindex[0 + (kk << 2)]));

			REP(16) data[kk] = _mm_srli_epi32(data[kk], nreadbits);

			pos += 4;

			//set nbits
			nbits = 32 - nreadbits;
		}
		else
		{
			//read 32 bits
			REP(16) data[kk] = _mm_set_epi32(*(uint*)(pos + vindex[3 + (kk << 2)]), *(uint*)(pos + vindex[2 + (kk << 2)]), *(uint*)(pos + vindex[1 + (kk << 2)]), *(uint*)(pos + vindex[0 + (kk << 2)]));

			pos += 4;

			//set nbits
			nbits = 32;
		}
	}

	__forceinline
		TARGETSSE void Read(__m512i* gid2)
	{
		__m128i* gid = (__m128i*)gid2;
		// if data is empty
		if (nbits < size) [[unlikely]]
		{
			// move nbits data to gid
			memcpy(gid, data, 16 * sizeof(__m128i));

			// remain number of bits to read
			int rbits = size - nbits;
			__m128i tmask = _mm_set1_epi32((1u << rbits) - 1u);

			// read 32 bits to data
			REP(16) data[kk] = _mm_set_epi32(*(uint*)(pos + vindex[3 + (kk << 2)]), *(uint*)(pos + vindex[2 + (kk << 2)]), *(uint*)(pos + vindex[1 + (kk << 2)]), *(uint*)(pos + vindex[0 + (kk << 2)]));

			// read rbits from data and concate to higher bits in gid
			REP(16) gid[kk] = _mm_or_si128(gid[kk], _mm_slli_epi32(_mm_and_si128(data[kk], tmask), nbits));

			//shift right
			REP(16) data[kk] = _mm_srli_epi32(data[kk], rbits);

			pos += 4;

			nbits = 32 - rbits;
		}
		else [[likely]]
		{
			//read size bits
			REP(16) gid[kk] = _mm_and_si128(data[kk], msize);

			//shift right
			REP(16) data[kk] = _mm_srli_epi32(data[kk], size);

			nbits -= size;
		}
	}
};
#endif

__forceinline TARGET512 double __mm512_reduce_add_pd(__m512d v0)
{
	__m256d v1 = _mm256_add_pd(_mm512_extractf64x4_pd(v0, 0), _mm512_extractf64x4_pd(v0, 1));
	__m128d v2 = _mm_add_pd(_mm256_extractf128_pd(v1, 0), _mm256_extractf128_pd(v1, 1));
	return simd_f64(v2, 0) + simd_f64(v2, 1);
}

__forceinline TARGET512 double __mm512_reduce_mul_pd(__m512d v0)
{
	__m256d v1 = _mm256_mul_pd(_mm512_extractf64x4_pd(v0, 0), _mm512_extractf64x4_pd(v0, 1));
	__m128d v2 = _mm_mul_pd(_mm256_extractf128_pd(v1, 0), _mm256_extractf128_pd(v1, 1));
	return simd_f64(v2, 0) * simd_f64(v2, 1);
}

__forceinline TARGET512 float __mm512_reduce_add_ps(__m512 v0)
{
	__m256 v1 = _mm256_mul_ps(_mm512_extractf32x8_ps(v0, 0), _mm512_extractf32x8_ps(v0, 1));
	__m128 v2 = _mm_mul_ps(_mm256_extractf128_ps(v1, 0), _mm256_extractf128_ps(v1, 1));
	__m128 v3 = _mm_mul_ps(v2, _mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(v2), 8)));
	return simd_f64(v3, 0) + simd_f64(v3, 1);
}

__forceinline TARGET512 float __mm512_reduce_mul_ps(__m512 v0)
{
	__m256 v1 = _mm256_add_ps(_mm512_extractf32x8_ps(v0, 0), _mm512_extractf32x8_ps(v0, 1));
	__m128 v2 = _mm_add_ps(_mm256_extractf128_ps(v1, 0), _mm256_extractf128_ps(v1, 1));
	__m128 v3 = _mm_add_ps(v2, _mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(v2), 8)));
	return simd_f64(v3, 0) * simd_f64(v3, 1);
}

/* Update a priori ancetral proportion for non-admix model */
template<typename REAL>
TARGET512 void BAYESIAN<REAL>::UpdateQNoAdmix512(int tid)
{
	if (tid == -1)
	{
		SetZero(Q, N * K);
		OpenLog((int64*)bufNK1, bufNK2, N * K * structure_nsubthread);

		//add priori probability
		double* buf1 = bufNK1, * buf2 = bufNK2;
		if (locpriori) for (int i = 0; i < N; ++i, buf1 += K, buf2 += K)
		{
			if (ainds[i]->vt == 0) continue;
			ChargeLog((int64*)buf1, buf2, Gamma + ainds[i]->popid * K, K);
		}

		//////////////////////////////////////////////////////////

		SetZero((int64*)l_atomic, 32);

		UpdateQNoAdmix512(0);

		//avoid thread-conflict
		for (int i = 1; i < structure_nsubthread; ++i)
			Add(bufNK1, bufNK1 + N * K * i, N * K);

		buf1 = bufNK1;
		REAL* q = Q;
		RNG<REAL> rng(seed + m, RNG_SALT_UPDATEQ);//checked
		for (int i = 0; i < N; ++i, buf1 += K, q += K)
		{
			if (ainds[i]->vt == 0) continue;
			ushort k2 = (ushort)rng.PolyLog(buf1, K);
			q[k2] = 1;
			Z[i] = k2;
		}
		return;
	}

	static __m512d maskoned = _mm512_set1_pd(1.0);
	static __m256  maskones256 = _mm256_set1_ps(1.0);
	static __m512i mask00 = _mm512_set1_epi32(0);
	static __m512i mask24 = _mm512_set1_epi32(0xFFFFFF);
	static __m512i mask16 = _mm512_set1_epi32(0xFFFF);
	static __m512i mask01 = _mm512_set1_epi32(1);
	static __m512i mask02 = _mm512_set1_epi32(2);
	static int PT_PLOIDYxNALLELES[150] = 									//Pattern index to ploidy level
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

				__m512i gtab[4];
				for (int64 k = 0, lt = l; k < 4; ++k, lt += 16)
					gtab[k] = _mm512_set_epi32(
						GetLocTabDiff(lt + 15), GetLocTabDiff(lt + 14), GetLocTabDiff(lt + 13), GetLocTabDiff(lt + 12),
						GetLocTabDiff(lt + 11), GetLocTabDiff(lt + 10), GetLocTabDiff(lt + 9), GetLocTabDiff(lt + 8),
						GetLocTabDiff(lt + 7), GetLocTabDiff(lt + 6), GetLocTabDiff(lt + 5), GetLocTabDiff(lt + 4),
						GetLocTabDiff(lt + 3), GetLocTabDiff(lt + 2), GetLocTabDiff(lt + 1), GetLocTabDiff(lt + 0));

				ushort lmask[4];
				byte* lmask1 = (byte*)lmask;
				uint64& lmask0 = *(uint64*)lmask;
				lmask0 = lend - l >= 64 ? 0xFFFFFFFFFFFFFFFF : (1ull << (lend - l)) - 1ull;

#ifdef __GNUC__
				GENO_READERSSE<REAL> rt(0, l, lend - l);
#else
				GENO_READER512<REAL> rt(0, l, lend - l);
#endif

				__m512i oindex[4]; __m256i* oindex2 = (__m256i*)oindex;
				__m512i oindex0 = _mm512_set1_epi64(allele_freq_offset[l]);

				REP(8) oindex2[kk] = _mm512_cvtepi64_epi32(_mm512_sub_epi64(_mm512_maskz_loadu_epi64(lmask1[kk], &allele_freq_offset[l + (kk << 3)]), oindex0));
				REP(4) oindex[kk] = _mm512_mask_mov_epi32(mask00, lmask[kk], oindex[kk]);

				__m512i gtaddr[4], gtlow[4], gtploidy[4], gtals[4], allele[4];
				__m512d freq[8];
				__m256i* allele2 = (__m256i*) & allele;
				__mmask16 typed[4];
				byte* typed2 = (byte*)typed;

				int64* buf1 = (int64*)bufNK1 + N * K * tid2; double* buf2 = bufNK2 + N * K * tid2;

				for (int i = 0; i < N; ++i, buf1 += K, buf2 += K)
				{
					rt.Read(gtaddr);

					REP(4) gtaddr[kk] = _mm512_add_epi32(gtab[kk], _mm512_slli_epi32(gtaddr[kk], 2));

					REP(4) gtlow[kk] = _mm512_i32gather_epi32(gtaddr[kk], (int*)gtab_base, 1);

					REP(4) gtploidy[kk] = _mm512_srli_epi32(gtlow[kk], 24);

					REP(4) gtlow[kk] = _mm512_and_epi32(gtlow[kk], mask24);

					REP(4) gtploidy[kk] = _mm512_mask_i32gather_epi32(mask00, lmask[kk], gtploidy[kk], PT_PLOIDYxNALLELES, sizeof(int));

					REP(4) gtals[kk] = _mm512_add_epi32(gtaddr[kk], gtlow[kk]);

					int maxv = maxploidy == minploidy ? maxploidy :
						_mm512_reduce_max_epi32(_mm512_max_epi32(
							_mm512_max_epi32(gtploidy[0], gtploidy[1]),
							_mm512_max_epi32(gtploidy[2], gtploidy[3])));

					__m512i ai = _mm512_set1_epi32(0);
					for (int a = 0; a < maxv; ++a, ai = _mm512_add_epi32(ai, mask01))
					{
						REP(4)
						{
							typed[kk] = _mm512_cmplt_epi32_mask(ai, gtploidy[kk]);

							allele[kk] = _mm512_mask_i32gather_epi32(mask00, typed[kk], gtals[kk], (void*)gtab_base, 1);

							gtals[kk] = _mm512_add_epi32(gtals[kk], mask02);

							allele[kk] = _mm512_and_epi32(allele[kk], mask16);

							allele[kk] = _mm512_add_epi32(allele[kk], oindex[kk]);
						}

						REAL* p2 = p;
						for (int k = 0; k < K; ++k, p2 += KT)
						{
							if constexpr (sizeof(REAL) == 8)
								REP(8) freq[kk] = _mm512_mask_i32gather_pd(maskoned, typed2[kk], allele2[kk], p2, sizeof(double));
							else
								REP(8) freq[kk] = _mm512_cvtps_pd(_mm256_mmask_i32gather_ps(maskones256, typed2[kk], allele2[kk], p2, sizeof(float)));

							REP(4) freq[kk] = _mm512_mul_pd(freq[kk], freq[kk + 4]);
							REP(2) freq[kk] = _mm512_mul_pd(freq[kk], freq[kk + 2]);
							REP(1) freq[kk] = _mm512_mul_pd(freq[kk], freq[kk + 1]);

							__m128d* freq2 = (__m128d*)freq;
							REP(2) freq2[kk] = _mm_mul_pd(freq2[kk], freq2[kk + 2]);
							REP(1) freq2[kk] = _mm_mul_pd(freq2[kk], freq2[kk + 1]);

							ChargeLog(buf1[k], buf2[k], simp_f64(freq, 0) * simp_f64(freq, 1));
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
TARGET512 void BAYESIAN<REAL>::UpdateZNoAdmix512(int tid)
{
	if (tid == -1)
	{
		SetZero(Mi, N * K);
		SetZero(Ni, K * KT);

		for (int i = 0; i < N; ++i)
			Mi[i * K + Z[i]] = ainds[i]->vt;

		//////////////////////////////////////////////////////////

		SetZero((int64*)l_atomic, 32);

		UpdateZNoAdmix512(0);

		return;
	}

	static __m512i maskff = _mm512_set1_epi32(0xFFFFFFFF);
	static __m512i mask00 = _mm512_set1_epi32(0);
	static __m512i mask24 = _mm512_set1_epi32(0xFFFFFF);
	static __m512i mask16 = _mm512_set1_epi32(0xFFFF);
	static __m512i mask01 = _mm512_set1_epi32(1);
	static __m512i mask02 = _mm512_set1_epi32(2);
	static int PT_PLOIDYxNALLELES[150] = 									//Pattern index to ploidy level
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

				__m512i gtab[4];
				for (int64 k = 0, lt = l; k < 4; ++k, lt += 16)
					gtab[k] = _mm512_set_epi32(
						GetLocTabDiff(lt + 15), GetLocTabDiff(lt + 14), GetLocTabDiff(lt + 13), GetLocTabDiff(lt + 12),
						GetLocTabDiff(lt + 11), GetLocTabDiff(lt + 10), GetLocTabDiff(lt + 9), GetLocTabDiff(lt + 8),
						GetLocTabDiff(lt + 7), GetLocTabDiff(lt + 6), GetLocTabDiff(lt + 5), GetLocTabDiff(lt + 4),
						GetLocTabDiff(lt + 3), GetLocTabDiff(lt + 2), GetLocTabDiff(lt + 1), GetLocTabDiff(lt + 0));

				ushort lmask[4];
				byte* lmask1 = (byte*)lmask;
				uint64& lmask0 = *(uint64*)lmask;
				lmask0 = lend - l >= 64 ? 0xFFFFFFFFFFFFFFFF : (1ull << (lend - l)) - 1ull;

#ifdef __GNUC__
				GENO_READERSSE<REAL> rt(0, l, lend - l);
#else
				GENO_READER512<REAL> rt(0, l, lend - l);
#endif

				__m512i oindex[4]; __m256i* oindex2 = (__m256i*)oindex;
				__m512i oindex0 = _mm512_set1_epi64(allele_freq_offset[l]);

				REP(8) oindex2[kk] = _mm512_cvtepi64_epi32(_mm512_sub_epi64(_mm512_maskz_loadu_epi64(lmask1[kk], &allele_freq_offset[l + (kk << 3)]), oindex0));
				REP(4) oindex[kk] = _mm512_mask_mov_epi32(mask00, lmask[kk], oindex[kk]);

				__m512i gtaddr[4], gtlow[4], gtploidy[4], gtals[4], allele[4], data[4];
				__mmask16 typed[4];
				int* allele2 = (int*)allele, * data2 = (int*)data;

				for (int i = 0; i < N; i++)
				{
					rt.Read(gtaddr);

					REP(4) gtaddr[kk] = _mm512_add_epi32(gtab[kk], _mm512_slli_epi32(gtaddr[kk], 2));

					REP(4) gtlow[kk] = _mm512_mask_i32gather_epi32(mask00, lmask[kk], gtaddr[kk], (int*)gtab_base, 1);

					REP(4) gtploidy[kk] = _mm512_srli_epi32(gtlow[kk], 24);

					REP(4) gtlow[kk] = _mm512_and_epi32(gtlow[kk], mask24);

					REP(4) gtploidy[kk] = _mm512_i32gather_epi32(gtploidy[kk], PT_PLOIDYxNALLELES, sizeof(int));

					REP(4) gtals[kk] = _mm512_add_epi32(gtaddr[kk], gtlow[kk]);

					int maxv = maxploidy == minploidy ? maxploidy :
						_mm512_reduce_max_epi32(_mm512_max_epi32(
							_mm512_max_epi32(gtploidy[0], gtploidy[1]),
							_mm512_max_epi32(gtploidy[2], gtploidy[3])));

					int* ni = Ni + Z[i] * KT + o;
					__m512i ai = _mm512_set1_epi32(0);
					for (int a = 0; a < maxv; ++a, ai = _mm512_add_epi32(ai, mask01))
					{
						REP(4)
						{
							typed[kk] = _mm512_cmplt_epi32_mask(ai, gtploidy[kk]);

							allele[kk] = _mm512_mask_i32gather_epi32(mask00, typed[kk], gtals[kk], (void*)gtab_base, 1);

							data[kk] = _mm512_mask_blend_epi32(typed[kk], mask00, maskff);

							allele[kk] = _mm512_and_epi32(allele[kk], mask16);

							allele[kk] = _mm512_add_epi32(allele[kk], oindex[kk]);

							gtals[kk] = _mm512_add_epi32(gtals[kk], mask02);
						}

						REP(64) ni[allele2[kk]] -= data2[kk];
					}
				}
			}
		}
	}
}

/* Update a priori ancetral proportion for non-admix model */
template<typename REAL>
template<bool fast_fp32>
TARGET512 void BAYESIAN<REAL>::UpdateQMetro512(int tid)
{
	if (tid == -1)
	{
		RNG<REAL> rng(seed + m, RNG_SALT_UPDATEQ);//checked
		REAL* bufi = (REAL*)bufNK1;
		REAL* q = NULL;

		for (int i = 0; i < N; ++i, bufi += K)
		{
			if (ainds[i]->vt == 0) continue;
			if (locpriori) rng.Dirichlet(bufi, AlphaLocal + ainds[i]->popid * K, K);
			else           rng.Dirichlet(bufi, Alpha, K);
		}

		OpenLog((int64*)bufN1, bufN2, N * structure_nsubthread);

		//////////////////////////////////////////////////////////

		SetZero((int64*)l_atomic, 32);

		g_fastsingle_val == 1 ? UpdateQMetro512<true >(0) : UpdateQMetro512<false>(0);

		//avoid thread-conflict
		for (int i = 1; i < structure_nsubthread; ++i)
			Add(bufN1, bufN1 + N * i, N);

		bufi = (REAL*)bufNK1; q = Q;
		for (int i = 0; i < N; ++i, q += K, bufi += K)
		{
			if (ainds[i]->vt == 0) continue;
			if (bufN1[i] >= NZERO || rng.Uniform() < exp(bufN1[i]))
				SetVal(q, bufi, K);
		}
		return;
	}

	static __m512d maskoned = _mm512_set1_pd(1.0);
	static __m512  maskones = _mm512_set1_ps(1.0);
	static __m512i mask00 = _mm512_set1_epi32(0);
	static __m512i mask24 = _mm512_set1_epi32(0xFFFFFF);
	static __m512i mask16 = _mm512_set1_epi32(0xFFFF);
	static __m512i mask01 = _mm512_set1_epi32(1);
	static __m512i mask02 = _mm512_set1_epi32(2);
	static int PT_PLOIDYxNALLELES[150] = 									//Pattern index to ploidy level
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

				__m512i gtab[4];
				for (int64 k = 0, lt = l; k < 4; ++k, lt += 16)
					gtab[k] = _mm512_set_epi32(
						GetLocTabDiff(lt + 15), GetLocTabDiff(lt + 14), GetLocTabDiff(lt + 13), GetLocTabDiff(lt + 12),
						GetLocTabDiff(lt + 11), GetLocTabDiff(lt + 10), GetLocTabDiff(lt + 9), GetLocTabDiff(lt + 8),
						GetLocTabDiff(lt + 7), GetLocTabDiff(lt + 6), GetLocTabDiff(lt + 5), GetLocTabDiff(lt + 4),
						GetLocTabDiff(lt + 3), GetLocTabDiff(lt + 2), GetLocTabDiff(lt + 1), GetLocTabDiff(lt + 0));

				ushort lmask[4];
				byte* lmask1 = (byte*)lmask;
				uint64& lmask0 = *(uint64*)lmask;
				lmask0 = lend - l >= 64 ? 0xFFFFFFFFFFFFFFFF : (1ull << (lend - l)) - 1ull;

#ifdef __GNUC__
				GENO_READERSSE<REAL> rt(0, l, lend - l);
#else
				GENO_READER512<REAL> rt(0, l, lend - l);
#endif

				__m512i oindex[4]; __m256i* oindex2 = (__m256i*)oindex;
				__m512i oindex0 = _mm512_set1_epi64(allele_freq_offset[l]);

				REP(8) oindex2[kk] = _mm512_cvtepi64_epi32(_mm512_sub_epi64(_mm512_maskz_loadu_epi64(lmask1[kk], &allele_freq_offset[l + (kk << 3)]), oindex0));
				REP(4) oindex[kk] = _mm512_mask_mov_epi32(mask00, lmask[kk], oindex[kk]);

				__m512i gtaddr[4], gtlow[4], gtploidy[4], gtals[4], allele[4];
				__m512d f1[8], f2[8];
				__m512 f1s[4], f2s[4];
				__m256i* allele2 = (__m256i*) & allele;
				__mmask16 typed[4];
				byte* typed2 = (byte*)typed;

				REAL* bufi = (REAL*)bufNK1, * q = Q;
				//avoid thread-conflict
				double* buf1 = bufN1 + N * tid2, * buf2 = bufN2 + N * tid2;

				for (int i = 0; i < N; ++i, q += K, bufi += K, buf1++, buf2++)
				{
					rt.Read(gtaddr);

					REP(4) gtaddr[kk] = _mm512_add_epi32(gtab[kk], _mm512_slli_epi32(gtaddr[kk], 2));

					REP(4) gtlow[kk] = _mm512_i32gather_epi32(gtaddr[kk], (int*)gtab_base, 1);

					REP(4) gtploidy[kk] = _mm512_srli_epi32(gtlow[kk], 24);

					REP(4) gtlow[kk] = _mm512_and_epi32(gtlow[kk], mask24);

					REP(4) gtploidy[kk] = _mm512_mask_i32gather_epi32(mask00, lmask[kk], gtploidy[kk], PT_PLOIDYxNALLELES, sizeof(int));

					REP(4) gtals[kk] = _mm512_add_epi32(gtaddr[kk], gtlow[kk]);

					int maxv = maxploidy == minploidy ? maxploidy :
						_mm512_reduce_max_epi32(_mm512_max_epi32(
							_mm512_max_epi32(gtploidy[0], gtploidy[1]),
							_mm512_max_epi32(gtploidy[2], gtploidy[3])));

					__m512i ai = _mm512_set1_epi32(0);

					for (int a = 0; a < maxv; ++a, ai = _mm512_add_epi32(ai, mask01))
					{
						REP(4)
						{
							typed[kk] = _mm512_cmplt_epi32_mask(ai, gtploidy[kk]);

							allele[kk] = _mm512_mask_i32gather_epi32(mask00, typed[kk], gtals[kk], (void*)gtab_base, 1);

							gtals[kk] = _mm512_add_epi32(gtals[kk], mask02);

							allele[kk] = _mm512_and_epi32(allele[kk], mask16);

							allele[kk] = _mm512_add_epi32(allele[kk], oindex[kk]);
						}

						if constexpr (sizeof(REAL) == 8 || !fast_fp32)
						{
							REP(8) f1[kk] = _mm512_setzero_pd();
							REP(8) f2[kk] = _mm512_setzero_pd();
						}
						else
						{
							REP(4) f1s[kk] = _mm512_setzero_ps();
							REP(4) f2s[kk] = _mm512_setzero_ps();
						}

						REAL* p2 = p;
						for (int k = 0; k < K; ++k, p2 += KT)
						{
							if constexpr (sizeof(REAL) == 8)
							{
								__m512d pp[8], qq = _mm512_set1_pd(q[k]), ii = _mm512_set1_pd(bufi[k]);

								REP(8) pp[kk] = _mm512_i32gather_pd(allele2[kk], p2, sizeof(double));
								REP(8) f1[kk] = _mm512_mask_add_pd(maskoned, typed2[kk], f1[kk], _mm512_mul_pd(pp[kk], ii));
								REP(8) f2[kk] = _mm512_mask_add_pd(maskoned, typed2[kk], f2[kk], _mm512_mul_pd(pp[kk], qq));
							}
							else if constexpr (fast_fp32)
							{
								__m512 pp[4], qq = _mm512_set1_ps(q[k]), ii = _mm512_set1_ps(bufi[k]);

								REP(4) pp[kk] = _mm512_i32gather_ps(allele[kk], p2, sizeof(float));
								REP(4) f1s[kk] = _mm512_mask_add_ps(maskones, typed[kk], f1s[kk], _mm512_mul_ps(pp[kk], ii));
								REP(4) f2s[kk] = _mm512_mask_add_ps(maskones, typed[kk], f2s[kk], _mm512_mul_ps(pp[kk], qq));
							}
							else
							{
								__m512d pp[8], qq = _mm512_set1_pd(q[k]), ii = _mm512_set1_pd(bufi[k]);
								__m512 ps[4]; __m256* ps2 = (__m256*)ps;

								REP(4) ps[kk] = _mm512_i32gather_ps(allele[kk], p2, sizeof(float));
								REP(8) pp[kk] = _mm512_cvtps_pd(ps2[kk]);
								REP(8) f1[kk] = _mm512_mask_add_pd(maskoned, typed2[kk], f1[kk], _mm512_mul_pd(pp[kk], ii));
								REP(8) f2[kk] = _mm512_mask_add_pd(maskoned, typed2[kk], f2[kk], _mm512_mul_pd(pp[kk], qq));
							}
						}

						if constexpr (sizeof(REAL) == 4 && fast_fp32)
						{
							__m256* f1t = (__m256*)f1s;
							REP(4) f1s[kk] = _mm512_div_ps(f1s[kk], f2s[kk]);
							//REP(4) f1s[kk] = _mm512_mask_mov_ps(maskones, typed[kk], f1s[kk]);
							REP(8) f1[kk] = _mm512_cvtps_pd(f1t[kk]);
						}
						else
						{
							REP(8) f1[kk] = _mm512_div_pd(f1[kk], f2[kk]);
							//REP(8) f1[kk] = _mm512_mask_mov_pd(maskoned, typed2[kk], f1[kk]);
						}

						for (int KK = sizeof(f1) / sizeof(f1[0]) / 2; KK >= 1; KK >>= 1)
							REP(KK) f1[kk] = _mm512_mul_pd(f1[kk], f1[kk + KK]);

						ChargeLog(*(int64*)buf1, *buf2, __mm512_reduce_mul_pd(f1[0]));
					}
				}
			}
		}

		//avoid thread-conflict
		CloseLog((int64*)bufN1 + N * tid2, bufN2 + N * tid2, N);
	}
} 

/* Update individual or allele origin when ancetral proportion */
template<typename REAL>
template<bool fast_fp32>
TARGET512 void BAYESIAN<REAL>::UpdateZAdmix512(int tid)
{
	if (tid == -1)
	{
		SetZero(Mi, N * K * structure_nsubthread);
		SetZero(Ni, K * KT);

		//////////////////////////////////////////////////////////

		SetZero((int64*)l_atomic, 32);

		g_fastsingle_val == 1 ? UpdateZAdmix512<true >(0) : UpdateZAdmix512<false>(0);

		//avoid thread-conflict
		for (int i = 1; i < structure_nsubthread; ++i)
			Add(Mi, Mi + N * K * i, N * K);

		return;
	}

	static __m512i maskKT_64 = _mm512_set1_epi64(KT);
	static __m512i mask00 = _mm512_set1_epi32(0);
	static __m256i mask00_256 = _mm256_set1_epi32(0);
	static __m512i mask24 = _mm512_set1_epi32(0xFFFFFF);
	static __m512i mask16 = _mm512_set1_epi32(0xFFFF);
	static __m512i mask01 = _mm512_set1_epi32(1);
	static __m512i mask02 = _mm512_set1_epi32(2);
	static __m512d minfreqd = _mm512_set1_pd(MIN_FREQ);
	static __m512 minfreqs = _mm512_set1_ps(MIN_FREQ);
	static int PT_PLOIDYxNALLELES[150] = 									//Pattern index to ploidy level
	{ 0, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	atomic<int> thread_counter = 0;
#pragma omp parallel num_threads(structure_nsubthread)
	{	
		int tid2 = thread_counter.fetch_add(1);

		//64 * K  * sizeof(double)
		__m512d* bufkd = (__m512d*)Align64((byte*)bufNK2 + (Max(N, 64) * K * sizeof(double) + 63) * tid2);
		__m512* bufks = (__m512*)Align64((byte*)bufNK2 + (Max(N, 64) * K * sizeof(double) + 63) * tid2);

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

				__m512i gtab[4];
				for (int64 k = 0, lt = l; k < 4; ++k, lt += 16)
					gtab[k] = _mm512_set_epi32(
						GetLocTabDiff(lt + 15), GetLocTabDiff(lt + 14), GetLocTabDiff(lt + 13), GetLocTabDiff(lt + 12),
						GetLocTabDiff(lt + 11), GetLocTabDiff(lt + 10), GetLocTabDiff(lt + 9), GetLocTabDiff(lt + 8),
						GetLocTabDiff(lt + 7), GetLocTabDiff(lt + 6), GetLocTabDiff(lt + 5), GetLocTabDiff(lt + 4),
						GetLocTabDiff(lt + 3), GetLocTabDiff(lt + 2), GetLocTabDiff(lt + 1), GetLocTabDiff(lt + 0));

				ushort lmask[4];
				byte* lmask1 = (byte*)lmask;
				uint64& lmask0 = *(uint64*)lmask;
				lmask0 = lend - l >= 64 ? 0xFFFFFFFFFFFFFFFF : (1ull << (lend - l)) - 1ull;

#ifdef __GNUC__
				GENO_READERSSE<REAL> rt(0, l, lend - l);
#else
				GENO_READER512<REAL> rt(0, l, lend - l);
#endif

				__m512i oindex[4]; __m256i* oindex2 = (__m256i*)oindex;
				__m512i oindex0 = _mm512_set1_epi64(allele_freq_offset[l]);

				REP(8) oindex2[kk] = _mm512_cvtepi64_epi32(_mm512_sub_epi64(_mm512_maskz_loadu_epi64(lmask1[kk], &allele_freq_offset[l + (kk << 3)]), oindex0));
				REP(4) oindex[kk] = _mm512_mask_mov_epi32(mask00, lmask[kk], oindex[kk]);

				__m512i gtaddr[4], gtlow[4], gtploidy[4], gtals[4], allele[4], k2[8];
				__m256i* allele2 = (__m256i*) & allele;
				__mmask16 typed[4];
				uint64& typed0 = *(uint64*)typed;
				byte* typed2 = (byte*)typed;
				uint64* k22 = (uint64*)&k2;

				RNG512<double> rngd; RNG512<float > rngs;

				if constexpr (sizeof(REAL) == 8 || !fast_fp32)
					new (&rngd) RNG512<double>(seed + m * L + l, RNG_SALT_UPDATEZ);
				else
					new (&rngs) RNG512<float >(seed + m * L + l, RNG_SALT_UPDATEZ);

				REAL* q = Q;
				//avoid thread-conflict
				int64* mi = Mi + N * K * tid2;

				for (int i = 0; i < N; ++i, q += K, mi += K)
				{
					rt.Read(gtaddr);

					REP(4) gtaddr[kk] = _mm512_add_epi32(gtab[kk], _mm512_slli_epi32(gtaddr[kk], 2));

					REP(4) gtlow[kk] = _mm512_i32gather_epi32(gtaddr[kk], (int*)gtab_base, 1);

					REP(4) gtploidy[kk] = _mm512_srli_epi32(gtlow[kk], 24);

					REP(4) gtlow[kk] = _mm512_and_epi32(gtlow[kk], mask24);

					REP(4) gtploidy[kk] = _mm512_mask_i32gather_epi32(mask00, lmask[kk], gtploidy[kk], PT_PLOIDYxNALLELES, sizeof(int));

					REP(4) gtals[kk] = _mm512_add_epi32(gtaddr[kk], gtlow[kk]);

					int maxv = maxploidy;

					__m512i ai = _mm512_set1_epi32(0);
					for (int a = 0; a < maxv; ++a, ai = _mm512_add_epi32(ai, mask01))
					{
						REP(4)
						{
							typed[kk] = _mm512_cmplt_epi32_mask(ai, gtploidy[kk]);

							allele[kk] = _mm512_mask_i32gather_epi32(mask00, typed[kk], gtals[kk], (void*)gtab_base, 1);

							gtals[kk] = _mm512_add_epi32(gtals[kk], mask02);

							allele[kk] = _mm512_and_epi32(allele[kk], mask16);

							allele[kk] = _mm512_add_epi32(allele[kk], oindex[kk]);
						}

						REAL* p2 = p;
						for (int k = 0; k < K; ++k, p2 += KT)
						{
							__m512d qd; __m512 qs;

							if constexpr (sizeof(REAL) == 8)
							{
								qd = _mm512_set1_pd(q[k]);

								REP(8) bufkd[kk + (k << 3)] = _mm512_i32gather_pd(allele2[kk], p2, sizeof(double));
								REP(8) bufkd[kk + (k << 3)] = _mm512_mul_pd(qd, bufkd[kk + (k << 3)]);
								REP(8) bufkd[kk + (k << 3)] = _mm512_mask_add_pd(_mm512_castsi512_pd(mask00), typed2[kk], minfreqd, bufkd[kk + (k << 3)]);
							}
							else if constexpr (fast_fp32)
							{
								qs = _mm512_set1_ps(q[k]);

								REP(4) bufks[kk + (k << 2)] = _mm512_i32gather_ps(allele[kk], p2, sizeof(float));
								REP(4) bufks[kk + (k << 2)] = _mm512_mul_ps(qs, bufks[kk + (k << 2)]);
								REP(4) bufks[kk + (k << 2)] = _mm512_mask_add_ps(_mm512_castsi512_ps(mask00), typed[kk], minfreqs, bufks[kk + (k << 2)]);
							}
							else
							{
								qd = _mm512_set1_pd(q[k]);

								REP(8) bufkd[kk + (k << 3)] = _mm512_cvtps_pd(_mm256_i32gather_ps(p2, allele2[kk], sizeof(float)));
								REP(8) bufkd[kk + (k << 3)] = _mm512_mul_pd(qd, bufkd[kk + (k << 3)]);
								REP(8) bufkd[kk + (k << 3)] = _mm512_mask_add_pd(_mm512_castsi512_pd(mask00), typed2[kk], minfreqd, bufkd[kk + (k << 3)]);
							}
						}

						//draw cluster for each allele copy
						if constexpr (sizeof(REAL) == 8 || !fast_fp32)
							rngd.Poly<64>(bufkd, K, k2);
						else
							rngs.Poly<64>(bufks, K, k2);

						//Update Mi
						uint64 typedx = typed0;
						REP(64) { mi[k22[kk]] += typedx & 1; typedx >>= 1; }

						REP(8) k2[kk] = _mm512_add_epi64(_mm512_mullo_epi64(k2[kk], maskKT_64), _mm512_cvtepi32_epi64(allele2[kk]));

						//ni[k2 * KT + als[a]]++;
						typedx = typed0;
						REP(64) { ni[k22[kk]] += typedx & 1;  typedx >>= 1; }
					}
				}
			}
		}
	}
}

/* Record updated MCMC parameters */
template<typename REAL>
template<bool isadmix, bool fast_fp32>
TARGET512 void BAYESIAN<REAL>::Record512(int tid)
{
	if (tid == -1) 
	{
		Add(MiSum, Mi, N * K);

		//////////////////////////////////////////////////////////

		SetZero((int64*)l_atomic, 32);
		
        switch ((int)binaryq * 10 + g_fastsingle_val)
        {
        case 01: BAYESIAN<REAL>::Record512<true , true >(0); break;
        case 02: BAYESIAN<REAL>::Record512<true , false>(0); break;
        case 11: BAYESIAN<REAL>::Record512<false, true >(0); break;
        case 12: BAYESIAN<REAL>::Record512<false, false>(0); break;
        }

		bufNK1[0] = Sum(bufNK1, structure_nsubthread);
		return;
	}

	static __m512d maskoned = _mm512_set1_pd(1.0);
	static __m256  maskones256 = _mm256_set1_ps(1.0);
	static __m512  maskones = _mm512_set1_ps(1.0);
	static __m512i mask00 = _mm512_set1_epi32(0);
	static __m512i mask24 = _mm512_set1_epi32(0xFFFFFF);
	static __m512i mask16 = _mm512_set1_epi32(0xFFFF);
	static __m512i mask01 = _mm512_set1_epi32(1);
	static __m512i mask02 = _mm512_set1_epi32(2);
	static int PT_PLOIDYxNALLELES[150] = 									//Pattern index to ploidy level
	{ 0, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	atomic<int> thread_counter = 0;
#pragma omp parallel num_threads(structure_nsubthread)
	{
		int tid2 = thread_counter.fetch_add(1);

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

				__m512i gtab[4];
				for (int64 k = 0, lt = l; k < 4; ++k, lt += 16)
					gtab[k] = _mm512_set_epi32(
						GetLocTabDiff(lt + 15), GetLocTabDiff(lt + 14), GetLocTabDiff(lt + 13), GetLocTabDiff(lt + 12),
						GetLocTabDiff(lt + 11), GetLocTabDiff(lt + 10), GetLocTabDiff(lt + 9), GetLocTabDiff(lt + 8),
						GetLocTabDiff(lt + 7), GetLocTabDiff(lt + 6), GetLocTabDiff(lt + 5), GetLocTabDiff(lt + 4),
						GetLocTabDiff(lt + 3), GetLocTabDiff(lt + 2), GetLocTabDiff(lt + 1), GetLocTabDiff(lt + 0));

				ushort lmask[4];
				byte* lmask1 = (byte*)lmask;
				uint64& lmask0 = *(uint64*)lmask;
				lmask0 = lend - l >= 64 ? 0xFFFFFFFFFFFFFFFF : (1ull << (lend - l)) - 1ull;

#ifdef __GNUC__
				GENO_READERSSE<REAL> rt(0, l, lend - l);
#else
				GENO_READER512<REAL> rt(0, l, lend - l);
#endif

				__m512i oindex[4]; __m256i* oindex2 = (__m256i*)oindex;
				__m512i oindex0 = _mm512_set1_epi64(allele_freq_offset[l]);

				REP(8) oindex2[kk] = _mm512_cvtepi64_epi32(_mm512_sub_epi64(_mm512_maskz_loadu_epi64(lmask1[kk], &allele_freq_offset[l + (kk << 3)]), oindex0));
				REP(4) oindex[kk] = _mm512_mask_mov_epi32(mask00, lmask[kk], oindex[kk]);

				__m512i gtaddr[4], gtlow[4], gtploidy[4], gtals[4], allele[4];
				__m512d f2[8]; __m512 f2s[4];
				__m256i* allele2 = (__m256i*)&allele;
				__mmask16 typed[4]; byte* typed2 = (byte*)typed;

				REAL* q = Q;

				for (int i = 0; i < N; ++i, q += K)
				{
					rt.Read(gtaddr);

					REP(4) gtaddr[kk] = _mm512_add_epi32(gtab[kk], _mm512_slli_epi32(gtaddr[kk], 2));

					REP(4) gtlow[kk] = _mm512_mask_i32gather_epi32(mask00, lmask[kk], gtaddr[kk], (int*)gtab_base, 1);

					REP(4) gtploidy[kk] = _mm512_srli_epi32(gtlow[kk], 24);

					REP(4) gtlow[kk] = _mm512_and_epi32(gtlow[kk], mask24);

					REP(4) gtploidy[kk] = _mm512_i32gather_epi32(gtploidy[kk], PT_PLOIDYxNALLELES, sizeof(int));

					REP(4) gtals[kk] = _mm512_add_epi32(gtaddr[kk], gtlow[kk]);

					int maxv = maxploidy == minploidy ? maxploidy :
						_mm512_reduce_max_epi32(_mm512_max_epi32(
							_mm512_max_epi32(gtploidy[0], gtploidy[1]),
							_mm512_max_epi32(gtploidy[2], gtploidy[3])));

					__m512i ai = _mm512_set1_epi32(0);

					for (int a = 0; a < maxv; ++a, ai = _mm512_add_epi32(ai, mask01))
					{
						REP(4) typed[kk] = _mm512_cmplt_epi32_mask(ai, gtploidy[kk]);

						REP(4) allele[kk] = _mm512_and_epi32(_mm512_mask_i32gather_epi32(mask00, typed[kk], gtals[kk], (void*)gtab_base, 1), mask16);

						REP(4) gtals[kk] = _mm512_add_epi32(gtals[kk], mask02);

						REP(4) allele[kk] = _mm512_add_epi32(allele[kk], oindex[kk]);

						if constexpr (sizeof(REAL) == 8 || !fast_fp32)
							REP(8) f2[kk] = _mm512_setzero_pd();
						else
							REP(4) f2s[kk] = _mm512_setzero_ps();

						if constexpr (!isadmix)
						{
							REAL* p2 = p + Z[i] * KT;

							if constexpr (sizeof(REAL) == 8)
								REP(8) f2[kk] = _mm512_mask_i32gather_pd(maskoned, typed2[kk], allele2[kk], p2, sizeof(double));
							else if constexpr (fast_fp32)
								REP(4) f2s[kk] = _mm512_mask_i32gather_ps(maskones, typed[kk], allele[kk], p2, sizeof(float));
							else
								REP(8) f2[kk] = _mm512_cvtps_pd(_mm256_mmask_i32gather_ps(maskones256, typed2[kk], allele2[kk], p2, sizeof(float)));
						}
						else
						{
							REAL* p2 = p;
							for (int k = 0; k < K; ++k, p2 += KT)
							{
								if constexpr (sizeof(REAL) == 8)
								{
									__m512d pp[8], qq = _mm512_set1_pd(q[k]);

									REP(8) pp[kk] = _mm512_mask_i32gather_pd(maskoned, typed2[kk], allele2[kk], p2, sizeof(double));
									//REP(8) f2[kk] = _mm512_mask_add_pd(maskoned, typed2[kk], f2[kk], _mm512_mul_pd(pp[kk], qq));
									REP(8) f2[kk] = _mm512_add_pd(f2[kk], _mm512_mul_pd(pp[kk], qq));
								}
								else if constexpr (fast_fp32)
								{
									__m512 pp[4], qq = _mm512_set1_ps(q[k]);

									REP(4) pp[kk] = _mm512_mask_i32gather_ps(maskones, typed[kk], allele[kk], p2, sizeof(float));
									//REP(4) f2s[kk] = _mm512_mask_add_ps(maskones, typed[kk], f2s[kk], _mm512_mul_ps(pp[kk], qq));
									REP(4) f2s[kk] = _mm512_add_ps(f2s[kk], _mm512_mul_ps(pp[kk], qq));
								}
								else
								{
									__m512d pp[8], qq = _mm512_set1_pd(q[k]);

									REP(8) pp[kk] = _mm512_cvtps_pd(_mm256_mmask_i32gather_ps(maskones256, typed2[kk], allele2[kk], p2, sizeof(float)));
									//REP(8) f2[kk] = _mm512_mask_add_pd(maskoned, typed2[kk], f2[kk], _mm512_mul_pd(pp[kk], qq));
									REP(8) f2[kk] = _mm512_add_pd(f2[kk], _mm512_mul_pd(pp[kk], qq));
								}
							}
						}

						if constexpr (sizeof(REAL) == 4 && fast_fp32)
						{
							__m256* f2t = (__m256*)f2s;
							REP(8) f2[kk] = _mm512_cvtps_pd(f2t[kk]);
						}

						REP(8) f2[kk] = _mm512_mask_mov_pd(maskoned, typed2[kk], f2[kk]);

						for (int KK = sizeof(f2) / sizeof(f2[0]) / 2; KK >= 1; KK >>= 1)
							REP(KK) f2[kk] = _mm512_mul_pd(f2[kk], f2[kk + KK]);

						ChargeLog(slog, prod, __mm512_reduce_mul_pd(f2[0]));
					}
				}
			}
		}

		CloseLog(slog, prod);
		bufNK1[tid2] = prod;
	}
}

#else

/* Update a priori ancetral proportion for non-admix model */
template<typename REAL>
TARGET512 void BAYESIAN<REAL>::UpdateQNoAdmix512(int tid) { }

/* Update individual or allele origin when ancetral proportion is binary */
template<typename REAL>
TARGET512 void BAYESIAN<REAL>::UpdateZNoAdmix512(int tid) { }

/* Update a priori ancetral proportion for non-admix model */
template<typename REAL>
template<bool fast_fp32>
TARGET512 void BAYESIAN<REAL>::UpdateQMetro512(int tid) { }

/* Update individual or allele origin when ancetral proportion */
template<typename REAL>
template<bool fast_fp32>
TARGET512 void BAYESIAN<REAL>::UpdateZAdmix512(int tid) { }

/* Record updated MCMC parameters */
template<typename REAL>
template<bool isadmix, bool fast_fp32>
TARGET512 void BAYESIAN<REAL>::Record512(int tid) { }

#endif
