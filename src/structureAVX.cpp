/* AVX Bayesian clustering functions */

#include "vcfpop.h"

template struct BAYESIAN<double>;
template struct BAYESIAN<float >;

template TARGETAVX void BAYESIAN<double>::UpdateQNoAdmixAVX(int tid);
template TARGETAVX void BAYESIAN<float >::UpdateQNoAdmixAVX(int tid);

template TARGETAVX void BAYESIAN<double>::UpdateZNoAdmixAVX(int tid);
template TARGETAVX void BAYESIAN<float >::UpdateZNoAdmixAVX(int tid);

template TARGETAVX void BAYESIAN<double>::UpdateQMetroAVX<true >(int tid);
template TARGETAVX void BAYESIAN<double>::UpdateQMetroAVX<false>(int tid);
template TARGETAVX void BAYESIAN<float >::UpdateQMetroAVX<true >(int tid);
template TARGETAVX void BAYESIAN<float >::UpdateQMetroAVX<false>(int tid);

template TARGETAVX void BAYESIAN<double>::UpdateZAdmixAVX<true >(int tid);
template TARGETAVX void BAYESIAN<double>::UpdateZAdmixAVX<false>(int tid);
template TARGETAVX void BAYESIAN<float >::UpdateZAdmixAVX<true >(int tid);
template TARGETAVX void BAYESIAN<float >::UpdateZAdmixAVX<false>(int tid);

template TARGETAVX void BAYESIAN<double>::RecordAVX<true , true >(int tid);
template TARGETAVX void BAYESIAN<double>::RecordAVX<true , false>(int tid);
template TARGETAVX void BAYESIAN<double>::RecordAVX<false, true >(int tid);
template TARGETAVX void BAYESIAN<double>::RecordAVX<false, false>(int tid);
template TARGETAVX void BAYESIAN<float >::RecordAVX<true , true >(int tid);
template TARGETAVX void BAYESIAN<float >::RecordAVX<true , false>(int tid);
template TARGETAVX void BAYESIAN<float >::RecordAVX<false, true >(int tid);
template TARGETAVX void BAYESIAN<float >::RecordAVX<false, false>(int tid);

#ifndef __aarch64__

#ifndef _GENO_READERAVX
template<typename REAL>
struct GENO_READERAVX
{
	__m256i vindex[8];						//Offset of 16 loci
	__m256i data[8];						//Readed bits
	__m256i msize;
	uint* pos;								//Current read pointer
	byte size;								//Number of bits a genotype id used
	byte nbits;								//Number of bits remaining in data

	TARGETAVX GENO_READERAVX()
	{

	}

	TARGETAVX GENO_READERAVX(int indid, int64 l, int64 num, BUCKET* bucket = NULL)
	{
		//set pos and size
		SetZero(this, 1);
		num = Min(num, 64);

		//set bucket from default bucket or assigned bucket
		if (bucket == NULL) bucket = &geno_bucket;

		LIST<OFFSET>& offset = bucket->offset;
		uint64 offset0 = offset[l].offset;
		pos = (uint*)(bucket->base_addr + offset0);

		int* vindex2 = (int*)vindex;
		OFFSET* offset1 = &offset.bucket[l];
		for (int i = 0; i < num; ++i)
			vindex2[i] = offset1[i].offset - offset0;

		size = bucket->offset[l].size;
		msize = _mm256_set1_epi32((1u << size) - 1u);

		if (indid)
		{
			//skip indid * size bits
			int nreadbits = indid * size;

			//pos move nreadbits / 32
			pos += nreadbits >> 5;

			//remain bits to read
			nreadbits &= 31;

			//read 32 bits to data and read remain bits from data
			REP(8) data[kk] = _mm256_srli_epi32(_mm256_i32gather_epi32((int*)pos, vindex[kk], 1), nreadbits);

			pos++;

			//set nbits
			nbits = 32 - nreadbits;
		}
		else
		{
			//read 32 bits
			REP(8) data[kk] = _mm256_i32gather_epi32((int*)pos, vindex[kk], 1);

			pos++;

			//set nbits
			nbits = 32;
		}
	}

	__forceinline
	TARGETAVX void Read(__m256i* gid)
	{
		// if data is empty
		if (nbits < size) [[unlikely]]
		{
			// remain number of bits to read
			int rbits = size - nbits;
			__m256i tmask = _mm256_set1_epi32((1u << rbits) - 1u);

			// move nbits data to gid
			REP(8) gid[kk] = data[kk];

			// read 32 bits to data
			REP(8) data[kk] = _mm256_i32gather_epi32((int*)pos, vindex[kk], 1);

			// read rbits from data and concate to higher bits in gid
			REP(8) gid[kk] = _mm256_or_si256(gid[kk], _mm256_slli_epi32(_mm256_and_si256(data[kk], tmask), nbits));

			//shift right
			REP(8) data[kk] = _mm256_srli_epi32(data[kk], rbits);

			pos++;

			nbits = 32 - rbits;
		}
		else [[likely]]
		{
			//read size bits
			REP(8) gid[kk] = _mm256_and_si256(data[kk], msize);

			//shift right
			REP(8) data[kk] = _mm256_srli_epi32(data[kk], size);

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
	TARGETSSE void Read(__m256i* gid2)
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

__forceinline TARGETAVX double _mm256_reduce_add_pd(__m256d v1)
{
	__m128d v2 = _mm_add_pd(_mm256_extractf128_pd(v1, 0), _mm256_extractf128_pd(v1, 1));
	return simd_f64(v2, 0) + simd_f64(v2, 1);
}

__forceinline TARGETAVX double _mm256_reduce_mul_pd(__m256d v1)
{
	__m128d v2 = _mm_mul_pd(_mm256_extractf128_pd(v1, 0), _mm256_extractf128_pd(v1, 1));
	return simd_f64(v2, 0) * simd_f64(v2, 1);
}

__forceinline TARGETAVX float _mm256_reduce_add_ps(__m256 v1)
{
	__m128 v2 = _mm_mul_ps(_mm256_extractf128_ps(v1, 0), _mm256_extractf128_ps(v1, 1));
	__m128 v3 = _mm_mul_ps(v2, _mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(v2), 8)));
	return simd_f64(v3, 0) + simd_f64(v3, 1);
}

__forceinline TARGETAVX float _mm256_reduce_mul_ps(__m256 v1)
{
	__m128 v2 = _mm_add_ps(_mm256_extractf128_ps(v1, 0), _mm256_extractf128_ps(v1, 1));
	__m128 v3 = _mm_add_ps(v2, _mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(v2), 8)));
	return simd_f64(v3, 0) * simd_f64(v3, 1);
}

__forceinline __m128i _mm256_cvtepi64_epi32x(__m256i a)
{
	return _mm_castps_si128(_mm_shuffle_ps(
		_mm_castpd_ps(_mm256_castpd256_pd128(_mm256_castsi256_pd(a))), 
		_mm_castpd_ps(_mm256_extractf128_pd(_mm256_castsi256_pd(a), 1)), 
		_MM_SHUFFLE(2, 0, 2, 0)));
}

/* Update a priori ancetral proportion for non-admix model */
template<typename REAL>
TARGETAVX void BAYESIAN<REAL>::UpdateQNoAdmixAVX(int tid)
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

		UpdateQNoAdmixAVX(0);

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

	static __m256d maskoned = _mm256_set1_pd(1.0);
	static __m128  maskones128 = _mm_set1_ps(1.0);
	static __m256  maskones = _mm256_set1_ps(1.0);
	static __m256i mask00 = _mm256_set1_epi32(0);
	static __m256i mask24 = _mm256_set1_epi32(0xFFFFFF);
	static __m256i mask16 = _mm256_set1_epi32(0xFFFF);
	static __m256i mask01 = _mm256_set1_epi32(1);
	static __m256i mask02 = _mm256_set1_epi32(2);
	static __m256i maskshift = _mm256_set_epi32(24, 25, 26, 27, 28, 29, 30, 31);
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

				__m256i gtab[8];
				for (int64 k = 0, lt = l; k < 8; ++k, lt += 8)
					gtab[k] = _mm256_set_epi32(
						GetLocTabDiff(lt + 7), GetLocTabDiff(lt + 6), GetLocTabDiff(lt + 5), GetLocTabDiff(lt + 4),
						GetLocTabDiff(lt + 3), GetLocTabDiff(lt + 2), GetLocTabDiff(lt + 1), GetLocTabDiff(lt + 0));

				__m256i lmask[8]; __m128i* lmask128 = (__m128i*)lmask;  __m256i lmask64[16];
				uint64 lmask0 = lend - l >= 64 ? 0xFFFFFFFFFFFFFFFF : (1ull << (lend - l)) - 1ull;
				byte* lmask1 = (byte*)&lmask0;

				REP(8) lmask[kk] = _mm256_sllv_epi32(_mm256_set1_epi32(lmask1[kk]), maskshift);

				REP(16) lmask64[kk] = _mm256_cvtepi32_epi64(lmask128[kk]);

#ifdef __GNUC__
				GENO_READERSSE<REAL> rt(0, l, lend - l);
#else
				GENO_READERAVX<REAL> rt(0, l, lend - l);
#endif

				__m256i oindex[8]; __m128i* oindex2 = (__m128i*)oindex;
				__m256i oindex0 = _mm256_set1_epi64x(allele_freq_offset[l]);

				REP(16) oindex2[kk] = _mm256_cvtepi64_epi32x(_mm256_sub_epi64(_mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(oindex0), _mm256_castsi256_pd(_mm256_maskload_epi64((int64*)&allele_freq_offset[l + (kk << 2)], lmask64[kk])), _mm256_castsi256_pd(lmask64[kk]))), oindex0));

				__m256i gtaddr[8], gtlow[8], gtploidy[8], gtals[8], allele[8], typed[8], typed64[16];
				__m256d freq[16];
				__m128i* allele2 = (__m128i*)allele;

				int64* buf1 = (int64*)bufNK1 + N * K * tid2; double* buf2 = bufNK2 + N * K * tid2;

				for (int i = 0; i < N; ++i, buf1 += K, buf2 += K)
				{
					rt.Read(gtaddr);

					REP(8) gtaddr[kk] = _mm256_add_epi32(gtab[kk], _mm256_slli_epi32(gtaddr[kk], 2));

					REP(8) gtlow[kk] = _mm256_i32gather_epi32((int*)gtab_base, gtaddr[kk], 1);

					REP(8) gtploidy[kk] = _mm256_srli_epi32(gtlow[kk], 24);

					REP(8) gtlow[kk] = _mm256_and_si256(gtlow[kk], mask24);

					REP(8) gtploidy[kk] = _mm256_mask_i32gather_epi32(mask00, PT_PLOIDYxNALLELES, gtploidy[kk], lmask[kk], sizeof(int));

					REP(8) gtals[kk] = _mm256_add_epi32(gtaddr[kk], gtlow[kk]);

					int maxv = maxploidy;
					if (maxploidy != minploidy)
					{
						__m256i maxv1 = _mm256_max_epi32(_mm256_max_epi32(_mm256_max_epi32(gtploidy[0], gtploidy[1]), _mm256_max_epi32(gtploidy[2], gtploidy[3])), _mm256_max_epi32(_mm256_max_epi32(gtploidy[4], gtploidy[5]), _mm256_max_epi32(gtploidy[6], gtploidy[7])));
						int* maxv2 = (int*)&maxv1;
						maxv = Max(Max(Max(maxv2[0], maxv2[1]), Max(maxv2[2], maxv2[3])), Max(Max(maxv2[4], maxv2[5]), Max(maxv2[6], maxv2[7])));
					}

					__m256i ai = _mm256_set1_epi32(0);
					for (int a = 0; a < maxv; ++a, ai = _mm256_add_epi32(ai, mask01))
					{
						REP(8)
						{
							typed[kk] = _mm256_cmpgt_epi32(gtploidy[kk], ai);

							if constexpr (sizeof(REAL) == 8)
							{
								typed64[0 + (kk << 1)] = _mm256_cvtepi32_epi64(_mm256_extracti128_si256(typed[kk], 0));
								typed64[1 + (kk << 1)] = _mm256_cvtepi32_epi64(_mm256_extracti128_si256(typed[kk], 1));
							}

							allele[kk] = _mm256_mask_i32gather_epi32(mask00, (int*)gtab_base, gtals[kk], typed[kk], 1);

							gtals[kk] = _mm256_add_epi32(gtals[kk], mask02);

							allele[kk] = _mm256_and_si256(allele[kk], mask16);

							allele[kk] = _mm256_add_epi32(allele[kk], oindex[kk]);
						}

						REAL* p2 = p;
						for (int k = 0; k < K; ++k, p2 += KT)
						{
							if constexpr (sizeof(REAL) == 8)
								REP(16) freq[kk] = _mm256_mask_i32gather_pd(maskoned, p2, allele2[kk], _mm256_castsi256_pd(typed64[kk]), sizeof(double));
							else
								REP(8)
							{
								__m256 v1 = _mm256_mask_i32gather_ps(maskones, p2, allele[kk], _mm256_castsi256_ps(typed[kk]), sizeof(float));
								freq[0 + (kk << 1)] = _mm256_cvtps_pd(_mm256_castps256_ps128(v1));
								freq[1 + (kk << 1)] = _mm256_cvtps_pd(_mm256_extractf128_ps(v1, 1));
							}

							REP(8) freq[kk] = _mm256_mul_pd(freq[kk], freq[kk + 8]);
							REP(4) freq[kk] = _mm256_mul_pd(freq[kk], freq[kk + 4]);
							REP(2) freq[kk] = _mm256_mul_pd(freq[kk], freq[kk + 2]);
							REP(1) freq[kk] = _mm256_mul_pd(freq[kk], freq[kk + 1]);

							__m128d* freq2 = (__m128d*)freq;
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
TARGETAVX void BAYESIAN<REAL>::UpdateZNoAdmixAVX(int tid)
{
	if (tid == -1)
	{
		SetZero(Mi, N * K); 
		SetZero(Ni, K * KT);

		for (int i = 0; i < N; ++i)
			Mi[i * K + Z[i]] = ainds[i]->vt;

		//////////////////////////////////////////////////////////

		SetZero((int64*)l_atomic, 32);

		UpdateZNoAdmixAVX(0);

		return;
	}

	static __m256i mask00 = _mm256_set1_epi32(0);
	static __m256i mask24 = _mm256_set1_epi32(0xFFFFFF);
	static __m256i mask16 = _mm256_set1_epi32(0xFFFF);
	static __m256i mask01 = _mm256_set1_epi32(1);
	static __m256i mask02 = _mm256_set1_epi32(2);
	static __m256i maskshift = _mm256_set_epi32(24, 25, 26, 27, 28, 29, 30, 31);
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

				__m256i gtab[8];
				for (int64 k = 0, lt = l; k < 8; ++k, lt += 8)
					gtab[k] = _mm256_set_epi32(
						GetLocTabDiff(lt + 7), GetLocTabDiff(lt + 6), GetLocTabDiff(lt + 5), GetLocTabDiff(lt + 4),
						GetLocTabDiff(lt + 3), GetLocTabDiff(lt + 2), GetLocTabDiff(lt + 1), GetLocTabDiff(lt + 0));

				__m256i lmask[8]; __m128i* lmask128 = (__m128i*)lmask;  __m256i lmask64[16];
				uint64 lmask0 = lend - l >= 64 ? 0xFFFFFFFFFFFFFFFF : (1ull << (lend - l)) - 1ull;
				byte* lmask1 = (byte*)&lmask0;

				REP(8) lmask[kk] = _mm256_sllv_epi32(_mm256_set1_epi32(lmask1[kk]), maskshift);

				REP(16) lmask64[kk] = _mm256_cvtepi32_epi64(lmask128[kk]);

#ifdef __GNUC__
				GENO_READERSSE<REAL> rt(0, l, lend - l);
#else
				GENO_READERAVX<REAL> rt(0, l, lend - l);
#endif

				__m256i oindex[8]; __m128i* oindex2 = (__m128i*)oindex;
				__m256i oindex0 = _mm256_set1_epi64x(allele_freq_offset[l]);

				REP(16) oindex2[kk] = _mm256_cvtepi64_epi32x(_mm256_sub_epi64(_mm256_castpd_si256(_mm256_blendv_pd(
					_mm256_castsi256_pd(oindex0),
					_mm256_castsi256_pd(_mm256_maskload_epi64((int64*)&allele_freq_offset[l + (kk << 2)], lmask64[kk])),
					_mm256_castsi256_pd(lmask64[kk]))), oindex0));

				__m256i gtaddr[8], gtlow[8], gtploidy[8], gtals[8], allele[8];
				__m256i typed[8];
				int* typed2 = (int*)typed, * allele2 = (int*)allele;

				for (int i = 0; i < N; i++)
				{
					rt.Read(gtaddr);

					REP(8) gtaddr[kk] = _mm256_add_epi32(gtab[kk], _mm256_slli_epi32(gtaddr[kk], 2));

					REP(8) gtlow[kk] = _mm256_mask_i32gather_epi32(mask00, (int*)gtab_base, gtaddr[kk], lmask[kk], 1);

					REP(8) gtploidy[kk] = _mm256_srli_epi32(gtlow[kk], 24);

					REP(8) gtlow[kk] = _mm256_and_si256(gtlow[kk], mask24);

					REP(8) gtploidy[kk] = _mm256_i32gather_epi32(PT_PLOIDYxNALLELES, gtploidy[kk], sizeof(int));

					REP(8) gtals[kk] = _mm256_add_epi32(gtaddr[kk], gtlow[kk]);

					int maxv = maxploidy;
					if (maxploidy != minploidy)
					{
						__m256i maxv1 = _mm256_max_epi32(_mm256_max_epi32(_mm256_max_epi32(gtploidy[0], gtploidy[1]), _mm256_max_epi32(gtploidy[2], gtploidy[3])), _mm256_max_epi32(_mm256_max_epi32(gtploidy[4], gtploidy[5]), _mm256_max_epi32(gtploidy[6], gtploidy[7])));
						int* maxv2 = (int*)&maxv1;
						maxv = Max(Max(Max(maxv2[0], maxv2[1]), Max(maxv2[2], maxv2[3])), Max(Max(maxv2[4], maxv2[5]), Max(maxv2[6], maxv2[7])));
					}

					int* ni = Ni + Z[i] * KT + o;
					__m256i ai = _mm256_set1_epi32(0);
					for (int a = 0; a < maxv; ++a, ai = _mm256_add_epi32(ai, mask01))
					{
						REP(8)
						{
							typed[kk] = _mm256_cmpgt_epi32(gtploidy[kk], ai);

							allele[kk] = _mm256_mask_i32gather_epi32(mask00, (int*)gtab_base, gtals[kk], typed[kk], 1);

							gtals[kk] = _mm256_add_epi32(gtals[kk], mask02);

							allele[kk] = _mm256_and_si256(allele[kk], mask16);

							allele[kk] = _mm256_add_epi32(allele[kk], oindex[kk]);
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
TARGETAVX void BAYESIAN<REAL>::UpdateQMetroAVX(int tid)
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

		g_fastsingle_val == 1 ? UpdateQMetroAVX<true >(0) : UpdateQMetroAVX<false>(0);

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

	static __m256d maskoned = _mm256_set1_pd(1.0);
	static __m256  maskones = _mm256_set1_ps(1.0);
	static __m256i mask00 = _mm256_set1_epi32(0);
	static __m256i maskff = _mm256_set1_epi32(0xFFFFFFFF);
	static __m256i mask24 = _mm256_set1_epi32(0xFFFFFF);
	static __m256i mask16 = _mm256_set1_epi32(0xFFFF);
	static __m256i mask01 = _mm256_set1_epi32(1);
	static __m256i mask02 = _mm256_set1_epi32(2);
	static __m256i maskshift = _mm256_set_epi32(24, 25, 26, 27, 28, 29, 30, 31);
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

				__m256i gtab[8];
				for (int64 k = 0, lt = l; k < 8; ++k, lt += 8)
					gtab[k] = _mm256_set_epi32(
						GetLocTabDiff(lt + 7), GetLocTabDiff(lt + 6), GetLocTabDiff(lt + 5), GetLocTabDiff(lt + 4),
						GetLocTabDiff(lt + 3), GetLocTabDiff(lt + 2), GetLocTabDiff(lt + 1), GetLocTabDiff(lt + 0));

				__m256i lmask[8]; __m128i* lmask128 = (__m128i*)lmask;  __m256i lmask64[16];
				uint64 lmask0 = lend - l >= 64 ? 0xFFFFFFFFFFFFFFFF : (1ull << (lend - l)) - 1ull;
				byte* lmask1 = (byte*)&lmask0;

				REP(8) lmask[kk] = _mm256_sllv_epi32(_mm256_set1_epi32(lmask1[kk]), maskshift);
				REP(8) lmask[kk] = _mm256_castps_si256(_mm256_blendv_ps(_mm256_castsi256_ps(mask00), _mm256_castsi256_ps(maskff), _mm256_castsi256_ps(lmask[kk])));
				REP(16) lmask64[kk] = _mm256_cvtepi32_epi64(lmask128[kk]);

#ifdef __GNUC__
				GENO_READERSSE<REAL> rt(0, l, lend - l);
#else
				GENO_READERAVX<REAL> rt(0, l, lend - l);
#endif

				__m256i oindex[8]; __m128i* oindex2 = (__m128i*)oindex;
				__m256i oindex0 = _mm256_set1_epi64x(allele_freq_offset[l]);

				REP(16) oindex2[kk] = _mm256_cvtepi64_epi32x(_mm256_sub_epi64(_mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(oindex0), _mm256_castsi256_pd(_mm256_maskload_epi64((int64*)&allele_freq_offset[l + (kk << 2)], lmask64[kk])), _mm256_castsi256_pd(lmask64[kk]))), oindex0));

				__m256i gtaddr[8], gtlow[8], gtploidy[8], gtals[8], allele[8], typed[8], typed64[16];
				__m128i* allele2 = (__m128i*)allele;
				__m256d f1[16], f2[16];
				__m256 f1s[8], f2s[8];

				REAL* bufi = (REAL*)bufNK1, * q = Q;
				//avoid thread-conflict
				double* buf1 = bufN1 + N * tid2, * buf2 = bufN2 + N * tid2;

				for (int i = 0; i < N; ++i, q += K, bufi += K, buf1++, buf2++)
				{
					rt.Read(gtaddr);

					REP(8) gtaddr[kk] = _mm256_add_epi32(gtab[kk], _mm256_slli_epi32(gtaddr[kk], 2));

					REP(8) gtlow[kk] = _mm256_i32gather_epi32((int*)gtab_base, gtaddr[kk], 1);

					REP(8) gtploidy[kk] = _mm256_srli_epi32(gtlow[kk], 24);

					REP(8) gtlow[kk] = _mm256_and_si256(gtlow[kk], mask24);

					REP(8) gtploidy[kk] = _mm256_mask_i32gather_epi32(mask00, PT_PLOIDYxNALLELES, gtploidy[kk], lmask[kk], sizeof(int));

					REP(8) gtals[kk] = _mm256_add_epi32(gtaddr[kk], gtlow[kk]);

					int maxv = maxploidy;
					if (maxploidy != minploidy)
					{
						__m256i maxv1 = _mm256_max_epi32(_mm256_max_epi32(_mm256_max_epi32(gtploidy[0], gtploidy[1]), _mm256_max_epi32(gtploidy[2], gtploidy[3])), _mm256_max_epi32(_mm256_max_epi32(gtploidy[4], gtploidy[5]), _mm256_max_epi32(gtploidy[6], gtploidy[7])));
						int* maxv2 = (int*)&maxv1;
						maxv = Max(Max(Max(maxv2[0], maxv2[1]), Max(maxv2[2], maxv2[3])), Max(Max(maxv2[4], maxv2[5]), Max(maxv2[6], maxv2[7])));
					}

					__m256i ai = _mm256_set1_epi32(0);
					for (int a = 0; a < maxv; ++a, ai = _mm256_add_epi32(ai, mask01))
					{
						REP(8)
						{
							typed[kk] = _mm256_cmpgt_epi32(gtploidy[kk], ai);

							typed64[0 + (kk << 1)] = _mm256_cvtepi32_epi64(_mm256_extracti128_si256(typed[kk], 0));
							typed64[1 + (kk << 1)] = _mm256_cvtepi32_epi64(_mm256_extracti128_si256(typed[kk], 1));

							allele[kk] = _mm256_mask_i32gather_epi32(mask00, (int*)gtab_base, gtals[kk], typed[kk], 1);

							gtals[kk] = _mm256_add_epi32(gtals[kk], mask02);

							allele[kk] = _mm256_and_si256(allele[kk], mask16);

							allele[kk] = _mm256_add_epi32(allele[kk], oindex[kk]);
						}

						if constexpr (sizeof(REAL) == 8 || !fast_fp32)
						{
							REP(16) f1[kk] = _mm256_setzero_pd();
							REP(16) f2[kk] = _mm256_setzero_pd();
						}
						else
						{
							REP(8) f1s[kk] = _mm256_setzero_ps();
							REP(8) f2s[kk] = _mm256_setzero_ps();
						}

						REAL* p2 = p;
						for (int k = 0; k < K; ++k, p2 += KT)
						{
							if constexpr (sizeof(REAL) == 8)
							{
								__m256d pp[16], v1, v2, qq = _mm256_set1_pd(q[k]), ii = _mm256_set1_pd(bufi[k]);
								REP(16) pp[kk] = _mm256_i32gather_pd(p2, allele2[kk], sizeof(double));
								REP(16)
								{
									//v1 = _mm256_mul_pd(pp[kk], ii);
									//v2 = _mm256_mul_pd(pp[kk], qq);
									//memory fence to diable fmadd
									//if (L == -1) [[unlikely]] { v1 = v2 = qq; }
									f1[kk] = _mm256_add_pd(f1[kk], _mm256_mul_pd(pp[kk], ii));
									f2[kk] = _mm256_add_pd(f2[kk], _mm256_mul_pd(pp[kk], qq));
								}
							}
							else if constexpr (fast_fp32)
							{
								__m256 pp[8], v1, v2, qq = _mm256_set1_ps(q[k]), ii = _mm256_set1_ps(bufi[k]);

								REP(8) pp[kk] = _mm256_i32gather_ps(p2, allele[kk], sizeof(float));
								REP(8)
								{
									//v1 = _mm256_mul_ps(pp[kk], ii);
									//v2 = _mm256_mul_ps(pp[kk], qq);
									//memory fence to diable fmadd
									//if (L == -1) [[unlikely]] { v1 = v2 = qq; }
									f1s[kk] = _mm256_add_ps(f1s[kk], _mm256_mul_ps(pp[kk], ii));
									f2s[kk] = _mm256_add_ps(f2s[kk], _mm256_mul_ps(pp[kk], qq));
								}
							}
							else
							{
								__m256d pp[16], v1, v2, qq = _mm256_set1_pd(q[k]), ii = _mm256_set1_pd(bufi[k]);

								REP(16) pp[kk] = _mm256_cvtps_pd(_mm_i32gather_ps(p2, allele2[kk], sizeof(float)));
								REP(16)
								{
									//v1 = _mm256_mul_pd(pp[kk], ii);
									//v2 = _mm256_mul_pd(pp[kk], qq);
									//memory fence to diable fmadd
									//if (L == -1) [[unlikely]] { v1 = v2 = qq; }
									f1[kk] = _mm256_add_pd(f1[kk], _mm256_mul_pd(pp[kk], ii));
									f2[kk] = _mm256_add_pd(f2[kk], _mm256_mul_pd(pp[kk], qq));
								}
							}
						}

						if constexpr (sizeof(REAL) == 4 && fast_fp32)
						{
							__m128* f1t = (__m128*)f1s;

							REP(8) f1s[kk] = _mm256_div_ps(f1s[kk], f2s[kk]);
							REP(8) f1s[kk] = _mm256_blendv_ps(maskones, f1s[kk], _mm256_castsi256_ps(typed[kk]));
							REP(16) f1[kk] = _mm256_cvtps_pd(f1t[kk]);
						}
						else
						{
							REP(16) f1[kk] = _mm256_div_pd(f1[kk], f2[kk]);
							REP(16) f1[kk] = _mm256_blendv_pd(maskoned, f1[kk], _mm256_castsi256_pd(typed64[kk]));
						}

						for (int KK = sizeof(f1) / sizeof(f1[0]) / 2; KK >= 1; KK >>= 1)
							REP(KK) f1[kk] = _mm256_mul_pd(f1[kk], f1[kk + KK]);

						ChargeLog(*(int64*)buf1, *buf2, _mm256_reduce_mul_pd(f1[0]));
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
TARGETAVX void BAYESIAN<REAL>::UpdateZAdmixAVX(int tid)
{
	if (tid == -1)
	{
		SetZero(Mi, N * K * structure_nsubthread);
		SetZero(Ni, K * KT);

		//////////////////////////////////////////////////////////

		SetZero((int64*)l_atomic, 32);

		g_fastsingle_val == 1 ? UpdateZAdmixAVX<true >(0) : UpdateZAdmixAVX<false>(0);

		//avoid thread-conflict
		for (int i = 1; i < structure_nsubthread; ++i)
			Add(Mi, Mi + N * K * i, N * K);

		return;
	}

	static __m256d maskoned = _mm256_set1_pd(1.0);
	static __m256  maskones = _mm256_set1_ps(1.0);
	static __m256i mask00 = _mm256_set1_epi32(0);
	static __m256i mask24 = _mm256_set1_epi32(0xFFFFFF);
	static __m256i mask16 = _mm256_set1_epi32(0xFFFF);
	static __m256i mask01 = _mm256_set1_epi32(1);
	static __m256i mask02 = _mm256_set1_epi32(2);
	static __m256i maskshift = _mm256_set_epi32(24, 25, 26, 27, 28, 29, 30, 31);
	static __m256d minfreqd = _mm256_set1_pd(MIN_FREQ);
	static __m256 minfreqs = _mm256_set1_ps(MIN_FREQ);
	static int PT_PLOIDYxNALLELES[150] = 									//Pattern index to ploidy level
	{ 0, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	atomic<int> thread_counter = 0;
#pragma omp parallel num_threads(structure_nsubthread)
	{
		int tid2 = thread_counter.fetch_add(1);

		//64 * K  * sizeof(double)
		__m256d* bufkd = (__m256d*)Align64((byte*)bufNK2 + (Max(N, 64) * K * sizeof(double) + 63) * tid2);
		__m256* bufks = (__m256*)Align64((byte*)bufNK2 + (Max(N, 64) * K * sizeof(double) + 63) * tid2);

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

				__m256i gtab[8];
				for (int64 k = 0, lt = l; k < 8; ++k, lt += 8)
					gtab[k] = _mm256_set_epi32(
						GetLocTabDiff(lt + 7), GetLocTabDiff(lt + 6), GetLocTabDiff(lt + 5), GetLocTabDiff(lt + 4),
						GetLocTabDiff(lt + 3), GetLocTabDiff(lt + 2), GetLocTabDiff(lt + 1), GetLocTabDiff(lt + 0));

				__m256i lmask[8]; __m128i* lmask128 = (__m128i*)lmask;  __m256i lmask64[16];
				uint64 lmask0 = lend - l >= 64 ? 0xFFFFFFFFFFFFFFFF : (1ull << (lend - l)) - 1ull;
				byte* lmask1 = (byte*)&lmask0;

				REP(8) lmask[kk] = _mm256_sllv_epi32(_mm256_set1_epi32(lmask1[kk]), maskshift);
				REP(16) lmask64[kk] = _mm256_cvtepi32_epi64(lmask128[kk]);

#ifdef __GNUC__
				GENO_READERSSE<REAL> rt(0, l, lend - l);
#else
				GENO_READERAVX<REAL> rt(0, l, lend - l);
#endif

				__m256i oindex[8], ot[16];
				__m128i* oindex2 = (__m128i*)oindex, * ot2 = (__m128i*)ot;
				__m256i oindex0 = _mm256_set1_epi64x(allele_freq_offset[l]);

				REP(16) ot[kk] = _mm256_maskload_epi64((int64*)&allele_freq_offset[l + (kk << 2)], lmask64[kk]);
				REP(16) ot[kk] = _mm256_castpd_si256(_mm256_blendv_pd(
					_mm256_castsi256_pd(oindex0),
					_mm256_castsi256_pd(ot[kk]),
					_mm256_castsi256_pd(lmask64[kk])));
				REP(16) ot[kk] = _mm256_sub_epi64(ot[kk], oindex0);
				REP(16) oindex2[kk] = _mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(ot2[0 + (kk << 1)]), _mm_castsi128_ps(ot2[1 + (kk << 1)]), _MM_SHUFFLE(2, 0, 2, 0)));

				__m256i gtaddr[8], gtlow[8], gtploidy[8], gtals[8], allele[8], allele64[16], k2[16], typed[8], typed64[16];
				__m128i* allele3 = (__m128i*)allele, * typed3 = (__m128i*)typed;
				int* typed2 = (int*)typed, * allele2 = (int*)allele;
				uint64* k22 = (uint64*)&k2;

				RNGAVX<double> rngd; RNGAVX<float > rngs;

				if constexpr (sizeof(REAL) == 8 || !fast_fp32)
					new (&rngd) RNGAVX<double>(seed + m * L + l, RNG_SALT_UPDATEZ);
				else
					new (&rngs) RNGAVX<float >(seed + m * L + l, RNG_SALT_UPDATEZ);

				REAL* q = Q;
				//avoid thread-conflict
				int64* mi = Mi + N * K * tid2;

				for (int i = 0; i < N; ++i, q += K, mi += K)
				{
					rt.Read(gtaddr);

					REP(8) gtaddr[kk] = _mm256_add_epi32(gtab[kk], _mm256_slli_epi32(gtaddr[kk], 2));

					REP(8) gtlow[kk] = _mm256_mask_i32gather_epi32(mask00, (int*)gtab_base, gtaddr[kk], lmask[kk], 1);

					REP(8) gtploidy[kk] = _mm256_srli_epi32(gtlow[kk], 24);

					REP(8) gtlow[kk] = _mm256_and_si256(gtlow[kk], mask24);

					REP(8) gtploidy[kk] = _mm256_i32gather_epi32(PT_PLOIDYxNALLELES, gtploidy[kk], sizeof(int));

					REP(8) gtals[kk] = _mm256_add_epi32(gtaddr[kk], gtlow[kk]);

					int maxv = maxploidy;

					__m256i ai = _mm256_set1_epi32(0);
					for (int a = 0; a < maxv; ++a, ai = _mm256_add_epi32(ai, mask01))
					{
						REP(8)
						{
							typed[kk] = _mm256_cmpgt_epi32(gtploidy[kk], ai);

							allele[kk] = _mm256_mask_i32gather_epi32(mask00, (int*)gtab_base, gtals[kk], typed[kk], 1);

							gtals[kk] = _mm256_add_epi32(gtals[kk], mask02);

							allele[kk] = _mm256_and_si256(allele[kk], mask16);

							allele[kk] = _mm256_add_epi32(allele[kk], oindex[kk]);
						}

						REP(16)
						{
							allele64[kk] = _mm256_cvtepi32_epi64(allele3[kk]);
							if constexpr (sizeof(REAL) == 8 || !fast_fp32)
								typed64[kk] = _mm256_cvtepi32_epi64(typed3[kk]);
						}

						REAL* p2 = p;
						for (int k = 0; k < K; ++k, p2 += KT)
						{
							__m256d qd; __m256 qs;

							if constexpr (sizeof(REAL) == 8)
							{
								qd = _mm256_set1_pd(q[k]);
								REP(16) bufkd[kk + (k << 4)] = _mm256_i32gather_pd(p2, allele3[kk], sizeof(double));
								REP(16) bufkd[kk + (k << 4)] = _mm256_mul_pd(bufkd[kk + (k << 4)], qd);
								REP(16) bufkd[kk + (k << 4)] = _mm256_add_pd(bufkd[kk + (k << 4)], minfreqd);
								REP(16) bufkd[kk + (k << 4)] = _mm256_and_pd(bufkd[kk + (k << 4)], _mm256_castsi256_pd(typed64[kk]));
							}
							else if constexpr (fast_fp32)
							{
								qs = _mm256_set1_ps(q[k]);
								REP(8) bufks[kk + (k << 3)] = _mm256_i32gather_ps(p2, allele[kk], sizeof(float));
								REP(8) bufks[kk + (k << 3)] = _mm256_mul_ps(bufks[kk + (k << 3)], qs);
								REP(8) bufks[kk + (k << 3)] = _mm256_add_ps(bufks[kk + (k << 3)], minfreqs);
								REP(8) bufks[kk + (k << 3)] = _mm256_and_ps(bufks[kk + (k << 3)], _mm256_castsi256_ps(typed[kk]));
							}
							else
							{
								qd = _mm256_set1_pd(q[k]);
								REP(16) bufkd[kk + (k << 4)] = _mm256_cvtps_pd(_mm_i32gather_ps(p2, allele3[kk], sizeof(float)));
								REP(16) bufkd[kk + (k << 4)] = _mm256_mul_pd(bufkd[kk + (k << 4)], qd);
								REP(16) bufkd[kk + (k << 4)] = _mm256_add_pd(bufkd[kk + (k << 4)], minfreqd);
								REP(16) bufkd[kk + (k << 4)] = _mm256_and_pd(bufkd[kk + (k << 4)], _mm256_castsi256_pd(typed64[kk]));
							}
						}

						//draw cluster for each allele copy
						if constexpr (sizeof(REAL) == 8 || !fast_fp32)
							rngd.Poly<64>(bufkd, K, k2);
						else
							rngs.Poly<64>(bufks, K, k2);

						//Update Mi
						REP(64) mi[k22[kk]] -= typed2[kk];

						REP(16) k2[kk] = _mm256_set_epi64x(k22[3 + (kk << 2)] * KT, k22[2 + (kk << 2)] * KT, k22[1 + (kk << 2)] * KT, k22[0 + (kk << 2)] * KT);

						REP(16) k2[kk] = _mm256_add_epi64(k2[kk], allele64[kk]);

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
TARGETAVX void BAYESIAN<REAL>::RecordAVX(int tid)
{
	if (tid == -1)
	{
		Add(MiSum, Mi, N * K);

		//////////////////////////////////////////////////////////

		SetZero((int64*)l_atomic, 32);

		switch ((int)binaryq * 10 + g_fastsingle_val)
		{
		case 01: BAYESIAN<REAL>::RecordAVX<true , true >(0); break;
		case 02: BAYESIAN<REAL>::RecordAVX<true , false>(0); break;
		case 11: BAYESIAN<REAL>::RecordAVX<false, true >(0); break;
		case 12: BAYESIAN<REAL>::RecordAVX<false, false>(0); break;
		}

		bufNK1[0] = Sum(bufNK1, structure_nsubthread);
		return;
	}

	static __m256d maskoned = _mm256_set1_pd(1.0);
	static __m256  maskones = _mm256_set1_ps(1.0);
	static __m128  maskones128 = _mm_set1_ps(1.0);
	static __m256i mask00 = _mm256_set1_epi32(0);
	static __m256i mask24 = _mm256_set1_epi32(0xFFFFFF);
	static __m256i mask16 = _mm256_set1_epi32(0xFFFF);
	static __m256i mask01 = _mm256_set1_epi32(1);
	static __m256i mask02 = _mm256_set1_epi32(2);
	static __m256i maskshift = _mm256_set_epi32(24, 25, 26, 27, 28, 29, 30, 31);
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

				__m256i gtab[8];
				for (int64 k = 0, lt = l; k < 8; ++k, lt += 8)
					gtab[k] = _mm256_set_epi32(
						GetLocTabDiff(lt + 7), GetLocTabDiff(lt + 6), GetLocTabDiff(lt + 5), GetLocTabDiff(lt + 4),
						GetLocTabDiff(lt + 3), GetLocTabDiff(lt + 2), GetLocTabDiff(lt + 1), GetLocTabDiff(lt + 0));

				__m256i lmask[8]; __m128i* lmask128 = (__m128i*)lmask;  __m256i lmask64[16];
				uint64 lmask0 = lend - l >= 64 ? 0xFFFFFFFFFFFFFFFF : (1ull << (lend - l)) - 1ull;
				byte* lmask1 = (byte*)&lmask0;

				REP(8) lmask[kk] = _mm256_sllv_epi32(_mm256_set1_epi32(lmask1[kk]), maskshift);
				REP(16) lmask64[kk] = _mm256_cvtepi32_epi64(lmask128[kk]);

#ifdef __GNUC__
				GENO_READERSSE<REAL> rt(0, l, lend - l);
#else
				GENO_READERAVX<REAL> rt(0, l, lend - l);
#endif

				__m256i oindex[8], ot[16];
				__m128i* oindex2 = (__m128i*)oindex, * ot2 = (__m128i*)ot;
				__m256i oindex0 = _mm256_set1_epi64x(allele_freq_offset[l]);

				REP(16) ot[kk] = _mm256_maskload_epi64((int64*)&allele_freq_offset[l + (kk << 2)], lmask64[kk]);
				REP(16) ot[kk] = _mm256_castpd_si256(_mm256_blendv_pd(
					_mm256_castsi256_pd(oindex0),
					_mm256_castsi256_pd(ot[kk]),
					_mm256_castsi256_pd(lmask64[kk])));
				REP(16) ot[kk] = _mm256_sub_epi64(ot[kk], oindex0);
				REP(16) oindex2[kk] = _mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(ot2[0 + (kk << 1)]), _mm_castsi128_ps(ot2[1 + (kk << 1)]), _MM_SHUFFLE(2, 0, 2, 0)));

				__m256i gtaddr[8], gtlow[8], gtploidy[8], gtals[8], allele[8], allele64[16], k2[16], typed[8], typed64[16];
				__m256d f2[16];
				__m256 f2s[8];
				__m128i* allele3 = (__m128i*)allele, * typed3 = (__m128i*)typed;

				REAL* q = Q;

				for (int i = 0; i < N; ++i, q += K)
				{
					rt.Read(gtaddr);

					REP(8) gtaddr[kk] = _mm256_add_epi32(gtab[kk], _mm256_slli_epi32(gtaddr[kk], 2));

					REP(8) gtlow[kk] = _mm256_mask_i32gather_epi32(mask00, (int*)gtab_base, gtaddr[kk], lmask[kk], 1);

					REP(8) gtploidy[kk] = _mm256_srli_epi32(gtlow[kk], 24);

					REP(8) gtlow[kk] = _mm256_and_si256(gtlow[kk], mask24);

					REP(8) gtploidy[kk] = _mm256_i32gather_epi32(PT_PLOIDYxNALLELES, gtploidy[kk], sizeof(int));

					REP(8) gtals[kk] = _mm256_add_epi32(gtaddr[kk], gtlow[kk]);

					int maxv = maxploidy;
					if (maxploidy != minploidy)
					{
						__m256i maxv1 = _mm256_max_epi32(_mm256_max_epi32(_mm256_max_epi32(gtploidy[0], gtploidy[1]), _mm256_max_epi32(gtploidy[2], gtploidy[3])), _mm256_max_epi32(_mm256_max_epi32(gtploidy[4], gtploidy[5]), _mm256_max_epi32(gtploidy[6], gtploidy[7])));
						int* maxv2 = (int*)&maxv1;
						maxv = Max(Max(Max(maxv2[0], maxv2[1]), Max(maxv2[2], maxv2[3])), Max(Max(maxv2[4], maxv2[5]), Max(maxv2[6], maxv2[7])));
					}

					__m256i ai = _mm256_set1_epi32(0);
					for (int a = 0; a < maxv; ++a, ai = _mm256_add_epi32(ai, mask01))
					{
						REP(8)
						{
							typed[kk] = _mm256_cmpgt_epi32(gtploidy[kk], ai);

							allele[kk] = _mm256_mask_i32gather_epi32(mask00, (int*)gtab_base, gtals[kk], typed[kk], 1);

							gtals[kk] = _mm256_add_epi32(gtals[kk], mask02);

							allele[kk] = _mm256_and_si256(allele[kk], mask16);

							allele[kk] = _mm256_add_epi32(allele[kk], oindex[kk]);
						}

						REP(16)
						{
							allele64[kk] = _mm256_cvtepi32_epi64(allele3[kk]);
							typed64[kk] = _mm256_cvtepi32_epi64(typed3[kk]);
						}

						if constexpr (sizeof(REAL) == 8 || !fast_fp32)
							REP(16) f2[kk] = _mm256_setzero_pd();
						else
							REP(8) f2s[kk] = _mm256_setzero_ps();

						if constexpr (!isadmix)
						{
							REAL* p2 = p + Z[i] * KT;

							if constexpr (sizeof(REAL) == 8)
								REP(16) f2[kk] = _mm256_mask_i32gather_pd(maskoned, p2, allele3[kk], _mm256_castsi256_pd(typed64[kk]), sizeof(double));
							else if constexpr (fast_fp32)
								REP(8) f2s[kk] = _mm256_mask_i32gather_ps(maskones, p2, allele[kk], _mm256_castsi256_ps(typed[kk]), sizeof(float));
							else
								REP(16) f2[kk] = _mm256_cvtps_pd(_mm_mask_i32gather_ps(maskones128, p2, allele3[kk], _mm_castsi128_ps(typed3[kk]), sizeof(float)));
						}
						else
						{
							REAL* p2 = p;
							for (int k = 0; k < K; ++k, p2 += KT)
							{
								if constexpr (sizeof(REAL) == 8)
								{
									__m256d pp[16], qq = _mm256_set1_pd(q[k]);

									REP(16) pp[kk] = _mm256_mask_i32gather_pd(maskoned, p2, allele3[kk], _mm256_castsi256_pd(typed64[kk]), sizeof(double));
									REP(16) f2[kk] = _mm256_add_pd(f2[kk], _mm256_mul_pd(pp[kk], qq));
								}
								else if constexpr (fast_fp32)
								{
									__m256 pp[8], qq = _mm256_set1_ps(q[k]);

									REP(8) pp[kk] = _mm256_mask_i32gather_ps(maskones, p2, allele[kk], _mm256_castsi256_ps(typed[kk]), sizeof(float));
									REP(8) f2s[kk] = _mm256_add_ps(f2s[kk], _mm256_mul_ps(pp[kk], qq));
								}
								else
								{
									__m256d pp[16], qq = _mm256_set1_pd(q[k]);

									REP(16) pp[kk] = _mm256_cvtps_pd(_mm_mask_i32gather_ps(maskones128, p2, allele3[kk], _mm_castsi128_ps(typed3[kk]), sizeof(float)));
									REP(16) f2[kk] = _mm256_add_pd(f2[kk], _mm256_mul_pd(pp[kk], qq));
								}
							}
						}

						if constexpr (sizeof(REAL) == 4 && fast_fp32)
						{
							__m128* f2t = (__m128*)f2s;
							REP(16) f2[kk] = _mm256_cvtps_pd(f2t[kk]);
						}

						REP(16) f2[kk] = _mm256_blendv_pd(maskoned, f2[kk], _mm256_castsi256_pd(typed64[kk]));

						for (int KK = sizeof(f2) / sizeof(f2[0]) / 2; KK >= 1; KK >>= 1)
							REP(KK) f2[kk] = _mm256_mul_pd(f2[kk], f2[kk + KK]);

						ChargeLog(slog, prod, _mm256_reduce_mul_pd(f2[0]));
					}
				}
			}
		}

		CloseLog(slog, prod);
		bufNK1[tid2] = prod;
	}
}

#else

template<typename REAL>
TARGETAVX void BAYESIAN<REAL>::UpdateQNoAdmixAVX(int tid) { }

template<typename REAL>
TARGETAVX void BAYESIAN<REAL>::UpdateZNoAdmixAVX(int tid) { }

template<typename REAL>
template<bool fast_fp32>
TARGETAVX void BAYESIAN<REAL>::UpdateQMetroAVX(int tid) { }

template<typename REAL>
template<bool fast_fp32>
TARGETAVX void BAYESIAN<REAL>::UpdateZAdmixAVX(int tid) { }

template<typename REAL>
template<bool isadmix, bool fast_fp32>
TARGETAVX void BAYESIAN<REAL>::RecordAVX(int tid) { }

#endif