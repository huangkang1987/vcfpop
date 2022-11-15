/* SSE Bayesian clustering functions */

#include "vcfpop.h"

template struct BAYESIAN<double>;
template struct BAYESIAN<float >;

template TARGETSSE void BAYESIAN<double>::UpdateQNoAdmixSSE(int tid);
template TARGETSSE void BAYESIAN<float >::UpdateQNoAdmixSSE(int tid);

template TARGETSSE void BAYESIAN<double>::UpdateZNoAdmixSSE(int tid);
template TARGETSSE void BAYESIAN<float >::UpdateZNoAdmixSSE(int tid);

template TARGETSSE void BAYESIAN<double>::UpdateQMetroSSE<true >(int tid);
template TARGETSSE void BAYESIAN<double>::UpdateQMetroSSE<false>(int tid);
template TARGETSSE void BAYESIAN<float >::UpdateQMetroSSE<true >(int tid);
template TARGETSSE void BAYESIAN<float >::UpdateQMetroSSE<false>(int tid);

template TARGETSSE void BAYESIAN<double>::UpdateZAdmixSSE<true >(int tid);
template TARGETSSE void BAYESIAN<double>::UpdateZAdmixSSE<false>(int tid);
template TARGETSSE void BAYESIAN<float >::UpdateZAdmixSSE<true >(int tid);
template TARGETSSE void BAYESIAN<float >::UpdateZAdmixSSE<false>(int tid);

template TARGETSSE void BAYESIAN<double>::RecordSSE<true , true >(int tid);
template TARGETSSE void BAYESIAN<double>::RecordSSE<true , false>(int tid);
template TARGETSSE void BAYESIAN<double>::RecordSSE<false, true >(int tid);
template TARGETSSE void BAYESIAN<double>::RecordSSE<false, false>(int tid);
template TARGETSSE void BAYESIAN<float >::RecordSSE<true , true >(int tid);
template TARGETSSE void BAYESIAN<float >::RecordSSE<true , false>(int tid);
template TARGETSSE void BAYESIAN<float >::RecordSSE<false, true >(int tid);
template TARGETSSE void BAYESIAN<float >::RecordSSE<false, false>(int tid);

#ifndef __aarch64__

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
	TARGETSSE void Read(__m128i* gid)
	{
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

__forceinline TARGETSSE double _mm_reduce_add_pd(__m128d v2)
{
	return simd_f64(v2, 0) + simd_f64(v2, 1);
}

__forceinline TARGETSSE double _mm_reduce_mul_pd(__m128d v2)
{
	return simd_f64(v2, 0) * simd_f64(v2, 1);
}

__forceinline TARGETSSE float _mm_reduce_add_ps(__m128 v2)
{
	__m128 v3 = _mm_add_ps(v2, _mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(v2), 8)));
	return simd_f32(v3, 0) + simd_f32(v3, 1);
}

__forceinline TARGETSSE float _mm_reduce_mul_ps(__m128 v2)
{
	__m128 v3 = _mm_mul_ps(v2, _mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(v2), 8)));
	return simd_f32(v3, 0) * simd_f32(v3, 1);
}

__forceinline TARGETSSE double _mm_reduce_add_psd(__m128 v2)
{
	__m128d v2b = _mm_add_pd(_mm_cvtps_pd(v2),
		_mm_cvtps_pd(_mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(v2), 8))));
	return simd_f64(v2b, 0) + simd_f64(v2b, 1);
}

__forceinline TARGETSSE double _mm_reduce_mul_psd(__m128 v2)
{
	__m128d v2b = _mm_mul_pd(_mm_cvtps_pd(v2),
		_mm_cvtps_pd(_mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(v2), 8))));
	return simd_f64(v2b, 0) * simd_f64(v2b, 1);
}

/* Update a priori ancetral proportion for non-admix model */
template<typename REAL>
TARGETSSE void BAYESIAN<REAL>::UpdateQNoAdmixSSE(int tid)
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

		UpdateQNoAdmixSSE(0);

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

	static __m128d maskoned = _mm_set1_pd(1.0);
	static __m128  maskones = _mm_set1_ps(1.0);
	static __m128i mask00 = _mm_set1_epi32(0);
	static __m128i maskff = _mm_set1_epi32(0xFFFFFFFF);
	static __m128i mask24 = _mm_set1_epi32(0xFFFFFF);
	static __m128i mask01 = _mm_set1_epi32(1);
	static __m128i mask02 = _mm_set1_epi32(2);
	static __m128i maskidx[16] = {
		_mm_set_epi32( 3,  2,  1,  0), _mm_set_epi32( 7,  6,  5,  4), _mm_set_epi32(11, 10,  9,  8), _mm_set_epi32(15, 14, 13, 12),
		_mm_set_epi32(19, 18, 17, 16), _mm_set_epi32(23, 22, 21, 20), _mm_set_epi32(27, 26, 25, 24), _mm_set_epi32(31, 30, 29, 28), 
		_mm_set_epi32(35, 34, 33, 32), _mm_set_epi32(39, 38, 37, 36), _mm_set_epi32(43, 42, 41, 40), _mm_set_epi32(47, 46, 45, 44), 
		_mm_set_epi32(51, 50, 49, 48), _mm_set_epi32(55, 54, 53, 52), _mm_set_epi32(59, 58, 57, 56), _mm_set_epi32(63, 62, 61, 60) };
	static int PT_PLOIDYxNALLELES[150] = 									//Pattern index to ploidy level
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

				__m128i gtab[16];
				for (int64 k = 0, lt = l; k < 16; ++k, lt += 4)
					gtab[k] = _mm_set_epi32(GetLocTabDiff(lt + 3), GetLocTabDiff(lt + 2), GetLocTabDiff(lt + 1), GetLocTabDiff(lt + 0));

				__m128i lmask[16], lmaskidx[16];
				uint64 lmask0 = lend - l >= 64 ? 0xFFFFFFFFFFFFFFFF : (1ull << (lend - l)) - 1ull, lmaskt = lmask0;
				REP(16)
				{
					lmask[kk] = _mm_castps_si128(_mm_blendv_ps(
						_mm_castsi128_ps(mask00),
						_mm_castsi128_ps(maskff),
						_mm_castsi128_ps(_mm_set_epi32(lmaskt << 28, lmaskt << 29, lmaskt << 30, lmaskt << 31))));
					lmaskt >>= 4;
					lmaskidx[kk] = _mm_and_si128(lmask[kk], maskidx[kk]);
				}

				GENO_READERSSE<REAL> rt(0, l, lend - l);

				__m128i oindex[16];
				uint64* toffset = allele_freq_offset + l;
				uint64 oindex01 = toffset[0];
				int* lmaskidx2 = (int*)lmaskidx;

				REP(16) oindex[kk] = _mm_set_epi32(toffset[lmaskidx2[3 + (kk << 2)]] - oindex01, toffset[lmaskidx2[2 + (kk << 2)]] - oindex01, toffset[lmaskidx2[1 + (kk << 2)]] - oindex01, toffset[lmaskidx2[0 + (kk << 2)]] - oindex01);

				__m128i gtaddr[16], gtlow[16], gtploidy[16], gtals[16], allele[16], typed[16], typed64[32];
				int* allele2 = (int*)allele, * gtploidy2 = (int*)gtploidy, * gtaddr2 = (int*)gtaddr;
				__m128d freq[32];

				int64* buf1 = (int64*)bufNK1 + N * K * tid; double* buf2 = bufNK2 + N * K * tid;

				for (int i = 0; i < N; ++i, buf1 += K, buf2 += K)
				{
					rt.Read(gtaddr);

					REP(16) gtaddr[kk] = _mm_add_epi32(gtab[kk], _mm_slli_epi32(gtaddr[kk], 2));

					REP(16) gtlow[kk] = _mm_set_epi32(*(uint*)(gtab_base + gtaddr2[3 + (kk << 2)]), *(uint*)(gtab_base + gtaddr2[2 + (kk << 2)]), *(uint*)(gtab_base + gtaddr2[1 + (kk << 2)]), *(uint*)(gtab_base + gtaddr2[0 + (kk << 2)]));

					REP(16) gtploidy[kk] = _mm_and_si128(lmask[kk], _mm_srli_epi32(gtlow[kk], 24));

					REP(16) gtlow[kk] = _mm_and_si128(gtlow[kk], mask24);

					REP(16) gtploidy[kk] = _mm_set_epi32(PT_PLOIDYxNALLELES[gtploidy2[3 + (kk << 2)]], PT_PLOIDYxNALLELES[gtploidy2[2 + (kk << 2)]], PT_PLOIDYxNALLELES[gtploidy2[1 + (kk << 2)]], PT_PLOIDYxNALLELES[gtploidy2[0 + (kk << 2)]]);

					REP(16) gtals[kk] = _mm_add_epi32(gtaddr[kk], gtlow[kk]);

					int maxv = maxploidy;
					if (maxploidy != minploidy)
					{
						__m128i maxv1 =
							_mm_max_epi32(
								_mm_max_epi32(
									_mm_max_epi32(
										_mm_max_epi32(gtploidy[0], gtploidy[1]),
										_mm_max_epi32(gtploidy[2], gtploidy[3])),
									_mm_max_epi32(
										_mm_max_epi32(gtploidy[4], gtploidy[5]),
										_mm_max_epi32(gtploidy[6], gtploidy[7]))),
								_mm_max_epi32(
									_mm_max_epi32(
										_mm_max_epi32(gtploidy[8], gtploidy[9]),
										_mm_max_epi32(gtploidy[10], gtploidy[11])),
									_mm_max_epi32(
										_mm_max_epi32(gtploidy[12], gtploidy[13]),
										_mm_max_epi32(gtploidy[14], gtploidy[15]))));
						int* maxv2 = (int*)&maxv1;
						maxv = Max(Max(maxv2[0], maxv2[1]), Max(maxv2[2], maxv2[3]));
					}

					__m128i ai = _mm_set1_epi32(0);
					for (int a = 0; a < maxv; ++a, ai = _mm_add_epi32(ai, mask01))
					{
						REP(16)
						{
							typed[kk] = _mm_cmpgt_epi32(gtploidy[kk], ai);

							if constexpr (sizeof(REAL) == 8)
							{
								typed64[0 + (kk << 1)] = _mm_cvtepi32_epi64(typed[kk]);
								typed64[1 + (kk << 1)] = _mm_cvtepi32_epi64(_mm_shuffle_epi32(typed[kk], _MM_SHUFFLE(1, 0, 3, 2)));
							}

							__m128i als3 = _mm_and_si128(gtals[kk], typed[kk]);

							allele[kk] = _mm_set_epi32(
								*(ushort*)(gtab_base + simd_i32(als3, 3)),
								*(ushort*)(gtab_base + simd_i32(als3, 2)),
								*(ushort*)(gtab_base + simd_i32(als3, 1)),
								*(ushort*)(gtab_base + simd_i32(als3, 0)));

							gtals[kk] = _mm_add_epi32(gtals[kk], mask02);

							allele[kk] = _mm_add_epi32(allele[kk], oindex[kk]);

							allele[kk] = _mm_and_si128(allele[kk], typed[kk]);
						}

						REAL* p2 = p;
						for (int k = 0; k < K; ++k, p2 += KT)
						{
							if constexpr (sizeof(REAL) == 8)
							{
								REP(32)
								{
									__m128d v1 = _mm_set_pd(p2[allele2[1 + (kk << 1)]], p2[allele2[0 + (kk << 1)]]);
									freq[kk] = _mm_blendv_pd(maskoned, v1, _mm_castsi128_pd(typed64[kk]));
								}
							}
							else
							{
								REP(16)
								{
									__m128 v1 = _mm_set_ps(p2[allele2[3 + (kk << 2)]], p2[allele2[2 + (kk << 2)]], p2[allele2[1 + (kk << 2)]], p2[allele2[0 + (kk << 2)]]);
									v1 = _mm_blendv_ps(maskones, v1, _mm_castsi128_ps(typed[kk]));
									freq[0 + (kk << 1)] = _mm_cvtps_pd(v1);
									freq[1 + (kk << 1)] = _mm_cvtps_pd(_mm_shuffle_ps(v1, v1, _MM_SHUFFLE(1, 0, 3, 2)));
								}
							}

							for (int K = sizeof(freq) / sizeof(freq[0]) / 2; K >= 1; K >>= 1)
								REP(K) freq[kk] = _mm_mul_pd(freq[kk], freq[kk + K]);

							ChargeLog(buf1[k], buf2[k], simp_f64(freq, 0) * simp_f64(freq, 1));
						}
					}
				}
			}
		}

		//avoid thread-conflict
		CloseLog((int64*)bufNK1 + N * K * tid, bufNK2 + N * K * tid, N * K);
	}
}

/* Update individual or allele origin when ancetral proportion is binary */
template<typename REAL>
TARGETSSE void BAYESIAN<REAL>::UpdateZNoAdmixSSE(int tid)
{
	if (tid == -1)
	{
		SetZero(Mi, N * K);
		SetZero(Ni, K * KT);

		for (int i = 0; i < N; ++i)
			Mi[i * K + Z[i]] = ainds[i]->vt;

		//////////////////////////////////////////////////////////

		SetZero((int64*)l_atomic, 32);

		UpdateZNoAdmixSSE(0);

		return;
	}

	static __m128i mask00 = _mm_set1_epi32(0);
	static __m128i maskff = _mm_set1_epi32(0xFFFFFFFF);
	static __m128i mask24 = _mm_set1_epi32(0xFFFFFF);
	static __m128i mask01 = _mm_set1_epi32(1);
	static __m128i mask02 = _mm_set1_epi32(2);
	static __m128i maskidx[16] = {
		_mm_set_epi32( 3,  2,  1,  0), _mm_set_epi32( 7,  6,  5,  4), _mm_set_epi32(11, 10,  9,  8), _mm_set_epi32(15, 14, 13, 12),
		_mm_set_epi32(19, 18, 17, 16), _mm_set_epi32(23, 22, 21, 20), _mm_set_epi32(27, 26, 25, 24), _mm_set_epi32(31, 30, 29, 28), 
		_mm_set_epi32(35, 34, 33, 32), _mm_set_epi32(39, 38, 37, 36), _mm_set_epi32(43, 42, 41, 40), _mm_set_epi32(47, 46, 45, 44), 
		_mm_set_epi32(51, 50, 49, 48), _mm_set_epi32(55, 54, 53, 52), _mm_set_epi32(59, 58, 57, 56), _mm_set_epi32(63, 62, 61, 60) };
	static int PT_PLOIDYxNALLELES[150] = 									//Pattern index to ploidy level
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
				int64 o = allele_freq_offset[l];
				int64 gtab_base = (int64)GetLoc(l).GetGtab();

				__m128i gtab[16];
				for (int64 k = 0, lt = l; k < 16; ++k, lt += 4)
					gtab[k] = _mm_set_epi32(GetLocTabDiff(lt + 3), GetLocTabDiff(lt + 2), GetLocTabDiff(lt + 1), GetLocTabDiff(lt + 0));

				__m128i lmask[16], lmaskidx[16];
				uint64 lmask0 = lend - l >= 64 ? 0xFFFFFFFFFFFFFFFF : (1ull << (lend - l)) - 1ull, lmaskt = lmask0;
				REP(16)
				{
					lmask[kk] = _mm_castps_si128(_mm_blendv_ps(
						_mm_castsi128_ps(mask00),
						_mm_castsi128_ps(maskff),
						_mm_castsi128_ps(_mm_set_epi32(lmaskt << 28, lmaskt << 29, lmaskt << 30, lmaskt << 31))));
					lmaskt >>= 4;
					lmaskidx[kk] = _mm_and_si128(lmask[kk], maskidx[kk]);
				}

				GENO_READERSSE<REAL> rt(0, l, lend - l);

				__m128i oindex[16];
				uint64* toffset = allele_freq_offset + l;
				uint64 oindex01 = toffset[0];
				int* lmaskidx2 = (int*)lmaskidx;

				REP(16) oindex[kk] = _mm_set_epi32(toffset[lmaskidx2[3 + (kk << 2)]] - oindex01, toffset[lmaskidx2[2 + (kk << 2)]] - oindex01, toffset[lmaskidx2[1 + (kk << 2)]] - oindex01, toffset[lmaskidx2[0 + (kk << 2)]] - oindex01);

				__m128i gtaddr[16], gtlow[16], gtploidy[16], gtals[16], allele[16], typed[16];
				int* typed2 = (int*)typed, * allele2 = (int*)allele, * gtploidy2 = (int*)gtploidy, * gtaddr2 = (int*)gtaddr;

				for (int i = 0; i < N; i++)
				{
					rt.Read(gtaddr);

					REP(16) gtaddr[kk] = _mm_add_epi32(gtab[kk], _mm_slli_epi32(gtaddr[kk], 2));

					REP(16) gtlow[kk] = _mm_set_epi32(*(uint*)(gtab_base + gtaddr2[3 + (kk << 2)]), *(uint*)(gtab_base + gtaddr2[2 + (kk << 2)]), *(uint*)(gtab_base + gtaddr2[1 + (kk << 2)]), *(uint*)(gtab_base + gtaddr2[0 + (kk << 2)]));

					REP(16) gtploidy[kk] = _mm_and_si128(lmask[kk], _mm_srli_epi32(gtlow[kk], 24));

					REP(16) gtlow[kk] = _mm_and_si128(gtlow[kk], mask24);

					REP(16) gtploidy[kk] = _mm_set_epi32(PT_PLOIDYxNALLELES[gtploidy2[3 + (kk << 2)]], PT_PLOIDYxNALLELES[gtploidy2[2 + (kk << 2)]], PT_PLOIDYxNALLELES[gtploidy2[1 + (kk << 2)]], PT_PLOIDYxNALLELES[gtploidy2[0 + (kk << 2)]]);

					REP(16) gtals[kk] = _mm_add_epi32(gtaddr[kk], gtlow[kk]);

					int maxv = maxploidy;
					if (maxploidy != minploidy)
					{
						__m128i maxv1 =
							_mm_max_epi32(
								_mm_max_epi32(
									_mm_max_epi32(
										_mm_max_epi32(gtploidy[0], gtploidy[1]),
										_mm_max_epi32(gtploidy[2], gtploidy[3])),
									_mm_max_epi32(
										_mm_max_epi32(gtploidy[4], gtploidy[5]),
										_mm_max_epi32(gtploidy[6], gtploidy[7]))),
								_mm_max_epi32(
									_mm_max_epi32(
										_mm_max_epi32(gtploidy[8], gtploidy[9]),
										_mm_max_epi32(gtploidy[10], gtploidy[11])),
									_mm_max_epi32(
										_mm_max_epi32(gtploidy[12], gtploidy[13]),
										_mm_max_epi32(gtploidy[14], gtploidy[15]))));
						int* maxv2 = (int*)&maxv1;
						maxv = Max(Max(maxv2[0], maxv2[1]), Max(maxv2[2], maxv2[3]));
					}

					int* ni = Ni + Z[i] * KT + o;
					__m128i ai = _mm_set1_epi32(0);
					for (int a = 0; a < maxv; ++a, ai = _mm_add_epi32(ai, mask01))
					{
						REP(16)
						{
							typed[kk] = _mm_cmpgt_epi32(gtploidy[kk], ai);

							__m128i als3 = _mm_and_si128(gtals[kk], typed[kk]);

							allele[kk] = _mm_set_epi32(
								*(ushort*)(gtab_base + simd_i32(als3, 3)),
								*(ushort*)(gtab_base + simd_i32(als3, 2)),
								*(ushort*)(gtab_base + simd_i32(als3, 1)),
								*(ushort*)(gtab_base + simd_i32(als3, 0)));

							allele[kk] = _mm_add_epi32(allele[kk], oindex[kk]);

							gtals[kk] = _mm_add_epi32(gtals[kk], mask02);

							allele[kk] = _mm_and_si128(allele[kk], typed[kk]);
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
TARGETSSE void BAYESIAN<REAL>::UpdateQMetroSSE(int tid)
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

		g_fastsingle_val == 1 ? UpdateQMetroSSE<true >(0) : UpdateQMetroSSE<false>(0);

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
	
	static __m128d maskoned = _mm_set1_pd(1.0);
	static __m128  maskones = _mm_set1_ps(1.0);
	static __m128i mask00 = _mm_set1_epi32(0);
	static __m128i maskff = _mm_set1_epi32(0xFFFFFFFF);
	static __m128i mask24 = _mm_set1_epi32(0xFFFFFF);
	static __m128i mask01 = _mm_set1_epi32(1);
	static __m128i mask02 = _mm_set1_epi32(2);
	static __m128i maskidx[16] = {
		_mm_set_epi32( 3,  2,  1,  0), _mm_set_epi32( 7,  6,  5,  4), _mm_set_epi32(11, 10,  9,  8), _mm_set_epi32(15, 14, 13, 12),
		_mm_set_epi32(19, 18, 17, 16), _mm_set_epi32(23, 22, 21, 20), _mm_set_epi32(27, 26, 25, 24), _mm_set_epi32(31, 30, 29, 28), 
		_mm_set_epi32(35, 34, 33, 32), _mm_set_epi32(39, 38, 37, 36), _mm_set_epi32(43, 42, 41, 40), _mm_set_epi32(47, 46, 45, 44), 
		_mm_set_epi32(51, 50, 49, 48), _mm_set_epi32(55, 54, 53, 52), _mm_set_epi32(59, 58, 57, 56), _mm_set_epi32(63, 62, 61, 60) };
	static int PT_PLOIDYxNALLELES[150] = 									//Pattern index to ploidy level
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

				__m128i gtab[16];
				for (int64 k = 0, lt = l; k < 16; ++k, lt += 4)
					gtab[k] = _mm_set_epi32(GetLocTabDiff(lt + 3), GetLocTabDiff(lt + 2), GetLocTabDiff(lt + 1), GetLocTabDiff(lt + 0));

				__m128i lmask[16], lmaskidx[16];
				uint64 lmask0 = lend - l >= 64 ? 0xFFFFFFFFFFFFFFFF : (1ull << (lend - l)) - 1ull, lmaskt = lmask0;
				REP(16)
				{
					lmask[kk] = _mm_castps_si128(_mm_blendv_ps(
						_mm_castsi128_ps(mask00),
						_mm_castsi128_ps(maskff),
						_mm_castsi128_ps(_mm_set_epi32(lmaskt << 28, lmaskt << 29, lmaskt << 30, lmaskt << 31))));
					lmaskt >>= 4;
					lmaskidx[kk] = _mm_and_si128(lmask[kk], maskidx[kk]);
				}

				GENO_READERSSE<REAL> rt(0, l, lend - l);

				__m128i oindex[16];
				uint64* toffset = allele_freq_offset + l;
				uint64 oindex01 = toffset[0];
				int* lmaskidx2 = (int*)lmaskidx;

				REP(16) oindex[kk] = _mm_set_epi32(toffset[lmaskidx2[3 + (kk << 2)]] - oindex01, toffset[lmaskidx2[2 + (kk << 2)]] - oindex01, toffset[lmaskidx2[1 + (kk << 2)]] - oindex01, toffset[lmaskidx2[0 + (kk << 2)]] - oindex01);

				__m128i gtaddr[16], gtlow[16], gtploidy[16], gtals[16], allele[16], typed[16], typed64[32];
				int* allele2 = (int*)allele, * gtploidy2 = (int*)gtploidy, * gtaddr2 = (int*)gtaddr;
				__m128d f1[32], f2[32];
				__m128 f1s[16], f2s[16];

				REAL* bufi = (REAL*)bufNK1, * q = Q;
				//avoid thread-conflict
				double* buf1 = bufN1 + N * tid, * buf2 = bufN2 + N * tid;

				for (int i = 0; i < N; ++i, q += K, bufi += K, buf1++, buf2++)
				{
					rt.Read(gtaddr);

					REP(16) gtaddr[kk] = _mm_add_epi32(gtab[kk], _mm_slli_epi32(gtaddr[kk], 2));

					REP(16) gtlow[kk] = _mm_set_epi32(*(uint*)(gtab_base + gtaddr2[3 + (kk << 2)]), *(uint*)(gtab_base + gtaddr2[2 + (kk << 2)]), *(uint*)(gtab_base + gtaddr2[1 + (kk << 2)]), *(uint*)(gtab_base + gtaddr2[0 + (kk << 2)]));

					REP(16) gtploidy[kk] = _mm_and_si128(lmask[kk], _mm_srli_epi32(gtlow[kk], 24));

					REP(16) gtlow[kk] = _mm_and_si128(gtlow[kk], mask24);

					REP(16) gtploidy[kk] = _mm_set_epi32(PT_PLOIDYxNALLELES[gtploidy2[3 + (kk << 2)]], PT_PLOIDYxNALLELES[gtploidy2[2 + (kk << 2)]], PT_PLOIDYxNALLELES[gtploidy2[1 + (kk << 2)]], PT_PLOIDYxNALLELES[gtploidy2[0 + (kk << 2)]]);

					REP(16) gtals[kk] = _mm_add_epi32(gtaddr[kk], gtlow[kk]);

					int maxv = maxploidy;
					if (maxploidy != minploidy)
					{
						__m128i maxv1 =
							_mm_max_epi32(
								_mm_max_epi32(
									_mm_max_epi32(
										_mm_max_epi32(gtploidy[0], gtploidy[1]),
										_mm_max_epi32(gtploidy[2], gtploidy[3])),
									_mm_max_epi32(
										_mm_max_epi32(gtploidy[4], gtploidy[5]),
										_mm_max_epi32(gtploidy[6], gtploidy[7]))),
								_mm_max_epi32(
									_mm_max_epi32(
										_mm_max_epi32(gtploidy[8], gtploidy[9]),
										_mm_max_epi32(gtploidy[10], gtploidy[11])),
									_mm_max_epi32(
										_mm_max_epi32(gtploidy[12], gtploidy[13]),
										_mm_max_epi32(gtploidy[14], gtploidy[15]))));
						int* maxv2 = (int*)&maxv1;
						maxv = Max(Max(maxv2[0], maxv2[1]), Max(maxv2[2], maxv2[3]));
					}

					__m128i ai = _mm_set1_epi32(0);
					for (int a = 0; a < maxv; ++a, ai = _mm_add_epi32(ai, mask01))
					{
						REP(16)
						{
							typed[kk] = _mm_cmpgt_epi32(gtploidy[kk], ai);

							if constexpr (sizeof(REAL) == 8 || !fast_fp32)
							{
								typed64[0 + (kk << 1)] = _mm_cvtepi32_epi64(typed[kk]);
								typed64[1 + (kk << 1)] = _mm_cvtepi32_epi64(_mm_shuffle_epi32(typed[kk], _MM_SHUFFLE(1, 0, 3, 2)));
							}

							__m128i als3 = _mm_and_si128(gtals[kk], typed[kk]);

							allele[kk] = _mm_set_epi32(
								*(ushort*)(gtab_base + simd_i32(als3, 3)),
								*(ushort*)(gtab_base + simd_i32(als3, 2)),
								*(ushort*)(gtab_base + simd_i32(als3, 1)),
								*(ushort*)(gtab_base + simd_i32(als3, 0)));

							gtals[kk] = _mm_add_epi32(gtals[kk], mask02);

							allele[kk] = _mm_add_epi32(allele[kk], oindex[kk]);

							allele[kk] = _mm_and_si128(allele[kk], typed[kk]);
						}

						if constexpr (sizeof(REAL) == 8 || !fast_fp32)
						{
							REP(32) f1[kk] = _mm_setzero_pd();
							REP(32) f2[kk] = _mm_setzero_pd();
						}
						else
						{
							REP(16) f1s[kk] = _mm_setzero_ps();
							REP(16) f2s[kk] = _mm_setzero_ps();
						}

						REAL* p2 = p;
						for (int k = 0; k < K; ++k, p2 += KT)
						{
							if constexpr (sizeof(REAL) == 8 || !fast_fp32)
							{
								__m128d pp[32], qq = _mm_set1_pd(q[k]), ii = _mm_set1_pd(bufi[k]);

								REP(32)
								{
									pp[kk] = _mm_set_pd(p2[allele2[1 + (kk << 1)]], p2[allele2[0 + (kk << 1)]]);
									f1[kk] = _mm_add_pd(f1[kk], _mm_mul_pd(pp[kk], ii));
									f2[kk] = _mm_add_pd(f2[kk], _mm_mul_pd(pp[kk], qq));
								}
							}
							else
							{
								__m128 pp[16], v1[16], v2[16], qq = _mm_set1_ps(q[k]), ii = _mm_set1_ps(bufi[k]);

								REP(16)
								{
									pp[kk] = _mm_set_ps(p2[allele2[3 + (kk << 2)]], p2[allele2[2 + (kk << 2)]], p2[allele2[1 + (kk << 2)]], p2[allele2[0 + (kk << 2)]]);
									v1[kk] = _mm_mul_ps(pp[kk], ii);
									v2[kk] = _mm_mul_ps(pp[kk], qq);
									f1s[kk] = _mm_add_ps(f1s[kk], v1[kk]);
									f2s[kk] = _mm_add_ps(f2s[kk], v2[kk]);
								}
							}
						}

						if constexpr (sizeof(REAL) == 4 && fast_fp32)
						{
							REP(16) f1s[kk] = _mm_div_ps(f1s[kk], f2s[kk]);
							REP(16) f1s[kk] = _mm_blendv_ps(maskones, f1s[kk], _mm_castsi128_ps(typed[kk]));
							REP(16)
							{
								f1[0 + (kk << 1)] = _mm_cvtps_pd(f1s[kk]);
								f1[1 + (kk << 1)] = _mm_cvtps_pd(_mm_shuffle_ps(f1s[kk], f1s[kk], _MM_SHUFFLE(1, 0, 3, 2)));
							}
						}
						else
						{
							REP(32)
							{
								//f2[kk] = _mm_blendv_pd(maskoned, f2[kk], _mm_castsi128_pd(typed64[kk]));
								f1[kk] = _mm_div_pd(f1[kk], f2[kk]);
								f1[kk] = _mm_blendv_pd(maskoned, f1[kk], _mm_castsi128_pd(typed64[kk]));
							}
						}

						for (int K = sizeof(f1) / sizeof(f1[0]) / 2; K >= 1; K >>= 1)
							REP(K) f1[kk] = _mm_mul_pd(f1[kk], f1[kk + K]);

						ChargeLog(*(int64*)buf1, *buf2, simp_f64(f1, 0) * simp_f64(f1, 1));
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
TARGETSSE void BAYESIAN<REAL>::UpdateZAdmixSSE(int tid)
{
	if (tid == -1)
	{
		SetZero(Mi, N * K * structure_nsubthread);
		SetZero(Ni, K * KT);

		//////////////////////////////////////////////////////////

		SetZero((int64*)l_atomic, 32);

		g_fastsingle_val == 1 ? UpdateZAdmixSSE<true >(0) : UpdateZAdmixSSE<false>(0);

		//avoid thread-conflict
		for (int i = 1; i < structure_nsubthread; ++i)
			Add(Mi, Mi + N * K * i, N * K);

		return;
	}

	static __m128i mask00 = _mm_set1_epi32(0);
	static __m128i maskff = _mm_set1_epi32(0xFFFFFFFF);
	static __m128i mask24 = _mm_set1_epi32(0xFFFFFF);
	static __m128i mask01 = _mm_set1_epi32(1);
	static __m128i mask02 = _mm_set1_epi32(2);
	static __m128i maskidx[16] = {
		_mm_set_epi32( 3,  2,  1,  0), _mm_set_epi32( 7,  6,  5,  4), _mm_set_epi32(11, 10,  9,  8), _mm_set_epi32(15, 14, 13, 12),
		_mm_set_epi32(19, 18, 17, 16), _mm_set_epi32(23, 22, 21, 20), _mm_set_epi32(27, 26, 25, 24), _mm_set_epi32(31, 30, 29, 28), 
		_mm_set_epi32(35, 34, 33, 32), _mm_set_epi32(39, 38, 37, 36), _mm_set_epi32(43, 42, 41, 40), _mm_set_epi32(47, 46, 45, 44), 
		_mm_set_epi32(51, 50, 49, 48), _mm_set_epi32(55, 54, 53, 52), _mm_set_epi32(59, 58, 57, 56), _mm_set_epi32(63, 62, 61, 60) };
	static __m128d minfreqd = _mm_set1_pd(MIN_FREQ);
	static __m128 minfreqs = _mm_set1_ps(MIN_FREQ);
	static int PT_PLOIDYxNALLELES[150] = 									//Pattern index to ploidy level
	{ 0, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	atomic<int> thread_counter = 0;
#pragma omp parallel num_threads(structure_nsubthread)
	{
		int tid = thread_counter.fetch_add(1);

		//64 * K  * sizeof(double)
		__m128d* bufkd = (__m128d*)Align64((byte*)bufNK2 + (Max(N, 64) * K * sizeof(double) + 63) * tid);
		__m128* bufks = (__m128*)Align64((byte*)bufNK2 + (Max(N, 64) * K * sizeof(double) + 63) * tid);

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

				__m128i gtab[16];
				for (int64 k = 0, lt = l; k < 16; ++k, lt += 4)
					gtab[k] = _mm_set_epi32(GetLocTabDiff(lt + 3), GetLocTabDiff(lt + 2), GetLocTabDiff(lt + 1), GetLocTabDiff(lt + 0));

				__m128i lmask[16], lmaskidx[16];
				uint64 lmask0 = lend - l >= 64 ? 0xFFFFFFFFFFFFFFFF : (1ull << (lend - l)) - 1ull, lmaskt = lmask0;
				REP(16)
				{
					lmask[kk] = _mm_castps_si128(_mm_blendv_ps(
						_mm_castsi128_ps(mask00),
						_mm_castsi128_ps(maskff),
						_mm_castsi128_ps(_mm_set_epi32(lmaskt << 28, lmaskt << 29, lmaskt << 30, lmaskt << 31))));
					lmaskt >>= 4;
					lmaskidx[kk] = _mm_and_si128(lmask[kk], maskidx[kk]);
				}

				GENO_READERSSE<REAL> rt(0, l, lend - l);

				__m128i oindex[16];
				uint64* toffset = allele_freq_offset + l;
				uint64 oindex01 = toffset[0];
				int* lmaskidx2 = (int*)lmaskidx;

				REP(16) oindex[kk] = _mm_set_epi32(toffset[lmaskidx2[3 + (kk << 2)]] - oindex01, toffset[lmaskidx2[2 + (kk << 2)]] - oindex01, toffset[lmaskidx2[1 + (kk << 2)]] - oindex01, toffset[lmaskidx2[0 + (kk << 2)]] - oindex01);

				__m128i gtaddr[16], gtlow[16], gtploidy[16], gtals[16], allele[16], allele64[32], k2[32], typed[16], typed64[32];
				int* typed2 = (int*)typed, * allele2 = (int*)allele, * gtploidy2 = (int*)gtploidy, * gtaddr2 = (int*)gtaddr;
				uint64* k22 = (uint64*)&k2;

				RNGSSE<double> rngd; RNGSSE<float > rngs;

				if constexpr (sizeof(REAL) == 8 || !fast_fp32)
					new (&rngd) RNGSSE<double>(seed + m * L + l, RNG_SALT_UPDATEZ);
				else
					new (&rngs) RNGSSE<float >(seed + m * L + l, RNG_SALT_UPDATEZ);

				REAL* q = Q;
				//avoid thread-conflict
				int64* mi = Mi + N * K * tid;

				for (int i = 0; i < N; ++i, q += K, mi += K)
				{
					rt.Read(gtaddr);

					REP(16) gtaddr[kk] = _mm_add_epi32(gtab[kk], _mm_slli_epi32(gtaddr[kk], 2));

					REP(16) gtlow[kk] = _mm_set_epi32(*(uint*)(gtab_base + gtaddr2[3 + (kk << 2)]), *(uint*)(gtab_base + gtaddr2[2 + (kk << 2)]), *(uint*)(gtab_base + gtaddr2[1 + (kk << 2)]), *(uint*)(gtab_base + gtaddr2[0 + (kk << 2)]));

					REP(16) gtploidy[kk] = _mm_and_si128(lmask[kk], _mm_srli_epi32(gtlow[kk], 24));

					REP(16) gtlow[kk] = _mm_and_si128(gtlow[kk], mask24);

					REP(16) gtploidy[kk] = _mm_set_epi32(PT_PLOIDYxNALLELES[gtploidy2[3 + (kk << 2)]], PT_PLOIDYxNALLELES[gtploidy2[2 + (kk << 2)]], PT_PLOIDYxNALLELES[gtploidy2[1 + (kk << 2)]], PT_PLOIDYxNALLELES[gtploidy2[0 + (kk << 2)]]);

					REP(16) gtals[kk] = _mm_add_epi32(gtaddr[kk], gtlow[kk]);

					int maxv = maxploidy;

					__m128i ai = _mm_set1_epi32(0);
					for (int a = 0; a < maxv; ++a, ai = _mm_add_epi32(ai, mask01))
					{
						REP(16)
						{
							typed[kk] = _mm_cmpgt_epi32(gtploidy[kk], ai);

							__m128i als3 = _mm_and_si128(gtals[kk], typed[kk]);

							allele[kk] = _mm_set_epi32(
								*(ushort*)(gtab_base + simd_i32(als3, 3)),
								*(ushort*)(gtab_base + simd_i32(als3, 2)),
								*(ushort*)(gtab_base + simd_i32(als3, 1)),
								*(ushort*)(gtab_base + simd_i32(als3, 0)));

							gtals[kk] = _mm_add_epi32(gtals[kk], mask02);

							allele[kk] = _mm_add_epi32(allele[kk], oindex[kk]);

							allele[kk] = _mm_and_si128(allele[kk], typed[kk]);
						}

						REP(16)
						{
							allele64[0 + (kk << 1)] = _mm_cvtepi32_epi64(allele[kk]);
							allele64[1 + (kk << 1)] = _mm_cvtepi32_epi64(_mm_shuffle_epi32(allele[kk], _MM_SHUFFLE(1, 0, 3, 2)));

							if constexpr (sizeof(REAL) == 8 || !fast_fp32)
							{
								typed64[0 + (kk << 1)] = _mm_cvtepi32_epi64(typed[kk]);
								typed64[1 + (kk << 1)] = _mm_cvtepi32_epi64(_mm_shuffle_epi32(typed[kk], _MM_SHUFFLE(1, 0, 3, 2)));
							}
						}

						REAL* p2 = p;
						for (int k = 0; k < K; ++k, p2 += KT)
						{
							__m128d qd; __m128 qs;

							if constexpr (sizeof(REAL) == 8)
							{
								qd = _mm_set1_pd(q[k]);

								REP(32) bufkd[kk + (k << 5)] = _mm_set_pd(p2[allele2[1 + (kk << 1)]], p2[allele2[0 + (kk << 1)]]);
								REP(32) bufkd[kk + (k << 5)] = _mm_mul_pd(bufkd[kk + (k << 5)], qd);
								REP(32) bufkd[kk + (k << 5)] = _mm_add_pd(bufkd[kk + (k << 5)], minfreqd);
								REP(32) bufkd[kk + (k << 5)] = _mm_and_pd(bufkd[kk + (k << 5)], _mm_castsi128_pd(typed64[kk]));
							}
							else if constexpr (fast_fp32)
							{
								qs = _mm_set1_ps(q[k]);

								REP(16) bufks[kk + (k << 4)] = _mm_set_ps(p2[allele2[3 + (kk << 2)]], p2[allele2[2 + (kk << 2)]], p2[allele2[1 + (kk << 2)]], p2[allele2[0 + (kk << 2)]]);
								REP(16) bufks[kk + (k << 4)] = _mm_mul_ps(bufks[kk + (k << 4)], qs);
								REP(16) bufks[kk + (k << 4)] = _mm_add_ps(bufks[kk + (k << 4)], minfreqs);
								REP(16) bufks[kk + (k << 4)] = _mm_and_ps(bufks[kk + (k << 4)], _mm_castsi128_ps(typed[kk]));
							}
							else
							{
								qd = _mm_set1_pd(q[k]);
								REP(16)
								{
									__m128 v1 = _mm_set_ps(p2[allele2[3 + (kk << 2)]], p2[allele2[2 + (kk << 2)]], p2[allele2[1 + (kk << 2)]], p2[allele2[0 + (kk << 2)]]);
									bufkd[0 + (kk << 1) + (k << 5)] = _mm_cvtps_pd(v1);
									bufkd[1 + (kk << 1) + (k << 5)] = _mm_cvtps_pd(_mm_shuffle_ps(v1, v1, _MM_SHUFFLE(1, 0, 3, 2)));
								}

								REP(32) bufkd[kk + (k << 5)] = _mm_mul_pd(bufkd[kk + (k << 5)], qd);
								REP(32) bufkd[kk + (k << 5)] = _mm_add_pd(bufkd[kk + (k << 5)], minfreqd);
								REP(32) bufkd[kk + (k << 5)] = _mm_and_pd(bufkd[kk + (k << 5)], _mm_castsi128_pd(typed64[kk]));
							}
						}

						//draw cluster for each allele copy
						if constexpr (sizeof(REAL) == 8 || !fast_fp32)
							rngd.Poly<64>(bufkd, K, k2);
						else
							rngs.Poly<64>(bufks, K, k2);

						//Update Mi
						REP(64) mi[k22[kk]] -= typed2[kk];

						REP(32) k2[kk] = _mm_set_epi64x(k22[1 + (kk << 1)] * KT, k22[0 + (kk << 1)] * KT);

						REP(32) k2[kk] = _mm_add_epi64(k2[kk], allele64[kk]);

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
TARGETSSE void BAYESIAN<REAL>::RecordSSE(int tid)
{
	if (tid == -1)
	{
		Add(MiSum, Mi, N * K);

		//////////////////////////////////////////////////////////

		SetZero((int64*)l_atomic, 32);

		switch ((int)binaryq * 10 + g_fastsingle_val)
		{
		case 01: BAYESIAN<REAL>::RecordSSE<true , true >(0); break;
		case 02: BAYESIAN<REAL>::RecordSSE<true , false>(0); break;
		case 11: BAYESIAN<REAL>::RecordSSE<false, true >(0); break;
		case 12: BAYESIAN<REAL>::RecordSSE<false, false>(0); break;
		}

		bufNK1[0] = Sum(bufNK1, structure_nsubthread);
		return;
	}

	static __m128d maskoned = _mm_set1_pd(1.0);
	static __m128i mask00 = _mm_set1_epi32(0);
	static __m128i maskff = _mm_set1_epi32(0xFFFFFFFF);
	static __m128i mask24 = _mm_set1_epi32(0xFFFFFF);
	static __m128i mask01 = _mm_set1_epi32(1);
	static __m128i mask02 = _mm_set1_epi32(2);
	static __m128i maskidx[16] = {
		_mm_set_epi32(3,  2,  1,  0), _mm_set_epi32(7,  6,  5,  4), _mm_set_epi32(11, 10,  9,  8), _mm_set_epi32(15, 14, 13, 12),
		_mm_set_epi32(19, 18, 17, 16), _mm_set_epi32(23, 22, 21, 20), _mm_set_epi32(27, 26, 25, 24), _mm_set_epi32(31, 30, 29, 28),
		_mm_set_epi32(35, 34, 33, 32), _mm_set_epi32(39, 38, 37, 36), _mm_set_epi32(43, 42, 41, 40), _mm_set_epi32(47, 46, 45, 44),
		_mm_set_epi32(51, 50, 49, 48), _mm_set_epi32(55, 54, 53, 52), _mm_set_epi32(59, 58, 57, 56), _mm_set_epi32(63, 62, 61, 60) };
	static int PT_PLOIDYxNALLELES[150] = 									//Pattern index to ploidy level
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

				__m128i gtab[16];
				for (int64 k = 0, lt = l; k < 16; ++k, lt += 4)
					gtab[k] = _mm_set_epi32(GetLocTabDiff(lt + 3), GetLocTabDiff(lt + 2), GetLocTabDiff(lt + 1), GetLocTabDiff(lt + 0));

				__m128i lmask[16], lmaskidx[16];
				uint64 lmask0 = lend - l >= 64 ? 0xFFFFFFFFFFFFFFFF : (1ull << (lend - l)) - 1ull, lmaskt = lmask0;
				REP(16)
				{
					lmask[kk] = _mm_castps_si128(_mm_blendv_ps(
						_mm_castsi128_ps(mask00),
						_mm_castsi128_ps(maskff),
						_mm_castsi128_ps(_mm_set_epi32(lmaskt << 28, lmaskt << 29, lmaskt << 30, lmaskt << 31))));
					lmaskt >>= 4;
					lmaskidx[kk] = _mm_and_si128(lmask[kk], maskidx[kk]);
				}

				GENO_READERSSE<REAL> rt(0, l, lend - l);

				__m128i oindex[16];
				uint64* toffset = allele_freq_offset + l;
				uint64 oindex01 = toffset[0];
				int* lmaskidx2 = (int*)lmaskidx;

				REP(16) oindex[kk] = _mm_set_epi32(toffset[lmaskidx2[3 + (kk << 2)]] - oindex01, toffset[lmaskidx2[2 + (kk << 2)]] - oindex01, toffset[lmaskidx2[1 + (kk << 2)]] - oindex01, toffset[lmaskidx2[0 + (kk << 2)]] - oindex01);

				__m128i gtaddr[16], gtlow[16], gtploidy[16], gtals[16], allele[16], allele64[32], typed[16], typed64[32];
				int* allele2 = (int*)allele, * gtploidy2 = (int*)gtploidy, * gtaddr2 = (int*)gtaddr;
				__m128d f2[32];
				__m128 f2s[16];

				REAL* q = Q;

				for (int i = 0; i < N; ++i, q += K)
				{
					rt.Read(gtaddr);

					REP(16) gtaddr[kk] = _mm_add_epi32(gtab[kk], _mm_slli_epi32(gtaddr[kk], 2));

					REP(16) gtlow[kk] = _mm_set_epi32(*(uint*)(gtab_base + gtaddr2[3 + (kk << 2)]), *(uint*)(gtab_base + gtaddr2[2 + (kk << 2)]), *(uint*)(gtab_base + gtaddr2[1 + (kk << 2)]), *(uint*)(gtab_base + gtaddr2[0 + (kk << 2)]));

					REP(16) gtploidy[kk] = _mm_and_si128(lmask[kk], _mm_srli_epi32(gtlow[kk], 24));

					REP(16) gtlow[kk] = _mm_and_si128(gtlow[kk], mask24);

					REP(16) gtploidy[kk] = _mm_set_epi32(PT_PLOIDYxNALLELES[gtploidy2[3 + (kk << 2)]], PT_PLOIDYxNALLELES[gtploidy2[2 + (kk << 2)]], PT_PLOIDYxNALLELES[gtploidy2[1 + (kk << 2)]], PT_PLOIDYxNALLELES[gtploidy2[0 + (kk << 2)]]);

					REP(16) gtals[kk] = _mm_add_epi32(gtaddr[kk], gtlow[kk]);

					int maxv = maxploidy;
					if (maxploidy != minploidy)
					{
						__m128i maxv1 =
							_mm_max_epi32(
								_mm_max_epi32(
									_mm_max_epi32(
										_mm_max_epi32(gtploidy[0], gtploidy[1]),
										_mm_max_epi32(gtploidy[2], gtploidy[3])),
									_mm_max_epi32(
										_mm_max_epi32(gtploidy[4], gtploidy[5]),
										_mm_max_epi32(gtploidy[6], gtploidy[7]))),
								_mm_max_epi32(
									_mm_max_epi32(
										_mm_max_epi32(gtploidy[8], gtploidy[9]),
										_mm_max_epi32(gtploidy[10], gtploidy[11])),
									_mm_max_epi32(
										_mm_max_epi32(gtploidy[12], gtploidy[13]),
										_mm_max_epi32(gtploidy[14], gtploidy[15]))));
						int* maxv2 = (int*)&maxv1;
						maxv = Max(Max(maxv2[0], maxv2[1]), Max(maxv2[2], maxv2[3]));
					}

					__m128i ai = _mm_set1_epi32(0);
					for (int a = 0; a < maxv; ++a, ai = _mm_add_epi32(ai, mask01))
					{
						REP(16)
						{
							typed[kk] = _mm_cmpgt_epi32(gtploidy[kk], ai);

							__m128i als3 = _mm_and_si128(gtals[kk], typed[kk]);

							allele[kk] = _mm_set_epi32(
								*(ushort*)(gtab_base + simd_i32(als3, 3)),
								*(ushort*)(gtab_base + simd_i32(als3, 2)),
								*(ushort*)(gtab_base + simd_i32(als3, 1)),
								*(ushort*)(gtab_base + simd_i32(als3, 0)));

							gtals[kk] = _mm_add_epi32(gtals[kk], mask02);

							allele[kk] = _mm_add_epi32(allele[kk], oindex[kk]);

							allele[kk] = _mm_and_si128(allele[kk], typed[kk]);
						}

						REP(16)
						{
							allele64[0 + (kk << 1)] = _mm_cvtepi32_epi64(allele[kk]);
							allele64[1 + (kk << 1)] = _mm_cvtepi32_epi64(_mm_shuffle_epi32(allele[kk], _MM_SHUFFLE(1, 0, 3, 2)));

							typed64[0 + (kk << 1)] = _mm_cvtepi32_epi64(typed[kk]);
							typed64[1 + (kk << 1)] = _mm_cvtepi32_epi64(_mm_shuffle_epi32(typed[kk], _MM_SHUFFLE(1, 0, 3, 2)));
						}

						if constexpr (sizeof(REAL) == 8 || !fast_fp32)
							REP(32) f2[kk] = _mm_setzero_pd();
						else
							REP(16) f2s[kk] = _mm_setzero_ps();

						if constexpr (!isadmix)
						{
							REAL* p2 = p + Z[i] * KT;

							if constexpr (sizeof(REAL) == 8 || !fast_fp32)
								REP(32) f2[kk] = _mm_set_pd(p2[allele2[1 + (kk << 1)]], p2[allele2[0 + (kk << 1)]]);
							else
								REP(16) f2s[kk] = _mm_set_ps(p2[allele2[3 + (kk << 2)]], p2[allele2[2 + (kk << 2)]], p2[allele2[1 + (kk << 2)]], p2[allele2[0 + (kk << 2)]]);
						}
						else
						{
							REAL* p2 = p;
							for (int k = 0; k < K; ++k, p2 += KT)
							{
								//disable fmadd
								if constexpr (sizeof(REAL) == 8 || !fast_fp32)
								{
									__m128d pp[32], qq = _mm_set1_pd(q[k]);

									REP(32) pp[kk] = _mm_set_pd(p2[allele2[1 + (kk << 1)]], p2[allele2[0 + (kk << 1)]]);
									REP(32) f2[kk] = _mm_add_pd(f2[kk], _mm_mul_pd(pp[kk], qq));
								}
								else
								{
									__m128 pp[16], qq = _mm_set1_ps(q[k]);

									REP(16) pp[kk] = _mm_set_ps(p2[allele2[3 + (kk << 2)]], p2[allele2[2 + (kk << 2)]], p2[allele2[1 + (kk << 2)]], p2[allele2[0 + (kk << 2)]]);
									REP(16) f2s[kk] = _mm_add_ps(f2s[kk], _mm_mul_ps(pp[kk], qq));
								}
							}
						}

						if constexpr (sizeof(REAL) == 4 && fast_fp32)
						{
							REP(16)
							{
								f2[0 + (kk << 1)] = _mm_cvtps_pd(f2s[kk]);
								f2[1 + (kk << 1)] = _mm_cvtps_pd(_mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(f2s[kk]), 8)));
							}
						}

						REP(32) f2[kk] = _mm_blendv_pd(maskoned, f2[kk], _mm_castsi128_pd(typed64[kk]));

						for (int K = sizeof(f2) / sizeof(f2[0]) / 2; K >= 1; K >>= 1)
							REP(K) f2[kk] = _mm_mul_pd(f2[kk], f2[kk + K]);

						ChargeLog(slog, prod, simp_f64(f2, 0) * simp_f64(f2, 1));
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
TARGETSSE void BAYESIAN<REAL>::UpdateQNoAdmixSSE(int tid) { }

template<typename REAL>
TARGETSSE void BAYESIAN<REAL>::UpdateZNoAdmixSSE(int tid) { }

template<typename REAL>
template<bool fast_fp32>
TARGETSSE void BAYESIAN<REAL>::UpdateQMetroSSE(int tid) { }

template<typename REAL>
template<bool fast_fp32>
TARGETSSE void BAYESIAN<REAL>::UpdateZAdmixSSE(int tid) { }

template<typename REAL>
template<bool isadmix, bool fast_fp32>
TARGETSSE void BAYESIAN<REAL>::RecordSSE(int tid) { }

#endif