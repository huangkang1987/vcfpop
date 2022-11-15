/* AVX512 Instruction Set Functions */
#include "vcfpop.h"

#ifndef __aarch64__

template TARGET512 void RNG512<double>::Poly<32>(__m512d* arr, int n, __m512i* re);
template TARGET512 void RNG512<double>::Poly<64>(__m512d* arr, int n, __m512i* re);
template TARGET512 void RNG512<float >::Poly<32>(__m512* arr, int n, __m512i* re);
template TARGET512 void RNG512<float >::Poly<64>(__m512* arr, int n, __m512i* re);

#ifndef _RNG512_FP64
/* Initialize rng */
TARGET512 RNG512<double>::RNG512()
{

}

/* Initialize rng */
TARGET512 RNG512<double>::RNG512(uint64 seed, uint64 salt)
{
	__m512i a[8], s, m;

	REP(8) { a[kk] = _mm512_set_epi64(seed + 7, seed + 6, seed + 5, seed + 4, seed + 3, seed + 2, seed + 1, seed + 0); seed += 8; }

	s = _mm512_set1_epi64(salt);
	m = _mm512_set1_epi32(0x5bd1e995);

	REP(8) a[kk] = _mm512_xor_si512(a[kk], _mm512_slli_epi64(_mm512_andnot_si512(a[kk], _mm512_set1_epi64(0xFFFFFFFF)), 32));

	s = _mm512_xor_si512(s, _mm512_slli_epi64(_mm512_andnot_si512(s, _mm512_set1_epi64(0xFFFFFFFF)), 32));

	// uint s = s ^ sizeof(uint);
	s = _mm512_xor_si512(s, _mm512_set1_epi32(sizeof(uint)));

	// a *= m;
	REP(8) a[kk] = _mm512_mullo_epi32(a[kk], m);//OK

	// a ^= a >> 24;
	REP(8) a[kk] = _mm512_xor_si512(a[kk], _mm512_srli_epi32(a[kk], 24));

	// a *= m;
	REP(8) a[kk] = _mm512_mullo_epi32(a[kk], m);//OK

	// s *= m;
	s = _mm512_mullo_epi32(s, m);//OK

	// a ^= s;
	REP(8) a[kk] = _mm512_xor_si512(a[kk], s);

	// a ^= a >> 13;
	REP(8) a[kk] = _mm512_xor_si512(a[kk], _mm512_srli_epi32(a[kk], 13));

	// a *= m;
	REP(8) a[kk] = _mm512_mullo_epi32(a[kk], m);//OK

	// a ^= a >> 15;
	REP(8) a[kk] = _mm512_xor_si512(a[kk], _mm512_srli_epi32(a[kk], 15));

	// original
	REP(8) x[kk] = _mm512_xor_si512(_mm512_set1_epi64(0x159A55E5075BCD15), a[kk]);

	REP(8) a[kk] = _mm512_slli_epi64(a[kk], 6);

	REP(8) y[kk] = _mm512_xor_si512(_mm512_set1_epi64(0x054913331F123BB5), a[kk]);
}

/* Draw a uniform distriubted real number, added minfreq outside */
template<int nbits>
TARGET512 void RNG512<double>::Poly(__m512d* arr, int n, __m512i* re)
{
	__m512d t[8], s[8];
	__m512d one = _mm512_set1_pd(1.0);
	__m512i mask1 = _mm512_set1_epi64(0x000FFFFFFFFFFFFF);
	__m512i mask2 = _mm512_set1_epi64(0x3FF0000000000000);
	__m512i* r = (__m512i*)t;

	REP(8) s[kk] = _mm512_setzero_pd();

	for (int i8 = 0; i8 < n * 8; i8 += 8)
		REP(8) s[kk] = _mm512_add_pd(s[kk], arr[kk + i8]);

	{
		__m512i a[8], b[8];

		REP(8) a[kk] = x[kk];

		REP(8) b[kk] = y[kk];

		REP(8) x[kk] = b[kk];

		REP(8) a[kk] = _mm512_xor_si512(a[kk], _mm512_slli_epi64(a[kk], 23));

		REP(8) a[kk] = _mm512_xor_si512(a[kk], _mm512_srli_epi64(a[kk], 18));

		REP(8) a[kk] = _mm512_xor_si512(a[kk], b[kk]);

		REP(8) a[kk] = _mm512_xor_si512(a[kk], _mm512_srli_epi64(b[kk], 5));

		REP(8) y[kk] = a[kk];

		REP(8) r[kk] = _mm512_add_epi64(a[kk], b[kk]);
	}

	REP(8) r[kk] = _mm512_and_si512(r[kk], mask1);

	REP(8) r[kk] = _mm512_or_si512(r[kk], mask2);

	REP(8) t[kk] = _mm512_sub_pd(t[kk], one);

	REP(8) t[kk] = _mm512_mul_pd(t[kk], s[kk]);

	__m512i midx[8], nidx[8], ninc[8];
	__mmask8 b[8], f[8] = { 0 };

	REP(8) midx[kk] = _mm512_set1_epi64(n - 1);
	REP(8) nidx[kk] = _mm512_setzero_si512();
	REP(8) ninc[kk] = _mm512_set1_epi64(1);

	for (int i8 = 0; i8 < n * 8; i8 += 8)
	{
		REP(8) b[kk] = _mm512_cmp_pd_mask(t[kk], arr[kk + i8], _CMP_LT_OS);

		REP(8) t[kk] = _mm512_sub_pd(t[kk], arr[kk + i8]);

		REP(8) b[kk] = (~f[kk]) & b[kk];

		REP(8) f[kk] = f[kk] | b[kk];

		REP(8) midx[kk] = _mm512_mask_blend_epi64(b[kk], midx[kk], nidx[kk]);//ok

		REP(8) nidx[kk] = _mm512_add_epi64(nidx[kk], ninc[kk]);
	}

	if constexpr (nbits == 32)
	{
		__m256i* re2 = (__m256i*)re;
		REP(8) re2[kk] = _mm512_cvtepi64_epi32(midx[kk]);
	}
	else
	{
		REP(8) re[kk] = midx[kk];
	}
}

#endif

#ifndef _RNG512_FP32
/* Initialize rng */
TARGET512 RNG512<float>::RNG512()
{

}

/* Initialize rng */
TARGET512 RNG512<float>::RNG512(uint64 seed, uint64 salt)
{
	__m512i a[4], s, m;
	REP(4) { a[kk] = _mm512_set_epi32(Mix(seed + 15), Mix(seed + 14), Mix(seed + 13), Mix(seed + 12), Mix(seed + 11), Mix(seed + 10), Mix(seed + 9), Mix(seed + 8), Mix(seed + 7), Mix(seed + 6), Mix(seed + 5), Mix(seed + 4), Mix(seed + 3), Mix(seed + 2), Mix(seed + 1), Mix(seed + 0)); seed += 16; }

	s = _mm512_set1_epi32(Mix(salt));
	m = _mm512_set1_epi32(0x5bd1e995);

	// uint s = s ^ sizeof(uint);
	s = _mm512_xor_si512(s, _mm512_set1_epi32(sizeof(uint)));

	// a *= m;
	REP(4) a[kk] = _mm512_mullo_epi32(a[kk], m);//OK

	// a ^= a >> 24;
	REP(4) a[kk] = _mm512_xor_si512(a[kk], _mm512_srli_epi32(a[kk], 24));

	// a *= m;
	REP(4) a[kk] = _mm512_mullo_epi32(a[kk], m);//OK

	// s *= m;
	s = _mm512_mullo_epi32(s, m);//OK

	// a ^= s;
	REP(4) a[kk] = _mm512_xor_si512(a[kk], s);

	// a ^= a >> 13;
	REP(4) a[kk] = _mm512_xor_si512(a[kk], _mm512_srli_epi32(a[kk], 13));

	// a *= m;
	REP(4) a[kk] = _mm512_mullo_epi32(a[kk], m);//OK

	// a ^= a >> 15;
	REP(4) a[kk] = _mm512_xor_si512(a[kk], _mm512_srli_epi32(a[kk], 15));

	// original
	REP(4) x[kk] = _mm512_xor_si512(_mm512_set1_epi32(0x075BCD15), a[kk]);

	REP(4) a[kk] = _mm512_slli_epi32(a[kk], 3);

	REP(4) y[kk] = _mm512_xor_si512(_mm512_set1_epi32(0x159A55E5), a[kk]);

	REP(4) a[kk] = _mm512_slli_epi32(a[kk], 3);

	REP(4) z[kk] = _mm512_xor_si512(_mm512_set1_epi32(0x1F123BB5), a[kk]);
}

/* Draw a uniform distriubted real number, added minfreq outside */
template<int nbits>
TARGET512 void RNG512<float>::Poly(__m512* arr, int n, __m512i* re)
{
	__m512 t[4], s[4]; __m512i u[16];
	__m512 one = _mm512_set1_ps(1.0f);
	__m512i mask1 = _mm512_set1_epi32(0x007FFFFF);
	__m512i mask2 = _mm512_set1_epi32(0x3F800000);
	__m512i* r = (__m512i*)t;

	REP(4) s[kk] = _mm512_setzero_ps();

	for (int i4 = 0; i4 < n * 4; i4 += 4)
		REP(4) s[kk] = _mm512_add_ps(s[kk], arr[kk + i4]);

	{
		//xorshift
		REP(4) u[kk] = _mm512_slli_epi32(x[kk], 16);
		REP(4) x[kk] = _mm512_xor_si512(x[kk], u[kk]);

		REP(4) u[kk] = _mm512_srli_epi32(x[kk], 5);
		REP(4) x[kk] = _mm512_xor_si512(x[kk], u[kk]);

		REP(4) u[kk] = _mm512_slli_epi32(x[kk], 1);
		REP(4) x[kk] = _mm512_xor_si512(x[kk], u[kk]);

		REP(4) u[kk] = x[kk];

		REP(4) x[kk] = y[kk];

		REP(4) y[kk] = z[kk];

		REP(4) z[kk] = _mm512_xor_si512(u[kk], x[kk]);

		REP(4) z[kk] = _mm512_xor_si512(z[kk], y[kk]);
	}

	REP(4) r[kk] = _mm512_and_si512(z[kk], mask1);

	REP(4) r[kk] = _mm512_or_si512(r[kk], mask2);

	REP(4) t[kk] = _mm512_sub_ps(t[kk], one);

	REP(4) t[kk] = _mm512_mul_ps(t[kk], s[kk]);

	__m512i midx[4] = { _mm512_set1_epi32(n - 1), _mm512_set1_epi32(n - 1), _mm512_set1_epi32(n - 1), _mm512_set1_epi32(n - 1) };
	__m512i nidx[4] = { _mm512_setzero_si512(), _mm512_setzero_si512(), _mm512_setzero_si512(), _mm512_setzero_si512() };
	__m512i ninc[4] = { _mm512_set1_epi32(1), _mm512_set1_epi32(1), _mm512_set1_epi32(1), _mm512_set1_epi32(1) };
	__mmask16 b[4], f[4] = { 0 };

	for (int i4 = 0; i4 < n * 4; i4 += 4)
	{
		REP(4) b[kk] = _mm512_cmp_ps_mask(t[kk], arr[kk + i4], _CMP_LT_OS);

		REP(4) t[kk] = _mm512_sub_ps(t[kk], arr[kk + i4]);

		REP(4) b[kk] = (~f[kk]) & b[kk];

		REP(4) f[kk] = f[kk] | b[kk];

		//GCC has a bug here, cannot correct move nidx to midx according to the mask b2
		REP(4) midx[kk] = _mm512_mask_mov_epi32(midx[kk], b[kk], nidx[kk]);

		REP(4) nidx[kk] = _mm512_add_epi32(nidx[kk], ninc[kk]);
	}

	if constexpr (nbits == 32)
	{
		REP(4) re[kk] = midx[kk];
	}
	else
	{
		__m256i* midx2 = (__m256i*)midx;
		REP(8) re[kk] = _mm512_cvtepi32_epi64(midx2[kk]);
	}
}

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
	__m256 v1 = _mm256_add_ps(_mm512_extractf32x8_ps(v0, 0), _mm512_extractf32x8_ps(v0, 1));
	__m128 v2 = _mm_add_ps(_mm256_extractf128_ps(v1, 0), _mm256_extractf128_ps(v1, 1));
	__m128 v3 = _mm_add_ps(v2, _mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(v2), 8)));
	return simd_f32(v3, 0) + simd_f32(v3, 1);
}

__forceinline TARGET512 float __mm512_reduce_mul_ps(__m512 v0)
{
	__m256 v1 = _mm256_mul_ps(_mm512_extractf32x8_ps(v0, 0), _mm512_extractf32x8_ps(v0, 1));
	__m128 v2 = _mm_mul_ps(_mm256_extractf128_ps(v1, 0), _mm256_extractf128_ps(v1, 1));
	__m128 v3 = _mm_mul_ps(v2, _mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(v2), 8)));
	return simd_f32(v3, 0) * simd_f32(v3, 1);
}

__forceinline TARGET512 double __mm512_reduce_add_psd(__m512 v0)
{
	__m512d v0b = _mm512_add_pd(
		_mm512_cvtps_pd(_mm512_extractf32x8_ps(v0, 0)),
		_mm512_cvtps_pd(_mm512_extractf32x8_ps(v0, 1)));
	__m256d v1 = _mm256_add_pd(_mm512_extractf64x4_pd(v0b, 0), _mm512_extractf64x4_pd(v0b, 1));
	__m128d v2 = _mm_add_pd(_mm256_extractf128_pd(v1, 0), _mm256_extractf128_pd(v1, 1));
	return simd_f64(v2, 0) + simd_f64(v2, 1);
}

__forceinline TARGET512 double __mm512_reduce_mul_psd(__m512 v0)
{
	__m512d v0b = _mm512_mul_pd(
		_mm512_cvtps_pd(_mm512_extractf32x8_ps(v0, 0)),
		_mm512_cvtps_pd(_mm512_extractf32x8_ps(v0, 1)));
	__m256d v1 = _mm256_mul_pd(_mm512_extractf64x4_pd(v0b, 0), _mm512_extractf64x4_pd(v0b, 1));
	__m128d v2 = _mm_mul_pd(_mm256_extractf128_pd(v1, 0), _mm256_extractf128_pd(v1, 1));
	return simd_f64(v2, 0) * simd_f64(v2, 1);
}

TARGET512 int64 GetMinIdx512(double* A, int64 n, double& val)
{
	constexpr int N = 2;
	int64 i = 0;
	val = DBL_MAX;
	uint64 idx = (uint64)-1;

	if (n >= N * sizeof(__m512d) / sizeof(double))
	{
		__m512d min1[N], a[N];
		__m512i midx[N], nidx[N], msep = _mm512_set1_epi64(N * sizeof(__m512d) / sizeof(double));
		REP(N) min1[kk] = _mm512_set1_pd(val);
		REP(N) midx[kk] = _mm512_set1_epi64(0xFFFFFFFFFFFFFFFF);
		REP(N) nidx[kk] = _mm512_set_epi64(7 + (kk << 3), 6 + (kk << 3), 5 + (kk << 3), 4 + (kk << 3), 3 + (kk << 3), 2 + (kk << 3), 1 + (kk << 3), 0 + (kk << 3));

		for (int64 l1 = n - N * sizeof(__m512d) / sizeof(double); i <= l1; i += N * sizeof(__m512d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm512_loadu_pd(A); A += sizeof(__m512d) / sizeof(double); }

			REP(N) midx[kk] = _mm512_mask_mov_epi64(midx[kk], _mm512_cmp_pd_mask(min1[kk], a[kk], _CMP_GT_OS), nidx[kk]);

			REP(N) min1[kk] = _mm512_min_pd(min1[kk], a[kk]);

			REP(N) nidx[kk] = _mm512_add_epi64(nidx[kk], msep);
		}

		for (int K = sizeof(min1) / sizeof(min1[0]) / 2; K >= 1; K >>= 1)
		{
			REP(K) midx[kk] = _mm512_mask_mov_epi64(midx[kk], _mm512_cmp_pd_mask(min1[kk], min1[kk + K], _CMP_GT_OS), midx[kk + K]);
			REP(K) min1[kk] = _mm512_min_pd(min1[kk], min1[kk + K]);
		}

		for (int64 j = 0; j < N * sizeof(__m512d) / sizeof(double); ++j)
		{
			if (simp_f64(min1, j) > val) continue;
			val = simp_f64(min1, j);
			idx = simp_u64(midx, j);
		}
	}

	for (; i < n; ++i, ++A)
	{
		if (*A > val) continue;
		val = *A;
		idx = i;
	}

	return (int64)idx;
}

TARGET512 int64 GetMinIdx512(float* A, int64 n, float& val)
{
	constexpr int N = 2;
	int64 i = 0;
	val = FLT_MAX;
	uint idx = (uint)-1;

	if (n >= N * sizeof(__m512) / sizeof(float))
	{
		__m512 min1[N], a[N];
		__m512i midx[N], nidx[N], msep = _mm512_set1_epi32(N * sizeof(__m512) / sizeof(float));
		REP(N) min1[kk] = _mm512_set1_ps(val);
		REP(N) midx[kk] = _mm512_set1_epi8((char)0xFF);
		REP(N) nidx[kk] = _mm512_set_epi32(15 + (kk << 4), 14 + (kk << 4), 13 + (kk << 4), 12 + (kk << 4), 11 + (kk << 4), 10 + (kk << 4), 9 + (kk << 4), 8 + (kk << 4), 7 + (kk << 4), 6 + (kk << 4), 5 + (kk << 4), 4 + (kk << 4), 3 + (kk << 4), 2 + (kk << 4), 1 + (kk << 4), 0 + (kk << 4));

		for (int64 l1 = n - N * sizeof(__m512) / sizeof(float); i <= l1; i += N * sizeof(__m512) / sizeof(float))
		{
			REP(N) { a[kk] = _mm512_loadu_ps(A); A += sizeof(__m512) / sizeof(float); }

			REP(N) midx[kk] = _mm512_mask_mov_epi32(midx[kk], _mm512_cmp_ps_mask(min1[kk], a[kk], _CMP_GT_OS), nidx[kk]);

			REP(N) min1[kk] = _mm512_min_ps(min1[kk], a[kk]);

			REP(N) nidx[kk] = _mm512_add_epi64(nidx[kk], msep);
		}

		for (int K = sizeof(min1) / sizeof(min1[0]) / 2; K >= 1; K >>= 1)
		{
			REP(K) midx[kk] = _mm512_mask_mov_epi32(midx[kk], _mm512_cmp_ps_mask(min1[kk], min1[kk + K], _CMP_GT_OS), midx[kk + K]);
			REP(K) min1[kk] = _mm512_min_ps(min1[kk], min1[kk + K]);
		}

		for (int64 j = 0; j < sizeof(__m512) / sizeof(float); ++j)
		{
			if (simp_f32(min1, j) > val) continue;
			val = simp_f32(min1, j);
			idx = simp_u32(midx, j);
		}
	}

	for (; i < n; ++i, ++A)
	{
		if (*A > val) continue;
		val = *A;
		idx = i;
	}

	return (int64)idx;
}

TARGET512 void GetMinMaxVal512(double* A, int64 n, double& minv, double& maxv)
{
	constexpr int N = 4;
	int64 i = 0;
	minv = DBL_MAX;
	maxv = -DBL_MAX;

	if (n >= N * sizeof(__m512d) / sizeof(double))
	{
		__m512d min1[N], max1[N], a[N];
		REP(N) min1[kk] = _mm512_set1_pd(minv);
		REP(N) max1[kk] = _mm512_set1_pd(maxv);

		for (int64 l1 = n - N * sizeof(__m512d) / sizeof(double); i <= l1; i += N * sizeof(__m512d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm512_loadu_pd(A); A += sizeof(__m512d) / sizeof(double); }

			REP(N)
			{
				min1[kk] = _mm512_min_pd(min1[kk], a[kk]);
				max1[kk] = _mm512_max_pd(max1[kk], a[kk]);
			}
		}

		for (int K = sizeof(max1) / sizeof(max1[0]) / 2; K >= 1; K >>= 1)
		{
			REP(K) min1[kk] = _mm512_min_pd(min1[kk], min1[kk + K]);
			REP(K) max1[kk] = _mm512_max_pd(max1[kk], max1[kk + K]);
		}

		minv = _mm512_reduce_min_pd(min1[0]);
		maxv = _mm512_reduce_max_pd(max1[0]);
	}

	for (; i < n; ++i, ++A)
	{
		if (*A < minv) minv = *A;
		if (*A > maxv) maxv = *A;
	}
}

TARGET512 void GetMinMaxVal512(float* A, int64 n, float& minv, float& maxv)
{
	constexpr int N = 4;
	int64 i = 0;
	minv = FLT_MAX;
	maxv = -FLT_MAX;

	if (n >= N * sizeof(__m512) / sizeof(float))
	{
		__m512 min1[N], max1[N], a[N];
		REP(N) min1[kk] = _mm512_set1_ps(minv);
		REP(N) max1[kk] = _mm512_set1_ps(maxv);

		for (int64 l1 = n - N * sizeof(__m512) / sizeof(float); i <= l1; i += N * sizeof(__m512) / sizeof(float))
		{
			REP(N) { a[kk] = _mm512_loadu_ps(A); A += sizeof(__m512) / sizeof(float); }

			REP(N)
			{
				min1[kk] = _mm512_min_ps(min1[kk], a[kk]);
				max1[kk] = _mm512_max_ps(max1[kk], a[kk]);
			}
		}

		for (int K = sizeof(max1) / sizeof(max1[0]) / 2; K >= 1; K >>= 1)
		{
			REP(K) min1[kk] = _mm512_min_ps(min1[kk], min1[kk + K]);
			REP(K) max1[kk] = _mm512_max_ps(max1[kk], max1[kk + K]);
		}

		minv = _mm512_reduce_min_ps(min1[0]);
		maxv = _mm512_reduce_max_ps(max1[0]);
	}

	for (; i < n; ++i, ++A)
	{
		if (*A < minv) minv = *A;
		if (*A > maxv) maxv = *A;
	}
}

TARGET512 double GetMaxVal512(double* A, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;
	double val = -DBL_MAX;

	if (n >= N * sizeof(__m512d) / sizeof(double))
	{
		__m512d max1[N];
		REP(N) max1[kk] = _mm512_set1_pd(val);

		for (int64 l1 = n - N * sizeof(__m512d) / sizeof(double); i <= l1; i += N * sizeof(__m512d) / sizeof(double))
		{
			REP(N)
			{
				max1[kk] = _mm512_max_pd(max1[kk], _mm512_loadu_pd(A)); 
				A += sizeof(__m512d) / sizeof(double);
			}
		}

		for (int K = sizeof(max1) / sizeof(max1[0]) / 2; K >= 1; K >>= 1)
			REP(K) max1[kk] = _mm512_max_pd(max1[kk], max1[kk + K]);

		val = _mm512_reduce_max_pd(max1[0]);
	}

	for (; i < n; ++i, ++A)
	{
		if (*A < val) continue;
		val = *A;
	}
	return val;
}

TARGET512 float GetMaxVal512(float* A, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;
	float val = -FLT_MAX;

	if (n >= N * sizeof(__m512) / sizeof(float))
	{
		__m512 max1[N];
		REP(N) max1[kk] = _mm512_set1_ps(val);

		for (int64 l1 = n - N * sizeof(__m512) / sizeof(float); i <= l1; i += N * sizeof(__m512) / sizeof(float))
		{
			REP(N) { max1[kk] = _mm512_max_ps(max1[kk], _mm512_loadu_ps(A)); A += sizeof(__m512) / sizeof(float); }
		}

		for (int K = sizeof(max1) / sizeof(max1[0]) / 2; K >= 1; K >>= 1)
			REP(K) max1[kk] = _mm512_max_ps(max1[kk], max1[kk + K]);

		__m128* max2 = (__m128*)max1;
		REP(2) max2[kk] = _mm_max_ps(max2[kk], max2[kk + 2]);
		REP(1) max2[kk] = _mm_max_ps(max2[kk], max2[kk + 1]);

		val = Max(Max(simp_f32(max1, 0), simp_f32(max1, 1)),
				  Max(simp_f32(max1, 2), simp_f32(max1, 3)));
	}

	for (; i < n; ++i, ++A)
	{
		if (*A < val) continue;
		val = *A;
	}

	return val;
}

TARGET512 double GetMaxVal512(double* A, int64 n, int64 sep)
{
	constexpr int N = 2;
	int64 i = 0;
	double val = -DBL_MAX;

	if (n >= N * sizeof(__m512d) / sizeof(double))
	{
		__m512i vindex = _mm512_mullo_epi64(_mm512_set_epi64(7, 6, 5, 4, 3, 2, 1, 0), _mm512_set1_epi64(sep));
		__m512d max1[N];
		REP(N) max1[kk] = _mm512_set1_pd(val);

		for (int64 l1 = n - N * sizeof(__m512d) / sizeof(double); i <= l1; i += N * sizeof(__m512d) / sizeof(double))
		{
			REP(N)
			{
				max1[kk] = _mm512_max_pd(max1[kk], _mm512_i64gather_pd(vindex, A, sizeof(double)));
				A += sizeof(__m512d) / sizeof(double) * sep;
			}
		}

		for (int K = sizeof(max1) / sizeof(max1[0]) / 2; K >= 1; K >>= 1)
			REP(K) max1[kk] = _mm512_max_pd(max1[kk], max1[kk + K]);

		val = _mm512_reduce_max_pd(max1[0]);
	}

	for (; i < n; ++i, ++A)
	{
		if (*A < val) continue;
		val = *A;
	}
	return val;
}

TARGET512 float GetMaxVal512(float* A, int64 n, int64 sep)
{
	constexpr int N = 2;
	int64 i = 0;
	float val = -FLT_MAX;

	if (n >= N * sizeof(__m256) / sizeof(float))
	{
		__m512i vindex = _mm512_mullo_epi64(_mm512_set1_epi64(sep), _mm512_set_epi64(7, 6, 5, 4, 3, 2, 1, 0));
		__m256 max1[N];
		REP(N) max1[kk] = _mm256_set1_ps(val);

		for (int64 l1 = n - N * sizeof(__m256) / sizeof(float); i <= l1; i += N * sizeof(__m256) / sizeof(float))
		{
			REP(N) 
			{ 
				max1[kk] = _mm256_max_ps(max1[kk], _mm512_i64gather_ps(vindex, A, sizeof(float))); 
				A += sizeof(__m256) / sizeof(float) * sep; 
			}
		}

		for (int K = sizeof(max1) / sizeof(max1[0]) / 2; K >= 1; K >>= 1)
			REP(K) max1[kk] = _mm256_max_ps(max1[kk], max1[kk + K]);

		__m128* max2 = (__m128*)max1;
		REP(1) max2[kk] = _mm_max_ps(max2[kk], max2[kk + 1]);

		val = Max(Max(simp_f32(max1, 0), simp_f32(max1, 1)),
				  Max(simp_f32(max1, 2), simp_f32(max1, 3)));
	}

	for (; i < n; ++i, A += sep)
	{
		if (*A < val) continue;
		val = *A;
	}

	return val;
}

TARGET512 double GetMinVal512(double* A, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;
	double val = DBL_MAX;

	if (n >= N * sizeof(__m512d) / sizeof(double))
	{
		__m512d min1[N];
		REP(N) min1[kk] = _mm512_set1_pd(val);

		for (int64 l1 = n - N * sizeof(__m512d) / sizeof(double); i <= l1; i += N * sizeof(__m512d) / sizeof(double))
		{
			REP(N)
			{
				min1[kk] = _mm512_min_pd(min1[kk], _mm512_loadu_pd(A));
				A += sizeof(__m512d) / sizeof(double);
			}
		}

		for (int K = sizeof(min1) / sizeof(min1[0]) / 2; K >= 1; K >>= 1)
			REP(K) min1[kk] = _mm512_min_pd(min1[kk], min1[kk + K]);

		val = _mm512_reduce_min_pd(min1[0]);
	}

	for (; i < n; ++i, ++A)
	{
		if (*A > val) continue;
		val = *A;
	}
	return val;
}

TARGET512 float GetMinVal512(float* A, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;
	float val = FLT_MAX;

	if (n >= N * sizeof(__m512) / sizeof(float))
	{
		__m512 min1[N];
		REP(N) min1[kk] = _mm512_set1_ps(val);

		for (int64 l1 = n - N * sizeof(__m512) / sizeof(float); i <= l1; i += N * sizeof(__m512) / sizeof(float))
		{
			REP(N) { min1[kk] = _mm512_min_ps(min1[kk], _mm512_loadu_ps(A)); A += sizeof(__m512) / sizeof(float); }
		}

		for (int K = sizeof(min1) / sizeof(min1[0]) / 2; K >= 1; K >>= 1)
			REP(K) min1[kk] = _mm512_min_ps(min1[kk], min1[kk + K]);

		__m128* min2 = (__m128*)min1;
		REP(2) min2[kk] = _mm_min_ps(min2[kk], min2[kk + 2]);
		REP(1) min2[kk] = _mm_min_ps(min2[kk], min2[kk + 1]);

		val = Min(Min(simp_f32(min1, 0), simp_f32(min1, 1)),
				  Min(simp_f32(min1, 2), simp_f32(min1, 3)));
	}

	for (; i < n; ++i, ++A)
	{
		if (*A > val) continue;
		val = *A;
	}
	return val;
}

TARGET512 int64 GetMinVal512(int64* A, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;
	int64 val = 0x7FFFFFFFFFFFFFFF;

	if (n >= N * sizeof(__m512i) / sizeof(int64))
	{
		__m512i min1[N];
		REP(N) min1[kk] = _mm512_set1_epi64(0x7FFFFFFFFFFFFFFF);

		for (int64 l1 = n - N * sizeof(__m512i) / sizeof(int64); i <= l1; i += N * sizeof(__m512i) / sizeof(int64))
		{
			REP(N) 
			{ 
				min1[kk] = _mm512_min_epi64(min1[kk], _mm512_loadu_si512(A)); 
				A += sizeof(__m512i) / sizeof(int64); 
			}
		}

		for (int K = sizeof(min1) / sizeof(min1[0]) / 2; K >= 1; K >>= 1)
			REP(K) min1[kk] = _mm512_min_epi64(min1[kk], min1[kk + K]);

		val = _mm512_reduce_min_epi64(min1[0]);
	}

	for (; i < n; ++i, ++A)
	{
		if (*A > val) continue;
		val = *A;
	}

	return val;
}

TARGET512 void SetVal512(uint* A, ushort* B, int64 n)
{
	int64 i = 0;

	if (n >= 32)
	{
		__m512i b;
		__m256i& b1 = *((__m256i*) & b + 0);
		__m256i& b2 = *((__m256i*) & b + 1);
		__m512i v1, v2;

		for (int64 l1 = n - 32; i <= l1; i += 32)
		{
			b = _mm512_loadu_si512(B); B += 32;

			v1 = _mm512_cvtepu16_epi32(b1);
			v2 = _mm512_cvtepu16_epi32(b2);

			_mm512_storeu_si512(A, v1); A += 16;
			_mm512_storeu_si512(A, v2); A += 16;
		}
	}

	for (; i < n; ++i)
		*A++ = *B++;
}

TARGET512 void AddExponent512(int64& slog, __m512d& val)
{
	__m512i& vv = *(__m512i*)&val;
	__m512i mask1 = _mm512_set1_epi64(0x7FF0000000000000);//3+4+4 = 11 bits, 1 + 11 + 52
	__m512i mask2 = _mm512_set1_epi64(0x800FFFFFFFFFFFFF);
	__m512i mask3 = _mm512_set1_epi64(0x3FF0000000000000);
	__m512i subv = _mm512_set1_epi64(1023);// 2^11 = 2048

	__m512i t = _mm512_sub_epi64(_mm512_srli_epi64(_mm512_and_si512(vv, mask1), 52), subv);

	slog += _mm512_reduce_add_epi64(t);
	vv = _mm512_or_si512(_mm512_and_si512(vv, mask2), mask3);
}

TARGET512 void AddExponent512(int64& slog, __m512& val)
{
	__m512i& vv = *(__m512i*)&val;
	__m512i mask1 = _mm512_set1_epi32(0x7F800000);        //3+4+1 = 8 bits, 1 + 8 + 23
	__m512i mask2 = _mm512_set1_epi32(0x807FFFFF);
	__m512i mask3 = _mm512_set1_epi32(0x3F800000);
	__m512i subv = _mm512_set1_epi32(127);// 2^8 = 256

	__m512i t = _mm512_sub_epi32(_mm512_srli_epi32(_mm512_and_si512(vv, mask1), 23), subv);

	slog += _mm512_reduce_add_epi32(t);
	vv = _mm512_or_si512(_mm512_and_si512(vv, mask2), mask3);
}

TARGET512 void ChargeLog512(int64& slog, double& prod, __m512d& val)
{
	AddExponent512(slog, val);

	prod = prod * __mm512_reduce_mul_pd(val);

	if (prod < DOUBLE_UNDERFLOW || prod > DOUBLE_OVERFLOW) [[unlikely]]
		AddExponent(slog, prod);
}

TARGET512 void ChargeLog512(int64& slog, double& prod, __m512& val)
{
	AddExponent512(slog, val);

	prod = prod * __mm512_reduce_mul_psd(val);

	if (prod < DOUBLE_UNDERFLOW || prod > DOUBLE_OVERFLOW) [[unlikely]]
		AddExponent(slog, prod);
}

TARGET512 double LogProd512(double* A, int64 n)
{
	return LogProdAVX(A, n);
	
	/*
	int64 i = 0;
	int64 slog = 0; double prod = 1;

	if (n >= 32)
	{
		__m512d pd = _mm512_set1_pd(1.0), dunder = _mm512_set1_pd(DOUBLE_UNDERFLOW), dover = _mm512_set1_pd(DOUBLE_OVERFLOW);
		__mmask8 flag;

		for (int64 l1 = n - 32; i <= l1; i += 32)
		{
			REP(2)
			{
				pd = _mm512_mul_pd(pd, _mm512_mul_pd(_mm512_loadu_pd(A), _mm512_loadu_pd(A + 8)));
				A += 16;
			}

			//
			flag = _mm512_cmp_pd_mask(pd, dunder, _CMP_LT_OS) | _mm512_cmp_pd_mask(dover, pd, _CMP_LT_OS);

			//if (flag) [[unlikely]]

			AddExponent512(slog, pd);
		}

		__m128d* pd1 = (__m128d*)&pd;
		REP(4) ChargeLogSSE(slog, prod, pd1[kk]);
	}

	for (; i < n; ++i, ++A)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);
	return prod;
	*/
}

TARGET512 double LogProd512(float* A, int64 n)
{
	return LogProdAVX(A, n);

	/*
	int64 i = 0;
	int64 slog = 0; double prod = 1;

	if (n >= 64)
	{
		__m512d pd = _mm512_set1_pd(1.0);

		for (int64 l1 = n - 64; i <= l1; i += 64)
		{
			REP(4)
			{
				pd = _mm512_mul_pd(pd, _mm512_mul_pd(
						_mm512_cvtps_pd(_mm256_loadu_ps(A)),
						_mm512_cvtps_pd(_mm256_loadu_ps(A + 8))));
				A += 16;
			}

			AddExponent512(slog, pd);
		}

		__m128d* pd1 = (__m128d*) & pd;
		REP(4) ChargeLogSSE(slog, prod, pd1[kk]);
	}

	for (; i < n; ++i, ++A)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);
	return prod;
	*/
}

TARGET512 float LogProd512x(float* A, int64 n)
{
	int64 i = 0;
	int64 slog = 0; double prod = 1;

	if (n >= 32)
	{
		__m512 pd = _mm512_set1_ps(1.0f);

		for (int64 l1 = n - 32; i <= l1; i += 32)
		{
			pd = _mm512_mul_ps(pd, _mm512_loadu_ps(A));
			pd = _mm512_mul_ps(pd, _mm512_loadu_ps(A + 16));

			A += 32;

			AddExponent512(slog, pd);
		}

		__m128* pd1 = (__m128*)&pd;
		REP(4) ChargeLogSSE(slog, prod, pd1[kk]);
	}

	for (; i < n; ++i, ++A)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);
	return prod;
}

TARGET512 double LogProd512(double* A, int64 n, int64 sep)
{
	int64 i = 0;
	int64 slog = 0; double prod = 1;

	if (n >= 8)
	{
		_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
		_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
		_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
		_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
		_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
		_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
		_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
		_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;

		__m512d pd = _mm512_set1_pd(1.0), dunder = _mm512_set1_pd(DOUBLE_UNDERFLOW), dover = _mm512_set1_pd(DOUBLE_OVERFLOW);
		__m512i vindex = _mm512_mullo_epi64(_mm512_set_epi64(-9, -10, -11, -12, -13, -14, -15, -16), _mm512_set1_epi64(sep));
		__mmask8 flag;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;

			pd = _mm512_mul_pd(pd, _mm512_i64gather_pd(vindex, A, sizeof(double)));

			flag = _mm512_cmp_pd_mask(pd, dunder, _CMP_LT_OS) | _mm512_cmp_pd_mask(dover, pd, _CMP_LT_OS);

			if (flag) [[unlikely]]
				AddExponent512(slog, pd);
		}

		__m128d* pd1 = (__m128d*)&pd;
		REP(4) ChargeLogSSE(slog, prod, pd1[kk]);
		A -= 8 * sep;
	}

	for (; i < n; ++i, A += sep)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);
	return prod;
}

TARGET512 double LogProd512(float* A, int64 n, int64 sep)
{
	return LogProdSSE(A, n, sep);
}

TARGET512 float LogProd512x(float* A, int64 n, int64 sep)
{
	return LogProdSSEx(A, n, sep);
	/*
	int64 i = 0;
	int64 slog = 0; double prod = 1;

	if (n >= 32)
	{
		__m512 pd = _mm512_set1_ps(1.0f);
		__m512i vindex = _mm512_mullo_epi32(_mm512_set1_epi32(sep), _mm512_set_epi32(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0)); xxx

		for (int64 l1 = n - 32; i <= l1; i += 32)
		{
			pd = _mm512_mul_ps(pd, _mm512_i32xxgather_ps(vindex, A, sizeof(float)));
			A += 16 * sep;

			pd = _mm512_mul_ps(pd, _mm512_i32xxgather_ps(vindex, A, sizeof(float)));
			A += 16 * sep;

			AddExponent512(slog, pd);
		}

		__m128* pd1 = (__m128*)&pd;
		REP(4) ChargeLogSSE(slog, prod, pd1[kk]);
	}

	for (; i < n; ++i, A += sep)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);
	return prod;
	*/
}

TARGET512 double LogProdDiv512(double* A, double* B, int64 n, int64 sep)
{
	int64 slog1 = 0; double prod1 = 1;
	int64 slog2 = 0; double prod2 = 1;

	for (int64 i = 0; i < n; ++i, A += sep, B += sep)
	{
		ChargeLog(slog1, prod1, *A);
		ChargeLog(slog2, prod2, *B);
	}

	CloseLog(slog1, prod1);
	CloseLog(slog2, prod2);
	return prod1 - prod2;

	/*
	int64 i = 0;
	int64 slog1 = 0; double prod1 = 1;
	int64 slog2 = 0; double prod2 = 1;

	if (n >= 64)
	{
		__m512d pd1 = _mm512_set1_pd(1.0), pd2 = _mm512_set1_pd(1.0);
		__m512i vindex = _mm512_mullo_epi64(_mm512_set_epi64(7, 6, 5, 4, 3, 2, 1, 0), _mm512_set1_epi64(sep));

		for (int64 l1 = n - 64; i <= l1; i += 64)
		{
			REP(8)
			{
				pd1 = _mm512_mul_pd(pd1, _mm512_i64gather_pd(vindex, A, sizeof(double))); A += sep * 8;
				pd2 = _mm512_mul_pd(pd2, _mm512_i64gather_pd(vindex, B, sizeof(double))); B += sep * 8;
			}

			AddExponent512(slog1, pd1);
			AddExponent512(slog2, pd2);
		}

		__m128d* pd11 = (__m128d*)&pd1;
		__m128d* pd22 = (__m128d*)&pd2;
		REP(4) ChargeLogSSE(slog1, prod1, pd11[kk]);
		REP(4) ChargeLogSSE(slog2, prod2, pd22[kk]);
	}

	for (; i < n; ++i, A += sep, B += sep)
	{
		ChargeLog(slog1, prod1, *A);
		ChargeLog(slog2, prod2, *B);
	}

	CloseLog(slog1, prod1);
	CloseLog(slog2, prod2);

	return prod1 - prod2;
	*/
}

TARGET512 double LogProdDiv512(float* A, float* B, int64 n, int64 sep)
{
	int64 i = 0;
	int64 slog1 = 0; double prod1 = 1;
	int64 slog2 = 0; double prod2 = 1;

	if (n >= 64)
	{
		__m512d pd1 = _mm512_set1_pd(1.0), pd2 = _mm512_set1_pd(1.0);

		for (int64 l1 = n - 64; i <= l1; i += 64)
		{
			REP(8)
			{
				pd1 = _mm512_mul_pd(pd1, _mm512_set_pd(A[7 * sep], A[6 * sep], A[5 * sep], A[4 * sep], A[3 * sep], A[2 * sep], A[1 * sep], A[0 * sep])); A += 8 * sep;
				pd2 = _mm512_mul_pd(pd2, _mm512_set_pd(B[7 * sep], B[6 * sep], B[5 * sep], B[4 * sep], B[3 * sep], B[2 * sep], B[1 * sep], B[0 * sep])); B += 8 * sep;
			}

			AddExponent512(slog1, pd1);
			AddExponent512(slog2, pd2);
		}

		__m128d* pd11 = (__m128d*)&pd1;
		__m128d* pd22 = (__m128d*)&pd2;

		REP(4) ChargeLogSSE(slog1, prod1, pd11[kk]);
		REP(4) ChargeLogSSE(slog2, prod2, pd22[kk]);
	}

	for (; i < n; ++i, A += sep, B += sep)
	{
		ChargeLog(slog1, prod1, *A);
		ChargeLog(slog2, prod2, *B);
	}

	CloseLog(slog1, prod1);
	CloseLog(slog2, prod2);

	return prod1 - prod2;
}

TARGET512 float LogProdDiv512x(float* A, float* B, int64 n, int64 sep)
{
	int64 i = 0;
	int64 slog1 = 0; double prod1 = 1;
	int64 slog2 = 0; double prod2 = 1;

	if (n >= 64)
	{
		__m512 pd1 = _mm512_set1_ps(1.0f), pd2 = _mm512_set1_ps(1.0f);

		for (int64 l1 = n - 64; i <= l1; i += 64)
		{
			REP(4)
			{
				pd1 = _mm512_mul_ps(pd1, _mm512_set_ps(
					A[15 * sep], A[14 * sep], A[13 * sep], A[12 * sep], A[11 * sep], A[10 * sep], A[9 * sep], A[8 * sep], A[7 * sep], A[6 * sep], A[5 * sep], A[4 * sep], A[3 * sep], A[2 * sep], A[1 * sep], A[0 * sep]));
				pd2 = _mm512_mul_ps(pd2, _mm512_set_ps(
					B[15 * sep], B[14 * sep], B[13 * sep], B[12 * sep], B[11 * sep], B[10 * sep], B[9 * sep], B[8 * sep], B[7 * sep], B[6 * sep], B[5 * sep], B[4 * sep], B[3 * sep], B[2 * sep], B[1 * sep], B[0 * sep]));
				A += 16 * sep; B += 16 * sep;
			}

			AddExponent512(slog1, pd1);
			AddExponent512(slog2, pd2);
		}

		__m128* pd11 = (__m128*)&pd1;
		__m128* pd22 = (__m128*)&pd2;

		REP(4) ChargeLogSSE(slog1, prod1, pd11[kk]);
		REP(4) ChargeLogSSE(slog2, prod2, pd22[kk]);
	}

	for (; i < n; ++i, A += sep, B += sep)
	{
		ChargeLog(slog1, prod1, *A);
		ChargeLog(slog2, prod2, *B);
	}

	CloseLog(slog1, prod1);
	CloseLog(slog2, prod2);

	return prod1 - prod2;
}

TARGET512 int64 CountNonZero512(byte* A, int64 n)
{
	int64 i = 0;
	uint64 r1 = 0, r2 = 0, r3 = 0, r4 = 0;

	if (n >= 256)
	{
		uint64 x1 = 0, x2 = 0, x3 = 0, x4 = 0;
		__m512i z = _mm512_setzero_si512();
		__m512i a1, a2, a3, a4;

		for (int64 l1 = n - 256; i <= l1; i += 256)
		{
			a1 = _mm512_loadu_si512(A); A += 64;
			a2 = _mm512_loadu_si512(A); A += 64;
			a3 = _mm512_loadu_si512(A); A += 64;
			a4 = _mm512_loadu_si512(A); A += 64;

			x1 = _mm512_cmpgt_epu8_mask(a1, z);
			x2 = _mm512_cmpgt_epu8_mask(a2, z);
			x3 = _mm512_cmpgt_epu8_mask(a3, z);
			x4 = _mm512_cmpgt_epu8_mask(a4, z);

			r1 += _mm_popcnt_u64(x1);
			r2 += _mm_popcnt_u64(x2);
			r3 += _mm_popcnt_u64(x3);
			r4 += _mm_popcnt_u64(x4);
		}

		r1 = r1 + r2 + r3 + r4;
	}

	for (; i < n; ++i, ++A)
		if (*A) r1++;

	return (int64)r1;
}

TARGET512 double Sum512(double* A, int64 n)
{
	constexpr int N = 2;
	int64 i = 0;
	double re = 0;

	if (n >= N * sizeof(__m512d) / sizeof(double))
	{
		__m512d s[N], a[N];
		REP(N) s[kk] = _mm512_setzero_pd();

		for (int64 l1 = n - N * sizeof(__m512d) / sizeof(double); i <= l1; i += N * sizeof(__m512d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm512_loadu_pd(A); A += sizeof(__m512d) / sizeof(double); }

			REP(N) s[kk] = _mm512_add_pd(s[kk], a[kk]);
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm512_add_pd(s[kk], s[kk + K]);

		re = __mm512_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *A++;
		re += v1;
	}

	return re;
}

TARGET512 double Sum512(float* A, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;
	double re = 0;

	if (n >= N * sizeof(__m256) / sizeof(float))
	{
		__m512d s[N], a[N];
		REP(N) s[kk] = _mm512_setzero_pd();

		for (int64 l1 = n - N * sizeof(__m256) / sizeof(float); i <= l1; i += N * sizeof(__m256) / sizeof(float))
		{
			REP(N) { a[kk] = _mm512_cvtps_pd(_mm256_loadu_ps(A)); A += sizeof(__m256) / sizeof(float); }

			REP(N) s[kk] = _mm512_add_pd(s[kk], a[kk]);
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm512_add_pd(s[kk], s[kk + K]);

		re = __mm512_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A++;
		re += v1;
	}

	return re;
}

TARGET512 float Sum512x(float* A, int64 n)
{
	constexpr int N = 2;
	int64 i = 0;
	float re = 0;

	if (n >= N * sizeof(__m512) / sizeof(float))
	{
		__m512 s[N], a[N];
		REP(N) s[kk] = _mm512_setzero_ps();

		for (int64 l1 = n - N * sizeof(__m512) / sizeof(float); i <= l1; i += N * sizeof(__m512) / sizeof(float))
		{
			REP(N) { a[kk] = _mm512_loadu_ps(A); A += sizeof(__m512) / sizeof(float); }

			REP(N) s[kk] = _mm512_add_ps(s[kk], a[kk]);
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm512_add_ps(s[kk], s[kk + K]);

		re = __mm512_reduce_add_ps(s[0]);
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A++;
		re += v1;
	}

	return re;
}

TARGET512 int64 Sum512(byte* A, int64 n)
{
	constexpr int N = 4;
	uint64 re = 0;
	int64 i = 0;

	if (n >= N * sizeof(__m512i) / sizeof(byte))
	{
		__m512i s[N], a[N], z = _mm512_setzero_si512();
		REP(N) s[kk] = _mm512_setzero_si512();

		for (int64 l1 = n - N * sizeof(__m512i) / sizeof(byte); i <= l1; i += N * sizeof(__m512i) / sizeof(byte))
		{
			REP(N) { a[kk] = _mm512_loadu_si512(A); A += sizeof(__m512i) / sizeof(byte); }

			REP(N) { s[kk] = _mm512_add_epi64(s[kk], _mm512_sad_epu8(a[kk], z)); }
		}

		s[0] = _mm512_add_epi64(_mm512_add_epi64(s[0], s[1]), _mm512_add_epi64(s[2], s[3]));
		re += _mm512_reduce_add_epi64(s[0]);
	}

	for (; i < n; ++i)
		re += *A++;

	return re;
}

TARGET512 double Sum512(double* A, int64 n, int64 sep)
{
	constexpr int N = 2;
	int64 i = 0;
	double re = 0;

	if (n >= N * sizeof(__m512d) / sizeof(double))
	{
		__m512d s[N], a[N];
		REP(N) s[kk] = _mm512_setzero_pd();
		__m512i vindex = _mm512_mullo_epi64(_mm512_set_epi64(7, 6, 5, 4, 3, 2, 1, 0), _mm512_set1_epi64(sep));

		for (int64 l1 = n - N * sizeof(__m512d) / sizeof(double); i <= l1; i += N * sizeof(__m512d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm512_i64gather_pd(vindex, A, sizeof(double)); A += sizeof(__m512d) / sizeof(double) * sep; }

			REP(N) s[kk] = _mm512_add_pd(s[kk], a[kk]);
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm512_add_pd(s[kk], s[kk + K]);

		re = __mm512_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, A += sep)
	{
		volatile double v1 = *A;
		re += v1;
	}

	return re;
}

TARGET512 double Sum512(float* A, int64 n, int64 sep)
{
	return SumSSE(A, n, sep);

	/*
	constexpr int N = 4;
	int64 i = 0;
	double re = 0;

	if (n >= N * sizeof(__m256) / sizeof(float))
	{
		__m512d s[N];
		__m256 a[N];
		REP(N) s[kk] = _mm512_setzero_pd();
		__m512i vindex = _mm512_mullo_epi64(_mm512_set1_epi64(sep), _mm512_set_epi64(7, 6, 5, 4, 3, 2, 1, 0));

		for (int64 l1 = n - N * sizeof(__m256) / sizeof(float); i <= l1; i += N * sizeof(__m256) / sizeof(float))
		{
			REP(N)
			{
				a[kk] = _mm512_i64gather_ps(vindex, A, sizeof(float));
				s[kk] = _mm512_add_pd(s[kk], _mm512_cvtps_pd(a[kk]));
				A += sep * sizeof(__m256) / sizeof(float);
			}
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm512_add_pd(s[kk], s[kk + K]);

		re = __mm512_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, A += sep)
	{
		volatile float v1 = *A;
		re += v1;
	}

	return re;
	*/
}

TARGET512 float Sum512x(float* A, int64 n, int64 sep)
{
	return SumSSEx(A, n, sep);
	/*

	constexpr int N = 2;
	int64 i = 0;
	float re = 0;

	if (n >= N * sizeof(__m512) / sizeof(float))
	{
		__m512 s[N];
		REP(N) s[kk] = _mm512_setzero_ps();
		__m512i vindex = _mm512_mullo_epi32(_mm512_set_epi32(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0), _mm512_set1_epi32(sep));

		for (int64 l1 = n - N * sizeof(__m512) / sizeof(float); i <= l1; i += N * sizeof(__m512) / sizeof(float))
		{
			REP(N) { s[kk] = _mm512_add_ps(s[kk], _mm512_i32xxgather_ps(vindex, A, sizeof(float))); A += sep * sizeof(__m512) / sizeof(float); }
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm512_add_ps(s[kk], s[kk + K]);

		re = __mm512_reduce_add_ps(s[0]);
	}

	for (; i < n; ++i, A += sep)
	{	
		volatile float v1 = *A;
		re += v1;
	}

	return re;
	*/
}

TARGET512 void Sum512(double* A, double** B, int64 k, int64 n)
{
	constexpr int N = 16;
	int64 i = 0;

	if (n >= N * sizeof(__m512d) / sizeof(double))
	{
		__m512d a[N];

		for (int64 l1 = n - N * sizeof(__m512d) / sizeof(double); i <= l1; i += N * sizeof(__m512d) / sizeof(double))
		{
			REP(N) a[kk] = _mm512_setzero_pd();

			for (int64 j = 0; j < k; ++j)
				REP(N) a[kk] = _mm512_add_pd(a[kk], _mm512_loadu_pd(&B[j][i + kk * sizeof(__m512d) / sizeof(double)]));

			REP(N) { _mm512_storeu_pd(A, a[kk]); A += sizeof(__m512d) / sizeof(double); }
		}
	}

	for (; i < n; ++i)
	{
		double Ai = 0;
		for (int64 j = 0; j < k; ++j)
			Ai += B[j][i];
		*A++ = Ai;
	}
}

TARGET512 void Sum512(float* A, float** B, int64 k, int64 n)
{
	constexpr int N = 16;
	int64 i = 0;

	if (n >= N * sizeof(__m256) / sizeof(float))
	{
		__m512d a[N];

		for (int64 l1 = n - N * sizeof(__m256) / sizeof(float); i <= l1; i += N * sizeof(__m256) / sizeof(float))
		{
			REP(N) a[kk] = _mm512_setzero_pd();

			for (int64 j = 0; j < k; ++j)
				REP(N) a[kk] = _mm512_add_pd(a[kk], _mm512_cvtps_pd(_mm256_loadu_ps(&B[j][i + kk * sizeof(__m256) / sizeof(float)])));

			REP(N) { _mm256_storeu_ps(A, _mm512_cvtpd_ps(a[kk])); A += sizeof(__m256) / sizeof(float); }
		}
	}

	for (; i < n; ++i)
	{
		double Ai = 0;
		for (int64 j = 0; j < k; ++j)
			Ai += B[j][i];
		*A++ = Ai;
	}
}

TARGET512 double Prod512(double* A, int64 n)
{
	constexpr int N = 2;
	int64 i = 0;
	double re = 1;

	if (n >= N * sizeof(__m512d) / sizeof(double))
	{
		__m512d s[N], a[N];
		REP(N) s[kk] = _mm512_set1_pd(1);

		for (int64 l1 = n - N * sizeof(__m512d) / sizeof(double); i <= l1; i += N * sizeof(__m512d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm512_loadu_pd(A); A += sizeof(__m512d) / sizeof(double); }

			REP(N) s[kk] = _mm512_mul_pd(s[kk], a[kk]);
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm512_mul_pd(s[kk], s[kk + K]);

		re = __mm512_reduce_mul_pd(s[0]);
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *A++;
		re *= v1;
	}

	return re;
}

TARGET512 double Prod512(float* A, int64 n)
{
	return ProdAVX(A, n);

	/*
	constexpr int N = 4;
	int64 i = 0;
	double re = 0;

	if (n >= N * sizeof(__m256) / sizeof(float))
	{
		__m512d s[N], a[N];
		REP(N) s[kk] = _mm512_set1_pd(1);

		for (int64 l1 = n - N * sizeof(__m256) / sizeof(float); i <= l1; i += N * sizeof(__m256) / sizeof(float))
		{
			REP(N) { a[kk] = _mm512_cvtps_pd(_mm256_loadu_ps(A)); A += sizeof(__m256) / sizeof(float); }

			REP(N) s[kk] = _mm512_mul_pd(s[kk], a[kk]);
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm512_mul_pd(s[kk], s[kk + K]);

		re = __mm512_reduce_mul_pd(s[0]);
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A++;
		re *= v1;
	}

	return re;
	*/
}

TARGET512 float Prod512x(float* A, int64 n)
{
	constexpr int N = 2;
	int64 i = 0;
	volatile float re = 1;

	if (n >= N * sizeof(__m512) / sizeof(float))
	{
		__m512 s[N], a[N];
		REP(N) s[kk] = _mm512_set1_ps(1);

		for (int64 l1 = n - N * sizeof(__m512) / sizeof(float); i <= l1; i += N * sizeof(__m512) / sizeof(float))
		{
			REP(N) { a[kk] = _mm512_loadu_ps(A); A += sizeof(__m512) / sizeof(float); }

			REP(N) s[kk] = _mm512_mul_ps(s[kk], a[kk]);
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm512_mul_ps(s[kk], s[kk + K]);

		re = __mm512_reduce_mul_ps(s[0]);
	}

	for (; i < n; ++i)
		re *= *A++;

	return re;
}

TARGET512 double Prod512(double* A, int64 n, int64 sep)
{
	return ProdSSE(A, n, sep);
}

TARGET512 double Prod512(float* A, int64 n, int64 sep)
{
	return ProdSSE(A, n, sep);
}

TARGET512 float Prod512x(float* A, int64 n, int64 sep)
{
	return ProdSSEx(A, n, sep);
}

TARGET512 double SumSquare512(double* A, int64 n)
{
	constexpr int N = 2;
	int64 i = 0;
	double re = 0;

	if (n >= N * sizeof(__m512d) / sizeof(double))
	{
		__m512d s[N], a[N];
		REP(N) s[kk] = _mm512_setzero_pd();

		for (int64 l1 = n - N * sizeof(__m512d) / sizeof(double); i <= l1; i += N * sizeof(__m512d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm512_loadu_pd(A); A += sizeof(__m512d) / sizeof(double); }

			REP(N) a[kk] = _mm512_mul_pd(a[kk], a[kk]);

			REP(N) s[kk] = _mm512_add_pd(s[kk], a[kk]);
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm512_add_pd(s[kk], s[kk + K]);

		re = __mm512_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, ++A)
	{
		volatile double v1 = *A * *A;
		re += v1;
	}

	return re;
}

TARGET512 double SumSquare512(float* A, int64 n)
{
	return SumSquareAVX(A, n);
	/*
	constexpr int N = 2;
	int64 i = 0;
	double re = 0;

	if (n >= N * sizeof(__m256) / sizeof(float))
	{
		__m512d s[N], a[N];
		REP(N) s[kk] = _mm512_setzero_pd();

		for (int64 l1 = n - N * sizeof(__m256) / sizeof(float); i <= l1; i += N * sizeof(__m256) / sizeof(float))
		{
			REP(N) { a[kk] = _mm512_cvtps_pd(_mm256_loadu_ps(A)); A += sizeof(__m256) / sizeof(float); }
			REP(N) a[kk] = _mm512_mul_pd(a[kk], a[kk]);
			REP(N) s[kk] = _mm512_add_pd(s[kk], a[kk]);
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm512_add_pd(s[kk], s[kk + K]);

		re = __mm512_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, ++A)
	{
		volatile double v1 = (double)*A * (double)*A;
		re += v1;
	}

	return re;
	*/
}

TARGET512 float SumSquare512x(float* A, int64 n)
{
	constexpr int N = 2;
	int64 i = 0;
	float re = 0;

	if (n >= N * sizeof(__m512) / sizeof(float))
	{
		__m512 s[N], a[N];
		REP(N) s[kk] = _mm512_setzero_ps();

		for (int64 l1 = n - N * sizeof(__m512) / sizeof(float); i <= l1; i += N * sizeof(__m512) / sizeof(float))
		{
			REP(N) { a[kk] = _mm512_loadu_ps(A); A += sizeof(__m512) / sizeof(float); }

			REP(N) a[kk] = _mm512_mul_ps(a[kk], a[kk]);

			REP(N) s[kk] = _mm512_add_ps(s[kk], a[kk]);
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm512_add_ps(s[kk], s[kk + K]);

		re = __mm512_reduce_add_ps(s[0]);
	}

	for (; i < n; ++i, ++A)
	{
		volatile float v1 = *A * *A;
		re += v1;
	}

	return re;
}

TARGET512 int64 SumSquare512(byte* A, int64 n)
{
	constexpr int N = 2;
	int64 i = 0;
	uint64 re = 0;

	if (n >= sizeof(__m512i) / sizeof(byte))
	{
		__m512i a, t[2], s = _mm512_setzero_si512();
		REP(2) t[kk] = _mm512_setzero_si512();
		__m256i* s2 = (__m256i*)&s;

		for (int64 l1 = n - N * sizeof(__m512i) / sizeof(byte); i <= l1; i += N * sizeof(__m512i) / sizeof(byte))
		{
			REP(N)
			{
				a = _mm512_loadu_si512((__m512i*)A);
				A += sizeof(__m512i) / sizeof(byte);
				s = _mm512_add_epi16(s, _mm512_maddubs_epi16(a, a));
			}

			if ((i & (sizeof(__m512i) / sizeof(byte) * 128 - 1)) == 0)
			{
				t[0] = _mm512_add_epi32(t[0], _mm512_cvtepi16_epi32(s2[0]));
				t[1] = _mm512_add_epi32(t[1], _mm512_cvtepi16_epi32(s2[1]));
				s = _mm512_setzero_si512();
			}
		}

		t[0] = _mm512_add_epi32(t[0], _mm512_cvtepi16_epi32(s2[0]));
		t[1] = _mm512_add_epi32(t[1], _mm512_cvtepi16_epi32(s2[1]));
		t[0] = _mm512_add_epi32(t[0], t[1]);

		re = _mm512_reduce_add_epi32(t[0]);
	}

	for (; i < n; ++i, ++A)
		re += *A * *A;

	return re;
}

TARGET512 void SumSumSquare512(double* A, int64 n, double& sum, double& sumsq)
{
	constexpr int N = 2;
	int64 i = 0;
	double re1 = 0, re2 = 0;

	if (n >= N * sizeof(__m512d) / sizeof(double))
	{
		__m512d s1[N], s2[N], a[N];
		REP(N) s1[kk] = s2[kk] = _mm512_setzero_pd();

		for (int64 l1 = n - N * sizeof(__m512d) / sizeof(double); i <= l1; i += N * sizeof(__m512d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm512_loadu_pd(A); A += sizeof(__m512d) / sizeof(double); }

			REP(N) s1[kk] = _mm512_add_pd(s1[kk], a[kk]);

			REP(N) a[kk] = _mm512_mul_pd(a[kk], a[kk]);

			REP(N) s2[kk] = _mm512_add_pd(s2[kk], a[kk]);
		}

		for (int K = sizeof(s1) / sizeof(s1[0]) / 2; K >= 1; K >>= 1)
		{
			REP(K) s1[kk] = _mm512_add_pd(s1[kk], s1[kk + K]);
			REP(K) s2[kk] = _mm512_add_pd(s2[kk], s2[kk + K]);
		}

		re1 = __mm512_reduce_add_pd(s1[0]);
		re2 = __mm512_reduce_add_pd(s2[0]);
	}

	for (; i < n; ++i, ++A)
	{
		volatile double v1 = *A;
		volatile double v2 = *A * *A;
		re1 += v1;
		re2 += v2;
	}

	sum = re1;
	sumsq = re2;
}

TARGET512 void SumSumSquare512(float* A, int64 n, double& sum, double& sumsq)
{
	constexpr int N = 4;
	int64 i = 0;
	double re1 = 0, re2 = 0;

	if (n >= N * sizeof(__m256) / sizeof(float))
	{
		__m512d s1[N], s2[N], a[N];
		REP(N) s1[kk] = s2[kk] = _mm512_setzero_pd();

		for (int64 l1 = n - N * sizeof(__m256) / sizeof(float); i <= l1; i += N * sizeof(__m256) / sizeof(float))
		{
			REP(N) { a[kk] = _mm512_cvtps_pd(_mm256_loadu_ps(A)); A += sizeof(__m256) / sizeof(float); }

			REP(N) s1[kk] = _mm512_add_pd(s1[kk], a[kk]);

			REP(N) a[kk] = _mm512_mul_pd(a[kk], a[kk]);

			REP(N) s2[kk] = _mm512_add_pd(s2[kk], a[kk]);
		}

		for (int K = sizeof(s1) / sizeof(s1[0]) / 2; K >= 1; K >>= 1)
		{
			REP(K) s1[kk] = _mm512_add_pd(s1[kk], s1[kk + K]);
			REP(K) s2[kk] = _mm512_add_pd(s2[kk], s2[kk + K]);
		}

		re1 = __mm512_reduce_add_pd(s1[0]);
		re2 = __mm512_reduce_add_pd(s2[0]);
	}

	for (; i < n; ++i, ++A)
	{
		volatile double v1 = (double)*A;
		volatile double v2 = (double)*A * (double)*A;
		re1 += v1;
		re2 += v2;
	}

	sum = re1;
	sumsq = re2;
}

TARGET512 double SumProdDiv512(double* A1, double* A2, double* B, int64 sep, int64 n)
{
	return SumProdDivAVX(A1, A2, B, sep, n);
}

TARGET512 double SumProdDiv512(double* A1, float* A2, float* B, int64 sep, int64 n)
{
	return SumProdDivSSE(A1, A2, B, sep, n);
}

TARGET512 double SumProdDiv512(float* A1, float* A2, float* B, int64 sep, int64 n)
{
	return SumProdDivSSE(A1, A2, B, sep, n);
}

TARGET512 float SumProdDiv512x(float* A1, float* A2, float* B, int64 sep, int64 n)
{
	return SumProdDivSSEx(A1, A2, B, sep, n);
}

TARGET512 double SumProd512(double* A, double* B, int64 sep, int64 n)
{
	return SumProdSSE(A, B, sep, n);
}

TARGET512 double SumProd512(float* A, float* B, int64 sep, int64 n)
{
	return SumProdSSE(A, B, sep, n);
}

TARGET512 float SumProd512x(float* A, float* B, int64 sep, int64 n)
{
	return SumProdSSEx(A, B, sep, n);
}

TARGET512 double SumProd512(double* A, double* B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * sizeof(__m512d) / sizeof(double))
	{
		__m512d s[N], a[N], b[N];
		REP(N) s[kk] = _mm512_setzero_pd();

		for (int64 l1 = n - N * sizeof(__m512d) / sizeof(double); i <= l1; i += N * sizeof(__m512d) / sizeof(double))
		{
			REP(N)
			{
				a[kk] = _mm512_loadu_pd(A); A += sizeof(__m512d) / sizeof(double);
				b[kk] = _mm512_loadu_pd(B); B += sizeof(__m512d) / sizeof(double);
			}

			REP(N) a[kk] = _mm512_mul_pd(a[kk], b[kk]);

			REP(N) s[kk] = _mm512_add_pd(s[kk], a[kk]);
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm512_add_pd(s[kk], s[kk + K]);

		re = __mm512_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, ++A, ++B)
	{
		volatile double v1 = *A * *B;
		re += v1;
	}

	return re;
}

TARGET512 double SumProd512(float* A, float* B, int64 n)
{
	return SumProdAVX(A, B, n);
	/*
	constexpr int N = 4;
	int64 i = 0;
	double re = 0;

	if (n >= N * sizeof(__m256) / sizeof(float))
	{
		__m512d s[N];
		REP(N) s[kk] = _mm512_setzero_pd();

		for (int64 l1 = n - N * sizeof(__m256) / sizeof(float); i <= l1; i += N * sizeof(__m256) / sizeof(float))
		{
			REP(N) 
			{
				s[kk] = _mm512_add_pd(s[kk], _mm512_mul_pd(_mm512_cvtps_pd(_mm256_loadu_ps(A)), _mm512_cvtps_pd(_mm256_loadu_ps(B))));
				A += sizeof(__m256) / sizeof(float); 
				B += sizeof(__m256) / sizeof(float);
			}
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm512_add_pd(s[kk], s[kk + K]);

		re = __mm512_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, ++A, ++B)
	{
		volatile double v1 = (double)*A * (double)*B;
		re += v1;
	}

	return re;
	*/
}

TARGET512 float SumProd512x(float* A, float* B, int64 n)
{
	constexpr int N = 2;
	int64 i = 0;
	volatile float re = 0;

	if (n >= N * sizeof(__m512) / sizeof(float))
	{
		__m512 s[N];
		REP(N) s[kk] = _mm512_setzero_ps();

		for (int64 l1 = n - N * sizeof(__m512) / sizeof(float); i <= l1; i += N * sizeof(__m512) / sizeof(float))
		{
			REP(N)
			{
				s[kk] = _mm512_add_ps(s[kk], _mm512_mul_ps(_mm512_loadu_ps(A), _mm512_loadu_ps(B)));
				A += sizeof(__m512) / sizeof(float);
				B += sizeof(__m512) / sizeof(float);
			}
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm512_add_ps(s[kk], s[kk + K]);

		re = __mm512_reduce_add_ps(s[0]);
	}

	for (; i < n; ++i, ++A, ++B)
	{
		volatile float v1 = *A * *B;
		re += v1;
	}

	return re;
}

TARGET512 void Add512(double* A, double* B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m512d) / sizeof(double))
	{
		__m512d a[N], b[N];

		for (int64 l1 = n - N * sizeof(__m512d) / sizeof(double); i <= l1; i += N * sizeof(__m512d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm512_loadu_pd(A); A += 8; }

			REP(N) { b[kk] = _mm512_loadu_pd(B); B += 8; }

			REP(N) a[kk] = _mm512_add_pd(a[kk], b[kk]);

			REP(N) _mm512_storeu_pd(A + (kk - N) * sizeof(__m512d) / sizeof(double), a[kk]);
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGET512 void Add512(float* A, float* B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m512) / sizeof(float))
	{
		__m512 a[N], b[N];

		for (int64 l1 = n - N * sizeof(__m512) / sizeof(float); i <= l1; i += N * sizeof(__m512) / sizeof(float))
		{
			REP(N) { a[kk] = _mm512_loadu_ps(A); A += sizeof(__m512) / sizeof(float); }

			REP(N) { b[kk] = _mm512_loadu_ps(B); B += sizeof(__m512) / sizeof(float); }

			REP(N) a[kk] = _mm512_add_ps(a[kk], b[kk]);

			REP(N) _mm512_storeu_ps(A + (kk - N) * sizeof(__m512) / sizeof(float), a[kk]);
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGET512 void Add512(int64* A, int64* B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m512i) / sizeof(int64))
	{
		__m512i a[N], b[N];

		for (int64 l1 = n - N * sizeof(__m512i) / sizeof(int64); i <= l1; i += N * sizeof(__m512i) / sizeof(int64))
		{
			REP(N) { a[kk] = _mm512_loadu_si512((__m512i*)A); A += sizeof(__m512i) / sizeof(int64); }

			REP(N) { b[kk] = _mm512_loadu_si512((__m512i*)B); B += sizeof(__m512i) / sizeof(int64); }

			REP(N) a[kk] = _mm512_add_epi64(a[kk], b[kk]);

			REP(N) _mm512_storeu_si512((__m512i*)(A + (kk - N) * sizeof(__m512i) / sizeof(int64)), a[kk]);
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGET512 void Add512(int* A, int* B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m512i) / sizeof(int))
	{
		__m512i a[N], b[N];

		for (int64 l1 = n - N * sizeof(__m512i) / sizeof(int); i <= l1; i += N * sizeof(__m512i) / sizeof(int))
		{
			REP(N) { a[kk] = _mm512_loadu_si512((__m512i*)A); A += sizeof(__m512i) / sizeof(int); }

			REP(N) { b[kk] = _mm512_loadu_si512((__m512i*)B); B += sizeof(__m512i) / sizeof(int); }

			REP(N) a[kk] = _mm512_add_epi32(a[kk], b[kk]);

			REP(N) _mm512_storeu_si512((__m512i*)(A + (kk - N) * sizeof(__m512i) / sizeof(int)), a[kk]);
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGET512 void Add512(double* A, double B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m512d) / sizeof(double))
	{
		__m512d b = _mm512_set1_pd(B), a[N];

		for (int64 l1 = n - N * sizeof(__m512d) / sizeof(double); i <= l1; i += N * sizeof(__m512d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm512_loadu_pd(A); A += sizeof(__m512d) / sizeof(double); }

			REP(N) a[kk] = _mm512_add_pd(a[kk], b);

			REP(N) _mm512_storeu_pd(A + (kk - N) * sizeof(__m512d) / sizeof(double), a[kk]);
		}
	}

	for (; i < n; ++i, A++)
		*A += B;
}

TARGET512 void Add512(float* A, float B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m512) / sizeof(float))
	{
		__m512 b = _mm512_set1_ps(B), a[N];

		for (int64 l1 = n - N * sizeof(__m512) / sizeof(float); i <= l1; i += N * sizeof(__m512) / sizeof(float))
		{
			REP(N) { a[kk] = _mm512_loadu_ps(A); A += sizeof(__m512) / sizeof(float); }

			REP(N) a[kk] = _mm512_add_ps(a[kk], b);

			REP(N) _mm512_storeu_ps(A + (kk - N) * sizeof(__m512) / sizeof(float), a[kk]);
		}
	}

	for (; i < n; ++i, A++, B)
		*A += B;
}

TARGET512 void Mul512(double* C, double* A, double* B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m512d) / sizeof(double))
	{
		__m512d a[N], b[N];

		for (int64 l1 = n - N * sizeof(__m512d) / sizeof(double); i <= l1; i += N * sizeof(__m512d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm512_loadu_pd(A); A += sizeof(__m512d) / sizeof(double); }

			REP(N) { b[kk] = _mm512_loadu_pd(B); B += sizeof(__m512d) / sizeof(double); }

			REP(N) a[kk] = _mm512_mul_pd(a[kk], b[kk]);

			REP(N) { _mm512_storeu_pd(C, a[kk]); C += sizeof(__m512d) / sizeof(double); }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *A++ * *B++;
		*C++ = v1;
	}
}

TARGET512 void Mul512(float* C, float* A, float* B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m512) / sizeof(float))
	{
		__m512 a[N], b[N];

		for (int64 l1 = n - N * sizeof(__m512) / sizeof(float); i <= l1; i += N * sizeof(__m512) / sizeof(float))
		{
			REP(N) { a[kk] = _mm512_loadu_ps(A); A += sizeof(__m512) / sizeof(float); }

			REP(N) { b[kk] = _mm512_loadu_ps(B); B += sizeof(__m512) / sizeof(float); }

			REP(N) a[kk] = _mm512_mul_ps(a[kk], b[kk]);

			REP(N) { _mm512_storeu_ps(C, a[kk]); C += sizeof(__m512) / sizeof(float); }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A++ * *B++;
		*C++ = v1;
	}
}

TARGET512 void Mul512(double* C, double* A, double B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m512d) / sizeof(double))
	{
		__m512d b = _mm512_set1_pd(B), a[N];

		for (int64 l1 = n - N * sizeof(__m512d) / sizeof(double); i <= l1; i += N * sizeof(__m512d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm512_loadu_pd(A); A += sizeof(__m512d) / sizeof(double); }

			REP(N) a[kk] = _mm512_mul_pd(a[kk], b);

			REP(N) { _mm512_storeu_pd(C, a[kk]); C += sizeof(__m512d) / sizeof(double); }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *A++ * B;
		*C++ = v1;
	}
}

TARGET512 void Mul512(float* C, float* A, float B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m512) / sizeof(float))
	{
		__m512 b = _mm512_set1_ps(B), a[N];

		for (int64 l1 = n - N * sizeof(__m512) / sizeof(float); i <= l1; i += N * sizeof(__m512) / sizeof(float))
		{
			REP(N) { a[kk] = _mm512_loadu_ps(A); A += sizeof(__m512) / sizeof(float); }

			REP(N) a[kk] = _mm512_mul_ps(a[kk], b);

			REP(N) { _mm512_storeu_ps(C, a[kk]); C += sizeof(__m512) / sizeof(float); }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A++ * B;
		*C++ = v1;
	}
}

TARGET512 void Mul512(double* A, double B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m512d) / sizeof(double))
	{
		__m512d b = _mm512_set1_pd(B), a[N];

		for (int64 l1 = n - N * sizeof(__m512d) / sizeof(double); i <= l1; i += N * sizeof(__m512d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm512_loadu_pd(A); A += sizeof(__m512d) / sizeof(double); }

			REP(N) a[kk] = _mm512_mul_pd(a[kk], b);

			REP(N) _mm512_storeu_pd(A + (kk - N) * sizeof(__m512d) / sizeof(double), a[kk]);
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *A * B;
		*A++ = v1;
	}
}

TARGET512 void Mul512(float* A, float B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m512) / sizeof(float))
	{
		__m512 b = _mm512_set1_ps(B), a[N];

		for (int64 l1 = n - N * sizeof(__m512) / sizeof(float); i <= l1; i += N * sizeof(__m512) / sizeof(float))
		{
			REP(N) { a[kk] = _mm512_loadu_ps(A); A += sizeof(__m512) / sizeof(float); }

			REP(N) a[kk] = _mm512_mul_ps(a[kk], b);

			REP(N) _mm512_storeu_ps(A + (kk - N) * sizeof(__m512) / sizeof(float), a[kk]);
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A * B;
		*A++ = v1;
	}
}

TARGET512 void AddProd512(double* C, double* A, double* B, int64 n)
{
	constexpr int N = 1;
	int64 i = 0;

	if (n >= N * sizeof(__m512d) / sizeof(double))
	{
		__m512d a[N];

		for (int64 l1 = n - N * sizeof(__m512d) / sizeof(double); i <= l1; i += N * sizeof(__m512d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm512_loadu_pd(A); A += sizeof(__m512d) / sizeof(double); }

			REP(N) { a[kk] = _mm512_mul_pd(a[kk], _mm512_loadu_pd(B)); B += sizeof(__m512d) / sizeof(double); }

			double* C2 = C;

			REP(N) { a[kk] = _mm512_add_pd(a[kk], _mm512_loadu_pd(C)); C += sizeof(__m512d) / sizeof(double); }

			REP(N) { _mm512_storeu_pd(C2, a[kk]); C2 += sizeof(__m512d) / sizeof(double); }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *A++ * *B++;
		*C++ += v1;
	}
}

TARGET512 void AddProd512(float* C, float* A, float* B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m512) / sizeof(float))
	{
		__m512 a[N];
		ushort maskff = 0 - (n > -1);

		for (int64 l1 = n - N * sizeof(__m512) / sizeof(float); i <= l1; i += N * sizeof(__m512) / sizeof(float))
		{
			REP(N) { a[kk] = _mm512_loadu_ps(A); A += sizeof(__m512) / sizeof(float); }

			REP(N) { a[kk] = _mm512_mul_ps(a[kk], _mm512_loadu_ps(B)); B += sizeof(__m512) / sizeof(float); }

			float* C2 = C;

			REP(N) { a[kk] = _mm512_mask_add_ps(a[kk], maskff, a[kk], _mm512_loadu_ps(C)); C += sizeof(__m512) / sizeof(float); }

			REP(N) { _mm512_storeu_ps(C2, a[kk]); C2 += sizeof(__m512) / sizeof(float); }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A++ * *B++;
		*C++ += v1;
	}
}

TARGET512 void AddProd512(double* C, double* A, double B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m512d) / sizeof(double))
	{
		__m512d a[N], b = _mm512_set1_pd(B);

		for (int64 l1 = n - N * sizeof(__m512d) / sizeof(double); i <= l1; i += N * sizeof(__m512d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm512_loadu_pd(A); A += sizeof(__m512d) / sizeof(double); }

			REP(N)  a[kk] = _mm512_mul_pd(a[kk], b);

			double* C2 = C;

			REP(N) { a[kk] = _mm512_add_pd(a[kk], _mm512_loadu_pd(C)); C += sizeof(__m512d) / sizeof(double); }

			REP(N) { _mm512_storeu_pd(C2, a[kk]); C2 += sizeof(__m512d) / sizeof(double); }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *A++ * B;
		*C++ += v1;
	}
}

TARGET512 void AddProd512(double* C, float* A, double B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m512d) / sizeof(double))
	{
		__m512d a[N], b = _mm512_set1_pd(B);

		for (int64 l1 = n - N * sizeof(__m512d) / sizeof(double); i <= l1; i += N * sizeof(__m512d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm512_cvtps_pd(_mm256_loadu_ps(A)); A += sizeof(__m256) / sizeof(float); }

			REP(N)  a[kk] = _mm512_mul_pd(a[kk], b);

			double* C2 = C;

			REP(N) { a[kk] = _mm512_add_pd(a[kk], _mm512_loadu_pd(C)); C += sizeof(__m512d) / sizeof(double); }

			REP(N) { _mm512_storeu_pd(C2, a[kk]); C2 += sizeof(__m512d) / sizeof(double); }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *A++ * B;
		*C++ += v1;
	}
}

TARGET512 void AddProd512(float* C, float* A, float B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;
	ushort maskff = 0 - (n > -1);

	if (n >= N * sizeof(__m512) / sizeof(float))
	{
		__m512 a[N], b = _mm512_set1_ps(B);

		for (int64 l1 = n - N * sizeof(__m512) / sizeof(float); i <= l1; i += N * sizeof(__m512) / sizeof(float))
		{
			REP(N) { a[kk] = _mm512_loadu_ps(A); A += sizeof(__m512) / sizeof(float); }

			REP(N)  a[kk] = _mm512_mul_ps(a[kk], b); 

			float* C2 = C;

			REP(N) { a[kk] = _mm512_mask_add_ps(a[kk], maskff, a[kk], _mm512_loadu_ps(C)); C += sizeof(__m512) / sizeof(float); }

			REP(N) { _mm512_storeu_ps(C2, a[kk]); C2 += sizeof(__m512) / sizeof(float); }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A++ * B;
		*C++ += v1;
	}
}

TARGET512 void Unify512(double* A, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;
	double invsum = 1.0 / (Sum512(A, n) + n * MIN_FREQ);

	if (n >= N * sizeof(__m512) / sizeof(double))
	{
		__m512d a[N], minf = _mm512_set1_pd(MIN_FREQ), invs = _mm512_set1_pd(invsum);

		for (int64 l1 = n - N * sizeof(__m512) / sizeof(double); i <= l1; i += N * sizeof(__m512) / sizeof(double))
		{
			double* A2 = A;

			REP(N) { a[kk] = _mm512_loadu_pd(A); A += sizeof(__m512) / sizeof(double); }

			REP(N) a[kk] = _mm512_add_pd(a[kk], minf);

			REP(N) a[kk] = _mm512_mul_pd(a[kk], invs);

			REP(N) { _mm512_storeu_pd(A2, a[kk]); A2 += sizeof(__m512) / sizeof(double); }
		}
	}

	for (; i < n; ++i, ++A)
	{
		volatile double v1 = (double)*A + MIN_FREQ;
		*A = v1 * invsum;
	}
}

TARGET512 void Unify512(float* A, int64 n)
{
	constexpr int N = 2;
	int64 i = 0;
	double invsum = 1.0 / (Sum512(A, n) + n * MIN_FREQ);

	if (n >= N * sizeof(__m256) / sizeof(float))
	{
		__m512d a[N], minf = _mm512_set1_pd(MIN_FREQ), invs = _mm512_set1_pd(invsum);

		for (int64 l1 = n - N * sizeof(__m256) / sizeof(float); i <= l1; i += N * sizeof(__m256) / sizeof(float))
		{
			float* A2 = A;

			REP(N) { a[kk] = _mm512_cvtps_pd(_mm256_loadu_ps(A)); A += sizeof(__m256) / sizeof(float); }

			REP(N) a[kk] = _mm512_add_pd(a[kk], minf);

			REP(N) a[kk] = _mm512_mul_pd(a[kk], invs);

			REP(N) { _mm256_storeu_ps(A2, _mm512_cvtpd_ps(a[kk])); A2 += sizeof(__m256) / sizeof(float); }
		}
	}

	for (; i < n; ++i, ++A)
	{
		volatile double v1 = (double)*A + MIN_FREQ;
		*A = v1 * invsum;
	}
}

TARGET512 char* StrNextIdx512(char* A, char val, int64 rep, int64 n)
{
	A++; n--;
	int64 i = 0;

	if (n >= 64)
	{
		__m512i v = _mm512_set1_epi8(val);

		for (int64 l1 = n - 64; i <= l1; i += 64, A += 64)
		{
			__m512i a = _mm512_loadu_si512(A);
			__mmask64 mask = _mm512_cmpeq_epu8_mask(a, v);
			int64 count = (int64)_mm_popcnt_u64(mask);

			if (rep > count)
			{
				rep -= count;
				continue;
			}
			else for (;; A++, mask >>= 1)
				if ((mask & 1) && !--rep)
					return A;
		}
	}

	for (; i < n; ++i, A++)
		if (*A == val && !--rep)
			return A;

	return NULL;
}

TARGET512 int64 CountChar512(char* A, char val, int64 n)
{
	uint64 re = 0;
	int64 i = 0;

	if (n >= 64)
	{
		__m512i v = _mm512_set1_epi8(val);

		for (int64 l1 = n - 64; i <= l1; i += 64)
		{
			re += _mm_popcnt_u64(_mm512_cmpeq_epu8_mask(_mm512_loadu_si512(A), v));
			A += 64;
		}
	}

	for (; i < n; ++i, A++)
		if (*A == val) re++;

	return (int64)re;
}

#endif
