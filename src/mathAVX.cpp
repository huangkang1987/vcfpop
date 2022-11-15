/* AVX Instruction Set Functions */

#include "vcfpop.h"

#ifndef __aarch64__

template TARGETAVX void RNGAVX<double>::Poly<32>(__m256d* arr, int n, __m256i* re);
template TARGETAVX void RNGAVX<double>::Poly<64>(__m256d* arr, int n, __m256i* re);
template TARGETAVX void RNGAVX<float >::Poly<32>(__m256*  arr, int n, __m256i* re);
template TARGETAVX void RNGAVX<float >::Poly<64>(__m256*  arr, int n, __m256i* re);

#ifndef _RNGAVX_FP64
/* Initialize rng */
TARGETAVX RNGAVX<double>::RNGAVX()
{

}

/* Initialize rng */
TARGETAVX RNGAVX<double>::RNGAVX(uint64 seed, uint64 salt)
{
	__m256i a[16], s, m;

	REP(16) { a[kk] = _mm256_set_epi64x(seed + 3, seed + 2, seed + 1, seed + 0); seed += 4; }

	s    = _mm256_set1_epi64x(salt);
	m    = _mm256_set1_epi32(0x5bd1e995);

	REP(16) a[kk] = _mm256_xor_si256(a[kk], _mm256_slli_epi64(_mm256_andnot_si256(a[kk], _mm256_set1_epi64x(0xFFFFFFFF)), 32));

	s = _mm256_xor_si256(s, _mm256_slli_epi64(_mm256_andnot_si256(s, _mm256_set1_epi64x(0xFFFFFFFF)), 32));

	// uint s = s ^ sizeof(uint);
	s = _mm256_xor_si256(s, _mm256_set1_epi32(sizeof(uint)));

	// a *= m;
	REP(16) a[kk] = _mm256_mullo_epi32(a[kk], m);//OK

	// a ^= a >> 24;
	REP(16) a[kk] = _mm256_xor_si256(a[kk], _mm256_srli_epi32(a[kk], 24));

	// a *= m;
	REP(16) a[kk] = _mm256_mullo_epi32(a[kk], m);//OK

	// s *= m;
	s = _mm256_mullo_epi32(s, m);//OK

	// a ^= s;
	REP(16) a[kk] = _mm256_xor_si256(a[kk], s);

	// a ^= a >> 13;
	REP(16) a[kk] = _mm256_xor_si256(a[kk], _mm256_srli_epi32(a[kk], 13));
	
	// a *= m;
	REP(16) a[kk] = _mm256_mullo_epi32(a[kk], m);//OK

	// a ^= a >> 15;
	REP(16) a[kk] = _mm256_xor_si256(a[kk], _mm256_srli_epi32(a[kk], 15));

	// original
	REP(16) x[kk] = _mm256_xor_si256(_mm256_set1_epi64x(0x159A55E5075BCD15), a[kk]);

	REP(16) a[kk] = _mm256_slli_epi64(a[kk], 6);

	REP(16) y[kk] = _mm256_xor_si256(_mm256_set1_epi64x(0x054913331F123BB5), a[kk]);
}

/* Draw a uniform distriubted real number */
template<int nbits>
TARGETAVX void RNGAVX<double>::Poly(__m256d* arr, int n, __m256i* re)
{
	__m256d t[16], s[16];
	__m256d one = _mm256_set1_pd(1.0);
	__m256i mask1 = _mm256_set1_epi64x(0x000FFFFFFFFFFFFF);
	__m256i mask2 = _mm256_set1_epi64x(0x3FF0000000000000);
	__m256i* r = (__m256i*)t;

	REP(16) s[kk] = _mm256_setzero_pd();

	for (int i16 = 0; i16 < n * 16; i16 += 16)
		REP(16) s[kk] = _mm256_add_pd(s[kk], arr[kk + i16]);

	{
		__m256i a[16], b[16];

		REP(16) a[kk] = x[kk];

		REP(16) b[kk] = y[kk];

		REP(16) x[kk] = b[kk];

		REP(16) a[kk] = _mm256_xor_si256(a[kk], _mm256_slli_epi64(a[kk], 23));

		REP(16) a[kk] = _mm256_xor_si256(a[kk], _mm256_srli_epi64(a[kk], 18));

		REP(16) a[kk] = _mm256_xor_si256(a[kk], b[kk]);

		REP(16) a[kk] = _mm256_xor_si256(a[kk], _mm256_srli_epi64(b[kk], 5));

		REP(16) y[kk] = a[kk];

		REP(16) r[kk] = _mm256_add_epi64(a[kk], b[kk]);
	}

	REP(16) r[kk] = _mm256_and_si256(r[kk], mask1);

	REP(16) r[kk] = _mm256_or_si256(r[kk], mask2);

	REP(16) t[kk] = _mm256_sub_pd(t[kk], one);

	REP(16) t[kk] = _mm256_mul_pd(t[kk], s[kk]);

	__m256i midx[16], nidx = _mm256_setzero_si256(), ninc = _mm256_set1_epi64x(1);
	__m256d f[16], b[16];

	REP(16) midx[kk] = _mm256_set1_epi64x(n - 1);

	REP(16) f[kk] = _mm256_setzero_pd();

	for (int i16 = 0; i16 < n * 16; i16 += 16)
	{
		REP(16) b[kk] = _mm256_cmp_pd(t[kk], arr[kk + i16], _CMP_LT_OS);

		t[ 0] = _mm256_sub_pd(t[ 0], arr[ 0 + i16]);
		t[ 1] = _mm256_sub_pd(t[ 1], arr[ 1 + i16]);
		t[ 2] = _mm256_sub_pd(t[ 2], arr[ 2 + i16]);
		t[ 3] = _mm256_sub_pd(t[ 3], arr[ 3 + i16]);

		REP(16) b[kk] = _mm256_andnot_pd(f[kk], b[kk]);

		t[ 4] = _mm256_sub_pd(t[ 4], arr[ 4 + i16]);
		t[ 5] = _mm256_sub_pd(t[ 5], arr[ 5 + i16]);
		t[ 6] = _mm256_sub_pd(t[ 6], arr[ 6 + i16]);
		t[ 7] = _mm256_sub_pd(t[ 7], arr[ 7 + i16]);

		REP(16) f[kk] = _mm256_or_pd(f[kk], b[kk]);
			
		t[ 8] = _mm256_sub_pd(t[ 8], arr[ 8 + i16]);
		t[ 9] = _mm256_sub_pd(t[ 9], arr[ 9 + i16]);
		t[11] = _mm256_sub_pd(t[11], arr[11 + i16]);
		t[10] = _mm256_sub_pd(t[10], arr[10 + i16]);

		REP(16) midx[kk] = _mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(midx[kk]), _mm256_castsi256_pd(nidx), b[kk]));//ok

		t[12] = _mm256_sub_pd(t[12], arr[12 + i16]);
		t[13] = _mm256_sub_pd(t[13], arr[13 + i16]);
		t[14] = _mm256_sub_pd(t[14], arr[14 + i16]);
		t[15] = _mm256_sub_pd(t[15], arr[15 + i16]);

		nidx = _mm256_add_epi64(nidx, ninc);
	}

	if constexpr (nbits == 32)
	{
		__m128i* midx2 = (__m128i*)midx, *re2 = (__m128i*)re;
		REP(16) re2[kk] = _mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(midx2[0 + (kk << 1)]), _mm_castsi128_ps(midx2[1 + (kk << 1)]), _MM_SHUFFLE(2, 0, 2, 0)));
	}
	else
	{
		REP(16) re[kk] = midx[kk];
	}
}
#endif

#ifndef _RNGAVX_FP32
/* Initialize rng */
TARGETAVX RNGAVX<float>::RNGAVX()
{
}

/* Initialize rng */
TARGETAVX RNGAVX<float>::RNGAVX(uint64 seed, uint64 salt)
{
	__m256i a[8], s, m;
	REP(8) { a[kk] = _mm256_set_epi32(Mix(seed + 7), Mix(seed + 6), Mix(seed + 5), Mix(seed + 4), Mix(seed + 3), Mix(seed + 2), Mix(seed + 1), Mix(seed + 0)); seed += 8; }
		
	s = _mm256_set1_epi32(Mix(salt));
	m = _mm256_set1_epi32(0x5bd1e995);

	// uint s = s ^ sizeof(uint);
	s = _mm256_xor_si256(s, _mm256_set1_epi32(sizeof(uint)));

	// a *= m;
	REP(8) a[kk] = _mm256_mullo_epi32(a[kk], m);//OK

	// a ^= a >> 24;
	REP(8) a[kk] = _mm256_xor_si256(a[kk], _mm256_srli_epi32(a[kk], 24));

	// a *= m;
	REP(8) a[kk] = _mm256_mullo_epi32(a[kk], m);//OK

	// s *= m;
	s = _mm256_mullo_epi32(s, m);//OK

	// a ^= s;
	REP(8) a[kk] = _mm256_xor_si256(a[kk], s);

	// a ^= a >> 13;
	REP(8) a[kk] = _mm256_xor_si256(a[kk], _mm256_srli_epi32(a[kk], 13));

	// a *= m;
	REP(8) a[kk] = _mm256_mullo_epi32(a[kk], m);//OK

	// a ^= a >> 15;
	REP(8) a[kk] = _mm256_xor_si256(a[kk], _mm256_srli_epi32(a[kk], 15));

	// original
	REP(8) x[kk] = _mm256_xor_si256(_mm256_set1_epi32(0x075BCD15), a[kk]);

	REP(8) a[kk] = _mm256_slli_epi32(a[kk], 3);

	REP(8) y[kk] = _mm256_xor_si256(_mm256_set1_epi32(0x159A55E5), a[kk]);

	REP(8) a[kk] = _mm256_slli_epi32(a[kk], 3);

	REP(8) z[kk] = _mm256_xor_si256(_mm256_set1_epi32(0x1F123BB5), a[kk]);
}

/* Draw a uniform distriubted real number */
template<int nbits>
TARGETAVX void RNGAVX<float>::Poly(__m256* arr, int n, __m256i* re)
{
	__m256 t[8], s[8]; __m256i u[8];
	__m256 one = _mm256_set1_ps(1.0f);
	__m256i mask1 = _mm256_set1_epi32(0x007FFFFF);
	__m256i mask2 = _mm256_set1_epi32(0x3F800000);
	__m256i* r = (__m256i*)t;

	REP(8) s[kk] = _mm256_setzero_ps();

	for (int i8 = 0; i8 < n * 8; i8 += 8)
		REP(8) s[kk] = _mm256_add_ps(s[kk], arr[kk + i8]);

	{
		//xorshift
		REP(8) u[kk] = _mm256_slli_epi32(x[kk], 16);
		REP(8) x[kk] = _mm256_xor_si256(x[kk], u[kk]);

		REP(8) u[kk] = _mm256_srli_epi32(x[kk], 5);
		REP(8) x[kk] = _mm256_xor_si256(x[kk], u[kk]);

		REP(8) u[kk] = _mm256_slli_epi32(x[kk], 1);
		REP(8) x[kk] = _mm256_xor_si256(x[kk], u[kk]);

		REP(8) u[kk] = x[kk];

		REP(8) x[kk] = y[kk];

		REP(8) y[kk] = z[kk];

		REP(8) z[kk] = _mm256_xor_si256(u[kk], x[kk]);

		REP(8) z[kk] = _mm256_xor_si256(z[kk], y[kk]);
	}

	REP(8) r[kk] = _mm256_and_si256(z[kk], mask1);

	REP(8) r[kk] = _mm256_or_si256(r[kk], mask2);

	REP(8) t[kk] = _mm256_sub_ps(t[kk], one);

	REP(8) t[kk] = _mm256_mul_ps(t[kk], s[kk]);

	__m256i midx[8], nidx = _mm256_setzero_si256(), ninc = _mm256_set1_epi32(1);
	__m256 f[8], b[8];
	REP(8) midx[kk] = _mm256_set1_epi32(n - 1);
	REP(8) f[kk] = _mm256_setzero_ps();

	for (int i8 = 0; i8 < n * 8; i8 += 8)
	{
		REP(8) b[kk] = _mm256_cmp_ps(t[kk], arr[kk + i8], _CMP_LT_OS);

		t[0] = _mm256_sub_ps(t[0], arr[0 + i8]);
		t[1] = _mm256_sub_ps(t[1], arr[1 + i8]);

		REP(8) b[kk] = _mm256_andnot_ps(f[kk], b[kk]);

		t[2] = _mm256_sub_ps(t[2], arr[2 + i8]);
		t[3] = _mm256_sub_ps(t[3], arr[3 + i8]);

		REP(8) f[kk] = _mm256_or_ps(f[kk], b[kk]);

		t[4] = _mm256_sub_ps(t[4], arr[4 + i8]);
		t[5] = _mm256_sub_ps(t[5], arr[5 + i8]);

		REP(8) midx[kk] = _mm256_castps_si256(_mm256_blendv_ps(_mm256_castsi256_ps(midx[kk]), _mm256_castsi256_ps(nidx), b[kk]));//ok

		t[6] = _mm256_sub_ps(t[6], arr[6 + i8]);
		t[7] = _mm256_sub_ps(t[7], arr[7 + i8]);

		nidx = _mm256_add_epi32(nidx, ninc);
	}

	if constexpr (nbits == 32)
	{
		REP(8) re[kk] = midx[kk];
	}
	else
	{
		__m128i* midx2 = (__m128i*)midx;
		REP(16) re[kk] = _mm256_cvtepi32_epi64(midx2[kk]);
	}
}
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
	__m128 v2 = _mm_add_ps(_mm256_extractf128_ps(v1, 0), _mm256_extractf128_ps(v1, 1));
	__m128 v3 = _mm_add_ps(v2, _mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(v2), 8)));
	return simd_f32(v3, 0) + simd_f32(v3, 1);
}

__forceinline TARGETAVX float _mm256_reduce_mul_ps(__m256 v1)
{
	__m128 v2 = _mm_mul_ps(_mm256_extractf128_ps(v1, 0), _mm256_extractf128_ps(v1, 1));
	__m128 v3 = _mm_mul_ps(v2, _mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(v2), 8)));
	return simd_f32(v3, 0) * simd_f32(v3, 1);
}

__forceinline TARGETAVX double _mm256_reduce_add_psd(__m256 v1)
{
	__m256d v1b = _mm256_add_pd(
		_mm256_cvtps_pd(_mm256_extractf128_ps(v1, 0)),
		_mm256_cvtps_pd(_mm256_extractf128_ps(v1, 1)));
	__m128d v2 = _mm_add_pd(_mm256_extractf128_pd(v1b, 0), _mm256_extractf128_pd(v1b, 1));
	return simd_f64(v2, 0) + simd_f64(v2, 1);
}

__forceinline TARGETAVX double _mm256_reduce_mul_psd(__m256 v1)
{
	__m256d v1b = _mm256_mul_pd(
		_mm256_cvtps_pd(_mm256_extractf128_ps(v1, 0)),
		_mm256_cvtps_pd(_mm256_extractf128_ps(v1, 1)));
	__m128d v2 = _mm_mul_pd(_mm256_extractf128_pd(v1b, 0), _mm256_extractf128_pd(v1b, 1));
	return simd_f64(v2, 0) * simd_f64(v2, 1);
}

TARGETAVX int64 GetMinIdxAVX(double* A, int64 n, double& val)
{
	constexpr int N = 4;
	int64 i = 0;
	val = DBL_MAX;
	uint64 idx = (uint64)-1;

	if (n >= N * sizeof(__m256d) / sizeof(double))
	{
		__m256d min1[N], f[N], a[N];
		__m256i midx[N], nidx[N], msep = _mm256_set1_epi64x(N * sizeof(__m256d) / sizeof(double));
		REP(N) min1[kk] = _mm256_set1_pd(val);
		REP(N) midx[kk] = _mm256_set1_epi64x(0xFFFFFFFFFFFFFFFF);
		REP(N) nidx[kk] = _mm256_set_epi64x(3 + (kk << 2), 2 + (kk << 2), 1 + (kk << 2), 0 + (kk << 2));

		for (int64 l1 = n - N * sizeof(__m256d) / sizeof(double); i <= l1; i += N * sizeof(__m256d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm256_loadu_pd(A); A += sizeof(__m256d) / sizeof(double); }

			REP(N) f[kk] = _mm256_cmp_pd(min1[kk], a[kk], _CMP_GT_OS);

			REP(N) min1[kk] = _mm256_min_pd(min1[kk], a[kk]);

			REP(N) midx[kk] = _mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(midx[kk]), _mm256_castsi256_pd(nidx[kk]), f[kk]));//ok

			REP(N) nidx[kk] = _mm256_add_epi64(nidx[kk], msep);
		}

		for (int K = sizeof(min1) / sizeof(min1[0]) / 2; K >= 1; K >>= 1)
		{
			REP(K) f[kk] = _mm256_cmp_pd(min1[kk], min1[kk + K], _CMP_GT_OS);
			REP(K) min1[kk] = _mm256_min_pd(min1[kk], min1[kk + K]);
			REP(K) midx[kk] = _mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(midx[kk]), _mm256_castsi256_pd(midx[kk + K]), f[kk]));
		}

		for (int64 j = 0; j < sizeof(__m256d) / sizeof(double); ++j)
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

TARGETAVX int64 GetMinIdxAVX(float* A, int64 n, float& val)
{
	constexpr int N = 4;
	int64 i = 0;
	val = FLT_MAX;
	uint idx = (uint)-1;

	if (n >= N * sizeof(__m256) / sizeof(float))
	{
		__m256 min1[N], f[N], a[N];
		__m256i midx[N], nidx[N], msep = _mm256_set1_epi32(N * sizeof(__m256) / sizeof(float));
		REP(N) min1[kk] = _mm256_set1_ps(val);
		REP(N) midx[kk] = _mm256_set1_epi8((char)0xFF);
		REP(N) nidx[kk] = _mm256_set_epi32(7 + (kk << 3), 6 + (kk << 3), 5 + (kk << 3), 4 + (kk << 3), 3 + (kk << 3), 2 + (kk << 3), 1 + (kk << 3), 0 + (kk << 3));

		for (int64 l1 = n - N * sizeof(__m256) / sizeof(float); i <= l1; i += N * sizeof(__m256) / sizeof(float))
		{
			REP(N) { a[kk] = _mm256_loadu_ps(A); A += sizeof(__m256) / sizeof(float); }

			REP(N) f[kk] = _mm256_cmp_ps(min1[kk], a[kk], _CMP_GT_OS);

			REP(N) min1[kk] = _mm256_min_ps(min1[kk], a[kk]);

			REP(N) midx[kk] = _mm256_castps_si256(_mm256_blendv_ps(_mm256_castsi256_ps(midx[kk]), _mm256_castsi256_ps(nidx[kk]), f[kk]));//ok

			REP(N) nidx[kk] = _mm256_add_epi32(nidx[kk], msep);
		}

		for (int K = sizeof(min1) / sizeof(min1[0]) / 2; K >= 1; K >>= 1)
		{
			REP(K) f[kk] = _mm256_cmp_ps(min1[kk], min1[kk + K], _CMP_GT_OS);
			REP(K) min1[kk] = _mm256_min_ps(min1[kk], min1[kk + K]);
			REP(K) midx[kk] = _mm256_castps_si256(_mm256_blendv_ps(_mm256_castsi256_ps(midx[kk]), _mm256_castsi256_ps(midx[kk + K]), f[kk]));
		}

		for (int64 j = 0; j < sizeof(__m256) / sizeof(float); ++j)
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

TARGETAVX void GetMinMaxValAVX(double* A, int64 n, double& minv, double& maxv)
{
	constexpr int N = 4;
	int64 i = 0;
	minv = DBL_MAX;
	maxv = -DBL_MAX;

	if (n >= N * sizeof(__m256d) / sizeof(double))
	{
		__m256d min1[N], max1[N], a[N];
		REP(N) min1[kk] = _mm256_set1_pd(minv);
		REP(N) max1[kk] = _mm256_set1_pd(maxv);

		for (int64 l1 = n - N * sizeof(__m256d) / sizeof(double); i <= l1; i += N * sizeof(__m256d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm256_loadu_pd(A); A += sizeof(__m256d) / sizeof(double); }

			REP(N)
			{
				min1[kk] = _mm256_min_pd(min1[kk], a[kk]);
				max1[kk] = _mm256_max_pd(max1[kk], a[kk]);
			}
		}

		for (int K = sizeof(max1) / sizeof(max1[0]) / 2; K >= 1; K >>= 1)
		{
			REP(K) min1[kk] = _mm256_min_pd(min1[kk], min1[kk + K]);
			REP(K) max1[kk] = _mm256_max_pd(max1[kk], max1[kk + K]);
		}

		minv = Min(Min(simp_f64(min1, 0), simp_f64(min1, 1)),
				   Min(simp_f64(min1, 2), simp_f64(min1, 3)));
		maxv = Max(Max(simp_f64(max1, 0), simp_f64(max1, 1)),
				   Max(simp_f64(max1, 2), simp_f64(max1, 3)));
	}

	for (; i < n; ++i, ++A)
	{
		if (*A < minv) minv = *A;
		if (*A > maxv) maxv = *A;
	}
}

TARGETAVX void GetMinMaxValAVX(float* A, int64 n, float& minv, float& maxv)
{
	constexpr int N = 4;
	int64 i = 0;
	minv = FLT_MAX;
	maxv = -FLT_MAX;

	if (n >= N * sizeof(__m256) / sizeof(float))
	{
		__m256 min1[N], max1[N], a[N];
		REP(N) min1[kk] = _mm256_set1_ps(minv);
		REP(N) max1[kk] = _mm256_set1_ps(maxv);

		for (int64 l1 = n - N * sizeof(__m256) / sizeof(float); i <= l1; i += N * sizeof(__m256) / sizeof(float))
		{
			REP(N) { a[kk] = _mm256_loadu_ps(A); A += sizeof(__m256) / sizeof(float); }

			REP(N)
			{
				min1[kk] = _mm256_min_ps(min1[kk], a[kk]);
				max1[kk] = _mm256_max_ps(max1[kk], a[kk]);
			}
		}

		for (int K = sizeof(max1) / sizeof(max1[0]) / 2; K >= 1; K >>= 1)
		{
			REP(K) min1[kk] = _mm256_min_ps(min1[kk], min1[kk + K]);
			REP(K) max1[kk] = _mm256_max_ps(max1[kk], max1[kk + K]);
		}

		__m128* min2 = (__m128*)min1, *max2 = (__m128*)max1;
		min2[0] = _mm_min_ps(min2[0], min2[1]);
		max2[0] = _mm_max_ps(max2[0], max2[1]);

		minv = Min(Min(simp_f32(min1, 0), simp_f32(min1, 1)),
				   Min(simp_f32(min1, 2), simp_f32(min1, 3)));
		maxv = Max(Max(simp_f32(max1, 0), simp_f32(max1, 1)),
				   Max(simp_f32(max1, 2), simp_f32(max1, 3)));
	}

	for (; i < n; ++i, ++A)
	{
		if (*A < minv) minv = *A;
		if (*A > maxv) maxv = *A;
	}
}

TARGETAVX double GetMaxValAVX(double* A, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;
	double val = -DBL_MAX;

	if (n >= N * sizeof(__m256d) / sizeof(double))
	{
		__m256d max1[N];
		REP(N) max1[kk] = _mm256_set1_pd(val);

		for (int64 l1 = n - N * sizeof(__m256d) / sizeof(double); i <= l1; i += N * sizeof(__m256d) / sizeof(double))
		{
			REP(N)
			{
				max1[kk] = _mm256_max_pd(max1[kk], _mm256_loadu_pd(A));
				A += sizeof(__m256d) / sizeof(double);
			}
		}

		for (int K = sizeof(max1) / sizeof(max1[0]) / 2; K >= 1; K >>= 1)
			REP(K) max1[kk] = _mm256_max_pd(max1[kk], max1[kk + K]);

		val = Max(Max(simp_f64(max1, 0), simp_f64(max1, 1)),
				  Max(simp_f64(max1, 2), simp_f64(max1, 3)));
	}

	for (; i < n; ++i, ++A)
	{
		if (*A < val) continue;
		val = *A;
	}

	return val;
}

TARGETAVX float GetMaxValAVX(float* A, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;
	float val = -FLT_MAX;

	if (n >= N * sizeof(__m256) / sizeof(float))
	{
		__m256 max1[N], a[N];
		REP(N) max1[kk] = _mm256_set1_ps(val);

		for (int64 l1 = n - N * sizeof(__m256) / sizeof(float); i <= l1; i += N * sizeof(__m256) / sizeof(float))
		{
			REP(N) { a[kk] = _mm256_loadu_ps(A); A += sizeof(__m256) / sizeof(float); }

			REP(N) { max1[kk] = _mm256_max_ps(max1[kk], a[kk]); }
		}

		for (int K = sizeof(max1) / sizeof(max1[0]) / 2; K >= 1; K >>= 1)
			REP(K) max1[kk] = _mm256_max_ps(max1[kk], max1[kk + K]);

		__m128* max2 = (__m128*)max1;
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

TARGETAVX double GetMaxValAVX(double* A, int64 n, int64 sep)
{
	int64 i = 0;
	double val = -DBL_MAX;

	if (n >= 8)
	{
		REP(8) { _mm_prefetch((const char*)A, _MM_HINT_T0); A += sep; }

		__m256i vindex = _mm256_set_epi64x(-9 * sep, -10 * sep, -11 * sep, -12 * sep);
		__m256d max1[2]; 
		REP(2) max1[kk] = _mm256_set1_pd(val);

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			REP(4) { _mm_prefetch((const char*)A, _MM_HINT_T0); A += sep; }

			max1[0] = _mm256_max_pd(max1[0], _mm256_i64gather_pd(A, vindex, sizeof(double)));

			REP(4) { _mm_prefetch((const char*)A, _MM_HINT_T0); A += sep; }

			max1[1] = _mm256_max_pd(max1[1], _mm256_i64gather_pd(A, vindex, sizeof(double)));
		}

		max1[0] = _mm256_max_pd(max1[0], max1[1]);

		val = Max(Max(simp_f64(max1, 0), simp_f64(max1, 1)),
				  Max(simp_f64(max1, 2), simp_f64(max1, 3)));
		
		A -= 8 * sep;
	}

	for (; i < n; ++i, A += sep)
	{
		if (*A < val) continue;
		val = *A;
	}

	return val;
}

TARGETAVX float GetMaxValAVX(float* A, int64 n, int64 sep)
{
	constexpr int N = 4;
	int64 i = 0;
	float val = -FLT_MAX;

	if (n >= sizeof(__m128) / sizeof(float))
	{
		__m256i vindex = _mm256_set_epi64x(3 * sep, 2 * sep, 1 * sep, 0 * sep);
		__m128 max1[N];
		REP(N) max1[kk] = _mm_set1_ps(val);

		for (int64 l1 = N * sizeof(__m128) / sizeof(float); i <= l1; i += N * sizeof(__m128) / sizeof(float))
		{
			REP(N) 
			{	
				max1[kk] = _mm_max_ps(max1[kk], _mm256_i64gather_ps(A, vindex, sizeof(float)));
				A += sizeof(__m128) / sizeof(float) * sep; 
			}
		}

		for (int K = sizeof(max1) / sizeof(max1[0]) / 2; K >= 1; K >>= 1)
			REP(K) max1[kk] = _mm_max_ps(max1[kk], max1[kk + K]);

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

TARGETAVX double GetMinValAVX(double* A, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;
	double val = DBL_MAX;

	if (n >= N * sizeof(__m256d) / sizeof(double))
	{
		__m256d min1[N];
		REP(N) min1[kk] = _mm256_set1_pd(val);

		for (int64 l1 = n - N * sizeof(__m256d) / sizeof(double); i <= l1; i += N * sizeof(__m256d) / sizeof(double))
		{
			REP(N)
			{
				min1[kk] = _mm256_min_pd(min1[kk], _mm256_loadu_pd(A));
				A += sizeof(__m256d) / sizeof(double);
			}
		}

		for (int K = sizeof(min1) / sizeof(min1[0]) / 2; K >= 1; K >>= 1)
			REP(K) min1[kk] = _mm256_min_pd(min1[kk], min1[kk + K]);

		val = Min(Min(simp_f64(min1, 0), simp_f64(min1, 1)),
				  Min(simp_f64(min1, 2), simp_f64(min1, 3)));
	}

	for (; i < n; ++i, ++A)
	{
		if (*A > val) continue;
		val = *A;
	}

	return val;
}

TARGETAVX float GetMinValAVX(float* A, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;
	float val = FLT_MAX;

	if (n >= N * sizeof(__m256) / sizeof(float))
	{
		__m256 min1[N], a[N];
		REP(N) min1[kk] = _mm256_set1_ps(val);

		for (int64 l1 = n - N * sizeof(__m256) / sizeof(float); i <= l1; i += N * sizeof(__m256) / sizeof(float))
		{
			REP(N) { a[kk] = _mm256_loadu_ps(A); A += sizeof(__m256) / sizeof(float); }

			REP(N) { min1[kk] = _mm256_min_ps(min1[kk], a[kk]); }
		}

		for (int K = sizeof(min1) / sizeof(min1[0]) / 2; K >= 1; K >>= 1)
			REP(K) min1[kk] = _mm256_min_ps(min1[kk], min1[kk + K]);

		__m128* min2 = (__m128*)min1;
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

TARGETAVX int64 GetMinValAVX(int64* A, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;
	int64 val = 0x7FFFFFFFFFFFFFFF;

	if (n >= N * sizeof(__m256i) / sizeof(int64))
	{
		__m256i min1[N], a[N], f[N];
		REP(N) min1[kk] = _mm256_set1_epi64x(0x7FFFFFFFFFFFFFFF);

		for (int64 l1 = n - N * sizeof(__m256i) / sizeof(int64); i <= l1; i += N * sizeof(__m256i) / sizeof(int64))
		{
			REP(N) { a[kk] = _mm256_loadu_si256((__m256i*)A); A += sizeof(__m256i) / sizeof(int64); }
			REP(N) f[kk] = _mm256_cmpgt_epi64(min1[kk], a[kk]); 
			REP(N)
			{
				min1[kk] =
				_mm256_castpd_si256(_mm256_blendv_pd(
					_mm256_castsi256_pd(min1[kk]),
					_mm256_castsi256_pd(a[kk]),
					_mm256_castsi256_pd(f[kk])));
			}
		}

		for (int K = sizeof(min1) / sizeof(min1[0]) / 2; K >= 1; K >>= 1)
			REP(K) min1[kk] =
			_mm256_castpd_si256(_mm256_blendv_pd(
				_mm256_castsi256_pd(min1[kk]),
				_mm256_castsi256_pd(min1[kk + K]),
				_mm256_castsi256_pd(_mm256_cmpgt_epi64(min1[kk], min1[kk + K]))));

		val = Min(Min(simp_i64(min1, 0), simp_i64(min1, 1)),
				  Min(simp_i64(min1, 2), simp_i64(min1, 3)));
	}

	for (; i < n; ++i, ++A)
	{
		if (*A > val) continue;
		val = *A;
	}

	return val;
}

TARGETAVX void SetValAVX(uint* A, ushort* B, int64 n)
{
	int64 i = 0;

	if (n >= 32)
	{
		__m256i a1, a2, a3, a4;
		__m128i b1, b2, b3, b4;

		for (int64 l1 = n - 32; i <= l1; i += 32)
		{
			b1 = _mm_loadu_si128((__m128i*)B); B += 8;
			b2 = _mm_loadu_si128((__m128i*)B); B += 8;
			b3 = _mm_loadu_si128((__m128i*)B); B += 8;
			b4 = _mm_loadu_si128((__m128i*)B); B += 8;

			a1 = _mm256_cvtepu16_epi32(b1);
			a2 = _mm256_cvtepu16_epi32(b2);
			a3 = _mm256_cvtepu16_epi32(b3);
			a4 = _mm256_cvtepu16_epi32(b4);

			_mm256_storeu_si256((__m256i*)A, a1); A += 8;
			_mm256_storeu_si256((__m256i*)A, a2); A += 8;
			_mm256_storeu_si256((__m256i*)A, a3); A += 8;
			_mm256_storeu_si256((__m256i*)A, a4); A += 8;
		}
	}

	for (; i < n; ++i)
		*A++ = *B++;
}

TARGETAVX void AddExponentAVX(int64& slog, __m256d& val)
{
	__m256i& vv = *(__m256i*)&val;
	__m256i mask1 = _mm256_set1_epi64x(0x7FF0000000000000);
	__m256i mask2 = _mm256_set1_epi64x(0x800FFFFFFFFFFFFF);
	__m256i mask3 = _mm256_set1_epi64x(0x3FF0000000000000);
	__m256i subv = _mm256_set1_epi64x(1023);

	__m256i t = _mm256_sub_epi64(_mm256_srli_epi64(_mm256_and_si256(vv, mask1), 52), subv);

	slog += (int)simd_u64(t, 0) + (int)simd_u64(t, 1) + (int)simd_u64(t, 2) + (int)simd_u64(t, 3);

	vv = _mm256_or_si256(_mm256_and_si256(vv, mask2), mask3);
}

TARGETAVX void AddExponentAVX(int64& slog, __m256& val)
{
	__m256i& vv = *(__m256i*)&val;
	__m256i mask1 = _mm256_set1_epi32(0x7F800000);
	__m256i mask2 = _mm256_set1_epi32(0x807FFFFF);
	__m256i mask3 = _mm256_set1_epi32(0x3F800000);
	__m256i subv = _mm256_set1_epi32(127);

	__m256i t = _mm256_sub_epi32(_mm256_srli_epi32(_mm256_and_si256(vv, mask1), 23), subv);

	t = _mm256_hadd_epi32(t, t);
	t = _mm256_hadd_epi32(t, t);
	slog += (int)simd_u32(t, 0) + (int)simd_u32(t, 4);

	vv = _mm256_or_si256(_mm256_and_si256(vv, mask2), mask3);
}

TARGETAVX void ChargeLogAVX(int64& slog, double& prod, __m256d& val)
{
	AddExponentAVX(slog, val);

	prod *= _mm256_reduce_mul_pd(val); 

	if (prod < DOUBLE_UNDERFLOW || prod > DOUBLE_OVERFLOW) [[unlikely]]
		AddExponent(slog, prod);
}

TARGETAVX void ChargeLogAVX(int64& slog, double& prod, __m256& val)
{
	AddExponentAVX(slog, val);

	prod *= _mm256_reduce_mul_psd(val);

	if (prod < DOUBLE_UNDERFLOW || prod > DOUBLE_OVERFLOW) [[unlikely]]
	AddExponent(slog, prod);
}

TARGETAVX double LogProdAVX(double* A, int64 n)
{
	int64 i = 0;
	int64 slog = 0; double prod = 1;
	
	if (n >= 32)
	{
		__m256d pd[2];
		REP(2) pd[kk] = _mm256_set1_pd(1.0);

		for (int64 l1 = n - 32; i <= l1; i += 32)
		{
			REP(2)
			{
				pd[0] = _mm256_mul_pd(pd[0], _mm256_mul_pd(_mm256_loadu_pd(A + 0), _mm256_loadu_pd(A + 8)));
				pd[1] = _mm256_mul_pd(pd[1], _mm256_mul_pd(_mm256_loadu_pd(A + 4), _mm256_loadu_pd(A + 12)));
				A += 16;
			}

			/*
			__m256d flag = _mm256_or_pd(
				_mm256_or_pd(_mm256_cmp_pd(pd[0], dunder, _CMP_LT_OS), _mm256_cmp_pd(dover, pd[0], _CMP_LT_OS)),
				_mm256_or_pd(_mm256_cmp_pd(pd[1], dunder, _CMP_LT_OS), _mm256_cmp_pd(dover, pd[1], _CMP_LT_OS)));

			if (_mm256_testz_pd(flag, flag)) [[likely]] continue;
			*/

			REP(2) AddExponentAVX(slog, pd[kk]);
		}

		__m128d* pd1 = (__m128d*)pd;
		REP(4) ChargeLogSSE(slog, prod, pd1[kk]);
	}

	for (; i < n; ++i, ++A)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);

	return prod;
}

TARGETAVX double LogProdAVX(float* A, int64 n)
{
	int64 i = 0;
	int64 slog = 0; double prod = 1;
	
	if (n >= 64)
	{
		__m256d pd[2];
		REP(2) pd[kk] = _mm256_set1_pd(1.0);

		for (int64 l1 = n - 64; i <= l1; i += 64)
		{
			REP(4)
			{
				pd[0] = _mm256_mul_pd(pd[0], _mm256_mul_pd(
					_mm256_cvtps_pd(_mm_loadu_ps(A)),
					_mm256_cvtps_pd(_mm_loadu_ps(A + 8))));
				pd[1] = _mm256_mul_pd(pd[1], _mm256_mul_pd(
					_mm256_cvtps_pd(_mm_loadu_ps(A + 4)),
					_mm256_cvtps_pd(_mm_loadu_ps(A + 12))));
				A += 16;
			}

			REP(2) AddExponentAVX(slog, pd[kk]);
		}

		__m128d* pd1 = (__m128d*)pd;
		REP(4) ChargeLogSSE(slog, prod, pd1[kk]);
	}

	for (; i < n; ++i, ++A)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);

	return prod;
}

TARGETAVX float LogProdAVXx(float* A, int64 n)
{
	int64 i = 0;
	int64 slog = 0; double prod = 1;

	if (n >= 32)
	{
		__m256 pd[2];
		REP(2) pd[kk] = _mm256_set1_ps(1.0f);

		for (int64 l1 = n - 32; i <= l1; i += 32)
		{
			pd[0] = _mm256_mul_ps(pd[0], _mm256_loadu_ps(A)); 
			pd[1] = _mm256_mul_ps(pd[1], _mm256_loadu_ps(A + 8));
			pd[0] = _mm256_mul_ps(pd[0], _mm256_loadu_ps(A + 16));
			pd[1] = _mm256_mul_ps(pd[1], _mm256_loadu_ps(A + 24));

			A += 32;
			REP(2) AddExponentAVX(slog, pd[kk]);
		}

		__m128* pd1 = (__m128*)pd;
		REP(4) ChargeLogSSE(slog, prod, pd1[kk]);
	}

	for (; i < n; ++i, ++A)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);

	return prod;
}

TARGETAVX double LogProdAVX(double* A, int64 n, int64 sep)
{
	int64 i = 0;
	int64 slog = 0; double prod = 1;

	if (n >= 8)
	{
		REP(8) { _mm_prefetch((const char*)A, _MM_HINT_T0); A += sep; }

		__m256d pd[2], f1, f2;
		REP(2) pd[kk] = _mm256_set1_pd(1.0);
		__m256d dunder = _mm256_set1_pd(DOUBLE_UNDERFLOW), dover = _mm256_set1_pd(DOUBLE_OVERFLOW);
		__m256i vindex = _mm256_set_epi64x(-13 * sep, -14 * sep, -15 * sep, -16 * sep);

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

			pd[0] = _mm256_mul_pd(pd[0], _mm256_i64gather_pd(A, vindex, sizeof(double)));
			f1 = _mm256_cmp_pd(pd[0], dunder, _CMP_LT_OS);
			f2 = _mm256_cmp_pd(dover, pd[0], _CMP_LT_OS);

			pd[1] = _mm256_mul_pd(pd[1], _mm256_i64gather_pd(A + sep * 4, vindex, sizeof(double)));
			f1 = _mm256_or_pd(f1, _mm256_cmp_pd(pd[0], dunder, _CMP_LT_OS));
			f2 = _mm256_or_pd(f2, _mm256_cmp_pd(dover, pd[0], _CMP_LT_OS));

			f1 = _mm256_or_pd(f1, f2);

			if (_mm256_testz_pd(f1, f1)) [[likely]] continue;

			AddExponentAVX(slog, pd[0]);
			AddExponentAVX(slog, pd[1]);
		}

		ChargeLogAVX(slog, prod, pd[0]);
		ChargeLogAVX(slog, prod, pd[1]);

		A -= 8 * sep;
	}

	for (; i < n; ++i, A += sep)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);

	return prod;
}

TARGETAVX double LogProdAVX(float* A, int64 n, int64 sep)
{
	return LogProdSSE(A, n, sep);
}

TARGETAVX float LogProdAVXx(float* A, int64 n, int64 sep)
{
	return LogProdSSEx(A, n, sep);

	/*
	int64 i = 0;
	int64 slog = 0; double prod = 1;

	if (n >= 32)
	{
		__m256 pd[2];
		REP(2) pd[kk] = _mm256_set1_ps(1.0f);
		__m256i vindex = _mm256_set_epi32(7 * sep, 6 * sep, 5 * sep, 4 * sep, 3 * sep, 2 * sep, 1 * sep, 0 * sep);

		for (int64 l1 = n - 32; i <= l1; i += 32)
		{
			pd[0] = _mm256_mul_ps(pd[0], _mm256_i32xxgather_ps(A, vindex, sizeof(float))); A += 8 * sep;
			pd[1] = _mm256_mul_ps(pd[1], _mm256_i32xxgather_ps(A, vindex, sizeof(float))); A += 8 * sep;

			pd[0] = _mm256_mul_ps(pd[0], _mm256_i32xxgather_ps(A, vindex, sizeof(float))); A += 8 * sep;
			pd[1] = _mm256_mul_ps(pd[1], _mm256_i32xxgather_ps(A, vindex, sizeof(float))); A += 8 * sep;

			REP(2) AddExponentAVX(slog, pd[kk]);
		}

		__m128* pd1 = (__m128*)pd;
		REP(4) ChargeLogSSE(slog, prod, pd1[kk]);
	}

	for (; i < n; ++i, A += sep)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);

	return prod;
	*/
}

TARGETAVX double LogProdDivAVX(double* A, double* B, int64 n, int64 sep)
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
		__m256d pd1[2], pd2[2];
		REP(2) pd1[kk] = pd2[kk] = _mm256_set1_pd(1.0);
		__m256i vindex = _mm256_set_epi64x(3 * sep, 2 * sep, 1 * sep, 0 * sep);

		for (int64 l1 = n - 64; i <= l1; i += 64)
		{
			REP(8)
			{
				pd1[0] = _mm256_mul_pd(pd1[0], _mm256_i64gather_pd(A, vindex, sizeof(double))); A += 4 * sep;

				pd2[0] = _mm256_mul_pd(pd2[0], _mm256_i64gather_pd(B, vindex, sizeof(double))); B += 4 * sep;

				pd1[1] = _mm256_mul_pd(pd1[1], _mm256_i64gather_pd(A, vindex, sizeof(double))); A += 4 * sep;

				pd2[1] = _mm256_mul_pd(pd2[1], _mm256_i64gather_pd(B, vindex, sizeof(double))); B += 4 * sep;
			}

			REP(2) AddExponentAVX(slog1, pd1[kk]);
			REP(2) AddExponentAVX(slog2, pd2[kk]);
		}

		__m128d* pd11 = (__m128d*)pd1;
		__m128d* pd22 = (__m128d*)pd2;
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

TARGETAVX double LogProdDivAVX(float* A, float* B, int64 n, int64 sep)
{
	int64 i = 0;
	int64 slog1 = 0; double prod1 = 1;
	int64 slog2 = 0; double prod2 = 1;

	if (n >= 64)
	{
		__m256d pd1[2], pd2[2];
		REP(2) pd1[kk] = pd2[kk] = _mm256_set1_pd(1.0);

		for (int64 l1 = n - 64; i <= l1; i += 64)
		{
			REP(8)
			{
				pd1[0] = _mm256_mul_pd(pd1[0], _mm256_set_pd(A[3 * sep], A[2 * sep], A[1 * sep], A[0 * sep])); A += sep * 4;
				pd2[0] = _mm256_mul_pd(pd2[0], _mm256_set_pd(B[3 * sep], B[2 * sep], B[1 * sep], B[0 * sep])); B += sep * 4;
				pd1[0] = _mm256_mul_pd(pd1[0], _mm256_set_pd(A[3 * sep], A[2 * sep], A[1 * sep], A[0 * sep])); A += sep * 4;
				pd2[1] = _mm256_mul_pd(pd2[1], _mm256_set_pd(B[3 * sep], B[2 * sep], B[1 * sep], B[0 * sep])); B += sep * 4;
			}

			REP(2) AddExponentAVX(slog1, pd1[kk]);
			REP(2) AddExponentAVX(slog2, pd2[kk]);
		}

		__m128d* pd11 = (__m128d*)pd1;
		__m128d* pd22 = (__m128d*)pd2;
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

TARGETAVX float LogProdDivAVXx(float* A, float* B, int64 n, int64 sep)
{
	int64 i = 0;
	int64 slog1 = 0; double prod1 = 1;
	int64 slog2 = 0; double prod2 = 1;

	if (n >= 64)
	{
		__m256 pd1[2], pd2[2];
		REP(2) pd1[kk] = pd2[kk] = _mm256_set1_ps(1.0f);

		__m256i vindex = _mm256_set_epi32(7 * sep, 6 * sep, 5 * sep, 4 * sep, 3 * sep, 2 * sep, 1 * sep, 0 * sep);

		for (int64 l1 = n - 64; i <= l1; i += 64)
		{
			REP(4)
			{
				pd1[0] = _mm256_mul_ps(pd1[0], _mm256_set_ps(A[7 * sep], A[6 * sep], A[5 * sep], A[4 * sep], A[3 * sep], A[2 * sep], A[1 * sep], A[0 * sep])); A += 8 * sep;
				pd2[0] = _mm256_mul_ps(pd2[0], _mm256_set_ps(B[7 * sep], B[6 * sep], B[5 * sep], B[4 * sep], B[3 * sep], B[2 * sep], B[1 * sep], B[0 * sep])); B += 8 * sep;
				pd1[1] = _mm256_mul_ps(pd1[1], _mm256_set_ps(A[7 * sep], A[6 * sep], A[5 * sep], A[4 * sep], A[3 * sep], A[2 * sep], A[1 * sep], A[0 * sep])); A += 8 * sep;
				pd2[1] = _mm256_mul_ps(pd2[1], _mm256_set_ps(B[7 * sep], B[6 * sep], B[5 * sep], B[4 * sep], B[3 * sep], B[2 * sep], B[1 * sep], B[0 * sep])); B += 8 * sep;
			}

			REP(2) AddExponentAVX(slog1, pd1[kk]);
			REP(2) AddExponentAVX(slog2, pd2[kk]);
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

TARGETAVX int64 CountNonZeroAVX(byte* A, int64 n)
{
	uint64 re = 0;
	int64 i = 0;

	if (n >= 32)
	{
		__m256i a = _mm256_setzero_si256(), z = _mm256_setzero_si256();

		for (int64 l1 = n - 32; i <= l1; i += 32)
		{
			a = _mm256_add_epi64(a, _mm256_sad_epu8(z, _mm256_sub_epi8(z, _mm256_cmpeq_epi8(_mm256_loadu_si256((__m256i*)A), z))));
			A += 32;
		}

		re = i - (simd_u64(a, 0) + simd_u64(a, 1) + simd_u64(a, 2) + simd_u64(a, 3));
	}

	for (; i < n; ++i, ++A)
		if (*A) re++;

	return (int64)re;
}

TARGETAVX double SumAVX(double* A, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;
	double re = 0;

	if (n >= N * sizeof(__m256d) / sizeof(double))
	{
		__m256d s[N], a[N];
		REP(N) s[kk] = _mm256_setzero_pd();

		for (int64 l1 = n - N * sizeof(__m256d) / sizeof(double); i <= l1; i += N * sizeof(__m256d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm256_loadu_pd(A); A += sizeof(__m256d) / sizeof(double); }

			REP(N) s[kk] = _mm256_add_pd(s[kk], a[kk]);
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm256_add_pd(s[kk], s[kk + K]);

		re = _mm256_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *A++;
		re += v1;
	}

	return re;
}

TARGETAVX double SumAVX(float* A, int64 n)
{
	constexpr int N = 8;
	int64 i = 0;
	double re = 0;

	if (n >= N * sizeof(__m128) / sizeof(float))
	{
		__m256d s[N], a[N];
		REP(N) s[kk] = _mm256_setzero_pd();

		for (int64 l1 = n - N * sizeof(__m128) / sizeof(float); i <= l1; i += N * sizeof(__m128) / sizeof(float))
		{
			REP(N) { a[kk] = _mm256_cvtps_pd(_mm_loadu_ps(A)); A += sizeof(__m128) / sizeof(float); }

			REP(N) s[kk] = _mm256_add_pd(s[kk], a[kk]);
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm256_add_pd(s[kk], s[kk + K]);
		
		re = _mm256_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A++;
		re += v1;
	}

	return re;
}

TARGETAVX float SumAVXx(float* A, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;
	float re = 0;

	if (n >= N * sizeof(__m256) / sizeof(float))
	{
		__m256 s[N], a[N];
		REP(N) s[kk] = _mm256_setzero_ps();

		for (int64 l1 = n - N * sizeof(__m256) / sizeof(float); i <= l1; i += N * sizeof(__m256) / sizeof(float))
		{
			REP(N) { a[kk] = _mm256_loadu_ps(A); A += sizeof(__m256) / sizeof(float); }

			REP(N) s[kk] = _mm256_add_ps(s[kk], a[kk]);
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm256_add_ps(s[kk], s[kk + K]);

		re = _mm256_reduce_add_ps(s[0]);
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A++;
		re += v1;
	}

	return re;
}

TARGETAVX int64 SumAVX(byte* A, int64 n)
{
	constexpr int N = 4;
	uint64 re = 0;
	int64 i = 0;

	if (n >= N * sizeof(__m256i) / sizeof(byte))
	{
		__m256i s[N], a[N], z = _mm256_setzero_si256();
		REP(N) s[kk] = _mm256_setzero_si256();

		for (int64 l1 = n - 128; i <= l1; i += 128)
		{
			REP(N) { a[kk] = _mm256_loadu_si256((__m256i*)A); A += sizeof(__m256i) / sizeof(byte); }

			REP(N) { s[kk] = _mm256_add_epi64(s[kk], _mm256_sad_epu8(a[kk], z)); }
		}

		s[0] = _mm256_add_epi64(_mm256_add_epi64(s[0], s[1]), _mm256_add_epi64(s[2], s[3]));

		re += simd_u64(s, 0) + simd_u64(s, 1) + simd_u64(s, 2) + simd_u64(s, 3);
	}

	for (; i < n; ++i)
		re += *A++;

	return re;
}

TARGETAVX double SumAVX(double* A, int64 n, int64 sep)
{
	constexpr int N = 4;
	int64 i = 0;
	double re = 0;

	if (n >= N * sizeof(__m256d) / sizeof(double))
	{
		__m256d s[N], a[N];
		REP(N) s[kk] = _mm256_setzero_pd();
		__m256i vindex = _mm256_set_epi64x(3 * sep, 2 * sep, 1 * sep, 0 * sep);

		for (int64 l1 = n - N * sizeof(__m256d) / sizeof(double); i <= l1; i += N * sizeof(__m256d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm256_i64gather_pd(A, vindex, sizeof(double)); A += sizeof(__m256d) / sizeof(double) * sep; }

			REP(N) s[kk] = _mm256_add_pd(s[kk], a[kk]);
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm256_add_pd(s[kk], s[kk + K]);

		re = _mm256_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, A += sep)
	{
		volatile double v1 = *A;
		re += v1;
	}

	return re;
}

TARGETAVX double SumAVX(float* A, int64 n, int64 sep)
{
	return SumSSE(A, n, sep);

	/*
	int64 i = 0;
	double re = 0;

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

		__m256d s1 = _mm256_setzero_pd();
		__m256i vindex = _mm256_set_epi32(-9 * sep, -10 * sep, -11 * sep, -12 * sep, -13 * sep, -14 * sep, -15 * sep, -16 * sep);
		__m256 a;
		__m128& a1 = *((__m128*) & a + 0);
		__m128& a2 = *((__m128*) & a + 1);

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
			a = _mm256_i32xxgather_ps(A, vindex, sizeof(float));
			s1 = _mm256_add_pd(s1, _mm256_add_pd(_mm256_cvtps_pd(a1), _mm256_cvtps_pd(a2)));
		}

		re = _mm256_reduce_add_pd(s1);

		A -= 8 * sep;
	}

	for (; i < n; ++i, A += sep)
	{
		volatile double v1 = *A;
		re += v1;
	}

	return re;
	*/
}

TARGETAVX float SumAVXx(float* A, int64 n, int64 sep)
{
	return SumSSEx(A, n, sep);

	/*
	constexpr int N = 4;
	int64 i = 0;
	float re = 0;

	if (n >= N * sizeof(__m256) / sizeof(float))
	{
		__m256 s[N];
		REP(N) s[kk] = _mm256_setzero_ps();
		__m256i vindex = _mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0);

		for (int64 l1 = n - N * sizeof(__m256) / sizeof(float); i <= l1; i += N * sizeof(__m256) / sizeof(float))
		{
			REP(N) { s[kk] = _mm256_add_ps(s[kk], _mm256_i32xxgather_ps(A, vindex, sizeof(float))); A += sep * sizeof(__m256) / sizeof(float); }
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm256_add_ps(s[kk], s[kk + K]);

		re = _mm256_reduce_add_ps(s[0]);
	}

	for (; i < n; ++i, A += sep)
	{
		volatile float v1 = *A;
		re += v1;
	}

	return re;
	*/
}

TARGETAVX void SumAVX(double* A, double** B, int64 k, int64 n)
{
	constexpr int N = 16;
	int64 i = 0;

	if (n >= N * sizeof(__m256d) / sizeof(double))
	{
		__m256d a[N];

		for (int64 l1 = n - N * sizeof(__m256d) / sizeof(double); i <= l1; i += N * sizeof(__m256d) / sizeof(double))
		{
			REP(N) a[kk] = _mm256_setzero_pd();

			for (int64 j = 0; j < k; ++j)
				REP(N) a[kk] = _mm256_add_pd(a[kk], _mm256_loadu_pd(&B[j][i + kk * sizeof(__m256d) / sizeof(double)]));

			REP(N) { _mm256_storeu_pd(A, a[kk]); A += sizeof(__m256d) / sizeof(double); }
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

TARGETAVX void SumAVX(float* A, float** B, int64 k, int64 n)
{
	constexpr int N = 16;
	int64 i = 0;

	if (n >= N * sizeof(__m128) / sizeof(float))
	{
		__m256d a[N];

		for (int64 l1 = n - N * sizeof(__m128) / sizeof(float); i <= l1; i += N * sizeof(__m128) / sizeof(float))
		{
			REP(N) a[kk] = _mm256_setzero_pd();

			for (int64 j = 0; j < k; ++j)
				REP(N) a[kk] = _mm256_add_pd(a[kk], _mm256_cvtps_pd(_mm_loadu_ps(&B[j][i + kk * sizeof(__m128) / sizeof(float)])));

			REP(N) { _mm_storeu_ps(A, _mm256_cvtpd_ps(a[kk])); A += sizeof(__m128) / sizeof(float); }
		}
	}

	for (; i < n; ++i)
	{
		double Ai = 0;
		for (int64 j = 0; j < k; ++j)
			Ai += B[j][i];
		*A++ = Ai;
	}
	
	/*
	int64 i = 0;

	if (n >= 16)
	{
		__m256 a1, a2;
		for (int64 l1 = n - 16; i <= l1; i += 16)
		{
			a1 = _mm256_setzero_ps();
			a2 = _mm256_setzero_ps();

			for (int64 j = 0; j < k; ++j)
			{
				a1 = _mm256_add_ps(a1, _mm256_loadu_ps(&B[j][i + 0]));
				a2 = _mm256_add_ps(a2, _mm256_loadu_ps(&B[j][i + 8]));
			}

			_mm256_storeu_ps(&A[i + 0], a1);
			_mm256_storeu_ps(&A[i + 8], a2);
		}
	}

	for (; i < n; ++i)
	{
		double Ai = 0;
		for (int64 j = 0; j < k; ++j)
			Ai += B[j][i];
		A[i] = Ai;
	}
	*/
}

TARGETAVX double ProdAVX(double* A, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;
	double re = 1;

	if (n >= N * sizeof(__m256d) / sizeof(double))
	{
		__m256d s[N], a[N];
		REP(N) s[kk] = _mm256_set1_pd(1);

		for (int64 l1 = n - N * sizeof(__m256d) / sizeof(double); i <= l1; i += N * sizeof(__m256d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm256_loadu_pd(A); A += sizeof(__m256d) / sizeof(double); }

			REP(N) s[kk] = _mm256_mul_pd(s[kk], a[kk]);
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm256_mul_pd(s[kk], s[kk + K]);

		re = _mm256_reduce_mul_pd(s[0]);
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *A++;
		re *= v1;
	}

	return re;
}

TARGETAVX double ProdAVX(float* A, int64 n)
{
	constexpr int N = 8;
	int64 i = 0;
	double re = 1;

	if (n >= N * sizeof(__m128) / sizeof(float))
	{
		__m256d s[N], a[N];
		REP(N) s[kk] = _mm256_set1_pd(1);

		for (int64 l1 = n - N * sizeof(__m128) / sizeof(float); i <= l1; i += N * sizeof(__m128) / sizeof(float))
		{
			REP(N) { a[kk] = _mm256_cvtps_pd(_mm_loadu_ps(A)); A += sizeof(__m128) / sizeof(float); }

			REP(N) s[kk] = _mm256_mul_pd(s[kk], a[kk]);
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm256_mul_pd(s[kk], s[kk + K]);

		re = _mm256_reduce_mul_pd(s[0]);
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A++;
		re *= v1;
	}

	return re;
}

TARGETAVX float ProdAVXx(float* A, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;
	volatile float re = 1;

	if (n >= N * sizeof(__m256) / sizeof(float))
	{
		__m256 s[N], a[N];
		REP(N) s[kk] = _mm256_set1_ps(1);

		for (int64 l1 = n - N * sizeof(__m256) / sizeof(float); i <= l1; i += N * sizeof(__m256) / sizeof(float))
		{
			REP(N) { a[kk] = _mm256_loadu_ps(A); A += sizeof(__m256) / sizeof(float); }

			REP(N) s[kk] = _mm256_mul_ps(s[kk], a[kk]);
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm256_mul_ps(s[kk], s[kk + K]);

		re = _mm256_reduce_mul_ps(s[0]);
	}

	for (; i < n; ++i)
		re *= *A++;

	return re;
}

TARGETAVX double ProdAVX(double* A, int64 n, int64 sep)
{
	return ProdSSE(A, n, sep);
}

TARGETAVX double ProdAVX(float* A, int64 n, int64 sep)
{
	return ProdSSE(A, n, sep);
}

TARGETAVX float ProdAVXx(float* A, int64 n, int64 sep)
{
	return ProdSSEx(A, n, sep);
}

TARGETAVX double SumSquareAVX(double* A, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;
	double re = 0;

	if (n >= N * sizeof(__m256d) / sizeof(double))
	{
		__m256d s[N], a[N];
		REP(N) s[kk] = _mm256_setzero_pd();

		for (int64 l1 = n - N * sizeof(__m256d) / sizeof(double); i <= l1; i += N * sizeof(__m256d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm256_loadu_pd(A); A += sizeof(__m256d) / sizeof(double); }

			REP(N) a[kk] = _mm256_mul_pd(a[kk], a[kk]);

			REP(N) s[kk] = _mm256_add_pd(s[kk], a[kk]);
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm256_add_pd(s[kk], s[kk + K]);

		re = _mm256_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, ++A)
	{
		volatile double v1 = *A * *A;
		re += v1;
	}

	return re;
}

TARGETAVX double SumSquareAVX(float* A, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;
	double re = 0;

	if (n >= N * sizeof(__m128) / sizeof(float))
	{
		__m256d s[N], a[N];
		REP(N) s[kk] = _mm256_setzero_pd();

		for (int64 l1 = n - N * sizeof(__m128) / sizeof(float); i <= l1; i += N * sizeof(__m128) / sizeof(float))
		{
			REP(N) { a[kk] = _mm256_cvtps_pd(_mm_loadu_ps(A)); A += sizeof(__m128) / sizeof(float); }
			REP(N) a[kk] = _mm256_mul_pd(a[kk], a[kk]);
			REP(N) s[kk] = _mm256_add_pd(s[kk], a[kk]);
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm256_add_pd(s[kk], s[kk + K]);

		re = _mm256_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, ++A)
	{
		volatile double v1 = (double)*A * (double)*A;
		re += v1;
	}

	return re;
}

TARGETAVX float SumSquareAVXx(float* A, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;
	float re = 0;

	if (n >= N * sizeof(__m256) / sizeof(float))
	{
		__m256 s[N], a[N];
		REP(N) s[kk] = _mm256_setzero_ps();

		for (int64 l1 = n - N * sizeof(__m256) / sizeof(float); i <= l1; i += N * sizeof(__m256) / sizeof(float))
		{
			REP(N) { a[kk] = _mm256_loadu_ps(A); A += sizeof(__m256) / sizeof(float); }

			REP(N) a[kk] = _mm256_mul_ps(a[kk], a[kk]);

			REP(N) s[kk] = _mm256_add_ps(s[kk], a[kk]);
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm256_add_ps(s[kk], s[kk + K]);

		re = _mm256_reduce_add_ps(s[0]);
	}

	for (; i < n; ++i, ++A)
	{
		volatile float v1 = *A * *A;
		re += v1;
	}

	return re;
}

TARGETAVX int64 SumSquareAVX(byte* A, int64 n)
{
	constexpr int N = 2;
	int64 i = 0;
	uint64 re = 0;

	if (n >= N * sizeof(__m256i) / sizeof(byte))
	{
		__m256i a, t[2], s = _mm256_setzero_si256();
		REP(2) t[kk] = _mm256_setzero_si256();
		__m128i* s2 = (__m128i*) & s;

		for (int64 l1 = n - N * sizeof(__m256i) / sizeof(byte); i <= l1; i += N * sizeof(__m256i) / sizeof(byte))
		{
			REP(N)
			{
				a = _mm256_loadu_si256((__m256i*)A);
				A += sizeof(__m256i) / sizeof(byte);
				s = _mm256_add_epi16(s, _mm256_maddubs_epi16(a, a));
			}

			if ((i & (sizeof(__m256i) / sizeof(byte) * 128 - 1)) == 0)
			{
				t[0] = _mm256_add_epi32(t[0], _mm256_cvtepi16_epi32(s2[0]));
				t[1] = _mm256_add_epi32(t[1], _mm256_cvtepi16_epi32(s2[1]));
				s = _mm256_setzero_si256();
			}
		}

		t[0] = _mm256_add_epi32(t[0], _mm256_cvtepi16_epi32(s2[0]));
		t[1] = _mm256_add_epi32(t[1], _mm256_cvtepi16_epi32(s2[1]));
		t[0] = _mm256_add_epi32(t[0], t[1]);

		re = simp_i32(t, 0) + simp_i32(t, 1) + simp_i32(t, 2) + simp_i32(t, 3) + 
			 simp_i32(t, 4) + simp_i32(t, 5) + simp_i32(t, 6) + simp_i32(t, 7);
	}

	for (; i < n; ++i, ++A)
		re += *A * *A;

	return re;
}

TARGETAVX void SumSumSquareAVX(double* A, int64 n, double& sum, double& sumsq)
{
	constexpr int N = 4;
	int64 i = 0;
	double re1 = 0, re2 = 0;

	if (n >= N * sizeof(__m256d) / sizeof(double))
	{
		__m256d s1[N], s2[N], a[N];
		REP(N) s1[kk] = s2[kk] = _mm256_setzero_pd();

		for (int64 l1 = n - N * sizeof(__m256d) / sizeof(double); i <= l1; i += N * sizeof(__m256d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm256_loadu_pd(A); A += sizeof(__m256d) / sizeof(double); }

			REP(N) s1[kk] = _mm256_add_pd(s1[kk], a[kk]);

			REP(N) a[kk] = _mm256_mul_pd(a[kk], a[kk]);

			REP(N) s2[kk] = _mm256_add_pd(s2[kk], a[kk]);
		}

		for (int K = sizeof(s1) / sizeof(s1[0]) / 2; K >= 1; K >>= 1)
		{
			REP(K) s1[kk] = _mm256_add_pd(s1[kk], s1[kk + K]);
			REP(K) s2[kk] = _mm256_add_pd(s2[kk], s2[kk + K]);
		}

		re1 = _mm256_reduce_add_pd(s1[0]);
		re2 = _mm256_reduce_add_pd(s2[0]);
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

TARGETAVX void SumSumSquareAVX(float* A, int64 n, double& sum, double& sumsq)
{
	constexpr int N = 8;
	int64 i = 0;
	double re1 = 0, re2 = 0;

	if (n >= N * sizeof(__m128) / sizeof(float))
	{
		__m256d s1[N], s2[N], a[N];
		REP(N) s1[kk] = s2[kk] = _mm256_setzero_pd();

		for (int64 l1 = n - N * sizeof(__m128) / sizeof(float); i <= l1; i += N * sizeof(__m128) / sizeof(float))
		{
			REP(N) { a[kk] = _mm256_cvtps_pd(_mm_loadu_ps(A)); A += sizeof(__m128) / sizeof(float); }

			REP(N) s1[kk] = _mm256_add_pd(s1[kk], a[kk]);

			REP(N) a[kk] = _mm256_mul_pd(a[kk], a[kk]);

			REP(N) s2[kk] = _mm256_add_pd(s2[kk], a[kk]);
		}

		for (int K = sizeof(s1) / sizeof(s1[0]) / 2; K >= 1; K >>= 1)
		{
			REP(K) s1[kk] = _mm256_add_pd(s1[kk], s1[kk + K]);
			REP(K) s2[kk] = _mm256_add_pd(s2[kk], s2[kk + K]);
		}

		re1 = _mm256_reduce_add_pd(s1[0]);
		re2 = _mm256_reduce_add_pd(s2[0]);
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

TARGETAVX double SumProdDivAVX(double* A1, double* A2, double* B, int64 sep, int64 n)
{
	int64 i = 0;
	volatile double re1 = 0, re2 = 0;

	if (n >= 8)
	{
		REP(8) { _mm_prefetch((const char*)B, _MM_HINT_T0); B += sep; }

		__m256d s1 = _mm256_setzero_pd(), s2 = _mm256_setzero_pd(), b;
		__m256i vindex = _mm256_set_epi64x(-9 * sep, -10 * sep, -11 * sep, -12 * sep);

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			REP(4) { _mm_prefetch((const char*)B, _MM_HINT_T0); B += sep; }

			b = _mm256_i64gather_pd(B, vindex, sizeof(double));

			s1 = _mm256_add_pd(s1, _mm256_mul_pd(_mm256_loadu_pd(A1), b)); A1 += 4;
			s2 = _mm256_add_pd(s2, _mm256_mul_pd(_mm256_loadu_pd(A2), b)); A2 += 4;

			REP(4) { _mm_prefetch((const char*)B, _MM_HINT_T0); B += sep; }

			b = _mm256_i64gather_pd(B, vindex, sizeof(double));

			s1 = _mm256_add_pd(s1, _mm256_mul_pd(_mm256_loadu_pd(A1), b)); A1 += 4;
			s2 = _mm256_add_pd(s2, _mm256_mul_pd(_mm256_loadu_pd(A2), b)); A2 += 4;
		}

		re1 = _mm256_reduce_add_pd(s1);
		re2 = _mm256_reduce_add_pd(s2);

		B -= 8 * sep;
	}

	for (; i < n; ++i, A1++, A2++, B += sep)
	{
		volatile double v1 = *A1 * *B;
		volatile double v2 = *A2 * *B;
		re1 += v1;
		re2 += v2;
	}

	return re1 / re2;
}

TARGETAVX double SumProdDivAVX(double* A1, float* A2, float* B, int64 sep, int64 n)
{
	return SumProdDivSSE(A1, A2, B, sep, n);

	/*
	int64 i = 0;
	double re1 = 0, re2 = 0;

	if (n >= 16)
	{
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;

		__m256d s1 = _mm256_setzero_pd(), s2 = _mm256_setzero_pd(), bb1, bb2;
		__m256i vindex = _mm256_set_epi32(-9 * sep, -10 * sep, -11 * sep, -12 * sep, -13 * sep, -14 * sep, -15 * sep, -16 * sep);
		__m256 a2, b;
		__m128& a21 = *((__m128*) & a2 + 0);
		__m128& a22 = *((__m128*) & a2 + 1);
		__m128& b1 = *((__m128*) & b + 0);
		__m128& b2 = *((__m128*) & b + 1);

		for (int64 l1 = n - 16; i <= l1; i += 16)
		{
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;

			b = _mm256_i32xxgather_ps(B, vindex, sizeof(float));
			a2 = _mm256_loadu_ps(A2); A2 += 8;
			bb1 = _mm256_cvtps_pd(b1);
			bb2 = _mm256_cvtps_pd(b2);

			s1 = _mm256_fmadd_pd(_mm256_loadu_pd(A1), bb1, s1); A1 += 4;
			s2 = _mm256_fmadd_pd(_mm256_cvtps_pd(a21), bb1, s2);
			s1 = _mm256_fmadd_pd(_mm256_loadu_pd(A1), bb2, s1); A1 += 4;
			s2 = _mm256_fmadd_pd(_mm256_cvtps_pd(a22), bb2, s2);

			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;

			b = _mm256_i32xxgather_ps(B, vindex, sizeof(float));
			a2 = _mm256_loadu_ps(A2); A2 += 8;
			bb1 = _mm256_cvtps_pd(b1);
			bb2 = _mm256_cvtps_pd(b2);

			s1 = _mm256_fmadd_pd(_mm256_loadu_pd(A1), bb1, s1); A1 += 4;
			s2 = _mm256_fmadd_pd(_mm256_cvtps_pd(a21), bb1, s2);
			s1 = _mm256_fmadd_pd(_mm256_loadu_pd(A1), bb2, s1); A1 += 4;
			s2 = _mm256_fmadd_pd(_mm256_cvtps_pd(a22), bb2, s2);
		}

		re1 = _mm256_reduce_add_pd(s1);
		re2 = _mm256_reduce_add_pd(s2);

		B -= 8 * sep;
	}

	for (; i < n; ++i, A1++, A2++, B += sep)
	{
		volatile double v1 = (double)*A1 * (double)*B;
		volatile double v2 = (double)*A2 * (double)*B;
		re1 += v1;
		re2 += v2;
	}

	return re1 / re2;
	*/
}

TARGETAVX double SumProdDivAVX(float* A1, float* A2, float* B, int64 sep, int64 n)
{
	return SumProdDivSSE(A1, A2, B, sep, n);

	/*
	int64 i = 0;
	double re1 = 0, re2 = 0;

	if (n >= 16)
	{
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;

		__m256d s1 = _mm256_setzero_pd(), s2 = _mm256_setzero_pd(), bb1, bb2;
		__m256i vindex = _mm256_set_epi32(-9 * sep, -10 * sep, -11 * sep, -12 * sep, -13 * sep, -14 * sep, -15 * sep, -16 * sep);
		__m256 a1, a2, b;
		__m128& a11 = *((__m128*)&a1 + 0);
		__m128& a12 = *((__m128*)&a1 + 1);
		__m128& a21 = *((__m128*)&a2 + 0);
		__m128& a22 = *((__m128*)&a2 + 1);
		__m128& b1 = *((__m128*)&b + 0);
		__m128& b2 = *((__m128*)&b + 1);

		for (int64 l1 = n - 16; i <= l1; i += 16)
		{
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;

			b = _mm256_i32xxgather_ps(B, vindex, sizeof(float));
			a1 = _mm256_loadu_ps(A1); A1 += 8;
			a2 = _mm256_loadu_ps(A2); A2 += 8;
			bb1 = _mm256_cvtps_pd(b1);
			bb2 = _mm256_cvtps_pd(b2);

			s1 = _mm256_add_pd(s1, _mm256_mul_pd(_mm256_cvtps_pd(a11), bb1));
			s2 = _mm256_add_pd(s2, _mm256_mul_pd(_mm256_cvtps_pd(a21), bb1));
			s1 = _mm256_add_pd(s1, _mm256_mul_pd(_mm256_cvtps_pd(a12), bb2));
			s2 = _mm256_add_pd(s2, _mm256_mul_pd(_mm256_cvtps_pd(a22), bb2));

			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;

			b = _mm256_i32xxgather_ps(B, vindex, sizeof(float));
			a1 = _mm256_loadu_ps(A1); A1 += 8;
			a2 = _mm256_loadu_ps(A2); A2 += 8;
			bb1 = _mm256_cvtps_pd(b1);
			bb2 = _mm256_cvtps_pd(b2);

			s1 = _mm256_add_pd(s1, _mm256_mul_pd(_mm256_cvtps_pd(a11), bb1));
			s2 = _mm256_add_pd(s2, _mm256_mul_pd(_mm256_cvtps_pd(a21), bb1));
			s1 = _mm256_add_pd(s1, _mm256_mul_pd(_mm256_cvtps_pd(a12), bb2));
			s2 = _mm256_add_pd(s2, _mm256_mul_pd(_mm256_cvtps_pd(a22), bb2));
		}

		re1 = _mm256_reduce_add_pd(s1);
		re2 = _mm256_reduce_add_pd(s2);

		B -= 8 * sep;
	}

	for (; i < n; ++i, A1++, A2++, B += sep)
	{
		volatile double v1 = (double)*A1 * (double)*B;
		volatile double v2 = (double)*A2 * (double)*B;
		re1 += v1;
		re2 += v2;
	}

	return re1 / re2;
	*/
}

TARGETAVX float SumProdDivAVXx(float* A1, float* A2, float* B, int64 sep, int64 n)
{
	constexpr int N = 2;
	int64 i = 0;
	volatile float re1 = 0, re2 = 0;

	if (n >= N * sizeof(__m256) / sizeof(float))
	{
		__m256 s1[N], s2[N], b;
		REP(N) s1[kk] = s2[kk] = _mm256_setzero_ps();

		for (int64 l1 = n - N * sizeof(__m256) / sizeof(float); i <= l1; i += N * sizeof(__m256) / sizeof(float))
		{
			REP(N)
			{
				b = _mm256_set_ps(B[7 * sep], B[6 * sep], B[5 * sep], B[4 * sep],
							      B[3 * sep], B[2 * sep], B[1 * sep], B[0 * sep]);
				B += sizeof(__m256) / sizeof(float) * sep;
				s1[kk] = _mm256_add_ps(s1[kk], _mm256_mul_ps(_mm256_loadu_ps(A1), b)); A1 += sizeof(__m256) / sizeof(float);
				s2[kk] = _mm256_add_ps(s2[kk], _mm256_mul_ps(_mm256_loadu_ps(A2), b)); A2 += sizeof(__m256) / sizeof(float);
			}
		}

		for (int K = sizeof(s1) / sizeof(s1[0]) / 2; K >= 1; K >>= 1)
		{
			REP(K) s1[kk] = _mm256_add_ps(s1[kk], s1[kk + K]);
			REP(K) s2[kk] = _mm256_add_ps(s2[kk], s2[kk + K]);
		}

		re1 = _mm256_reduce_add_ps(s1[0]);
		re2 = _mm256_reduce_add_ps(s2[0]);
	}

	for (; i < n; ++i, A1++, A2++, B += sep)
	{
		volatile float v1 = *A1 * *B;
		volatile float v2 = *A2 * *B;
		re1 += v1;
		re2 += v2;
	}

	return re1 / re2;
}

TARGETAVX double SumProdAVX(double* A, double* B, int64 sep, int64 n)
{
	return SumProdSSE(A, B, sep, n);
}

TARGETAVX double SumProdAVX(float* A, float* B, int64 sep, int64 n)
{
	return SumProdSSE(A, B, sep, n);
}

TARGETAVX float SumProdAVXx(float* A, float* B, int64 sep, int64 n)
{
	return SumProdSSEx(A, B, sep, n);
}

TARGETAVX double SumProdAVX(double* A, double* B, int64 n)
{
	constexpr int N = 8;
	int64 i = 0;
	double re = 0;

	if (n >= N * sizeof(__m256d) / sizeof(double))
	{
		__m256d s[N], a[N];
		REP(N) s[kk] = _mm256_setzero_pd();

		for (int64 l1 = n - N * sizeof(__m256d) / sizeof(double); i <= l1; i += N * sizeof(__m256d) / sizeof(double))
		{
			REP(N)
			{
				a[kk] = _mm256_mul_pd(_mm256_loadu_pd(A), _mm256_loadu_pd(B));
				A += sizeof(__m256d) / sizeof(double); 
				B += sizeof(__m256d) / sizeof(double);
			}

			REP(N) s[kk] = _mm256_add_pd(s[kk], a[kk]);
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm256_add_pd(s[kk], s[kk + K]);

		re = _mm256_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, ++A, ++B)
	{
		volatile double v1 = *A * *B;
		re += v1;
	}

	return re;
}

TARGETAVX double SumProdAVX(float* A, float* B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;
	double re = 0;

	if (n >= N * sizeof(__m128) / sizeof(float))
	{
		__m256d s[N];
		REP(N) s[kk] = _mm256_setzero_pd();

		for (int64 l1 = n - N * sizeof(__m128) / sizeof(float); i <= l1; i += N * sizeof(__m128) / sizeof(float))
		{
			REP(N)
			{
				s[kk] = _mm256_add_pd(s[kk], _mm256_mul_pd(_mm256_cvtps_pd(_mm_loadu_ps(A)), _mm256_cvtps_pd(_mm_loadu_ps(B))));
				A += sizeof(__m128) / sizeof(float);
				B += sizeof(__m128) / sizeof(float);
			}
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm256_add_pd(s[kk], s[kk + K]);

		re = _mm256_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, ++A, ++B)
	{
		volatile double v1 = (double)*A * (double)*B;
		re += v1;
	}

	return re;
}

TARGETAVX float SumProdAVXx(float* A, float* B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;
	volatile float re = 0;

	if (n >= N * sizeof(__m256) / sizeof(float))
	{
		__m256 s[N];
		REP(N) s[kk] = _mm256_setzero_ps();

		for (int64 l1 = n - N * sizeof(__m256) / sizeof(float); i <= l1; i += N * sizeof(__m256) / sizeof(float))
		{
			REP(N)
			{
				s[kk] = _mm256_add_ps(s[kk], _mm256_mul_ps(_mm256_loadu_ps(A), _mm256_loadu_ps(B)));
				A += sizeof(__m256) / sizeof(float);
				B += sizeof(__m256) / sizeof(float);
			}
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm256_add_ps(s[kk], s[kk + K]);

		re = _mm256_reduce_add_ps(s[0]);
	}

	for (; i < n; ++i, ++A, ++B)
	{
		volatile float v1 = *A * *B;
		re += v1;
	}

	return re;
}

TARGETAVX void AddAVX(double* A, double* B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m256d) / sizeof(double))
	{
		__m256d a[N], b[N];

		for (int64 l1 = n - N * sizeof(__m256d) / sizeof(double); i <= l1; i += N * sizeof(__m256d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm256_loadu_pd(A); A += 4; }

			REP(N) { b[kk] = _mm256_loadu_pd(B); B += 4; }

			REP(N) a[kk] = _mm256_add_pd(a[kk], b[kk]);

			REP(N) _mm256_storeu_pd(A + (kk - N) * sizeof(__m256d) / sizeof(double), a[kk]);
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGETAVX void AddAVX(float* A, float* B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m256) / sizeof(float))
	{
		__m256 a[N], b[N];

		for (int64 l1 = n - N * sizeof(__m256) / sizeof(float); i <= l1; i += N * sizeof(__m256) / sizeof(float))
		{
			REP(N) { a[kk] = _mm256_loadu_ps(A); A += sizeof(__m256) / sizeof(float); }
	
			REP(N) { b[kk] = _mm256_loadu_ps(B); B += sizeof(__m256) / sizeof(float); }

			REP(N) a[kk] = _mm256_add_ps(a[kk], b[kk]);

			REP(N) _mm256_storeu_ps(A + (kk - N) * sizeof(__m256) / sizeof(float), a[kk]);
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGETAVX void AddAVX(int64* A, int64* B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m256i) / sizeof(int64))
	{
		__m256i a[N], b[N];

		for (int64 l1 = n - N * sizeof(__m256i) / sizeof(int64); i <= l1; i += N * sizeof(__m256i) / sizeof(int64))
		{
			REP(N) { a[kk] = _mm256_loadu_si256((__m256i*)A); A += sizeof(__m256i) / sizeof(int64); }

			REP(N) { b[kk] = _mm256_loadu_si256((__m256i*)B); B += sizeof(__m256i) / sizeof(int64); }

			REP(N) a[kk] = _mm256_add_epi64(a[kk], b[kk]);

			REP(N) _mm256_storeu_si256((__m256i*)(A + (kk - N) * sizeof(__m256i) / sizeof(int64)), a[kk]);
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGETAVX void AddAVX(int* A, int* B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m256i) / sizeof(int))
	{
		__m256i a[N], b[N];

		for (int64 l1 = n - N * sizeof(__m256i) / sizeof(int); i <= l1; i += N * sizeof(__m256i) / sizeof(int))
		{
			REP(N) { a[kk] = _mm256_loadu_si256((__m256i*)A); A += 8; }

			REP(N) { b[kk] = _mm256_loadu_si256((__m256i*)B); B += 8; }

			REP(N) a[kk] = _mm256_add_epi32(a[kk], b[kk]);

			REP(N) _mm256_storeu_si256((__m256i*)(A + (kk - N) * sizeof(__m256i) / sizeof(int)), a[kk]);
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGETAVX void AddAVX(double* A, double B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m256d) / sizeof(double))
	{
		__m256d b = _mm256_set1_pd(B), a[N];

		for (int64 l1 = n - N * sizeof(__m256d) / sizeof(double); i <= l1; i += N * sizeof(__m256d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm256_loadu_pd(A); A += sizeof(__m256d) / sizeof(double); }

			REP(N) a[kk] = _mm256_add_pd(a[kk], b);

			REP(N) _mm256_storeu_pd(A + (kk - N) * sizeof(__m256d) / sizeof(double), a[kk]);
		}
	}

	for (; i < n; ++i, A++)
		*A += B;
}

TARGETAVX void AddAVX(float* A, float B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m256) / sizeof(float))
	{
		__m256 b = _mm256_set1_ps(B), a[N];

		for (int64 l1 = n - N * sizeof(__m256) / sizeof(float); i <= l1; i += N * sizeof(__m256) / sizeof(float))
		{
			REP(N) { a[kk] = _mm256_loadu_ps(A); A += sizeof(__m256) / sizeof(float); }

			REP(N) a[kk] = _mm256_add_ps(a[kk], b);

			REP(N) _mm256_storeu_ps(A + (kk - N) * sizeof(__m256) / sizeof(float), a[kk]);
		}
	}

	for (; i < n; ++i, A++)
		*A += B;
}

TARGETAVX void MulAVX(double* C, double* A, double* B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m256d) / sizeof(double))
	{
		__m256d a[N], b[N];

		for (int64 l1 = n - N * sizeof(__m256d) / sizeof(double); i <= l1; i += N * sizeof(__m256d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm256_loadu_pd(A); A += sizeof(__m256d) / sizeof(double); }

			REP(N) { b[kk] = _mm256_loadu_pd(B); B += sizeof(__m256d) / sizeof(double); }

			REP(N) a[kk] = _mm256_mul_pd(a[kk], b[kk]);

			REP(N) { _mm256_storeu_pd(C, a[kk]); C += sizeof(__m256d) / sizeof(double); }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *A++ * *B++;
		*C++ = v1;
	}
}

TARGETAVX void MulAVX(float* C, float* A, float* B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m256) / sizeof(float))
	{
		__m256 a[N], b[N];

		for (int64 l1 = n - N * sizeof(__m256) / sizeof(float); i <= l1; i += N * sizeof(__m256) / sizeof(float))
		{
			REP(N) { a[kk] = _mm256_loadu_ps(A); A += sizeof(__m256) / sizeof(float); }

			REP(N) { b[kk] = _mm256_loadu_ps(B); B += sizeof(__m256) / sizeof(float); }

			REP(N) a[kk] = _mm256_mul_ps(a[kk], b[kk]);

			REP(N) { _mm256_storeu_ps(C, a[kk]); C += sizeof(__m256) / sizeof(float); }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A++ * *B++;
		*C++ = v1;
	}
}

TARGETAVX void MulAVX(double* C, double* A, double B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m256d) / sizeof(double))
	{
		__m256d b = _mm256_set1_pd(B), a[N];

		for (int64 l1 = n - N * sizeof(__m256d) / sizeof(double); i <= l1; i += N * sizeof(__m256d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm256_loadu_pd(A); A += sizeof(__m256d) / sizeof(double); }

			REP(N) a[kk] = _mm256_mul_pd(a[kk], b);

			REP(N) { _mm256_storeu_pd(C, a[kk]); C += sizeof(__m256d) / sizeof(double); }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *A++ * B;
		*C++ = v1;
	}
}

TARGETAVX void MulAVX(float* C, float* A, float B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m256) / sizeof(float))
	{
		__m256 b = _mm256_set1_ps(B), a[N];

		for (int64 l1 = n - N * sizeof(__m256) / sizeof(float); i <= l1; i += N * sizeof(__m256) / sizeof(float))
		{
			REP(N) { a[kk] = _mm256_loadu_ps(A); A += sizeof(__m256) / sizeof(float); }

			REP(N) a[kk] = _mm256_mul_ps(a[kk], b);

			REP(N) { _mm256_storeu_ps(C, a[kk]); C += sizeof(__m256) / sizeof(float); }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A++ * B;
		*C++ = v1;
	}
}

TARGETAVX void MulAVX(double* A, double B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m256d) / sizeof(double))
	{
		__m256d b = _mm256_set1_pd(B), a[N];

		for (int64 l1 = n - N * sizeof(__m256d) / sizeof(double); i <= l1; i += N * sizeof(__m256d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm256_loadu_pd(A); A += sizeof(__m256d) / sizeof(double); }

			REP(N) a[kk] = _mm256_mul_pd(a[kk], b);

			REP(N) _mm256_storeu_pd(A + (kk - N) * sizeof(__m256d) / sizeof(double), a[kk]);
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *A * B;
		*A++ = v1;
	}
}

TARGETAVX void MulAVX(float* A, float B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m256) / sizeof(float))
	{
		__m256 b = _mm256_set1_ps(B), a[N];

		for (int64 l1 = n - N * sizeof(__m256) / sizeof(float); i <= l1; i += N * sizeof(__m256) / sizeof(float))
		{
			REP(N) { a[kk] = _mm256_loadu_ps(A); A += sizeof(__m256) / sizeof(float); }

			REP(N) a[kk] = _mm256_mul_ps(a[kk], b);

			REP(N) _mm256_storeu_ps(A + (kk - N) * sizeof(__m256) / sizeof(float), a[kk]);
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A * B;
		*A++ = v1;
	}
}

TARGETAVX void AddProdAVX(double* C, double* A, double* B, int64 n)
{
	constexpr int N = 2;

	int64 i = 0;

	if (n >= N * sizeof(__m256d) / sizeof(double))
	{
		__m256d a[N];

		for (int64 l1 = n - N * sizeof(__m256d) / sizeof(double); i <= l1; i += N * sizeof(__m256d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm256_loadu_pd(A); A += sizeof(__m256d) / sizeof(double); }

			REP(N) { a[kk] = _mm256_mul_pd(a[kk], _mm256_loadu_pd(B)); B += sizeof(__m256d) / sizeof(double); }

			double* C2 = C; 

			REP(N) { a[kk] = _mm256_add_pd(a[kk], _mm256_loadu_pd(C)); C += sizeof(__m256d) / sizeof(double); }

			REP(N) { _mm256_storeu_pd(C2, a[kk]); C2 += sizeof(__m256d) / sizeof(double); }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *A++ * *B++;
		*C++ += v1;
	}
}

TARGETAVX void AddProdAVX(float* C, float* A, float* B, int64 n)
{
	return AddProdSSE(C, A, B, n);

	/*
	constexpr int N = 4;
	int64 i = 0;

	if (n >= 2 * N * sizeof(__m256) / sizeof(float))
	{
		__m256 a[N];
		byte maskff = 0 - (n > -1);

		for (int64 l1 = n - N * sizeof(__m256) / sizeof(float); i <= l1; i += N * sizeof(__m256) / sizeof(float))
		{
			REP(N) { a[kk] = _mm256_loadu_ps(A); A += sizeof(__m256) / sizeof(float); }

			REP(N) { a[kk] = _mm256_mul_ps(a[kk], _mm256_loadu_ps(B)); B += sizeof(__m256) / sizeof(float); }

			float* C2 = C;

			REP(N) { a[kk] = _mm256_add_ps(a[kk], _mm256_loadu_ps(C)); C += sizeof(__m256) / sizeof(float); }

			REP(N) { _mm256_storeu_ps(C2, a[kk]); C2 += sizeof(__m256) / sizeof(float); }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A++ * *B++;
		*C++ += v1;
	}
	*/
}

TARGETAVX void AddProdAVX(double* C, double* A, double B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= 2 * N * sizeof(__m256d) / sizeof(double))
	{
		__m256d a[N], b = _mm256_set1_pd(B);

		for (int64 l1 = n - N * sizeof(__m256d) / sizeof(double); i <= l1; i += N * sizeof(__m256d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm256_loadu_pd(A); A += sizeof(__m256d) / sizeof(double); }

			REP(N)  a[kk] = _mm256_mul_pd(a[kk], b);

			double* C2 = C;

			REP(N) { a[kk] = _mm256_add_pd(a[kk], _mm256_loadu_pd(C)); C += sizeof(__m256d) / sizeof(double); }

			REP(N) { _mm256_storeu_pd(C2, a[kk]); C2 += sizeof(__m256d) / sizeof(double); }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *A++ * B;
		*C++ += v1;
	}
}

TARGETAVX void AddProdAVX(double* C, float* A, double B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= 2 * N * sizeof(__m256d) / sizeof(double))
	{
		__m256d a[N], b = _mm256_set1_pd(B);

		for (int64 l1 = n - N * sizeof(__m256d) / sizeof(double); i <= l1; i += N * sizeof(__m256d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm256_cvtps_pd(_mm_loadu_ps(A)); A += sizeof(__m128) / sizeof(float); }

			REP(N)  a[kk] = _mm256_mul_pd(a[kk], b);

			double* C2 = C;

			REP(N) { a[kk] = _mm256_add_pd(a[kk], _mm256_loadu_pd(C)); C += sizeof(__m256d) / sizeof(double); }

			REP(N) { _mm256_storeu_pd(C2, a[kk]); C2 += sizeof(__m256d) / sizeof(double); }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *A++ * B;
		*C++ += v1;
	}
}

TARGETAVX void AddProdAVX(float* C, float* A, float B, int64 n)
{
	return AddProdSSE(C, A, B, n);
	/*
	constexpr int N = 4;
	int64 i = 0;

	if (n >= 2 * N * sizeof(__m256) / sizeof(float))
	{
		__m256 a[N], b = _mm256_set1_ps(B);

		for (int64 l1 = n - N * sizeof(__m256) / sizeof(float); i <= l1; i += N * sizeof(__m256) / sizeof(float))
		{
			REP(N) { a[kk] = _mm256_loadu_ps(A); A += sizeof(__m256) / sizeof(float); }

			REP(N)  a[kk] = _mm256_mul_ps(a[kk], b);

			float* C2 = C;

			REP(N) { a[kk] = _mm256_add_ps(a[kk], _mm256_loadu_ps(C)); C += sizeof(__m256) / sizeof(float); }

			REP(N) { _mm256_storeu_ps(C2, a[kk]); C2 += sizeof(__m256) / sizeof(float); }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A++ * B;
		*C++ += v1;
	}
	*/
}

TARGETAVX void UnifyAVX(double* A, int64 n)
{
	constexpr int N = 8;
	int64 i = 0;
	double invsum = 1.0 / (SumAVX(A, n) + n * MIN_FREQ);

	if (n >= N * sizeof(__m256) / sizeof(double))
	{
		__m256d a[N], minf = _mm256_set1_pd(MIN_FREQ), invs = _mm256_set1_pd(invsum);

		for (int64 l1 = n - N * sizeof(__m256) / sizeof(double); i <= l1; i += N * sizeof(__m256) / sizeof(double))
		{
			double* A2 = A;

			REP(N) { a[kk] = _mm256_loadu_pd(A); A += sizeof(__m256) / sizeof(double); }

			REP(N) a[kk] = _mm256_add_pd(a[kk], minf);

			REP(N) a[kk] = _mm256_mul_pd(a[kk], invs);

			REP(N) { _mm256_storeu_pd(A2, a[kk]); A2 += sizeof(__m256) / sizeof(double); }
		}
	}

	for (; i < n; ++i, ++A)
	{
		volatile double v1 = (double)*A + MIN_FREQ;
		*A = v1 * invsum;
	}
}

TARGETAVX void UnifyAVX(float* A, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;
	double invsum = 1.0 / (SumAVX(A, n) + n * MIN_FREQ);

	if (n >= N * sizeof(__m128) / sizeof(float))
	{
		__m256d a[N], minf = _mm256_set1_pd(MIN_FREQ), invs = _mm256_set1_pd(invsum);

		for (int64 l1 = n - N * sizeof(__m128) / sizeof(float); i <= l1; i += N * sizeof(__m128) / sizeof(float))
		{
			float* A2 = A;

			REP(N) { a[kk] = _mm256_cvtps_pd(_mm_loadu_ps(A)); A += sizeof(__m128) / sizeof(float); }

			REP(N) a[kk] = _mm256_add_pd(a[kk], minf);

			REP(N) a[kk] = _mm256_mul_pd(a[kk], invs);

			REP(N) { _mm_storeu_ps(A2, _mm256_cvtpd_ps(a[kk])); A2 += sizeof(__m128) / sizeof(float); }
		}
	}

	for (; i < n; ++i, ++A)
	{
		volatile double v1 = (double)*A + MIN_FREQ;
		*A = v1 * invsum;
	}
}

TARGETAVX char* StrNextIdxAVX(char* A, char val, int64 rep, int64 n)
{
	A++; n--;
	int64 i = 0;

	if (n >= 32)
	{
		__m256i r = _mm256_setzero_si256(), v = _mm256_set1_epi8(val), o = _mm256_setzero_si256();

		for (int64 l1 = n - 32; i <= l1; )
		{
			__m256i a = _mm256_loadu_si256((__m256i*)A);
			r = _mm256_cmpeq_epi8(a, v);

			if (_mm256_testz_si256(r, r)) { A += 32; i += 32; continue; }

			r = _mm256_sad_epu8(o, _mm256_sub_epi8(o, r)); //4 int 64, each 8 bytes

			for (int64 j = 0; j < 4; ++j, A += 8, i += 8)
			{
				if (rep > (int64)simd_u64(r, j))
					rep -= simd_u64(r, j);
				else for (;; A++)
					if (*A == val && !--rep) return A;
			}
		}
	}

	for (; i < n; ++i, A++)
		if (*A == val && !--rep)
			return A;

	return NULL;
}

TARGETAVX int64 CountCharAVX(char* A, char val, int64 n)
{
	uint64 re = 0;
	int64 i = 0;

	if (n >= 32)
	{
		__m256i r = _mm256_setzero_si256(), v = _mm256_set1_epi8(val), o = _mm256_setzero_si256();

		for (int64 l1 = n - 32; i <= l1; i += 32)
		{
			r = _mm256_add_epi64(r, _mm256_sad_epu8(o, _mm256_sub_epi8(o, _mm256_cmpeq_epi8(_mm256_loadu_si256((__m256i*)A), v))));
			A += 32;
		}

		re = simd_u64(r, 0) + simd_u64(r, 1) + simd_u64(r, 2) + simd_u64(r, 3);
	}

	for (; i < n; ++i, A++)
		if (*A == val) re++;

	return (int64)re;
}

#endif