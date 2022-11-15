/* SSE Instruction Set Functions */

#include "vcfpop.h"

#ifndef __aarch64__

template TARGETSSE void RNGSSE<double>::Poly<32>(__m128d* arr, int n, __m128i* re);
template TARGETSSE void RNGSSE<double>::Poly<64>(__m128d* arr, int n, __m128i* re);
template TARGETSSE void RNGSSE<float >::Poly<32>(__m128*  arr, int n, __m128i* re);
template TARGETSSE void RNGSSE<float >::Poly<64>(__m128*  arr, int n, __m128i* re);

#ifndef _RNGSSE_FP64
/* Initialize rng */
TARGETSSE RNGSSE<double>::RNGSSE()
{

}

/* Initialize rng */
TARGETSSE RNGSSE<double>::RNGSSE(uint64 seed, uint64 salt)
{
	__m128i a[32], s, m;

	REP(32) { a[kk] = _mm_set_epi64x(seed + 1, seed + 0); seed += 2; }

	s = _mm_set1_epi64x(salt);
	m = _mm_set1_epi32(0x5bd1e995);

	REP(32) a[kk] = _mm_xor_si128(a[kk], _mm_slli_epi64(_mm_andnot_si128(a[kk], _mm_set1_epi64x(0xFFFFFFFF)), 32));

	s             = _mm_xor_si128(s    , _mm_slli_epi64(_mm_andnot_si128(s    , _mm_set1_epi64x(0xFFFFFFFF)), 32));

	// uint s = s ^ sizeof(uint);
	s = _mm_xor_si128(s, _mm_set1_epi32(sizeof(uint)));

	// a *= m;
	REP(32) a[kk] = _mm_mullo_epi32(a[kk], m);

	// a ^= a >> 24;
	REP(32) a[kk] = _mm_xor_si128(a[kk], _mm_srli_epi32(a[kk], 24));

	// a *= m;
	REP(32) a[kk] = _mm_mullo_epi32(a[kk], m);

	// s *= m;
	s = _mm_mullo_epi32(s, m);

	// a ^= s;
	REP(32) a[kk] = _mm_xor_si128(a[kk], s);

	// a ^= a >> 13;
	REP(32) a[kk] = _mm_xor_si128(a[kk], _mm_srli_epi32(a[kk], 13));

	// a *= m;
	REP(32) a[kk] = _mm_mullo_epi32(a[kk], m);

	// a ^= a >> 15;
	REP(32) a[kk] = _mm_xor_si128(a[kk], _mm_srli_epi32(a[kk], 15));

	// original
	REP(32) x[kk] = _mm_xor_si128(_mm_set1_epi64x(0x159A55E5075BCD15), a[kk]);

	REP(32) a[kk] = _mm_slli_epi64(a[kk], 6);

	REP(32) y[kk] = _mm_xor_si128(_mm_set1_epi64x(0x054913331F123BB5), a[kk]);
}

/* Draw a uniform distriubted real number */
template<int nbits>
TARGETSSE void RNGSSE<double>::Poly(__m128d* arr, int n, __m128i* re)
{
	__m128d t[32], s[32];
	__m128d one = _mm_set1_pd(1.0);
	__m128i mask1 = _mm_set1_epi64x(0x000FFFFFFFFFFFFF);
	__m128i mask2 = _mm_set1_epi64x(0x3FF0000000000000);
	__m128i* r = (__m128i*)t;

	REP(32) s[kk] = _mm_setzero_pd();

	for (int i32 = 0; i32 < n * 32; i32 += 32)
		REP(32) s[kk] = _mm_add_pd(s[kk], arr[kk + i32]);

	{
		__m128i a[32], b[32];

		REP(32) a[kk] = x[kk];

		REP(32) b[kk] = y[kk];

		REP(32) x[kk] = b[kk];

		REP(32) a[kk] = _mm_xor_si128(a[kk], _mm_slli_epi64(a[kk], 23));

		REP(32) a[kk] = _mm_xor_si128(a[kk], _mm_srli_epi64(a[kk], 18));

		REP(32) a[kk] = _mm_xor_si128(a[kk], b[kk]);

		REP(32) a[kk] = _mm_xor_si128(a[kk], _mm_srli_epi64(b[kk], 5));

		REP(32) y[kk] = a[kk];

		REP(32) r[kk] = _mm_add_epi64(a[kk], b[kk]);
	}

	REP(32) r[kk] = _mm_and_si128(r[kk], mask1);

	REP(32) r[kk] = _mm_or_si128(r[kk], mask2);

	REP(32) t[kk] = _mm_sub_pd(t[kk], one);

	REP(32) t[kk] = _mm_mul_pd(t[kk], s[kk]);

	__m128i midx[32], nidx = _mm_setzero_si128(), ninc = _mm_set1_epi64x(1);
	__m128d f[32], b[32];
	REP(32) midx[kk] = _mm_set1_epi64x(n - 1);
	REP(32) f[kk] = _mm_setzero_pd();

	for (int i32 = 0; i32 < n * 32; i32 += 32)
	{
		REP(32) b[kk] = _mm_cmplt_pd(t[kk], arr[kk + i32]);

		REP(8) t[ 0 + kk] = _mm_sub_pd(t[ 0 + kk], arr[ 0 + kk + i32]);

		REP(32) b[kk] = _mm_andnot_pd(f[kk], b[kk]);
			
		REP(8) t[ 8 + kk] = _mm_sub_pd(t[ 8 + kk], arr[ 8 + kk + i32]);

		REP(32) f[kk] = _mm_or_pd(f[kk], b[kk]);
			
		REP(8) t[16 + kk] = _mm_sub_pd(t[16 + kk], arr[16 + kk + i32]);

		REP(32) midx[kk] = _mm_castpd_si128(_mm_blendv_pd(_mm_castsi128_pd(midx[kk]), _mm_castsi128_pd(nidx), b[kk]));//ok
			
		REP(8) t[24 + kk] = _mm_sub_pd(t[24 + kk], arr[24 + kk + i32]);

		nidx = _mm_add_epi64(nidx, ninc);
	}

	if constexpr (nbits == 32)
	{
		REP(16) re[kk] = _mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(midx[0 + (kk << 1)]), _mm_castsi128_ps(midx[1 + (kk << 1)]), _MM_SHUFFLE(2, 0, 2, 0)));
	}
	else
	{
		REP(32) re[kk] = midx[kk];
	}
}
#endif

#ifndef _RNGSSE_FP32
/* Initialize rng */
TARGETSSE RNGSSE<float>::RNGSSE()
{

}

/* Initialize rng */
TARGETSSE RNGSSE<float>::RNGSSE(uint64 seed, uint64 salt)
{	
	__m128i a[16], s, m;
	REP(16) { a[kk] = _mm_set_epi32(Mix(seed +  3), Mix(seed +  2), Mix(seed + 1), Mix(seed + 0)); seed += 4; }

	s = _mm_set1_epi32(Mix(salt));
	m = _mm_set1_epi32(0x5bd1e995);

	// uint s = s ^ sizeof(uint);
	s = _mm_xor_si128(s, _mm_set1_epi32(sizeof(uint)));

	// a *= m;
	REP(16) a[kk] = _mm_mullo_epi32(a[kk], m);

	// a ^= a >> 24;
	REP(16) a[kk] = _mm_xor_si128(a[kk], _mm_srli_epi32(a[kk], 24));

	// a *= m;
	REP(16) a[kk] = _mm_mullo_epi32(a[kk], m);

	// s *= m;
	s = _mm_mullo_epi32(s, m);

	// a ^= s;
	REP(16) a[kk] = _mm_xor_si128(a[kk], s);

	// a ^= a >> 13;
	REP(16) a[kk] = _mm_xor_si128(a[kk], _mm_srli_epi32(a[kk], 13));

	// a *= m;
	REP(16) a[kk] = _mm_mullo_epi32(a[kk], m);

	// a ^= a >> 15;
	REP(16) a[kk] = _mm_xor_si128(a[kk], _mm_srli_epi32(a[kk], 15));

	// original
	REP(16) x[kk] = _mm_xor_si128(_mm_set1_epi32(0x075BCD15), a[kk]);

	REP(16) a[kk] = _mm_slli_epi32(a[kk], 3);

	REP(16) y[kk] = _mm_xor_si128(_mm_set1_epi32(0x159A55E5), a[kk]);

	REP(16) a[kk] = _mm_slli_epi32(a[kk], 3);

	REP(16) z[kk] = _mm_xor_si128(_mm_set1_epi32(0x1F123BB5), a[kk]);
}

/* Draw a uniform distriubted real number */
template<int nbits>
TARGETSSE void RNGSSE<float>::Poly(__m128* arr, int n, __m128i* re)
{
	__m128 t[16], s[16]; __m128i u[16];
	__m128 one = _mm_set1_ps(1.0f);
	__m128i mask1 = _mm_set1_epi32(0x007FFFFF);
	__m128i mask2 = _mm_set1_epi32(0x3F800000);
	__m128i* r = (__m128i*)t;

	REP(16) s[kk] = _mm_setzero_ps();

	for (int i16 = 0; i16 < n * 16; i16 += 16)
		REP(16) s[kk] = _mm_add_ps(s[kk], arr[kk + i16]);

	{
		//xorshift
		REP(16) u[kk] = _mm_slli_epi32(x[kk], 16);
		REP(16) x[kk] = _mm_xor_si128(x[kk], u[kk]);

		REP(16) u[kk] = _mm_srli_epi32(x[kk], 5);
		REP(16) x[kk] = _mm_xor_si128(x[kk], u[kk]);

		REP(16) u[kk] = _mm_slli_epi32(x[kk], 1);
		REP(16) x[kk] = _mm_xor_si128(x[kk], u[kk]);

		REP(16) u[kk] = x[kk];

		REP(16) x[kk] = y[kk];

		REP(16) y[kk] = z[kk];

		REP(16) z[kk] = _mm_xor_si128(u[kk], x[kk]);

		REP(16) z[kk] = _mm_xor_si128(z[kk], y[kk]);
	}

	REP(16) r[kk] = _mm_and_si128(z[kk], mask1);

	REP(16) r[kk] = _mm_or_si128(r[kk], mask2);

	REP(16) t[kk] = _mm_sub_ps(t[kk], one);

	REP(16) t[kk] = _mm_mul_ps(t[kk], s[kk]);

	__m128i midx[16], nidx = _mm_setzero_si128(), ninc = _mm_set1_epi32(1);
	__m128 f[16], b[16];
	REP(16) midx[kk] = _mm_set1_epi32(n - 1);
	REP(16) f[kk] = _mm_setzero_ps();

	for (int i16 = 0; i16 < n * 16; i16 += 16)
	{
		REP(16) b[kk] = _mm_cmplt_ps(t[kk], arr[kk + i16]);

		REP(8) t[ 0 + kk] = _mm_sub_ps(t[ 0 + kk], arr[ 0 + kk + i16]);

		REP(16) b[kk] = _mm_andnot_ps(f[kk], b[kk]);
		
		REP(8) t[ 8 + kk] = _mm_sub_ps(t[ 8 + kk], arr[ 8 + kk + i16]);

		REP(16)
		{
			f[kk] = _mm_or_ps(f[kk], b[kk]);
			midx[kk] = _mm_castps_si128(_mm_blendv_ps(_mm_castsi128_ps(midx[kk]), _mm_castsi128_ps(nidx), b[kk]));//ok
		}

		nidx = _mm_add_epi32(nidx, ninc);
	}

	if constexpr (nbits == 32)
	{
		REP(16) re[kk] = midx[kk];
	}
	else
	{
		REP(16)
		{
			re[0 + (kk << 1)] = _mm_cvtepi32_epi64(midx[kk]);
			re[1 + (kk << 1)] = _mm_cvtepi32_epi64(_mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(midx[kk]), _mm_castsi128_ps(midx[kk]), _MM_SHUFFLE(1, 0, 3, 2))));
		}
	}
}
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

TARGETSSE int64 GetMinIdxSSE(double* A, int64 n, double& val)
{
	constexpr int N = 4;
	int64 i = 0;
	val = DBL_MAX;
	uint64 idx = (uint64)-1;

	if (n >= N * sizeof(__m128d) / sizeof(double))
	{
		__m128d min1[N], f[N], a[N];
		__m128i midx[N], nidx[N], msep = _mm_set1_epi64x(N * sizeof(__m128d) / sizeof(double));
		REP(N) min1[kk] = _mm_set1_pd(val);
		REP(N) midx[kk] = _mm_set1_epi64x(0xFFFFFFFFFFFFFFFF);
		REP(N) nidx[kk] = _mm_set_epi64x(1 + (kk << 1), 0 + (kk << 1));

		for (int64 l1 = n - N * sizeof(__m128d) / sizeof(double); i <= l1; i += N * sizeof(__m128d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm_loadu_pd(A); A += sizeof(__m128d) / sizeof(double); }

			REP(N) f[kk] = _mm_cmpgt_pd(min1[kk], a[kk]);

			REP(N) min1[kk] = _mm_blendv_pd(min1[kk], a[kk], f[kk]);//ok

			REP(N) midx[kk] = _mm_castpd_si128(_mm_blendv_pd(_mm_castsi128_pd(midx[kk]), _mm_castsi128_pd(nidx[kk]), f[kk]));//ok

			REP(N) nidx[kk] = _mm_add_epi64(nidx[kk], msep);
		}

		for (int K = sizeof(min1) / sizeof(min1[0]) / 2; K >= 1; K >>= 1)
		{
			REP(K) f[kk] = _mm_cmpgt_pd(min1[kk], min1[kk + K]);
			REP(K) min1[kk] = _mm_blendv_pd(min1[kk], min1[kk + K], f[kk]);
			REP(K) midx[kk] = _mm_castpd_si128(_mm_blendv_pd(_mm_castsi128_pd(midx[kk]), _mm_castsi128_pd(midx[kk + K]), f[kk]));
		}

		for (int64 j = 0; j < sizeof(__m128d) / sizeof(double); ++j)
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

TARGETSSE int64 GetMinIdxSSE(float* A, int64 n, float& val)
{
	constexpr int N = 4;
	int64 i = 0;
	val = FLT_MAX;
	uint idx = (uint)-1;

	if (n >= N * sizeof(__m128) / sizeof(float))
	{
		__m128 min1[N], f[N], a[N];
		__m128i midx[N], nidx[N], msep = _mm_set1_epi32(N * sizeof(__m128) / sizeof(float));
		REP(N) min1[kk] = _mm_set1_ps(val);
		REP(N) midx[kk] = _mm_set1_epi32(0xFFFFFFFF);
		REP(N) nidx[kk] = _mm_set_epi32(3 + (kk << 2), 2 + (kk << 2), 1 + (kk << 2), 0 + (kk << 2));

		for (int64 l1 = n - N * sizeof(__m128) / sizeof(float); i <= l1; i += N * sizeof(__m128) / sizeof(float))
		{
			REP(N) { a[kk] = _mm_loadu_ps(A); A += sizeof(__m128) / sizeof(float); }

			REP(N) f[kk] = _mm_cmpgt_ps(min1[kk], a[kk]);

			REP(N) min1[kk] = _mm_blendv_ps(min1[kk], a[kk], f[kk]);//ok

			REP(N) midx[kk] = _mm_castps_si128(_mm_blendv_ps(_mm_castsi128_ps(midx[kk]), _mm_castsi128_ps(nidx[kk]), f[kk]));

			REP(N) nidx[kk] = _mm_add_epi32(nidx[kk], msep);
		}

		for (int K = sizeof(min1) / sizeof(min1[0]) / 2; K >= 1; K >>= 1)
		{
			REP(K) f[kk] = _mm_cmpgt_ps(min1[kk], min1[kk + K]);
			REP(K) min1[kk] = _mm_blendv_ps(min1[kk], min1[kk + K], f[kk]);
			REP(K) midx[kk] = _mm_castps_si128(_mm_blendv_ps(_mm_castsi128_ps(midx[kk]), _mm_castsi128_ps(midx[kk + K]), f[kk]));
		}

		for (int64 j = 0; j < sizeof(__m128) / sizeof(float); ++j)
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

TARGETSSE void GetMinMaxValSSE(double* A, int64 n, double& minv, double& maxv)
{
	constexpr int N = 4;
	int64 i = 0;
	minv = DBL_MAX;
	maxv = -DBL_MAX;

	if (n >= N * sizeof(__m128d) / sizeof(double))
	{
		__m128d min1[N], max1[N], a[N];
		REP(N) min1[kk] = _mm_set1_pd(minv);
		REP(N) max1[kk] = _mm_set1_pd(maxv);

		for (int64 l1 = n - N * sizeof(__m128d) / sizeof(double); i <= l1; i += N * sizeof(__m128d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm_loadu_pd(A); A += sizeof(__m128d) / sizeof(double); }

			REP(N)
			{
				min1[kk] = _mm_min_pd(min1[kk], a[kk]);
				max1[kk] = _mm_max_pd(max1[kk], a[kk]);
			}
		}

		for (int K = sizeof(max1) / sizeof(max1[0]) / 2; K >= 1; K >>= 1)
		{
			REP(K) min1[kk] = _mm_min_pd(min1[kk], min1[kk + K]);
			REP(K) max1[kk] = _mm_max_pd(max1[kk], max1[kk + K]);
		}

		minv = Min(simp_f64(min1, 0), simp_f64(min1, 1));
		maxv = Max(simp_f64(max1, 0), simp_f64(max1, 1));
	}

	for (; i < n; ++i, ++A)
	{
		if (*A < minv) minv = *A;
		if (*A > maxv) maxv = *A;
	}
}

TARGETSSE void GetMinMaxValSSE(float* A, int64 n, float& minv, float& maxv)
{
	constexpr int N = 4;
	int64 i = 0;
	minv = FLT_MAX;
	maxv = -FLT_MAX;

	if (n >= N * sizeof(__m128) / sizeof(float))
	{
		__m128 min1[N], max1[N], a[N];
		REP(N) min1[kk] = _mm_set1_ps(minv);
		REP(N) max1[kk] = _mm_set1_ps(maxv);

		for (int64 l1 = n - N * sizeof(__m128) / sizeof(float); i <= l1; i += N * sizeof(__m128) / sizeof(float))
		{
			REP(N) { a[kk] = _mm_loadu_ps(A); A += sizeof(__m128) / sizeof(float); }

			REP(N)
			{
				min1[kk] = _mm_min_ps(min1[kk], a[kk]);
				max1[kk] = _mm_max_ps(max1[kk], a[kk]);
			}
		}

		for (int K = sizeof(max1) / sizeof(max1[0]) / 2; K >= 1; K >>= 1)
		{
			REP(K) min1[kk] = _mm_min_ps(min1[kk], min1[kk + K]);
			REP(K) max1[kk] = _mm_max_ps(max1[kk], max1[kk + K]);
		}

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

TARGETSSE double GetMaxValSSE(double* A, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;
	double val = -DBL_MAX;

	if (n >= N * sizeof(__m128d) / sizeof(double))
	{
		__m128d max1[N];
		REP(N) max1[kk] = _mm_set1_pd(val);

		for (int64 l1 = n - N * sizeof(__m128d) / sizeof(double); i <= l1; i += N * sizeof(__m128d) / sizeof(double))
		{
			REP(N)
			{
				max1[kk] = _mm_max_pd(max1[kk], _mm_loadu_pd(A));
				A += sizeof(__m128d) / sizeof(double);
			}
		}

		for (int K = sizeof(max1) / sizeof(max1[0]) / 2; K >= 1; K >>= 1)
			REP(K) max1[kk] = _mm_max_pd(max1[kk], max1[kk + K]);

		val = Max(simp_f64(max1, 0), simp_f64(max1, 1));
	}

	for (; i < n; ++i, ++A)
	{
		if (*A < val) continue;
		val = *A;
	}

	return val;
}

TARGETSSE float GetMaxValSSE(float* A, int64 n)
{
	constexpr int N = 8;
	int64 i = 0;
	float val = -FLT_MAX;

	if (n >= N * sizeof(__m128) / sizeof(float))
	{
		__m128 max1[N], a[N];
		REP(N) max1[kk] = _mm_set1_ps(val);

		for (int64 l1 = n - N * sizeof(__m128) / sizeof(float); i <= l1; i += N * sizeof(__m128) / sizeof(float))
		{
			REP(N) { a[kk] = _mm_loadu_ps(A); A += sizeof(__m128) / sizeof(float); }

			REP(N) { max1[kk] = _mm_max_ps(max1[kk], a[kk]); }
		}

		for (int K = sizeof(max1) / sizeof(max1[0]) / 2; K >= 1; K >>= 1)
			REP(K) max1[kk] = _mm_max_ps(max1[kk], max1[kk + K]);

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

TARGETSSE double GetMaxValSSE(double* A, int64 n, int64 sep)
{
	int64 i = 0;
	double val = -DBL_MAX;

	if (n >= 8)
	{
		__m128d max1[2];
		REP(2) max1[kk] = _mm_set1_pd(val);
		REP(4) { _mm_prefetch((const char*)A, _MM_HINT_T0); A += sep; }

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			REP(4) { _mm_prefetch((const char*)A, _MM_HINT_T0); A += sep; }
			max1[0] = _mm_max_pd(max1[0], _mm_set_pd(A[-7 * sep], A[-8 * sep]));
			max1[1] = _mm_max_pd(max1[1], _mm_set_pd(A[-5 * sep], A[-6 * sep]));

			REP(4) { _mm_prefetch((const char*)A, _MM_HINT_T0); A += sep; }
			max1[0] = _mm_max_pd(max1[0], _mm_set_pd(A[-7 * sep], A[-8 * sep]));
			max1[1] = _mm_max_pd(max1[1], _mm_set_pd(A[-5 * sep], A[-6 * sep]));
		}

		max1[0] = _mm_max_pd(max1[0], max1[1]);

		val = Max(simp_f64(max1, 0), simp_f64(max1, 1));

		A -= 4 * sep;
	}

	for (; i < n; ++i, A += sep)
	{
		if (*A < val) continue;
		val = *A;
	}

	return val;
}

TARGETSSE float GetMaxValSSE(float* A, int64 n, int64 sep)
{
	constexpr int N = 4;
	int64 i = 0;
	float val = -FLT_MAX;

	if (n >= N * sizeof(__m128) / sizeof(float))
	{
		__m128 max1[N];
		REP(N) max1[kk] = _mm_set1_ps(val);

		for (int64 l1 = n - N * sizeof(__m128) / sizeof(float); i <= l1; i += N * sizeof(__m128) / sizeof(float))
		{
			REP(N)
			{
				max1[kk] = _mm_max_ps(max1[kk], _mm_set_ps(A[3 * sep], A[2 * sep], A[1 * sep], A[0 * sep]));
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

TARGETSSE double GetMinValSSE(double* A, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;
	double val = DBL_MAX;

	if (n >= N * sizeof(__m128d) / sizeof(double))
	{
		__m128d min1[N];
		REP(N) min1[kk] = _mm_set1_pd(val);

		for (int64 l1 = n - N * sizeof(__m128d) / sizeof(double); i <= l1; i += N * sizeof(__m128d) / sizeof(double))
		{
			REP(N)
			{
				min1[kk] = _mm_min_pd(min1[kk], _mm_loadu_pd(A));
				A += sizeof(__m128d) / sizeof(double);
			}
		}

		for (int K = sizeof(min1) / sizeof(min1[0]) / 2; K >= 1; K >>= 1)
			REP(K) min1[kk] = _mm_min_pd(min1[kk], min1[kk + K]);

		val = Min(simp_f64(min1, 0), simp_f64(min1, 1));
	}

	for (; i < n; ++i, ++A)
	{
		if (*A > val) continue;
		val = *A;
	}

	return val;
}

TARGETSSE float GetMinValSSE(float* A, int64 n)
{
	constexpr int N = 8;
	int64 i = 0;
	float val = FLT_MAX;

	if (n >= N * sizeof(__m128) / sizeof(float))
	{
		__m128 min1[N], a[N];
		REP(N) min1[kk] = _mm_set1_ps(val);

		for (int64 l1 = n - N * sizeof(__m128) / sizeof(float); i <= l1; i += N * sizeof(__m128) / sizeof(float))
		{
			REP(N) { a[kk] = _mm_loadu_ps(A); A += sizeof(__m128) / sizeof(float); }

			REP(N) { min1[kk] = _mm_min_ps(min1[kk], a[kk]); }
		}

		for (int K = sizeof(min1) / sizeof(min1[0]) / 2; K >= 1; K >>= 1)
			REP(K) min1[kk] = _mm_min_ps(min1[kk], min1[kk + K]);

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

TARGETSSE int64 GetMinValSSE(int64* A, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;
	int64 val = 0x7FFFFFFFFFFFFFFF;

	if (n >= N * sizeof(__m128i) / sizeof(int64))
	{
		__m128i min1[N], a[N], f[N];
		REP(N) min1[kk] = _mm_set1_epi64x(0x7FFFFFFFFFFFFFFF);

		for (int64 l1 = n - N * sizeof(__m128i) / sizeof(int64); i <= l1; i += N * sizeof(__m128i) / sizeof(int64))
		{
			REP(N) { a[kk] = _mm_loadu_si128((__m128i*)A); A += sizeof(__m128i) / sizeof(int64); }
			
			REP(N) f[kk] = _mm_cmpgt_epi64(min1[kk], a[kk]);

			REP(N) min1[kk] = _mm_castpd_si128(_mm_blendv_pd( _mm_castsi128_pd(min1[kk]), _mm_castsi128_pd(a[kk]), _mm_castsi128_pd(f[kk])));
		}

		for (int K = sizeof(min1) / sizeof(min1[0]) / 2; K >= 1; K >>= 1)
			REP(K) min1[kk] = _mm_castpd_si128(_mm_blendv_pd( _mm_castsi128_pd(min1[kk]), _mm_castsi128_pd(min1[kk + K]), _mm_castsi128_pd(_mm_cmpgt_epi64(min1[kk], min1[kk + K]))));

		val = Min(simp_i64(min1, 0), simp_i64(min1, 1));
	}

	for (; i < n; ++i, ++A)
	{
		if (*A > val) continue;
		val = *A;
	}

	return val;
}

TARGETSSE void SetValSSE(uint* A, ushort* B, int64 n)
{
	int64 i = 0;

	if (n >= 8)
	{
		__m128i b;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			b = _mm_loadu_si128((__m128i*)B);  B += 8;

			_mm_storeu_si128((__m128i*)A, _mm_cvtepu16_epi32(b)); A += 4;

			b = _mm_srli_si128(b, 8);
			_mm_storeu_si128((__m128i*)A, _mm_cvtepu16_epi32(b)); A += 4;
		}
	}

	for (; i < n; ++i)
		*A++ = *B++;
}

TARGETSSE void AddExponentSSE(int64& slog, __m128d& val)
{
	__m128i& vv = *(__m128i*)&val;
	__m128i mask1 = _mm_set1_epi64x(0x7FF0000000000000);
	__m128i mask2 = _mm_set1_epi64x(0x800FFFFFFFFFFFFF);
	__m128i mask3 = _mm_set1_epi64x(0x3FF0000000000000);
	__m128i subv = _mm_set1_epi64x(1023);

	__m128i t = _mm_sub_epi64(_mm_srli_epi64(_mm_and_si128(vv, mask1), 52), subv);

	slog += (int)simd_u64(t, 0) + (int)simd_u64(t, 1);

	vv = _mm_or_si128(_mm_and_si128(vv, mask2), mask3);
}

TARGETSSE void AddExponentSSE(int64& slog, __m128& val)
{
	__m128i& vv = *(__m128i*) & val;
	__m128i mask1 = _mm_set1_epi32(0x7F800000);
	__m128i mask2 = _mm_set1_epi32(0x807FFFFF);
	__m128i mask3 = _mm_set1_epi32(0x3F800000);
	__m128i subv = _mm_set1_epi32(127);

	__m128i t = _mm_sub_epi32(_mm_srli_epi32(_mm_and_si128(vv, mask1), 23), subv);

	slog += (int)simd_u32(t, 0) + (int)simd_u32(t, 1) + (int)simd_u32(t, 2) + (int)simd_u32(t, 3);

	vv = _mm_or_si128(_mm_and_si128(vv, mask2), mask3);
}

TARGETSSE void ChargeLogSSE(int64& slog, double& prod, __m128d& val)
{
	AddExponentSSE(slog, val);

	prod *= _mm_reduce_mul_pd(val);

	if (prod < DOUBLE_UNDERFLOW || prod > DOUBLE_OVERFLOW) [[unlikely]]
		AddExponent(slog, prod);
}

TARGETSSE void ChargeLogSSE(int64& slog, double& prod, __m128& val)
{
	AddExponentSSE(slog, val);

	prod *= _mm_reduce_mul_psd(val);

	if (prod < DOUBLE_UNDERFLOW || prod > DOUBLE_OVERFLOW) [[unlikely]]
		AddExponent(slog, prod);
}

TARGETSSE double LogProdSSE(double* A, int64 n)
{
	int64 i = 0;
	int64 slog = 0; double prod = 1;

	if (n >= 32)
	{
		__m128d pd[4];
		REP(4) pd[kk] = _mm_set1_pd(1.0);

		for (int64 l1 = n - 32; i <= l1; i += 32)
		{
			REP(2)
			{
				pd[0] = _mm_mul_pd(pd[0], _mm_mul_pd(_mm_loadu_pd(A + 0), _mm_loadu_pd(A + 8)));
				//f1 = _mm_cmplt_pd(pd[0], dunder);
				//f2 = _mm_cmplt_pd(dover, pd[0]);

				pd[1] = _mm_mul_pd(pd[1], _mm_mul_pd(_mm_loadu_pd(A + 2), _mm_loadu_pd(A + 10)));
				//f1 = _mm_or_pd(f1, _mm_cmplt_pd(pd[1], dunder));
				//f2 = _mm_or_pd(f2, _mm_cmplt_pd(dover, pd[1]));

				pd[2] = _mm_mul_pd(pd[2], _mm_mul_pd(_mm_loadu_pd(A + 4), _mm_loadu_pd(A + 12)));
				//f1 = _mm_or_pd(f1, _mm_cmplt_pd(pd[2], dunder));
				//f2 = _mm_or_pd(f2, _mm_cmplt_pd(dover, pd[2]));

				pd[3] = _mm_mul_pd(pd[3], _mm_mul_pd(_mm_loadu_pd(A + 6), _mm_loadu_pd(A + 14)));
				//f1 = _mm_or_pd(f1, _mm_cmplt_pd(pd[3], dunder));
				//f2 = _mm_or_pd(f2, _mm_cmplt_pd(dover, pd[3]));

				A += 16;
			}

			//f1 = _mm_or_pd(f1, f2);
			//if (_mm_testz_pd(f1, f1)) [[likely]] continue;

			REP(4) AddExponentSSE(slog, pd[kk]);
		}

		REP(4) ChargeLogSSE(slog, prod, pd[kk]);
	}

	for (; i < n; ++i, ++A)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);
	return prod;
}

TARGETSSE double LogProdSSE(float* A, int64 n)
{
	int64 i = 0;
	int64 slog = 0; double prod = 1;

	if (n >= 64)
	{
		__m128d pd[4];
		REP(4) pd[kk] = _mm_set1_pd(1.0);

		for (int64 l1 = n - 64; i <= l1; i += 64)
		{
			REP(4)
			{
				pd[0] = _mm_mul_pd(pd[0], _mm_mul_pd(
							_mm_cvtps_pd(_mm_loadu_ps(A + 0)),
							_mm_cvtps_pd(_mm_loadu_ps(A + 8))));
				pd[1] = _mm_mul_pd(pd[1], _mm_mul_pd(
							_mm_cvtps_pd(_mm_loadu_ps(A + 2)),
							_mm_cvtps_pd(_mm_loadu_ps(A + 10))));
				pd[2] = _mm_mul_pd(pd[2], _mm_mul_pd(
							_mm_cvtps_pd(_mm_loadu_ps(A + 4)),
							_mm_cvtps_pd(_mm_loadu_ps(A + 12))));
				pd[3] = _mm_mul_pd(pd[3], _mm_mul_pd(
							_mm_cvtps_pd(_mm_loadu_ps(A + 6)),
							_mm_cvtps_pd(_mm_loadu_ps(A + 14))));
				A += 16;
			}

			REP(4) AddExponentSSE(slog, pd[kk]);
		}

		REP(4) ChargeLogSSE(slog, prod, pd[kk]);
	}

	for (; i < n; ++i, ++A)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);
	return prod;
}

TARGETSSE float LogProdSSEx(float* A, int64 n)
{
	int64 i = 0;
	int64 slog = 0; double prod = 1;

	if (n >= 32)
	{
		__m128 pd[4];
		REP(4) pd[kk] = _mm_set1_ps(1.0f);

		for (int64 l1 = n - 32; i <= l1; i += 32)
		{
			pd[0] = _mm_mul_ps(pd[0], _mm_loadu_ps(A));
			pd[1] = _mm_mul_ps(pd[1], _mm_loadu_ps(A + 4));
			pd[2] = _mm_mul_ps(pd[2], _mm_loadu_ps(A + 8));
			pd[3] = _mm_mul_ps(pd[3], _mm_loadu_ps(A + 12));
			pd[0] = _mm_mul_ps(pd[0], _mm_loadu_ps(A + 16));
			pd[1] = _mm_mul_ps(pd[1], _mm_loadu_ps(A + 20));
			pd[2] = _mm_mul_ps(pd[2], _mm_loadu_ps(A + 24));
			pd[3] = _mm_mul_ps(pd[3], _mm_loadu_ps(A + 28));

			A += 32;
			REP(4) AddExponentSSE(slog, pd[kk]);
		}

		REP(4) ChargeLogSSE(slog, prod, pd[kk]);
	}

	for (; i < n; ++i, ++A)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);
	return prod;
}

TARGETSSE double LogProdSSE(double* A, int64 n, int64 sep)
{
	int64 i = 0;
	int64 slog = 0; double prod = 1;

	if (n >= 8)
	{
		REP(4) { _mm_prefetch((const char*)A, _MM_HINT_T0); A += sep; }
		__m128d f1, f2, pd[4], v1[4], dunder = _mm_set1_pd(DOUBLE_UNDERFLOW), dover = _mm_set1_pd(DOUBLE_OVERFLOW);
		REP(4) pd[kk] = _mm_set1_pd(1.0);

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			REP(4)
			{
				_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
				_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
				v1[kk] = _mm_set_pd(A[-5 * sep], A[-6 * sep]);
			}

			f1 = _mm_setzero_pd();
			f2 = _mm_setzero_pd();

			REP(4)
			{
				pd[kk] = _mm_mul_pd(pd[kk], v1[kk]);

				f1 = _mm_or_pd(f1, _mm_cmplt_pd(pd[kk], dunder));
				f2 = _mm_or_pd(f2, _mm_cmplt_pd(dover, pd[kk]));
			}

			f1 = _mm_or_pd(f1, f2);

			if (_mm_testz_si128(_mm_castpd_si128(f1), _mm_castpd_si128(f1))) [[likely]] continue;

			REP(4) AddExponentSSE(slog, pd[kk]);
		}

		REP(4) ChargeLogSSE(slog, prod, pd[kk]);
		A -= 4 * sep;
	}

	for (; i < n; ++i, A += sep)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);
	return prod;
}

TARGETSSE double LogProdSSE(float* A, int64 n, int64 sep)
{
	int64 i = 0;
	int64 slog = 0; double prod = 1;

	if (n >= 32)
	{
		REP(4) { _mm_prefetch((const char*)A, _MM_HINT_T0); A += sep; }

		__m128d pd[2];
		__m128 a1, a2;
		REP(2) pd[kk] = _mm_set1_pd(1.0);

		for (int64 l1 = n - 32; i <= l1; i += 32)
		{
			REP(8)
			{
				_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
				_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
				_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
				_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;

				a1 = _mm_set_ps(A[-5 * sep], A[-6 * sep], A[-7 * sep], A[-8 * sep]);
				a2 = _mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(a1), 8));

				pd[0] = _mm_mul_pd(pd[0], _mm_cvtps_pd(a1));
				pd[1] = _mm_mul_pd(pd[1], _mm_cvtps_pd(a2));
			}

			REP(2) AddExponentSSE(slog, pd[kk]);
		}

		REP(2) ChargeLogSSE(slog, prod, pd[kk]);
		A -= 4 * sep;
	}

	for (; i < n; ++i, A += sep)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);
	return prod;
}

TARGETSSE float LogProdSSEx(float* A, int64 n, int64 sep)
{
	int64 i = 0;
	int64 slog = 0; double prod = 1;

	if (n >= 32)
	{
		__m128 pd[4];
		REP(4) pd[kk] = _mm_set1_ps(1.0f);

		for (int64 l1 = n - 32; i <= l1; i += 32)
		{
			pd[0] = _mm_mul_ps(pd[0], _mm_set_ps(A[3 * sep], A[2 * sep], A[1 * sep], A[0 * sep])); A += 4 * sep;
			pd[1] = _mm_mul_ps(pd[1], _mm_set_ps(A[3 * sep], A[2 * sep], A[1 * sep], A[0 * sep])); A += 4 * sep;
			pd[2] = _mm_mul_ps(pd[2], _mm_set_ps(A[3 * sep], A[2 * sep], A[1 * sep], A[0 * sep])); A += 4 * sep;
			pd[3] = _mm_mul_ps(pd[3], _mm_set_ps(A[3 * sep], A[2 * sep], A[1 * sep], A[0 * sep])); A += 4 * sep;

			pd[0] = _mm_mul_ps(pd[0], _mm_set_ps(A[3 * sep], A[2 * sep], A[1 * sep], A[0 * sep])); A += 4 * sep;
			pd[1] = _mm_mul_ps(pd[1], _mm_set_ps(A[3 * sep], A[2 * sep], A[1 * sep], A[0 * sep])); A += 4 * sep;
			pd[2] = _mm_mul_ps(pd[2], _mm_set_ps(A[3 * sep], A[2 * sep], A[1 * sep], A[0 * sep])); A += 4 * sep;
			pd[3] = _mm_mul_ps(pd[3], _mm_set_ps(A[3 * sep], A[2 * sep], A[1 * sep], A[0 * sep])); A += 4 * sep;

			REP(4) AddExponentSSE(slog, pd[kk]);
		}

		REP(4) ChargeLogSSE(slog, prod, pd[kk]);
	}

	for (; i < n; ++i, A += sep)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);
	return prod;
}

TARGETSSE double LogProdDivSSE(double* A, double* B, int64 n, int64 sep)
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
		__m128d pd1[4], pd2[4];
		REP(4) pd1[kk] = pd2[kk] = _mm_set1_pd(1.0);

		for (int64 l1 = n - 64; i <= l1; i += 64)
		{
			REP(8)
			{
				pd1[0] = _mm_mul_pd(pd1[0], _mm_set_pd(A[1 * sep], A[0 * sep])); A += sep * 2;
				pd2[0] = _mm_mul_pd(pd2[0], _mm_set_pd(B[1 * sep], B[0 * sep])); B += sep * 2;

				pd1[1] = _mm_mul_pd(pd1[1], _mm_set_pd(A[1 * sep], A[0 * sep])); A += sep * 2;
				pd2[1] = _mm_mul_pd(pd2[1], _mm_set_pd(B[1 * sep], B[0 * sep])); B += sep * 2;

				pd1[2] = _mm_mul_pd(pd1[2], _mm_set_pd(A[1 * sep], A[0 * sep])); A += sep * 2;
				pd2[2] = _mm_mul_pd(pd2[2], _mm_set_pd(B[1 * sep], B[0 * sep])); B += sep * 2;

				pd1[3] = _mm_mul_pd(pd1[3], _mm_set_pd(A[1 * sep], A[0 * sep])); A += sep * 2;
				pd2[3] = _mm_mul_pd(pd2[3], _mm_set_pd(B[1 * sep], B[0 * sep])); B += sep * 2;
			}

			REP(4) AddExponentSSE(slog1, pd1[kk]);
			REP(4) AddExponentSSE(slog2, pd2[kk]);
		}

		REP(4) ChargeLogSSE(slog1, prod1, pd1[kk]);
		REP(4) ChargeLogSSE(slog2, prod2, pd2[kk]);
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

TARGETSSE double LogProdDivSSE(float* A, float* B, int64 n, int64 sep)
{
	int64 i = 0;
	int64 slog1 = 0; double prod1 = 1;
	int64 slog2 = 0; double prod2 = 1;

	if (n >= 64)
	{
		__m128d pd1[4], pd2[4];
		REP(4) pd1[kk] = pd2[kk] = _mm_set1_pd(1.0);

		for (int64 l1 = n - 64; i <= l1; i += 64)
		{
			REP(8)
			{
				pd1[0] = _mm_mul_pd(pd1[0], _mm_set_pd(A[1 * sep], A[0 * sep])); A += sep * 2;
				pd2[0] = _mm_mul_pd(pd2[0], _mm_set_pd(B[1 * sep], B[0 * sep])); B += sep * 2;

				pd1[1] = _mm_mul_pd(pd1[1], _mm_set_pd(A[1 * sep], A[0 * sep])); A += sep * 2;
				pd2[1] = _mm_mul_pd(pd2[1], _mm_set_pd(B[1 * sep], B[0 * sep])); B += sep * 2;

				pd1[2] = _mm_mul_pd(pd1[2], _mm_set_pd(A[1 * sep], A[0 * sep])); A += sep * 2;
				pd2[2] = _mm_mul_pd(pd2[2], _mm_set_pd(B[1 * sep], B[0 * sep])); B += sep * 2;

				pd1[3] = _mm_mul_pd(pd1[3], _mm_set_pd(A[1 * sep], A[0 * sep])); A += sep * 2;
				pd2[3] = _mm_mul_pd(pd2[3], _mm_set_pd(B[1 * sep], B[0 * sep])); B += sep * 2;
			}

			REP(4) AddExponentSSE(slog1, pd1[kk]);
			REP(4) AddExponentSSE(slog2, pd2[kk]);
		}

		REP(4) ChargeLogSSE(slog1, prod1, pd1[kk]);
		REP(4) ChargeLogSSE(slog2, prod2, pd2[kk]);
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

TARGETSSE float LogProdDivSSEx(float* A, float* B, int64 n, int64 sep)
{
	int64 i = 0;
	int64 slog1 = 0; double prod1 = 1;
	int64 slog2 = 0; double prod2 = 1;

	if (n >= 64)
	{
		__m128 pd1[4], pd2[4];
		REP(4) pd1[kk] = pd2[kk] = _mm_set1_ps(1.0f);

		for (int64 l1 = n - 64; i <= l1; i += 64)
		{
			REP(4)
			{
				pd1[0] = _mm_mul_ps(pd1[0], _mm_set_ps(A[3 * sep], A[2 * sep], A[1 * sep], A[0 * sep])); A += 4 * sep;
				pd2[0] = _mm_mul_ps(pd2[0], _mm_set_ps(B[3 * sep], B[2 * sep], B[1 * sep], B[0 * sep])); B += 4 * sep;

				pd1[1] = _mm_mul_ps(pd1[1], _mm_set_ps(A[3 * sep], A[2 * sep], A[1 * sep], A[0 * sep])); A += 4 * sep;
				pd2[1] = _mm_mul_ps(pd2[1], _mm_set_ps(B[3 * sep], B[2 * sep], B[1 * sep], B[0 * sep])); B += 4 * sep;

				pd1[2] = _mm_mul_ps(pd1[2], _mm_set_ps(A[3 * sep], A[2 * sep], A[1 * sep], A[0 * sep])); A += 4 * sep;
				pd2[2] = _mm_mul_ps(pd2[2], _mm_set_ps(B[3 * sep], B[2 * sep], B[1 * sep], B[0 * sep])); B += 4 * sep;

				pd1[3] = _mm_mul_ps(pd1[3], _mm_set_ps(A[3 * sep], A[2 * sep], A[1 * sep], A[0 * sep])); A += 4 * sep;
				pd2[3] = _mm_mul_ps(pd2[3], _mm_set_ps(B[3 * sep], B[2 * sep], B[1 * sep], B[0 * sep])); B += 4 * sep;
			}

			REP(4) AddExponentSSE(slog1, pd1[kk]);
			REP(4) AddExponentSSE(slog2, pd2[kk]);
		}

		REP(4) ChargeLogSSE(slog1, prod1, pd1[kk]);
		REP(4) ChargeLogSSE(slog2, prod2, pd2[kk]);
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

TARGETSSE int64 CountNonZeroSSE(byte* A, int64 n)
{
	uint64 re = 0;
	int64 i = 0;

	if (n >= 16)
	{
		__m128i a = _mm_setzero_si128(), z = _mm_setzero_si128();

		for (int64 l1 = n - 16; i <= l1; i += 16)
		{
			a = _mm_add_epi64(a, _mm_sad_epu8(z, _mm_sub_epi8(z, _mm_cmpeq_epi8(_mm_loadu_si128((__m128i*)A), z))));
			A += 16;
		}

		re = i - (simd_u64(a, 0) + simd_u64(a, 1));
	}

	for (; i < n; ++i, ++A)
		if (*A) re++;

	return (int64)re;
}

TARGETSSE double SumSSE(double* A, int64 n)
{
	constexpr int N = 8;
	int64 i = 0;
	double re = 0;

	if (n >= N * sizeof(__m128d) / sizeof(double))
	{
		__m128d s[N], a[N];
		REP(N) s[kk] = _mm_setzero_pd();

		for (int64 l1 = n - N * sizeof(__m128d) / sizeof(double); i <= l1; i += N * sizeof(__m128d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm_loadu_pd(A); A += sizeof(__m128d) / sizeof(double); }

			REP(N) s[kk] = _mm_add_pd(s[kk], a[kk]);
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm_add_pd(s[kk], s[kk + K]);

		re = _mm_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *A++;
		re += v1;
	}

	return re;
}

TARGETSSE double SumSSE(float* A, int64 n)
{
	constexpr int N = 16;
	int64 i = 0;
	double re = 0;

	if (n >= N / 2 * sizeof(__m128) / sizeof(float))
	{
		__m128d s[N];
		REP(N) s[kk] = _mm_setzero_pd();

		for (int64 l1 = n - N / 2 * sizeof(__m128) / sizeof(float); i <= l1; i += N / 2 * sizeof(__m128) / sizeof(float))
		{
			REP(N / 2)
			{
				__m128 v1 = _mm_loadu_ps(A); A += sizeof(__m128) / sizeof(float);
				__m128d a1 = _mm_cvtps_pd(v1);
				__m128d a2 = _mm_cvtps_pd(_mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(v1), 8)));

				s[0 + (kk << 1)] = _mm_add_pd(s[0 + (kk << 1)], a1);
				s[1 + (kk << 1)] = _mm_add_pd(s[1 + (kk << 1)], a2);
			}
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm_add_pd(s[kk], s[kk + K]);

		re = _mm_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A++;
		re += v1;
	}

	return re;
}

TARGETSSE float SumSSEx(float* A, int64 n)
{
	constexpr int N = 8;
	int64 i = 0;
	float re = 0;

	if (n >= N * sizeof(__m128) / sizeof(float))
	{
		__m128 s[N], a[N];
		REP(N) s[kk] = _mm_setzero_ps();

		for (int64 l1 = n - N * sizeof(__m128) / sizeof(float); i <= l1; i += N * sizeof(__m128) / sizeof(float))
		{
			REP(N) { a[kk] = _mm_loadu_ps(A); A += sizeof(__m128) / sizeof(float); }

			REP(N) s[kk] = _mm_add_ps(s[kk], a[kk]);
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm_add_ps(s[kk], s[kk + K]);

		re = _mm_reduce_add_ps(s[0]);
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A++;
		re += v1;
	}

	return re;
}

TARGETSSE int64 SumSSE(byte* A, int64 n)
{
	constexpr int N = 4;
	uint64 re = 0;
	int64 i = 0;

	if (n >= N * sizeof(__m128i) / sizeof(byte))
	{
		__m128i s[N], a[N], z = _mm_setzero_si128();
		REP(N) s[kk] = _mm_setzero_si128();

		for (int64 l1 = n - N * sizeof(__m128i) / sizeof(byte); i <= l1; i += N * sizeof(__m128i) / sizeof(byte))
		{
			REP(N) { a[kk] = _mm_loadu_si128((__m128i*)A); A += sizeof(__m128i) / sizeof(byte); }

			REP(N) s[kk] = _mm_add_epi64(s[kk], _mm_sad_epu8(a[kk], z));
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm_add_epi64(s[kk], s[kk + K]);

		re += simd_u64(s, 0) + simd_u64(s, 1);
	}

	for (; i < n; ++i)
		re += *A++;

	return re;
}

TARGETSSE double SumSSE(double* A, int64 n, int64 sep)
{
	constexpr int N = 8;
	int64 i = 0;
	double re = 0;

	if (n >= N * sizeof(__m128d) / sizeof(double))
	{
		REP(4) _mm_prefetch((const char*)(A + sep * kk), _MM_HINT_T0);

		__m128d s[N], a[N];
		REP(N) s[kk] = _mm_setzero_pd();

		for (int64 l1 = n - N * sizeof(__m128d) / sizeof(double); i <= l1; i += N * sizeof(__m128d) / sizeof(double))
		{
			REP(N * 2) { _mm_prefetch((const char*)(A + sep * kk), _MM_HINT_T0); }

			REP(N) { a[kk] = _mm_set_pd(A[sep], A[0]); A += sep * sizeof(__m128d) / sizeof(double); }

			REP(N) s[kk] = _mm_add_pd(s[kk], a[kk]);
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm_add_pd(s[kk], s[kk + K]);

		re = _mm_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, A += sep)
	{
		volatile double v1 = *A;
		re += v1;
	}

	return re;
}

TARGETSSE double SumSSE(float* A, int64 n, int64 sep)
{
	int64 i = 0;
	double re = 0;

	if (n >= 4)
	{
		REP(4) { _mm_prefetch((const char*)A, _MM_HINT_T0); A += sep; }

		__m128d s1 = _mm_setzero_pd(), s2 = _mm_setzero_pd();

		for (int64 l1 = n - 4; i <= l1; i += 4)
		{
			REP(4) { _mm_prefetch((const char*)A, _MM_HINT_T0); A += sep; }

			s1 = _mm_add_pd(s1, _mm_set_pd(A[-7 * sep], A[-8 * sep]));
			s2 = _mm_add_pd(s2, _mm_set_pd(A[-5 * sep], A[-6 * sep]));
		}

		s1 = _mm_add_pd(s1, s2);
		re = _mm_reduce_add_pd(s1);
		A -= 4 * sep;
	}

	for (; i < n; ++i, A += sep)
	{
		volatile double v1 = *A;
		re += v1;
	}

	return re;
}

TARGETSSE float SumSSEx(float* A, int64 n, int64 sep)
{
	int64 i = 0;
	float re = 0;

	if (n >= 8)
	{
		__m128 s1 = _mm_setzero_ps(), s2 = _mm_setzero_ps();

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			s1 = _mm_add_ps(s1, _mm_set_ps(A[3 * sep], A[2 * sep], A[1 * sep], A[0 * sep]));
			s2 = _mm_add_ps(s2, _mm_set_ps(A[7 * sep], A[6 * sep], A[5 * sep], A[4 * sep]));
			A += 8 * sep;
		}

		s1 = _mm_add_ps(s1, s2);
		re = _mm_reduce_add_ps(s1);
	}

	for (; i < n; ++i, A += sep)
	{
		volatile float v1 = *A;
		re += v1;
	}

	return re;

	/*
	constexpr int N = 8;
	int64 i = 0;
	float re = 0;

	if (n >= N * sizeof(__m128) / sizeof(float))
	{
		__m128 s[N];
		REP(N) s[kk] = _mm_setzero_ps();

		for (int64 l1 = n - N * sizeof(__m128) / sizeof(float); i <= l1; i += N * sizeof(__m128) / sizeof(float))
		{
			REP(N) { s[kk] = _mm_add_ps(s[kk], _mm_set_ps(A[3 * sep], A[2 * sep], A[1 * sep], A[0 * sep])); A += sep * sizeof(__m128) / sizeof(float); }
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm_add_ps(s[kk], s[kk + K]);

		re = _mm_reduce_add_ps(s[0]);
	}

	for (; i < n; ++i, A += sep)
	{
		volatile float v1 = *A;
		re += v1;
	}

	return re;
	*/
}

TARGETSSE void SumSSE(double* A, double** B, int64 k, int64 n)
{
	constexpr int N = 16;
	int64 i = 0;

	if (n >= N * sizeof(__m128d) / sizeof(double))
	{
		__m128d a[N];

		for (int64 l1 = n - N * sizeof(__m128d) / sizeof(double); i <= l1; i += N * sizeof(__m128d) / sizeof(double))
		{
			REP(N) a[kk] = _mm_setzero_pd();

			for (int64 j = 0; j < k; ++j)
				REP(N) a[kk] = _mm_add_pd(a[kk], _mm_loadu_pd(&B[j][i + kk * sizeof(__m128d) / sizeof(double)]));

			REP(N) { _mm_storeu_pd(A, a[kk]); A += sizeof(__m128d) / sizeof(double); }
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

TARGETSSE void SumSSE(float* A, float** B, int64 k, int64 n)
{
	constexpr int N = 16;
	int64 i = 0;

	if (n >= N * sizeof(__m128) / sizeof(float))
	{
		__m128d a[N * 2];

		for (int64 l1 = n - N * sizeof(__m128) / sizeof(float); i <= l1; i += N * sizeof(__m128) / sizeof(float))
		{
			REP(N * 2) a[kk] = _mm_setzero_pd();
			
			for (int64 j = 0; j < k; ++j)
				REP(N)
				{
					__m128 v1 = _mm_loadu_ps(&B[j][i + kk * sizeof(__m128) / sizeof(float)]);
					a[0 + (kk << 1)] = _mm_add_pd(a[0 + (kk << 1)], _mm_cvtps_pd(v1));
					a[1 + (kk << 1)] = _mm_add_pd(a[1 + (kk << 1)], _mm_cvtps_pd(_mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(v1), 8))));
				}

			REP(N)
			{
				_mm_storeu_ps(A, _mm_shuffle_ps(
					_mm_cvtpd_ps(a[0 + (kk << 1)]), 
					_mm_cvtpd_ps(a[1 + (kk << 1)]), _MM_SHUFFLE(1, 0, 1, 0)));
				
				A += sizeof(__m128) / sizeof(float); 
			}
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
		for (int64 l1 = n - 16; i <= l1; i += 16)
		{
			__m128 a1 = _mm_setzero_ps();
			__m128 a2 = _mm_setzero_ps();
			__m128 a3 = _mm_setzero_ps();
			__m128 a4 = _mm_setzero_ps();

			for (int64 j = 0; j < k; ++j)
			{
				a1 = _mm_add_ps(a1, _mm_loadu_ps(&B[j][i +  0]));
				a2 = _mm_add_ps(a2, _mm_loadu_ps(&B[j][i +  4]));
				a3 = _mm_add_ps(a3, _mm_loadu_ps(&B[j][i +  8]));
				a4 = _mm_add_ps(a4, _mm_loadu_ps(&B[j][i + 12]));
			}

			_mm_storeu_ps(&A[i +  0], a1);
			_mm_storeu_ps(&A[i +  4], a2);
			_mm_storeu_ps(&A[i +  8], a3);
			_mm_storeu_ps(&A[i + 12], a4);
		}
	}

	for (; i < n; ++i)
	{
		float Ai = 0;
		for (int64 j = 0; j < k; ++j)
			Ai += B[j][i];
		A[i] = Ai;
	}
	*/
}

TARGETSSE double ProdSSE(double* A, int64 n)
{
	constexpr int N = 8;
	int64 i = 0;
	double re = 1;

	if (n >= N * sizeof(__m128d) / sizeof(double))
	{
		__m128d s[N], a[N];
		REP(N) s[kk] = _mm_set1_pd(1);

		for (int64 l1 = n - N * sizeof(__m128d) / sizeof(double); i <= l1; i += N * sizeof(__m128d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm_loadu_pd(A); A += sizeof(__m128d) / sizeof(double); }

			REP(N) s[kk] = _mm_mul_pd(s[kk], a[kk]);
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm_mul_pd(s[kk], s[kk + K]);

		re = _mm_reduce_mul_pd(s[0]);
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *A++;
		re *= v1;
	}

	return re;
}

TARGETSSE double ProdSSE(float* A, int64 n)
{
	constexpr int N = 16;
	int64 i = 0;
	double re = 1;

	if (n >= N / 2 * sizeof(__m128) / sizeof(float))
	{
		__m128d s[N];
		REP(N) s[kk] = _mm_set1_pd(1);

		for (int64 l1 = n - N / 2 * sizeof(__m128) / sizeof(float); i <= l1; i += N / 2 * sizeof(__m128) / sizeof(float))
		{
			REP(N / 2)
			{
				__m128 v1 = _mm_loadu_ps(A); A += sizeof(__m128) / sizeof(float);
				__m128d a1 = _mm_cvtps_pd(v1);
				__m128d a2 = _mm_cvtps_pd(_mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(v1), 8)));

				s[0 + (kk << 1)] = _mm_mul_pd(s[0 + (kk << 1)], a1);
				s[1 + (kk << 1)] = _mm_mul_pd(s[1 + (kk << 1)], a2);
			}
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm_mul_pd(s[kk], s[kk + K]);

		re = _mm_reduce_mul_pd(s[0]);
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A++;
		re *= v1;
	}

	return re;
}

TARGETSSE float ProdSSEx(float* A, int64 n)
{
	constexpr int N = 8;
	int64 i = 0;
	volatile float re = 1;

	if (n >= N * sizeof(__m128) / sizeof(float))
	{
		__m128 s[N], a[N];
		REP(N) s[kk] = _mm_set1_ps(1);

		for (int64 l1 = n - N * sizeof(__m128) / sizeof(float); i <= l1; i += N * sizeof(__m128) / sizeof(float))
		{
			REP(N) { a[kk] = _mm_loadu_ps(A); A += sizeof(__m128) / sizeof(float); }

			REP(N) s[kk] = _mm_mul_ps(s[kk], a[kk]);
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm_mul_ps(s[kk], s[kk + K]);

		re = _mm_reduce_mul_ps(s[0]);
	}

	for (; i < n; ++i)
		re *= *A++;

	return re;
}

TARGETSSE double ProdSSE(double* A, int64 n, int64 sep)
{
	int64 i = 0;
	double re = 1;

	if (n >= 8)
	{
		REP(8) { _mm_prefetch((const char*)A, _MM_HINT_T0); A += sep; }

		__m128d pd1 = _mm_set1_pd(1.0), pd2 = _mm_set1_pd(1.0);
		__m128d pd3 = _mm_set1_pd(1.0), pd4 = _mm_set1_pd(1.0);

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			REP(2) { _mm_prefetch((const char*)A, _MM_HINT_T0); A += sep; }
			pd1 = _mm_mul_pd(pd1, _mm_set_pd(A[-9 * sep], A[-10 * sep]));

			REP(2) { _mm_prefetch((const char*)A, _MM_HINT_T0); A += sep; }
			pd2 = _mm_mul_pd(pd2, _mm_set_pd(A[-9 * sep], A[-10 * sep]));

			REP(2) { _mm_prefetch((const char*)A, _MM_HINT_T0); A += sep; }
			pd3 = _mm_mul_pd(pd3, _mm_set_pd(A[-9 * sep], A[-10 * sep]));

			REP(2) { _mm_prefetch((const char*)A, _MM_HINT_T0); A += sep; }
			pd4 = _mm_mul_pd(pd4, _mm_set_pd(A[-9 * sep], A[-10 * sep]));
		}

		pd1 = _mm_mul_pd(_mm_mul_pd(pd1, pd2), _mm_mul_pd(pd3, pd4));
		re = _mm_reduce_mul_pd(pd1);
		A -= 8 * sep;
	}

	for (; i < n; ++i, A += sep)
	{
		volatile double v1 = *A;
		re *= v1;
	}

	return re;
}

TARGETSSE double ProdSSE(float* A, int64 n, int64 sep)
{
	int64 i = 0;
	double re = 1;

	if (n >= 8)
	{
		REP(8) { _mm_prefetch((const char*)A, _MM_HINT_T0); A += sep; }

		__m128d pd1 = _mm_set1_pd(1.0), pd2 = _mm_set1_pd(1.0);
		__m128d pd3 = _mm_set1_pd(1.0), pd4 = _mm_set1_pd(1.0);

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			REP(2) { _mm_prefetch((const char*)A, _MM_HINT_T0); A += sep; }
			pd1 = _mm_mul_pd(pd1, _mm_set_pd(A[-9 * sep], A[-10 * sep]));

			REP(2) { _mm_prefetch((const char*)A, _MM_HINT_T0); A += sep; }
			pd2 = _mm_mul_pd(pd2, _mm_set_pd(A[-9 * sep], A[-10 * sep]));

			REP(2) { _mm_prefetch((const char*)A, _MM_HINT_T0); A += sep; }
			pd3 = _mm_mul_pd(pd3, _mm_set_pd(A[-9 * sep], A[-10 * sep]));

			REP(2) { _mm_prefetch((const char*)A, _MM_HINT_T0); A += sep; }
			pd4 = _mm_mul_pd(pd4, _mm_set_pd(A[-9 * sep], A[-10 * sep]));
		}

		pd1 = _mm_mul_pd(_mm_mul_pd(pd1, pd2), _mm_mul_pd(pd3, pd4));
		re = _mm_reduce_mul_pd(pd1);
		A -= 8 * sep;
	}

	for (; i < n; ++i, A += sep)
	{
		volatile double v1 = *A;
		re *= v1;
	}

	return re;
}

TARGETSSE float ProdSSEx(float* A, int64 n, int64 sep)
{
	constexpr int N = 4;
	int64 i = 0;
	volatile float re = 1;

	if (n >= N * sizeof(__m128) / sizeof(float))
	{
		__m128 pd[4];
		REP(N) pd[kk] = _mm_set1_ps(1.0f);

		for (int64 l1 = n - N * sizeof(__m128) / sizeof(float); i <= l1; i += N * sizeof(__m128) / sizeof(float))
		{
			REP(N)
			{
				pd[kk] = _mm_mul_ps(pd[kk], _mm_set_ps(A[3 * sep], A[2 * sep], A[1 * sep], A[0 * sep]));
				A += 4 * sep;
			}
		}

		for (int K = sizeof(pd) / sizeof(pd[0]) / 2; K >= 1; K >>= 1)
			REP(K) pd[kk] = _mm_mul_ps(pd[kk], pd[kk + K]);

		re = _mm_reduce_mul_ps(pd[0]);
	}

	for (; i < n; ++i, A += sep)
	{
		volatile float v1 = *A;
		re *= v1;
	}

	return re;
}

TARGETSSE double SumSquareSSE(double* A, int64 n)
{
	constexpr int N = 8;
	int64 i = 0;
	double re = 0;

	if (n >= N * sizeof(__m128d) / sizeof(double))
	{
		__m128d s[N], a[N];
		REP(N) s[kk] = _mm_setzero_pd();

		for (int64 l1 = n - N * sizeof(__m128d) / sizeof(double); i <= l1; i += N * sizeof(__m128d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm_loadu_pd(A); A += sizeof(__m128d) / sizeof(double); }

			REP(N) a[kk] = _mm_mul_pd(a[kk], a[kk]);

			REP(N) s[kk] = _mm_add_pd(s[kk], a[kk]);
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm_add_pd(s[kk], s[kk + K]);

		re = _mm_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, ++A)
	{
		volatile double v1 = *A * *A;
		re += v1;
	}

	return re;
}

TARGETSSE double SumSquareSSE(float* A, int64 n)
{
	constexpr int N = 8;
	int64 i = 0;
	double re = 0;

	if (n >= N / 2 * sizeof(__m128) / sizeof(float))
	{
		__m128d s[N];
		REP(N) s[kk] = _mm_setzero_pd();

		for (int64 l1 = n - N / 2 * sizeof(__m128) / sizeof(float); i <= l1; i += N / 2 * sizeof(__m128) / sizeof(float))
		{
			REP(N / 2)
			{
				__m128 v1 = _mm_loadu_ps(A); A += sizeof(__m128) / sizeof(float);
				__m128d a1 = _mm_cvtps_pd(v1);
				__m128d a2 = _mm_cvtps_pd(_mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(v1), 8)));

				a1 = _mm_mul_pd(a1, a1);
				a2 = _mm_mul_pd(a2, a2);

				s[0 + (kk << 1)] = _mm_add_pd(s[0 + (kk << 1)], a1);
				s[1 + (kk << 1)] = _mm_add_pd(s[1 + (kk << 1)], a2);
			}
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm_add_pd(s[kk], s[kk + K]);

		re = _mm_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, ++A)
	{
		volatile double v1 = (double)*A * (double)*A;
		re += v1;
	}

	return re;
}

TARGETSSE float SumSquareSSEx(float* A, int64 n)
{
	constexpr int N = 8;
	int64 i = 0;
	float re = 0;

	if (n >= N * sizeof(__m128) / sizeof(float))
	{
		__m128 s[N], a[N];
		REP(N) s[kk] = _mm_setzero_ps();

		for (int64 l1 = n - N * sizeof(__m128) / sizeof(float); i <= l1; i += N * sizeof(__m128) / sizeof(float))
		{
			REP(N) { a[kk] = _mm_loadu_ps(A); A += sizeof(__m128) / sizeof(float); }

			REP(N) a[kk] = _mm_mul_ps(a[kk], a[kk]);

			REP(N) s[kk] = _mm_add_ps(s[kk], a[kk]);
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm_add_ps(s[kk], s[kk + K]);

		re = _mm_reduce_add_ps(s[0]);
	}

	for (; i < n; ++i, ++A)
	{
		volatile float v1 = *A * *A;
		re += v1;
	}

	return re;
}

TARGETSSE int64 SumSquareSSE(byte* A, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;
	uint64 re = 0;

	if (n >= N * sizeof(__m128i) / sizeof(byte))
	{
		__m128i a, t[2], s = _mm_setzero_si128();
		REP(2) t[kk] = _mm_setzero_si128();

		for (int64 l1 = n - N * sizeof(__m128i) / sizeof(byte); i <= l1; i += N * sizeof(__m128i) / sizeof(byte))
		{
			REP(N)
			{
				a = _mm_loadu_si128((__m128i*)A);
				A += sizeof(__m128i) / sizeof(byte);
				s = _mm_add_epi16(s, _mm_maddubs_epi16(a, a));
			}

			if ((i & (sizeof(__m128i) / sizeof(byte) * 128 - 1)) == 0)
			{
				t[0] = _mm_add_epi32(t[0], _mm_cvtepi16_epi32(s));
				t[1] = _mm_add_epi32(t[1], _mm_cvtepi16_epi32(_mm_srli_si128(s, 8)));
				s = _mm_setzero_si128();
			}
		}

		t[0] = _mm_add_epi32(t[0], _mm_cvtepi16_epi32(s));
		t[1] = _mm_add_epi32(t[1], _mm_cvtepi16_epi32(_mm_srli_si128(s, 8)));
		t[0] = _mm_add_epi32(t[0], t[1]);

		re = simp_i32(t, 0) + simp_i32(t, 1) + simp_i32(t, 2) + simp_i32(t, 3);
	}

	for (; i < n; ++i, ++A)
		re += *A * *A;

	return re;
}

TARGETSSE void SumSumSquareSSE(double* A, int64 n, double& sum, double& sumsq)
{
	constexpr int N = 8;
	int64 i = 0;
	double re1 = 0, re2 = 0;

	if (n >= N * sizeof(__m128d) / sizeof(double))
	{
		__m128d s1[N], s2[N], a[N];
		REP(N) s1[kk] = s2[kk] = _mm_setzero_pd();

		for (int64 l1 = n - N * sizeof(__m128d) / sizeof(double); i <= l1; i += N * sizeof(__m128d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm_loadu_pd(A); A += sizeof(__m128d) / sizeof(double); }

			REP(N) s1[kk] = _mm_add_pd(s1[kk], a[kk]);

			REP(N) a[kk] = _mm_mul_pd(a[kk], a[kk]);

			REP(N) s2[kk] = _mm_add_pd(s2[kk], a[kk]);
		}

		for (int K = sizeof(s1) / sizeof(s1[0]) / 2; K >= 1; K >>= 1)
		{
			REP(K) s1[kk] = _mm_add_pd(s1[kk], s1[kk + K]);
			REP(K) s2[kk] = _mm_add_pd(s2[kk], s2[kk + K]);
		}

		re1 = _mm_reduce_add_pd(s1[0]);
		re2 = _mm_reduce_add_pd(s2[0]);
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

TARGETSSE void SumSumSquareSSE(float* A, int64 n, double& sum, double& sumsq)
{
	constexpr int N = 16;
	int64 i = 0;
	double re1 = 0, re2 = 0;

	if (n >= N / 2 * sizeof(__m128) / sizeof(float))
	{
		__m128d s1[N], s2[N];
		REP(N) s1[kk] = s2[kk] = _mm_setzero_pd();

		for (int64 l1 = n - N / 2 * sizeof(__m128) / sizeof(float); i <= l1; i += N / 2 * sizeof(__m128) / sizeof(float))
		{
			REP(N / 2)
			{
				__m128 v1 = _mm_loadu_ps(A); A += sizeof(__m128) / sizeof(float);
				__m128d a1 = _mm_cvtps_pd(v1);
				__m128d a2 = _mm_cvtps_pd(_mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(v1), 8)));

				s1[0 + (kk << 1)] = _mm_add_pd(s1[0 + (kk << 1)], a1);
				s1[1 + (kk << 1)] = _mm_add_pd(s1[1 + (kk << 1)], a2);

				a1 = _mm_mul_pd(a1, a1);
				a2 = _mm_mul_pd(a2, a2);

				s2[0 + (kk << 1)] = _mm_add_pd(s2[0 + (kk << 1)], a1);
				s2[1 + (kk << 1)] = _mm_add_pd(s2[1 + (kk << 1)], a2);
			}
		}

		for (int K = sizeof(s1) / sizeof(s1[0]) / 2; K >= 1; K >>= 1)
		{
			REP(K) s1[kk] = _mm_add_pd(s1[kk], s1[kk + K]);
			REP(K) s2[kk] = _mm_add_pd(s2[kk], s2[kk + K]);
		}

		re1 = _mm_reduce_add_pd(s1[0]);
		re2 = _mm_reduce_add_pd(s2[0]);
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

TARGETSSE double SumProdDivSSE(double* A1, double* A2, double* B, int64 sep, int64 n)
{
	constexpr int N = 2;
	int64 i = 0;
	volatile double re1 = 0, re2 = 0;

	if (n >= N * sizeof(__m128) / sizeof(float))
	{
		REP(8) { _mm_prefetch((const char*)B, _MM_HINT_T0); B += sep; }

		__m128d s1[2], s2[2], b;
		REP(2) s1[kk] = s2[kk] = _mm_setzero_pd();

		for (int64 l1 = n - N * sizeof(__m128) / sizeof(float); i <= l1; i += N * sizeof(__m128) / sizeof(float))
		{
			REP(N)
			{
				_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
				_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
				b = _mm_set_pd(B[-9 * sep], B[-10 * sep]);
				s1[0] = _mm_add_pd(s1[0], _mm_mul_pd(_mm_loadu_pd(A1), b)); A1 += 2;
				s2[0] = _mm_add_pd(s2[0], _mm_mul_pd(_mm_loadu_pd(A2), b)); A2 += 2;

				_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
				_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
				b = _mm_set_pd(B[-9 * sep], B[-10 * sep]);
				s1[1] = _mm_add_pd(s1[1], _mm_mul_pd(_mm_loadu_pd(A1), b)); A1 += 2;
				s2[1] = _mm_add_pd(s2[1], _mm_mul_pd(_mm_loadu_pd(A2), b)); A2 += 2;
			}
		}

		s1[0] = _mm_add_pd(s1[0], s1[1]);
		s2[0] = _mm_add_pd(s2[0], s2[1]);

		re1 = _mm_reduce_add_pd(s1[0]);
		re2 = _mm_reduce_add_pd(s2[0]);

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

TARGETSSE double SumProdDivSSE(double* A1, float* A2, float* B, int64 sep, int64 n)
{
	constexpr int N = 2;
	int64 i = 0;
	volatile double re1 = 0, re2 = 0;

	if (n >= N * sizeof(__m128) / sizeof(float))
	{
		REP(6) { _mm_prefetch((const char*)B, _MM_HINT_T0); B += sep; }

		__m128 a2[N * 2];
		__m128d a1[N * 2], b[N * 2], s1[N * 2], s2[N * 2];
		REP(N * 2) s1[kk] = s2[kk] = _mm_setzero_pd();

		for (int64 l1 = n - N * sizeof(__m128) / sizeof(float); i <= l1; i += N * sizeof(__m128) / sizeof(float))
		{
			REP(N)
			{
				_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
				_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
				a1[0 + (kk << 1)] = _mm_loadu_pd(A1); A1 += 2;
				a2[0 + (kk << 1)] = _mm_loadu_ps(A2); A2 += 4;
				b[0 + (kk << 1)] = _mm_set_pd(B[-7 * sep], B[-8 * sep]);

				_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
				_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
				a1[1 + (kk << 1)] = _mm_loadu_pd(A1); A1 += 2;
				a2[1 + (kk << 1)] = _mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(a2[0 + (kk << 1)]), 8));
				b[1 + (kk << 1)] = _mm_set_pd(B[-7 * sep], B[-8 * sep]);

				s1[0 + (kk << 1)] = _mm_add_pd(s1[0 + (kk << 1)], _mm_mul_pd(a1[0 + (kk << 1)], b[0 + (kk << 1)]));
				s2[0 + (kk << 1)] = _mm_add_pd(s2[0 + (kk << 1)], _mm_mul_pd(_mm_cvtps_pd(a2[0 + (kk << 1)]), b[0 + (kk << 1)]));

				s1[1 + (kk << 1)] = _mm_add_pd(s1[1 + (kk << 1)], _mm_mul_pd(a1[1 + (kk << 1)], b[1 + (kk << 1)]));
				s2[1 + (kk << 1)] = _mm_add_pd(s2[1 + (kk << 1)], _mm_mul_pd(_mm_cvtps_pd(a2[1 + (kk << 1)]), b[1 + (kk << 1)]));
			}
		}

		for (int K = sizeof(s1) / sizeof(s1[0]) / 2; K >= 1; K >>= 1)
		{
			REP(K) s1[kk] = _mm_add_pd(s1[kk], s1[kk + K]);
			REP(K) s2[kk] = _mm_add_pd(s2[kk], s2[kk + K]);
		}

		re1 = _mm_reduce_add_pd(s1[0]);
		re2 = _mm_reduce_add_pd(s2[0]);

		B -= 6 * sep;
	}

	for (; i < n; ++i, A1++, A2++, B += sep)
	{
		volatile double v1 = (double)*A1 * (double)*B;
		volatile double v2 = (double)*A2 * (double)*B;
		re1 += v1;
		re2 += v2;
	}

	return re1 / re2;
}

TARGETSSE double SumProdDivSSE(float* A1, float* A2, float* B, int64 sep, int64 n)
{
	constexpr int N = 2;
	int64 i = 0;
	volatile double re1 = 0, re2 = 0;

	if (n >= N * sizeof(__m128) / sizeof(float))
	{
		REP(6) { _mm_prefetch((const char*)B, _MM_HINT_T0); B += sep; }

		__m128d s1[N * 2], s2[N * 2], b[N * 2];
		REP(N * 2) s1[kk] = s2[kk] = _mm_setzero_pd();

		for (int64 l1 = n - N * sizeof(__m128) / sizeof(float); i <= l1; i += N * sizeof(__m128) / sizeof(float))
		{
			REP(N)
			{
				_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
				_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
				b[0 + (kk << 1)] = _mm_set_pd(B[-7 * sep], B[-8 * sep]);

				__m128 a11 = _mm_loadu_ps(A1); A1 += sizeof(__m128) / sizeof(float);
				__m128 a21 = _mm_loadu_ps(A2); A2 += sizeof(__m128) / sizeof(float);

				_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
				_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
				b[1 + (kk << 1)] = _mm_set_pd(B[-7 * sep], B[-8 * sep]);

				__m128 a12 = _mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(a11), 8));
				__m128 a22 = _mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(a21), 8));

				s1[0 + (kk << 1)] = _mm_add_pd(s1[0 + (kk << 1)], _mm_mul_pd(_mm_cvtps_pd(a11), b[0 + (kk << 1)]));
				s2[0 + (kk << 1)] = _mm_add_pd(s2[0 + (kk << 1)], _mm_mul_pd(_mm_cvtps_pd(a21), b[0 + (kk << 1)]));
				s1[1 + (kk << 1)] = _mm_add_pd(s1[1 + (kk << 1)], _mm_mul_pd(_mm_cvtps_pd(a12), b[1 + (kk << 1)]));
				s2[1 + (kk << 1)] = _mm_add_pd(s2[1 + (kk << 1)], _mm_mul_pd(_mm_cvtps_pd(a22), b[1 + (kk << 1)])); 
			}
		}

		for (int K = sizeof(s1) / sizeof(s1[0]) / 2; K >= 1; K >>= 1)
		{
			REP(K) s1[kk] = _mm_add_pd(s1[kk], s1[kk + K]);
			REP(K) s2[kk] = _mm_add_pd(s2[kk], s2[kk + K]);
		}

		re1 = _mm_reduce_add_pd(s1[0]);
		re2 = _mm_reduce_add_pd(s2[0]);

		B -= 6 * sep;
	}

	for (; i < n; ++i, A1++, A2++, B += sep)
	{
		volatile double v1 = (double)*A1 * (double)*B;
		volatile double v2 = (double)*A2 * (double)*B;
		re1 += v1;
		re2 += v2;
	}

	return re1 / re2;
}

TARGETSSE float SumProdDivSSEx(float* A1, float* A2, float* B, int64 sep, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;
	volatile float re1 = 0, re2 = 0;

	if (n >= N * sizeof(__m128) / sizeof(float))
	{
		__m128 s1[N], s2[N], b;
		REP(N) s1[kk] = s2[kk] = _mm_setzero_ps();

		for (int64 l1 = n - N * sizeof(__m128) / sizeof(float); i <= l1; i += N * sizeof(__m128) / sizeof(float))
		{
			REP(N)
			{
				b = _mm_set_ps(B[3 * sep], B[2 * sep], B[1 * sep], B[0 * sep]); B += sizeof(__m128) / sizeof(float) * sep;
				s1[kk] = _mm_add_ps(s1[kk], _mm_mul_ps(_mm_loadu_ps(A1), b)); A1 += sizeof(__m128) / sizeof(float);
				s2[kk] = _mm_add_ps(s2[kk], _mm_mul_ps(_mm_loadu_ps(A2), b)); A2 += sizeof(__m128) / sizeof(float);
			}
		}

		for (int K = sizeof(s1) / sizeof(s1[0]) / 2; K >= 1; K >>= 1)
		{
			REP(K) s1[kk] = _mm_add_ps(s1[kk], s1[kk + K]);
			REP(K) s2[kk] = _mm_add_ps(s2[kk], s2[kk + K]);
		}

		re1 = _mm_reduce_add_ps(s1[0]);
		re2 = _mm_reduce_add_ps(s2[0]);
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

TARGETSSE double SumProdSSE(double* A, double* B, int64 sep, int64 n)
{
	int64 i = 0;
	double re = 0;

	if (n >= 2)
	{
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;

		__m128d s = _mm_setzero_pd();

		for (int64 l1 = n - 2; i <= l1; i += 2)
		{
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			s = _mm_add_pd(s, _mm_mul_pd(_mm_loadu_pd(A), _mm_set_pd(B[-3 * sep], B[-4 * sep])));
			A += 2;
		}

		re = _mm_reduce_add_pd(s); 
		B -= 2 * sep;
	}

	for (; i < n; ++i, ++A, ++B)
	{
		volatile double v1 = *A * *B;
		re += v1;
	}

	return re;
}

TARGETSSE double SumProdSSE(float* A, float* B, int64 sep, int64 n)
{
	int64 i = 0;
	double re = 0;

	if (n >= 4)
	{
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;

		__m128d s1 = _mm_setzero_pd(), s2 = _mm_setzero_pd();
		__m128 a;

		for (int64 l1 = n - 4; i <= l1; i += 4)
		{
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			a = _mm_loadu_ps(A); A += 4;
			s1 = _mm_add_pd(s1, _mm_mul_pd(_mm_cvtps_pd(a), _mm_set_pd(B[-5 * sep], B[-6 * sep])));

			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			a = _mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(a), 8));
			s2 = _mm_add_pd(s2, _mm_mul_pd(_mm_cvtps_pd(a), _mm_set_pd(B[-5 * sep], B[-6 * sep])));
		}

		s1 = _mm_add_pd(s1, s2);
		re = _mm_reduce_add_pd(s1); 

		B -= 4 * sep;
	}

	for (; i < n; ++i, A++, B += sep)
	{
		volatile double v1 = (double)*A * (double)*B;
		re += v1;
	}

	return re;
}

TARGETSSE float SumProdSSEx(float* A, float* B, int64 sep, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;
	volatile float re = 0;

	if (n >= N * sizeof(__m128) / sizeof(float))
	{
		REP(N) { _mm_prefetch((const char*)B, _MM_HINT_T0); B += sep; }

		__m128 s[N];
		REP(N) s[kk] = _mm_setzero_ps();

		for (int64 l1 = n - N * sizeof(__m128) / sizeof(float); i <= l1; i += N * sizeof(__m128) / sizeof(float))
		{
			REP(N)
			{
				_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
				_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
				_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
				_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;

				s[kk] = _mm_add_ps(s[kk], _mm_mul_ps(_mm_loadu_ps(A), _mm_set_ps(B[-5 * sep], B[-6 * sep], B[-7 * sep], B[-8 * sep])));
				A += sizeof(__m128) / sizeof(float);
			}
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm_add_ps(s[kk], s[kk + K]);

		re = _mm_reduce_add_ps(s[0]);

		B -= 4 * sep;
	}

	for (; i < n; ++i, A++, B += sep)
	{
		volatile float v1 = *A * *B;
		re += v1;
	}

	return re;
}

TARGETSSE double SumProdSSE(double* A, double* B, int64 n)
{
	constexpr int N = 16;
	int64 i = 0;
	double re = 0;

	if (n >= N * sizeof(__m128d) / sizeof(double))
	{
		__m128d s[N];
		REP(N) s[kk] = _mm_setzero_pd();

		for (int64 l1 = n - N * sizeof(__m128d) / sizeof(double); i <= l1; i += N * sizeof(__m128d) / sizeof(double))
		{
			REP(N)
			{
				s[kk] = _mm_add_pd(s[kk], _mm_mul_pd(_mm_loadu_pd(A), _mm_loadu_pd(B)));
				A += sizeof(__m128d) / sizeof(double);
				B += sizeof(__m128d) / sizeof(double);
			}
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm_add_pd(s[kk], s[kk + K]);

		re = _mm_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, ++A, ++B)
	{
		volatile double v1 = *A * *B;
		re += v1;
	}

	return re;
}

TARGETSSE double SumProdSSE(float* A, float* B, int64 n)
{
	constexpr int N = 8;
	int64 i = 0;
	double re = 0;

	if (n >= N / 2 * sizeof(__m128) / sizeof(float))
	{
		__m128d s[N];
		REP(N) s[kk] = _mm_setzero_pd();

		for (int64 l1 = n - N / 2 * sizeof(__m128) / sizeof(float); i <= l1; i += N / 2 * sizeof(__m128) / sizeof(float))
		{
			REP(N / 2)
			{
				__m128 v1 = _mm_loadu_ps(A), v2 = _mm_loadu_ps(B);
				A += sizeof(__m128) / sizeof(float);
				B += sizeof(__m128) / sizeof(float);
				s[0 + (kk << 1)] = _mm_add_pd(s[0 + (kk << 1)], _mm_mul_pd(_mm_cvtps_pd(v1), _mm_cvtps_pd(v2)));

				v1 = _mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(v1), 8));
				v2 = _mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(v2), 8));
				s[1 + (kk << 1)] = _mm_add_pd(s[1 + (kk << 1)], _mm_mul_pd(_mm_cvtps_pd(v1), _mm_cvtps_pd(v2)));
			}
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm_add_pd(s[kk], s[kk + K]);

		re = _mm_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, ++A, ++B)
	{
		volatile double v1 = (double)*A * (double)*B;
		re += v1;
	}

	return re;
}

TARGETSSE float SumProdSSEx(float* A, float* B, int64 n)
{
	constexpr int N = 8;
	int64 i = 0;
	volatile float re = 0;

	if (n >= N * sizeof(__m128) / sizeof(float))
	{
		__m128 s[N];
		REP(N) s[kk] = _mm_setzero_ps();

		for (int64 l1 = n - N * sizeof(__m128) / sizeof(float); i <= l1; i += N * sizeof(__m128) / sizeof(float))
		{
			REP(N)
			{
				s[kk] = _mm_add_ps(s[kk], _mm_mul_ps(_mm_loadu_ps(A), _mm_loadu_ps(B)));
				A += sizeof(__m128) / sizeof(float);
				B += sizeof(__m128) / sizeof(float);
			}
		}

		for (int K = sizeof(s) / sizeof(s[0]) / 2; K >= 1; K >>= 1)
			REP(K) s[kk] = _mm_add_ps(s[kk], s[kk + K]);

		re = _mm_reduce_add_ps(s[0]);
	}

	for (; i < n; ++i, ++A, ++B)
	{
		volatile float v1 = *A * *B;
		re += v1;
	}

	return re;
}

TARGETSSE void AddSSE(double* A, double* B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m128d) / sizeof(double))
	{
		__m128d a[N], b[N];

		for (int64 l1 = n - N * sizeof(__m128d) / sizeof(double); i <= l1; i += N * sizeof(__m128d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm_loadu_pd(A); A += sizeof(__m128d) / sizeof(double); }

			REP(N) { b[kk] = _mm_loadu_pd(B); B += sizeof(__m128d) / sizeof(double); }

			REP(N) a[kk] = _mm_add_pd(a[kk], b[kk]);

			REP(N) _mm_storeu_pd(A + (kk - N) * sizeof(__m128d) / sizeof(double), a[kk]);
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGETSSE void AddSSE(float* A, float* B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m128) / sizeof(float))
	{
		__m128 a[N], b[N];

		for (int64 l1 = n - N * sizeof(__m128) / sizeof(float); i <= l1; i += N * sizeof(__m128) / sizeof(float))
		{
			REP(N) { a[kk] = _mm_loadu_ps(A); A += sizeof(__m128) / sizeof(float); }

			REP(N) { b[kk] = _mm_loadu_ps(B); B += sizeof(__m128) / sizeof(float); }

			REP(N) a[kk] = _mm_add_ps(a[kk], b[kk]);

			REP(N) _mm_storeu_ps(A + (kk - N) * sizeof(__m128) / sizeof(float), a[kk]);
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGETSSE void AddSSE(int64* A, int64* B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m128i) / sizeof(int64))
	{
		__m128i a[N], b[N];

		for (int64 l1 = n - N * sizeof(__m128i) / sizeof(int64); i <= l1; i += N * sizeof(__m128i) / sizeof(int64))
		{
			REP(N) { a[kk] = _mm_loadu_si128((__m128i*)A); A += sizeof(__m128i) / sizeof(int64); }

			REP(N) { b[kk] = _mm_loadu_si128((__m128i*)B); B += sizeof(__m128i) / sizeof(int64); }

			REP(N) a[kk] = _mm_add_epi64(a[kk], b[kk]);

			REP(N) _mm_storeu_si128((__m128i*)(A + (kk - N) * sizeof(__m128i) / sizeof(int64)), a[kk]);
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGETSSE void AddSSE(int* A, int* B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m128i) / sizeof(int))
	{
		__m128i a[N], b[N];

		for (int64 l1 = n - N * sizeof(__m128i) / sizeof(int); i <= l1; i += N * sizeof(__m128i) / sizeof(int))
		{
			REP(N) { a[kk] = _mm_loadu_si128((__m128i*)A); A += sizeof(__m128i) / sizeof(int); }

			REP(N) { b[kk] = _mm_loadu_si128((__m128i*)B); B += sizeof(__m128i) / sizeof(int); }

			REP(N) a[kk] = _mm_add_epi32(a[kk], b[kk]);

			REP(N) _mm_storeu_si128((__m128i*)(A + (kk - N) * sizeof(__m128i) / sizeof(int)), a[kk]);
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGETSSE void AddSSE(double* A, double B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m128d) / sizeof(double))
	{
		__m128d b = _mm_set1_pd(B), a[N];

		for (int64 l1 = n - N * sizeof(__m128d) / sizeof(double); i <= l1; i += N * sizeof(__m128d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm_loadu_pd(A); A += sizeof(__m128d) / sizeof(double); }

			REP(N) a[kk] = _mm_add_pd(a[kk], b);

			REP(N) _mm_storeu_pd(A + (kk - N) * sizeof(__m128d) / sizeof(double), a[kk]);
		}
	}

	for (; i < n; ++i, A++)
		*A += B;
}

TARGETSSE void AddSSE(float* A, float B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m128) / sizeof(float))
	{
		__m128 b = _mm_set1_ps(B), a[N];

		for (int64 l1 = n - N * sizeof(__m128) / sizeof(float); i <= l1; i += N * sizeof(__m128) / sizeof(float))
		{
			REP(N) { a[kk] = _mm_loadu_ps(A); A += sizeof(__m128) / sizeof(float); }

			REP(N) a[kk] = _mm_add_ps(a[kk], b);

			REP(N) _mm_storeu_ps(A + (kk - N) * sizeof(__m128) / sizeof(float), a[kk]);
		}
	}

	for (; i < n; ++i, A++)
		*A += B;
}

TARGETSSE void MulSSE(double* C, double* A, double* B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m128d) / sizeof(double))
	{
		__m128d a[N], b[N];

		for (int64 l1 = n - N * sizeof(__m128d) / sizeof(double); i <= l1; i += N * sizeof(__m128d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm_loadu_pd(A); A += sizeof(__m128d) / sizeof(double); }

			REP(N) { b[kk] = _mm_loadu_pd(B); B += sizeof(__m128d) / sizeof(double); }

			REP(N) a[kk] = _mm_mul_pd(a[kk], b[kk]);

			REP(N) { _mm_storeu_pd(C, a[kk]); C += sizeof(__m128d) / sizeof(double); }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *A++ * *B++;
		*C++ = v1;
	}
}

TARGETSSE void MulSSE(float* C, float* A, float* B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m128) / sizeof(float))
	{
		__m128 a[N], b[N];

		for (int64 l1 = n - N * sizeof(__m128) / sizeof(float); i <= l1; i += N * sizeof(__m128) / sizeof(float))
		{
			REP(N) { a[kk] = _mm_loadu_ps(A); A += sizeof(__m128) / sizeof(float); }

			REP(N) { b[kk] = _mm_loadu_ps(B); B += sizeof(__m128) / sizeof(float); }

			REP(N) a[kk] = _mm_mul_ps(a[kk], b[kk]);

			REP(N) { _mm_storeu_ps(C, a[kk]); C += sizeof(__m128) / sizeof(float); }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A++ * *B++;
		*C++ = v1;
	}
}

TARGETSSE void MulSSE(double* C, double* A, double B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m128d) / sizeof(double))
	{
		__m128d b = _mm_set1_pd(B), a[N];

		for (int64 l1 = n - N * sizeof(__m128d) / sizeof(double); i <= l1; i += N * sizeof(__m128d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm_loadu_pd(A); A += sizeof(__m128d) / sizeof(double); }

			REP(N) a[kk] = _mm_mul_pd(a[kk], b);

			REP(N) { _mm_storeu_pd(C, a[kk]); C += sizeof(__m128d) / sizeof(double); }
		}
	}
	
	for (; i < n; ++i)
	{
		volatile double v1 = *A++ * B;
		*C++ = v1;
	}
}

TARGETSSE void MulSSE(float* C, float* A, float B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m128) / sizeof(float))
	{
		__m128 b = _mm_set1_ps(B), a[N];

		for (int64 l1 = n - N * sizeof(__m128) / sizeof(float); i <= l1; i += N * sizeof(__m128) / sizeof(float))
		{
			REP(N) { a[kk] = _mm_loadu_ps(A); A += sizeof(__m128) / sizeof(float); }

			REP(N) a[kk] = _mm_mul_ps(a[kk], b);

			REP(N) { _mm_storeu_ps(C, a[kk]); C += sizeof(__m128) / sizeof(float); }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A++ * B;
		*C++ = v1;
	}
}

TARGETSSE void MulSSE(double* A, double B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m128d) / sizeof(double))
	{
		__m128d b = _mm_set1_pd(B), a[N];

		for (int64 l1 = n - N * sizeof(__m128d) / sizeof(double); i <= l1; i += N * sizeof(__m128d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm_loadu_pd(A); A += sizeof(__m128d) / sizeof(double); }

			REP(N) a[kk] = _mm_mul_pd(a[kk], b);

			REP(N) _mm_storeu_pd(A + (kk - N) * sizeof(__m128d) / sizeof(double), a[kk]);
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *A * B;
		*A++ = v1;
	}
}

TARGETSSE void MulSSE(float* A, float B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m128) / sizeof(float))
	{
		__m128 b = _mm_set1_ps(B), a[N];

		for (int64 l1 = n - N * sizeof(__m128) / sizeof(float); i <= l1; i += N * sizeof(__m128) / sizeof(float))
		{
			REP(N) { a[kk] = _mm_loadu_ps(A); A += sizeof(__m128) / sizeof(float); }

			REP(N) a[kk] = _mm_mul_ps(a[kk], b);

			REP(N) _mm_storeu_ps(A + (kk - N) * sizeof(__m128) / sizeof(float), a[kk]);
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A * B;
		*A++ = v1;
	}
}

TARGETSSE void AddProdSSE(double* C, double* A, double* B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(__m128d) / sizeof(double))
	{
		__m128d a[N];

		for (int64 l1 = n - N * sizeof(__m128d) / sizeof(double); i <= l1; i += N * sizeof(__m128d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm_loadu_pd(A); A += sizeof(__m128d) / sizeof(double); }

			REP(N) { a[kk] = _mm_mul_pd(a[kk], _mm_loadu_pd(B)); B += sizeof(__m128d) / sizeof(double); }

			double* C2 = C;

			REP(N) { a[kk] = _mm_add_pd(a[kk], _mm_loadu_pd(C)); C += sizeof(__m128d) / sizeof(double); }

			REP(N) { _mm_storeu_pd(C2, a[kk]); C2 += sizeof(__m128d) / sizeof(double); }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *A++ * *B++;
		*C++ += v1;
	}
}

TARGETSSE void AddProdSSE(float* C, float* A, float* B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= 4 * N * sizeof(__m128) / sizeof(float))
	{
		__m128 a[N];

		for (int64 l1 = n - N * sizeof(__m128) / sizeof(float); i <= l1; i += N * sizeof(__m128) / sizeof(float))
		{
			REP(N) { a[kk] = _mm_loadu_ps(A); A += sizeof(__m128) / sizeof(float); }

			REP(N) { a[kk] = _mm_mul_ps(a[kk], _mm_loadu_ps(B)); B += sizeof(__m128) / sizeof(float); }

			float* C2 = C;

			REP(N) { a[kk] = _mm_add_ps(a[kk], _mm_loadu_ps(C)); C += sizeof(__m128) / sizeof(float); }

			REP(N) { _mm_storeu_ps(C2, a[kk]); C2 += sizeof(__m128) / sizeof(float); }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A++ * *B++;
		*C++ += v1;
	}
}

TARGETSSE void AddProdSSE(double* C, double* A, double B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= 4 * N * sizeof(__m128d) / sizeof(double))
	{
		__m128d a[N], b = _mm_set1_pd(B);

		for (int64 l1 = n - N * sizeof(__m128d) / sizeof(double); i <= l1; i += N * sizeof(__m128d) / sizeof(double))
		{
			REP(N) { a[kk] = _mm_loadu_pd(A); A += sizeof(__m128d) / sizeof(double); }

			REP(N)  a[kk] = _mm_mul_pd(a[kk], b);

			double* C2 = C;

			REP(N) { a[kk] = _mm_add_pd(a[kk], _mm_loadu_pd(C)); C += sizeof(__m128d) / sizeof(double); }

			REP(N) { _mm_storeu_pd(C2, a[kk]); C2 += sizeof(__m128d) / sizeof(double); }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *A++ * B;
		*C++ += v1;
	}
}

TARGETSSE void AddProdSSE(double* C, float* A, double B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= 4 * N * sizeof(__m128d) / sizeof(double))
	{
		__m128d a[N], b = _mm_set1_pd(B);

		for (int64 l1 = n - N * sizeof(__m128d) / sizeof(double); i <= l1; i += N * sizeof(__m128d) / sizeof(double))
		{
			REP(N / 2) 
			{
				__m128 v1 = _mm_loadu_ps(A); A += sizeof(__m128) / sizeof(float);
				
				a[0 + (kk << 1)] = _mm_cvtps_pd(v1); 
				a[1 + (kk << 1)] = _mm_cvtps_pd(_mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(v1), 8)));
			}

			REP(N)  a[kk] = _mm_mul_pd(a[kk], b);

			double* C2 = C;

			REP(N) { a[kk] = _mm_add_pd(a[kk], _mm_loadu_pd(C)); C += sizeof(__m128d) / sizeof(double); }

			REP(N) { _mm_storeu_pd(C2, a[kk]); C2 += sizeof(__m128d) / sizeof(double); }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *A++ * B;
		*C++ += v1;
	}
}

TARGETSSE void AddProdSSE(float* C, float* A, float B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= 4 * N * sizeof(__m128) / sizeof(float))
	{
		__m128 a[N], b = _mm_set1_ps(B);

		for (int64 l1 = n - N * sizeof(__m128) / sizeof(float); i <= l1; i += N * sizeof(__m128) / sizeof(float))
		{
			REP(N) { a[kk] = _mm_loadu_ps(A); A += sizeof(__m128) / sizeof(float); }

			REP(N)  a[kk] = _mm_mul_ps(a[kk], b);

			float* C2 = C;

			REP(N) { a[kk] = _mm_add_ps(a[kk], _mm_loadu_ps(C)); C += sizeof(__m128) / sizeof(float); }

			REP(N) { _mm_storeu_ps(C2, a[kk]); C2 += sizeof(__m128) / sizeof(float); }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A++ * B;
		*C++ += v1;
	}
}

TARGETSSE void UnifySSE(double* A, int64 n)
{
	constexpr int N = 16;
	int64 i = 0;
	double invsum = 1.0 / (SumSSE(A, n) + n * MIN_FREQ);

	if (n >= N * sizeof(__m128) / sizeof(double))
	{
		__m128d a[N], minf = _mm_set1_pd(MIN_FREQ), invs = _mm_set1_pd(invsum);

		for (int64 l1 = n - N * sizeof(__m128) / sizeof(double); i <= l1; i += N * sizeof(__m128) / sizeof(double))
		{
			double* A2 = A;

			REP(N) { a[kk] = _mm_loadu_pd(A); A += sizeof(__m128) / sizeof(double); }

			REP(N) a[kk] = _mm_add_pd(a[kk], minf);

			REP(N) a[kk] = _mm_mul_pd(a[kk], invs);

			REP(N) { _mm_storeu_pd(A2, a[kk]); A2 += sizeof(__m128) / sizeof(double); }
		}
	}

	for (; i < n; ++i, ++A)
	{
		volatile double v1 = (double)*A + MIN_FREQ;
		*A = v1 * invsum;
	}
}

TARGETSSE void UnifySSE(float* A, int64 n)
{
	constexpr int N = 8;
	int64 i = 0;
	double invsum = 1.0 / (SumSSE(A, n) + n * MIN_FREQ);

	if (n >= N / 2 * sizeof(__m128) / sizeof(float))
	{
		__m128d a[N], minf = _mm_set1_pd(MIN_FREQ), invs = _mm_set1_pd(invsum);

		for (int64 l1 = n - N / 2 * sizeof(__m128) / sizeof(float); i <= l1; i += N / 2 * sizeof(__m128) / sizeof(float))
		{
			float* A2 = A;

			REP(N) { a[kk] = _mm_cvtps_pd(_mm_loadu_ps(A)); A += sizeof(__m128) / sizeof(float) / 2; }

			REP(N) a[kk] = _mm_add_pd(a[kk], minf);

			REP(N) a[kk] = _mm_mul_pd(a[kk], invs);

			REP(N / 2) 
			{ 
				_mm_storeu_ps(A2, _mm_shuffle_ps(
					_mm_cvtpd_ps(a[0 + (kk << 1)]),
					_mm_cvtpd_ps(a[1 + (kk << 1)]), _MM_SHUFFLE(1, 0, 1, 0)));

				A2 += sizeof(__m128) / sizeof(float); 
			}
		}
	}

	for (; i < n; ++i, ++A)
	{
		volatile double v1 = (double)*A + MIN_FREQ;
		*A = v1 * invsum;
	}
}

TARGETSSE char* StrNextIdxSSE(char* A, char val, int64 rep, int64 n)
{
	A++; n--;
	int64 i = 0;

	if (n >= 16)
	{
		__m128i r = _mm_setzero_si128(), v = _mm_set1_epi8(val), o = _mm_setzero_si128();

		for (int64 l1 = n - 16; i <= l1; )
		{
			__m128i a = _mm_loadu_si128((__m128i*)A);
			r = _mm_cmpeq_epi8(a, v);

			if (_mm_testz_si128(r, r)) { A += 16; i += 16; continue; }

			r = _mm_sad_epu8(o, _mm_sub_epi8(o, r)); //4 int 64, each 8 bytes

			for (int64 j = 0; j < 2; ++j, A += 8, i += 8)
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

TARGETSSE int64 CountCharSSE(char* A, char val, int64 n)
{
	uint64 re = 0;
	int64 i = 0;

	if (n >= 16)
	{
		__m128i r = _mm_setzero_si128(), v = _mm_set1_epi8(val), o = _mm_setzero_si128();

		for (int64 l1 = n - 16; i <= l1; i += 16)
		{
			r = _mm_add_epi64(r, _mm_sad_epu8(o, _mm_sub_epi8(o, _mm_cmpeq_epi8(_mm_loadu_si128((__m128i*)A), v))));
			A += 16;
		}

		re = simd_u64(r, 0) + simd_u64(r, 1);
	}

	for (; i < n; ++i, A++)
		if (*A == val) re++;

	return (int64)re;
}

#endif