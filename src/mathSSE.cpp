/* SSE Instruction Set Functions */

#include "vcfpop.h"

#ifndef __aarch64__

template struct RNGSSE<double>;
template struct RNGSSE<float >;
template TARGETSSE void RNGSSE<double>::Integer<uint  >(uint  * re, int64 n, uint   minv, uint   maxv);
template TARGETSSE void RNGSSE<double>::Integer<uint64>(uint64* re, int64 n, uint64 minv, uint64 maxv);
template TARGETSSE void RNGSSE<float >::Integer<uint  >(uint  * re, int64 n, uint   minv, uint   maxv);
template TARGETSSE void RNGSSE<float >::Integer<uint64>(uint64* re, int64 n, uint64 minv, uint64 maxv);

#ifndef _RNGSSE_FP64
#define XS 32

/* Initialize rng */
TARGETSSE RNGSSE<double>::RNGSSE()
{

}

/* Initialize rng */
TARGETSSE RNGSSE<double>::RNGSSE(uint64 seed, uint64 salt)
{
	__m128i a[XS], s, m;

	UNROLL(XS) { a[kk] = _mm_set_epi64x(seed + 1, seed + 0); seed += 2; }

	s = _mm_set1_epi64x(salt);
	m = _mm_set1_epi32(0x5bd1e995);

	UNROLL(XS) a[kk] = _mm_xor_si128(a[kk], _mm_slli_epi64(_mm_andnot_si128(a[kk], _mm_set1_epi64x(0xFFFFFFFF)), 32));

	s             = _mm_xor_si128(s    , _mm_slli_epi64(_mm_andnot_si128(s    , _mm_set1_epi64x(0xFFFFFFFF)), 32));

	// uint s = s ^ 4;
	s = _mm_xor_si128(s, _mm_set1_epi32(4));

	// a *= m;
	UNROLL(XS) a[kk] = _mm_mullo_epi32(a[kk], m);

	// a ^= a >> 24;
	UNROLL(XS) a[kk] = _mm_xor_si128(a[kk], _mm_srli_epi32(a[kk], 24));

	// a *= m;
	UNROLL(XS) a[kk] = _mm_mullo_epi32(a[kk], m);

	// s *= m;
	s = _mm_mullo_epi32(s, m);

	// a ^= s;
	UNROLL(XS) a[kk] = _mm_xor_si128(a[kk], s);

	// a ^= a >> 13;
	UNROLL(XS) a[kk] = _mm_xor_si128(a[kk], _mm_srli_epi32(a[kk], 13));

	// a *= m;
	UNROLL(XS) a[kk] = _mm_mullo_epi32(a[kk], m);

	// a ^= a >> 15;
	UNROLL(XS) a[kk] = _mm_xor_si128(a[kk], _mm_srli_epi32(a[kk], 15));

	// original
	UNROLL(XS) x[kk] = _mm_xor_si128(_mm_set1_epi64x(0x159A55E5075BCD15), a[kk]);

	UNROLL(XS) a[kk] = _mm_slli_epi64(a[kk], 6);

	UNROLL(XS) y[kk] = _mm_xor_si128(_mm_set1_epi64x(0x054913331F123BB5), a[kk]);
}

/* Draw 64 64-bit integers in [0,n), 64*n frequencies are in arr */
TARGETSSE void RNGSSE<double>::Poly(__m128d* arr, int n, __m128i* re)
{
	__m128d t[XS], s[XS];
	__m128d one = _mm_set1_pd(1.0);
	__m128i mask1 = _mm_set1_epi64x(0x000FFFFFFFFFFFFF);
	__m128i mask2 = _mm_set1_epi64x(0x3FF0000000000000);
	__m128i* r = (__m128i*)t;

	UNROLL(XS) s[kk] = _mm_setzero_pd();

	for (int i = 0; i < n * XS; i += XS)
		UNROLL(XS) s[kk] = _mm_add_pd(s[kk], arr[kk + i]);

	XorShift();

	UNROLL(XS) r[kk] = _mm_add_epi64(x[kk], y[kk]);

	UNROLL(XS) r[kk] = _mm_and_si128(r[kk], mask1);

	UNROLL(XS) r[kk] = _mm_or_si128(r[kk], mask2);

	UNROLL(XS) t[kk] = _mm_sub_pd(t[kk], one);

	UNROLL(XS) t[kk] = _mm_mul_pd(t[kk], s[kk]);

	__m128i midx[XS], nidx = _mm_setzero_si128(), ninc = _mm_set1_epi64x(1);
	__m128d f[XS], b[XS];
	UNROLL(XS) midx[kk] = _mm_set1_epi64x(n - 1);
	UNROLL(XS) f[kk] = _mm_setzero_pd();

	for (int i = 0; i < n * XS; i += XS)
	{
		UNROLL(XS) b[kk] = _mm_cmplt_pd(t[kk], arr[kk + i]);

		UNROLL(XS) t[kk] = _mm_sub_pd(t[kk], arr[kk + i]);

		UNROLL(XS) b[kk] = _mm_andnot_pd(f[kk], b[kk]);

		UNROLL(XS) f[kk] = _mm_or_pd(f[kk], b[kk]);
			
		UNROLL(XS) midx[kk] = _mm_castpd_si128(_mm_blendv_pd(_mm_castsi128_pd(midx[kk]), _mm_castsi128_pd(nidx), b[kk]));//ok

		nidx = _mm_add_epi64(nidx, ninc);
	}

	UNROLL(XS) re[kk] = midx[kk];
}

/* Draw uniform distriubted intergers */
TARGETSSE void RNGSSE<double>::XorShift()
{
	__m128i a[XS], b[XS];

	UNROLL(XS) a[kk] = x[kk];

	UNROLL(XS) b[kk] = y[kk];

	UNROLL(XS) x[kk] = b[kk];

	UNROLL(XS) a[kk] = _mm_xor_si128(a[kk], _mm_slli_epi64(a[kk], 23));

	UNROLL(XS) a[kk] = _mm_xor_si128(a[kk], _mm_srli_epi64(a[kk], 18));

	UNROLL(XS) a[kk] = _mm_xor_si128(a[kk], b[kk]);

	UNROLL(XS) a[kk] = _mm_xor_si128(a[kk], _mm_srli_epi64(b[kk], 5));

	UNROLL(XS) y[kk] = a[kk];
}

/* Draw uniform distriubted integers */
template<typename INT>
TARGETSSE void RNGSSE<double>::Integer(INT* re, int64 n, INT minv, INT maxv)
{
	constexpr int xesize = sizeof(x) / sizeof(INT);
	int64 i = 0;
	INT modv = maxv - minv;
	INT* rei = re;

	for (; i <= n - xesize; i += xesize)
	{
		XorShift();
		UNROLL(XS) { _mm_storeu_si128((__m128i*)rei, _mm_add_epi64(x[kk], y[kk])); rei += E128B / sizeof(INT); };
	}

	if (i != n)
	{
		__m128i re2[XS];
		XorShift();
		UNROLL(XS) re2[kk] = _mm_add_epi64(x[kk], y[kk]);
		SetVal((INT*)rei, (INT*)re2, n - i);
	}

	if (maxv != (INT)-1 || minv != 0)
	{
		for (i = 0; i < n; ++i)
			re[i] = re[i] % modv + minv;
	}
}

/* Draw uniform distriubted real numbers */
TARGETSSE void RNGSSE<double>::Uniform(double* re, int n, double minv, double maxv)
{
	constexpr int xesize = sizeof(x) / sizeof(double);
	int i = 0;
	double range = maxv - minv;

	__m128i mask1 = _mm_set1_epi64x(0x000FFFFFFFFFFFFF);
	__m128i mask2 = _mm_set1_epi64x(0x3FF0000000000000);
	__m128d v1 = _mm_set1_pd(minv - range);
	__m128d v2 = _mm_set1_pd(range);

	if (range == 1.0)
	{
		for (; i <= n - xesize; i += xesize)
		{
			XorShift();
			UNROLL(XS) { _mm_storeu_pd(re, _mm_add_pd(_mm_castsi128_pd(_mm_or_si128(_mm_and_si128(_mm_add_epi64(x[kk], y[kk]), mask1), mask2)), v1)); re += E128D; }
		}
	}
	else
	{
		for (; i <= n - xesize; i += xesize)
		{
			XorShift();
			UNROLL(XS) { _mm_storeu_pd(re, _mm_fmaddx_pd(_mm_castsi128_pd(_mm_or_si128(_mm_and_si128(_mm_add_epi64(x[kk], y[kk]), mask1), mask2)), v2, v1)); re += E128D; }
		}
	}

	if (i != n)
	{
		__m128d ref[XS];
		XorShift();
		UNROLL(XS) ref[kk] = _mm_fmaddx_pd(_mm_castsi128_pd(_mm_or_si128(_mm_and_si128(_mm_add_epi64(x[kk], y[kk]), mask1), mask2)), v2, v1);
		SetVal((double*)re, (double*)ref, n - i);
	}
}

/* Draw uniform distriubted real numbers */
TARGETSSE void RNGSSE<double>::Normal(double* re, int n, double mean, double sd)
{
	constexpr int xhsize = XS / 2;
	constexpr int xesize = XS * E128D;

	int i = 0;

	__m128i mask1 = _mm_set1_epi64x(0x000FFFFFFFFFFFFF);
	__m128i mask2 = _mm_set1_epi64x(0x3FF0000000000000);
	__m128d v1 = _mm_set1_pd(-1);
	__m128d min_freq = _mm_set1_pd(MIN_FREQ);
	__m128d pi2 = _mm_set1_pd(2.0 * M_PI);
	__m128d mu = _mm_set1_pd(mean);
	__m128d s = _mm_set1_pd(sd);

	for (; i <= n - xesize; i += xesize)
	{
		XorShift();
		UNROLL(XS) _mm_storeu_pd(re + E128D * kk, _mm_add_pd(_mm_castsi128_pd(_mm_or_si128(_mm_and_si128(_mm_add_epi64(x[kk], y[kk]), mask1), mask2)), v1));

		__m128d u1, u2, u3, u4;
		for (int j = 0; j < xhsize; ++j)
		{
			u1 = _mm_max_pd(_mm_loadu_pd(re + j * E128D), min_freq);
			u2 = _mm_mul_pd(_mm_loadu_pd(re + j * E128D + xhsize * E128D), pi2);

			UNROLL(E128D) simd_f64(u1, kk) = sqrt(-2.0 * log(simd_f64(u1, kk)));
			UNROLL(E128D) simd_f64(u3, kk) = cos(simd_f64(u2, kk));
			UNROLL(E128D) simd_f64(u4, kk) = sin(simd_f64(u2, kk));

			_mm_storeu_pd(re + j * E128D, _mm_mul_pd(u1, u3));
			_mm_storeu_pd(re + j * E128D + xhsize * E128D, _mm_mul_pd(u1, u4));
		}
			
		if (sd != 1 || mean != 0)
			UNROLL(XS) { _mm_storeu_pd(re, _mm_fmaddx_pd(_mm_loadu_pd(re), s, mu)); re += E128D; }
		else
			re += XS * E128D;
	}

	if (i != n)
	{
		double ref[XS * E128D];
		
		XorShift();
		UNROLL(XS) _mm_storeu_pd(ref + E128D * kk, _mm_add_pd(_mm_castsi128_pd(_mm_or_si128(_mm_and_si128(_mm_add_epi64(x[kk], y[kk]), mask1), mask2)), v1));

		__m128d u1, u2, u3, u4;
		for (int j = 0; j < xhsize; ++j)
		{
			u1 = _mm_max_pd(_mm_loadu_pd(ref + j * E128D), min_freq);
			u2 = _mm_mul_pd(_mm_loadu_pd(ref + j * E128D + xhsize * E128D), pi2);

			UNROLL(E128D) simd_f64(u1, kk) = sqrt(-2.0 * log(simd_f64(u1, kk)));
			UNROLL(E128D) simd_f64(u3, kk) = cos(simd_f64(u2, kk));
			UNROLL(E128D) simd_f64(u4, kk) = sin(simd_f64(u2, kk));

			_mm_storeu_pd(ref + j * E128D, _mm_mul_pd(u1, u3));
			_mm_storeu_pd(ref + j * E128D + xhsize * E128D, _mm_mul_pd(u1, u4));
		}
			
		if (sd != 1 || mean != 0)
			UNROLL(XS) _mm_storeu_pd(ref + kk * E128D, _mm_fmaddx_pd(_mm_loadu_pd(ref), s, mu)); 

		SetVal((double*)re, (double*)ref, n - i);
	}
}
#endif

#ifndef _RNGSSE_FP32
#define XS 16
#define XS2 32

/* Initialize rng */
TARGETSSE RNGSSE<float>::RNGSSE()
{

}

/* Initialize rng */
TARGETSSE RNGSSE<float>::RNGSSE(uint64 seed, uint64 salt)
{
	__m128i a[XS], s, m;
	UNROLL(XS) { a[kk] = _mm_set_epi32(Mix(seed +  3), Mix(seed +  2), Mix(seed + 1), Mix(seed + 0)); seed += 4; }

	s = _mm_set1_epi32(Mix(salt));
	m = _mm_set1_epi32(0x5bd1e995);

	// uint s = s ^ 4;
	s = _mm_xor_si128(s, _mm_set1_epi32(4));

	// a *= m;
	UNROLL(XS) a[kk] = _mm_mullo_epi32(a[kk], m);

	// a ^= a >> 24;
	UNROLL(XS) a[kk] = _mm_xor_si128(a[kk], _mm_srli_epi32(a[kk], 24));

	// a *= m;
	UNROLL(XS) a[kk] = _mm_mullo_epi32(a[kk], m);

	// s *= m;
	s = _mm_mullo_epi32(s, m);

	// a ^= s;
	UNROLL(XS) a[kk] = _mm_xor_si128(a[kk], s);

	// a ^= a >> 13;
	UNROLL(XS) a[kk] = _mm_xor_si128(a[kk], _mm_srli_epi32(a[kk], 13));

	// a *= m;
	UNROLL(XS) a[kk] = _mm_mullo_epi32(a[kk], m);

	// a ^= a >> 15;
	UNROLL(XS) a[kk] = _mm_xor_si128(a[kk], _mm_srli_epi32(a[kk], 15));

	// original
	UNROLL(XS) x[kk] = _mm_xor_si128(_mm_set1_epi32(0x075BCD15), a[kk]);

	UNROLL(XS) a[kk] = _mm_slli_epi32(a[kk], 3);

	UNROLL(XS) y[kk] = _mm_xor_si128(_mm_set1_epi32(0x159A55E5), a[kk]);

	UNROLL(XS) a[kk] = _mm_slli_epi32(a[kk], 3);

	UNROLL(XS) z[kk] = _mm_xor_si128(_mm_set1_epi32(0x1F123BB5), a[kk]);
}

/* Draw 64 64-bit integers in [0,n), 64*n frequencies are in arr */
TARGETSSE void RNGSSE<float>::Poly(__m128* arr, int n, __m128i* re)
{
	__m128 t[XS], s[XS];
	__m128 one = _mm_set1_ps(1.0f);
	__m128i mask1 = _mm_set1_epi32(0x007FFFFF);
	__m128i mask2 = _mm_set1_epi32(0x3F800000);
	__m128i* r = (__m128i*)t;

	UNROLL(XS) s[kk] = _mm_setzero_ps();

	for (int i = 0; i < n * XS; i += XS)
		UNROLL(XS) s[kk] = _mm_add_ps(s[kk], arr[kk + i]);

	XorShift();

	UNROLL(XS) r[kk] = _mm_and_si128(z[kk], mask1);

	UNROLL(XS) r[kk] = _mm_or_si128(r[kk], mask2);

	UNROLL(XS) t[kk] = _mm_sub_ps(t[kk], one);

	UNROLL(XS) t[kk] = _mm_mul_ps(t[kk], s[kk]);

	__m128i midx[XS], nidx = _mm_setzero_si128(), ninc = _mm_set1_epi32(1);
	__m128 f[XS], b[XS];
	UNROLL(XS) midx[kk] = _mm_set1_epi32(n - 1);
	UNROLL(XS) f[kk] = _mm_setzero_ps();

	for (int i = 0; i < n * XS; i += XS)
	{
		UNROLL(XS) b[kk] = _mm_cmplt_ps(t[kk], arr[kk + i]);

		UNROLL(XS) t[kk] = _mm_sub_ps(t[kk], arr[kk + i]);

		UNROLL(XS) b[kk] = _mm_andnot_ps(f[kk], b[kk]);

		UNROLL(XS) f[kk] = _mm_or_ps(f[kk], b[kk]);

		UNROLL(XS) midx[kk] = _mm_castps_si128(_mm_blendv_ps(_mm_castsi128_ps(midx[kk]), _mm_castsi128_ps(nidx), b[kk]));//ok

		nidx = _mm_add_epi32(nidx, ninc);
	}

	UNROLL(XS)
	{
		re[0 + (kk << 1)] = _mm_cvtepi32_epi64(midx[kk]);
		re[1 + (kk << 1)] = _mm_cvtepi32_epi64(_mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(midx[kk]), _mm_castsi128_ps(midx[kk]), _MM_SHUFFLE(1, 0, 3, 2))));
	}
}

/* Draw uniform distriubted intergers */
TARGETSSE void RNGSSE<float>::XorShift()
{
	__m128i u[XS];

	UNROLL(XS) u[kk] = _mm_slli_epi32(x[kk], 16);
	UNROLL(XS) x[kk] = _mm_xor_si128(x[kk], u[kk]);

	UNROLL(XS) u[kk] = _mm_srli_epi32(x[kk], 5);
	UNROLL(XS) x[kk] = _mm_xor_si128(x[kk], u[kk]);

	UNROLL(XS) u[kk] = _mm_slli_epi32(x[kk], 1);
	UNROLL(XS) x[kk] = _mm_xor_si128(x[kk], u[kk]);

	UNROLL(XS) u[kk] = x[kk];

	UNROLL(XS) x[kk] = y[kk];

	UNROLL(XS) y[kk] = z[kk];

	UNROLL(XS) z[kk] = _mm_xor_si128(u[kk], x[kk]);

	UNROLL(XS) z[kk] = _mm_xor_si128(z[kk], y[kk]);
}

/* Draw uniform distriubted integers */
template<typename INT>
TARGETSSE void RNGSSE<float>::Integer(INT* re, int64 n, INT minv, INT maxv)
{
	constexpr int xesize = sizeof(x) / sizeof(INT);
	int64 i = 0;
	INT modv = maxv - minv;
	INT* rei = re;

	for (; i <= n - xesize; i += xesize)
	{
		XorShift();
		UNROLL(XS) { _mm_storeu_si128((__m128i*)rei, z[kk]); rei += E128B / sizeof(INT); }
	}

	if (i != n)
	{
		XorShift();
		SetVal((INT*)rei, (INT*)z, n - i);
	}

	if (maxv != (INT)-1 || minv != 0)
	{
		for (i = 0; i < n; ++i)
			re[i] = re[i] % modv + minv;
	}
}

/* Draw uniform distriubted real numbers */
TARGETSSE void RNGSSE<float>::Uniform(float* re, int n, float minv, float maxv)
{
	constexpr int xesize = sizeof(x) / sizeof(float);
	int i = 0;
	float range = maxv - minv;

	__m128i mask1 = _mm_set1_epi32(0x007FFFFF);
	__m128i mask2 = _mm_set1_epi32(0x3F800000);
	__m128 v1 = _mm_set1_ps(minv - range);
	__m128 v2 = _mm_set1_ps(range);

	if (range == 1.0)
	{
		for (; i <= n - xesize; i += xesize)
		{
			XorShift();
			UNROLL(XS) { _mm_storeu_ps(re, _mm_add_ps(_mm_castsi128_ps(_mm_or_si128(_mm_and_si128(z[kk], mask1), mask2)), v1)); re += E128F; }
		}
	}
	else
	{
		for (; i <= n - xesize; i += xesize)
		{
			XorShift();
			UNROLL(XS) { _mm_storeu_ps(re, _mm_fmaddx_ps(_mm_castsi128_ps(_mm_or_si128(_mm_and_si128(z[kk], mask1), mask2)), v2, v1)); re += E128F; }
		}
	}

	if (i != n)
	{
		__m128 ref[XS];
		XorShift();
		UNROLL(XS) ref[kk] = _mm_fmaddx_ps(_mm_castsi128_ps(_mm_or_si128(_mm_and_si128(z[kk], mask1), mask2)), v2, v1);
		SetVal((float*)re, (float*)ref, n - i);
	}
}

/* Draw uniform distriubted real numbers */
TARGETSSE void RNGSSE<float>::Normal(float* re, int n, float mean, float sd)
{
	constexpr int xhsize = XS / 2;
	constexpr int xesize = XS * E128F;

	int i = 0;

	__m128i mask1 = _mm_set1_epi32(0x007FFFFF);
	__m128i mask2 = _mm_set1_epi32(0x3F800000);
	__m128 v1 = _mm_set1_ps(-1);
	__m128 min_freq = _mm_set1_ps((float)MIN_FREQ);
	__m128 pi2 = _mm_set1_ps((float)(2.0 * M_PI));
	__m128 mu = _mm_set1_ps(mean);
	__m128 s = _mm_set1_ps(sd);

	for (; i <= n - xesize; i += xesize)
	{
		XorShift();
		UNROLL(XS) _mm_storeu_ps(re + kk * E128F, _mm_add_ps(_mm_castsi128_ps(_mm_or_si128(_mm_and_si128(z[kk], mask1), mask2)), v1));

		__m128 u1, u2, u3, u4;
		for (int j = 0; j < xhsize; ++j)
		{
			u1 = _mm_max_ps(_mm_loadu_ps(re + j * E128F), min_freq);
			u2 = _mm_mul_ps(_mm_loadu_ps(re + j * E128F + xhsize * E128F), pi2);

			UNROLL(E128F) simd_f32(u1, kk) = sqrt(-2.0 * log(simd_f32(u1, kk)));
			UNROLL(E128F) simd_f32(u3, kk) = cos(simd_f32(u2, kk));
			UNROLL(E128F) simd_f32(u4, kk) = sin(simd_f32(u2, kk));
				
			_mm_storeu_ps(re + j * E128F, _mm_mul_ps(u1, u3));
			_mm_storeu_ps(re + j * E128F + xhsize * E128F, _mm_mul_ps(u1, u4));
		}
			
		if (sd != 1 || mean != 0)
			UNROLL(XS) { _mm_storeu_ps(re, _mm_fmaddx_ps(_mm_loadu_ps(re), s, mu)); re += E128F; }
		else
			re += XS * E128F;
	}

	if (i != n)
	{
		float ref[XS * E128F];
		
		XorShift();
		UNROLL(XS) _mm_storeu_ps(ref + kk * E128F, _mm_add_ps(_mm_castsi128_ps(_mm_or_si128(_mm_and_si128(z[kk], mask1), mask2)), v1));

		__m128 u1, u2, u3, u4;
		for (int j = 0; j < xhsize; ++j)
		{
			u1 = _mm_max_ps(_mm_loadu_ps(ref + j * E128F), min_freq);
			u2 = _mm_mul_ps(_mm_loadu_ps(ref + j * E128F + xhsize * E128F), pi2);

			UNROLL(E128F) simd_f32(u1, kk) = sqrt(-2.0 * log(simd_f32(u1, kk)));
			UNROLL(E128F) simd_f32(u3, kk) = cos(simd_f32(u2, kk));
			UNROLL(E128F) simd_f32(u4, kk) = sin(simd_f32(u2, kk));
				
			_mm_storeu_ps(ref + j * E128F, _mm_mul_ps(u1, u3));
			_mm_storeu_ps(ref + j * E128F + xhsize * E128F, _mm_mul_ps(u1, u4));
		}
			
		if (sd != 1 || mean != 0)
			UNROLL(XS) _mm_storeu_ps(ref + kk * E128F, _mm_fmaddx_ps(_mm_loadu_ps(re), s, mu)); 

		SetVal((float*)re, (float*)ref, n - i);
	}
}
#endif

TARGETSSE int64 GetMinIdxSSE(double* A, int64 n, double& val)
{
#define N 4
	int64 i = 0;
	val = DBL_MAX;
	uint64 idx = (uint64)-1;

	if (n >= N * E128D)
	{
		__m128d min1[N], f[N];
		__m128i midx[N], nidx[N], msep = _mm_set1_epi64x(N * E128D);
		UNROLL(N) min1[kk] = _mm_set1_pd(val);
		UNROLL(N) midx[kk] = _mm_set1_epi64x(0xFFFFFFFFFFFFFFFF);
		UNROLL(N) nidx[kk] = _mm_set_epi64x(1 + (kk << 1), 0 + (kk << 1));

		for (int64 l1 = n - N * E128D; i <= l1; i += N * E128D)
		{
			UNROLL(N) 
			{	
				f[kk] = _mm_cmpgt_pd(min1[kk], _mm_loadu_pd(A)); 
				min1[kk] = _mm_blendv_pd(min1[kk], _mm_loadu_pd(A), f[kk]); 
				A += E128D; 
			}

			UNROLL(N) midx[kk] = _mm_castpd_si128(_mm_blendv_pd(_mm_castsi128_pd(midx[kk]), _mm_castsi128_pd(nidx[kk]), f[kk]));//ok

			UNROLL(N) nidx[kk] = _mm_add_epi64(nidx[kk], msep);
		}

		REDUCE(min1)
		{
			f[kk] = _mm_cmpgt_pd(min1[kk], min1[kk + KK]);
			min1[kk] = _mm_min_pd(min1[kk], min1[kk + KK]);
			midx[kk] = _mm_castpd_si128(_mm_blendv_pd(_mm_castsi128_pd(midx[kk]), _mm_castsi128_pd(midx[kk + KK]), f[kk]));
		}

		for (int64 j = 0; j < E128D; ++j)
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
#define N 4
	int64 i = 0;
	val = FLT_MAX;
	uint idx = (uint)-1;

	if (n >= N * E128F)
	{
		__m128 min1[N], f[N];
		__m128i midx[N], nidx[N], msep = _mm_set1_epi32(N * E128F);
		UNROLL(N) min1[kk] = _mm_set1_ps(val);
		UNROLL(N) midx[kk] = _mm_set1_epi32(0xFFFFFFFF);
		UNROLL(N) nidx[kk] = _mm_set_epi32(3 + (kk << 2), 2 + (kk << 2), 1 + (kk << 2), 0 + (kk << 2));

		for (int64 l1 = n - N * E128F; i <= l1; i += N * E128F)
		{
			UNROLL(N) 
			{ 
				f[kk] = _mm_cmpgt_ps(min1[kk], _mm_loadu_ps(A));
				min1[kk] = _mm_blendv_ps(min1[kk], _mm_loadu_ps(A), f[kk]); 
				A += E128F; 
			}

			UNROLL(N) midx[kk] = _mm_castps_si128(_mm_blendv_ps(_mm_castsi128_ps(midx[kk]), _mm_castsi128_ps(nidx[kk]), f[kk]));

			UNROLL(N) nidx[kk] = _mm_add_epi32(nidx[kk], msep);
		}

		REDUCE(min1)
		{
			f[kk] = _mm_cmpgt_ps(min1[kk], min1[kk + KK]);
			min1[kk] = _mm_blendv_ps(min1[kk], min1[kk + KK], f[kk]);
			midx[kk] = _mm_castps_si128(_mm_blendv_ps(_mm_castsi128_ps(midx[kk]), _mm_castsi128_ps(midx[kk + KK]), f[kk]));
		}

		for (int64 j = 0; j < E128F; ++j)
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
#define N 4
	int64 i = 0;
	minv = DBL_MAX;
	maxv = -DBL_MAX;

	if (n >= N * E128D)
	{
		__m128d min1[N], max1[N];
		UNROLL(N) min1[kk] = _mm_set1_pd(minv);
		UNROLL(N) max1[kk] = _mm_set1_pd(maxv);

		for (int64 l1 = n - N * E128D; i <= l1; i += N * E128D)
		{
			UNROLL(N)
			{
				min1[kk] = _mm_min_pd(min1[kk], _mm_loadu_pd(A));
				max1[kk] = _mm_max_pd(max1[kk], _mm_loadu_pd(A));
				A += E128D;
			}
		}

		REDUCE(min1)
		{
			min1[kk] = _mm_min_pd(min1[kk], min1[kk + KK]);
			max1[kk] = _mm_max_pd(max1[kk], max1[kk + KK]);
		}
		
		minv = _mm_reduce_min_pd(min1[0]);
		maxv = _mm_reduce_max_pd(max1[0]);
	}

	for (; i < n; ++i, ++A)
	{
		minv = std::min(minv, *A);
		maxv = std::max(maxv, *A);
	}
}

TARGETSSE void GetMinMaxValSSE(float* A, int64 n, float& minv, float& maxv)
{
#define N 4
	int64 i = 0;
	minv = FLT_MAX;
	maxv = -FLT_MAX;

	if (n >= N * E128F)
	{
		__m128 min1[N], max1[N];
		UNROLL(N) min1[kk] = _mm_set1_ps(minv);
		UNROLL(N) max1[kk] = _mm_set1_ps(maxv);

		for (int64 l1 = n - N * E128F; i <= l1; i += N * E128F)
		{
			UNROLL(N)
			{
				min1[kk] = _mm_min_ps(min1[kk], _mm_loadu_ps(A));
				max1[kk] = _mm_max_ps(max1[kk], _mm_loadu_ps(A));
				A += E128F;
			}
		}

		REDUCE(min1)
		{
			min1[kk] = _mm_min_ps(min1[kk], min1[kk + KK]);
			max1[kk] = _mm_max_ps(max1[kk], max1[kk + KK]);
		}
		
		minv = _mm_reduce_min_ps(min1[0]);
		maxv = _mm_reduce_max_ps(max1[0]);
	}

	for (; i < n; ++i, ++A)
	{
		minv = std::min(minv, *A);
		maxv = std::max(maxv, *A);
	}
}

TARGETSSE double GetMaxValSSE(double* A, int64 n)
{
#define N 4
	int64 i = 0;
	double val = -DBL_MAX;

	if (n >= N * E128D)
	{
		__m128d max1[N];
		UNROLL(N) max1[kk] = _mm_set1_pd(val);

		for (int64 l1 = n - N * E128D; i <= l1; i += N * E128D)
		{
			UNROLL(N) { max1[kk] = _mm_max_pd(max1[kk], _mm_loadu_pd(A)); A += E128D; }
		}

		REDUCE(max1) max1[kk] = _mm_max_pd(max1[kk], max1[kk + KK]);
		
		val = _mm_reduce_max_pd(max1[0]);
	}

	for (; i < n; ++i, ++A)
	{
		val = std::max(val, *A);
	}

	return val;
}

TARGETSSE float GetMaxValSSE(float* A, int64 n)
{
#define N 8
	int64 i = 0;
	float val = -FLT_MAX;

	if (n >= N * E128F)
	{
		__m128 max1[N];
		UNROLL(N) max1[kk] = _mm_set1_ps(val);

		for (int64 l1 = n - N * E128F; i <= l1; i += N * E128F)
		{
			UNROLL(N) { max1[kk] = _mm_max_ps(max1[kk], _mm_loadu_ps(A)); A += E128F; }
		}

		REDUCE(max1) max1[kk] = _mm_max_ps(max1[kk], max1[kk + KK]);
		
		val = _mm_reduce_max_ps(max1[0]);
	}

	for (; i < n; ++i, ++A)
	{
		val = std::max(val, *A);
	}

	return val;
}

TARGETSSE double GetMaxValSSE(double* A, int64 n, int64 sep)
{
	//suboptimal to compile
	{
		double val = -DBL_MAX;
		UNROLLHEAD(4)
		for (int64 i = 0; i < n; ++i, A += sep)
			val = std::max(val, *A);
		return val;
	}

#define N 4
	int64 i = 0;
	double val = -DBL_MAX;

	if (n >= N * E128D)
	{
		__m128d max1[N];
		UNROLL(N) max1[kk] = _mm_set1_pd(val);
		
		for (int64 l1 = n - N * E128D; i <= l1; i += N * E128D)
		{
			UNROLL(N)
			{
				max1[kk] = _mm_max_pd(max1[kk], _mm_set_pd(A[1 * sep], A[0 * sep]));
				A += E128D * sep;
			}
		}

		REDUCE(max1) max1[kk] = _mm_max_pd(max1[kk], max1[kk + KK]);

		val = _mm_reduce_max_pd(max1[0]);
	}

	for (; i < n; ++i, A += sep)
	{
		val = std::max(val, *A);
	}

	return val;
}

TARGETSSE float GetMaxValSSE(float* A, int64 n, int64 sep)
{
	//suboptimal to compile
	{
		float val = -FLT_MAX;
		UNROLLHEAD(4)
		for (int64 i = 0; i < n; ++i, A += sep)
			val = std::max(val, *A);
		return val;
	}

#define N 4
	int64 i = 0;
	float val = -FLT_MAX;

	if (n >= N * E128F)
	{
		__m128 max1[N];
		UNROLL(N) max1[kk] = _mm_set1_ps(val);

		for (int64 l1 = n - N * E128F; i <= l1; i += N * E128F)
		{
			UNROLL(N)
			{
				max1[kk] = _mm_max_ps(max1[kk], _mm_set_ps(A[3 * sep], A[2 * sep], A[1 * sep], A[0 * sep]));
				A += E128F * sep;
			}
		}

		REDUCE(max1) max1[kk] = _mm_max_ps(max1[kk], max1[kk + KK]);

		val = _mm_reduce_max_ps(max1[0]);
	}

	for (; i < n; ++i, A += sep)
	{
		val = std::max(val, *A);
	}

	return val;
}

TARGETSSE double GetMinValSSE(double* A, int64 n)
{
#define N 4
	int64 i = 0;
	double val = DBL_MAX;

	if (n >= N * E128D)
	{
		__m128d min1[N];
		UNROLL(N) min1[kk] = _mm_set1_pd(val);

		for (int64 l1 = n - N * E128D; i <= l1; i += N * E128D)
		{
			UNROLL(N) { min1[kk] = _mm_min_pd(min1[kk], _mm_loadu_pd(A)); A += E128D; }
		}

		REDUCE(min1) min1[kk] = _mm_min_pd(min1[kk], min1[kk + KK]);
		
		val = _mm_reduce_min_pd(min1[0]);
	}

	for (; i < n; ++i, ++A)
	{
		val = std::min(val, *A);
	}

	return val;
}

TARGETSSE float GetMinValSSE(float* A, int64 n)
{
#define N 8
	int64 i = 0;
	float val = FLT_MAX;

	if (n >= N * E128F)
	{
		__m128 min1[N], a[N];
		UNROLL(N) min1[kk] = _mm_set1_ps(val);

		for (int64 l1 = n - N * E128F; i <= l1; i += N * E128F)
		{
			UNROLL(N) { min1[kk] = _mm_min_ps(min1[kk], _mm_loadu_ps(A)); A += E128F; }
		}

		REDUCE(min1) min1[kk] = _mm_min_ps(min1[kk], min1[kk + KK]);
		
		val = _mm_reduce_min_ps(min1[0]);
	}

	for (; i < n; ++i, ++A)
	{
		val = std::min(val, *A);
	}

	return val;
}

TARGETSSE int64 GetMinValSSE(int64* A, int64 n)
{
#define N 4
	int64 i = 0;
	int64 val = 0x7FFFFFFFFFFFFFFF;

	if (n >= N * E128D)
	{
		__m128i min1[N];
		UNROLL(N) min1[kk] = _mm_set1_epi64x(0x7FFFFFFFFFFFFFFF);

		for (int64 l1 = n - N * E128D; i <= l1; i += N * E128D)
		{
			UNROLL(N) 
			{ 
				min1[kk] = _mm_castpd_si128(_mm_blendv_pd(_mm_castsi128_pd(min1[kk]), _mm_castsi128_pd(_mm_loadu_si128((__m128i*)A)), _mm_castsi128_pd(_mm_cmpgt_epi64(min1[kk], _mm_loadu_si128((__m128i*)A))))); 

				A += E128D; 
			}
		}

		REDUCE(min1) min1[kk] = _mm_castpd_si128(_mm_blendv_pd( _mm_castsi128_pd(min1[kk]), _mm_castsi128_pd(min1[kk + KK]), _mm_castsi128_pd(_mm_cmpgt_epi64(min1[kk], min1[kk + KK]))));
		
		val = _mm_reduce_min_epi64(min1[0]);
	}

	for (; i < n; ++i, ++A)
	{
		val = std::min(val, *A);
	}

	return val;
}

TARGETSSE void SetValSSE(uint* A, ushort* B, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E128F)
	{
		for (int64 l1 = n - N * E128F; i <= l1; i += N * E128F)
		{
			UNROLL(N) { _mm_storeu_si128((__m128i*)A, _mm_cvtepu16_epi32(_mm_loadu_si64(B))); A += E128F; B += E128F;}
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
#define N 4
	int64 i = 0;
	int64 slog = 0; double prod = 1;

	if (n >= N * E512D)
	{
		__m128d pd[E512_128];
		UNROLL(E512_128) pd[kk] = _mm_set1_pd(1.0);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128)
			{ pd[kk] = _mm_mul_pd(pd[kk], _mm_loadu_pd(A));  A += E128D; }

			UNROLL(E512_128) AddExponentSSE(slog, pd[kk]);
		}

		__m128d* pd1 = (__m128d*)pd;
		UNROLL(E512_128) ChargeLogSSE(slog, prod, pd1[kk]);
	}

	for (; i < n; ++i, ++A)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);
	return prod;
}

TARGETSSE double LogProdSSE(float* A, int64 n)
{
#define N 4
	int64 i = 0;
	int64 slog = 0; double prod = 1;

	if (n >= N * E512D)
	{
		__m128d pd[E512_128];
		UNROLL(E512_128) pd[kk] = _mm_set1_pd(1.0);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128)
			{ pd[kk] = _mm_mul_pd(pd[kk], _mm_cvtps_pd(_mm_castsi128_ps(_mm_loadu_si64(A)))); A += E128D; }
				
			UNROLL(E512_128) AddExponentSSE(slog, pd[kk]);
		}
		
		__m128d* pd1 = (__m128d*)pd;
		UNROLL(E512_128) ChargeLogSSE(slog, prod, pd1[kk]);
	}

	for (; i < n; ++i, ++A)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);
	return prod;
}

TARGETSSE double LogProdSSE(double* A, int64 n, int64 sep)
{
	//suboptimal to compile
	{
		int64 slog = 0; double prod = 1;
		UNROLLHEAD(4)
		for (int64 i = 0; i < n; ++i, A += sep)
			ChargeLog(slog, prod, *A);

		CloseLog(slog, prod);
		return prod;
	}

#define N 4
	int64 i = 0;
	int64 slog = 0; double prod = 1;

	if (n >= N * E512D)
	{
		__m128d pd[E512_128];
		UNROLL(E512_128) pd[kk] = _mm_set1_pd(1.0);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128) 
			{ pd[kk] = _mm_mul_pd(pd[kk], _mm_set_pd(A[1 * sep], A[0 * sep])); A += E128D * sep; }

			UNROLL(E512_128) AddExponentSSE(slog, pd[kk]);
		}

		UNROLL(E512_128) ChargeLogSSE(slog, prod, pd[kk]);
	}

	for (; i < n; ++i, A += sep)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);

	return prod;
}

TARGETSSE double LogProdSSE(float* A, int64 n, int64 sep)
{
	//suboptimal to compile
	{
		int64 slog = 0; double prod = 1;
		UNROLLHEAD(4)
		for (int64 i = 0; i < n; ++i, A += sep)
			ChargeLog(slog, prod, *A);

		CloseLog(slog, prod);
		return prod;
	}

#define N 4
	int64 i = 0;
	int64 slog = 0; double prod = 1;

	if (n >= N * E128F)
	{
		__m128d pd1 = _mm_set1_pd(1.0), pd2 = _mm_set1_pd(1.0);
		__m128 a1, a2;

		for (int64 l1 = n - N * E128F; i <= l1; i += N * E128F)
		{
			UNROLL(N)
			{
				a1 = _mm_set_ps(A[3 * sep], A[2 * sep], A[1 * sep], A[0 * sep]); A += E128F * sep;
				a2 = _mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(a1), 8));
				pd1 = _mm_mul_pd(pd1, _mm_cvtps_pd(a1));
				pd2 = _mm_mul_pd(pd2, _mm_cvtps_pd(a2));
			}

			AddExponentSSE(slog, pd1);
			AddExponentSSE(slog, pd2);
		}

		ChargeLogSSE(slog, prod, pd1);
		ChargeLogSSE(slog, prod, pd2);
	}

	for (; i < n; ++i, A += sep)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);
	return prod;
}

TARGETSSE double LogProdDivSSE(double* A, double* B, int64 n, int64 sep)
{
	//suboptimal to compile
	{
		int64 slog1 = 0; double prod1 = 1;
		int64 slog2 = 0; double prod2 = 1;
		UNROLLHEAD(4)
		for (int64 i = 0; i < n; ++i, A += sep, B += sep)
		{
			ChargeLog(slog1, prod1, *A);
			ChargeLog(slog2, prod2, *B);
		}

		CloseLog(slog1, prod1);
		CloseLog(slog2, prod2);
		return prod1 - prod2;
	}
}

TARGETSSE double LogProdDivSSE(float* A, float* B, int64 n, int64 sep)
{
	//suboptimal to compile
	{
		int64 slog1 = 0; double prod1 = 1;
		int64 slog2 = 0; double prod2 = 1;
		UNROLLHEAD(4)
		for (int64 i = 0; i < n; ++i, A += sep, B += sep)
		{
			ChargeLog(slog1, prod1, *A);
			ChargeLog(slog2, prod2, *B);
		}

		CloseLog(slog1, prod1);
		CloseLog(slog2, prod2);
		return prod1 - prod2;
	}

#define N 8
	int64 i = 0;
	int64 slog1 = 0; double prod1 = 1;
	int64 slog2 = 0; double prod2 = 1;

	if (n >= N * E512D)
	{
		__m128d pd1[E512_128], pd2[E512_128];
		UNROLL(E512_128) pd1[kk] = pd2[kk] = _mm_set1_pd(1.0);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128)
			{
				pd1[kk] = _mm_mul_pd(pd1[kk], _mm_set_pd(A[1 * sep], A[0 * sep])); A += sep * E128D;
				pd2[kk] = _mm_mul_pd(pd2[kk], _mm_set_pd(B[1 * sep], B[0 * sep])); B += sep * E128D;
			}

			UNROLL(E512_128) AddExponentSSE(slog1, pd1[kk]);
			UNROLL(E512_128) AddExponentSSE(slog2, pd2[kk]);
		}

		UNROLL(E512_128) ChargeLogSSE(slog1, prod1, pd1[kk]);
		UNROLL(E512_128) ChargeLogSSE(slog2, prod2, pd2[kk]);
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
	
	if (n >= E128B)
	{
		__m128i a = _mm_setzero_si128(), z = _mm_setzero_si128(), o = _mm_set1_epi8(0x01);

		for (int64 l1 = n - E128B; i <= l1; i += E128B)
		{
			a = _mm_add_epi64(a, _mm_sad_epu8(z, _mm_min_epu8(_mm_loadu_si128((__m128i*)A), o))); A += E128B;
		}

		re = _mm_reduce_add_epi64(a);
	}

	for (; i < n; ++i, ++A)
		if (*A) re++;

	return (int64)re;
}

TARGETSSE double SumSSE(double* A, int64 n)
{
#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		__m128d s[E512_128] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128)
			{ s[kk] = _mm_add_pd(s[kk], _mm_loadu_pd(A)); A += E128D; }
		}

		REDUCE(s) s[kk] = _mm_add_pd(s[kk], s[kk + KK]);

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
#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		__m128d s[E512_128] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128)
			{ s[kk] = _mm_add_pd(s[kk], _mm_set_pd(A[1], A[0])); A += 2; }
		}

		REDUCE(s) s[kk] = _mm_add_pd(s[kk], s[kk + KK]);

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
#define N 4
	int64 i = 0;
	volatile float re = 0;

	if (n >= N * E512F)
	{
		__m128 s[E512_128] = { 0 };

		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) UNROLL(E512_128)
			{ s[kk] = _mm_add_ps(s[kk], _mm_loadu_ps(A)); A += E128F; }
		}

		REDUCE(s) s[kk] = _mm_add_ps(s[kk], s[kk + KK]);

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
#define N 4
	int64 re = 0;
	int64 i = 0;

	if (n >= N * E128B)
	{
		__m128i s = _mm_setzero_si128(), z = _mm_setzero_si128();

		for (int64 l1 = n - N * E128B; i <= l1; i += N * E128B)
		{
			UNROLL(N)
			{ s = _mm_add_epi64(s, _mm_sad_epu8(_mm_loadu_si128((__m128i*)A), z)); A += E128B; }
		}
		
		re += _mm_reduce_add_epi64(s);
	}

	for (; i < n; ++i)
		re += *A++;

	return re;
}

TARGETSSE double SumSSE(double* A, int64 n, int64 sep)
{
	//suboptimal to compile
	{
		double re = 0;
		UNROLLHEAD(4)
		for (int64 i = 0; i < n; ++i, A += sep)
			re += *A;
		return re;
	}

#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		__m128d s[E512_128] = { 0 };
		UNROLL(E128D) _mm_prefetch((const char*)&A[kk * sep], _MM_HINT_T0);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128)
			{
				UNROLL(E128D) _mm_prefetch((const char*)&A[E128D * sep + kk * sep], _MM_HINT_T0);
				s[kk] = _mm_add_pd(s[kk], _mm_set_pd(A[1 * sep], A[0 * sep]));
				A += E128D * sep;
			}
		}

		REDUCE(s) s[kk] = _mm_add_pd(s[kk], s[kk + KK]);

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
	//suboptimal to compile
	{
		double re = 0;
		UNROLLHEAD(4)
		for (int64 i = 0; i < n; ++i, A += sep)
			re += *A;
		return re;
	}

#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		__m128d s[E512_128] = { 0 };
		UNROLL(E128D) _mm_prefetch((const char*)&A[kk * sep], _MM_HINT_T0);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128)
			{
				UNROLL(E128D) _mm_prefetch((const char*)&A[E128D * sep + kk * sep], _MM_HINT_T0);
				s[kk] = _mm_add_pd(s[kk], _mm_set_pd(A[1 * sep], A[0 * sep])); 
				A += E128D * sep;
			}
		}

		REDUCE(s) s[kk] = _mm_add_pd(s[kk], s[kk + KK]);

		re = _mm_reduce_add_pd(s[0]);
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
	//suboptimal to compile
	{
		float re = 0;
		UNROLLHEAD(4)
		for (int64 i = 0; i < n; ++i, A += sep)
			re += *A;
		return re;
	}

#define N 4
	int64 i = 0;
	volatile float re = 0;

	if (n >= N * E512F)
	{
		__m128 s[E512_128] = { 0 };
		UNROLL(E128F) _mm_prefetch((const char*)&A[kk * sep], _MM_HINT_T0);

		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) UNROLL(E512_128)
			{
				UNROLL(E128F) _mm_prefetch((const char*)&A[E128F * sep + kk * sep], _MM_HINT_T0);
				s[kk] = _mm_add_ps(s[kk], _mm_set_ps(A[3 * sep], A[2 * sep], A[1 * sep], A[0 * sep]));
				A += E128F * sep; 
			}
		}

		REDUCE(s) s[kk] = _mm_add_ps(s[kk], s[kk + KK]);

		re = _mm_reduce_add_ps(s[0]);
	}

	for (; i < n; ++i, A += sep)
	{
		volatile float v1 = *A;
		re += v1;
	}

	return re;
}

TARGETSSE double ProdSSE(double* A, int64 n)
{
#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		__m128d s[E512_128];
		UNROLL(E512_128) s[kk] = _mm_set1_pd(1);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128)
			{ s[kk] = _mm_mul_pd(s[kk], _mm_loadu_pd(A)); A += E128D; }
		}

		REDUCE(s) s[kk] = _mm_mul_pd(s[kk], s[kk + KK]);

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
#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		__m128d s[E512_128];
		UNROLL(E512_128) s[kk] = _mm_set1_pd(1);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128)
			{ s[kk] = _mm_mul_pd(s[kk], _mm_set_pd(A[1], A[0])); A += E128D; }
		}

		REDUCE(s) s[kk] = _mm_mul_pd(s[kk], s[kk + KK]);

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
#define N 4
	int64 i = 0;
	volatile float re = 0;

	if (n >= N * E512F)
	{
		__m128 s[E512_128];
		UNROLL(E512_128) s[kk] = _mm_set1_ps(1);

		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) UNROLL(E512_128) 
			{ s[kk] = _mm_mul_ps(s[kk], _mm_loadu_ps(A)); A += E128F; }
		}

		REDUCE(s) s[kk] = _mm_mul_ps(s[kk], s[kk + KK]);

		re = _mm_reduce_mul_ps(s[0]);
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A++;
		re *= v1;
	}

	return re;
}

TARGETSSE double ProdSSE(double* A, int64 n, int64 sep)
{
	//suboptimal to compile
	{
		double re = 1;
		UNROLLHEAD(4)
		for (int64 i = 0; i < n; ++i, A += sep)
			re *= *A;
		return re;
	}

#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		__m128d s[E512_128];
		UNROLL(E512_128) s[kk] = _mm_set1_pd(1);
		UNROLL(E128D) _mm_prefetch((const char*)&A[kk * sep], _MM_HINT_T0);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128)
			{
				UNROLL(E128D) _mm_prefetch((const char*)&A[E128D * sep + kk * sep], _MM_HINT_T0);
				s[kk] = _mm_mul_pd(s[kk], _mm_set_pd(A[1 * sep], A[0 * sep]));
				A += E128D * sep;
			}
		}

		REDUCE(s) s[kk] = _mm_mul_pd(s[kk], s[kk + KK]);

		re = _mm_reduce_mul_pd(s[0]);
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
	//suboptimal to compile
	{
		double re = 1;
		UNROLLHEAD(4)
		for (int64 i = 0; i < n; ++i, A += sep)
			re *= *A;
		return re;
	}

#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		__m128d s[E512_128];
		UNROLL(E512_128) s[kk] = _mm_set1_pd(1);
		UNROLL(E128D) _mm_prefetch((const char*)&A[kk * sep], _MM_HINT_T0);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128)
			{
				UNROLL(E128D) _mm_prefetch((const char*)&A[E128D * sep + kk * sep], _MM_HINT_T0);
				s[kk] = _mm_mul_pd(s[kk], _mm_set_pd(A[1 * sep], A[0 * sep]));
				A += E128D * sep;
			}
		}

		REDUCE(s) s[kk] = _mm_mul_pd(s[kk], s[kk + KK]);

		re = _mm_reduce_mul_pd(s[0]);
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
	//suboptimal to compile
	{
		float re = 1;
		UNROLLHEAD(4)
		for (int64 i = 0; i < n; ++i, A += sep)
			re *= *A;
		return re;
	}

#define N 4
	int64 i = 0;
	volatile float re = 0;

	if (n >= N * E512F)
	{
		__m128 s[E512_128];
		UNROLL(E512_128) s[kk] = _mm_set1_ps(1);
		UNROLL(E128F) _mm_prefetch((const char*)&A[kk * sep], _MM_HINT_T0);

		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) UNROLL(E512_128)
			{
				UNROLL(E128F) _mm_prefetch((const char*)&A[E128F * sep + kk * sep], _MM_HINT_T0);
				s[kk] = _mm_mul_ps(s[kk], _mm_set_ps(A[3 * sep], A[2 * sep], A[1 * sep], A[0 * sep]));
				A += E128F * sep;
			}
		}

		REDUCE(s) s[kk] = _mm_mul_ps(s[kk], s[kk + KK]);

		re = _mm_reduce_mul_ps(s[0]);
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
#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		__m128d s[E512_128] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128)
			{ s[kk] = _mm_fmaddx_pd(_mm_loadu_pd(A), _mm_loadu_pd(A), s[kk]); A += E128D; }
		}

		REDUCE(s) s[kk] = _mm_add_pd(s[kk], s[kk + KK]);

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
#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		__m128d s[E512_128] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128)
			{
				__m128d v1 = _mm_cvtps_pd(_mm_castsi128_ps(_mm_loadu_si64(A))); A += E128D;
				s[kk] = _mm_fmaddx_pd(v1, v1, s[kk]);
			}
		}

		REDUCE(s) s[kk] = _mm_add_pd(s[kk], s[kk + KK]);

		re = _mm_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, ++A)
	{
		volatile double v1 = (double)*A * (double)*A;
		re += v1;
	}

	return re;
}

TARGETSSE int64 SumSquareSSE(byte* A, int64 n)
{
#define N 4
	int64 i = 0;
	uint64 re = 0;

	if (n >= N * E128B)
	{
		__m128i t = _mm_setzero_si128(), s = _mm_setzero_si128();

		for (int64 l1 = n - N * E128B; i <= l1; i += N * E128B)
		{
			UNROLL(N) { s = _mm_add_epi16(s, _mm_maddubs_epi16(_mm_loadu_si128((__m128i*)A), _mm_loadu_si128((__m128i*)A))); A += E128B; }

			if ((i & (E128B * 128 - 1)) == 0) [[unlikely]]
			{
				UNROLL(E512_128) { t = _mm_add_epi64(t, _mm_cvtepi16_epi64(s)); s = _mm_srli_si128(s, 4); }
			}
		}
		
		UNROLL(E512_128) { t = _mm_add_epi64(t, _mm_cvtepi16_epi64(s)); s = _mm_srli_si128(s, 4); }
		re = _mm_reduce_add_epi64(t);
	}

	for (; i < n; ++i, ++A)
		re += *A * *A;

	return re;
}

TARGETSSE void SumSumSquareSSE(double* A, int64 n, double& sum, double& sumsq)
{
#define N 4
	int64 i = 0;
	volatile double re1 = 0, re2 = 0;

	if (n >= N * E512D)
	{
		__m128d s1[E512_128] = { 0 }, s2[E512_128] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128)
			{
				__m128d v1 = _mm_loadu_pd(A); A += E128D;

				s1[kk] = _mm_add_pd(v1, s1[kk]);

				s2[kk] = _mm_fmaddx_pd(v1, v1, s2[kk]);
			}
		}

		REDUCE(s1)
		{
			s1[kk] = _mm_add_pd(s1[kk], s1[kk + KK]);
			s2[kk] = _mm_add_pd(s2[kk], s2[kk + KK]);
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
#define N 4
	int64 i = 0;
	volatile double re1 = 0, re2 = 0;

	if (n >= N * E512D)
	{
		__m128d s1[E512_128] = { 0 }, s2[E512_128] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128)
			{
				__m128d v1 = _mm_cvtps_pd(_mm_castsi128_ps(_mm_loadu_si64(A))); A += E128D;

				s1[kk] = _mm_add_pd(v1, s1[kk]);

				s2[kk] = _mm_fmaddx_pd(v1, v1, s2[kk]);
			}
		}

		REDUCE(s1)
		{
			s1[kk] = _mm_add_pd(s1[kk], s1[kk + KK]);
			s2[kk] = _mm_add_pd(s2[kk], s2[kk + KK]);
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
	//suboptimal to compile
	{
		double re1 = 0, re2 = 0;
		UNROLLHEAD(4)
		for (int64 i = 0; i < n; ++i, A1++, A2++, B += sep)
		{
			volatile double v1 = (double)*A1 * (double)*B;
			volatile double v2 = (double)*A2 * (double)*B;
			re1 += v1;
			re2 += v2;
		}
		return re1 / re2;
	}

#define N 2
	int64 i = 0;
	volatile double re1 = 0, re2 = 0;

	if (n >= N * E512D)
	{
		UNROLL(E128D) { _mm_prefetch((const char*)B, _MM_HINT_T0); B += sep; }

		__m128d s1[E512_128] = { 0 }, s2[E512_128] = { 0 };
		
		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128)
			{
				UNROLL(E128D) { _mm_prefetch((const char*)B, _MM_HINT_T0); B += sep; }

				__m128d b = _mm_set_pd(B[-3 * sep], B[-4 * sep]);
				
				s1[kk] = _mm_fmaddx_pd(b, _mm_loadu_pd(A1), s1[kk]); A1 += E128D;

				s2[kk] = _mm_fmaddx_pd(b, _mm_loadu_pd(A2), s2[kk]); A2 += E128D;
			}
		}

		REDUCE(s1)
		{
			s1[kk] = _mm_add_pd(s1[kk], s1[kk + KK]);
			s2[kk] = _mm_add_pd(s2[kk], s2[kk + KK]);
		}

		re1 = _mm_reduce_add_pd(s1[0]);
		re2 = _mm_reduce_add_pd(s2[0]);

		B -= E128D * sep;
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

TARGETSSE double SumProdDivSSE(double* A1, float* A2, float* B, int64 sep, int64 n)
{
	//suboptimal to compile
	{
		double re1 = 0, re2 = 0;
		UNROLLHEAD(4)
		for (int64 i = 0; i < n; ++i, A1++, A2++, B += sep)
		{
			volatile double v1 = (double)*A1 * (double)*B;
			volatile double v2 = (double)*A2 * (double)*B;
			re1 += v1;
			re2 += v2;
		}
		return re1 / re2;
	}

#define N 2
	int64 i = 0;
	volatile double re1 = 0, re2 = 0;

	if (n >= N * E512D)
	{
		UNROLL(E128D) { _mm_prefetch((const char*)B, _MM_HINT_T0); B += sep; }

		__m128d s1[E512_128] = { 0 }, s2[E512_128] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128)
			{
				UNROLL(E128D) { _mm_prefetch((const char*)B, _MM_HINT_T0); B += sep; }

				__m128d b = _mm_set_pd(B[-3 * sep], B[-4 * sep]);

				s1[kk] = _mm_fmaddx_pd(_mm_loadu_pd(A1), b, s1[kk]); A1 += E128D;

				s2[kk] = _mm_fmaddx_pd(_mm_cvtps_pd(_mm_castsi128_ps(_mm_loadu_si64(A2))), b, s2[kk]); A2 += E128D;
			}
		}

		REDUCE(s1)
		{
			s1[kk] = _mm_add_pd(s1[kk], s1[kk + KK]);
			s2[kk] = _mm_add_pd(s2[kk], s2[kk + KK]);
		}

		re1 = _mm_reduce_add_pd(s1[0]);
		re2 = _mm_reduce_add_pd(s2[0]);

		B -= E128D * sep;
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
	//suboptimal to compile
	{
		double re1 = 0, re2 = 0;
		UNROLLHEAD(4)
		for (int64 i = 0; i < n; ++i, A1++, A2++, B += sep)
		{
			volatile double v1 = (double)*A1 * (double)*B;
			volatile double v2 = (double)*A2 * (double)*B;
			re1 += v1;
			re2 += v2;
		}
		return re1 / re2;
	}

#define N 2
	int64 i = 0;
	volatile double re1 = 0, re2 = 0;

	if (n >= N * E512D)
	{
		UNROLL(E128D) { _mm_prefetch((const char*)B, _MM_HINT_T0); B += sep; }

		__m128d s1[E512_128] = { 0 }, s2[E512_128] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128)
			{
				UNROLL(E128D) { _mm_prefetch((const char*)B, _MM_HINT_T0); B += sep; }

				__m128d b = _mm_set_pd(B[-3 * sep], B[-4 * sep]);

				s1[kk] = _mm_fmaddx_pd(_mm_cvtps_pd(_mm_castsi128_ps(_mm_loadu_si64(A1))), b, s1[kk]); A1 += E128D;

				s2[kk] = _mm_fmaddx_pd(_mm_cvtps_pd(_mm_castsi128_ps(_mm_loadu_si64(A2))), b, s2[kk]); A2 += E128D;
			}
		}

		REDUCE(s1)
		{
			s1[kk] = _mm_add_pd(s1[kk], s1[kk + KK]);
			s2[kk] = _mm_add_pd(s2[kk], s2[kk + KK]);
		}

		re1 = _mm_reduce_add_pd(s1[0]);
		re2 = _mm_reduce_add_pd(s2[0]);

		B -= E128D * sep;
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
	//suboptimal to compile
	{
		float re1 = 0, re2 = 0;
		UNROLLHEAD(4)
		for (int64 i = 0; i < n; ++i, A1++, A2++, B += sep)
		{
			volatile float v1 = (float)*A1 * (float)*B;
			volatile float v2 = (float)*A2 * (float)*B;
			re1 += v1;
			re2 += v2;
		}
		return re1 / re2;
	}
}

TARGETSSE double SumProdSSE(double* A, double* B, int64 sep, int64 n)
{
	//suboptimal to compile
	{
		double re = 0;
		UNROLLHEAD(4)
		for (int64 i = 0; i < n; ++i, A++, B += sep)
			re += (double)*A * (double)*B;
		return re;
	}

#define N 2
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		UNROLL(E128D) { _mm_prefetch((const char*)B, _MM_HINT_T0); B += sep; }

		__m128d s[E512_128] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N)
			{
				UNROLL(E512_128)
				{
					UNROLL(E128D) { _mm_prefetch((const char*)B, _MM_HINT_T0); B += sep; }
					
					s[kk] = _mm_fmaddx_pd(_mm_loadu_pd(A), _mm_set_pd(B[-4 * sep], B[-3 * sep]), s[kk]); A += E128D;
				}
			}
		}

		REDUCE(s) s[kk] = _mm_add_pd(s[kk], s[kk + KK]); 

		re = _mm_reduce_add_pd(s[0]);

		B -= E128D * sep;
	}

	for (; i < n; ++i, ++A, B += sep)
	{
		volatile double v1 = *A * *B;
		re += v1;
	}

	return re;
}

TARGETSSE double SumProdSSE(float* A, float* B, int64 sep, int64 n)
{
	//suboptimal to compile
	{
		double re = 0;
		UNROLLHEAD(4)
		for (int64 i = 0; i < n; ++i, A++, B += sep)
			re += (double)*A * (double)*B;
		return re;
	}

#define N 2
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		UNROLL(E128D) { _mm_prefetch((const char*)B, _MM_HINT_T0); B += sep; }

		__m128d s[E512_128] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128)
			{
				UNROLL(E128D) { _mm_prefetch((const char*)B, _MM_HINT_T0); B += sep; }
				
				s[kk] = _mm_fmaddx_pd(_mm_cvtps_pd(_mm_castsi128_ps(_mm_loadu_si64(A))), _mm_set_pd(B[-3 * sep], B[-4 * sep]), s[kk]); 
				
				A += E128D;
			}
		}

		REDUCE(s) s[kk] = _mm_add_pd(s[kk], s[kk + KK]);

		re = _mm_reduce_add_pd(s[0]);

		B -= E128D * sep;
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
	//suboptimal to compile
	{
		float re = 0;
		UNROLLHEAD(4)
		for (int64 i = 0; i < n; ++i, A++, B += sep)
		{
			volatile float v1 = *A * *B;
			re += v1;
		}
		return re;
	}
}

TARGETSSE double SumProdSSE(double* A, double* B, int64 n)
{
#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		__m128d s[E512_128] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128) 
			{ s[kk] = _mm_fmaddx_pd(_mm_loadu_pd(A), _mm_loadu_pd(B), s[kk]); A += E128D; B += E128D; }
		}

		REDUCE(s) s[kk] = _mm_add_pd(s[kk], s[kk + KK]);

		re = _mm_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, ++A, ++B)
	{
		volatile double v1 = (double)*A * (double)*B;
		re += v1;
	}

	return re;
}

TARGETSSE double SumProdSSE(float* A, float* B, int64 n)
{
#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		__m128d s[E512_128] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128) 
			{ s[kk] = _mm_fmaddx_pd(_mm_cvtps_pd(_mm_castsi128_ps(_mm_loadu_si64(A))), _mm_cvtps_pd(_mm_castsi128_ps(_mm_loadu_si64(B))), s[kk]); A += E128D; B += E128D; }
		}

		REDUCE(s) s[kk] = _mm_add_pd(s[kk], s[kk + KK]);

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
#define N 4
	int64 i = 0;
	volatile float re = 0;

	if (n >= N * E512F)
	{
		__m128 s[E512_128] = { 0 };

		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) UNROLL(E512_128) 
			{ s[kk] = _mm_fmaddx_ps(_mm_loadu_ps(A), _mm_loadu_ps(B), s[kk]); A += E128F; B += E128F; }
		}

		REDUCE(s) s[kk] = _mm_add_ps(s[kk], s[kk + KK]);

		re = _mm_reduce_add_ps(s[0]);
	}

	for (; i < n; ++i, ++A, ++B)
	{
		volatile float v1 = (float)*A * (float)*B;
		re += v1;
	}

	return re;
}

TARGETSSE double SumProdSSE(double* A, double* B, double* C, int64 n)
{
#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		__m128d s[E512_128] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128) 
			{ s[kk] = _mm_fmaddx_pd(_mm_mul_pd(_mm_loadu_pd(A), _mm_loadu_pd(B)), _mm_loadu_pd(C), s[kk]); A += E128D; B += E128D; C += E128D; }
		}

		REDUCE(s) s[kk] = _mm_add_pd(s[kk], s[kk + KK]);

		re = _mm_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, ++A, ++B, ++C)
	{
		volatile double v1 = (double)*A * (double)*B * (double)*C;
		re += v1;
	}

	return re;
}

TARGETSSE float SumProdSSE(float* A, float* B, float* C, int64 n)
{
#define N 4
	int64 i = 0;
	volatile float re = 0;

	if (n >= N * E512F)
	{
		__m128 s[E512_128] = { 0 };

		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) UNROLL(E512_128)
			{ s[kk] = _mm_fmaddx_ps(_mm_mul_ps(_mm_loadu_ps(A), _mm_loadu_ps(B)), _mm_loadu_ps(C), s[kk]); A += E128F; B += E128F; C += E128F; }
		}

		REDUCE(s) s[kk] = _mm_add_ps(s[kk], s[kk + KK]);

		re = _mm_reduce_add_ps(s[0]);
	}

	for (; i < n; ++i, ++A, ++B, ++C)
	{
		volatile float v1 = (float)*A * (float)*B * (float)*C;
		re += v1;
	}

	return re;
}

TARGETSSE double SumSqProdSSE(double* A, double* B, int64 n)
{
#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		__m128d s[E512_128] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128) 
			{ s[kk] = _mm_fmaddx_pd(_mm_mul_pd(_mm_loadu_pd(A), _mm_loadu_pd(A)), _mm_loadu_pd(B), s[kk]); A += E128D; B += E128D; }
		}

		REDUCE(s) s[kk] = _mm_add_pd(s[kk], s[kk + KK]);

		re = _mm_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, ++A, ++B)
	{
		volatile double v1 = (double)*A * (double)*A * (double)*B;
		re += v1;
	}

	return re;
}

TARGETSSE float SumSqProdSSE(float* A, float* B, int64 n)
{
#define N 4
	int64 i = 0;
	volatile float re = 0;

	if (n >= N * E512F)
	{
		__m128 s[E512_128] = { 0 };

		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) UNROLL(E512_128) 
			{ s[kk] = _mm_fmaddx_ps(_mm_mul_ps(_mm_loadu_ps(A), _mm_loadu_ps(A)), _mm_loadu_ps(B), s[kk]); A += E128F; B += E128F; }
		}

		REDUCE(s) s[kk] = _mm_add_ps(s[kk], s[kk + KK]);

		re = _mm_reduce_add_ps(s[0]);
	}

	for (; i < n; ++i, ++A, ++B)
	{
		volatile float v1 = (float)*A * (float)*A * (float)*B;
		re += v1;
	}

	return re;
}

TARGETSSE void AddSSE(double* A, double* B, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E512D)
	{
		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128) 
			{ _mm_storeu_pd(A, _mm_add_pd(_mm_loadu_pd(A), _mm_loadu_pd(B))); A += E128D; B += E128D; }
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGETSSE void AddSSE(float* A, float* B, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E512F)
	{
		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) UNROLL(E512_128) 
			{ _mm_storeu_ps(A, _mm_add_ps(_mm_loadu_ps(A), _mm_loadu_ps(B))); A += E128F; B += E128F; }
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGETSSE void AddSSE(int64* A, int64* B, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E128D)
	{
		for (int64 l1 = n - N * E128D; i <= l1; i += N * E128D)
		{
			UNROLL(N) { _mm_storeu_si128((__m128i*)A, _mm_add_epi64(_mm_loadu_si128((__m128i*)A), _mm_loadu_si128((__m128i*)B))); A += E128D; B += E128D; }
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGETSSE void AddSSE(int* A, int* B, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E128F)
	{
		for (int64 l1 = n - N * E128F; i <= l1; i += N * E128F)
		{
			UNROLL(N) { _mm_storeu_si128((__m128i*)A, _mm_add_epi32(_mm_loadu_si128((__m128i*)A), _mm_loadu_si128((__m128i*)B))); A += E128F; B += E128F; }
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGETSSE void AddSSE(int* A, int B, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E128F)
	{
		__m128i b = _mm_set1_epi32(B);

		for (int64 l1 = n - N * E128F; i <= l1; i += N * E128F)
		{
			UNROLL(N) { _mm_storeu_si128((__m128i*)A, _mm_add_epi32(_mm_loadu_si128((__m128i*)A), b)); A += E128F; }
		}
	}

	for (; i < n; ++i, A++)
		*A += B;
}

TARGETSSE void AddSSE(double* A, double B, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E128D)
	{
		__m128d b = _mm_set1_pd(B);

		for (int64 l1 = n - N * E128D; i <= l1; i += N * E128D)
		{
			UNROLL(N) {_mm_storeu_pd(A, _mm_add_pd(_mm_loadu_pd(A), b)); A += E128D; }
		}
	}

	for (; i < n; ++i, A++)
		*A += B;
}

TARGETSSE void AddSSE(float* A, float B, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E128F)
	{
		__m128 b = _mm_set1_ps(B);

		for (int64 l1 = n - N * E128F; i <= l1; i += N * E128F)
		{
			UNROLL(N) { _mm_storeu_ps(A, _mm_add_ps(_mm_loadu_ps(A), b)); A += E128F; }
		}
	}

	for (; i < n; ++i, A++)
		*A += B;
}

TARGETSSE void MulSSE(double* A, double* B, double* C, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E128D)
	{
		for (int64 l1 = n - N * E128D; i <= l1; i += N * E128D)
		{
			UNROLL(N) { _mm_storeu_pd(A, _mm_mul_pd(_mm_loadu_pd(B), _mm_loadu_pd(C))); A += E128D; B += E128D; C += E128D; }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *B++ * *C++;
		*A++ = v1;
	}
}

TARGETSSE void MulSSE(float* A, float* B, float* C, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E128F)
	{
		for (int64 l1 = n - N * E128F; i <= l1; i += N * E128F)
		{
			UNROLL(N) { _mm_storeu_ps(A, _mm_mul_ps(_mm_loadu_ps(B), _mm_loadu_ps(C))); A += E128F; B += E128F; C += E128F; }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *B++ * *C++;
		*A++ = v1;
	}
}

TARGETSSE void MulSSE(double* A, double* B, double C, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E128D)
	{
		__m128d c = _mm_set1_pd(C);

		for (int64 l1 = n - N * E128D; i <= l1; i += N * E128D)
		{
			UNROLL(N) { _mm_storeu_pd(A, _mm_mul_pd(_mm_loadu_pd(B), c)); A += E128D; B += E128D; }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *B++ * C;
		*A++ = v1;
	}
}

TARGETSSE void MulSSE(float* A, float* B, float C, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E128F)
	{
		__m128 c = _mm_set1_ps(C);

		for (int64 l1 = n - N * E128F; i <= l1; i += N * E128F)
		{
			UNROLL(N) { _mm_storeu_ps(A, _mm_mul_ps(_mm_loadu_ps(B), c)); A += E128F; B += E128F;}
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *B++ * C;
		*A++ = v1;
	}
}

TARGETSSE void MulSSE(double* A, double B, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E128D)
	{
		__m128d b = _mm_set1_pd(B);

		for (int64 l1 = n - N * E128D; i <= l1; i += N * E128D)
		{
			UNROLL(N) { _mm_storeu_pd(A, _mm_mul_pd(_mm_loadu_pd(A), b)); A += E128D; }
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
#define N 4
	int64 i = 0;

	if (n >= N * E128F)
	{
		__m128 b = _mm_set1_ps(B);

		for (int64 l1 = n - N * E128F; i <= l1; i += N * E128F)
		{
			UNROLL(N) { _mm_storeu_ps(A,_mm_mul_ps(_mm_loadu_ps(A), b)); A += E128F; }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A * B;
		*A++ = v1;
	}
}

TARGETSSE void DivSSE(double* A, double B, double* C, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E128D)
	{
		__m128d b = _mm_set1_pd(B);

		for (int64 l1 = n - N * E128D; i <= l1; i += N * E128D)
		{
			UNROLL(N) { _mm_storeu_pd(A, _mm_div_pd(b, _mm_loadu_pd(C))); A += E128D; C += E128D; }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = B / *C++;
		*A++ = v1;
	}
}

TARGETSSE void DivSSE(float* A, float B, float* C, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E128F)
	{
		__m128 b = _mm_set1_ps(B);

		for (int64 l1 = n - N * E128F; i <= l1; i += N * E128F)
		{
			UNROLL(N) { _mm_storeu_ps(A, _mm_div_ps(b, _mm_loadu_ps(C))); A += E128F; C += E128F; }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = B / *C++;
		*A++ = v1;
	}
}

TARGETSSE void DivSSE(double* A, double* B, double* C, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E128D)
	{
		for (int64 l1 = n - N * E128D; i <= l1; i += N * E128D)
		{
			UNROLL(N) { _mm_storeu_pd(A, _mm_div_pd(_mm_loadu_pd(B), _mm_loadu_pd(C))); A += E128D; B += E128D; C += E128D; }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *B++ / *C++;
		*A++ = v1;
	}
}

TARGETSSE void DivSSE(float* A, float* B, float* C, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E128F)
	{
		for (int64 l1 = n - N * E128F; i <= l1; i += N * E128F)
		{
			UNROLL(N) { _mm_storeu_ps(A, _mm_div_ps(_mm_loadu_ps(B), _mm_loadu_ps(C))); A += E128F; B += E128F; C += E128F; }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *B++ / *C++;
		*A++ = v1;
	}
}

TARGETSSE void AddProdSSE(double* A, double* B, double* C, int64 n)
{
#define N 1
	int64 i = 0;

	if (n >= N * E128D)
	{
		for (int64 l1 = n - N * E128D; i <= l1; i += N * E128D)
		{
			UNROLL(N) { _mm_storeu_pd(A, _mm_fmaddx_pd(_mm_loadu_pd(B), _mm_loadu_pd(C), _mm_loadu_pd(A))); A += E128D; B += E128D; C += E128D; }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *B++ * *C++;
		*A++ += v1;
	}
}

TARGETSSE void AddProdSSE(float* A, float* B, float* C, int64 n)
{
#define N 1
	int64 i = 0;

	if (n >= N * E128F)
	{
		for (int64 l1 = n - N * E128F; i <= l1; i += N * E128F)
		{
			UNROLL(N) { _mm_storeu_ps(A, _mm_fmaddx_ps(_mm_loadu_ps(B), _mm_loadu_ps(C), _mm_loadu_ps(A))); A += E128F; B += E128F; C += E128F; }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *B++ * *C++;
		*A++ += v1;
	}
}

TARGETSSE void AddProdSSE(double* A, double* B, double C, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E128D)
	{
		__m128d b[N], c = _mm_set1_pd(C);

		for (int64 l1 = n - N * E128D; i <= l1; i += N * E128D)
		{
			UNROLL(N) { b[kk] = _mm_loadu_pd(B); B += E128D; }

			UNROLL(N)  b[kk] = _mm_mul_pd(b[kk], c);

			double* A2 = A;

			UNROLL(N) { b[kk] = _mm_add_pd(b[kk], _mm_loadu_pd(A)); A += E128D; }

			UNROLL(N) { _mm_storeu_pd(A2, b[kk]); A2 += E128D; }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *B++ * C;
		*A++ += v1;
	}
}

TARGETSSE void AddProdSSE(double* A, float* B, double C, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E128D)
	{
		__m128d b[N], c = _mm_set1_pd(C);

		for (int64 l1 = n - N * E128D; i <= l1; i += N * E128D)
		{
			__m128 v1 = _mm_loadu_ps(B); B += E128F;
			b[0] = _mm_cvtps_pd(v1);
			b[1] = _mm_cvtps_pd(_mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(v1), 8)));

			__m128 v2 = _mm_loadu_ps(B); B += E128F;
			b[2] = _mm_cvtps_pd(v2);
			b[3] = _mm_cvtps_pd(_mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(v2), 8)));

			UNROLL(N) b[kk] = _mm_mul_pd(b[kk], c);

			double* A2 = A;

			UNROLL(N) { b[kk] = _mm_add_pd(b[kk], _mm_loadu_pd(A)); A += E128D; }

			UNROLL(N) { _mm_storeu_pd(A2, b[kk]); A2 += E128D; }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *B++ * C;
		*A++ += v1;
	}
}

TARGETSSE void AddProdSSE(float* A, float* B, float C, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E128F)
	{
		__m128 b[N], c = _mm_set1_ps(C);

		for (int64 l1 = n - N * E128F; i <= l1; i += N * E128F)
		{
			UNROLL(N) { b[kk] = _mm_loadu_ps(B); B += E128F; }

			UNROLL(N)  b[kk] = _mm_mul_ps(b[kk], c);

			float* A2 = A;

			UNROLL(N) { b[kk] = _mm_add_ps(b[kk], _mm_loadu_ps(A)); A += E128F; }

			UNROLL(N) { _mm_storeu_ps(A2, b[kk]); A2 += E128F; }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *B++ * C;
		*A++ += v1;
	}
}

TARGETSSE void UnifySSE(double* A, int64 n)
{
#define N 4
	int64 i = 0;
	double invsum = 1.0 / (SumSSE(A, n) + n * MIN_FREQ);
	
	if (n >= E128D)
	{
		__m128d b = _mm_set1_pd(invsum), c = _mm_set1_pd(MIN_FREQ * invsum);
		
		UNROLLHEAD(N)
		for (int64 l1 = n - E128D; i <= l1; i += E128D)
		{
			_mm_storeu_pd(A,_mm_fmaddx_pd(_mm_loadu_pd(A), b, c)); 
			A += E128D;
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
#define N 4
	int64 i = 0;
	double invsum = 1.0 / (SumSSE(A, n) + n * MIN_FREQ);

	if (n >= E128F)
	{
		__m128d b = _mm_set1_pd(invsum), c = _mm_set1_pd(MIN_FREQ * invsum);
		
		UNROLLHEAD(N)
		for (int64 l1 = n - E128F; i <= l1; i += E128F)
		{
			_mm_storeu_ps(A, _mm_shuffle_ps(
						_mm_cvtpd_ps(_mm_fmaddx_pd(_mm_set_pd(A[1], A[0]), b, c)), 
						_mm_cvtpd_ps(_mm_fmaddx_pd(_mm_set_pd(A[3], A[2]), b, c)), 
						_MM_SHUFFLE(1, 0, 1, 0)));
			A += E128F;
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
#define N 4
	A++; n--;
	int64 i = 0;

	if (n >= E512B)
	{
		__m128i r[N], v = _mm_set1_epi8(val);

		for (int64 l1 = n - E512B; i <= l1; i += E512B)
		{
			char* Ab = A;
			__m128i rz = _mm_setzero_si128();

			UNROLL(N) { rz = _mm_or_si128(rz, r[kk] = _mm_cmpeq_epi8(_mm_loadu_si128((__m128i*)A), v)); A += E128B; }

			if (_mm_testz_si128(rz, rz)) continue;

			uint64 mask; ushort* mask2 = (ushort*) & mask;
			UNROLL(N) mask2[kk] = (ushort)_mm_movemask_epi8(r[kk]);

			int count = _mm_popcnt_u64(mask);
			
			if (rep > count) [[likely]]
			{
				rep -= count;
				continue;
			}
			else for (;; Ab++, mask >>= 1)
				if ((mask & 1) && !--rep)
					return Ab;
		}
	}

	for (; i < n; ++i, A++)
		if (*A == val && !--rep)
			return A;

	return NULL;
}

TARGETSSE int64 CountCharSSE(char* A, char val, int64 n)
{
#define N 4
	uint64 re = 0;
	int64 i = 0;

	if (n >= E128B)
	{
		__m128i r = _mm_setzero_si128(), v = _mm_set1_epi8(val), o = _mm_setzero_si128();
		
		UNROLLHEAD(N)
		for (int64 l1 = n - E128B; i <= l1; i += E128B)
		{
			r = _mm_add_epi64(r, _mm_sad_epu8(o, _mm_cmpeq_epi8(_mm_loadu_si128((__m128i*)A), v)));
			A += E128B;
		}

		re = _mm_reduce_add_epi64(r) / 0xFF;
	}

	for (; i < n; ++i, A++)
		if (*A == val) re++;

	return (int64)re;
}

TARGETSSE void DiagQuadFormSSE(double* res1, double* A, double* D, int64 m, int64 n)
{
#define N 3

#define DECLARE			double* pA1 = (A + n * i), *pA2 = (A + n * j), *pD = D; __m128d a1[N], a2[N], r[N][N] = { 0 }
#define ALOAD1(ii)		a1[ii] = _mm_loadu_pd(pA1 + n * ii)
#define ALOAD2(ii)		a2[ii] = _mm_loadu_pd(pA2 + n * ii)
#define ALOAD3(ii)		a1[ii] = a2[ii] = _mm_loadu_pd(pA1 + n * ii)
#define ADMUL(ii)		a1[ii] = _mm_mul_pd(a1[ii], _mm_loadu_pd(pD))
#define FMADD(ii,jj)	r[ii][jj] = _mm_fmaddx_pd(a1[ii], a2[jj], r[ii][jj])
#define RDUADD(ii,jj)	res1[(i+ii) * m + (j+jj)] = _mm_reduce_add_pd(r[ii][jj])
#define REMAIN(ii,jj)	res1[(i+ii) * m + (j+jj)] += A[(i+ii) * n + k] * A[(j+jj) * n + k] * D[k]
#define FINAL(ii,jj)	res1[(i+ii) + (j+jj) * m] = res1[(i+ii) * m + (j+jj)]
	
    int64 i = 0, j = 0, k = 0;
    for (i = 0; i + N <= m; i += N)
	{
        for (j = 0; j < i; j += N)
		{
			DECLARE;

            for (k = 0; k + E128D <= n; k += E128D, pA1 += E128D, pA2 += E128D, pD += E128D) 
			{ LOOP(ALOAD1); LOOP(ALOAD2); LOOP(ADMUL); LOOP2(FMADD); }
			
			LOOP2(RDUADD);
			
			VECTORIZE
            for (; k < n; k++) 
			{ LOOP2(REMAIN); }

			LOOP2(FINAL);
        }
		
		for (; j <= i; j += N)
		{
			DECLARE;

            for (k = 0; k + E128D <= n; k += E128D, pA1 += E128D, pD += E128D) 
			{ LOOP(ALOAD3); LOOP(ADMUL); LOOP3(FMADD); }

			LOOP3(RDUADD);

			VECTORIZE
            for (; k < n; k++) 
			{ LOOP3(REMAIN); }

			LOOP3(FINAL);
		}
    }
	
    for (; i < m; i += N)
	{
		int64 Na = std::min(m - i, (int64)N);

        for (j = 0; j < i; j += N)
		{
			DECLARE;

            for (k = 0; k + E128D <= n; k += E128D, pA1 += E128D, pA2 += E128D, pD += E128D) 
			{ LOOPNa(ALOAD1); LOOP(ALOAD2); LOOPNa(ADMUL); LOOP2Na(FMADD); }
			
			LOOP2Na(RDUADD);
			
			VECTORIZE
            for (; k < n; k++) 
			{ LOOP2Na(REMAIN); }

			LOOP2Na(FINAL);
        }
		
		for (; j <= i; j += N)
		{
			DECLARE;

            for (k = 0; k + E128D <= n; k += E128D, pA1 += E128D, pD += E128D) 
			{ LOOPNa(ALOAD3); LOOP(ADMUL); LOOP3Na(FMADD); }

			LOOP3Na(RDUADD);

			VECTORIZE
            for (; k < n; k++) 
			{ LOOP3Na(REMAIN); }

			LOOP3Na(FINAL);
		}
    }

#undef DECLARE
#undef ALOAD1
#undef ALOAD2
#undef ALOAD3
#undef ADMUL
#undef FMADD
#undef RDUADD
#undef REMAIN
#undef FINAL
}

TARGETSSE void DiagQuadFormSSE(double* res2, double* A, double* D, double* B, int64 m, int64 n)
{
#define N 8

#define DECLARE			double* pA = (A + n * i), *pB = B, *pD = D; __m128d bd, r[N] = { 0 }
#define BDMUL			bd = _mm_mul_pd(_mm_loadu_pd(pB), _mm_loadu_pd(pD)); pB += E128D; pD += E128D
#define FMADD(ii)		r[ii] = _mm_fmaddx_pd(_mm_loadu_pd(pA + n * ii), bd, r[ii]); 
#define RDUADD(ii)		res2[(i+ii)] = _mm_reduce_add_pd(r[ii])
#define FINAL(ii)		res2[(i+ii)] += A[(i+ii) * n + k] * B[k] * D[k]

    int64 i = 0, k = 0;
    for (i = 0; i + N <= m; i += N)
	{
		DECLARE;

        for (k = 0; k + E128D <= n; k += E128D, pA += E128D) 
		{ BDMUL; LOOP(FMADD); }
			
		LOOP(RDUADD);

		VECTORIZE
        for (; k < n; k++) 
		{ LOOP(FINAL); }
	}
	
    for (; i < m; i += N)
	{
		int64 Na = std::min(m - i, (int64)N);
			
		DECLARE;

        for (k = 0; k + E128D <= n; k += E128D, pA += E128D) 
		{ BDMUL; LOOPNa(FMADD); }
			
		LOOPNa(RDUADD);
			
		VECTORIZE
        for (; k < n; k++) 
		{ LOOPNa(FINAL); }
	}

#undef DECLARE
#undef BDMUL
#undef FMADD
#undef RDUADD
#undef FINAL
}

TARGETSSE void DiagQuadFormSSE(double* res3, double* B, double* D, int64 n)
{
#define N 4
	__m128d s[4] = { 0 };

	int64 k = 0;
    for (; k + N * E128D <= n; k += N * E128D) 
	{
		UNROLL(N) 
		{
			s[kk] = _mm_fmaddx_pd(_mm_mul_pd(_mm_loadu_pd(B), _mm_loadu_pd(B)), _mm_loadu_pd(D), s[kk]); 
			B += E128D;
			D += E128D;
		}
	}

	REDUCE(s) s[kk] = _mm_add_pd(s[kk], s[kk + KK]);

	res3[0] = _mm_reduce_add_pd(s[0]);
			
	VECTORIZE
    for (; k < n; k++, B++, D++) 
		res3[0] += B[0] * B[0] * D[0];
}

TARGETSSE void MatrixMulSSE(double* res, double* A, double* B, int64 m, int64 n, int64 p)
{
#define N 3
	
#define DECLARE			double* pA = (A + n * i), *pB = (B + n * j); __m128d a[N], b[N], r[N][N] = { 0 }
#define ALOAD(kk)		a[kk] = _mm_loadu_pd(pA + n * kk)
#define BLOAD(kk)		b[kk] = _mm_loadu_pd(pB + n * kk)
#define FMADD(ii,jj)	r[ii][jj] = _mm_fmaddx_pd(a[ii], b[jj], r[ii][jj])
#define RDUADD(ii,jj)	res[(i+ii) + (j+jj) * m] = _mm_reduce_add_pd(r[ii][jj])
#define REMAIN(ii,jj)	res[(i+ii) + (j+jj) * m] += A[(i+ii) * n + k] * B[(j+jj) * n + k]

    int64 i = 0, j = 0, k = 0;
    for (i = 0; i + N <= m; i += N)
	{
        for (j = 0; j + N <= p; j += N)
		{
			DECLARE;

            for (k = 0; k + E128D <= n; k += E128D, pA += E128D, pB += E128D) 
			{ LOOP(ALOAD); LOOP(BLOAD); LOOP2(FMADD); }
			
			LOOP2(RDUADD);
			
			VECTORIZE
            for (; k < n; k++) 
			{ LOOP2(REMAIN); }
        }
		
		for (; j < p; j += N)
		{
			int64 Nb = std::min(p - j, (int64)N);
			
			DECLARE;

            for (k = 0; k + E128D <= n; k += E128D, pA += E128D, pB += E128D) 
			{ LOOP(ALOAD); LOOPNb(BLOAD); LOOP2Nb(FMADD); }
			
			LOOP2Nb(RDUADD);

			VECTORIZE
            for (; k < n; k++) 
			{ LOOP2Nb(REMAIN); }

		}
    }
	
    for (; i < m; i += N)
	{
		int64 Na = std::min(m - i, (int64)N);
        for (j = 0; j + N <= p; j += N)
		{
			DECLARE;

            for (k = 0; k + E128D <= n; k += E128D, pA += E128D, pB += E128D) 
			{ LOOPNa(ALOAD); LOOP(BLOAD); LOOP2Na(FMADD); }
			
			LOOP2Na(RDUADD);
			
			VECTORIZE
            for (; k < n; k++) 
			{ LOOP2Na(REMAIN); }
        }

		
		for (; j < p; j += N)
		{
			int64 Nb = std::min(p - j, (int64)N);
			
			DECLARE;

            for (k = 0; k + E128D <= n; k += E128D, pA += E128D, pB += E128D) 
			{ LOOPNa(ALOAD); LOOPNb(BLOAD); LOOP2NaNb(FMADD); }
			
			LOOP2NaNb(RDUADD);
			
			VECTORIZE
            for (; k < n; k++) 
			{ LOOP2NaNb(REMAIN); }

		}
    }

#undef DECLARE
#undef ALOAD
#undef BLOAD
#undef FMADD
#undef RDUADD
#undef FINAL
}

TARGETSSE void DiagQuadFormSSE(float* res1, float* A, float* D, int64 m, int64 n)
{
#define N 3

#define DECLARE			float* pA1 = (A + n * i), *pA2 = (A + n * j), *pD = D; __m128 a1[N], a2[N], r[N][N] = { 0 }
#define ALOAD1(ii)		a1[ii] = _mm_loadu_ps(pA1 + n * ii)
#define ALOAD2(ii)		a2[ii] = _mm_loadu_ps(pA2 + n * ii)
#define ALOAD3(ii)		a1[ii] = a2[ii] = _mm_loadu_ps(pA1 + n * ii)
#define ADMUL(ii)		a1[ii] = _mm_mul_ps(a1[ii], _mm_loadu_ps(pD))
#define FMADD(ii,jj)	r[ii][jj] = _mm_fmaddx_ps(a1[ii], a2[jj], r[ii][jj])
#define RDUADD(ii,jj)	res1[(i+ii) * m + (j+jj)] = _mm_reduce_add_ps(r[ii][jj])
#define REMAIN(ii,jj)	res1[(i+ii) * m + (j+jj)] += A[(i+ii) * n + k] * A[(j+jj) * n + k] * D[k]
#define FINAL(ii,jj)	res1[(i+ii) + (j+jj) * m] = res1[(i+ii) * m + (j+jj)]
	
    int64 i = 0, j = 0, k = 0;
    for (i = 0; i + N <= m; i += N)
	{
        for (j = 0; j < i; j += N)
		{
			DECLARE;

            for (k = 0; k + E128F <= n; k += E128F, pA1 += E128F, pA2 += E128F, pD += E128F) 
			{ LOOP(ALOAD1); LOOP(ALOAD2); LOOP(ADMUL); LOOP2(FMADD); }
			
			LOOP2(RDUADD);
			
			VECTORIZE
            for (; k < n; k++) 
			{ LOOP2(REMAIN); }

			LOOP2(FINAL);
        }
		
		for (; j <= i; j += N)
		{
			DECLARE;

            for (k = 0; k + E128F <= n; k += E128F, pA1 += E128F, pD += E128F) 
			{ LOOP(ALOAD3); LOOP(ADMUL); LOOP3(FMADD); }

			LOOP3(RDUADD);

			VECTORIZE
            for (; k < n; k++) 
			{ LOOP3(REMAIN); }

			LOOP3(FINAL);
		}
    }
	
    for (; i < m; i += N)
	{
		int64 Na = std::min(m - i, (int64)N);

        for (j = 0; j < i; j += N)
		{
			DECLARE;

            for (k = 0; k + E128F <= n; k += E128F, pA1 += E128F, pA2 += E128F, pD += E128F) 
			{ LOOPNa(ALOAD1); LOOP(ALOAD2); LOOPNa(ADMUL); LOOP2Na(FMADD); }
			
			LOOP2Na(RDUADD);
			
			VECTORIZE
            for (; k < n; k++) 
			{ LOOP2Na(REMAIN); }

			LOOP2Na(FINAL);
        }
		
		for (; j <= i; j += N)
		{
			DECLARE;

            for (k = 0; k + E128F <= n; k += E128F, pA1 += E128F, pD += E128F) 
			{ LOOPNa(ALOAD3); LOOP(ADMUL); LOOP3Na(FMADD); }

			LOOP3Na(RDUADD);

			VECTORIZE
            for (; k < n; k++) 
			{ LOOP3Na(REMAIN); }

			LOOP3Na(FINAL);
		}
    }

#undef DECLARE
#undef ALOAD1
#undef ALOAD2
#undef ALOAD3
#undef ADMUL
#undef FMADD
#undef RDUADD
#undef REMAIN
#undef FINAL
}

TARGETSSE void DiagQuadFormSSE(float* res2, float* A, float* D, float* B, int64 m, int64 n)
{
#define N 8

#define DECLARE			float* pA = (A + n * i), *pB = B, *pD = D; __m128 bd, r[N] = { 0 }
#define BDMUL			bd = _mm_mul_ps(_mm_loadu_ps(pB), _mm_loadu_ps(pD)); pB += E128F; pD += E128F
#define FMADD(ii)		r[ii] = _mm_fmaddx_ps(_mm_loadu_ps(pA + n * ii), bd, r[ii]); 
#define RDUADD(ii)		res2[(i+ii)] = _mm_reduce_add_ps(r[ii])
#define FINAL(ii)		res2[(i+ii)] += A[(i+ii) * n + k] * B[k] * D[k]

    int64 i = 0, k = 0;
    for (i = 0; i + N <= m; i += N)
	{
		DECLARE;

        for (k = 0; k + E128F <= n; k += E128F, pA += E128F) 
		{ BDMUL; LOOP(FMADD); }
			
		LOOP(RDUADD);

		VECTORIZE
        for (; k < n; k++) 
		{ LOOP(FINAL); }
	}
	
    for (; i < m; i += N)
	{
		int64 Na = std::min(m - i, (int64)N);
			
		DECLARE;

        for (k = 0; k + E128F <= n; k += E128F, pA += E128F) 
		{ BDMUL; LOOPNa(FMADD); }
			
		LOOPNa(RDUADD);
			
		VECTORIZE
        for (; k < n; k++) 
		{ LOOPNa(FINAL); }
	}

#undef DECLARE
#undef BDMUL
#undef FMADD
#undef RDUADD
#undef FINAL
}

TARGETSSE void DiagQuadFormSSE(float* res3, float* B, float* D, int64 n)
{
#define N 4
	__m128 s[4] = { 0 };

	int64 k = 0;
    for (; k + N * E128F <= n; k += N * E128F) 
	{
		UNROLL(N) 
		{
			s[kk] = _mm_fmaddx_ps(_mm_mul_ps(_mm_loadu_ps(B), _mm_loadu_ps(B)), _mm_loadu_ps(D), s[kk]); 
			B += E128F;
			D += E128F;
		}
	}

	REDUCE(s) s[kk] = _mm_add_ps(s[kk], s[kk + KK]);

	res3[0] = _mm_reduce_add_ps(s[0]);
			
	VECTORIZE
    for (; k < n; k++, B++, D++) 
		res3[0] += B[0] * B[0] * D[0];
}

TARGETSSE void MatrixMulSSE(float* res, float* A, float* B, int64 m, int64 n, int64 p)
{
#define N 3
	
#define DECLARE			float* pA = (A + n * i), *pB = (B + n * j); __m128 a[N], b[N], r[N][N] = { 0 }
#define ALOAD(kk)		a[kk] = _mm_loadu_ps(pA + n * kk)
#define BLOAD(kk)		b[kk] = _mm_loadu_ps(pB + n * kk)
#define FMADD(ii,jj)	r[ii][jj] = _mm_fmaddx_ps(a[ii], b[jj], r[ii][jj])
#define RDUADD(ii,jj)	res[(i+ii) + (j+jj) * m] = _mm_reduce_add_ps(r[ii][jj])
#define REMAIN(ii,jj)	res[(i+ii) + (j+jj) * m] += A[(i+ii) * n + k] * B[(j+jj) * n + k]

    int64 i = 0, j = 0, k = 0;
    for (i = 0; i + N <= m; i += N)
	{
        for (j = 0; j + N <= p; j += N)
		{
			DECLARE;

            for (k = 0; k + E128F <= n; k += E128F, pA += E128F, pB += E128F) 
			{ LOOP(ALOAD); LOOP(BLOAD); LOOP2(FMADD); }
			
			LOOP2(RDUADD);
			
			VECTORIZE
            for (; k < n; k++) 
			{ LOOP2(REMAIN); }
        }
		
		for (; j < p; j += N)
		{
			int64 Nb = std::min(p - j, (int64)N);
			
			DECLARE;

            for (k = 0; k + E128F <= n; k += E128F, pA += E128F, pB += E128F) 
			{ LOOP(ALOAD); LOOPNb(BLOAD); LOOP2Nb(FMADD); }
			
			LOOP2Nb(RDUADD);

			VECTORIZE
            for (; k < n; k++) 
			{ LOOP2Nb(REMAIN); }

		}
    }
	
    for (; i < m; i += N)
	{
		int64 Na = std::min(m - i, (int64)N);
        for (j = 0; j + N <= p; j += N)
		{
			DECLARE;

            for (k = 0; k + E128F <= n; k += E128F, pA += E128F, pB += E128F) 
			{ LOOPNa(ALOAD); LOOP(BLOAD); LOOP2Na(FMADD); }
			
			LOOP2Na(RDUADD);
			
			VECTORIZE
            for (; k < n; k++) 
			{ LOOP2Na(REMAIN); }
        }

		
		for (; j < p; j += N)
		{
			int64 Nb = std::min(p - j, (int64)N);
			
			DECLARE;

            for (k = 0; k + E128F <= n; k += E128F, pA += E128F, pB += E128F) 
			{ LOOPNa(ALOAD); LOOPNb(BLOAD); LOOP2NaNb(FMADD); }
			
			LOOP2NaNb(RDUADD);
			
			VECTORIZE
            for (; k < n; k++) 
			{ LOOP2NaNb(REMAIN); }

		}
    }

#undef DECLARE
#undef ALOAD
#undef BLOAD
#undef FMADD
#undef RDUADD
#undef FINAL
}

#endif
