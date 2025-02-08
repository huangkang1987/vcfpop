/* AVX Instruction Set Functions */

#include "vcfpop.h"

#ifndef __aarch64__

template struct RNGAVX<double>;
template struct RNGAVX<float >;
template TARGETAVX void RNGAVX<double>::Integer<uint  >(uint  * re, int64 n, uint   minv, uint   maxv);
template TARGETAVX void RNGAVX<double>::Integer<uint64>(uint64* re, int64 n, uint64 minv, uint64 maxv);
template TARGETAVX void RNGAVX<float >::Integer<uint  >(uint  * re, int64 n, uint   minv, uint   maxv);
template TARGETAVX void RNGAVX<float >::Integer<uint64>(uint64* re, int64 n, uint64 minv, uint64 maxv);

#ifndef _RNGAVX_FP64
#define XS 16

/* Initialize rng */
TARGETAVX RNGAVX<double>::RNGAVX()
{

}

/* Initialize rng */
TARGETAVX RNGAVX<double>::RNGAVX(uint64 seed, uint64 salt)
{
	__m256i a[XS], s, m;

	UNROLL(XS) { a[kk] = _mm256_set_epi64x(seed + 3, seed + 2, seed + 1, seed + 0); seed += 4; }

	s    = _mm256_set1_epi64x(salt);
	m    = _mm256_set1_epi32(0x5bd1e995);

	UNROLL(XS) a[kk] = _mm256_xor_si256(a[kk], _mm256_slli_epi64(_mm256_andnot_si256(a[kk], _mm256_set1_epi64x(0xFFFFFFFF)), 32));

	s = _mm256_xor_si256(s, _mm256_slli_epi64(_mm256_andnot_si256(s, _mm256_set1_epi64x(0xFFFFFFFF)), 32));

	// uint s = s ^ 4;
	s = _mm256_xor_si256(s, _mm256_set1_epi32(4));

	// a *= m;
	UNROLL(XS) a[kk] = _mm256_mullo_epi32(a[kk], m);//OK

	// a ^= a >> 24;
	UNROLL(XS) a[kk] = _mm256_xor_si256(a[kk], _mm256_srli_epi32(a[kk], 24));

	// a *= m;
	UNROLL(XS) a[kk] = _mm256_mullo_epi32(a[kk], m);//OK

	// s *= m;
	s = _mm256_mullo_epi32(s, m);//OK

	// a ^= s;
	UNROLL(XS) a[kk] = _mm256_xor_si256(a[kk], s);

	// a ^= a >> 13;
	UNROLL(XS) a[kk] = _mm256_xor_si256(a[kk], _mm256_srli_epi32(a[kk], 13));
	
	// a *= m;
	UNROLL(XS) a[kk] = _mm256_mullo_epi32(a[kk], m);//OK

	// a ^= a >> 15;
	UNROLL(XS) a[kk] = _mm256_xor_si256(a[kk], _mm256_srli_epi32(a[kk], 15));

	// original
	UNROLL(XS) x[kk] = _mm256_xor_si256(_mm256_set1_epi64x(0x159A55E5075BCD15), a[kk]);

	UNROLL(XS) a[kk] = _mm256_slli_epi64(a[kk], 6);

	UNROLL(XS) y[kk] = _mm256_xor_si256(_mm256_set1_epi64x(0x054913331F123BB5), a[kk]);
}

/* Draw 64 64-bit integers in [0,n), 64*n frequencies are in arr */
TARGETAVX void RNGAVX<double>::Poly(__m256d* arr, int n, __m256i* re)
{
	__m256d t[XS], s[XS];
	__m256d one = _mm256_set1_pd(1.0);
	__m256i mask1 = _mm256_set1_epi64x(0x000FFFFFFFFFFFFF);
	__m256i mask2 = _mm256_set1_epi64x(0x3FF0000000000000);
	__m256i* r = (__m256i*)t;

	UNROLL(XS) s[kk] = _mm256_setzero_pd();

	for (int i = 0; i < n * XS; i += XS)
		UNROLL(XS) s[kk] = _mm256_add_pd(s[kk], arr[kk + i]);

	XorShift();

	UNROLL(XS) r[kk] = _mm256_add_epi64(x[kk], y[kk]);

	UNROLL(XS) r[kk] = _mm256_and_si256(r[kk], mask1);

	UNROLL(XS) r[kk] = _mm256_or_si256(r[kk], mask2);

	UNROLL(XS) t[kk] = _mm256_sub_pd(t[kk], one);

	UNROLL(XS) t[kk] = _mm256_mul_pd(t[kk], s[kk]);

	__m256i midx[XS], nidx = _mm256_setzero_si256(), ninc = _mm256_set1_epi64x(1);
	__m256d f[XS], b[XS];

	UNROLL(XS) midx[kk] = _mm256_set1_epi64x(n - 1);

	UNROLL(XS) f[kk] = _mm256_setzero_pd();

	for (int i = 0; i < n * XS; i += XS)
	{
		UNROLL(XS) b[kk] = _mm256_cmp_pd(t[kk], arr[kk + i], _CMP_LT_OS);

		UNROLL(XS) t[kk] = _mm256_sub_pd(t[kk], arr[kk + i]);

		UNROLL(XS) b[kk] = _mm256_andnot_pd(f[kk], b[kk]);

		UNROLL(XS) f[kk] = _mm256_or_pd(f[kk], b[kk]);

		UNROLL(XS) midx[kk] = _mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(midx[kk]), _mm256_castsi256_pd(nidx), b[kk]));//ok

		nidx = _mm256_add_epi64(nidx, ninc);
	}

	UNROLL(XS) re[kk] = midx[kk];
}

/* Draw uniform distriubted intergers */
TARGET void RNGAVX<double>::XorShift()
{
	__m256i a[XS], b[XS];

	UNROLL(XS) a[kk] = x[kk];

	UNROLL(XS) b[kk] = y[kk];

	UNROLL(XS) x[kk] = b[kk];

	UNROLL(XS) a[kk] = _mm256_xor_si256(a[kk], _mm256_slli_epi64(a[kk], 23));

	UNROLL(XS) a[kk] = _mm256_xor_si256(a[kk], _mm256_srli_epi64(a[kk], 18));

	UNROLL(XS) a[kk] = _mm256_xor_si256(a[kk], b[kk]);

	UNROLL(XS) a[kk] = _mm256_xor_si256(a[kk], _mm256_srli_epi64(b[kk], 5));

	UNROLL(XS) y[kk] = a[kk];
}

/* Draw uniform distriubted integers */
template<typename INT>
TARGET void RNGAVX<double>::Integer(INT* re, int64 n, INT minv, INT maxv)
{
	constexpr int xesize = sizeof(x) / sizeof(INT);
	int64 i = 0;
	INT modv = maxv - minv;
	INT* rei = re;

	for (; i <= n - xesize; i += xesize)
	{
		XorShift();
		UNROLL(XS) { _mm256_storeu_si256((__m256i*)rei, _mm256_add_epi64(x[kk], y[kk])); rei += E256B / sizeof(INT); };
	}

	if (i != n)
	{
		__m256i re2[XS];
		XorShift();
		UNROLL(XS) re2[kk] = _mm256_add_epi64(x[kk], y[kk]);
		SetVal((INT*)rei, (INT*)re2, n - i);
	}

	if (maxv != (INT)-1 || minv != 0)
	{
		for (i = 0; i < n; ++i)
			re[i] = re[i] % modv + minv;
	}
}

/* Draw uniform distriubted real numbers */
TARGET void RNGAVX<double>::Uniform(double* re, int n, double minv, double maxv)
{
	constexpr int xesize = sizeof(x) / sizeof(double);
	int i = 0;
	double range = maxv - minv;

	__m256i mask1 = _mm256_set1_epi64x(0x000FFFFFFFFFFFFF);
	__m256i mask2 = _mm256_set1_epi64x(0x3FF0000000000000);
	__m256d v1 = _mm256_set1_pd(minv - range);
	__m256d v2 = _mm256_set1_pd(range);

	if (range == 1.0)
	{
		for (; i <= n - xesize; i += xesize)
		{
			XorShift();
			UNROLL(XS) { _mm256_storeu_pd(re, _mm256_add_pd(_mm256_castsi256_pd(_mm256_or_si256(_mm256_and_si256(_mm256_add_epi64(x[kk], y[kk]), mask1), mask2)), v1)); re += E256D; }
		}
	}
	else
	{
		for (; i <= n - xesize; i += xesize)
		{
			XorShift();
			UNROLL(XS) { _mm256_storeu_pd(re, _mm256_fmaddx_pd(_mm256_castsi256_pd(_mm256_or_si256(_mm256_and_si256(_mm256_add_epi64(x[kk], y[kk]), mask1), mask2)), v2, v1)); re += E256D; }
		}
	}

	if (i != n)
	{
		__m256d ref[XS];
		XorShift();
		UNROLL(XS) ref[kk] = _mm256_fmaddx_pd(_mm256_castsi256_pd(_mm256_or_si256(_mm256_and_si256(_mm256_add_epi64(x[kk], y[kk]), mask1), mask2)), v2, v1);
		SetVal((double*)re, (double*)ref, n - i);
	}
}

/* Draw uniform distriubted real numbers */
TARGET void RNGAVX<double>::Normal(double* re, int n, double mean, double sd)
{
	constexpr int xhsize = XS / 2;
	constexpr int xesize = XS * E256D;

	int i = 0;

	__m256i mask1 = _mm256_set1_epi64x(0x000FFFFFFFFFFFFF);
	__m256i mask2 = _mm256_set1_epi64x(0x3FF0000000000000);
	__m256d v1 = _mm256_set1_pd(-1);
	__m256d min_freq = _mm256_set1_pd(MIN_FREQ);
	__m256d pi2 = _mm256_set1_pd(2.0 * M_PI);
	__m256d mu = _mm256_set1_pd(mean);
	__m256d s = _mm256_set1_pd(sd);

	for (; i <= n - xesize; i += xesize)
	{
		XorShift();
		UNROLL(XS) _mm256_storeu_pd(re + E256D * kk, _mm256_add_pd(_mm256_castsi256_pd(_mm256_or_si256(_mm256_and_si256(_mm256_add_epi64(x[kk], y[kk]), mask1), mask2)), v1));

		__m256d u1, u2, u3, u4;
		for (int j = 0; j < xhsize; ++j)
		{
			u1 = _mm256_max_pd(_mm256_loadu_pd(re + j * E256D), min_freq);
			u2 = _mm256_mul_pd(_mm256_loadu_pd(re + j * E256D + xhsize * E256D), pi2);

			UNROLL(E256D) simd_f64(u1, kk) = sqrt(-2.0 * log(simd_f64(u1, kk)));
			UNROLL(E256D) simd_f64(u3, kk) = cos(simd_f64(u2, kk));
			UNROLL(E256D) simd_f64(u4, kk) = sin(simd_f64(u2, kk));

			_mm256_storeu_pd(re + j * E256D, _mm256_mul_pd(u1, u3));
			_mm256_storeu_pd(re + j * E256D + xhsize * E256D, _mm256_mul_pd(u1, u4));
		}
		
		if (sd != 1 || mean != 0)
			UNROLL(XS) { _mm256_storeu_pd(re, _mm256_fmaddx_pd(_mm256_loadu_pd(re), s, mu)); re += E256D; }
		else
			re += XS * E256D;
	}

	if (i != n)
	{
		double ref[XS * E256D];

		XorShift();
		UNROLL(XS) _mm256_storeu_pd(ref + E256D * kk, _mm256_add_pd(_mm256_castsi256_pd(_mm256_or_si256(_mm256_and_si256(_mm256_add_epi64(x[kk], y[kk]), mask1), mask2)), v1));

		__m256d u1, u2, u3, u4;
		for (int j = 0; j < xhsize; ++j)
		{
			u1 = _mm256_max_pd(_mm256_loadu_pd(ref + j * E256D), min_freq);
			u2 = _mm256_mul_pd(_mm256_loadu_pd(ref + j * E256D + xhsize * E256D), pi2);

			UNROLL(E256D) simd_f64(u1, kk) = sqrt(-2.0 * log(simd_f64(u1, kk)));
			UNROLL(E256D) simd_f64(u3, kk) = cos(simd_f64(u2, kk));
			UNROLL(E256D) simd_f64(u4, kk) = sin(simd_f64(u2, kk));

			_mm256_storeu_pd(ref + j * E256D, _mm256_mul_pd(u1, u3));
			_mm256_storeu_pd(ref + j * E256D + xhsize * E256D, _mm256_mul_pd(u1, u4));
		}
		
		if (sd != 1 || mean != 0)
			UNROLL(XS) _mm256_storeu_pd(ref + kk * E256D, _mm256_fmaddx_pd(_mm256_loadu_pd(re), s, mu)); 

		SetVal((double*)re, (double*)ref, n - i);
	}
}
#endif

#ifndef _RNGAVX_FP32
#define XS 8
#define XS2 16
/* Initialize rng */
TARGETAVX RNGAVX<float>::RNGAVX()
{
}

/* Initialize rng */
TARGETAVX RNGAVX<float>::RNGAVX(uint64 seed, uint64 salt)
{
	__m256i a[XS], s, m;
	UNROLL(XS) { a[kk] = _mm256_set_epi32(Mix(seed + 7), Mix(seed + 6), Mix(seed + 5), Mix(seed + 4), Mix(seed + 3), Mix(seed + 2), Mix(seed + 1), Mix(seed + 0)); seed += 8; }
		
	s = _mm256_set1_epi32(Mix(salt));
	m = _mm256_set1_epi32(0x5bd1e995);

	// uint s = s ^ 4;
	s = _mm256_xor_si256(s, _mm256_set1_epi32(4));

	// a *= m;
	UNROLL(XS) a[kk] = _mm256_mullo_epi32(a[kk], m);//OK

	// a ^= a >> 24;
	UNROLL(XS) a[kk] = _mm256_xor_si256(a[kk], _mm256_srli_epi32(a[kk], 24));

	// a *= m;
	UNROLL(XS) a[kk] = _mm256_mullo_epi32(a[kk], m);//OK

	// s *= m;
	s = _mm256_mullo_epi32(s, m);//OK

	// a ^= s;
	UNROLL(XS) a[kk] = _mm256_xor_si256(a[kk], s);

	// a ^= a >> 13;
	UNROLL(XS) a[kk] = _mm256_xor_si256(a[kk], _mm256_srli_epi32(a[kk], 13));

	// a *= m;
	UNROLL(XS) a[kk] = _mm256_mullo_epi32(a[kk], m);//OK

	// a ^= a >> 15;
	UNROLL(XS) a[kk] = _mm256_xor_si256(a[kk], _mm256_srli_epi32(a[kk], 15));

	// original
	UNROLL(XS) x[kk] = _mm256_xor_si256(_mm256_set1_epi32(0x075BCD15), a[kk]);

	UNROLL(XS) a[kk] = _mm256_slli_epi32(a[kk], 3);

	UNROLL(XS) y[kk] = _mm256_xor_si256(_mm256_set1_epi32(0x159A55E5), a[kk]);

	UNROLL(XS) a[kk] = _mm256_slli_epi32(a[kk], 3);

	UNROLL(XS) z[kk] = _mm256_xor_si256(_mm256_set1_epi32(0x1F123BB5), a[kk]);
}

/* Draw 64 64-bit integers in [0,n), 64*n frequencies are in arr */
TARGETAVX void RNGAVX<float>::Poly(__m256* arr, int n, __m256i* re)
{
	__m256 t[XS], s[XS];
	__m256 one = _mm256_set1_ps(1.0f);
	__m256i mask1 = _mm256_set1_epi32(0x007FFFFF);
	__m256i mask2 = _mm256_set1_epi32(0x3F800000);
	__m256i* r = (__m256i*)t;

	UNROLL(XS) s[kk] = _mm256_setzero_ps();

	for (int i = 0; i < n * XS; i += XS)
		UNROLL(XS) s[kk] = _mm256_add_ps(s[kk], arr[kk + i]);

	XorShift();

	UNROLL(XS) r[kk] = _mm256_and_si256(z[kk], mask1);

	UNROLL(XS) r[kk] = _mm256_or_si256(r[kk], mask2);

	UNROLL(XS) t[kk] = _mm256_sub_ps(t[kk], one);

	UNROLL(XS) t[kk] = _mm256_mul_ps(t[kk], s[kk]);

	__m256i midx[XS], nidx = _mm256_setzero_si256(), ninc = _mm256_set1_epi32(1);
	__m256 f[XS], b[XS];
	UNROLL(XS) midx[kk] = _mm256_set1_epi32(n - 1);
	UNROLL(XS) f[kk] = _mm256_setzero_ps();

	for (int i = 0; i < n * XS; i += XS)
	{
		UNROLL(XS) b[kk] = _mm256_cmp_ps(t[kk], arr[kk + i], _CMP_LT_OS);

		UNROLL(XS) t[kk] = _mm256_sub_ps(t[kk], arr[kk + i]);

		UNROLL(XS) b[kk] = _mm256_andnot_ps(f[kk], b[kk]);

		UNROLL(XS) f[kk] = _mm256_or_ps(f[kk], b[kk]);

		UNROLL(XS) midx[kk] = _mm256_castps_si256(_mm256_blendv_ps(_mm256_castsi256_ps(midx[kk]), _mm256_castsi256_ps(nidx), b[kk]));//ok

		nidx = _mm256_add_epi32(nidx, ninc);
	}

	__m128i* midx2 = (__m128i*)midx;
	UNROLL(XS2) re[kk] = _mm256_cvtepi32_epi64(midx2[kk]);
}

/* Draw uniform distriubted intergers */
TARGET void RNGAVX<float>::XorShift()
{
	__m256i u[XS];

	UNROLL(XS) u[kk] = _mm256_slli_epi32(x[kk], 16);
	UNROLL(XS) x[kk] = _mm256_xor_si256(x[kk], u[kk]);

	UNROLL(XS) u[kk] = _mm256_srli_epi32(x[kk], 5);
	UNROLL(XS) x[kk] = _mm256_xor_si256(x[kk], u[kk]);

	UNROLL(XS) u[kk] = _mm256_slli_epi32(x[kk], 1);
	UNROLL(XS) x[kk] = _mm256_xor_si256(x[kk], u[kk]);

	UNROLL(XS) u[kk] = x[kk];

	UNROLL(XS) x[kk] = y[kk];

	UNROLL(XS) y[kk] = z[kk];

	UNROLL(XS) z[kk] = _mm256_xor_si256(u[kk], x[kk]);

	UNROLL(XS) z[kk] = _mm256_xor_si256(z[kk], y[kk]);
}

/* Draw uniform distriubted integers */
template<typename INT>
TARGET void RNGAVX<float>::Integer(INT* re, int64 n, INT minv, INT maxv)
{
	constexpr int xesize = sizeof(x) / sizeof(INT);
	int64 i = 0;
	INT modv = maxv - minv;
	INT* rei = re;

	for (; i <= n - xesize; i += xesize)
	{
		XorShift();
		UNROLL(XS) { _mm256_storeu_si256((__m256i*)rei, z[kk]); rei += E256B / sizeof(INT); }
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
TARGET void RNGAVX<float>::Uniform(float* re, int n, float minv, float maxv)
{
	constexpr int xesize = sizeof(x) / sizeof(float);

	int i = 0;
	float range = maxv - minv;

	__m256i mask1 = _mm256_set1_epi32(0x007FFFFF);
	__m256i mask2 = _mm256_set1_epi32(0x3F800000);
	__m256 v1 = _mm256_set1_ps(minv - range);
	__m256 v2 = _mm256_set1_ps(range);

	if (range == 1.0)
	{
		for (; i <= n - xesize; i += xesize)
		{
			XorShift();
			UNROLL(XS) { _mm256_storeu_ps(re, _mm256_add_ps(_mm256_castsi256_ps(_mm256_or_si256(_mm256_and_si256(z[kk], mask1), mask2)), v1)); re += E256F; }
		}
	}
	else
	{
		for (; i <= n - xesize; i += xesize)
		{
			XorShift();
			UNROLL(XS) { _mm256_storeu_ps(re, _mm256_fmaddx_ps(_mm256_castsi256_ps(_mm256_or_si256(_mm256_and_si256(z[kk], mask1), mask2)), v2, v1)); re += E256F; }
		}
	}

	if (i != n)
	{
		__m256 ref[XS];
		XorShift();
		UNROLL(XS) ref[kk] = _mm256_fmaddx_ps(_mm256_castsi256_ps(_mm256_or_si256(_mm256_and_si256(z[kk], mask1), mask2)), v2, v1);
		SetVal((float*)re, (float*)ref, n - i);
	}
}

/* Draw uniform distriubted real numbers */
TARGET void RNGAVX<float>::Normal(float* re, int n, float mean, float sd)
{
	constexpr int xhsize = XS / 2;
	constexpr int xesize = XS * E256F;

	int i = 0;

	__m256i mask1 = _mm256_set1_epi32(0x007FFFFF);
	__m256i mask2 = _mm256_set1_epi32(0x3F800000);
	__m256 v1 = _mm256_set1_ps(- 1);
	__m256 min_freq = _mm256_set1_ps((float)MIN_FREQ);
	__m256 pi2 = _mm256_set1_ps((float)(2.0 * M_PI));
	__m256 mu = _mm256_set1_ps(mean);
	__m256 s  = _mm256_set1_ps(sd);

	for (; i <= n - xesize; i += xesize)
	{
		XorShift();
		UNROLL(XS) _mm256_storeu_ps(re + kk * E256F, _mm256_add_ps(_mm256_castsi256_ps(_mm256_or_si256(_mm256_and_si256(z[kk], mask1), mask2)), v1));

		__m256 u1, u2, u3, u4;
		for (int j = 0; j < xhsize; ++j)
		{
			u1 = _mm256_max_ps(_mm256_loadu_ps(re + j * E256F), min_freq);
			u2 = _mm256_mul_ps(_mm256_loadu_ps(re + j * E256F + xhsize * E256F), pi2);

			UNROLL(E256F) simd_f32(u1, kk) = sqrt(-2.0 * log(simd_f32(u1, kk)));
			UNROLL(E256F) simd_f32(u3, kk) = cos(simd_f32(u2, kk));
			UNROLL(E256F) simd_f32(u4, kk) = sin(simd_f32(u2, kk));
				
			_mm256_storeu_ps(re + j * E256F, _mm256_mul_ps(u1, u3));
			_mm256_storeu_ps(re + j * E256F + xhsize * E256F, _mm256_mul_ps(u1, u4));
		}
			
		if (sd != 1 || mean != 0)
			UNROLL(XS) { _mm256_storeu_ps(re, _mm256_fmaddx_ps(_mm256_loadu_ps(re), s, mu)); re += E256F; }
		else
			re += XS * E256F;
	}

	if (i != n)
	{
		float ref[XS * E256F];

		XorShift();
		UNROLL(XS) _mm256_storeu_ps(ref + kk * E256F, _mm256_add_ps(_mm256_castsi256_ps(_mm256_or_si256(_mm256_and_si256(z[kk], mask1), mask2)), v1));
		
		__m256 u1, u2, u3, u4;
		for (int j = 0; j < xhsize; ++j)
		{
			u1 = _mm256_max_ps(_mm256_loadu_ps(ref + j * E256F), min_freq);
			u2 = _mm256_mul_ps(_mm256_loadu_ps(ref + j * E256F + xhsize * E256F), pi2);

			UNROLL(E256F) simd_f32(u1, kk) = sqrt(-2.0 * log(simd_f32(u1, kk)));
			UNROLL(E256F) simd_f32(u3, kk) = cos(simd_f32(u2, kk));
			UNROLL(E256F) simd_f32(u4, kk) = sin(simd_f32(u2, kk));

			_mm256_storeu_ps(&ref[j * E256F], _mm256_mul_ps(u1, u3));
			_mm256_storeu_ps(&ref[j * E256F + xhsize * E256F], _mm256_mul_ps(u1, u4));
		}
			
		if (sd != 1 || mean != 0)
			UNROLL(XS) _mm256_storeu_ps(ref + kk * E256F, _mm256_fmaddx_ps(_mm256_loadu_ps(re), s, mu)); 

		SetVal((float*)re, (float*)ref, n - i);
	}
}
#endif

TARGETAVX int64 GetMinIdxAVX(double* A, int64 n, double& val)
{
#define N 4
	int64 i = 0;
	val = DBL_MAX;
	uint64 idx = (uint64)-1;

	if (n >= N * E256D)
	{
		__m256d min1[N], f[N];
		__m256i midx[N], nidx[N], msep = _mm256_set1_epi64x(N * E256D);
		UNROLL(N) min1[kk] = _mm256_set1_pd(val);
		UNROLL(N) midx[kk] = _mm256_set1_epi64x(0xFFFFFFFFFFFFFFFF);
		UNROLL(N) nidx[kk] = _mm256_set_epi64x(3 + (kk << 2), 2 + (kk << 2), 1 + (kk << 2), 0 + (kk << 2));

		for (int64 l1 = n - N * E256D; i <= l1; i += N * E256D)
		{
			UNROLL(N) 
			{
				f[kk] = _mm256_cmp_pd(min1[kk], _mm256_loadu_pd(A ), _CMP_GT_OS);
				min1[kk] = _mm256_min_pd(min1[kk], _mm256_loadu_pd(A)); 
				A += E256D;
			}

			UNROLL(N) midx[kk] = _mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(midx[kk]), _mm256_castsi256_pd(nidx[kk]), f[kk]));

			UNROLL(N) nidx[kk] = _mm256_add_epi64(nidx[kk], msep);
		}

		REDUCE(min1)
		{
			f[kk] = _mm256_cmp_pd(min1[kk], min1[kk + KK], _CMP_GT_OS);
			min1[kk] = _mm256_min_pd(min1[kk], min1[kk + KK]);
			midx[kk] = _mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(midx[kk]), _mm256_castsi256_pd(midx[kk + KK]), f[kk]));
		}

		for (int64 j = 0; j < E256D; ++j)
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
#define N 4
	int64 i = 0;
	val = FLT_MAX;
	uint idx = (uint)-1;

	if (n >= N * E256F)
	{
		__m256 min1[N], f[N];
		__m256i midx[N], nidx[N], msep = _mm256_set1_epi32(N * E256F);
		UNROLL(N) min1[kk] = _mm256_set1_ps(val);
		UNROLL(N) midx[kk] = _mm256_set1_epi8((char)0xFF);
		UNROLL(N) nidx[kk] = _mm256_set_epi32(7 + (kk << 3), 6 + (kk << 3), 5 + (kk << 3), 4 + (kk << 3), 3 + (kk << 3), 2 + (kk << 3), 1 + (kk << 3), 0 + (kk << 3));

		for (int64 l1 = n - N * E256F; i <= l1; i += N * E256F)
		{
			UNROLL(N) 
			{
				f[kk] = _mm256_cmp_ps(min1[kk], _mm256_loadu_ps(A), _CMP_GT_OS);
				min1[kk] = _mm256_min_ps(min1[kk], _mm256_loadu_ps(A));
				A += E256F; 
			}

			UNROLL(N) midx[kk] = _mm256_castps_si256(_mm256_blendv_ps(_mm256_castsi256_ps(midx[kk]), _mm256_castsi256_ps(nidx[kk]), f[kk]));

			UNROLL(N) nidx[kk] = _mm256_add_epi32(nidx[kk], msep);
		}

		REDUCE(min1)
		{
			f[kk] = _mm256_cmp_ps(min1[kk], min1[kk + KK], _CMP_GT_OS);
			min1[kk] = _mm256_min_ps(min1[kk], min1[kk + KK]);
			midx[kk] = _mm256_castps_si256(_mm256_blendv_ps(_mm256_castsi256_ps(midx[kk]), _mm256_castsi256_ps(midx[kk + KK]), f[kk]));
		}

		for (int64 j = 0; j < E256F; ++j)
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
#define N 4
	int64 i = 0;
	minv = DBL_MAX;
	maxv = -DBL_MAX;

	if (n >= N * E256D)
	{
		__m256d min1[N], max1[N];
		UNROLL(N) min1[kk] = _mm256_set1_pd(minv);
		UNROLL(N) max1[kk] = _mm256_set1_pd(maxv);

		for (int64 l1 = n - N * E256D; i <= l1; i += N * E256D)
		{
			UNROLL(N)
			{
				min1[kk] = _mm256_min_pd(min1[kk], _mm256_loadu_pd(A));
				max1[kk] = _mm256_max_pd(max1[kk], _mm256_loadu_pd(A));
				A += E256D;
			}
		}

		REDUCE(min1)
		{
			min1[kk] = _mm256_min_pd(min1[kk], min1[kk + KK]);
			max1[kk] = _mm256_max_pd(max1[kk], max1[kk + KK]);
		}
		
		minv = _mm256_reduce_min_pd(min1[0]);
		maxv = _mm256_reduce_max_pd(max1[0]);
	}

	for (; i < n; ++i, ++A)
	{
		minv = std::min(minv, *A);
		maxv = std::max(maxv, *A);
	}
}

TARGETAVX void GetMinMaxValAVX(float* A, int64 n, float& minv, float& maxv)
{
#define N 4
	int64 i = 0;
	minv = FLT_MAX;
	maxv = -FLT_MAX;

	if (n >= N * E256F)
	{
		__m256 min1[N], max1[N];
		UNROLL(N) min1[kk] = _mm256_set1_ps(minv);
		UNROLL(N) max1[kk] = _mm256_set1_ps(maxv);

		for (int64 l1 = n - N * E256F; i <= l1; i += N * E256F)
		{
			UNROLL(N)
			{
				min1[kk] = _mm256_min_ps(min1[kk], _mm256_loadu_ps(A));
				max1[kk] = _mm256_max_ps(max1[kk], _mm256_loadu_ps(A));
				A += E256F;
			}
		}

		REDUCE(min1)
		{
			min1[kk] = _mm256_min_ps(min1[kk], min1[kk + KK]);
			max1[kk] = _mm256_max_ps(max1[kk], max1[kk + KK]);
		}
		
		minv = _mm256_reduce_min_ps(min1[0]);
		maxv = _mm256_reduce_max_ps(max1[0]);
	}

	for (; i < n; ++i, ++A)
	{
		minv = std::min(minv, *A);
		maxv = std::max(maxv, *A);
	}
}

TARGETAVX double GetMaxValAVX(double* A, int64 n)
{
#define N 4
	int64 i = 0;
	double val = -DBL_MAX;

	if (n >= N * E256D)
	{
		__m256d max1[N];
		UNROLL(N) max1[kk] = _mm256_set1_pd(val);

		for (int64 l1 = n - N * E256D; i <= l1; i += N * E256D)
		{
			UNROLL(N) { max1[kk] = _mm256_max_pd(max1[kk], _mm256_loadu_pd(A)); A += E256D; }
		}

		REDUCE(max1) max1[kk] = _mm256_max_pd(max1[kk], max1[kk + KK]);
		
		val = _mm256_reduce_max_pd(max1[0]);
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
#define N 4
	int64 i = 0;
	float val = -FLT_MAX;

	if (n >= N * E256F)
	{
		__m256 max1[N], a[N];
		UNROLL(N) max1[kk] = _mm256_set1_ps(val);

		for (int64 l1 = n - N * E256F; i <= l1; i += N * E256F)
		{
			UNROLL(N) { max1[kk] = _mm256_max_ps(max1[kk], _mm256_loadu_ps(A)); A += E256F; }
		}

		REDUCE(max1) max1[kk] = _mm256_max_ps(max1[kk], max1[kk + KK]);
		
		val = _mm256_reduce_max_ps(max1[0]);
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
	//suboptimal to compile
	{
		double val = -DBL_MAX;
		UNROLLHEAD(4)
		for (int64 i = 0; i < n; ++i, A += sep)
			val = std::max(val, *A);
		return val;
	}

#define N 2
	int64 i = 0;
	double val = -DBL_MAX;

	if (n >= N * E256D)
	{
		__m256d max1 = _mm256_set1_pd(val);
		__m256i vindex = _mm256_set_epi64x(3 * sep, 2 * sep, 1 * sep, 0 * sep);

		for (int64 l1 = n - N * E256D; i <= l1; i += N * E256D)
		{
			UNROLL(N)
			{
				max1 = _mm256_max_pd(max1, _mm256_i64gather_pd(A, vindex, sizeof(double)));
				
				A += E256D * sep;
			}
		}

		val = _mm256_reduce_max_pd(max1);
	}

	for (; i < n; ++i, ++A)
	{
		val = std::max(val, *A);
	}

	return val;
}

TARGETAVX float GetMaxValAVX(float* A, int64 n, int64 sep)
{
	//suboptimal to compile
	{
		float val = -FLT_MAX;
		UNROLLHEAD(4)
		for (int64 i = 0; i < n; ++i, A += sep)
			val = std::max(val, *A);
		return val;
	}

#define N 2
	int64 i = 0;
	float val = -FLT_MAX;

	if (n >= N * E256F)
	{
		__m256 max1 = _mm256_set1_ps(val);
		__m256i vindex = _mm256_set_epi32(7 * sep, 6 * sep, 5 * sep, 4 * sep, 3 * sep, 2 * sep, 1 * sep, 0 * sep);

		for (int64 l1 = n - N * E256F; i <= l1; i += N * E256F)
		{
			UNROLL(N)
			{
				max1 = _mm256_max_ps(max1, _mm256_i32gather_ps(A, vindex, sizeof(float)));
				
				A += E256F * sep;
			}
		}

		val = _mm256_reduce_max_ps(max1);
	}

	for (; i < n; ++i, ++A)
	{
		val = std::max(val, *A);
	}

	return val;
}

TARGETAVX double GetMinValAVX(double* A, int64 n)
{
#define N 4
	int64 i = 0;
	double val = DBL_MAX;

	if (n >= N * E256D)
	{
		__m256d min1[N];
		UNROLL(N) min1[kk] = _mm256_set1_pd(val);

		for (int64 l1 = n - N * E256D; i <= l1; i += N * E256D)
		{
			UNROLL(N)
			{
				min1[kk] = _mm256_min_pd(min1[kk], _mm256_loadu_pd(A));
				A += E256D;
			}
		}

		REDUCE(min1) min1[kk] = _mm256_min_pd(min1[kk], min1[kk + KK]);

		val = _mm256_reduce_min_pd(min1[0]);
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
#define N 4
	int64 i = 0;
	float val = FLT_MAX;

	if (n >= N * E256F)
	{
		__m256 min1[N], a[N];
		UNROLL(N) min1[kk] = _mm256_set1_ps(val);

		for (int64 l1 = n - N * E256F; i <= l1; i += N * E256F)
		{
			UNROLL(N) { a[kk] = _mm256_loadu_ps(A); A += E256F; }

			UNROLL(N) { min1[kk] = _mm256_min_ps(min1[kk], a[kk]); }
		}

		REDUCE(min1) min1[kk] = _mm256_min_ps(min1[kk], min1[kk + KK]);
		
		val = _mm256_reduce_min_ps(min1[0]);
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
#define N 4
	int64 i = 0;
	int64 val = 0x7FFFFFFFFFFFFFFF;

	if (n >= N * E256D)
	{
		__m256i min1[N];
		UNROLL(N) min1[kk] = _mm256_set1_epi64x(0x7FFFFFFFFFFFFFFF);

		for (int64 l1 = n - N * E256D; i <= l1; i += N * E256D)
		{
			UNROLL(N) 
			{
				min1[kk] =_mm256_castpd_si256(_mm256_blendv_pd(
					_mm256_castsi256_pd(min1[kk]),
					_mm256_castsi256_pd(_mm256_loadu_si256((__m256i*)A)),
					_mm256_castsi256_pd(_mm256_cmpgt_epi64(min1[kk], _mm256_loadu_si256((__m256i*)A)))));
				A += E256D;
			}
		}

		REDUCE(min1) min1[kk] = _mm256_castpd_si256(_mm256_blendv_pd(
				_mm256_castsi256_pd(min1[kk]),
				_mm256_castsi256_pd(min1[kk + KK]),
				_mm256_castsi256_pd(_mm256_cmpgt_epi64(min1[kk], min1[kk + KK]))));
		
		val = _mm256_reduce_min_epi64(min1[0]);
	}

	for (; i < n; ++i, ++A)
	{
		val = std::min(val, *A);
	}

	return val;
}

TARGETAVX void SetValAVX(uint* A, ushort* B, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E256F)
	{
		for (int64 l1 = n - N * E256F; i <= l1; i += N * E256F)
		{
			UNROLL(N) { _mm256_storeu_si256((__m256i*)A, _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i*)B))); A += E256F; B += E256F; }
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
#define N 4
	int64 i = 0;
	int64 slog = 0; double prod = 1;
	
	if (n >= N * E512D)
	{
		__m256d pd[E512_256];
		UNROLL(E512_256) pd[kk] = _mm256_set1_pd(1.0);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_256)
			{ pd[kk] = _mm256_mul_pd(pd[kk], _mm256_loadu_pd(A)); A += E256D; }

			UNROLL(E512_256) AddExponentAVX(slog, pd[kk]);
		}

		__m128d* pd1 = (__m128d*)pd;
		UNROLL(E512_128) ChargeLogSSE(slog, prod, pd1[kk]);
	}

	for (; i < n; ++i, ++A)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);

	return prod;
}

TARGETAVX double LogProdAVX(float* A, int64 n)
{
#define N 4
	int64 i = 0;
	int64 slog = 0; double prod = 1;
	
	if (n >= N * E512D)
	{
		__m256d pd[E512_256];
		UNROLL(E512_256) pd[kk] = _mm256_set1_pd(1.0);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_256)
			{ pd[kk] = _mm256_mul_pd(pd[kk], _mm256_cvtps_pd(_mm_loadu_ps(A))); A += E128F; }
				
			UNROLL(E512_256) AddExponentAVX(slog, pd[kk]);
		}

		__m128d* pd1 = (__m128d*)pd;
		UNROLL(E512_128) ChargeLogSSE(slog, prod, pd1[kk]);
	}

	for (; i < n; ++i, ++A)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);

	return prod;
}

TARGETAVX double LogProdAVX(double* A, int64 n, int64 sep)
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
		__m256d pd[E512_256];
		UNROLL(E512_256) pd[kk] = _mm256_set1_pd(1.0);
		__m256i vindex = _mm256_set_epi64x(3 * sep, 2 * sep, 1 * sep, 0 * sep);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_256) 
			{ pd[kk] = _mm256_mul_pd(pd[kk], _mm256_i64gather_pd(A, vindex, sizeof(double))); A += E256D * sep; }

			UNROLL(E512_256) AddExponentAVX(slog, pd[kk]);
		}
		
		__m128d* pd1 = (__m128d*)&pd;
		UNROLL(E512_128) ChargeLogSSE(slog, prod, pd1[kk]);
	}

	for (; i < n; ++i, A += sep)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);

	return prod;
}

TARGETAVX double LogProdAVX(float* A, int64 n, int64 sep)
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
}

TARGETAVX double LogProdDivAVX(double* A, double* B, int64 n, int64 sep)
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

TARGETAVX double LogProdDivAVX(float* A, float* B, int64 n, int64 sep)
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

TARGETAVX int64 CountNonZeroAVX(byte* A, int64 n)
{
	uint64 re = 0;
	int64 i = 0;

	if (n >= E256B)
	{
		__m256i a = _mm256_setzero_si256(), z = _mm256_setzero_si256(), o = _mm256_set1_epi8(0x01);

		for (int64 l1 = n - E256B; i <= l1; i += E256B)
		{
			a = _mm256_add_epi64(a, _mm256_sad_epu8(z, _mm256_min_epu8(_mm256_loadu_si256((__m256i*)A), o)));
			A += E256B;
		}

		re = _mm256_reduce_add_epi64(a);
	}

	for (; i < n; ++i, ++A)
		if (*A) re++;

	return (int64)re;
}

TARGETAVX double SumAVX(double* A, int64 n)
{
#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		__m256d s[E512_256] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_256) 
			{ s[kk] = _mm256_add_pd(s[kk], _mm256_loadu_pd(A)); A += E256D; }
		}

		REDUCE(s) s[kk] = _mm256_add_pd(s[kk], s[kk + KK]);

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
#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		__m256d s[E512_256] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_256)
			{ s[kk] = _mm256_add_pd(s[kk], _mm256_cvtps_pd(_mm_loadu_ps(A))); A += E128F; }
		}

		REDUCE(s) s[kk] = _mm256_add_pd(s[kk], s[kk + KK]);
		
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
#define N 4
	int64 i = 0;
	volatile float re = 0;

	if (n >= N * E512F)
	{
		__m256 s[E512_256] = { 0 };

		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) UNROLL(E512_256) 
			{ s[kk] = _mm256_add_ps(s[kk], _mm256_loadu_ps(A)); A += E256F; }
		}

		REDUCE(s) s[kk] = _mm256_add_ps(s[kk], s[kk + KK]);

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
#define N 4
	int64 re = 0;
	int64 i = 0;

	if (n >= N * E256B)
	{
		__m256i s = _mm256_setzero_si256(), z = _mm256_setzero_si256();

		for (int64 l1 = n - N * E256B; i <= l1; i += N * E256B)
		{
			UNROLL(N) { s = _mm256_add_epi64(s, _mm256_sad_epu8(_mm256_loadu_si256((__m256i*)A), z)); A += E256B;}
		}

		re += _mm256_reduce_add_epi64(s);
	}

	for (; i < n; ++i)
		re += *A++;

	return re;
}

TARGETAVX double SumAVX(double* A, int64 n, int64 sep)
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
		__m256d s[E512_256] = { 0 };
		__m256i vindex = _mm256_set_epi64x(3 * sep, 2 * sep, 1 * sep, 0 * sep);
		UNROLL(E256D) _mm_prefetch((const char*)&A[kk * sep], _MM_HINT_T0);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_256)
			{
				UNROLL(E256D) _mm_prefetch((const char*)&A[E256D * sep + kk * sep], _MM_HINT_T0);
				s[kk] = _mm256_add_pd(s[kk], _mm256_i64gather_pd(A, vindex, sizeof(double)));
				A += E256D * sep;
			}
		}

		REDUCE(s) s[kk] = _mm256_add_pd(s[kk], s[kk + KK]);

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
		__m256d s[E512_256] = { 0 };
		__m128i vindex = _mm_set_epi32(3 * sep, 2 * sep, 1 * sep, 0 * sep);
		UNROLL(E256D) _mm_prefetch((const char*)&A[kk * sep], _MM_HINT_T0);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_256)
			{
				UNROLL(E256D) _mm_prefetch((const char*)&A[E256D * sep + kk * sep], _MM_HINT_T0);
				s[kk] = _mm256_add_pd(s[kk], _mm256_cvtps_pd(_mm_i32gather_ps(A, vindex, sizeof(float))));
				A += E256D * sep;
			}
		}

		REDUCE(s) s[kk] = _mm256_add_pd(s[kk], s[kk + KK]);

		re = _mm256_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, A += sep)
	{
		volatile double v1 = *A;
		re += v1;
	}

	return re;
}

TARGETAVX float SumAVXx(float* A, int64 n, int64 sep)
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
		__m256 s[E512_256] = { 0 };
		__m256i vindex = _mm256_set_epi32(7 * sep, 6 * sep, 5 * sep, 4 * sep, 3 * sep, 2 * sep, 1 * sep, 0 * sep);
		UNROLL(E256F) _mm_prefetch((const char*)&A[kk * sep], _MM_HINT_T0);

		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) UNROLL(E512_256)
			{
				UNROLL(E256F) _mm_prefetch((const char*)&A[E256F * sep + kk * sep], _MM_HINT_T0);
				s[kk] = _mm256_add_ps(s[kk], _mm256_i32gather_ps(A, vindex, sizeof(float)));
				A += E256F * sep;
			}
		}

		REDUCE(s) s[kk] = _mm256_add_ps(s[kk], s[kk + KK]);

		re = _mm256_reduce_add_ps(s[0]);
	}

	for (; i < n; ++i, A += sep)
	{
		volatile float v1 = *A;
		re += v1;
	}

	return re;
}

TARGETAVX double ProdAVX(double* A, int64 n)
{
#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		__m256d s[E512_256];
		UNROLL(E512_256) s[kk] = _mm256_set1_pd(1);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_256)
			{ s[kk] = _mm256_mul_pd(s[kk], _mm256_loadu_pd(A)); A += E256D; }
		}

		REDUCE(s) s[kk] = _mm256_mul_pd(s[kk], s[kk + KK]);

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
#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		__m256d s[E512_256];
		UNROLL(E512_256) s[kk] = _mm256_set1_pd(1);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_256)
			{ s[kk] = _mm256_mul_pd(s[kk], _mm256_cvtps_pd(_mm_loadu_ps(A))); A += E128F; }
		}

		REDUCE(s) s[kk] = _mm256_mul_pd(s[kk], s[kk + KK]);

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
#define N 4
	int64 i = 0;
	volatile float re = 0;

	if (n >= N * E512F)
	{
		__m256 s[E512_256];
		UNROLL(E512_256) s[kk] = _mm256_set1_ps(1);

		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) UNROLL(E512_256)
			{ s[kk] = _mm256_mul_ps(s[kk], _mm256_loadu_ps(A)); A += E256F; }
		}

		REDUCE(s) s[kk] = _mm256_mul_ps(s[kk], s[kk + KK]);

		re = _mm256_reduce_mul_ps(s[0]);
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A++;
		re *= v1;
	}

	return re;
}

TARGETAVX double ProdAVX(double* A, int64 n, int64 sep)
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
		__m256d s[E512_256];
		UNROLL(E512_256) s[kk] = _mm256_set1_pd(1);
		__m256i vindex = _mm256_set_epi64x(3 * sep, 2 * sep, 1 * sep, 0 * sep);
		UNROLL(E256D) _mm_prefetch((const char*)&A[kk * sep], _MM_HINT_T0);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_256)
			{
				UNROLL(E256D) _mm_prefetch((const char*)&A[E256D * sep + kk * sep], _MM_HINT_T0);
				s[kk] = _mm256_mul_pd(s[kk], _mm256_i64gather_pd(A, vindex, sizeof(double)));
				A += E256D * sep;
			}
		}

		REDUCE(s) s[kk] = _mm256_mul_pd(s[kk], s[kk + KK]);

		re = _mm256_reduce_mul_pd(s[0]);
	}

	for (; i < n; ++i, A += sep)
	{
		volatile double v1 = *A;
		re *= v1;
	}

	return re;
}

TARGETAVX double ProdAVX(float* A, int64 n, int64 sep)
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
		__m256d s[E512_256];
		UNROLL(E512_256) s[kk] = _mm256_set1_pd(1);
		__m128i vindex = _mm_set_epi32(3 * sep, 2 * sep, 1 * sep, 0 * sep);
		UNROLL(E256D) _mm_prefetch((const char*)&A[kk * sep], _MM_HINT_T0);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_256)
			{
				UNROLL(E256D) _mm_prefetch((const char*)&A[E256D * sep + kk * sep], _MM_HINT_T0);
				s[kk] = _mm256_mul_pd(s[kk], _mm256_cvtps_pd(_mm_i32gather_ps(A, vindex, sizeof(float))));
				A += E256D * sep;
			}
		}

		REDUCE(s) s[kk] = _mm256_mul_pd(s[kk], s[kk + KK]);

		re = _mm256_reduce_mul_pd(s[0]);
	}

	for (; i < n; ++i, A += sep)
	{
		volatile double v1 = *A;
		re *= v1;
	}

	return re;
}

TARGETAVX float ProdAVXx(float* A, int64 n, int64 sep)
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
		__m256 s[E512_256];
		UNROLL(E512_256) s[kk] = _mm256_set1_ps(1);
		__m256i vindex = _mm256_set_epi32(7 * sep, 6 * sep, 5 * sep, 4 * sep, 3 * sep, 2 * sep, 1 * sep, 0 * sep);
		UNROLL(E256F) _mm_prefetch((const char*)&A[kk * sep], _MM_HINT_T0);

		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) UNROLL(E512_256)
			{
				UNROLL(E256F) _mm_prefetch((const char*)&A[E256F * sep + kk * sep], _MM_HINT_T0);
				s[kk] = _mm256_mul_ps(s[kk], _mm256_i32gather_ps(A, vindex, sizeof(float)));
				A += E256F * sep;
			}
		}

		REDUCE(s) s[kk] = _mm256_mul_ps(s[kk], s[kk + KK]);

		re = _mm256_reduce_mul_ps(s[0]);
	}

	for (; i < n; ++i, A += sep)
	{
		volatile float v1 = *A;
		re *= v1;
	}

	return re;
}

TARGETAVX double SumSquareAVX(double* A, int64 n)
{
#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		__m256d s[E512_256] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_256) 
			{ s[kk] = _mm256_fmaddx_pd(_mm256_loadu_pd(A), _mm256_loadu_pd(A), s[kk]); A += E256D; }
		}

		REDUCE(s) s[kk] = _mm256_add_pd(s[kk], s[kk + KK]);

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
#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		__m256d s[E512_256] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_256)
			{
				__m256d v1 = _mm256_cvtps_pd(_mm_loadu_ps(A)); A += E128F;
				s[kk] = _mm256_fmaddx_pd(v1, v1, s[kk]);
			}
		}

		REDUCE(s) s[kk] = _mm256_add_pd(s[kk], s[kk + KK]);

		re = _mm256_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, ++A)
	{
		volatile double v1 = (double)*A * (double)*A;
		re += v1;
	}

	return re;
}

TARGETAVX int64 SumSquareAVX(byte* A, int64 n)
{
#define N 2
	int64 i = 0;
	uint64 re = 0;

	if (n >= N * E256B)
	{
		__m256i t = _mm256_setzero_si256(), s[2] = { 0 };
		__m128i* s2 = (__m128i*)&s;

		for (int64 l1 = n - N * E256B; i <= l1; i += N * E256B)
		{
			UNROLL(N) { s[0] = _mm256_add_epi16(s[0], _mm256_maddubs_epi16(_mm256_loadu_si256((__m256i*)A), _mm256_loadu_si256((__m256i*)A))); A += E256B; }

			if ((i & (E256B * 128 - 1)) == 0) [[unlikely]]
			{
				s[1] = _mm256_srli_si256(s[0], 8);
				UNROLL(E512_128) t = _mm256_add_epi64(t, _mm256_cvtepi16_epi64(s2[kk]));
				s[0] = _mm256_setzero_si256();
			}
		}
		
		s[1] = _mm256_srli_si256(s[0], 8);
		UNROLL(E512_128) t = _mm256_add_epi64(t, _mm256_cvtepi16_epi64(s2[kk]));
		re = _mm256_reduce_add_epi64(t);
	}

	for (; i < n; ++i, ++A)
		re += *A * *A;

	return re;
}

TARGETAVX void SumSumSquareAVX(double* A, int64 n, double& sum, double& sumsq)
{
#define N 4
	int64 i = 0;
	volatile double re1 = 0, re2 = 0;

	if (n >= N * E512D)
	{
		__m256d s1[E512_256] = { 0 }, s2[E512_256] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_256)
			{
				__m256d v1 = _mm256_loadu_pd(A); A += E256D;

				s1[kk] = _mm256_add_pd(v1, s1[kk]);

				s2[kk] = _mm256_fmaddx_pd(v1, v1, s2[kk]);
			}
		}

		REDUCE(s1)
		{
			s1[kk] = _mm256_add_pd(s1[kk], s1[kk + KK]);
			s2[kk] = _mm256_add_pd(s2[kk], s2[kk + KK]);
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
#define N 4
	int64 i = 0;
	volatile double re1 = 0, re2 = 0;

	if (n >= N * E512D)
	{
		__m256d s1[E512_256] = { 0 }, s2[E512_256] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_256)
			{
				__m256d v1 = _mm256_cvtps_pd(_mm_loadu_ps(A)); A += E128F;

				s1[kk] = _mm256_add_pd(v1, s1[kk]);

				s2[kk] = _mm256_fmaddx_pd(v1, v1, s2[kk]);
			}
		}

		REDUCE(s1)
		{
			s1[kk] = _mm256_add_pd(s1[kk], s1[kk + KK]);
			s2[kk] = _mm256_add_pd(s2[kk], s2[kk + KK]);
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
		UNROLL(E256D) { _mm_prefetch((const char*)B, _MM_HINT_T0); B += sep; }

		__m256d s1[E512_256] = { 0 }, s2[E512_256] = { 0 };
		__m256i vindex = _mm256_set_epi64x(-5 * sep, -6 * sep, -7 * sep, -8 * sep);
		
		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_256)
			{
				UNROLL(E256D) { _mm_prefetch((const char*)B, _MM_HINT_T0); B += sep; }

				__m256d b = _mm256_i64gather_pd(B, vindex, sizeof(double));
				
				s1[kk] = _mm256_fmaddx_pd(b, _mm256_loadu_pd(A1), s1[kk]); A1 += E256D;

				s2[kk] = _mm256_fmaddx_pd(b, _mm256_loadu_pd(A2), s2[kk]); A2 += E256D;
			}
		}

		REDUCE(s1)
		{
			s1[kk] = _mm256_add_pd(s1[kk], s1[kk + KK]);
			s2[kk] = _mm256_add_pd(s2[kk], s2[kk + KK]);
		}

		re1 = _mm256_reduce_add_pd(s1[0]);
		re2 = _mm256_reduce_add_pd(s2[0]);

		B -= E256D * sep;
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

TARGETAVX double SumProdDivAVX(double* A1, float* A2, float* B, int64 sep, int64 n)
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
		UNROLL(E256D) { _mm_prefetch((const char*)B, _MM_HINT_T0); B += sep; }

		__m256d s1[E512_256] = { 0 }, s2[E512_256] = { 0 };
		__m128i vindex = _mm_set_epi32(-5 * sep, -6 * sep, -7 * sep, -8 * sep);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_256)
			{
				UNROLL(E256D) { _mm_prefetch((const char*)B, _MM_HINT_T0); B += sep; }

				__m256d b = _mm256_cvtps_pd(_mm_i32gather_ps(B, vindex, sizeof(float)));

				s1[kk] = _mm256_fmaddx_pd(_mm256_loadu_pd(A1), b, s1[kk]); A1 += E256D; 

				s2[kk] = _mm256_fmaddx_pd(_mm256_cvtps_pd(_mm_loadu_ps(A2)), b, s2[kk]); A2 += E256D; 
			}
		}

		REDUCE(s1)
		{
			s1[kk] = _mm256_add_pd(s1[kk], s1[kk + KK]);
			s2[kk] = _mm256_add_pd(s2[kk], s2[kk + KK]);
		}

		re1 = _mm256_reduce_add_pd(s1[0]);
		re2 = _mm256_reduce_add_pd(s2[0]);

		B -= E256D * sep;
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

TARGETAVX double SumProdDivAVX(float* A1, float* A2, float* B, int64 sep, int64 n)
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
		UNROLL(E256D) { _mm_prefetch((const char*)B, _MM_HINT_T0); B += sep; }

		__m256d s1[E512_256] = { 0 }, s2[E512_256] = { 0 };
		__m128i vindex = _mm_set_epi32(-5 * sep, -6 * sep, -7 * sep, -8 * sep);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_256)
			{
				UNROLL(E256D) { _mm_prefetch((const char*)B, _MM_HINT_T0); B += sep; }

				__m256d b = _mm256_cvtps_pd(_mm_i32gather_ps(B, vindex, sizeof(float)));

				s1[kk] = _mm256_fmaddx_pd(_mm256_cvtps_pd(_mm_loadu_ps(A1)), b, s1[kk]); A1 += E256D;

				s2[kk] = _mm256_fmaddx_pd(_mm256_cvtps_pd(_mm_loadu_ps(A2)), b, s2[kk]); A1 += E256D;
			}
		}

		REDUCE(s1)
		{
			s1[kk] = _mm256_add_pd(s1[kk], s1[kk + KK]);
			s2[kk] = _mm256_add_pd(s2[kk], s2[kk + KK]);
		}

		re1 = _mm256_reduce_add_pd(s1[0]);
		re2 = _mm256_reduce_add_pd(s2[0]);

		B -= E256D * sep;
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

TARGETAVX float SumProdDivAVXx(float* A1, float* A2, float* B, int64 sep, int64 n)
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

TARGETAVX double SumProdAVX(double* A, double* B, int64 sep, int64 n)
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
		UNROLL(E512D) { _mm_prefetch((const char*)B, _MM_HINT_T0); B += sep; }

		__m256d s[E512_256] = { 0 };
		__m256i vindex = _mm256_set_epi64x(-9 * sep, -10 * sep, -11 * sep, -12 * sep);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N)
			{
				UNROLL(E512_256)
				{
					UNROLL(E256D) { _mm_prefetch((const char*)B, _MM_HINT_T0); B += sep; }
					
					s[kk] = _mm256_fmaddx_pd(_mm256_loadu_pd(A), _mm256_i64gather_pd(B, vindex, sizeof(double)), s[kk]); A += E256D;
				} 
			}
		}

		REDUCE(s) s[kk] = _mm256_add_pd(s[kk], s[kk + KK]); 

		re = _mm256_reduce_add_pd(s[0]); 

		B -= E512D * sep;
	}

	for (; i < n; ++i, ++A, B += sep)
	{
		volatile double v1 = *A * *B;
		re += v1;
	}

	return re;
}

TARGETAVX double SumProdAVX(float* A, float* B, int64 sep, int64 n)
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
		UNROLL(E256D) { _mm_prefetch((const char*)B, _MM_HINT_T0); B += sep; }

		__m256d s[E512_256] = { 0 };
		__m128i vindex = _mm_set_epi32(-5 * sep, -6 * sep, -7 * sep, -8 * sep);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_256)
			{
				UNROLL(E256D) { _mm_prefetch((const char*)B, _MM_HINT_T0); B += sep; }

				s[kk] = _mm256_fmaddx_pd(_mm256_cvtps_pd(_mm_loadu_ps(A)), _mm256_cvtps_pd(_mm_i32gather_ps(B, vindex, sizeof(float))), s[kk]);
				
				A += E256D;
			}
		}

		REDUCE(s) s[kk] = _mm256_add_pd(s[kk], s[kk + KK]);

		re = _mm256_reduce_add_pd(s[0]);

		B -= E128F * sep;
	}

	for (; i < n; ++i, A++, B += sep)
	{
		volatile double v1 = (double)*A * (double)*B;
		re += v1;
	}

	return re;
}

TARGETAVX float SumProdAVXx(float* A, float* B, int64 sep, int64 n)
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

TARGETAVX double SumProdAVX(double* A, double* B, int64 n)
{
#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		__m256d s[E512_256] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_256) 
			{ s[kk] = _mm256_fmaddx_pd(_mm256_loadu_pd(A), _mm256_loadu_pd(B), s[kk]); A += E256D; B += E256D; }
		}

		REDUCE(s) s[kk] = _mm256_add_pd(s[kk], s[kk + KK]);

		re = _mm256_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, ++A, ++B)
	{
		volatile double v1 = (double)*A * (double)*B;
		re += v1;
	}

	return re;
}

TARGETAVX double SumProdAVX(float* A, float* B, int64 n)
{
#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		__m256d s[E512_256] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_256) 
			{ s[kk] = _mm256_fmaddx_pd(_mm256_cvtps_pd(_mm_loadu_ps(A)), _mm256_cvtps_pd(_mm_loadu_ps(B)), s[kk]); A += E256D; B += E256D; }
		}

		REDUCE(s) s[kk] = _mm256_add_pd(s[kk], s[kk + KK]);

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
#define N 4
	int64 i = 0;
	volatile float re = 0;

	if (n >= N * E512F)
	{
		__m256 s[E512_256] = { 0 };

		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) UNROLL(E512_256) 
			{ s[kk] = _mm256_fmaddx_ps(_mm256_loadu_ps(A), _mm256_loadu_ps(B), s[kk]); A += E256F; B += E256F; }
		}

		REDUCE(s) s[kk] = _mm256_add_ps(s[kk], s[kk + KK]);

		re = _mm256_reduce_add_ps(s[0]);
	}

	for (; i < n; ++i, ++A, ++B)
	{
		volatile float v1 = (float)*A * (float)*B;
		re += v1;
	}

	return re;
}

TARGETAVX double SumProdAVX(double* A, double* B, double* C, int64 n)
{
#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		__m256d s[E512_256] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_256)
			{ s[kk] = _mm256_fmaddx_pd(_mm256_mul_pd(_mm256_loadu_pd(A), _mm256_loadu_pd(B)), _mm256_loadu_pd(C), s[kk]); A += E256D; B += E256D; C += E256D; }
		}

		REDUCE(s) s[kk] = _mm256_add_pd(s[kk], s[kk + KK]);

		re = _mm256_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, ++A, ++B, ++C)
	{
		volatile double v1 = (double)*A * (double)*B * (double)*C;
		re += v1;
	}

	return re;
}

TARGETAVX float SumProdAVX(float* A, float* B, float* C, int64 n)
{
#define N 4
	int64 i = 0;
	volatile float re = 0;

	if (n >= N * E512F)
	{
		__m256 s[E512_256] = { 0 };

		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) UNROLL(E512_256) 
			{ s[kk] = _mm256_fmaddx_ps(_mm256_mul_ps(_mm256_loadu_ps(A), _mm256_loadu_ps(B)), _mm256_loadu_ps(C), s[kk]); A += E256F; B += E256F; C += E256F; }
		}

		REDUCE(s) s[kk] = _mm256_add_ps(s[kk], s[kk + KK]);

		re = _mm256_reduce_add_ps(s[0]);
	}

	for (; i < n; ++i, ++A, ++B, ++C)
	{
		volatile float v1 = (float)*A * (float)*B * (float)*C;
		re += v1;
	}

	return re;
}

TARGETAVX double SumSqProdAVX(double* A, double* B, int64 n)
{
#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		__m256d s[E512_256] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_256) 
			{ s[kk] = _mm256_fmaddx_pd(_mm256_mul_pd(_mm256_loadu_pd(A), _mm256_loadu_pd(A)), _mm256_loadu_pd(B), s[kk]); A += E256D; B += E256D; }
		}

		REDUCE(s) s[kk] = _mm256_add_pd(s[kk], s[kk + KK]);

		re = _mm256_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, ++A, ++B)
	{
		volatile double v1 = (double)*A * (double)*A * (double)*B;
		re += v1;
	}

	return re;
}

TARGETAVX float SumSqProdAVX(float* A, float* B, int64 n)
{
#define N 4
	int64 i = 0;
	volatile float re = 0;

	if (n >= N * E512F)
	{
		__m256 s[E512_256] = { 0 };

		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) UNROLL(E512_256) 
			{ s[kk] = _mm256_fmaddx_ps(_mm256_mul_ps(_mm256_loadu_ps(A), _mm256_loadu_ps(A)), _mm256_loadu_ps(B), s[kk]); A += E256F; B += E256F; }
		}

		REDUCE(s) s[kk] = _mm256_add_ps(s[kk], s[kk + KK]);

		re = _mm256_reduce_add_ps(s[0]);
	}

	for (; i < n; ++i, ++A, ++B)
	{
		volatile float v1 = (float)*A * (float)*A * (float)*B;
		re += v1;
	}

	return re;
}

TARGETAVX void AddAVX(double* A, double* B, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E256D)
	{
		for (int64 l1 = n - N * E256D; i <= l1; i += N * E256D)
		{
			UNROLL(N) { _mm256_storeu_pd(A, _mm256_add_pd(_mm256_loadu_pd(A), _mm256_loadu_pd(B))); A += E256D; B += E256D; }
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGETAVX void AddAVX(float* A, float* B, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E256F)
	{
		for (int64 l1 = n - N * E256F; i <= l1; i += N * E256F)
		{
			UNROLL(N) { _mm256_storeu_ps(A, _mm256_add_ps(_mm256_loadu_ps(A), _mm256_loadu_ps(B))); A += E256F; B += E256F; }
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGETAVX void AddAVX(int64* A, int64* B, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E256D)
	{
		for (int64 l1 = n - N * E256D; i <= l1; i += N * E256D)
		{
			UNROLL(N) { _mm256_storeu_si256((__m256i*)A, _mm256_add_epi64(_mm256_loadu_si256((__m256i*)A), _mm256_loadu_si256((__m256i*)B))); A += E256D; B += E256D; }
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGETAVX void AddAVX(int* A, int* B, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E256F)
	{
		for (int64 l1 = n - N * E256F; i <= l1; i += N * E256F)
		{
			UNROLL(N) { _mm256_storeu_si256((__m256i*)A, _mm256_add_epi32(_mm256_loadu_si256((__m256i*)A), _mm256_loadu_si256((__m256i*)B))); A += E256F; B += E256F; }
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGETAVX void AddAVX(int* A, int B, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E256F)
	{
		__m256i b = _mm256_set1_epi32(B);

		for (int64 l1 = n - N * E256F; i <= l1; i += N * E256F)
		{
			UNROLL(N) { _mm256_storeu_si256((__m256i*)A, _mm256_add_epi32(_mm256_loadu_si256((__m256i*)A), b)); A += E256F; }
		}
	}

	for (; i < n; ++i, A++)
		*A += B;
}

TARGETAVX void AddAVX(double* A, double B, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E256D)
	{
		__m256d b = _mm256_set1_pd(B);

		for (int64 l1 = n - N * E256D; i <= l1; i += N * E256D)
		{
			UNROLL(N) { _mm256_storeu_pd(A, _mm256_add_pd(_mm256_loadu_pd(A), b)); A += E256D; }
		}
	}

	for (; i < n; ++i, A++)
		*A += B;
}

TARGETAVX void AddAVX(float* A, float B, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E256F)
	{
		__m256 b = _mm256_set1_ps(B);

		for (int64 l1 = n - N * E256F; i <= l1; i += N * E256F)
		{
			UNROLL(N) { _mm256_storeu_ps(A, _mm256_add_ps(_mm256_loadu_ps(A), b)); A += E256F; }
		}
	}

	for (; i < n; ++i, A++)
		*A += B;
}

TARGETAVX void MulAVX(double* A, double* B, double* C, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E256D)
	{
		for (int64 l1 = n - N * E256D; i <= l1; i += N * E256D)
		{
			UNROLL(N) { _mm256_storeu_pd(A, _mm256_mul_pd(_mm256_loadu_pd(B), _mm256_loadu_pd(C))); A += E256D; B += E256D; C += E256D; }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *B++ * *C++;
		*A++ = v1;
	}
}

TARGETAVX void MulAVX(float* A, float* B, float* C, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E256F)
	{
		for (int64 l1 = n - N * E256F; i <= l1; i += N * E256F)
		{
			UNROLL(N) { _mm256_storeu_ps(A, _mm256_mul_ps(_mm256_loadu_ps(B), _mm256_loadu_ps(C))); A += E256F; B += E256F; C += E256F; }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *B++ * *C++;
		*A++ = v1;
	}
}

TARGETAVX void MulAVX(double* A, double* B, double C, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E256D)
	{
		__m256d c = _mm256_set1_pd(C);

		for (int64 l1 = n - N * E256D; i <= l1; i += N * E256D)
		{
			UNROLL(N) { _mm256_storeu_pd(A, _mm256_mul_pd(_mm256_loadu_pd(B), c)); A += E256D; B += E256D; }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *B++ * C;
		*A++ = v1;
	}
}

TARGETAVX void MulAVX(float* A, float* B, float C, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E256F)
	{
		__m256 c = _mm256_set1_ps(C);

		for (int64 l1 = n - N * E256F; i <= l1; i += N * E256F)
		{
			UNROLL(N) { _mm256_storeu_ps(A, _mm256_mul_ps(_mm256_loadu_ps(B), c)); A += E256F; B += E256F; }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *B++ * C;
		*A++ = v1;
	}
}

TARGETAVX void MulAVX(double* A, double B, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E256D)
	{
		__m256d b = _mm256_set1_pd(B);

		for (int64 l1 = n - N * E256D; i <= l1; i += N * E256D)
		{
			UNROLL(N) { _mm256_storeu_pd(A, _mm256_mul_pd(_mm256_loadu_pd(A), b)); A += E256D; }
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
#define N 4
	int64 i = 0;

	if (n >= N * E256F)
	{
		__m256 b = _mm256_set1_ps(B);

		for (int64 l1 = n - N * E256F; i <= l1; i += N * E256F)
		{
			UNROLL(N) { _mm256_storeu_ps(A, _mm256_mul_ps(_mm256_loadu_ps(A), b)); A += E256F; }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A * B;
		*A++ = v1;
	}
}

TARGETAVX void DivAVX(double* A, double B, double* C, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E256D)
	{
		__m256d b = _mm256_set1_pd(B);

		for (int64 l1 = n - N * E256D; i <= l1; i += N * E256D)
		{
			UNROLL(N) { _mm256_storeu_pd(A, _mm256_div_pd(b, _mm256_loadu_pd(C))); A += E256D; C += E256D; }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = B / *C++;
		*A++ = v1;
	}
}

TARGETAVX void DivAVX(float* A, float B, float* C, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E256F)
	{
		__m256 b = _mm256_set1_ps(B);

		for (int64 l1 = n - N * E256F; i <= l1; i += N * E256F)
		{
			UNROLL(N) { _mm256_storeu_ps(A, _mm256_div_ps(b, _mm256_loadu_ps(C))); A += E256F; C += E256F; }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = B / *C++;
		*A++ = v1;
	}
}

TARGETAVX void DivAVX(double* A, double* B, double* C, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E256D)
	{
		for (int64 l1 = n - N * E256D; i <= l1; i += N * E256D)
		{
			UNROLL(N) { _mm256_storeu_pd(A, _mm256_div_pd(_mm256_loadu_pd(B), _mm256_loadu_pd(C))); A += E256D; B += E256D; C += E256D; }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *B++ / *C++;
		*A++ = v1;
	}
}

TARGETAVX void DivAVX(float* A, float* B, float* C, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E256F)
	{
		for (int64 l1 = n - N * E256F; i <= l1; i += N * E256F)
		{
			UNROLL(N) { _mm256_storeu_ps(A, _mm256_div_ps(_mm256_loadu_ps(B), _mm256_loadu_ps(C))); A += E256F; B += E256F; C += E256F; }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *B++ / *C++;
		*A++ = v1;
	}
}

TARGETAVX void AddProdAVX(double* A, double* B, double* C, int64 n)
{
#define N 1
	int64 i = 0;

	if (n >= N * E256D)
	{
		for (int64 l1 = n - N * E256D; i <= l1; i += N * E256D)
		{
			UNROLL(N) { _mm256_storeu_pd(A, _mm256_fmaddx_pd(_mm256_loadu_pd(B), _mm256_loadu_pd(C), _mm256_loadu_pd(A))); A += E256D; B += E256D; C += E256D; }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *B++ * *C++;
		*A++ += v1;
	}
}

TARGETAVX void AddProdAVX(float* A, float* B, float* C, int64 n)
{
#define N 1
	int64 i = 0;

	if (n >= N * E256F)
	{
		for (int64 l1 = n - N * E256F; i <= l1; i += N * E256F)
		{
			UNROLL(N) { _mm256_storeu_ps(A, _mm256_fmaddx_ps(_mm256_loadu_ps(B), _mm256_loadu_ps(C), _mm256_loadu_ps(A))); A += E256F; B += E256F; C += E256F; }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *B++ * *C++;
		*A++ += v1;
	}
}

TARGETAVX void AddProdAVX(double* A, double* B, double C, int64 n)
{
#define N 1
	int64 i = 0;

	if (n >= N * E256D)
	{
		__m256d c = _mm256_set1_pd(C);

		for (int64 l1 = n - N * E256D; i <= l1; i += N * E256D)
		{
			UNROLL(N) { _mm256_storeu_pd(A, _mm256_fmaddx_pd(_mm256_loadu_pd(B), c, _mm256_loadu_pd(A))); A += E256D; B += E256D; }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *B++ * C;
		*A++ += v1;
	}
}

TARGETAVX void AddProdAVX(double* A, float* B, double C, int64 n)
{
#define N 1
	int64 i = 0;

	if (n >= N * E256D)
	{
		__m256d c = _mm256_set1_pd(C);

		for (int64 l1 = n - N * E256D; i <= l1; i += N * E256D)
		{
			UNROLL(N) { _mm256_storeu_pd(A, _mm256_fmaddx_pd(_mm256_cvtps_pd(_mm_loadu_ps(B)), c, _mm256_loadu_pd(A))); A += E256D; B += E256D; }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *B++ * C;
		*A++ += v1;
	}
}

TARGETAVX void AddProdAVX(float* A, float* B, float C, int64 n)
{
#define N 1
	int64 i = 0;

	if (n >= N * E256F)
	{
		__m256 c = _mm256_set1_ps(C);

		for (int64 l1 = n - N * E256F; i <= l1; i += N * E256F)
		{
			UNROLL(N) { _mm256_storeu_ps(A, _mm256_fmaddx_ps(_mm256_loadu_ps(B), c, _mm256_loadu_ps(A))); A += E256F; B += E256F; }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *B++ * C;
		*A++ += v1;
	}
}

TARGETAVX void UnifyAVX(double* A, int64 n)
{
#define N 4
	int64 i = 0;
	double invsum = 1.0 / (SumAVX(A, n) + n * MIN_FREQ);

	if (n >= E256D)
	{
		__m256d b = _mm256_set1_pd(invsum), c = _mm256_set1_pd(MIN_FREQ * invsum);
		
		UNROLLHEAD(N)
		for (int64 l1 = n - E256D; i <= l1; i += E256D)
		{
			_mm256_storeu_pd(A, _mm256_fmaddx_pd(_mm256_loadu_pd(A), b, c)); 
			A += E256D;
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
#define N 4
	int64 i = 0;
	double invsum = 1.0 / (SumAVX(A, n) + n * MIN_FREQ);

	if (n >= E256D)
	{
		__m256d b = _mm256_set1_pd(invsum), c = _mm256_set1_pd(MIN_FREQ * invsum);

		UNROLLHEAD(N)
		for (int64 l1 = n - E256D; i <= l1; i += E256D)
		{
			_mm_storeu_ps(A, _mm256_cvtpd_ps(_mm256_fmaddx_pd(_mm256_cvtps_pd(_mm_loadu_ps(A)), b, c))); 
			A += E256D;
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
#define N 2
	A++; n--;
	int64 i = 0;

	if (n >= E512B)
	{
		__m256i r[N], v = _mm256_set1_epi8(val), o = _mm256_setzero_si256();

		for (int64 l1 = n - E512B; i <= l1; i += E512B)
		{
			char* Ab = A;
			UNROLL(N) { r[kk] = _mm256_cmpeq_epi8(_mm256_loadu_si256((__m256i*)A), v); A += E256B; }

			__m256i rz = _mm256_or_si256(r[0], r[1]);

			if (_mm256_testz_si256(rz, rz)) continue;

			uint64 mask; uint* mask2 = (uint*) & mask;
			UNROLL(N) mask2[kk] = _mm256_movemask_epi8(r[kk]);
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

TARGETAVX int64 CountCharAVX(char* A, char val, int64 n)
{
#define N 4
	uint64 re = 0;
	int64 i = 0;

	if (n >= E256B)
	{
		__m256i r = _mm256_setzero_si256(), v = _mm256_set1_epi8(val), o = _mm256_setzero_si256();
		
		UNROLLHEAD(N)
		for (int64 l1 = n - E256B; i <= l1; i += E256B)
		{
			r = _mm256_add_epi64(r, _mm256_sad_epu8(o, _mm256_cmpeq_epi8(_mm256_loadu_si256((__m256i*)A), v)));
			A += E256B;
		}

		re = _mm256_reduce_add_epi64(r) / 0xFF;
	}

	for (; i < n; ++i, A++)
		if (*A == val) re++;

	return (int64)re;
}

TARGETAVX void DiagQuadFormAVX(double* res1, double* A, double* D, int64 m, int64 n)
{
#define N 3

#define DECLARE			double* pA1 = (A + n * i), *pA2 = (A + n * j), *pD = D; __m256d a1[N], a2[N], r[N][N] = { 0 }
#define ALOAD1(ii)		a1[ii] = _mm256_loadu_pd(pA1 + n * ii)
#define ALOAD2(ii)		a2[ii] = _mm256_loadu_pd(pA2 + n * ii)
#define ALOAD3(ii)		a1[ii] = a2[ii] = _mm256_loadu_pd(pA1 + n * ii)
#define ADMUL(ii)		a1[ii] = _mm256_mul_pd(a1[ii], _mm256_loadu_pd(pD))
#define FMADD(ii,jj)	r[ii][jj] = _mm256_fmaddx_pd(a1[ii], a2[jj], r[ii][jj])
#define RDUADD(ii,jj)	res1[(i+ii) * m + (j+jj)] = _mm256_reduce_add_pd(r[ii][jj])
#define REMAIN(ii,jj)	res1[(i+ii) * m + (j+jj)] += A[(i+ii) * n + k] * A[(j+jj) * n + k] * D[k]
#define FINAL(ii,jj)	res1[(i+ii) + (j+jj) * m] = res1[(i+ii) * m + (j+jj)]
	
    int64 i = 0, j = 0, k = 0;
    for (i = 0; i + N <= m; i += N)
	{
        for (j = 0; j < i; j += N)
		{
			DECLARE;

            for (k = 0; k + E256D <= n; k += E256D, pA1 += E256D, pA2 += E256D, pD += E256D) 
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

            for (k = 0; k + E256D <= n; k += E256D, pA1 += E256D, pD += E256D) 
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

            for (k = 0; k + E256D <= n; k += E256D, pA1 += E256D, pA2 += E256D, pD += E256D) 
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

            for (k = 0; k + E256D <= n; k += E256D, pA1 += E256D, pD += E256D) 
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

TARGETAVX void DiagQuadFormAVX(double* res2, double* A, double* D, double* B, int64 m, int64 n)
{
#define N 8

#define DECLARE			double* pA = (A + n * i), *pB = B, *pD = D; __m256d bd, r[N] = { 0 }
#define BDMUL			bd = _mm256_mul_pd(_mm256_loadu_pd(pB), _mm256_loadu_pd(pD)); pB += E256D; pD += E256D
#define FMADD(ii)		r[ii] = _mm256_fmaddx_pd(_mm256_loadu_pd(pA + n * ii), bd, r[ii]); 
#define RDUADD(ii)		res2[(i+ii)] = _mm256_reduce_add_pd(r[ii])
#define FINAL(ii)		res2[(i+ii)] += A[(i+ii) * n + k] * B[k] * D[k]

    int64 i = 0, k = 0;
    for (i = 0; i + N <= m; i += N)
	{
		DECLARE;

        for (k = 0; k + E256D <= n; k += E256D, pA += E256D) 
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

        for (k = 0; k + E256D <= n; k += E256D, pA += E256D) 
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

TARGETAVX void DiagQuadFormAVX(double* res3, double* B, double* D, int64 n)
{
#define N 4
	__m256d s[4] = { 0 };

	int64 k = 0;
    for (; k + N * E256D <= n; k += N * E256D) 
	{
		UNROLL(N) 
		{
			s[kk] = _mm256_fmaddx_pd(_mm256_mul_pd(_mm256_loadu_pd(B), _mm256_loadu_pd(B)), _mm256_loadu_pd(D), s[kk]); 
			B += E256D;
			D += E256D;
		}
	}

	REDUCE(s) s[kk] = _mm256_add_pd(s[kk], s[kk + KK]);

	res3[0] = _mm256_reduce_add_pd(s[0]);
			
	VECTORIZE
    for (; k < n; k++, B++, D++) 
		res3[0] += B[0] * B[0] * D[0];
}

TARGETAVX void MatrixMulAVX(double* res, double* A, double* B, int64 m, int64 n, int64 p)
{
#define N 3
	
#define DECLARE			double* pA = (A + n * i), *pB = (B + n * j); __m256d a[N], b[N], r[N][N] = { 0 }
#define ALOAD(kk)		a[kk] = _mm256_loadu_pd(pA + n * kk)
#define BLOAD(kk)		b[kk] = _mm256_loadu_pd(pB + n * kk)
#define FMADD(ii,jj)	r[ii][jj] = _mm256_fmaddx_pd(a[ii], b[jj], r[ii][jj])
#define RDUADD(ii,jj)	res[(i+ii) + (j+jj) * m] = _mm256_reduce_add_pd(r[ii][jj])
#define REMAIN(ii,jj)	res[(i+ii) + (j+jj) * m] += A[(i+ii) * n + k] * B[(j+jj) * n + k]

    int64 i = 0, j = 0, k = 0;
    for (i = 0; i + N <= m; i += N)
	{
        for (j = 0; j + N <= p; j += N)
		{
			DECLARE;

            for (k = 0; k + E256D <= n; k += E256D, pA += E256D, pB += E256D) 
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

            for (k = 0; k + E256D <= n; k += E256D, pA += E256D, pB += E256D) 
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

            for (k = 0; k + E256D <= n; k += E256D, pA += E256D, pB += E256D) 
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

            for (k = 0; k + E256D <= n; k += E256D, pA += E256D, pB += E256D) 
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

TARGETAVX void DiagQuadFormAVX(float* res1, float* A, float* D, int64 m, int64 n)
{
#define N 3

#define DECLARE			float* pA1 = (A + n * i), *pA2 = (A + n * j), *pD = D; __m256 a1[N], a2[N], r[N][N] = { 0 }
#define ALOAD1(ii)		a1[ii] = _mm256_loadu_ps(pA1 + n * ii)
#define ALOAD2(ii)		a2[ii] = _mm256_loadu_ps(pA2 + n * ii)
#define ALOAD3(ii)		a1[ii] = a2[ii] = _mm256_loadu_ps(pA1 + n * ii)
#define ADMUL(ii)		a1[ii] = _mm256_mul_ps(a1[ii], _mm256_loadu_ps(pD))
#define FMADD(ii,jj)	r[ii][jj] = _mm256_fmaddx_ps(a1[ii], a2[jj], r[ii][jj])
#define RDUADD(ii,jj)	res1[(i+ii) * m + (j+jj)] = _mm256_reduce_add_ps(r[ii][jj])
#define REMAIN(ii,jj)	res1[(i+ii) * m + (j+jj)] += A[(i+ii) * n + k] * A[(j+jj) * n + k] * D[k]
#define FINAL(ii,jj)	res1[(i+ii) + (j+jj) * m] = res1[(i+ii) * m + (j+jj)]
	
    int64 i = 0, j = 0, k = 0;
    for (i = 0; i + N <= m; i += N)
	{
        for (j = 0; j < i; j += N)
		{
			DECLARE;

            for (k = 0; k + E256F <= n; k += E256F, pA1 += E256F, pA2 += E256F, pD += E256F) 
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

            for (k = 0; k + E256F <= n; k += E256F, pA1 += E256F, pD += E256F) 
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

            for (k = 0; k + E256F <= n; k += E256F, pA1 += E256F, pA2 += E256F, pD += E256F) 
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

            for (k = 0; k + E256F <= n; k += E256F, pA1 += E256F, pD += E256F) 
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

TARGETAVX void DiagQuadFormAVX(float* res2, float* A, float* D, float* B, int64 m, int64 n)
{
#define N 8

#define DECLARE			float* pA = (A + n * i), *pB = B, *pD = D; __m256 bd, r[N] = { 0 }
#define BDMUL			bd = _mm256_mul_ps(_mm256_loadu_ps(pB), _mm256_loadu_ps(pD)); pB += E256F; pD += E256F
#define FMADD(ii)		r[ii] = _mm256_fmaddx_ps(_mm256_loadu_ps(pA + n * ii), bd, r[ii]); 
#define RDUADD(ii)		res2[(i+ii)] = _mm256_reduce_add_ps(r[ii])
#define FINAL(ii)		res2[(i+ii)] += A[(i+ii) * n + k] * B[k] * D[k]

    int64 i = 0, k = 0;
    for (i = 0; i + N <= m; i += N)
	{
		DECLARE;

        for (k = 0; k + E256F <= n; k += E256F, pA += E256F) 
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

        for (k = 0; k + E256F <= n; k += E256F, pA += E256F) 
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

TARGETAVX void DiagQuadFormAVX(float* res3, float* B, float* D, int64 n)
{
#define N 4
	__m256 s[4] = { 0 };

	int64 k = 0;
    for (; k + N * E256F <= n; k += N * E256F) 
	{
		UNROLL(N) 
		{
			s[kk] = _mm256_fmaddx_ps(_mm256_mul_ps(_mm256_loadu_ps(B), _mm256_loadu_ps(B)), _mm256_loadu_ps(D), s[kk]); 
			B += E256F; 
			D += E256F; 
		}
	}

	REDUCE(s) s[kk] = _mm256_add_ps(s[kk], s[kk + KK]);

	res3[0] = _mm256_reduce_add_ps(s[0]);
			
	VECTORIZE
    for (; k < n; k++, B++, D++) 
		res3[0] += B[0] * B[0] * D[0];
}

TARGETAVX void MatrixMulAVX(float* res, float* A, float* B, int64 m, int64 n, int64 p)
{
#define N 3
	
#define DECLARE			float* pA = (A + n * i), *pB = (B + n * j); __m256 a[N], b[N], r[N][N] = { 0 }
#define ALOAD(kk)		a[kk] = _mm256_loadu_ps(pA + n * kk)
#define BLOAD(kk)		b[kk] = _mm256_loadu_ps(pB + n * kk)
#define FMADD(ii,jj)	r[ii][jj] = _mm256_fmaddx_ps(a[ii], b[jj], r[ii][jj])
#define RDUADD(ii,jj)	res[(i+ii) + (j+jj) * m] = _mm256_reduce_add_ps(r[ii][jj])
#define REMAIN(ii,jj)	res[(i+ii) + (j+jj) * m] += A[(i+ii) * n + k] * B[(j+jj) * n + k]

    int64 i = 0, j = 0, k = 0;
    for (i = 0; i + N <= m; i += N)
	{
        for (j = 0; j + N <= p; j += N)
		{
			DECLARE;

            for (k = 0; k + E256F <= n; k += E256F, pA += E256F, pB += E256F) 
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

            for (k = 0; k + E256F <= n; k += E256F, pA += E256F, pB += E256F) 
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

            for (k = 0; k + E256F <= n; k += E256F, pA += E256F, pB += E256F) 
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

            for (k = 0; k + E256F <= n; k += E256F, pA += E256F, pB += E256F) 
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
