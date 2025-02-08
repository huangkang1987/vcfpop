/* AVX512 Instruction Set Functions */

#include "vcfpop.h"

#ifndef __aarch64__

template struct RNG512<double>;
template struct RNG512<float >;
template TARGET512 void RNG512<double>::Integer<uint  >(uint  * re, int64 n, uint   minv, uint   maxv);
template TARGET512 void RNG512<double>::Integer<uint64>(uint64* re, int64 n, uint64 minv, uint64 maxv);
template TARGET512 void RNG512<float >::Integer<uint  >(uint  * re, int64 n, uint   minv, uint   maxv);
template TARGET512 void RNG512<float >::Integer<uint64>(uint64* re, int64 n, uint64 minv, uint64 maxv);

#ifndef _RNG512_FP64
#define XS 8
#define XS2 16

/* Initialize rng */
TARGET512 RNG512<double>::RNG512()
{

}

/* Initialize rng */
TARGET512 RNG512<double>::RNG512(uint64 seed, uint64 salt)
{
	__m512i a[XS], s, m;

	UNROLL(XS) { a[kk] = _mm512_set_epi64(seed + 7, seed + 6, seed + 5, seed + 4, seed + 3, seed + 2, seed + 1, seed + 0); seed += 8; }

	s = _mm512_set1_epi64(salt);
	m = _mm512_set1_epi32(0x5bd1e995);

	UNROLL(XS) a[kk] = _mm512_xor_si512(a[kk], _mm512_slli_epi64(_mm512_andnot_si512(a[kk], _mm512_set1_epi64(0xFFFFFFFF)), 32));

	s = _mm512_xor_si512(s, _mm512_slli_epi64(_mm512_andnot_si512(s, _mm512_set1_epi64(0xFFFFFFFF)), 32));

	// uint s = s ^ 4;
	s = _mm512_xor_si512(s, _mm512_set1_epi32(4));

	// a *= m;
	UNROLL(XS) a[kk] = _mm512_mullo_epi32(a[kk], m);//OK

	// a ^= a >> 24;
	UNROLL(XS) a[kk] = _mm512_xor_si512(a[kk], _mm512_srli_epi32(a[kk], 24));

	// a *= m;
	UNROLL(XS) a[kk] = _mm512_mullo_epi32(a[kk], m);//OK

	// s *= m;
	s = _mm512_mullo_epi32(s, m);//OK

	// a ^= s;
	UNROLL(XS) a[kk] = _mm512_xor_si512(a[kk], s);

	// a ^= a >> 13;
	UNROLL(XS) a[kk] = _mm512_xor_si512(a[kk], _mm512_srli_epi32(a[kk], 13));

	// a *= m;
	UNROLL(XS) a[kk] = _mm512_mullo_epi32(a[kk], m);//OK

	// a ^= a >> 15;
	UNROLL(XS) a[kk] = _mm512_xor_si512(a[kk], _mm512_srli_epi32(a[kk], 15));

	// original
	UNROLL(XS) x[kk] = _mm512_xor_si512(_mm512_set1_epi64(0x159A55E5005BCD15), a[kk]);

	UNROLL(XS) a[kk] = _mm512_slli_epi64(a[kk], 6);

	UNROLL(XS) y[kk] = _mm512_xor_si512(_mm512_set1_epi64(0x054913331F123BB5), a[kk]);
}

/* Draw 64 64-bit integers in [0,n), 64*n frequencies are in arr */
TARGET512 void RNG512<double>::Poly(__m512d* arr, int n, __m512i* re)
{
	__m512d t[XS], s[XS];
	__m512d one = _mm512_set1_pd(1.0);
	__m512i mask1 = _mm512_set1_epi64(0x000FFFFFFFFFFFFF);
	__m512i mask2 = _mm512_set1_epi64(0x3FF0000000000000);
	__m512i* r = (__m512i*)t;

	UNROLL(XS) s[kk] = _mm512_setzero_pd();

	for (int i = 0; i < n * XS; i += XS)
		UNROLL(XS) s[kk] = _mm512_add_pd(s[kk], arr[kk + i]);

	XorShift();

	UNROLL(XS) r[kk] = _mm512_add_epi64(x[kk], y[kk]);

	UNROLL(XS) r[kk] = _mm512_and_si512(r[kk], mask1);

	UNROLL(XS) r[kk] = _mm512_or_si512(r[kk], mask2);

	UNROLL(XS) t[kk] = _mm512_sub_pd(t[kk], one);

	UNROLL(XS) t[kk] = _mm512_mul_pd(t[kk], s[kk]);

	__m512i midx[XS], nidx[XS], ninc[XS];
	__mmask8 b[XS], f[XS] = { 0 };

	UNROLL(XS) midx[kk] = _mm512_set1_epi64(n - 1);
	UNROLL(XS) nidx[kk] = _mm512_setzero_si512();
	UNROLL(XS) ninc[kk] = _mm512_set1_epi64(1);

	for (int i = 0; i < n * XS; i += XS)
	{
		UNROLL(XS) b[kk] = _mm512_cmp_pd_mask(t[kk], arr[kk + i], _CMP_LT_OS);

		UNROLL(XS) t[kk] = _mm512_sub_pd(t[kk], arr[kk + i]);

		UNROLL(XS) b[kk] = (~f[kk]) & b[kk];

		UNROLL(XS) f[kk] = f[kk] | b[kk];

		UNROLL(XS) midx[kk] = _mm512_mask_blend_epi64(b[kk], midx[kk], nidx[kk]);//ok

		UNROLL(XS) nidx[kk] = _mm512_add_epi64(nidx[kk], ninc[kk]);
	}

	UNROLL(XS) re[kk] = midx[kk];
}

/* Draw uniform distriubted intergers */
TARGET512 void RNG512<double>::XorShift()
{
	__m512i a[XS], b[XS];

	UNROLL(XS) a[kk] = x[kk];

	UNROLL(XS) b[kk] = y[kk];

	UNROLL(XS) x[kk] = b[kk];

	UNROLL(XS) a[kk] = _mm512_xor_si512(a[kk], _mm512_slli_epi64(a[kk], 23));

	UNROLL(XS) a[kk] = _mm512_xor_si512(a[kk], _mm512_srli_epi64(a[kk], 18));

	UNROLL(XS) a[kk] = _mm512_xor_si512(a[kk], b[kk]);

	UNROLL(XS) a[kk] = _mm512_xor_si512(a[kk], _mm512_srli_epi64(b[kk], 5));

	UNROLL(XS) y[kk] = a[kk];
}

/* Draw uniform distriubted integers */
template<typename INT>
TARGET512 void RNG512<double>::Integer(INT* re, int64 n, INT minv, INT maxv)
{
	constexpr int xesize = sizeof(x) / sizeof(INT);
	int64 i = 0;
	INT modv = maxv - minv;
	INT* rei = re;

	for (; i <= n - xesize; i += xesize)
	{
		XorShift();
		UNROLL(XS) { _mm512_storeu_si512((__m512i*)rei, _mm512_add_epi64(x[kk], y[kk])); rei += E512B / sizeof(INT); };
	}

	if (i != n)
	{
		__m512i re2[XS];
		XorShift();
		UNROLL(XS) re2[kk] = _mm512_add_epi64(x[kk], y[kk]);
		SetVal((INT*)rei, (INT*)re2, n - i);
	}

	if (maxv != (INT)-1 || minv != 0)
	{
		for (i = 0; i < n; ++i)
			re[i] = re[i] % modv + minv;
	}
}

/* Draw uniform distriubted real numbers */
TARGET512 void RNG512<double>::Uniform(double* re, int n, double minv, double maxv)
{
	constexpr int xesize = sizeof(x) / sizeof(double);
	int i = 0;
	double range = maxv - minv;

	__m512i mask1 = _mm512_set1_epi64(0x000FFFFFFFFFFFFF);
	__m512i mask2 = _mm512_set1_epi64(0x3FF0000000000000);
	__m512d v1 = _mm512_set1_pd(minv - range);
	__m512d v2 = _mm512_set1_pd(range);

	if (range == 1.0)
	{
		for (; i <= n - xesize; i += xesize)
		{
			XorShift();
			UNROLL(XS) { _mm512_storeu_pd(re, _mm512_add_pd(_mm512_castsi512_pd(_mm512_or_si512(_mm512_and_si512(_mm512_add_epi64(x[kk], y[kk]), mask1), mask2)), v1)); re += E512D; }
		}
	}
	else
	{
		for (; i <= n - xesize; i += xesize)
		{
			XorShift();
			UNROLL(XS) { _mm512_storeu_pd(re, _mm512_fmaddx_pd(_mm512_castsi512_pd(_mm512_or_si512(_mm512_and_si512(_mm512_add_epi64(x[kk], y[kk]), mask1), mask2)), v2, v1)); re += E512D; }
		}
	}

	if (i != n)
	{
		__m512d ref[XS];
		XorShift();
		UNROLL(XS) ref[kk] = _mm512_fmaddx_pd(_mm512_castsi512_pd(_mm512_or_si512(_mm512_and_si512(_mm512_add_epi64(x[kk], y[kk]), mask1), mask2)), v2, v1);
		SetVal((double*)re, (double*)ref, n - i);
	}
}

/* Draw uniform distriubted real numbers */
TARGET512 void RNG512<double>::Normal(double* re, int n, double mean, double sd)
{
	constexpr int xhsize = XS / 2;
	constexpr int xesize = XS * E512D;

	int i = 0;

	__m512i mask1 = _mm512_set1_epi64(0x000FFFFFFFFFFFFF);
	__m512i mask2 = _mm512_set1_epi64(0x3FF0000000000000);
	__m512d v1 = _mm512_set1_pd(-1);
	__m512d min_freq = _mm512_set1_pd(MIN_FREQ);
	__m512d pi2 = _mm512_set1_pd(2.0 * M_PI);
	__m512d mu = _mm512_set1_pd(mean);
	__m512d s = _mm512_set1_pd(sd);

	for (; i <= n - xesize; i += xesize)
	{
		XorShift();
		UNROLL(XS) _mm512_storeu_pd(re + E512D * kk, _mm512_add_pd(_mm512_castsi512_pd(_mm512_or_si512(_mm512_and_si512(_mm512_add_epi64(x[kk], y[kk]), mask1), mask2)), v1));

		__m512d u1, u2, u3, u4;
		for (int j = 0; j < xhsize; ++j)
		{
			u1 = _mm512_max_pd(_mm512_loadu_pd(re + j * E512D), min_freq);
			u2 = _mm512_mul_pd(_mm512_loadu_pd(re + j * E512D + xhsize * E512D), pi2);

			UNROLL(E512D) simd_f64(u1, kk) = sqrt(-2.0 * log(simd_f64(u1, kk)));
			UNROLL(E512D) simd_f64(u3, kk) = cos(simd_f64(u2, kk));
			UNROLL(E512D) simd_f64(u4, kk) = sin(simd_f64(u2, kk));
				
			_mm512_storeu_pd(re + j * E512D, _mm512_mul_pd(u1, u3));
			_mm512_storeu_pd(re + j * E512D + xhsize * E512D, _mm512_mul_pd(u1, u4));
		}
			
		if (sd != 1 || mean != 0)
			UNROLL(XS) { _mm512_storeu_pd(re, _mm512_fmaddx_pd(_mm512_loadu_pd(re), s, mu)); re += E512D; }
		else
			re += XS * E512D;
	}

	if (i != n)
	{
		double ref[XS * E512D];

		XorShift();
		UNROLL(XS) _mm512_storeu_pd(ref + E512D * kk, _mm512_add_pd(_mm512_castsi512_pd(_mm512_or_si512(_mm512_and_si512(_mm512_add_epi64(x[kk], y[kk]), mask1), mask2)), v1));

		__m512d u1, u2, u3, u4;
		for (int j = 0; j < xhsize; ++j)
		{
			u1 = _mm512_max_pd(_mm512_loadu_pd(ref + j * E512D), min_freq);
			u2 = _mm512_mul_pd(_mm512_loadu_pd(ref + j * E512D + xhsize * E512D), pi2);

			UNROLL(E512D) simd_f64(u1, kk) = sqrt(-2.0 * log(simd_f64(u1, kk)));
			UNROLL(E512D) simd_f64(u3, kk) = cos(simd_f64(u2, kk));
			UNROLL(E512D) simd_f64(u4, kk) = sin(simd_f64(u2, kk));
			
			_mm512_storeu_pd(ref + j * E512D, _mm512_mul_pd(u1, u3));
			_mm512_storeu_pd(ref + j * E512D + xhsize * E512D, _mm512_mul_pd(u1, u4));
		}
		
		if (sd != 1 || mean != 0)
			UNROLL(XS) _mm512_storeu_pd(ref + kk * E512D, _mm512_fmaddx_pd(_mm512_loadu_pd(re), s, mu)); 

		SetVal((double*)re, (double*)ref, n - i);
	}
}
#endif

#ifndef _RNG512_FP32
#define XS 4
#define XS2 8

/* Initialize rng */
TARGET512 RNG512<float>::RNG512()
{

}

/* Initialize rng */
TARGET512 RNG512<float>::RNG512(uint64 seed, uint64 salt)
{
	__m512i a[XS], s, m;
	UNROLL(XS) { a[kk] = _mm512_set_epi32(Mix(seed + 15), Mix(seed + 14), Mix(seed + 13), Mix(seed + 12), Mix(seed + 11), Mix(seed + 10), Mix(seed + 9), Mix(seed + 8), Mix(seed + 7), Mix(seed + 6), Mix(seed + 5), Mix(seed + 4), Mix(seed + 3), Mix(seed + 2), Mix(seed + 1), Mix(seed + 0)); seed += 16; }

	s = _mm512_set1_epi32(Mix(salt));
	m = _mm512_set1_epi32(0x5bd1e995);

	// uint s = s ^ 4;
	s = _mm512_xor_si512(s, _mm512_set1_epi32(4));

	// a *= m;
	UNROLL(XS) a[kk] = _mm512_mullo_epi32(a[kk], m);//OK

	// a ^= a >> 24;
	UNROLL(XS) a[kk] = _mm512_xor_si512(a[kk], _mm512_srli_epi32(a[kk], 24));

	// a *= m;
	UNROLL(XS) a[kk] = _mm512_mullo_epi32(a[kk], m);//OK

	// s *= m;
	s = _mm512_mullo_epi32(s, m);//OK

	// a ^= s;
	UNROLL(XS) a[kk] = _mm512_xor_si512(a[kk], s);

	// a ^= a >> 13;
	UNROLL(XS) a[kk] = _mm512_xor_si512(a[kk], _mm512_srli_epi32(a[kk], 13));

	// a *= m;
	UNROLL(XS) a[kk] = _mm512_mullo_epi32(a[kk], m);//OK

	// a ^= a >> 15;
	UNROLL(XS) a[kk] = _mm512_xor_si512(a[kk], _mm512_srli_epi32(a[kk], 15));

	// original
	UNROLL(XS) x[kk] = _mm512_xor_si512(_mm512_set1_epi32(0x005BCD15), a[kk]);

	UNROLL(XS) a[kk] = _mm512_slli_epi32(a[kk], 3);

	UNROLL(XS) y[kk] = _mm512_xor_si512(_mm512_set1_epi32(0x159A55E5), a[kk]);

	UNROLL(XS) a[kk] = _mm512_slli_epi32(a[kk], 3);

	UNROLL(XS) z[kk] = _mm512_xor_si512(_mm512_set1_epi32(0x1F123BB5), a[kk]);
}

/* Draw 64 64-bit integers in [0,n), 64*n frequencies are in arr */
TARGET512 void RNG512<float>::Poly(__m512* arr, int n, __m512i* re)
{
	__m512 t[XS], s[XS]; 
	__m512 one = _mm512_set1_ps(1.0f);
	__m512i mask1 = _mm512_set1_epi32(0x000FFFFF);
	__m512i mask2 = _mm512_set1_epi32(0x3F800000);
	__m512i* r = (__m512i*)t;

	UNROLL(XS) s[kk] = _mm512_setzero_ps();

	for (int i = 0; i < n * XS; i += XS)
		UNROLL(XS) s[kk] = _mm512_add_ps(s[kk], arr[kk + i]);

	XorShift();

	UNROLL(XS) r[kk] = _mm512_and_si512(z[kk], mask1);

	UNROLL(XS) r[kk] = _mm512_or_si512(r[kk], mask2);

	UNROLL(XS) t[kk] = _mm512_sub_ps(t[kk], one);

	UNROLL(XS) t[kk] = _mm512_mul_ps(t[kk], s[kk]);

	__m512i midx[XS] = { _mm512_set1_epi32(n - 1), _mm512_set1_epi32(n - 1), _mm512_set1_epi32(n - 1), _mm512_set1_epi32(n - 1) };
	__m512i nidx[XS] = { _mm512_setzero_si512(), _mm512_setzero_si512(), _mm512_setzero_si512(), _mm512_setzero_si512() };
	__m512i ninc[XS] = { _mm512_set1_epi32(1), _mm512_set1_epi32(1), _mm512_set1_epi32(1), _mm512_set1_epi32(1) };
	__mmask16 b[XS], f[XS] = { 0 };

	for (int i = 0; i < n * XS; i += XS)
	{
		UNROLL(XS) b[kk] = _mm512_cmp_ps_mask(t[kk], arr[kk + i], _CMP_LT_OS);

		UNROLL(XS) t[kk] = _mm512_sub_ps(t[kk], arr[kk + i]);

		UNROLL(XS) b[kk] = (~f[kk]) & b[kk];

		UNROLL(XS) f[kk] = f[kk] | b[kk];

		//GCC has a bug here, cannot correct move nidx to midx according to the mask b2
		UNROLL(XS) midx[kk] = _mm512_mask_mov_epi32(midx[kk], b[kk], nidx[kk]);

		UNROLL(XS) nidx[kk] = _mm512_add_epi32(nidx[kk], ninc[kk]);
	}

	__m256i* midx2 = (__m256i*)midx;
	UNROLL(XS2) re[kk] = _mm512_cvtepi32_epi64(midx2[kk]);
}

/* Draw uniform distriubted intergers */
TARGET512 void RNG512<float>::XorShift()
{
	__m512i u[XS];

	UNROLL(XS) u[kk] = _mm512_slli_epi32(x[kk], 16);
	UNROLL(XS) x[kk] = _mm512_xor_si512(x[kk], u[kk]);

	UNROLL(XS) u[kk] = _mm512_srli_epi32(x[kk], 5);
	UNROLL(XS) x[kk] = _mm512_xor_si512(x[kk], u[kk]);

	UNROLL(XS) u[kk] = _mm512_slli_epi32(x[kk], 1);
	UNROLL(XS) x[kk] = _mm512_xor_si512(x[kk], u[kk]);

	UNROLL(XS) u[kk] = x[kk];

	UNROLL(XS) x[kk] = y[kk];

	UNROLL(XS) y[kk] = z[kk];

	UNROLL(XS) z[kk] = _mm512_xor_si512(u[kk], x[kk]);

	UNROLL(XS) z[kk] = _mm512_xor_si512(z[kk], y[kk]);
}

/* Draw uniform distriubted integers */
template<typename INT>
TARGET512 void RNG512<float>::Integer(INT* re, int64 n, INT minv, INT maxv)
{
	constexpr int xesize = sizeof(x) / sizeof(INT);
	int64 i = 0;
	INT modv = maxv - minv;
	INT* rei = re;

	for (; i <= n - xesize; i += xesize)
	{
		XorShift();
		UNROLL(XS) { _mm512_storeu_si512((__m512i*)rei, z[kk]); rei += E512B / sizeof(INT); }
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
TARGET512 void RNG512<float>::Uniform(float* re, int n, float minv, float maxv)
{
	constexpr int xesize = sizeof(x) / sizeof(float);

	int i = 0;
	float range = maxv - minv;

	__m512i mask1 = _mm512_set1_epi32(0x000FFFFF);
	__m512i mask2 = _mm512_set1_epi32(0x3F800000);
	__m512 v1 = _mm512_set1_ps(minv - range);
	__m512 v2 = _mm512_set1_ps(range);

	if (range == 1.0)
	{
		for (; i <= n - xesize; i += xesize)
		{
			XorShift();
			UNROLL(XS) { _mm512_storeu_ps(re, _mm512_add_ps(_mm512_castsi512_ps(_mm512_or_si512(_mm512_and_si512(z[kk], mask1), mask2)), v1)); re += E512F; }
		}
	}
	else
	{
		for (; i <= n - xesize; i += xesize)
		{
			XorShift();
			UNROLL(XS) { _mm512_storeu_ps(re, _mm512_fmaddx_ps(_mm512_castsi512_ps(_mm512_or_si512(_mm512_and_si512(z[kk], mask1), mask2)), v2, v1)); re += E512F; }
		}
	}

	if (i != n)
	{
		__m512 ref[XS];
		XorShift();
		UNROLL(XS) ref[kk] = _mm512_fmaddx_ps(_mm512_castsi512_ps(_mm512_or_si512(_mm512_and_si512(z[kk], mask1), mask2)), v2, v1);
		SetVal((float*)re, (float*)ref, n - i);
	}
}

/* Draw uniform distriubted real numbers */
TARGET512 void RNG512<float>::Normal(float* re, int n, float mean, float sd)
{
	constexpr int xhsize = XS / 2;
	constexpr int xesize = XS * E512F;

	int i = 0;

	__m512i mask1 = _mm512_set1_epi32(0x000FFFFF);
	__m512i mask2 = _mm512_set1_epi32(0x3F800000);
	__m512 v1 = _mm512_set1_ps(-1);
	__m512 min_freq = _mm512_set1_ps((float)MIN_FREQ);
	__m512 pi2 = _mm512_set1_ps((float)(2.0 * M_PI));
	__m512 mu = _mm512_set1_ps(mean);
	__m512 s = _mm512_set1_ps(sd);

	for (; i <= n - xesize; i += xesize)
	{
		XorShift();
		UNROLL(XS) _mm512_storeu_ps(re + kk * E256F, _mm512_add_ps(_mm512_castsi512_ps(_mm512_or_si512(_mm512_and_si512(z[kk], mask1), mask2)), v1));
		
		__m512 u1, u2, u3, u4;
		for (int j = 0; j < xhsize; ++j)
		{
			u1 = _mm512_max_ps(_mm512_loadu_ps(re + j * E512F), min_freq);
			u2 = _mm512_mul_ps(_mm512_loadu_ps(re + j * E512F + xhsize * E512F), pi2);

			UNROLL(E512F) simd_f32(u1, kk) = sqrt(-2.0 * log(simd_f32(u1, kk)));
			UNROLL(E512F) simd_f32(u3, kk) = cos(simd_f32(u2, kk));
			UNROLL(E512F) simd_f32(u4, kk) = sin(simd_f32(u2, kk));
			
			_mm512_storeu_ps(re + j * E512F, _mm512_mul_ps(u1, u3));
			_mm512_storeu_ps(re + j * E512F + xhsize * E512F, _mm512_mul_ps(u1, u4));
		}
		
		if (sd != 1 || mean != 0)
			UNROLL(XS) { _mm512_storeu_ps(re, _mm512_fmaddx_ps(_mm512_loadu_ps(re), s, mu)); re += E512F; }
		else
			re += XS * E512F;
	}

	if (i != n)
	{
		float ref[XS * E512F];

		XorShift();
		UNROLL(XS) _mm512_storeu_ps(ref + kk * E512F, _mm512_add_ps(_mm512_castsi512_ps(_mm512_or_si512(_mm512_and_si512(z[kk], mask1), mask2)), v1));

		__m512 u1, u2, u3, u4;
		for (int j = 0; j < xhsize; ++j)
		{
			u1 = _mm512_max_ps(_mm512_loadu_ps(ref + j * E512F), min_freq);
			u2 = _mm512_mul_ps(_mm512_loadu_ps(ref + j * E512F + xhsize * E512F), pi2);

			UNROLL(E512F) simd_f32(u1, kk) = sqrt(-2.0 * log(simd_f32(u1, kk)));
			UNROLL(E512F) simd_f32(u3, kk) = cos(simd_f32(u2, kk));
			UNROLL(E512F) simd_f32(u4, kk) = sin(simd_f32(u2, kk));
			
			_mm512_storeu_ps(&ref[j * E512F], _mm512_mul_ps(u1, u3));
			_mm512_storeu_ps(&ref[j * E512F + xhsize * E512F], _mm512_mul_ps(u1, u4));
		}
		
		if (sd != 1 || mean != 0)
			UNROLL(XS) _mm512_storeu_ps(ref + kk * E512F, _mm512_fmaddx_ps(_mm512_loadu_ps(re), s, mu)); 

		SetVal((float*)re, (float*)ref, n - i);
	}
}
#endif

TARGET512 int64 GetMinIdx512(double* A, int64 n, double& val)
{
#define N 2
	int64 i = 0;
	val = DBL_MAX;
	uint64 idx = (uint64)-1;

	if (n >= N * E512D)
	{
		__m512d min1[N];
		__m512i midx[N], nidx[N], msep = _mm512_set1_epi64(N * E512D);
		UNROLL(N) min1[kk] = _mm512_set1_pd(val);
		UNROLL(N) midx[kk] = _mm512_set1_epi64(0xFFFFFFFFFFFFFFFF);
		UNROLL(N) nidx[kk] = _mm512_set_epi64(7 + (kk << 3), 6 + (kk << 3), 5 + (kk << 3), 4 + (kk << 3), 3 + (kk << 3), 2 + (kk << 3), 1 + (kk << 3), 0 + (kk << 3));

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) 
			{
				midx[kk] = _mm512_mask_mov_epi64(midx[kk], _mm512_cmp_pd_mask(min1[kk], _mm512_loadu_pd(A), _CMP_GT_OS), nidx[kk]);
				min1[kk] = _mm512_min_pd(min1[kk], _mm512_loadu_pd(A));
				A += E512D;
			}

			UNROLL(N) nidx[kk] = _mm512_add_epi64(nidx[kk], msep);
		}

		REDUCE(min1)
		{
			midx[kk] = _mm512_mask_mov_epi64(midx[kk], _mm512_cmp_pd_mask(min1[kk], min1[kk + KK], _CMP_GT_OS), midx[kk + KK]);
			min1[kk] = _mm512_min_pd(min1[kk], min1[kk + KK]);
		}

		for (int64 j = 0; j < E512D; ++j)
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
#define N 2
	int64 i = 0;
	val = FLT_MAX;
	uint idx = (uint)-1;

	if (n >= N * E512F)
	{
		__m512 min1[N];
		__m512i midx[N], nidx[N], msep = _mm512_set1_epi32(N * E512F);
		UNROLL(N) min1[kk] = _mm512_set1_ps(val);
		UNROLL(N) midx[kk] = _mm512_set1_epi8((char)0xFF);
		UNROLL(N) nidx[kk] = _mm512_set_epi32(15 + (kk << 4), 14 + (kk << 4), 13 + (kk << 4), 12 + (kk << 4), 11 + (kk << 4), 10 + (kk << 4), 9 + (kk << 4), 8 + (kk << 4), 7 + (kk << 4), 6 + (kk << 4), 5 + (kk << 4), 4 + (kk << 4), 3 + (kk << 4), 2 + (kk << 4), 1 + (kk << 4), 0 + (kk << 4));

		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) 
			{
				midx[kk] = _mm512_mask_mov_epi32(midx[kk], _mm512_cmp_ps_mask(min1[kk], _mm512_loadu_ps(A), _CMP_GT_OS), nidx[kk]);
				min1[kk] = _mm512_min_ps(min1[kk], _mm512_loadu_ps(A));
				A += E512F;
			}

			UNROLL(N) nidx[kk] = _mm512_add_epi64(nidx[kk], msep);
		}

		REDUCE(min1)
		{
			midx[kk] = _mm512_mask_mov_epi32(midx[kk], _mm512_cmp_ps_mask(min1[kk], min1[kk + KK], _CMP_GT_OS), midx[kk + KK]);
			min1[kk] = _mm512_min_ps(min1[kk], min1[kk + KK]);
		}

		for (int64 j = 0; j < E512F; ++j)
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
#define N 4
	int64 i = 0;
	minv = DBL_MAX;
	maxv = -DBL_MAX;

	if (n >= N * E512D)
	{
		__m512d min1[N], max1[N];
		UNROLL(N) min1[kk] = _mm512_set1_pd(minv);
		UNROLL(N) max1[kk] = _mm512_set1_pd(maxv);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N)
			{
				min1[kk] = _mm512_min_pd(min1[kk], _mm512_loadu_pd(A));
				max1[kk] = _mm512_max_pd(max1[kk], _mm512_loadu_pd(A));
				A += E512D;
			}
		}

		REDUCE(min1)
		{
			min1[kk] = _mm512_min_pd(min1[kk], min1[kk + KK]);
			max1[kk] = _mm512_max_pd(max1[kk], max1[kk + KK]);
		}

		minv = _mm512_reduce_min_pd(min1[0]);
		maxv = _mm512_reduce_max_pd(max1[0]);
	}

	for (; i < n; ++i, ++A)
	{
		minv = std::min(minv, *A);
		maxv = std::max(maxv, *A);
	}
}

TARGET512 void GetMinMaxVal512(float* A, int64 n, float& minv, float& maxv)
{
#define N 4
	int64 i = 0;
	minv = FLT_MAX;
	maxv = -FLT_MAX;

	if (n >= N * E512F)
	{
		__m512 min1[N], max1[N];
		UNROLL(N) min1[kk] = _mm512_set1_ps(minv);
		UNROLL(N) max1[kk] = _mm512_set1_ps(maxv);

		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N)
			{
				min1[kk] = _mm512_min_ps(min1[kk], _mm512_loadu_ps(A));
				max1[kk] = _mm512_max_ps(max1[kk], _mm512_loadu_ps(A));
				A += E512F;
			}
		}

		REDUCE(min1)
		{
			min1[kk] = _mm512_min_ps(min1[kk], min1[kk + KK]);
			max1[kk] = _mm512_max_ps(max1[kk], max1[kk + KK]);
		}

		minv = _mm512_reduce_min_ps(min1[0]);
		maxv = _mm512_reduce_max_ps(max1[0]);
	}

	for (; i < n; ++i, ++A)
	{
		minv = std::min(minv, *A);
		maxv = std::max(maxv, *A);
	}
}

TARGET512 double GetMaxVal512(double* A, int64 n)
{
#define N 4
	int64 i = 0;
	double val = -DBL_MAX;

	if (n >= N * E512D)
	{
		__m512d max1[N];
		UNROLL(N) max1[kk] = _mm512_set1_pd(val);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) { max1[kk] = _mm512_max_pd(max1[kk], _mm512_loadu_pd(A)); A += E512D; }
		}

		REDUCE(max1) max1[kk] = _mm512_max_pd(max1[kk], max1[kk + KK]);
		
		val = _mm512_reduce_max_pd(max1[0]);
	}

	for (; i < n; ++i, ++A)
	{
		val = std::max(val, *A);
	}
	return val;
}

TARGET512 float GetMaxVal512(float* A, int64 n)
{
#define N 4
	int64 i = 0;
	float val = -FLT_MAX;

	if (n >= N * E512F)
	{
		__m512 max1[N];
		UNROLL(N) max1[kk] = _mm512_set1_ps(val);

		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) { max1[kk] = _mm512_max_ps(max1[kk], _mm512_loadu_ps(A)); A += E512F; }
		}

		REDUCE(max1) max1[kk] = _mm512_max_ps(max1[kk], max1[kk + KK]);
		
		val = _mm512_reduce_max_ps(max1[0]);
	}

	for (; i < n; ++i, ++A)
	{
		val = std::max(val, *A);
	}

	return val;
}

TARGET512 double GetMaxVal512(double* A, int64 n, int64 sep)
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

	if (n >= N * E512D)
	{
		__m512d max1 = _mm512_set1_pd(val);
		__m512i vindex = _mm512_set_epi64(7 * sep, 6 * sep, 5 * sep, 4 * sep, 3 * sep, 2 * sep, 1 * sep, 0 * sep);
		volatile __m512d t0 = _mm512_i64gather_pd(vindex, A, sizeof(double));

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N)
			{
				volatile __m512d t1 = _mm512_i64gather_pd(vindex, A + E512D * sep, sizeof(double));
				
				max1 = _mm512_max_pd(max1, _mm512_i64gather_pd(vindex, A, sizeof(double)));
				
				A += E512D * sep;
			}
		}

		val = _mm512_reduce_max_pd(max1);
	}

	for (; i < n; ++i, ++A)
	{
		val = std::max(val, *A);
	}

	return val;
}

TARGET512 float GetMaxVal512(float* A, int64 n, int64 sep)
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

	if (n >= N * E512F)
	{
		__m512 max1 = _mm512_set1_ps(val);
		__m512i vindex = _mm512_set_epi32(15 * sep, 14 * sep, 13 * sep, 12 * sep, 11 * sep, 10 * sep, 9 * sep, 8 * sep, 7 * sep, 6 * sep, 5 * sep, 4 * sep, 3 * sep, 2 * sep, 1 * sep, 0 * sep);
		//volatile __m512 t0 = _mm512_i32gather_ps(vindex, A, sizeof(float));

		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N)
			{
				//volatile __m512 t1 = _mm512_i32gather_ps(vindex, A + E512F * sep, sizeof(float));
				
				max1 = _mm512_max_ps(max1, _mm512_i32gather_ps(vindex, A, sizeof(float)));
				
				A += E512F * sep;
			}
		}

		val = _mm512_reduce_max_ps(max1);
	}

	for (; i < n; ++i, ++A)
	{
		val = std::max(val, *A);
	}

	return val;
}

TARGET512 double GetMinVal512(double* A, int64 n)
{
#define N 4
	int64 i = 0;
	double val = DBL_MAX;

	if (n >= N * E512D)
	{
		__m512d min1[N];
		UNROLL(N) min1[kk] = _mm512_set1_pd(val);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) { min1[kk] = _mm512_min_pd(min1[kk], _mm512_loadu_pd(A)); A += E512D; }
		}

		REDUCE(min1) min1[kk] = _mm512_min_pd(min1[kk], min1[kk + KK]);
		
		val = _mm512_reduce_min_pd(min1[0]);
	}

	for (; i < n; ++i, ++A)
	{
		val = std::min(val, *A);
	}
	return val;
}

TARGET512 float GetMinVal512(float* A, int64 n)
{
#define N 4
	int64 i = 0;
	float val = FLT_MAX;

	if (n >= N * E512F)
	{
		__m512 min1[N];
		UNROLL(N) min1[kk] = _mm512_set1_ps(val);

		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) { min1[kk] = _mm512_min_ps(min1[kk], _mm512_loadu_ps(A)); A += E512F; }
		}

		REDUCE(min1) min1[kk] = _mm512_min_ps(min1[kk], min1[kk + KK]);
		
		val = _mm512_reduce_min_ps(min1[0]);
	}

	for (; i < n; ++i, ++A)
	{
		val = std::min(val, *A);
	}

	return val;
}

TARGET512 int64 GetMinVal512(int64* A, int64 n)
{
#define N 4
	int64 i = 0;
	int64 val = 0x7FFFFFFFFFFFFFFF;

	if (n >= N * E512D)
	{
		__m512i min1[N];
		UNROLL(N) min1[kk] = _mm512_set1_epi64(0x7FFFFFFFFFFFFFFF);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) { min1[kk] = _mm512_min_epi64(min1[kk], _mm512_loadu_si512(A)); A += E512D; }
		}

		REDUCE(min1) min1[kk] = _mm512_min_epi64(min1[kk], min1[kk + KK]);
		
		val = _mm512_reduce_min_epi64(min1[0]);
	}

	for (; i < n; ++i, ++A)
	{
		val = std::min(val, *A);
	}
	return val;
}

TARGET512 void SetVal512(uint* A, ushort* B, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E512F)
	{
		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) { _mm512_storeu_si512(A, _mm512_cvtepu16_epi32(_mm256_loadu_si256((__m256i*)B))); A += E512F; B += E512F; }
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
#define N 4
	int64 i = 0;
	int64 slog = 0; double prod = 1;

	if (n >= N * E512D)
	{
		__m512d pd[E512_512];
		UNROLL(E512_512) pd[kk] = _mm512_set1_pd(1.0);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_512)
			{ pd[kk] = _mm512_mul_pd(pd[kk], _mm512_loadu_pd(A)); A += E512D; }

			UNROLL(E512_512) AddExponent512(slog, pd[kk]);
		}

		__m128d* pd1 = (__m128d*)pd;
		UNROLL(E512_128) ChargeLogSSE(slog, prod, pd1[kk]);
	}

	for (; i < n; ++i, ++A)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);
	return prod;
}

TARGET512 double LogProd512(float* A, int64 n)
{
#define N 4
	int64 i = 0;
	int64 slog = 0; double prod = 1;

	if (n >= N * E512D)
	{
		__m512d pd[E512_512];
		UNROLL(E512_512) pd[kk] = _mm512_set1_pd(1.0);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_512) 
			{ pd[kk] = _mm512_mul_pd(pd[kk], _mm512_cvtps_pd(_mm256_loadu_ps(A))); A += E512D; }

			UNROLL(E512_512) AddExponent512(slog, pd[kk]);
		}

		__m128d* pd1 = (__m128d*)pd;
		UNROLL(E512_128) ChargeLogSSE(slog, prod, pd1[kk]);
	}

	for (; i < n; ++i, ++A)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);
	return prod;
}

TARGET512 double LogProd512(double* A, int64 n, int64 sep)
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
		__m512d pd[E512_512];
		UNROLL(E512_512) pd[kk] = _mm512_set1_pd(1.0);
		__m512i vindex = _mm512_mullo_epi64(_mm512_sub_epi64(_mm512_set_epi64(7, 6, 5, 4, 3, 2, 1, 0), _mm512_set1_epi64(0)), _mm512_set1_epi64(sep));

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_512) 
			{ pd[kk] = _mm512_mul_pd(pd[kk], _mm512_i64gather_pd(vindex, A, sizeof(double))); A += E512D * sep; }

			UNROLL(E512_512) AddExponent512(slog, pd[kk]);
		}

		__m128d* pd1 = (__m128d*)pd;
		UNROLL(E512_128) ChargeLogSSE(slog, prod, pd1[kk]);
	}

	for (; i < n; ++i, A += sep)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);

	return prod;
}

TARGET512 double LogProd512(float* A, int64 n, int64 sep)
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

TARGET512 double LogProdDiv512(double* A, double* B, int64 n, int64 sep)
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

TARGET512 double LogProdDiv512(float* A, float* B, int64 n, int64 sep)
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

TARGET512 int64 CountNonZero512(byte* A, int64 n)
{
	int64 i = 0;
	uint64 re = 0;

	if (n >= E512B)
	{
		__m512i a = _mm512_setzero_si512(), z = _mm512_setzero_si512(), o = _mm512_set1_epi8(0x01);

		for (int64 l1 = n - E512B; i <= l1; i += E512B)
		{
			a = _mm512_add_epi64(a, _mm512_sad_epu8(z, _mm512_min_epu8(o, _mm512_loadu_si512(A))));
			A += E512B;
		}
		
		re += _mm512_reduce_add_epi64(a);
	}

	for (; i < n; ++i, ++A)
		if (*A) re++;

	return (int64)re;
}

TARGET512 double Sum512(double* A, int64 n)
{
#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		__m512d s[E512_512] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_512) 
			{ s[kk] = _mm512_add_pd(s[kk], _mm512_loadu_pd(A)); A += E512D; }
		}

		REDUCE(s) s[kk] = _mm512_add_pd(s[kk], s[kk + KK]);

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
#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		__m512d s[E512_512] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_512)
			{ s[kk] = _mm512_add_pd(s[kk], _mm512_cvtps_pd(_mm256_loadu_ps(A))); A += E512D; }
		}

		REDUCE(s) s[kk] = _mm512_add_pd(s[kk], s[kk + KK]);

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
#define N 4
	int64 i = 0;
	volatile float re = 0;

	if (n >= N * E512F)
	{
		__m512 s[E512_512] = { 0 };

		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N)  UNROLL(E512_512) 
			{ s[kk] = _mm512_add_ps(s[kk], _mm512_loadu_ps(A)); A += E512F; }
		}

		REDUCE(s) s[kk] = _mm512_add_ps(s[kk], s[kk + KK]);

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
#define N 4
	int64 re = 0;
	int64 i = 0;

	if (n >= N * E512B)
	{
		__m512i s = _mm512_setzero_si512(), z = _mm512_setzero_si512();

		for (int64 l1 = n - N * E512B; i <= l1; i += N * E512B)
		{
			UNROLL(N) { s = _mm512_add_epi64(s, _mm512_sad_epu8(_mm512_loadu_si512(A), z)); A += E512B; }
		}
		
		re += _mm512_reduce_add_epi64(s);
	}

	for (; i < n; ++i)
		re += *A++;

	return re;
}

TARGET512 double Sum512(double* A, int64 n, int64 sep)
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
		__m512d s[E512_512] = { 0 };
		__m512i vindex = _mm512_set_epi64(7 * sep, 6 * sep, 5 * sep, 4 * sep, 3 * sep, 2 * sep, 1 * sep, 0 * sep);
		UNROLL(E512D) _mm_prefetch((const char*)&A[kk * sep], _MM_HINT_T0);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_512)
			{
				UNROLL(E512D) _mm_prefetch((const char*)&A[E512D * sep + kk * sep], _MM_HINT_T0);
				s[kk] = _mm512_add_pd(s[kk], _mm512_i64gather_pd(vindex, A, sizeof(double)));
				A += E512D * sep;
			}
		}

		REDUCE(s) s[kk] = _mm512_add_pd(s[kk], s[kk + KK]);

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
		__m512d s[E512_512] = { 0 };
		__m256i vindex = _mm256_set_epi32(7 * sep, 6 * sep, 5 * sep, 4 * sep, 3 * sep, 2 * sep, 1 * sep, 0 * sep);
		UNROLL(E512D) _mm_prefetch((const char*)&A[kk * sep], _MM_HINT_T0);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_512)
			{
				UNROLL(E512D) _mm_prefetch((const char*)&A[E512D * sep + kk * sep], _MM_HINT_T0);
				s[kk] = _mm512_add_pd(s[kk], _mm512_cvtps_pd(_mm256_i32gather_ps(A, vindex, sizeof(float))));
				A += E512D * sep;
			}
		}

		REDUCE(s) s[kk] = _mm512_add_pd(s[kk], s[kk + KK]);

		re = __mm512_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, A += sep)
	{
		volatile double v1 = *A;
		re += v1;
	}

	return re;
}

TARGET512 float Sum512x(float* A, int64 n, int64 sep)
{
	//suboptimal to compile
	{
		float re = 0;
		UNROLLHEAD(4)
		for (int64 i = 0; i < n; ++i, A += sep)
			re += *A;
		return re;
	}

	//suboptimal to AVX
#define N 4
	int64 i = 0;
	volatile float re = 0;

	if (n >= N * E512F)
	{
		__m512 s[E512_512] = { 0 };
		__m512i vindex = _mm512_set_epi32(15 * sep, 14 * sep, 13 * sep, 12 * sep, 11 * sep, 10 * sep, 9 * sep, 8 * sep, 7 * sep, 6 * sep, 5 * sep, 4 * sep, 3 * sep, 2 * sep, 1 * sep, 0 * sep);
		UNROLL(E512F) _mm_prefetch((const char*)&A[kk * sep], _MM_HINT_T0);

		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) UNROLL(E512_512)
			{
				UNROLL(E512F) _mm_prefetch((const char*)&A[E512F * sep + kk * sep], _MM_HINT_T0);
				s[kk] = _mm512_add_ps(s[kk], _mm512_i32gather_ps(vindex, A, sizeof(float)));
				A += E512F * sep;
			}
		}

		REDUCE(s) s[kk] = _mm512_add_ps(s[kk], s[kk + KK]);

		re = __mm512_reduce_add_ps(s[0]);
	}

	for (; i < n; ++i, A += sep)
	{
		volatile float v1 = *A;
		re += v1;
	}

	return re;
}

TARGET512 double Prod512(double* A, int64 n)
{
#define N 4
	int64 i = 0;
	volatile double re = 1;

	if (n >= N * E512D)
	{
		__m512d s[E512_512];
		UNROLL(E512_512) s[kk] = _mm512_set1_pd(1);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_512)
			{ s[kk] = _mm512_mul_pd(s[kk], _mm512_loadu_pd(A)); A += E512D; }
		}

		REDUCE(s) s[kk] = _mm512_mul_pd(s[kk], s[kk + KK]);

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
#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		__m512d s[E512_512];
		UNROLL(E512_512) s[kk] = _mm512_set1_pd(1);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_512)
			{ s[kk] = _mm512_mul_pd(s[kk], _mm512_cvtps_pd(_mm256_loadu_ps(A))); A += E512D; }
		}

		REDUCE(s) s[kk] = _mm512_mul_pd(s[kk], s[kk + KK]);

		re = __mm512_reduce_mul_pd(s[0]);
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A++;
		re *= v1;
	}

	return re;
}

TARGET512 float Prod512x(float* A, int64 n)
{
#define N 4
	int64 i = 0;
	volatile float re = 0;

	if (n >= N * E512F)
	{
		__m512 s[E512_512];
		UNROLL(E512_512) s[kk] = _mm512_set1_ps(1);

		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) UNROLL(E512_512)
			{ s[kk] = _mm512_mul_ps(s[kk], _mm512_loadu_ps(A)); A += E512F; }
		}

		REDUCE(s) s[kk] = _mm512_mul_ps(s[kk], s[kk + KK]);

		re = __mm512_reduce_mul_ps(s[0]);
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A++;
		re *= v1;
	}

	return re;
}

TARGET512 double Prod512(double* A, int64 n, int64 sep)
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
		__m512d s[E512_512];
		UNROLL(E512_512) s[kk] = _mm512_set1_pd(1);
		__m512i vindex = _mm512_set_epi64(7 * sep, 6 * sep, 5 * sep, 4 * sep, 3 * sep, 2 * sep, 1 * sep, 0 * sep);
		UNROLL(E512D) _mm_prefetch((const char*)&A[kk * sep], _MM_HINT_T0);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_512)
			{
				UNROLL(E512D) _mm_prefetch((const char*)&A[E512D * sep + kk * sep], _MM_HINT_T0);
				s[kk] = _mm512_mul_pd(s[kk], _mm512_i64gather_pd(vindex, A, sizeof(double)));
				A += E512D * sep;
			}
		}

		REDUCE(s) s[kk] = _mm512_mul_pd(s[kk], s[kk + KK]);

		re = __mm512_reduce_mul_pd(s[0]);
	}

	for (; i < n; ++i, A += sep)
	{
		volatile double v1 = *A;
		re *= v1;
	}

	return re;
}

TARGET512 double Prod512(float* A, int64 n, int64 sep)
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
		__m512d s[E512_512];
		UNROLL(E512_512) s[kk] = _mm512_set1_pd(1);
		__m256i vindex = _mm256_set_epi32(7 * sep, 6 * sep, 5 * sep, 4 * sep, 3 * sep, 2 * sep, 1 * sep, 0 * sep);
		UNROLL(E512D) _mm_prefetch((const char*)&A[kk * sep], _MM_HINT_T0);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_512)
			{
				UNROLL(E512D) _mm_prefetch((const char*)&A[E512D * sep + kk * sep], _MM_HINT_T0);
				s[kk] = _mm512_mul_pd(s[kk], _mm512_cvtps_pd(_mm256_i32gather_ps(A, vindex, sizeof(float))));
				A += E512D * sep;
			}
		}

		REDUCE(s) s[kk] = _mm512_mul_pd(s[kk], s[kk + KK]);

		re = __mm512_reduce_mul_pd(s[0]);
	}

	for (; i < n; ++i, A += sep)
	{
		volatile double v1 = *A;
		re *= v1;
	}

	return re;
}

TARGET512 float Prod512x(float* A, int64 n, int64 sep)
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
		__m512 s[E512_512];
		UNROLL(E512_512) s[kk] = _mm512_set1_ps(1);
		__m512i vindex = _mm512_set_epi32(15 * sep, 14 * sep, 13 * sep, 12 * sep, 11 * sep, 10 * sep, 9 * sep, 8 * sep, 7 * sep, 6 * sep, 5 * sep, 4 * sep, 3 * sep, 2 * sep, 1 * sep, 0 * sep);
		UNROLL(E512F) _mm_prefetch((const char*)&A[kk * sep], _MM_HINT_T0);

		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) UNROLL(E512_512)
			{
				UNROLL(E512F) _mm_prefetch((const char*)&A[E512F * sep + kk * sep], _MM_HINT_T0);
				s[kk] = _mm512_mul_ps(s[kk], _mm512_i32gather_ps(vindex, A, sizeof(float)));
				A += E512F * sep;
			}
		}

		REDUCE(s) s[kk] = _mm512_mul_ps(s[kk], s[kk + KK]);

		re = __mm512_reduce_mul_ps(s[0]);
	}

	for (; i < n; ++i, A += sep)
	{
		volatile float v1 = *A;
		re *= v1;
	}

	return re;
}

TARGET512 double SumSquare512(double* A, int64 n)
{
#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		__m512d s[E512_512] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_512) 
			{ s[kk] = _mm512_fmaddx_pd(_mm512_loadu_pd(A), _mm512_loadu_pd(A), s[kk]); A += E512D; }
		}

		REDUCE(s) s[kk] = _mm512_add_pd(s[kk], s[kk + KK]);

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
#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		__m512d s[E512_512] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_512)
			{
				__m512d v1  = _mm512_cvtps_pd(_mm256_loadu_ps(A)); A += E512D;
				s[kk] = _mm512_fmaddx_pd(v1, v1, s[kk]);
			}
		}

		REDUCE(s) s[kk] = _mm512_add_pd(s[kk], s[kk + KK]);

		re = __mm512_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, ++A)
	{
		volatile double v1 = (double)*A * (double)*A;
		re += v1;
	}

	return re;
}

TARGET512 int64 SumSquare512(byte* A, int64 n)
{
#define N 2
	int64 i = 0;
	uint64 re = 0;

	if (n >= N * E512B)
	{
		__m512i t = _mm512_setzero_si512(), s = _mm512_setzero_si512();
		__m128i* s2 = (__m128i*)&s;

		for (int64 l1 = n - N * E512B; i <= l1; i += N * E512B)
		{
			UNROLL(N) { s = _mm512_add_epi16(s, _mm512_maddubs_epi16(_mm512_loadu_si512(A), _mm512_loadu_si512(A))); A += E512B; }
			
			if ((i & (E512B * 128 - 1)) == 0) [[unlikely]]
			{
				UNROLL(E512_128) t = _mm512_add_epi64(t, _mm512_cvtepi16_epi64(s2[kk]));
				s = _mm512_setzero_si512();
			}
		}
		
		UNROLL(E512_128) t = _mm512_add_epi64(t, _mm512_cvtepi16_epi64(s2[kk]));
		re = _mm512_reduce_add_epi64(t);
	}

	for (; i < n; ++i, ++A)
		re += *A * *A;

	return re;
}

TARGET512 void SumSumSquare512(double* A, int64 n, double& sum, double& sumsq)
{
#define N 4
	int64 i = 0;
	volatile double re1 = 0, re2 = 0;

	if (n >= N * E512D)
	{
		__m512d s1[E512_512] = { 0 }, s2[E512_512] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_512)
			{
				__m512d v1 = _mm512_loadu_pd(A); A += E512D;

				s1[kk] = _mm512_add_pd(v1, s1[kk]);

				s2[kk] = _mm512_fmaddx_pd(v1, v1, s2[kk]);
			}
		}

		REDUCE(s1)
		{
			s1[kk] = _mm512_add_pd(s1[kk], s1[kk + KK]);
			s2[kk] = _mm512_add_pd(s2[kk], s2[kk + KK]);
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
#define N 4
	int64 i = 0;
	volatile double re1 = 0, re2 = 0;

	if (n >= N * E512D)
	{
		__m512d s1[E512_512] = { 0 }, s2[E512_512] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_512)
			{ 
				__m512d v1 = _mm512_cvtps_pd(_mm256_loadu_ps(A)); A += E512D;

				s1[kk] = _mm512_add_pd(v1, s1[kk]);

				s2[kk] = _mm512_fmaddx_pd(v1, v1, s2[kk]);
			}
		}

		REDUCE(s1)
		{
			s1[kk] = _mm512_add_pd(s1[kk], s1[kk + KK]);
			s2[kk] = _mm512_add_pd(s2[kk], s2[kk + KK]);
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
		UNROLL(E512D) { _mm_prefetch((const char*)B, _MM_HINT_T0); B += sep; }

		__m512d s1[E512_512] = { 0 }, s2[E512_512] = { 0 };
		__m512i vindex = _mm512_set_epi64(-9 * sep, -10 * sep, -11 * sep, -12 * sep, -13 * sep, -14 * sep, -15 * sep, -16 * sep);
		
		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_512)
			{
				UNROLL(E512D) { _mm_prefetch((const char*)B, _MM_HINT_T0); B += sep; }

				__m512d b = _mm512_i64gather_pd(vindex, B, sizeof(double));
				
				s1[kk] = _mm512_fmaddx_pd(b, _mm512_loadu_pd(A1), s1[kk]); A1 += E512D;

				s2[kk] = _mm512_fmaddx_pd(b, _mm512_loadu_pd(A2), s2[kk]); A1 += E512D;
			}
		}

		re1 = __mm512_reduce_add_pd(s1[0]);
		re2 = __mm512_reduce_add_pd(s2[0]);

		B -= E512D * sep;
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

TARGET512 double SumProdDiv512(double* A1, float* A2, float* B, int64 sep, int64 n)
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
		UNROLL(E512D) { _mm_prefetch((const char*)B, _MM_HINT_T0); B += sep; }

		__m512d s1[E512_512] = { 0 }, s2[E512_512] = { 0 };
		__m256i vindex = _mm256_set_epi32(-9 * sep, -10 * sep, -11 * sep, -12 * sep, -13 * sep, -14 * sep, -15 * sep, -16 * sep);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_512)
			{
				UNROLL(E512D) { _mm_prefetch((const char*)B, _MM_HINT_T0); B += sep; }

				__m512d b = _mm512_cvtps_pd(_mm256_i32gather_ps(B, vindex, sizeof(float)));

				s1[0] = _mm512_fmaddx_pd(_mm512_loadu_pd(A1), b, s1[kk]); A1 += E512D; 

				s2[0] = _mm512_fmaddx_pd(_mm512_cvtps_pd(_mm256_loadu_ps(A2)), b, s2[kk]); A2 += E512D; 
			}
		}

		re1 = __mm512_reduce_add_pd(s1[0]);
		re2 = __mm512_reduce_add_pd(s2[0]);

		B -= E512D * sep;
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

TARGET512 double SumProdDiv512(float* A1, float* A2, float* B, int64 sep, int64 n)
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
		UNROLL(E512D) { _mm_prefetch((const char*)B, _MM_HINT_T0); B += sep; }

		__m512d s1[E512_512] = { 0 }, s2[E512_512] = { 0 };
		__m256i vindex = _mm256_set_epi32(-9 * sep, -10 * sep, -11 * sep, -12 * sep, -13 * sep, -14 * sep, -15 * sep, -16 * sep);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_512)
			{
				UNROLL(E512D) { _mm_prefetch((const char*)B, _MM_HINT_T0); B += sep; }

				__m512d b = _mm512_cvtps_pd(_mm256_i32gather_ps(B, vindex, sizeof(float)));

				s1[0] = _mm512_fmaddx_pd(_mm512_cvtps_pd(_mm256_loadu_ps(A1)), b, s1[0]); A1 += E512D;

				s2[0] = _mm512_fmaddx_pd(_mm512_cvtps_pd(_mm256_loadu_ps(A2)), b, s2[0]); A2 += E512D;
			}
		}

		re1 = __mm512_reduce_add_pd(s1[0]);
		re2 = __mm512_reduce_add_pd(s2[0]);

		B -= E512D * sep;
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

TARGET512 float SumProdDiv512x(float* A1, float* A2, float* B, int64 sep, int64 n)
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

TARGET512 double SumProd512(double* A, double* B, int64 sep, int64 n)
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
		__m512d s = _mm512_setzero_pd();
		__m512i vindex = _mm512_set_epi64(7 * sep, 6 * sep, 5 * sep, 4 * sep, 3 * sep, 2 * sep, 1 * sep, 0 * sep);
		volatile __m512d t0 = _mm512_i64gather_pd(vindex, B, sizeof(double));

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N)
			{
				volatile __m512d t1 = _mm512_i64gather_pd(vindex, B + E512D * sep, sizeof(double));
				
				s = _mm512_fmaddx_pd(_mm512_loadu_pd(A), _mm512_i64gather_pd(vindex, B, sizeof(double)), s); A += E512D;
				
				B += E512D * sep;
			}
		}

		re = __mm512_reduce_add_pd(s); 
	}

	for (; i < n; ++i, ++A, B += sep)
	{
		volatile double v1 = *A * *B;
		re += v1;
	}

	return re;
}

TARGET512 double SumProd512(float* A, float* B, int64 sep, int64 n)
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
		UNROLL(E256F) { _mm_prefetch((const char*)B, _MM_HINT_T0); B += sep; }

		__m512d s[E512_512] = { 0 };
		__m256i vindex = _mm256_set_epi32(-9 * sep, -10 * sep, -11 * sep, -12 * sep, -13 * sep, -14 * sep, -15 * sep, -16 * sep);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_512)
			{
				UNROLL(E256F) { _mm_prefetch((const char*)B, _MM_HINT_T0); B += sep; }
				
				s[kk] = _mm512_fmaddx_pd(_mm512_cvtps_pd(_mm256_loadu_ps(A)), _mm512_cvtps_pd(_mm256_i32gather_ps(B, vindex, sizeof(float))), s[kk]);
				
				A += E512D;
			}
		}

		REDUCE(s) s[kk] = _mm512_add_pd(s[kk], s[kk + KK]);

		re = __mm512_reduce_add_pd(s[0]);

		B -= E256F * sep;
	}

	for (; i < n; ++i, A++, B += sep)
	{
		volatile double v1 = (double)*A * (double)*B;
		re += v1;
	}

	return re;
}

TARGET512 float SumProd512x(float* A, float* B, int64 sep, int64 n)
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

TARGET512 double SumProd512(double* A, double* B, int64 n)
{
#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		__m512d s[E512_512] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N)UNROLL(E512_512) 
			{ s[kk] = _mm512_fmaddx_pd(_mm512_loadu_pd(A), _mm512_loadu_pd(B), s[kk]); A += E512D; B += E512D; }
		}

		REDUCE(s) s[kk] = _mm512_add_pd(s[kk], s[kk + KK]);

		re = __mm512_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, ++A, ++B)
	{
		volatile double v1 = (double)*A * (double)*B;
		re += v1;
	}

	return re;
}

TARGET512 double SumProd512(float* A, float* B, int64 n)
{
#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		__m512d s[E512_512] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_512) 
			{ s[kk] = _mm512_fmaddx_pd(_mm512_cvtps_pd(_mm256_loadu_ps(A)), _mm512_cvtps_pd(_mm256_loadu_ps(B)), s[kk]); A += E512D; B += E512D; }
		}

		REDUCE(s) s[kk] = _mm512_add_pd(s[kk], s[kk + KK]);

		re = __mm512_reduce_add_pd(s[0]);
	}
	
	for (; i < n; ++i, ++A, ++B)
	{
		volatile double v1 = (double)*A * (double)*B;
		re += v1;
	}

	return re;
}

TARGET512 float SumProd512x(float* A, float* B, int64 n)
{
#define N 4
	int64 i = 0;
	volatile float re = 0;

	if (n >= N * E512F)
	{
		__m512 s[E512_512] = { 0 };

		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) UNROLL(E512_512) 
			{ s[kk] = _mm512_fmaddx_ps(_mm512_loadu_ps(A), _mm512_loadu_ps(B), s[kk]); A += E512F; B += E512F; }
		}

		REDUCE(s) s[kk] = _mm512_add_ps(s[kk], s[kk + KK]);

		re = __mm512_reduce_add_ps(s[0]);
	}

	for (; i < n; ++i, ++A, ++B)
	{
		volatile float v1 = (float)*A * (float)*B;
		re += v1;
	}

	return re;
}

TARGET512 double SumProd512(double* A, double* B, double* C, int64 n)
{
#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		__m512d s[E512_512] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_512)
			{ s[kk] = _mm512_fmaddx_pd(_mm512_mul_pd(_mm512_loadu_pd(A), _mm512_loadu_pd(B)), _mm512_loadu_pd(C), s[kk]); A += E512D; B += E512D; C += E512D; }
		}

		REDUCE(s) s[kk] = _mm512_add_pd(s[kk], s[kk + KK]);

		re = __mm512_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, ++A, ++B, ++C)
	{
		volatile double v1 = (double)*A * (double)*B * (double)*C;
		re += v1;
	}

	return re;
}

TARGET512 float SumProd512(float* A, float* B, float* C, int64 n)
{
#define N 4
	int64 i = 0;
	volatile float re = 0;

	if (n >= N * E512F)
	{
		__m512 s[E512_512] = { 0 };

		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) UNROLL(E512_512) 
			{ s[kk] = _mm512_fmaddx_ps(_mm512_mul_ps(_mm512_loadu_ps(A), _mm512_loadu_ps(B)), _mm512_loadu_ps(C), s[kk]); A += E512F; B += E512F; C += E512F; }
		}

		REDUCE(s) s[kk] = _mm512_add_ps(s[kk], s[kk + KK]);

		re = __mm512_reduce_add_ps(s[0]);
	}

	for (; i < n; ++i, ++A, ++B, ++C)
	{
		volatile float v1 = (float)*A * (float)*B * (float)*C;
		re += v1;
	}

	return re;
}

TARGET512 double SumSqProd512(double* A, double* B, int64 n)
{
#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		__m512d s[E512_512] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_512) 
			{ s[kk] = _mm512_fmaddx_pd(_mm512_mul_pd(_mm512_loadu_pd(A), _mm512_loadu_pd(A)), _mm512_loadu_pd(B), s[kk]); A += E512D; B += E512D; }
		}

		REDUCE(s) s[kk] = _mm512_add_pd(s[kk], s[kk + KK]);

		re = __mm512_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, ++A, ++B)
	{
		volatile double v1 = (double)*A * (double)*A * (double)*B;
		re += v1;
	}

	return re;
}

TARGET512 float SumSqProd512(float* A, float* B, int64 n)
{
#define N 4
	int64 i = 0;
	volatile float re = 0;

	if (n >= N * E512F)
	{
		__m512 s[E512_512] = { 0 };

		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) UNROLL(E512_512)
			{ s[kk] = _mm512_fmaddx_ps(_mm512_mul_ps(_mm512_loadu_ps(A), _mm512_loadu_ps(A)), _mm512_loadu_ps(B), s[kk]); A += E512F; B += E512F; }
		}

		REDUCE(s) s[kk] = _mm512_add_ps(s[kk], s[kk + KK]);

		re = __mm512_reduce_add_ps(s[0]);
	}

	for (; i < n; ++i, ++A, ++B)
	{
		volatile float v1 = (float)*A * (float)*A * (float)*B;
		re += v1;
	}

	return re;
}

TARGET512 void Add512(double* A, double* B, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E512D)
	{
		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) { _mm512_storeu_pd(A, _mm512_add_pd(_mm512_loadu_pd(A), _mm512_loadu_pd(B))); A += E512D; B += E512D; }
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGET512 void Add512(float* A, float* B, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E512F)
	{
		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) { _mm512_storeu_ps(A, _mm512_add_ps(_mm512_loadu_ps(A), _mm512_loadu_ps(B))); A += E512F; B += E512F; }
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGET512 void Add512(int64* A, int64* B, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E512D)
	{
		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) { _mm512_storeu_si512(A, _mm512_add_epi64(_mm512_loadu_si512(A), _mm512_loadu_si512(B))); A += E512D; B += E512D; }
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGET512 void Add512(int* A, int* B, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E512F)
	{
		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) { _mm512_storeu_si512(A, _mm512_add_epi32(_mm512_loadu_si512(A), _mm512_loadu_si512(B))); A += E512F; B += E512F; }
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGET512 void Add512(int* A, int B, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E512F)
	{
		__m512i b = _mm512_set1_epi32(B);

		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) { _mm512_storeu_si512(A, _mm512_add_epi32(_mm512_loadu_si512(A), b)); A += E512F; }
		}
	}

	for (; i < n; ++i, A++)
		*A += B;
}

TARGET512 void Add512(double* A, double B, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E512D)
	{
		__m512d b = _mm512_set1_pd(B);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) { _mm512_storeu_pd(A, _mm512_add_pd(_mm512_loadu_pd(A), b)); A += E512D; }
		}
	}

	for (; i < n; ++i, A++)
		*A += B;
}

TARGET512 void Add512(float* A, float B, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E512F)
	{
		__m512 b = _mm512_set1_ps(B);

		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) { _mm512_storeu_ps(A, _mm512_add_ps(_mm512_loadu_ps(A), b)); A += E512F; }
		}
	}

	for (; i < n; ++i, A++)
		*A += B;
}

TARGET512 void Mul512(double* A, double* B, double* C, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E512D)
	{
		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) { _mm512_storeu_pd(A, _mm512_mul_pd(_mm512_loadu_pd(B), _mm512_loadu_pd(C))); A += E512D; B += E512D; C += E512D; }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *B++ * *C++;
		*A++ = v1;
	}
}

TARGET512 void Mul512(float* A, float* B, float* C, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E512F)
	{
		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) { _mm512_storeu_ps(A, _mm512_mul_ps(_mm512_loadu_ps(B), _mm512_loadu_ps(C))); A += E512F; B += E512F; C += E512F; }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *B++ * *C++;
		*A++ = v1;
	}
}

TARGET512 void Mul512(double* A, double* B, double C, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E512D)
	{
		__m512d c = _mm512_set1_pd(C);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) { _mm512_storeu_pd(A, _mm512_mul_pd(_mm512_loadu_pd(B), c)); A += E512D; B += E512D; }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *B++ * C;
		*A++ = v1;
	}
}

TARGET512 void Mul512(float* A, float* B, float C, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E512F)
	{
		__m512 c = _mm512_set1_ps(C);

		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) { _mm512_storeu_ps(A, _mm512_mul_ps(_mm512_loadu_ps(B), c)); A += E512F; B += E512F; }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *B++ * C;
		*A++ = v1;
	}
}

TARGET512 void Mul512(double* A, double B, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E512D)
	{
		__m512d b = _mm512_set1_pd(B);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) { _mm512_storeu_pd(A, _mm512_mul_pd(_mm512_loadu_pd(A), b)); A += E512D; }
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
#define N 4
	int64 i = 0;

	if (n >= N * E512F)
	{
		__m512 b = _mm512_set1_ps(B);

		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) { _mm512_storeu_ps(A, _mm512_mul_ps(_mm512_loadu_ps(A), b)); A += E512F; }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A * B;
		*A++ = v1;
	}
}

TARGET512 void Div512(double* A, double B, double* C, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E512D)
	{
		__m512d b = _mm512_set1_pd(B);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) { _mm512_storeu_pd(A, _mm512_div_pd(b, _mm512_loadu_pd(C))); A += E512D; C += E512D; }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = B / *C++;
		*A++ = v1;
	}
}

TARGET512 void Div512(float* A, float B, float* C, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E512F)
	{
		__m512 b = _mm512_set1_ps(B);

		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) { _mm512_storeu_ps(A, _mm512_div_ps(b, _mm512_loadu_ps(C))); A += E512F; C += E512F; }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = B / *C++;
		*A++ = v1;
	}
}

TARGET512 void Div512(double* A, double* B, double* C, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E512D)
	{
		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) { _mm512_storeu_pd(A, _mm512_div_pd(_mm512_loadu_pd(B), _mm512_loadu_pd(C))); A += E512D; B += E512D; C += E512D; }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *B++ / *C++;
		*A++ = v1;
	}
}

TARGET512 void Div512(float* A, float* B, float* C, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E512F)
	{
		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) { _mm512_storeu_ps(A, _mm512_div_ps(_mm512_loadu_ps(B), _mm512_loadu_ps(C))); A += E512F; B += E512F; C += E512F; }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *B++ / *C++;
		*A++ = v1;
	}
}

TARGET512 void AddProd512(double* A, double* B, double* C, int64 n)
{
#define N 1
	int64 i = 0;

	if (n >= N * E512D)
	{
		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) { _mm512_storeu_pd(A, _mm512_fmaddx_pd(_mm512_loadu_pd(B), _mm512_loadu_pd(C), _mm512_loadu_pd(A))); A += E512D; B += E512D; C += E512D; }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *B++ * *C++;
		*A++ += v1;
	}
}

TARGET512 void AddProd512(float* A, float* B, float* C, int64 n)
{
#define N 1
	int64 i = 0;

	if (n >= N * E512F)
	{
		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) { _mm512_storeu_ps(A, _mm512_fmaddx_ps(_mm512_loadu_ps(B), _mm512_loadu_ps(C), _mm512_loadu_ps(A))); A += E512F; B += E512F; C += E512F; }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *B++ * *C++;
		*A++ += v1;
	}
}

TARGET512 void AddProd512(double* A, double* B, double C, int64 n)
{
#define N 1
	int64 i = 0;

	if (n >= N * E512D)
	{
		__m512d c = _mm512_set1_pd(C);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) { _mm512_storeu_pd(A, _mm512_fmaddx_pd(_mm512_loadu_pd(B), c, _mm512_loadu_pd(A))); A += E512D; B += E512D; }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *B++ * C;
		*A++ += v1;
	}
}

TARGET512 void AddProd512(double* A, float* B, double C, int64 n)
{
#define N 1
	int64 i = 0;

	if (n >= N * E512D)
	{
		__m512d c = _mm512_set1_pd(C);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) { _mm512_storeu_pd(A, _mm512_fmaddx_pd(_mm512_cvtps_pd(_mm256_loadu_ps(B)), c, _mm512_loadu_pd(A))); A += E512D; B += E512D; }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *B++ * C;
		*A++ += v1;
	}
}

TARGET512 void AddProd512(float* A, float* B, float C, int64 n)
{
#define N 1
	int64 i = 0;

	if (n >= N * E512F)
	{
		__m512 c = _mm512_set1_ps(C);

		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) { _mm512_storeu_ps(A, _mm512_fmaddx_ps(_mm512_loadu_ps(B), c, _mm512_loadu_ps(A))); A += E512F; B += E512F; }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *B++ * C;
		*A++ += v1;
	}
}

TARGET512 void Unify512(double* A, int64 n)
{
#define N 4
	int64 i = 0;
	double invsum = 1.0 / (Sum512(A, n) + n * MIN_FREQ);

	if (n >= E512D)
	{
		__m512d b = _mm512_set1_pd(invsum), c = _mm512_set1_pd(MIN_FREQ * invsum);

		UNROLLHEAD(N)
		for (int64 l1 = n - E512D; i <= l1; i += E512D)
		{
			_mm512_storeu_pd(A, _mm512_fmaddx_pd(_mm512_loadu_pd(A), b, c)); 
			A += E512D;
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
#define N 4
	int64 i = 0;
	double invsum = 1.0 / (Sum512(A, n) + n * MIN_FREQ);

	if (n >= E512D)
	{
		__m512d b = _mm512_set1_pd(invsum), c = _mm512_set1_pd(MIN_FREQ * invsum);

		UNROLLHEAD(N)
		for (int64 l1 = n - E512D; i <= l1; i += E512D)
		{
			_mm256_storeu_ps(A, _mm512_cvtpd_ps(_mm512_fmaddx_pd(_mm512_cvtps_pd(_mm256_loadu_ps(A)), b, c))); 
			A += E512D;
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

	if (n >= E512B)
	{
		__m512i v = _mm512_set1_epi8(val);

		UNROLL(N)
		for (int64 l1 = n - E512B; i <= l1; i += E512B)
		{
			char* Ab = A;
			__mmask64 mask = _mm512_cmpeq_epu8_mask(_mm512_loadu_si512(A), v); A += E512B;

			int64 count = (int64)_mm_popcnt_u64(mask);

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

TARGET512 int64 CountChar512(char* A, char val, int64 n)
{
#define N 4
	uint64 re = 0;
	int64 i = 0;

	if (n >= E512B)
	{
		__m512i v = _mm512_set1_epi8(val);
		
		UNROLLHEAD(N)
		for (int64 l1 = n - E512B; i <= l1; i += E512B)
		{
			re += _mm_popcnt_u64(_mm512_cmpeq_epu8_mask(_mm512_loadu_si512(A), v));
			A += E512B;
		}
	}

	for (; i < n; ++i, A++)
		if (*A == val) re++;

	return (int64)re;
}

TARGET512 void DiagQuadForm512x(double* res1, double* A, double* D, int64 m, int64 n)
{
#define NA 2
#define NA8 16
#define NB8 8

    for (int64 i = 0; i < m; i += NA8)
	{
        for (int64 j = 0; j + NB8 <= i + NA8; j += NB8)
		{
			double* pA1 = A + i, *pA2 = A + j, *pD = D; 
			__m512d r[NB8][NA] = { 0 };

            for (int64 k = 0; k < n; k++, pA1 += m, pA2 += m, pD ++) 
			{
				 _mm_prefetch((const char*)(pA1 + m), _MM_HINT_T0);
				 _mm_prefetch((const char*)(pA2 + m), _MM_HINT_T0);

				__m512d d = _mm512_set1_pd(*pD);
				__m512d ad[NA];

				UNROLL(NA) ad[kk] = _mm512_mul_pd(_mm512_loadu_pd(pA1 + kk * 8), d);

				UNROLL(NB8) 
				{
					__m512d a2 = _mm512_set1_pd(pA2[kk]);
					UNROLLHEAD(NA) for (int ii = 0; ii < NA; ++ii)
						r[kk][ii] = _mm512_fmadd_pd(a2, ad[ii], r[kk][ii]);
				}
			}

			for (int jj = 0; jj < NB8; ++jj)
				for (int ii = 0; ii < NA8; ++ii)
					res1[i + ii + (j + jj) * m] = res1[j + jj + (i + ii) * m] = *((double*)&r[jj][0] + ii);
        }
    }
}

TARGET512 void DiagQuadForm512x(float* res1, float* A, float* D, int64 m, int64 n)
{
#define NA 1
#define NB 1
#define NA16 16
#define NB16 16

    for (int64 i = 0; i < m; i += NA16)
	{
        for (int64 j = 0; j + NB16 <= i + NA16; j += NB16)
		{
			float* pA1 = A + i, *pA2 = A + j, *pD = D; 
			__m512 r[NB16][NA] = { 0 };

            for (int64 k = 0; k < n; k++, pA1 += m, pA2 += m, pD ++) 
			{
				 _mm_prefetch((const char*)(pA1 + m), _MM_HINT_T0);
				 _mm_prefetch((const char*)(pA2 + m), _MM_HINT_T0);

				__m512 d = _mm512_set1_ps(*pD);
				__m512 ad[NA];

				UNROLL(NA) ad[kk] = _mm512_mul_ps(_mm512_loadu_ps(pA1 + kk * E512F), d);

				UNROLL(NB16) 
				{
					__m512 a2 = _mm512_set1_ps(pA2[kk]);
					UNROLLHEAD(NA) for (int ii = 0; ii < NA; ++ii)
						r[kk][ii] = _mm512_fmadd_ps(a2, ad[ii], r[kk][ii]);
				}
			}

			for (int jj = 0; jj < NB16; ++jj)
				for (int ii = 0; ii < NA16; ++ii)
					res1[i + ii + (j + jj) * m] = res1[j + jj + (i + ii) * m] = *((float*)&r[jj][0] + ii);
        }
    }
}

TARGET512 void DiagQuadForm512(double* res1, double* A, double* D, int64 m, int64 n)
{
#define N 4
#define DECLARE			double* pA1 = (A + n * i), *pA2 = (A + n * j), *pD = D; __m512d a1[N], a2[N], r[N][N] = { 0 }
#define ALOAD1(ii)		a1[ii] = _mm512_loadu_pd(pA1 + n * ii)
#define ALOAD2(ii)		a2[ii] = _mm512_loadu_pd(pA2 + n * ii)
#define ALOAD3(ii)		a1[ii] = a2[ii] = _mm512_loadu_pd(pA1 + n * ii)
#define ADMUL(ii)		a1[ii] = _mm512_mul_pd(a1[ii], _mm512_loadu_pd(pD))
#define FMADD(ii,jj)	r[ii][jj] = _mm512_fmaddx_pd(a1[ii], a2[jj], r[ii][jj])
#define RDUADD(ii,jj)	res1[(i+ii) * m + (j+jj)] = _mm512_reduce_add_pd(r[ii][jj])
#define REMAIN(ii,jj)	res1[(i+ii) * m + (j+jj)] += A[(i+ii) * n + k] * A[(j+jj) * n + k] * D[k]
#define FINAL(ii,jj)	res1[(i+ii) + (j+jj) * m] = res1[(i+ii) * m + (j+jj)]
	
    int64 i = 0, j = 0, k = 0;
    for (i = 0; i + N <= m; i += N)
	{
        for (j = 0; j < i; j += N)
		{
			DECLARE;

            for (k = 0; k + E512D <= n; k += E512D, pA1 += E512D, pA2 += E512D, pD += E512D) 
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

            for (k = 0; k + E512D <= n; k += E512D, pA1 += E512D, pD += E512D) 
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

            for (k = 0; k + E512D <= n; k += E512D, pA1 += E512D, pA2 += E512D, pD += E512D) 
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

            for (k = 0; k + E512D <= n; k += E512D, pA1 += E512D, pD += E512D) 
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

TARGET512 void DiagQuadForm512(double* res2, double* A, double* D, double* B, int64 m, int64 n)
{
#define N 8
#define DECLARE			double* pA = (A + n * i), *pB = B, *pD = D; __m512d bd, r[N] = { 0 }
#define BDMUL			bd = _mm512_mul_pd(_mm512_loadu_pd(pB), _mm512_loadu_pd(pD)); pB += E512D; pD += E512D
#define FMADD(ii)		r[ii] = _mm512_fmaddx_pd(_mm512_loadu_pd(pA + n * ii), bd, r[ii]); 
#define RDUADD(ii)		res2[(i+ii)] = _mm512_reduce_add_pd(r[ii])
#define FINAL(ii)		res2[(i+ii)] += A[(i+ii) * n + k] * B[k] * D[k]

    int64 i = 0, k = 0;
    for (i = 0; i + N <= m; i += N)
	{
		DECLARE;

        for (k = 0; k + E512D <= n; k += E512D, pA += E512D) 
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

        for (k = 0; k + E512D <= n; k += E512D, pA += E512D) 
		{ BDMUL; LOOPNa(FMADD); }
			
		LOOPNa(RDUADD);
			
		VECTORIZE
        for (; k < n; k++) 
		{ LOOPNa(FINAL); }
	}

#undef DECLARE
#undef BDMUL
#undef FMADD
#undef REDUCE2
#undef FINAL
}

TARGET512 void DiagQuadForm512(double* res3, double* B, double* D, int64 n)
{
#define N 4
	__m512d s[N] = { 0 };

	int64 k = 0;
    for (; k + N * E512D <= n; k += N * E512D) 
	{
		UNROLL(N)
		{
			s[kk] = _mm512_fmaddx_pd(_mm512_mul_pd(_mm512_loadu_pd(B), _mm512_loadu_pd(B)), _mm512_loadu_pd(D), s[kk]); 
			B += E512D;
			D += E512D;
		}
	}

	REDUCE(s) s[kk] = _mm512_add_pd(s[kk], s[kk + KK]);
			
	res3[0] = _mm512_reduce_add_pd(s[0]);
			
	VECTORIZE
    for (; k < n; k++, B++, D++) 
		res3[0] += B[0] * B[0] * D[0];
}

TARGET512 void MatrixMul512(double* res, double* A, double* B, int64 m, int64 n, int64 p)
{
#define N 4
#define DECLARE			double* pA = (A + n * i), *pB = (B + n * j); __m512d a[N], b[N], r[N][N] = { 0 }
#define ALOAD(kk)		a[kk] = _mm512_loadu_pd(pA + n * kk)
#define BLOAD(kk)		b[kk] = _mm512_loadu_pd(pB + n * kk)
#define FMADD(ii,jj)	r[ii][jj] = _mm512_fmaddx_pd(a[ii], b[jj], r[ii][jj])
#define RDUADD(ii,jj)	res[(i+ii) + (j+jj) * m] = _mm512_reduce_add_pd(r[ii][jj])
#define REMAIN(ii,jj)	res[(i+ii) + (j+jj) * m] += A[(i+ii) * n + k] * B[(j+jj) * n + k]

    int64 i = 0, j = 0, k = 0;
    for (i = 0; i + N <= m; i += N)
	{
        for (j = 0; j + N <= p; j += N)
		{
			DECLARE;

            for (k = 0; k + E512D <= n; k += E512D, pA += E512D, pB += E512D) 
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

            for (k = 0; k + E512D <= n; k += E512D, pA += E512D, pB += E512D) 
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

            for (k = 0; k + E512D <= n; k += E512D, pA += E512D, pB += E512D) 
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

            for (k = 0; k + E512D <= n; k += E512D, pA += E512D, pB += E512D) 
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

TARGET512 void DiagQuadForm512(float* res1, float* A, float* D, int64 m, int64 n)
{
#define N 4
#define DECLARE			float* pA1 = (A + n * i), *pA2 = (A + n * j), *pD = D; __m512 a1[N], a2[N], r[N][N] = { 0 }
#define ALOAD1(ii)		a1[ii] = _mm512_loadu_ps(pA1 + n * ii)
#define ALOAD2(ii)		a2[ii] = _mm512_loadu_ps(pA2 + n * ii)
#define ALOAD3(ii)		a1[ii] = a2[ii] = _mm512_loadu_ps(pA1 + n * ii)
#define ADMUL(ii)		a1[ii] = _mm512_mul_ps(a1[ii], _mm512_loadu_ps(pD))
#define FMADD(ii,jj)	r[ii][jj] = _mm512_fmaddx_ps(a1[ii], a2[jj], r[ii][jj])
#define RDUADD(ii,jj)	res1[(i+ii) * m + (j+jj)] = _mm512_reduce_add_ps(r[ii][jj])
#define REMAIN(ii,jj)	res1[(i+ii) * m + (j+jj)] += A[(i+ii) * n + k] * A[(j+jj) * n + k] * D[k]
#define FINAL(ii,jj)	res1[(i+ii) + (j+jj) * m] = res1[(i+ii) * m + (j+jj)]
	
    int64 i = 0, j = 0, k = 0;
    for (i = 0; i + N <= m; i += N)
	{
        for (j = 0; j < i; j += N)
		{
			DECLARE;

            for (k = 0; k + E512F <= n; k += E512F, pA1 += E512F, pA2 += E512F, pD += E512F) 
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

            for (k = 0; k + E512F <= n; k += E512F, pA1 += E512F, pD += E512F) 
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

            for (k = 0; k + E512F <= n; k += E512F, pA1 += E512F, pA2 += E512F, pD += E512F) 
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

            for (k = 0; k + E512F <= n; k += E512F, pA1 += E512F, pD += E512F) 
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

TARGET512 void DiagQuadForm512(float* res2, float* A, float* D, float* B, int64 m, int64 n)
{
#define N 8
#define DECLARE			float* pA = (A + n * i), *pB = B, *pD = D; __m512 bd, r[N] = { 0 }
#define BDMUL			bd = _mm512_mul_ps(_mm512_loadu_ps(pB), _mm512_loadu_ps(pD)); pB += E512F; pD += E512F
#define FMADD(ii)		r[ii] = _mm512_fmaddx_ps(_mm512_loadu_ps(pA + n * ii), bd, r[ii]); 
#define RDUADD(ii)		res2[(i+ii)] = _mm512_reduce_add_ps(r[ii])
#define FINAL(ii)		res2[(i+ii)] += A[(i+ii) * n + k] * B[k] * D[k]

    int64 i = 0, k = 0;
    for (i = 0; i + N <= m; i += N)
	{
		DECLARE;

        for (k = 0; k + E512F <= n; k += E512F, pA += E512F) 
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

        for (k = 0; k + E512F <= n; k += E512F, pA += E512F) 
		{ BDMUL; LOOPNa(FMADD); }
			
		LOOPNa(RDUADD);
			
		VECTORIZE
        for (; k < n; k++) 
		{ LOOPNa(FINAL); }
	}

#undef DECLARE
#undef BDMUL
#undef FMADD
#undef REDUCE2
#undef FINAL
}

TARGET512 void DiagQuadForm512(float* res3, float* B, float* D, int64 n)
{
#define N 4
	__m512 s[N] = { 0 };

	int64 k = 0;
    for (; k + N * E512F <= n; k += N * E512F) 
	{
		UNROLL(N)
		{
			s[kk] = _mm512_fmaddx_ps(_mm512_mul_ps(_mm512_loadu_ps(B), _mm512_loadu_ps(B)), _mm512_loadu_ps(D), s[kk]); 
			B += E512F; 
			D += E512F; 
		}
	}

	REDUCE(s) s[kk] = _mm512_add_ps(s[kk], s[kk + KK]);
			
	res3[0] = _mm512_reduce_add_ps(s[0]);
			
	VECTORIZE
    for (; k < n; k++, B++, D++) 
		res3[0] += B[0] * B[0] * D[0];
}

TARGET512 void MatrixMul512(float* res, float* A, float* B, int64 m, int64 n, int64 p)
{
#define N 4
#define DECLARE			float* pA = (A + n * i), *pB = (B + n * j); __m512 a[N], b[N], r[N][N] = { 0 }
#define ALOAD(kk)		a[kk] = _mm512_loadu_ps(pA + n * kk)
#define BLOAD(kk)		b[kk] = _mm512_loadu_ps(pB + n * kk)
#define FMADD(ii,jj)	r[ii][jj] = _mm512_fmaddx_ps(a[ii], b[jj], r[ii][jj])
#define RDUADD(ii,jj)	res[(i+ii) + (j+jj) * m] = _mm512_reduce_add_ps(r[ii][jj])
#define REMAIN(ii,jj)	res[(i+ii) + (j+jj) * m] += A[(i+ii) * n + k] * B[(j+jj) * n + k]

    int64 i = 0, j = 0, k = 0;
    for (i = 0; i + N <= m; i += N)
	{
        for (j = 0; j + N <= p; j += N)
		{
			DECLARE;

            for (k = 0; k + E512F <= n; k += E512F, pA += E512F, pB += E512F) 
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

            for (k = 0; k + E512F <= n; k += E512F, pA += E512F, pB += E512F) 
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

            for (k = 0; k + E512F <= n; k += E512F, pA += E512F, pB += E512F) 
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

            for (k = 0; k + E512F <= n; k += E512F, pA += E512F, pB += E512F) 
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
