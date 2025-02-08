/* NEON Instruction Set Functions */

#include "vcfpop.h"

#ifdef __aarch64__

template struct RNGNEO<double>;
template struct RNGNEO<float >;
template TARGETNEO void RNGNEO<double>::Integer<uint  >(uint  * re, int64 n, uint   minv, uint   maxv);
template TARGETNEO void RNGNEO<double>::Integer<uint64>(uint64* re, int64 n, uint64 minv, uint64 maxv);
template TARGETNEO void RNGNEO<float >::Integer<uint  >(uint  * re, int64 n, uint   minv, uint   maxv);
template TARGETNEO void RNGNEO<float >::Integer<uint64>(uint64* re, int64 n, uint64 minv, uint64 maxv);

#ifndef _RNGNEO_FP64
#define XS 32

/* Initialize rng */
TARGETNEO RNGNEO<double>::RNGNEO()
{

}

/* Initialize rng */
TARGETNEO RNGNEO<double>::RNGNEO(uint64 seed, uint64 salt)
{
	uint64x2_t a[XS], s, m;

	UNROLL(XS) { a[kk] = vld1q_u64(((uint64[]) { seed, seed + 1 })); seed += 2; }

	s = vdupq_n_u64(salt);
	m = vdupq_n_u32(0x5bd1e995);

	UNROLL(XS) a[kk] = veorq_u64(a[kk], vshlq_n_u64(vmvnq_u32(a[kk]), 32));

	s             = veorq_u64(s    , vshlq_n_u64(vmvnq_u32(s    ), 32));

	// uint s = s ^ 4;
	s = veorq_u32(s, vdupq_n_u32(4));

	// a *= m;
	UNROLL(XS) a[kk] = vmulq_u32(a[kk], m);

	// a ^= a >> 24;
	UNROLL(XS) a[kk] = veorq_u64(a[kk], vshrq_n_u32(a[kk], 24));

	// a *= m;
	UNROLL(XS) a[kk] = vmulq_u32(a[kk], m);

	// s *= m;
	s = vmulq_u32(s, m);

	// a ^= s;
	UNROLL(XS) a[kk] = veorq_u64(a[kk], s);

	// a ^= a >> 13;
	UNROLL(XS) a[kk] = veorq_u64(a[kk], vshrq_n_u32(a[kk], 13));

	// a *= m;
	UNROLL(XS) a[kk] = vmulq_u32(a[kk], m);

	// a ^= a >> 15;
	UNROLL(XS) a[kk] = veorq_u64(a[kk], vshrq_n_u32(a[kk], 15));

	// original
	UNROLL(XS) x[kk] = veorq_u64(vdupq_n_u64(0x159A55E5075BCD15), a[kk]);

	UNROLL(XS) a[kk] = vshlq_n_u64(a[kk], 6);

	UNROLL(XS) y[kk] = veorq_u64(vdupq_n_u64(0x054913331F123BB5), a[kk]);
}

/* Draw 64 64-bit integers in [0,n), 64*n frequencies are in arr */
TARGETNEO void RNGNEO<double>::Poly(float64x2_t* arr, int n, void* re)
{
	float64x2_t t[XS], s[XS];
	float64x2_t one = vdupq_n_f64(1.0);
	uint64x2_t mask1 = vdupq_n_u64(0x000FFFFFFFFFFFFF);
	uint64x2_t mask2 = vdupq_n_u64(0x3FF0000000000000);
	uint64x2_t* r = (uint64x2_t*)t; uint64x2_t* re64 = (uint64x2_t*)re; uint32x4_t* re32 = (uint32x4_t*)re;

	UNROLL(XS) s[kk] = vdupq_n_f64(0);

	for (int i = 0; i < n * XS; i += XS)
		UNROLL(XS) s[kk] = vaddq_f64(s[kk], arr[kk + i]);

	XorShift();

	UNROLL(XS) r[kk] = vaddq_s64(x[kk], y[kk]);

	UNROLL(XS) r[kk] = vandq_u64(r[kk], mask1);
	
	UNROLL(XS) r[kk] = vorrq_u64(r[kk], mask2);

	UNROLL(XS) t[kk] = vsubq_f64(t[kk], one);

	UNROLL(XS) t[kk] = vmulq_f64(t[kk], s[kk]);

	uint64x2_t midx[XS], nidx = vdupq_n_u64(0), ninc = vdupq_n_u64(1);
	uint64x2_t f[XS], b[XS];
	UNROLL(XS) midx[kk] = vdupq_n_u64(n - 1);
	UNROLL(XS) f[kk] = vdupq_n_f64(0);

	for (int i = 0; i < n * XS; i += XS)
	{
		UNROLL(XS) b[kk] = vcltq_f64(t[kk], arr[kk + i]);

		UNROLL(XS) t[kk] = vsubq_f64(t[kk], arr[kk + i]);

		UNROLL(XS) b[kk] = vbicq_u64(b[kk], f[kk]);

		UNROLL(XS) f[kk] = vorrq_u64(f[kk], b[kk]);

		UNROLL(XS) midx[kk] = vbslq_s64(b[kk], nidx, midx[kk]);

		nidx = vaddq_s64(nidx, ninc);
	}

	UNROLL(XS) re64[kk] = midx[kk];
}

/* Draw uniform distriubted intergers */
TARGETNEO void RNGNEO<double>::XorShift()
{
	uint64x2_t a[XS], b[XS];

	UNROLL(XS) a[kk] = x[kk];

	UNROLL(XS) b[kk] = y[kk];

	UNROLL(XS) x[kk] = b[kk];

	UNROLL(XS) a[kk] = veorq_u64(a[kk], vshlq_n_u64(a[kk], 23));

	UNROLL(XS) a[kk] = veorq_u64(a[kk], vshlq_n_u64(a[kk], 18));

	UNROLL(XS) a[kk] = veorq_u64(a[kk], b[kk]);

	UNROLL(XS) a[kk] = veorq_u64(a[kk], vshlq_n_u64(b[kk], 5));

	UNROLL(XS) y[kk] = a[kk];
}

/* Draw uniform distriubted integers */
template<typename INT>
TARGETNEO void RNGNEO<double>::Integer(INT* re, int64 n, INT minv, INT maxv)
{
	constexpr int xesize = sizeof(x) / sizeof(INT);
	int64 i = 0;
	INT modv = maxv - minv;
	INT* rei = re;

	for (; i <= n - xesize; i += xesize)
	{
		XorShift();
		UNROLL(XS) { vst1q_u64((uint64*)rei, vaddq_u64(x[kk], y[kk])); rei += E128B / sizeof(INT); };
	}

	if (i != n)
	{
		uint64x2_t re2[XS];
		XorShift();
		UNROLL(XS) re2[kk] = vaddq_u64(x[kk], y[kk]);
		SetVal((INT*)rei, (INT*)re2, n - i);
	}

	if (maxv != (INT)-1 || minv != 0)
	{
		for (i = 0; i < n; ++i)
			re[i] = re[i] % modv + minv;
	}
}

/* Draw uniform distriubted real numbers */
TARGETNEO void RNGNEO<double>::Uniform(double* re, int n, double minv, double maxv)
{
	constexpr int xesize = sizeof(x) / sizeof(double);
	int i = 0;
	double range = maxv - minv;

	uint64x2_t mask1 = vdupq_n_u64(0x000FFFFFFFFFFFFF);
	uint64x2_t mask2 = vdupq_n_u64(0x3FF0000000000000);
	float64x2_t v1 = vdupq_n_f64(minv - range);
	float64x2_t v2 = vdupq_n_f64(range);

	if (range == 1.0)
	{
		for (; i <= n - xesize; i += xesize)
		{
			XorShift();
			UNROLL(XS) { vst1q_f64(re, vaddq_f64(vreinterpretq_f64_u64(vorrq_u64(vandq_u64(vaddq_u64(x[kk], y[kk]), mask1), mask2)), v1)); re += E128D; }
		}
	}
	else
	{
		for (; i <= n - xesize; i += xesize)
		{
			XorShift();
			UNROLL(XS) { vst1q_f64(re, _neo_fmaddx_pd(vreinterpretq_f64_u64(vorrq_u64(vandq_u64(vaddq_u64(x[kk], y[kk]), mask1), mask2)), v2, v1)); re += E128D; }
		}
	}

	if (i != n)
	{
		float64x2_t ref[XS];
		XorShift();
		UNROLL(XS) ref[kk] = _neo_fmaddx_pd(vreinterpretq_f64_u64(vorrq_u64(vandq_u64(vaddq_u64(x[kk], y[kk]), mask1), mask2)), v2, v1);
		SetVal((double*)re, (double*)ref, n - i);
	}
}

/* Draw uniform distriubted real numbers */
TARGETNEO void RNGNEO<double>::Normal(double* re, int n, double mean, double sd)
{
	constexpr int xhsize = XS / 2;
	constexpr int xesize = XS * E128D;

	int i = 0;

	uint64x2_t mask1 = vdupq_n_u64(0x000FFFFFFFFFFFFF);
	uint64x2_t mask2 = vdupq_n_u64(0x3FF0000000000000);
	float64x2_t v1 = vdupq_n_f64(-1);
	float64x2_t min_freq = vdupq_n_f64(MIN_FREQ);
	float64x2_t pi2 = vdupq_n_f64(2.0 * M_PI);
	float64x2_t mu = vdupq_n_f64(mean);
	float64x2_t s = vdupq_n_f64(sd);

	for (; i <= n - xesize; i += xesize)
	{
		XorShift();
		UNROLL(XS) vst1q_f64(re + E128D * kk, vaddq_f64(vreinterpretq_f64_u64(vorrq_u64(vandq_u64(vaddq_u64(x[kk], y[kk]), mask1), mask2)), v1));

		float64x2_t u1, u2, u3, u4;
		for (int j = 0; j < xhsize; ++j)
		{
			u1 = vmaxq_f64(vld1q_f64(re + j * E128D), min_freq);
			u2 = vmulq_f64(vld1q_f64(re + j * E128D + xhsize * E128D), pi2);

			UNROLL(E128D) simd_f64(u1, kk) = sqrt(-2.0 * log(simd_f64(u1, kk)));
			UNROLL(E128D) simd_f64(u3, kk) = cos(simd_f64(u2, kk));
			UNROLL(E128D) simd_f64(u4, kk) = sin(simd_f64(u2, kk));
				
			vst1q_f64(re + j * E128D, vmulq_f64(u1, u3));
			vst1q_f64(re + j * E128D + xhsize * E128D, vmulq_f64(u1, u4));
		}
	}
	
	if (sd != 1 || mean != 0)
		UNROLL(XS) { vst1q_f64(re, _neo_fmaddx_pd(vld1q_f64(re), s, mu)); re += E128D; }
	else
		re += XS * E128D;

	if (i != n)
	{
		double ref[XS * E128D];

		XorShift();
		UNROLL(XS) vst1q_f64(ref + E128D * kk, vaddq_f64(vreinterpretq_f64_u64(vorrq_u64(vandq_u64(vaddq_u64(x[kk], y[kk]), mask1), mask2)), v1));

		float64x2_t u1, u2, u3, u4;
		for (int j = 0; j < xhsize; ++j)
		{
			u1 = vmaxq_f64(vld1q_f64(ref + j * E128D), min_freq);
			u2 = vmulq_f64(vld1q_f64(ref + j * E128D + xhsize * E128D), pi2);

			UNROLL(E128D) simd_f64(u1, kk) = sqrt(-2.0 * log(simd_f64(u1, kk)));
			UNROLL(E128D) simd_f64(u3, kk) = cos(simd_f64(u2, kk));
			UNROLL(E128D) simd_f64(u4, kk) = sin(simd_f64(u2, kk));

			vst1q_f64(ref + j * E128D, vmulq_f64(u1, u3));
			vst1q_f64(ref + j * E128D + xhsize * E128D, vmulq_f64(u1, u4));
		}
		
		if (sd != 1 || mean != 0)
			UNROLL(XS) vst1q_f64(ref + kk * E128D, _neo_fmaddx_pd(vld1q_f64(ref), s, mu)); 

		SetVal((double*)re, (double*)ref, n - i);
	}
}
#endif

#ifndef _RNGNEO_FP32
#define XS 16
#define XS2 32

/* Initialize rng */
TARGETNEO RNGNEO<float>::RNGNEO()
{

}

/* Initialize rng */
TARGETNEO RNGNEO<float>::RNGNEO(uint64 seed, uint64 salt)
{
	uint32x4_t a[XS], s, m;
	UNROLL(XS) { a[kk] = vld1q_u32(((uint[]) { Mix(seed + 0), Mix(seed + 1), Mix(seed + 2), Mix(seed + 3) })); seed += 4; }

	s = vdupq_n_u32(Mix(salt));
	m = vdupq_n_u32(0x5bd1e995);

	// uint s = s ^ 4;
	s = veorq_u32(s, vdupq_n_u32(4));

	// a *= m;
	UNROLL(XS) a[kk] = vmulq_u32(a[kk], m);

	// a ^= a >> 24;
	UNROLL(XS) a[kk] = veorq_u32(a[kk], vshrq_n_u32(a[kk], 24));

	// a *= m;
	UNROLL(XS) a[kk] = vmulq_u32(a[kk], m);

	// s *= m;
	s = vmulq_u32(s, m);

	// a ^= s;
	UNROLL(XS) a[kk] = veorq_u32(a[kk], s);

	// a ^= a >> 13;
	UNROLL(XS) a[kk] = veorq_u32(a[kk], vshrq_n_u32(a[kk], 13));

	// a *= m;
	UNROLL(XS) a[kk] = vmulq_u32(a[kk], m);

	// a ^= a >> 15;
	UNROLL(XS) a[kk] = veorq_u32(a[kk], vshrq_n_u32(a[kk], 15));

	// original
	UNROLL(XS) x[kk] = veorq_u32(vdupq_n_u32(0x075BCD15), a[kk]);

	UNROLL(XS) a[kk] = vshlq_n_u32(a[kk], 3);

	UNROLL(XS) y[kk] = veorq_u32(vdupq_n_u32(0x159A55E5), a[kk]);

	UNROLL(XS) a[kk] = vshlq_n_u32(a[kk], 3);

	UNROLL(XS) z[kk] = veorq_u32(vdupq_n_u32(0x1F123BB5), a[kk]);
}

/* Draw 64 64-bit integers in [0,n), 64*n frequencies are in arr */
TARGETNEO void RNGNEO<float>::Poly(float32x4_t* arr, int n, void* re)
{
	float32x4_t t[XS], s[XS]; 
	float32x4_t one = vdupq_n_f32(1.0f);
	uint32x4_t mask1 = vdupq_n_u32(0x007FFFFF);
	uint32x4_t mask2 = vdupq_n_u32(0x3F800000);
	uint32x4_t* r = (uint32x4_t*)t; uint64x2_t* re64 = (uint64x2_t*)re; uint32x4_t* re32 = (uint32x4_t*)re;

	UNROLL(XS) s[kk] = vdupq_n_f32(0);

	for (int i = 0; i < n * XS; i += XS)
		UNROLL(XS) s[kk] = vaddq_f32(s[kk], arr[kk + i]);

	XorShift();

	UNROLL(XS) r[kk] = vandq_u32(z[kk], mask1);

	UNROLL(XS) r[kk] = vorrq_u32(r[kk], mask2);

	UNROLL(XS) t[kk] = vsubq_f32(t[kk], one);

	UNROLL(XS) t[kk] = vmulq_f32(t[kk], s[kk]);

	uint32x4_t midx[XS], nidx = vdupq_n_u32(0), ninc = vdupq_n_u32(1);
	uint32x4_t f[XS], b[XS];
	UNROLL(XS) midx[kk] = vdupq_n_u32(n - 1);
	UNROLL(XS) f[kk] = vdupq_n_u32(0);

	for (int i = 0; i < n * XS; i += XS)
	{
		UNROLL(XS) b[kk] = vcltq_f32(t[kk], arr[kk + i]);

		UNROLL(XS) t[kk] = vsubq_f32(t[kk], arr[kk + i]);

		UNROLL(XS) b[kk] = vbicq_u32(b[kk], f[kk]);

		UNROLL(XS) f[kk] = vorrq_u32(f[kk], b[kk]);
		
		UNROLL(XS) midx[kk] = vbslq_s32(b[kk], nidx, midx[kk]);

		nidx = vaddq_s32(nidx, ninc);
	}

	UNROLL(XS) 
	{
		re64[0 + (kk << 1)] = vld1q_u64(((uint64[]) { vgetq_lane_u32(midx[kk], 0), vgetq_lane_u32(midx[kk], 1) }));
		re64[1 + (kk << 1)] = vld1q_u64(((uint64[]) { vgetq_lane_u32(midx[kk], 2), vgetq_lane_u32(midx[kk], 3) }));
	}
}

/* Draw uniform distriubted intergers */
TARGETNEO void RNGNEO<float>::XorShift()
{
	uint32x4_t u[XS];

	UNROLL(XS) u[kk] = vshlq_n_u32(x[kk], 16);
	UNROLL(XS) x[kk] = veorq_u32(x[kk], u[kk]);

	UNROLL(XS) u[kk] = vshrq_n_u32(x[kk], 5);
	UNROLL(XS) x[kk] = veorq_u32(x[kk], u[kk]);

	UNROLL(XS) u[kk] = vshlq_n_u32(x[kk], 1);
	UNROLL(XS) x[kk] = veorq_u32(x[kk], u[kk]);

	UNROLL(XS) u[kk] = x[kk];

	UNROLL(XS) x[kk] = y[kk];

	UNROLL(XS) y[kk] = z[kk];

	UNROLL(XS) z[kk] = veorq_u32(u[kk], x[kk]);

	UNROLL(XS) z[kk] = veorq_u32(z[kk], y[kk]);
}

/* Draw uniform distriubted integers */
template<typename INT>
TARGETNEO void RNGNEO<float>::Integer(INT* re, int64 n, INT minv, INT maxv)
{
	constexpr int xesize = sizeof(x) / sizeof(INT);
	int64 i = 0;
	INT modv = maxv - minv;
	INT* rei = re;

	for (; i <= n - xesize; i += xesize)
	{
		XorShift();
		UNROLL(XS) { vst1q_u32((uint*)rei, z[kk]); rei += E128B / sizeof(INT); }
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
TARGETNEO void RNGNEO<float>::Uniform(float* re, int n, float minv, float maxv)
{
	constexpr int xesize = sizeof(x) / sizeof(float);
	int i = 0;
	float range = maxv - minv;

	uint32x4_t mask1 = vdupq_n_u32(0x007FFFFF);
	uint32x4_t mask2 = vdupq_n_u32(0x3F800000);
	float32x4_t v1 = vdupq_n_f32(minv - range);
	float32x4_t v2 = vdupq_n_f32(range);

	if (range == 1.0)
	{
		for (; i <= n - xesize; i += xesize)
		{
			XorShift();
			UNROLL(XS) { vst1q_f32(re, vaddq_f32(vreinterpretq_f32_u32(vorrq_u32(vandq_u32(z[kk], mask1), mask2)), v1)); re += E128F; }
		}
	}
	else
	{
		for (; i <= n - xesize; i += xesize)
		{
			XorShift();
			UNROLL(XS) { vst1q_f32(re, _neo_fmaddx_ps(vreinterpretq_f32_u32(vorrq_u32(vandq_u32(z[kk], mask1), mask2)), v2, v1)); re += E128F; }
		}
	}

	if (i != n)
	{
		float32x4_t ref[XS];
		XorShift();
		UNROLL(XS) ref[kk] = _neo_fmaddx_ps(vreinterpretq_f32_u32(vorrq_u32(vandq_u32(z[kk], mask1), mask2)), v2, v1);
		SetVal((float*)re, (float*)ref, n - i);
	}
}

/* Draw uniform distriubted real numbers */
TARGETNEO void RNGNEO<float>::Normal(float* re, int n, float mean, float sd)
{
	constexpr int xhsize = XS / 2;
	constexpr int xesize = XS * E128F;

	int i = 0;

	uint32x4_t mask1 = vdupq_n_u32(0x007FFFFF);
	uint32x4_t mask2 = vdupq_n_u32(0x3F800000);
	float32x4_t v1 = vdupq_n_f32(-1);
	float32x4_t min_freq = vdupq_n_f32((float)MIN_FREQ);
	float32x4_t pi2 = vdupq_n_f32((float)(2.0 * M_PI));
	float32x4_t mu = vdupq_n_f32(mean);
	float32x4_t s = vdupq_n_f32(sd);

	for (; i <= n - xesize; i += xesize)
	{
		XorShift();
		UNROLL(XS) vst1q_f32(re + kk * E128F, vaddq_f32(vreinterpretq_f32_u32(vorrq_u32(vandq_u32(z[kk], mask1), mask2)), v1));

		float32x4_t u1, u2, u3, u4;
		for (int j = 0; j < xhsize; ++j)
		{
			u1 = vmaxq_f32(vld1q_f32(re + j * E128F), min_freq);
			u2 = vmulq_f32(vld1q_f32(re + j * E128F + xhsize * E128F), pi2);

			UNROLL(E128F) simd_f32(u1, kk) = sqrt(-2.0 * log(simd_f32(u1, kk)));
			UNROLL(E128F) simd_f32(u3, kk) = cos(simd_f32(u2, kk));
			UNROLL(E128F) simd_f32(u4, kk) = sin(simd_f32(u2, kk));
			
			vst1q_f32(re + j * E128F, vmulq_f32(u1, u3));
			vst1q_f32(re + j * E128F + xhsize * E128F, vmulq_f32(u1, u4));
		}
			
		if (sd != 1 || mean != 0)
			UNROLL(XS) { vst1q_f32(re, _neo_fmaddx_ps(vld1q_f32(re), s, mu)); re += E128F; }
		else
			re += XS * E128F;
	}

	if (i != n)
	{
		float ref[XS * E128F];

		XorShift();
		UNROLL(XS) vst1q_f32(ref + kk * E128F, vaddq_f32(vreinterpretq_f32_u32(vorrq_u32(vandq_u32(z[kk], mask1), mask2)), v1));

		float32x4_t u1, u2, u3, u4;
		for (int j = 0; j < xhsize; ++j)
		{
			u1 = vmaxq_f32(vld1q_f32(ref + j * E128F), min_freq);
			u2 = vmulq_f32(vld1q_f32(ref + j * E128F + xhsize * E128F), pi2);

			UNROLL(E128F) simd_f32(u1, kk) = sqrt(-2.0 * log(simd_f32(u1, kk)));
			UNROLL(E128F) simd_f32(u3, kk) = cos(simd_f32(u2, kk));
			UNROLL(E128F) simd_f32(u4, kk) = sin(simd_f32(u2, kk));
			
			vst1q_f32(ref + j * E128F, vmulq_f32(u1, u3));
			vst1q_f32(ref + j * E128F + xhsize * E128F, vmulq_f32(u1, u4));
		}
			
		if (sd != 1 || mean != 0)
			UNROLL(XS) vst1q_f32(ref + kk * E128F, _neo_fmaddx_ps(vld1q_f32(re), s, mu)); 

		SetVal((float*)re, (float*)ref, n - i);
	}
}
#endif

TARGETNEO int64 GetMinIdxNEO(double* A, int64 n, double& val)
{
#define N 4
	int64 i = 0;
	val = DBL_MAX;
	uint64 idx = (uint64)-1;

	if (n >= N * E128D)
	{
		float64x2_t min1[N];
		uint64x2_t f, midx[N], nidx[N], msep = vdupq_n_u64(8);
		UNROLL(N) min1[kk] = vdupq_n_f64(val);
		UNROLL(N) midx[kk] = vdupq_n_u64(0xFFFFFFFFFFFFFFFF);
		UNROLL(N) nidx[kk] = vld1q_u64(((uint64[]) { 0ull + (kk << 1), 1ull + (kk << 1) }));

		for (int64 l1 = n - N * E128D; i <= l1; i += N * E128D)
		{
			UNROLL(N) 
			{
				f = vcgtq_f64(min1[kk], vld1q_f64(A));
				min1[kk] = vbslq_f64(f, vld1q_f64(A), min1[kk]);
                midx[kk] = vbslq_u64(f, nidx[kk], midx[kk]);
                nidx[kk] = vaddq_u64(nidx[kk], msep);
				A += E128D;
			}
		}

		REDUCE(min1)
		{
			f = vcgtq_f64(min1[kk], min1[kk + KK]);
			min1[kk] = vbslq_f64(f, min1[kk + KK], min1[kk]);
			midx[kk] = vbslq_u64(f, midx[kk + KK], midx[kk]);
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

TARGETNEO int64 GetMinIdxNEO(float* A, int64 n, float& val)
{
#define N 4
	int64 i = 0;
	val = FLT_MAX;
	uint idx = (uint)-1;

	if (n >= N * E128F)
	{
		float32x4_t min1[N];
		uint32x4_t f[N], midx[N], nidx[N], msep = vdupq_n_u32(N * E128F);
		UNROLL(N) min1[kk] = vdupq_n_f32(val);
		UNROLL(N) midx[kk] = vdupq_n_u32(0xFFFFFFFF);
		UNROLL(N) nidx[kk] = vld1q_u32(((uint[]) { 0u + (kk << 2), 1u + (kk << 2), 2u + (kk << 2), 3u + (kk << 2) }));

		for (int64 l1 = n - N * E128F; i <= l1; i += N * E128F)
		{
			UNROLL(N) 
			{
				f[kk] = vcgtq_f32(min1[kk], vld1q_f32(A));
				min1[kk] = vbslq_f32(f[kk], vld1q_f32(A), min1[kk]);
				A += E128F; 
			}

			UNROLL(N) midx[kk] = vbslq_s32(f[kk], nidx[kk], midx[kk]);

			UNROLL(N) nidx[kk] = vaddq_s32(nidx[kk], msep);
		}

		REDUCE(min1)
		{
			f[kk] = vcgtq_f32(min1[kk], min1[kk + KK]);
			min1[kk] = vbslq_s32(f[kk], min1[kk + KK], min1[kk]);
			midx[kk] = vbslq_s32(f[kk], midx[kk + KK], midx[kk]);
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

TARGETNEO void GetMinMaxValNEO(double* A, int64 n, double& minv, double& maxv)
{
#define N 4
	int64 i = 0;
	minv = DBL_MAX;
	maxv = -DBL_MAX;

	if (n >= N * E128D)
	{
		float64x2_t min1[N], max1[N];
		UNROLL(N) min1[kk] = vdupq_n_f64(minv);
		UNROLL(N) max1[kk] = vdupq_n_f64(maxv);

		for (int64 l1 = n - N * E128D; i <= l1; i += N * E128D)
		{
			UNROLL(N)
			{
				min1[kk] = vminq_f64(min1[kk], vld1q_f64(A));
				max1[kk] = vmaxq_f64(max1[kk], vld1q_f64(A));
				A += E128D;
			}
		}

		REDUCE(min1)
		{
			min1[kk] = vminq_f64(min1[kk], min1[kk + KK]);
			max1[kk] = vmaxq_f64(max1[kk], max1[kk + KK]);
		}
		
		minv = _neo_reduce_min_pd(min1[0]);
		maxv = _neo_reduce_max_pd(max1[0]);
	}

	for (; i < n; ++i, ++A)
	{
		minv = std::min(minv, *A);
		maxv = std::max(maxv, *A);
	}
}

TARGETNEO void GetMinMaxValNEO(float* A, int64 n, float& minv, float& maxv)
{
#define N 4
	int64 i = 0;
	minv = FLT_MAX;
	maxv = -FLT_MAX;

	if (n >= N * E128F)
	{
		float32x4_t min1[N], max1[N];
		UNROLL(N) min1[kk] = vdupq_n_f32(minv);
		UNROLL(N) max1[kk] = vdupq_n_f32(maxv);

		for (int64 l1 = n - N * E128F; i <= l1; i += N * E128F)
		{
			UNROLL(N)
			{
				min1[kk] = vminq_f32(min1[kk], vld1q_f32(A));
				max1[kk] = vmaxq_f32(max1[kk], vld1q_f32(A));
				A += E128F;
			}
		}

		REDUCE(min1)
		{
			min1[kk] = vminq_f32(min1[kk], min1[kk + KK]);
			max1[kk] = vmaxq_f32(max1[kk], max1[kk + KK]);
		}
		
		minv = _neo_reduce_min_ps(min1[0]);
		maxv = _neo_reduce_max_ps(max1[0]);
	}

	for (; i < n; ++i, ++A)
	{
		minv = std::min(minv, *A);
		maxv = std::max(maxv, *A);
	}
}

TARGETNEO double GetMaxValNEO(double* A, int64 n)
{
#define N 4
	int64 i = 0;
	double val = -DBL_MAX;

	if (n >= N * E128D)
	{
		float64x2_t max1[N];
		UNROLL(N) max1[kk] = vdupq_n_f64(val);

		for (int64 l1 = n - N * E128D; i <= l1; i += N * E128D)
		{
			UNROLL(N) { max1[kk] = vmaxq_f64(max1[kk], vld1q_f64(A)); A += E128D; }
		}

		REDUCE(max1) max1[kk] = vmaxq_f64(max1[kk], max1[kk + KK]);
		
		val = _neo_reduce_max_pd(max1[0]);
	}

	for (; i < n; ++i, ++A)
	{
		val = std::max(val, *A);
	}

	return val;
}

TARGETNEO float GetMaxValNEO(float* A, int64 n)
{
#define N 8
	int64 i = 0;
	float val = -FLT_MAX;

	if (n >= N * E128F)
	{
		float32x4_t max1[N];
		UNROLL(N) max1[kk] = vdupq_n_f32(val);

		for (int64 l1 = n - N * E128F; i <= l1; i += N * E128F)
		{
			UNROLL(N) { max1[kk] = vmaxq_f32(max1[kk], vld1q_f32(A)); A += E128F; }
		}

		REDUCE(max1) max1[kk] = vmaxq_f32(max1[kk], max1[kk + KK]);
		
		val = _neo_reduce_max_ps(max1[0]);
	}

	for (; i < n; ++i, ++A)
	{
		val = std::max(val, *A);
	}

	return val;
}

TARGETNEO double GetMaxValNEO(double* A, int64 n, int64 sep)
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
		float64x2_t max1[N];
		UNROLL(N) max1[kk] = vdupq_n_f64(val);

		for (int64 l1 = n - N * E128D; i <= l1; i += N * E128D)
		{
			UNROLL(N)
			{
				max1[kk] = vmaxq_f64(max1[kk], vld1q_f64(((double[]) { A[0 * sep], A[1 * sep], A[2 * sep], A[3 * sep] })));
				A += E128D * sep;
			}
		}

		REDUCE(max1) max1[kk] = vmaxq_f64(max1[kk], max1[kk + KK]);
		
		val = _neo_reduce_max_epi64(max1[0]);
	}

	for (; i < n; ++i, A += sep)
	{
		val = std::max(val, *A);
	}

	return val;
}

TARGETNEO float GetMaxValNEO(float* A, int64 n, int64 sep)
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
		float32x4_t max1[N];
		UNROLL(N) max1[kk] = vdupq_n_f32(val);

		for (int64 l1 = n - N * E128F; i <= l1; i += N * E128F)
		{
			UNROLL(N)
			{
				max1[kk] = vmaxq_f32(max1[kk], vld1q_f32(((float[]) { A[0 * sep], A[1 * sep], A[2 * sep], A[3 * sep] })));
				A += E128F * sep;
			}
		}

		REDUCE(max1) max1[kk] = vmaxq_f32(max1[kk], max1[kk + KK]);
		
		val = _neo_reduce_max_ps(max1[0]);
	}

	for (; i < n; ++i, A += sep)
	{
		val = std::min(val, *A);
	}

	return val;
}

TARGETNEO double GetMinValNEO(double* A, int64 n)
{
#define N 4
	int64 i = 0;
	double val = DBL_MAX;

	if (n >= N * E128D)
	{
		float64x2_t min1[N], a[N];
		UNROLL(N) min1[kk] = vdupq_n_f64(val);

		for (int64 l1 = n - N * E128D; i <= l1; i += N * E128D)
		{
			UNROLL(N) { min1[kk] = vminq_f64(min1[kk], vld1q_f64(A)); A += E128D; }
		}

		REDUCE(min1) min1[kk] = vminq_f64(min1[kk], min1[kk + KK]);
		
		val = _neo_reduce_min_pd(min1[0]);
	}

	for (; i < n; ++i, ++A)
	{
		val = std::min(val, *A);
	}

	return val;
}

TARGETNEO float GetMinValNEO(float* A, int64 n)
{
#define N 8
	int64 i = 0;
	float val = FLT_MAX;

	if (n >= N * E128F)
	{
		float32x4_t min1[N], a[N];
		UNROLL(N) min1[kk] = vdupq_n_f32(val);

		for (int64 l1 = n - N * E128F; i <= l1; i += N * E128F)
		{
			UNROLL(N) { min1[kk] = vminq_f32(min1[kk], vld1q_f32(A)); A += E128F; }
		}

		REDUCE(min1) min1[kk] = vminq_f32(min1[kk], min1[kk + KK]);
		
		val = _neo_reduce_min_ps(min1[0]);
	}

	for (; i < n; ++i, ++A)
	{
		val = std::min(val, *A);
	}

	return val;
}

TARGETNEO int64 GetMinValNEO(int64* A, int64 n)
{
#define N 4
	int64 i = 0;
	int64 val = 0x7FFFFFFFFFFFFFFF;

	if (n >= N * E128D)
	{
		int64x2_t min1[N];
		UNROLL(N) min1[kk] = vdupq_n_u64(0x7FFFFFFFFFFFFFFF);

		for (int64 l1 = n - N * E128D; i <= l1; i += N * E128D)
		{
			UNROLL(N) 
			{
				min1[kk] = vbslq_s64(vcgtq_s64(min1[kk], vld1q_s64(A)), vld1q_s64(A), min1[kk]); 

				A += E128D; 
			}
		}

		REDUCE(min1) min1[kk] = vbslq_s64(vcgtq_s64(min1[kk], min1[kk + KK]), min1[kk + KK], min1[kk]);

		val = _neo_reduce_min_epi64(min1[0]);
	}

	for (; i < n; ++i, ++A)
	{
		val = std::min(val, *A);
	}

	return val;
}

TARGETNEO void SetValNEO(uint* A, ushort* B, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E128F)
	{
		for (int64 l1 = n - N * E128F; i <= l1; i += N * E128F)
		{
			UNROLL(N) { vst1q_u32(A, vmovl_u16(vld1_u16(B))); A += E128F; B += E128F;}
		}
	}

	for (; i < n; ++i)
		*A++ = *B++;
}

TARGETNEO void AddExponentNEO(int64& slog, float64x2_t& val)
{
	uint64x2_t& vv = *(uint64x2_t*)&val;
	uint64x2_t mask1 = vdupq_n_u64(0x7FF0000000000000);
	uint64x2_t mask2 = vdupq_n_u64(0x800FFFFFFFFFFFFF);
	uint64x2_t mask3 = vdupq_n_u64(0x3FF0000000000000);
	uint64x2_t subv = vdupq_n_u64(1023);

	uint64x2_t t = vsubq_s64(vshrq_n_u64(vandq_u64(vv, mask1), 52), subv);

	slog += (int)(vgetq_lane_s64(t, 0) + vgetq_lane_s64(t, 1));

	vv = vorrq_u64(vandq_u64(vv, mask2), mask3);
}

TARGETNEO void AddExponentNEO(int64& slog, float32x4_t& val)
{
	uint32x4_t& vv = *(uint32x4_t*)&val;
	uint32x4_t mask1 = vdupq_n_u32(0x7F800000);
	uint32x4_t mask2 = vdupq_n_u32(0x807FFFFF);
	uint32x4_t mask3 = vdupq_n_u32(0x3F800000);
	uint32x4_t subv = vdupq_n_u32(127);

	uint32x4_t t = vsubq_s32(vshrq_n_u32(vandq_u32(vv, mask1), 23), subv);

	slog += (int)vgetq_lane_s32(t, 0) + (int)vgetq_lane_s32(t, 1) + (int)vgetq_lane_s32(t, 2) + (int)vgetq_lane_s32(t, 3);

	vv = vorrq_u32(vandq_u32(vv, mask2), mask3);
}

TARGETNEO void ChargeLogNEO(int64& slog, double& prod, float64x2_t& val)
{
	AddExponentNEO(slog, val);

	prod *= _neo_reduce_mul_pd(val);

	if (prod < DOUBLE_UNDERFLOW || prod > DOUBLE_OVERFLOW) [[unlikely]]
		AddExponent(slog, prod);
}

TARGETNEO void ChargeLogNEO(int64& slog, double& prod, float32x4_t& val)
{
	AddExponentNEO(slog, val);

	prod *= _neo_reduce_mul_psd(val);

	if (prod < DOUBLE_UNDERFLOW || prod > DOUBLE_OVERFLOW) [[unlikely]]
		AddExponent(slog, prod);
}

TARGETNEO double LogProdNEO(double* A, int64 n)
{
#define N 4
	int64 i = 0;
	int64 slog = 0; double prod = 1;

	if (n >= N * E512D)
	{
		float64x2_t pd[E512_128];
		UNROLL(E512_128) pd[kk] = vdupq_n_f64(1.0);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128)
			{ pd[kk] = vmulq_f64(pd[kk], vld1q_f64(A));  A += E128D; }

			UNROLL(E512_128) AddExponentNEO(slog, pd[kk]);
		}
		
		float64x2_t* pd1 = (float64x2_t*)pd;
		UNROLL(E512_128) ChargeLogNEO(slog, prod, pd1[kk]);
	}

	for (; i < n; ++i, ++A)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);
	return prod;
}

TARGETNEO double LogProdNEO(float* A, int64 n)
{
#define N 4
	int64 i = 0;
	int64 slog = 0; double prod = 1;

	if (n >= N * E512D)
	{
        float64x2_t pd[E512_128];
		UNROLL(E512_128) pd[kk] = vdupq_n_f64(1.0);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128)
			{ pd[kk] = vmulq_f64(pd[kk], vcvt_f64_f32(vld1_f32(A))); A += E128D; }

			UNROLL(E512_128) AddExponentNEO(slog, pd[kk]);
		}
		
		float64x2_t* pd1 = (float64x2_t*)pd;
		UNROLL(E512_128) ChargeLogNEO(slog, prod, pd1[kk]);
	}

	for (; i < n; ++i, ++A)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);
	return prod;
}

TARGETNEO double LogProdNEO(double* A, int64 n, int64 sep)
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
		float64x2_t pd[E512_128];
		UNROLL(E512_128) pd[kk] = vdupq_n_f64(1.0);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128) 
			{ pd[kk] = vmulq_f64(pd[kk], vld1q_f64(((double[]) { A[0 * sep], A[1 * sep] }))); A += E128D * sep; }
			
			UNROLL(E512_128) AddExponentNEO(slog, pd[kk]);
		}
		
		UNROLL(E512_128) ChargeLogNEO(slog, prod, pd[kk]);
	}

	for (; i < n; ++i, A += sep)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);

	return prod;
}

TARGETNEO double LogProdNEO(float* A, int64 n, int64 sep)
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
		float64x2_t pd1 = vdupq_n_f64(1.0f), pd2 = vdupq_n_f64(1.0f);
		float32x4_t a1;

		for (int64 l1 = n - N * E128F; i <= l1; i += N * E128F)
		{
			UNROLL(N)
			{
				a1 = vld1q_f32(((float[]) { A[0 * sep], A[1 * sep], A[2 * sep], A[3 * sep] })); A += E128F * sep;

				pd1 = vmulq_f64(pd1, vcvt_f64_f32(vget_low_f32 (a1)));
				pd2 = vmulq_f64(pd2, vcvt_f64_f32(vget_high_f32(a1)));
			}

			AddExponentNEO(slog, pd1);
			AddExponentNEO(slog, pd2);
		}

		ChargeLogNEO(slog, prod, pd1);
		ChargeLogNEO(slog, prod, pd2);
	}

	for (; i < n; ++i, A += sep)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);
	return prod;
}

TARGETNEO double LogProdDivNEO(double* A, double* B, int64 n, int64 sep)
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

TARGETNEO double LogProdDivNEO(float* A, float* B, int64 n, int64 sep)
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
		float64x2_t pd1[E512_128], pd2[E512_128];
		UNROLL(E512_128) pd1[kk] = pd2[kk] = vdupq_n_f64(1.0);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128)
			{
				pd1[kk] = vmulq_f64(pd1[kk], vld1q_f64(((double[]) { A[0 * sep], A[1 * sep] }))); A += sep * 2;
				pd2[kk] = vmulq_f64(pd2[kk], vld1q_f64(((double[]) { B[0 * sep], B[1 * sep] }))); B += sep * 2;
			}

			UNROLL(E512_128) AddExponentNEO(slog1, pd1[kk]);
			UNROLL(E512_128) AddExponentNEO(slog2, pd2[kk]);
		}

		UNROLL(E512_128) ChargeLogNEO(slog1, prod1, pd1[kk]);
		UNROLL(E512_128) ChargeLogNEO(slog2, prod2, pd2[kk]);
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

TARGETNEO int64 CountNonZeroNEO(byte* A, int64 n)
{
	uint64 re = 0;
	int64 i = 0;

	if (n >= E128B)
	{
		uint8x16_t z = vdupq_n_u8(0), v;
		uint64 s1 = 0, s2 = 0;
		byte* Ab = A;

		for (int64 l1 = n - E128B; i <= l1; i += E128B)
		{
			v = vceqq_u8(vld1q_u8(A), z); A += E128B;
			s1 += __builtin_popcountll(vgetq_lane_u64((uint64x2_t)v, 0));
			s2 += __builtin_popcountll(vgetq_lane_u64((uint64x2_t)v, 1));
		}

		re = (A - Ab) - ((s1 + s2) >> 3);
	}

	for (; i < n; ++i, ++A)
		if (*A) re++;

	return (int64)re;
}

TARGETNEO double SumNEO(double* A, int64 n)
{
#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		float64x2_t s[E512_128] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128)
			{ s[kk] = vaddq_f64(s[kk], vld1q_f64(A)); A += E128D; }
		}

		REDUCE(s) s[kk] = vaddq_f64(s[kk], s[kk + KK]);

		re = _neo_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *A++;
		re += v1;
	}

	return re;
}

TARGETNEO double SumNEO(float* A, int64 n)
{
#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		float64x2_t s[E512_128] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128)
			{ s[kk] = vaddq_f64(s[kk], vld1q_f64(((double[]) { A[0], A[1] }))); A += 2; }
		}

		REDUCE(s) s[kk] = vaddq_f64(s[kk], s[kk + KK]);

		re = _neo_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A++;
		re += v1;
	}

	return re;
}

TARGETNEO float SumNEOx(float* A, int64 n)
{
#define N 4
	int64 i = 0;
	volatile float re = 0;

	if (n >= N * E512F)
	{
		float32x4_t s[E512_128] = { 0 };

		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) UNROLL(E512_128)
			{ s[kk] = vaddq_f32(s[kk], vld1q_f32(A)); A += E128F; }
		}

		REDUCE(s) s[kk] = vaddq_f32(s[kk], s[kk + KK]);

		re = _neo_reduce_add_ps(s[0]);
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A++;
		re += v1;
	}

	return re;
}

TARGETNEO int64 SumNEO(byte* A, int64 n)
{
#define N 4
	int64 re = 0;
	int64 i = 0;

	if (n >= N * E128B)
	{
		uint64x2_t s = vdupq_n_u64(0);

		for (int64 l1 = n - N * E128B; i <= l1; i += N * E128B)
		{
			UNROLL(N) 
			{ s = vaddq_u64(s, vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(vld1q_u8(A))))); A += E128B; }
		}

		re += _neo_reduce_add_epi64(s);
	}

	for (; i < n; ++i)
		re += *A++;

	return re;
}

TARGETNEO double SumNEO(double* A, int64 n, int64 sep)
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
		float64x2_t s[E512_128] = { 0 };
		UNROLL(E128D) __builtin_prefetch(&A[kk * sep], 0);
					   
		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128)
			{
				UNROLL(E128D) __builtin_prefetch(&A[E128D * sep + kk * sep], 0);
				s[kk] = vaddq_f64(s[kk], vld1q_f64(((double[]) { A[0 * sep], A[1 * sep] })));
				A += E128D * sep;
			}
		}

		REDUCE(s) s[kk] = vaddq_f64(s[kk], s[kk + KK]);

		re = _neo_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, A += sep)
	{
		volatile double v1 = *A;
		re += v1;
	}

	return re;
}

TARGETNEO double SumNEO(float* A, int64 n, int64 sep)
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
		float64x2_t s[E512_128] = { 0 };
		UNROLL(E128D) __builtin_prefetch(&A[kk * sep], 0);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128)
			{
				UNROLL(E128D) __builtin_prefetch(&A[E128D * sep + kk * sep], 0);
				s[kk] = vaddq_f64(s[kk], vld1q_f64(((double[]) { A[0 * sep], A[1 * sep] })));
				A += E128D * sep;
			}
		}

		REDUCE(s) s[kk] = vaddq_f64(s[kk], s[kk + KK]);

		re = _neo_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, A += sep)
	{
		volatile double v1 = *A;
		re += v1;
	}

	return re;
}

TARGETNEO float SumNEOx(float* A, int64 n, int64 sep)
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
		float32x4_t s[E512_128] = { 0 };
		UNROLL(E128F) __builtin_prefetch(&A[kk * sep], 0);

		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) UNROLL(E512_128)
			{
				UNROLL(E128F) __builtin_prefetch(&A[E128F * sep + kk * sep], 0);
				s[kk] = vaddq_f32(s[kk], vld1q_f32(((float[]) { A[0 * sep], A[1 * sep], A[2 * sep], A[3 * sep] })));
				A += E128F * sep;
			}
		}

		REDUCE(s) s[kk] = vaddq_f32(s[kk], s[kk + KK]);

		re = _neo_reduce_add_ps(s[0]);
	}

	for (; i < n; ++i, A += sep)
	{
		volatile float v1 = *A;
		re += v1;
	}

	return re;
}

TARGETNEO double ProdNEO(double* A, int64 n)
{
#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		float64x2_t s[E512_128];
		UNROLL(E512_128) s[kk] = vdupq_n_f64(1);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128)
			{ s[kk] = vmulq_f64(s[kk], vld1q_f64(A)); A += E128D; }
		}

		REDUCE(s) s[kk] = vmulq_f64(s[kk], s[kk + KK]);

		re = _neo_reduce_mul_pd(s[0]);
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *A++;
		re *= v1;
	}

	return re;
}

TARGETNEO double ProdNEO(float* A, int64 n)
{
#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		float64x2_t s[E512_128];
		UNROLL(E512_128) s[kk] = vdupq_n_f64(1);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128)
			{ s[kk] = vmulq_f64(s[kk], vld1q_f64(((double[]) { A[0], A[1] }))); A += E128D; }
		}

		REDUCE(s) s[kk] = vmulq_f64(s[kk], s[kk + KK]);

		re = _neo_reduce_mul_pd(s[0]);
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A++;
		re *= v1;
	}

	return re;
}

TARGETNEO float ProdNEOx(float* A, int64 n)
{
#define N 4
	int64 i = 0;
	volatile float re = 0;

	if (n >= N * E512F)
	{
		float32x4_t s[E512_128];
		UNROLL(E512_128) s[kk] = vdupq_n_f32(1);

		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) UNROLL(E512_128) 
			{ s[kk] = vmulq_f32(s[kk], vld1q_f32(A)); A += E128F; }
		}

		REDUCE(s) s[kk] = vmulq_f32(s[kk], s[kk + KK]);

		re = _neo_reduce_mul_ps(s[0]);
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A++;
		re *= v1;
	}

	return re;
}

TARGETNEO double ProdNEO(double* A, int64 n, int64 sep)
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
		float64x2_t s[E512_128];
		UNROLL(E512_128) s[kk] = vdupq_n_f64(1);
		UNROLL(E128D) __builtin_prefetch(&A[kk * sep], 0);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128)
			{
				UNROLL(E128D) __builtin_prefetch(&A[E128D * sep + kk * sep], 0);
				s[kk] = vmulq_f64(s[kk], vld1q_f64(((double[]) { A[0 * sep], A[1 * sep] })));
				A += E128D * sep;
			}
		}

		REDUCE(s) s[kk] = vmulq_f64(s[kk], s[kk + KK]);

		re = _neo_reduce_mul_pd(s[0]);
	}

	for (; i < n; ++i, A += sep)
	{
		volatile double v1 = *A;
		re *= v1;
	}

	return re;
}

TARGETNEO double ProdNEO(float* A, int64 n, int64 sep)
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
		float64x2_t s[E512_128];
		UNROLL(E512_128) s[kk] = vdupq_n_f64(1);
		UNROLL(E128D) __builtin_prefetch(&A[kk * sep], 0);

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128)
			{
				UNROLL(E128D) __builtin_prefetch(&A[E128D * sep + kk * sep], 0);
				s[kk] = vmulq_f64(s[kk], vld1q_f64(((double[]) { A[0 * sep], A[1 * sep] })));
				A += E128D * sep;
			}
		}

		REDUCE(s) s[kk] = vmulq_f64(s[kk], s[kk + KK]);

		re = _neo_reduce_mul_pd(s[0]);
	}

	for (; i < n; ++i, A += sep)
	{
		volatile double v1 = *A;
		re *= v1;
	}

	return re;
}

TARGETNEO float ProdNEOx(float* A, int64 n, int64 sep)
{
	//suboptimal to compile
	{
		float re = 1;
		UNROLLHEAD(4)
		for (int64 i = 0; i < n; ++i, A += sep)
			re *= *A;
		return re;
	}

#define N 2
	int64 i = 0;
	volatile float re = 0;

	if (n >= N * E512F)
	{
		float32x4_t s[E512_128];
		UNROLL(E512_128) s[kk] = vdupq_n_f32(1);
		UNROLL(E128F) __builtin_prefetch(&A[kk * sep], 0);

		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) UNROLL(E512_128)
			{
				UNROLL(E128F) __builtin_prefetch(&A[E128F * sep + kk * sep], 0);
				s[kk] = vmulq_f32(s[kk], vld1q_f32(((float[]) { A[0 * sep], A[1 * sep], A[2 * sep], A[3 * sep] })));
				A += E128F * sep;
			}
		}

		REDUCE(s) s[kk] = vmulq_f32(s[kk], s[kk + KK]);

		re = _neo_reduce_mul_ps(s[0]);
	}

	for (; i < n; ++i, A += sep)
	{
		volatile float v1 = *A;
		re *= v1;
	}

	return re;
}

TARGETNEO double SumSquareNEO(double* A, int64 n)
{
#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		float64x2_t s[E512_128] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128)
			{ s[kk] = _neo_fmaddx_pd(vld1q_f64(A), vld1q_f64(A), s[kk]); A += E128D; }
		}

		REDUCE(s) s[kk] = vaddq_f64(s[kk], s[kk + KK]);

		re = _neo_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, ++A)
	{
		volatile double v1 = *A * *A;
		re += v1;
	}

	return re;
}

TARGETNEO double SumSquareNEO(float* A, int64 n)
{
#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		float64x2_t s[E512_128] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128)
			{
				float64x2_t v1 = vcvt_f64_f32(vld1_f32(A)); A += E128D;
				s[kk] = _neo_fmaddx_pd(v1, v1, s[kk]); 
			}
		}

		REDUCE(s) s[kk] = vaddq_f64(s[kk], s[kk + KK]);

		re = _neo_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, ++A)
	{
		volatile double v1 = (double)*A * (double)*A;
		re += v1;
	}

	return re;
}

TARGETNEO int64 SumSquareNEO(byte* A, int64 n)
{
#define N 4
	int64 i = 0;
	uint64 re = 0;

	if (n >= N * E128B)
	{
		uint64x2_t t = vdupq_n_u64(0); 
		uint16x8_t s = vdupq_n_u16(0);

		for (int64 l1 = n - N * E128B; i <= l1; i += N * E128B)
		{
			UNROLL(N) { s = vaddq_u16(s, vpaddlq_u8(vmulq_u8(vld1q_u8(A), vld1q_u8(A)))); A += E128B; }
			
			if ((i & (E128B * 128 - 1)) == 0) [[unlikely]]
			{
				t = vaddq_u64(t, vpaddlq_u32(vpaddlq_u16(s)));
				s = vdupq_n_u16(0);
			}
		}
		
		t = vaddq_u64(t, vpaddlq_u32(vpaddlq_u16(s)));
		re = _neo_reduce_add_epi64(t);
	}

	for (; i < n; ++i, ++A)
		re += *A * *A;

	return re;
}

TARGETNEO void SumSumSquareNEO(double* A, int64 n, double& sum, double& sumsq)
{
#define N 4
	int64 i = 0;
	volatile double re1 = 0, re2 = 0;

	if (n >= N * E512D)
	{
		float64x2_t s1[E512_128] = { 0 }, s2[E512_128] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128)
			{
				float64x2_t v1 = vld1q_f64(A); A += E128D;

				s1[kk] = vaddq_f64(v1, s1[kk]);

				s2[kk] = _neo_fmaddx_pd(v1, v1, s2[kk]);
			}
		}

		REDUCE(s1)
		{
			s1[kk] = vaddq_f64(s1[kk], s1[kk + KK]);
			s2[kk] = vaddq_f64(s2[kk], s2[kk + KK]);
		}

		re1 = _neo_reduce_add_pd(s1[0]);
		re2 = _neo_reduce_add_pd(s2[0]);
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

TARGETNEO void SumSumSquareNEO(float* A, int64 n, double& sum, double& sumsq)
{
#define N 4
	int64 i = 0;
	volatile double re1 = 0, re2 = 0;

	if (n >= N * E512D)
	{
		float64x2_t s1[E512_128] = { 0 }, s2[E512_128] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128)
			{
				float64x2_t v1 = vcvt_f64_f32(vld1_f32(A)); A += E128D;

				s1[kk] = vaddq_f64(v1, s1[kk]);

				s2[kk] = _neo_fmaddx_pd(v1, v1, s2[kk]);
			}
		}

		REDUCE(s1)
		{
			s1[kk] = vaddq_f64(s1[kk], s1[kk + KK]);
			s2[kk] = vaddq_f64(s2[kk], s2[kk + KK]);
		}

		re1 = _neo_reduce_add_pd(s1[0]);
		re2 = _neo_reduce_add_pd(s2[0]);
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

TARGETNEO double SumProdDivNEO(double* A1, double* A2, double* B, int64 sep, int64 n)
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
		UNROLL(E128D) { __builtin_prefetch(B, 0); B += sep; }

		float64x2_t s1[E512_128] = { 0 }, s2[E512_128] = { 0 };
		
		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128)
			{
				UNROLL(E128D) { __builtin_prefetch(B, 0); B += sep; }

				float64x2_t b = vld1q_f64(((double[]) { B[-4 * sep], B[-3 * sep] }));
				
				s1[kk] = vaddq_f64(vmulq_f64(b, vld1q_f64(A1)), s1[kk]); A1 += E128D;

				s2[kk] = vaddq_f64(vmulq_f64(b, vld1q_f64(A2)), s2[kk]); A2 += E128D;
			}
		}

		REDUCE(s1)
		{
			s1[kk] = vaddq_f64(s1[kk], s1[kk + KK]);
			s2[kk] = vaddq_f64(s2[kk], s2[kk + KK]);
		}

		re1 = _neo_reduce_add_pd(s1[0]);
		re2 = _neo_reduce_add_pd(s2[0]);

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

TARGETNEO double SumProdDivNEO(double* A1, float* A2, float* B, int64 sep, int64 n)
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
		UNROLL(E128D) { __builtin_prefetch((const char*)B, 0); B += sep; }

		float64x2_t s1[E512_128] = { 0 }, s2[E512_128] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128)
			{
				UNROLL(E128D) { __builtin_prefetch((const char*)B, 0); B += sep; }

				float64x2_t b = vld1q_f64(((double[]) { B[-4 * sep], B[-3 * sep] }));

				s1[kk] = vaddq_f64(vmulq_f64(vld1q_f64(A1), b), s1[kk]); A1 += E128D;

				s2[kk] = vaddq_f64(vmulq_f64(vcvt_f64_f32(vld1_f32(A2)), b), s2[kk]); A2 += E128D;
			}
		}

		REDUCE(s1)
		{
			s1[kk] = vaddq_f64(s1[kk], s1[kk + KK]);
			s2[kk] = vaddq_f64(s2[kk], s2[kk + KK]);
		}

		re1 = _neo_reduce_add_pd(s1[0]);
		re2 = _neo_reduce_add_pd(s2[0]);

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

TARGETNEO double SumProdDivNEO(float* A1, float* A2, float* B, int64 sep, int64 n)
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
		UNROLL(E128D) { __builtin_prefetch(B, 0); B += sep; }
		
		float64x2_t s1[E512_128] = { 0 }, s2[E512_128] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128)
			{
				UNROLL(E128D) { __builtin_prefetch(B, 0); B += sep; }

				float64x2_t b = vld1q_f64(((double[]) { B[-4 * sep], B[-3 * sep] }));
				
				s1[kk] = _neo_fmaddx_pd(vcvt_f64_f32(vld1_f32(A1)), b, s1[kk]); A1 += E128D;

				s2[kk] = _neo_fmaddx_pd(vcvt_f64_f32(vld1_f32(A2)), b, s2[kk]); A2 += E128D;
			}
		}

		REDUCE(s1)
		{
			s1[kk] = vaddq_f64(s1[kk], s1[kk + KK]);
			s2[kk] = vaddq_f64(s2[kk], s2[kk + KK]);
		}

		re1 = _neo_reduce_add_pd(s1[0]);
		re2 = _neo_reduce_add_pd(s2[0]);

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

TARGETNEO float SumProdDivNEOx(float* A1, float* A2, float* B, int64 sep, int64 n)
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

TARGETNEO double SumProdNEO(double* A, double* B, int64 sep, int64 n)
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
		UNROLL(E128D) { __builtin_prefetch(B, 0); B += sep; }
		
		float64x2_t s[E512_128] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N)
			{
				UNROLL(E512_128)
				{
					UNROLL(E128D) { __builtin_prefetch(B, 0); B += sep; }

					s[kk] = _neo_fmaddx_pd(vld1q_f64(A), vld1q_f64(((double[]) { B[-4 * sep], B[-3 * sep] })), s[kk]);  A += E128D;
				}
			}
		}

		REDUCE(s) s[kk] = vaddq_f64(s[kk], s[kk + KK]); 

		re = _neo_reduce_add_pd(s[0]);
        
		B -= E128D * sep;
	}

	for (; i < n; ++i, ++A, ++B)
	{
		volatile double v1 = *A * *B;
		re += v1;
	}

	return re;
}

TARGETNEO double SumProdNEO(float* A, float* B, int64 sep, int64 n)
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
		UNROLL(E128D) { __builtin_prefetch(B, 0); B += sep; }
		
		float64x2_t s[E512_128] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128)
			{
				UNROLL(E128D) { __builtin_prefetch(B, 0); B += sep; }

				s[kk] = _neo_fmaddx_pd(vcvt_f64_f32(vld1_f32(A)), vld1q_f64(((double[]) { B[-4 * sep], B[-3 * sep] })), s[kk]);
			
				A += E128D;
			}
		}

		REDUCE(s) s[kk] = vaddq_f64(s[kk], s[kk + KK]);

		re = _neo_reduce_add_pd(s[0]);

		B -= E128D * sep;
	}

	for (; i < n; ++i, A++, B += sep)
	{
		volatile double v1 = (double)*A * (double)*B;
		re += v1;
	}

	return re;
}

TARGETNEO float SumProdNEOx(float* A, float* B, int64 sep, int64 n)
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

TARGETNEO double SumProdNEO(double* A, double* B, int64 n)
{
#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		float64x2_t s[E512_128] = { 0 } ;

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128) 
			{ s[kk] = _neo_fmaddx_pd(vld1q_f64(A), vld1q_f64(B), s[kk]); A += E128D; B += E128D; }
		}

		REDUCE(s) s[kk] = vaddq_f64(s[kk], s[kk + KK]);

		re = _neo_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, ++A, ++B)
	{
		volatile double v1 = (double)*A * (double)*B;
		re += v1;
	}

	return re;
}

TARGETNEO double SumProdNEO(float* A, float* B, int64 n)
{
#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		float64x2_t s[E512_128] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128) 
			{ s[kk] = _neo_fmaddx_pd(vcvt_f64_f32(vld1_f32(A)), vcvt_f64_f32(vld1_f32(B)), s[kk]); A += E128D; B += E128D; }
		}

		REDUCE(s) s[kk] = vaddq_f64(s[kk], s[kk + KK]);

		re = _neo_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, ++A, ++B)
	{
		volatile double v1 = (double)*A * (double)*B;
		re += v1;
	}

	return re;
}

TARGETNEO float SumProdNEOx(float* A, float* B, int64 n)
{
#define N 4
	int64 i = 0;
	volatile float re = 0;

	if (n >= N * E512F)
	{
		float32x4_t s[E512_128] = { 0 };

		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) UNROLL(E512_128) 
			{ s[kk] = _neo_fmaddx_ps(vld1q_f32(A), vld1q_f32(B), s[kk]); A += E128F; B += E128F; }
		}

		REDUCE(s) s[kk] = vaddq_f32(s[kk], s[kk + KK]);

		re = _neo_reduce_add_ps(s[0]);
	}

	for (; i < n; ++i, ++A, ++B)
	{
		volatile float v1 = (float)*A * (float)*B;
		re += v1;
	}

	return re;
}

TARGETNEO double SumProdNEO(double* A, double* B, double* C, int64 n)
{
#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		float64x2_t s[E512_128] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128) 
			{ s[kk] = _neo_fmaddx_pd(vmulq_f64(vld1q_f64(A), vld1q_f64(B)), vld1q_f64(C), s[kk]); A += E128D; B += E128D; C += E128D; }
		}

		REDUCE(s) s[kk] = vaddq_f64(s[kk], s[kk + KK]);

		re = _neo_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, ++A, ++B, ++C)
	{
		volatile double v1 = (double)*A * (double)*B * (double)*C;
		re += v1;
	}

	return re;
}

TARGETNEO float SumProdNEO(float* A, float* B, float* C, int64 n)
{
#define N 4
	int64 i = 0;
	volatile float re = 0;

	if (n >= N * E512F)
	{
		float32x4_t s[E512_128] = { 0 };

		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) UNROLL(E512_128)
			{ s[kk] = _neo_fmaddx_ps(vmulq_f32(vld1q_f32(A), vld1q_f32(B)), vld1q_f32(C), s[kk]); A += E128F; B += E128F; C += E128F; }
		}

		REDUCE(s) s[kk] = vaddq_f32(s[kk], s[kk + KK]);

		re = _neo_reduce_add_ps(s[0]);
	}

	for (; i < n; ++i, ++A, ++B, ++C)
	{
		volatile float v1 = (float)*A * (float)*B * (float)*C;
		re += v1;
	}

	return re;
}

TARGETNEO double SumSqProdNEO(double* A, double* B, int64 n)
{
#define N 4
	int64 i = 0;
	volatile double re = 0;

	if (n >= N * E512D)
	{
		float64x2_t s[E512_128] = { 0 };

		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128)
			{ s[kk] = _neo_fmaddx_pd(vmulq_f64(vld1q_f64(A), vld1q_f64(A)), vld1q_f64(B), s[kk]); A += E128D; B += E128D; }
		}

		REDUCE(s) s[kk] = vaddq_f64(s[kk], s[kk + KK]);

		re = _neo_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, ++A, ++B)
	{
		volatile double v1 = (double)*A * (double)*A * (double)*B;
		re += v1;
	}

	return re;
}

TARGETNEO float SumSqProdNEO(float* A, float* B, int64 n)
{
#define N 4
	int64 i = 0;
	volatile float re = 0;

	if (n >= N * E512F)
	{
		float32x4_t s[E512_128] = { 0 };

		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) UNROLL(E512_128) 
			{ s[kk] = _neo_fmaddx_ps(vmulq_f32(vld1q_f32(A), vld1q_f32(A)), vld1q_f32(B), s[kk]); A += E128F; B += E128F; }
		}

		REDUCE(s) s[kk] = vaddq_f32(s[kk], s[kk + KK]);

		re = _neo_reduce_add_ps(s[0]);
	}

	for (; i < n; ++i, ++A, ++B)
	{
		volatile float v1 = (float)*A * (float)*A * (float)*B;
		re += v1;
	}

	return re;
}

TARGETNEO void AddNEO(double* A, double* B, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E512D)
	{
		for (int64 l1 = n - N * E512D; i <= l1; i += N * E512D)
		{
			UNROLL(N) UNROLL(E512_128) 
			{ vst1q_f64(A, vaddq_f64(vld1q_f64(A), vld1q_f64(B))); A += E128D; B += E128D; }
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGETNEO void AddNEO(float* A, float* B, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E512F)
	{
		for (int64 l1 = n - N * E512F; i <= l1; i += N * E512F)
		{
			UNROLL(N) UNROLL(E512_128) 
			{ vst1q_f32(A, vaddq_f32(vld1q_f32(A), vld1q_f32(B))); A += E128F; B += E128F; }
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGETNEO void AddNEO(int64* A, int64* B, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E128D)
	{
		for (int64 l1 = n - N * E128D; i <= l1; i += N * E128D)
		{
			UNROLL(N) { vst1q_s64(A, vaddq_s64(vld1q_s64(A), vld1q_s64(B))); A += E128D; B += E128D; }
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGETNEO void AddNEO(int* A, int* B, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E128F)
	{
		for (int64 l1 = n - N * E128F; i <= l1; i += N * E128F)
		{
			UNROLL(N) { vst1q_s32(A, vaddq_s32(vld1q_s32(A), vld1q_s32(B)));  A += E128F; B += E128F; }
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGETNEO void AddNEO(int* A, int B, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E128F)
	{
		int32x4_t b = b = vdupq_n_s32(B);

		for (int64 l1 = n - N * E128F; i <= l1; i += N * E128F)
		{
			UNROLL(N) { vst1q_s32(A, vaddq_s32(vld1q_s32(A), b)); A += E128F; }
		}
	}

	for (; i < n; ++i, A++)
		*A += B;
}

TARGETNEO void AddNEO(double* A, double B, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E128D)
	{
		float64x2_t b = vdupq_n_f64(B);

		for (int64 l1 = n - N * E128D; i <= l1; i += N * E128D)
		{
			UNROLL(N) { vst1q_f64(A, vaddq_f64(vld1q_f64(A), b)); A += E128D; }
		}
	}

	for (; i < n; ++i, A++)
		*A += B;
}

TARGETNEO void AddNEO(float* A, float B, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E128F)
	{
		float32x4_t b = vdupq_n_f32(B);

		for (int64 l1 = n - N * E128F; i <= l1; i += N * E128F)
		{
			UNROLL(N) { vst1q_f32(A, vaddq_f32(vld1q_f32(A), b)); A += E128F; }
		}
	}

	for (; i < n; ++i, A++)
		*A += B;
}

TARGETNEO void MulNEO(double* A, double* B, double* C, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E128D)
	{
		for (int64 l1 = n - N * E128D; i <= l1; i += N * E128D)
		{
			UNROLL(N) { vst1q_f64(A, vmulq_f64(vld1q_f64(B), vld1q_f64(C))); A += E128D; B += E128D; C += E128D; }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *B++ * *C++;
		*A++ = v1;
	}
}

TARGETNEO void MulNEO(float* A, float* B, float* C, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E128F)
	{
		for (int64 l1 = n - N * E128F; i <= l1; i += N * E128F)
		{
			UNROLL(N) { vst1q_f32(A, vmulq_f32(vld1q_f32(B), vld1q_f32(C))); A += E128F; B += E128F; C += E128F; }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *B++ * *C++;
		*A++ = v1;
	}
}

TARGETNEO void MulNEO(double* A, double* B, double C, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E128D)
	{
		float64x2_t c = vdupq_n_f64(C);

		for (int64 l1 = n - N * E128D; i <= l1; i += N * E128D)
		{
			UNROLL(N) { vst1q_f64(A, vmulq_f64(vld1q_f64(B), c)); A += E128D; B += E128D; }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *B++ * C;
		*A++ = v1;
	}
}

TARGETNEO void MulNEO(float* A, float* B, float C, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E128F)
	{
		float32x4_t c = vdupq_n_f32(C);

		for (int64 l1 = n - N * E128F; i <= l1; i += N * E128F)
		{
			UNROLL(N) { vst1q_f32(A, vmulq_f32(vld1q_f32(B), c)); A += E128F; B += E128F;}
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *B++ * C;
		*A++ = v1;
	}
}

TARGETNEO void MulNEO(double* A, double B, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E128D)
	{
		float64x2_t b = vdupq_n_f64(B);

		for (int64 l1 = n - N * E128D; i <= l1; i += N * E128D)
		{
			UNROLL(N) { vst1q_f64(A, vmulq_f64(vld1q_f64(A), b)); A += E128D; }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *A * B;
		*A++ = v1;
	}
}

TARGETNEO void MulNEO(float* A, float B, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E128F)
	{
		float32x4_t b = vdupq_n_f32(B);

		for (int64 l1 = n - N * E128F; i <= l1; i += N * E128F)
		{
			UNROLL(N) { vst1q_f32(A, vmulq_f32(vld1q_f32(A), b)); A += E128F; }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A * B;
		*A++ = v1;
	}
}

TARGETNEO void DivNEO(double* A, double B, double* C, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E128D)
	{
		float64x2_t b = vdupq_n_f64(B);

		for (int64 l1 = n - N * E128D; i <= l1; i += N * E128D)
		{
			UNROLL(N) { vst1q_f64(A, vdivq_f64(b, vld1q_f64(C))); A += E128D; C += E128D; }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = B / *C++;
		*A++ = v1;
	}
}

TARGETNEO void DivNEO(float* A, float B, float* C, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E128F)
	{
		float32x4_t b = vdupq_n_f32(B);

		for (int64 l1 = n - N * E128F; i <= l1; i += N * E128F)
		{
			UNROLL(N) { vst1q_f32(A, vdivq_f32(b, vld1q_f32(C))); A += E128F; C += E128F; }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = B / *C++;
		*A++ = v1;
	}
}

TARGETNEO void DivNEO(double* A, double* B, double* C, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E128D)
	{
		for (int64 l1 = n - N * E128D; i <= l1; i += N * E128D)
		{
			UNROLL(N) { vst1q_f64(A, vdivq_f64(vld1q_f64(B), vld1q_f64(C))); A += E128D; B += E128D; C += E128D; }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *B++ / *C++;
		*A++ = v1;
	}
}

TARGETNEO void DivNEO(float* A, float* B, float* C, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E128F)
	{
		for (int64 l1 = n - N * E128F; i <= l1; i += N * E128F)
		{
			UNROLL(N) { vst1q_f32(A, vdivq_f32(vld1q_f32(B), vld1q_f32(C))); A += E128F; B += E128F; C += E128F; }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *B++ / *C++;
		*A++ = v1;
	}
}

TARGETNEO void AddProdNEO(double* A, double* B, double* C, int64 n)
{
#define N 1
	int64 i = 0;

	if (n >= N * E128D)
	{
		for (int64 l1 = n - N * E128D; i <= l1; i += N * E128D)
		{
			UNROLL(N) { vst1q_f64(A, _neo_fmaddx_pd(vld1q_f64(B), vld1q_f64(C), vld1q_f64(A))); A += E128D; B += E128D; C += E128D; }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *B++ * *C++;
		*A++ += v1;
	}
}

TARGETNEO void AddProdNEO(float* A, float* B, float* C, int64 n)
{
#define N 1
	int64 i = 0;

	if (n >= N * E128F)
	{
		for (int64 l1 = n - N * E128F; i <= l1; i += N * E128F)
		{
			UNROLL(N) { vst1q_f32(A, _neo_fmaddx_ps(vld1q_f32(B), vld1q_f32(C), vld1q_f32(A))); A += E128F; B += E128F; C += E128F; }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *B++ * *C++;
		*A++ += v1;
	}
}

TARGETNEO void AddProdNEO(double* A, double* B, double C, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E128D)
	{
		float64x2_t b[N], c = vdupq_n_f64(C);

		for (int64 l1 = n - N * E128D; i <= l1; i += N * E128D)
		{
			UNROLL(N) { b[kk] = vld1q_f64(B); B += E128D; }

			UNROLL(N)  b[kk] = vmulq_f64(b[kk], c);

			double* A2 = A;

			UNROLL(N) { b[kk] = vaddq_f64(b[kk], vld1q_f64(A)); A += E128D; }

			UNROLL(N) { vst1q_f64(A2, b[kk]); A2 += E128D; }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *B++ * C;
		*A++ += v1;
	}
}

TARGETNEO void AddProdNEO(double* A, float* B, double C, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E128D)
	{
		float64x2_t b[N], c = vdupq_n_f64(C);

		for (int64 l1 = n - N * E128D; i <= l1; i += N * E128D)
		{
			float32x4_t v1 = vld1q_f32(B); B += E128F;
			b[0] = vcvt_f64_f32(vget_low_f32(v1));
			b[1] = vcvt_f64_f32(vget_high_f32(v1));

			float32x4_t v2 = vld1q_f32(B); B += E128F;
			b[2] = vcvt_f64_f32(vget_low_f32(v2));
			b[3] = vcvt_f64_f32(vget_high_f32(v2));

			UNROLL(N) b[kk] = vmulq_f64(b[kk], c);

			double* A2 = A;

			UNROLL(N) { b[kk] = vaddq_f64(b[kk], vld1q_f64(A)); A += E128D; }

			UNROLL(N) { vst1q_f64(A2, b[kk]); A2 += E128D; }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *B++ * C;
		*A++ += v1;
	}
}

TARGETNEO void AddProdNEO(float* A, float* B, float C, int64 n)
{
#define N 4
	int64 i = 0;

	if (n >= N * E128F)
	{
		float32x4_t b[N], c = vdupq_n_f32(C);

		for (int64 l1 = n - N * E128F; i <= l1; i += N * E128F)
		{
			UNROLL(N) { b[kk] = vld1q_f32(B); B += E128F; }

			UNROLL(N)  b[kk] = vmulq_f32(b[kk], c);

			float* A2 = A;

			UNROLL(N) { b[kk] = vaddq_f32(b[kk], vld1q_f32(A)); A += E128F; }

			UNROLL(N) { vst1q_f32(A2, b[kk]); A2 += E128F; }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *B++ * C;
		*A++ += v1;
	}
}

TARGETNEO void UnifyNEO(double* A, int64 n)
{
#define N 4
	int64 i = 0;
	double invsum = 1.0 / (SumNEO(A, n) + n * MIN_FREQ);

	if (n >= E128D)
	{
		float64x2_t b = vdupq_n_f64(invsum), c = vdupq_n_f64(MIN_FREQ * invsum);
		
		UNROLLHEAD(N)
		for (int64 l1 = n - E128D; i <= l1; i += E128D)
		{
			vst1q_f64(A, _neo_fmaddx_pd(vld1q_f64(A), b, c)); 
			A += E128D;
		}
	}

	for (; i < n; ++i, ++A)
	{
		volatile double v1 = (double)*A + MIN_FREQ;
		*A = v1 * invsum;
	}
}

TARGETNEO void UnifyNEO(float* A, int64 n)
{
#define N 4
	int64 i = 0;
	double invsum = 1.0 / (SumNEO(A, n) + n * MIN_FREQ);

	if (n >= E128F)
	{
		float64x2_t b = vdupq_n_f64(invsum), c = vdupq_n_f64(MIN_FREQ * invsum);
		
		UNROLLHEAD(N)
		for (int64 l1 = n - E128F; i <= l1; i += E128F)
		{
			vst1q_f32(A, vcombine_f32(
					vcvt_f32_f64(_neo_fmaddx_pd(vcvt_f64_f32(vld1_f32(A)), b, c)), 
					vcvt_f32_f64(_neo_fmaddx_pd(vcvt_f64_f32(vld1_f32(A + E128D)), b, c))));
			A += E128F;
		}
	}

	for (; i < n; ++i, ++A)
	{
		volatile double v1 = (double)*A + MIN_FREQ;
		*A = v1 * invsum;
	}
}

TARGETNEO char* StrNextIdxNEO(char* A, char val, int64 rep, int64 n)
{
#define N 4
	A++; n--;
	int64 i = 0;

	if (n >= E512B)
	{
		uint8x16_t v = vdupq_n_u8(val), o = vdupq_n_u8(0);

		for (int64 l1 = n - E512B; i <= l1; i += E512B)
		{
			char* Ab = A;
			int count = 0;

			UNROLL(N) { count += vaddvq_u8(vsubq_u8(o,vceqq_u8(vld1q_u8((byte*)A), v))); A += E128B; }

			if (rep > count) [[likely]]
			{
				rep -= count;
				continue;
			}
			else for (;; Ab++)
				if (*Ab == val && !--rep)
					return Ab;
		}
	}

	for (; i < n; ++i, A++)
		if (*A == val && !--rep)
			return A;

	return NULL;
}

TARGETNEO int64 CountCharNEO(char* A, char val, int64 n)
{
#define N 4
	uint64 re = 0;
	int64 i = 0;

	if (n >= E128B)
	{
		uint8x16_t v = vdupq_n_u8(val), o = vdupq_n_u8(0);
		
		UNROLLHEAD(N)
		for (int64 l1 = n - E128B; i <= l1; i += E128B)
		{
			re += vaddvq_u8(vsubq_u8(o, vceqq_u8(vld1q_u8((byte*)A), v)));
			A += E128B;
		}
	}

	for (; i < n; ++i, A++)
		if (*A == val) re++;

	return (int64)re;
}

TARGETNEO void DiagQuadFormNEO(double* res1, double* A, double* D, int64 m, int64 n)
{
#define N 3

#define DECLARE			double* pA1 = (A + n * i), *pA2 = (A + n * j), *pD = D; float64x2_t a1[N], a2[N], r[N][N] = { 0 }
#define ALOAD1(ii)		a1[ii] = vld1q_f64(pA1 + n * ii)
#define ALOAD2(ii)		a2[ii] = vld1q_f64(pA2 + n * ii)
#define ALOAD3(ii)		a1[ii] = a2[ii] = vld1q_f64(pA1 + n * ii)
#define ADMUL(ii)		a1[ii] = vmulq_f64(a1[ii], vld1q_f64(pD))
#define FMADD(ii,jj)	r[ii][jj] = _neo_fmaddx_pd(a1[ii], a2[jj], r[ii][jj])
#define RDUADD(ii,jj)	res1[(i+ii) * m + (j+jj)] = _neo_reduce_add_pd(r[ii][jj])
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

TARGETNEO void DiagQuadFormNEO(double* res2, double* A, double* D, double* B, int64 m, int64 n)
{
#define N 8

#define DECLARE			double* pA = (A + n * i), *pB = B, *pD = D; float64x2_t bd, r[N] = { 0 }
#define BDMUL			bd = vmulq_f64(vld1q_f64(pB), vld1q_f64(pD)); pB += E128D; pD += E128D
#define FMADD(ii)		r[ii] = _neo_fmaddx_pd(vld1q_f64(pA + n * ii), bd, r[ii]); 
#define RDUADD(ii)		res2[(i+ii)] = _neo_reduce_add_pd(r[ii])
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

TARGETNEO void DiagQuadFormNEO(double* res3, double* B, double* D, int64 n)
{
#define N 4
    double re = 0;
	float64x2_t s[N] = { 0 };

	int64 k = 0;
    for (; k + N * E128D <= n; k += N * E128D)
	{
		UNROLL(N) 
		{
			s[kk] = _neo_fmaddx_pd(vmulq_f64(vld1q_f64(B), vld1q_f64(B)), vld1q_f64(D), s[kk]);
			B += E128D;
			D += E128D;
		}
	}

	REDUCE(s) s[kk] = vaddq_f64(s[kk], s[kk + KK]);
    
    res3[0] = _neo_reduce_add_pd(s[0]);
            
    VECTORIZE
    for (; k < n; k++, B++, D++)
        res3[0] += B[0] * B[0] * D[0];
}

TARGETNEO void MatrixMulNEO(double* res, double* A, double* B, int64 m, int64 n, int64 p)
{
    //suboptimal to compile
    {
        Mat<double> a(A, n, m, false, true);
        Mat<double> b(B, n, p, false, true);
        Mat<double> r(res, m, p, false, true);

        r = a.t() * b;
        return;
    }
    
#define N 3
	
#define DECLARE			double* pA = (A + n * i), *pB = (B + n * j); float64x2_t a[N], b[N], r[N][N] = { 0 }
#define ALOAD(kk)		a[kk] = vld1q_f64(pA + n * kk)
#define BLOAD(kk)		b[kk] = vld1q_f64(pB + n * kk)
#define FMADD(ii,jj)	r[ii][jj] = _neo_fmaddx_pd(a[ii], b[jj], r[ii][jj])
#define RDUADD(ii,jj)	res[(i+ii) + (j+jj) * m] = _neo_reduce_add_pd(r[ii][jj])
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

TARGETNEO void DiagQuadFormNEO(float* res1, float* A, float* D, int64 m, int64 n)
{
#define N 3

#define DECLARE			float* pA1 = (A + n * i), *pA2 = (A + n * j), *pD = D; float32x4_t a1[N], a2[N], r[N][N] = { 0 }
#define ALOAD1(ii)		a1[ii] = vld1q_f32(pA1 + n * ii)
#define ALOAD2(ii)		a2[ii] = vld1q_f32(pA2 + n * ii)
#define ALOAD3(ii)		a1[ii] = a2[ii] = vld1q_f32(pA1 + n * ii)
#define ADMUL(ii)		a1[ii] = vmulq_f32(a1[ii], vld1q_f32(pD))
#define FMADD(ii,jj)	r[ii][jj] = _neo_fmaddx_ps(a1[ii], a2[jj], r[ii][jj])
#define RDUADD(ii,jj)	res1[(i+ii) * m + (j+jj)] = _neo_reduce_add_ps(r[ii][jj])
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

TARGETNEO void DiagQuadFormNEO(float* res2, float* A, float* D, float* B, int64 m, int64 n)
{
#define N 8

#define DECLARE			float* pA = (A + n * i), *pB = B, *pD = D; float32x4_t bd, r[N] = { 0 }
#define BDMUL			bd = vmulq_f32(vld1q_f32(pB), vld1q_f32(pD)); pB += E128F; pD += E128F
#define FMADD(ii)		r[ii] = _neo_fmaddx_ps(vld1q_f32(pA + n * ii), bd, r[ii]); 
#define RDUADD(ii)		res2[(i+ii)] = _neo_reduce_add_ps(r[ii])
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

TARGETNEO void DiagQuadFormNEO(float* res3, float* B, float* D, int64 n)
{
#define N 4
	float32x4_t s[N] = { 0 };

	int64 k = 0;
    for (; k + N * E128F <= n; k += N * E128F) 
	{
		UNROLL(N) 
		{
			s[kk] = _neo_fmaddx_ps(vmulq_f32(vld1q_f32(B), vld1q_f32(B)), vld1q_f32(D), s[kk]); 
			B += E128F;
			D += E128F;
		}
	}

	REDUCE(s) s[kk] = vaddq_f32(s[kk], s[kk + KK]);
    
    res3[0] = _neo_reduce_add_ps(s[0]);
            
    VECTORIZE
    for (; k < n; k++, B++, D++)
        res3[0] += B[0] * B[0] * D[0];
}

TARGETNEO void MatrixMulNEO(float* res, float* A, float* B, int64 m, int64 n, int64 p)
{
    //suboptimal to compile
    {
        Mat<float> a(A, n, m, false, true);
        Mat<float> b(B, n, p, false, true);
        Mat<float> r(res, m, p, false, true);

        r = a.t() * b;
        return;
    }
    
#define N 3
	
#define DECLARE			float* pA = (A + n * i), *pB = (B + n * j); float32x4_t a[N], b[N], r[N][N] = { 0 }
#define ALOAD(kk)		a[kk] = vld1q_f32(pA + n * kk)
#define BLOAD(kk)		b[kk] = vld1q_f32(pB + n * kk)
#define FMADD(ii,jj)	r[ii][jj] = _neo_fmaddx_ps(a[ii], b[jj], r[ii][jj])
#define RDUADD(ii,jj)	res[(i+ii) + (j+jj) * m] = _neo_reduce_add_ps(r[ii][jj])
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
