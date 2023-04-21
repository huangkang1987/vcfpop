/* NEON Instruction Set Functions */

#include "vcfpop.h"

#ifdef __aarch64__

template TARGETNEO void RNGNEO<double>::Poly<32>(float64x2_t* arr, int n, void* re);
template TARGETNEO void RNGNEO<double>::Poly<64>(float64x2_t* arr, int n, void* re);
template TARGETNEO void RNGNEO<float >::Poly<32>(float32x4_t* arr, int n, void* re);
template TARGETNEO void RNGNEO<float >::Poly<64>(float32x4_t* arr, int n, void* re);

#ifndef _RNGNEO_FP64
/* Initialize rng */
TARGETNEO RNGNEO<double>::RNGNEO()
{

}

/* Initialize rng */
TARGETNEO RNGNEO<double>::RNGNEO(uint64 seed, uint64 salt)
{
	uint64x2_t a[32], s, m;

	REP(32) { a[kk] = vld1q_u64(((uint64[]) { seed, seed + 1 })); seed += 2; }

	s = vdupq_n_u64(salt);
	m = vdupq_n_u32(0x5bd1e995);

	REP(32) a[kk] = veorq_u64(a[kk], vshlq_n_u64(vmvnq_u32(a[kk]), 32));

	s             = veorq_u64(s    , vshlq_n_u64(vmvnq_u32(s    ), 32));

	// uint s = s ^ sizeof(uint);
	s = veorq_u32(s, vdupq_n_u32(sizeof(uint)));

	// a *= m;
	REP(32) a[kk] = vmulq_u32(a[kk], m);

	// a ^= a >> 24;
	REP(32) a[kk] = veorq_u64(a[kk], vshrq_n_u32(a[kk], 24));

	// a *= m;
	REP(32) a[kk] = vmulq_u32(a[kk], m);

	// s *= m;
	s = vmulq_u32(s, m);

	// a ^= s;
	REP(32) a[kk] = veorq_u64(a[kk], s);

	// a ^= a >> 13;
	REP(32) a[kk] = veorq_u64(a[kk], vshrq_n_u32(a[kk], 13));

	// a *= m;
	REP(32) a[kk] = vmulq_u32(a[kk], m);

	// a ^= a >> 15;
	REP(32) a[kk] = veorq_u64(a[kk], vshrq_n_u32(a[kk], 15));

	// original
	REP(32) x[kk] = veorq_u64(vdupq_n_u64(0x159A55E5075BCD15), a[kk]);

	REP(32) a[kk] = vshlq_n_u64(a[kk], 6);

	REP(32) y[kk] = veorq_u64(vdupq_n_u64(0x054913331F123BB5), a[kk]);
}

/* Draw a uniform distriubted real number */
template<int nbits>
TARGETNEO void RNGNEO<double>::Poly(float64x2_t* arr, int n, void* re)
{
	float64x2_t t[32], s[32];
	float64x2_t one = vdupq_n_f64(1.0);
	uint64x2_t mask1 = vdupq_n_u64(0x000FFFFFFFFFFFFF);
	uint64x2_t mask2 = vdupq_n_u64(0x3FF0000000000000);
	uint64x2_t* r = (uint64x2_t*)t; uint64x2_t* re64 = (uint64x2_t*)re; uint32x4_t* re32 = (uint32x4_t*)re;

	REP(32) s[kk] = vdupq_n_f64(0);

	for (int i32 = 0; i32 < n * 32; i32 += 32)
		REP(32) s[kk] = vaddq_f64(s[kk], arr[kk + i32]);

	{
		uint64x2_t a[32], b[32];

		REP(32) a[kk] = x[kk];

		REP(32) b[kk] = y[kk];

		REP(32) x[kk] = b[kk];

		REP(32) a[kk] = veorq_u64(a[kk], vshlq_n_u64(a[kk], 23));

		REP(32) a[kk] = veorq_u64(a[kk], vshrq_n_u64(a[kk], 18));

		REP(32) a[kk] = veorq_u64(a[kk], b[kk]);

		REP(32) a[kk] = veorq_u64(a[kk], vshrq_n_u64(b[kk], 5));

		REP(32) y[kk] = a[kk];

		REP(32) r[kk] = vaddq_s64(a[kk], b[kk]);
	}

	REP(32) r[kk] = vandq_u64(r[kk], mask1);
	
	REP(32) r[kk] = vorrq_u64(r[kk], mask2);

	REP(32) t[kk] = vsubq_f64(t[kk], one);

	REP(32) t[kk] = vmulq_f64(t[kk], s[kk]);

	uint64x2_t midx[32], nidx = vdupq_n_u64(0), ninc = vdupq_n_u64(1);
	uint64x2_t f[32], b[32];
	REP(32) midx[kk] = vdupq_n_u64(n - 1);
	REP(32) f[kk] = vdupq_n_f64(0);

	for (int i32 = 0; i32 < n * 32; i32 += 32)
	{
		REP(32) b[kk] = vcltq_f64(t[kk], arr[kk + i32]);

		REP(8) t[ 0 + kk] = vsubq_f64(t[ 0 + kk], arr[ 0 + kk + i32]);

		REP(32) b[kk] = vbicq_u64(b[kk], f[kk]);

		REP(8) t[ 8 + kk] = vsubq_f64(t[ 8 + kk], arr[ 8 + kk + i32]);

		REP(32) f[kk] = vorrq_u64(f[kk], b[kk]);

		REP(8) t[16 + kk] = vsubq_f64(t[16 + kk], arr[16 + kk + i32]);

		REP(32) midx[kk] = vbslq_s64(b[kk], nidx, midx[kk]);

		REP(8) t[24 + kk] = vsubq_f64(t[24 + kk], arr[24 + kk + i32]);

		nidx = vaddq_s64(nidx, ninc);
	}

	if constexpr (nbits == 32)
	{
		REP(16) re32[kk] = vld1q_u32(((uint[]) { (uint)vgetq_lane_u64(midx[0 + (kk << 1)], 0), (uint)vgetq_lane_u64(midx[0 + (kk << 1)], 1), (uint)vgetq_lane_u64(midx[1 + (kk << 1)], 0), (uint)vgetq_lane_u64(midx[1 + (kk << 1)], 1) }));
	}
	else
	{
		REP(32) re64[kk] = midx[kk];
	}
}
#endif

#ifndef _RNGNEO_FP32
/* Initialize rng */
TARGETNEO RNGNEO<float>::RNGNEO()
{

}

/* Initialize rng */
TARGETNEO RNGNEO<float>::RNGNEO(uint64 seed, uint64 salt)
{
	uint32x4_t a[16], s, m;
	REP(16) { a[kk] = vld1q_u32(((uint[]) { Mix(seed + 0), Mix(seed + 1), Mix(seed + 2), Mix(seed + 3) })); seed += 4; }

	s = vdupq_n_u32(Mix(salt));
	m = vdupq_n_u32(0x5bd1e995);

	// uint s = s ^ sizeof(uint);
	s = veorq_u32(s, vdupq_n_u32(sizeof(uint)));

	// a *= m;
	REP(16) a[kk] = vmulq_u32(a[kk], m);

	// a ^= a >> 24;
	REP(16) a[kk] = veorq_u32(a[kk], vshrq_n_u32(a[kk], 24));

	// a *= m;
	REP(16) a[kk] = vmulq_u32(a[kk], m);

	// s *= m;
	s = vmulq_u32(s, m);

	// a ^= s;
	REP(16) a[kk] = veorq_u32(a[kk], s);

	// a ^= a >> 13;
	REP(16) a[kk] = veorq_u32(a[kk], vshrq_n_u32(a[kk], 13));

	// a *= m;
	REP(16) a[kk] = vmulq_u32(a[kk], m);

	// a ^= a >> 15;
	REP(16) a[kk] = veorq_u32(a[kk], vshrq_n_u32(a[kk], 15));

	// original
	REP(16) x[kk] = veorq_u32(vdupq_n_u32(0x075BCD15), a[kk]);

	REP(16) a[kk] = vshlq_n_u32(a[kk], 3);

	REP(16) y[kk] = veorq_u32(vdupq_n_u32(0x159A55E5), a[kk]);

	REP(16) a[kk] = vshlq_n_u32(a[kk], 3);

	REP(16) z[kk] = veorq_u32(vdupq_n_u32(0x1F123BB5), a[kk]);
}

/* Draw a uniform distriubted real number */
template<int nbits>
TARGETNEO void RNGNEO<float>::Poly(float32x4_t* arr, int n, void* re)
{
	float32x4_t t[16], s[16]; uint32x4_t u[16];
	float32x4_t one = vdupq_n_f32(1.0f);
	uint32x4_t mask1 = vdupq_n_u32(0x007FFFFF);
	uint32x4_t mask2 = vdupq_n_u32(0x3F800000);
	uint32x4_t* r = (uint32x4_t*)t; uint64x2_t* re64 = (uint64x2_t*)re; uint32x4_t* re32 = (uint32x4_t*)re;

	REP(16) s[kk] = vdupq_n_f32(0);

	for (int i16 = 0; i16 < n * 16; i16 += 16)
		REP(16) s[kk] = vaddq_f32(s[kk], arr[kk + i16]);

	{
		//xorshit 
		REP(16) u[kk] = vshlq_n_u32(x[kk], 16);
		REP(16) x[kk] = veorq_u32(x[kk], u[kk]);

		REP(16) u[kk] = vshrq_n_u32(x[kk], 5);
		REP(16) x[kk] = veorq_u32(x[kk], u[kk]);

		REP(16) u[kk] = vshlq_n_u32(x[kk], 1);
		REP(16) x[kk] = veorq_u32(x[kk], u[kk]);

		REP(16) u[kk] = x[kk];

		REP(16) x[kk] = y[kk];

		REP(16) y[kk] = z[kk];

		REP(16) z[kk] = veorq_u32(u[kk], x[kk]);

		REP(16) z[kk] = veorq_u32(z[kk], y[kk]);
	}

	REP(16) r[kk] = vandq_u32(z[kk], mask1);

	REP(16) r[kk] = vorrq_u32(r[kk], mask2);

	REP(16) t[kk] = vsubq_f32(t[kk], one);

	REP(16) t[kk] = vmulq_f32(t[kk], s[kk]);

	uint32x4_t midx[16], nidx = vdupq_n_u32(0), ninc = vdupq_n_u32(1);
	uint32x4_t f[16], b[16];
	REP(16) midx[kk] = vdupq_n_u32(n - 1);
	REP(16) f[kk] = vdupq_n_u32(0);

	for (int i16 = 0; i16 < n * 16; i16 += 16)
	{
		REP(16) b[kk] = vcltq_f32(t[kk], arr[kk + i16]);

		REP(8) t[0 + kk] = vsubq_f32(t[0 + kk], arr[0 + kk + i16]);

		REP(16) b[kk] = vbicq_u32(b[kk], f[kk]);

		REP(8) t[8 + kk] = vsubq_f32(t[8 + kk], arr[8 + kk + i16]);

		REP(16)
		{
			f[kk] = vorrq_u32(f[kk], b[kk]);
			midx[kk] = vbslq_s32(b[kk], nidx, midx[kk]);
		}

		nidx = vaddq_s32(nidx, ninc);
	}

	if constexpr (nbits == 32)
	{
		REP(16) re32[kk] = midx[kk];
	}
	else
	{
		REP(16) 
		{
			re64[0 + (kk << 1)] = vld1q_u64(((uint64[]) { vgetq_lane_u32(midx[kk], 0), vgetq_lane_u32(midx[kk], 1) }));
			re64[1 + (kk << 1)] = vld1q_u64(((uint64[]) { vgetq_lane_u32(midx[kk], 2), vgetq_lane_u32(midx[kk], 3) }));
		}
	}
}
#endif

__forceinline TARGETNEO double _neo_reduce_add_pd(float64x2_t v2)
{
	return vgetq_lane_f64(v2, 0) + vgetq_lane_f64(v2, 1);
}

__forceinline TARGETNEO double _neo_reduce_mul_pd(float64x2_t v2)
{
	return vgetq_lane_f64(v2, 0) * vgetq_lane_f64(v2, 1);
}

__forceinline TARGETNEO float _neo_reduce_add_ps(float32x4_t v2)
{
	volatile float a1 = vgetq_lane_f32(v2, 0) + vgetq_lane_f32(v2, 2), a2 = vgetq_lane_f32(v2, 1) + vgetq_lane_f32(v2, 3);
	return a1 + a2;
}

__forceinline TARGETNEO float _neo_reduce_mul_ps(float32x4_t v2)
{
	volatile float a1 = vgetq_lane_f32(v2, 0) * vgetq_lane_f32(v2, 2), a2 = vgetq_lane_f32(v2, 1) * vgetq_lane_f32(v2, 3);
	return a1 * a2;
}

__forceinline TARGETNEO double _neo_reduce_add_psd(float32x4_t v2)
{
	volatile double a1 = (double)vgetq_lane_f32(v2, 0) + (double)vgetq_lane_f32(v2, 2);
	volatile double a2 = (double)vgetq_lane_f32(v2, 1) + (double)vgetq_lane_f32(v2, 3);
	return a1 + a2;
}

__forceinline TARGETNEO double _neo_reduce_mul_psd(float32x4_t v2)
{
	volatile double a1 = (double)vgetq_lane_f32(v2, 0) * (double)vgetq_lane_f32(v2, 2);
	volatile double a2 = (double)vgetq_lane_f32(v2, 1) * (double)vgetq_lane_f32(v2, 3);
	return a1 * a2;
}

TARGETNEO int64 GetMinIdxNEO(double* A, int64 n, double& val)
{
	constexpr int N = 4;
	int64 i = 0;
	val = DBL_MAX;
	uint64 idx = (uint64)-1;

	if (n >= N * sizeof(float64x2_t) / sizeof(double))
	{
		float64x2_t min1[N], f[N], a[N];
		uint64x2_t midx[N], nidx[N], msep = vdupq_n_u64(8);
		REP(N) min1[kk] = vdupq_n_f64(val);
		REP(N) midx[kk] = vdupq_n_u64(0xFFFFFFFFFFFFFFFF);
		REP(N) nidx[kk] = vld1q_u64(((uint64[]) { 0ull + (kk << 1), 1ull + (kk << 1) }));

		for (int64 l1 = n - N * sizeof(float64x2_t) / sizeof(double); i <= l1; i += N * sizeof(float64x2_t) / sizeof(double))
		{
			REP(N) { a[kk] = vld1q_f64(A); A += sizeof(float64x2_t) / sizeof(double); }

			REP(N) f[kk] = vcgtq_f64(min1[kk], a[kk]);

			REP(N) min1[kk] = vbslq_s64(f[kk], a[kk], min1[kk]);

			REP(N) midx[kk] = vbslq_s64(f[kk], nidx[kk], midx[kk]);

			REP(N) nidx[kk] = vaddq_s64(nidx[kk], msep);
		}

		for (int KK = sizeof(min1) / sizeof(min1[0]) / 2; KK >= 1; KK >>= 1)
		{
			REP(KK) f[kk] = vcgtq_f64(min1[kk], min1[kk + KK]);
			REP(KK) min1[kk] = vbslq_s64(f[kk], min1[kk + KK], min1[kk]);
			REP(KK) midx[kk] = vbslq_s64(f[kk], midx[kk + KK], midx[kk]);
		}

		if (vgetq_lane_f64(min1[0], 0) < val)
		{
			val = vgetq_lane_f64(min1[0], 0);
			idx = vgetq_lane_u64(midx[0], 0);
		}

		if (vgetq_lane_f64(min1[0], 1) < val)
		{
			val = vgetq_lane_f64(min1[0], 1);
			idx = vgetq_lane_u64(midx[0], 1);
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
	constexpr int N = 4;
	int64 i = 0;
	val = FLT_MAX;
	uint idx = (uint)-1;

	if (n >= N * sizeof(float32x4_t) / sizeof(float))
	{
		float32x4_t min1[N], f[N], a[N];
		uint32x4_t midx[N], nidx[N], msep = vdupq_n_u32(N * sizeof(float32x4_t) / sizeof(float));
		REP(N) min1[kk] = vdupq_n_f32(val);
		REP(N) midx[kk] = vdupq_n_u32(0xFFFFFFFF);
		REP(N) nidx[kk] = vld1q_u32(((uint[]) { 0u + (kk << 2), 1u + (kk << 2), 2u + (kk << 2), 3u + (kk << 2) }));

		for (int64 l1 = n - N * sizeof(float32x4_t) / sizeof(float); i <= l1; i += N * sizeof(float32x4_t) / sizeof(float))
		{
			REP(N) { a[kk] = vld1q_f32(A); A += sizeof(float32x4_t) / sizeof(float); }

			REP(N) f[kk] = vcgtq_f32(min1[kk], a[kk]);

			REP(N) min1[kk] = vbslq_s32(f[kk], a[kk], min1[kk]);

			REP(N) midx[kk] = vbslq_s32(f[kk], nidx[kk], midx[kk]);

			REP(N) nidx[kk] = vaddq_s32(nidx[kk], msep);
		}

		for (int KK = sizeof(min1) / sizeof(min1[0]) / 2; KK >= 1; KK >>= 1)
		{
			REP(KK) f[kk] = vcgtq_f32(min1[kk], min1[kk + KK]);
			REP(KK) min1[kk] = vbslq_s32(f[kk], min1[kk + KK], min1[kk]);
			REP(KK) midx[kk] = vbslq_s32(f[kk], midx[kk + KK], midx[kk]);
		}

		if (vgetq_lane_f32(min1[0], 0) < val)
		{
			val = vgetq_lane_f32(min1[0], 0);
			idx = vgetq_lane_u32(midx[0], 0);
		}

		if (vgetq_lane_f32(min1[0], 1) < val)
		{
			val = vgetq_lane_f32(min1[0], 1);
			idx = vgetq_lane_u32(midx[0], 1);
		}

		if (vgetq_lane_f32(min1[0], 2) < val)
		{
			val = vgetq_lane_f32(min1[0], 2);
			idx = vgetq_lane_u32(midx[0], 2);
		}

		if (vgetq_lane_f32(min1[0], 3) < val)
		{
			val = vgetq_lane_f32(min1[0], 3);
			idx = vgetq_lane_u32(midx[0], 3);
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
	constexpr int N = 4;
	int64 i = 0;
	minv = DBL_MAX;
	maxv = -DBL_MAX;

	if (n >= N * sizeof(float64x2_t) / sizeof(double))
	{
		float64x2_t min1[N], max1[N], a[N];
		REP(N) min1[kk] = vdupq_n_f64(minv);
		REP(N) max1[kk] = vdupq_n_f64(maxv);

		for (int64 l1 = n - N * sizeof(float64x2_t) / sizeof(double); i <= l1; i += N * sizeof(float64x2_t) / sizeof(double))
		{
			REP(N) { a[kk] = vld1q_f64(A); A += sizeof(float64x2_t) / sizeof(double); }

			REP(N)
			{
				min1[kk] = vminq_f64(min1[kk], a[kk]);
				max1[kk] = vmaxq_f64(max1[kk], a[kk]);
			}
		}

		for (int KK = sizeof(max1) / sizeof(max1[0]) / 2; KK >= 1; KK >>= 1)
		{
			REP(KK) min1[kk] = vminq_f64(min1[kk], min1[kk + KK]);
			REP(KK) max1[kk] = vmaxq_f64(max1[kk], max1[kk + KK]);
		}

		minv = Min(vgetq_lane_f64(min1[0], 0), vgetq_lane_f64(min1[0], 1));
		maxv = Max(vgetq_lane_f64(max1[0], 0), vgetq_lane_f64(max1[0], 1));
	}

	for (; i < n; ++i, ++A)
	{
		if (*A < minv) minv = *A;
		if (*A > maxv) maxv = *A;
	}
}

TARGETNEO void GetMinMaxValNEO(float* A, int64 n, float& minv, float& maxv)
{
	constexpr int N = 4;
	int64 i = 0;
	minv = FLT_MAX;
	maxv = -FLT_MAX;

	if (n >= N * sizeof(float32x4_t) / sizeof(float))
	{
		float32x4_t min1[N], max1[N], a[N];
		REP(N) min1[kk] = vdupq_n_f32(minv);
		REP(N) max1[kk] = vdupq_n_f32(maxv);

		for (int64 l1 = n - N * sizeof(float32x4_t) / sizeof(float); i <= l1; i += N * sizeof(float32x4_t) / sizeof(float))
		{
			REP(N) { a[kk] = vld1q_f32(A); A += sizeof(float32x4_t) / sizeof(float); }

			REP(N)
			{
				min1[kk] = vminq_f32(min1[kk], a[kk]);
				max1[kk] = vmaxq_f32(max1[kk], a[kk]);
			}
		}

		for (int KK = sizeof(max1) / sizeof(max1[0]) / 2; KK >= 1; KK >>= 1)
		{
			REP(KK) min1[kk] = vminq_f32(min1[kk], min1[kk + KK]);
			REP(KK) max1[kk] = vmaxq_f32(max1[kk], max1[kk + KK]);
		}

		minv = Min(Min(vgetq_lane_f32(min1[0], 0), vgetq_lane_f32(min1[0], 1)),
				   Min(vgetq_lane_f32(min1[0], 2), vgetq_lane_f32(min1[0], 3)));
		maxv = Max(Max(vgetq_lane_f32(max1[0], 0), vgetq_lane_f32(max1[0], 1)),
				   Max(vgetq_lane_f32(max1[0], 2), vgetq_lane_f32(max1[0], 3)));
	}

	for (; i < n; ++i, ++A)
	{
		if (*A < minv) minv = *A;
		if (*A > maxv) maxv = *A;
	}
}

TARGETNEO double GetMaxValNEO(double* A, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;
	double val = -DBL_MAX;

	if (n >= N * sizeof(float64x2_t) / sizeof(double))
	{
		float64x2_t max1[N];
		REP(N) max1[kk] = vdupq_n_f64(val);

		for (int64 l1 = n - N * sizeof(float64x2_t) / sizeof(double); i <= l1; i += N * sizeof(float64x2_t) / sizeof(double))
		{
			REP(N)
			{
				max1[kk] = vmaxq_f64(max1[kk], vld1q_f64(A));
				A += sizeof(float64x2_t) / sizeof(double);
			}
		}

		for (int KK = sizeof(max1) / sizeof(max1[0]) / 2; KK >= 1; KK >>= 1)
			REP(KK) max1[kk] = vmaxq_f64(max1[kk], max1[kk + KK]);

		val = Max(vgetq_lane_f64(max1[0], 0), vgetq_lane_f64(max1[0], 1));
	}

	for (; i < n; ++i, ++A)
	{
		if (*A < val) continue;
		val = *A;
	}

	return val;
}

TARGETNEO float GetMaxValNEO(float* A, int64 n)
{
	constexpr int N = 8;
	int64 i = 0;
	float val = -FLT_MAX;

	if (n >= N * sizeof(float32x4_t) / sizeof(float))
	{
		float32x4_t max1[N], a[N];
		REP(N) max1[kk] = vdupq_n_f32(val);

		for (int64 l1 = n - N * sizeof(float32x4_t) / sizeof(float); i <= l1; i += N * sizeof(float32x4_t) / sizeof(float))
		{
			REP(N) { a[kk] = vld1q_f32(A); A += sizeof(float32x4_t) / sizeof(float); }

			REP(N) { max1[kk] = vmaxq_f32(max1[kk], a[kk]); }
		}

		for (int KK = sizeof(max1) / sizeof(max1[0]) / 2; KK >= 1; KK >>= 1)
			REP(KK) max1[kk] = vmaxq_f32(max1[kk], max1[kk + KK]);

		val = Max(Max(vgetq_lane_f32(max1[0], 0), vgetq_lane_f32(max1[0], 1)),
				  Max(vgetq_lane_f32(max1[0], 2), vgetq_lane_f32(max1[0], 3)));
	}

	for (; i < n; ++i, ++A)
	{
		if (*A < val) continue;
		val = *A;
	}

	return val;
}

TARGETNEO double GetMaxValNEO(double* A, int64 n, int64 sep)
{
	int64 i = 0;
	double val = -DBL_MAX;

	if (n >= 8)
	{
		float64x2_t max1[2];
		REP(2) max1[kk] = vdupq_n_f64(val);
        REP(4) { __builtin_prefetch(A, 0); A += sep; }

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
            REP(4) { __builtin_prefetch(A, 0); A += sep; }
            max1[0] = vmaxq_f64(max1[0], vld1q_f64(((double[]) { A[-8 * sep], A[-7 * sep] } )));
            max1[1] = vmaxq_f64(max1[1], vld1q_f64(((double[]) { A[-6 * sep], A[-5 * sep] } )));
            
            REP(4) { __builtin_prefetch(A, 0); A += sep; }
            max1[0] = vmaxq_f64(max1[0], vld1q_f64(((double[]) { A[-8 * sep], A[-7 * sep] } )));
            max1[1] = vmaxq_f64(max1[1], vld1q_f64(((double[]) { A[-6 * sep], A[-5 * sep] } )));
		}

		max1[0] = vmaxq_f64(max1[0], max1[1]);

		val = Max(vgetq_lane_f64(max1[0], 0), vgetq_lane_f64(max1[0], 1));
        
        A -= 4 * sep;
	}

	for (; i < n; ++i, A += sep)
	{
		if (*A < val) continue;
		val = *A;
	}

	return val;
}

TARGETNEO float GetMaxValNEO(float* A, int64 n, int64 sep)
{
	constexpr int N = 4;
	int64 i = 0;
	float val = -FLT_MAX;

	if (n >= N * sizeof(float32x4_t) / sizeof(float))
	{
		float32x4_t max1[N];
		REP(N) max1[kk] = vdupq_n_f32(val);

		for (int64 l1 = n - N * sizeof(float32x4_t) / sizeof(float); i <= l1; i += N * sizeof(float32x4_t) / sizeof(float))
		{
			REP(N)
			{
				max1[kk] = vmaxq_f32(max1[kk], vld1q_f32(((float[]) { A[0 * sep], A[1 * sep], A[2 * sep], A[3 * sep] })));
				A += sizeof(float32x4_t) / sizeof(float) * sep;
			}
		}

		for (int KK = sizeof(max1) / sizeof(max1[0]) / 2; KK >= 1; KK >>= 1)
			REP(KK) max1[kk] = vmaxq_f32(max1[kk], max1[kk + KK]);

		val = Max(Max(vgetq_lane_f32(max1[0], 0), vgetq_lane_f32(max1[0], 1)),
				  Max(vgetq_lane_f32(max1[0], 2), vgetq_lane_f32(max1[0], 3)));
	}

	for (; i < n; ++i, A += sep)
	{
		if (*A < val) continue;
		val = *A;
	}

	return val;
}

TARGETNEO double GetMinValNEO(double* A, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;
	double val = DBL_MAX;

	if (n >= N * sizeof(float64x2_t) / sizeof(double))
	{
		float64x2_t min1[N], a[N];
		REP(N) min1[kk] = vdupq_n_f64(val);

		for (int64 l1 = n - N * sizeof(float64x2_t) / sizeof(double); i <= l1; i += N * sizeof(float64x2_t) / sizeof(double))
		{
			REP(N)
			{
				min1[kk] = vminq_f64(min1[kk], vld1q_f64(A));
				A += sizeof(float64x2_t) / sizeof(double);
			}
		}

		for (int KK = sizeof(min1) / sizeof(min1[0]) / 2; KK >= 1; KK >>= 1)
			REP(KK) min1[kk] = vminq_f64(min1[kk], min1[kk + KK]);

		val = Min(vgetq_lane_f64(min1[0], 0), vgetq_lane_f64(min1[0], 1));
	}

	for (; i < n; ++i, ++A)
	{
		if (*A > val) continue;
		val = *A;
	}

	return val;
}

TARGETNEO float GetMinValNEO(float* A, int64 n)
{
	constexpr int N = 8;
	int64 i = 0;
	float val = FLT_MAX;

	if (n >= N * sizeof(float32x4_t) / sizeof(float))
	{
		float32x4_t min1[N], a[N];
		REP(N) min1[kk] = vdupq_n_f32(val);

		for (int64 l1 = n - N * sizeof(float32x4_t) / sizeof(float); i <= l1; i += N * sizeof(float32x4_t) / sizeof(float))
		{
			REP(N) { a[kk] = vld1q_f32(A); A += sizeof(float32x4_t) / sizeof(float); }

			REP(N) { min1[kk] = vminq_f32(min1[kk], a[kk]); }
		}

		for (int KK = sizeof(min1) / sizeof(min1[0]) / 2; KK >= 1; KK >>= 1)
			REP(KK) min1[kk] = vminq_f32(min1[kk], min1[kk + KK]);

		val = Min(Min(vgetq_lane_f32(min1[0], 0), vgetq_lane_f32(min1[0], 1)),
				  Min(vgetq_lane_f32(min1[0], 2), vgetq_lane_f32(min1[0], 3)));
	}

	for (; i < n; ++i, ++A)
	{
		if (*A > val) continue;
		val = *A;
	}

	return val;
}

TARGETNEO int64 GetMinValNEO(int64* A, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;
	int64 val = 0x7FFFFFFFFFFFFFFF;

	if (n >= N * sizeof(uint64x2_t) / sizeof(int64))
	{
		uint64x2_t min1[N], a[N], f[N];
		REP(N) min1[kk] = vdupq_n_u64(0x7FFFFFFFFFFFFFFF);

		for (int64 l1 = n - N * sizeof(uint64x2_t) / sizeof(int64); i <= l1; i += N * sizeof(uint64x2_t) / sizeof(int64))
		{
			REP(N) { a[kk] = vld1q_s64(A); A += sizeof(uint64x2_t) / sizeof(int64); }
			
			REP(N) f[kk] = vcgtq_s64(min1[kk], a[kk]);

			REP(N) min1[kk] = vbslq_s64(f[kk], a[kk], min1[kk]);
		}

		for (int KK = sizeof(min1) / sizeof(min1[0]) / 2; KK >= 1; KK >>= 1)
			REP(KK) min1[kk] = vbslq_s64(vcgtq_s64(min1[kk], min1[kk + KK]), min1[kk + KK], min1[kk]);

		val = Min(vgetq_lane_s64(min1[0], 0), vgetq_lane_s64(min1[0], 1));
	}

	for (; i < n; ++i, ++A)
	{
		if (*A > val) continue;
		val = *A;
	}

	return val;
}

TARGETNEO void SetValNEO(uint* A, ushort* B, int64 n)
{
	int64 i = 0;

	if (n >= 8)
	{
		uint16x8_t b;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			b = vld1q_u16(B);  B += 8;

			vst1q_u32(A, vmovl_u16(vget_low_u16 (b))); A += 4;

			vst1q_u32(A, vmovl_u16(vget_high_u16(b))); A += 4;
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
	int64 i = 0;
	int64 slog = 0; double prod = 1;

	if (n >= 32)
	{
		float64x2_t pd[4], dunder = vdupq_n_f64(DOUBLE_UNDERFLOW), dover = vdupq_n_f64(DOUBLE_OVERFLOW);
		REP(4) pd[kk] = vdupq_n_f64(1.0);

		for (int64 l1 = n - 32; i <= l1; i += 32)
		{
			REP(2)
			{
				pd[0] = vmulq_f64(pd[0], vmulq_f64(vld1q_f64(A + 0), vld1q_f64(A + 8)));
				pd[1] = vmulq_f64(pd[1], vmulq_f64(vld1q_f64(A + 2), vld1q_f64(A + 10)));
				pd[2] = vmulq_f64(pd[2], vmulq_f64(vld1q_f64(A + 4), vld1q_f64(A + 12)));
				pd[3] = vmulq_f64(pd[3], vmulq_f64(vld1q_f64(A + 6), vld1q_f64(A + 14)));

				A += 16;
			}

			REP(4) AddExponentNEO(slog, pd[kk]);
		}

		REP(4) ChargeLogNEO(slog, prod, pd[kk]);
	}

	for (; i < n; ++i, ++A)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);
	return prod;
}

TARGETNEO double LogProdNEO(float* A, int64 n)
{
	int64 i = 0;
	int64 slog = 0; double prod = 1;

	if (n >= 64)
	{
		float32x4_t a[4];
        float64x2_t pd[4];
		REP(4) pd[kk] = vdupq_n_f64(1.0);

		for (int64 l1 = n - 64; i <= l1; i += 64)
		{
			REP(4)
			{
				a[0] = vld1q_f32(A); A += 4;
				a[1] = vld1q_f32(A); A += 4;
				a[2] = vld1q_f32(A); A += 4;
				a[3] = vld1q_f32(A); A += 4;

				pd[0] = vmulq_f64(pd[0], vmulq_f64(vcvt_f64_f32(vget_low_f32 (a[0])), vcvt_f64_f32(vget_low_f32 (a[2]))));
				pd[1] = vmulq_f64(pd[1], vmulq_f64(vcvt_f64_f32(vget_high_f32(a[0])), vcvt_f64_f32(vget_high_f32(a[2]))));
				pd[2] = vmulq_f64(pd[2], vmulq_f64(vcvt_f64_f32(vget_low_f32 (a[1])), vcvt_f64_f32(vget_low_f32 (a[3]))));
				pd[3] = vmulq_f64(pd[3], vmulq_f64(vcvt_f64_f32(vget_high_f32(a[1])), vcvt_f64_f32(vget_high_f32(a[3]))));
			}

			REP(4) AddExponentNEO(slog, pd[kk]);
		}

		REP(4) ChargeLogNEO(slog, prod, pd[kk]);
	}

	for (; i < n; ++i, ++A)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);
	return prod;
}

TARGETNEO float LogProdNEOx(float* A, int64 n)
{
	int64 i = 0;
	int64 slog = 0; double prod = 1;

	if (n >= 32)
	{
		float32x4_t pd[4];
		REP(4) pd[kk] = vdupq_n_f32(1.0f);

		for (int64 l1 = n - 32; i <= l1; i += 32)
		{
			pd[0] = vmulq_f32(pd[0], vld1q_f32(A +  0));
			pd[1] = vmulq_f32(pd[1], vld1q_f32(A +  4));
			pd[2] = vmulq_f32(pd[2], vld1q_f32(A +  8));
			pd[3] = vmulq_f32(pd[3], vld1q_f32(A + 12));
			pd[0] = vmulq_f32(pd[0], vld1q_f32(A + 16));
			pd[1] = vmulq_f32(pd[1], vld1q_f32(A + 20));
			pd[2] = vmulq_f32(pd[2], vld1q_f32(A + 24));
			pd[3] = vmulq_f32(pd[3], vld1q_f32(A + 28));

			A += 32;
			REP(4) AddExponentNEO(slog, pd[kk]);
		}

		REP(4) ChargeLogNEO(slog, prod, pd[kk]);
	}

	for (; i < n; ++i, ++A)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);
	return prod;
}

TARGETNEO double LogProdNEO(double* A, int64 n, int64 sep)
{
	int64 i = 0;
	int64 slog = 0; double prod = 1;

	if (n >= 8)
	{
		REP(4) { __builtin_prefetch(A, 0); A += sep; }
		float64x2_t pd[4], v1[4], dunder = vdupq_n_f64(DOUBLE_UNDERFLOW), dover = vdupq_n_f64(DOUBLE_OVERFLOW);
		uint64x2_t f1, f2;
		REP(4) pd[kk] = vdupq_n_f64(1.0);

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			REP(4)
			{
				__builtin_prefetch(A, 0); A += sep;
				__builtin_prefetch(A, 0); A += sep;
				v1[kk] = vld1q_f64(((double[]) { A[-6 * sep], A[-5 * sep] }));
			}

			f1 = vdupq_n_f64(0);
			f2 = vdupq_n_f64(0);

			REP(4)
			{
				pd[kk] = vmulq_f64(pd[kk], v1[kk]);

				f1 = vorrq_u64(f1, vcltq_f64(pd[kk], dunder));
				f2 = vorrq_u64(f2, vcltq_f64(dover, pd[kk]));
			}

			f1 = vorrq_u64(f1, f2);

			if (vgetq_lane_u64(f1, 0) | vgetq_lane_u64(f1, 1)) [[unlikely]]
				REP(4) AddExponentNEO(slog, pd[kk]);
		}

		REP(4) ChargeLogNEO(slog, prod, pd[kk]);
		A -= 4 * sep;
	}

	for (; i < n; ++i, A += sep)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);
	return prod;
}

TARGETNEO double LogProdNEO(float* A, int64 n, int64 sep)
{
	int64 i = 0;
	int64 slog = 0; double prod = 1;

	if (n >= 32)
	{
		REP(4) { __builtin_prefetch(A, 0); A += sep; }

		float64x2_t pd[2];
		float32x4_t a1;
		REP(2) pd[kk] = vdupq_n_f64(1.0f);

		for (int64 l1 = n - 32; i <= l1; i += 32)
		{
			REP(8)
			{
				__builtin_prefetch(A, 0); A += sep;
				__builtin_prefetch(A, 0); A += sep;
				__builtin_prefetch(A, 0); A += sep;
				__builtin_prefetch(A, 0); A += sep;

				a1 = vld1q_f32(((float[]) { A[-8 * sep], A[-7 * sep], A[-6 * sep], A[-5 * sep] }));

				pd[0] = vmulq_f64(pd[0], vcvt_f64_f32(vget_low_f32 (a1)));
				pd[1] = vmulq_f64(pd[1], vcvt_f64_f32(vget_high_f32(a1)));
			}

			REP(2) AddExponentNEO(slog, pd[kk]);
		}

		REP(2) ChargeLogNEO(slog, prod, pd[kk]);
		A -= 4 * sep;
	}

	for (; i < n; ++i, A += sep)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);
	return prod;
}

TARGETNEO float LogProdNEOx(float* A, int64 n, int64 sep)
{
	int64 i = 0;
	int64 slog = 0; double prod = 1;

	if (n >= 32)
	{
		float32x4_t pd[4];
		REP(4) pd[kk] = vdupq_n_f32(1.0f);

		for (int64 l1 = n - 32; i <= l1; i += 32)
		{
			pd[0] = vmulq_f32(pd[0], vld1q_f32(((float[]) { A[0 * sep], A[1 * sep], A[2 * sep], A[3 * sep] }))); A += 4 * sep;
			pd[1] = vmulq_f32(pd[1], vld1q_f32(((float[]) { A[0 * sep], A[1 * sep], A[2 * sep], A[3 * sep] }))); A += 4 * sep;
			pd[2] = vmulq_f32(pd[2], vld1q_f32(((float[]) { A[0 * sep], A[1 * sep], A[2 * sep], A[3 * sep] }))); A += 4 * sep;
			pd[3] = vmulq_f32(pd[3], vld1q_f32(((float[]) { A[0 * sep], A[1 * sep], A[2 * sep], A[3 * sep] }))); A += 4 * sep;

			pd[0] = vmulq_f32(pd[0], vld1q_f32(((float[]) { A[0 * sep], A[1 * sep], A[2 * sep], A[3 * sep] }))); A += 4 * sep;
			pd[1] = vmulq_f32(pd[1], vld1q_f32(((float[]) { A[0 * sep], A[1 * sep], A[2 * sep], A[3 * sep] }))); A += 4 * sep;
			pd[2] = vmulq_f32(pd[2], vld1q_f32(((float[]) { A[0 * sep], A[1 * sep], A[2 * sep], A[3 * sep] }))); A += 4 * sep;
			pd[3] = vmulq_f32(pd[3], vld1q_f32(((float[]) { A[0 * sep], A[1 * sep], A[2 * sep], A[3 * sep] }))); A += 4 * sep;

			REP(4) AddExponentNEO(slog, pd[kk]);
		}

		REP(4) ChargeLogNEO(slog, prod, pd[kk]);
	}

	for (; i < n; ++i, A += sep)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);
	return prod;
}

TARGETNEO double LogProdDivNEO(double* A, double* B, int64 n, int64 sep)
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
	int64 slog = 0; double prod = 1;

	if (n >= 4)
	{
		__builtin_prefetch(A, 0); A += sep;
		__builtin_prefetch(A, 0); A += sep;
		__builtin_prefetch(A, 0); A += sep;
		__builtin_prefetch(A, 0); A += sep;
		__builtin_prefetch(B, 0); B += sep;
		__builtin_prefetch(B, 0); B += sep;
		__builtin_prefetch(B, 0); B += sep;
		__builtin_prefetch(B, 0); B += sep;

		float64x2_t pd = vdupq_n_f64(1.0), dunder = vdupq_n_f64(DOUBLE_UNDERFLOW), dover = vdupq_n_f64(DOUBLE_OVERFLOW);
		uint64x2_t flag;

		for (int64 l1 = n - 4; i <= l1; i += 4)
		{
			__builtin_prefetch(A, 0); A += sep;
			__builtin_prefetch(A, 0); A += sep;
			__builtin_prefetch(A, 0); A += sep;
			__builtin_prefetch(A, 0); A += sep;
			__builtin_prefetch(B, 0); B += sep;
			__builtin_prefetch(B, 0); B += sep;
			__builtin_prefetch(B, 0); B += sep;
			__builtin_prefetch(B, 0); B += sep;
			
			double a[4] = { A[-8 * sep], A[-7 * sep], A[-6 * sep], A[-5 * sep] };
			double b[4] = { B[-4 * sep], B[-3 * sep], B[-6 * sep], B[-5 * sep] };

			pd = vmulq_f64(pd, vdivq_f64(
				vmulq_f64(vld1q_f64(a), vld1q_f64(a + 2)), 
				vmulq_f64(vld1q_f64(b), vld1q_f64(b + 2))));

			flag = vorrq_u64(vcltq_f64(pd, dunder), vcltq_f64(dover, pd));

			if (vgetq_lane_u64(flag, 0) | vgetq_lane_u64(flag, 1)) [[unlikely]]
				AddExponentNEO(slog, pd);
		}

		ChargeLogNEO(slog, prod, pd);
		A -= 4 * sep; 
		B -= 4 * sep;
	}

	for (; i < n; ++i, A += sep, B += sep)
		ChargeLog(slog, prod, *A / *B);

	CloseLog(slog, prod);
	return prod;
	*/
}

TARGETNEO double LogProdDivNEO(float* A, float* B, int64 n, int64 sep)
{
	int64 i = 0;
	int64 slog1 = 0; double prod1 = 1;
	int64 slog2 = 0; double prod2 = 1;

	if (n >= 64)
	{
		float64x2_t pd1[4], pd2[4];
		REP(4) pd1[kk] = pd2[kk] = vdupq_n_f64(1.0);

		for (int64 l1 = n - 64; i <= l1; i += 64)
		{
			REP(8)
			{
				pd1[0] = vmulq_f64(pd1[0], vld1q_f64(((double[]) { A[0 * sep], A[1 * sep] }))); A += sep * 2;
				pd2[0] = vmulq_f64(pd2[0], vld1q_f64(((double[]) { B[0 * sep], B[1 * sep] }))); B += sep * 2;

				pd1[1] = vmulq_f64(pd1[1], vld1q_f64(((double[]) { A[0 * sep], A[1 * sep] }))); A += sep * 2;
				pd2[1] = vmulq_f64(pd2[1], vld1q_f64(((double[]) { B[0 * sep], B[1 * sep] }))); B += sep * 2;

				pd1[2] = vmulq_f64(pd1[2], vld1q_f64(((double[]) { A[0 * sep], A[1 * sep] }))); A += sep * 2;
				pd2[2] = vmulq_f64(pd2[2], vld1q_f64(((double[]) { B[0 * sep], B[1 * sep] }))); B += sep * 2;

				pd1[3] = vmulq_f64(pd1[3], vld1q_f64(((double[]) { A[0 * sep], A[1 * sep] }))); A += sep * 2;
				pd2[3] = vmulq_f64(pd2[3], vld1q_f64(((double[]) { B[0 * sep], B[1 * sep] }))); B += sep * 2;
			}

			REP(4) AddExponentNEO(slog1, pd1[kk]);
			REP(4) AddExponentNEO(slog2, pd2[kk]);
		}

		REP(4) ChargeLogNEO(slog1, prod1, pd1[kk]);
		REP(4) ChargeLogNEO(slog2, prod2, pd2[kk]);
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

TARGETNEO float LogProdDivNEOx(float* A, float* B, int64 n, int64 sep)
{
	int64 i = 0;
	int64 slog1 = 0; double prod1 = 1;
	int64 slog2 = 0; double prod2 = 1;

	if (n >= 64)
	{
		float32x4_t pd1[4], pd2[4];
		REP(4) pd1[kk] = pd2[kk] = vdupq_n_f32(1.0f);

		for (int64 l1 = n - 64; i <= l1; i += 64)
		{
			REP(4)
			{
				pd1[0] = vmulq_f32(pd1[0], vld1q_f32(((float[]) { A[0 * sep], A[1 * sep], A[2 * sep], A[3 * sep] }))); A += 4 * sep;
				pd2[0] = vmulq_f32(pd2[0], vld1q_f32(((float[]) { B[0 * sep], B[1 * sep], B[2 * sep], B[3 * sep] }))); B += 4 * sep;

				pd1[1] = vmulq_f32(pd1[1], vld1q_f32(((float[]) { A[0 * sep], A[1 * sep], A[2 * sep], A[3 * sep] }))); A += 4 * sep;
				pd2[1] = vmulq_f32(pd2[1], vld1q_f32(((float[]) { B[0 * sep], B[1 * sep], B[2 * sep], B[3 * sep] }))); B += 4 * sep;

				pd1[2] = vmulq_f32(pd1[2], vld1q_f32(((float[]) { A[0 * sep], A[1 * sep], A[2 * sep], A[3 * sep] }))); A += 4 * sep;
				pd2[2] = vmulq_f32(pd2[2], vld1q_f32(((float[]) { B[0 * sep], B[1 * sep], B[2 * sep], B[3 * sep] }))); B += 4 * sep;

				pd1[3] = vmulq_f32(pd1[3], vld1q_f32(((float[]) { A[0 * sep], A[1 * sep], A[2 * sep], A[3 * sep] }))); A += 4 * sep;
				pd2[3] = vmulq_f32(pd2[3], vld1q_f32(((float[]) { B[0 * sep], B[1 * sep], B[2 * sep], B[3 * sep] }))); B += 4 * sep;
			}

			REP(4) AddExponentNEO(slog1, pd1[kk]);
			REP(4) AddExponentNEO(slog2, pd2[kk]);
		}

		REP(4) ChargeLogNEO(slog1, prod1, pd1[kk]);
		REP(4) ChargeLogNEO(slog2, prod2, pd2[kk]);
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

	if (n >= 16)
	{
		uint8x16_t z = vdupq_n_u8(0), v;
		uint64 s1 = 0, s2 = 0;
		byte* Ab = A;

		for (int64 l1 = n - 16; i <= l1; i += 16)
		{
			v = vld1q_u8(A);
			v = vceqq_u8(v, z);
			s1 += __builtin_popcountll(vgetq_lane_u64((uint64x2_t)v, 0));
			s2 += __builtin_popcountll(vgetq_lane_u64((uint64x2_t)v, 1));
			A += 16;
		}

		re = (A - Ab) - ((s1 + s2) >> 3);
	}

	for (; i < n; ++i, ++A)
		if (*A) re++;

	return (int64)re;
}

TARGETNEO double SumNEO(double* A, int64 n)
{
	constexpr int N = 8;
	int64 i = 0;
	double re = 0;

	if (n >= N * sizeof(float64x2_t) / sizeof(double))
	{
		float64x2_t s[N], a[N];
		REP(N) s[kk] = vdupq_n_f64(0);

		for (int64 l1 = n - N * sizeof(float64x2_t) / sizeof(double); i <= l1; i += N * sizeof(float64x2_t) / sizeof(double))
		{
			REP(N) { a[kk] = vld1q_f64(A); A += sizeof(float64x2_t) / sizeof(double); }

			REP(N) s[kk] = vaddq_f64(s[kk], a[kk]);
		}

		for (int KK = sizeof(s) / sizeof(s[0]) / 2; KK >= 1; KK >>= 1)
			REP(KK) s[kk] = vaddq_f64(s[kk], s[kk + KK]);

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
	constexpr int N = 16;
	int64 i = 0;
	double re = 0;

	if (n >= N / 2 * sizeof(float32x4_t) / sizeof(float))
	{
		float64x2_t s[N];
		REP(N) s[kk] = vdupq_n_f64(0);

		for (int64 l1 = n - N / 2 * sizeof(float32x4_t) / sizeof(float); i <= l1; i += N / 2 * sizeof(float32x4_t) / sizeof(float))
		{
			REP(N / 2)
			{
				float32x4_t a = vld1q_f32(A); A += sizeof(float32x4_t) / sizeof(float);
				s[0 + (kk << 1)] = vaddq_f64(s[0 + (kk << 1)], vcvt_f64_f32(vget_low_f32 (a)));
				s[1 + (kk << 1)] = vaddq_f64(s[1 + (kk << 1)], vcvt_f64_f32(vget_high_f32(a)));
			}
		}

		for (int KK = sizeof(s) / sizeof(s[0]) / 2; KK >= 1; KK >>= 1)
			REP(KK) s[kk] = vaddq_f64(s[kk], s[kk + KK]);

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
	constexpr int N = 8;
	int64 i = 0;
	float re = 0;

	if (n >= N * sizeof(float32x4_t) / sizeof(float))
	{
		float32x4_t s[N], a[N];
		REP(N) s[kk] = vdupq_n_f32(0);

		for (int64 l1 = n - N * sizeof(float32x4_t) / sizeof(float); i <= l1; i += N * sizeof(float32x4_t) / sizeof(float))
		{
			REP(N) { a[kk] = vld1q_f32(A); A += sizeof(float32x4_t) / sizeof(float); }

			REP(N) s[kk] = vaddq_f32(s[kk], a[kk]);
		}

		for (int KK = sizeof(s) / sizeof(s[0]) / 2; KK >= 1; KK >>= 1)
			REP(KK) s[kk] = vaddq_f32(s[kk], s[kk + KK]);

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
	constexpr int N = 4;
	uint64 re = 0;
	int64 i = 0;

	if (n >= N * sizeof(uint64x2_t) / sizeof(byte))
	{
		uint64x2_t s[N], a[N], z = vdupq_n_u64(0);
		REP(N) s[kk] = vdupq_n_u64(0);

		for (int64 l1 = n - N * sizeof(uint64x2_t) / sizeof(byte); i <= l1; i += N * sizeof(uint64x2_t) / sizeof(byte))
		{
			REP(N) { a[kk] = vld1q_u8(A); A += sizeof(uint64x2_t) / sizeof(byte); }

			REP(N) s[kk] = vaddq_u64(s[kk], vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(a[kk]))));
		}

		for (int KK = sizeof(s) / sizeof(s[0]) / 2; KK >= 1; KK >>= 1)
			REP(KK) s[kk] = vaddq_u64(s[kk], s[kk + KK]);

		re += vgetq_lane_u64(s[0], 0) + vgetq_lane_u64(s[0], 1);
	}

	for (; i < n; ++i)
		re += *A++;

	return re;
}

TARGETNEO double SumNEO(double* A, int64 n, int64 sep)
{
	constexpr int N = 8;
	int64 i = 0;
	double re = 0;

	if (n >= N * sizeof(float64x2_t) / sizeof(double))
	{
		REP(4) __builtin_prefetch(A + sep * kk, 0);

		float64x2_t s[N], a[N];
		REP(N) s[kk] = vdupq_n_f64(0);

		for (int64 l1 = n - N * sizeof(float64x2_t) / sizeof(double); i <= l1; i += N * sizeof(float64x2_t) / sizeof(double))
		{
			REP(N * 2) __builtin_prefetch(A + sep * kk, 0);

			REP(N) { a[kk] = vld1q_f64(((double[]) { A[0], A[sep] })); A += sep * sizeof(float64x2_t) / sizeof(double); }

			REP(N) s[kk] = vaddq_f64(s[kk], a[kk]);
		}

		for (int KK = sizeof(s) / sizeof(s[0]) / 2; KK >= 1; KK >>= 1)
			REP(KK) s[kk] = vaddq_f64(s[kk], s[kk + KK]);

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
	int64 i = 0;
	double re = 0;

	if (n >= 4)
	{
		REP(4) { __builtin_prefetch(A, 0); A += sep; }

		float64x2_t s1 = vdupq_n_f64(0), s2 = vdupq_n_f64(0);

		for (int64 l1 = n - 4; i <= l1; i += 4)
		{
			REP(4) { __builtin_prefetch(A, 0); A += sep; }

			s1 = vaddq_f64(s1, vld1q_f64(((double[]) { A[-8 * sep], A[-7 * sep] })));
			s2 = vaddq_f64(s2, vld1q_f64(((double[]) { A[-6 * sep], A[-5 * sep] })));
		}

		s1 = vaddq_f64(s1, s2);
		re = _neo_reduce_add_pd(s1);
		A -= 4 * sep;
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
	int64 i = 0;
	float re = 0;

	if (n >= 8)
	{
		float32x4_t s1 = vdupq_n_f64(0), s2 = vdupq_n_f64(0);

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			s1 = vaddq_f32(s1, vld1q_f32(((float[]) { A[0 * sep], A[1 * sep], A[2 * sep], A[3 * sep] })));
			s2 = vaddq_f32(s2, vld1q_f32(((float[]) { A[4 * sep], A[5 * sep], A[6 * sep], A[7 * sep] })));
			A += 8 * sep;
		}

		s1 = vaddq_f32(s1, s2);
		re = _neo_reduce_add_ps(s1);
	}

	for (; i < n; ++i, A += sep)
	{
		volatile float v1 = *A;
		re += v1;
	}

	return re;
}

TARGETNEO void SumNEO(double* A, double** B, int64 k, int64 n)
{
	constexpr int N = 16;
	int64 i = 0;

	if (n >= N * sizeof(float64x2_t) / sizeof(double))
	{
		float64x2_t a[N];

		for (int64 l1 = n - N * sizeof(float64x2_t) / sizeof(double); i <= l1; i += N * sizeof(float64x2_t) / sizeof(double))
		{
			REP(N) a[kk] = vdupq_n_f64(0);

			for (int64 j = 0; j < k; ++j)
				REP(N) a[kk] = vaddq_f64(a[kk], vld1q_f64(&B[j][i + kk * sizeof(float64x2_t) / sizeof(double)]));

			REP(N) { vst1q_f64(A, a[kk]); A += sizeof(float64x2_t) / sizeof(double); }
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

TARGETNEO void SumNEO(float* A, float** B, int64 k, int64 n)
{
	constexpr int N = 16;
	int64 i = 0;

	if (n >= N * sizeof(float32x4_t) / sizeof(float))
	{
		float64x2_t a[N * 2];

		for (int64 l1 = n - N * sizeof(float32x4_t) / sizeof(float); i <= l1; i += N * sizeof(float32x4_t) / sizeof(float))
		{
			REP(N * 2) a[kk] = vdupq_n_f64(0);

			for (int64 j = 0; j < k; ++j)
				REP(N)
				{
					float64x2_t b = vld1q_f32(&B[j][i + kk * sizeof(float32x4_t) / sizeof(float)]);
					a[0 + (kk << 1)] = vaddq_f64(a[0 + (kk << 1)], vcvt_f64_f32(vget_low_f32(b)));
					a[1 + (kk << 1)] = vaddq_f64(a[1 + (kk << 1)], vcvt_f64_f32(vget_high_f32(b)));
				}

			REP(N)
			{
				vst1q_f32(A, vcombine_f32(
					vcvt_f32_f64(a[0 + (kk << 1)]),
					vcvt_f32_f64(a[1 + (kk << 1)])));

				A += sizeof(float32x4_t) / sizeof(float);
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
}

TARGETNEO double ProdNEO(double* A, int64 n)
{
	constexpr int N = 8;
	int64 i = 0;
	double re = 1;

	if (n >= N * sizeof(float64x2_t) / sizeof(double))
	{
		float64x2_t s[N], a[N];
		REP(N) s[kk] = vdupq_n_f64(1);

		for (int64 l1 = n - N * sizeof(float64x2_t) / sizeof(double); i <= l1; i += N * sizeof(float64x2_t) / sizeof(double))
		{
			REP(N) { a[kk] = vld1q_f64(A); A += sizeof(float64x2_t) / sizeof(double); }

			REP(N) s[kk] = vmulq_f64(s[kk], a[kk]);
		}

		for (int KK = sizeof(s) / sizeof(s[0]) / 2; KK >= 1; KK >>= 1)
			REP(KK) s[kk] = vmulq_f64(s[kk], s[kk + KK]);

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
	constexpr int N = 16;
	int64 i = 0;
	double re = 1;

	if (n >= N / 2 * sizeof(float32x4_t) / sizeof(float))
	{
		float64x2_t s[N];
		REP(N) s[kk] = vdupq_n_f64(1);

		for (int64 l1 = n - N / 2 * sizeof(float32x4_t) / sizeof(float); i <= l1; i += N / 2 * sizeof(float32x4_t) / sizeof(float))
		{
			REP(N / 2)
			{
				float32x4_t v1 = vld1q_f32(A); A += sizeof(float32x4_t) / sizeof(float);
				float64x2_t a1 = vcvt_f64_f32(vget_low_f32 (v1));
				float64x2_t a2 = vcvt_f64_f32(vget_high_f32(v1));

				s[0 + (kk << 1)] = vmulq_f64(s[0 + (kk << 1)], a1);
				s[1 + (kk << 1)] = vmulq_f64(s[1 + (kk << 1)], a2);
			}
		}

		for (int KK = sizeof(s) / sizeof(s[0]) / 2; KK >= 1; KK >>= 1)
			REP(KK) s[kk] = vmulq_f64(s[kk], s[kk + KK]);

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
	constexpr int N = 8;
	int64 i = 0;
	float re = 1;

	if (n >= N * sizeof(float32x4_t) / sizeof(float))
	{
		float32x4_t s[N], a[N];
		REP(N) s[kk] = vdupq_n_f32(1);

		for (int64 l1 = n - N * sizeof(float32x4_t) / sizeof(float); i <= l1; i += N * sizeof(float32x4_t) / sizeof(float))
		{
			REP(N) { a[kk] = vld1q_f32(A); A += sizeof(float32x4_t) / sizeof(float); }

			REP(N) s[kk] = vmulq_f32(s[kk], a[kk]);
		}

		for (int KK = sizeof(s) / sizeof(s[0]) / 2; KK >= 1; KK >>= 1)
			REP(KK) s[kk] = vmulq_f32(s[kk], s[kk + KK]);

		re = _neo_reduce_mul_ps(s[0]);
	}

	for (; i < n; ++i)
		re *= *A++;

	return re;
}

TARGETNEO double ProdNEO(double* A, int64 n, int64 sep)
{
	int64 i = 0;
	double re = 1;

	if (n >= 8)
	{
		REP(8) { __builtin_prefetch(A, 0); A += sep; }

		float64x2_t pd1 = vdupq_n_f64(1.0), pd2 = vdupq_n_f64(1.0);
		float64x2_t pd3 = vdupq_n_f64(1.0), pd4 = vdupq_n_f64(1.0);

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			REP(2) { __builtin_prefetch(A, 0); A += sep; }
			pd1 = vmulq_f64(pd1, vld1q_f64(((double[]) { A[-10 * sep], A[-9 * sep] })));

			REP(2) { __builtin_prefetch(A, 0); A += sep; }
			pd2 = vmulq_f64(pd2, vld1q_f64(((double[]) { A[-10 * sep], A[-9 * sep] })));

			REP(2) { __builtin_prefetch(A, 0); A += sep; }
			pd3 = vmulq_f64(pd3, vld1q_f64(((double[]) { A[-10 * sep], A[-9 * sep] })));

			REP(2) { __builtin_prefetch(A, 0); A += sep; }
			pd4 = vmulq_f64(pd4, vld1q_f64(((double[]) { A[-10 * sep], A[-9 * sep] })));
		}

		pd1 = vmulq_f64(vmulq_f64(pd1, pd2), vmulq_f64(pd3, pd4));
		re = _neo_reduce_mul_pd(pd1);
		A -= 8 * sep;
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
	int64 i = 0;
	double re = 1;

	if (n >= 8)
	{
		REP(8) { __builtin_prefetch(A, 0); A += sep; }

		float64x2_t pd1 = vdupq_n_f64(1.0), pd2 = vdupq_n_f64(1.0);
		float64x2_t pd3 = vdupq_n_f64(1.0), pd4 = vdupq_n_f64(1.0);

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			REP(2) { __builtin_prefetch(A, 0); A += sep; }
			pd1 = vmulq_f64(pd1, vld1q_f64(((double[]) { A[-10 * sep], A[-9 * sep] })));

			REP(2) { __builtin_prefetch(A, 0); A += sep; }
			pd2 = vmulq_f64(pd2, vld1q_f64(((double[]) { A[-10 * sep], A[-9 * sep] })));

			REP(2) { __builtin_prefetch(A, 0); A += sep; }
			pd3 = vmulq_f64(pd3, vld1q_f64(((double[]) { A[-10 * sep], A[-9 * sep] })));

			REP(2) { __builtin_prefetch(A, 0); A += sep; }
			pd4 = vmulq_f64(pd4, vld1q_f64(((double[]) { A[-10 * sep], A[-9 * sep] })));
		}

		pd1 = vmulq_f64(vmulq_f64(pd1, pd2), vmulq_f64(pd3, pd4));
		re = _neo_reduce_mul_pd(pd1);
		A -= 8 * sep;
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
	constexpr int N = 4;
	int64 i = 0;
	volatile float re = 1;

	if (n >= N * sizeof(float32x4_t) / sizeof(float))
	{
		float32x4_t pd[4];
		REP(N) pd[kk] = vdupq_n_f32(1.0f);

		for (int64 l1 = n - N * sizeof(float32x4_t) / sizeof(float); i <= l1; i += N * sizeof(float32x4_t) / sizeof(float))
		{
			REP(N)
			{
				pd[kk] = vmulq_f32(pd[kk], vld1q_f32(((float[]) { A[0 * sep], A[1 * sep], A[2 * sep], A[3 * sep] })));
				A += 4 * sep;
			}
		}

		for (int KK = sizeof(pd) / sizeof(pd[0]) / 2; KK >= 1; KK >>= 1)
			REP(KK) pd[kk] = vmulq_f32(pd[kk], pd[kk + KK]);

		re = _neo_reduce_mul_ps(pd[0]);
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
	constexpr int N = 8;
	int64 i = 0;
	double re = 0;

	if (n >= N * sizeof(float64x2_t) / sizeof(double))
	{
		float64x2_t s[N], a[N];
		REP(N) s[kk] = vdupq_n_f64(0);

		for (int64 l1 = n - N * sizeof(float64x2_t) / sizeof(double); i <= l1; i += N * sizeof(float64x2_t) / sizeof(double))
		{
			REP(N) { a[kk] = vld1q_f64(A); A += sizeof(float64x2_t) / sizeof(double); }

			REP(N) a[kk] = vmulq_f64(a[kk], a[kk]);

			REP(N) s[kk] = vaddq_f64(s[kk], a[kk]);
		}

		for (int KK = sizeof(s) / sizeof(s[0]) / 2; KK >= 1; KK >>= 1)
			REP(KK) s[kk] = vaddq_f64(s[kk], s[kk + KK]);

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
	constexpr int N = 8;
	int64 i = 0;
	double re = 0;

	if (n >= N / 2 * sizeof(float32x4_t) / sizeof(float))
	{
		float64x2_t s[N];
		REP(N) s[kk] = vdupq_n_f64(0);

		for (int64 l1 = n - N / 2 * sizeof(float32x4_t) / sizeof(float); i <= l1; i += N / 2 * sizeof(float32x4_t) / sizeof(float))
		{
			REP(N / 2)
			{
				float32x4_t v1 = vld1q_f32(A); A += sizeof(float32x4_t) / sizeof(float);
				float64x2_t a1 = vcvt_f64_f32(vget_low_f32 (v1));
				float64x2_t a2 = vcvt_f64_f32(vget_high_f32(v1));

				a1 = vmulq_f64(a1, a1);
				a2 = vmulq_f64(a2, a2);

				s[0 + (kk << 1)] = vaddq_f64(s[0 + (kk << 1)], a1);
				s[1 + (kk << 1)] = vaddq_f64(s[1 + (kk << 1)], a2);
			}
		}

		for (int KK = sizeof(s) / sizeof(s[0]) / 2; KK >= 1; KK >>= 1)
			REP(KK) s[kk] = vaddq_f64(s[kk], s[kk + KK]);

		re = _neo_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, ++A)
	{
		volatile double v1 = (double)*A * (double)*A;
		re += v1;
	}

	return re;
}

TARGETNEO float SumSquareNEOx(float* A, int64 n)
{
	constexpr int N = 8;
	int64 i = 0;
	float re = 0;

	if (n >= N * sizeof(float32x4_t) / sizeof(float))
	{
		float32x4_t s[N], a[N];
		REP(N) s[kk] = vdupq_n_f32(0);

		for (int64 l1 = n - N * sizeof(float32x4_t) / sizeof(float); i <= l1; i += N * sizeof(float32x4_t) / sizeof(float))
		{
			REP(N) { a[kk] = vld1q_f32(A); A += sizeof(float32x4_t) / sizeof(float); }

			REP(N) a[kk] = vmulq_f32(a[kk], a[kk]);

			REP(N) s[kk] = vaddq_f32(s[kk], a[kk]);
		}

		for (int KK = sizeof(s) / sizeof(s[0]) / 2; KK >= 1; KK >>= 1)
			REP(KK) s[kk] = vaddq_f32(s[kk], s[kk + KK]);

		re = _neo_reduce_add_ps(s[0]);
	}

	for (; i < n; ++i, ++A)
	{
		volatile float v1 = *A * *A;
		re += v1;
	}

	return re;
}

TARGETNEO int64 SumSquareNEO(byte* A, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;
	uint64 re = 0;

	if (n >= 64)
	{
		uint64x2_t s1 = vdupq_n_u64(0), s2 = vdupq_n_u64(0);
		uint8x16_t a[4];

		for (int64 l1 = n - 64; i <= l1; i += 64)
		{
			REP(4) { a[kk] = vld1q_u8(A); A += 16; }

			REP(4) a[kk] = vmulq_u8(a[kk], a[kk]);

			a[0] = vaddq_u8(a[0], a[1]);
			a[2] = vaddq_u8(a[2], a[3]);

			s1 = vaddq_u64(s1, vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(a[0]))));
			s2 = vaddq_u64(s2, vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(a[2]))));
		}
		
		s1 = vaddq_u64(s1, s2);
		re += vgetq_lane_u64(s1, 0) + vgetq_lane_u64(s1, 1);
	}

	for (; i < n; ++i, ++A)
		re += *A * *A;

	return re;
}

TARGETNEO void SumSumSquareNEO(double* A, int64 n, double& sum, double& sumsq)
{
	constexpr int N = 8;
	int64 i = 0;
	double re1 = 0, re2 = 0;


	if (n >= N * sizeof(float64x2_t) / sizeof(double))
	{
		float64x2_t s1[N], s2[N], a[N];
		REP(N) s1[kk] = s2[kk] = vdupq_n_f64(0);

		for (int64 l1 = n - N * sizeof(float64x2_t) / sizeof(double); i <= l1; i += N * sizeof(float64x2_t) / sizeof(double))
		{
			REP(N) { a[kk] = vld1q_f64(A); A += sizeof(float64x2_t) / sizeof(double); }

			REP(N) s1[kk] = vaddq_f64(s1[kk], a[kk]);

			REP(N) a[kk] = vmulq_f64(a[kk], a[kk]);

			REP(N) s2[kk] = vaddq_f64(s2[kk], a[kk]);
		}

		for (int KK = sizeof(s1) / sizeof(s1[0]) / 2; KK >= 1; KK >>= 1)
		{
			REP(KK) s1[kk] = vaddq_f64(s1[kk], s1[kk + KK]);
			REP(KK) s2[kk] = vaddq_f64(s2[kk], s2[kk + KK]);
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
	constexpr int N = 16;
	int64 i = 0;
	double re1 = 0, re2 = 0;

	if (n >= N / 2 * sizeof(float32x4_t) / sizeof(float))
	{
		float64x2_t s1[N], s2[N];
		REP(N) s1[kk] = s2[kk] = vdupq_n_f64(0);

		for (int64 l1 = n - N / 2 * sizeof(float32x4_t) / sizeof(float); i <= l1; i += N / 2 * sizeof(float32x4_t) / sizeof(float))
		{
			REP(N / 2)
			{
				float32x4_t v1 = vld1q_f32(A); A += sizeof(float32x4_t) / sizeof(float);
				float64x2_t a1 = vcvt_f64_f32(vget_low_f32(v1));
				float64x2_t a2 = vcvt_f64_f32(vget_high_f32(v1));

				s1[0 + (kk << 1)] = vaddq_f64(s1[0 + (kk << 1)], a1);
				s1[1 + (kk << 1)] = vaddq_f64(s1[1 + (kk << 1)], a2);

				a1 = vmulq_f64(a1, a1);
				a2 = vmulq_f64(a2, a2);

				s2[0 + (kk << 1)] = vaddq_f64(s2[0 + (kk << 1)], a1);
				s2[1 + (kk << 1)] = vaddq_f64(s2[1 + (kk << 1)], a2);
			}
		}

		for (int KK = sizeof(s1) / sizeof(s1[0]) / 2; KK >= 1; KK >>= 1)
		{
			REP(KK) s1[kk] = vaddq_f64(s1[kk], s1[kk + KK]);
			REP(KK) s2[kk] = vaddq_f64(s2[kk], s2[kk + KK]);
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
	constexpr int N = 2;
	int64 i = 0;
	volatile double re1 = 0, re2 = 0;

	if (n >= N * sizeof(float32x4_t) / sizeof(float))
	{
		REP(8) { __builtin_prefetch(B, 0); B += sep; }

		float64x2_t s1[2], s2[2], b;
		REP(2) s1[kk] = s2[kk] = vdupq_n_f64(0);

		for (int64 l1 = n - N * sizeof(float32x4_t) / sizeof(float); i <= l1; i += N * sizeof(float32x4_t) / sizeof(float))
		{
			REP(N)
			{
				__builtin_prefetch(B, 0); B += sep;
				__builtin_prefetch(B, 0); B += sep;
				b = vld1q_f64(((double[]) { B[-10 * sep], B[-9 * sep] }));
				s1[0] = vaddq_f64(s1[0], vmulq_f64(vld1q_f64(A1), b)); A1 += 2;
				s2[0] = vaddq_f64(s2[0], vmulq_f64(vld1q_f64(A2), b)); A2 += 2;

				__builtin_prefetch(B, 0); B += sep;
				__builtin_prefetch(B, 0); B += sep;
				b = vld1q_f64(((double[]) { B[-10 * sep], B[-9 * sep] }));
				s1[1] = vaddq_f64(s1[1], vmulq_f64(vld1q_f64(A1), b)); A1 += 2;
				s2[1] = vaddq_f64(s2[1], vmulq_f64(vld1q_f64(A2), b)); A2 += 2;
			}
		}

		s1[0] = vaddq_f64(s1[0], s1[1]);
		s2[0] = vaddq_f64(s2[0], s2[1]);

		re1 = _neo_reduce_add_pd(s1[0]);
		re2 = _neo_reduce_add_pd(s2[0]);

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

TARGETNEO double SumProdDivNEO(double* A1, float* A2, float* B, int64 sep, int64 n)
{
	constexpr int N = 2;
	int64 i = 0;
	double re1 = 0, re2 = 0;

	if (n >= N * sizeof(float32x4_t) / sizeof(float))
	{
		REP(6) { __builtin_prefetch(B, 0); B += sep; }

		float32x4_t a2;
		float64x2_t a1[N * 2], b[N * 2], s1[N * 2], s2[N * 2];
		REP(N * 2) s1[kk] = s2[kk] = vdupq_n_f64(0);

		for (int64 l1 = n - N * sizeof(float32x4_t) / sizeof(float); i <= l1; i += N * sizeof(float32x4_t) / sizeof(float))
		{
			REP(N)
			{
				__builtin_prefetch(B, 0); B += sep; 
				__builtin_prefetch(B, 0); B += sep;
				a1[0 + (kk << 1)] = vld1q_f64(A1); A1 += 2;
				a2 = vld1q_f32(A2); A2 += 4;
				b[0 + (kk << 1)] = vld1q_f64(((double[]) { B[-8 * sep], B[-7 * sep] }));

				__builtin_prefetch(B, 0); B += sep;
				__builtin_prefetch(B, 0); B += sep;
				a1[1 + (kk << 1)] = vld1q_f64(A1); A1 += 2;
				//a2
				b[1 + (kk << 1)] = vld1q_f64(((double[]) { B[-8 * sep], B[-7 * sep] }));

				s1[0 + (kk << 1)] = vaddq_f64(s1[0 + (kk << 1)], vmulq_f64(a1[0 + (kk << 1)], b[0 + (kk << 1)]));
				s2[0 + (kk << 1)] = vaddq_f64(s2[0 + (kk << 1)], vmulq_f64(vcvt_f64_f32(vget_low_f32 (a2)), b[0 + (kk << 1)]));

				s1[1 + (kk << 1)] = vaddq_f64(s1[1 + (kk << 1)], vmulq_f64(a1[1 + (kk << 1)], b[1 + (kk << 1)]));
				s2[1 + (kk << 1)] = vaddq_f64(s2[1 + (kk << 1)], vmulq_f64(vcvt_f64_f32(vget_high_f32(a2)), b[1 + (kk << 1)]));
			}
		}

		for (int KK = sizeof(s1) / sizeof(s1[0]) / 2; KK >= 1; KK >>= 1)
		{
			REP(KK) s1[kk] = vaddq_f64(s1[kk], s1[kk + KK]);
			REP(KK) s2[kk] = vaddq_f64(s2[kk], s2[kk + KK]);
		}

		re1 = _neo_reduce_add_pd(s1[0]);
		re2 = _neo_reduce_add_pd(s2[0]);

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

TARGETNEO double SumProdDivNEO(float* A1, float* A2, float* B, int64 sep, int64 n)
{
	constexpr int N = 2;
	int64 i = 0;
	volatile double re1 = 0, re2 = 0;

	if (n >= N * sizeof(float32x4_t) / sizeof(float))
	{
		REP(6) { __builtin_prefetch(B, 0); B += sep; }

		float64x2_t s1[N * 2], s2[N * 2], b[N * 2];
		REP(N * 2) s1[kk] = s2[kk] = vdupq_n_f64(0);

		for (int64 l1 = n - N * sizeof(float32x4_t) / sizeof(float); i <= l1; i += N * sizeof(float32x4_t) / sizeof(float))
		{
			REP(N)
			{
				__builtin_prefetch(B, 0); B += sep;
				__builtin_prefetch(B, 0); B += sep;

				b[0 + (kk << 1)] = vld1q_f64(((double[]) { B[-8 * sep], B[-7 * sep] }));

				float32x4_t a1 = vld1q_f32(A1); A1 += sizeof(float32x4_t) / sizeof(float);
				float32x4_t a2 = vld1q_f32(A2); A2 += sizeof(float32x4_t) / sizeof(float);

				__builtin_prefetch(B, 0); B += sep;
				__builtin_prefetch(B, 0); B += sep; 
				b[1 + (kk << 1)] = vld1q_f64(((double[]) { B[-8 * sep], B[-7 * sep] }));

				//a1
				//a2

				s1[0 + (kk << 1)] = vaddq_f64(s1[0 + (kk << 1)], vmulq_f64(vcvt_f64_f32(vget_low_f32 (a1)), b[0 + (kk << 1)]));
				s2[0 + (kk << 1)] = vaddq_f64(s2[0 + (kk << 1)], vmulq_f64(vcvt_f64_f32(vget_low_f32 (a2)), b[0 + (kk << 1)]));
				s1[1 + (kk << 1)] = vaddq_f64(s1[1 + (kk << 1)], vmulq_f64(vcvt_f64_f32(vget_high_f32(a1)), b[1 + (kk << 1)]));
				s2[1 + (kk << 1)] = vaddq_f64(s2[1 + (kk << 1)], vmulq_f64(vcvt_f64_f32(vget_high_f32(a2)), b[1 + (kk << 1)]));
			}
		}

		for (int KK = sizeof(s1) / sizeof(s1[0]) / 2; KK >= 1; KK >>= 1)
		{
			REP(KK) s1[kk] = vaddq_f64(s1[kk], s1[kk + KK]);
			REP(KK) s2[kk] = vaddq_f64(s2[kk], s2[kk + KK]);
		}

		re1 = _neo_reduce_add_pd(s1[0]);
		re2 = _neo_reduce_add_pd(s2[0]);

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

TARGETNEO float SumProdDivNEOx(float* A1, float* A2, float* B, int64 sep, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;
	volatile float re1 = 0, re2 = 0;

	if (n >= N * sizeof(float32x4_t) / sizeof(float))
	{
		float32x4_t s1[N], s2[N], b;
		REP(N) s1[kk] = s2[kk] = vdupq_n_f32(0);

		for (int64 l1 = n - N * sizeof(float32x4_t) / sizeof(float); i <= l1; i += N * sizeof(float32x4_t) / sizeof(float))
		{
			REP(N)
			{
				b = vld1q_f32(((float[]) { B[0 * sep], B[1 * sep], B[2 * sep], B[3 * sep] })); B += sizeof(float32x4_t) / sizeof(float) * sep;
				s1[kk] = vaddq_f32(s1[kk], vmulq_f32(vld1q_f32(A1), b)); A1 += sizeof(float32x4_t) / sizeof(float);
				s2[kk] = vaddq_f32(s2[kk], vmulq_f32(vld1q_f32(A2), b)); A2 += sizeof(float32x4_t) / sizeof(float);
			}
		}

		for (int KK = sizeof(s1) / sizeof(s1[0]) / 2; KK >= 1; KK >>= 1)
		{
			REP(KK) s1[kk] = vaddq_f32(s1[kk], s1[kk + KK]);
			REP(KK) s2[kk] = vaddq_f32(s2[kk], s2[kk + KK]);
		}

		re1 = _neo_reduce_add_ps(s1[0]);
		re2 = _neo_reduce_add_ps(s2[0]);
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

TARGETNEO double SumProdNEO(double* A, double* B, int64 sep, int64 n)
{
	int64 i = 0;
	double re = 0;

	if (n >= 2)
	{
		__builtin_prefetch(B, 0); B += sep;
		__builtin_prefetch(B, 0); B += sep;

		float64x2_t s = vdupq_n_f64(0);

		for (int64 l1 = n - 2; i <= l1; i += 2)
		{
			__builtin_prefetch(B, 0); B += sep;
			__builtin_prefetch(B, 0); B += sep;
			s = vaddq_f64(s, vmulq_f64(vld1q_f64(A), vld1q_f64(((double[]) { B[-4 * sep], B[-3 * sep] }))));
			A += 2;
		}

		re = _neo_reduce_add_pd(s);
		B -= 2 * sep;
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
	int64 i = 0;
	double re = 0;

	if (n >= 4)
	{
		__builtin_prefetch(B, 0); B += sep;
		__builtin_prefetch(B, 0); B += sep;
		__builtin_prefetch(B, 0); B += sep;
		__builtin_prefetch(B, 0); B += sep;

		float64x2_t s1 = vdupq_n_f64(0), s2 = vdupq_n_u64(0);
		float32x4_t a;

		for (int64 l1 = n - 4; i <= l1; i += 4)
		{
			__builtin_prefetch(B, 0); B += sep;
			__builtin_prefetch(B, 0); B += sep;
			a = vld1q_f32(A); A += 4;
			s1 = vaddq_f64(s1, vmulq_f64(vcvt_f64_f32(vget_low_f32 (a)), vld1q_f64(((double[]) { B[-6 * sep], B[-5 * sep] }))));

			__builtin_prefetch(B, 0); B += sep;
			__builtin_prefetch(B, 0); B += sep;
			//a
			s2 = vaddq_f64(s2, vmulq_f64(vcvt_f64_f32(vget_high_f32(a)), vld1q_f64(((double[]) { B[-6 * sep], B[-5 * sep] }))));
		}

		s1 = vaddq_f64(s1, s2);
		re = _neo_reduce_add_pd(s1);

		B -= 4 * sep;
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
	constexpr int N = 4;
	int64 i = 0;
	volatile float re = 0;

	if (n >= N * sizeof(float32x4_t) / sizeof(float))
	{
		REP(N) { __builtin_prefetch(B, 0); B += sep; }

		float32x4_t s[N];
		REP(N) s[kk] = vdupq_n_f32(0);

		for (int64 l1 = n - N * sizeof(float32x4_t) / sizeof(float); i <= l1; i += N * sizeof(float32x4_t) / sizeof(float))
		{
			REP(N)
			{
				__builtin_prefetch(B, 0); B += sep;
				__builtin_prefetch(B, 0); B += sep;
				__builtin_prefetch(B, 0); B += sep;
				__builtin_prefetch(B, 0); B += sep;

				s[kk] = vaddq_f32(s[kk], vmulq_f32(vld1q_f32(A), vld1q_f32(((float[]) { B[-8 * sep], B[-7 * sep], B[-6 * sep], B[-5 * sep] }))));
				A += sizeof(float32x4_t) / sizeof(float);
			}
		}

		for (int KK = sizeof(s) / sizeof(s[0]) / 2; KK >= 1; KK >>= 1)
			REP(KK) s[kk] = vaddq_f32(s[kk], s[kk + KK]);

		re = _neo_reduce_add_ps(s[0]);

		B -= 4 * sep;
	}

	for (; i < n; ++i, A++, B += sep)
	{
		volatile float v1 = *A * *B;
		re += v1;
	}

	return re;
}

TARGETNEO double SumProdNEO(double* A, double* B, int64 n)
{
	constexpr int N = 16;
	int64 i = 0;
	double re = 0;

	if (n >= N * sizeof(float64x2_t) / sizeof(double))
	{
		float64x2_t s[N];
		REP(N) s[kk] = vdupq_n_f64(0);

		for (int64 l1 = n - N * sizeof(float64x2_t) / sizeof(double); i <= l1; i += N * sizeof(float64x2_t) / sizeof(double))
		{
			REP(N)
			{
				s[kk] = vaddq_f64(s[kk], vmulq_f64(vld1q_f64(A), vld1q_f64(B)));
				A += sizeof(float64x2_t) / sizeof(double);
				B += sizeof(float64x2_t) / sizeof(double);
			}
		}

		for (int KK = sizeof(s) / sizeof(s[0]) / 2; KK >= 1; KK >>= 1)
			REP(KK) s[kk] = vaddq_f64(s[kk], s[kk + KK]);

		re = _neo_reduce_add_pd(s[0]);
	}

	for (; i < n; ++i, ++A, ++B)
	{
		volatile double v1 = *A * *B;
		re += v1;
	}

	return re;
}

TARGETNEO double SumProdNEO(float* A, float* B, int64 n)
{
	constexpr int N = 8;
	int64 i = 0;
	double re = 0;

	if (n >= N / 2 * sizeof(float32x4_t) / sizeof(float))
	{
		float64x2_t s[N];
		REP(N) s[kk] = vdupq_n_f64(0);

		for (int64 l1 = n - N / 2 * sizeof(float32x4_t) / sizeof(float); i <= l1; i += N / 2 * sizeof(float32x4_t) / sizeof(float))
		{
			REP(N / 2)
			{
				float32x4_t v1 = vld1q_f32(A), v2 = vld1q_f32(B);
				A += sizeof(float32x4_t) / sizeof(float);
				B += sizeof(float32x4_t) / sizeof(float);

				s[0 + (kk << 1)] = vaddq_f64(s[0 + (kk << 1)], vmulq_f64(vcvt_f64_f32(vget_low_f32 (v1)), vcvt_f64_f32(vget_low_f32 (v2))));

				s[1 + (kk << 1)] = vaddq_f64(s[1 + (kk << 1)], vmulq_f64(vcvt_f64_f32(vget_high_f32(v1)), vcvt_f64_f32(vget_high_f32(v2))));
			}
		}

		for (int KK = sizeof(s) / sizeof(s[0]) / 2; KK >= 1; KK >>= 1)
			REP(KK) s[kk] = vaddq_f64(s[kk], s[kk + KK]);

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
	constexpr int N = 8;
	int64 i = 0;
	float re = 0;

	if (n >= N * sizeof(float32x4_t) / sizeof(float))
	{
		float32x4_t s[N];
		REP(N) s[kk] = vdupq_n_f32(0);

		for (int64 l1 = n - N * sizeof(float32x4_t) / sizeof(float); i <= l1; i += N * sizeof(float32x4_t) / sizeof(float))
		{
			REP(N)
			{
				s[kk] = vaddq_f32(s[kk], vmulq_f32(vld1q_f32(A), vld1q_f32(B)));
				A += sizeof(float32x4_t) / sizeof(float);
				B += sizeof(float32x4_t) / sizeof(float);
			}
		}

		for (int KK = sizeof(s) / sizeof(s[0]) / 2; KK >= 1; KK >>= 1)
			REP(KK) s[kk] = vaddq_f32(s[kk], s[kk + KK]);

		re = _neo_reduce_add_ps(s[0]);
	}

	for (; i < n; ++i, ++A, ++B)
	{
		volatile float v1 = *A * *B;
		re += v1;
	}

	return re;
}

TARGETNEO void AddNEO(double* A, double* B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(float64x2_t) / sizeof(double))
	{
		float64x2_t a[N], b[N];

		for (int64 l1 = n - N * sizeof(float64x2_t) / sizeof(double); i <= l1; i += N * sizeof(float64x2_t) / sizeof(double))
		{
			REP(4) { a[kk] = vld1q_f64(A); A += sizeof(float64x2_t) / sizeof(double); }

			REP(4) { b[kk] = vld1q_f64(B); B += sizeof(float64x2_t) / sizeof(double); }

			REP(4) a[kk] = vaddq_f64(a[kk], b[kk]);

			REP(4) vst1q_f64(A + (kk - N) * sizeof(float64x2_t) / sizeof(double), a[kk]);
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGETNEO void AddNEO(float* A, float* B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(float32x4_t) / sizeof(float))
	{
		float32x4_t a[N], b[N];

		for (int64 l1 = n - N * sizeof(float32x4_t) / sizeof(float); i <= l1; i += N * sizeof(float32x4_t) / sizeof(float))
		{
			REP(N) { a[kk] = vld1q_f32(A); A += sizeof(float32x4_t) / sizeof(float); }

			REP(N) { b[kk] = vld1q_f32(B); B += sizeof(float32x4_t) / sizeof(float); }

			REP(N) a[kk] = vaddq_f32(a[kk], b[kk]);

			REP(N) vst1q_f32(A + (kk - N) * sizeof(float32x4_t) / sizeof(float), a[kk]);
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGETNEO void AddNEO(int64* A, int64* B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(int64x2_t) / sizeof(int64))
	{
		int64x2_t a[N], b[N];

		for (int64 l1 = n - N * sizeof(int64x2_t) / sizeof(int64); i <= l1; i += N * sizeof(int64x2_t) / sizeof(int64))
		{
			REP(N) { a[kk] = vld1q_s64(A); A += sizeof(int64x2_t) / sizeof(int64); }

			REP(N) { b[kk] = vld1q_s64(B); B += sizeof(int64x2_t) / sizeof(int64); }

			REP(N) a[kk] = vaddq_s64(a[kk], b[kk]);

			REP(N) vst1q_s64(A + (kk - N) * sizeof(int64x2_t) / sizeof(int64), a[kk]);
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGETNEO void AddNEO(int* A, int* B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(int32x4_t) / sizeof(int))
	{
		int32x4_t a[N], b[N];

		for (int64 l1 = n - N * sizeof(int32x4_t) / sizeof(int); i <= l1; i += N * sizeof(int32x4_t) / sizeof(int))
		{
			REP(N) { a[kk] = vld1q_s32(A); A += sizeof(int32x4_t) / sizeof(int); }

			REP(N) { b[kk] = vld1q_s32(B); B += sizeof(int32x4_t) / sizeof(int); }

			REP(N) a[kk] = vaddq_s32(a[kk], b[kk]);

			REP(N) vst1q_s32(A + (kk - N) * sizeof(__m128i) / sizeof(int), a[kk]);
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGETNEO void AddNEO(int* A, int B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(int32x4_t) / sizeof(int))
	{
		int32x4_t a[N], b = vdupq_n_s32(B);

		for (int64 l1 = n - N * sizeof(int32x4_t) / sizeof(int); i <= l1; i += N * sizeof(int32x4_t) / sizeof(int))
		{
			REP(N) { a[kk] = vld1q_s32(A); A += sizeof(int32x4_t) / sizeof(int); }

			REP(N) a[kk] = vaddq_s32(a[kk], b);

			REP(N) vst1q_s32(A + (kk - N) * sizeof(__m128i) / sizeof(int), a[kk]);
		}
	}

	for (; i < n; ++i, A++)
		*A += B;
}

TARGETNEO void AddNEO(double* A, double B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(float64x2_t) / sizeof(double))
	{
		float64x2_t b = vdupq_n_f64(B), a[N];

		for (int64 l1 = n - N * sizeof(float64x2_t) / sizeof(double); i <= l1; i += N * sizeof(float64x2_t) / sizeof(double))
		{
			REP(N) { a[kk] = vld1q_f64(A); A += sizeof(float64x2_t) / sizeof(double); }

			REP(N) a[kk] = vaddq_f64(a[kk], b);

			REP(N) vst1q_f64(A + (kk - N) * sizeof(float64x2_t) / sizeof(double), a[kk]);
		}
	}

	for (; i < n; ++i, A++)
		*A += B;
}

TARGETNEO void AddNEO(float* A, float B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(float32x4_t) / sizeof(float))
	{
		float32x4_t b = vdupq_n_f32(B), a[N];

		for (int64 l1 = n - N * sizeof(float32x4_t) / sizeof(float); i <= l1; i += N * sizeof(float32x4_t) / sizeof(float))
		{
			REP(N) { a[kk] = vld1q_f32(A); A += sizeof(float32x4_t) / sizeof(float); }

			REP(N) a[kk] = vaddq_f32(a[kk], b);

			REP(N) vst1q_f32(A + (kk - N) * sizeof(float32x4_t) / sizeof(float), a[kk]);
		}
	}

	for (; i < n; ++i, A++)
		*A += B;
}

TARGETNEO void MulNEO(double* C, double* A, double* B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(float64x2_t) / sizeof(double))
	{
		float64x2_t a[N], b[N];

		for (int64 l1 = n - N * sizeof(float64x2_t) / sizeof(double); i <= l1; i += N * sizeof(float64x2_t) / sizeof(double))
		{
			REP(N) { a[kk] = vld1q_f64(A); A += sizeof(float64x2_t) / sizeof(double); }

			REP(N) { b[kk] = vld1q_f64(B); B += sizeof(float64x2_t) / sizeof(double); }

			REP(N) a[kk] = vmulq_f64(a[kk], b[kk]);

			REP(N) { vst1q_f64(C, a[kk]); C += sizeof(float64x2_t) / sizeof(double); }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *A++ * *B++;
		*C++ = v1;
	}
}

TARGETNEO void MulNEO(float* C, float* A, float* B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(float32x4_t) / sizeof(float))
	{
		float32x4_t a[N], b[N];

		for (int64 l1 = n - N * sizeof(float32x4_t) / sizeof(float); i <= l1; i += N * sizeof(float32x4_t) / sizeof(float))
		{
			REP(N) { a[kk] = vld1q_f32(A); A += sizeof(float32x4_t) / sizeof(float); }

			REP(N) { b[kk] = vld1q_f32(B); B += sizeof(float32x4_t) / sizeof(float); }

			REP(N) a[kk] = vmulq_f32(a[kk], b[kk]);

			REP(N) { vst1q_f32(C, a[kk]); C += sizeof(float32x4_t) / sizeof(float); }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A++ * *B++;
		*C++ = v1;
	}
}

TARGETNEO void MulNEO(double* C, double* A, double B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(float64x2_t) / sizeof(double))
	{
		float64x2_t b = vdupq_n_f64(B), a[N];

		for (int64 l1 = n - N * sizeof(float64x2_t) / sizeof(double); i <= l1; i += N * sizeof(float64x2_t) / sizeof(double))
		{
			REP(N) { a[kk] = vld1q_f64(A); A += sizeof(float64x2_t) / sizeof(double); }

			REP(N) a[kk] = vmulq_f64(a[kk], b);

			REP(N) { vst1q_f64(C, a[kk]); C += sizeof(float64x2_t) / sizeof(double); }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *A++ * B;
		*C++ = v1;
	}
}

TARGETNEO void MulNEO(float* C, float* A, float B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(float32x4_t) / sizeof(float))
	{
		float32x4_t b = vdupq_n_f32(B), a[N];

		for (int64 l1 = n - N * sizeof(float32x4_t) / sizeof(float); i <= l1; i += N * sizeof(float32x4_t) / sizeof(float))
		{
			REP(N) { a[kk] = vld1q_f32(A); A += sizeof(float32x4_t) / sizeof(float); }

			REP(N) a[kk] = vmulq_f32(a[kk], b);

			REP(N) { vst1q_f32(C, a[kk]); C += sizeof(float32x4_t) / sizeof(float); }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A++ * B;
		*C++ = v1;
	}
}

TARGETNEO void MulNEO(double* A, double B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(float64x2_t) / sizeof(double))
	{
		float64x2_t b = vdupq_n_f64(B), a[N];

		for (int64 l1 = n - N * sizeof(float64x2_t) / sizeof(double); i <= l1; i += N * sizeof(float64x2_t) / sizeof(double))
		{
			REP(N) { a[kk] = vld1q_f64(A); A += sizeof(float64x2_t) / sizeof(double); }

			REP(N) a[kk] = vmulq_f64(a[kk], b);

			REP(N) vst1q_f64(A + (kk - N) * sizeof(float64x2_t) / sizeof(double), a[kk]);
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
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(float32x4_t) / sizeof(float))
	{
		float32x4_t b = vdupq_n_f32(B), a[N];

		for (int64 l1 = n - N * sizeof(float32x4_t) / sizeof(float); i <= l1; i += N * sizeof(float32x4_t) / sizeof(float))
		{
			REP(N) { a[kk] = vld1q_f32(A); A += sizeof(float32x4_t) / sizeof(float); }

			REP(N) a[kk] = vmulq_f32(a[kk], b);

			REP(N) vst1q_f32(A + (kk - N) * sizeof(float32x4_t) / sizeof(float), a[kk]);
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A * B;
		*A++ = v1;
	}
}

TARGETNEO void AddProdNEO(double* C, double* A, double* B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= N * sizeof(float64x2_t) / sizeof(double))
	{
		float64x2_t a[N];

		for (int64 l1 = n - N * sizeof(float64x2_t) / sizeof(double); i <= l1; i += N * sizeof(float64x2_t) / sizeof(double))
		{
			REP(N) { a[kk] = vld1q_f64(A); A += sizeof(float64x2_t) / sizeof(double); }

			REP(N) { a[kk] = vmulq_f64(a[kk], vld1q_f64(B)); B += sizeof(float64x2_t) / sizeof(double); }

			double* C2 = C;

			REP(N) { a[kk] = vaddq_f64(a[kk], vld1q_f64(C)); C += sizeof(float64x2_t) / sizeof(double); }

			REP(N) { vst1q_f64(C2, a[kk]); C2 += sizeof(float64x2_t) / sizeof(double); }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *A++ * *B++;
		*C++ += v1;
	}
}

TARGETNEO void AddProdNEO(float* C, float* A, float* B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= 4 * N * sizeof(float32x4_t) / sizeof(float))
	{
		float32x4_t a[N];

		for (int64 l1 = n - N * sizeof(float32x4_t) / sizeof(float); i <= l1; i += N * sizeof(float32x4_t) / sizeof(float))
		{
			REP(N) { a[kk] = vld1q_f32(A); A += sizeof(float32x4_t) / sizeof(float); }

			REP(N) { a[kk] = vmulq_f32(a[kk], vld1q_f32(B)); B += sizeof(float32x4_t) / sizeof(float); }

			float* C2 = C;

			REP(N) { a[kk] = vaddq_f32(a[kk], vld1q_f32(C)); C += sizeof(float32x4_t) / sizeof(float); }

			REP(N) { vst1q_f32(C2, a[kk]); C2 += sizeof(float32x4_t) / sizeof(float); }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A++ * *B++;
		*C++ += v1;
	}
}

TARGETNEO void AddProdNEO(double* C, double* A, double B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= 4 * N * sizeof(float64x2_t) / sizeof(double))
	{
		float64x2_t a[N], b = vdupq_n_f64(B);

		for (int64 l1 = n - N * sizeof(float64x2_t) / sizeof(double); i <= l1; i += N * sizeof(float64x2_t) / sizeof(double))
		{
			REP(N) { a[kk] = vld1q_f64(A); A += sizeof(float64x2_t) / sizeof(double); }

			REP(N)  a[kk] = vmulq_f64(a[kk], b);

			double* C2 = C;

			REP(N) { a[kk] = vaddq_f64(a[kk], vld1q_f64(C)); C += sizeof(float64x2_t) / sizeof(double); }

			REP(N) { vst1q_f64(C2, a[kk]); C2 += sizeof(float64x2_t) / sizeof(double); }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *A++ * B;
		*C++ += v1;
	}
}

TARGETNEO void AddProdNEO(double* C, float* A, double B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= 4 * N * sizeof(float64x2_t) / sizeof(double))
	{
		float64x2_t a[N], b = vdupq_n_f64(B);

		for (int64 l1 = n - N * sizeof(float64x2_t) / sizeof(double); i <= l1; i += N * sizeof(float64x2_t) / sizeof(double))
		{
			REP(N / 2)
			{
				float32x4_t v1 = vld1q_f32(A); A += sizeof(float32x4_t) / sizeof(float);

				a[0 + (kk << 1)] = vcvt_f64_f32(vget_low_f32 (v1));
				a[1 + (kk << 1)] = vcvt_f64_f32(vget_high_f32(v1));
			}

			REP(N)  a[kk] = vmulq_f64(a[kk], b);

			double* C2 = C;

			REP(N) { a[kk] = vaddq_f64(a[kk], vld1q_f64(C)); C += sizeof(float64x2_t) / sizeof(double); }

			REP(N) { vst1q_f64(C2, a[kk]); C2 += sizeof(float64x2_t) / sizeof(double); }
		}
	}

	for (; i < n; ++i)
	{
		volatile double v1 = *A++ * B;
		*C++ += v1;
	}
}

TARGETNEO void AddProdNEO(float* C, float* A, float B, int64 n)
{
	constexpr int N = 4;
	int64 i = 0;

	if (n >= 4 * N * sizeof(float32x4_t) / sizeof(float))
	{
		float32x4_t a[N], b = vdupq_n_f32(B);

		for (int64 l1 = n - N * sizeof(float32x4_t) / sizeof(float); i <= l1; i += N * sizeof(float32x4_t) / sizeof(float))
		{
			REP(N) { a[kk] = vld1q_f32(A); A += sizeof(float32x4_t) / sizeof(float); }

			REP(N)  a[kk] = vmulq_f32(a[kk], b);

			float* C2 = C;

			REP(N) { a[kk] = vaddq_f32(a[kk], vld1q_f32(C)); C += sizeof(float32x4_t) / sizeof(float); }

			REP(N) { vst1q_f32(C2, a[kk]); C2 += sizeof(float32x4_t) / sizeof(float); }
		}
	}

	for (; i < n; ++i)
	{
		volatile float v1 = *A++ * B;
		*C++ += v1;
	}
}

TARGETNEO void UnifyNEO(double* A, int64 n)
{
	constexpr int N = 16;
	int64 i = 0;
	double invsum = 1.0 / (SumNEO(A, n) + n * MIN_FREQ);

	if (n >= N * sizeof(float32x4_t) / sizeof(double))
	{
		float64x2_t a[N], minf = vdupq_n_f64(MIN_FREQ), invs = vdupq_n_f64(invsum);

		for (int64 l1 = n - N * sizeof(float32x4_t) / sizeof(double); i <= l1; i += N * sizeof(float32x4_t) / sizeof(double))
		{
			double* A2 = A;

			REP(N) { a[kk] = vld1q_f64(A); A += sizeof(float32x4_t) / sizeof(double); }

			REP(N) a[kk] = vaddq_f64(a[kk], minf);

			REP(N) a[kk] = vmulq_f64(a[kk], invs);

			REP(N) { vst1q_f64(A2, a[kk]); A2 += sizeof(float32x4_t) / sizeof(double); }
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
	constexpr int N = 8;
	int64 i = 0;
	double invsum = 1.0 / (SumNEO(A, n) + n * MIN_FREQ);

	if (n >= N / 2 * sizeof(float32x4_t) / sizeof(float))
	{
		float64x2_t a[N], minf = vdupq_n_f64(MIN_FREQ), invs = vdupq_n_f64(invsum);

		for (int64 l1 = n - N / 2 * sizeof(float32x4_t) / sizeof(float); i <= l1; i += N / 2 * sizeof(float32x4_t) / sizeof(float))
		{
			float* A2 = A;

			REP(N / 2) 
			{
				float32x4_t aa = vld1q_f32(A); A += sizeof(float32x4_t) / sizeof(float);
				a[0 + (kk << 1)] = vcvt_f64_f32(vget_low_f32 (aa));
				a[1 + (kk << 1)] = vcvt_f64_f32(vget_high_f32(aa));
			}

			REP(N) a[kk] = vaddq_f64(a[kk], minf);

			REP(N) a[kk] = vmulq_f64(a[kk], invs);

			REP(N / 2)
			{
				vst1q_f32(A2, vcombine_f32(
					vcvt_f32_f64(a[0 + (kk << 1)]),
					vcvt_f32_f64(a[1 + (kk << 1)])));

				A2 += sizeof(float32x4_t) / sizeof(float);
			}
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
	A++; n--;
	int64 i = 0;

	if (n >= 16)
	{
		uint8x16_t r = vdupq_n_u8(0), v = vdupq_n_u8(val), o = vdupq_n_u8(0), a;
		uint64 s1, s2;

		for (int64 l1 = n - 16; i <= l1; )
		{
			a = vld1q_u8((byte*)A);

			r = vceqq_u8(a, v);
			
			if (!(vgetq_lane_u64((uint64x2_t)r, 0) | vgetq_lane_u64((uint64x2_t)r, 1))) { A += 16; i += 16; continue; }
			
			r = vsubq_u8(o, r);

			//r = _mm_sad_epu8(o, r);
			s1 = __builtin_popcountll(vgetq_lane_u64((uint64x2_t)r, 0));
			s2 = __builtin_popcountll(vgetq_lane_u64((uint64x2_t)r, 1));

			if (rep > s1) 
				rep -= s1;
			else for (;; A++) 
				if (*A == val && !--rep) return A;
			A += 8;

			if (rep > s2) 
				rep -= s2;
			else for (;; A++) 
				if (*A == val && !--rep) return A;
			A += 8;
		}
	}

	for (; i < n; ++i, A++)
		if (*A == val && !--rep)
			return A;

	return NULL;
}

TARGETNEO int64 CountCharNEO(char* A, char val, int64 n)
{
	uint64 re = 0;
	int64 i = 0;

	if (n >= 16)
	{
		uint8x16_t v = vdupq_n_u8(val), o = vdupq_n_u8(0), a;
		uint64 s1 = 0, s2 = 0;

		for (int64 l1 = n - 16; i <= l1; i += 16)
		{
			a = vld1q_u8((byte*)A);
			A += 16;
			a = vceqq_u8(a, v);
			a = vsubq_u8(o, a);
			s1 += __builtin_popcountll(vgetq_lane_u64((uint64x2_t)a, 0));
			s2 += __builtin_popcountll(vgetq_lane_u64((uint64x2_t)a, 1));
		}

		re = s1 + s2;
	}

	for (; i < n; ++i, A++)
		if (*A == val) re++;

	return (int64)re;
}

#endif
