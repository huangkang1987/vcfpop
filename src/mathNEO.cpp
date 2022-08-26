/* NEON Instruction Set Functions */

#include "vcfpop.h"

#ifdef __aarch64__

#ifndef _RNGNEO
/* Initialize rng */
TARGETNEO RNGNEO::RNGNEO()
{
}

/* Initialize rng */
TARGETNEO RNGNEO::RNGNEO(uint64 s)
{
	uint64 ss[8] = {
		Hash64ULong(s + 0), Hash64ULong(s + 1),
		Hash64ULong(s + 2), Hash64ULong(s + 3),
		Hash64ULong(s + 4), Hash64ULong(s + 5),
		Hash64ULong(s + 6), Hash64ULong(s + 7),
	};

	uint64x2_t* seed = (uint64x2_t*)ss;

	uint64x2_t a1 = vdupq_n_u64(0x159A55E5075BCD15);
	uint64x2_t a2 = vdupq_n_u64(0x054913331F123BB5);

	x[0] = veorq_u64(a1, seed[0]);
	x[1] = veorq_u64(a1, seed[1]);
	x[2] = veorq_u64(a1, seed[2]);
	x[3] = veorq_u64(a1, seed[3]);

	seed[0] = vshlq_n_u64(seed[0], 6);
	seed[1] = vshlq_n_u64(seed[1], 6);
	seed[2] = vshlq_n_u64(seed[2], 6);
	seed[3] = vshlq_n_u64(seed[3], 6);

	y[0] = veorq_u64(a2, seed[0]);
	y[1] = veorq_u64(a2, seed[1]);
	y[2] = veorq_u64(a2, seed[2]);
	y[3] = veorq_u64(a2, seed[3]);
}

/* Draw a uniform distriubted interger */
TARGETNEO void RNGNEO::XorShift128p(uint64x2_t* re)
{
	uint64x2_t a[4], b[4];

	a[0] = x[0];
	a[1] = x[1];
	a[2] = x[2];
	a[3] = x[3];

	b[0] = y[0];
	b[1] = y[1];
	b[2] = y[2];
	b[3] = y[3];

	x[0] = b[0];
	x[1] = b[1];
	x[2] = b[2];
	x[3] = b[3];

	a[0] = veorq_u64(a[0], vshlq_n_u64(a[0], 23));
	a[1] = veorq_u64(a[1], vshlq_n_u64(a[1], 23));
	a[2] = veorq_u64(a[2], vshlq_n_u64(a[2], 23));
	a[3] = veorq_u64(a[3], vshlq_n_u64(a[3], 23));

	a[0] = veorq_u64(a[0], vshrq_n_u64(a[0], 18));
	a[1] = veorq_u64(a[1], vshrq_n_u64(a[1], 18));
	a[2] = veorq_u64(a[2], vshrq_n_u64(a[2], 18));
	a[3] = veorq_u64(a[3], vshrq_n_u64(a[3], 18));

	a[0] = veorq_u64(a[0], b[0]);
	a[1] = veorq_u64(a[1], b[1]);
	a[2] = veorq_u64(a[2], b[2]);
	a[3] = veorq_u64(a[3], b[3]);

	a[0] = veorq_u64(a[0], vshrq_n_u64(b[0], 5));
	a[1] = veorq_u64(a[1], vshrq_n_u64(b[1], 5));
	a[2] = veorq_u64(a[2], vshrq_n_u64(b[2], 5));
	a[3] = veorq_u64(a[3], vshrq_n_u64(b[3], 5));

	y[0] = a[0];
	y[1] = a[1];
	y[2] = a[2];
	y[3] = a[3];

	re[0] = vaddq_s64(a[0], b[0]);
	re[1] = vaddq_s64(a[1], b[1]);
	re[2] = vaddq_s64(a[2], b[2]);
	re[3] = vaddq_s64(a[3], b[3]);
}

/* Draw a uniform distriubted real number */
TARGETNEO void RNGNEO::Uniform(float64x2_t* re)
{
	float64x2_t one = vdupq_n_f64(1.0);
	uint64x2_t mask1 = vdupq_n_u64(0x000FFFFFFFFFFFFF);
	uint64x2_t mask2 = vdupq_n_u64(0x3FF0000000000000);
	uint64x2_t* r = (uint64x2_t*)re;

	XorShift128p(r);

	r[0] = vorrq_u64(vandq_u64(r[0], mask1), mask2);
	r[1] = vorrq_u64(vandq_u64(r[1], mask1), mask2);
	r[2] = vorrq_u64(vandq_u64(r[2], mask1), mask2);
	r[3] = vorrq_u64(vandq_u64(r[3], mask1), mask2);

	re[0] = vsubq_f64(re[0], one);
	re[1] = vsubq_f64(re[1], one);
	re[2] = vsubq_f64(re[2], one);
	re[3] = vsubq_f64(re[3], one);
}

/* Draw a uniform distriubted real number */
TARGETNEO void RNGNEO::Poly(float64x2_t* a, float64x2_t* s, int n, uint64x2_t* re)
{
	float64x2_t t[4];

	Uniform(t);

	float64x2_t t1 = vmulq_f64(t[0], s[0]);
	float64x2_t t2 = vmulq_f64(t[1], s[1]);
	float64x2_t t3 = vmulq_f64(t[2], s[2]);
	float64x2_t t4 = vmulq_f64(t[3], s[3]);

	uint64x2_t f1 = vdupq_n_u64(0);
	uint64x2_t f2 = vdupq_n_u64(0);
	uint64x2_t f3 = vdupq_n_u64(0);
	uint64x2_t f4 = vdupq_n_u64(0);

	uint64x2_t midx1 = vdupq_n_u64(n - 1);
	uint64x2_t midx2 = vdupq_n_u64(n - 1);
	uint64x2_t midx3 = vdupq_n_u64(n - 1);
	uint64x2_t midx4 = vdupq_n_u64(n - 1);

	uint64x2_t nidx = vdupq_n_u64(0);
	uint64x2_t ninc = vdupq_n_u64(1);

	uint64x2_t b1, b2, b3, b4;
	float64x2_t v1, v2, v3, v4;

	for (int i = 0; i < n; ++i)
	{
		v1 = a[i * 4 + 0];
		v2 = a[i * 4 + 1];
		v3 = a[i * 4 + 2];
		v4 = a[i * 4 + 3];

		b1 = vcltq_f64(t1, v1);
		b2 = vcltq_f64(t2, v2);
		b3 = vcltq_f64(t3, v3);
		b4 = vcltq_f64(t4, v4);

		t1 = vsubq_f64(t1, v1);
		t2 = vsubq_f64(t2, v2);
		t3 = vsubq_f64(t3, v3);
		t4 = vsubq_f64(t4, v4);

		b1 = vbicq_u64(b1, f1);
		b2 = vbicq_u64(b2, f2);
		b3 = vbicq_u64(b3, f3);
		b4 = vbicq_u64(b4, f4);

		f1 = vorrq_u64(f1, b1);
		f2 = vorrq_u64(f2, b2);
		f3 = vorrq_u64(f3, b3);
		f4 = vorrq_u64(f4, b4);

		midx1 = vbslq_s64(b1, nidx, midx1);
		midx2 = vbslq_s64(b2, nidx, midx2);
		midx3 = vbslq_s64(b3, nidx, midx3);
		midx4 = vbslq_s64(b4, nidx, midx4);

		nidx = vaddq_s64(nidx, ninc);
	}

	re[0] = midx1;
	re[1] = midx2;
	re[2] = midx3;
	re[3] = midx4;
}

/* Draw a polynormial distriubted integer with propoirtions in natural logarithm */
TARGETNEO void RNGNEO::PolyLog(float64x2_t* a, int n, uint64x2_t* re)
{
	//proportional polynomial distribution, will overwrite a
	float64x2_t max1 = vdupq_n_f64(-1e300);
	float64x2_t max2 = vdupq_n_f64(-1e300);
	float64x2_t max3 = vdupq_n_f64(-1e300);
	float64x2_t max4 = vdupq_n_f64(-1e300);

	float64x2_t s1 = vdupq_n_f64(MIN_FREQ * n);
	float64x2_t s2 = vdupq_n_f64(MIN_FREQ * n);
	float64x2_t s3 = vdupq_n_f64(MIN_FREQ * n);
	float64x2_t s4 = vdupq_n_f64(MIN_FREQ * n);

	double* af = (double*)a;

	for (int i = 0; i < n; ++i)
	{
		max1 = vmaxq_f64(max1, a[i * 4 + 0]);
		max2 = vmaxq_f64(max2, a[i * 4 + 1]);
		max3 = vmaxq_f64(max3, a[i * 4 + 2]);
		max4 = vmaxq_f64(max4, a[i * 4 + 3]);
	}

	for (int i = 0; i < n; ++i)
	{
		a[i * 4 + 0] = vsubq_f64(a[i * 4 + 0], max1);
		a[i * 4 + 1] = vsubq_f64(a[i * 4 + 1], max2);
		a[i * 4 + 2] = vsubq_f64(a[i * 4 + 2], max3);
		a[i * 4 + 3] = vsubq_f64(a[i * 4 + 3], max4);

		af[i * 8 + 0] = (af[i * 8 + 0] < -23) ? MIN_FREQ : exp(af[i * 8 + 0]);
		af[i * 8 + 1] = (af[i * 8 + 1] < -23) ? MIN_FREQ : exp(af[i * 8 + 1]);
		af[i * 8 + 2] = (af[i * 8 + 2] < -23) ? MIN_FREQ : exp(af[i * 8 + 2]);
		af[i * 8 + 3] = (af[i * 8 + 3] < -23) ? MIN_FREQ : exp(af[i * 8 + 3]);
		af[i * 8 + 4] = (af[i * 8 + 4] < -23) ? MIN_FREQ : exp(af[i * 8 + 4]);
		af[i * 8 + 5] = (af[i * 8 + 5] < -23) ? MIN_FREQ : exp(af[i * 8 + 5]);
		af[i * 8 + 6] = (af[i * 8 + 6] < -23) ? MIN_FREQ : exp(af[i * 8 + 6]);
		af[i * 8 + 7] = (af[i * 8 + 7] < -23) ? MIN_FREQ : exp(af[i * 8 + 7]);

		s1 = vaddq_f64(s1, a[i * 4 + 0]);
		s2 = vaddq_f64(s2, a[i * 4 + 1]);
		s3 = vaddq_f64(s3, a[i * 4 + 2]);
		s4 = vaddq_f64(s4, a[i * 4 + 3]);
	}

	float64x2_t t[4];

	Uniform(t);

	float64x2_t t1 = vmulq_f64(t[0], s1);
	float64x2_t t2 = vmulq_f64(t[1], s2);
	float64x2_t t3 = vmulq_f64(t[2], s3);
	float64x2_t t4 = vmulq_f64(t[3], s4);

	uint64x2_t f1 = vdupq_n_u64(0);
	uint64x2_t f2 = vdupq_n_u64(0);
	uint64x2_t f3 = vdupq_n_u64(0);
	uint64x2_t f4 = vdupq_n_u64(0);

	uint64x2_t midx1 = vdupq_n_u64(n - 1);
	uint64x2_t midx2 = vdupq_n_u64(n - 1);
	uint64x2_t midx3 = vdupq_n_u64(n - 1);
	uint64x2_t midx4 = vdupq_n_u64(n - 1);

	uint64x2_t nidx = vdupq_n_u64(0);
	uint64x2_t ninc = vdupq_n_u64(1);

	uint64x2_t b1, b2, b3, b4;
	float64x2_t v1, v2, v3, v4;

	for (int i = 0; i < n; ++i)
	{
		v1 = a[i * 4 + 0];
		v2 = a[i * 4 + 1];
		v3 = a[i * 4 + 2];
		v4 = a[i * 4 + 3];

		b1 = vcltq_f64(t1, v1);
		b2 = vcltq_f64(t2, v2);
		b3 = vcltq_f64(t3, v3);
		b4 = vcltq_f64(t4, v4);

		t1 = vsubq_f64(t1, v1);
		t2 = vsubq_f64(t2, v2);
		t3 = vsubq_f64(t3, v3);
		t4 = vsubq_f64(t4, v4);

		b1 = vbicq_u64(b1, f1);
		b2 = vbicq_u64(b2, f2);
		b3 = vbicq_u64(b3, f3);
		b4 = vbicq_u64(b4, f4);

		f1 = vorrq_u64(f1, b1);
		f2 = vorrq_u64(f2, b2);
		f3 = vorrq_u64(f3, b3);
		f4 = vorrq_u64(f4, b4);

		midx1 = vbslq_s64(b1, nidx, midx1);
		midx2 = vbslq_s64(b2, nidx, midx2);
		midx3 = vbslq_s64(b3, nidx, midx3);
		midx4 = vbslq_s64(b4, nidx, midx4);

		nidx = vaddq_s64(nidx, ninc);
	}

	re[0] = midx1;
	re[1] = midx2;
	re[2] = midx3;
	re[3] = midx4;
}
#endif

TARGETNEO int64 GetMinIdxNEO(double* A, int64 n, double& val)
{
	int64 i = 0;
	val = 1e300;
	uint64 idx = (uint64)-1;

	if (n >= 8)
	{
		float64x2_t min1 = vdupq_n_f64(val);
		float64x2_t min2 = vdupq_n_f64(val);
		float64x2_t min3 = vdupq_n_f64(val);
		float64x2_t min4 = vdupq_n_f64(val);

		uint64x2_t midx1 = vdupq_n_u64(0xFFFFFFFFFFFFFFFF);
		uint64x2_t midx2 = vdupq_n_u64(0xFFFFFFFFFFFFFFFF);
		uint64x2_t midx3 = vdupq_n_u64(0xFFFFFFFFFFFFFFFF);
		uint64x2_t midx4 = vdupq_n_u64(0xFFFFFFFFFFFFFFFF);

		uint64 ni[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };

		uint64x2_t nidx1 = vld1q_u64(&ni[0]);
		uint64x2_t nidx2 = vld1q_u64(&ni[2]);
		uint64x2_t nidx3 = vld1q_u64(&ni[4]);
		uint64x2_t nidx4 = vld1q_u64(&ni[6]);

		uint64x2_t msep = vdupq_n_u64(8);

		float64x2_t f1, f2, f3, f4;
		float64x2_t v1, v2, v3, v4;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			v1 = vld1q_f64(A); A += 2;
			v2 = vld1q_f64(A); A += 2;
			v3 = vld1q_f64(A); A += 2;
			v4 = vld1q_f64(A); A += 2;

			f1 = vcgtq_f64(min1, v1);
			f2 = vcgtq_f64(min2, v2);
			f3 = vcgtq_f64(min3, v3);
			f4 = vcgtq_f64(min4, v4);

			min1 = vbslq_s64(f1, v1, min1);
			min2 = vbslq_s64(f2, v2, min2);
			min3 = vbslq_s64(f3, v3, min3);
			min4 = vbslq_s64(f4, v4, min4);

			midx1 = vbslq_s64(f1, nidx1, midx1);
			midx2 = vbslq_s64(f2, nidx2, midx2);
			midx3 = vbslq_s64(f3, nidx3, midx3);
			midx4 = vbslq_s64(f4, nidx4, midx4);

			nidx1 = vaddq_s64(nidx1, msep);
			nidx2 = vaddq_s64(nidx2, msep);
			nidx3 = vaddq_s64(nidx3, msep);
			nidx4 = vaddq_s64(nidx4, msep);
		}

		f1 = vcgtq_f64(min1, min2);
		f3 = vcgtq_f64(min3, min4);

		min1 = vbslq_s64(f1, min2, min1);
		min3 = vbslq_s64(f3, min4, min3);

		midx1 = vbslq_s64(f1, midx2, midx1);
		midx3 = vbslq_s64(f3, midx4, midx3);

		f1 = vcgtq_f64(min1, min3);
		min1 = vbslq_s64(f1, min3, min1);
		midx1 = vbslq_s64(f1, midx3, midx1);

		if (vgetq_lane_f64(min1, 0) < val)
		{
			val = vgetq_lane_f64(min1, 0);
			idx = vgetq_lane_u64(midx1, 0);
		}

		if (vgetq_lane_f64(min1, 1) < val)
		{
			val = vgetq_lane_f64(min1, 1);
			idx = vgetq_lane_u64(midx1, 1);
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
	int64 i = 0;
	minv = 1e300;
	maxv = -1e300;

	if (n >= 8)
	{
		float64x2_t max1 = vdupq_n_f64(maxv);
		float64x2_t max2 = vdupq_n_f64(maxv);
		float64x2_t max3 = vdupq_n_f64(maxv);
		float64x2_t max4 = vdupq_n_f64(maxv);
		float64x2_t min1 = vdupq_n_f64(minv);
		float64x2_t min2 = vdupq_n_f64(minv);
		float64x2_t min3 = vdupq_n_f64(minv);
		float64x2_t min4 = vdupq_n_f64(minv);
		float64x2_t v1, v2, v3, v4;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			v1 = vld1q_f64(A); A += 2;
			v2 = vld1q_f64(A); A += 2;
			v3 = vld1q_f64(A); A += 2;
			v4 = vld1q_f64(A); A += 2;

			min1 = vminq_f64(min1, v1);
			min1 = vminq_f64(min1, v1);
			min2 = vminq_f64(min2, v2);
			min3 = vminq_f64(min3, v3);
			min4 = vminq_f64(min4, v4);

			max1 = vmaxq_f64(max1, v1);
			max2 = vmaxq_f64(max2, v2);
			max3 = vmaxq_f64(max3, v3);
			max4 = vmaxq_f64(max4, v4);
		}

		min1 = vminq_f64(vminq_f64(min1, min2), vminq_f64(min3, min4));
		max1 = vmaxq_f64(vmaxq_f64(max1, max2), vmaxq_f64(max3, max4));

		if (vgetq_lane_f64(min1, 0) < minv) minv = vgetq_lane_f64(min1, 0);
		if (vgetq_lane_f64(max1, 0) > maxv) maxv = vgetq_lane_f64(max1, 0);
		if (vgetq_lane_f64(min1, 1) < minv) minv = vgetq_lane_f64(min1, 1);
		if (vgetq_lane_f64(max1, 1) > maxv) maxv = vgetq_lane_f64(max1, 1);
	}

	for (; i < n; ++i, ++A)
	{
		if (*A < minv) minv = *A;
		if (*A > maxv) maxv = *A;
	}
}

TARGETNEO double GetMaxValNEO(double* A, int64 n)
{
	int64 i = 0;
	double val = -1e300;

	if (n >= 8)
	{
		float64x2_t max1 = vdupq_n_f64(val);
		float64x2_t max2 = vdupq_n_f64(val);
		float64x2_t max3 = vdupq_n_f64(val);
		float64x2_t max4 = vdupq_n_f64(val);
		float64x2_t v1, v2, v3, v4;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			v1 = vld1q_f64(A); A += 2;
			v2 = vld1q_f64(A); A += 2;
			v3 = vld1q_f64(A); A += 2;
			v4 = vld1q_f64(A); A += 2;

			max1 = vmaxq_f64(max1, v1);
			max2 = vmaxq_f64(max2, v2);
			max3 = vmaxq_f64(max3, v3);
			max4 = vmaxq_f64(max4, v4);
		}

		max1 = vmaxq_f64(vmaxq_f64(max1, max2), vmaxq_f64(max3, max4));

		if (vgetq_lane_f64(max1, 0) > val) val = vgetq_lane_f64(max1, 0);
		if (vgetq_lane_f64(max1, 1) > val) val = vgetq_lane_f64(max1, 1);
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
	double val = -1e300;

	if (n >= 8)
	{
		float64x2_t max1 = vdupq_n_f64(val);
		float64x2_t max2 = vdupq_n_f64(val); 
		float64x2_t max3 = vdupq_n_f64(val);
		float64x2_t max4 = vdupq_n_f64(val);
		float64x2_t v1, v2, v3, v4;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			double t[8] = {
				A[0 * sep],  A[1 * sep], A[2 * sep],  A[3 * sep],
				A[1 * sep],  A[5 * sep], A[6 * sep],  A[7 * sep],
			};

			v1 = vld1q_f64(&t[0]);
			v2 = vld1q_f64(&t[2]);
			v3 = vld1q_f64(&t[4]);
			v4 = vld1q_f64(&t[6]);

			A += sep * 8;
			max1 = vmaxq_f64(max1, v1);
			max2 = vmaxq_f64(max2, v2);
			max3 = vmaxq_f64(max3, v3);
			max4 = vmaxq_f64(max4, v4);
		}

		max1 = vmaxq_f64(vmaxq_f64(max1, max2), vmaxq_f64(max3, max4));

		if (vgetq_lane_f64(max1, 0) > val) val = vgetq_lane_f64(max1, 0);
		if (vgetq_lane_f64(max1, 1) > val) val = vgetq_lane_f64(max1, 1);
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
	int64 i = 0;
	double val = 1e300;

	if (n >= 8)
	{
		float64x2_t min1 = vdupq_n_f64(val);
		float64x2_t min2 = vdupq_n_f64(val);
		float64x2_t min3 = vdupq_n_f64(val);
		float64x2_t min4 = vdupq_n_f64(val);
		float64x2_t v1, v2, v3, v4;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			v1 = vld1q_f64(A); A += 2;
			v2 = vld1q_f64(A); A += 2;
			v3 = vld1q_f64(A); A += 2;
			v4 = vld1q_f64(A); A += 2;

			min1 = vminq_f64(min1, v1);
			min2 = vminq_f64(min2, v2);
			min3 = vminq_f64(min3, v3);
			min4 = vminq_f64(min4, v4);
		}

		min1 = vminq_f64(vminq_f64(min1, min2), vminq_f64(min3, min4));

		if (vgetq_lane_f64(min1, 0) < val) val = vgetq_lane_f64(min1, 0);
		if (vgetq_lane_f64(min1, 1) < val) val = vgetq_lane_f64(min1, 1);

		for (; i < n; ++i, ++A)
		{
			if (*A > val) continue;
			val = *A;
		}
	}

	return val;
}

TARGETNEO int64 GetMinValNEO(int64* A, int64 n)
{
	int64 i = 0;
	int64 val = 0x7FFFFFFFFFFFFFFF;

	if (n >= 8)
	{
		uint64x2_t min1 = vdupq_n_u64(0x7FFFFFFFFFFFFFFF);
		uint64x2_t min2 = vdupq_n_u64(0x7FFFFFFFFFFFFFFF);
		uint64x2_t min3 = vdupq_n_u64(0x7FFFFFFFFFFFFFFF);
		uint64x2_t min4 = vdupq_n_u64(0x7FFFFFFFFFFFFFFF);
		uint64x2_t v1, v2, v3, v4;
		uint64x2_t f1, f2, f3, f4;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			v1 = vld1q_s64(A); A += 2;
			v2 = vld1q_s64(A); A += 2;
			v3 = vld1q_s64(A); A += 2;
			v4 = vld1q_s64(A); A += 2;

			f1 = vcgtq_s64(min1, v1);
			f2 = vcgtq_s64(min2, v2);
			f3 = vcgtq_s64(min3, v3);
			f4 = vcgtq_s64(min4, v4);

			min1 = vbslq_s64(f1, v1, min1);
			min2 = vbslq_s64(f2, v2, min2);
			min3 = vbslq_s64(f3, v3, min3);
			min4 = vbslq_s64(f4, v4, min4);
		}

		f1 = vcgtq_s64(min1, min2);
		f3 = vcgtq_s64(min3, min4);

		min1 = vbslq_s64(f1, min2, min1);
		min3 = vbslq_s64(f3, min4, min3);

		f1 = vcgtq_s64(min1, min3);
		min1 = vbslq_s64(f1, min3, min1);

		if (vgetq_lane_s64(min1, 0) < val) val = vgetq_lane_s64(min1, 0);
		if (vgetq_lane_s64(min1, 1) < val) val = vgetq_lane_s64(min1, 1);
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

TARGETNEO void ChargeLogNEO(int64& slog, double& prod, float64x2_t& val)
{
	AddExponentNEO(slog, val);
	prod = prod * vgetq_lane_f64(val, 0) * vgetq_lane_f64(val, 1);

	if (prod < DOUBLE_UNDERFLOW || prod > DOUBLE_OVERFLOW) [[unlikely]]
		AddExponent(slog, prod);
}

TARGETNEO double LogProdNEO(double* A, int64 n)
{
	int64 i = 0;
	int64 slog = 0; double prod = 1;

	if (n >= 2)
	{
		float64x2_t pd = vdupq_n_f64(1.0), dunder = vdupq_n_f64(DOUBLE_UNDERFLOW), dover = vdupq_n_f64(DOUBLE_OVERFLOW);
		uint64x2_t flag;

		for (int64 l1 = n - 2; i <= l1; i += 2)
		{
			pd = vmulq_f64(pd, vld1q_f64(A)); A += 2;

			flag = vorrq_u64(vcltq_f64(pd, dunder), vcltq_f64(dover, pd));

			if (vgetq_lane_u64(flag, 0) | vgetq_lane_u64(flag, 1)) [[unlikely]]
				AddExponentNEO(slog, pd);
		}

		ChargeLogNEO(slog, prod, pd);
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

	if (n >= 2)
	{
		__builtin_prefetch(A, 0); A += sep;
		__builtin_prefetch(A, 0); A += sep;

		float64x2_t pd = vdupq_n_f64(1.0), dunder = vdupq_n_f64(DOUBLE_UNDERFLOW), dover = vdupq_n_f64(DOUBLE_OVERFLOW);
		uint64x2_t flag;

		for (int64 l1 = n - 2; i <= l1; i += 2)
		{
			__builtin_prefetch(A, 0); A += sep;
			__builtin_prefetch(A, 0); A += sep;

			double a[2] = { A[-4 * sep], A[-3 * sep] };

			pd = vmulq_f64(pd, vld1q_f64(a));

			flag = vorrq_u64(vcltq_f64(pd, dunder), vcltq_f64(dover, pd));

			if (vgetq_lane_u64(flag, 0) | vgetq_lane_u64(flag, 1)) [[unlikely]]
				AddExponentNEO(slog, pd);
		}

		ChargeLogNEO(slog, prod, pd);
		A -= sep * 2;
	}

	for (; i < n; ++i, A += sep)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);
	return prod;
}

TARGETNEO double LogProdDivNEO(double* A, double* B, int64 n, int64 sep)
{
	int64 i = 0;
	int64 slog = 0; double prod = 1;

	if (n >= 2)
	{
		__builtin_prefetch(A, 0); A += sep;
		__builtin_prefetch(A, 0); A += sep;
		__builtin_prefetch(B, 0); B += sep;
		__builtin_prefetch(B, 0); B += sep;

		float64x2_t pd = vdupq_n_f64(1.0), dunder = vdupq_n_f64(DOUBLE_UNDERFLOW), dover = vdupq_n_f64(DOUBLE_OVERFLOW);
		uint64x2_t flag;

		for (int64 l1 = n - 2; i <= l1; i += 2)
		{
			__builtin_prefetch(A, 0); A += sep;
			__builtin_prefetch(A, 0); A += sep;
			__builtin_prefetch(B, 0); B += sep;
			__builtin_prefetch(B, 0); B += sep;
			
			double a[2] = { A[-4 * sep], A[-3 * sep] };
			double b[2] = { B[-4 * sep], B[-3 * sep] };

			pd = vmulq_f64(pd, vdivq_f64(vld1q_f64(a), vld1q_f64(b)));

			flag = vorrq_u64(vcltq_f64(pd, dunder), vcltq_f64(dover, pd));

			if (vgetq_lane_u64(flag, 0) | vgetq_lane_u64(flag, 1)) [[unlikely]]
				AddExponentNEO(slog, pd);
		}

		ChargeLogNEO(slog, prod, pd);
		A -= sep * 2; 
		B -= sep * 2;
	}

	for (; i < n; ++i, A += sep, B += sep)
		ChargeLog(slog, prod, *A / *B);

	CloseLog(slog, prod);
	return prod;
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
	int64 i = 0;
	double re = 0;

	if (n >= 8)
	{
		float64x2_t s1 = (float64x2_t)vdupq_n_u64(0);
		float64x2_t s2 = (float64x2_t)vdupq_n_u64(0);
		float64x2_t s3 = (float64x2_t)vdupq_n_u64(0);
		float64x2_t s4 = (float64x2_t)vdupq_n_u64(0);
		float64x2_t v1, v2, v3, v4;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			v1 = vld1q_f64(A); A += 2;
			v2 = vld1q_f64(A); A += 2;
			v3 = vld1q_f64(A); A += 2;
			v4 = vld1q_f64(A); A += 2;

			s1 = vaddq_f64(s1, v1);
			s2 = vaddq_f64(s2, v2);
			s3 = vaddq_f64(s3, v3);
			s4 = vaddq_f64(s4, v4);
		}

		s1 = vaddq_f64(vaddq_f64(s1, s2), vaddq_f64(s3, s4));
		re = vgetq_lane_f64(s1, 0) + vgetq_lane_f64(s1, 1);
	}

	for (; i < n; ++i)
		re += *A++;

	return re;
}

TARGETNEO int64 SumNEO(byte* A, int64 n)
{
	uint64 re = 0;
	int64 i = 0;

	if (n >= 64)
	{
		uint8x16_t v1, v2, v3, v4;
		uint64x2_t s1 = vdupq_n_u64(0);

		for (int64 l1 = n - 64; i <= l1; i += 64)
		{
			v1 = vld1q_u8(A); A += 16;
			v2 = vld1q_u8(A); A += 16;
			v3 = vld1q_u8(A); A += 16;
			v4 = vld1q_u8(A); A += 16;

			v1 = vaddq_u8(v1, v2);
			v3 = vaddq_u8(v3, v4);
			v1 = vaddq_u8(v1, v3);

			s1 = vaddq_u64(s1, vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(v1))));
		}

		re += vgetq_lane_u64(s1, 0) + vgetq_lane_u64(s1, 1);
	}

	for (; i < n; ++i)
		re += *A++;

	return re;
}

TARGETNEO double SumNEO(double* A, int64 n, int64 sep)
{
	int64 i = 0;
	double re = 0;

	if (n >= 4)
	{
		__builtin_prefetch(A, 0); A += sep;
		__builtin_prefetch(A, 0); A += sep;
		__builtin_prefetch(A, 0); A += sep;
		__builtin_prefetch(A, 0); A += sep;

		float64x2_t s1 = (float64x2_t)vdupq_n_u64(0), s2 = (float64x2_t)vdupq_n_u64(0);

		for (int64 l1 = n - 4; i <= l1; i += 4)
		{
			__builtin_prefetch(A, 0); A += sep;
			__builtin_prefetch(A, 0); A += sep;
			__builtin_prefetch(A, 0); A += sep;
			__builtin_prefetch(A, 0); A += sep;

			double t[8] = { A[-8 * sep], A[-7 * sep], A[-6 * sep], A[-5 * sep] };

			s1 = vaddq_f64(s1, vld1q_f64(&t[0]));

			s2 = vaddq_f64(s2, vld1q_f64(&t[2]));
		}

		s1 = vaddq_f64(s1, s2);
		//s1 = _mm_hadd_pd(s1, s1);
		re = vgetq_lane_f64(s1, 0) + vgetq_lane_f64(s1, 1);
		A -= sep * 4;
	}

	for (; i < n; ++i, A += sep)
		re += *A;

	return re;
}

TARGETNEO void SumNEO(double* A, double** B, int64 k, int64 n)
{
	int64 i = 0;

	if (n >= 8)
	{
		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			float64x2_t a1 = (float64x2_t)vdupq_n_u64(0);
			float64x2_t a2 = (float64x2_t)vdupq_n_u64(0);
			float64x2_t a3 = (float64x2_t)vdupq_n_u64(0);
			float64x2_t a4 = (float64x2_t)vdupq_n_u64(0);

			for (int64 j = 0; j < k; ++j)
			{
				a1 = vaddq_f64(a1, vld1q_f64(&B[j][i + 0]));
				a2 = vaddq_f64(a2, vld1q_f64(&B[j][i + 2]));
				a3 = vaddq_f64(a3, vld1q_f64(&B[j][i + 4]));
				a4 = vaddq_f64(a4, vld1q_f64(&B[j][i + 6]));
			}

			vst1q_f64(&A[i + 0], a1);
			vst1q_f64(&A[i + 2], a2);
			vst1q_f64(&A[i + 4], a3);
			vst1q_f64(&A[i + 6], a4);
		}
	}

	for (; i < n; ++i)
	{
		A[i] = 0;
		for (int64 j = 0; j < k; ++j)
			A[i] += B[j][i];
	}
}

TARGETNEO double ProdNEO(double* A, int64 n)
{
	int64 i = 0;
	double re = 1;

	if (n >= 16)
	{
		float64x2_t pd1 = vdupq_n_f64(1.0);
		float64x2_t pd2 = vdupq_n_f64(1.0);
		float64x2_t pd3 = vdupq_n_f64(1.0);
		float64x2_t pd4 = vdupq_n_f64(1.0);
		float64x2_t pd5 = vdupq_n_f64(1.0);
		float64x2_t pd6 = vdupq_n_f64(1.0);
		float64x2_t pd7 = vdupq_n_f64(1.0);
		float64x2_t pd8 = vdupq_n_f64(1.0);

		for (int64 l1 = n - 16; i <= l1; i += 16)
		{
			pd1 = vmulq_f64(pd1, vld1q_f64(A));
			A += 2;

			pd2 = vmulq_f64(pd2, vld1q_f64(A));
			A += 2;

			pd3 = vmulq_f64(pd3, vld1q_f64(A));
			A += 2;

			pd4 = vmulq_f64(pd4, vld1q_f64(A));
			A += 2;

			pd5 = vmulq_f64(pd5, vld1q_f64(A));
			A += 2;

			pd6 = vmulq_f64(pd6, vld1q_f64(A));
			A += 2;

			pd7 = vmulq_f64(pd7, vld1q_f64(A));
			A += 2;

			pd8 = vmulq_f64(pd8, vld1q_f64(A));
			A += 2;
		}

		pd1 = vmulq_f64(pd1, pd2);
		pd3 = vmulq_f64(pd3, pd4);
		pd5 = vmulq_f64(pd5, pd6);
		pd7 = vmulq_f64(pd7, pd8);

		pd1 = vmulq_f64(pd1, pd3);
		pd5 = vmulq_f64(pd5, pd7);

		pd1 = vmulq_f64(pd1, pd5);

		re = vgetq_lane_f64(pd1, 0) * vgetq_lane_f64(pd1, 1);
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
		__builtin_prefetch(A, 0); A += sep;
		__builtin_prefetch(A, 0); A += sep;
		__builtin_prefetch(A, 0); A += sep;
		__builtin_prefetch(A, 0); A += sep;
		__builtin_prefetch(A, 0); A += sep;
		__builtin_prefetch(A, 0); A += sep;
		__builtin_prefetch(A, 0); A += sep;
		__builtin_prefetch(A, 0); A += sep;

		float64x2_t pd1 = vdupq_n_f64(1.0), pd2 = vdupq_n_f64(1.0);
		float64x2_t pd3 = vdupq_n_f64(1.0), pd4 = vdupq_n_f64(1.0);

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			__builtin_prefetch(A, 0); A += sep;
			__builtin_prefetch(A, 0); A += sep;
			__builtin_prefetch(A, 0); A += sep;
			__builtin_prefetch(A, 0); A += sep;
			__builtin_prefetch(A, 0); A += sep;
			__builtin_prefetch(A, 0); A += sep;
			__builtin_prefetch(A, 0); A += sep;
			__builtin_prefetch(A, 0); A += sep;

			double t[8] = {  
				A[-16 * sep], A[-15 * sep], A[-14 * sep], A[-13 * sep],
				A[-12 * sep], A[-11 * sep], A[-10 * sep], A[ -9 * sep]
			};

			pd1 = vmulq_f64(pd1, vld1q_f64(&t[0]));
			pd2 = vmulq_f64(pd2, vld1q_f64(&t[2]));
			pd3 = vmulq_f64(pd3, vld1q_f64(&t[4]));
			pd4 = vmulq_f64(pd4, vld1q_f64(&t[6]));
		}

		pd1 = vmulq_f64(vmulq_f64(pd1, pd2), vmulq_f64(pd3, pd4));
		re = vgetq_lane_f64(pd1, 0) * vgetq_lane_f64(pd1, 1);
		A -= sep * 8;
	}

	for (; i < n; ++i, A += sep)
		re *= *A;

	return re;
}

TARGETNEO double SumSquareNEO(double* A, int64 n)
{
	int64 i = 0;
	double re = 0;

	if (n >= 8)
	{
		float64x2_t s1 = (float64x2_t)vdupq_n_u64(0);
		float64x2_t s2 = (float64x2_t)vdupq_n_u64(0);
		float64x2_t s3 = (float64x2_t)vdupq_n_u64(0);
		float64x2_t s4 = (float64x2_t)vdupq_n_u64(0);
		float64x2_t a1, a2, a3, a4;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			a1 = vld1q_f64(A); A += 2;
			a2 = vld1q_f64(A); A += 2;
			a3 = vld1q_f64(A); A += 2;
			a4 = vld1q_f64(A); A += 2;

			s1 = vaddq_f64(s1, vmulq_f64(a1, a1));
			s2 = vaddq_f64(s2, vmulq_f64(a2, a2));
			s3 = vaddq_f64(s3, vmulq_f64(a3, a3));
			s4 = vaddq_f64(s4, vmulq_f64(a4, a4));
		}

		s1 = vaddq_f64(vaddq_f64(s1, s2), vaddq_f64(s3, s4));
		//s1 = _mm_hadd_pd(s1, s1);
		re = vgetq_lane_f64(s1, 0) + vgetq_lane_f64(s1, 1);
	}

	for (; i < n; ++i, ++A)
		re += *A * *A;

	return re;
}

TARGETNEO int64 SumSquareNEO(byte* A, int64 n)
{
	int64 i = 0;
	uint64 re = 0;

	if (n >= 64)
	{
		uint64x2_t s1 = vdupq_n_u64(0), s2 = vdupq_n_u64(0);
		uint8x16_t v1, v2, v3, v4;

		for (int64 l1 = n - 64; i <= l1; i += 64)
		{
			v1 = vld1q_u8(A); A += 16;
			v2 = vld1q_u8(A); A += 16;
			v3 = vld1q_u8(A); A += 16;
			v4 = vld1q_u8(A); A += 16;

			v1 = vmulq_u8(v1, v1);
			v2 = vmulq_u8(v2, v2);
			v3 = vmulq_u8(v3, v3);
			v4 = vmulq_u8(v4, v4);

			v1 = vaddq_u8(v1, v2);
			v3 = vaddq_u8(v3, v4);

			s1 = vaddq_u64(s1, vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(v1))));
			s2 = vaddq_u64(s2, vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(v3))));
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
	int64 i = 0;
	sum = sumsq = 0;

	if (n >= 8)
	{
		float64x2_t s1 = (float64x2_t)vdupq_n_u64(0), sq1 = (float64x2_t)vdupq_n_u64(0);
		float64x2_t s2 = (float64x2_t)vdupq_n_u64(0), sq2 = (float64x2_t)vdupq_n_u64(0);
		float64x2_t a;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			a = vld1q_f64(A); A += 2;
			s1 = vaddq_f64(s1, a);
			sq1 = vaddq_f64(sq1, vmulq_f64(a, a));

			a = vld1q_f64(A); A += 2;
			s2 = vaddq_f64(s2, a);
			sq2 = vaddq_f64(sq2, vmulq_f64(a, a));

			a = vld1q_f64(A); A += 2;
			s1 = vaddq_f64(s1, a);
			sq1 = vaddq_f64(sq1, vmulq_f64(a, a));

			a = vld1q_f64(A); A += 2;
			s2 = vaddq_f64(s2, a);
			sq2 = vaddq_f64(sq2, vmulq_f64(a, a));
		}

		s1 = vaddq_f64(s1, s2);
		sq1 = vaddq_f64(sq1, sq2);

		//s1 = _mm_hadd_pd(s1, s1);
		sum = vgetq_lane_f64(s1, 0) + vgetq_lane_f64(s1, 1);
		//sq1 = _mm_hadd_pd(sq1, sq1);
		sumsq = vgetq_lane_f64(sq1, 0) + vgetq_lane_f64(sq1, 1);
	}

	for (; i < n; ++i, ++A)
	{
		sum += *A;
		sumsq += *A * *A;
	}
}

/* re = Sum(A1[i++] * B[j += sep]) / Sum(A2[i++] * B[j += sep]) */
TARGETNEO double SumProdDivNEO(double* A1, double* A2, double* B, int64 sep, int64 n)
{
	int64 i = 0;
	double re1 = 0, re2 = 0;

	if (n >= 4)
	{
		__builtin_prefetch(B, 0); B += sep;
		__builtin_prefetch(B, 0); B += sep;
		__builtin_prefetch(B, 0); B += sep;
		__builtin_prefetch(B, 0); B += sep;
		__builtin_prefetch(B, 0); B += sep;
		__builtin_prefetch(B, 0); B += sep;
		__builtin_prefetch(B, 0); B += sep;
		__builtin_prefetch(B, 0); B += sep;

		float64x2_t s1 = (float64x2_t)vdupq_n_u64(0), s2 = (float64x2_t)vdupq_n_u64(0), b;

		for (int64 l1 = n - 4; i <= l1; i += 4)
		{
			__builtin_prefetch(B, 0); B += sep;
			__builtin_prefetch(B, 0); B += sep;

			double t1[2] = { B[-10 * sep], B[-9 * sep] };
			b = vld1q_f64(t1);

			s1 = vaddq_f64(s1, vmulq_f64(vld1q_f64(A1), b)); A1 += 2;
			s2 = vaddq_f64(s2, vmulq_f64(vld1q_f64(A2), b)); A2 += 2;

			__builtin_prefetch(B, 0); B += sep;
			__builtin_prefetch(B, 0); B += sep;

			double t2[2] = { B[-10 * sep], B[-9 * sep] };
			b = vld1q_f64(t2);

			s1 = vaddq_f64(s1, vmulq_f64(vld1q_f64(A1), b)); A1 += 2;
			s2 = vaddq_f64(s2, vmulq_f64(vld1q_f64(A2), b)); A2 += 2;
		}

		//s1 = _mm_hadd_pd(s1, s1);
		//s2 = _mm_hadd_pd(s2, s2);
		re1 = vgetq_lane_f64(s1, 0) + vgetq_lane_f64(s1, 1);
		re2 = vgetq_lane_f64(s2, 0) + vgetq_lane_f64(s2, 1);
		B -= sep * 8;
	}

	for (; i < n; ++i, A1++, A2++, B += sep)
	{
		re1 += *A1 * *B;
		re2 += *A2 * *B;
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

		float64x2_t s = (float64x2_t)vdupq_n_u64(0);

		for (int64 l1 = n - 2; i <= l1; i += 2)
		{
			__builtin_prefetch(B, 0); B += sep;
			__builtin_prefetch(B, 0); B += sep;
			
			double b[2] = { B[-4 * sep], B[-3 * sep] };
			s = vaddq_f64(s, vmulq_f64(vld1q_f64(A), vld1q_f64(b)));
			A += 2;
		}

		//s = _mm_hadd_pd(s, s);
		re = vgetq_lane_f64(s, 0) + vgetq_lane_f64(s, 1);
		B -= sep * 2;
	}

	for (; i < n; ++i, A++, B += sep)
		re += *A * *B;

	return re;
}

TARGETNEO double SumProdNEO(double* A, double* B, int64 n)
{
	int64 i = 0;
	double re = 0;

	if (n >= 8)
	{
		float64x2_t s1 = (float64x2_t)vdupq_n_u64(0);
		float64x2_t s2 = (float64x2_t)vdupq_n_u64(0);
		float64x2_t s3 = (float64x2_t)vdupq_n_u64(0);
		float64x2_t s4 = (float64x2_t)vdupq_n_u64(0);
		float64x2_t a1, a2, a3, a4;
		float64x2_t b1, b2, b3, b4;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			a1 = vld1q_f64(A); A += 2;
			a2 = vld1q_f64(A); A += 2;
			a3 = vld1q_f64(A); A += 2;
			a4 = vld1q_f64(A); A += 2;

			b1 = vld1q_f64(B); B += 2;
			b2 = vld1q_f64(B); B += 2;
			b3 = vld1q_f64(B); B += 2;
			b4 = vld1q_f64(B); B += 2;

			s1 = vaddq_f64(s1, vmulq_f64(a1, b1));
			s2 = vaddq_f64(s2, vmulq_f64(a2, b2));
			s3 = vaddq_f64(s3, vmulq_f64(a3, b3));
			s4 = vaddq_f64(s4, vmulq_f64(a4, b4));
		}

		s1 = vaddq_f64(vaddq_f64(s1, s2), vaddq_f64(s3, s4));
		//s1 = _mm_hadd_pd(s1, s1);
		re = vgetq_lane_f64(s1, 0) + vgetq_lane_f64(s1, 1);
	}

	for (; i < n; ++i, ++A, ++B)
		re += *A * *B;

	return re;
}

TARGETNEO void AddNEO(double* A, double* B, int64 n)
{
	int64 i = 0;

	if (n >= 8)
	{
		float64x2_t a1, a2, a3, a4;
		float64x2_t b1, b2, b3, b4;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			a1 = vld1q_f64(A); A += 2;
			a2 = vld1q_f64(A); A += 2;
			a3 = vld1q_f64(A); A += 2;
			a4 = vld1q_f64(A); A += 2;

			b1 = vld1q_f64(B); B += 2;
			b2 = vld1q_f64(B); B += 2;
			b3 = vld1q_f64(B); B += 2;
			b4 = vld1q_f64(B); B += 2;

			a1 = vaddq_f64(a1, b1);
			a2 = vaddq_f64(a2, b2);
			a3 = vaddq_f64(a3, b3);
			a4 = vaddq_f64(a4, b4);

			vst1q_f64(A - 8, a1);
			vst1q_f64(A - 6, a2);
			vst1q_f64(A - 4, a3);
			vst1q_f64(A - 2, a4);
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGETNEO void AddNEO(int64* A, int64* B, int64 n)
{
	int64 i = 0;

	if (n >= 8)
	{
		int64x2_t a1, a2, a3, a4;
		int64x2_t b1, b2, b3, b4;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			a1 = vld1q_s64(A); A += 2;
			a2 = vld1q_s64(A); A += 2;
			a3 = vld1q_s64(A); A += 2;
			a4 = vld1q_s64(A); A += 2;

			b1 = vld1q_s64(B); B += 2;
			b2 = vld1q_s64(B); B += 2;
			b3 = vld1q_s64(B); B += 2;
			b4 = vld1q_s64(B); B += 2;

			a1 = vaddq_s64(a1, b1);
			a2 = vaddq_s64(a2, b2);
			a3 = vaddq_s64(a3, b3);
			a4 = vaddq_s64(a4, b4);

			vst1q_s64(&A[-8], a1);
			vst1q_s64(&A[-6], a2);
			vst1q_s64(&A[-4], a3);
			vst1q_s64(&A[-2], a4);
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGETNEO void AddNEO(int* A, int* B, int64 n)
{
	int64 i = 0;

	if (n >= 16)
	{
		int32x4_t a1, a2, a3, a4;
		int32x4_t b1, b2, b3, b4;

		for (int64 l1 = n - 16; i <= l1; i += 16)
		{
			a1 = vld1q_s32(A); A += 4;
			a2 = vld1q_s32(A); A += 4;
			a3 = vld1q_s32(A); A += 4;
			a4 = vld1q_s32(A); A += 4;

			b1 = vld1q_s32(B); B += 4;
			b2 = vld1q_s32(B); B += 4;
			b3 = vld1q_s32(B); B += 4;
			b4 = vld1q_s32(B); B += 4;

			a1 = vaddq_s32(a1, b1);
			a2 = vaddq_s32(a2, b2);
			a3 = vaddq_s32(a3, b3);
			a4 = vaddq_s32(a4, b4);

			vst1q_s32(&A[-16], a1);
			vst1q_s32(&A[-12], a2);
			vst1q_s32(&A[ -8], a3);
			vst1q_s32(&A[ -4], a4);
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGETNEO void AddNEO(double* A, double B, int64 n)
{
	int64 i = 0;

	if (n >= 8)
	{
		float64x2_t b = vdupq_n_f64(B);
		float64x2_t a1, a2, a3, a4;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			a1 = vld1q_f64(A); A += 2;
			a2 = vld1q_f64(A); A += 2;
			a3 = vld1q_f64(A); A += 2;
			a4 = vld1q_f64(A); A += 2;

			a1 = vaddq_f64(a1, b);
			a2 = vaddq_f64(a2, b);
			a3 = vaddq_f64(a3, b);
			a4 = vaddq_f64(a4, b);

			vst1q_f64(&A[-8], a1);
			vst1q_f64(&A[-6], a2);
			vst1q_f64(&A[-4], a3);
			vst1q_f64(&A[-2], a4);
		}
	}

	for (; i < n; ++i, A++)
		*A += B;
}

TARGETNEO void MulNEO(double* C, double* A, double* B, int64 n)
{
	int64 i = 0;

	if (n >= 8)
	{
		float64x2_t a1, a2, a3, a4;
		float64x2_t b1, b2, b3, b4;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			a1 = vld1q_f64(A); A += 2;
			a2 = vld1q_f64(A); A += 2;
			a3 = vld1q_f64(A); A += 2;
			a4 = vld1q_f64(A); A += 2;

			b1 = vld1q_f64(B); B += 2;
			b2 = vld1q_f64(B); B += 2;
			b3 = vld1q_f64(B); B += 2;
			b4 = vld1q_f64(B); B += 2;

			a1 = vmulq_f64(a1, b1);
			a2 = vmulq_f64(a2, b2);
			a3 = vmulq_f64(a3, b3);
			a4 = vmulq_f64(a4, b4);

			vst1q_f64(C, a1); C += 2;
			vst1q_f64(C, a2); C += 2;
			vst1q_f64(C, a3); C += 2;
			vst1q_f64(C, a4); C += 2;
		}
	}

	for (; i < n; ++i)
		*C++ = *A++ * *B++;
}

TARGETNEO void MulNEO(double* C, double* A, double B, int64 n)
{
	int64 i = 0;

	if (n >= 8)
	{
		float64x2_t b = vdupq_n_f64(B);
		float64x2_t a1, a2, a3, a4;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			a1 = vld1q_f64(A); A += 2;
			a2 = vld1q_f64(A); A += 2;
			a3 = vld1q_f64(A); A += 2;
			a4 = vld1q_f64(A); A += 2;

			a1 = vmulq_f64(a1, b);
			a2 = vmulq_f64(a2, b);
			a3 = vmulq_f64(a3, b);
			a4 = vmulq_f64(a4, b);

			vst1q_f64(C, a1); C += 2;
			vst1q_f64(C, a2); C += 2;
			vst1q_f64(C, a3); C += 2;
			vst1q_f64(C, a4); C += 2;
		}
	}

	for (; i < n; ++i)
		*C++ = *A++ * B;
}

TARGETNEO void MulNEO(double* A, double B, int64 n)
{
	int64 i = 0;

	if (n >= 8)
	{
		float64x2_t b = vdupq_n_f64(B);
		float64x2_t a1, a2, a3, a4;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			a1 = vld1q_f64(A); A += 2;
			a2 = vld1q_f64(A); A += 2;
			a3 = vld1q_f64(A); A += 2;
			a4 = vld1q_f64(A); A += 2;

			a1 = vmulq_f64(a1, b);
			a2 = vmulq_f64(a2, b);
			a3 = vmulq_f64(a3, b);
			a4 = vmulq_f64(a4, b);

			vst1q_f64(A - 8, a1);
			vst1q_f64(A - 6, a2);
			vst1q_f64(A - 4, a3);
			vst1q_f64(A - 2, a4);
		}
	}

	for (; i < n; ++i)
		*A++ *= B;
}

TARGETNEO void AddProdNEO(double* C, double* A, double* B, int64 n)
{
	int64 i = 0;

	if (n >= 8)
	{
		float64x2_t a1, a2, a3, a4;
		float64x2_t b1, b2, b3, b4;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			a1 = vld1q_f64(A); A += 2;
			a2 = vld1q_f64(A); A += 2;
			a3 = vld1q_f64(A); A += 2;
			a4 = vld1q_f64(A); A += 2;

			b1 = vld1q_f64(B); B += 2;
			b2 = vld1q_f64(B); B += 2;
			b3 = vld1q_f64(B); B += 2;
			b4 = vld1q_f64(B); B += 2;

			a1 = vmulq_f64(a1, b1);
			a2 = vmulq_f64(a2, b2);
			a3 = vmulq_f64(a3, b3);
			a4 = vmulq_f64(a4, b4);

			b1 = vld1q_f64(C); C += 2;
			b2 = vld1q_f64(C); C += 2;
			b3 = vld1q_f64(C); C += 2;
			b4 = vld1q_f64(C); C += 2;

			a1 = vaddq_f64(a1, b1);
			a2 = vaddq_f64(a2, b2);
			a3 = vaddq_f64(a3, b3);
			a4 = vaddq_f64(a4, b4);

			vst1q_f64(&C[-8], a1);
			vst1q_f64(&C[-6], a2);
			vst1q_f64(&C[-4], a3);
			vst1q_f64(&C[-2], a4);
		}
	}

	for (; i < n; ++i)
		*C++ += *A++ * *B++;
}

TARGETNEO void AddProdNEO(double* C, double* A, double B, int64 n)
{
	int64 i = 0;

	if (n >= 2)
	{
		float64x2_t b = vdupq_n_f64(B);

		for (int64 l1 = n - 2; i <= l1; i += 2)
		{
			vst1q_f64(C, vaddq_f64(vld1q_f64(C), vmulq_f64(vld1q_f64(A), b)));
			A += 2; C += 2;
		}
	}

	for (; i < n; ++i)
		*C++ += *A++ * B;
}

TARGETNEO void UnifyNEO(double* A, int64 n)
{
	int64 i = 0;
	double invsum = 1.0 / (SumNEO(A, n) + n * MIN_FREQ);

	if (n >= 8)
	{
		float64x2_t minv = vdupq_n_f64(MIN_FREQ), invs = vdupq_n_f64(invsum);

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			vst1q_f64(A, vmulq_f64(invs, vaddq_f64(vld1q_f64(A), minv))); A += 2;
			vst1q_f64(A, vmulq_f64(invs, vaddq_f64(vld1q_f64(A), minv))); A += 2;
			vst1q_f64(A, vmulq_f64(invs, vaddq_f64(vld1q_f64(A), minv))); A += 2;
			vst1q_f64(A, vmulq_f64(invs, vaddq_f64(vld1q_f64(A), minv))); A += 2;
		}
	}

	for (; i < n; ++i, ++A)
		*A = (*A + MIN_FREQ) * invsum;
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
