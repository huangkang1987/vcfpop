/* SSE Instruction Set Functions */

#include "vcfpop.h"

#ifndef __aarch64__

#ifndef _RNGSSE
/* Initialize rng */
TARGETSSE RNGSSE::RNGSSE()
{
}

/* Initialize rng */
TARGETSSE RNGSSE::RNGSSE(uint64 s)
{
	__m128i seed[4];

	seed[0] = _mm_set_epi64x(Hash64ULong(s + 1), Hash64ULong(s + 0));
	seed[1] = _mm_set_epi64x(Hash64ULong(s + 3), Hash64ULong(s + 2));
	seed[2] = _mm_set_epi64x(Hash64ULong(s + 5), Hash64ULong(s + 4));
	seed[3] = _mm_set_epi64x(Hash64ULong(s + 7), Hash64ULong(s + 6));

	x[0] = _mm_xor_si128(_mm_set1_epi64x(0x159A55E5075BCD15), seed[0]); 
	x[1] = _mm_xor_si128(_mm_set1_epi64x(0x159A55E5075BCD15), seed[1]);
	x[2] = _mm_xor_si128(_mm_set1_epi64x(0x159A55E5075BCD15), seed[2]);
	x[3] = _mm_xor_si128(_mm_set1_epi64x(0x159A55E5075BCD15), seed[3]);

	seed[0] = _mm_slli_epi64(seed[0], 6);
	seed[1] = _mm_slli_epi64(seed[1], 6);
	seed[2] = _mm_slli_epi64(seed[2], 6);
	seed[3] = _mm_slli_epi64(seed[3], 6);

	y[0] = _mm_xor_si128(_mm_set1_epi64x(0x054913331F123BB5), seed[0]); 
	y[1] = _mm_xor_si128(_mm_set1_epi64x(0x054913331F123BB5), seed[1]); 
	y[2] = _mm_xor_si128(_mm_set1_epi64x(0x054913331F123BB5), seed[2]); 
	y[3] = _mm_xor_si128(_mm_set1_epi64x(0x054913331F123BB5), seed[3]); 
}

/* Draw a uniform distriubted interger */
TARGETSSE void RNGSSE::XorShift128p(__m128i* re)
{
	__m128i a[4], b[4];

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

	a[0] = _mm_xor_si128(a[0], _mm_slli_epi64(a[0], 23));
	a[1] = _mm_xor_si128(a[1], _mm_slli_epi64(a[1], 23));
	a[2] = _mm_xor_si128(a[2], _mm_slli_epi64(a[2], 23));
	a[3] = _mm_xor_si128(a[3], _mm_slli_epi64(a[3], 23));

	a[0] = _mm_xor_si128(a[0], _mm_srli_epi64(a[0], 18));
	a[1] = _mm_xor_si128(a[1], _mm_srli_epi64(a[1], 18));
	a[2] = _mm_xor_si128(a[2], _mm_srli_epi64(a[2], 18));
	a[3] = _mm_xor_si128(a[3], _mm_srli_epi64(a[3], 18));

	a[0] = _mm_xor_si128(a[0], b[0]);
	a[1] = _mm_xor_si128(a[1], b[1]);
	a[2] = _mm_xor_si128(a[2], b[2]);
	a[3] = _mm_xor_si128(a[3], b[3]);

	a[0] = _mm_xor_si128(a[0], _mm_srli_epi64(b[0], 5));
	a[1] = _mm_xor_si128(a[1], _mm_srli_epi64(b[1], 5));
	a[2] = _mm_xor_si128(a[2], _mm_srli_epi64(b[2], 5));
	a[3] = _mm_xor_si128(a[3], _mm_srli_epi64(b[3], 5));

	y[0] = a[0];
	y[1] = a[1];
	y[2] = a[2];
	y[3] = a[3];

	re[0] = _mm_add_epi64(a[0], b[0]);
	re[1] = _mm_add_epi64(a[1], b[1]);
	re[2] = _mm_add_epi64(a[2], b[2]);
	re[3] = _mm_add_epi64(a[3], b[3]);
}

/* Draw a uniform distriubted real number */
TARGETSSE void RNGSSE::Uniform(__m128d* re)
{
	__m128d one = _mm_set1_pd(1.0);
	__m128i mask1 = _mm_set1_epi64x(0x000FFFFFFFFFFFFF);
	__m128i mask2 = _mm_set1_epi64x(0x3FF0000000000000);
	__m128i* r = (__m128i*)re;

	XorShift128p(r);

	r[0] = _mm_or_si128(_mm_and_si128(r[0], mask1), mask2);
	r[1] = _mm_or_si128(_mm_and_si128(r[1], mask1), mask2);
	r[2] = _mm_or_si128(_mm_and_si128(r[2], mask1), mask2);
	r[3] = _mm_or_si128(_mm_and_si128(r[3], mask1), mask2);

	re[0] = _mm_sub_pd(re[0], one);
	re[1] = _mm_sub_pd(re[1], one);
	re[2] = _mm_sub_pd(re[2], one);
	re[3] = _mm_sub_pd(re[3], one);
}

/* Draw a uniform distriubted real number */
TARGETSSE void RNGSSE::Poly(__m128d* a, __m128d* s, int n, __m128i* re)
{
	__m128d t[4];

	Uniform(t);

	__m128d t1 = _mm_mul_pd(t[0], s[0]);
	__m128d t2 = _mm_mul_pd(t[1], s[1]);
	__m128d t3 = _mm_mul_pd(t[2], s[2]);
	__m128d t4 = _mm_mul_pd(t[3], s[3]);
	__m128d f1 = _mm_set1_pd(0);
	__m128d f2 = f1;
	__m128d f3 = f1;
	__m128d f4 = f1;
	__m128i midx1 = _mm_set1_epi64x(n - 1);
	__m128i midx2 = midx1;
	__m128i midx3 = midx1;
	__m128i midx4 = midx1;
	__m128i nidx = _mm_setzero_si128();
	__m128i ninc = _mm_set1_epi64x(1);
	__m128d b1, b2, b3, b4;
	__m128d v1, v2, v3, v4;

	for (int i = 0; i < n; ++i)
	{
		v1 = a[i * 4 + 0];
		v2 = a[i * 4 + 1];
		v3 = a[i * 4 + 2];
		v4 = a[i * 4 + 3];

		b1 = _mm_cmplt_pd(t1, v1);
		b2 = _mm_cmplt_pd(t2, v2);
		b3 = _mm_cmplt_pd(t3, v3);
		b4 = _mm_cmplt_pd(t4, v4);

		t1 = _mm_sub_pd(t1, v1);
		t2 = _mm_sub_pd(t2, v2);
		t3 = _mm_sub_pd(t3, v3);
		t4 = _mm_sub_pd(t4, v4);

		b1 = _mm_andnot_pd(f1, b1);
		b2 = _mm_andnot_pd(f2, b2);
		b3 = _mm_andnot_pd(f3, b3);
		b4 = _mm_andnot_pd(f4, b4);

		f1 = _mm_or_pd(f1, b1);
		f2 = _mm_or_pd(f2, b2);
		f3 = _mm_or_pd(f3, b3);
		f4 = _mm_or_pd(f4, b4);

		midx1 = _mm_castpd_si128(_mm_blendv_pd(_mm_castsi128_pd(midx1), _mm_castsi128_pd(nidx), b1));//ok
		midx2 = _mm_castpd_si128(_mm_blendv_pd(_mm_castsi128_pd(midx2), _mm_castsi128_pd(nidx), b2));
		midx3 = _mm_castpd_si128(_mm_blendv_pd(_mm_castsi128_pd(midx3), _mm_castsi128_pd(nidx), b3));
		midx4 = _mm_castpd_si128(_mm_blendv_pd(_mm_castsi128_pd(midx4), _mm_castsi128_pd(nidx), b4));

		nidx = _mm_add_epi64(nidx, ninc);
	}

	re[0] = midx1;
	re[1] = midx2;
	re[2] = midx3;
	re[3] = midx4;
}

/* Draw a polynormial distriubted integer with propoirtions in natural logarithm */
TARGETSSE void RNGSSE::PolyLog(__m128d* a, int n, __m128i* re)
{
	//proportional polynomial distribution, will overwrite a
	__m128d maxv[4] = { _mm_set1_pd(-1e300), _mm_set1_pd(-1e300), _mm_set1_pd(-1e300), _mm_set1_pd(-1e300) };
	__m128d s[4] = { _mm_set1_pd(MIN_FREQ * n), _mm_set1_pd(MIN_FREQ * n), _mm_set1_pd(MIN_FREQ * n), _mm_set1_pd(MIN_FREQ * n) };
	double* af = (double*)a;

	for (int i = 0; i < n; ++i)
	{
		maxv[0] = _mm_max_pd(maxv[0], a[i * 4 + 0]);
		maxv[1] = _mm_max_pd(maxv[1], a[i * 4 + 1]);
		maxv[2] = _mm_max_pd(maxv[2], a[i * 4 + 2]);
		maxv[3] = _mm_max_pd(maxv[3], a[i * 4 + 3]);
	}

	for (int i = 0; i < n; ++i)
	{
		a[i * 4 + 0] = _mm_sub_pd(a[i * 4 + 0], maxv[0]);
		a[i * 4 + 1] = _mm_sub_pd(a[i * 4 + 1], maxv[1]);
		a[i * 4 + 2] = _mm_sub_pd(a[i * 4 + 2], maxv[2]);
		a[i * 4 + 3] = _mm_sub_pd(a[i * 4 + 3], maxv[3]);

		af[i * 8 + 0] = (af[i * 8 + 0] < -23) ? MIN_FREQ : exp(af[i * 8 + 0]);
		af[i * 8 + 1] = (af[i * 8 + 1] < -23) ? MIN_FREQ : exp(af[i * 8 + 1]);
		af[i * 8 + 2] = (af[i * 8 + 2] < -23) ? MIN_FREQ : exp(af[i * 8 + 2]);
		af[i * 8 + 3] = (af[i * 8 + 3] < -23) ? MIN_FREQ : exp(af[i * 8 + 3]);
		af[i * 8 + 4] = (af[i * 8 + 4] < -23) ? MIN_FREQ : exp(af[i * 8 + 4]);
		af[i * 8 + 5] = (af[i * 8 + 5] < -23) ? MIN_FREQ : exp(af[i * 8 + 5]);
		af[i * 8 + 6] = (af[i * 8 + 6] < -23) ? MIN_FREQ : exp(af[i * 8 + 6]);
		af[i * 8 + 7] = (af[i * 8 + 7] < -23) ? MIN_FREQ : exp(af[i * 8 + 7]);

		s[0] = _mm_add_pd(s[0], a[i * 4 + 0]);
		s[1] = _mm_add_pd(s[1], a[i * 4 + 1]);
		s[2] = _mm_add_pd(s[2], a[i * 4 + 2]);
		s[3] = _mm_add_pd(s[3], a[i * 4 + 3]);
	}

	__m128d t[4];

	Uniform(t);

	__m128d t1 = _mm_mul_pd(t[0], s[0]);
	__m128d t2 = _mm_mul_pd(t[1], s[1]);
	__m128d t3 = _mm_mul_pd(t[2], s[2]);
	__m128d t4 = _mm_mul_pd(t[3], s[3]);
	__m128d f1 = _mm_set1_pd(0);
	__m128d f2 = f1;
	__m128d f3 = f1;
	__m128d f4 = f1;
	__m128i midx1 = _mm_set1_epi64x(n - 1);
	__m128i midx2 = midx1;
	__m128i midx3 = midx1;
	__m128i midx4 = midx1;
	__m128i nidx = _mm_setzero_si128();
	__m128i ninc = _mm_set1_epi64x(1);
	__m128d b1, b2, b3, b4;
	__m128d v1, v2, v3, v4;

	for (int i = 0; i < n; ++i)
	{
		v1 = a[i * 4 + 0];
		v2 = a[i * 4 + 1];
		v3 = a[i * 4 + 2];
		v4 = a[i * 4 + 3];

		b1 = _mm_cmplt_pd(t1, v1);
		b2 = _mm_cmplt_pd(t2, v2);
		b3 = _mm_cmplt_pd(t3, v3);
		b4 = _mm_cmplt_pd(t4, v4);

		t1 = _mm_sub_pd(t1, v1);
		t2 = _mm_sub_pd(t2, v2);
		t3 = _mm_sub_pd(t3, v3);
		t4 = _mm_sub_pd(t4, v4);

		b1 = _mm_andnot_pd(f1, b1);
		b2 = _mm_andnot_pd(f2, b2);
		b3 = _mm_andnot_pd(f3, b3);
		b4 = _mm_andnot_pd(f4, b4);

		f1 = _mm_or_pd(f1, b1);
		f2 = _mm_or_pd(f2, b2);
		f3 = _mm_or_pd(f3, b3);
		f4 = _mm_or_pd(f4, b4);

		midx1 = _mm_castpd_si128(_mm_blendv_pd(_mm_castsi128_pd(midx1), _mm_castsi128_pd(nidx), b1));//ok
		midx2 = _mm_castpd_si128(_mm_blendv_pd(_mm_castsi128_pd(midx2), _mm_castsi128_pd(nidx), b2));
		midx3 = _mm_castpd_si128(_mm_blendv_pd(_mm_castsi128_pd(midx3), _mm_castsi128_pd(nidx), b3));
		midx4 = _mm_castpd_si128(_mm_blendv_pd(_mm_castsi128_pd(midx4), _mm_castsi128_pd(nidx), b4));

		nidx = _mm_add_epi64(nidx, ninc);
	}

	re[0] = midx1;
	re[1] = midx2;
	re[2] = midx3;
	re[3] = midx4;
}
#endif

TARGETSSE int64 GetMinIdxSSE(double* A, int64 n, double& val)
{
	int64 i = 0;
	val = 1e300;
	uint64 idx = (uint64)-1;

	if (n >= 8)
	{
		__m128d min1 = _mm_set1_pd(val);
		__m128d min2 = _mm_set1_pd(val);
		__m128d min3 = _mm_set1_pd(val);
		__m128d min4 = _mm_set1_pd(val);
		__m128i midx1 = _mm_set1_epi8((char)0xFF);
		__m128i midx2 = _mm_set1_epi8((char)0xFF);
		__m128i midx3 = _mm_set1_epi8((char)0xFF);
		__m128i midx4 = _mm_set1_epi8((char)0xFF);
		__m128i nidx1 = _mm_set_epi64x(1, 0);
		__m128i nidx2 = _mm_set_epi64x(3, 2);
		__m128i nidx3 = _mm_set_epi64x(5, 4);
		__m128i nidx4 = _mm_set_epi64x(7, 6);
		__m128i msep = _mm_set1_epi64x(8);
		__m128d f1, f2, f3, f4;
		__m128d v1, v2, v3, v4;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			v1 = _mm_loadu_pd(A); A += 2;
			v2 = _mm_loadu_pd(A); A += 2;
			v3 = _mm_loadu_pd(A); A += 2;
			v4 = _mm_loadu_pd(A); A += 2;

			f1 = _mm_cmpgt_pd(min1, v1);
			f2 = _mm_cmpgt_pd(min2, v2);
			f3 = _mm_cmpgt_pd(min3, v3);
			f4 = _mm_cmpgt_pd(min4, v4);

			min1 = _mm_blendv_pd(min1, v1, f1);//ok
			min2 = _mm_blendv_pd(min2, v2, f2);
			min3 = _mm_blendv_pd(min3, v3, f3);
			min4 = _mm_blendv_pd(min4, v4, f4);

			midx1 = _mm_castpd_si128(_mm_blendv_pd(_mm_castsi128_pd(midx1), _mm_castsi128_pd(nidx1), f1));//ok
			midx2 = _mm_castpd_si128(_mm_blendv_pd(_mm_castsi128_pd(midx2), _mm_castsi128_pd(nidx2), f2));
			midx3 = _mm_castpd_si128(_mm_blendv_pd(_mm_castsi128_pd(midx3), _mm_castsi128_pd(nidx3), f3));
			midx4 = _mm_castpd_si128(_mm_blendv_pd(_mm_castsi128_pd(midx4), _mm_castsi128_pd(nidx4), f4));

			nidx1 = _mm_add_epi64(nidx1, msep);
			nidx2 = _mm_add_epi64(nidx2, msep);
			nidx3 = _mm_add_epi64(nidx3, msep);
			nidx4 = _mm_add_epi64(nidx4, msep);
		}

		f1 = _mm_cmpgt_pd(min1, min2);
		f3 = _mm_cmpgt_pd(min3, min4);

		min1 = _mm_blendv_pd(min1, min2, f1);
		min3 = _mm_blendv_pd(min3, min4, f3);

		midx1 = _mm_castpd_si128(_mm_blendv_pd(_mm_castsi128_pd(midx1), _mm_castsi128_pd(midx2), f1));
		midx3 = _mm_castpd_si128(_mm_blendv_pd(_mm_castsi128_pd(midx3), _mm_castsi128_pd(midx4), f3));

		f1 = _mm_cmpgt_pd(min1, min3);
		min1 = _mm_blendv_pd(min1, min3, f1);
		midx1 = _mm_castpd_si128(_mm_blendv_pd(_mm_castsi128_pd(midx1), _mm_castsi128_pd(midx3), f1));

		for (int64 j = 0; j < 2; ++j)
		{
			if (simd_f64(min1, j) > val) continue;
			val = simd_f64(min1, j);
			idx = simd_u64(midx1, j);
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
	int64 i = 0;
	minv = 1e300;
	maxv = -1e300;

	if (n >= 8)
	{
		__m128d max1 = _mm_set1_pd(maxv);
		__m128d max2 = max1;
		__m128d max3 = max1;
		__m128d max4 = max1;
		__m128d min1 = _mm_set1_pd(minv);
		__m128d min2 = min1;
		__m128d min3 = min1;
		__m128d min4 = min1;
		__m128d v1, v2, v3, v4;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			v1 = _mm_loadu_pd(A); A += 2;
			v2 = _mm_loadu_pd(A); A += 2;
			v3 = _mm_loadu_pd(A); A += 2;
			v4 = _mm_loadu_pd(A); A += 2;

			min1 = _mm_min_pd(min1, v1);
			min1 = _mm_min_pd(min1, v1);
			min2 = _mm_min_pd(min2, v2);
			min3 = _mm_min_pd(min3, v3);
			min4 = _mm_min_pd(min4, v4);

			max1 = _mm_max_pd(max1, v1);
			max2 = _mm_max_pd(max2, v2);
			max3 = _mm_max_pd(max3, v3);
			max4 = _mm_max_pd(max4, v4);
		}

		min1 = _mm_min_pd(_mm_min_pd(min1, min2), _mm_min_pd(min3, min4));
		max1 = _mm_max_pd(_mm_max_pd(max1, max2), _mm_max_pd(max3, max4));

		for (int64 j = 0; j < 2; ++j)
		{
			if (simd_f64(min1, j) < minv) minv = simd_f64(min1, j);
			if (simd_f64(max1, j) > maxv) maxv = simd_f64(max1, j);
		}
	}

	for (; i < n; ++i, ++A)
	{
		if (*A < minv) minv = *A;
		if (*A > maxv) maxv = *A;
	}
}

TARGETSSE double GetMaxValSSE(double* A, int64 n)
{
	int64 i = 0;
	double val = -1e300;

	if (n >= 8)
	{
		__m128d max1 = _mm_set1_pd(val);
		__m128d max2 = _mm_set1_pd(val);
		__m128d max3 = _mm_set1_pd(val);
		__m128d max4 = _mm_set1_pd(val);
		__m128d v1, v2, v3, v4;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			v1 = _mm_loadu_pd(A); A += 2;
			v2 = _mm_loadu_pd(A); A += 2;
			v3 = _mm_loadu_pd(A); A += 2;
			v4 = _mm_loadu_pd(A); A += 2;

			max1 = _mm_max_pd(max1, v1);
			max2 = _mm_max_pd(max2, v2);
			max3 = _mm_max_pd(max3, v3);
			max4 = _mm_max_pd(max4, v4);
		}

		max1 = _mm_max_pd(_mm_max_pd(max1, max2), _mm_max_pd(max3, max4));

		for (int64 j = 0; j < 2; ++j)
		{
			if (simd_f64(max1, j) < val) continue;
			val = simd_f64(max1, j);
		}
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
	double val = -1e300;

	if (n >= 8)
	{
		__m128d max1 = _mm_set1_pd(val);
		__m128d max2 = _mm_set1_pd(val);
		__m128d max3 = _mm_set1_pd(val);
		__m128d max4 = _mm_set1_pd(val);
		__m128d v1, v2, v3, v4;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			v1 = _mm_set_pd(A[sep], A[0]); A += sep * 2;
			v2 = _mm_set_pd(A[sep], A[0]); A += sep * 2;
			v3 = _mm_set_pd(A[sep], A[0]); A += sep * 2;
			v4 = _mm_set_pd(A[sep], A[0]); A += sep * 2;

			max1 = _mm_max_pd(max1, v1);
			max2 = _mm_max_pd(max2, v2);
			max3 = _mm_max_pd(max3, v3);
			max4 = _mm_max_pd(max4, v4);
		}

		max1 = _mm_max_pd(_mm_max_pd(max1, max2), _mm_max_pd(max3, max4));

		for (int64 j = 0; j < 2; ++j)
		{
			if (simd_f64(max1, j) < val) continue;
			val = simd_f64(max1, j);
		}
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
	int64 i = 0;
	double val = 1e300;

	if (n >= 8)
	{
		__m128d min1 = _mm_set1_pd(val);
		__m128d min2 = _mm_set1_pd(val);
		__m128d min3 = _mm_set1_pd(val);
		__m128d min4 = _mm_set1_pd(val);
		__m128d v1, v2, v3, v4;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			v1 = _mm_loadu_pd(A); A += 2;
			v2 = _mm_loadu_pd(A); A += 2;
			v3 = _mm_loadu_pd(A); A += 2;
			v4 = _mm_loadu_pd(A); A += 2;

			min1 = _mm_min_pd(min1, v1);
			min2 = _mm_min_pd(min2, v2);
			min3 = _mm_min_pd(min3, v3);
			min4 = _mm_min_pd(min4, v4);
		}

		min1 = _mm_min_pd(_mm_min_pd(min1, min2), _mm_min_pd(min3, min4));

		for (int64 j = 0; j < 2; ++j)
		{
			if (simd_f64(min1, j) > val) continue;
			val = simd_f64(min1, j);
		}

		for (; i < n; ++i, ++A)
		{
			if (*A > val) continue;
			val = *A;
		}
	}

	return val;
}

TARGETSSE int64 GetMinValSSE(int64* A, int64 n)
{
	int64 i = 0;
	int64 val = 0x7FFFFFFFFFFFFFFF;

	if (n >= 8)
	{
		__m128i min1 = _mm_set1_epi64x(0x7FFFFFFFFFFFFFFF);
		__m128i min2 = _mm_set1_epi64x(0x7FFFFFFFFFFFFFFF);
		__m128i min3 = _mm_set1_epi64x(0x7FFFFFFFFFFFFFFF);
		__m128i min4 = _mm_set1_epi64x(0x7FFFFFFFFFFFFFFF);
		__m128i v1, v2, v3, v4;
		__m128i f1, f2, f3, f4;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			v1 = _mm_loadu_si128((const __m128i*)A); A += 2;
			v2 = _mm_loadu_si128((const __m128i*)A); A += 2;
			v3 = _mm_loadu_si128((const __m128i*)A); A += 2;
			v4 = _mm_loadu_si128((const __m128i*)A); A += 2;

			f1 = _mm_cmpgt_epi64(min1, v1);
			f2 = _mm_cmpgt_epi64(min2, v2);
			f3 = _mm_cmpgt_epi64(min3, v3);
			f4 = _mm_cmpgt_epi64(min4, v4);

			min1 = _mm_castpd_si128(_mm_blendv_pd(_mm_castsi128_pd(min1), _mm_castsi128_pd(v1), _mm_castsi128_pd(f1)));
			min2 = _mm_castpd_si128(_mm_blendv_pd(_mm_castsi128_pd(min2), _mm_castsi128_pd(v2), _mm_castsi128_pd(f2)));
			min3 = _mm_castpd_si128(_mm_blendv_pd(_mm_castsi128_pd(min3), _mm_castsi128_pd(v3), _mm_castsi128_pd(f3)));
			min4 = _mm_castpd_si128(_mm_blendv_pd(_mm_castsi128_pd(min4), _mm_castsi128_pd(v4), _mm_castsi128_pd(f4)));
		}

		f1 = _mm_cmpgt_epi64(min1, min2);
		f3 = _mm_cmpgt_epi64(min3, min4);

		min1 = _mm_castpd_si128(_mm_blendv_pd(_mm_castsi128_pd(min1), _mm_castsi128_pd(min2), _mm_castsi128_pd(f1)));
		min3 = _mm_castpd_si128(_mm_blendv_pd(_mm_castsi128_pd(min3), _mm_castsi128_pd(min4), _mm_castsi128_pd(f3)));

		f1 = _mm_cmpgt_epi64(min1, min3);
		min1 = _mm_castpd_si128(_mm_blendv_pd(_mm_castsi128_pd(min1), _mm_castsi128_pd(min3), _mm_castsi128_pd(f1)));

		for (int64 j = 0; j < 2; ++j)
		{
			if (simd_u64(min1, j) > val) continue;
			val = simd_u64(min1, j);
		}
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

TARGETSSE void ChargeLogSSE(int64& slog, double& prod, __m128d& val)
{
	AddExponentSSE(slog, val);
	prod = prod * simd_f64(val, 0) * simd_f64(val, 1);

	if (prod < DOUBLE_UNDERFLOW || prod > DOUBLE_OVERFLOW) [[unlikely]]
		AddExponent(slog, prod);
}

TARGETSSE double LogProdSSE(double* A, int64 n)
{
	int64 i = 0;
	int64 slog = 0; double prod = 1;

	if (n >= 2)
	{
		__m128d pd = _mm_set1_pd(1.0), dunder = _mm_set1_pd(DOUBLE_UNDERFLOW), dover = _mm_set1_pd(DOUBLE_OVERFLOW);
		__m128i flag, mask = _mm_set1_epi64x(0xFFFFFFFFFFFFFFFF);

		for (int64 l1 = n - 2; i <= l1; i += 2)
		{
			pd = _mm_mul_pd(pd, _mm_loadu_pd(A)); A += 2;
			flag = _mm_castpd_si128(_mm_or_pd(_mm_cmplt_pd(pd, dunder), _mm_cmplt_pd(dover, pd)));

			if (_mm_test_all_zeros(mask, flag)) [[likely]] continue;

			AddExponentSSE(slog, pd);
		}

		ChargeLogSSE(slog, prod, pd);
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

	if (n >= 2)
	{
		_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
		_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;

		__m128d pd = _mm_set1_pd(1.0), dunder = _mm_set1_pd(DOUBLE_UNDERFLOW), dover = _mm_set1_pd(DOUBLE_OVERFLOW);
		__m128i flag, mask = _mm_set1_epi64x(0xFFFFFFFFFFFFFFFF);

		for (int64 l1 = n - 2; i <= l1; i += 2)
		{
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;

			pd = _mm_mul_pd(pd, _mm_set_pd(A[-3 * sep], A[-4 * sep]));
			flag = _mm_castpd_si128(_mm_or_pd(_mm_cmplt_pd(pd, dunder), _mm_cmplt_pd(dover, pd)));

			if (_mm_test_all_zeros(mask, flag)) [[likely]] continue;

			AddExponentSSE(slog, pd);
		}

		ChargeLogSSE(slog, prod, pd);
		A -= sep * 2;
	}

	for (; i < n; ++i, A += sep)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);
	return prod;
}

TARGETSSE double LogProdDivSSE(double* A, double* B, int64 n, int64 sep)
{
	int64 i = 0;
	int64 slog = 0; double prod = 1;

	if (n >= 2)
	{
		_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
		_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;

		__m128d pd = _mm_set1_pd(1.0), dunder = _mm_set1_pd(DOUBLE_UNDERFLOW), dover = _mm_set1_pd(DOUBLE_OVERFLOW);
		__m128i flag, mask = _mm_set1_epi64x(0xFFFFFFFFFFFFFFFF);

		for (int64 l1 = n - 2; i <= l1; i += 2)
		{
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;

			pd = _mm_mul_pd(pd, _mm_div_pd(_mm_set_pd(A[-3 * sep], A[-4 * sep]), _mm_set_pd(B[-3 * sep], B[-4 * sep])));
			flag = _mm_castpd_si128(_mm_or_pd(_mm_cmplt_pd(pd, dunder), _mm_cmplt_pd(dover, pd)));

			if (_mm_test_all_zeros(mask, flag)) [[likely]] continue;

			AddExponentSSE(slog, pd);
		}

		ChargeLogSSE(slog, prod, pd);
		A -= sep * 2; 
		B -= sep * 2;
	}

	for (; i < n; ++i, A += sep, B += sep)
		ChargeLog(slog, prod, *A / *B);

	CloseLog(slog, prod);
	return prod;
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
	int64 i = 0;
	double re = 0;

	if (n >= 8)
	{
		__m128d s1 = _mm_setzero_pd(), s2 = _mm_setzero_pd(), s3 = _mm_setzero_pd(), s4 = _mm_setzero_pd();
		__m128d v1, v2, v3, v4;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			v1 = _mm_loadu_pd(A); A += 2;
			v2 = _mm_loadu_pd(A); A += 2;
			v3 = _mm_loadu_pd(A); A += 2;
			v4 = _mm_loadu_pd(A); A += 2;

			s1 = _mm_add_pd(s1, v1);
			s2 = _mm_add_pd(s2, v2);
			s3 = _mm_add_pd(s3, v3);
			s4 = _mm_add_pd(s4, v4);
		}

		s1 = _mm_add_pd(_mm_add_pd(s1, s2), _mm_add_pd(s3, s4));
		s1 = _mm_hadd_pd(s1, s1);
		re = simd_f64(s1, 0);
	}

	for (; i < n; ++i)
		re += *A++;

	return re;
}

TARGETSSE int64 SumSSE(byte* A, int64 n)
{
	uint64 re = 0;
	int64 i = 0;

	if (n >= 64)
	{
		__m128i s1 = _mm_setzero_si128();
		__m128i s2 = _mm_setzero_si128();
		__m128i s3 = _mm_setzero_si128();
		__m128i s4 = _mm_setzero_si128();
		__m128i z = _mm_setzero_si128();
		__m128i v1, v2, v3, v4;

		for (int64 l1 = n - 64; i <= l1; i += 64)
		{
			v1 = _mm_loadu_si128((__m128i*)A); A += 16;
			v2 = _mm_loadu_si128((__m128i*)A); A += 16;
			v3 = _mm_loadu_si128((__m128i*)A); A += 16;
			v4 = _mm_loadu_si128((__m128i*)A); A += 16;

			s1 = _mm_add_epi64(s1, _mm_sad_epu8(v1, z));
			s2 = _mm_add_epi64(s2, _mm_sad_epu8(v2, z));
			s3 = _mm_add_epi64(s3, _mm_sad_epu8(v3, z));
			s4 = _mm_add_epi64(s4, _mm_sad_epu8(v4, z));
		}

		s1 = _mm_add_epi64(_mm_add_epi64(s1, s2), _mm_add_epi64(s3, s4));
		re += simd_u64(s1, 0) + simd_u64(s1, 1);
	}

	for (; i < n; ++i)
		re += *A++;

	return re;
}

TARGETSSE double SumSSE(double* A, int64 n, int64 sep)
{
	int64 i = 0;
	double re = 0;

	if (n >= 4)
	{
		_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
		_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
		_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
		_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;

		__m128d s1 = _mm_setzero_pd(), s2 = _mm_setzero_pd();

		for (int64 l1 = n - 4; i <= l1; i += 4)
		{
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;

			s1 = _mm_add_pd(s1, _mm_set_pd(A[-7 * sep], A[-8 * sep]));
			s2 = _mm_add_pd(s2, _mm_set_pd(A[-5 * sep], A[-6 * sep]));
		}

		s1 = _mm_add_pd(s1, s2);
		//s1 = _mm_hadd_pd(s1, s1);
		re = simd_f64(s1, 0) + simd_f64(s1, 1);
		A -= sep * 4;
	}

	for (; i < n; ++i, A += sep)
		re += *A;

	return re;
}

TARGETSSE void SumSSE(double* A, double** B, int64 k, int64 n)
{
	int64 i = 0;

	if (n >= 8)
	{
		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			__m128d a1 = _mm_setzero_pd();
			__m128d a2 = _mm_setzero_pd();
			__m128d a3 = _mm_setzero_pd();
			__m128d a4 = _mm_setzero_pd();

			for (int64 j = 0; j < k; ++j)
			{
				a1 = _mm_add_pd(a1, _mm_loadu_pd(&B[j][i + 0]));
				a2 = _mm_add_pd(a2, _mm_loadu_pd(&B[j][i + 2]));
				a3 = _mm_add_pd(a3, _mm_loadu_pd(&B[j][i + 4]));
				a4 = _mm_add_pd(a4, _mm_loadu_pd(&B[j][i + 6]));
			}

			_mm_storeu_pd(&A[i + 0], a1);
			_mm_storeu_pd(&A[i + 2], a2);
			_mm_storeu_pd(&A[i + 4], a3);
			_mm_storeu_pd(&A[i + 6], a4);
		}
	}

	for (; i < n; ++i)
	{
		A[i] = 0;
		for (int64 j = 0; j < k; ++j)
			A[i] += B[j][i];
	}
}

TARGETSSE double ProdSSE(double* A, int64 n)
{
	int64 i = 0;
	double re = 1;

	if (n >= 16)
	{
		__m128d pd1 = _mm_set1_pd(1.0);
		__m128d pd2 = _mm_set1_pd(1.0);
		__m128d pd3 = _mm_set1_pd(1.0);
		__m128d pd4 = _mm_set1_pd(1.0);
		__m128d pd5 = _mm_set1_pd(1.0);
		__m128d pd6 = _mm_set1_pd(1.0);
		__m128d pd7 = _mm_set1_pd(1.0);
		__m128d pd8 = _mm_set1_pd(1.0);

		for (int64 l1 = n - 16; i <= l1; i += 16)
		{
			pd1 = _mm_mul_pd(pd1, _mm_loadu_pd(A));
			A += 2;

			pd2 = _mm_mul_pd(pd2, _mm_loadu_pd(A));
			A += 2;

			pd3 = _mm_mul_pd(pd3, _mm_loadu_pd(A));
			A += 2;

			pd4 = _mm_mul_pd(pd4, _mm_loadu_pd(A));
			A += 2;

			pd5 = _mm_mul_pd(pd5, _mm_loadu_pd(A));
			A += 2;

			pd6 = _mm_mul_pd(pd6, _mm_loadu_pd(A));
			A += 2;

			pd7 = _mm_mul_pd(pd7, _mm_loadu_pd(A));
			A += 2;

			pd8 = _mm_mul_pd(pd8, _mm_loadu_pd(A));
			A += 2;
		}

		pd1 = _mm_mul_pd(pd1, pd2);
		pd3 = _mm_mul_pd(pd3, pd4);
		pd5 = _mm_mul_pd(pd5, pd6);
		pd7 = _mm_mul_pd(pd7, pd8);

		pd1 = _mm_mul_pd(pd1, pd3);
		pd5 = _mm_mul_pd(pd5, pd7);

		pd1 = _mm_mul_pd(pd1, pd5);

		re = simd_f64(pd1, 0) * simd_f64(pd1, 1);
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
		_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
		_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
		_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
		_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;

		__m128d pd1 = _mm_set1_pd(1.0), pd2 = _mm_set1_pd(1.0);
		__m128d pd3 = _mm_set1_pd(1.0), pd4 = _mm_set1_pd(1.0);

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			pd1 = _mm_mul_pd(pd1, _mm_set_pd(A[-5 * sep], A[-6 * sep]));

			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			pd2 = _mm_mul_pd(pd2, _mm_set_pd(A[-5 * sep], A[-6 * sep]));

			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			pd3 = _mm_mul_pd(pd3, _mm_set_pd(A[-5 * sep], A[-6 * sep]));

			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			pd4 = _mm_mul_pd(pd4, _mm_set_pd(A[-5 * sep], A[-6 * sep]));
		}

		pd1 = _mm_mul_pd(_mm_mul_pd(pd1, pd2), _mm_mul_pd(pd3, pd4));
		re = simd_f64(pd1, 0) * simd_f64(pd1, 1);
		A -= sep * 4;
	}

	for (; i < n; ++i, A += sep)
		re *= *A;

	return re;
}

TARGETSSE double SumSquareSSE(double* A, int64 n)
{
	int64 i = 0;
	double re = 0;

	if (n >= 8)
	{
		__m128d s1 = _mm_setzero_pd();
		__m128d s2 = _mm_setzero_pd();
		__m128d s3 = _mm_setzero_pd();
		__m128d s4 = _mm_setzero_pd();
		__m128d a1, a2, a3, a4;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			a1 = _mm_loadu_pd(A); A += 2;
			a2 = _mm_loadu_pd(A); A += 2;
			a3 = _mm_loadu_pd(A); A += 2;
			a4 = _mm_loadu_pd(A); A += 2;

			s1 = _mm_add_pd(s1, _mm_mul_pd(a1, a1));
			s2 = _mm_add_pd(s2, _mm_mul_pd(a2, a2));
			s3 = _mm_add_pd(s3, _mm_mul_pd(a3, a3));
			s4 = _mm_add_pd(s4, _mm_mul_pd(a4, a4));
		}

		s1 = _mm_add_pd(_mm_add_pd(s1, s2), _mm_add_pd(s3, s4));
		s1 = _mm_hadd_pd(s1, s1);
		re = simd_f64(s1, 0);
	}

	for (; i < n; ++i, ++A)
		re += *A * *A;

	return re;
}

TARGETSSE int64 SumSquareSSE(byte* A, int64 n)
{
	int64 i = 0;
	uint64 re = 0;

	if (n >= 64)
	{
		__m128i s1 = _mm_setzero_si128();
		__m128i s2 = _mm_setzero_si128();
		__m128i s3 = _mm_setzero_si128();
		__m128i s4 = _mm_setzero_si128();
		__m128i a1, a2, a3, a4;

		for (int64 l1 = n - 64; i <= l1; i += 64)
		{
			a1 = _mm_loadu_si128((__m128i*)A); A += 16;
			a2 = _mm_loadu_si128((__m128i*)A); A += 16;
			a3 = _mm_loadu_si128((__m128i*)A); A += 16;
			a4 = _mm_loadu_si128((__m128i*)A); A += 16;

			a1 = _mm_maddubs_epi16(a1, a1);
			a2 = _mm_maddubs_epi16(a2, a2);
			a3 = _mm_maddubs_epi16(a3, a3);
			a4 = _mm_maddubs_epi16(a4, a4);

			s1 = _mm_add_epi32(s1, _mm_add_epi32(_mm_cvtepu8_epi16(a1),
				_mm_cvtepu8_epi16(_mm_srli_si128(a1, 8))));//should use _mm_srli_si128
			s2 = _mm_add_epi32(s2, _mm_add_epi32(_mm_cvtepu8_epi16(a2),
				_mm_cvtepu8_epi16(_mm_srli_si128(a2, 8))));//should use _mm_srli_si128
			s3 = _mm_add_epi32(s3, _mm_add_epi32(_mm_cvtepu8_epi16(a3),
				_mm_cvtepu8_epi16(_mm_srli_si128(a3, 8))));//should use _mm_srli_si128
			s4 = _mm_add_epi32(s4, _mm_add_epi32(_mm_cvtepu8_epi16(a4),
				_mm_cvtepu8_epi16(_mm_srli_si128(a4, 8))));//should use _mm_srli_si128
		}

		s1 = _mm_add_epi32(_mm_add_epi32(s1, s2), _mm_add_epi32(s3, s4));
		s1 = _mm_hadd_epi32(s1, s1);
		s1 = _mm_hadd_epi32(s1, s1);
		re += simd_u32(s1, 0);
	}

	for (; i < n; ++i, ++A)
		re += *A * *A;

	return re;
}

TARGETSSE void SumSumSquareSSE(double* A, int64 n, double& sum, double& sumsq)
{
	int64 i = 0;
	sum = sumsq = 0;

	if (n >= 8)
	{
		__m128d s1 = _mm_setzero_pd(), sq1 = _mm_setzero_pd();
		__m128d s2 = _mm_setzero_pd(), sq2 = _mm_setzero_pd();
		__m128d a;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			a = _mm_loadu_pd(A); A += 2;
			s1 = _mm_add_pd(s1, a);
			sq1 = _mm_add_pd(sq1, _mm_mul_pd(a, a));

			a = _mm_loadu_pd(A); A += 2;
			s2 = _mm_add_pd(s2, a);
			sq2 = _mm_add_pd(sq2, _mm_mul_pd(a, a));

			a = _mm_loadu_pd(A); A += 2;
			s1 = _mm_add_pd(s1, a);
			sq1 = _mm_add_pd(sq1, _mm_mul_pd(a, a));

			a = _mm_loadu_pd(A); A += 2;
			s2 = _mm_add_pd(s2, a);
			sq2 = _mm_add_pd(sq2, _mm_mul_pd(a, a));
		}

		s1 = _mm_add_pd(s1, s2);
		sq1 = _mm_add_pd(sq1, sq2);

		s1 = _mm_hadd_pd(s1, s1);
		sum = simd_f64(s1, 0);
		sq1 = _mm_hadd_pd(sq1, sq1);
		sumsq = simd_f64(sq1, 0);
	}

	for (; i < n; ++i, ++A)
	{
		sum += *A;
		sumsq += *A * *A;
	}
}

/* re = Sum(A1[i++] * B[j += sep]) / Sum(A2[i++] * B[j += sep]) */
TARGETSSE double SumProdDivSSE(double* A1, double* A2, double* B, int64 sep, int64 n)
{
	int64 i = 0;
	double re1 = 0, re2 = 0;

	if (n >= 4)
	{
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;

		__m128d s1 = _mm_setzero_pd(), s2 = _mm_setzero_pd(), b;

		for (int64 l1 = n - 4; i <= l1; i += 4)
		{
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			b = _mm_set_pd(B[-9 * sep], B[-10 * sep]);
			s1 = _mm_add_pd(s1, _mm_mul_pd(_mm_loadu_pd(A1), b)); A1 += 2;
			s2 = _mm_add_pd(s2, _mm_mul_pd(_mm_loadu_pd(A2), b)); A2 += 2;

			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			b = _mm_set_pd(B[-9 * sep], B[-10 * sep]);
			s1 = _mm_add_pd(s1, _mm_mul_pd(_mm_loadu_pd(A1), b)); A1 += 2;
			s2 = _mm_add_pd(s2, _mm_mul_pd(_mm_loadu_pd(A2), b)); A2 += 2;
		}

		s1 = _mm_hadd_pd(s1, s1);
		s2 = _mm_hadd_pd(s2, s2);
		re1 = simd_f64(s1, 0);
		re2 = simd_f64(s2, 0);
		B -= sep * 8;
	}

	for (; i < n; ++i, A1++, A2++, B += sep)
	{
		re1 += *A1 * *B;
		re2 += *A2 * *B;
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

		s = _mm_hadd_pd(s, s);
		re = simd_f64(s, 0);
		B -= sep * 2;
	}

	for (; i < n; ++i, A++, B += sep)
		re += *A * *B;

	return re;
}

TARGETSSE double SumProdSSE(double* A, double* B, int64 n)
{
	int64 i = 0;
	double re = 0;

	if (n >= 8)
	{
		__m128d s1 = _mm_setzero_pd();
		__m128d s2 = _mm_setzero_pd();
		__m128d s3 = _mm_setzero_pd();
		__m128d s4 = _mm_setzero_pd();
		__m128d a1, a2, a3, a4;
		__m128d b1, b2, b3, b4;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			a1 = _mm_loadu_pd(A); A += 2;
			a2 = _mm_loadu_pd(A); A += 2;
			a3 = _mm_loadu_pd(A); A += 2;
			a4 = _mm_loadu_pd(A); A += 2;

			b1 = _mm_loadu_pd(B); B += 2;
			b2 = _mm_loadu_pd(B); B += 2;
			b3 = _mm_loadu_pd(B); B += 2;
			b4 = _mm_loadu_pd(B); B += 2;

			s1 = _mm_add_pd(s1, _mm_mul_pd(a1, b1));
			s2 = _mm_add_pd(s2, _mm_mul_pd(a2, b2));
			s3 = _mm_add_pd(s3, _mm_mul_pd(a3, b3));
			s4 = _mm_add_pd(s4, _mm_mul_pd(a4, b4));
		}

		s1 = _mm_add_pd(_mm_add_pd(s1, s2), _mm_add_pd(s3, s4));
		s1 = _mm_hadd_pd(s1, s1);
		re = simd_f64(s1, 0);
	}

	for (; i < n; ++i, ++A, ++B)
		re += *A * *B;

	return re;
}

TARGETSSE void AddSSE(double* A, double* B, int64 n)
{
	int64 i = 0;

	if (n >= 8)
	{
		__m128d a1, a2, a3, a4;
		__m128d b1, b2, b3, b4;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			a1 = _mm_loadu_pd(A); A += 2;
			a2 = _mm_loadu_pd(A); A += 2;
			a3 = _mm_loadu_pd(A); A += 2;
			a4 = _mm_loadu_pd(A); A += 2;

			b1 = _mm_loadu_pd(B); B += 2;
			b2 = _mm_loadu_pd(B); B += 2;
			b3 = _mm_loadu_pd(B); B += 2;
			b4 = _mm_loadu_pd(B); B += 2;

			a1 = _mm_add_pd(a1, b1);
			a2 = _mm_add_pd(a2, b2);
			a3 = _mm_add_pd(a3, b3);
			a4 = _mm_add_pd(a4, b4);

			_mm_storeu_pd(A - 8, a1);
			_mm_storeu_pd(A - 6, a2);
			_mm_storeu_pd(A - 4, a3);
			_mm_storeu_pd(A - 2, a4);
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGETSSE void AddSSE(int64* A, int64* B, int64 n)
{
	int64 i = 0;

	if (n >= 8)
	{
		__m128i a1, a2, a3, a4;
		__m128i b1, b2, b3, b4;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			a1 = _mm_loadu_si128((__m128i*)A); A += 2;
			a2 = _mm_loadu_si128((__m128i*)A); A += 2;
			a3 = _mm_loadu_si128((__m128i*)A); A += 2;
			a4 = _mm_loadu_si128((__m128i*)A); A += 2;

			b1 = _mm_loadu_si128((__m128i*)B); B += 2;
			b2 = _mm_loadu_si128((__m128i*)B); B += 2;
			b3 = _mm_loadu_si128((__m128i*)B); B += 2;
			b4 = _mm_loadu_si128((__m128i*)B); B += 2;

			a1 = _mm_add_epi64(a1, b1);
			a2 = _mm_add_epi64(a2, b2);
			a3 = _mm_add_epi64(a3, b3);
			a4 = _mm_add_epi64(a4, b4);

			_mm_storeu_si128((__m128i*)(A - 8), a1);
			_mm_storeu_si128((__m128i*)(A - 6), a2);
			_mm_storeu_si128((__m128i*)(A - 4), a3);
			_mm_storeu_si128((__m128i*)(A - 2), a4);
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGETSSE void AddSSE(int* A, int* B, int64 n)
{
	int64 i = 0;

	if (n >= 16)
	{
		__m128i a1, a2, a3, a4;
		__m128i b1, b2, b3, b4;

		for (int64 l1 = n - 16; i <= l1; i += 16)
		{
			a1 = _mm_loadu_si128((__m128i*)A); A += 4;
			a2 = _mm_loadu_si128((__m128i*)A); A += 4;
			a3 = _mm_loadu_si128((__m128i*)A); A += 4;
			a4 = _mm_loadu_si128((__m128i*)A); A += 4;

			b1 = _mm_loadu_si128((__m128i*)B); B += 4;
			b2 = _mm_loadu_si128((__m128i*)B); B += 4;
			b3 = _mm_loadu_si128((__m128i*)B); B += 4;
			b4 = _mm_loadu_si128((__m128i*)B); B += 4;

			a1 = _mm_add_epi32(a1, b1);
			a2 = _mm_add_epi32(a2, b2);
			a3 = _mm_add_epi32(a3, b3);
			a4 = _mm_add_epi32(a4, b4);

			_mm_storeu_si128((__m128i*)(A - 16), a1);
			_mm_storeu_si128((__m128i*)(A - 12), a2);
			_mm_storeu_si128((__m128i*)(A - 8), a3);
			_mm_storeu_si128((__m128i*)(A - 4), a4);
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGETSSE void AddSSE(double* A, double B, int64 n)
{
	int64 i = 0;

	if (n >= 8)
	{
		__m128d b = _mm_set1_pd(B);
		__m128d a1, a2, a3, a4;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			a1 = _mm_loadu_pd(A); A += 2;
			a2 = _mm_loadu_pd(A); A += 2;
			a3 = _mm_loadu_pd(A); A += 2;
			a4 = _mm_loadu_pd(A); A += 2;

			a1 = _mm_add_pd(a1, b);
			a2 = _mm_add_pd(a2, b);
			a3 = _mm_add_pd(a3, b);
			a4 = _mm_add_pd(a4, b);

			_mm_storeu_pd(A - 8, a1);
			_mm_storeu_pd(A - 6, a2);
			_mm_storeu_pd(A - 4, a3);
			_mm_storeu_pd(A - 2, a4);
		}
	}

	for (; i < n; ++i, A++)
		*A += B;
}

TARGETSSE void MulSSE(double* C, double* A, double* B, int64 n)
{
	int64 i = 0;

	if (n >= 8)
	{
		__m128d a1, a2, a3, a4;
		__m128d b1, b2, b3, b4;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			a1 = _mm_loadu_pd(A); A += 2;
			a2 = _mm_loadu_pd(A); A += 2;
			a3 = _mm_loadu_pd(A); A += 2;
			a4 = _mm_loadu_pd(A); A += 2;

			b1 = _mm_loadu_pd(B); B += 2;
			b2 = _mm_loadu_pd(B); B += 2;
			b3 = _mm_loadu_pd(B); B += 2;
			b4 = _mm_loadu_pd(B); B += 2;

			a1 = _mm_mul_pd(a1, b1);
			a2 = _mm_mul_pd(a2, b2);
			a3 = _mm_mul_pd(a3, b3);
			a4 = _mm_mul_pd(a4, b4);

			_mm_storeu_pd(C, a1); C += 2;
			_mm_storeu_pd(C, a2); C += 2;
			_mm_storeu_pd(C, a3); C += 2;
			_mm_storeu_pd(C, a4); C += 2;
		}
	}

	for (; i < n; ++i)
		*C++ = *A++ * *B++;
}

TARGETSSE void MulSSE(double* C, double* A, double B, int64 n)
{
	int64 i = 0;

	if (n >= 8)
	{
		__m128d b = _mm_set1_pd(B);
		__m128d a1, a2, a3, a4;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			a1 = _mm_loadu_pd(A); A += 2;
			a2 = _mm_loadu_pd(A); A += 2;
			a3 = _mm_loadu_pd(A); A += 2;
			a4 = _mm_loadu_pd(A); A += 2;

			a1 = _mm_mul_pd(a1, b);
			a2 = _mm_mul_pd(a2, b);
			a3 = _mm_mul_pd(a3, b);
			a4 = _mm_mul_pd(a4, b);

			_mm_storeu_pd(C, a1); C += 2;
			_mm_storeu_pd(C, a2); C += 2;
			_mm_storeu_pd(C, a3); C += 2;
			_mm_storeu_pd(C, a4); C += 2;
		}
	}

	for (; i < n; ++i)
		*C++ = *A++ * B;
}

TARGETSSE void MulSSE(double* A, double B, int64 n)
{
	int64 i = 0;

	if (n >= 8)
	{
		__m128d b = _mm_set1_pd(B);
		__m128d a1, a2, a3, a4;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			a1 = _mm_loadu_pd(A); A += 2;
			a2 = _mm_loadu_pd(A); A += 2;
			a3 = _mm_loadu_pd(A); A += 2;
			a4 = _mm_loadu_pd(A); A += 2;

			a1 = _mm_mul_pd(a1, b);
			a2 = _mm_mul_pd(a2, b);
			a3 = _mm_mul_pd(a3, b);
			a4 = _mm_mul_pd(a4, b);

			_mm_storeu_pd(A - 8, a1);
			_mm_storeu_pd(A - 6, a2);
			_mm_storeu_pd(A - 4, a3);
			_mm_storeu_pd(A - 2, a4);
		}
	}

	for (; i < n; ++i)
		*A++ *= B;
}

TARGETSSE void AddProdSSE(double* C, double* A, double* B, int64 n)
{
	int64 i = 0;

	if (n >= 8)
	{
		__m128d a1, a2, a3, a4;
		__m128d b1, b2, b3, b4;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			a1 = _mm_loadu_pd(A); A += 2;
			a2 = _mm_loadu_pd(A); A += 2;
			a3 = _mm_loadu_pd(A); A += 2;
			a4 = _mm_loadu_pd(A); A += 2;

			b1 = _mm_loadu_pd(B); B += 2;
			b2 = _mm_loadu_pd(B); B += 2;
			b3 = _mm_loadu_pd(B); B += 2;
			b4 = _mm_loadu_pd(B); B += 2;

			a1 = _mm_mul_pd(a1, b1);
			a2 = _mm_mul_pd(a2, b2);
			a3 = _mm_mul_pd(a3, b3);
			a4 = _mm_mul_pd(a4, b4);

			b1 = _mm_loadu_pd(C); C += 2;
			b2 = _mm_loadu_pd(C); C += 2;
			b3 = _mm_loadu_pd(C); C += 2;
			b4 = _mm_loadu_pd(C); C += 2;

			a1 = _mm_add_pd(a1, b1);
			a2 = _mm_add_pd(a2, b2);
			a3 = _mm_add_pd(a3, b3);
			a4 = _mm_add_pd(a4, b4);

			_mm_storeu_pd(C - 8, a1);
			_mm_storeu_pd(C - 6, a2);
			_mm_storeu_pd(C - 4, a3);
			_mm_storeu_pd(C - 2, a4);
		}
	}

	for (; i < n; ++i)
		*C++ += *A++ * *B++;
}

TARGETSSE void AddProdSSE(double* C, double* A, double B, int64 n)
{
	int64 i = 0;

	if (n >= 2)
	{
		__m128d b = _mm_set1_pd(B);

		for (int64 l1 = n - 2; i <= l1; i += 2)
		{
			_mm_storeu_pd(C, _mm_add_pd(_mm_loadu_pd(C), _mm_mul_pd(_mm_loadu_pd(A), b)));
			A += 2; C += 2;
		}
	}

	for (; i < n; ++i)
		*C++ += *A++ * B;
}

TARGETSSE void UnifySSE(double* A, int64 n)
{
	int64 i = 0;
	double invsum = 1.0 / (SumSSE(A, n) + n * MIN_FREQ);

	if (n >= 8)
	{
		__m128d minv = _mm_set1_pd(MIN_FREQ), invs = _mm_set1_pd(invsum);

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			_mm_storeu_pd(A, _mm_mul_pd(invs, _mm_add_pd(_mm_loadu_pd(A), minv))); A += 2;
			_mm_storeu_pd(A, _mm_mul_pd(invs, _mm_add_pd(_mm_loadu_pd(A), minv))); A += 2;
			_mm_storeu_pd(A, _mm_mul_pd(invs, _mm_add_pd(_mm_loadu_pd(A), minv))); A += 2;
			_mm_storeu_pd(A, _mm_mul_pd(invs, _mm_add_pd(_mm_loadu_pd(A), minv))); A += 2;
		}
	}

	for (; i < n; ++i, ++A)
		*A = (*A + MIN_FREQ) * invsum;
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