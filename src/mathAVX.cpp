/* AVX Instruction Set Functions */

#include "vcfpop.h"

#ifndef __aarch64__

#ifndef _RNGAVX
/* Initialize rng */
TARGETAVX RNGAVX::RNGAVX()
{
}

/* Initialize rng */
TARGETAVX RNGAVX::RNGAVX(uint64 s)
{
	__m256i seed[2];

	seed[0] = _mm256_set_epi64x(Hash64ULong(s + 3), Hash64ULong(s + 2), Hash64ULong(s + 1), Hash64ULong(s + 0));
	seed[1] = _mm256_set_epi64x(Hash64ULong(s + 7), Hash64ULong(s + 6), Hash64ULong(s + 5), Hash64ULong(s + 4));
	
	x[0] = _mm256_xor_si256(_mm256_set1_epi64x(0x159A55E5075BCD15), seed[0]);
	x[1] = _mm256_xor_si256(_mm256_set1_epi64x(0x159A55E5075BCD15), seed[1]);

	seed[0] = _mm256_slli_epi64(seed[0], 6);
	seed[1] = _mm256_slli_epi64(seed[1], 6);

	y[0] = _mm256_xor_si256(_mm256_set1_epi64x(0x054913331F123BB5), seed[0]);
	y[1] = _mm256_xor_si256(_mm256_set1_epi64x(0x054913331F123BB5), seed[1]);
}

/* Draw a uniform distriubted interger */
TARGETAVX void RNGAVX::XorShift128p(__m256i* re)
{
	__m256i a[4], b[4];

	a[0] = x[0];
	a[1] = x[1];

	b[0] = y[0];
	b[1] = y[1];

	x[0] = b[0];
	x[1] = b[1];

	a[0] = _mm256_xor_si256(a[0], _mm256_slli_epi64(a[0], 23));
	a[1] = _mm256_xor_si256(a[1], _mm256_slli_epi64(a[1], 23));

	a[0] = _mm256_xor_si256(a[0], _mm256_srli_epi64(a[0], 18));
	a[1] = _mm256_xor_si256(a[1], _mm256_srli_epi64(a[1], 18));

	a[0] = _mm256_xor_si256(a[0], b[0]);
	a[1] = _mm256_xor_si256(a[1], b[1]);

	a[0] = _mm256_xor_si256(a[0], _mm256_srli_epi64(b[0], 5));
	a[1] = _mm256_xor_si256(a[1], _mm256_srli_epi64(b[1], 5));

	y[0] = a[0];
	y[1] = a[1];

	re[0] = _mm256_add_epi64(a[0], b[0]);
	re[1] = _mm256_add_epi64(a[1], b[1]);
}

/* Draw a uniform distriubted real number */
TARGETAVX void RNGAVX::Uniform(__m256d* re)
{
	__m256d one = _mm256_set1_pd(1.0);
	__m256i mask1 = _mm256_set1_epi64x(0x000FFFFFFFFFFFFF);
	__m256i mask2 = _mm256_set1_epi64x(0x3FF0000000000000);
	__m256i* r = (__m256i*)re;

	XorShift128p(r);

	r[0] = _mm256_or_si256(_mm256_and_si256(r[0], mask1), mask2);
	r[1] = _mm256_or_si256(_mm256_and_si256(r[1], mask1), mask2);

	re[0] = _mm256_sub_pd(re[0], one);
	re[1] = _mm256_sub_pd(re[1], one);
}

/* Draw a uniform distriubted real number */
TARGETAVX void RNGAVX::Poly(__m256d* a, __m256d* s, int n, __m256i* re)
{
	__m256d t[2];

	Uniform(t);

	__m256d t1 = _mm256_mul_pd(t[0], s[0]);
	__m256d t2 = _mm256_mul_pd(t[1], s[1]);
	__m256d f1 = _mm256_set1_pd(0);
	__m256d f2 = f1;
	__m256i midx1 = _mm256_set1_epi64x(n - 1);
	__m256i midx2 = midx1;
	__m256i nidx = _mm256_setzero_si256();
	__m256i ninc = _mm256_set1_epi64x(1);
	__m256d b1, b2;
	__m256d v1, v2;

	for (int i = 0; i < n; ++i)
	{
		v1 = a[i * 2 + 0];
		v2 = a[i * 2 + 1];

		b1 = _mm256_cmp_pd(t1, v1, _CMP_LT_OS);
		b2 = _mm256_cmp_pd(t2, v2, _CMP_LT_OS);

		t1 = _mm256_sub_pd(t1, v1);
		t2 = _mm256_sub_pd(t2, v2);

		b1 = _mm256_andnot_pd(f1, b1);
		b2 = _mm256_andnot_pd(f2, b2);

		f1 = _mm256_or_pd(f1, b1);
		f2 = _mm256_or_pd(f2, b2);

		midx1 = _mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(midx1), _mm256_castsi256_pd(nidx), b1));//ok
		midx2 = _mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(midx2), _mm256_castsi256_pd(nidx), b2));

		nidx = _mm256_add_epi64(nidx, ninc);
	}

	re[0] = midx1;
	re[1] = midx2;
}

/* Draw a polynormial distriubted integer with propoirtions in natural logarithm */
TARGETAVX void RNGAVX::PolyLog(__m256d* a, int n, __m256i* re)
{
	//proportional polynomial distribution, will overwrite a
	__m256d maxval[2] = { _mm256_set1_pd(-1e300), _mm256_set1_pd(-1e300) };
	__m256d s[2] = { _mm256_set1_pd(MIN_FREQ * n), _mm256_set1_pd(MIN_FREQ * n) };
	double* af = (double*)a;

	for (int i = 0; i < n; ++i)
	{
		maxval[0] = _mm256_max_pd(maxval[0], a[i * 2 + 0]);
		maxval[1] = _mm256_max_pd(maxval[1], a[i * 2 + 1]);
	}

	for (int i = 0; i < n; ++i)
	{
		a[i * 2 + 0] = _mm256_sub_pd(a[i * 2 + 0], maxval[0]);
		a[i * 2 + 1] = _mm256_sub_pd(a[i * 2 + 1], maxval[1]);

		af[i * 8 + 0] = (af[i * 8 + 0] < -23) ? MIN_FREQ : exp(af[i * 8 + 0]);
		af[i * 8 + 1] = (af[i * 8 + 1] < -23) ? MIN_FREQ : exp(af[i * 8 + 1]);
		af[i * 8 + 2] = (af[i * 8 + 2] < -23) ? MIN_FREQ : exp(af[i * 8 + 2]);
		af[i * 8 + 3] = (af[i * 8 + 3] < -23) ? MIN_FREQ : exp(af[i * 8 + 3]);
		af[i * 8 + 4] = (af[i * 8 + 4] < -23) ? MIN_FREQ : exp(af[i * 8 + 4]);
		af[i * 8 + 5] = (af[i * 8 + 5] < -23) ? MIN_FREQ : exp(af[i * 8 + 5]);
		af[i * 8 + 6] = (af[i * 8 + 6] < -23) ? MIN_FREQ : exp(af[i * 8 + 6]);
		af[i * 8 + 7] = (af[i * 8 + 7] < -23) ? MIN_FREQ : exp(af[i * 8 + 7]);

		s[0] = _mm256_add_pd(s[0], a[i * 2 + 0]);
		s[1] = _mm256_add_pd(s[1], a[i * 2 + 1]);
	}

	__m256d t[2];

	Uniform(t);

	__m256d t1 = _mm256_mul_pd(t[0], s[0]);
	__m256d t2 = _mm256_mul_pd(t[1], s[1]);
	__m256d f1 = _mm256_set1_pd(0);
	__m256d f2 = f1;
	__m256i midx1 = _mm256_set1_epi64x(n - 1);
	__m256i midx2 = midx1;
	__m256i nidx = _mm256_setzero_si256();
	__m256i ninc = _mm256_set1_epi64x(1);
	__m256d b1, b2;
	__m256d v1, v2;

	for (int i = 0; i < n; ++i)
	{
		v1 = a[i * 2 + 0];
		v2 = a[i * 2 + 1];

		b1 = _mm256_cmp_pd(t1, v1, _CMP_LT_OS);
		b2 = _mm256_cmp_pd(t2, v2, _CMP_LT_OS);

		t1 = _mm256_sub_pd(t1, v1);
		t2 = _mm256_sub_pd(t2, v2);

		b1 = _mm256_andnot_pd(f1, b1);
		b2 = _mm256_andnot_pd(f2, b2);

		f1 = _mm256_or_pd(f1, b1);
		f2 = _mm256_or_pd(f2, b2);

		midx1 = _mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(midx1), _mm256_castsi256_pd(nidx), b1));//ok
		midx2 = _mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(midx2), _mm256_castsi256_pd(nidx), b2));

		nidx = _mm256_add_epi64(nidx, ninc);
	}

	re[0] = midx1;
	re[1] = midx2;
}
#endif

TARGETAVX int64 GetMinIdxAVX(double* A, int64 n, double& val)
{
	int64 i = 0;
	val = 1e300;
	uint64 idx = (uint64)-1;

	if (n >= 16)
	{
		__m256d min1 = _mm256_set1_pd(val);
		__m256d min2 = _mm256_set1_pd(val);
		__m256d min3 = _mm256_set1_pd(val);
		__m256d min4 = _mm256_set1_pd(val);
		__m256i midx1 = _mm256_set1_epi8((char)0xFF);
		__m256i midx2 = _mm256_set1_epi8((char)0xFF);
		__m256i midx3 = _mm256_set1_epi8((char)0xFF);
		__m256i midx4 = _mm256_set1_epi8((char)0xFF);
		__m256i nidx1 = _mm256_set_epi64x(3, 2, 1, 0);
		__m256i nidx2 = _mm256_set_epi64x(7, 6, 5, 4);
		__m256i nidx3 = _mm256_set_epi64x(11, 10, 9, 8);
		__m256i nidx4 = _mm256_set_epi64x(15, 14, 13, 12);
		__m256i msep = _mm256_set1_epi64x(16);
		__m256d f1, f2, f3, f4;
		__m256d v1, v2, v3, v4;

		for (int64 l1 = n - 16; i <= l1; i += 16)
		{
			v1 = _mm256_loadu_pd(A); A += 4;
			v2 = _mm256_loadu_pd(A); A += 4;
			v3 = _mm256_loadu_pd(A); A += 4;
			v4 = _mm256_loadu_pd(A); A += 4;

			f1 = _mm256_cmp_pd(min1, v1, _CMP_GT_OS);
			f2 = _mm256_cmp_pd(min2, v2, _CMP_GT_OS);
			f3 = _mm256_cmp_pd(min3, v3, _CMP_GT_OS);
			f4 = _mm256_cmp_pd(min4, v4, _CMP_GT_OS);

			min1 = _mm256_min_pd(min1, v1);
			min2 = _mm256_min_pd(min2, v2);
			min3 = _mm256_min_pd(min3, v3);
			min4 = _mm256_min_pd(min4, v4);

			midx1 = _mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(midx1), _mm256_castsi256_pd(nidx1), f1));//ok
			midx2 = _mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(midx2), _mm256_castsi256_pd(nidx2), f2));
			midx3 = _mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(midx3), _mm256_castsi256_pd(nidx3), f3));
			midx4 = _mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(midx4), _mm256_castsi256_pd(nidx4), f4));

			nidx1 = _mm256_add_epi64(nidx1, msep);
			nidx2 = _mm256_add_epi64(nidx2, msep);
			nidx3 = _mm256_add_epi64(nidx3, msep);
			nidx4 = _mm256_add_epi64(nidx4, msep);
		}

		f1 = _mm256_cmp_pd(min1, min2, _CMP_GT_OS);
		f3 = _mm256_cmp_pd(min3, min4, _CMP_GT_OS);

		min1 = _mm256_min_pd(min1, min2);
		min3 = _mm256_min_pd(min3, min4);

		midx1 = _mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(midx1), _mm256_castsi256_pd(midx2), f1));
		midx3 = _mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(midx3), _mm256_castsi256_pd(midx4), f3));

		f1 = _mm256_cmp_pd(min1, min3, _CMP_GT_OS);
		min1 = _mm256_min_pd(min1, min3);
		midx1 = _mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(midx1), _mm256_castsi256_pd(midx3), f1));

		for (int64 j = 0; j < 4; ++j)
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

TARGETAVX void GetMinMaxValAVX(double* A, int64 n, double& minv, double& maxv)
{
	int64 i = 0;
	minv = 1e300;
	maxv = -1e300;

	if (n >= 16)
	{
		__m256d min1 = _mm256_set1_pd(minv);
		__m256d min2 = _mm256_set1_pd(minv);
		__m256d min3 = _mm256_set1_pd(minv);
		__m256d min4 = _mm256_set1_pd(minv);
		__m256d max1 = _mm256_set1_pd(maxv);
		__m256d max2 = _mm256_set1_pd(maxv);
		__m256d max3 = _mm256_set1_pd(maxv);
		__m256d max4 = _mm256_set1_pd(maxv);
		__m256d v1, v2, v3, v4;

		for (int64 l1 = n - 16; i <= l1; i += 16)
		{
			v1 = _mm256_loadu_pd(A); A += 4;
			v2 = _mm256_loadu_pd(A); A += 4;
			v3 = _mm256_loadu_pd(A); A += 4;
			v4 = _mm256_loadu_pd(A); A += 4;

			min1 = _mm256_min_pd(min1, v1);
			min2 = _mm256_min_pd(min2, v2);
			min3 = _mm256_min_pd(min3, v3);
			min4 = _mm256_min_pd(min4, v4);

			max1 = _mm256_max_pd(max1, v1);
			max2 = _mm256_max_pd(max2, v2);
			max3 = _mm256_max_pd(max3, v3);
			max4 = _mm256_max_pd(max4, v4);
		}

		min1 = _mm256_min_pd(_mm256_min_pd(min1, min2), _mm256_min_pd(min3, min4));
		max1 = _mm256_max_pd(_mm256_max_pd(max1, max2), _mm256_max_pd(max3, max4));

		for (int64 j = 0; j < 4; ++j)
		{
			if (simd_f64(min1, j) < minv)  minv = simd_f64(min1, j);
			if (simd_f64(max1, j) > maxv)  maxv = simd_f64(max1, j);
		}
	}

	for (; i < n; ++i, ++A)
	{
		if (*A < minv) minv = *A;
		if (*A > maxv) maxv = *A;
	}
}

TARGETAVX double GetMaxValAVX(double* A, int64 n)
{
	int64 i = 0;
	double val = -1e300;

	if (n >= 16)
	{
		__m256d max1 = _mm256_set1_pd(val);
		__m256d max2 = _mm256_set1_pd(val);
		__m256d max3 = _mm256_set1_pd(val);
		__m256d max4 = _mm256_set1_pd(val);
		__m256d v1, v2, v3, v4;

		for (int64 l1 = n - 16; i <= l1; i += 16)
		{
			v1 = _mm256_loadu_pd(A); A += 4;
			v2 = _mm256_loadu_pd(A); A += 4;
			v3 = _mm256_loadu_pd(A); A += 4;
			v4 = _mm256_loadu_pd(A); A += 4;

			max1 = _mm256_max_pd(max1, v1);
			max2 = _mm256_max_pd(max2, v2);
			max3 = _mm256_max_pd(max3, v3);
			max4 = _mm256_max_pd(max4, v4);
		}

		max1 = _mm256_max_pd(_mm256_max_pd(max1, max2), _mm256_max_pd(max3, max4));

		for (int64 j = 0; j < 4; ++j)
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

TARGETAVX double GetMaxValAVX(double* A, int64 n, int64 sep)
{
	int64 i = 0;
	double val = -1e300;

	if (n >= 16)
	{
		_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
		_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
		_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
		_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
		_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
		_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
		_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
		_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;

		__m256i vindex = _mm256_set_epi64x(-9 * sep, -10 * sep, -11 * sep, -12 * sep);
		__m256d max1 = _mm256_set1_pd(val), max2 = _mm256_set1_pd(val);

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;

			max1 = _mm256_max_pd(max1, _mm256_i64gather_pd(A, vindex, sizeof(double)));

			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;

			max2 = _mm256_max_pd(max2, _mm256_i64gather_pd(A, vindex, sizeof(double)));
		}

		max1 = _mm256_max_pd(max1, max2);

		for (int64 j = 0; j < 4; ++j)
		{
			if (simd_f64(max1, j) < val) continue;
			val = simd_f64(max1, j);
		}
		
		A -= sep * 8;
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
	int64 i = 0;
	double val = 1e300;

	if (n >= 16)
	{
		__m256d min1 = _mm256_set1_pd(val);
		__m256d min2 = _mm256_set1_pd(val);
		__m256d min3 = _mm256_set1_pd(val);
		__m256d min4 = _mm256_set1_pd(val);
		__m256d v1, v2, v3, v4;

		for (int64 l1 = n - 16; i <= l1; i += 16)
		{
			v1 = _mm256_loadu_pd(A); A += 4;
			v2 = _mm256_loadu_pd(A); A += 4;
			v3 = _mm256_loadu_pd(A); A += 4;
			v4 = _mm256_loadu_pd(A); A += 4;

			min1 = _mm256_min_pd(min1, v1);
			min2 = _mm256_min_pd(min2, v2);
			min3 = _mm256_min_pd(min3, v3);
			min4 = _mm256_min_pd(min4, v4);
		}

		min1 = _mm256_min_pd(_mm256_min_pd(min1, min2), _mm256_min_pd(min3, min4));

		for (int64 j = 0; j < 4; ++j)
		{
			if (simd_f64(min1, j) > val) continue;
			val = simd_f64(min1, j);
		}
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
	int64 i = 0;
	int64 val = 0x7FFFFFFFFFFFFFFF;

	if (n >= 16)
	{
		__m256i min1 = _mm256_set1_epi64x(0x7FFFFFFFFFFFFFFF);
		__m256i min2 = _mm256_set1_epi64x(0x7FFFFFFFFFFFFFFF);
		__m256i min3 = _mm256_set1_epi64x(0x7FFFFFFFFFFFFFFF);
		__m256i min4 = _mm256_set1_epi64x(0x7FFFFFFFFFFFFFFF);
		__m256i v1, v2, v3, v4;
		__m256i f1, f2, f3, f4;

		for (int64 l1 = n - 16; i <= l1; i += 16)
		{
			v1 = _mm256_loadu_si256((__m256i*)A); A += 4;
			v2 = _mm256_loadu_si256((__m256i*)A); A += 4;
			v3 = _mm256_loadu_si256((__m256i*)A); A += 4;
			v4 = _mm256_loadu_si256((__m256i*)A); A += 4;

			f1 = _mm256_cmpgt_epi64(min1, v1);
			f2 = _mm256_cmpgt_epi64(min2, v2);
			f3 = _mm256_cmpgt_epi64(min3, v3);
			f4 = _mm256_cmpgt_epi64(min4, v4);

			min1 = _mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(min1), _mm256_castsi256_pd(v1), _mm256_castsi256_pd(f1)));//ok
			min2 = _mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(min2), _mm256_castsi256_pd(v2), _mm256_castsi256_pd(f2)));
			min3 = _mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(min3), _mm256_castsi256_pd(v3), _mm256_castsi256_pd(f3)));
			min4 = _mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(min4), _mm256_castsi256_pd(v4), _mm256_castsi256_pd(f4)));
		}

		f1 = _mm256_cmpgt_epi64(min1, min2);
		f3 = _mm256_cmpgt_epi64(min3, min4);

		min1 = _mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(min1), _mm256_castsi256_pd(min2), _mm256_castsi256_pd(f1)));
		min3 = _mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(min3), _mm256_castsi256_pd(min4), _mm256_castsi256_pd(f3)));

		f1 = _mm256_cmpgt_epi64(min1, min3);
		min1 = _mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(min1), _mm256_castsi256_pd(min3), _mm256_castsi256_pd(f1)));

		for (int64 j = 0; j < 4; ++j)
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

TARGETAVX void ChargeLogAVX(int64& slog, double& prod, __m256d& val)
{
	AddExponentAVX(slog, val);
	prod = prod * simd_f64(val, 0) * simd_f64(val, 1) * simd_f64(val, 2) * simd_f64(val, 3);

	if (prod < DOUBLE_UNDERFLOW || prod > DOUBLE_OVERFLOW) [[unlikely]]
		AddExponent(slog, prod);
}

TARGETAVX double LogProdAVX(double* A, int64 n)
{
	int64 i = 0;
	int64 slog = 0; double prod = 1;
	
	if (n >= 4)
	{
		__m256d dunder = _mm256_set1_pd(DOUBLE_UNDERFLOW), dover = _mm256_set1_pd(DOUBLE_OVERFLOW);
		__m256d pd = _mm256_set1_pd(1.0);
		__m256i maskff = _mm256_set1_epi64x(0xFFFFFFFFFFFFFFFF);

		for (int64 l1 = n - 4; i <= l1; i += 4)
		{
			pd = _mm256_mul_pd(pd, _mm256_loadu_pd(A)); A += 4;

			if (_mm256_testz_si256(
				_mm256_castpd_si256(_mm256_or_pd(_mm256_cmp_pd(pd, dunder, _CMP_LT_OS), _mm256_cmp_pd(dover, pd, _CMP_LT_OS))),
				maskff)) [[likely]] continue;

			//100% slow if add exponent in each loop
			AddExponentAVX(slog, pd);
		}
		ChargeLogAVX(slog, prod, pd);
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

	if (n >= 4)
	{
		_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
		_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
		_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
		_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;

		__m256d pd = _mm256_set1_pd(1.0);
		__m256d dunder = _mm256_set1_pd(DOUBLE_UNDERFLOW), dover = _mm256_set1_pd(DOUBLE_OVERFLOW);
		__m256i vindex = _mm256_set_epi64x(-5 * sep, -6 * sep, -7 * sep, -8 * sep), mask = _mm256_set1_epi64x(0xFFFFFFFFFFFFFFFF);
		__m256i maskff = _mm256_set1_epi64x(0xFFFFFFFFFFFFFFFF);

		for (int64 l1 = n - 4; i <= l1; i += 4)
		{
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;

			pd = _mm256_mul_pd(pd, _mm256_i64gather_pd(A, vindex, sizeof(double)));

			if (_mm256_testz_si256(
				_mm256_castpd_si256(_mm256_or_pd(_mm256_cmp_pd(pd, dunder, _CMP_LT_OS), _mm256_cmp_pd(dover, pd, _CMP_LT_OS))),
				maskff)) [[likely]] continue;

			AddExponentAVX(slog, pd);
		}

		ChargeLogAVX(slog, prod, pd);

		A -= sep * 4;
	}

	for (; i < n; ++i, A += sep)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);

	return prod;
}

TARGETAVX double LogProdDivAVX(double* A, double* B, int64 n, int64 sep)
{
	int64 i = 0;
	int64 slog = 0; double prod = 1;

	if (n >= 4)
	{
		_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
		_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
		_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
		_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;

		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;

		__m256d pd = _mm256_set1_pd(1.0);
		__m256d dunder = _mm256_set1_pd(DOUBLE_UNDERFLOW), dover = _mm256_set1_pd(DOUBLE_OVERFLOW);
		__m256i vindex = _mm256_set_epi64x(-5 * sep, -6 * sep, -7 * sep, -8 * sep), mask = _mm256_set1_epi64x(0xFFFFFFFFFFFFFFFF);
		__m256i maskff = _mm256_set1_epi64x(0xFFFFFFFFFFFFFFFF);

		for (int64 l1 = n - 4; i <= l1; i += 4)
		{
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;

			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;

			pd = _mm256_mul_pd(pd, _mm256_div_pd(_mm256_i64gather_pd(A, vindex, sizeof(double)), _mm256_i64gather_pd(B, vindex, sizeof(double))));

			if (_mm256_testz_si256(
				_mm256_castpd_si256(_mm256_or_pd(_mm256_cmp_pd(pd, dunder, _CMP_LT_OS), _mm256_cmp_pd(dover, pd, _CMP_LT_OS))),
				maskff)) [[likely]] continue;

			AddExponentAVX(slog, pd);
		}

		ChargeLogAVX(slog, prod, pd);
		A -= sep * 4;
		B -= sep * 4;
	}

	for (; i < n; ++i, A += sep, B += sep)
		ChargeLog(slog, prod, *A / *B);

	CloseLog(slog, prod);

	return prod;
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
	int64 i = 0;
	double re = 0;

	if (n >= 16)
	{
		__m256d s1 = _mm256_setzero_pd();
		__m256d s2 = _mm256_setzero_pd();
		__m256d s3 = _mm256_setzero_pd();
		__m256d s4 = _mm256_setzero_pd();
		__m256d v1, v2, v3, v4;

		for (int64 l1 = n - 16; i <= l1; i += 16)
		{
			v1 = _mm256_loadu_pd(A); A += 4;
			v2 = _mm256_loadu_pd(A); A += 4;
			v3 = _mm256_loadu_pd(A); A += 4;
			v4 = _mm256_loadu_pd(A); A += 4;

			s1 = _mm256_add_pd(s1, v1);
			s2 = _mm256_add_pd(s2, v2);
			s3 = _mm256_add_pd(s3, v3);
			s4 = _mm256_add_pd(s4, v4);
		}

		s1 = _mm256_add_pd(_mm256_add_pd(s1, s2), _mm256_add_pd(s3, s4));
		s1 = _mm256_hadd_pd(s1, s1);
		re = simd_f64(s1, 0) + simd_f64(s1, 2);
	}

	for (; i < n; ++i)
		re += *A++;

	return re;
}

TARGETAVX int64 SumAVX(byte* A, int64 n)
{
	uint64 re = 0;
	int64 i = 0;

	if (n >= 128)
	{
		__m256i s1 = _mm256_setzero_si256();
		__m256i s2 = _mm256_setzero_si256();
		__m256i s3 = _mm256_setzero_si256();
		__m256i s4 = _mm256_setzero_si256();
		__m256i z = _mm256_setzero_si256();
		__m256i v1, v2, v3, v4;

		for (int64 l1 = n - 128; i <= l1; i += 128)
		{
			v1 = _mm256_loadu_si256((__m256i*)A); A += 32;
			v2 = _mm256_loadu_si256((__m256i*)A); A += 32;
			v3 = _mm256_loadu_si256((__m256i*)A); A += 32;
			v4 = _mm256_loadu_si256((__m256i*)A); A += 32;

			s1 = _mm256_add_epi64(s1, _mm256_sad_epu8(v1, z));
			s2 = _mm256_add_epi64(s2, _mm256_sad_epu8(v2, z));
			s3 = _mm256_add_epi64(s3, _mm256_sad_epu8(v3, z));
			s4 = _mm256_add_epi64(s4, _mm256_sad_epu8(v4, z));
		}

		s1 = _mm256_add_epi64(_mm256_add_epi64(s1, s2), _mm256_add_epi64(s3, s4));
		re += simd_u64(s1, 0) + simd_u64(s1, 1) + simd_u64(s1, 2) + simd_u64(s1, 3);
	}

	for (; i < n; ++i)
		re += *A++;

	return re;
}

TARGETAVX double SumAVX(double* A, int64 n, int64 sep)
{
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

		__m256d s1 = _mm256_setzero_pd(), s2 = _mm256_setzero_pd();
		__m256i vindex = _mm256_set_epi64x(-9 * sep, -10 * sep, -11 * sep, -12 * sep);

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			s1 = _mm256_add_pd(s1, _mm256_i64gather_pd(A, vindex, sizeof(double)));

			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			s2 = _mm256_add_pd(s2, _mm256_i64gather_pd(A, vindex, sizeof(double)));
		}

		s1 = _mm256_add_pd(s1, s2);
		s1 = _mm256_hadd_pd(s1, s1);
		re = simd_f64(s1, 0) + simd_f64(s1, 2); 
		A -= sep * 8;
	}

	for (; i < n; ++i, A += sep)
		re += *A;

	return re;
}

TARGETAVX void SumAVX(double* A, double** B, int64 k, int64 n)
{
	int64 i = 0;

	if (n >= 8)
	{
		__m256d a1, a2;
		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			a1 = _mm256_setzero_pd();
			a2 = _mm256_setzero_pd();

			for (int64 j = 0; j < k; ++j)
			{
				a1 = _mm256_add_pd(a1, _mm256_loadu_pd(&B[j][i + 0]));
				a2 = _mm256_add_pd(a2, _mm256_loadu_pd(&B[j][i + 4]));
			}

			_mm256_storeu_pd(&A[i + 0], a1);
			_mm256_storeu_pd(&A[i + 4], a2);
		}
	}

	for (; i < n; ++i)
	{
		A[i] = 0;
		for (int64 j = 0; j < k; ++j)
			A[i] += B[j][i];
	}
}

TARGETAVX double ProdAVX(double* A, int64 n)
{
	int64 i = 0;
	double re = 1;

	if (n >= 16)
	{
		__m256d pd1 = _mm256_set1_pd(1.0);
		__m256d pd2 = _mm256_set1_pd(1.0);
		__m256d pd3 = _mm256_set1_pd(1.0);
		__m256d pd4 = _mm256_set1_pd(1.0);
		__m256d v1, v2, v3, v4;

		for (int64 l1 = n - 16; i <= l1; i += 16)
		{
			v1 = _mm256_loadu_pd(A); A += 4;
			v2 = _mm256_loadu_pd(A); A += 4;
			v3 = _mm256_loadu_pd(A); A += 4;
			v4 = _mm256_loadu_pd(A); A += 4;

			pd1 = _mm256_mul_pd(pd1, v1);
			pd2 = _mm256_mul_pd(pd2, v2);
			pd3 = _mm256_mul_pd(pd3, v3);
			pd4 = _mm256_mul_pd(pd4, v4);
		}

		pd1 = _mm256_mul_pd(_mm256_mul_pd(pd1, pd2), _mm256_mul_pd(pd3, pd4));
		re = simd_f64(pd1, 0) * simd_f64(pd1, 1) * simd_f64(pd1, 2) * simd_f64(pd1, 3);
	}

	for (; i < n; ++i)
		re *= *A++;

	return re;
}

TARGETAVX double ProdAVX(double* A, int64 n, int64 sep)
{
	int64 i = 0;
	double re = 1;

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

		__m256d pd1 = _mm256_set1_pd(1.0), pd2 = _mm256_set1_pd(1.0);
		__m256i vindex = _mm256_set_epi64x(-9 * sep, -10 * sep, -11 * sep, -12 * sep);

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			pd1 = _mm256_mul_pd(pd1, _mm256_i64gather_pd(A, vindex, sizeof(double)));

			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			_mm_prefetch((const char*)A, _MM_HINT_T0); A += sep;
			pd2 = _mm256_mul_pd(pd2, _mm256_i64gather_pd(A, vindex, sizeof(double)));
		}

		pd1 = _mm256_mul_pd(pd1, pd2);
		re = simd_f64(pd1, 0) * simd_f64(pd1, 1) * simd_f64(pd1, 2) * simd_f64(pd1, 3);
		A -= sep * 8;
	}

	for (; i < n; ++i, A += sep)
		re *= *A;

	return re;
}

TARGETAVX double SumSquareAVX(double* A, int64 n)
{
	int64 i = 0;
	double re = 0;

	if (n >= 16)
	{
		__m256d s1 = _mm256_setzero_pd();
		__m256d s2 = _mm256_setzero_pd();
		__m256d s3 = _mm256_setzero_pd();
		__m256d s4 = _mm256_setzero_pd();
		__m256d a1, a2, a3, a4;

		for (int64 l1 = n - 16; i <= l1; i += 16)
		{
			a1 = _mm256_loadu_pd(A); A += 4;
			a2 = _mm256_loadu_pd(A); A += 4;
			a3 = _mm256_loadu_pd(A); A += 4;
			a4 = _mm256_loadu_pd(A); A += 4;

			s1 = _mm256_fmadd_pd(a1, a1, s1);
			s2 = _mm256_fmadd_pd(a2, a2, s2);
			s3 = _mm256_fmadd_pd(a3, a3, s3);
			s4 = _mm256_fmadd_pd(a4, a4, s4);
		}

		s1 = _mm256_add_pd(_mm256_add_pd(s1, s2), _mm256_add_pd(s3, s4));
		s1 = _mm256_hadd_pd(s1, s1);
		re = simd_f64(s1, 0) + simd_f64(s1, 2);
	}

	for (; i < n; ++i, ++A)
		re += *A * *A;

	return re;
}

TARGETAVX int64 SumSquareAVX(byte* A, int64 n)
{
	int64 i = 0;
	uint64 re = 0;

	if (n >= 128)
	{
		__m256i a1, a2, a3, a4;
		__m256i s1 = _mm256_setzero_si256();
		__m256i s2 = _mm256_setzero_si256();
		__m256i s3 = _mm256_setzero_si256();
		__m256i s4 = _mm256_setzero_si256();

		for (int64 l1 = n - 128; i <= l1; i += 128)
		{
			a1 = _mm256_loadu_si256((__m256i*)A); A += 32;
			a2 = _mm256_loadu_si256((__m256i*)A); A += 32;
			a3 = _mm256_loadu_si256((__m256i*)A); A += 32;
			a4 = _mm256_loadu_si256((__m256i*)A); A += 32;

			a1 = _mm256_maddubs_epi16(a1, a1);
			a2 = _mm256_maddubs_epi16(a2, a2);
			a3 = _mm256_maddubs_epi16(a3, a3);
			a4 = _mm256_maddubs_epi16(a4, a4);

			s1 = _mm256_add_epi32(s1, _mm256_add_epi32(_mm256_cvtepu16_epi32(_mm256_castsi256_si128(a1)), _mm256_cvtepu16_epi32(_mm256_extractf128_si256(a1, 1))));
			s2 = _mm256_add_epi32(s2, _mm256_add_epi32(_mm256_cvtepu16_epi32(_mm256_castsi256_si128(a2)), _mm256_cvtepu16_epi32(_mm256_extractf128_si256(a2, 1))));
			s3 = _mm256_add_epi32(s3, _mm256_add_epi32(_mm256_cvtepu16_epi32(_mm256_castsi256_si128(a3)), _mm256_cvtepu16_epi32(_mm256_extractf128_si256(a3, 1))));
			s4 = _mm256_add_epi32(s4, _mm256_add_epi32(_mm256_cvtepu16_epi32(_mm256_castsi256_si128(a4)), _mm256_cvtepu16_epi32(_mm256_extractf128_si256(a4, 1))));
		}

		s1 = _mm256_add_epi32(_mm256_add_epi32(s1, s2), _mm256_add_epi32(s3, s4));
		s1 = _mm256_hadd_epi32(s1, s1);
		s1 = _mm256_hadd_epi32(s1, s1);
		re += simd_u32(s1, 0) + simd_u32(s1, 4);
	}

	for (; i < n; ++i, ++A)
		re += *A * *A;

	return re;
}

TARGETAVX void SumSumSquareAVX(double* A, int64 n, double& sum, double& sumsq)
{
	int64 i = 0;
	sum = sumsq = 0;

	if (n >= 8)
	{
		__m256d s1 = _mm256_setzero_pd(), sq1 = _mm256_setzero_pd();
		__m256d s2 = _mm256_setzero_pd(), sq2 = _mm256_setzero_pd();
		__m256d a1, a2;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			a1 = _mm256_loadu_pd(A); A += 4;
			s1 = _mm256_add_pd(s1, a1);
			sq1 = _mm256_fmadd_pd(a1, a1, sq1);

			a2 = _mm256_loadu_pd(A); A += 4;
			s2 = _mm256_add_pd(s2, a2);
			sq2 = _mm256_fmadd_pd(a2, a2, sq2);
		}

		s1 = _mm256_add_pd(s1, s2);
		sq1 = _mm256_add_pd(sq1, sq2);

		s1 = _mm256_hadd_pd(s1, s1);
		sum = simd_f64(s1, 0) + simd_f64(s1, 2);
		sq1 = _mm256_hadd_pd(sq1, sq1);
		sumsq = simd_f64(sq1, 0) + simd_f64(sq1, 2);
	}

	for (; i < n; ++i, ++A)
	{
		sum += *A;
		sumsq += *A * *A;
	}
}

TARGETAVX double SumProdDivAVX(double* A1, double* A2, double* B, int64 sep, int64 n)
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

		__m256d s1 = _mm256_setzero_pd(), s2 = _mm256_setzero_pd(), b;
		__m256i vindex = _mm256_set_epi64x(-9 * sep, -10 * sep, -11 * sep, -12 * sep);

		for (int64 l1 = n - 4; i <= l1; i += 4)
		{
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;

			b = _mm256_i64gather_pd(B, vindex, sizeof(double));

			s1 = _mm256_fmadd_pd(_mm256_loadu_pd(A1), b, s1); A1 += 4;
			s2 = _mm256_fmadd_pd(_mm256_loadu_pd(A2), b, s2); A2 += 4;
		}

		s1 = _mm256_hadd_pd(s1, s1);  re1 = simd_f64(s1, 0) + simd_f64(s1, 2);
		s2 = _mm256_hadd_pd(s2, s2);  re2 = simd_f64(s2, 0) + simd_f64(s2, 2);
		B -= sep * 8;
	}

	for (; i < n; ++i, A1++, A2++, B += sep)
	{
		re1 += *A1 * *B;
		re2 += *A2 * *B;
	}

	return re1 / re2;
}

TARGETAVX double SumProdAVX(double* A, double* B, int64 sep, int64 n)
{
	int64 i = 0;
	double re = 0;

	if (n >= 4)
	{
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;

		__m256d s = _mm256_setzero_pd();
		__m256i vindex = _mm256_set_epi64x(-5 * sep, -6 * sep, -7 * sep, -8 * sep);

		for (int64 l1 = n - 4; i <= l1; i += 4)
		{
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;

			s = _mm256_fmadd_pd(_mm256_loadu_pd(A), _mm256_i64gather_pd(B, vindex, sizeof(double)), s);
			A += 4;
		}

		s = _mm256_hadd_pd(s, s);
		re = simd_f64(s, 0) + simd_f64(s, 2);
		B -= sep * 4;
	}

	for (; i < n; ++i, A++, B += sep)
		re += *A * *B;

	return re;
}

TARGETAVX double SumProdAVX(double* A, double* B, int64 n)
{
	int64 i = 0;
	double re = 0;

	if (n >= 16)
	{
		__m256d s1 = _mm256_setzero_pd();
		__m256d s2 = _mm256_setzero_pd();
		__m256d s3 = _mm256_setzero_pd();
		__m256d s4 = _mm256_setzero_pd();
		__m256d a1, a2, a3, a4;
		__m256d b1, b2, b3, b4;

		for (int64 l1 = n - 16; i <= l1; i += 16)
		{
			a1 = _mm256_loadu_pd(A); A += 4;
			a2 = _mm256_loadu_pd(A); A += 4;
			a3 = _mm256_loadu_pd(A); A += 4;
			a4 = _mm256_loadu_pd(A); A += 4;

			b1 = _mm256_loadu_pd(B); B += 4;
			b2 = _mm256_loadu_pd(B); B += 4;
			b3 = _mm256_loadu_pd(B); B += 4;
			b4 = _mm256_loadu_pd(B); B += 4;

			s1 = _mm256_fmadd_pd(a1, b1, s1);
			s2 = _mm256_fmadd_pd(a2, b2, s2);
			s3 = _mm256_fmadd_pd(a3, b3, s3);
			s4 = _mm256_fmadd_pd(a4, b4, s4);
		}

		s1 = _mm256_add_pd(_mm256_add_pd(s1, s2), _mm256_add_pd(s3, s4));
		s1 = _mm256_hadd_pd(s1, s1);
		re = simd_f64(s1, 0) + simd_f64(s1, 2);
	}

	for (; i < n; ++i, ++A, ++B)
		re += *A * *B;

	return re;
}

TARGETAVX void AddAVX(double* A, double* B, int64 n)
{
	int64 i = 0;

	if (n >= 16)
	{
		__m256d a1, a2, a3, a4;
		__m256d b1, b2, b3, b4;

		for (int64 l1 = n - 16; i <= l1; i += 16)
		{
			a1 = _mm256_loadu_pd(A); A += 4;
			a2 = _mm256_loadu_pd(A); A += 4;
			a3 = _mm256_loadu_pd(A); A += 4;
			a4 = _mm256_loadu_pd(A); A += 4;

			b1 = _mm256_loadu_pd(B); B += 4;
			b2 = _mm256_loadu_pd(B); B += 4;
			b3 = _mm256_loadu_pd(B); B += 4;
			b4 = _mm256_loadu_pd(B); B += 4;

			_mm256_storeu_pd(A - 16, _mm256_add_pd(a1, b1));
			_mm256_storeu_pd(A - 12, _mm256_add_pd(a2, b2));
			_mm256_storeu_pd(A - 8, _mm256_add_pd(a3, b3));
			_mm256_storeu_pd(A - 4, _mm256_add_pd(a4, b4));
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGETAVX void AddAVX(int64* A, int64* B, int64 n)
{
	int64 i = 0;

	if (n >= 16)
	{
		__m256i a1, a2, a3, a4;
		__m256i b1, b2, b3, b4;

		for (int64 l1 = n - 16; i <= l1; i += 16)
		{
			a1 = _mm256_loadu_si256((__m256i*)A); A += 4;
			a2 = _mm256_loadu_si256((__m256i*)A); A += 4;
			a3 = _mm256_loadu_si256((__m256i*)A); A += 4;
			a4 = _mm256_loadu_si256((__m256i*)A); A += 4;

			b1 = _mm256_loadu_si256((__m256i*)B); B += 4;
			b2 = _mm256_loadu_si256((__m256i*)B); B += 4;
			b3 = _mm256_loadu_si256((__m256i*)B); B += 4;
			b4 = _mm256_loadu_si256((__m256i*)B); B += 4;

			a1 = _mm256_add_epi64(a1, b1);
			a2 = _mm256_add_epi64(a2, b2);
			a3 = _mm256_add_epi64(a3, b3);
			a4 = _mm256_add_epi64(a4, b4);

			_mm256_storeu_si256((__m256i*)(A - 16), a1);
			_mm256_storeu_si256((__m256i*)(A - 12), a2);
			_mm256_storeu_si256((__m256i*)(A - 8), a3);
			_mm256_storeu_si256((__m256i*)(A - 4), a4);
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGETAVX void AddAVX(int* A, int* B, int64 n)
{
	int64 i = 0;

	if (n >= 32)
	{
		__m256i a1, a2, a3, a4;
		__m256i b1, b2, b3, b4;

		for (int64 l1 = n - 32; i <= l1; i += 32)
		{
			a1 = _mm256_loadu_si256((__m256i*)A); A += 8;
			a2 = _mm256_loadu_si256((__m256i*)A); A += 8;
			a3 = _mm256_loadu_si256((__m256i*)A); A += 8;
			a4 = _mm256_loadu_si256((__m256i*)A); A += 8;

			b1 = _mm256_loadu_si256((__m256i*)B); B += 8;
			b2 = _mm256_loadu_si256((__m256i*)B); B += 8;
			b3 = _mm256_loadu_si256((__m256i*)B); B += 8;
			b4 = _mm256_loadu_si256((__m256i*)B); B += 8;

			a1 = _mm256_add_epi32(a1, b1);
			a2 = _mm256_add_epi32(a2, b2);
			a3 = _mm256_add_epi32(a3, b3);
			a4 = _mm256_add_epi32(a4, b4);

			_mm256_storeu_si256((__m256i*)(A - 32), a1);
			_mm256_storeu_si256((__m256i*)(A - 24), a2);
			_mm256_storeu_si256((__m256i*)(A - 16), a3);
			_mm256_storeu_si256((__m256i*)(A - 8), a4);
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGETAVX void AddAVX(double* A, double B, int64 n)
{
	int64 i = 0;

	if (n >= 16)
	{
		__m256d b = _mm256_set1_pd(B);
		__m256d a1, a2, a3, a4;

		for (int64 l1 = n - 16; i <= l1; i += 16)
		{
			a1 = _mm256_loadu_pd(A); A += 4;
			a2 = _mm256_loadu_pd(A); A += 4;
			a3 = _mm256_loadu_pd(A); A += 4;
			a4 = _mm256_loadu_pd(A); A += 4;

			a1 = _mm256_add_pd(a1, b);
			a2 = _mm256_add_pd(a2, b);
			a3 = _mm256_add_pd(a3, b);
			a4 = _mm256_add_pd(a4, b);

			_mm256_storeu_pd(A - 16, a1);
			_mm256_storeu_pd(A - 12, a2);
			_mm256_storeu_pd(A - 8, a3);
			_mm256_storeu_pd(A - 4, a4);
		}
	}

	for (; i < n; ++i, A++)
		*A += B;
}

TARGETAVX void MulAVX(double* C, double* A, double* B, int64 n)
{
	int64 i = 0;

	if (n >= 16)
	{
		__m256d a1, a2, a3, a4;
		__m256d b1, b2, b3, b4;

		for (int64 l1 = n - 16; i <= l1; i += 16)
		{
			a1 = _mm256_loadu_pd(A); A += 4;
			a2 = _mm256_loadu_pd(A); A += 4;
			a3 = _mm256_loadu_pd(A); A += 4;
			a4 = _mm256_loadu_pd(A); A += 4;

			b1 = _mm256_loadu_pd(B); B += 4;
			b2 = _mm256_loadu_pd(B); B += 4;
			b3 = _mm256_loadu_pd(B); B += 4;
			b4 = _mm256_loadu_pd(B); B += 4;

			a1 = _mm256_mul_pd(a1, b1);
			a2 = _mm256_mul_pd(a2, b2);
			a3 = _mm256_mul_pd(a3, b3);
			a4 = _mm256_mul_pd(a4, b4);

			_mm256_storeu_pd(C, a1); C += 4;
			_mm256_storeu_pd(C, a2); C += 4;
			_mm256_storeu_pd(C, a3); C += 4;
			_mm256_storeu_pd(C, a4); C += 4;
		}
	}

	for (; i < n; ++i)
		*C++ = *A++ * *B++;
}

TARGETAVX void MulAVX(double* C, double* A, double B, int64 n)
{
	int64 i = 0;

	if (n >= 16)
	{
		__m256d b = _mm256_set1_pd(B);
		__m256d a1, a2, a3, a4;

		for (int64 l1 = n - 16; i <= l1; i += 16)
		{
			a1 = _mm256_loadu_pd(A); A += 4;
			a2 = _mm256_loadu_pd(A); A += 4;
			a3 = _mm256_loadu_pd(A); A += 4;
			a4 = _mm256_loadu_pd(A); A += 4;

			a1 = _mm256_mul_pd(a1, b);
			a2 = _mm256_mul_pd(a2, b);
			a3 = _mm256_mul_pd(a3, b);
			a4 = _mm256_mul_pd(a4, b);

			_mm256_storeu_pd(C, a1); C += 4;
			_mm256_storeu_pd(C, a2); C += 4;
			_mm256_storeu_pd(C, a3); C += 4;
			_mm256_storeu_pd(C, a4); C += 4;
		}
	}

	for (; i < n; ++i)
		*C++ = *A++ * B;
}

TARGETAVX void MulAVX(double* A, double B, int64 n)
{
	int64 i = 0;

	if (n >= 16)
	{
		__m256d b = _mm256_set1_pd(B);
		__m256d a1, a2, a3, a4;

		for (int64 l1 = n - 16; i <= l1; i += 16)
		{
			a1 = _mm256_loadu_pd(A); A += 4;
			a2 = _mm256_loadu_pd(A); A += 4;
			a3 = _mm256_loadu_pd(A); A += 4;
			a4 = _mm256_loadu_pd(A); A += 4;

			a1 = _mm256_mul_pd(a1, b);
			a2 = _mm256_mul_pd(a2, b);
			a3 = _mm256_mul_pd(a3, b);
			a4 = _mm256_mul_pd(a4, b);

			_mm256_storeu_pd(A - 16, a1);
			_mm256_storeu_pd(A - 12, a2);
			_mm256_storeu_pd(A - 8, a3);
			_mm256_storeu_pd(A - 4, a4);
		}
	}

	for (; i < n; ++i)
		*A++ *= B;
}

TARGETAVX void AddProdAVX(double* C, double* A, double* B, int64 n)
{
	int64 i = 0;

	if (n >= 4)
	{
		for (int64 l1 = n - 4; i <= l1; i += 4)
		{
			_mm256_storeu_pd(C, _mm256_fmadd_pd(_mm256_loadu_pd(A), _mm256_loadu_pd(B), _mm256_loadu_pd(C)));
			A += 4; B += 4; C += 4;
		}
	}

	for (; i < n; ++i)
		*C++ += *A++ * *B++;
}

TARGETAVX void AddProdAVX(double* C, double* A, double B, int64 n)
{
	int64 i = 0;

	if (n >= 4)
	{
		__m256d b = _mm256_set1_pd(B);

		for (int64 l1 = n - 4; i <= l1; i += 4)
		{
			_mm256_storeu_pd(C, _mm256_fmadd_pd(_mm256_loadu_pd(A), b, _mm256_loadu_pd(C)));
			A += 4; C += 4;
		}
	}

	for (; i < n; ++i)
		*C++ += *A++ * B;
}

TARGETAVX void UnifyAVX(double* A, int64 n)
{
	int64 i = 0;
	double invsum = 1.0 / (SumAVX(A, n) + n * MIN_FREQ);

	if (n >= 8)
	{
		__m256d minv = _mm256_set1_pd(invsum * MIN_FREQ), invs = _mm256_set1_pd(invsum);

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			_mm256_storeu_pd(A, _mm256_fmadd_pd(_mm256_loadu_pd(A), invs, minv));
			A += 4;

			_mm256_storeu_pd(A, _mm256_fmadd_pd(_mm256_loadu_pd(A), invs, minv));
			A += 4;
		}
	}

	for (; i < n; ++i, ++A)
		*A = (*A + MIN_FREQ) * invsum;
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