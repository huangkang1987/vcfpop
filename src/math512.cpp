/* AVX512 Instruction Set Functions */
#include "vcfpop.h"

#ifndef __aarch64__

#ifndef _RNG512
/* Initialize rng */
TARGET512 RNG512::RNG512()
{

}

/* Initialize rng */
TARGET512 RNG512::RNG512(uint64 s)
{
	__m512i seed = _mm512_set_epi64(
		Hash64ULong(s + 7), Hash64ULong(s + 6), Hash64ULong(s + 5), Hash64ULong(s + 4),
		Hash64ULong(s + 3), Hash64ULong(s + 2), Hash64ULong(s + 1), Hash64ULong(s + 0));

	x = _mm512_xor_si512(_mm512_set1_epi64(0x159A55E5075BCD15), seed);

	seed = _mm512_slli_epi64(seed, 6);

	y = _mm512_xor_si512(_mm512_set1_epi64(0x054913331F123BB5), seed);
}

/* Draw a uniform distriubted interger */
TARGET512 void RNG512::XorShift128p(__m512i& re)
{
	__m512i a, b;

	a = x;
	b = y;
	x = b;

	a = _mm512_xor_si512(a, _mm512_slli_epi64(a, 23));
	a = _mm512_xor_si512(a, _mm512_srli_epi64(a, 18));
	a = _mm512_xor_si512(a, b);
	a = _mm512_xor_si512(a, _mm512_srli_epi64(b, 5));
	y = a;
	re = _mm512_add_epi64(a, b);
}

/* Draw a uniform distriubted real number */
TARGET512 void RNG512::Uniform(__m512d& re)
{
	__m512d one = _mm512_set1_pd(1.0);
	__m512i mask1 = _mm512_set1_epi64(0x000FFFFFFFFFFFFF);
	__m512i mask2 = _mm512_set1_epi64(0x3FF0000000000000);
	__m512i& r = *(__m512i*)&re;

	XorShift128p(r);

	r = _mm512_or_si512(_mm512_and_si512(r, mask1), mask2);

	re = _mm512_sub_pd(re, one);
}

/* Draw a uniform distriubted real number */
TARGET512 void RNG512::Poly(__m512d* a, __m512d& s, int n, __m512i& re)
{
	__m512d t;

	Uniform(t);

	t = _mm512_mul_pd(t, s);
	__mmask8 f = 0;
	__m512i midx = _mm512_set1_epi64(n - 1);
	__m512i nidx = _mm512_setzero_si512();
	__m512i ninc = _mm512_set1_epi64(1);
	__mmask8 b;
	__m512d v;

	for (int i = 0; i < n; ++i)
	{
		v = a[i];

		b = _mm512_cmp_pd_mask(t, v, _CMP_LT_OS);

		t = _mm512_sub_pd(t, v);

		b = (~f) & b;

		f = f | b;

		midx = _mm512_mask_blend_epi64(b, midx, nidx);//ok

		nidx = _mm512_add_epi64(nidx, ninc);
	}

	re = midx;
}

/* Draw a polynormial distriubted integer with propoirtions in natural logarithm */
TARGET512 void RNG512::PolyLog(__m512d* a, int n, __m512i& re)
{
	//proportional polynomial distribution, will overwrite a
	__m512d maxval = _mm512_set1_pd(-1e300);
	__m512d s = _mm512_set1_pd(MIN_FREQ * n);
	__m512d minfreq = _mm512_set1_pd(MIN_FREQ);
	double* af = (double*)a;

	for (int i = 0; i < n; ++i)
		maxval = _mm512_max_pd(maxval, a[i]);

	for (int i = 0; i < n; ++i)
	{
		a[i] = _mm512_sub_pd(a[i], maxval);

		af[i * 8 + 0] = (af[i * 8 + 0] < -23) ? MIN_FREQ : exp(af[i * 8 + 0]);
		af[i * 8 + 1] = (af[i * 8 + 1] < -23) ? MIN_FREQ : exp(af[i * 8 + 1]);
		af[i * 8 + 2] = (af[i * 8 + 2] < -23) ? MIN_FREQ : exp(af[i * 8 + 2]);
		af[i * 8 + 3] = (af[i * 8 + 3] < -23) ? MIN_FREQ : exp(af[i * 8 + 3]);
		af[i * 8 + 4] = (af[i * 8 + 4] < -23) ? MIN_FREQ : exp(af[i * 8 + 4]);
		af[i * 8 + 5] = (af[i * 8 + 5] < -23) ? MIN_FREQ : exp(af[i * 8 + 5]);
		af[i * 8 + 6] = (af[i * 8 + 6] < -23) ? MIN_FREQ : exp(af[i * 8 + 6]);
		af[i * 8 + 7] = (af[i * 8 + 7] < -23) ? MIN_FREQ : exp(af[i * 8 + 7]);

		s = _mm512_add_pd(s, a[i]);
	}

	__m512d t;
	Uniform(t);

	t = _mm512_mul_pd(t, s);
	__mmask8 f = 0;
	__m512i midx = _mm512_set1_epi64(n - 1);
	__m512i nidx = _mm512_setzero_si512();
	__m512i ninc = _mm512_set1_epi64(1);
	__mmask8 b;
	__m512d v;

	for (int i = 0; i < n; ++i)
	{
		v = a[i];

		b = _mm512_cmp_pd_mask(t, v, _CMP_LT_OS);

		t = _mm512_sub_pd(t, v);

		b = (~f) & b;

		f = f | b;

		midx = _mm512_mask_blend_epi64(b, midx, nidx);//ok

		nidx = _mm512_add_epi64(nidx, ninc);
	}

	re = midx;
}
#endif

TARGET512 int64 GetMinIdx512(double* A, int64 n, double& val)
{
	int64 i = 0;
	val = 1e300;
	uint64 idx = (uint64)-1;

	if (n >= 16)
	{
		__m512d min1 = _mm512_set1_pd(val), min2 = _mm512_set1_pd(val);
		__m512i midx1 = _mm512_set1_epi8((char)0xFF);
		__m512i midx2 = _mm512_set1_epi8((char)0xFF);
		__m512i nidx1 = _mm512_set_epi64(7, 6, 5, 4, 3, 2, 1, 0);
		__m512i nidx2 = _mm512_set_epi64(15, 14, 13, 12, 11, 10, 9, 8);
		__m512i msep = _mm512_set1_epi64(16);
		__m512d ma1, ma2;

		for (int64 l1 = n - 16; i <= l1; i += 16)
		{
			ma1 = _mm512_loadu_pd(A); A += 8;
			ma2 = _mm512_loadu_pd(A); A += 8;

			midx1 = _mm512_mask_and_epi64(midx1, _mm512_cmp_pd_mask(min1, ma1, _CMP_GT_OS), nidx1, nidx1);
			midx2 = _mm512_mask_and_epi64(midx2, _mm512_cmp_pd_mask(min2, ma2, _CMP_GT_OS), nidx2, nidx2);

			min1 = _mm512_min_pd(min1, ma1);
			min2 = _mm512_min_pd(min2, ma2);

			nidx1 = _mm512_add_epi64(nidx1, msep);
			nidx2 = _mm512_add_epi64(nidx2, msep);
		}

		midx1 = _mm512_mask_and_epi64(midx1, _mm512_cmp_pd_mask(min1, min2, _CMP_GT_OS), midx2, midx2);
		min1 = _mm512_min_pd(min1, min2);

		for (int64 j = 0; j < 8; ++j)
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

TARGET512 void GetMinMaxVal512(double* A, int64 n, double& minv, double& maxv)
{
	int64 i = 0;
	minv = 1e300;
	maxv = -1e300;

	if (n >= 32)
	{
		__m512d min1 = _mm512_set1_pd(minv);
		__m512d min2 = _mm512_set1_pd(minv);
		__m512d min3 = _mm512_set1_pd(minv);
		__m512d min4 = _mm512_set1_pd(minv);
		__m512d max1 = _mm512_set1_pd(maxv);
		__m512d max2 = _mm512_set1_pd(maxv);
		__m512d max3 = _mm512_set1_pd(maxv);
		__m512d max4 = _mm512_set1_pd(maxv);
		__m512d v1, v2, v3, v4;

		for (int64 l1 = n - 32; i <= l1; i += 32)
		{
			v1 = _mm512_loadu_pd(A); A += 8;
			v2 = _mm512_loadu_pd(A); A += 8;
			v3 = _mm512_loadu_pd(A); A += 8;
			v4 = _mm512_loadu_pd(A); A += 8;

			min1 = _mm512_min_pd(min1, v1);
			min2 = _mm512_min_pd(min2, v2);
			min3 = _mm512_min_pd(min3, v3);
			min4 = _mm512_min_pd(min4, v4);

			max1 = _mm512_max_pd(max1, v1);
			max2 = _mm512_max_pd(max2, v2);
			max3 = _mm512_max_pd(max3, v3);
			max4 = _mm512_max_pd(max4, v4);
		}

		min1 = _mm512_min_pd(_mm512_min_pd(min1, min2), _mm512_min_pd(min3, min4));
		max1 = _mm512_max_pd(_mm512_max_pd(max1, max2), _mm512_max_pd(max3, max4));

		for (int64 j = 0; j < 8; ++j)
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

TARGET512 double GetMaxVal512(double* A, int64 n)
{
	int64 i = 0;
	double val = -1e300;

	if (n >= 32)
	{
		__m512d max1 = _mm512_set1_pd(val);
		__m512d max2 = _mm512_set1_pd(val);
		__m512d max3 = _mm512_set1_pd(val);
		__m512d max4 = _mm512_set1_pd(val);
		__m512d v1, v2, v3, v4;

		for (int64 l1 = n - 32; i <= l1; i += 32)
		{
			v1 = _mm512_loadu_pd(A); A += 8;
			v2 = _mm512_loadu_pd(A); A += 8;
			v3 = _mm512_loadu_pd(A); A += 8;
			v4 = _mm512_loadu_pd(A); A += 8;

			max1 = _mm512_max_pd(max1, v1);
			max2 = _mm512_max_pd(max2, v2);
			max3 = _mm512_max_pd(max3, v3);
			max4 = _mm512_max_pd(max4, v4);
		}

		max1 = _mm512_max_pd(_mm512_max_pd(max1, max2), _mm512_max_pd(max3, max4));

		for (int64 j = 0; j < 8; ++j)
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

TARGET512 double GetMaxVal512(double* A, int64 n, int64 sep)
{
	int64 i = 0;
	double val = -1e300;

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

		__m512i vindex = _mm512_set_epi64(-9 * sep, -10 * sep, -11 * sep, -12 * sep, -13 * sep, -14 * sep, -15 * sep, -16 * sep);
		__m512d max1 = _mm512_set1_pd(val);

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

			max1 = _mm512_max_pd(max1, _mm512_i64gather_pd(vindex, A, sizeof(double)));
		}

		for (int64 j = 0; j < 8; ++j)
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

TARGET512 double GetMinVal512(double* A, int64 n)
{
	int64 i = 0;
	double val = 1e300;

	if (n >= 32)
	{
		__m512d min1 = _mm512_set1_pd(val);
		__m512d min2 = _mm512_set1_pd(val);
		__m512d min3 = _mm512_set1_pd(val);
		__m512d min4 = _mm512_set1_pd(val);
		__m512d v1, v2, v3, v4;

		for (int64 l1 = n - 32; i <= l1; i += 32)
		{
			v1 = _mm512_loadu_pd(A); A += 8;
			v2 = _mm512_loadu_pd(A); A += 8;
			v3 = _mm512_loadu_pd(A); A += 8;
			v4 = _mm512_loadu_pd(A); A += 8;

			min1 = _mm512_min_pd(min1, v1);
			min2 = _mm512_min_pd(min2, v2);
			min3 = _mm512_min_pd(min3, v3);
			min4 = _mm512_min_pd(min4, v4);
		}

		min1 = _mm512_min_pd(_mm512_min_pd(min1, min2), _mm512_min_pd(min3, min4));

		for (int64 j = 0; j < 8; ++j)
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

TARGET512 int64 GetMinVal512(int64* A, int64 n)
{
	int64 i = 0;
	int64 val = 0x7FFFFFFFFFFFFFFF;

	if (n >= 32)
	{
		__m512i min1 = _mm512_set1_epi64(0x7FFFFFFFFFFFFFFF);
		__m512i min2 = _mm512_set1_epi64(0x7FFFFFFFFFFFFFFF);
		__m512i min3 = _mm512_set1_epi64(0x7FFFFFFFFFFFFFFF);
		__m512i min4 = _mm512_set1_epi64(0x7FFFFFFFFFFFFFFF);
		__m512i v1, v2, v3, v4;

		for (int64 l1 = n - 32; i <= l1; i += 32)
		{
			v1 = _mm512_loadu_si512(A); A += 8;
			v2 = _mm512_loadu_si512(A); A += 8;
			v3 = _mm512_loadu_si512(A); A += 8;
			v4 = _mm512_loadu_si512(A); A += 8;

			min1 = _mm512_min_epi64(min1, v1);
			min2 = _mm512_min_epi64(min2, v2);
			min3 = _mm512_min_epi64(min3, v3);
			min4 = _mm512_min_epi64(min4, v4);
		}

		min1 = _mm512_min_epi64(_mm512_min_epi64(min1, min2), _mm512_min_epi64(min3, min4));

		for (int64 j = 0; j < 8; ++j)
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

TARGET512 void SetVal512(uint* A, ushort* B, int64 n)
{
	int64 i = 0;

	if (n >= 64)
	{
		__m512i V1, V3;
		__m256i& v1 = *((__m256i*) & V1 + 0);
		__m256i& v2 = *((__m256i*) & V1 + 1);
		__m256i& v3 = *((__m256i*) & V3 + 0);
		__m256i& v4 = *((__m256i*) & V3 + 1);
		__m512i x1, x2, x3, x4;

		for (int64 l1 = n - 64; i <= l1; i += 64)
		{
			V1 = _mm512_loadu_si512(B); B += 32;
			V3 = _mm512_loadu_si512(B); B += 32;

			x1 = _mm512_cvtepu16_epi32(v1);
			x2 = _mm512_cvtepu16_epi32(v2);
			x3 = _mm512_cvtepu16_epi32(v3);
			x4 = _mm512_cvtepu16_epi32(v4);

			_mm512_storeu_si512(A, x1); A += 16;
			_mm512_storeu_si512(A, x2); A += 16;
			_mm512_storeu_si512(A, x3); A += 16;
			_mm512_storeu_si512(A, x4); A += 16;
		}
	}
	
	for (; i < n; ++i)
		*A++ = *B++;
}

TARGET512 void AddExponent512(int64& slog, __m512d& val)
{
	__m512i& vv = *(__m512i*)&val;
	__m512i mask1 = _mm512_set1_epi64(0x7FF0000000000000);
	__m512i mask2 = _mm512_set1_epi64(0x800FFFFFFFFFFFFF);
	__m512i mask3 = _mm512_set1_epi64(0x3FF0000000000000);
	__m512i subv = _mm512_set1_epi64(1023);

	__m512i t = _mm512_sub_epi64(_mm512_srli_epi64(_mm512_and_si512(vv, mask1), 52), subv);

	slog += _mm512_reduce_add_epi64(t);
	vv = _mm512_or_si512(_mm512_and_si512(vv, mask2), mask3);
}

TARGET512 void ChargeLog512(int64& slog, double& prod, __m512d& val)
{
	AddExponent512(slog, val);
	prod = prod * _mm512_reduce_mul_pd(val);

	if (prod < DOUBLE_UNDERFLOW || prod > DOUBLE_OVERFLOW) [[unlikely]]
		AddExponent(slog, prod);
}

TARGET512 double LogProd512(double* A, int64 n)
{
	int64 i = 0;
	int64 slog = 0; double prod = 1;

	if (n >= 8)
	{
		__m512d pd = _mm512_set1_pd(1.0), dunder = _mm512_set1_pd(DOUBLE_UNDERFLOW), dover = _mm512_set1_pd(DOUBLE_OVERFLOW);
		__mmask8 flag;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			pd = _mm512_mul_pd(pd, _mm512_loadu_pd(A)); A += 8;
			flag = _mm512_cmp_pd_mask(pd, dunder, _CMP_LT_OS) | _mm512_cmp_pd_mask(dover, pd, _CMP_LT_OS);

			if (flag) [[unlikely]]
				AddExponent512(slog, pd);
		}

		ChargeLog512(slog, prod, pd);
	}

	for (; i < n; ++i, ++A)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);
	return prod;
}

TARGET512 double LogProd512(double* A, int64 n, int64 sep)
{
	int64 i = 0;
	int64 slog = 0; double prod = 1;

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

		__m512d pd = _mm512_set1_pd(1.0), dunder = _mm512_set1_pd(DOUBLE_UNDERFLOW), dover = _mm512_set1_pd(DOUBLE_OVERFLOW);
		__m512i vindex = _mm512_set_epi64(-9 * sep, -10 * sep, -11 * sep, -12 * sep, -13 * sep, -14 * sep, -15 * sep, -16 * sep);
		__mmask8 flag;

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

			pd = _mm512_mul_pd(pd, _mm512_i64gather_pd(vindex, A, sizeof(double)));
			flag = _mm512_cmp_pd_mask(pd, dunder, _CMP_LT_OS) | _mm512_cmp_pd_mask(dover, pd, _CMP_LT_OS);

			if (flag) [[unlikely]]
				AddExponent512(slog, pd);
		}

		ChargeLog512(slog, prod, pd);
		A -= sep * 8;
	}

	for (; i < n; ++i, A += sep)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);
	return prod;
}

TARGET512 double LogProdDiv512(double* A, double* B, int64 n, int64 sep)
{
	int64 i = 0;
	int64 slog = 0; double prod = 1;

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

		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;

		__m512d pd = _mm512_set1_pd(1.0), dunder = _mm512_set1_pd(DOUBLE_UNDERFLOW), dover = _mm512_set1_pd(DOUBLE_OVERFLOW);
		__m512i vindex = _mm512_set_epi64(-9 * sep, -10 * sep, -11 * sep, -12 * sep, -13 * sep, -14 * sep, -15 * sep, -16 * sep);
		__mmask8 flag;

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

			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;

			pd = _mm512_mul_pd(pd, _mm512_div_pd(_mm512_i64gather_pd(vindex, A, sizeof(double)), _mm512_i64gather_pd(vindex, B, sizeof(double))));
			flag = _mm512_cmp_pd_mask(pd, dunder, _CMP_LT_OS) | _mm512_cmp_pd_mask(dover, pd, _CMP_LT_OS);

			if (flag) [[unlikely]]
				AddExponent512(slog, pd);
		}

		ChargeLog512(slog, prod, pd);
		A -= sep * 8;
		B -= sep * 8;
	}

	for (; i < n; ++i, A += sep, B += sep)
		ChargeLog(slog, prod, *A / *B);

	CloseLog(slog, prod);
	return prod;
}

TARGET512 int64 CountNonZero512(byte* A, int64 n)
{
	int64 i = 0;
	uint64 r1 = 0, r2 = 0, r3 = 0, r4 = 0;

	if (n >= 256)
	{
		uint64 x1 = 0, x2 = 0, x3 = 0, x4 = 0;
		__m512i z = _mm512_setzero_si512();
		__m512i v1, v2, v3, v4;

		for (int64 l1 = n - 256; i <= l1; i += 256)
		{
			v1 = _mm512_loadu_si512(A); A += 64;
			v2 = _mm512_loadu_si512(A); A += 64;
			v3 = _mm512_loadu_si512(A); A += 64;
			v4 = _mm512_loadu_si512(A); A += 64;

			x1 = _mm512_cmpgt_epu8_mask(v1, z);
			x2 = _mm512_cmpgt_epu8_mask(v2, z);
			x3 = _mm512_cmpgt_epu8_mask(v3, z);
			x4 = _mm512_cmpgt_epu8_mask(v4, z);

			r1 += _mm_popcnt_u64(x1);
			r2 += _mm_popcnt_u64(x2);
			r3 += _mm_popcnt_u64(x3);
			r4 += _mm_popcnt_u64(x4);
		}

		r1 = r1 + r2 + r3 + r4;
	}

	for (; i < n; ++i, ++A)
		if (*A) r1++;

	return (int64)r1;
}

TARGET512 double Sum512(double* A, int64 n)
{
	int64 i = 0;
	double re = 0;

	if (n >= 32)
	{
		__m512d s1 = _mm512_setzero_pd();
		__m512d s2 = _mm512_setzero_pd();
		__m512d s3 = _mm512_setzero_pd();
		__m512d s4 = _mm512_setzero_pd();
		__m512d v1, v2, v3, v4;

		for (int64 l1 = n - 32; i <= l1; i += 32)
		{
			v1 = _mm512_loadu_pd(A); A += 8;
			v2 = _mm512_loadu_pd(A); A += 8;
			v3 = _mm512_loadu_pd(A); A += 8;
			v4 = _mm512_loadu_pd(A); A += 8;

			s1 = _mm512_add_pd(s1, v1);
			s2 = _mm512_add_pd(s2, v2);
			s3 = _mm512_add_pd(s3, v3);
			s4 = _mm512_add_pd(s4, v4);
		}

		s1 = _mm512_add_pd(_mm512_add_pd(s1, s2), _mm512_add_pd(s3, s4));
		re = _mm512_reduce_add_pd(s1);
	}

	for (; i < n; ++i)
		re += *A++;

	return re;
}

TARGET512 int64 Sum512(byte* A, int64 n)
{
	uint64 re = 0;
	int64 i = 0;

	if (n >= 256)
	{
		__m512i s1 = _mm512_setzero_si512();
		__m512i s2 = _mm512_setzero_si512();
		__m512i s3 = _mm512_setzero_si512();
		__m512i s4 = _mm512_setzero_si512();
		__m512i z = _mm512_setzero_si512();
		__m512i v1, v2, v3, v4;

		for (int64 l1 = n - 256; i <= l1; i += 256)
		{
			v1 = _mm512_loadu_si512(A); A += 64;
			v2 = _mm512_loadu_si512(A); A += 64;
			v3 = _mm512_loadu_si512(A); A += 64;
			v4 = _mm512_loadu_si512(A); A += 64;

			s1 = _mm512_add_epi64(s1, _mm512_sad_epu8(v1, z));
			s2 = _mm512_add_epi64(s2, _mm512_sad_epu8(v2, z));
			s3 = _mm512_add_epi64(s3, _mm512_sad_epu8(v3, z));
			s4 = _mm512_add_epi64(s4, _mm512_sad_epu8(v4, z));
		}

		s1 = _mm512_add_epi64(_mm512_add_epi64(s1, s2), _mm512_add_epi64(s3, s4));
		re += _mm512_reduce_add_epi64(s1);
	}

	for (; i < n; ++i)
		re += *A++;

	return re;
}

TARGET512 double Sum512(double* A, int64 n, int64 sep)
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

		__m512d s = _mm512_setzero_pd();
		__m512i vindex = _mm512_set_epi64(-9 * sep, -10 * sep, -11 * sep, -12 * sep, -13 * sep, -14 * sep, -15 * sep, -16 * sep);

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

			s = _mm512_add_pd(s, _mm512_i64gather_pd(vindex, A, sizeof(double)));
		}

		re = _mm512_reduce_add_pd(s);
		A -= sep * 8;
	}

	for (; i < n; ++i, A += sep)
		re += *A;

	return re;
}

TARGET512 void Sum512(double* A, double** B, int64 k, int64 n)
{
	int64 i = 0;

	if (n >= 8)
	{
		__m512d a;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			a = _mm512_setzero_pd();

			for (int64 j = 0; j < k; ++j)
				a = _mm512_add_pd(a, _mm512_loadu_pd(&B[j][i + 0]));

			_mm512_storeu_pd(&A[i + 0], a);
		}
	}

	for (; i < n; ++i)
	{
		A[i] = 0;
		for (int64 j = 0; j < k; ++j)
			A[i] += B[j][i];
	}
}

TARGET512 double Prod512(double* A, int64 n)
{
	int64 i = 0;
	double re = 1;

	if (n >= 32)
	{
		__m512d pd1 = _mm512_set1_pd(1.0);
		__m512d pd2 = _mm512_set1_pd(1.0);
		__m512d pd3 = _mm512_set1_pd(1.0);
		__m512d pd4 = _mm512_set1_pd(1.0);
		__m512d v1, v2, v3, v4;

		for (int64 l1 = n - 32; i <= l1; i += 32)
		{
			v1 = _mm512_loadu_pd(A); A += 8;
			v2 = _mm512_loadu_pd(A); A += 8;
			v3 = _mm512_loadu_pd(A); A += 8;
			v4 = _mm512_loadu_pd(A); A += 8;

			pd1 = _mm512_mul_pd(pd1, v1);
			pd2 = _mm512_mul_pd(pd2, v2);
			pd3 = _mm512_mul_pd(pd3, v3);
			pd4 = _mm512_mul_pd(pd4, v4);
		}

		pd1 = _mm512_mul_pd(_mm512_mul_pd(pd1, pd2), _mm512_mul_pd(pd3, pd4));
		re = _mm512_reduce_mul_pd(pd1);
	}

	for (; i < n; ++i)
		re *= *A++;

	return re;
}

TARGET512 double Prod512(double* A, int64 n, int64 sep)
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

		__m512d pd = _mm512_set1_pd(1.0);
		__m512i vindex = _mm512_set_epi64(-9 * sep, -10 * sep, -11 * sep, -12 * sep, -13 * sep, -14 * sep, -15 * sep, -16 * sep);

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

			pd = _mm512_mul_pd(pd, _mm512_i64gather_pd(vindex, A, sizeof(double)));
		}

		re = _mm512_reduce_mul_pd(pd);
		A -= sep * 8;
	}

	for (; i < n; ++i, A += sep)
		re *= *A;

	return re;
}

TARGET512 double SumSquare512(double* A, int64 n)
{
	int64 i = 0;
	double re = 0;

	if (n >= 32)
	{
		__m512d s1 = _mm512_setzero_pd();
		__m512d s2 = _mm512_setzero_pd();
		__m512d s3 = _mm512_setzero_pd();
		__m512d s4 = _mm512_setzero_pd();
		__m512d v1, v2, v3, v4;

		for (int64 l1 = n - 32; i <= l1; i += 32)
		{
			v1 = _mm512_loadu_pd(A); A += 8;
			v2 = _mm512_loadu_pd(A); A += 8;
			v3 = _mm512_loadu_pd(A); A += 8;
			v4 = _mm512_loadu_pd(A); A += 8;

			s1 = _mm512_fmadd_pd(v1, v1, s1);
			s2 = _mm512_fmadd_pd(v2, v2, s2);
			s3 = _mm512_fmadd_pd(v3, v3, s3);
			s4 = _mm512_fmadd_pd(v4, v4, s4);
		}

		s1 = _mm512_add_pd(_mm512_add_pd(s1, s2), _mm512_add_pd(s3, s4));
		re = _mm512_reduce_add_pd(s1);
	}

	for (; i < n; ++i, ++A)
		re += *A * *A;

	return re;
}

TARGET512 int64 SumSquare512(byte* A, int64 n)
{
	int64 i = 0;
	uint64 re = 0;

	if (n >= 128)
	{
		__m512i A1, A3;
		__m512i s1 = _mm512_setzero_si512();
		__m512i s2 = _mm512_setzero_si512();
		__m512i s3 = _mm512_setzero_si512();
		__m512i s4 = _mm512_setzero_si512();
		__m256i& a1 = simd_u256(A1, 0);
		__m256i& a2 = simd_u256(A1, 1);
		__m256i& a3 = simd_u256(A3, 0);
		__m256i& a4 = simd_u256(A3, 1);

		for (int64 l1 = n - 128; i <= l1; i += 128)
		{
			A1 = _mm512_loadu_si512((__m512i*)A); A += 64;
			A3 = _mm512_loadu_si512((__m512i*)A); A += 64;

			A1 = _mm512_maddubs_epi16(A1, A1);
			A3 = _mm512_maddubs_epi16(A3, A3);

			s1 = _mm512_add_epi32(s1, _mm512_cvtepi16_epi32(a1));
			s2 = _mm512_add_epi32(s2, _mm512_cvtepi16_epi32(a2));
			s3 = _mm512_add_epi32(s3, _mm512_cvtepi16_epi32(a3));
			s4 = _mm512_add_epi32(s4, _mm512_cvtepi16_epi32(a4));
		}

		s1 = _mm512_add_epi32(_mm512_add_epi32(s1, s2), _mm512_add_epi32(s3, s4));
		re = _mm512_reduce_add_epi32(s1);
	}

	for (; i < n; ++i, ++A)
		re += *A * *A;

	return re;
}

TARGET512 void SumSumSquare512(double* A, int64 n, double& sum, double& sumsq)
{
	int64 i = 0;
	sum = sumsq = 0;

	if (n >= 16)
	{
		__m512d s1 = _mm512_setzero_pd(), sq1 = _mm512_setzero_pd();
		__m512d s2 = _mm512_setzero_pd(), sq2 = _mm512_setzero_pd();
		__m512d a1, a2;

		for (int64 l1 = n - 16; i <= l1; i += 16)
		{
			a1 = _mm512_loadu_pd(A); A += 8;
			s1 = _mm512_add_pd(a1, s1);
			sq1 = _mm512_fmadd_pd(a1, a1, sq1);

			a2 = _mm512_loadu_pd(A); A += 8;
			s2 = _mm512_add_pd(a2, s2);
			sq2 = _mm512_fmadd_pd(a2, a2, sq2);
		}

		s1 = _mm512_add_pd(s1, s2);
		sq1 = _mm512_add_pd(sq1, sq2);

		sum = _mm512_reduce_add_pd(s1);
		sumsq = _mm512_reduce_add_pd(sq1);
	}

	for (; i < n; ++i, ++A)
	{
		sum += *A;
		sumsq += *A * *A;
	}
}

TARGET512 double SumProdDiv512(double* A1, double* A2, double* B, int64 sep, int64 n)
{
	int64 i = 0;
	double re1 = 0, re2 = 0;

	if (n >= 8)
	{
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;

		__m512d s1 = _mm512_setzero_pd(), s2 = _mm512_setzero_pd(), b;
		__m512i vindex = _mm512_set_epi64(-9 * sep, -10 * sep, -11 * sep, -12 * sep, -13 * sep, -14 * sep, -15 * sep, -16 * sep);
		__m512d v1, v2;

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;

			b = _mm512_i64gather_pd(vindex, B, sizeof(double));

			v1 = _mm512_loadu_pd(A1); A1 += 8;
			v2 = _mm512_loadu_pd(A2); A2 += 8;

			s1 = _mm512_fmadd_pd(v1, b, s1); 
			s2 = _mm512_fmadd_pd(v2, b, s2);
		}

		re1 = _mm512_reduce_add_pd(s1);
		re2 = _mm512_reduce_add_pd(s2);
		B -= sep * 8;
	}

	for (; i < n; ++i, A1++, A2++, B += sep)
	{
		re1 += *A1 * *B;
		re2 += *A2 * *B;
	}
	return re1 / re2;
}

TARGET512 double SumProd512(double* A, double* B, int64 sep, int64 n)
{
	int64 i = 0;
	double re = 0;

	if (n >= 8)
	{
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
		_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;

		__m512d s = _mm512_setzero_pd();
		__m512i vindex = _mm512_set_epi64(-9 * sep, -10 * sep, -11 * sep, -12 * sep, -13 * sep, -14 * sep, -15 * sep, -16 * sep);

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;
			_mm_prefetch((const char*)B, _MM_HINT_T0); B += sep;

			s = _mm512_fmadd_pd(_mm512_loadu_pd(A), _mm512_i64gather_pd(vindex, B, sizeof(double)), s); A += 8;
		}

		re = _mm512_reduce_add_pd(s); 
		B -= sep * 8;
	}

	for (; i < n; ++i, A++, B += sep)
		re += *A * *B;

	return re;
}

TARGET512 double SumProd512(double* A, double* B, int64 n)
{
	int64 i = 0;
	double re = 0;

	if (n >= 32)
	{
		__m512d s1 = _mm512_setzero_pd();
		__m512d s2 = _mm512_setzero_pd();
		__m512d s3 = _mm512_setzero_pd();
		__m512d s4 = _mm512_setzero_pd();
		__m512d a1, a2, a3, a4;
		__m512d b1, b2, b3, b4;

		for (int64 l1 = n - 32; i <= l1; i += 32)
		{
			a1 = _mm512_loadu_pd(A); A += 8;
			a2 = _mm512_loadu_pd(A); A += 8;
			a3 = _mm512_loadu_pd(A); A += 8;
			a4 = _mm512_loadu_pd(A); A += 8;

			b1 = _mm512_loadu_pd(B); B += 8;
			b2 = _mm512_loadu_pd(B); B += 8;
			b3 = _mm512_loadu_pd(B); B += 8;
			b4 = _mm512_loadu_pd(B); B += 8;

			s1 = _mm512_fmadd_pd(a1, b1, s1);
			s2 = _mm512_fmadd_pd(a2, b2, s2);
			s3 = _mm512_fmadd_pd(a3, b3, s3);
			s4 = _mm512_fmadd_pd(a4, b4, s4);
		}

		s1 = _mm512_add_pd(_mm512_add_pd(s1, s2), _mm512_add_pd(s3, s4));
		re = _mm512_reduce_add_pd(s1);
	}

	for (; i < n; ++i, ++A, ++B)
		re += *A * *B;

	return re;
}

TARGET512 void Add512(double* A, double* B, int64 n)
{
	int64 i = 0;

	if (n >= 32)
	{
		__m512d a1, a2, a3, a4;
		__m512d b1, b2, b3, b4;

		for (int64 l1 = n - 32; i <= l1; i += 32)
		{
			a1 = _mm512_loadu_pd(A); A += 8;
			a2 = _mm512_loadu_pd(A); A += 8;
			a3 = _mm512_loadu_pd(A); A += 8;
			a4 = _mm512_loadu_pd(A); A += 8;

			b1 = _mm512_loadu_pd(B); B += 8;
			b2 = _mm512_loadu_pd(B); B += 8;
			b3 = _mm512_loadu_pd(B); B += 8;
			b4 = _mm512_loadu_pd(B); B += 8;

			a1 = _mm512_add_pd(a1, b1);
			a2 = _mm512_add_pd(a2, b2);
			a3 = _mm512_add_pd(a3, b3);
			a4 = _mm512_add_pd(a4, b4);

			_mm512_storeu_pd(A - 32, a1);
			_mm512_storeu_pd(A - 24, a2);
			_mm512_storeu_pd(A - 16, a3);
			_mm512_storeu_pd(A - 8, a4);
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGET512 void Add512(int64* A, int64* B, int64 n)
{
	int64 i = 0;

	if (n >= 32)
	{
		__m512i a1, a2, a3, a4;
		__m512i b1, b2, b3, b4;

		for (int64 l1 = n - 32; i <= l1; i += 32)
		{
			a1 = _mm512_loadu_si512((__m512i*)A); A += 8;
			a2 = _mm512_loadu_si512((__m512i*)A); A += 8;
			a3 = _mm512_loadu_si512((__m512i*)A); A += 8;
			a4 = _mm512_loadu_si512((__m512i*)A); A += 8;

			b1 = _mm512_loadu_si512((__m512i*)B); B += 8;
			b2 = _mm512_loadu_si512((__m512i*)B); B += 8;
			b3 = _mm512_loadu_si512((__m512i*)B); B += 8;
			b4 = _mm512_loadu_si512((__m512i*)B); B += 8;

			a1 = _mm512_add_epi64(a1, b1);
			a2 = _mm512_add_epi64(a2, b2);
			a3 = _mm512_add_epi64(a3, b3);
			a4 = _mm512_add_epi64(a4, b4);

			_mm512_storeu_si512((__m512i*)(A - 32), a1);
			_mm512_storeu_si512((__m512i*)(A - 24), a2);
			_mm512_storeu_si512((__m512i*)(A - 16), a3);
			_mm512_storeu_si512((__m512i*)(A - 8), a4);
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGET512 void Add512(int* A, int* B, int64 n)
{
	int64 i = 0;

	if (n >= 64)
	{
		__m512i a1, a2, a3, a4;
		__m512i b1, b2, b3, b4;

		for (int64 l1 = n - 64; i <= l1; i += 64)
		{
			a1 = _mm512_loadu_si512((__m512i*)A); A += 16;
			a2 = _mm512_loadu_si512((__m512i*)A); A += 16;
			a3 = _mm512_loadu_si512((__m512i*)A); A += 16;
			a4 = _mm512_loadu_si512((__m512i*)A); A += 16;

			b1 = _mm512_loadu_si512((__m512i*)B); B += 16;
			b2 = _mm512_loadu_si512((__m512i*)B); B += 16;
			b3 = _mm512_loadu_si512((__m512i*)B); B += 16;
			b4 = _mm512_loadu_si512((__m512i*)B); B += 16;

			a1 = _mm512_add_epi32(a1, b1);
			a2 = _mm512_add_epi32(a2, b2);
			a3 = _mm512_add_epi32(a3, b3);
			a4 = _mm512_add_epi32(a4, b4);

			_mm512_storeu_si512((__m512i*)(A - 64), a1);
			_mm512_storeu_si512((__m512i*)(A - 48), a2);
			_mm512_storeu_si512((__m512i*)(A - 32), a3);
			_mm512_storeu_si512((__m512i*)(A - 16), a4);
		}
	}

	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGET512 void Add512(double* A, double B, int64 n)
{
	int64 i = 0;

	if (n >= 32)
	{
		__m512d b = _mm512_set1_pd(B);
		__m512d a1, a2, a3, a4;

		for (int64 l1 = n - 32; i <= l1; i += 32)
		{
			a1 = _mm512_loadu_pd(A); A += 8;
			a2 = _mm512_loadu_pd(A); A += 8;
			a3 = _mm512_loadu_pd(A); A += 8;
			a4 = _mm512_loadu_pd(A); A += 8;

			a1 = _mm512_add_pd(a1, b);
			a2 = _mm512_add_pd(a2, b);
			a3 = _mm512_add_pd(a3, b);
			a4 = _mm512_add_pd(a4, b);

			_mm512_storeu_pd(A - 32, a1);
			_mm512_storeu_pd(A - 24, a2);
			_mm512_storeu_pd(A - 16, a3);
			_mm512_storeu_pd(A - 8, a4);
		}
	}

	for (; i < n; ++i, A++)
		*A += B;
}

TARGET512 void Mul512(double* C, double* A, double* B, int64 n)
{
	int64 i = 0;

	if (n >= 32)
	{
		__m512d a1, a2, a3, a4;
		__m512d b1, b2, b3, b4;

		for (int64 l1 = n - 32; i <= l1; i += 32)
		{
			a1 = _mm512_loadu_pd(A); A += 8;
			a2 = _mm512_loadu_pd(A); A += 8;
			a3 = _mm512_loadu_pd(A); A += 8;
			a4 = _mm512_loadu_pd(A); A += 8;

			b1 = _mm512_loadu_pd(B); B += 8;
			b2 = _mm512_loadu_pd(B); B += 8;
			b3 = _mm512_loadu_pd(B); B += 8;
			b4 = _mm512_loadu_pd(B); B += 8;

			a1 = _mm512_mul_pd(a1, b1);
			a2 = _mm512_mul_pd(a2, b2);
			a3 = _mm512_mul_pd(a3, b3);
			a4 = _mm512_mul_pd(a4, b4);

			_mm512_storeu_pd(C, a1); C += 8;
			_mm512_storeu_pd(C, a2); C += 8;
			_mm512_storeu_pd(C, a3); C += 8;
			_mm512_storeu_pd(C, a4); C += 8;
		}
	}

	for (; i < n; ++i)
		*C++ = *A++ * *B++;
}

TARGET512 void Mul512(double* C, double* A, double B, int64 n)
{
	int64 i = 0;

	if (n >= 32)
	{
		__m512d b = _mm512_set1_pd(B);
		__m512d a1, a2, a3, a4;

		for (int64 l1 = n - 32; i <= l1; i += 32)
		{
			a1 = _mm512_loadu_pd(A); A += 8;
			a2 = _mm512_loadu_pd(A); A += 8;
			a3 = _mm512_loadu_pd(A); A += 8;
			a4 = _mm512_loadu_pd(A); A += 8;

			a1 = _mm512_mul_pd(a1, b);
			a2 = _mm512_mul_pd(a2, b);
			a3 = _mm512_mul_pd(a3, b);
			a4 = _mm512_mul_pd(a4, b);

			_mm512_storeu_pd(C, a1); C += 8;
			_mm512_storeu_pd(C, a2); C += 8;
			_mm512_storeu_pd(C, a3); C += 8;
			_mm512_storeu_pd(C, a4); C += 8;
		}
	}

	for (; i < n; ++i)
		*C++ = *A++ * B;
}

TARGET512 void Mul512(double* A, double B, int64 n)
{
	int64 i = 0;

	if (n >= 32)
	{
		__m512d b = _mm512_set1_pd(B);
		__m512d a1, a2, a3, a4;

		for (int64 l1 = n - 32; i <= l1; i += 32)
		{
			a1 = _mm512_loadu_pd(A); A += 8;
			a2 = _mm512_loadu_pd(A); A += 8;
			a3 = _mm512_loadu_pd(A); A += 8;
			a4 = _mm512_loadu_pd(A); A += 8;

			a1 = _mm512_mul_pd(a1, b);
			a2 = _mm512_mul_pd(a2, b);
			a3 = _mm512_mul_pd(a3, b);
			a4 = _mm512_mul_pd(a4, b);

			_mm512_storeu_pd(A - 32, a1);
			_mm512_storeu_pd(A - 24, a2);
			_mm512_storeu_pd(A - 16, a3);
			_mm512_storeu_pd(A -  8, a4);
		}
	}

	for (; i < n; ++i)
		*A++ *= B;
}

TARGET512 void AddProd512(double* C, double* A, double* B, int64 n)
{
	int64 i = 0;

	if (n >= 8)
	{
		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			_mm512_storeu_pd(C, _mm512_fmadd_pd(_mm512_loadu_pd(A), _mm512_loadu_pd(B), _mm512_loadu_pd(C)));
			A += 8; B += 8; C += 8;
		}
	}

	for (; i < n; ++i)
		*C++ += *A++ * *B++;
}

TARGET512 void AddProd512(double* C, double* A, double B, int64 n)
{
	int64 i = 0;

	if (n >= 8)
	{
		__m512d b = _mm512_set1_pd(B);

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			_mm512_storeu_pd(C, _mm512_fmadd_pd(_mm512_loadu_pd(A), b, _mm512_loadu_pd(C)));
			A += 8; C += 8;
		}
	}

	for (; i < n; ++i)
		*C++ += *A++ * B;
}

TARGET512 void Unify512(double* A, int64 n)
{
	int64 i = 0;
	double invsum = 1.0 / (Sum512(A, n) + n * MIN_FREQ);

	if (n >= 8)
	{
		__m512d minv = _mm512_set1_pd(invsum * MIN_FREQ), invs = _mm512_set1_pd(invsum);

		for (int64 l1 = n - 8; i <= l1; i += 8)
		{
			_mm512_storeu_pd(A, _mm512_fmadd_pd(_mm512_loadu_pd(A), invs, minv));
			A += 8;
		}
	}

	for (; i < n; ++i, ++A)
		*A = (*A + MIN_FREQ) * invsum;
}

TARGET512 char* StrNextIdx512(char* A, char val, int64 rep, int64 n)
{
	A++; n--;
	int64 i = 0;

	if (n >= 64)
	{
		__m512i v = _mm512_set1_epi8(val);

		for (int64 l1 = n - 64; i <= l1; i += 64, A += 64)
		{
			__m512i a = _mm512_loadu_si512(A); 
			__mmask64 mask = _mm512_cmpeq_epu8_mask(a, v);
			int64 count = (int64)_mm_popcnt_u64(mask);

			if (rep > count)
			{
				rep -= count;
				continue;
			}
			else for (;; A++, mask >>= 1)
				if ((mask & 1) && !--rep)
					return A;
		}
	}

	for (; i < n; ++i, A++)
		if (*A == val && !--rep)
			return A;

	return NULL;
}

TARGET512 int64 CountChar512(char* A, char val, int64 n)
{
	uint64 re = 0;
	int64 i = 0;

	if (n >= 64)
	{
		__m512i v = _mm512_set1_epi8(val);

		for (int64 l1 = n - 64; i <= l1; i += 64)
		{
			re += _mm_popcnt_u64(_mm512_cmpeq_epu8_mask(_mm512_loadu_si512(A), v));
			A += 64;
		}
	}

	for (; i < n; ++i, A++)
		if (*A == val) re++;
	
	return (int64)re;
}

#endif