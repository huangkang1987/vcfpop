/* AVX512 Instruction Set Functions */

#pragma once
#include "vcfpop.h"

#ifndef __aarch64__

#pragma pack(push, 1)

template<typename REAL>
struct RNG512
{

};

template<>
struct RNG512<double>
{
	__m512i x[8];
	__m512i y[8];

	/* Initialize rng */
	TARGET512 RNG512();

	/* Initialize rng */
	TARGET512 RNG512(uint64 s, uint64 salt);

	/* Draw 64 64-bit integers in [0,n), 64*n frequencies are in arr */
	TARGET512 void Poly(__m512d* arr, int n, __m512i* re);

	/* Draw uniform distriubted intergers */
	TARGET512 void XorShift();

	/* Draw uniform distriubted integers */
	template<typename INT>
	TARGET512 void Integer(INT* re, int64 n, INT minv = 0, INT maxv = (INT)-1);

	/* Draw uniform distriubted real numbers */
	TARGET512 void Uniform(double* re, int n, double minv = 0, double maxv = 1);

	/* Draw uniform distriubted real numbers */
	TARGET512 void Normal(double* re, int n, double mean = 0, double sd = 1);
};

template<>
struct RNG512<float>
{
	__m512i x[4];
	__m512i y[4];
	__m512i z[4];

	/* Initialize rng */
	TARGET512 RNG512();

	/* Initialize rng */
	TARGET512 RNG512(uint64 s, uint64 salt);

	/* Draw 64 64-bit integers in [0,n), 64*n frequencies are in arr */
	TARGET512 void Poly(__m512* arr, int n, __m512i* re);

	/* Draw uniform distriubted intergers */
	TARGET512 void XorShift();

	/* Draw uniform distriubted intergers */
	template<typename INT>
	TARGET512 void Integer(INT* re, int64 n, INT minv = 0, INT maxv = (INT)-1);

	/* Draw uniform distriubted real numbers */
	TARGET512 void Uniform(float* re, int n, float minv = 0, float maxv = 1);

	/* Draw normal distriubted real numbers */
	TARGET512 void Normal(float* re, int n, float mean = 0, float sd = 1);
};

#pragma pack(pop)

static forceinline TARGET512 double __mm512_reduce_add_pd(__m512d v1)
{
	__m256d v2 = _mm256_add_pd(_mm512_extractf64x4_pd(v1, 0), _mm512_extractf64x4_pd(v1, 1));
	__m128d v3 = _mm_add_pd(_mm256_extractf128_pd(v2, 0), _mm256_extractf128_pd(v2, 1));
	return simd_f64(v3, 0) + simd_f64(v3, 1);
}

static forceinline TARGET512 double __mm512_reduce_mul_pd(__m512d v1)
{
	__m256d v2 = _mm256_mul_pd(_mm512_extractf64x4_pd(v1, 0), _mm512_extractf64x4_pd(v1, 1));
	__m128d v3 = _mm_mul_pd(_mm256_extractf128_pd(v2, 0), _mm256_extractf128_pd(v2, 1));
	return simd_f64(v3, 0) * simd_f64(v3, 1);
}

static forceinline TARGET512 float __mm512_reduce_add_ps(__m512 v1)
{
	__m256 v2 = _mm256_add_ps(_mm512_extractf32x8_ps(v1, 0), _mm512_extractf32x8_ps(v1, 1));
	__m128 v3 = _mm_add_ps(_mm256_extractf128_ps(v2, 0), _mm256_extractf128_ps(v2, 1));
	__m128 v4 = _mm_add_ps(v3, _mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(v3), 8)));
	return simd_f32(v4, 0) + simd_f32(v4, 1);
}

static forceinline TARGET512 float __mm512_reduce_mul_ps(__m512 v1)
{
	__m256 v2 = _mm256_mul_ps(_mm512_extractf32x8_ps(v1, 0), _mm512_extractf32x8_ps(v1, 1));
	__m128 v3 = _mm_mul_ps(_mm256_extractf128_ps(v2, 0), _mm256_extractf128_ps(v2, 1));
	__m128 v4 = _mm_mul_ps(v3, _mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(v3), 8)));
	return simd_f32(v4, 0) * simd_f32(v4, 1);
}

static forceinline TARGET512 double __mm512_reduce_add_psd(__m512 v1)
{
	__m512d v1b = _mm512_add_pd(
		_mm512_cvtps_pd(_mm512_extractf32x8_ps(v1, 0)),
		_mm512_cvtps_pd(_mm512_extractf32x8_ps(v1, 1)));
	__m256d v2 = _mm256_add_pd(_mm512_extractf64x4_pd(v1b, 0), _mm512_extractf64x4_pd(v1b, 1));
	__m128d v3 = _mm_add_pd(_mm256_extractf128_pd(v2, 0), _mm256_extractf128_pd(v2, 1));
	return simd_f64(v3, 0) + simd_f64(v3, 1);
}

static forceinline TARGET512 double __mm512_reduce_mul_psd(__m512 v1)
{
	__m512d v1b = _mm512_mul_pd(
		_mm512_cvtps_pd(_mm512_extractf32x8_ps(v1, 0)),
		_mm512_cvtps_pd(_mm512_extractf32x8_ps(v1, 1)));
	__m256d v2 = _mm256_mul_pd(_mm512_extractf64x4_pd(v1b, 0), _mm512_extractf64x4_pd(v1b, 1));
	__m128d v3 = _mm_mul_pd(_mm256_extractf128_pd(v2, 0), _mm256_extractf128_pd(v2, 1));
	return simd_f64(v3, 0) * simd_f64(v3, 1);
}

TARGET512 int64 GetMinIdx512(double* A, int64 n, double& val);

TARGET512 int64 GetMinIdx512(float* A, int64 n, float& val);

TARGET512 void GetMinMaxVal512(double* A, int64 n, double& minv, double& maxv);

TARGET512 void GetMinMaxVal512(float* A, int64 n, float& minv, float& maxv);

TARGET512 double GetMaxVal512(double* A, int64 n);

TARGET512 float GetMaxVal512(float* A, int64 n);

TARGET512 double GetMaxVal512(double* A, int64 n, int64 sep);

TARGET512 float GetMaxVal512(float* A, int64 n, int64 sep);

TARGET512 double GetMinVal512(double* A, int64 n);

TARGET512 float GetMinVal512(float* A, int64 n);

TARGET512 int64 GetMinVal512(int64* A, int64 n);

TARGET512 void SetVal512(uint* a, ushort* b, int64 n);

TARGET512 void AddExponent512(int64& slog, __m512d& val);

TARGET512 void AddExponent512(int64& slog, __m512& val);

TARGET512 void ChargeLog512(int64& slog, double& prod, __m512d& val);

TARGET512 void ChargeLog512(int64& slog, double& prod, __m512& val);

TARGET512 double LogProd512(double* A, int64 n);

TARGET512 double LogProd512(float* A, int64 n);

TARGET512 double LogProd512(double* A, int64 n, int64 sep);

TARGET512 double LogProd512(float* A, int64 n, int64 sep);

TARGET512 double LogProdDiv512(double* A, double* B, int64 n, int64 sep);

TARGET512 double LogProdDiv512(float* A, float* B, int64 n, int64 sep);

TARGET512 int64 CountNonZero512(byte* A, int64 n);

TARGET512 double Sum512(double* A, int64 n);

TARGET512 double Sum512(float* A, int64 n);

TARGET512 int64 Sum512(byte* A, int64 n);

TARGET512 double Sum512(double* A, int64 n, int64 sep);

TARGET512 double Sum512(float* A, int64 n, int64 sep);

TARGET512 double Prod512(double* A, int64 n);

TARGET512 double Prod512(float* A, int64 n);

TARGET512 double Prod512(double* A, int64 n, int64 sep);

TARGET512 double Prod512(float* A, int64 n, int64 sep);

TARGET512 double SumSquare512(double* A, int64 n);

TARGET512 double SumSquare512(float* A, int64 n);

TARGET512 int64 SumSquare512(byte* A, int64 n);

TARGET512 void SumSumSquare512(double* A, int64 n, double& sum, double& sumsq);

TARGET512 void SumSumSquare512(float* A, int64 n, double& sum, double& sumsq);

TARGET512 double SumProdDiv512(double* A1, double* A2, double* B, int64 sep, int64 n);

TARGET512 double SumProdDiv512(double* A1, float* A2, float* B, int64 sep, int64 n);

TARGET512 double SumProdDiv512(float* A1, float* A2, float* B, int64 sep, int64 n);

TARGET512 float SumProdDiv512x(float* A1, float* A2, float* B, int64 sep, int64 n);

TARGET512 double SumProd512(double* A, double* B, int64 sep, int64 n);

TARGET512 double SumProd512(float* A, float* B, int64 sep, int64 n);

TARGET512 double SumProd512(double* A, double* B, int64 n);

TARGET512 double SumProd512(double* A, double* B, double* C, int64 n);

TARGET512 double SumProd512(float* A, float* B, int64 n);

TARGET512 float SumProd512(float* A, float* B, float* C, int64 n);

TARGET512 double SumSqProd512(double* A, double* B, int64 n);

TARGET512 float SumSqProd512(float* A, float* B, int64 n);

TARGET512 void Add512(double* A, double* B, int64 n);

TARGET512 void Add512(float* A, float* B, int64 n);

TARGET512 void Add512(int64* A, int64* B, int64 n);

TARGET512 void Add512(int* A, int* B, int64 n);

TARGET512 void Add512(int* A, int B, int64 n);

TARGET512 void Add512(double* A, double B, int64 n);

TARGET512 void Add512(float* A, float B, int64 n);

TARGET512 void Mul512(double* A, double* B, double* C, int64 n);

TARGET512 void Mul512(float* A, float* B, float* C, int64 n);

TARGET512 void Mul512(double* A, double* B, double C, int64 n);

TARGET512 void Mul512(float* A, float* B, float C, int64 n);

TARGET512 void Mul512(double* A, double B, int64 n);

TARGET512 void Mul512(float* A, float B, int64 n);

TARGET512 void Div512(double* A, double B, double* C, int64 n);

TARGET512 void Div512(float* A, float B, float* C, int64 n);

TARGET512 void Div512(double* A, double* B, double* C, int64 n);

TARGET512 void Div512(float* A, float* B, float* C, int64 n);

TARGET512 void AddProd512(double* A, double* B, double* C, int64 n);

TARGET512 void AddProd512(float* A, float* B, float* C, int64 n);

TARGET512 void AddProd512(double* A, double* B, double C, int64 n);

TARGET512 void AddProd512(double* A, float* B, double C, int64 n);

TARGET512 void AddProd512(float* A, float* B, float C, int64 n);

TARGET512 void Unify512(double* A, int64 n);

TARGET512 void Unify512(float* A, int64 n);

TARGET512 char* StrNextIdx512(char* A, char val, int64 rep, int64 n);

TARGET512 int64 CountChar512(char* A, char val, int64 n);

TARGET512 float Sum512x(float* A, int64 n);

TARGET512 float Sum512x(float* A, int64 n, int64 sep);

TARGET512 float Prod512x(float* A, int64 n);

TARGET512 float Prod512x(float* A, int64 n, int64 sep);

TARGET512 float SumProd512x(float* A, float* B, int64 sep, int64 n);

TARGET512 float SumProd512x(float* A, float* B, int64 n);

TARGET512 void DiagQuadForm512(double* res, double* A, double* D, int64 m, int64 n);

TARGET512 void DiagQuadForm512(float* res, float* A, float* D, int64 m, int64 n);

TARGET512 void DiagQuadForm512(double* res, double* A, double* D, double* B, int64 m, int64 n);

TARGET512 void DiagQuadForm512(float* res, float* A, float* D, float* B, int64 m, int64 n);

TARGET512 void DiagQuadForm512(double* res, double* A, double* D, int64 n);

TARGET512 void DiagQuadForm512(float* res, float* A, float* D, int64 n);

TARGET512 void MatrixMul512(double* res, double* A, double* B, int64 m, int64 n, int64 p);

TARGET512 void MatrixMul512(float* res, float* A, float* B, int64 m, int64 n, int64 p);

#endif