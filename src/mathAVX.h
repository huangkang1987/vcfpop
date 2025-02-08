/* AVX Instruction Set Functions */

#pragma once
#include "vcfpop.h"

#ifndef __aarch64__

#pragma pack(push, 1)

template<typename REAL>
struct RNGAVX
{

};

template<>
struct RNGAVX<double>
{
	__m256i x[16];
	__m256i y[16];

	/* Initialize rng */
	TARGETAVX RNGAVX();

	/* Initialize rng */
	TARGETAVX RNGAVX(uint64 s, uint64 salt);

	/* Draw 64 64-bit integers in [0,n), 64*n frequencies are in arr */
	TARGETAVX void Poly(__m256d* arr, int n, __m256i* re);

	/* Draw uniform distriubted intergers */
	TARGETAVX void XorShift();

	/* Draw uniform distriubted intergers */
	template<typename INT>
	TARGETAVX void Integer(INT* re, int64 n, INT minv = 0, INT maxv = (INT)-1);

	/* Draw uniform distriubted real numbers */
	TARGETAVX void Uniform(double* re, int n, double minv = 0, double maxv = 1);

	/* Draw normal distriubted real numbers */
	TARGETAVX void Normal(double* re, int n, double mean = 0, double sd = 1);
};

template<>
struct RNGAVX<float>
{
	__m256i x[8];
	__m256i y[8];
	__m256i z[8];

	/* Initialize rng */
	TARGETAVX RNGAVX();

	/* Initialize rng */
	TARGETAVX RNGAVX(uint64 s, uint64 salt);

	/* Draw 64 64-bit integers in [0,n), 64*n frequencies are in arr */
	TARGETAVX void Poly(__m256* arr, int n, __m256i* re);

	/* Draw uniform distriubted intergers */
	TARGETAVX void XorShift();

	/* Draw uniform distriubted intergers */
	template<typename INT>
	TARGETAVX void Integer(INT* re, int64 n, INT minv = 0, INT maxv = (INT)-1);

	/* Draw uniform distriubted real numbers */
	TARGETAVX void Uniform(float* re, int n, float minv = 0, float maxv = 1);

	/* Draw normal distriubted real numbers */
	TARGETAVX void Normal(float* re, int n, float mean = 0, float sd = 1);
};

#pragma pack(pop)

static forceinline TARGETAVX int64 _mm256_reduce_add_epi64(__m256i v1)
{
	return simd_i64(v1, 0) + simd_i64(v1, 1) + simd_i64(v1, 2) + simd_i64(v1, 3);
}

static forceinline TARGETAVX double _mm256_reduce_add_pd(__m256d v1)
{
	__m128d v2 = _mm_add_pd(_mm256_extractf128_pd(v1, 0), _mm256_extractf128_pd(v1, 1));
	return simd_f64(v2, 0) + simd_f64(v2, 1);
}

static forceinline TARGETAVX double _mm256_reduce_mul_pd(__m256d v1)
{
	__m128d v2 = _mm_mul_pd(_mm256_extractf128_pd(v1, 0), _mm256_extractf128_pd(v1, 1));
	return simd_f64(v2, 0) * simd_f64(v2, 1);
}

static forceinline TARGETAVX float _mm256_reduce_add_ps(__m256 v1)
{
	__m128 v2 = _mm_add_ps(_mm256_extractf128_ps(v1, 0), _mm256_extractf128_ps(v1, 1));
	__m128 v3 = _mm_add_ps(v2, _mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(v2), 8)));
	return simd_f32(v3, 0) + simd_f32(v3, 1);
}

static forceinline TARGETAVX float _mm256_reduce_mul_ps(__m256 v1)
{
	__m128 v2 = _mm_mul_ps(_mm256_extractf128_ps(v1, 0), _mm256_extractf128_ps(v1, 1));
	__m128 v3 = _mm_mul_ps(v2, _mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(v2), 8)));
	return simd_f32(v3, 0) * simd_f32(v3, 1);
}

static forceinline TARGETAVX double _mm256_reduce_add_psd(__m256 v1)
{
	__m256d v1b = _mm256_add_pd(
		_mm256_cvtps_pd(_mm256_extractf128_ps(v1, 0)),
		_mm256_cvtps_pd(_mm256_extractf128_ps(v1, 1)));
	__m128d v2 = _mm_add_pd(_mm256_extractf128_pd(v1b, 0), _mm256_extractf128_pd(v1b, 1));
	return simd_f64(v2, 0) + simd_f64(v2, 1);
}

static forceinline TARGETAVX double _mm256_reduce_mul_psd(__m256 v1)
{
	__m256d v1b = _mm256_mul_pd(
		_mm256_cvtps_pd(_mm256_extractf128_ps(v1, 0)),
		_mm256_cvtps_pd(_mm256_extractf128_ps(v1, 1)));
	__m128d v2 = _mm_mul_pd(_mm256_extractf128_pd(v1b, 0), _mm256_extractf128_pd(v1b, 1));
	return simd_f64(v2, 0) * simd_f64(v2, 1);
}

static forceinline TARGETAVX int64 _mm256_reduce_min_epi64(__m256i v1) 
{
    return std::min(std::min(simd_i64(v1, 0), simd_i64(v1, 1)),
		            std::min(simd_i64(v1, 2), simd_i64(v1, 3)));
}

static forceinline TARGETAVX int _mm256_reduce_min_epi32(__m256i v1) 
{
    __m256i v2 = _mm256_min_epi32(v1, _mm256_permute2f128_si256(v1, v1, 0x01));
    v2 = _mm256_min_epi32(v2, _mm256_shuffle_epi32(v2, _MM_SHUFFLE(2, 3, 0, 1)));
    v2 = _mm256_min_epi32(v2, _mm256_shuffle_epi32(v2, _MM_SHUFFLE(1, 0, 3, 2)));
    return simd_i32(v2, 0);
}

static forceinline TARGETAVX double _mm256_reduce_min_pd(__m256d v1) 
{
    __m256d v2 = _mm256_min_pd(v1, _mm256_permute2f128_pd(v1, v1, 0x01));
    v2 = _mm256_min_pd(v2, _mm256_permute_pd(v2, 0b0101));
    return simd_f64(v2, 0);
}

static forceinline TARGETAVX float _mm256_reduce_min_ps(__m256 v1) 
{
    __m256 v2 = _mm256_min_ps(v1, _mm256_permute2f128_ps(v1, v1, 0x01));
    v2 = _mm256_min_ps(v2, _mm256_shuffle_ps(v2, v2, _MM_SHUFFLE(2, 3, 0, 1)));
    v2 = _mm256_min_ps(v2, _mm256_shuffle_ps(v2, v2, _MM_SHUFFLE(1, 0, 3, 2)));
    return simd_f32(v2, 0);
}

static forceinline TARGETAVX int64 _mm256_reduce_max_epi64(__m256i v1) 
{
    return std::max(std::max(simd_i64(v1, 0), simd_i64(v1, 1)),
		            std::max(simd_i64(v1, 2), simd_i64(v1, 3)));
}

static forceinline TARGETAVX int _mm256_reduce_max_epi32(__m256i v1) 
{
    __m256i v2 = _mm256_max_epi32(v1, _mm256_permute2f128_si256(v1, v1, 0x01));
    v2 = _mm256_max_epi32(v2, _mm256_shuffle_epi32(v2, _MM_SHUFFLE(2, 3, 0, 1)));
    v2 = _mm256_max_epi32(v2, _mm256_shuffle_epi32(v2, _MM_SHUFFLE(1, 0, 3, 2)));
    return simd_i32(v2, 0);
}

static forceinline TARGETAVX double _mm256_reduce_max_pd(__m256d v1) 
{
    __m256d v2 = _mm256_max_pd(v1, _mm256_permute2f128_pd(v1, v1, 0x01));
    v2 = _mm256_max_pd(v2, _mm256_permute_pd(v2, 0b0101));
    return simd_f64(v2, 0);
}

static forceinline TARGETAVX float _mm256_reduce_max_ps(__m256 v1) 
{
    __m256 v2 = _mm256_max_ps(v1, _mm256_permute2f128_ps(v1, v1, 0x01));
    v2 = _mm256_max_ps(v2, _mm256_shuffle_ps(v2, v2, _MM_SHUFFLE(2, 3, 0, 1)));
    v2 = _mm256_max_ps(v2, _mm256_shuffle_ps(v2, v2, _MM_SHUFFLE(1, 0, 3, 2)));
    return simd_f32(v2, 0);
}

static forceinline __m128i _mm256_cvtepi64_epi32x(__m256i a)
{
	return _mm_castps_si128(_mm_shuffle_ps(
		_mm_castpd_ps(_mm256_castpd256_pd128(_mm256_castsi256_pd(a))), 
		_mm_castpd_ps(_mm256_extractf128_pd(_mm256_castsi256_pd(a), 1)), 
		_MM_SHUFFLE(2, 0, 2, 0)));
}

TARGETAVX int64 GetMinIdxAVX(double* A, int64 n, double& val);

TARGETAVX int64 GetMinIdxAVX(float* A, int64 n, float& val);

TARGETAVX void GetMinMaxValAVX(double* A, int64 n, double& minv, double& maxv);

TARGETAVX void GetMinMaxValAVX(float* A, int64 n, float& minv, float& maxv);

TARGETAVX double GetMaxValAVX(double* A, int64 n);

TARGETAVX float GetMaxValAVX(float* A, int64 n);

TARGETAVX double GetMaxValAVX(double* A, int64 n, int64 sep);

TARGETAVX float GetMaxValAVX(float* A, int64 n, int64 sep);

TARGETAVX double GetMinValAVX(double* A, int64 n);

TARGETAVX float GetMinValAVX(float* A, int64 n);

TARGETAVX int64 GetMinValAVX(int64* A, int64 n);

TARGETAVX void SetValAVX(uint* a, ushort* b, int64 n);

TARGETAVX void AddExponentAVX(int64& slog, __m256d& val);

TARGETAVX void AddExponentAVX(int64& slog, __m256& val);

TARGETAVX void ChargeLogAVX(int64& slog, double& prod, __m256d& val);

TARGETAVX void ChargeLogAVX(int64& slog, double& prod, __m256& val);

TARGETAVX double LogProdAVX(double* A, int64 n);

TARGETAVX double LogProdAVX(float* A, int64 n);

TARGETAVX double LogProdAVX(double* A, int64 n, int64 sep);

TARGETAVX double LogProdAVX(float* A, int64 n, int64 sep);

TARGETAVX double LogProdDivAVX(double* A, double* B, int64 n, int64 sep);

TARGETAVX double LogProdDivAVX(float* A, float* B, int64 n, int64 sep);

TARGETAVX int64 CountNonZeroAVX(byte* A, int64 n);

TARGETAVX double SumAVX(double* A, int64 n);

TARGETAVX double SumAVX(float* A, int64 n);

TARGETAVX int64 SumAVX(byte* A, int64 n);

TARGETAVX double SumAVX(double* A, int64 n, int64 sep);

TARGETAVX double SumAVX(float* A, int64 n, int64 sep);

TARGETAVX double ProdAVX(double* A, int64 n);

TARGETAVX double ProdAVX(float* A, int64 n);

TARGETAVX double ProdAVX(double* A, int64 n, int64 sep);

TARGETAVX double ProdAVX(float* A, int64 n, int64 sep);

TARGETAVX double SumSquareAVX(double* A, int64 n);

TARGETAVX double SumSquareAVX(float* A, int64 n);

TARGETAVX int64 SumSquareAVX(byte* A, int64 n);

TARGETAVX void SumSumSquareAVX(double* A, int64 n, double& sum, double& sumsq);

TARGETAVX void SumSumSquareAVX(float* A, int64 n, double& sum, double& sumsq);

TARGETAVX double SumProdDivAVX(double* A1, double* A2, double* B, int64 sep, int64 n);

TARGETAVX double SumProdDivAVX(double* A1, float* A2, float* B, int64 sep, int64 n);

TARGETAVX double SumProdDivAVX(float* A1, float* A2, float* B, int64 sep, int64 n);

TARGETAVX float SumProdDivAVXx(float* A1, float* A2, float* B, int64 sep, int64 n);

TARGETAVX double SumProdAVX(double* A, double* B, int64 sep, int64 n);

TARGETAVX double SumProdAVX(float* A, float* B, int64 sep, int64 n);

TARGETAVX double SumProdAVX(double* A, double* B, int64 n);

TARGETAVX double SumProdAVX(double* A, double* B, double* C, int64 n);

TARGETAVX double SumProdAVX(float* A, float* B, int64 n);

TARGETAVX float SumProdAVX(float* A, float* B, float* C, int64 n);

TARGETAVX double SumSqProdAVX(double* A, double* B, int64 n);

TARGETAVX float SumSqProdAVX(float* A, float* B, int64 n);

TARGETAVX void AddAVX(double* A, double* B, int64 n);

TARGETAVX void AddAVX(float* A, float* B, int64 n);

TARGETAVX void AddAVX(int64* A, int64* B, int64 n);

TARGETAVX void AddAVX(int* A, int* B, int64 n);

TARGETAVX void AddAVX(int* A, int B, int64 n);

TARGETAVX void AddAVX(double* A, double B, int64 n);

TARGETAVX void AddAVX(float* A, float B, int64 n);

TARGETAVX void MulAVX(double* A, double* B, double* C, int64 n);

TARGETAVX void MulAVX(float* A, float* B, float* C, int64 n);

TARGETAVX void MulAVX(double* A, double* B, double C, int64 n);

TARGETAVX void MulAVX(float* A, float* B, float C, int64 n);

TARGETAVX void MulAVX(double* A, double B, int64 n);

TARGETAVX void MulAVX(float* A, float B, int64 n);

TARGETAVX void DivAVX(double* A, double B, double* C, int64 n);

TARGETAVX void DivAVX(float* A, float B, float* C, int64 n);

TARGETAVX void DivAVX(double* A, double* B, double* C, int64 n);

TARGETAVX void DivAVX(float* A, float* B, float* C, int64 n);

TARGETAVX void AddProdAVX(double* A, double* B, double* C, int64 n);

TARGETAVX void AddProdAVX(float* A, float* B, float* C, int64 n);

TARGETAVX void AddProdAVX(double* A, double* B, double C, int64 n);

TARGETAVX void AddProdAVX(double* A, float* B, double C, int64 n);

TARGETAVX void AddProdAVX(float* A, float* B, float C, int64 n);

TARGETAVX void UnifyAVX(double* A, int64 n);

TARGETAVX void UnifyAVX(float* A, int64 n);

TARGETAVX char* StrNextIdxAVX(char* A, char val, int64 rep, int64 n);

TARGETAVX int64 CountCharAVX(char* A, char val, int64 n);

TARGETAVX float SumAVXx(float* A, int64 n);

TARGETAVX float SumAVXx(float* A, int64 n, int64 sep);

TARGETAVX float ProdAVXx(float* A, int64 n);

TARGETAVX float ProdAVXx(float* A, int64 n, int64 sep);

TARGETAVX float SumProdAVXx(float* A, float* B, int64 sep, int64 n);

TARGETAVX float SumProdAVXx(float* A, float* B, int64 n);

TARGETAVX void DiagQuadFormAVX(double* res, double* A, double* D, int64 m, int64 n);

TARGETAVX void DiagQuadFormAVX(float* res, float* A, float* D, int64 m, int64 n);

TARGETAVX void DiagQuadFormAVX(double* res, double* A, double* D, double* B, int64 m, int64 n);

TARGETAVX void DiagQuadFormAVX(float* res, float* A, float* D, float* B, int64 m, int64 n);

TARGETAVX void DiagQuadFormAVX(double* res, double* A, double* D, int64 n);

TARGETAVX void DiagQuadFormAVX(float* res, float* A, float* D, int64 n);

TARGETAVX void MatrixMulAVX(double* res, double* A, double* B, int64 m, int64 n, int64 p);

TARGETAVX void MatrixMulAVX(float* res, float* A, float* B, int64 m, int64 n, int64 p);

#endif