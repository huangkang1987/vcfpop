/* SSE Instruction Set Functions */

#pragma once
#include "vcfpop.h"

#ifndef __aarch64__

#pragma pack(push, 1)

template<typename REAL>
struct RNGSSE
{

};

template<>
struct RNGSSE<double>
{
	__m128i x[32];
	__m128i y[32];

	/* Initialize rng */
	TARGETSSE RNGSSE();

	/* Initialize rng */
	TARGETSSE RNGSSE(uint64 s, uint64 salt);

	/* Draw 64 64-bit integers in [0,n), 64*n frequencies are in arr */
	TARGETSSE void Poly(__m128d* arr, int n, __m128i* re);

	/* Draw uniform distriubted intergers */
	TARGETSSE void XorShift();

	/* Draw uniform distriubted intergers */
	template<typename INT>
	TARGETSSE void Integer(INT* re, int64 n, INT minv = 0, INT maxv = (INT)-1);

	/* Draw uniform distriubted real numbers */
	TARGETSSE void Uniform(double* re, int n, double minv = 0, double maxv = 1);

	/* Draw normal distriubted real numbers */
	TARGETSSE void Normal(double* re, int n, double mean = 0, double sd = 1);
};

template<>
struct RNGSSE<float>
{
	__m128i x[16];
	__m128i y[16];
	__m128i z[16];

	/* Initialize rng */
	TARGETSSE RNGSSE();

	/* Initialize rng */
	TARGETSSE RNGSSE(uint64 s, uint64 salt);

	/* Draw 64 64-bit integers in [0,n), 64*n frequencies are in arr */
	TARGETSSE void Poly(__m128* arr, int n, __m128i* re);

	/* Draw uniform distriubted intergers */
	TARGETSSE void XorShift();

	/* Draw uniform distriubted intergers */
	template<typename INT>
	TARGETSSE void Integer(INT* re, int64 n, INT minv = 0, INT maxv = (INT)-1);

	/* Draw uniform distriubted real numbers */
	TARGETSSE void Uniform(float* re, int n, float minv = 0, float maxv = 1);

	/* Draw normal distriubted real numbers */
	TARGETSSE void Normal(float* re, int n, float mean = 0, float sd = 1);
};

#pragma pack(pop)

static forceinline TARGETSSE int64 _mm_reduce_add_epi64(__m128i v1)
{
	return simd_i64(v1, 0) + simd_i64(v1, 1);
}

static forceinline TARGETSSE double _mm_reduce_add_pd(__m128d v1)
{
	return simd_f64(v1, 0) + simd_f64(v1, 1);
}

static forceinline TARGETSSE double _mm_reduce_mul_pd(__m128d v1)
{
	return simd_f64(v1, 0) * simd_f64(v1, 1);
}

static forceinline TARGETSSE float _mm_reduce_add_ps(__m128 v1)
{
	__m128 v2 = _mm_add_ps(v1, _mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(v1), 8)));
	return simd_f32(v2, 0) + simd_f32(v2, 1);
}

static forceinline TARGETSSE float _mm_reduce_mul_ps(__m128 v1)
{
	__m128 v2 = _mm_mul_ps(v1, _mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(v1), 8)));
	return simd_f32(v2, 0) * simd_f32(v2, 1);
}

static forceinline TARGETSSE double _mm_reduce_add_psd(__m128 v1)
{
	__m128d v1b = _mm_add_pd(_mm_cvtps_pd(v1),
				 _mm_cvtps_pd(_mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(v1), 8))));
	return simd_f64(v1b, 0) + simd_f64(v1b, 1);
}

static forceinline TARGETSSE double _mm_reduce_mul_psd(__m128 v1)
{
	__m128d v1b = _mm_mul_pd(_mm_cvtps_pd(v1),
				 _mm_cvtps_pd(_mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(v1), 8))));
	return simd_f64(v1b, 0) * simd_f64(v1b, 1);
}

static forceinline TARGETSSE int64 _mm_reduce_min_epi64(__m128i v1) 
{
    return std::min(simd_i64(v1, 0), simd_i64(v1, 1));
}

static forceinline TARGETSSE int _mm_reduce_min_epi32(__m128i v1) 
{
	__m128i v2 = _mm_min_epi32(v1, _mm_shuffle_epi32(v1, _MM_SHUFFLE(2, 3, 0, 1)));
            v2 = _mm_min_epi32(v2, _mm_shuffle_epi32(v2, _MM_SHUFFLE(1, 0, 3, 2))); 
	return simd_i32(v2, 0);
}

static forceinline TARGETSSE double _mm_reduce_min_pd(__m128d v1) 
{
	__m128d v2 = _mm_min_pd(v1, _mm_shuffle_pd(v1, v1, 0b01));
	return simd_f64(v2, 0);
}

static forceinline TARGETSSE float _mm_reduce_min_ps(__m128 v1) 
{
	__m128 v2 = _mm_min_ps(v1, _mm_shuffle_ps(v1, v1, _MM_SHUFFLE(2, 3, 0, 1)));
           v2 = _mm_min_ps(v2, _mm_shuffle_ps(v2, v2, _MM_SHUFFLE(1, 0, 3, 2))); 
	return simd_f32(v2, 0);
}

static forceinline TARGETSSE int64 _mm_reduce_max_epi64(__m128i v1) 
{
    return std::max(simd_i64(v1, 0), simd_i64(v1, 1));
}

static forceinline TARGETSSE int _mm_reduce_max_epi32(__m128i v1) 
{
	__m128i v2 = _mm_max_epi32(v1, _mm_shuffle_epi32(v1, _MM_SHUFFLE(2, 3, 0, 1)));
            v2 = _mm_max_epi32(v2, _mm_shuffle_epi32(v2, _MM_SHUFFLE(1, 0, 3, 2))); 
	return simd_i32(v2, 0);
}

static forceinline TARGETSSE double _mm_reduce_max_pd(__m128d v1) 
{
	__m128d v2 = _mm_max_pd(v1, _mm_shuffle_pd(v1, v1, 0b01));
	return simd_f64(v2, 0);
}

static forceinline TARGETSSE float _mm_reduce_max_ps(__m128 v1) 
{
	__m128 v2 = _mm_max_ps(v1, _mm_shuffle_ps(v1, v1, _MM_SHUFFLE(2, 3, 0, 1)));
           v2 = _mm_max_ps(v2, _mm_shuffle_ps(v2, v2, _MM_SHUFFLE(1, 0, 3, 2))); 
	return simd_f32(v2, 0);
}

TARGETSSE int64 GetMinIdxSSE(double* A, int64 n, double& val);

TARGETSSE int64 GetMinIdxSSE(float* A, int64 n, float& val);

TARGETSSE void GetMinMaxValSSE(double* A, int64 n, double& minv, double& maxv);

TARGETSSE void GetMinMaxValSSE(float* A, int64 n, float& minv, float& maxv);

TARGETSSE double GetMaxValSSE(double* A, int64 n);

TARGETSSE float GetMaxValSSE(float* A, int64 n);

TARGETSSE double GetMaxValSSE(double* A, int64 n, int64 sep);

TARGETSSE float GetMaxValSSE(float* A, int64 n, int64 sep);

TARGETSSE double GetMinValSSE(double* A, int64 n);

TARGETSSE float GetMinValSSE(float* A, int64 n);

TARGETSSE int64 GetMinValSSE(int64* A, int64 n);

TARGETSSE void SetValSSE(uint* a, ushort* b, int64 n);

TARGETSSE void SetValSSE(uint* A, ushort* B, int64 n);

TARGETSSE void AddExponentSSE(int64& slog, __m128d& val);

TARGETSSE void AddExponentSSE(int64& slog, __m128& val);

TARGETSSE void ChargeLogSSE(int64& slog, double& prod, __m128d& val);

TARGETSSE void ChargeLogSSE(int64& slog, double& prod, __m128& val);

TARGETSSE double LogProdSSE(double* A, int64 n);

TARGETSSE double LogProdSSE(float* A, int64 n);

TARGETSSE double LogProdSSE(double* A, int64 n, int64 sep);

TARGETSSE double LogProdSSE(float* A, int64 n, int64 sep);

TARGETSSE double LogProdDivSSE(double* A, double* B, int64 n, int64 sep);

TARGETSSE double LogProdDivSSE(float* A, float* B, int64 n, int64 sep);

TARGETSSE int64 CountNonZeroSSE(byte* A, int64 n);

TARGETSSE double SumSSE(double* A, int64 n);

TARGETSSE double SumSSE(float* A, int64 n);

TARGETSSE int64 SumSSE(byte* A, int64 n);

TARGETSSE double SumSSE(double* A, int64 n, int64 sep);

TARGETSSE double SumSSE(float* A, int64 n, int64 sep);

TARGETSSE double ProdSSE(double* A, int64 n);

TARGETSSE double ProdSSE(float* A, int64 n);

TARGETSSE double ProdSSE(double* A, int64 n, int64 sep);

TARGETSSE double ProdSSE(float* A, int64 n, int64 sep);

TARGETSSE double SumSquareSSE(double* A, int64 n);

TARGETSSE double SumSquareSSE(float* A, int64 n);

TARGETSSE int64 SumSquareSSE(byte* A, int64 n);

TARGETSSE void SumSumSquareSSE(double* A, int64 n, double& sum, double& sumsq);

TARGETSSE void SumSumSquareSSE(float* A, int64 n, double& sum, double& sumsq);

TARGETSSE double SumProdDivSSE(double* A1, double* A2, double* B, int64 sep, int64 n);

TARGETSSE double SumProdDivSSE(double* A1, float* A2, float* B, int64 sep, int64 n);

TARGETSSE double SumProdDivSSE(float* A1, float* A2, float* B, int64 sep, int64 n);

TARGETSSE float SumProdDivSSEx(float* A1, float* A2, float* B, int64 sep, int64 n);

TARGETSSE double SumProdSSE(double* A, double* B, int64 sep, int64 n);

TARGETSSE double SumProdSSE(float* A, float* B, int64 sep, int64 n);

TARGETSSE double SumProdSSE(double* A, double* B, int64 n);

TARGETSSE double SumProdSSE(double* A, double* B, double* C, int64 n);

TARGETSSE double SumProdSSE(float* A, float* B, int64 n);

TARGETSSE float SumProdSSE(float* A, float* B, float* C, int64 n);

TARGETSSE double SumSqProdSSE(double* A, double* B, int64 n);

TARGETSSE float SumSqProdSSE(float* A, float* B, int64 n);

TARGETSSE void AddSSE(double* A, double* B, int64 n);

TARGETSSE void AddSSE(float* A, float* B, int64 n);

TARGETSSE void AddSSE(int64* A, int64* B, int64 n);

TARGETSSE void AddSSE(int* A, int* B, int64 n);

TARGETSSE void AddSSE(int* A, int B, int64 n);

TARGETSSE void AddSSE(double* A, double B, int64 n);

TARGETSSE void AddSSE(float* A, float B, int64 n);

TARGETSSE void MulSSE(double* A, double* B, double* C, int64 n);

TARGETSSE void MulSSE(float* A, float* B, float* C, int64 n);

TARGETSSE void MulSSE(double* A, double* B, double C, int64 n);

TARGETSSE void MulSSE(float* A, float* B, float C, int64 n);

TARGETSSE void MulSSE(double* A, double B, int64 n);

TARGETSSE void MulSSE(float* A, float B, int64 n);

TARGETSSE void DivSSE(double* A, double B, double* C, int64 n);

TARGETSSE void DivSSE(float* A, float B, float* C, int64 n);

TARGETSSE void DivSSE(double* A, double* B, double* C, int64 n);

TARGETSSE void DivSSE(float* A, float* B, float* C, int64 n);

TARGETSSE void AddProdSSE(double* A, double* B, double* C, int64 n);

TARGETSSE void AddProdSSE(float* A, float* B, float* C, int64 n);

TARGETSSE void AddProdSSE(double* A, double* B, double C, int64 n);

TARGETSSE void AddProdSSE(double* A, float* B, double C, int64 n);

TARGETSSE void AddProdSSE(float* A, float* B, float C, int64 n);

TARGETSSE void UnifySSE(double* A, int64 n);

TARGETSSE void UnifySSE(float* A, int64 n);

TARGETSSE char* StrNextIdxSSE(char* A, char val, int64 rep, int64 n);

TARGETSSE int64 CountCharSSE(char* A, char val, int64 n);

TARGETSSE float SumSSEx(float* A, int64 n);

TARGETSSE float SumSSEx(float* A, int64 n, int64 sep);

TARGETSSE float ProdSSEx(float* A, int64 n);

TARGETSSE float ProdSSEx(float* A, int64 n, int64 sep);

TARGETSSE float SumProdSSEx(float* A, float* B, int64 sep, int64 n);

TARGETSSE float SumProdSSEx(float* A, float* B, int64 n);

TARGETSSE void DiagQuadFormSSE(double* res, double* A, double* D, int64 m, int64 n);

TARGETSSE void DiagQuadFormSSE(float* res, float* A, float* D, int64 m, int64 n);

TARGETSSE void DiagQuadFormSSE(double* res, double* A, double* D, double* B, int64 m, int64 n);

TARGETSSE void DiagQuadFormSSE(float* res, float* A, float* D, float* B, int64 m, int64 n);

TARGETSSE void DiagQuadFormSSE(double* res, double* A, double* D, int64 n);

TARGETSSE void DiagQuadFormSSE(float* res, float* A, float* D, int64 n);

TARGETSSE void MatrixMulSSE(double* res, double* A, double* B, int64 m, int64 n, int64 p);

TARGETSSE void MatrixMulSSE(float* res, float* A, float* B, int64 m, int64 n, int64 p);

#endif