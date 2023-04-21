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

	/* Draw a uniform distriubted real number */
	template<int nbits>
	TARGETSSE void Poly(__m128d* arr, int n, __m128i* re);
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

	/* Draw a uniform distriubted real number */
	template<int nbits>
	TARGETSSE void Poly(__m128* arr, int n, __m128i* re);
};

#pragma pack(pop)

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

TARGETSSE void SumSSE(double* A, double** B, int64 k, int64 n);

TARGETSSE void SumSSE(float* A, float** B, int64 k, int64 n);

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

TARGETSSE double SumProdSSE(double* A, double* B, int64 sep, int64 n);

TARGETSSE double SumProdSSE(float* A, float* B, int64 sep, int64 n);

TARGETSSE double SumProdSSE(double* A, double* B, int64 n);

TARGETSSE double SumProdSSE(float* A, float* B, int64 n);

TARGETSSE void AddSSE(double* A, double* B, int64 n);

TARGETSSE void AddSSE(float* A, float* B, int64 n);

TARGETSSE void AddSSE(int64* A, int64* B, int64 n);

TARGETSSE void AddSSE(int* A, int* B, int64 n);

TARGETSSE void AddSSE(int* A, int B, int64 n);

TARGETSSE void AddSSE(double* A, double B, int64 n);

TARGETSSE void AddSSE(float* A, float B, int64 n);

TARGETSSE void MulSSE(double* C, double* A, double* B, int64 n);

TARGETSSE void MulSSE(float* C, float* A, float* B, int64 n);

TARGETSSE void MulSSE(double* C, double* A, double B, int64 n);

TARGETSSE void MulSSE(float* C, float* A, float B, int64 n);

TARGETSSE void MulSSE(double* A, double B, int64 n);

TARGETSSE void MulSSE(float* A, float B, int64 n);

TARGETSSE void AddProdSSE(double* C, double* A, double* B, int64 n);

TARGETSSE void AddProdSSE(float* C, float* A, float* B, int64 n);

TARGETSSE void AddProdSSE(double* C, double* A, double B, int64 n);

TARGETSSE void AddProdSSE(double* C, float* A, double B, int64 n);

TARGETSSE void AddProdSSE(float* C, float* A, float B, int64 n);

TARGETSSE void UnifySSE(double* A, int64 n);

TARGETSSE void UnifySSE(float* A, int64 n);

TARGETSSE char* StrNextIdxSSE(char* A, char val, int64 rep, int64 n);

TARGETSSE int64 CountCharSSE(char* A, char val, int64 n);

TARGETSSE float LogProdSSEx(float* A, int64 n);

TARGETSSE float LogProdSSEx(float* A, int64 n, int64 sep);

TARGETSSE float LogProdDivSSEx(float* A, float* B, int64 n, int64 sep);

TARGETSSE float SumSSEx(float* A, int64 n);

TARGETSSE float SumSSEx(float* A, int64 n, int64 sep);

TARGETSSE float ProdSSEx(float* A, int64 n);

TARGETSSE float ProdSSEx(float* A, int64 n, int64 sep);

TARGETSSE float SumSquareSSEx(float* A, int64 n);

TARGETSSE float SumProdDivSSEx(float* A1, float* A2, float* B, int64 sep, int64 n);

TARGETSSE float SumProdSSEx(float* A, float* B, int64 sep, int64 n);

TARGETSSE float SumProdSSEx(float* A, float* B, int64 n);

#endif