/* SSE Instruction Set Functions */

#pragma once
#include "vcfpop.h"

#ifndef __aarch64__

#pragma pack(push, 1)

struct RNGSSE
{
public:
	__m128i x[4];
	__m128i y[4];

	/* Initialize rng */
	TARGETSSE RNGSSE();

	/* Initialize rng */
	TARGETSSE RNGSSE(uint64 s);

	/* Draw a uniform distriubted interger */
	TARGETSSE void XorShift128p(__m128i* re);

	/* Draw a uniform distriubted real number */
	TARGETSSE void Uniform(__m128d* re);

	/* Draw a uniform distriubted real number */
	TARGETSSE void Poly(__m128d* a, __m128d* s, int n, __m128i* re);

	/* Draw a polynormial distriubted integer with propoirtions in natural logarithm */
	TARGETSSE void PolyLog(__m128d* a, int n, __m128i* re);
};

#pragma pack(pop)

TARGETSSE int64 GetMinIdxSSE(double* A, int64 n, double& val);

TARGETSSE void GetMinMaxValSSE(double* A, int64 n, double& minv, double& maxv);

TARGETSSE double GetMaxValSSE(double* A, int64 n);

TARGETSSE double GetMaxValSSE(double* A, int64 n, int64 sep);

TARGETSSE double GetMinValSSE(double* A, int64 n);

TARGETSSE int64 GetMinValSSE(int64* A, int64 n);

TARGETSSE void SetValSSE(uint* a, ushort* b, int64 n);

TARGETSSE void SetValSSE(uint* A, ushort* B, int64 n);

TARGETSSE void AddExponentSSE(int64& slog, __m128d& val);

TARGETSSE void ChargeLogSSE(int64& slog, double& prod, __m128d& val);

TARGETSSE double LogProdSSE(double* A, int64 n);

TARGETSSE double LogProdSSE(double* A, int64 n, int64 sep);

TARGETSSE double LogProdDivSSE(double* A, double* B, int64 n, int64 sep);

TARGETSSE int64 CountNonZeroSSE(byte* A, int64 n);

TARGETSSE double SumSSE(double* A, int64 n);

TARGETSSE int64 SumSSE(byte* A, int64 n);

TARGETSSE double SumSSE(double* A, int64 n, int64 sep);

TARGETSSE void SumSSE(double* A, double** B, int64 k, int64 n);

TARGETSSE double ProdSSE(double* A, int64 n);

TARGETSSE double ProdSSE(double* A, int64 n, int64 sep);

TARGETSSE double SumSquareSSE(double* A, int64 n);

TARGETSSE int64 SumSquareSSE(byte* A, int64 n);

TARGETSSE void SumSumSquareSSE(double* A, int64 n, double& sum, double& sumsq);

TARGETSSE double SumProdDivSSE(double* A1, double* A2, double* B, int64 sep, int64 n);

TARGETSSE double SumProdSSE(double* A, double* B, int64 sep, int64 n);

TARGETSSE double SumProdSSE(double* A, double* B, int64 n);

TARGETSSE void AddSSE(double* A, double* B, int64 n);

TARGETSSE void AddSSE(int64* A, int64* B, int64 n);

TARGETSSE void AddSSE(int* A, int* B, int64 n);

TARGETSSE void AddSSE(double* A, double B, int64 n);

TARGETSSE void MulSSE(double* C, double* A, double* B, int64 n);

TARGETSSE void MulSSE(double* C, double* A, double B, int64 n);

TARGETSSE void MulSSE(double* A, double B, int64 n);

TARGETSSE void AddProdSSE(double* C, double* A, double* B, int64 n);

TARGETSSE void AddProdSSE(double* C, double* A, double B, int64 n);

TARGETSSE void UnifySSE(double* A, int64 n);

TARGETSSE char* StrNextIdxSSE(char* A, char val, int64 rep, int64 n);

TARGETSSE int64 CountCharSSE(char* A, char val, int64 n);

#endif