/* AVX Instruction Set Functions */

#pragma once
#include "vcfpop.h"

#ifndef __aarch64__

#pragma pack(push, 1)

struct RNGAVX
{
public:
	__m256i x[2];
	__m256i y[2];

	/* Initialize rng */
	TARGETAVX RNGAVX();

	/* Initialize rng */
	TARGETAVX RNGAVX(uint64 s);

	/* Draw a uniform distriubted interger */
	TARGETAVX void XorShift128p(__m256i* re);

	/* Draw a uniform distriubted real number */
	TARGETAVX void Uniform(__m256d* re);

	/* Draw a uniform distriubted real number */
	TARGETAVX void Poly(__m256d* a, __m256d* s, int n, __m256i* re);

	TARGETAVX void PolyLog(__m256d* a, int n, __m256i* re);
};

#pragma pack(pop)

TARGETAVX int64 GetMinIdxAVX(double* A, int64 n, double& val);

TARGETAVX void GetMinMaxValAVX(double* A, int64 n, double& minv, double& maxv);

TARGETAVX double GetMaxValAVX(double* A, int64 n);

TARGETAVX double GetMaxValAVX(double* A, int64 n, int64 sep);

TARGETAVX double GetMinValAVX(double* A, int64 n);

TARGETAVX int64 GetMinValAVX(int64* A, int64 n);

TARGETAVX void SetValAVX(uint* a, ushort* b, int64 n);

TARGETAVX void AddExponentAVX(int64& slog, __m256d& val);

TARGETAVX void ChargeLogAVX(int64& slog, double& prod, __m256d& val);

TARGETAVX double LogProdAVX(double* A, int64 n);

TARGETAVX double LogProdAVX(double* A, int64 n, int64 sep);

TARGETAVX double LogProdDivAVX(double* A, double* B, int64 n, int64 sep);

TARGETAVX int64 CountNonZeroAVX(byte* A, int64 n);

TARGETAVX double SumAVX(double* A, int64 n);

TARGETAVX int64 SumAVX(byte* A, int64 n);

TARGETAVX double SumAVX(double* A, int64 n, int64 sep);

TARGETAVX void SumAVX(double* A, double** B, int64 k, int64 n);

TARGETAVX double ProdAVX(double* A, int64 n);

TARGETAVX double ProdAVX(double* A, int64 n, int64 sep);

TARGETAVX double SumSquareAVX(double* A, int64 n);

TARGETAVX int64 SumSquareAVX(byte* A, int64 n);

TARGETAVX void SumSumSquareAVX(double* A, int64 n, double& sum, double& sumsq);

TARGETAVX double SumProdDivAVX(double* A1, double* A2, double* B, int64 sep, int64 n);

TARGETAVX double SumProdAVX(double* A, double* B, int64 sep, int64 n);

TARGETAVX double SumProdAVX(double* A, double* B, int64 n);

TARGETAVX void AddAVX(double* A, double* B, int64 n);

TARGETAVX void AddAVX(int64* A, int64* B, int64 n);

TARGETAVX void AddAVX(int* A, int* B, int64 n);

TARGETAVX void AddAVX(double* A, double B, int64 n);

TARGETAVX void MulAVX(double* C, double* A, double* B, int64 n);

TARGETAVX void MulAVX(double* C, double* A, double B, int64 n);

TARGETAVX void MulAVX(double* A, double B, int64 n);

TARGETAVX void AddProdAVX(double* C, double* A, double* B, int64 n);

TARGETAVX void AddProdAVX(double* C, double* A, double B, int64 n);

TARGETAVX void UnifyAVX(double* A, int64 n);

TARGETAVX char* StrNextIdxAVX(char* A, char val, int64 rep, int64 n);

TARGETAVX int64 CountCharAVX(char* A, char val, int64 n);

#endif