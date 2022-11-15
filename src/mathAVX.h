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

	/* Draw a uniform distriubted real number */
	template<int nbits>
	TARGETAVX void Poly(__m256d* arr, int n, __m256i* re);

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

	/* Draw a uniform distriubted real number */
	template<int nbits>
	TARGETAVX void Poly(__m256* arr, int n, __m256i* re);
};

#pragma pack(pop)

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

TARGETAVX void SumAVX(double* A, double** B, int64 k, int64 n);

TARGETAVX void SumAVX(float* A, float** B, int64 k, int64 n);

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

TARGETAVX double SumProdAVX(double* A, double* B, int64 sep, int64 n);

TARGETAVX double SumProdAVX(float* A, float* B, int64 sep, int64 n);

TARGETAVX double SumProdAVX(double* A, double* B, int64 n);

TARGETAVX double SumProdAVX(float* A, float* B, int64 n);

TARGETAVX void AddAVX(double* A, double* B, int64 n);

TARGETAVX void AddAVX(float* A, float* B, int64 n);

TARGETAVX void AddAVX(int64* A, int64* B, int64 n);

TARGETAVX void AddAVX(int* A, int* B, int64 n);

TARGETAVX void AddAVX(double* A, double B, int64 n);

TARGETAVX void AddAVX(float* A, float B, int64 n);

TARGETAVX void MulAVX(double* C, double* A, double* B, int64 n);

TARGETAVX void MulAVX(float* C, float* A, float* B, int64 n);

TARGETAVX void MulAVX(double* C, double* A, double B, int64 n);

TARGETAVX void MulAVX(float* C, float* A, float B, int64 n);

TARGETAVX void MulAVX(double* A, double B, int64 n);

TARGETAVX void MulAVX(float* A, float B, int64 n);

TARGETAVX void AddProdAVX(double* C, double* A, double* B, int64 n);

TARGETAVX void AddProdAVX(float* C, float* A, float* B, int64 n);

TARGETAVX void AddProdAVX(double* C, double* A, double B, int64 n);

TARGETAVX void AddProdAVX(double* C, float* A, double B, int64 n);

TARGETAVX void AddProdAVX(float* C, float* A, float B, int64 n);

TARGETAVX void UnifyAVX(double* A, int64 n);

TARGETAVX void UnifyAVX(float* A, int64 n);

TARGETAVX char* StrNextIdxAVX(char* A, char val, int64 rep, int64 n);

TARGETAVX int64 CountCharAVX(char* A, char val, int64 n);

TARGETAVX float LogProdAVXx(float* A, int64 n);

TARGETAVX float LogProdAVXx(float* A, int64 n, int64 sep);

TARGETAVX float LogProdDivAVXx(float* A, float* B, int64 n, int64 sep);

TARGETAVX float SumAVXx(float* A, int64 n);

TARGETAVX float SumAVXx(float* A, int64 n, int64 sep);

TARGETAVX float ProdAVXx(float* A, int64 n);

TARGETAVX float ProdAVXx(float* A, int64 n, int64 sep);

TARGETAVX float SumSquareAVXx(float* A, int64 n);

TARGETAVX float SumProdDivAVXx(float* A1, float* A2, float* B, int64 sep, int64 n);

TARGETAVX float SumProdAVXx(float* A, float* B, int64 sep, int64 n);

TARGETAVX float SumProdAVXx(float* A, float* B, int64 n);

#endif