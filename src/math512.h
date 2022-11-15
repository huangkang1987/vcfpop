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

	/* Draw a uniform distriubted real number */
	template<int nbits>
	TARGET512 void Poly(__m512d* arr, int n, __m512i* re);
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

	/* Draw a uniform distriubted real number */
	template<int nbits>
	TARGET512 void Poly(__m512* arr, int n, __m512i* re);
};

#pragma pack(pop)

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

TARGET512 void Sum512(double* A, double** B, int64 k, int64 n);

TARGET512 void Sum512(float* A, float** B, int64 k, int64 n);

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

TARGET512 double SumProd512(double* A, double* B, int64 sep, int64 n);

TARGET512 double SumProd512(float* A, float* B, int64 sep, int64 n);

TARGET512 double SumProd512(double* A, double* B, int64 n);

TARGET512 double SumProd512(float* A, float* B, int64 n);

TARGET512 void Add512(double* A, double* B, int64 n);

TARGET512 void Add512(float* A, float* B, int64 n);

TARGET512 void Add512(int64* A, int64* B, int64 n);

TARGET512 void Add512(int* A, int* B, int64 n);

TARGET512 void Add512(double* A, double B, int64 n);

TARGET512 void Add512(float* A, float B, int64 n);

TARGET512 void Mul512(double* C, double* A, double* B, int64 n);

TARGET512 void Mul512(float* C, float* A, float* B, int64 n);

TARGET512 void Mul512(double* C, double* A, double B, int64 n);

TARGET512 void Mul512(float* C, float* A, float B, int64 n);

TARGET512 void Mul512(double* A, double B, int64 n);

TARGET512 void Mul512(float* A, float B, int64 n);

TARGET512 void AddProd512(double* C, double* A, double* B, int64 n);

TARGET512 void AddProd512(float* C, float* A, float* B, int64 n);

TARGET512 void AddProd512(double* C, double* A, double B, int64 n);

TARGET512 void AddProd512(double* C, float* A, double B, int64 n);

TARGET512 void AddProd512(float* C, float* A, float B, int64 n);

TARGET512 void Unify512(double* A, int64 n);

TARGET512 void Unify512(float* A, int64 n);

TARGET512 char* StrNextIdx512(char* A, char val, int64 rep, int64 n);

TARGET512 int64 CountChar512(char* A, char val, int64 n);

TARGET512 float LogProd512x(float* A, int64 n);

TARGET512 float LogProd512x(float* A, int64 n, int64 sep);

TARGET512 float LogProdDiv512x(float* A, float* B, int64 n, int64 sep);

TARGET512 float Sum512x(float* A, int64 n);

TARGET512 float Sum512x(float* A, int64 n, int64 sep);

TARGET512 float Prod512x(float* A, int64 n);

TARGET512 float Prod512x(float* A, int64 n, int64 sep);

TARGET512 float SumSquare512x(float* A, int64 n);

TARGET512 float SumProdDiv512x(float* A1, float* A2, float* B, int64 sep, int64 n);

TARGET512 float SumProd512x(float* A, float* B, int64 sep, int64 n);

TARGET512 float SumProd512x(float* A, float* B, int64 n); 

#endif