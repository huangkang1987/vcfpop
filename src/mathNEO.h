/* NEO Instruction Set Functions */

#pragma once
#include "vcfpop.h"

#ifdef __aarch64__

#pragma pack(push, 1)

template<typename REAL>
struct RNGNEO
{

};

template<>
struct RNGNEO<double>
{
	uint64x2_t x[32];
	uint64x2_t y[32];

	/* Initialize rng */
	TARGETNEO RNGNEO();

	/* Initialize rng */
	TARGETNEO RNGNEO(uint64 s, uint64 salt);

	/* Draw a uniform distriubted real number */
	template<int nbits>
	TARGETNEO void Poly(float64x2_t* arr, int n, void* re);
};

template<>
struct RNGNEO<float>
{
	uint32x4_t x[16];
	uint32x4_t y[16];
	uint32x4_t z[16];

	/* Initialize rng */
	TARGETNEO RNGNEO();

	/* Initialize rng */
	TARGETNEO RNGNEO(uint64 s, uint64 salt);

	/* Draw a uniform distriubted real number */
	template<int nbits>
	TARGETNEO void Poly(float32x4_t* arr, int n, void* re);
};

#pragma pack(pop)

TARGETNEO int64 GetMinIdxNEO(double* A, int64 n, double& val);

TARGETNEO int64 GetMinIdxNEO(float* A, int64 n, float& val);

TARGETNEO void GetMinMaxValNEO(double* A, int64 n, double& minv, double& maxv);

TARGETNEO void GetMinMaxValNEO(float* A, int64 n, float& minv, float& maxv);

TARGETNEO double GetMaxValNEO(double* A, int64 n);

TARGETNEO float GetMaxValNEO(float* A, int64 n);

TARGETNEO double GetMaxValNEO(double* A, int64 n, int64 sep);

TARGETNEO float GetMaxValNEO(float* A, int64 n, int64 sep);

TARGETNEO double GetMinValNEO(double* A, int64 n);

TARGETNEO float GetMinValNEO(float* A, int64 n);

TARGETNEO int64 GetMinValNEO(int64* A, int64 n);

TARGETNEO void SetValNEO(uint* a, ushort* b, int64 n);

TARGETNEO void SetValNEO(uint* A, ushort* B, int64 n);

TARGETNEO void AddExponentNEO(int64& slog, float64x2_t& val);

TARGETNEO void AddExponentNEO(int64& slog, float32x4_t& val);

TARGETNEO void ChargeLogNEO(int64& slog, double& prod, float64x2_t& val);

TARGETNEO void ChargeLogNEO(int64& slog, double& prod, float32x4_t& val);

TARGETNEO double LogProdNEO(double* A, int64 n);

TARGETNEO double LogProdNEO(float* A, int64 n);

TARGETNEO double LogProdNEO(double* A, int64 n, int64 sep);

TARGETNEO double LogProdNEO(float* A, int64 n, int64 sep);

TARGETNEO double LogProdDivNEO(double* A, double* B, int64 n, int64 sep);

TARGETNEO double LogProdDivNEO(float* A, float* B, int64 n, int64 sep);

TARGETNEO int64 CountNonZeroNEO(byte* A, int64 n);

TARGETNEO double SumNEO(double* A, int64 n);

TARGETNEO double SumNEO(float* A, int64 n);

TARGETNEO int64 SumNEO(byte* A, int64 n);

TARGETNEO double SumNEO(double* A, int64 n, int64 sep);

TARGETNEO double SumNEO(float* A, int64 n, int64 sep);

TARGETNEO void SumNEO(double* A, double** B, int64 k, int64 n);

TARGETNEO void SumNEO(float* A, float** B, int64 k, int64 n);

TARGETNEO double ProdNEO(double* A, int64 n);

TARGETNEO double ProdNEO(float* A, int64 n);

TARGETNEO double ProdNEO(double* A, int64 n, int64 sep);

TARGETNEO double ProdNEO(float* A, int64 n, int64 sep);

TARGETNEO double SumSquareNEO(double* A, int64 n);

TARGETNEO double SumSquareNEO(float* A, int64 n);

TARGETNEO int64 SumSquareNEO(byte* A, int64 n);

TARGETNEO void SumSumSquareNEO(double* A, int64 n, double& sum, double& sumsq);

TARGETNEO void SumSumSquareNEO(float* A, int64 n, double& sum, double& sumsq);

TARGETNEO double SumProdDivNEO(double* A1, double* A2, double* B, int64 sep, int64 n);

TARGETNEO double SumProdDivNEO(double* A1, float* A2, float* B, int64 sep, int64 n);

TARGETNEO double SumProdDivNEO(float* A1, float* A2, float* B, int64 sep, int64 n);

TARGETNEO double SumProdNEO(double* A, double* B, int64 sep, int64 n);

TARGETNEO double SumProdNEO(float* A, float* B, int64 sep, int64 n);

TARGETNEO double SumProdNEO(double* A, double* B, int64 n);

TARGETNEO double SumProdNEO(float* A, float* B, int64 n);

TARGETNEO void AddNEO(double* A, double* B, int64 n);

TARGETNEO void AddNEO(float* A, float* B, int64 n);

TARGETNEO void AddNEO(int64* A, int64* B, int64 n);

TARGETNEO void AddNEO(int* A, int* B, int64 n);

TARGETNEO void AddNEO(int* A, int B, int64 n);

TARGETNEO void AddNEO(double* A, double B, int64 n);

TARGETNEO void AddNEO(float* A, float B, int64 n);

TARGETNEO void MulNEO(double* C, double* A, double* B, int64 n);

TARGETNEO void MulNEO(float* C, float* A, float* B, int64 n);

TARGETNEO void MulNEO(double* C, double* A, double B, int64 n);

TARGETNEO void MulNEO(float* C, float* A, float B, int64 n);

TARGETNEO void MulNEO(double* A, double B, int64 n);

TARGETNEO void MulNEO(float* A, float B, int64 n);

TARGETNEO void AddProdNEO(double* C, double* A, double* B, int64 n);

TARGETNEO void AddProdNEO(float* C, float* A, float* B, int64 n);

TARGETNEO void AddProdNEO(double* C, double* A, double B, int64 n);

TARGETNEO void AddProdNEO(double* C, float* A, double B, int64 n);

TARGETNEO void AddProdNEO(float* C, float* A, float B, int64 n);

TARGETNEO void UnifyNEO(double* A, int64 n);

TARGETNEO void UnifyNEO(float* A, int64 n);

TARGETNEO char* StrNextIdxNEO(char* A, char val, int64 rep, int64 n);

TARGETNEO int64 CountCharNEO(char* A, char val, int64 n);

TARGETNEO float LogProdNEOx(float* A, int64 n);

TARGETNEO float LogProdNEOx(float* A, int64 n, int64 sep);

TARGETNEO float LogProdDivNEOx(float* A, float* B, int64 n, int64 sep);

TARGETNEO float SumNEOx(float* A, int64 n);

TARGETNEO float SumNEOx(float* A, int64 n, int64 sep);

TARGETNEO float ProdNEOx(float* A, int64 n);

TARGETNEO float ProdNEOx(float* A, int64 n, int64 sep);

TARGETNEO float SumSquareNEOx(float* A, int64 n);

TARGETNEO float SumProdDivNEOx(float* A1, float* A2, float* B, int64 sep, int64 n);

TARGETNEO float SumProdNEOx(float* A, float* B, int64 sep, int64 n);

TARGETNEO float SumProdNEOx(float* A, float* B, int64 n);

#endif