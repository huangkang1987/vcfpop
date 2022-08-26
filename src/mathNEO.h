/* NEO Instruction Set Functions */

#pragma once
#include "vcfpop.h"

#ifdef __aarch64__

#pragma pack(push, 1)

struct RNGNEO
{
public:
	uint64x2_t x[4];
	uint64x2_t y[4];

	/* Initialize rng */
	TARGETNEO RNGNEO();

	/* Initialize rng */
	TARGETNEO RNGNEO(uint64 s);

	/* Draw a uniform distriubted interger */
	TARGETNEO void XorShift128p(uint64x2_t* re);

	/* Draw a uniform distriubted real number */
	TARGETNEO void Uniform(float64x2_t* re);

	/* Draw a uniform distriubted real number */
	TARGETNEO void Poly(float64x2_t* a, float64x2_t* s, int n, uint64x2_t* re);

	/* Draw a polynormial distriubted integer with propoirtions in natural logarithm */
	TARGETNEO void PolyLog(float64x2_t* a, int n, uint64x2_t* re);
};

#pragma pack(pop)

TARGETNEO int64 GetMinIdxNEO(double* A, int64 n, double& val);

TARGETNEO void GetMinMaxValNEO(double* A, int64 n, double& minv, double& maxv);

TARGETNEO double GetMaxValNEO(double* A, int64 n);

TARGETNEO double GetMaxValNEO(double* A, int64 n, int64 sep);

TARGETNEO double GetMinValNEO(double* A, int64 n);

TARGETNEO int64 GetMinValNEO(int64* A, int64 n);

TARGETNEO void SetValNEO(uint* a, ushort* b, int64 n);

TARGETNEO void SetValNEO(uint* A, ushort* B, int64 n);

TARGETNEO void AddExponentNEO(int64& slog, float64x2_t& val);

TARGETNEO void ChargeLogNEO(int64& slog, double& prod, float64x2_t& val);

TARGETNEO double LogProdNEO(double* A, int64 n);

TARGETNEO double LogProdNEO(double* A, int64 n, int64 sep);

TARGETNEO double LogProdDivNEO(double* A, double* B, int64 n, int64 sep);

TARGETNEO int64 CountNonZeroNEO(byte* A, int64 n);

TARGETNEO double SumNEO(double* A, int64 n);

TARGETNEO int64 SumNEO(byte* A, int64 n);

TARGETNEO double SumNEO(double* A, int64 n, int64 sep);

TARGETNEO void SumNEO(double* A, double** B, int64 k, int64 n);

TARGETNEO double ProdNEO(double* A, int64 n);

TARGETNEO double ProdNEO(double* A, int64 n, int64 sep);

TARGETNEO double SumSquareNEO(double* A, int64 n);

TARGETNEO int64 SumSquareNEO(byte* A, int64 n);

TARGETNEO void SumSumSquareNEO(double* A, int64 n, double& sum, double& sumsq);

TARGETNEO double SumProdDivNEO(double* A1, double* A2, double* B, int64 sep, int64 n);

TARGETNEO double SumProdNEO(double* A, double* B, int64 sep, int64 n);

TARGETNEO double SumProdNEO(double* A, double* B, int64 n);

TARGETNEO void AddNEO(double* A, double* B, int64 n);

TARGETNEO void AddNEO(int64* A, int64* B, int64 n);

TARGETNEO void AddNEO(int* A, int* B, int64 n);

TARGETNEO void AddNEO(double* A, double B, int64 n);

TARGETNEO void MulNEO(double* C, double* A, double* B, int64 n);

TARGETNEO void MulNEO(double* C, double* A, double B, int64 n);

TARGETNEO void MulNEO(double* A, double B, int64 n);

TARGETNEO void AddProdNEO(double* C, double* A, double* B, int64 n);

TARGETNEO void AddProdNEO(double* C, double* A, double B, int64 n);

TARGETNEO void UnifyNEO(double* A, int64 n);

TARGETNEO char* StrNextIdxNEO(char* A, char val, int64 rep, int64 n);

TARGETNEO int64 CountCharNEO(char* A, char val, int64 n);

#endif