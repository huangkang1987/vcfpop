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

	/* Draw 64 64-bit integers in [0,n), 64*n frequencies are in arr */
	TARGETNEO void Poly(float64x2_t* arr, int n, void* re);

	/* Draw uniform distriubted intergers */
	TARGETNEO void XorShift();

	/* Draw uniform distriubted intergers */
	template<typename INT>
	TARGETNEO void Integer(INT* re, int64 n, INT minv = 0, INT maxv = (INT)-1);

	/* Draw uniform distriubted real numbers */
	TARGETNEO void Uniform(double* re, int n, double minv = 0, double maxv = 1);

	/* Draw normal distriubted real numbers */
	TARGETNEO void Normal(double* re, int n, double mean = 0, double sd = 1);
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

	/* Draw 64 64-bit integers in [0,n), 64*n frequencies are in arr */
	TARGETNEO void Poly(float32x4_t* arr, int n, void* re);

	/* Draw uniform distriubted intergers */
	TARGETNEO void XorShift();

	/* Draw uniform distriubted intergers */
	template<typename INT>
	TARGETNEO void Integer(INT* re, int64 n, INT minv = 0, INT maxv = (INT)-1);

	/* Draw uniform distriubted real numbers */
	TARGETNEO void Uniform(float* re, int n, float minv = 0, float maxv = 1);

	/* Draw normal distriubted real numbers */
	TARGETNEO void Normal(float* re, int n, float mean = 0, float sd = 1);
};

#pragma pack(pop)

static forceinline TARGETNEO int64 _neo_reduce_add_epi64(uint64x2_t v1)
{
	return vgetq_lane_s64(v1, 0) + vgetq_lane_s64(v1, 1);
}

static forceinline TARGETNEO double _neo_reduce_add_pd(float64x2_t v1)
{
	return vgetq_lane_f64(v1, 0) + vgetq_lane_f64(v1, 1);
}

static forceinline TARGETNEO double _neo_reduce_mul_pd(float64x2_t v1)
{
	return vgetq_lane_f64(v1, 0) * vgetq_lane_f64(v1, 1);
}

static forceinline TARGETNEO float _neo_reduce_add_ps(float32x4_t v1)
{
	volatile float a1 = vgetq_lane_f32(v1, 0) + vgetq_lane_f32(v1, 2), a2 = vgetq_lane_f32(v1, 1) + vgetq_lane_f32(v1, 3);
	return a1 + a2;
}

static forceinline TARGETNEO float _neo_reduce_mul_ps(float32x4_t v1)
{
	volatile float a1 = vgetq_lane_f32(v1, 0) * vgetq_lane_f32(v1, 2), a2 = vgetq_lane_f32(v1, 1) * vgetq_lane_f32(v1, 3);
	return a1 * a2;
}

static forceinline TARGETNEO double _neo_reduce_add_psd(float32x4_t v1)
{
	volatile double a1 = (double)vgetq_lane_f32(v1, 0) + (double)vgetq_lane_f32(v1, 2);
	volatile double a2 = (double)vgetq_lane_f32(v1, 1) + (double)vgetq_lane_f32(v1, 3);
	return a1 + a2;
}

static forceinline TARGETNEO double _neo_reduce_mul_psd(float32x4_t v1)
{
	volatile double a1 = (double)vgetq_lane_f32(v1, 0) * (double)vgetq_lane_f32(v1, 2);
	volatile double a2 = (double)vgetq_lane_f32(v1, 1) * (double)vgetq_lane_f32(v1, 3);
	return a1 * a2;
}

static forceinline TARGETNEO int64 _neo_reduce_min_epi64(uint64x2_t v1)
{
    return std::min(vgetq_lane_s64(v1, 0), vgetq_lane_s64(v1, 1));
}

static forceinline TARGETNEO int _neo_reduce_min_epi32(int32x4_t v1)
{
	return std::min(std::min(vgetq_lane_s32(v1, 0), vgetq_lane_s32(v1, 1)),
			        std::min(vgetq_lane_s32(v1, 2), vgetq_lane_s32(v1, 3)));
}

static forceinline TARGETNEO double _neo_reduce_min_pd(float64x2_t v1) 
{
	return std::min(vgetq_lane_f64(v1, 0), vgetq_lane_f64(v1, 1));
}

static forceinline TARGETNEO float _neo_reduce_min_ps(float32x4_t v1)
{
	return std::min(std::min(vgetq_lane_f32(v1, 0), vgetq_lane_f32(v1, 1)),
			        std::min(vgetq_lane_f32(v1, 2), vgetq_lane_f32(v1, 3)));
}

static forceinline TARGETNEO int64 _neo_reduce_max_epi64(uint64x2_t v1)
{
    return std::max(vgetq_lane_s64(v1, 0), vgetq_lane_s64(v1, 1));
}

static forceinline TARGETNEO int _neo_reduce_max_epi32(int32x4_t v1)
{
	return std::max(std::max(vgetq_lane_s32(v1, 0), vgetq_lane_s32(v1, 1)),
			        std::max(vgetq_lane_s32(v1, 2), vgetq_lane_s32(v1, 3)));
}

static forceinline TARGETNEO double _neo_reduce_max_pd(float64x2_t v1)
{
	return std::max(vgetq_lane_f64(v1, 0), vgetq_lane_f64(v1, 1));
}

static forceinline TARGETNEO float _neo_reduce_max_ps(float32x4_t v1)
{
	return std::max(std::max(vgetq_lane_f32(v1, 0), vgetq_lane_f32(v1, 1)),
			        std::max(vgetq_lane_f32(v1, 2), vgetq_lane_f32(v1, 3)));
}

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

TARGETNEO float SumProdDivNEOx(float* A1, float* A2, float* B, int64 sep, int64 n);

TARGETNEO double SumProdNEO(double* A, double* B, int64 sep, int64 n);

TARGETNEO double SumProdNEO(float* A, float* B, int64 sep, int64 n);

TARGETNEO double SumProdNEO(double* A, double* B, int64 n);

TARGETNEO double SumProdNEO(double* A, double* B, double* C, int64 n);

TARGETNEO double SumProdNEO(float* A, float* B, int64 n);

TARGETNEO float SumProdNEO(float* A, float* B, float* C, int64 n);

TARGETNEO double SumSqProdNEO(double* A, double* B, int64 n);

TARGETNEO float SumSqProdNEO(float* A, float* B, int64 n);

TARGETNEO void AddNEO(double* A, double* B, int64 n);

TARGETNEO void AddNEO(float* A, float* B, int64 n);

TARGETNEO void AddNEO(int64* A, int64* B, int64 n);

TARGETNEO void AddNEO(int* A, int* B, int64 n);

TARGETNEO void AddNEO(int* A, int B, int64 n);

TARGETNEO void AddNEO(double* A, double B, int64 n);

TARGETNEO void AddNEO(float* A, float B, int64 n);

TARGETNEO void MulNEO(double* A, double* B, double* C, int64 n);

TARGETNEO void MulNEO(float* A, float* B, float* C, int64 n);

TARGETNEO void MulNEO(double* A, double* B, double C, int64 n);

TARGETNEO void MulNEO(float* A, float* B, float C, int64 n);

TARGETNEO void MulNEO(double* A, double B, int64 n);

TARGETNEO void MulNEO(float* A, float B, int64 n);

TARGETNEO void DivNEO(double* A, double B, double* C, int64 n);

TARGETNEO void DivNEO(float* A, float B, float* C, int64 n);

TARGETNEO void DivNEO(double* A, double* B, double* C, int64 n);

TARGETNEO void DivNEO(float* A, float* B, float* C, int64 n);

TARGETNEO void AddProdNEO(double* A, double* B, double* C, int64 n);

TARGETNEO void AddProdNEO(float* A, float* B, float* C, int64 n);

TARGETNEO void AddProdNEO(double* A, double* B, double C, int64 n);

TARGETNEO void AddProdNEO(double* A, float* B, double C, int64 n);

TARGETNEO void AddProdNEO(float* A, float* B, float C, int64 n);

TARGETNEO void UnifyNEO(double* A, int64 n);

TARGETNEO void UnifyNEO(float* A, int64 n);

TARGETNEO char* StrNextIdxNEO(char* A, char val, int64 rep, int64 n);

TARGETNEO int64 CountCharNEO(char* A, char val, int64 n);

TARGETNEO float SumNEOx(float* A, int64 n);

TARGETNEO float SumNEOx(float* A, int64 n, int64 sep);

TARGETNEO float ProdNEOx(float* A, int64 n);

TARGETNEO float ProdNEOx(float* A, int64 n, int64 sep);

TARGETNEO float SumProdNEOx(float* A, float* B, int64 sep, int64 n);

TARGETNEO float SumProdNEOx(float* A, float* B, int64 n);

TARGETNEO void DiagQuadFormNEO(double* res, double* A, double* D, int64 m, int64 n);

TARGETNEO void DiagQuadFormNEO(float* res, float* A, float* D, int64 m, int64 n);

TARGETNEO void DiagQuadFormNEO(double* res, double* A, double* D, double* B, int64 m, int64 n);

TARGETNEO void DiagQuadFormNEO(float* res, float* A, float* D, float* B, int64 m, int64 n);

TARGETNEO void DiagQuadFormNEO(double* res, double* A, double* D, int64 n);

TARGETNEO void DiagQuadFormNEO(float* res, float* A, float* D, int64 n);

TARGETNEO void MatrixMulNEO(double* res, double* A, double* B, int64 m, int64 n, int64 p);

TARGETNEO void MatrixMulNEO(float* res, float* A, float* B, int64 m, int64 n, int64 p);

#endif
