/* Math  Functions */

#include "vcfpop.h"

/* Get 1D lower triangular index */
TARGET int GetLowerTriangularId(int id1, int id2, int n)
{
	int minid = Min(id1, id2);
	return Max(id1, id2) - ((minid * (1 + minid - n - n)) >> 1);
}

/* Is an erroneous real number */
TARGET bool IsError(double x)
{
	return isnan(x) || isinf(x);
}

/* Is an erroneous real number */
TARGET bool IsError(float x)
{
	return isnan(x) || isinf(x);
}

/* Is a normal real number */
TARGET bool IsNormal(double x)
{
	return !IsError(x);
}

/* Is a normal real number */
TARGET bool IsNormal(float x)
{
	return !IsError(x);
}

/* Find the index of the mimumum element */
TARGET int64 GetMinIdx(double* A, int64 n, double& val)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return GetMinIdxNEO(A, n, val);
#else
	case 4: return GetMinIdx512(A, n, val);
	case 3: return GetMinIdxAVX(A, n, val);
	case 2: return GetMinIdxSSE(A, n, val);
#endif
	}

	val = DBL_MAX;
	int64 idx = -1;
	for (int64 i = 0; i < n; ++i)
	{
		if (A[i] > val) continue;
		val = A[i];
		idx = i;
	}
	return idx;
}

/* Find the index of the mimumum element */
TARGET int64 GetMinIdx(float* A, int64 n, float& val)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return GetMinIdxNEO(A, n, val);
#else
	case 4: return GetMinIdx512(A, n, val);
	case 3: return GetMinIdxAVX(A, n, val);
	case 2: return GetMinIdxSSE(A, n, val);
#endif
	}

	val = FLT_MAX;
	int64 idx = -1;
	for (int64 i = 0; i < n; ++i)
	{
		if (A[i] > val) continue;
		val = A[i];
		idx = i;
	}
	return idx;
}

/* Find maximum and minimum element of A */
TARGET void GetMinMaxVal(double* A, int64 n, double& minv, double& maxv)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return GetMinMaxValNEO(A, n, minv, maxv);
#else
	case 4: return GetMinMaxVal512(A, n, minv, maxv);
	case 3: return GetMinMaxValAVX(A, n, minv, maxv);
	case 2: return GetMinMaxValSSE(A, n, minv, maxv);
#endif
	}
	
	minv = DBL_MAX;
	maxv = -DBL_MAX;
	for (int64 i = 0; i < n; ++i)
	{
		if (A[i] < minv) minv = A[i];
		if (A[i] > maxv) maxv = A[i];
	}
}

/* Find maximum and minimum element of A */
TARGET void GetMinMaxVal(float* A, int64 n, float& minv, float& maxv)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return GetMinMaxValNEO(A, n, minv, maxv);
#else
	case 4: return GetMinMaxVal512(A, n, minv, maxv);
	case 3: return GetMinMaxValAVX(A, n, minv, maxv);
	case 2: return GetMinMaxValSSE(A, n, minv, maxv);
#endif
	}
	
	minv = FLT_MAX;
	maxv = -FLT_MAX;
	for (int64 i = 0; i < n; ++i)
	{
		if (A[i] < minv) minv = A[i];
		if (A[i] > maxv) maxv = A[i];
	}
}

/* Find minimum element of A */
TARGET double GetMinVal(double* A, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return GetMinValNEO(A, n);
#else
	case 4: return GetMinVal512(A, n);
	case 3: return GetMinValAVX(A, n);
	case 2: return GetMinValSSE(A, n);
#endif
	}
	
	double val = DBL_MAX;
	for (int64 i = 0; i < n; ++i)
	{
		if (A[i] > val) continue;
		val = A[i];
	}
	return val;
}

/* Find minimum element of A */
TARGET float GetMinVal(float* A, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return GetMinValNEO(A, n);
#else
	case 4: return GetMinVal512(A, n);
	case 3: return GetMinValAVX(A, n);
	case 2: return GetMinValSSE(A, n);
#endif
	}
	
	float val = FLT_MAX;
	for (int64 i = 0; i < n; ++i)
	{
		if (A[i] > val) continue;
		val = A[i];
	}
	return val;
}

/* Find minimum element of A */
TARGET int64 GetMinVal(int64* A, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return GetMinValNEO(A, n);
#else
	case 4: return GetMinVal512(A, n);
	case 3: return GetMinValAVX(A, n);
	case 2: return GetMinValSSE(A, n);
#endif
	}

	int64 val = 0x7FFFFFFFFFFFFFFF;
	for (int64 i = 0; i < n; ++i)
	{
		if (A[i] > val) continue;
		val = A[i];
	}
	return val;
}

/* Find Maximum element of A */
TARGET double GetMaxVal(double* A, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return GetMaxValNEO(A, n);
#else
	case 4: return GetMaxVal512(A, n);
	case 3: return GetMaxValAVX(A, n);
	case 2: return GetMaxValSSE(A, n);
#endif
	}
	
	double val = -DBL_MAX;
	for (int64 i = 0; i < n; ++i)
	{
		if (A[i] < val) continue;
		val = A[i];
	}
	return val;
}

/* Find Maximum element of A */
TARGET float GetMaxVal(float* A, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return GetMaxValNEO(A, n);
#else
	case 4: return GetMaxVal512(A, n);
	case 3: return GetMaxValAVX(A, n);
	case 2: return GetMaxValSSE(A, n);
#endif
	}
	
	float val = -FLT_MAX;
	for (int64 i = 0; i < n; ++i)
	{
		if (A[i] < val) continue;
		val = A[i];
	}
	return val;
}

/* Find Maximum element of A */
TARGET double GetMaxVal(double* A, int64 n, int64 sep)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return GetMaxValNEO(A, n, sep);
#else
	case 4: return GetMaxVal512(A, n, sep);
	case 3: return GetMaxValAVX(A, n, sep);
	case 2: return GetMaxValSSE(A, n, sep);
#endif
	}

	double val = -DBL_MAX;
	for (int64 i = 0; i < n; ++i, A += sep)
	{
		if (*A < val) continue;
		val = *A;
	}
	return val;
}

/* Find Maximum element of A */
TARGET float GetMaxVal(float* A, int64 n, int64 sep)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return GetMaxValNEO(A, n, sep);
#else
	case 4: return GetMaxVal512(A, n, sep);
	case 3: return GetMaxValAVX(A, n, sep);
	case 2: return GetMaxValSSE(A, n, sep);
#endif
	}

	float val = -FLT_MAX;
	for (int64 i = 0; i < n; ++i, A += sep)
	{
		if (*A < val) continue;
		val = *A;
	}
	return val;
}

/* A[i] = B[i] */
TARGET void SetVal(uint* A, ushort* B, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SetValNEO(A, B, n);
#else
	case 4: return SetVal512(A, B, n);
	case 3: return SetValAVX(A, B, n);
	case 2: return SetValSSE(A, B, n);
#endif
	}
	
	for (int64 i = 0; i < n; ++i)
		A[i] = B[i];
}

/* log(prod(A[i++])) */
TARGET double LogProd(double* A, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return LogProdNEO(A, n);
#else
	case 4: return LogProd512(A, n);
	case 3: return LogProdAVX(A, n);
	case 2: return LogProdSSE(A, n);
#endif
	}

	int64 slog = 0; double prod = 1;
	for (int64 i = 0; i < n; ++i)
		ChargeLog(slog, prod, A[i]);

	CloseLog(slog, prod);
	return prod;
}

/* log(prod(A[i++])) */
TARGET double LogProd(float* A, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return LogProdNEO(A, n);
#else
	case 4: return LogProd512(A, n);
	case 3: return LogProdAVX(A, n);
	case 2: return LogProdSSE(A, n);
#endif
	}

	int64 slog = 0; double prod = 1;
	for (int64 i = 0; i < n; ++i)
		ChargeLog(slog, prod, A[i]);

	CloseLog(slog, prod);
	return prod;
}

/* log(prod(A[i++])) */
TARGET float LogProdx(float* A, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return LogProdNEOx(A, n);
#else
	case 4: return LogProd512x(A, n);
	case 3: return LogProdAVXx(A, n);
	case 2: return LogProdSSEx(A, n);
#endif
	}

	int64 slog = 0; double prod = 1;
	for (int64 i = 0; i < n; ++i)
		ChargeLog(slog, prod, A[i]);

	CloseLog(slog, prod);
	return prod;
}

/* log(prod(A[i += sep])) */
TARGET double LogProd(double* A, int64 n, int64 sep)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return LogProdNEO(A, n, sep);
#else
	case 4: return LogProd512(A, n, sep);
	case 3: return LogProdAVX(A, n, sep);
	case 2: return LogProdSSE(A, n, sep);
#endif
	}

	int64 slog = 0; double prod = 1;
	for (int64 i = 0; i < n; ++i, A += sep)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);
	return prod;
}

/* log(prod(A[i += sep])) */
TARGET double LogProd(float* A, int64 n, int64 sep)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return LogProdNEO(A, n, sep);
#else
	case 4: return LogProd512(A, n, sep);
	case 3: return LogProdAVX(A, n, sep);
	case 2: return LogProdSSE(A, n, sep);
#endif
	}

	int64 slog = 0; double prod = 1;
	for (int64 i = 0; i < n; ++i, A += sep)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);
	return prod;
}

/* log(prod(A[i += sep])) */
TARGET float LogProdx(float* A, int64 n, int64 sep)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return LogProdNEOx(A, n, sep);
#else
	case 4: return LogProd512x(A, n, sep);
	case 3: return LogProdAVXx(A, n, sep);
	case 2: return LogProdSSEx(A, n, sep);
#endif
	}

	int64 slog = 0; double prod = 1;
	for (int64 i = 0; i < n; ++i, A += sep)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);
	return prod;
}

/* log(prod(A[i += sep] / B[i += sep])) */
TARGET double LogProdDiv(double* A, double* B, int64 n, int64 sep)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return LogProdDivNEO(A, B, n, sep);
#else
	case 4: return LogProdDiv512(A, B, n, sep);
	case 3: return LogProdDivAVX(A, B, n, sep);
	case 2: return LogProdDivSSE(A, B, n, sep);
#endif
	}

	int64 slog1 = 0; double prod1 = 1;
	int64 slog2 = 0; double prod2 = 1;
	for (int64 i = 0; i < n; ++i, A += sep, B += sep)
	{
		ChargeLog(slog1, prod1, *A);
		ChargeLog(slog2, prod2, *B);
	}

	CloseLog(slog1, prod1);
	CloseLog(slog2, prod2);
	return prod1 - prod2;
}

/* log(prod(A[i += sep] / B[i += sep])) */
TARGET double LogProdDiv(float* A, float* B, int64 n, int64 sep)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return LogProdDivNEO(A, B, n, sep);
#else
	case 4: return LogProdDiv512(A, B, n, sep);
	case 3: return LogProdDivAVX(A, B, n, sep);
	case 2: return LogProdDivSSE(A, B, n, sep);
#endif
	}

	int64 slog1 = 0, slog2 = 0; 
	double prod1 = 1, prod2 = 1;
	for (int64 i = 0; i < n; ++i, A += sep, B += sep)
	{
		ChargeLog(slog1, prod1, *A);
		ChargeLog(slog2, prod2, *B);
	}

	CloseLog(slog1, prod1);
	CloseLog(slog2, prod2);

	return prod1 - prod2;
}

/* log(prod(A[i += sep] / B[i += sep])) */
TARGET float LogProdDivx(float* A, float* B, int64 n, int64 sep)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return LogProdDivNEOx(A, B, n, sep);
#else
	case 4: return LogProdDiv512x(A, B, n, sep);
	case 3: return LogProdDivAVXx(A, B, n, sep);
	case 2: return LogProdDivSSEx(A, B, n, sep);
#endif
	}

	int64 slog1 = 0; double prod1 = 1;
	int64 slog2 = 0; double prod2 = 1;

	for (int64 i = 0; i < n; ++i, A += sep, B += sep)
	{
		ChargeLog(slog1, prod1, *A);
		ChargeLog(slog2, prod2, *B);
	}

	CloseLog(slog1, prod1);
	CloseLog(slog2, prod2);

	return prod1 - prod2;
}

/* Natural logarithm with bounds */
TARGET double MyLog(double val)
{
	if (val < 1e-300)
		return -6.90775527898214000E+02;
	else if (val > 1e300)
		return +6.90775527898214000E+02;
	return log(val);
}

/* Charge a value to fast sum log */
TARGET void OpenLog(int64& slog, double& prod)
{
	slog = 0;
	prod = 1;
}

/* Charge a value to fast sum log */
TARGET void OpenLog(int64* slog, double* prod, int64 n)
{
	SetZero(slog, n);
	SetVal(prod, 1.0, n);
}

/* Add exponent to slog2 */
TARGET void AddExponent(int64& slog2, double& val)
{
	int64& vv = *(int64*)&val;

	// add the exponent to slog2
	slog2 += ((vv & 0x7FF0000000000000) >> 52) - 1023;

	// set the exponent of val to 0
	vv = (vv & 0x800FFFFFFFFFFFFF) | 0x3FF0000000000000;
}

/* Add exponent to slog2 */
TARGET void AddExponent(int64& slog2, float& val)
{
	int& vv = *(int*)&val;

	// add the exponent to slog2
	slog2 += ((vv & 0x7F800000) >> 23) - 127;

	// set the exponent of val to 0
	vv = (vv & 0x807FFFFF) | 0x3F800000;
}

/* Charge a value to fast sum log */
TARGET void ChargeLog(int64& slog, double& prod, double val)
{
	if (val < DOUBLE_UNDERFLOW || val > DOUBLE_OVERFLOW) [[unlikely]]
		AddExponent(slog, val);

	prod *= val;

	if (prod < DOUBLE_UNDERFLOW || prod > DOUBLE_OVERFLOW) [[unlikely]]
		AddExponent(slog, prod);
}

/* Charge a value to fast sum log */
TARGET void ChargeLog(int64& slog, double& prod, float val)
{
	prod *= val;

	if (prod < DOUBLE_UNDERFLOW || prod > DOUBLE_OVERFLOW) [[unlikely]]
		AddExponent(slog, prod);
}

/* Charge a value to fast sum log */
TARGET void ChargeLog(int64* slog, double* prod, double* val, int64 n, int64 sep)
{
	for (int64 i = 0; i < n; ++i, val += sep)
		ChargeLog(slog[i], prod[i], *val);
}

/* Charge a value to fast sum log */
TARGET void ChargeLog(int64* slog, double* prod, float* val, int64 n, int64 sep)
{
	for (int64 i = 0; i < n; ++i, val += sep)
		ChargeLog(slog[i], prod[i], *val);
}

/* Charge a value to fast sum log */
TARGET void ChargeLog(int64* slog, double* prod, double* val, int64 n)
{
	for (int64 i = 0; i < n; ++i)
		ChargeLog(slog[i], prod[i], val[i]);
}

/* Charge a value to fast sum log */
TARGET void ChargeLog(int64* slog, double* prod, float* val, int64 n)
{
	for (int64 i = 0; i < n; ++i)
		ChargeLog(slog[i], prod[i], val[i]);
}

/////////////////////////////////////////////

/* Add exponent to slog2 */
TARGET void AddExponentAtomic(int64& slog2, double& val)
{
	uint64& vv = *(uint64*)&val;
	atomic<int64>& slog_atomic = *(atomic<int64>*) & slog2;

	uint64 ov, nv;
	for (ov = vv, nv = (vv & 0x800FFFFFFFFFFFFF) | 0x3FF0000000000000;
#if defined(__clang__) || defined(__GNUC__)
		!__sync_bool_compare_and_swap(&vv, ov, nv);
#else
		ov != InterlockedCompareExchange(&vv, nv, ov);
#endif
		ov = vv, nv = (vv & 0x800FFFFFFFFFFFFF) | 0x3FF0000000000000);

	// add the exponent to slog2
	slog_atomic += ((ov & 0x7FF0000000000000) >> 52) - 1023;
}

/* Add exponent to slog2 */
TARGET void AddExponentAtomic(int64& slog2, float& val)
{
	uint& vv = *(uint*)&val;
	atomic<int64>& slog_atomic = *(atomic<int64>*) & slog2;

	uint ov, nv;
	for (ov = vv, nv = (ov & 0x807FFFFF) | 0x3F800000;
#if defined(__clang__) || defined(__GNUC__)
		!__sync_bool_compare_and_swap(&vv, ov, nv);
#else
		 ov != InterlockedCompareExchange(&vv, nv, ov);
#endif
		 ov = vv, nv = (ov & 0x807FFFFF) | 0x3F800000);


	// add the exponent to slog2
	slog_atomic += ((ov & 0x7F800000) >> 23) - 127;
}

/* Charge a value to fast sum log */
TARGET void ChargeLogAtomic(int64& slog, double& prod, double val)
{
	if (val < DOUBLE_UNDERFLOW || val > DOUBLE_OVERFLOW) [[unlikely]]
		AddExponentAtomic(slog, val);

	AtomicMulFloat(prod, val);

	if (prod < DOUBLE_UNDERFLOW || prod > DOUBLE_OVERFLOW) [[unlikely]]
		AddExponentAtomic(slog, prod);
}

/* Charge a value to fast sum log */
TARGET void ChargeLogAtomic(int64& slog, double& prod, float val)
{
	AtomicMulFloat(prod, val);

	if (prod < DOUBLE_UNDERFLOW || prod > DOUBLE_OVERFLOW) [[unlikely]]
		AddExponentAtomic(slog, prod);
}

/* Charge a value to fast sum log */
TARGET void ChargeLogAtomic(int64* slog, double* prod, double* val, int64 n, int64 sep)
{
	for (int64 i = 0; i < n; ++i, val += sep)
		ChargeLogAtomic(slog[i], prod[i], *val);
}

/* Charge a value to fast sum log */
TARGET void ChargeLogAtomic(int64* slog, double* prod, float* val, int64 n, int64 sep)
{
	for (int64 i = 0; i < n; ++i, val += sep)
		ChargeLogAtomic(slog[i], prod[i], *val);
}

/* Charge a value to fast sum log */
TARGET void ChargeLogAtomic(int64* slog, double* prod, double* val, int64 n)
{
	for (int64 i = 0; i < n; ++i)
		ChargeLogAtomic(slog[i], prod[i], val[i]);
}

/* Charge a value to fast sum log */
TARGET void ChargeLogAtomic(int64* slog, double* prod, float* val, int64 n)
{
	for (int64 i = 0; i < n; ++i)
		ChargeLogAtomic(slog[i], prod[i], val[i]);
}

/////////////////////////////////////////////

/* Charge a value to fast sum log, convert slog to double */
TARGET void CloseLog(int64& slog, double& prod)
{
	double& slog2 = *(double*)&slog;
	prod = slog2 = (slog + log2(prod)) * 0.693147180559945;
}

/* Finalize fast sum log */
TARGET void CloseLog(int64* slog, double* prod, int64 n)
{
	for (int64 i = 0; i < n; ++i)
		CloseLog(slog[i], prod[i]);
}

/* Count non-zero elements */
TARGET int64 CountNonZero(byte* A, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return CountNonZeroNEO(A, n);
#else
	case 4: return CountNonZero512(A, n);
	case 3: return CountNonZeroAVX(A, n);
	case 2: return CountNonZeroSSE(A, n);
#endif
	}

	int64 re = 0;
	for (int64 i = 0; i < n; ++i)
		if (A[i]) re++;
	return re;
}

/* Sum of A */
TARGET double Sum(double* A, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumNEO(A, n);
#else
	case 4: return Sum512(A, n);
	case 3: return SumAVX(A, n);
	case 2: return SumSSE(A, n);
#endif
	}

	double re = 0;
	for (int64 i = 0; i < n; ++i)
		re += A[i];
	return re;
}

/* Sum of A */
TARGET double Sum(float* A, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumNEO(A, n);
#else
	case 4: return Sum512(A, n);
	case 3: return SumAVX(A, n);
	case 2: return SumSSE(A, n);
#endif
	}

	double re = 0;
	for (int64 i = 0; i < n; ++i)
		re += A[i];
	return re;
}

/* Sum of A */
TARGET float Sumx(float* A, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumNEOx(A, n);
#else
	case 4: return Sum512x(A, n);
	case 3: return SumAVXx(A, n);
	case 2: return SumSSEx(A, n);
#endif
	}

	float re = 0;
	for (int64 i = 0; i < n; ++i)
		re += A[i];
	return re;
}

/* Sum of A */
TARGET int64 Sum(byte* A, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumNEO(A, n);
#else
	case 4: return Sum512(A, n);
	case 3: return SumAVX(A, n);
	case 2: return SumSSE(A, n);
#endif
	}

	uint64 re = 0;
	for (int64 i = 0; i < n; ++i)
		re += A[i];
	return re;
}

/* re += A[i += sep] */
TARGET double Sum(double* A, int64 n, int64 sep)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumNEO(A, n, sep);
#else
	case 4: return Sum512(A, n, sep);
	case 3: return SumAVX(A, n, sep);
	case 2: return SumSSE(A, n, sep);
#endif
	}

	double re = 0;
	for (int64 i = 0; i < n; ++i, A += sep)
		re += *A;
	return re;
}

/* re += A[i += sep] */
TARGET double Sum(float* A, int64 n, int64 sep)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumNEO(A, n, sep);
#else
	case 4: return Sum512(A, n, sep);
	case 3: return SumAVX(A, n, sep);
	case 2: return SumSSE(A, n, sep);
#endif
	}

	double re = 0;
	for (int64 i = 0; i < n; ++i, A += sep)
		re += *A;
	return re;
}

/* re += A[i += sep] */
TARGET float Sumx(float* A, int64 n, int64 sep)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumNEOx(A, n, sep);
#else
	case 4: return Sum512x(A, n, sep);
	case 3: return SumAVXx(A, n, sep);
	case 2: return SumSSEx(A, n, sep);
#endif
	}

	float re = 0;
	for (int64 i = 0; i < n; ++i, A += sep)
		re += *A;
	return re;
}

/* A[i] = B[0][i] + ... + B[k][i] */
TARGET void Sum(double* A, double** B, int64 k, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumNEO(A, B, k, n);
#else
	case 4: return Sum512(A, B, k, n);
	case 3: return SumAVX(A, B, k, n);
	case 2: return SumSSE(A, B, k, n);
#endif
	}

	for (int64 i = 0; i < n; ++i)
	{
		double Ai = 0;
		for (int64 j = 0; j < k; ++j)
			Ai += B[j][i];
		A[i] = Ai;
	}
}

/* A[i] = B[0][i] + ... + B[k][i] */
TARGET void Sum(float* A, float** B, int64 k, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumNEO(A, B, k, n);
#else
	case 4: return Sum512(A, B, k, n);
	case 3: return SumAVX(A, B, k, n);
	case 2: return SumSSE(A, B, k, n);
#endif
	}

	for (int64 i = 0; i < n; ++i)
	{
		double Ai = 0;
		for (int64 j = 0; j < k; ++j)
			Ai += B[j][i];
		A[i] = Ai;
	}
}

/* Product of A */
TARGET double Prod(double* A, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return ProdNEO(A, n);
#else
	case 4: return Prod512(A, n);
	case 3: return ProdAVX(A, n);
	case 2: return ProdSSE(A, n);
#endif
	}
	
	double re = 1;
	for (int64 i = 0; i < n; ++i)
		re *= A[i];
	return re;
}

/* Product of A */
TARGET double Prod(float* A, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return ProdNEO(A, n);
#else
	case 4: return Prod512(A, n);
	case 3: return ProdAVX(A, n);
	case 2: return ProdSSE(A, n);
#endif
	}
	
	double re = 1;
	for (int64 i = 0; i < n; ++i)
		re *= A[i];
	return re;
}

/* Product of A */
TARGET float Prodx(float* A, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return ProdNEOx(A, n);
#else
	case 4: return Prod512x(A, n);
	case 3: return ProdAVXx(A, n);
	case 2: return ProdSSEx(A, n);
#endif
	}
	
	float re = 1;
	for (int64 i = 0; i < n; ++i)
		re *= A[i];
	return re;
}

/* re *= A[i += sep] */
TARGET double Prod(double* A, int64 n, int64 sep)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return ProdNEO(A, n, sep);
#else
	case 4: return Prod512(A, n, sep);
	case 3: return ProdAVX(A, n, sep);
	case 2: return ProdSSE(A, n, sep);
#endif
	}
	
	double re = 1;
	for (int64 i = 0; i < n; ++i)
		re *= A[i * sep];
	return re;
}

/* re *= A[i += sep] */
TARGET double Prod(float* A, int64 n, int64 sep)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return ProdNEO(A, n, sep);
#else
	case 4: return Prod512(A, n, sep);
	case 3: return ProdAVX(A, n, sep);
	case 2: return ProdSSE(A, n, sep);
#endif
	}
	
	double re = 1;
	for (int64 i = 0; i < n; ++i)
		re *= A[i * sep];
	return re;
}

/* re *= A[i += sep] */
TARGET float Prodx(float* A, int64 n, int64 sep)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return ProdNEOx(A, n, sep);
#else
	case 4: return Prod512x(A, n, sep);
	case 3: return ProdAVXx(A, n, sep);
	case 2: return ProdSSEx(A, n, sep);
#endif
	}
	
	float re = 1;
	for (int64 i = 0; i < n; ++i)
		re *= A[i * sep];
	return re;
}

/* Sum of squared A */
TARGET double SumSquare(double* A, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumSquareNEO(A, n);
#else
	case 4: return SumSquare512(A, n);
	case 3: return SumSquareAVX(A, n);
	case 2: return SumSquareSSE(A, n);
#endif
	}
	
	double re = 0;
	for (int64 i = 0; i < n; ++i)
		re += A[i] * A[i];
	return re;
}

/* Sum of squared A */
TARGET double SumSquare(float* A, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumSquareNEO(A, n);
#else
	case 4: return SumSquare512(A, n);
	case 3: return SumSquareAVX(A, n);
	case 2: return SumSquareSSE(A, n);
#endif
	}
	
	double re = 0;
	for (int64 i = 0; i < n; ++i)
		re += (double)A[i] * (double)A[i];
	return re;
}

/* Sum of squared A */
TARGET float SumSquarex(float* A, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumSquareNEOx(A, n);
#else
	case 4: return SumSquare512x(A, n);
	case 3: return SumSquareAVXx(A, n);
	case 2: return SumSquareSSEx(A, n);
#endif
	}
	
	float re = 0;
	for (int64 i = 0; i < n; ++i)
		re += A[i] * A[i];
	return re;
}

/* Sum of squared A */
TARGET int64 SumSquare(byte* A, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumSquareNEO(A, n);
#else
	case 4: return SumSquare512(A, n);
	case 3: return SumSquareAVX(A, n);
	case 2: return SumSquareSSE(A, n);
#endif
	}
	
	uint64 re = 0;
	for (int64 i = 0; i < n; ++i)
		re += A[i] * A[i];
	return re;
}

/* Sum of A and sum of squared A */
TARGET void SumSumSquare(double* A, int64 n, double& sum, double& sumsq)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumSumSquareNEO(A, n, sum, sumsq);
#else
	case 4: return SumSumSquare512(A, n, sum, sumsq);
	case 3: return SumSumSquareAVX(A, n, sum, sumsq);
	case 2: return SumSumSquareSSE(A, n, sum, sumsq);
#endif
	}
	
	sum = sumsq = 0;
	for (int64 i = 0; i < n; ++i)
	{
		sum += A[i];
		sumsq += A[i] * A[i];
	}
}

/* Sum of A and sum of squared A */
TARGET void SumSumSquare(float* A, int64 n, double& sum, double& sumsq)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumSumSquareNEO(A, n, sum, sumsq);
#else
	case 4: return SumSumSquare512(A, n, sum, sumsq);
	case 3: return SumSumSquareAVX(A, n, sum, sumsq);
	case 2: return SumSumSquareSSE(A, n, sum, sumsq);
#endif
	}
	
	sum = sumsq = 0;
	for (int64 i = 0; i < n; ++i)
	{
		sum += (double)A[i];
		sumsq += (double)A[i] * (double)A[i];
	}
}

/* re = sum(A1[i++] * B[j += sep]) / Sum(A2[i++] * B[j += sep]) */
TARGET double SumProdDiv(double* A1, double* A2, double* B, int64 sep, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumProdDivNEO(A1, A2, B, sep, n);
#else
	case 4: return SumProdDiv512(A1, A2, B, sep, n);
	case 3: return SumProdDivAVX(A1, A2, B, sep, n);
	case 2: return SumProdDivSSE(A1, A2, B, sep, n);
#endif
	}

	double re1 = 0, re2 = 0;
	for (int64 i = 0; i < n; ++i, A1++, A2++, B += sep)
	{
		volatile double v1 = (double)*A1 * (double)*B;
		volatile double v2 = (double)*A2 * (double)*B;
		re1 += v1;
		re2 += v2;
	}
	return re1 / re2;
}

/* re = sum(A1[i++] * B[j += sep]) / Sum(A2[i++] * B[j += sep]) */
TARGET double SumProdDiv(double* A1, float* A2, float* B, int64 sep, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumProdDivNEO(A1, A2, B, sep, n);
#else
	case 4: return SumProdDiv512(A1, A2, B, sep, n);
	case 3: return SumProdDivAVX(A1, A2, B, sep, n);
	case 2: return SumProdDivSSE(A1, A2, B, sep, n);
#endif
	}

	double re1 = 0, re2 = 0;
	for (int64 i = 0; i < n; ++i, A1++, A2++, B += sep)
	{
		volatile double v1 = (double)*A1 * (double)*B;
		volatile double v2 = (double)*A2 * (double)*B;
		re1 += v1;
		re2 += v2;
	}
	return re1 / re2;
}

/* re = sum(A1[i++] * B[j += sep]) / Sum(A2[i++] * B[j += sep]) */
TARGET double SumProdDiv(float* A1, float* A2, float* B, int64 sep, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumProdDivNEO(A1, A2, B, sep, n);
#else
	case 4: return SumProdDiv512(A1, A2, B, sep, n);
	case 3: return SumProdDivAVX(A1, A2, B, sep, n);
	case 2: return SumProdDivSSE(A1, A2, B, sep, n);
#endif
	}

	double re1 = 0, re2 = 0;
	for (int64 i = 0; i < n; ++i, A1++, A2++, B += sep)
	{
		volatile double v1 = (double)*A1 * (double)*B;
		volatile double v2 = (double)*A2 * (double)*B;
		re1 += v1;
		re2 += v2;
	}
	return re1 / re2;
}

/* re = sum(A1[i++] * B[j += sep]) / Sum(A2[i++] * B[j += sep]) */
TARGET float SumProdDivx(float* A1, float* A2, float* B, int64 sep, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumProdDivNEOx(A1, A2, B, sep, n);
#else
	case 4: return SumProdDiv512x(A1, A2, B, sep, n);
	case 3: return SumProdDivAVXx(A1, A2, B, sep, n);
	case 2: return SumProdDivSSEx(A1, A2, B, sep, n);
#endif
	}

	float re1 = 0, re2 = 0;
	for (int64 i = 0; i < n; ++i, A1++, A2++, B += sep)
	{
		volatile float v1 = *A1 * *B;
		volatile float v2 = *A2 * *B;
		re1 += v1;
		re2 += v2;
	}
	return re1 / re2;
}

/* re = sum(A[i++] * B[j += sep]) */
TARGET double SumProd(double* A, double* B, int64 sep, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumProdNEO(A, B, sep, n);
#else
	case 4: return SumProd512(A, B, sep, n);
	case 3: return SumProdAVX(A, B, sep, n);
	case 2: return SumProdSSE(A, B, sep, n);
#endif
	}

	double re = 0;
	for (int64 i = 0; i < n; ++i, A++, B += sep)
		re += *A * *B;
	return re;
}

/* re = sum(A[i++] * B[j += sep]) */
TARGET double SumProd(float* A, float* B, int64 sep, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumProdNEO(A, B, sep, n);
#else
	case 4: return SumProd512(A, B, sep, n);
	case 3: return SumProdAVX(A, B, sep, n);
	case 2: return SumProdSSE(A, B, sep, n);
#endif
	}

	double re = 0;
	for (int64 i = 0; i < n; ++i, A++, B += sep)
		re += (double)*A * (double)*B;
	return re;
}

/* re = sum(A[i++] * B[j += sep]) */
TARGET float SumProdx(float* A, float* B, int64 sep, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumProdNEOx(A, B, sep, n);
#else
	case 4: return SumProd512x(A, B, sep, n);
	case 3: return SumProdAVXx(A, B, sep, n);
	case 2: return SumProdSSEx(A, B, sep, n);
#endif
	}

	float re = 0;
	for (int64 i = 0; i < n; ++i, A++, B += sep)
	{
		volatile float v1 = *A * *B;
		re += v1;
	}
	return re;
}

/* re = sum(A[i] * B[i]) */
TARGET double SumProd(double* A, double* B, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumProdNEO(A, B, n);
#else
	case 4: return SumProd512(A, B, n);
	case 3: return SumProdAVX(A, B, n);
	case 2: return SumProdSSE(A, B, n);
#endif
	}
	
	double re = 0;
	for (int64 i = 0; i < n; ++i)
		re += A[i] * B[i];
	return re;
}

/* re = sum(A[i] * B[i]) */
TARGET double SumProd(float* A, float* B, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumProdNEO(A, B, n);
#else
	case 4: return SumProd512(A, B, n);
	case 3: return SumProdAVX(A, B, n);
	case 2: return SumProdSSE(A, B, n);
#endif
	}
	
	double re = 0;
	for (int64 i = 0; i < n; ++i)
		re += (double)A[i] * (double)B[i];
	return re;
}

/* re = sum(A[i] * B[i]) */
TARGET float SumProdx(float* A, float* B, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumProdNEOx(A, B, n);
#else
	case 4: return SumProd512x(A, B, n);
	case 3: return SumProdAVXx(A, B, n);
	case 2: return SumProdSSEx(A, B, n);
#endif
	}
	
	float re = 0;
	for (int64 i = 0; i < n; ++i)
	{
		volatile float v1 = A[i] * B[i];
		re += v1;
	}
	return re;
}

/* Add B into A, A[i] += B[i] */
TARGET void Add(double* A, double* B, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return AddNEO(A, B, n);
#else
	case 4: return Add512(A, B, n);
	case 3: return AddAVX(A, B, n);
	case 2: return AddSSE(A, B, n);
#endif
	}
	
	for (int64 i = 0; i < n; ++i)
		A[i] += B[i];
}

/* Add B into A, A[i] += B[i] */
TARGET void Add(float* A, float* B, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return AddNEO(A, B, n);
#else
	case 4: return Add512(A, B, n);
	case 3: return AddAVX(A, B, n);
	case 2: return AddSSE(A, B, n);
#endif
	}

	for (int64 i = 0; i < n; ++i)
		A[i] += B[i];
}

/* Add B into A, A[i] += B[i] */
TARGET void Add(int64* A, int64* B, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return AddNEO(A, B, n);
#else
	case 4: return Add512(A, B, n);
	case 3: return AddAVX(A, B, n);
	case 2: return AddSSE(A, B, n);
#endif
	}

	for (int64 i = 0; i < n; ++i)
		A[i] += B[i];
}

/* Add B into A, A[i] += B[i] */
TARGET void Add(int* A, int* B, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return AddNEO(A, B, n);
#else
	case 4: return Add512(A, B, n);
	case 3: return AddAVX(A, B, n);
	case 2: return AddSSE(A, B, n);
#endif
	}

	for (int64 i = 0; i < n; ++i)
		A[i] += B[i];
}

/* Add B into A, A[i] += B */
TARGET void Add(double* A, double B, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return AddNEO(A, B, n);
#else
	case 4: return Add512(A, B, n);
	case 3: return AddAVX(A, B, n);
	case 2: return AddSSE(A, B, n);
#endif
	}

	for (int64 i = 0; i < n; ++i)
		A[i] += B;
}

/* Add B into A, A[i] += B */
TARGET void Add(float* A, float B, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return AddNEO(A, B, n);
#else
	case 4: return Add512(A, B, n);
	case 3: return AddAVX(A, B, n);
	case 2: return AddSSE(A, B, n);
#endif
	}

	for (int64 i = 0; i < n; ++i)
		A[i] += B;
}

/* C[i] = A[i] * B[i] */
TARGET void Mul(double* C, double* A, double* B, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return MulNEO(C, A, B, n);
#else
	case 4: return Mul512(C, A, B, n);
	case 3: return MulAVX(C, A, B, n);
	case 2: return MulSSE(C, A, B, n);
#endif
	}

	for (int64 i = 0; i < n; ++i)
		C[i] = A[i] * B[i];
}

/* C[i] = A[i] * B[i] */
TARGET void Mul(float* C, float* A, float* B, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return MulNEO(C, A, B, n);
#else
	case 4: return Mul512(C, A, B, n);
	case 3: return MulAVX(C, A, B, n);
	case 2: return MulSSE(C, A, B, n);
#endif
	}

	for (int64 i = 0; i < n; ++i)
		C[i] = A[i] * B[i];
}

/* C[i] = A[i] * B */
TARGET void Mul(double* C, double* A, double B, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return MulNEO(C, A, B, n);
#else
	case 4: return Mul512(C, A, B, n);
	case 3: return MulAVX(C, A, B, n);
	case 2: return MulSSE(C, A, B, n);
#endif
	}

	for (int64 i = 0; i < n; ++i)
		C[i] = A[i] * B;
}

/* C[i] = A[i] * B */
TARGET void Mul(float* C, float* A, float B, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return MulNEO(C, A, B, n);
#else
	case 4: return Mul512(C, A, B, n);
	case 3: return MulAVX(C, A, B, n);
	case 2: return MulSSE(C, A, B, n);
#endif
	}

	for (int64 i = 0; i < n; ++i)
		C[i] = A[i] * B;
}

/* A[i] *= B */
TARGET void Mul(double* A, double B, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return MulNEO(A, B, n);
#else
	case 4: return Mul512(A, B, n);
	case 3: return MulAVX(A, B, n);
	case 2: return MulSSE(A, B, n);
#endif
	}
	
	for (int64 i = 0; i < n; ++i)
		A[i] *= B;
}

/* A[i] *= B */
TARGET void Mul(float* A, float B, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return MulNEO(A, B, n);
#else
	case 4: return Mul512(A, B, n);
	case 3: return MulAVX(A, B, n);
	case 2: return MulSSE(A, B, n);
#endif
	}
	
	for (int64 i = 0; i < n; ++i)
		A[i] *= B;
}

/* C[i] += A[i] * B[i] */
TARGET void AddProd(double* C, double* A, double* B, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return AddProdNEO(C, A, B, n);
#else
	case 4: return AddProd512(C, A, B, n);
	case 3: return AddProdAVX(C, A, B, n);
	case 2: return AddProdSSE(C, A, B, n);
#endif
	}

	for (int64 i = 0; i < n; ++i)
		C[i] += A[i] * B[i];
}

/* C[i] += A[i] * B[i] */
TARGET void AddProd(float* C, float* A, float* B, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return AddProdNEO(C, A, B, n);
#else
	case 4: return AddProd512(C, A, B, n);
	case 3: return AddProdAVX(C, A, B, n);
	case 2: return AddProdSSE(C, A, B, n);
#endif
	}

	for (int64 i = 0; i < n; ++i)
		C[i] += A[i] * B[i];
}

/* C[i] += A[i] * B */
TARGET void AddProd(double* C, double* A, double B, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return AddProdNEO(C, A, B, n);
#else
	case 4: return AddProd512(C, A, B, n);
	case 3: return AddProdAVX(C, A, B, n);
	case 2: return AddProdSSE(C, A, B, n);
#endif
	}

	for (int64 i = 0; i < n; ++i)
		C[i] += A[i] * B;
}

/* C[i] += A[i] * B */
TARGET void AddProd(double* C, float* A, double B, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return AddProdNEO(C, A, B, n);
#else
	case 4: return AddProd512(C, A, B, n);
	case 3: return AddProdAVX(C, A, B, n);
	case 2: return AddProdSSE(C, A, B, n);
#endif
	}

	for (int64 i = 0; i < n; ++i)
		C[i] += A[i] * B;
}

/* C[i] += A[i] * B */
TARGET void AddProd(float* C, float* A, float B, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return AddProdNEO(C, A, B, n);
#else
	case 4: return AddProd512(C, A, B, n);
	case 3: return AddProdAVX(C, A, B, n);
	case 2: return AddProdSSE(C, A, B, n);
#endif
	}

	for (int64 i = 0; i < n; ++i)
		C[i] += A[i] * B;
}

/* Set the sum of A to one */
TARGET void Unify(double* A, int64 n)
{
	switch (SIMD_TYPE)
	//A[i] = A[i] * invs + MIN_FREQ * invs
	{
#ifdef __aarch64__
	case 2: return UnifyNEO(A, n);
#else
	case 4: return Unify512(A, n);
	case 3: return UnifyAVX(A, n);
	case 2: return UnifySSE(A, n);
#endif
	}

	double invsum = 1.0 / (Sum(A, n) + n * MIN_FREQ);
	for (int64 i = 0; i < n; ++i)
		A[i] = (A[i] + MIN_FREQ) * invsum;
}

/* Set the sum of A to one */
TARGET void Unify(float* A, int64 n)
{
	switch (SIMD_TYPE)
		//A[i] = A[i] * invs + MIN_FREQ * invs
	{
#ifdef __aarch64__
	case 2: return UnifyNEO(A, n);
#else
	case 4: return Unify512(A, n);
	case 3: return UnifyAVX(A, n);
	case 2: return UnifySSE(A, n);
#endif
	}

	double invsum = 1.0 / (Sum(A, n) + n * MIN_FREQ);
	for (int64 i = 0; i < n; ++i)
		A[i] = ((double)A[i] + MIN_FREQ) * invsum;
}

/* Find next position of val in string A*/
TARGET char* StrNextIdx(char* A, char val, int64 rep, int64 n)
{
	if (!n) return NULL;
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return StrNextIdxNEO(A, val, rep, n);
#else
	case 4: return StrNextIdx512(A, val, rep, n);
	case 3: return StrNextIdxAVX(A, val, rep, n);
	case 2: return StrNextIdxSSE(A, val, rep, n);
#endif
	}

	A++; n--;
	for (int64 i = 0; i < n; ++i, A++)
		if (*A == val && !--rep)
			return A;
	return NULL;
}

/* Count val in string A */
TARGET int64 CountChar(char* A, char val, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return CountCharNEO(A, val, n);
#else
	case 4: return CountChar512(A, val, n);
	case 3: return CountCharAVX(A, val, n);
	case 2: return CountCharSSE(A, val, n);
#endif
	}

	int64 re = 0;
	for (int64 i = 0; i < n; ++i)
		if (A[i] == val) 
			re++;
	return re;
}

/* A[i] = B */
TARGET void SetVal(double* A, double B, int64 n, int64 sep)
{
	for (int64 i = 0; i < n; ++i)
		A[i * sep] = B;
}

/* A[i] = B */
TARGET void SetVal(float* A, float B, int64 n, int64 sep)
{
	for (int64 i = 0; i < n; ++i)
		A[i * sep] = B;
}

/* A[i] = B */
TARGET void SetVal(float* A, double* B, int64 n)
{
	for (int64 i = 0; i < n; ++i)
		A[i] = B[i];
}

/* re = sum(A[i] * B[j] * alen[i * k + j]) */
TARGET double SumProdSMM(ushort* alen, double* A, double* B, int k)
{
	double re = 0;
	for (int i = 0; i < k; ++i)
		for (int j = 0; j < k; ++j)
			re += A[i] * B[j] * ((int)alen[i] - (int)alen[j]) * ((int)alen[i] - (int)alen[j]);
	return re;
}

/* re = sum(A[i] * B[j] * alen[i * k + j]) */
TARGET double SumProdSMM(ushort* alen, float* A, float* B, int k)
{
	double re = 0;
	for (int i = 0; i < k; ++i)
		for (int j = 0; j < k; ++j)
			re += (double)A[i] * (double)B[j] * ((int)alen[i] - (int)alen[j]) * ((int)alen[i] - (int)alen[j]);
	return re;
}

/* re += sum(A[i] * alen[i * k]) */
TARGET double SumProdSMM(ushort* alen, double* A, ushort j, int k)
{
	double re = 0;
	for (int i = 0, alenj = alen[j]; i < k; ++i)
		re += A[i] * ((int)alen[i] - alenj) * ((int)alen[i] - alenj);
	return re;
}

/* re += sum(A[i] * alen[i * k]) */
TARGET double SumProdSMM(ushort* alen, float* A, ushort j, int k)
{
	float re = 0;
	for (int i = 0, alenj = alen[j]; i < k; ++i)
		re += (double)A[i] * ((int)alen[i] - alenj) * ((int)alen[i] - alenj);
	return re;
}

/* Count number of non-zero elements */
TARGET int64 CountNonZero(double* A, int64 n)
{
	//allele freq
	int64 count = 0;
	for (int64 i = 0; i < n; ++i)
		if (A[i]) count++;
	return count;
}

/* Set the sum of A to one */
TARGET void Unify(double* A, int64 m, int64 n)
{
	for (int64 i = 0; i < m; ++i)
	{
		Unify(A, n);
		A += n;
	}
}

/* Set the sum of A to one */
TARGET void Unify(float* A, int64 m, int64 n)
{
	for (int64 i = 0; i < m; ++i)
	{
		Unify(A, n);
		A += n;
	}
}

/* Set the sum of A to one, convert int64 to double */
TARGET void UnifyInt64ToDouble(int64* A, int64 m, int64 n)
{
	double* Af = (double*)A;

	for (int64 i = 0, end = m * n; i < end; ++i)
		Af[i] = (double)A[i];

	for (int64 i = 0; i < m; ++i)
	{
		Unify(Af, n);
		Af += n;
	}
}

/* Calculate SSWP in AMOVA */
TARGET double SSP(double* p, int64 k, int64 nhap, bool isiam, ushort* alen2)
{
	//freq array, without missing alleles
	if (!k || !nhap) return 0;
	if (isiam)
	{
		double s1 = (double)nhap, s2 = SumSquare(p, k) * nhap * nhap;
		if (!s1) return 0;
		return (s1 * s1 - s2) / (double)(2 * s1);
	}
	else
	{
		alen2 += k;
		double re = 0;
		for (int64 i = 0; i < k; ++i)
			for (int64 j = 0; j < i; ++j)
				re += p[i] * p[j] * alen2[i * k + j];
		return re * nhap;
	}
}

/* Calculate SSWP in AMOVA */
TARGET double SSP(float* p, int64 k, int64 nhap, bool isiam, ushort* alen2)
{
	//freq array, without missing alleles
	if (!k || !nhap) return 0;
	if (isiam)
	{
		double s1 = nhap, s2 = SumSquare(p, k) * nhap * nhap;
		if (!s1) return 0;
		return (s1 * s1 - s2) / (float)(2 * s1);
	}
	else
	{
		alen2 += k;
		double re = 0;
		for (int64 i = 0; i < k; ++i)
			for (int64 j = 0; j < i; ++j)
				re += p[i] * p[j] * alen2[i * k + j];
		return re * nhap;
	}
}

/* Calculate SSTOT in AMOVA */
TARGET double SSC(double* a, int64 k, bool isiam, ushort* alen2)
{
	//allele count array, without missing alleles
	if (!k) return 0;
	double s1 = 0, s2 = 0;
	if (isiam)
	{
		SumSumSquare(a, k, s1, s2);
		if (!s1) return 0;
		return (s1 * s1 - s2) / (double)(2 * s1);
	}
	else
	{
		alen2 += k;
		double re = 0, nt = 0;
		for (int64 i = 0; i < k; ++i)
		{
			nt += a[i];
			for (int64 j = 0; j < i; ++j)
				re += a[i] * a[j] * alen2[i * k + j];
		}
		return nt > 0 ? re / nt : 0;
	}
}

/* Calculate SSTOT in AMOVA */
TARGET double SSC(float* a, int64 k, bool isiam, ushort* alen2)
{
	//allele count array, without missing alleles
	if (!k) return 0;
	double s1 = 0, s2 = 0;
	if (isiam)
	{
		SumSumSquare(a, k, s1, s2);
		if (!s1) return 0;
		return (s1 * s1 - s2) / (2 * s1);
	}
	else
	{
		alen2 += k;
		double re = 0, nt = 0;
		for (int64 i = 0; i < k; ++i)
		{
			nt += a[i];
			for (int64 j = 0; j < i; ++j)
				re += a[i] * a[j] * alen2[i * k + j];
		}
		return nt > 0 ? re / nt : 0;
	}
}

/* Add a value to sum and add count */
TARGET void ChargeSum(double val, double& mean, int& count)
{
	if (IsNormal(val))
	{
		mean += val;
		count++;
	}
}

/* Add a value to sum and add count */
TARGET void ChargeSum(float val, float& mean, int& count)
{
	if (IsNormal(val))
	{
		mean += val;
		count++;
	}
}

/* https://github.com/ygalanter/CyberGeeks/blob/master/src/math.c */
TARGET double MyRInt(double x)
{
	double t = floor(fabs(x) + 0.5);
	return (x < 0.0) ? -t : t;
}

/* Core function of cosine */
TARGET double CosCore(double x)
{
	double x2 = x * x;
	double x4 = x2 * x2;
	double x8 = x4 * x4;
	return (-2.7236370439787708e-7 * x2 + 2.4799852696610628e-5) * x8 +
		   (-1.3888885054799695e-3 * x2 + 4.1666666636943683e-2) * x4 +
		   (-4.9999999999963024e-1 * x2 + 1.0000000000000000e+0);
}

/* Core function of sine */
TARGET double SinCore(double x)
{
	double x2 = x * x;
	double x4 = x2 * x2;
	return ((2.7181216275479732e-6 * x2 - 1.9839312269456257e-4) * x4 +
		   (8.3333293048425631e-3 * x2 - 1.6666666640797048e-1)) * x2 * x + x;
}

/* Core function of arc sine */
TARGET double ArcSinCore(double x)
{
	double x2 = x * x;
	double x4 = x2 * x2;
	double x8 = x4 * x4;
	return (((4.5334220547132049e-2 * x2 - 1.1226216762576600e-2) * x4 +
		   (2.6334281471361822e-2 * x2 + 2.0596336163223834e-2)) * x8 +
		   (3.0582043602875735e-2 * x2 + 4.4630538556294605e-2) * x4 +
		   (7.5000364034134126e-2 * x2 + 1.6666666300567365e-1)) * x2 * x + x;
}

/* Sine function */
TARGET double MySin(double x)
{
	double q = MyRInt(x * 6.3661977236758138e-1); 
	int quadrant = (int)q;
	double t = x - q * 1.5707963267923333e+00;
	t = t - q * 2.5633441515945189e-12;
	t = quadrant & 1 ? CosCore(t) : SinCore(t);
	return (quadrant & 2) ? -t : t;
}

/* Cosine function */
TARGET double MyCos(double x)
{
	return MySin(x + (M_PI / 2));
}

/* ArcCosine function */
TARGET double MyArcCos(double x)
{
	double xa = abs(x);
	double t = xa > 0.5625 ? 
		2.0 * ArcSinCore(sqrt(0.5 * (1.0 - xa))) : 
		1.5707963267948966 - ArcSinCore(xa);
	return (x < 0.0) ? (3.1415926535897932 - t) : t;
}

/* ArcSine function */
TARGET double MyArcSin(double x)
{
	return (M_PI / 2) - MyArcCos(x);
}

/* Tangent function */
TARGET double MyTan(double x)
{
	return MySin(x) / MyCos(x);
}
