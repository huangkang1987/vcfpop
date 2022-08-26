/* Math  Functions */

#include "vcfpop.h"

/* Get square root of val, account some slightly negative input */
TARGET double MySqrt(double val)
{
	if (val > 0)
		return sqrt(val);
	else if (val > -MIN_FREQ)
		return 0;
	else
		return NA;
}

/* Get 1D lower triangular index */
TARGET int GetLowerTriangularId(int id1, int id2, int n)
{
	int minid = Min(id1, id2);
	return Max(id1, id2) - ((minid * (1 + minid - n - n)) >> 1);
}

/* Power of x with a integer index */
TARGET double IntegerPower(double x, int index)
{
	switch (index)
	{
	case 0: return 1;
	case 1: return x;
	case 2: return x * x;
	case 3: return x * x * x;
	case 4:
	{
		double x2 = x * x;
		return x2 * x2;
	}
	case 5:
	{
		double x2 = x * x;
		return x2 * x2 * x;
	}
	case 6:
	{
		double x3 = x * x * x;
		return x3 * x3;
	}
	case 7:
	{
		double x3 = x * x * x;
		return x3 * x3 * x;
	}
	case 8:
	{
		double x2 = x * x, x4 = x2 * x2;
		return x4 * x4;
	}
	case 9:
	{
		double x3 = x * x * x;
		return x3 * x3 * x3;
	}
	case 10:
	{
		double x2 = x * x, x4 = x2 * x2;
		return x4 * x4 * x2;
	}
	default: return pow(x, (double)index);
	}
}

/* Is an erroneous real number */
TARGET bool IsError(double x)
{
	return isnan(x) || isinf(x);
}

/* Is a normal real number */
TARGET bool IsNormal(double x)
{
	return !IsError(x);
}

/* Find the index of the mimumum element */
TARGETMMX int64 GetMinIdx(double* A, int64 n, double& val)
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

	val = 1e300;
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
TARGETMMX void GetMinMaxVal(double* A, int64 n, double& minv, double& maxv)
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
	
	minv = 1e300;
	maxv = -1e300;
	for (int64 i = 0; i < n; ++i)
	{
		if (A[i] < minv) minv = A[i];
		if (A[i] > maxv) maxv = A[i];
	}
}

/* Find minimum element of A */
TARGETMMX double GetMinVal(double* A, int64 n)
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
	
	double val = 1e300;
	for (int64 i = 0; i < n; ++i)
	{
		if (A[i] > val) continue;
		val = A[i];
	}
	return val;
}

/* Find minimum element of A */
TARGETMMX int64 GetMinVal(int64* A, int64 n)
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
TARGETMMX double GetMaxVal(double* A, int64 n)
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
	
	double val = -1e300;
	for (int64 i = 0; i < n; ++i)
	{
		if (A[i] < val) continue;
		val = A[i];
	}
	return val;
}

/* Find Maximum element of A */
TARGETMMX double GetMaxVal(double* A, int64 n, int64 sep)
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

	double val = -1e300;
	for (int64 i = 0; i < n; ++i, A += sep)
	{
		if (*A < val) continue;
		val = *A;
	}
	return val;
}

/* A[i] = B[i] */
TARGETMMX void SetVal(uint* A, ushort* B, int64 n)
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
TARGETMMX double LogProd(double* A, int64 n)
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

/* log(prod(A[i += sep])) */
TARGETMMX double LogProd(double* A, int64 n, int64 sep)
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

/* log(prod(A[i += sep] / B[i += sep])) */
TARGETMMX double LogProdDiv(double* A, double* B, int64 n, int64 sep)
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

	int64 slog = 0; double prod = 1;
	for (int64 i = 0; i < n; ++i, A += sep, B += sep)
		ChargeLog(slog, prod, *A / *B);

	CloseLog(slog, prod);
	return prod;
}

/* Count non-zero elements */
TARGETMMX int64 CountNonZero(byte* A, int64 n)
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
TARGETMMX double Sum(double* A, int64 n)
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
TARGETMMX int64 Sum(byte* A, int64 n)
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
TARGETMMX double Sum(double* A, int64 n, int64 sep)
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

/* A[i] = B[0][i] + ... + B[k][i] */
TARGETMMX void Sum(double* A, double** B, int64 k, int64 n)
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

	SetZero(A, n);
	for (int64 j = 0; j < k; ++j)
		for (int64 i = 0; i < n; ++i)
			A[i] += B[j][i];
}

/* Product of A */
TARGETMMX double Prod(double* A, int64 n)
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

/* re *= A[i += sep] */
TARGETMMX double Prod(double* A, int64 n, int64 sep)
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

/* Sum of squared A */
TARGETMMX double SumSquare(double* A, int64 n)
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
TARGETMMX int64 SumSquare(byte* A, int64 n)
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
TARGETMMX void SumSumSquare(double* A, int64 n, double& sum, double& sumsq)
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

/* re = sum(A1[i++] * B[j += sep]) / Sum(A2[i++] * B[j += sep]) */
TARGETMMX double SumProdDiv(double* A1, double* A2, double* B, int64 sep, int64 n)
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
	for (int64 i = 0; i < n; ++i)
	{
		re1 += A1[i] * B[sep * i];
		re2 += A2[i] * B[sep * i];
	}
	return re1 / re2;
}

/* re = sum(A[i++] * B[j += sep]) */
TARGETMMX double SumProd(double* A, double* B, int64 sep, int64 n)
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
	for (int64 i = 0; i < n; ++i)
		re += A[i] * B[i * sep];
	return re;
}

/* re = sum(A[i] * B[i]) */
TARGETMMX double SumProd(double* A, double* B, int64 n)
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

/* Add B into A, A[i] += B[i] */
TARGETMMX void Add(double* A, double* B, int64 n)
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
TARGETMMX void Add(int64* A, int64* B, int64 n)
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
TARGETMMX void Add(int* A, int* B, int64 n)
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
TARGETMMX void Add(double* A, double B, int64 n)
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
TARGETMMX void Mul(double* C, double* A, double* B, int64 n)
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
TARGETMMX void Mul(double* C, double* A, double B, int64 n)
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
TARGETMMX void Mul(double* A, double B, int64 n)
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
TARGETMMX void AddProd(double* C, double* A, double* B, int64 n)
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
TARGETMMX void AddProd(double* C, double* A, double B, int64 n)
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
TARGETMMX void Unify(double* A, int64 n)
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

/* Find next position of val in string A*/
TARGETMMX char* StrNextIdx(char* A, char val, int64 rep, int64 n)
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
TARGETMMX int64 CountChar(char* A, char val, int64 n)
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

/* re = sum(A[i] * B[j] * alen[i * k + j]) */
TARGET double SumProdSMM(ushort* alen, double* A, double* B, int k)
{
	double re = 0;
	for (int i = 0; i < k; ++i)
		for (int j = 0; j < k; ++j)
			re += A[i] * B[j] * ((int)alen[i] - (int)alen[j]) * ((int)alen[i] - (int)alen[j]);
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

/* Natural logarithm with bounds */
TARGET double MyLog(double val)
{
	if (val < 1e-300)
		return -6.90775527898214000E+02;
	else if (val > 1e300)
		return +6.90775527898214000E+02;
	return log(val);
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
