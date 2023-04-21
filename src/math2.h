/* Math  Functions */

#pragma once
#include "vcfpop.h"

#define Max(a,b)            (((a) > (b)) ? (a) : (b))
#define Min(a,b)            (((a) < (b)) ? (a) : (b))

/* Get square root of val, account some slightly negative input */
template<typename REAL>
TARGET REAL MySqrt(REAL val)
{
	if (val > 0)
		return sqrt(val);
	else if (val > -MIN_FREQ)
		return 0;
	else
		return NAN;
}

/* Power of x with a integer index */
template<typename REAL>
TARGET REAL IntegerPower(REAL x, int index)
{
	switch (index)
	{
	case 0: return 1;
	case 1: return x;
	case 2: return x * x;
	case 3: return x * x * x;
	case 4:
	{
		REAL x2 = x * x;
		return x2 * x2;
	}
	case 5:
	{
		REAL x2 = x * x;
		return x2 * x2 * x;
	}
	case 6:
	{
		REAL x3 = x * x * x;
		return x3 * x3;
	}
	case 7:
	{
		REAL x3 = x * x * x;
		return x3 * x3 * x;
	}
	case 8:
	{
		REAL x2 = x * x, x4 = x2 * x2;
		return x4 * x4;
	}
	case 9:
	{
		REAL x3 = x * x * x;
		return x3 * x3 * x3;
	}
	case 10:
	{
		REAL x2 = x * x, x4 = x2 * x2;
		return x4 * x4 * x2;
	}
	default: return pow(x, (REAL)index);
	}
}

/* Get 1D lower triangular index */
TARGET int GetLowerTriangularId(int id1, int id2, int n);

/* Is an erroneous real number */
TARGET bool IsError(double x);

/* Is an erroneous real number */
TARGET bool IsError(float x);

/* Is a normal real number */
TARGET bool IsNormal(double x);

/* Is a normal real number */
TARGET bool IsNormal(float x);

/* Find the index of the mimumum element */
TARGET int64 GetMinIdx(double* A, int64 n, double& val);

/* Find the index of the mimumum element */
TARGET int64 GetMinIdx(float* A, int64 n, float& val);

/* Find maximum and minimum element of A */
TARGET void GetMinMaxVal(double* A, int64 n, double& minv, double& maxv);

/* Find maximum and minimum element of A */
TARGET void GetMinMaxVal(float* A, int64 n, float& minv, float& maxv);

/* Find minimum element of A */
TARGET double GetMinVal(double* A, int64 n);

/* Find minimum element of A */
TARGET float GetMinVal(float* A, int64 n);

/* Find minimum element of A */
TARGET int64 GetMinVal(int64* A, int64 n);

/* Find Maximum element of A */
TARGET double GetMaxVal(double* A, int64 n);

/* Find Maximum element of A */
TARGET float GetMaxVal(float* A, int64 n);

/* Find Maximum element of A */
TARGET double GetMaxVal(double* A, int64 n, int64 sep);

/* Find Maximum element of A */
TARGET float GetMaxVal(float* A, int64 n, int64 sep);

/* A[i] = B[i] */
TARGET void SetVal(uint* A, ushort* B, int64 n);

/* log(prod(A[i++])) */
TARGET double LogProd(double* A, int64 n);

/* log(prod(A[i++])) */
TARGET double LogProd(float* A, int64 n);

/* log(prod(A[i += sep])) */
TARGET double LogProd(double* A, int64 n, int64 sep);

/* log(prod(A[i += sep])) */
TARGET double LogProd(float* A, int64 n, int64 sep);

/* log(prod(A[i += sep] / B[i += sep])) */
TARGET double LogProdDiv(double* A, double* B, int64 n, int64 sep);

/* log(prod(A[i += sep] / B[i += sep])) */
TARGET double LogProdDiv(float* A, float* B, int64 n, int64 sep);

/* Count non-zero elements */
TARGET int64 CountNonZero(byte* A, int64 n);

/* Sum of A */
TARGET double Sum(double* A, int64 n);

/* Sum of A */
TARGET double Sum(float* A, int64 n);

/* Sum of A */
TARGET int64 Sum(byte* A, int64 n);

/* re += A[i += sep] */
TARGET double Sum(double* A, int64 n, int64 sep);

/* re += A[i += sep] */
TARGET double Sum(float* A, int64 n, int64 sep);

/* A[i] = B[0][i] + ... + B[k][i] */
TARGET void Sum(double* A, double** B, int64 k, int64 n);

/* A[i] = B[0][i] + ... + B[k][i] */
TARGET void Sum(float* A, float** B, int64 k, int64 n);

/* Product of A */
TARGET double Prod(double* A, int64 n);

/* Product of A */
TARGET double Prod(float* A, int64 n);

/* re *= A[i += sep] */
TARGET double Prod(double* A, int64 n, int64 sep);

/* re *= A[i += sep] */
TARGET double Prod(float* A, int64 n, int64 sep);

/* Sum of squared A */
TARGET double SumSquare(double* A, int64 n);

/* Sum of squared A */
TARGET double SumSquare(float* A, int64 n);

/* Sum of squared A */
TARGET int64 SumSquare(byte* A, int64 n);

/* Sum of A and sum of squared A */
TARGET void SumSumSquare(double* A, int64 n, double& sum, double& sumsq);

/* Sum of A and sum of squared A */
TARGET void SumSumSquare(float* A, int64 n, double& sum, double& sumsq);

/* re = sum(A1[i++] * B[j += sep]) / Sum(A2[i++] * B[j += sep]) */
TARGET double SumProdDiv(double* A1, double* A2, double* B, int64 sep, int64 n);

/* re = sum(A1[i++] * B[j += sep]) / Sum(A2[i++] * B[j += sep]) */
TARGET double SumProdDiv(double* A1, float* A2, float* B, int64 sep, int64 n);

/* re = sum(A1[i++] * B[j += sep]) / Sum(A2[i++] * B[j += sep]) */
TARGET double SumProdDiv(float* A1, float* A2, float* B, int64 sep, int64 n);

/* re = sum(A[i++] * B[j += sep]) */
TARGET double SumProd(double* A, double* B, int64 sep, int64 n);

/* re = sum(A[i++] * B[j += sep]) */
TARGET double SumProd(float* A, float* B, int64 sep, int64 n);

/* re = sum(A[i] * B[i]) */
TARGET double SumProd(double* A, double* B, int64 n);

/* re = sum(A[i] * B[i]) */
TARGET double SumProd(float* A, float* B, int64 n);

/* Add B into A, A[i] += B[i] */
TARGET void Add(double* A, double* B, int64 n);

/* Add B into A, A[i] += B[i] */
TARGET void Add(float* A, float* B, int64 n);

/* Add B into A, A[i] += B[i] */
TARGET void Add(int64* A, int64* B, int64 n);

/* Add B into A, A[i] += B[i] */
TARGET void Add(int* A, int* B, int64 n);

/* Add B into A, A[i] += B[i] */
TARGET void Add(int* A, int B, int64 n);

/* Add B into A, A[i] += B */
TARGET void Add(double* A, double B, int64 n);

/* Add B into A, A[i] += B */
TARGET void Add(float* A, float B, int64 n);

/* C[i] = A[i] * B[i] */
TARGET void Mul(double* C, double* A, double* B, int64 n);

/* C[i] = A[i] * B[i] */
TARGET void Mul(float* C, float* A, float* B, int64 n);

/* C[i] = A[i] * B */
TARGET void Mul(double* C, double* A, double B, int64 n);

/* C[i] = A[i] * B */
TARGET void Mul(float* C, float* A, float B, int64 n);

/* A[i] *= B */
TARGET void Mul(double* A, double B, int64 n);

/* A[i] *= B */
TARGET void Mul(float* A, float B, int64 n);

/* C[i] += A[i] * B[i] */
TARGET void AddProd(double* C, double* A, double* B, int64 n);

/* C[i] += A[i] * B[i] */
TARGET void AddProd(float* C, float* A, float* B, int64 n);

/* C[i] += A[i] * B */
TARGET void AddProd(double* C, double* A, double B, int64 n);

/* C[i] += A[i] * B */
TARGET void AddProd(double* C, float* A, double B, int64 n);

/* C[i] += A[i] * B */
TARGET void AddProd(float* C, float* A, float B, int64 n);

/* Set the sum of A to one */
TARGET void Unify(double* A, int64 n);

/* Set the sum of A to one */
TARGET void Unify(float* A, int64 n);

/* Find next position of val in string A*/
TARGET char* StrNextIdx(char* A, char val, int64 rep, int64 n);

/* Count val in string A */
TARGET int64 CountChar(char* A, char val, int64 n);

TARGET float LogProdx(float* A, int64 n);

TARGET float LogProdx(float* A, int64 n, int64 sep);

TARGET float LogProdDivx(float* A, float* B, int64 n, int64 sep);

TARGET float Sumx(float* A, int64 n);

TARGET float Sumx(float* A, int64 n, int64 sep);

TARGET float Prodx(float* A, int64 n);

TARGET float Prodx(float* A, int64 n, int64 sep);

TARGET float SumSquarex(float* A, int64 n);

TARGET float SumProdDivx(float* A1, float* A2, float* B, int64 sep, int64 n);

TARGET float SumProdx(float* A, float* B, int64 sep, int64 n);

TARGET float SumProdx(float* A, float* B, int64 n);

/* Find the index of the mimumum element */
template <typename T>
TARGET uint64 GetMinID(T* val, int64 n)
{
	uint64 id = 0;
	for (int64 i = 1; i < n; ++i)
		if (val[id] > val[i])
			id = i;
	return id;
}

/* Find the index of the maximum element */
template <typename T>
TARGET uint64 GetMaxID(T* val, int64 n)
{
	uint64 id = 0;
	for (int64 i = 1; i < n; ++i)
		if (val[id] < val[i])
			id = i;
	return id;
}

/* Ceil of log2(v) */
TARGET inline int CeilLog2(int64 v)
{
	if (v <= 1) return 1;

#ifdef __aarch64__
	return 64 - (int)__builtin_clzll((uint64)(v - 1));
#else
	return 64 - (int)_lzcnt_u64((uint64)(v - 1)); 
#endif
}

/* Ceil of log2(v) */
TARGET inline int CeilLog2(int v)
{
	if (v <= 1) return 1;

#ifdef __aarch64__
	return 32 - (int)__builtin_clz((uint)(v - 1));
#else
	return 32 - (int)_lzcnt_u32((uint)(v - 1));
#endif
}

/* Ceil of log10(v) */
TARGET inline uint CeilLog10(uint64 v)
{
	if (v >= 10000000000000000) return 17;
	if (v >= 1000000000000000) return 16;
	if (v >= 100000000000000) return 15;
	if (v >= 10000000000000) return 14;
	if (v >= 1000000000000) return 13;
	if (v >= 100000000000) return 12;
	if (v >= 10000000000) return 11;
	if (v >= 1000000000) return 10;
	if (v >= 100000000) return 9;
	if (v >= 10000000) return 8;
	if (v >= 1000000) return 7;
	if (v >= 100000) return 6;
	if (v >= 10000) return 5;
	if (v >= 1000) return 4;
	if (v >= 100) return 3;
	if (v >= 10) return 2;
	return 1;
}

/* Set all bytes to 0xFF */
template <typename T>
TARGET void SetFF(T* A, int64 n)
{
	memset(A, 0xFF, sizeof(T) * n);
}

/* Set all bytes to 0x00 */
template <typename T>
TARGET void SetZero(T* A, int64 n)
{
	memset(A, 0, sizeof(T) * n);
}

/* A[i] = B[i] */
template<typename T>
TARGET void SetVal(T* A, T* B, int64 n)
{
	memmove(A, B, sizeof(T) * n);
}

/* A[i] <- 0 if A[i] < 0 */
template<typename T>
TARGET void Truncate(T* A, int64 n)
{
	for (int i = 0; i < n; ++i)
		A[i] = A[i] > 0 ? A[i] : 0;
}

/* A[i] = B */
TARGET void SetVal(double* A, double B, int64 n, int64 sep);

/* A[i] = B */
TARGET void SetVal(double* A, double B, int64 n, int64 sep);

/* A[i] = B */
TARGET void SetVal(float* A, float B, int64 n, int64 sep);

/* A[i] = B */
TARGET void SetVal(float* A, double* B, int64 n);

/* A[i] = B */
template<typename T, typename T3>
TARGET void SetVal(T* A, T B, T3 n)
{
	if constexpr (sizeof(T) == 1)
	{
		memset(A, B, n);
		return;
	}
	if constexpr (sizeof(T) == 2)
	{
		wmemset(A, B, n);
		return;
	}

	for (T3 i = 0; i < n; ++i)
		*A++ = B;
}

/* re = sum(A[i] * B[j] * alen[i * k + j]) */
TARGET double SumProdSMM(ushort* alen, double* A, double* B, int k);

/* re = sum(A[i] * B[j] * alen[i * k + j]) */
TARGET double SumProdSMM(ushort* alen, float* A, float* B, int k);

/* re += sum(A[i] * alen[i * k]) */
TARGET double SumProdSMM(ushort* alen, double* A, ushort j, int k);

/* re += sum(A[i] * alen[i * k]) */
TARGET double SumProdSMM(ushort* alen, float* A, ushort j, int k);

/* Natural logarithm with bounds */
TARGET double MyLog(double val);

/* Charge a value to fast sum log */
TARGET void OpenLog(int64& slog, double& prod);

/* Charge a value to fast sum log */
TARGET void OpenLog(int64* slog, double* prod, int64 n);

/* Add exponent to slog2 */
TARGET void AddExponent(int64& slog2, double& val);

/* Add exponent to slog2 */
TARGET void AddExponent(int64& slog2, float& val);

/* Charge a value to fast sum log */
TARGET void ChargeLog(int64& slog, double& prod, double val);

/* Charge a value to fast sum log */
TARGET void ChargeLog(int64& slog, double& prod, float val);

/* Charge a value to fast sum log */
TARGET void ChargeLog(int64* slog, double* prod, double* val, int64 n, int64 sep);

/* Charge a value to fast sum log */
TARGET void ChargeLog(int64* slog, double* prod, float* val, int64 n, int64 sep);

/* Charge a value to fast sum log */
TARGET void ChargeLog(int64* slog, double* prod, double* val, int64 n);

/* Charge a value to fast sum log */
TARGET void ChargeLog(int64* slog, double* prod, float* val, int64 n);

/////////////////////////////////////////////

/* Add exponent to slog2 */
TARGET void AddExponentAtomic(int64& slog2, double& val);

/* Add exponent to slog2 */
TARGET void AddExponentAtomic(int64& slog2, float& val);

/* Charge a value to fast sum log */
TARGET void ChargeLogAtomic(int64& slog, double& prod, double val);

/* Charge a value to fast sum log */
TARGET void ChargeLogAtomic(int64& slog, double& prod, float val);

/* Charge a value to fast sum log */
TARGET void ChargeLogAtomic(int64* slog, double* prod, double* val, int64 n, int64 sep);

/* Charge a value to fast sum log */
TARGET void ChargeLogAtomic(int64* slog, double* prod, float* val, int64 n, int64 sep);

/* Charge a value to fast sum log */
TARGET void ChargeLogAtomic(int64* slog, double* prod, double* val, int64 n);

/* Charge a value to fast sum log */
TARGET void ChargeLogAtomic(int64* slog, double* prod, float* val, int64 n);

/////////////////////////////////////////////

/* Charge a value to fast sum log, convert slog to double */
TARGET void CloseLog(int64& slog, double& prod);

/* Finalize fast sum log */
TARGET void CloseLog(int64* slog, double* prod, int64 n);

/* Count number of non-zero elements */
TARGET int64 CountNonZero(double* A, int64 n);

/* Set the sum of A to one */
TARGET void Unify(double* A, int64 m, int64 n);

/* Set the sum of A to one */
TARGET void Unify(float* A, int64 m, int64 n);

/* Set the sum of A to one, convert int64 to double */
TARGET void UnifyInt64ToDouble(int64* A, int64 m, int64 n);

/* Calculate SSWP in AMOVA */
TARGET double SSP(double* p, int64 k, int64 nhap, bool isiam, ushort* alen2);

/* Calculate SSWP in AMOVA */
TARGET double SSP(float* p, int64 k, int64 nhap, bool isiam, ushort* alen2);

/* Calculate SSTOT in AMOVA */
TARGET double SSC(double* a, int64 k, bool isiam, ushort* alen2);

/* Calculate SSTOT in AMOVA */
TARGET double SSC(float* a, int64 k, bool isiam, ushort* alen2);

/* Calculate SSTOT in AMOVA 
template <typename T>
TARGET double SSA(T* a, int v, int64 k, bool isiam, ushort* alen2)
{
	//allele array, without missing alleles
	if (!v) return 0;
	int d = 0;
	if (isiam)
	{
		for (int i = 0; i < v; ++i)
			for (int j = i + 1; j < v; ++j)
				if (a[i] != a[j])
					d++;
	}
	else
	{
		alen2 += k;
		for (int i = 0; i < v; ++i)
			for (int j = i + 1; j < v; ++j)
				d += alen2[a[i] * k + a[j]];
	}
	return d / (double)v;
}
*/

/* Swap values for two variables */
template<typename T>
TARGET void Swap(T& a, T& b)
{
	T c = a;
	a = b;
	b = c;
}

/* Bubble Sort */
template <typename T>
TARGET void Sort(T* d, int64 n)
{
	for (int64 i = 0; i < n; ++i)
		for (int64 j = i + 1; j < n; ++j)
			if (d[i] > d[j])
			{
				T ta = d[i];
				d[i] = d[j];
				d[j] = ta;
			}
}

/* Quick Sort */
template <typename T>
TARGET void QuickSort(T* arr, int64 left, int64 right)
{
	int64 i = left, j = right;
	T pivot = arr[(left + right) >> 1];

	while (left < j || i < right)
	{
		while (arr[i] < pivot) i++;
		while (arr[j] > pivot) j--;
		if (i <= j)
			Swap(arr[i++], arr[j--]);
		if (i > j)
		{
			if (left < j)
				QuickSort(arr, left, j);
			if (i < right)
				QuickSort(arr, i, right);
			return;
		}
	}
}

/* Add a value to sum and add weight */
TARGET void ChargeWeight(double val, double weight, double& numerator, double& denominator);

/* Add a value to sum and add weight */
TARGET void ChargeWeight(float val, float weight, float& numerator, float& denominator);

/* Add a value to sum and add count */
TARGET void ChargeSum(double val, double& mean, int& count);

/* Add a value to sum and add count */
TARGET void ChargeSum(float val, float& mean, int& count);

/* Core function of cosine */
TARGET double CosCore(double x);

/* Core function of sine */
TARGET double SinCore(double x);

/* Core function of arc sine */
TARGET double ArcSinCore(double x);

/* Sine function */
TARGET double MySin(double x);

/* Cosine function */
TARGET double MyCos(double x);

/* ArcCosine function */
TARGET double MyArcCos(double x);

/* ArcSine function */
TARGET double MyArcSin(double x);

/* Tangent function */
TARGET double MyTan(double x);
