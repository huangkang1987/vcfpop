/* Math  Functions */

#pragma once
#include "vcfpop.h"

#define Max(a,b)            (((a) > (b)) ? (a) : (b))
#define Min(a,b)            (((a) < (b)) ? (a) : (b))

/* Get square root of val, account some slightly negative input */
TARGET double MySqrt(double val);

/* Get 1D lower triangular index */
TARGET int GetLowerTriangularId(int id1, int id2, int n);

/* Power of x with a integer index */
TARGET double IntegerPower(double x, int index);

/* Is an erroneous real number */
TARGET bool IsError(double x);

/* Is a normal real number */
TARGET bool IsNormal(double x);

/* Find the index of the mimumum element */
TARGETMMX int64 GetMinIdx(double* A, int64 n, double& val);

/* Find maximum and minimum element of A */
TARGETMMX void GetMinMaxVal(double* A, int64 n, double& minv, double& maxv);

/* Find minimum element of A */
TARGETMMX double GetMinVal(double* A, int64 n);

/* Find minimum element of A */
TARGETMMX int64 GetMinVal(int64* A, int64 n);

/* Find Maximum element of A */
TARGETMMX double GetMaxVal(double* A, int64 n);

/* Find Maximum element of A */
TARGETMMX double GetMaxVal(double* A, int64 n, int64 sep);

/* A[i] = B[i] */
TARGETMMX void SetVal(uint* A, ushort* B, int64 n);

/* log(prod(A[i++])) */
TARGETMMX double LogProd(double* A, int64 n);

/* log(prod(A[i += sep])) */
TARGETMMX double LogProd(double* A, int64 n, int64 sep);

/* log(prod(A[i += sep] / B[i += sep])) */
TARGETMMX double LogProdDiv(double* A, double* B, int64 n, int64 sep);

/* Count non-zero elements */
TARGETMMX int64 CountNonZero(byte* A, int64 n);

/* Sum of A */
TARGETMMX double Sum(double* A, int64 n);

/* Sum of A */
TARGETMMX int64 Sum(byte* A, int64 n);

/* re += A[i += sep] */
TARGETMMX double Sum(double* A, int64 n, int64 sep);

/* A[i] = B[0][i] + ... + B[k][i] */
TARGETMMX void Sum(double* A, double** B, int64 k, int64 n);

/* Product of A */
TARGETMMX double Prod(double* A, int64 n);

/* re *= A[i += sep] */
TARGETMMX double Prod(double* A, int64 n, int64 sep);

/* Sum of squared A */
TARGETMMX double SumSquare(double* A, int64 n);

/* Sum of squared A */
TARGETMMX int64 SumSquare(byte* A, int64 n);

/* Sum of A and sum of squared A */
TARGETMMX void SumSumSquare(double* A, int64 n, double& sum, double& sumsq);

/* re = sum(A1[i++] * B[j += sep]) / Sum(A2[i++] * B[j += sep]) */
TARGETMMX double SumProdDiv(double* A1, double* A2, double* B, int64 sep, int64 n);

/* re = sum(A[i++] * B[j += sep]) */
TARGETMMX double SumProd(double* A, double* B, int64 sep, int64 n);

/* re = sum(A[i] * B[i]) */
TARGETMMX double SumProd(double* A, double* B, int64 n);

/* Add B into A, A[i] += B[i] */
TARGETMMX void Add(double* A, double* B, int64 n);

/* Add B into A, A[i] += B[i] */
TARGETMMX void Add(int64* A, int64* B, int64 n);

/* Add B into A, A[i] += B[i] */
TARGETMMX void Add(int* A, int* B, int64 n);

/* Add B into A, A[i] += B */
TARGETMMX void Add(double* A, double B, int64 n);

/* C[i] = A[i] * B[i] */
TARGETMMX void Mul(double* C, double* A, double* B, int64 n);

/* C[i] = A[i] * B */
TARGETMMX void Mul(double* C, double* A, double B, int64 n);

/* A[i] *= B */
TARGETMMX void Mul(double* A, double B, int64 n);

/* C[i] += A[i] * B[i] */
TARGETMMX void AddProd(double* C, double* A, double* B, int64 n);

/* C[i] += A[i] * B */
TARGETMMX void AddProd(double* C, double* A, double B, int64 n);

/* Set the sum of A to one */
TARGETMMX void Unify(double* A, int64 n);

/* Find next position of val in string A*/
TARGETMMX char* StrNextIdx(char* A, char val, int64 rep, int64 n);

/* Count val in string A */
TARGETMMX int64 CountChar(char* A, char val, int64 n);

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
template<typename T>
TARGET void SetVal(T* A, T B, int64 n)
{
	for (int64 i = 0; i < n; ++i)
		*A++ = B;
}

/* A[i] = B */
template<typename T>
TARGET void SetVal(T* A, T B, int n)
{
	for (int i = 0; i < n; ++i)
		*A++ = B;
}

/* re = sum(A[i] * B[j] * alen[i * k + j]) */
TARGET double SumProdSMM(ushort* alen, double* A, double* B, int k);

/* re += sum(A[i] * alen[i * k]) */
TARGET double SumProdSMM(ushort* alen, double* A, ushort j, int k);

/* Natural logarithm with bounds */
TARGET double MyLog(double val);

/* Charge a value to fast sum log */
template<typename T>
TARGET void OpenLog(T& slog, double& prod)
{
	slog = 0;
	prod = 1;
}

/* Charge a value to fast sum log */
template<typename T>
TARGET void OpenLog(T* slog, double* prod, int64 n)
{
	SetZero(slog, n);
	SetVal(prod, 1.0, n);
}

/* Add exponent to slog2 */
template<typename T>
TARGET void AddExponent(T& slog2, double& val)
{
	int64& vv = *(int64*)&val;

	// add the exponent to slog2
	slog2 += ((vv & 0x7FF0000000000000) >> 52) - 1023;

	// set the exponent of val to 0
	vv = (vv & 0x800FFFFFFFFFFFFF) | 0x3FF0000000000000;
}

/* Charge a value to fast sum log */
template<typename T>
TARGET void ChargeLog(T& slog, double& prod, double val)
{
	if (val < DOUBLE_UNDERFLOW || val > DOUBLE_OVERFLOW) [[unlikely]]
		AddExponent(slog, val);

	prod *= val;

	if (prod < DOUBLE_UNDERFLOW || prod > DOUBLE_OVERFLOW) [[unlikely]]
		AddExponent(slog, prod);
}

/* Charge a value to fast sum log, convert slog to double */
template<typename T>
TARGET void CloseLog(T& slog, double& prod)
{
	double& slog2 = *(double*)&slog;
	prod = slog2 = (slog + log2(prod)) * 0.693147180559945;
}

/* Charge a value to fast sum log */
template<typename T>
TARGET void ChargeLog(T* slog, double* prod, double* val, int64 n, int64 sep)
{
	for (int64 i = 0; i < n; ++i, val += sep)
		ChargeLog(slog[i], prod[i], *val);
}

/* Charge a value to fast sum log */
template<typename T>
TARGET void ChargeLog(T* slog, double* prod, double* val, int64 n)
{
	for (int64 i = 0; i < n; ++i)
		ChargeLog(slog[i], prod[i], val[i]);
}

/* Finalize fast sum log */
template<typename T>
TARGET void CloseLog(T* slog, double* prod, int64 n)
{
	for (int64 i = 0; i < n; ++i)
		CloseLog(slog[i], prod[i]);
}

/* Count number of non-zero elements */
TARGET int64 CountNonZero(double* A, int64 n);

/* Set the sum of A to one */
TARGET void Unify(double* A, int64 m, int64 n);

/* Set the sum of A to one, convert int64 to double */
TARGET void UnifyInt64ToDouble(int64* A, int64 m, int64 n);

/* Calculate SSWP in AMOVA */
TARGET double SSP(double* p, int64 k, int64 nhap, bool isiam, ushort* alen2);

/* Calculate SSTOT in AMOVA */
TARGET double SSC(double* a, int64 k, bool isiam, ushort* alen2);

/* Calculate SSTOT in AMOVA */
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

/* Add a value to sum and add count */
template <typename T>
TARGET void ChargeSum(T val, double& mean, int& count)
{
	if (IsNormal(val))
	{
		mean += val;
		count++;
	}
}
