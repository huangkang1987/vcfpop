/* Math  Functions */

#pragma once
#include "vcfpop.h"

#pragma pack(push, 1)

template<typename REAL>
struct RNGSIMD
{
	byte data[1024];

	/* Initialize rng */
	TARGETSIMD RNGSIMD();

	/* Initialize rng */
	TARGETSIMD RNGSIMD(uint64 s, uint64 salt);

	/* Draw 64 64-bit integers in [0,n), 64*n frequencies are in arr */
	TARGETSIMD void Poly(REAL* arr, int n, int64* re);

	/* Draw uniform distriubted intergers */
	TARGETSIMD void XorShift();

	/* Draw uniform distriubted integers */
	template<typename INT>
	TARGETSIMD void Integer(INT* re, int64 n, INT minv = 0, INT maxv = (INT)-1);

	/* Draw uniform distriubted real numbers */
	TARGETSIMD void Uniform(REAL* re, int n, REAL minv = 0, REAL maxv = 1);

	/* Draw uniform distriubted real numbers */
	TARGETSIMD void Normal(REAL* re, int n, REAL mean = 0, REAL sd = 1);
};

#pragma pack(pop)

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

/* re = sum(A1[i++] * B[j += sep]) / Sum(A2[i++] * B[j += sep]) */
TARGET float SumProdDivx(float* A1, float* A2, float* B, int64 sep, int64 n);

/* re = sum(A[i++] * B[j += sep]) */
TARGET double SumProd(double* A, double* B, int64 sep, int64 n);

/* re = sum(A[i++] * B[j += sep]) */
TARGET double SumProd(float* A, float* B, int64 sep, int64 n);

/* re = sum(A[i] * B[i]) */
TARGET double SumProd(double* A, double* B, int64 n);

/* re = sum(A[i] * B[i]) */
TARGET double SumProd(float* A, float* B, int64 n);

/* re = sum(A[i] * B[i] * C[i]) */
TARGET double SumProd(double* A, double* B, double* C, int64 n);

/* re = sum(A[i] * B[i] * C[i]) */
TARGET float SumProd(float* A, float* B, float* C, int64 n);

/* re = sum(A[i] * A[i] * B[i]) */
TARGET double SumSqProd(double* A, double* B, int64 n);

/* re = sum(A[i] * A[i] * B[i]) */
TARGET float SumSqProd(float* A, float* B, int64 n);

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

/* A[i] = B[i] * C[i] */
TARGET void Mul(double* A, double* B, double* C, int64 n);

/* A[i] = B[i] * C[i] */
TARGET void Mul(float* A, float* B, float* C, int64 n);

/* A[i] = B[i] * C */
TARGET void Mul(double* A, double* B, double C, int64 n);

/* A[i] = B[i] * C */
TARGET void Mul(float* A, float* B, float C, int64 n);

/* A[i] *= B */
TARGET void Mul(double* A, double B, int64 n);

/* A[i] *= B */
TARGET void Mul(float* A, float B, int64 n);

/* A[i] = B / C[i] */
TARGET void Div(double* A, double B, double* C, int64 n);

/* A[i] = B / C[i] */
TARGET void Div(float* A, float B, float* C, int64 n);

/* A[i] = B[i] / C[i] */
TARGET void Div(double* A, double* B, double* C, int64 n);

/* A[i] = B[i] / C[i] */
TARGET void Div(float* A, float* B, float* C, int64 n);

/* A[i] += B[i] * C[i] */
TARGET void AddProd(double* A, double* B, double* C, int64 n);

/* A[i] += B[i] * C[i] */
TARGET void AddProd(float* A, float* B, float* C, int64 n);

/* A[i] += B[i] * C */
TARGET void AddProd(double* A, double* B, double C, int64 nx);

/* A[i] += B[i] * C */
TARGET void AddProd(double* A, float* B, double C, int64 n);

/* A[i] += B[i] * C */
TARGET void AddProd(float* A, float* B, float C, int64 n);

/* Set the sum of A to one */
TARGET void Unify(double* A, int64 n);

/* Set the sum of A to one */
TARGET void Unify(float* A, int64 n);

/* Find next position of val in string A*/
TARGET char* StrNextIdx(char* A, char val, int64 rep, int64 n);

/* Count val in string A */
TARGET int64 CountChar(char* A, char val, int64 n);

/* Sum of A */
TARGET float Sumx(float* A, int64 n);

/* re += A[i += sep] */
TARGET float Sumx(float* A, int64 n, int64 sep);

/* Product of A */
TARGET float Prodx(float* A, int64 n);

/* re *= A[i += sep] */
TARGET float Prodx(float* A, int64 n, int64 sep);

/* re = sum(A[i] * B[i]) */
TARGET float SumProdx(float* A, float* B, int64 n);

/* re = sum(A[i] * B[i]) */
TARGET float SumProdx(float* A, float* B, int64 n, int64 sep);

/* Quadratic form A D A' with D being a diagonal matrix, A is m*n, D is n*n, ColMajor */
TARGET void DiagQuadForm(double* res, double* A, double* D, int64 m, int64 n);

/* Quadratic form A D A' with D being a diagonal matrix, A is m*n, D is n*n, ColMajor */
TARGET void DiagQuadForm(float* res, float* A, float* D, int64 m, int64 n);

/* Quadratic form A D B with D being a diagonal matrix, A is m*n, D is n*n, B is n*1, ColMajor */
TARGET void DiagQuadForm(double* res, double* A, double* D, double* B, int64 m, int64 n);

/* Quadratic form A D B with D being a diagonal matrix, A is m*n, D is n*n, B is n*1, ColMajor */
TARGET void DiagQuadForm(float* res, float* A, float* D, float* B, int64 m, int64 n);

/* Quadratic form A D A' with D being a diagonal matrix, A is 1*n, D is n*n, ColMajor */
TARGET void DiagQuadForm(double* res, double* A, double* D, int64 n);

/* Quadratic form A D A' with D being a diagonal matrix, A is 1*n, D is n*n, ColMajor */
TARGET void DiagQuadForm(float* res, float* A, float* D, int64 n);

/* Matrix Muplification for A'B, A is n*m, B is n*p, ColMajor */
TARGET void MatrixMul(double* res, double* A, double* B, int64 m, int64 n, int64 p);

/* Matrix Muplification for A'B, A is n*m, B is n*p, ColMajor */
TARGET void MatrixMul(float* res, float* A, float* B, int64 m, int64 n, int64 p);

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

/* Get 1D lower triangular index */
static forceinline TARGET int GetLowerTriangularId(int id1, int id2, int n)
{
	int minid = std::min(id1, id2);
	return std::max(id1, id2) - ((minid * (1 + minid - n - n)) >> 1);
}

/* Is an erroneous real number */
static forceinline TARGET bool IsError(double x)
{
	return isnan(x) || isinf(x);
}

/* Is an erroneous real number */
static forceinline TARGET bool IsError(float x)
{
	return isnan(x) || isinf(x);
}

/* Is a normal real number */
static forceinline TARGET bool IsNormal(double x)
{
	return !IsError(x);
}

/* Is a normal real number */
static forceinline TARGET bool IsNormal(float x)
{
	return !IsError(x);
}

/* Ceil of log2(v) */
static forceinline TARGET int CeilLog2(int64 v)
{
	if (v <= 1) return 1;

#if defined(__aarch64__)
	return 64 - (int)__builtin_clzll((uint64)(v - 1));
#else
	return 64 - (int)_lzcnt_u64((uint64)(v - 1)); 
#endif
}

/* Ceil of log2(v) */
static forceinline TARGET int CeilLog2(int v)
{
	if (v <= 1) return 1;

#if defined(__aarch64__)
	return 32 - (int)__builtin_clz((uint)(v - 1));
#else
	return 32 - (int)_lzcnt_u32((uint)(v - 1));
#endif
}

/* Ceil of log10(v) */
static forceinline TARGET uint CeilLog10(uint64 v)
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

/* A[i] = B */
static forceinline TARGET void SetVal(double* A, double B, int64 n, int64 sep)
{
	for (int64 i = 0; i < n; ++i)
		A[i * sep] = B;
}

/* A[i] = B */
static forceinline TARGET void SetVal(float* A, float B, int64 n, int64 sep)
{
	for (int64 i = 0; i < n; ++i)
		A[i * sep] = B;
}

/* A[i] = B */
static forceinline TARGET void SetVal(float* A, double* B, int64 n)
{
	for (int64 i = 0; i < n; ++i)
		A[i] = B[i];
}

/* re = sum(A[i] * B[j] * alen[i * k + j]) */
static forceinline TARGET double SumProdSMM(ushort* alen, double* A, double* B, int k)
{
	double re = 0;
	for (int i = 0; i < k; ++i)
		for (int j = 0; j < k; ++j)
			re += A[i] * B[j] * ((int)alen[i] - (int)alen[j]) * ((int)alen[i] - (int)alen[j]);
	return re;
}

/* re = sum(A[i] * B[j] * alen[i * k + j]) */
static forceinline TARGET double SumProdSMM(ushort* alen, float* A, float* B, int k)
{
	double re = 0;
	for (int i = 0; i < k; ++i)
		for (int j = 0; j < k; ++j)
			re += (double)A[i] * (double)B[j] * ((int)alen[i] - (int)alen[j]) * ((int)alen[i] - (int)alen[j]);
	return re;
}

/* re += sum(A[i] * alen[i * k]) */
static forceinline TARGET double SumProdSMM(ushort* alen, double* A, ushort j, int k)
{
	double re = 0;
	for (int i = 0, alenj = alen[j]; i < k; ++i)
		re += A[i] * ((int)alen[i] - alenj) * ((int)alen[i] - alenj);
	return re;
}

/* re += sum(A[i] * alen[i * k]) */
static forceinline TARGET double SumProdSMM(ushort* alen, float* A, ushort j, int k)
{
	float re = 0;
	for (int i = 0, alenj = alen[j]; i < k; ++i)
		re += (double)A[i] * ((int)alen[i] - alenj) * ((int)alen[i] - alenj);
	return re;
}

/* Natural logarithm with bounds */
static forceinline TARGET double MyLog(double val)
{
	if (val < 1e-300)
		return -6.90775527898214000E+02;
	else if (val > 1e300)
		return +6.90775527898214000E+02;
	return log(val);
}

/* Atomic multiply val to ref */
template<typename T>
static forceinline TARGET void AtomicMulFloat(volatile T& ref, T val)
{
	if constexpr (std::is_same_v<T, double>)
	{
		for (double ov = ref, nv = ov * val;
#if defined(__clang__) || defined(__GNUC__)
			!__sync_bool_compare_and_swap((uint64*)&ref, *(uint64*)&ov, *(uint64*)&nv);
#else
			* (uint64*)&ov != InterlockedCompareExchange((uint64*)&ref, *(uint64*)&nv, *(uint64*)&ov);
#endif
			ov = ref, nv = ov * val);
	}
	if constexpr (std::is_same_v<T, float>)
	{
		for (float ov = ref, nv = ov * val;
#if defined(__clang__) || defined(__GNUC__)
			!__sync_bool_compare_and_swap((uint*)&ref, *(uint*)&ov, *(uint*)&nv);
#else
			* (uint*)&ov != InterlockedCompareExchange((uint*)&ref, *(uint*)&nv, *(uint*)&ov);
#endif
			ov = ref, nv = ov * val);
	}
}

/* Atomic add val to ref */
template<typename T>
static forceinline TARGET void AtomicAddFloat(volatile T& ref, T val)
{
	if constexpr (std::is_same_v<T, double>)
	{
		for (double ov = ref, nv = ov + val;
#if defined(__clang__) || defined(__GNUC__)
			!__sync_bool_compare_and_swap((uint64*)&ref, *(uint64*)&ov, *(uint64*)&nv);
#else
			* (uint64*)&ov != InterlockedCompareExchange((uint64*)&ref, *(uint64*)&nv, *(uint64*)&ov);
#endif
			ov = ref, nv = ov + val);
	}
	if constexpr (std::is_same_v<T, float>)
	{
		for (float ov = ref, nv = ov + val;
#if defined(__clang__) || defined(__GNUC__)
			!__sync_bool_compare_and_swap((uint*)&ref, *(uint*)&ov, *(uint*)&nv);
#else
			* (uint*)&ov != InterlockedCompareExchange((uint*)&ref, *(uint*)&nv, *(uint*)&ov);
#endif
			ov = ref, nv = ov + val);
	}
}

/* Atomic add val to ref */
template<typename T>
static forceinline TARGET void AtomicAddFloat(volatile T* ref, T* val, int len)
{
	if constexpr (std::is_same_v<T, double>)
	{
		for (int i = 0; i < len; ++i)
			for (double ov = ref[i], nv = ov + val[i];
#if defined(__clang__) || defined(__GNUC__)
				!__sync_bool_compare_and_swap((uint64*)&ref[i], *(uint64*)&ov, *(uint64*)&nv);
#else
				* (uint64*)&ov != InterlockedCompareExchange((uint64*)&ref[i], *(uint64*)&nv, *(uint64*)&ov);
#endif
				ov = ref[i], nv = ov + val[i]);
	}
	if constexpr (std::is_same_v<T, float>)
	{
		for (int i = 0; i < len; ++i)
			for (float ov = ref[i], nv = ov + val[i];
#if defined(__clang__) || defined(__GNUC__)
				!__sync_bool_compare_and_swap((uint*)&ref[i], *(uint*)&ov, *(uint*)&nv);
#else
				* (uint*)&ov != InterlockedCompareExchange((uint*)&ref[i], *(uint*)&nv, *(uint*)&ov);
#endif
				ov = ref[i], nv = ov + val[i]);
	}
}

/* Atomic set max */
template<typename T>
static forceinline TARGET void AtomicMax(atomic<T>& ref, T value)
{
	for (T prev_value = ref; 
		prev_value < value && !ref.compare_exchange_weak(prev_value, value); 
		prev_value = ref);
}

/* Atomic set min */
template<typename T>
static forceinline TARGET void AtomicMin(atomic<T>& ref, T value)
{
	for (T prev_value = ref; 
		prev_value > value && !ref.compare_exchange_weak(prev_value, value);
		prev_value = ref);
}

/* Charge a value to fast sum log */
static forceinline TARGET void OpenLog(int64& slog, double& prod)
{
	slog = 0;
	prod = 1;
}

/* Charge a value to fast sum log */
static forceinline TARGET void OpenLog(int64* slog, double* prod, int64 n)
{
	SetZero(slog, n);
	SetVal(prod, 1.0, n);
}

/* Add exponent to slog2 */
static forceinline TARGET void AddExponent(int64& slog2, double& val)
{
	int64& vv = *(int64*)&val;

	// add the exponent to slog2
	slog2 += ((vv & 0x7FF0000000000000) >> 52) - 1023;

	// set the exponent of val to 0
	vv = (vv & 0x800FFFFFFFFFFFFF) | 0x3FF0000000000000;
}

/* Add exponent to slog2 */
static forceinline TARGET void AddExponent(int64& slog2, float& val)
{
	int& vv = *(int*)&val;

	// add the exponent to slog2
	slog2 += ((vv & 0x7F800000) >> 23) - 127;

	// set the exponent of val to 0
	vv = (vv & 0x807FFFFF) | 0x3F800000;
}

/* Charge a value to fast sum log */
static forceinline TARGET void ChargeLog(int64& slog, double& prod, double val)
{
	if (val < DOUBLE_UNDERFLOW || val > DOUBLE_OVERFLOW) [[unlikely]]
		AddExponent(slog, val);

	prod *= val;

	if (prod < DOUBLE_UNDERFLOW || prod > DOUBLE_OVERFLOW) [[unlikely]]
		AddExponent(slog, prod);
}

/* Charge a value to fast sum log */
static forceinline TARGET void ChargeLog(int64& slog, double& prod, float val)
{
	prod *= val;

	if (prod < DOUBLE_UNDERFLOW || prod > DOUBLE_OVERFLOW) [[unlikely]]
		AddExponent(slog, prod);
}

/* Charge a value to fast sum log */
static forceinline TARGET void ChargeLog(int64* slog, double* prod, double* val, int64 n, int64 sep)
{
	for (int64 i = 0; i < n; ++i, val += sep)
		ChargeLog(slog[i], prod[i], *val);
}

/* Charge a value to fast sum log */
static forceinline TARGET void ChargeLog(int64* slog, double* prod, float* val, int64 n, int64 sep)
{
	for (int64 i = 0; i < n; ++i, val += sep)
		ChargeLog(slog[i], prod[i], *val);
}

/* Charge a value to fast sum log */
static forceinline TARGET void ChargeLog(int64* slog, double* prod, double* val, int64 n)
{
	for (int64 i = 0; i < n; ++i)
		ChargeLog(slog[i], prod[i], val[i]);
}

/* Charge a value to fast sum log */
static forceinline TARGET void ChargeLog(int64* slog, double* prod, float* val, int64 n)
{
	for (int64 i = 0; i < n; ++i)
		ChargeLog(slog[i], prod[i], val[i]);
}

/* Add exponent to slog2 */
static forceinline TARGET void AddExponentAtomic(int64& slog2, double& val)
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
static forceinline TARGET void AddExponentAtomic(int64& slog2, float& val)
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
static forceinline TARGET void ChargeLogAtomic(int64& slog, double& prod, double val)
{
	if (val < DOUBLE_UNDERFLOW || val > DOUBLE_OVERFLOW) [[unlikely]]
		AddExponentAtomic(slog, val);

	AtomicMulFloat(prod, val);

	if (prod < DOUBLE_UNDERFLOW || prod > DOUBLE_OVERFLOW) [[unlikely]]
		AddExponentAtomic(slog, prod);
}

/* Charge a value to fast sum log */
static forceinline TARGET void ChargeLogAtomic(int64& slog, double& prod, float val)
{
	AtomicMulFloat(prod, (double)val);

	if (prod < DOUBLE_UNDERFLOW || prod > DOUBLE_OVERFLOW) [[unlikely]]
		AddExponentAtomic(slog, prod);
}

/* Charge a value to fast sum log */
static forceinline TARGET void ChargeLogAtomic(int64* slog, double* prod, double* val, int64 n, int64 sep)
{
	for (int64 i = 0; i < n; ++i, val += sep)
		ChargeLogAtomic(slog[i], prod[i], *val);
}

/* Charge a value to fast sum log */
static forceinline TARGET void ChargeLogAtomic(int64* slog, double* prod, float* val, int64 n, int64 sep)
{
	for (int64 i = 0; i < n; ++i, val += sep)
		ChargeLogAtomic(slog[i], prod[i], *val);
}

/* Charge a value to fast sum log */
static forceinline TARGET void ChargeLogAtomic(int64* slog, double* prod, double* val, int64 n)
{
	for (int64 i = 0; i < n; ++i)
		ChargeLogAtomic(slog[i], prod[i], val[i]);
}

/* Charge a value to fast sum log */
static forceinline TARGET void ChargeLogAtomic(int64* slog, double* prod, float* val, int64 n)
{
	for (int64 i = 0; i < n; ++i)
		ChargeLogAtomic(slog[i], prod[i], val[i]);
}

/* Charge a value to fast sum log, convert slog to double */
static forceinline TARGET void CloseLog(int64& slog, double& prod)
{
	double& slog2 = *(double*)&slog;
	prod = slog2 = (slog + log2(prod)) * 0.693147180559945;
}

/* Finalize fast sum log */
static forceinline TARGET void CloseLog(int64* slog, double* prod, int64 n)
{
	for (int64 i = 0; i < n; ++i)
		CloseLog(slog[i], prod[i]);
}

/* Count number of non-zero elements */
static forceinline TARGET int64 CountNonZero(double* A, int64 n)
{
	//allele freq
	int64 count = 0;
	for (int64 i = 0; i < n; ++i)
		if (A[i]) count++;
	return count;
}

/* Set the sum of A to one */
static forceinline TARGET void Unify(double* A, int64 m, int64 n)
{
	for (int64 i = 0; i < m; ++i)
	{
		Unify(A, n);
		A += n;
	}
}

/* Set the sum of A to one */
static forceinline TARGET void Unify(float* A, int64 m, int64 n)
{
	for (int64 i = 0; i < m; ++i)
	{
		Unify(A, n);
		A += n;
	}
}

/* Set the sum of A to one, convert int64 to double */
static forceinline TARGET void UnifyInt64ToDouble(int64* A, int64 m, int64 n)
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
static forceinline TARGET double SSP(double* p, int64 k, int64 nhap, bool isiam, ushort* alen2)
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
static forceinline TARGET double SSP(float* p, int64 k, int64 nhap, bool isiam, ushort* alen2)
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
static forceinline TARGET double SSC(double* a, int64 k, bool isiam, ushort* alen2)
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
static forceinline TARGET double SSC(float* a, int64 k, bool isiam, ushort* alen2)
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

/* Add a value to sum and add weight */
static forceinline TARGET void ChargeWeight(double val, double weight, double& numerator, double& denominator)
{
	double v1 = val * weight;
	if (IsNormal(v1))
	{
		numerator += v1;
		denominator += weight;
	}
}

/* Add a value to sum and add weight */
static forceinline TARGET void ChargeWeight(float val, float weight, float& numerator, float& denominator)
{
	float v1 = val * weight;
	if (IsNormal(v1))
	{
		numerator += v1;
		denominator += weight;
	}
}

/* Add a value to sum and add count */
static forceinline TARGET void ChargeSum(double val, double& mean, int& count)
{
	if (IsNormal(val))
	{
		mean += val;
		count++;
	}
}

/* Add a value to sum and add count */
static forceinline TARGET void ChargeSum(float val, float& mean, int& count)
{
	if (IsNormal(val))
	{
		mean += val;
		count++;
	}
}

/* Core function of cosine */
static forceinline TARGET double CosCore(double x)
{
	double x2 = x * x;
	double x4 = x2 * x2;
	double x8 = x4 * x4;
	return (-2.7236370439787708e-7 * x2 + 2.4799852696610628e-5) * x8 +
		   (-1.3888885054799695e-3 * x2 + 4.1666666636943683e-2) * x4 +
		   (-4.9999999999963024e-1 * x2 + 1.0000000000000000e+0);
}

/* Core function of sine */
static forceinline TARGET double SinCore(double x)
{
	double x2 = x * x;
	double x4 = x2 * x2;
	return ((2.7181216275479732e-6 * x2 - 1.9839312269456257e-4) * x4 +
		   (8.3333293048425631e-3 * x2 - 1.6666666640797048e-1)) * x2 * x + x;
}

/* Core function of arc sine */
static forceinline TARGET double ArcSinCore(double x)
{
	double x2 = x * x;
	double x4 = x2 * x2;
	double x8 = x4 * x4;
	return (((4.5334220547132049e-2 * x2 - 1.1226216762576600e-2) * x4 +
		   (2.6334281471361822e-2 * x2 + 2.0596336163223834e-2)) * x8 +
		   (3.0582043602875735e-2 * x2 + 4.4630538556294605e-2) * x4 +
		   (7.5000364034134126e-2 * x2 + 1.6666666300567365e-1)) * x2 * x + x;
}

/* https://github.com/ygalanter/CyberGeeks/blob/master/src/math.c */
static forceinline TARGET double MyRInt(double x)
{
	double t = floor(fabs(x) + 0.5);
	return (x < 0.0) ? -t : t;
}

/* Sine function */
static forceinline TARGET double MySin(double x)
{
	double q = MyRInt(x * 6.3661977236758138e-1); 
	int quadrant = (int)q;
	double t = x - q * 1.5707963267923333e+00;
	t = t - q * 2.5633441515945189e-12;
	t = quadrant & 1 ? CosCore(t) : SinCore(t);
	return (quadrant & 2) ? -t : t;
}

/* Cosine function */
static forceinline TARGET double MyCos(double x)
{
	return MySin(x + (M_PI / 2));
}

/* ArcCosine function */
static forceinline TARGET double MyArcCos(double x)
{
	double xa = abs(x);
	double t = xa > 0.5625 ? 
		2.0 * ArcSinCore(sqrt(0.5 * (1.0 - xa))) : 
		1.5707963267948966 - ArcSinCore(xa);
	return (x < 0.0) ? (3.1415926535897932 - t) : t;
}

/* ArcSine function */
static forceinline TARGET double MyArcSin(double x)
{
	return (M_PI / 2) - MyArcCos(x);
}

/* Tangent function */
static forceinline TARGET double MyTan(double x)
{
	return MySin(x) / MyCos(x);
}