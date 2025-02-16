#pragma once

#define WINDOWS_IGNORE_PACKING_MISMATCH
#define _CRT_SECURE_NO_WARNINGS

#ifdef __arm64__
    #define __aarch64__
#endif

#ifdef __aarch64__
	#include <arm_neon.h>
#else
	#include <immintrin.h>
#endif

#include <stdexcept>
#include <chrono>
#include <string>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cfloat>
#include <cmath>
#include <math.h>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <new>
#include <map>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <utility>
#include <atomic> 
#include <shared_mutex>
#include <mutex> 
#include <thread>
#include <omp.h>

#define ARMA_ALLOW_FAKE_GCC
#define ARMA_DONT_USE_OPENMP
#define ARMA_DONT_USE_WRAPPER
#define ARMA_NO_DEBUG
#define ARMA_WARN_LEVEL 0
#define ARMA_DONT_PRINT_ERRORS
#include <armadillo>

/* armadillo provide these functions */
#ifndef __APPLE__
    extern "C" void openblas_set_num_threads(int num_threads);
    extern "C" int openblas_get_num_threads();
    extern "C" int openblas_get_parallel();
    extern "C" char* openblas_get_config();
#else
    static void openblas_set_num_threads(int num_threads){ return; }
    static int openblas_get_num_threads(){ return 0; }
    static int openblas_get_parallel(){ return 0; }
    static char* openblas_get_config(){ return NULL; }
#endif

#ifndef __CUDA__
using namespace std::filesystem;
namespace fs = std::filesystem;
#endif

using namespace arma;
using std::pair;
using std::mutex;
using std::shared_mutex;
using std::unordered_map;
using std::map;
using std::vector;
using std::string;
using std::thread;
using std::atomic;
using std::atomic_flag;

template<typename K, typename V>
using umap = std::unordered_map<K, V>;

#define rcol arma::Col<REAL>
#define rrow arma::Row<REAL>
#define rmat arma::Mat<REAL>
#define urow arma::Row<arma::u64>
#define ucol arma::Col<arma::u64>
#define umat arma::Mat<arma::u64>
#define irow arma::Row<arma::s64>
#define icol arma::Col<arma::s64>
#define imat arma::Mat<arma::s64>

#define Z_LARGE64
#define ZLIB_WINAPI
#define LOCKFREE_

typedef unsigned char byte;
typedef unsigned int uint;
typedef unsigned short ushort;
typedef unsigned long long uint64;
typedef long long int64;
typedef std::chrono::time_point<std::chrono::steady_clock> timepoint;

#ifdef OVER4GLOC
	typedef int64 LOCN;
#else
	typedef uint LOCN;
#endif

#ifdef __APPLE__
	#include <mach-o/dyld.h>
	#include <mach/mach.h>
	#include <libproc.h>
	#define fseeko64 fseeko
	#define ftello64 ftello
	#define off64_t long long
#endif

#ifdef _WIN64
	#pragma warning(disable:26451)
	#pragma warning(disable:26495)
	#pragma warning(disable:4706)
	#pragma warning(disable:4996)
	#pragma warning(disable:4366)
	#pragma warning(disable:4310)
	#pragma warning(disable:4244)
	#pragma warning(disable:4503)
	#pragma warning(disable:4819)
	#pragma warning(disable:4702)
	#pragma warning(disable:6386)
	#pragma warning(disable:6387)
	#pragma warning(disable:6031)
	#pragma warning(disable:4661)

	#define TARGET    
	#define TARGETSSE 
	#define TARGETAVX 
	#define TARGET512 
	#define TARGETNEO 
	#define TARGETSIMD 

	#define fseeko64 _fseeki64
	#define ftello64 _ftelli64 

	#include <windows.h>
	#include <psapi.h>
	#pragma comment (lib,"psapi.lib")
	#include <zlib.h>
	#pragma comment (lib,"zlibstat.lib")
	#include <io.h>

	#ifdef LOCKFREE
		typedef atomic_flag LOCK;
	#else
		typedef CRITICAL_SECTION LOCK;
	#endif

	#define PATH_DELIM '\\'
	#define PATH_DELIM_REVERSE '/'
	#define PATH_DOUBLE_DELIM "\\\\"
	#define PATH_DOUBLE_DELIM_REVERSE "//"
#else
	#ifdef __aarch64__
		#define TARGET
		#define TARGETSSE
		#define TARGETAVX
		#define TARGET512
		#define TARGETSIMD 
		#define TARGETNEO  __attribute__((__target__("neon")))
	#else
		#define TARGET     
		#define TARGETSSE 
		#define TARGETAVX  __attribute__((__target__("avx"))) __attribute__((__target__("avx2")))
		#define TARGET512  __attribute__((__target__("avx"))) __attribute__((__target__("avx2"))) __attribute__((__target__("avx512f"))) __attribute__((__target__("avx512bw"))) __attribute__((__target__("avx512dq"))) __attribute__((__target__("avx512vl")))
		#define TARGETSIMD 
		#define TARGETNEO
	#endif

	//#include <segvcatch.h>
	#include <sys/types.h>
	#include <sys/mman.h>
	#include <fcntl.h>
	#include <unistd.h>
	#include <stdarg.h>
	#include <dirent.h>
	#include <zlib.h>
	#pragma comment (lib,"libz.a")
	#include <dlfcn.h>

	#ifdef LOCKFREE
		typedef atomic_flag LOCK;
	#else
		typedef pthread_mutex_t LOCK;
	#endif

	#define Sleep(x) usleep((x)*1000)
	#define PATH_DELIM '/'
	#define PATH_DELIM_REVERSE '\\'
	#define PATH_DOUBLE_DELIM "//"
	#define PATH_DOUBLE_DELIM_REVERSE "\\\\"
#endif

#undef min
#undef max

#if defined(__clang__) || defined(__GNUC__)
	/*
	typedef float		__m128  __attribute__((__vector_size__(16), __aligned__(16)));
	typedef double		__m128d __attribute__((__vector_size__(16), __aligned__(16)));
	typedef long long	__m128i __attribute__((__vector_size__(16), __aligned__(16)));

	typedef float		__m256  __attribute__((__vector_size__(32), __aligned__(32)));
	typedef double		__m256d __attribute__((__vector_size__(32), __aligned__(32)));
	typedef long long	__m256i __attribute__((__vector_size__(32), __aligned__(32)));

	typedef float		__m512  __attribute__((__vector_size__(64), __aligned__(64)));
	typedef double		__m512d __attribute__((__vector_size__(64), __aligned__(64)));
	typedef long long	__m512i __attribute__((__vector_size__(64), __aligned__(64)));
	*/
	#define VLA_NEW(name,type,size) type name##VLA[size]; type* name = (type*)name##VLA
	#define VLA_DELETE(name)
	#define forceinline __attribute__((always_inline))
#else
	#define VLA_NEW(name,type,size) type* name = new type[size]
	#define VLA_DELETE(name) { if (name) delete[] (name); (name) = NULL; }
	#define forceinline __forceinline
#endif
	
#ifdef _FMA
	#define _mm512_fmaddx_pd(a,b,c)		_mm512_fmadd_pd(a,b,c)
	#define _mm512_fmaddx_ps(a,b,c)		_mm512_fmadd_ps(a,b,c)
	#define _mm256_fmaddx_pd(a,b,c)		_mm256_fmadd_pd(a,b,c)
	#define _mm256_fmaddx_ps(a,b,c)		_mm256_fmadd_ps(a,b,c)
	#define _mm_fmaddx_pd(a,b,c)		_mm_add_pd(_mm_mul_pd(a,b),c)
	#define _mm_fmaddx_ps(a,b,c)		_mm_add_ps(_mm_mul_ps(a,b),c)
	#define _neo_fmaddx_pd(a,b,c)		vfmaq_f64(a,b,c)
	#define _neo_fmaddx_ps(a,b,c)		vfmaq_f32(a,b,c)
#else
	#define _mm512_fmaddx_pd(a,b,c)		_mm512_add_pd(_mm512_mul_pd(a,b),c)
	#define _mm512_fmaddx_ps(a,b,c)		_mm512_add_ps(_mm512_mul_ps(a,b),c)
	#define _mm256_fmaddx_pd(a,b,c)		_mm256_add_pd(_mm256_mul_pd(a,b),c)
	#define _mm256_fmaddx_ps(a,b,c)		_mm256_add_ps(_mm256_mul_ps(a,b),c)
	#define _mm_fmaddx_pd(a,b,c)		_mm_add_pd(_mm_mul_pd(a,b),c)
	#define _mm_fmaddx_ps(a,b,c)		_mm_add_ps(_mm_mul_ps(a,b),c)
    #define _neo_fmaddx_pd(a,b,c)		vaddq_f64(vmulq_f64(a,b),c)
	#define _neo_fmaddx_ps(a,b,c)		vaddq_f32(vmulq_f32(a,b),c)
#endif

#define simd_f64(x,y)		((double*)&(x))[y]
#define simd_f32(x,y)		((float*)&(x))[y]
#define simd_u256(x,y)		((__m256i*)&(x))[y]
#define simd_u128(x,y)		((__m128i*)&(x))[y]
#define simd_u64(x,y)		((uint64*)&(x))[y]
#define simd_u32(x,y)		((uint*)&(x))[y]
#define simd_u16(x,y)		((ushort*)&(x))[y]
#define simd_u8(x,y)		((byte*)&(x))[y]
#define simd_i256(x,y)		((__m256i*)&(x))[y]
#define simd_i128(x,y)		((__m128i*)&(x))[y]
#define simd_i64(x,y)		((int64*)&(x))[y]
#define simd_i32(x,y)		((int*)&(x))[y]
#define simd_i16(x,y)		((short*)&(x))[y]
#define simd_i8(x,y)		((char*)&(x))[y]

#define simp_f64(x,y)		((double*)(x))[y]
#define simp_f32(x,y)		((float*)(x))[y]
#define simp_u256(x,y)		((__m256i*)(x))[y]
#define simp_u128(x,y)		((__m128i*)(x))[y]
#define simp_u64(x,y)		((uint64*)(x))[y]
#define simp_u32(x,y)		((uint*)(x))[y]
#define simp_u16(x,y)		((ushort*)(x))[y]
#define simp_u8(x,y)		((byte*)(x))[y]
#define simp_i256(x,y)		((__m256i*)(x))[y]
#define simp_i128(x,y)		((__m128i*)(x))[y]
#define simp_i64(x,y)		((int64*)(x))[y]
#define simp_i32(x,y)		((int*)(x))[y]
#define simp_i16(x,y)		((short*)(x))[y]
#define simp_i8(x,y)		((char*)(x))[y]
	
#define E512_512 1
#define E512_256 2
#define E512_128 4
#define E256_128 2

#define E512D 8
#define E512F 16
#define E512S 32
#define E512B 64

#define E256D 4
#define E256F 8
#define E256S 16
#define E256B 32

#define E128D 2
#define E128F 4
#define E128S 8
#define E128B 16

#define GetLoc(x)			(slocus[(x)])				//(useslocus ? slocus[l] : locus[l])
#define GetLocPos(x)		(locus_pos[(x)])			//(useslocus ? locus_pos[l] : locus[l].pos)
#define GetLocId(x)			(locus_id[(x)])				//(useslocus ? locus_id[l] : locus[l].id)
#define GetLocK(x)			((x) < lend ? slocus[(x)].k : 0)
#define GetLocTabDiff(x)	((x) < lend ? (int64)slocus[(x)].GetGtab() - gtab_base : 0)	

#define REDUCE(x) for (int KK = sizeof(x) / sizeof(x[0]) / 2; KK >= 1; KK >>= 1) REP(KK)

#if defined(__clang__)
	#define STRINGIFY2(x) #x
	#define STRINGIFY1(x) STRINGIFY2(x)
	#define UNROLL2(x) unroll x
	#define UNROLL1(N)			_Pragma(STRINGIFY1(UNROLL2(N))) for (int kk = 0; kk < N; ++kk)

	#define UNROLL(N)			UNROLL1(N)
	#define UNROLLHEAD(N)		_Pragma(STRINGIFY1(UNROLL2(N)))
	#define VECTORIZE			_Pragma("clang loop vectorize(enable)")
#elif  defined(__GNUC__)
	#define STRINGIFY2(x) #x
	#define STRINGIFY1(x) STRINGIFY2(x)
	#define UNROLL2(x) GCC unroll x
	#define UNROLL1(N)			_Pragma(STRINGIFY1(UNROLL2(N))) for (int kk = 0; kk < N; ++kk)

	#define UNROLL(N)			UNROLL1(N)
	#define UNROLLHEAD(N)		_Pragma(STRINGIFY1(UNROLL2(N)))
	#define VECTORIZE			_Pragma("GCC ivdep")
#elif  defined(__CUDACC__)
	#define STRINGIFY2(x) #x
	#define STRINGIFY1(x) STRINGIFY2(x)
	#define UNROLL2(x) unroll x
	#define UNROLL1(N)			_Pragma(STRINGIFY1(UNROLL2(N))) for (int kk = 0; kk < N; ++kk)

	#define UNROLL(N)			UNROLL1(N)
	#define UNROLLHEAD(N)		_Pragma(STRINGIFY1(UNROLL2(N)))
	#define VECTORIZE			
#else
	#define UNROLL(x)			for (int kk = 0; kk < (x); ++kk)
	#define UNROLLHEAD(N)		
	#define VECTORIZE			_Pragma("loop(ivdep)")
#endif

#define REP_A1(expr) expr(0);
#define REP_A2(expr) expr(0);expr(1);
#define REP_A3(expr) expr(0);expr(1);expr(2);
#define REP_A4(expr) expr(0);expr(1);expr(2);expr(3);
#define REP_A5(expr) expr(0);expr(1);expr(2);expr(3);expr(4);
#define REP_A6(expr) expr(0);expr(1);expr(2);expr(3);expr(4);expr(5);
#define REP_A7(expr) expr(0);expr(1);expr(2);expr(3);expr(4);expr(5);expr(6);
#define REP_A8(expr) expr(0);expr(1);expr(2);expr(3);expr(4);expr(5);expr(6);expr(7);
#define REP_A9(expr) expr(0);expr(1);expr(2);expr(3);expr(4);expr(5);expr(6);expr(7);expr(8);
#define REP_A10(expr) expr(0);expr(1);expr(2);expr(3);expr(4);expr(5);expr(6);expr(7);expr(8);expr(9);
#define REP_A11(expr) expr(0);expr(1);expr(2);expr(3);expr(4);expr(5);expr(6);expr(7);expr(8);expr(9);expr(10);
#define REP_A12(expr) expr(0);expr(1);expr(2);expr(3);expr(4);expr(5);expr(6);expr(7);expr(8);expr(9);expr(10);expr(11);
#define REP_A13(expr) expr(0);expr(1);expr(2);expr(3);expr(4);expr(5);expr(6);expr(7);expr(8);expr(9);expr(10);expr(11);expr(12);
#define REP_A14(expr) expr(0);expr(1);expr(2);expr(3);expr(4);expr(5);expr(6);expr(7);expr(8);expr(9);expr(10);expr(11);expr(12);expr(13);
#define REP_A15(expr) expr(0);expr(1);expr(2);expr(3);expr(4);expr(5);expr(6);expr(7);expr(8);expr(9);expr(10);expr(11);expr(12);expr(13);expr(14);
#define REP_A16(expr) expr(0);expr(1);expr(2);expr(3);expr(4);expr(5);expr(6);expr(7);expr(8);expr(9);expr(10);expr(11);expr(12);expr(13);expr(14);expr(15);

#define REP_B1(expr,x) expr(x,0);
#define REP_B2(expr,x) expr(x,0);expr(x,1);
#define REP_B3(expr,x) expr(x,0);expr(x,1);expr(x,2);
#define REP_B4(expr,x) expr(x,0);expr(x,1);expr(x,2);expr(x,3);
#define REP_B5(expr,x) expr(x,0);expr(x,1);expr(x,2);expr(x,3);expr(x,4);
#define REP_B6(expr,x) expr(x,0);expr(x,1);expr(x,2);expr(x,3);expr(x,4);expr(x,5);
#define REP_B7(expr,x) expr(x,0);expr(x,1);expr(x,2);expr(x,3);expr(x,4);expr(x,5);expr(x,6);
#define REP_B8(expr,x) expr(x,0);expr(x,1);expr(x,2);expr(x,3);expr(x,4);expr(x,5);expr(x,6);expr(x,7);
#define REP_B9(expr,x) expr(x,0);expr(x,1);expr(x,2);expr(x,3);expr(x,4);expr(x,5);expr(x,6);expr(x,7);expr(x,8);
#define REP_B10(expr,x) expr(x,0);expr(x,1);expr(x,2);expr(x,3);expr(x,4);expr(x,5);expr(x,6);expr(x,7);expr(x,8);expr(x,9);
#define REP_B11(expr,x) expr(x,0);expr(x,1);expr(x,2);expr(x,3);expr(x,4);expr(x,5);expr(x,6);expr(x,7);expr(x,8);expr(x,9);expr(x,10);
#define REP_B12(expr,x) expr(x,0);expr(x,1);expr(x,2);expr(x,3);expr(x,4);expr(x,5);expr(x,6);expr(x,7);expr(x,8);expr(x,9);expr(x,10);expr(x,11);
#define REP_B13(expr,x) expr(x,0);expr(x,1);expr(x,2);expr(x,3);expr(x,4);expr(x,5);expr(x,6);expr(x,7);expr(x,8);expr(x,9);expr(x,10);expr(x,11);expr(x,12);
#define REP_B14(expr,x) expr(x,0);expr(x,1);expr(x,2);expr(x,3);expr(x,4);expr(x,5);expr(x,6);expr(x,7);expr(x,8);expr(x,9);expr(x,10);expr(x,11);expr(x,12);expr(x,13);
#define REP_B15(expr,x) expr(x,0);expr(x,1);expr(x,2);expr(x,3);expr(x,4);expr(x,5);expr(x,6);expr(x,7);expr(x,8);expr(x,9);expr(x,10);expr(x,11);expr(x,12);expr(x,13);expr(x,14);
#define REP_B16(expr,x) expr(x,0);expr(x,1);expr(x,2);expr(x,3);expr(x,4);expr(x,5);expr(x,6);expr(x,7);expr(x,8);expr(x,9);expr(x,10);expr(x,11);expr(x,12);expr(x,13);expr(x,14);expr(x,15);

#define REP_C1(expr,x) expr(0,x);
#define REP_C2(expr,x) expr(0,x);expr(1,x);
#define REP_C3(expr,x) expr(0,x);expr(1,x);expr(2,x);
#define REP_C4(expr,x) expr(0,x);expr(1,x);expr(2,x);expr(3,x);
#define REP_C5(expr,x) expr(0,x);expr(1,x);expr(2,x);expr(3,x);expr(4,x);
#define REP_C6(expr,x) expr(0,x);expr(1,x);expr(2,x);expr(3,x);expr(4,x);expr(5,x);
#define REP_C7(expr,x) expr(0,x);expr(1,x);expr(2,x);expr(3,x);expr(4,x);expr(5,x);expr(6,x);
#define REP_C8(expr,x) expr(0,x);expr(1,x);expr(2,x);expr(3,x);expr(4,x);expr(5,x);expr(6,x);expr(7,x);
#define REP_C9(expr,x) expr(0,x);expr(1,x);expr(2,x);expr(3,x);expr(4,x);expr(5,x);expr(6,x);expr(7,x);expr(8,x);
#define REP_C10(expr,x) expr(0,x);expr(1,x);expr(2,x);expr(3,x);expr(4,x);expr(5,x);expr(6,x);expr(7,x);expr(8,x);expr(9,x);
#define REP_C11(expr,x) expr(0,x);expr(1,x);expr(2,x);expr(3,x);expr(4,x);expr(5,x);expr(6,x);expr(7,x);expr(8,x);expr(9,x);expr(10,x);
#define REP_C12(expr,x) expr(0,x);expr(1,x);expr(2,x);expr(3,x);expr(4,x);expr(5,x);expr(6,x);expr(7,x);expr(8,x);expr(9,x);expr(10,x);expr(11,x);
#define REP_C13(expr,x) expr(0,x);expr(1,x);expr(2,x);expr(3,x);expr(4,x);expr(5,x);expr(6,x);expr(7,x);expr(8,x);expr(9,x);expr(10,x);expr(11,x);expr(12,x);
#define REP_C14(expr,x) expr(0,x);expr(1,x);expr(2,x);expr(3,x);expr(4,x);expr(5,x);expr(6,x);expr(7,x);expr(8,x);expr(9,x);expr(10,x);expr(11,x);expr(12,x);expr(13,x);
#define REP_C15(expr,x) expr(0,x);expr(1,x);expr(2,x);expr(3,x);expr(4,x);expr(5,x);expr(6,x);expr(7,x);expr(8,x);expr(9,x);expr(10,x);expr(11,x);expr(12,x);expr(13,x);expr(14,x);
#define REP_C16(expr,x) expr(0,x);expr(1,x);expr(2,x);expr(3,x);expr(4,x);expr(5,x);expr(6,x);expr(7,x);expr(8,x);expr(9,x);expr(10,x);expr(11,x);expr(12,x);expr(13,x);expr(14,x);expr(15,x);

#define REP_D1(fun,expr) fun(expr,0);
#define REP_D2(fun,expr) fun(expr,0);fun(expr,1);
#define REP_D3(fun,expr) fun(expr,0);fun(expr,1);fun(expr,2);
#define REP_D4(fun,expr) fun(expr,0);fun(expr,1);fun(expr,2);fun(expr,3);
#define REP_D5(fun,expr) fun(expr,0);fun(expr,1);fun(expr,2);fun(expr,3);fun(expr,4);
#define REP_D6(fun,expr) fun(expr,0);fun(expr,1);fun(expr,2);fun(expr,3);fun(expr,4);fun(expr,5);
#define REP_D7(fun,expr) fun(expr,0);fun(expr,1);fun(expr,2);fun(expr,3);fun(expr,4);fun(expr,5);fun(expr,6);
#define REP_D8(fun,expr) fun(expr,0);fun(expr,1);fun(expr,2);fun(expr,3);fun(expr,4);fun(expr,5);fun(expr,6);fun(expr,7);
#define REP_D9(fun,expr) fun(expr,0);fun(expr,1);fun(expr,2);fun(expr,3);fun(expr,4);fun(expr,5);fun(expr,6);fun(expr,7);fun(expr,8);
#define REP_D10(fun,expr) fun(expr,0);fun(expr,1);fun(expr,2);fun(expr,3);fun(expr,4);fun(expr,5);fun(expr,6);fun(expr,7);fun(expr,8);fun(expr,9);
#define REP_D11(fun,expr) fun(expr,0);fun(expr,1);fun(expr,2);fun(expr,3);fun(expr,4);fun(expr,5);fun(expr,6);fun(expr,7);fun(expr,8);fun(expr,9);fun(expr,10);
#define REP_D12(fun,expr) fun(expr,0);fun(expr,1);fun(expr,2);fun(expr,3);fun(expr,4);fun(expr,5);fun(expr,6);fun(expr,7);fun(expr,8);fun(expr,9);fun(expr,10);fun(expr,11);
#define REP_D13(fun,expr) fun(expr,0);fun(expr,1);fun(expr,2);fun(expr,3);fun(expr,4);fun(expr,5);fun(expr,6);fun(expr,7);fun(expr,8);fun(expr,9);fun(expr,10);fun(expr,11);fun(expr,12);
#define REP_D14(fun,expr) fun(expr,0);fun(expr,1);fun(expr,2);fun(expr,3);fun(expr,4);fun(expr,5);fun(expr,6);fun(expr,7);fun(expr,8);fun(expr,9);fun(expr,10);fun(expr,11);fun(expr,12);fun(expr,13);
#define REP_D15(fun,expr) fun(expr,0);fun(expr,1);fun(expr,2);fun(expr,3);fun(expr,4);fun(expr,5);fun(expr,6);fun(expr,7);fun(expr,8);fun(expr,9);fun(expr,10);fun(expr,11);fun(expr,12);fun(expr,13);fun(expr,14);
#define REP_D16(fun,expr) fun(expr,0);fun(expr,1);fun(expr,2);fun(expr,3);fun(expr,4);fun(expr,5);fun(expr,6);fun(expr,7);fun(expr,8);fun(expr,9);fun(expr,10);fun(expr,11);fun(expr,12);fun(expr,13);fun(expr,14);fun(expr,15);

#define REP_E1(fun,expr) fun##1(expr,0);
#define REP_E2(fun,expr) fun##1(expr,0);fun##2(expr,1);
#define REP_E3(fun,expr) fun##1(expr,0);fun##2(expr,1);fun##3(expr,2);
#define REP_E4(fun,expr) fun##1(expr,0);fun##2(expr,1);fun##3(expr,2);fun##4(expr,3);
#define REP_E5(fun,expr) fun##1(expr,0);fun##2(expr,1);fun##3(expr,2);fun##4(expr,3);fun##5(expr,4);
#define REP_E6(fun,expr) fun##1(expr,0);fun##2(expr,1);fun##3(expr,2);fun##4(expr,3);fun##5(expr,4);fun##6(expr,5);
#define REP_E7(fun,expr) fun##1(expr,0);fun##2(expr,1);fun##3(expr,2);fun##4(expr,3);fun##5(expr,4);fun##6(expr,5);fun##7(expr,6);
#define REP_E8(fun,expr) fun##1(expr,0);fun##2(expr,1);fun##3(expr,2);fun##4(expr,3);fun##5(expr,4);fun##6(expr,5);fun##7(expr,6);fun##8(expr,7);
#define REP_E9(fun,expr) fun##1(expr,0);fun##2(expr,1);fun##3(expr,2);fun##4(expr,3);fun##5(expr,4);fun##6(expr,5);fun##7(expr,6);fun##8(expr,7);fun##9(expr,8);
#define REP_E10(fun,expr) fun##1(expr,0);fun##2(expr,1);fun##3(expr,2);fun##4(expr,3);fun##5(expr,4);fun##6(expr,5);fun##7(expr,6);fun##8(expr,7);fun##9(expr,8);fun##10(expr,9);
#define REP_E11(fun,expr) fun##1(expr,0);fun##2(expr,1);fun##3(expr,2);fun##4(expr,3);fun##5(expr,4);fun##6(expr,5);fun##7(expr,6);fun##8(expr,7);fun##9(expr,8);fun##10(expr,9);fun##11(expr,10);
#define REP_E12(fun,expr) fun##1(expr,0);fun##2(expr,1);fun##3(expr,2);fun##4(expr,3);fun##5(expr,4);fun##6(expr,5);fun##7(expr,6);fun##8(expr,7);fun##9(expr,8);fun##10(expr,9);fun##11(expr,10);fun##12(expr,11);
#define REP_E13(fun,expr) fun##1(expr,0);fun##2(expr,1);fun##3(expr,2);fun##4(expr,3);fun##5(expr,4);fun##6(expr,5);fun##7(expr,6);fun##8(expr,7);fun##9(expr,8);fun##10(expr,9);fun##11(expr,10);fun##12(expr,11);fun##13(expr,12);
#define REP_E14(fun,expr) fun##1(expr,0);fun##2(expr,1);fun##3(expr,2);fun##4(expr,3);fun##5(expr,4);fun##6(expr,5);fun##7(expr,6);fun##8(expr,7);fun##9(expr,8);fun##10(expr,9);fun##11(expr,10);fun##12(expr,11);fun##13(expr,12);fun##14(expr,13);
#define REP_E15(fun,expr) fun##1(expr,0);fun##2(expr,1);fun##3(expr,2);fun##4(expr,3);fun##5(expr,4);fun##6(expr,5);fun##7(expr,6);fun##8(expr,7);fun##9(expr,8);fun##10(expr,9);fun##11(expr,10);fun##12(expr,11);fun##13(expr,12);fun##14(expr,13);fun##15(expr,14);
#define REP_E16(fun,expr) fun##1(expr,0);fun##2(expr,1);fun##3(expr,2);fun##4(expr,3);fun##5(expr,4);fun##6(expr,5);fun##7(expr,6);fun##8(expr,7);fun##9(expr,8);fun##10(expr,9);fun##11(expr,10);fun##12(expr,11);fun##13(expr,12);fun##14(expr,13);fun##15(expr,14);fun##16(expr,15);

#define REP_F1(expr) 
#define REP_F2(expr) case 1:expr(0);
#define REP_F3(expr) case 2:expr(1);case 1:expr(0);
#define REP_F4(expr) case 3:expr(2);case 2:expr(1);case 1:expr(0);
#define REP_F5(expr) case 4:expr(3);case 3:expr(2);case 2:expr(1);case 1:expr(0);
#define REP_F6(expr) case 5:expr(4);case 4:expr(3);case 3:expr(2);case 2:expr(1);case 1:expr(0);
#define REP_F7(expr) case 6:expr(5);case 5:expr(4);case 4:expr(3);case 3:expr(2);case 2:expr(1);case 1:expr(0);
#define REP_F8(expr) case 7:expr(6);case 6:expr(5);case 5:expr(4);case 4:expr(3);case 3:expr(2);case 2:expr(1);case 1:expr(0);
#define REP_F9(expr) case 8:expr(7);case 7:expr(6);case 6:expr(5);case 5:expr(4);case 4:expr(3);case 3:expr(2);case 2:expr(1);case 1:expr(0);
#define REP_F10(expr) case 9:expr(8);case 8:expr(7);case 7:expr(6);case 6:expr(5);case 5:expr(4);case 4:expr(3);case 3:expr(2);case 2:expr(1);case 1:expr(0);
#define REP_F11(expr) case 10:expr(9);case 9:expr(8);case 8:expr(7);case 7:expr(6);case 6:expr(5);case 5:expr(4);case 4:expr(3);case 3:expr(2);case 2:expr(1);case 1:expr(0);
#define REP_F12(expr) case 11:expr(10);case 10:expr(9);case 9:expr(8);case 8:expr(7);case 7:expr(6);case 6:expr(5);case 5:expr(4);case 4:expr(3);case 3:expr(2);case 2:expr(1);case 1:expr(0);
#define REP_F13(expr) case 12:expr(11);case 11:expr(10);case 10:expr(9);case 9:expr(8);case 8:expr(7);case 7:expr(6);case 6:expr(5);case 5:expr(4);case 4:expr(3);case 3:expr(2);case 2:expr(1);case 1:expr(0);
#define REP_F14(expr) case 13:expr(12);case 12:expr(11);case 11:expr(10);case 10:expr(9);case 9:expr(8);case 8:expr(7);case 7:expr(6);case 6:expr(5);case 5:expr(4);case 4:expr(3);case 3:expr(2);case 2:expr(1);case 1:expr(0);
#define REP_F15(expr) case 14:expr(13);case 13:expr(12);case 12:expr(11);case 11:expr(10);case 10:expr(9);case 9:expr(8);case 8:expr(7);case 7:expr(6);case 6:expr(5);case 5:expr(4);case 4:expr(3);case 3:expr(2);case 2:expr(1);case 1:expr(0);
#define REP_F16(expr) case 15:expr(14);case 14:expr(13);case 13:expr(12);case 12:expr(11);case 11:expr(10);case 10:expr(9);case 9:expr(8);case 8:expr(7);case 7:expr(6);case 6:expr(5);case 5:expr(4);case 4:expr(3);case 3:expr(2);case 2:expr(1);case 1:expr(0);

#define REP_G1(fun,expr)  
#define REP_G2(fun,expr)  case 1:fun(expr,0);
#define REP_G3(fun,expr)  case 2:fun(expr,1);case 1:fun(expr,0);
#define REP_G4(fun,expr)  case 3:fun(expr,2);case 2:fun(expr,1);case 1:fun(expr,0);
#define REP_G5(fun,expr)  case 4:fun(expr,3);case 3:fun(expr,2);case 2:fun(expr,1);case 1:fun(expr,0);
#define REP_G6(fun,expr)  case 5:fun(expr,4);case 4:fun(expr,3);case 3:fun(expr,2);case 2:fun(expr,1);case 1:fun(expr,0);
#define REP_G7(fun,expr)  case 6:fun(expr,5);case 5:fun(expr,4);case 4:fun(expr,3);case 3:fun(expr,2);case 2:fun(expr,1);case 1:fun(expr,0);
#define REP_G8(fun,expr)  case 7:fun(expr,6);case 6:fun(expr,5);case 5:fun(expr,4);case 4:fun(expr,3);case 3:fun(expr,2);case 2:fun(expr,1);case 1:fun(expr,0);
#define REP_G9(fun,expr)  case 8:fun(expr,7);case 7:fun(expr,6);case 6:fun(expr,5);case 5:fun(expr,4);case 4:fun(expr,3);case 3:fun(expr,2);case 2:fun(expr,1);case 1:fun(expr,0);
#define REP_G10(fun,expr)  case 9:fun(expr,8);case 8:fun(expr,7);case 7:fun(expr,6);case 6:fun(expr,5);case 5:fun(expr,4);case 4:fun(expr,3);case 3:fun(expr,2);case 2:fun(expr,1);case 1:fun(expr,0);
#define REP_G11(fun,expr)  case 10:fun(expr,9);case 9:fun(expr,8);case 8:fun(expr,7);case 7:fun(expr,6);case 6:fun(expr,5);case 5:fun(expr,4);case 4:fun(expr,3);case 3:fun(expr,2);case 2:fun(expr,1);case 1:fun(expr,0);
#define REP_G12(fun,expr)  case 11:fun(expr,10);case 10:fun(expr,9);case 9:fun(expr,8);case 8:fun(expr,7);case 7:fun(expr,6);case 6:fun(expr,5);case 5:fun(expr,4);case 4:fun(expr,3);case 3:fun(expr,2);case 2:fun(expr,1);case 1:fun(expr,0);
#define REP_G13(fun,expr)  case 12:fun(expr,11);case 11:fun(expr,10);case 10:fun(expr,9);case 9:fun(expr,8);case 8:fun(expr,7);case 7:fun(expr,6);case 6:fun(expr,5);case 5:fun(expr,4);case 4:fun(expr,3);case 3:fun(expr,2);case 2:fun(expr,1);case 1:fun(expr,0);
#define REP_G14(fun,expr)  case 13:fun(expr,12);case 12:fun(expr,11);case 11:fun(expr,10);case 10:fun(expr,9);case 9:fun(expr,8);case 8:fun(expr,7);case 7:fun(expr,6);case 6:fun(expr,5);case 5:fun(expr,4);case 4:fun(expr,3);case 3:fun(expr,2);case 2:fun(expr,1);case 1:fun(expr,0);
#define REP_G15(fun,expr)  case 14:fun(expr,13);case 13:fun(expr,12);case 12:fun(expr,11);case 11:fun(expr,10);case 10:fun(expr,9);case 9:fun(expr,8);case 8:fun(expr,7);case 7:fun(expr,6);case 6:fun(expr,5);case 5:fun(expr,4);case 4:fun(expr,3);case 3:fun(expr,2);case 2:fun(expr,1);case 1:fun(expr,0);
#define REP_G16(fun,expr)  case 15:fun(expr,14);case 14:fun(expr,13);case 13:fun(expr,12);case 12:fun(expr,11);case 11:fun(expr,10);case 10:fun(expr,9);case 9:fun(expr,8);case 8:fun(expr,7);case 7:fun(expr,6);case 6:fun(expr,5);case 5:fun(expr,4);case 4:fun(expr,3);case 3:fun(expr,2);case 2:fun(expr,1);case 1:fun(expr,0);

#define REP_H1(fun,expr)  break;
#define REP_H2(fun,expr)  case 1:fun##1(expr,0);break;
#define REP_H3(fun,expr)  case 2:fun##2(expr,1);case 1:fun##1(expr,0);break;
#define REP_H4(fun,expr)  case 3:fun##3(expr,2);case 2:fun##2(expr,1);case 1:fun##1(expr,0);break;
#define REP_H5(fun,expr)  case 4:fun##4(expr,3);case 3:fun##3(expr,2);case 2:fun##2(expr,1);case 1:fun##1(expr,0);break;
#define REP_H6(fun,expr)  case 5:fun##5(expr,4);case 4:fun##4(expr,3);case 3:fun##3(expr,2);case 2:fun##2(expr,1);case 1:fun##1(expr,0);break;
#define REP_H7(fun,expr)  case 6:fun##6(expr,5);case 5:fun##5(expr,4);case 4:fun##4(expr,3);case 3:fun##3(expr,2);case 2:fun##2(expr,1);case 1:fun##1(expr,0);break;
#define REP_H8(fun,expr)  case 7:fun##7(expr,6);case 6:fun##6(expr,5);case 5:fun##5(expr,4);case 4:fun##4(expr,3);case 3:fun##3(expr,2);case 2:fun##2(expr,1);case 1:fun##1(expr,0);break;
#define REP_H9(fun,expr)  case 8:fun##8(expr,7);case 7:fun##7(expr,6);case 6:fun##6(expr,5);case 5:fun##5(expr,4);case 4:fun##4(expr,3);case 3:fun##3(expr,2);case 2:fun##2(expr,1);case 1:fun##1(expr,0);break;
#define REP_H10(fun,expr)  case 9:fun##9(expr,8);case 8:fun##8(expr,7);case 7:fun##7(expr,6);case 6:fun##6(expr,5);case 5:fun##5(expr,4);case 4:fun##4(expr,3);case 3:fun##3(expr,2);case 2:fun##2(expr,1);case 1:fun##1(expr,0);break;
#define REP_H11(fun,expr)  case 10:fun##10(expr,9);case 9:fun##9(expr,8);case 8:fun##8(expr,7);case 7:fun##7(expr,6);case 6:fun##6(expr,5);case 5:fun##5(expr,4);case 4:fun##4(expr,3);case 3:fun##3(expr,2);case 2:fun##2(expr,1);case 1:fun##1(expr,0);break;
#define REP_H12(fun,expr)  case 11:fun##11(expr,10);case 10:fun##10(expr,9);case 9:fun##9(expr,8);case 8:fun##8(expr,7);case 7:fun##7(expr,6);case 6:fun##6(expr,5);case 5:fun##5(expr,4);case 4:fun##4(expr,3);case 3:fun##3(expr,2);case 2:fun##2(expr,1);case 1:fun##1(expr,0);break;
#define REP_H13(fun,expr)  case 12:fun##12(expr,11);case 11:fun##11(expr,10);case 10:fun##10(expr,9);case 9:fun##9(expr,8);case 8:fun##8(expr,7);case 7:fun##7(expr,6);case 6:fun##6(expr,5);case 5:fun##5(expr,4);case 4:fun##4(expr,3);case 3:fun##3(expr,2);case 2:fun##2(expr,1);case 1:fun##1(expr,0);break;
#define REP_H14(fun,expr)  case 13:fun##13(expr,12);case 12:fun##12(expr,11);case 11:fun##11(expr,10);case 10:fun##10(expr,9);case 9:fun##9(expr,8);case 8:fun##8(expr,7);case 7:fun##7(expr,6);case 6:fun##6(expr,5);case 5:fun##5(expr,4);case 4:fun##4(expr,3);case 3:fun##3(expr,2);case 2:fun##2(expr,1);case 1:fun##1(expr,0);break;
#define REP_H15(fun,expr)  case 14:fun##14(expr,13);case 13:fun##13(expr,12);case 12:fun##12(expr,11);case 11:fun##11(expr,10);case 10:fun##10(expr,9);case 9:fun##9(expr,8);case 8:fun##8(expr,7);case 7:fun##7(expr,6);case 6:fun##6(expr,5);case 5:fun##5(expr,4);case 4:fun##4(expr,3);case 3:fun##3(expr,2);case 2:fun##2(expr,1);case 1:fun##1(expr,0);break;
#define REP_H16(fun,expr)  case 15:fun##15(expr,14);case 14:fun##14(expr,13);case 13:fun##13(expr,12);case 12:fun##12(expr,11);case 11:fun##11(expr,10);case 10:fun##10(expr,9);case 9:fun##9(expr,8);case 8:fun##8(expr,7);case 7:fun##7(expr,6);case 6:fun##6(expr,5);case 5:fun##5(expr,4);case 4:fun##4(expr,3);case 3:fun##3(expr,2);case 2:fun##2(expr,1);case 1:fun##1(expr,0);break;

#define REP_I1(fun,expr,x)  break;
#define REP_I2(fun,expr,x)  case 0x##x##1:fun##x(expr,0);break;
#define REP_I3(fun,expr,x)  case 0x##x##2:fun##x(expr,1);case 0x##x##1:fun##x(expr,0);break;
#define REP_I4(fun,expr,x)  case 0x##x##3:fun##x(expr,2);case 0x##x##2:fun##x(expr,1);case 0x##x##1:fun##x(expr,0);break;
#define REP_I5(fun,expr,x)  case 0x##x##4:fun##x(expr,3);case 0x##x##3:fun##x(expr,2);case 0x##x##2:fun##x(expr,1);case 0x##x##1:fun##x(expr,0);break;
#define REP_I6(fun,expr,x)  case 0x##x##5:fun##x(expr,4);case 0x##x##4:fun##x(expr,3);case 0x##x##3:fun##x(expr,2);case 0x##x##2:fun##x(expr,1);case 0x##x##1:fun##x(expr,0);break;
#define REP_I7(fun,expr,x)  case 0x##x##6:fun##x(expr,5);case 0x##x##5:fun##x(expr,4);case 0x##x##4:fun##x(expr,3);case 0x##x##3:fun##x(expr,2);case 0x##x##2:fun##x(expr,1);case 0x##x##1:fun##x(expr,0);break;
#define REP_I8(fun,expr,x)  case 0x##x##7:fun##x(expr,6);case 0x##x##6:fun##x(expr,5);case 0x##x##5:fun##x(expr,4);case 0x##x##4:fun##x(expr,3);case 0x##x##3:fun##x(expr,2);case 0x##x##2:fun##x(expr,1);case 0x##x##1:fun##x(expr,0);break;
#define REP_I9(fun,expr,x)  case 0x##x##8:fun##x(expr,7);case 0x##x##7:fun##x(expr,6);case 0x##x##6:fun##x(expr,5);case 0x##x##5:fun##x(expr,4);case 0x##x##4:fun##x(expr,3);case 0x##x##3:fun##x(expr,2);case 0x##x##2:fun##x(expr,1);case 0x##x##1:fun##x(expr,0);break;
#define REP_I10(fun,expr,x)  case 0x##x##9:fun##x(expr,8);case 0x##x##8:fun##x(expr,7);case 0x##x##7:fun##x(expr,6);case 0x##x##6:fun##x(expr,5);case 0x##x##5:fun##x(expr,4);case 0x##x##4:fun##x(expr,3);case 0x##x##3:fun##x(expr,2);case 0x##x##2:fun##x(expr,1);case 0x##x##1:fun##x(expr,0);break;
#define REP_I11(fun,expr,x)  case 0x##x##10:fun##x(expr,9);case 0x##x##9:fun##x(expr,8);case 0x##x##8:fun##x(expr,7);case 0x##x##7:fun##x(expr,6);case 0x##x##6:fun##x(expr,5);case 0x##x##5:fun##x(expr,4);case 0x##x##4:fun##x(expr,3);case 0x##x##3:fun##x(expr,2);case 0x##x##2:fun##x(expr,1);case 0x##x##1:fun##x(expr,0);break;
#define REP_I12(fun,expr,x)  case 0x##x##11:fun##x(expr,10);case 0x##x##10:fun##x(expr,9);case 0x##x##9:fun##x(expr,8);case 0x##x##8:fun##x(expr,7);case 0x##x##7:fun##x(expr,6);case 0x##x##6:fun##x(expr,5);case 0x##x##5:fun##x(expr,4);case 0x##x##4:fun##x(expr,3);case 0x##x##3:fun##x(expr,2);case 0x##x##2:fun##x(expr,1);case 0x##x##1:fun##x(expr,0);break;
#define REP_I13(fun,expr,x)  case 0x##x##12:fun##x(expr,11);case 0x##x##11:fun##x(expr,10);case 0x##x##10:fun##x(expr,9);case 0x##x##9:fun##x(expr,8);case 0x##x##8:fun##x(expr,7);case 0x##x##7:fun##x(expr,6);case 0x##x##6:fun##x(expr,5);case 0x##x##5:fun##x(expr,4);case 0x##x##4:fun##x(expr,3);case 0x##x##3:fun##x(expr,2);case 0x##x##2:fun##x(expr,1);case 0x##x##1:fun##x(expr,0);break;
#define REP_I14(fun,expr,x)  case 0x##x##13:fun##x(expr,12);case 0x##x##12:fun##x(expr,11);case 0x##x##11:fun##x(expr,10);case 0x##x##10:fun##x(expr,9);case 0x##x##9:fun##x(expr,8);case 0x##x##8:fun##x(expr,7);case 0x##x##7:fun##x(expr,6);case 0x##x##6:fun##x(expr,5);case 0x##x##5:fun##x(expr,4);case 0x##x##4:fun##x(expr,3);case 0x##x##3:fun##x(expr,2);case 0x##x##2:fun##x(expr,1);case 0x##x##1:fun##x(expr,0);break;
#define REP_I15(fun,expr,x)  case 0x##x##14:fun##x(expr,13);case 0x##x##13:fun##x(expr,12);case 0x##x##12:fun##x(expr,11);case 0x##x##11:fun##x(expr,10);case 0x##x##10:fun##x(expr,9);case 0x##x##9:fun##x(expr,8);case 0x##x##8:fun##x(expr,7);case 0x##x##7:fun##x(expr,6);case 0x##x##6:fun##x(expr,5);case 0x##x##5:fun##x(expr,4);case 0x##x##4:fun##x(expr,3);case 0x##x##3:fun##x(expr,2);case 0x##x##2:fun##x(expr,1);case 0x##x##1:fun##x(expr,0);break;
#define REP_I16(fun,expr,x)  case 0x##x##15:fun##x(expr,14);case 0x##x##14:fun##x(expr,13);case 0x##x##13:fun##x(expr,12);case 0x##x##12:fun##x(expr,11);case 0x##x##11:fun##x(expr,10);case 0x##x##10:fun##x(expr,9);case 0x##x##9:fun##x(expr,8);case 0x##x##8:fun##x(expr,7);case 0x##x##7:fun##x(expr,6);case 0x##x##6:fun##x(expr,5);case 0x##x##5:fun##x(expr,4);case 0x##x##4:fun##x(expr,3);case 0x##x##3:fun##x(expr,2);case 0x##x##2:fun##x(expr,1);case 0x##x##1:fun##x(expr,0);break;

#define REP_J1(fun1,fun2,expr)  
#define REP_J2(fun1,fun2,expr)  fun1(fun2,expr,1);
#define REP_J3(fun1,fun2,expr)  fun1(fun2,expr,2);fun1(fun2,expr,1);
#define REP_J4(fun1,fun2,expr)  fun1(fun2,expr,3);fun1(fun2,expr,2);fun1(fun2,expr,1);
#define REP_J5(fun1,fun2,expr)  fun1(fun2,expr,4);fun1(fun2,expr,3);fun1(fun2,expr,2);fun1(fun2,expr,1);
#define REP_J6(fun1,fun2,expr)  fun1(fun2,expr,5);fun1(fun2,expr,4);fun1(fun2,expr,3);fun1(fun2,expr,2);fun1(fun2,expr,1);
#define REP_J7(fun1,fun2,expr)  fun1(fun2,expr,6);fun1(fun2,expr,5);fun1(fun2,expr,4);fun1(fun2,expr,3);fun1(fun2,expr,2);fun1(fun2,expr,1);
#define REP_J8(fun1,fun2,expr)  fun1(fun2,expr,7);fun1(fun2,expr,6);fun1(fun2,expr,5);fun1(fun2,expr,4);fun1(fun2,expr,3);fun1(fun2,expr,2);fun1(fun2,expr,1);
#define REP_J9(fun1,fun2,expr)  fun1(fun2,expr,8);fun1(fun2,expr,7);fun1(fun2,expr,6);fun1(fun2,expr,5);fun1(fun2,expr,4);fun1(fun2,expr,3);fun1(fun2,expr,2);fun1(fun2,expr,1);
#define REP_J10(fun1,fun2,expr)  fun1(fun2,expr,9);fun1(fun2,expr,8);fun1(fun2,expr,7);fun1(fun2,expr,6);fun1(fun2,expr,5);fun1(fun2,expr,4);fun1(fun2,expr,3);fun1(fun2,expr,2);fun1(fun2,expr,1);
#define REP_J11(fun1,fun2,expr)  fun1(fun2,expr,10);fun1(fun2,expr,9);fun1(fun2,expr,8);fun1(fun2,expr,7);fun1(fun2,expr,6);fun1(fun2,expr,5);fun1(fun2,expr,4);fun1(fun2,expr,3);fun1(fun2,expr,2);fun1(fun2,expr,1);
#define REP_J12(fun1,fun2,expr)  fun1(fun2,expr,11);fun1(fun2,expr,10);fun1(fun2,expr,9);fun1(fun2,expr,8);fun1(fun2,expr,7);fun1(fun2,expr,6);fun1(fun2,expr,5);fun1(fun2,expr,4);fun1(fun2,expr,3);fun1(fun2,expr,2);fun1(fun2,expr,1);
#define REP_J13(fun1,fun2,expr)  fun1(fun2,expr,12);fun1(fun2,expr,11);fun1(fun2,expr,10);fun1(fun2,expr,9);fun1(fun2,expr,8);fun1(fun2,expr,7);fun1(fun2,expr,6);fun1(fun2,expr,5);fun1(fun2,expr,4);fun1(fun2,expr,3);fun1(fun2,expr,2);fun1(fun2,expr,1);
#define REP_J14(fun1,fun2,expr)  fun1(fun2,expr,13);fun1(fun2,expr,12);fun1(fun2,expr,11);fun1(fun2,expr,10);fun1(fun2,expr,9);fun1(fun2,expr,8);fun1(fun2,expr,7);fun1(fun2,expr,6);fun1(fun2,expr,5);fun1(fun2,expr,4);fun1(fun2,expr,3);fun1(fun2,expr,2);fun1(fun2,expr,1);
#define REP_J15(fun1,fun2,expr)  fun1(fun2,expr,14);fun1(fun2,expr,13);fun1(fun2,expr,12);fun1(fun2,expr,11);fun1(fun2,expr,10);fun1(fun2,expr,9);fun1(fun2,expr,8);fun1(fun2,expr,7);fun1(fun2,expr,6);fun1(fun2,expr,5);fun1(fun2,expr,4);fun1(fun2,expr,3);fun1(fun2,expr,2);fun1(fun2,expr,1);
#define REP_J16(fun1,fun2,expr)  fun1(fun2,expr,15);fun1(fun2,expr,14);fun1(fun2,expr,13);fun1(fun2,expr,12);fun1(fun2,expr,11);fun1(fun2,expr,10);fun1(fun2,expr,9);fun1(fun2,expr,8);fun1(fun2,expr,7);fun1(fun2,expr,6);fun1(fun2,expr,5);fun1(fun2,expr,4);fun1(fun2,expr,3);fun1(fun2,expr,2);fun1(fun2,expr,1);

#define LOOP_helper2(n,expr)		REP_A##n(expr)
#define LOOP2_helper2(n,expr)		REP_D##n(REP_B##n,expr)
#define LOOP3_helper2(n,expr)		REP_E##n(REP_B,expr)
#define LOOPNa_helper2(n,expr)		switch (Na) { REP_F##n(expr); }
#define LOOPNb_helper2(n,expr)		switch (Nb) { REP_F##n(expr); }
#define LOOPNk_helper2(n,expr)      switch (Nk) { REP_F##n(expr); }
#define LOOP2Na_helper2(n,expr)		switch (Na) { REP_G##n(REP_B##n,expr); }
#define LOOP2Nb_helper2(n,expr)		switch (Nb) { REP_G##n(REP_C##n,expr); }
#define LOOP3Na_helper2(n,expr)		switch (Na) { REP_H##n(REP_B,expr); }
#define LOOP2NaNb_helper2(n,expr)	switch((Nb << 4) | Na)  { REP_J##n(REP_I##n,REP_B,expr); }

#define LOOP_helper1(n,expr)		LOOP_helper2(n,expr)
#define LOOP2_helper1(n,expr)		LOOP2_helper2(n,expr)
#define LOOP3_helper1(n,expr)		LOOP3_helper2(n,expr)
#define LOOPNa_helper1(n,expr)		LOOPNa_helper2(n,expr)
#define LOOPNb_helper1(n,expr)		LOOPNb_helper2(n,expr)
#define LOOPNk_helper1(n,expr)      LOOPNk_helper2(n,expr)
#define LOOP2Na_helper1(n,expr)		LOOP2Na_helper2(n,expr)
#define LOOP2Nb_helper1(n,expr)		LOOP2Nb_helper2(n,expr)
#define LOOP3Na_helper1(n,expr)		LOOP3Na_helper2(n,expr)
#define LOOP2NaNb_helper1(n,expr)	LOOP2NaNb_helper2(n,expr)

#define LOOP(expr)					LOOP_helper1(N,expr)
#define LOOP2(expr)					LOOP2_helper1(N,expr)
#define LOOP3(expr)					LOOP3_helper1(N,expr)
#define LOOPNa(expr)				LOOPNa_helper1(N,expr)
#define LOOPNb(expr)				LOOPNb_helper1(N,expr)
#define LOOPNk(expr)                LOOPNk_helper1(N,expr)
#define LOOP2Na(expr)				LOOP2Na_helper1(N,expr)
#define LOOP2Nb(expr)				LOOP2Nb_helper1(N,expr)
#define LOOP3Na(expr)				LOOP3Na_helper1(N,expr)
#define LOOP2NaNb(expr)				LOOP2NaNb_helper1(N,expr)

#define REP(x)				for (int kk = 0; kk < (x); ++kk)
#define DEL(x)				{ if (x) delete[] (x); (x) = NULL; }
#define FREE(x)				{ if (x) free(x); (x) = NULL; }
#define FREECUDA(x)			{ if (x) FreeHostCUDA(x); (x) = NULL; }

#ifdef HASH64
	typedef uint64 HASH;
#else
	typedef uint HASH;
#endif

#define THREADH(x) \
		TARGET void x(int id);\
		TARGET void x##In(void)

#define THREAD(x) \
		TARGET void x(int id)\
		{\
			threadid = (uint)(uint64)id;\
			x##In();\
		}\
		TARGET void x##In(void)

#define THREAD2H(x) \
		template<typename REAL> TARGET void x(int id);\
		template<typename REAL> TARGET void x##In(void)

#define THREAD2(x) \
		template TARGET void x<double>(int id); \
		template TARGET void x<float >(int id); \
		template TARGET void x##In<double>(void); \
		template TARGET void x##In<float >(void); \
		template<typename REAL> \
		TARGET void x(int id)\
		{\
			threadid = (uint)id;\
			x##In<REAL>();\
		}\
		template<typename REAL> \
		TARGET void x##In(void)

/* Lock-free thread competition */

//state_lock[ii % NBUF] == ii * 4 + (0 for available, 1 for thread processing, 3 for thread processed)

#define THREAD_BEGIN    while (ii >= progress1 + NBUF) Sleep(SLEEP_TIME_TINY); \
						int64 avail_val = ii * 4; \
						if (state_lock[ii % NBUF].compare_exchange_strong(avail_val, avail_val + 1)){

#define THREAD_END      state_lock[ii % NBUF] = ii * 4 + 3;}

#define GUARD_BEGIN     while (state_lock[ii % NBUF] != ii * 4 + 3) Sleep(SLEEP_TIME_TINY);

#define GUARD_END		state_lock[ii % NBUF] = (ii + NBUF) * 4;

#define THREAD_BEGIN2   while (ii >= progress1 + NBUF || ii >= progress2 + NBUF) Sleep(SLEEP_TIME_TINY); \
						int64 avail_val = ii * 4; \
						if (state_lock[ii % NBUF].compare_exchange_strong(avail_val, avail_val + 1)){

#define THREAD_END2     state_lock[ii % NBUF] = ii * 4 + 2;}

#define GUARD_BEGIN2    while (state_lock[ii % NBUF] >> 1 != ii * 2 + 1) Sleep(SLEEP_TIME_TINY);

#define GUARD_END2		int64 test_val = ii * 4 + 4, next_val = (ii + NBUF) * 4; \
						state_lock[ii % NBUF]++; \
						state_lock[ii % NBUF].compare_exchange_strong(test_val, next_val);

#define sfprintf(f,...) {sprintf(f,__VA_ARGS__); f += strlen(f);}

/* Class and Struct names*/
struct BCFHEADER;
struct VCF;
struct MEMORY;
struct GENOTYPE;
struct SLOCUS;
struct LOCUS;
template<typename REAL> struct POP;
template<typename REAL> struct IND;
struct STRUCTURE_RUNINFO;
template <typename T> struct LIST;
template <typename T, typename T2> struct TABLE_ENTRY;
template <typename T, typename T2> struct TABLE;
struct Huang2015ENTRY;
template<typename REAL> struct RNGSSE;
template<typename REAL> struct RNGAVX;
template<typename REAL> struct RNG512;
template<typename REAL> struct RNGNEO;
struct HAPLO_DUMMY_LOCUS;
struct FILEINFO;

#include "dre.h"
#include "ml.h"
#include "ml2.h"
#include "mlbin.h"
#include "mom.h"
#include "mom2.h"
#include "global.h"
#include "hash.h"
#include "parameters.h"
#include "math2.h"
#include "mathNEO.h"
#include "mathSSE.h"
#include "mathAVX.h"
#include "math512.h"
#include "misc.h"
#include "file.h"
#include "string2.h"
#include "statistics.h"
#include "matrix.h"

#include "spa.h"
#include "ploidyinfer.h"
#include "function.h"
#include "load.h"
#include "filter.h"
#include "diversity.h"
#include "slide.h"
#include "decay.h"
#include "block.h"
#include "gwas.h"
#include "haplotype.h"
#include "conversion.h"
#include "indstat.h"
#include "dist.h"
#include "pcoa.h"
#include "clustering.h"
#include "diff.h"
#include "kinship.h"
#include "relatedness.h"
#include "amova.h"
#include "popas.h"
#include "structure.h"
#include "structureNEO.h"
#include "structureSSE.h"
#include "structureAVX.h"
#include "structure512.h"
#include "structureCUDA.h"
#include "ad.h"
#include "menu.h"

TARGET int main(int _argc, char** _argv); 
