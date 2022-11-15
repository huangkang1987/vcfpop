#pragma once

#define WINDOWS_IGNORE_PACKING_MISMATCH
#define _CRT_SECURE_NO_WARNINGS

#ifdef testarm64
	//test code
	#define __aarch64__
	//#define __ARM_FP
	#define __ARM_NEON
	#define __ARM_ACLE
	#define _M_HYBRID_X86_ARM64
	//#define __ARM_FEATURE_BF16
	//#define __LITTLE_ENDIAN__
	#undef _WIN64
#endif

#ifdef __aarch64__
	#include <arm_neon.h>
	#include <arm_acle.h>
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
#include <vector>
#include <algorithm>
#include <utility>
#include <atomic> 
#include <shared_mutex>
#include <mutex> 
#include <thread>
#include <omp.h>

#ifndef __CUDA__
	using namespace std::filesystem;
#endif


using std::pair;
using std::mutex;
using std::shared_mutex;
using std::map;
using std::vector;
using std::string;
using std::thread;
using std::atomic;
using std::atomic_flag;

//#define EIGEN_DONT_ALIGN_STATICALLY
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

	#define TARGET    
	#define TARGETSSE 
	#define TARGETAVX 
	#define TARGET512 
	#define TARGETNEO 

	#define fseeko64 _fseeki64
	#define ftello64 _ftelli64 

	#include <windows.h>
	#include <psapi.h>
	#pragma comment (lib,"psapi.lib")
	#include "zlib/zlib.h"  
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
	#if __aarch64__
		#define TARGET
		#define TARGETSSE
		#define TARGETAVX
		#define TARGET512
		#define TARGETNEO  __attribute__((__target__("neon")))
	#else
		#define TARGET     
		#define TARGETSSE 
		#define TARGETAVX  __attribute__((__target__("avx"))) __attribute__((__target__("avx2")))
		#define TARGET512  __attribute__((__target__("avx"))) __attribute__((__target__("avx2"))) __attribute__((__target__("avx512f"))) __attribute__((__target__("avx512bw"))) __attribute__((__target__("avx512dq"))) __attribute__((__target__("avx512vl")))
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


#if defined(__clang__) || defined(__GNUC__)
	typedef float		__m128  __attribute__((__vector_size__(16), __aligned__(16)));
	typedef double		__m128d __attribute__((__vector_size__(16), __aligned__(16)));
	typedef long long	__m128i __attribute__((__vector_size__(16), __aligned__(16)));

	typedef float		__m256  __attribute__((__vector_size__(32), __aligned__(32)));
	typedef double		__m256d __attribute__((__vector_size__(32), __aligned__(32)));
	typedef long long	__m256i __attribute__((__vector_size__(32), __aligned__(32)));

	typedef float		__m512  __attribute__((__vector_size__(64), __aligned__(64)));
	typedef double		__m512d __attribute__((__vector_size__(64), __aligned__(64)));
	typedef long long	__m512i __attribute__((__vector_size__(64), __aligned__(64)));

	#define VLA_NEW(name,type,size) type name##VLA[size]; type* name = (type*)name##VLA
	#define VLA_DELETE(name)
	#define __forceinline __inline__ __attribute__((always_inline))
#else
	#define VLA_NEW(name,type,size) type* name = new type[size]
	#define VLA_DELETE(name) delete[] name
#endif

#define simd_f64(x,y)	((double*)&(x))[y]
#define simd_f32(x,y)	((float*)&(x))[y]
#define simd_u256(x,y)	((__m256i*)&(x))[y]
#define simd_u128(x,y)	((__m128i*)&(x))[y]
#define simd_u64(x,y)	((uint64*)&(x))[y]
#define simd_u32(x,y)	((uint*)&(x))[y]
#define simd_u16(x,y)	((ushort*)&(x))[y]
#define simd_u8(x,y)	((byte*)&(x))[y]
#define simd_i256(x,y)	((__m256i*)&(x))[y]
#define simd_i128(x,y)	((__m128i*)&(x))[y]
#define simd_i64(x,y)	((int64*)&(x))[y]
#define simd_i32(x,y)	((int*)&(x))[y]
#define simd_i16(x,y)	((short*)&(x))[y]
#define simd_i8(x,y)	((char*)&(x))[y]

#define simp_f64(x,y)	((double*)(x))[y]
#define simp_f32(x,y)	((float*)(x))[y]
#define simp_u256(x,y)	((__m256i*)(x))[y]
#define simp_u128(x,y)	((__m128i*)(x))[y]
#define simp_u64(x,y)	((uint64*)(x))[y]
#define simp_u32(x,y)	((uint*)(x))[y]
#define simp_u16(x,y)	((ushort*)(x))[y]
#define simp_u8(x,y)	((byte*)(x))[y]
#define simp_i256(x,y)	((__m256i*)(x))[y]
#define simp_i128(x,y)	((__m128i*)(x))[y]
#define simp_i64(x,y)	((int64*)(x))[y]
#define simp_i32(x,y)	((int*)(x))[y]
#define simp_i16(x,y)	((short*)(x))[y]
#define simp_i8(x,y)	((char*)(x))[y]

#define GetLoc(x)		(slocus[(x)])				//(useslocus ? slocus[l] : locus[l])
#define GetLocPos(x)	(locus_pos[(x)])			//(useslocus ? locus_pos[l] : locus[l].pos)
#define GetLocId(x)		(locus_id[(x)])				//(useslocus ? locus_id[l] : locus[l].id)
#define GetLocK(x)			((x) < lend ? slocus[(x)].k : 0)
#define GetLocTabDiff(x)	((x) < lend ? (int64)slocus[(x)].GetGtab() - gtab_base : 0)	
#define REP(x)				for (int kk = 0; kk < (x); ++kk)

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
#include "math2.h"
#include "mathNEO.h"
#include "mathSSE.h"
#include "mathAVX.h"
#include "math512.h"
#include "parameters.h"
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
