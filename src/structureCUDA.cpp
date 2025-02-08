/* CUDA Bayesian clustering functions */
#define __CUDA__

#include "vcfpop.h"

template TARGET void BAYESIAN<double>::UpdateQMetroCUDA();
template TARGET void BAYESIAN<float >::UpdateQMetroCUDA();

template TARGET void BAYESIAN<double>::UpdateQNoAdmixCUDA();
template TARGET void BAYESIAN<float >::UpdateQNoAdmixCUDA();

template TARGET void BAYESIAN<double>::UpdateZAdmixCUDA();
template TARGET void BAYESIAN<float >::UpdateZAdmixCUDA();

template TARGET void BAYESIAN<double>::UpdateZNoAdmixCUDA();
template TARGET void BAYESIAN<float >::UpdateZNoAdmixCUDA();

template TARGET void BAYESIAN<double>::RecordCUDA();
template TARGET void BAYESIAN<float >::RecordCUDA();

TARGET void* Malloc(uint64 size)
{
    void* re = NULL;
    while (re == NULL)
    {
        re = malloc(size);
        if (re) break;
        printf("Fail to allocate %0.3f Gib memory, increase virtual memory and retry (this will not terminate this process).", size / 1024.0 / 1024.0 / 1024.0);
        Pause();
    } 
    return re;
}

#ifdef CUDA

#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <cublas_v2.h>

//cublas.so
decltype(&cublasCreate_v2) cublasCreateA;
decltype(&cublasDestroy_v2) cublasDestroyA;
decltype(&cublasSgemm_v2) cublasSgemmA;
decltype(&cublasDgemm_v2) cublasDgemmA;

//cusolver.so
decltype(&cusolverDnCreate) cusolverDnCreateA;
decltype(&cusolverDnDestroy) cusolverDnDestroyA;
decltype(&cusolverDnSetStream) cusolverDnSetStreamA;
decltype(&cusolverDnDgesvd) cusolverDnDgesvdA;
decltype(&cusolverDnSgesvd) cusolverDnSgesvdA;
decltype(&cusolverDnDsyevd) cusolverDnDsyevdA;
decltype(&cusolverDnSsyevd) cusolverDnSsyevdA;
decltype(&cusolverDnDgesvd_bufferSize) cusolverDnDgesvd_bufferSizeA;
decltype(&cusolverDnSgesvd_bufferSize) cusolverDnSgesvd_bufferSizeA;
decltype(&cusolverDnDsyevd_bufferSize) cusolverDnDsyevd_bufferSizeA;
decltype(&cusolverDnSsyevd_bufferSize) cusolverDnSsyevd_bufferSizeA;

TARGETCUDA void InitCUDA()
{
    string cudaPath = std::getenv("CUDA_PATH");
#ifdef _WIN64
    #define LoadLib(x) LoadLibraryA(x)
    #define LoadFunc(x,y) GetProcAddress(x,y)
    #define LibExt ".dll"
    HMODULE hblas = NULL, hsolver = NULL;
#else
    #define LoadLib(x) dlopen(x, RTLD_LAZY)
    #define LoadFunc(x,y) dlsym(x,y)
    #define LibExt ".so"
    void* hblas = NULL, *hsolver = NULL;
#endif
    
    for (int i = 19; i >= 10; --i)
    {
        string cublas = cudaPath + PATH_DELIM + "bin" + PATH_DELIM + "cublas64_" + std::to_string(i) + LibExt;
        string cusolver = cudaPath + PATH_DELIM + "bin" + PATH_DELIM + "cusolver64_" + std::to_string(i) + LibExt;
        if (!hblas && FileExists(cublas.c_str()))
            hblas = LoadLib(cublas.c_str());
        if (!hsolver && FileExists(cusolver.c_str()))
            hsolver = LoadLib(cusolver.c_str());
    }
	
    if (!hblas) Exit("\nError: cublas.dll is not found!\n");
    if (!hsolver) Exit("\nError: cusolver.dll is not found!\n");

    cublasCreateA  = (decltype(cublasCreateA)) LoadFunc(hblas, "cublasCreate_v2");
    cublasDestroyA = (decltype(cublasDestroyA))LoadFunc(hblas, "cublasDestroy_v2");
    cublasSgemmA   = (decltype(cublasSgemmA))  LoadFunc(hblas, "cublasSgemm_v2");
    cublasDgemmA   = (decltype(cublasDgemmA))  LoadFunc(hblas, "cublasDgemm_v2");
    
    cusolverDnCreateA            = (decltype(cusolverDnCreateA))           LoadFunc(hsolver, "cusolverDnCreate");
    cusolverDnDestroyA           = (decltype(cusolverDnDestroyA))          LoadFunc(hsolver, "cusolverDnDestroy");
    cusolverDnSetStreamA         = (decltype(cusolverDnSetStreamA))        LoadFunc(hsolver, "cusolverDnSetStream");
    cusolverDnDgesvdA            = (decltype(cusolverDnDgesvdA))           LoadFunc(hsolver, "cusolverDnDgesvd");
    cusolverDnSgesvdA            = (decltype(cusolverDnSgesvdA))           LoadFunc(hsolver, "cusolverDnSgesvd");
    cusolverDnDsyevdA            = (decltype(cusolverDnDsyevdA))           LoadFunc(hsolver, "cusolverDnDsyevd");
    cusolverDnSsyevdA            = (decltype(cusolverDnSsyevdA))           LoadFunc(hsolver, "cusolverDnSsyevdA");
    cusolverDnDgesvd_bufferSizeA = (decltype(cusolverDnDgesvd_bufferSizeA))LoadFunc(hsolver, "cusolverDnDgesvd_bufferSize");
    cusolverDnSgesvd_bufferSizeA = (decltype(cusolverDnSgesvd_bufferSizeA))LoadFunc(hsolver, "cusolverDnSgesvd_bufferSize");
    cusolverDnDsyevd_bufferSizeA = (decltype(cusolverDnDsyevd_bufferSizeA))LoadFunc(hsolver, "cusolverDnDsyevd_bufferSize");
    cusolverDnSsyevd_bufferSizeA = (decltype(cusolverDnSsyevd_bufferSizeA))LoadFunc(hsolver, "cusolverDnSsyevd_bufferSize");
}

#define BLOCK_DIM (64)
#define checkErrors(val) Check((val), #val, __FILE__, __LINE__)

// CUDA API error checking
#define CUDA_CHECK(err)                                             \
  do {                                                              \
    cudaError_t err_ = (err);                                       \
    if (err_ != cudaSuccess) {                                      \
      printf("CUDA error %d at %s:%d\n", err_, __FILE__, __LINE__); \
      throw std::runtime_error("CUDA error");                       \
    }                                                               \
  } while (0)

// cusolver API error checking
#define CUSOLVER_CHECK(err)                                             \
  do {                                                                  \
    cusolverStatus_t err_ = (err);                                      \
    if (err_ != CUSOLVER_STATUS_SUCCESS) {                              \
      printf("cusolver error %d at %s:%d\n", err_, __FILE__, __LINE__); \
      throw std::runtime_error("cusolver error");                       \
    }                                                                   \
  } while (0)

// cublas API error checking
#define CUBLAS_CHECK(err)                                                                          \
    do {                                                                                           \
        cublasStatus_t err_ = (err);                                                               \
        if (err_ != CUBLAS_STATUS_SUCCESS) {                                                       \
            std::printf("cublas error %d at %s:%d\n", err_, __FILE__, __LINE__);                   \
            throw std::runtime_error("cublas error");                                              \
        }                                                                                          \
    } while (0)

template <typename T>
void Check(T result, char const* const func, const char* const file, int const line)
{
    if (result)
    {
        printf("CUDA error at %s:%d code=%d(%s) \"%s\" \n",
            file,
            line,
            static_cast<unsigned int>(result),
            cudaGetErrorName(result),
            func);
        exit(EXIT_FAILURE); 
    }
}

/* Get square root of val, account some slightly negative input */
template<typename REAL>
__device__ REAL MySqrtDevice(REAL val)
{
    if (val > 0)
        return sqrt(val);
    else if (val > -MIN_FREQ)
        return 0;
    else
        return NAN;
}

/* https://github.com/ygalanter/CyberGeeks/blob/master/src/math.c */
__device__ double MyRIntDevice(double x)
{
    double t = floor(fabs(x) + 0.5);
    return (x < 0.0) ? -t : t;
}

/* Core function of cosine */
__device__ double CosCoreDevice(double x)
{
    double x2 = x * x;
    double x4 = x2 * x2;
    double x8 = x4 * x4;
    return (-2.7236370439787708e-7 * x2 + 2.4799852696610628e-5) * x8 +
        (-1.3888885054799695e-3 * x2 + 4.1666666636943683e-2) * x4 +
        (-4.9999999999963024e-1 * x2 + 1.0000000000000000e+0);
}

/* Core function of sine */
__device__ double SinCoreDevice(double x)
{
    double x2 = x * x;
    double x4 = x2 * x2;
    return ((2.7181216275479732e-6 * x2 - 1.9839312269456257e-4) * x4 +
        (8.3333293048425631e-3 * x2 - 1.6666666640797048e-1)) * x2 * x + x;
}

/* Core function of arc sine */
__device__ double ArcSinCoreDevice(double x)
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
__device__ double MySinDevice(double x)
{
    double q = MyRIntDevice(x * 6.3661977236758138e-1);
    int quadrant = (int)q;
    double t = x - q * 1.5707963267923333e+00;
    t = t - q * 2.5633441515945189e-12;
    t = quadrant & 1 ? CosCoreDevice(t) : SinCoreDevice(t);
    return (quadrant & 2) ? -t : t;
}

/* Cosine function */
__device__ double MyCosDevice(double x)
{
    return MySinDevice(x + (M_PI / 2));
}

/* ArcCosine function */
__device__ double MyArcCosDevice(double x)
{
    double xa = abs(x);
    double t = xa > 0.5625 ?
        2.0 * ArcSinCoreDevice(sqrt(0.5 * (1.0 - xa))) :
        1.5707963267948966 - ArcSinCoreDevice(xa);
    return (x < 0.0) ? (3.1415926535897932 - t) : t;
}

/* ArcSine function */
__device__ double MyArcSinDevice(double x)
{
    return (M_PI / 2) - MyArcCosDevice(x);
}

/* Tangent function */
__device__ double MyTanDevice(double x)
{
    return MySinDevice(x) / MyCosDevice(x);
}

/*=============================== DATA STRUCTURE ===============================*/
#pragma pack(push, 1)
struct GENO_READER_CUDA
{
    uint64* pos;							//Current read pointer
    uint64 data;							//Readed bits
    int size;								//Number of bits a genotype id used
    int nbits;								//Number of bits remaining in data

    /* Initialize reader */
    __device__ GENO_READER_CUDA(uint64 offset);

    /* Get id of next ind (order by indid) */
    __device__ int Read();
};

struct GENOTYPE_CUDA
{
    uint v1;
    /*
    //4 bytes
    uint offset : 24;						//this pointer + offset * 2 is address of alleles array
    //ushort* alleles;				        //Alleles copy in ascending order [ploidy]
                                            //Unique alleles [nalleles] order by dosage descending
    uint patternid : 8;						//Pattern index
    */

    __device__ int Nalleles();

    __device__ int Ploidy();

    __device__ ushort* GetAlleleArray();
};

struct SLOCUS_CUDA
{
    //small locus, 12 bytes
    uint64 v1;
    uint   v2;

    /*
    uint64 bits1 : 48;						//Genotype table[ngeno]
                                            //alen[k] {for non-vcf SMM distance}
                                            //chrom \0 name \0 {(allele identifiers \0)[k] for vcf/bcf)}
                                            //genotype alleles[gasize]
    uint64 k : 16;							//Number of alleles

    uint   ngeno : 22;						//Number of genotypes
    uint   flag_pass : 1;					//0 pass filter
    uint   flag_alen : 1;					//Has alen table?
    uint   pes_model : 8;					//1 for RCS, 2 for PRCS, 3 for CES, 4+ for PES
    */

    __device__ GENOTYPE_CUDA* GetGtab();

    __device__ int GetK();

    __device__ int GetNgeno();

};

struct MEMORY_CUDA
{
    int KT;
    int N;
    int maxploidy;
    int minploidy;
    int pad;
    int64 L;

    MEMORY_CUDA* device_mem;        //GPU memory block
    byte* g_bucket;                 //address of genotype bucket
    uint64* loc_addr;               //address of genotype data in genotype_bucket
    int64* allele_freq_offset;                //the index of 1st allele among all alleles at all locus
    SLOCUS_CUDA* sloc;              //slocus in GPU memory
};

template<typename REAL>
struct RNG_CUDA
{

};

template<>
struct RNG_CUDA<double>
{
    uint64 x;
    uint64 y;

    /* Initialize rng */
    __device__ RNG_CUDA();

    /* Initialize rng */
    __device__ RNG_CUDA(uint64 seed, uint64 salt);

    /* Draw a uniform distriubted interger */
    __device__ uint64 XorShift();

    /* Draw a polynormial distriubted integer */
    __device__ int Poly(double* a, int n, bool typed);

    /*
    
    double U1;
    double U2;

    // Draw a vector from Dirichlet distribution D(a1 + b1, a2 + b2, ...)
    template<typename T1, typename T2, typename T3>
    __device__ void Dirichlet(T1* res, T2* a, T3* b, int n);

    // Draw a real number from gamma distribution
    __device__ double Gamma(double alpha, double beta = 1);

    // Draw a normal distriubted real number
    __device__ double Normal();

    // Draw a uniform distriubted real number
    __device__ double Uniform();

    */
};

template<>
struct RNG_CUDA<float >
{
    uint x;
    uint y;
    uint z;

    /* Initialize rng */
    __device__ RNG_CUDA();

    /* Initialize rng */
    __device__ RNG_CUDA(uint64 seed, uint64 salt);

    /* Draw a uniform distriubted interger */
    __device__ uint XorShift();

    /* Draw a polynormial distriubted integer */
    __device__ int Poly(float* a, int n, bool typed);

    /*
    
    double U1;
    double U2;

    //Draw a vector from Dirichlet distribution D(a1 + b1, a2 + b2, ...)
    template<typename T1, typename T2, typename T3>
    __device__ void Dirichlet(T1* res, T2* a, T3* b, int n);

    //Draw a real number from gamma distribution
    __device__ double Gamma(double alpha, double beta = 1);

    //Draw a normal distriubted real number
    __device__ double Normal();

    //Draw a uniform distriubted real number
    __device__ float Uniform(); 

    */
};

#pragma pack(pop)

/*=============================== GLOVAL VARIABLES ===============================*/

MEMORY_CUDA* cuda_mem;

thread_local cudaStream_t GPUstream;
thread_local cudaDeviceProp GPUprop;
__device__ cudaDeviceProp GPUprop_device;
__device__ MEMORY_CUDA cuda_mem_device;
__device__ ushort missing_array_CUDA[N_MAX_PLOIDY] = { 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF };
__device__ int PT_NALLELES[150] = 									//Pattern index to number of alleles
{ 0, 1, 1, 2, 1, 2, 3, 1, 2, 2, 3, 4, 1, 2, 2, 3, 3, 4, 5, 1, 2, 2, 3, 2, 3, 4, 3, 4, 5, 6, 1, 2, 2, 3, 2, 3, 4, 3, 3, 4, 5, 4, 5, 6, 7, 1, 2, 2, 3, 2, 3, 4, 2, 3, 3, 4, 5, 3, 4, 4, 5, 6, 4, 5, 6, 7, 8, 1, 2, 2, 3, 2, 3, 4, 2, 3, 3, 4, 5, 3, 3, 4, 4, 5, 6, 3, 4, 5, 4, 5, 6, 7, 5, 6, 7, 8, 9, 1, 2, 2, 3, 2, 3, 4, 2, 3, 3, 4, 5, 2, 3, 3, 4, 4, 5, 6, 3, 4, 3, 4, 5, 4, 5, 6, 7, 4, 4, 5, 6, 5, 6, 7, 8, 5, 6, 7, 8, 9, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
__device__ int PT_PLOIDY[150] = 									//Pattern index to ploidy level
{ 0, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };

/*=============================== MEMBER FUNCTIONS ===============================*/

#ifndef _GENOTYPE_CUDA
    __device__ int GENOTYPE_CUDA::Nalleles()
    {
        return PT_NALLELES[v1 >> 24];
    }

    __device__ int GENOTYPE_CUDA::Ploidy()
    {
        return PT_PLOIDY[v1 >> 24];
    }

    __device__ ushort* GENOTYPE_CUDA::GetAlleleArray()
    {
        uint offset = v1 & 0xFFFFFF;
        switch (offset)
        {
        case 0:
            return NULL;
        case 0xFFFFFF:
            return missing_array_CUDA;
        default:
            return (ushort*)((byte*)this + offset);
        }
    }

    __device__ GENOTYPE_CUDA* SLOCUS_CUDA::GetGtab()
    {
        return (GENOTYPE_CUDA*)(v1 & 0xFFFFFFFFFFFF);
    }
#endif

#ifndef _SLOCUS_CUDA
    __device__ int SLOCUS_CUDA::GetK()
    {
        return v1 >> 48;
    }

    __device__ int SLOCUS_CUDA::GetNgeno()
    {
        return v1 & 0x3FFFFF;
    }
#endif

#ifndef _GENO_READER_CUDA
    __device__ GENO_READER_CUDA::GENO_READER_CUDA(uint64 offset)
    {
        int shift = (offset & 7) << 3;
        pos = (uint64*)(offset & 0xFFFFFFFFFFF8);
        size = offset >> 48;
        data = (*pos++) >> shift;
        nbits = 64 - shift;
    }

    __device__ int GENO_READER_CUDA::Read()
    {
        // if data is empty
        if (nbits < size)
        {
            // move nbits data to gid
            int gid = (int)data;

            // remain number of bits to read
            int rbits = size - nbits;

            // read 64 bits to data
            data = *pos++;

            // read rbits from data and concate to higher bits in gid
            gid |= ((uint)data & ((1u << rbits) - 1u)) << nbits;

            //shift right
            data >>= rbits;
            nbits = 64 - rbits;
            return gid;
        }
        else
        {
            //read size bits
            int gid = (int)data & ((1u << size) - 1u);

            //shift right
            data >>= size;
            nbits -= size;
            return gid;
        }
    }
#endif

#ifndef _HASH_CUDA

    /* 32 bit Integer Hashing */
    __device__ uint MurmurHash2_CUDA(uint data, uint seed)
    {
        uint m = 0x5bd1e995;
        uint s = seed ^ sizeof(uint);

        uint a = data;

        a *= m;
        a ^= a >> 24;
        a *= m;
        s *= m;
        a ^= s;

        a ^= a >> 13;
        a *= m;
        a ^= a >> 15;

        return a;
    }

    /* 64 bit Integer Hashing */
    __device__ uint64 MurmurHash64_CUDA(uint64 data, uint64 seed)
    {
        uint* d = (uint*)&data;
        uint* s = (uint*)&seed;

        d[1] = (~d[0]) ^ d[1];
        s[1] = (~s[0]) ^ s[1];

        return ((uint64)MurmurHash2_CUDA(d[1], s[1]) << 32) |
               ((uint64)MurmurHash2_CUDA(d[0], s[0]));
    }

    /* 64 bit Integer Hashing */
    __device__ uint MurmurHash32_CUDA(uint64 data, uint64 seed)
    {
        uint a32 = (~(uint)data) ^ ((uint)(data >> 32));
        uint s32 = (~(uint)seed) ^ ((uint)(seed >> 32));

        return MurmurHash2_CUDA(a32, s32);
    }

    /* Mix high and low 32 bits */
    __device__ uint Mix_CUDA(uint64 x)
    {
        return (~(uint)x) ^ ((uint)(x >> 32));
    }
#endif

 #ifndef _RNG_CUDA_FP64
    __device__ RNG_CUDA<double>::RNG_CUDA()
    {

    }

    /* Initialize rng */
    __device__ RNG_CUDA<double>::RNG_CUDA(uint64 s, uint64 salt)
    {
        s = MurmurHash64_CUDA(s, salt);

        x = 0x159A55E5075BCD15 ^ (s);      //123456789, 362436069
        y = 0x054913331F123BB5 ^ (s << 6); //521288629, 88675123
    }

    /* Draw a uniform distriubted interger */
    __device__ uint64 RNG_CUDA<double>::XorShift()
    {
        //XorShift
        uint64 a = x, b = y;

        x = b;
        a ^= a << 23;
        a ^= a >> 18;
        a ^= b;
        a ^= b >> 5;
        y = a;

        return a + b;
    }

    /* Draw a polynormial distriubted integer */
    __device__ int RNG_CUDA<double>::Poly(double* a, int n, bool typed)
    {
        //row unify
        double s = 0;

        if (typed)
        {
            volatile double v1 = (double)MIN_FREQ;
            for (int i = 0; i < n; ++i)
            {
                a[i] += v1;
                s += a[i];
            }
        }

        uint64 u = XorShift();

        if (typed)
        {
            uint64 r = 0x3FF0000000000000;
            double& re = *(double*)&r;
            r |= u & 0x000FFFFFFFFFFFFF;

            double t = (re - 1.0) * s;

            for (int i = 0; i < n; ++i)
            {
                if (t < a[i]) return i;
                t -= a[i];
            }
        }
        return n - 1;
    }

    /*
    // Draw a vector from Dirichlet distribution D(a1 + b1, a2 + b2, ...)
    template<typename T1, typename T2, typename T3>
    __device__ void RNG_CUDA<double>::Dirichlet(T1* res, T2* a, T3* b, int n)
    {
        //Dirichlet distribution
        double s = 0;
        for (int i = 0; i < n; ++i)
        {
            double v = Gamma((double)a[i] + (double)b[i]);
            res[i] = v;
            s += v;
        }
        Mul(res, (T1)(1.0 / s), n);
    }

    // Draw a real number from gamma distributio
    __device__ double RNG_CUDA<double>::Gamma(double alpha, double beta)
    {
        //gamma distribution
        if (alpha < 1)
        {
            //bug fixed on 20220816 to keep code sequence
            volatile double v1 = Gamma(1.0 + alpha, beta);
            volatile double v2 = pow(Uniform(), 1.0 / alpha);
            return v1 * v2;
        }
        double t, v, u;
        double d = alpha - 0.333333333333333;
        double c = 1.0 / (sqrt(d) * 3.0);

        for (;;)
        {
            do
            {
                t = Normal();
                volatile double v1 = c * t;
                v = 1.0 + v1;
            } while (v <= 0);

            v = v * v * v;
            u = Uniform();

            if (u < 1.0 - 0.0331 * t * t * t * t) break;
            if (log(u) < 0.5 * t * t + d * (1.0 - v + log(v))) break;
        }
        return beta * d * v;
    }

    // Draw a normal distriubted real number
    __device__ double RNG_CUDA<double>::Normal()
    {
        //normal distribution
        if (*(uint*)&U1 != 0)
        {
            double re = U1 * MySinDevice(U2);
            *(uint*)&U1 = 0;
            return re;
        }

        volatile double v1 = Uniform();
        volatile double v2 = Uniform();

        U1 = MySqrtDevice(-2.0 * log(v1));
        U2 = 2.0 * M_PI * v2;
        return U1 * MyCosDevice(U2);
    }

    // Draw a uniform distriubted real number
    __device__ double RNG_CUDA<double>::Uniform()
    {
        uint64 u = XorShift(), r = 0x3FF0000000000000;
        double& re = *(double*)&r;
        r |= u & 0x000FFFFFFFFFFFFF;
        return re - 1.0;
    }

    */

#endif

#ifndef _RNG_CUDA_FP32
    __device__ RNG_CUDA<float>::RNG_CUDA()
    {

    }

    /* Initialize rng */
    __device__ RNG_CUDA<float>::RNG_CUDA(uint64 seed, uint64 salt)
    {
        uint s = MurmurHash32_CUDA(seed, salt);

        x = 0x075BCD15 ^ (s);
        y = 0x159A55E5 ^ (s << 3);
        z = 0x1F123BB5 ^ (s << 6);
    }

    /* Draw a uniform distriubted interger */
    __device__ uint RNG_CUDA<float>::XorShift()
    {
        uint t;
        x ^= x << 16;
        x ^= x >> 5;
        x ^= x << 1;
        t = x;
        x = y;
        y = z;
        z = t ^ x ^ y;
        return z;
    }

    /* Draw a polynormial distriubted integer */
    __device__ int RNG_CUDA<float>::Poly(float* a, int n, bool typed)
    {
        //row unify
        float s = 0;

        if (typed)
        {
            volatile float v1 = (float)MIN_FREQ;
            for (int i = 0; i < n; ++i)
            {
                a[i] += v1;
                s += a[i];
            }
        }

        uint u = XorShift();

        if (typed)
        {
            uint r = 0x3F800000;
            float& re = *(float*)&r;
            r |= u & 0x007FFFFF;

            float t = (re - 1.0f) * s;

            for (int i = 0; i < n; ++i)
            {
                if (t < a[i]) return i;
                t -= a[i];
            }
        }
        return n - 1;
    }
    
    /*
    // Draw a vector from Dirichlet distribution D(a1 + b1, a2 + b2, ...)
    template<typename T1, typename T2, typename T3>
    __device__ void RNG_CUDA<float>::Dirichlet(T1* res, T2* a, T3* b, int n)
    {
        //Dirichlet distribution
        double s = 0;
        for (T3 i = 0; i < n; ++i)
        {
            double v = Gamma((double)a[i] + (double)b[i]);
            res[i] = (T1)v;
            s += v;
        }
        Mul(res, (T1)(1.0 / s), n);
    }

    // Draw a real number from gamma distribution
    __device__ double RNG_CUDA<float>::Gamma(double alpha, double beta)
    {
        //gamma distribution
        if (alpha < 1)
        {
            //bug fixed on 20220816 to keep code sequence
            volatile double v1 = Gamma(1.0 + alpha, beta);
            volatile double v2 = pow(Uniform(), 1.0 / alpha);
            return v1 * v2;
        }
        double t, v, u;
        double d = alpha - 0.333333333333333;
        double c = 1.0 / (sqrt(d) * 3.0);

        for (;;)
        {
            do
            {
                t = Normal();
                v = 1.0 + c * t;
            } while (v <= 0);

            v = v * v * v;
            u = Uniform();

            if (u < 1.0 - 0.0331 * t * t * t * t) break;
            if (log(u) < 0.5 * t * t + d * (1.0 - v + log(v))) break;
        }
        return beta * d * v;
    }

    // Draw a normal distriubted real number
    __device__ double RNG_CUDA<float>::Normal()
    {
        //normal distribution
        if (*(uint*)&U1 != 0)
        {
            double re = U1 * MySinDevice(U2);
            *(uint*)&U1 = 0;
            return re;
        }

        volatile double v1 = Uniform();
        volatile double v2 = Uniform();
        U1 = MySqrtDevice(-2.0 * log(std::max(MIN_FREQ, v1)));
        U2 = 2.0 * M_PI * v2;
        return U1 * MyCosDevice(U2);
    }

    // Draw a uniform distriubted real number
    __device__ float RNG_CUDA<float>::Uniform()
    {
        uint u = XorShift(), r = 0x3F800000;
        float& re = *(float*)&r;
        r |= u & 0x007FFFFF;
        return re - 1.0f;
    }

    */
#endif


/*=============================== EXPORT FUNCTIONS ===============================*/

TARGETCUDA void WaitCUDA()
{
    checkErrors(cudaStreamSynchronize(GPUstream));
    checkErrors(cudaGetLastError());
}

TARGETCUDA void* MallocHostCUDA(uint64 size)
{
    if (nGPU == 0) return Malloc(size);

    void* re = NULL;
    while (re == NULL)
    {
        cudaMallocHost((void**)&re, size);
        if (re) break;
        printf("Fail to allocate %0.3f Gib memory, increase virtual memory and retry (this will not terminate this process).", size / 1024.0 / 1024.0 / 1024.0);
        Pause();
    }
    return re;
}

TARGETCUDA void FreeHostCUDA(void* addr)
{
    if (nGPU == 0) return free(addr);
    checkErrors(cudaFreeHost(addr));
}

TARGETCUDA void* MallocDeviceCUDA(uint64 size)
{
    void* re = NULL;
    while (re == NULL)
    {
        cudaMalloc((void**)&re, size);
        if (re) break;
        printf("Fail to allocate %0.3f Gib memory, increase virtual memory and retry (this will not terminate this process).", size / 1024.0 / 1024.0 / 1024.0);
        Pause();
    }
    return re;
}

TARGETCUDA void FreeCUDA(void* addr)
{
    checkErrors(cudaFree(addr));
}

TARGETCUDA void MemsetCUDA(void* addr, int val, uint64 size)
{
    checkErrors(cudaMemsetAsync(addr, val, size, GPUstream));
}

TARGETCUDA void MemcpyCUDA(void* dst, void* src, uint64 size, bool todevice)
{
    if (todevice)
        cudaMemcpyAsync(dst, src, size, cudaMemcpyHostToDevice, GPUstream);
    else
        cudaMemcpyAsync(dst, src, size, cudaMemcpyDeviceToHost, GPUstream);
}

__global__ void SetDeviceMemory(MEMORY_CUDA* device_mem)
{
    cuda_mem_device = *device_mem;
}

TARGETCUDA void FreeStructureMemory(int devID)
{
    FreeCUDA(cuda_mem[devID].device_mem);
}

TARGETCUDA void CopyStructureMemory(int devID)
{
    checkErrors(cudaSetDevice(devID));
    checkErrors(cudaStreamCreate(&GPUstream));

    MEMORY_CUDA& cmem = cuda_mem[devID];
    cmem.KT = KT;
    cmem.N = nind;
    cmem.maxploidy = maxploidy;
    cmem.minploidy = minploidy;
    cmem.pad = 0;
    cmem.L = nloc;

    int64 gbsize = 0;
    // calculate genotype bucket size
    gbsize += Align(geno_bucket.coffset.load(), sizeof(uint64)); 

    // calculate genotype data size
    int64 gssize = 0;
    for (int64 l = 0; l < nloc; ++l)
    {
        int gsize = (int)GetLoc(l).ngeno;
        int gasize = (int)GetLoc(l).GetGenoAlleleSize();//genotype allele array
        gssize += gsize * sizeof(GENOTYPE) + gasize * sizeof(ushort);
    }

    // calculate total usage
    int64 tsize = 0;
    /*MEMORY*/              tsize += Align16(sizeof(MEMORY_CUDA));
    /*GBUCKET*/             tsize += Align16(gbsize);
    /*loc_addr*/            tsize += Align16(nloc * sizeof(uint64));
    /*allele_freq_offset*/  tsize += Align16(nloc * sizeof(int64));
    /*slocus*/              tsize += Align16(nloc * sizeof(SLOCUS_CUDA));
    /*geno allele*/         tsize += Align16(gssize);

    // alloc memory
    byte* p = (byte*)MallocDeviceCUDA(tsize);
    /*MEMORY*/              cmem.device_mem             = (MEMORY_CUDA*)p;      p += Align16(sizeof(MEMORY_CUDA));
    /*GBUCKET*/             cmem.g_bucket               = p;                    p += Align16(gbsize);
    /*loc_addr*/            cmem.loc_addr               = (uint64*)p;           p += Align16(nloc * sizeof(uint64));
    /*allele_freq_offset*/  cmem.allele_freq_offset     = (int64*)p;            p += Align16(nloc * sizeof(int64));
    /*slocus*/              cmem.sloc                   = (SLOCUS_CUDA*)p;      p += Align16(nloc * sizeof(SLOCUS_CUDA));

    // copy memory
    /*MEMORY*/          MemcpyCUDA(cmem.device_mem, &cmem, sizeof(MEMORY_CUDA), true);
    /*GBUCKET*/         MemcpyCUDA(cmem.g_bucket, geno_bucket.base_addr, geno_bucket.coffset.load(), true);

    /*loc_addr*/
    {
        OFFSET* offset = new OFFSET[nloc];
        byte* base_addr = cmem.g_bucket;
        for (int64 l = 0; l < nloc; ++l)
            offset[l] = OFFSET{ (uint64)base_addr + geno_bucket.offset[l].offset, geno_bucket.offset[l].size };
        MemcpyCUDA(cmem.loc_addr, offset, nloc * sizeof(uint64), true);
        WaitCUDA();
        DEL(offset);
    }

    /*allele_freq_offset*/
    {
        MemcpyCUDA(cmem.allele_freq_offset, allele_freq_offset, nloc * sizeof(int64), true);
        WaitCUDA();
    }

    /*slocus*/
    {
        SLOCUS_CUDA* sloc1 = new SLOCUS_CUDA[nloc];
        SLOCUS* sloc2 = (SLOCUS*)sloc1;
        memcpy(sloc1, slocus, nloc * sizeof(SLOCUS_CUDA));
        byte* locbufo = new byte[gssize];

        /*geno allele*/
#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
        for (int id = 0; id < g_nthread_val; ++id)
        {
            byte* p2 = p, * locbuf = locbufo;
            for (int64 l = 0; l < nloc; ++l)
            {
                int gsize = (int)GetLoc(l).ngeno;
                int gasize = (int)GetLoc(l).GetGenoAlleleSize();//genotype allele array
                int gssize1 = gsize * sizeof(GENOTYPE) + gasize * sizeof(ushort);

                if (l % g_nthread_val != id)
                {
                    locbuf += gssize1;
                    p2 += gssize1;
                    continue;
                }

                sloc2[l].bits1 = (uint64)p2;
                sloc2[l].flag_alen = 0;

                GENOTYPE* gtab2 = (GENOTYPE*)locbuf;
                ushort* gatab2 = (ushort*)(locbuf + gsize * sizeof(GENOTYPE));

                memcpy(gtab2, GetLoc(l).GetGtab(), gsize * sizeof(GENOTYPE));
                memcpy(gatab2, GetLoc(l).GetGenoAlleleArray(), gasize * sizeof(ushort));

                for (int gi = 0; gi < gsize; ++gi)
                {
                    GENOTYPE& gt2 = gtab2[gi];
                    if (gt2.Nalleles())
                    {
                        gt2.SetAlleleArray(gatab2);
                        gatab2 += gt2.Nalleles() + gt2.Ploidy();
                    }
                    else
                        gt2.SetAlleleArray((ushort*)-1);
                }

                locbuf += gssize1;
                p2 += gssize1;
            }
        }

        MemcpyCUDA(p, locbufo, gssize, true);
        MemcpyCUDA(cmem.sloc, sloc1, nloc * sizeof(SLOCUS_CUDA), true);

        SetDeviceMemory<<<1, 1, 0, GPUstream>>> (cmem.device_mem);
        WaitCUDA();

        DEL(locbufo);
        DEL(sloc1);
        checkErrors(cudaStreamDestroy(GPUstream));
    }
}

TARGETCUDA int GetDeviceCountCUDA()
{
    int devNum = 0;
    if (cudaGetDeviceCount(&devNum))
        devNum = 0;
    return devNum;
}

__global__ void SetDeviceProperty(cudaDeviceProp* mem_device)
{
    GPUprop_device = *mem_device;
}

TARGETCUDA void CreateStreamCUDA(int devID)
{
    checkErrors(cudaSetDevice(devID));
    checkErrors(cudaStreamCreate(&GPUstream));
    checkErrors(cudaGetDeviceProperties(&GPUprop, devID));

    cudaDeviceProp* mem_device = (cudaDeviceProp*)MallocDeviceCUDA(sizeof(cudaDeviceProp));
    MemcpyCUDA(mem_device, &GPUprop, sizeof(cudaDeviceProp), true);
    SetDeviceProperty<<<1, 1, 0, GPUstream>>> (mem_device);
    FreeCUDA(mem_device);
}

TARGETCUDA void DestroyStreamCUDA()
{
    checkErrors(cudaStreamDestroy(GPUstream));
}

TARGETCUDA int GetNumCores(cudaDeviceProp& devProp)
{
    int cores = 0;
    int mp = devProp.multiProcessorCount;
    switch (devProp.major) {
    case 2: // Fermi
        if (devProp.minor == 1) cores = mp * 48;
        else cores = mp * 32;
        break;
    case 3: // Kepler
        cores = mp * 192;
        break;
    case 5: // Maxwell
        cores = mp * 128;
        break;
    case 6: // Pascal
        if ((devProp.minor == 1) || (devProp.minor == 2)) cores = mp * 128;
        else if (devProp.minor == 0) cores = mp * 64;
        else printf("Unknown device type\n");
        break;
    case 7: // Volta and Turing
        if ((devProp.minor == 0) || (devProp.minor == 5)) cores = mp * 64;
        else printf("Unknown device type\n");
        break;
    case 8: // Ampere
        if (devProp.minor == 0) cores = mp * 64;
        else if (devProp.minor == 6) cores = mp * 128;
        else printf("Unknown device type\n");
        break;
    default:
        printf("Unknown device type\n");
        break;
    }
    return cores;
}

TARGETCUDA void AllocMemoryCUDA()
{
    cuda_mem = new MEMORY_CUDA[nGPU];
}

TARGETCUDA void FreeMemoryCUDA()
{
    DEL(cuda_mem);
}

TARGETCUDA void ShowDevicesCUDA()
{
    InitCUDA();
    printf("List of CUDA devices:");
    for (int dev = 0; dev < nGPU; ++dev)
    {
        cudaDeviceProp GPUprop2;
        cudaGetDeviceProperties(&GPUprop2, dev);
        //byte* t = (byte*)GPUprop2.uuid.bytes;
        //UUID = %02x%02x%02x%02x-%02x%02x-%02x%02x-%02x%02x-%02x%02x%02x%02x%02x%02x, 
        printf("\n  %d: %s, v%d.%d, #cores = %d, Memory = %0.1f Gib",
            dev, GPUprop2.name, GPUprop2.major, GPUprop2.minor,
            GetNumCores(GPUprop2),
            //t[ 0], t[ 1], t[ 2], t[ 3], 
            //t[ 4], t[ 5], t[ 6], t[ 7], 
            //t[ 8], t[ 9], t[10], t[11], 
            //t[12], t[13], t[14], t[15], 
            GPUprop2.totalGlobalMem / (1073741824.0));
    }
    printf("\n\n");
}

TARGETCUDA void ResetDeviceCUDA()
{
    for (int dev = 0; dev < nGPU; ++dev)
    {
        checkErrors(cudaSetDevice(dev));
        checkErrors(cudaDeviceReset());
    }
}

TARGETCUDA void Eig64CUDA(double* A, double* U, double* V, int64 n)
{
    cusolverDnHandle_t handle = nullptr;
    cudaStream_t stream = nullptr;
    double* dA = nullptr, * dV = nullptr, * dwork = nullptr;
    int* dinfo = nullptr, lwork = 0, info = 0;
    cusolverStatus_t cusolverStatus = CUSOLVER_STATUS_SUCCESS;
    cudaError_t cudaStatus = cudaSuccess;

    // Step 1: Create cusolver handle and bind a stream
    cusolverStatus = cusolverDnCreateA(&handle);
    if (cusolverStatus != CUSOLVER_STATUS_SUCCESS)
        Exit("\nError: Failed to create cusolver handle.\n");

    cudaStatus = cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
    if (cudaStatus != cudaSuccess)
        Exit("\nError: Failed to create CUDA stream.\n");
    cusolverDnSetStreamA(handle, stream);

    // Step 2: Allocate device memory for A, eigenvalues (V), and info
    size_t dA_size = n * n * sizeof(double);
    size_t dV_size = n * sizeof(double);
    size_t dinfo_size = sizeof(int);

    cudaStatus = cudaMalloc((void**)&dA, dA_size + dV_size + dinfo_size);
    if (cudaStatus != cudaSuccess)
        Exit("\nError: Failed to allocate device memory.\n");

    dV = dA + n * n;
    dinfo = (int*)(dV + n);

    // Step 3: Copy matrix A to the device
    cudaStatus = cudaMemcpy(dA, A, dA_size, cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess)
        Exit("\nError: Failed to copy matrix A to device.\n");

    // Step 4: Query working space for eigenvalue decomposition
    cusolverStatus = cusolverDnDsyevd_bufferSizeA(handle, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER, n, dA, n, dV, &lwork);
    if (cusolverStatus != CUSOLVER_STATUS_SUCCESS)
        Exit("\nError: Failed to query workspace size for eigenvalue decomposition.\n");

    cudaStatus = cudaMalloc((void**)&dwork, lwork * sizeof(double));
    if (cudaStatus != cudaSuccess)
        Exit("\nError: Failed to allocate device memory for workspace.\n");

    // Step 5: Perform eigenvalue decomposition
    cusolverStatus = cusolverDnDsyevdA(handle, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER, n, dA, n, dV, dwork, lwork, dinfo);
    if (cusolverStatus != CUSOLVER_STATUS_SUCCESS)
        Exit("\nError: Failed to perform eigenvalue decomposition.\n");

    // Step 6: Check for errors in the decomposition
    cudaStatus = cudaMemcpy(&info, dinfo, dinfo_size, cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess)
        Exit("\nError: Failed to copy info from device to host.\n");

    if (info != 0)
        Exit("\nError: Eigenvalue decomposition failed. Info = %d\n", info);

    // Step 7: Copy eigenvalues and eigenvectors back to host
    cudaStatus = cudaMemcpyAsync(V, dV, dV_size, cudaMemcpyDeviceToHost, stream);
    if (cudaStatus != cudaSuccess)
        Exit("\nError: Failed to copy eigenvalues to host.\n");

    cudaStatus = cudaMemcpyAsync(U, dA, dA_size, cudaMemcpyDeviceToHost, stream);
    if (cudaStatus != cudaSuccess)
        Exit("\nError: Failed to copy eigenvectors to host.\n");

    // Synchronize the stream to ensure all operations are complete
    cudaStreamSynchronize(stream);

    // Step 8: Free resources
    if (dA) cudaFree(dA);
    if (dwork) cudaFree(dwork);
    if (handle) cusolverDnDestroyA(handle);
    if (stream) cudaStreamDestroy(stream);
}

TARGETCUDA void Eig32CUDA(float* A, float* U, float* V, int64 n)
{
    cusolverDnHandle_t handle = nullptr;
    cudaStream_t stream = nullptr;
    float* dA = nullptr, * dV = nullptr, * dwork = nullptr;
    int* dinfo = nullptr, lwork = 0, info = 0;
    cusolverStatus_t cusolverStatus = CUSOLVER_STATUS_SUCCESS;
    cudaError_t cudaStatus = cudaSuccess;

    // Step 1: Create cusolver handle and bind a stream
    cusolverStatus = cusolverDnCreateA(&handle);
    if (cusolverStatus != CUSOLVER_STATUS_SUCCESS)
        Exit("\nError: Failed to create cusolver handle.\n");

    cudaStatus = cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
    if (cudaStatus != cudaSuccess)
        Exit("\nError: Failed to create CUDA stream.\n");
    cusolverDnSetStreamA(handle, stream);

    // Step 2: Allocate device memory for A, eigenvalues (V), and info
    size_t dA_size = n * n * sizeof(float);
    size_t dV_size = n * sizeof(float);
    size_t dinfo_size = sizeof(int);

    cudaStatus = cudaMalloc((void**)&dA, dA_size + dV_size + dinfo_size);
    if (cudaStatus != cudaSuccess)
        Exit("\nError: Failed to allocate device memory.\n");

    dV = dA + n * n;
    dinfo = (int*)(dV + n);

    // Step 3: Copy matrix A to the device
    cudaStatus = cudaMemcpy(dA, A, dA_size, cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess)
        Exit("\nError: Failed to copy matrix A to device.\n");

    // Step 4: Query working space for eigenvalue decomposition
    cusolverStatus = cusolverDnSsyevd_bufferSizeA(handle, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER, n, dA, n, dV, &lwork);
    if (cusolverStatus != CUSOLVER_STATUS_SUCCESS)
        Exit("\nError: Failed to query workspace size for eigenvalue decomposition.\n");

    cudaStatus = cudaMalloc((void**)&dwork, lwork * sizeof(float));
    if (cudaStatus != cudaSuccess)
        Exit("\nError: Failed to allocate device memory for workspace.\n");

    // Step 5: Perform eigenvalue decomposition
    cusolverStatus = cusolverDnSsyevdA(handle, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER, n, dA, n, dV, dwork, lwork, dinfo);
    if (cusolverStatus != CUSOLVER_STATUS_SUCCESS)
        Exit("\nError: Failed to perform eigenvalue decomposition.\n");

    // Step 6: Check for errors in the decomposition
    cudaStatus = cudaMemcpy(&info, dinfo, dinfo_size, cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess)
        Exit("\nError: Failed to copy info from device to host.\n");

    if (info != 0)
        Exit("\nError: Eigenvalue decomposition failed. Info = %d\n", info);

    // Step 7: Copy eigenvalues and eigenvectors back to host
    cudaStatus = cudaMemcpyAsync(V, dV, dV_size, cudaMemcpyDeviceToHost, stream);
    if (cudaStatus != cudaSuccess)
        Exit("\nError: Failed to copy eigenvalues to host.\n");

    cudaStatus = cudaMemcpyAsync(U, dA, dA_size, cudaMemcpyDeviceToHost, stream);
    if (cudaStatus != cudaSuccess)
        Exit("\nError: Failed to copy eigenvectors to host.\n");

    // Synchronize the stream to ensure all operations are complete
    cudaStreamSynchronize(stream);

    // Step 8: Free resources
    if (dA) cudaFree(dA);
    if (dwork) cudaFree(dwork);
    if (handle) cusolverDnDestroyA(handle);
    if (stream) cudaStreamDestroy(stream);
}

TARGETCUDA void Svd64CUDA(double* A, double* U, double* S, double* VT, int64 m, int64 n)
{
    cusolverDnHandle_t handle = NULL;
    cudaStream_t stream = NULL;

    int64 mn = std::min(m, n);
    double* dA = nullptr, *dS = nullptr, *dU = nullptr, *dV = nullptr, *dwork = nullptr;
    int* dinfo = nullptr, lwork = 0, info = 0;

    // Step 1: Create cusolver handle and bind a stream
    cusolverDnCreateA(&handle);
    cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
    cusolverDnSetStreamA(handle, stream);

    // Step 2: Query working space of SVD
    cusolverDnDgesvd_bufferSizeA(handle, m, n, &lwork);

    // Step 3: Allocate device memory and copy A to device
    cudaMalloc((void**)(&dA), sizeof(double) * (m * n + m * mn + mn + mn * n + lwork) + sizeof(int));
    dU = dA + m * n;
    dS = dU + m * mn;
    dV = dS + mn;
    dwork = dV + mn * n;
    dinfo = (int*)(dwork + lwork);

    cudaMemcpyAsync(dA, A, sizeof(double) * m * n, cudaMemcpyHostToDevice, stream);
    cudaMemsetAsync(dinfo, 0, sizeof(int), stream); // Initialize dinfo to zero

    // Step 4: Compute SVD
    cusolverStatus_t t1 = cusolverDnDgesvdA(handle, 'S', 'S', m, n, dA, m, dS, dU, m, dV, n, dwork, lwork, NULL, dinfo);
    if (t1 != CUSOLVER_STATUS_SUCCESS) Exit("\nError: cusolverDnDgesvd failed.\n");

    // Step 5: Copy results back to host
    cudaMemcpyAsync(U, dU, sizeof(double) * m * mn, cudaMemcpyDeviceToHost, stream);
    cudaMemcpyAsync(S, dS, sizeof(double) * mn, cudaMemcpyDeviceToHost, stream);
    cudaMemcpyAsync(VT, dV, sizeof(double) * mn * n, cudaMemcpyDeviceToHost, stream);
    cudaMemcpyAsync(&info, dinfo, sizeof(int), cudaMemcpyDeviceToHost, stream);

    // Synchronize the stream to ensure all operations are complete
    cudaStreamSynchronize(stream);
    if (info != 0) Exit("\nError: CUDA failed to perform singular value decomposition.\n");

    // Step 6: Free resources
    cudaFree(dA);
    cusolverDnDestroyA(handle);
    cudaStreamDestroy(stream);
}

TARGETCUDA void Svd32CUDA(float* A, float* U, float* S, float* VT, int64 m, int64 n)
{
    cusolverDnHandle_t handle = NULL;
    cudaStream_t stream = NULL;

    int64 mn = std::min(m, n);
    float* dA = nullptr, *dS = nullptr, *dU = nullptr, *dV = nullptr, *dwork = nullptr;
    int* dinfo = nullptr, lwork = 0, info = 0;

    // Step 1: Create cusolver handle and bind a stream
    cusolverDnCreateA(&handle);
    cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
    cusolverDnSetStreamA(handle, stream);

    // Step 2: Query working space of SVD
    cusolverDnSgesvd_bufferSizeA(handle, m, n, &lwork);

    // Step 3: Allocate device memory and copy A to device
    cudaMalloc((void**)(&dA), sizeof(float) * (m * n + m * mn + mn + mn * n + lwork) + sizeof(int));
    dU = dA + m * n;
    dS = dU + m * mn;
    dV = dS + mn;
    dwork = dV + mn * n;
    dinfo = (int*)(dwork + lwork);

    cudaMemcpyAsync(dA, A, sizeof(float) * m * n, cudaMemcpyHostToDevice, stream);
    cudaMemsetAsync(dinfo, 0, sizeof(int), stream); // Initialize dinfo to zero

    // Step 4: Compute SVD
    cusolverStatus_t t1 = cusolverDnSgesvdA(handle, 'S', 'S', m, n, dA, m, dS, dU, m, dV, n, dwork, lwork, NULL, dinfo);
    if (t1 != CUSOLVER_STATUS_SUCCESS) Exit("\nError: cusolverDnSgesvd failed.\n");

    // Step 5: Copy results back to host
    cudaMemcpyAsync(U, dU, sizeof(float) * m * mn, cudaMemcpyDeviceToHost, stream);
    cudaMemcpyAsync(S, dS, sizeof(float) * mn, cudaMemcpyDeviceToHost, stream);
    cudaMemcpyAsync(VT, dV, sizeof(float) * mn * n, cudaMemcpyDeviceToHost, stream);
    cudaMemcpyAsync(&info, dinfo, sizeof(int), cudaMemcpyDeviceToHost, stream);

    // Synchronize the stream to ensure all operations are complete
    cudaStreamSynchronize(stream);
    if (info != 0) Exit("\nError: CUDA failed to perform singular value decomposition.\n");

    // Step 6: Free resources
    cudaFree(dA);
    cusolverDnDestroyA(handle);
    cudaStreamDestroy(stream);
}

TARGETCUDA void MatrixMul64CUDA(double* A, double* B, double* res, int64 m, int64 n, int64 p, bool Atrans, bool Btrans)
{
    // Allocate device memory for matrices A, B, and res
    double* dA = nullptr, * dB = nullptr, * dres = nullptr;
    cublasHandle_t handle = nullptr;
    cublasStatus_t cublasStatus = CUBLAS_STATUS_SUCCESS;
    cudaError_t cudaStatus = cudaSuccess;

    // Step 1: Allocate device memory
    size_t dA_size = m * n * sizeof(double);
    size_t dB_size = n * p * sizeof(double);
    size_t dres_size = m * p * sizeof(double);

    cudaStatus = cudaMalloc((void**)&dA, dA_size + dB_size + dres_size);
    if (cudaStatus != cudaSuccess)
        Exit("\nError: Failed to allocate device memory for dA.\n");
    
    dB = dA + m * n;
    dres = dB + n * p;

    // Step 2: Copy matrices A and B from the host to the device
    cudaStatus = cudaMemcpy(dA, A, dA_size, cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess)
        Exit("\nError: Failed to copy matrix A to device.\n");

    cudaStatus = cudaMemcpy(dB, B, dB_size, cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess)
        Exit("\nError: Failed to copy matrix B to device.\n");

    // Step 3: Create a cuBLAS handle
    cublasStatus = cublasCreateA(&handle);
    if (cublasStatus != CUBLAS_STATUS_SUCCESS)
        Exit("\nError: Failed to create cuBLAS handle.\n");

    // Step 4: Perform matrix multiplication using cuBLAS
    double alpha = 1.0, beta = 0.0;
    cublasStatus = cublasDgemmA(
        handle,
        Atrans ? CUBLAS_OP_T : CUBLAS_OP_N,  // Transpose A if Atrans is true
        Btrans ? CUBLAS_OP_T : CUBLAS_OP_N,  // Transpose B if Btrans is true
        m,                                   // Number of rows in A (or A^T)
        p,                                   // Number of columns in B (or B^T)
        n,                                   // Number of columns in A (or A^T) and rows in B (or B^T)
        &alpha,                              // Scalar multiplier for A * B
        dA,                                  // Pointer to matrix A on device
        Atrans ? n : m,                      // Leading dimension of A (or A^T)
        dB,                                  // Pointer to matrix B on device
        Btrans ? p : n,                      // Leading dimension of B (or B^T)
        &beta,                               // Scalar multiplier for result matrix
        dres,                                // Pointer to result matrix on device
        m                                    // Leading dimension of result matrix
    );

    if (cublasStatus != CUBLAS_STATUS_SUCCESS)
        Exit("\nError: cuBLAS failed to perform matrix multiplication.\n");

    // Step 5: Copy the result matrix from device to host
    cudaStatus = cudaMemcpy(res, dres, dres_size, cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess)
        Exit("\nError: Failed to copy result matrix from device to host.\n");

    // Step 6: Free resources
    if (dA) cudaFree(dA);
    if (handle) cublasDestroyA(handle);
}

TARGETCUDA void MatrixMul32CUDA(float* A, float* B, float* res, int64 m, int64 n, int64 p, bool Atrans, bool Btrans)
{
    // Allocate device memory for matrices A, B, and res
    float* dA = nullptr, * dB = nullptr, * dres = nullptr;
    cublasHandle_t handle = nullptr;
    cublasStatus_t cublasStatus = CUBLAS_STATUS_SUCCESS;
    cudaError_t cudaStatus = cudaSuccess;

    // Step 1: Allocate device memory
    size_t dA_size = m * n * sizeof(float);
    size_t dB_size = n * p * sizeof(float);
    size_t dres_size = m * p * sizeof(float);

    cudaStatus = cudaMalloc((void**)&dA, dA_size + dB_size + dres_size);
    if (cudaStatus != cudaSuccess)
        Exit("\nError: Failed to allocate device memory for dA.\n");
    
    dB = dA + m * n;
    dres = dB + n * p;

    // Step 2: Copy matrices A and B from the host to the device
    cudaStatus = cudaMemcpy(dA, A, dA_size, cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess)
        Exit("\nError: Failed to copy matrix A to device.\n");

    cudaStatus = cudaMemcpy(dB, B, dB_size, cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess)
        Exit("\nError: Failed to copy matrix B to device.\n");

    // Step 3: Create a cuBLAS handle
    cublasStatus = cublasCreateA(&handle);
    if (cublasStatus != CUBLAS_STATUS_SUCCESS)
        Exit("\nError: Failed to create cuBLAS handle.\n");

    // Step 4: Perform matrix multiplication using cuBLAS
    float alpha = 1.0, beta = 0.0;
    cublasStatus = cublasSgemmA(
        handle,
        Atrans ? CUBLAS_OP_T : CUBLAS_OP_N,  // Transpose A if Atrans is true
        Btrans ? CUBLAS_OP_T : CUBLAS_OP_N,  // Transpose B if Btrans is true
        m,                                   // Number of rows in A (or A^T)
        p,                                   // Number of columns in B (or B^T)
        n,                                   // Number of columns in A (or A^T) and rows in B (or B^T)
        &alpha,                              // Scalar multiplier for A * B
        dA,                                  // Pointer to matrix A on device
        Atrans ? n : m,                      // Leading dimension of A (or A^T)
        dB,                                  // Pointer to matrix B on device
        Btrans ? p : n,                      // Leading dimension of B (or B^T)
        &beta,                               // Scalar multiplier for result matrix
        dres,                                // Pointer to result matrix on device
        m                                    // Leading dimension of result matrix
    );

    if (cublasStatus != CUBLAS_STATUS_SUCCESS)
        Exit("\nError: cuBLAS failed to perform matrix multiplication.\n");

    // Step 5: Copy the result matrix from device to host
    cudaStatus = cudaMemcpy(res, dres, dres_size, cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess)
        Exit("\nError: Failed to copy result matrix from device to host.\n");

    // Step 6: Free resources
    if (dA) cudaFree(dA);
    if (handle) cublasDestroyA(handle);
}

/*=============================== DEVICE FUNCTIONS ===============================*/

template<typename REAL>
__device__ double SumProdDiv(REAL* A1, REAL* A2, REAL* B, int64 sep, int n)
{
    double re1 = 0, re2 = 0;
    for (int i = 0; i < n; ++i, A1++, A2++, B += sep)
    {
        volatile double v1 = (double)* A1 * (double)*B, v2 = (double)*A2 * (double)*B;
        __threadfence_block();
        re1 += v1;
        re2 += v2;
    }
    return re1 / re2;
}

__device__ float SumProdDivx(float* A1, float* A2, float* B, int64 sep, int n)
{
    float re1 = 0, re2 = 0;
    for (int i = 0; i < n; ++i, A1++, A2++, B += sep)
    {
        volatile float v1 = *A1 * *B, v2 = *A2 * *B;
        __threadfence_block();
        re1 += v1;
        re2 += v2;
    }
    return re1 / re2;
}

template<typename REAL>
__device__ double SumProd(REAL* A, REAL* B, int64 sep, int n)
{
    double re = 0;
    for (int i = 0; i < n; ++i, A++, B += sep)
    {
        volatile double v1 = (double)*A * (double)*B;
        __threadfence_block();
        re += v1;
    }
    return re;
}

__device__ float SumProdx(float* A, float* B, int64 sep, int n)
{
    float re = 0;
    for (int i = 0; i < n; ++i, A++, B += sep)
    {
        volatile float v1 = *A * *B;
        __threadfence_block();
        re += v1;
    }
    return re;
}

__device__ int ReduceMax(volatile int* A)
{
    int tid = threadIdx.x;
    __syncthreads();
    if (tid < 32) A[tid] = max(A[tid], A[tid + 32]);
    if (tid < 16) A[tid] = max(A[tid], A[tid + 16]);
    if (tid <  8) A[tid] = max(A[tid], A[tid +  8]);
    if (tid <  4) A[tid] = max(A[tid], A[tid +  4]);
    if (tid <  2) A[tid] = max(A[tid], A[tid +  2]);
    if (tid <  1) A[tid] = max(A[tid], A[tid +  1]);
    __syncthreads();
    return A[0];
}

__device__ double ReduceSum(volatile double* A)
{
    int tid = threadIdx.x;
    __syncthreads();
    if (tid < 32) A[tid] += A[tid + 32];
    if (tid < 16) A[tid] += A[tid + 16];
    if (tid <  8) A[tid] += A[tid +  8];
    if (tid <  4) A[tid] += A[tid +  4];
    if (tid <  2) A[tid] += A[tid +  2];
    if (tid <  1) A[tid] += A[tid +  1];
    __syncthreads();
    return A[0];
}

__device__ double ReduceProd(volatile double* A)
{
    int tid = threadIdx.x;
    __syncthreads();
    if (tid < 32) A[tid] *= A[tid + 32];
    if (tid < 16) A[tid] *= A[tid + 16];
    if (tid <  8) A[tid] *= A[tid +  8];
    if (tid <  4) A[tid] *= A[tid +  4];
    if (tid <  2) A[tid] *= A[tid +  2];
    if (tid <  1) A[tid] *= A[tid +  1];
    __syncthreads();
    return A[0];
}

#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600

#else
__device__ double atomicAdd(double* address, double val)
{
    uint64* addr = (uint64*)address;
    uint64 old = *addr, assumed;
    do
    {
        assumed = old;
        old = atomicCAS(addr, assumed, 
            __double_as_longlong(val + __longlong_as_double(assumed)));
    } while (assumed != old);

    return __longlong_as_double(old);
}

#endif

template<typename REAL>
__device__ REAL atomicMul(REAL* address, REAL val)
{
    if constexpr (std::is_same_v<REAL, double>)
    {
        uint64* addr = (uint64*)address;
        uint64 old = *addr, assumed;
        do
        {
            assumed = old;
            old = atomicCAS(addr, assumed, __double_as_longlong(val * __longlong_as_double(assumed)));
        } while (assumed != old);

        return __longlong_as_double(old);
    }
    else
    {
        int* address_as_int = (int*)address;
        int old = *address_as_int, assumed;
        do
        {
            assumed = old;
            old = atomicCAS(address_as_int, assumed, __float_as_int(val * __float_as_int(assumed)));
        } while (assumed != old); return __int_as_float(old);
    }
}

template<typename REAL>
__device__ void AddExponent1AtomicDevice(int64& slog2, REAL& val)
{
    if constexpr (std::is_same_v<REAL, double>)
    {
        int64& vv = *(int64*)&val;

        // add the exponent to slog2
        //slog2 += ((vv & 0x7FF0000000000000) >> 52) - 1023;
        atomicAdd((uint64*)&slog2, ((vv & 0x7FF0000000000000) >> 52) - 1023);

        // set the exponent of val to 0
        vv = (vv & 0x800FFFFFFFFFFFFF) | 0x3FF0000000000000;
    }
    else
    {
        int& vv = *(int*)&val;

        // add the exponent to slog2
        //slog2 += ((vv & 0x7FF0000000000000) >> 52) - 1023;
        atomicAdd((uint*)&slog2, ((vv & 0x7F800000) >> 23) - 127);

        // set the exponent of val to 0
        vv = (vv & 0x807FFFFF) | 0x3F800000;
    }
}

__device__ void AddExponent2AtomicDevice(int64& slog2, double& val)
{
    //fetch value and set to 1.0
    uint64 vv = atomicExch((uint64*)&val, 0x3FF0000000000000);

    //the thread get value
    if (vv != 0x3FF0000000000000)
    {
        //int64& vv = *(int64*)&val;

        // add the exponent to slog2
        //slog2 += ((vv & 0x7FF0000000000000) >> 52) - 1023;
        atomicAdd((uint64*)&slog2, ((vv & 0x7FF0000000000000) >> 52) - 1023);

        // set the exponent of val to 0
        vv = (vv & 0x800FFFFFFFFFFFFF) | 0x3FF0000000000000;

        atomicMul(&val, __longlong_as_double(vv));
    }
}

template<typename REAL>
__device__ void ChargeLogAtomicDevice(int64& slog, double& prod, REAL val)
{
    if constexpr (std::is_same_v<REAL, double>)
        if (val < DOUBLE_UNDERFLOW || val > DOUBLE_OVERFLOW) [[unlikely]]
            AddExponent1AtomicDevice(slog, val);

    //prod *= val;
    atomicMul(&prod, (double)val);

    if (prod < DOUBLE_UNDERFLOW || prod > DOUBLE_OVERFLOW) [[unlikely]]
        AddExponent2AtomicDevice(slog, prod);
}

template<typename REAL>
__device__ void ChargeLogAtomicDevice(int64* slog, double* prod, REAL* val, int64 n, int64 sep)
{
    for (int64 i = 0; i < n; ++i, val += sep)
        ChargeLogAtomicDevice(slog[i], prod[i], *val);
}

template<typename REAL>
__device__ void AddExponentDevice(int64& slog2, REAL& val)
{
    if constexpr (std::is_same_v<REAL, double>)
    {
        int64& vv = *(int64*)&val;

        // add the exponent to slog2
        slog2 += ((vv & 0x7FF0000000000000) >> 52) - 1023;

        // set the exponent of val to 0
        vv = (vv & 0x800FFFFFFFFFFFFF) | 0x3FF0000000000000;
    }
    else
    {
        int& vv = *(int*)&val;

        // add the exponent to slog2
        slog2 += ((vv & 0x7F800000) >> 23) - 127;

        // set the exponent of val to 0
        vv = (vv & 0x807FFFFF) | 0x3F800000;
    }
}

template<typename REAL>
__device__ void ChargeLogDevice(int64& slog, double& prod, REAL val)
{
    if constexpr (std::is_same_v<REAL, double>)
        if (val < DOUBLE_UNDERFLOW || val > DOUBLE_OVERFLOW) [[unlikely]]
            AddExponentDevice(slog, val);

    prod *= val;

    if (prod < DOUBLE_UNDERFLOW || prod > DOUBLE_OVERFLOW) [[unlikely]]
        AddExponentDevice(slog, prod);
}

template<typename REAL>
__device__ void ChargeLogDevice(int64* slog, double* prod, REAL* val, int64 n, int64 sep)
{
    for (int64 i = 0; i < n; ++i, val += sep)
        ChargeLogDevice(slog[i], prod[i], *val);
}

__device__ void OpenLogDevice(int64& slog, double& prod)
{
    slog = 0;
    prod = 1.0;
}

__device__ void OpenLogDevice(int64* slog, double* prod, int n)
{
    for (int i = 0; i < n; ++i)
        OpenLogDevice(slog[i], prod[i]);
}

__device__ void CloseLogDevice(int64& slog, double& prod)
{
    double& slog2 = *(double*)&slog;
    prod = slog2 = (slog + log2(prod)) * 0.693147180559945;
}

__device__ void CloseLogDevice(int64* slog, double* prod, int n)
{
    for (int i = 0; i < n; ++i)
        CloseLogDevice(slog[i], prod[i]);
}

/*=============================== FUNCTIONS ===============================*/

/*
template<typename REAL, bool fmodel, bool fsame>
__global__ void UpdatePCUDAKernel(BAYESIAN<REAL>* bayes, int64 m)
{
    uint64 l = blockDim.x * blockIdx.x + threadIdx.x;
    if (l >= cuda_mem_device.L) return;

    SLOCUS_CUDA* slocus = cuda_mem_device.sloc;
    int KT = cuda_mem_device.KT;
    int K = bayes->K;
    int k2 = GetLoc(l).GetK();
    int* ni = bayes->Ni_CUDA + cuda_mem_device.allele_freq_offset[l];
    REAL* p0 = bayes->Freq_CUDA + cuda_mem_device.allele_freq_offset[l];
    RNG_CUDA<REAL>rng (m * cuda_mem_device.L + l + bayes->seed, RNG_SALT_UPDATEP);

    if constexpr (fmodel)
    {
        REAL* f = bayes->f_CUDA, *pa = bayes->FreqA_CUDA + cuda_mem_device.allele_freq_offset[l];
        if constexpr (fsame)
        {
            REAL* p = p0, f0 = f[0];

            //Mul(ClusterFreq(0), AncestralFreqA, f[0], KT);
            for (int kk = 0; kk < k2; ++kk)
            {
                REAL pkk = p[kk] = pa[kk] * f0;

                p += KT;
                for (int k = 1; k < K; ++k, p += KT)
                {
                    //SetVal(ClusterFreq(k), ClusterFreq(0), KT);
                    p[kk] = pkk;
                }
            }
        }
        else
        {
            REAL* p = p0;
            for (int k = 0; k < K; ++k, p += KT)
            {
                //Mul(ClusterFreq(k), AncestralFreqA, f[k], KT);
                for (int kk = 0; kk < k2; ++kk)
                    p[kk] = pa[kk] * f[k];
            }
        }
    }
    else
    {
        REAL* p = p0, *Lambda = bayes->Lambda_CUDA;
        for (int k = 0; k < K; ++k, p += KT)
        {
            //SetVal(ClusterFreq(k), (REAL)Lambda[k], KT);
            REAL Lambdak = (REAL)Lambda[k];
            for (int kk = 0; kk < k2; ++kk)
                p[kk] = Lambdak;
        }
    }

    REAL* p = p0;
    for (int k = 0; k < K; ++k, ni += KT, p += KT)
        rng.Dirichlet(p, p, ni, k2);
}

template<typename REAL>
TARGET void BAYESIAN<REAL>::UpdatePCUDA()
{
    int nthread = BLOCK_DIM;
    int nblock = (nloc + BLOCK_DIM - 1) / BLOCK_DIM;

    if (sizeof(REAL) == 8)
    {
        if (fmodel && fsame)
            UpdatePCUDAKernel<double, true , true > <<<nblock, nthread, 0, GPUstream >>> (bayes_CUDA, m);
        else if (fmodel)
            UpdatePCUDAKernel<double, true , false> <<<nblock, nthread, 0, GPUstream >>> (bayes_CUDA, m);
        else
            UpdatePCUDAKernel<double, false, false> <<<nblock, nthread, 0, GPUstream >>> (bayes_CUDA, m);
    }
    else
    {
        if (fmodel && fsame)
            UpdatePCUDAKernel<float , true , true > <<<nblock, nthread, 0, GPUstream >>> (bayes_CUDA, m);
        else if (fmodel)
            UpdatePCUDAKernel<float , true , false> <<<nblock, nthread, 0, GPUstream >>> (bayes_CUDA, m);
        else
            UpdatePCUDAKernel<float , false, false> <<<nblock, nthread, 0, GPUstream >>> (bayes_CUDA, m);
    }

    WaitCUDA();
}
*/

template<typename REAL>
__global__ void UpdateZNoAdmixCUDAKernel(BAYESIAN<REAL>* bayes)
{
    uint64 l = blockDim.x * blockIdx.x + threadIdx.x;
    int64 L = cuda_mem_device.L;
    if (l >= L) return;

    int KT = cuda_mem_device.KT;
    int N = cuda_mem_device.N;
    ushort* Z = bayes->Z_CUDA;

    SLOCUS_CUDA* slocus = cuda_mem_device.sloc;
    int* Ni = bayes->Ni_CUDA + cuda_mem_device.allele_freq_offset[l];
    GENO_READER_CUDA rt(cuda_mem_device.loc_addr[l]);
    GENOTYPE_CUDA* gtab = GetLoc(l).GetGtab();

    for (int i = 0; i < N; ++i)
    {
        GENOTYPE_CUDA& gt = gtab[rt.Read()];
        if (gt.Nalleles() == 0) continue;
        ushort* als = gt.GetAlleleArray();
        int* ni = Ni + Z[i] * KT;

        for (int a = 0, vi = gt.Ploidy(); a < vi; ++a)
            ni[als[a]] ++;
    }
}

template<typename REAL>
TARGET void BAYESIAN<REAL>::UpdateZNoAdmixCUDA()
{
    SetZero(Mi, N * K);
    SetZero(Ni, K * KT);

    // Z in already updated in update Q
    //Count number of allele copies in individual i and cluster k
    for (int i = 0; i < N; ++i)
        Mi[i * K + Z[i]] = ainds<REAL>[i]->vt;

    int nthread = BLOCK_DIM;
    int nblock = (nloc + BLOCK_DIM - 1) / BLOCK_DIM;

    MemcpyCUDA(Z_CUDA, Z, N * sizeof(ushort), true);
    MemsetCUDA(Ni_CUDA, 0, K * KT * sizeof(int));

    UpdateZNoAdmixCUDAKernel <<<nblock, nthread, 0, GPUstream >>> (bayes_CUDA);

    MemcpyCUDA(Ni, Ni_CUDA, K * KT * sizeof(int), false);

    WaitCUDA();
}

//////////////////////////////////////////////////////////////////////////////////

template<typename REAL, bool fast_fp32>
__global__ void UpdateZAdmixCUDAKernel(BAYESIAN<REAL>* bayes, int64 m, bool useshare)
{
    extern __shared__ byte sdata[];

    uint64 l = blockDim.x * blockIdx.x + threadIdx.x;
    if (l >= cuda_mem_device.L) return;

    int KT = cuda_mem_device.KT;
    int N = cuda_mem_device.N;

    ///////////////////////////////////////////////////////////////////

    SLOCUS_CUDA* slocus = cuda_mem_device.sloc;
    int64* allele_freq_offset = cuda_mem_device.allele_freq_offset;
    int* Ni = bayes->Ni_CUDA + allele_freq_offset[l];
    REAL* q = bayes->Q_CUDA;
    uint64* mi = (uint64*)bayes->Mi_CUDA;
    REAL* Freq = bayes->Freq_CUDA;
    int K = bayes->K;
    RNG_CUDA<double> rngd; RNG_CUDA<float> rngs;
    int maxploidy = cuda_mem_device.maxploidy;

    if constexpr (std::is_same_v<REAL, double> || !fast_fp32)
        new (&rngd) RNG_CUDA<double>(m * cuda_mem_device.L + l + bayes->seed, RNG_SALT_UPDATEZ);
    else
        new (&rngs) RNG_CUDA<float >(m * cuda_mem_device.L + l + bayes->seed, RNG_SALT_UPDATEZ);

    GENO_READER_CUDA rt(cuda_mem_device.loc_addr[l]);
    GENOTYPE_CUDA* gtab = GetLoc(l).GetGtab();

    double gbuf[256];
    double* bufkd = useshare ? 
        (double*)(sdata + sizeof(double) * threadIdx.x * K) :
        (double*)gbuf;
    float* bufks = (float*)bufkd;

    for (int i = 0; i < N; ++i, q += K, mi += K)
    {
        GENOTYPE_CUDA& gt = gtab[rt.Read()];
        int ploidy = 0; 
        ushort* als;
        if (gt.Nalleles() > 0)
        {
            ploidy = gt.Ploidy();
            als = gt.GetAlleleArray();
        }

        for (int a = 0; a < maxploidy; ++a)
        {
            bool typed = a < ploidy;

            //The same to previous allele, do not update multinomial prob to save time
            if (typed && (a == 0 || als[a] != als[a - 1]))
                for (int k = 0; k < K; ++k)
                {
                    if constexpr (std::is_same_v<REAL, double> || !fast_fp32)
                        bufkd[k] = (double)q[k] * (double)ClusterAlleleFreq(k, l, als[a]);
                    else
                        bufks[k] = q[k] * ClusterAlleleFreq(k, l, als[a]);
                }

            ushort k2;
            //draw cluster for each allele copy
            if constexpr (std::is_same_v<REAL, double> || !fast_fp32)
                k2 = (ushort)rngd.Poly(bufkd, K, typed);
            else
                k2 = (ushort)rngs.Poly(bufks, K, typed);

            //Update Mi, NI
            if (typed)
            {
                atomicAdd(&mi[k2], 1);
                Ni[k2 * KT + als[a]]++;
            }
        }
    }
}

template<typename REAL>
TARGET void BAYESIAN<REAL>::UpdateZAdmixCUDA()
{
    SetZero(Mi, N * K);
    SetZero(Ni, K * KT);

    int nthread = 64;
    int block_share_size = sizeof(double) * K * nthread;
    bool useshare = GPUprop.sharedMemPerBlock > block_share_size;
    block_share_size = useshare ? block_share_size : 0;
    int nblock = (nloc + nthread - 1) / nthread;

    MemcpyCUDA(Z_CUDA, Z, N * sizeof(ushort), true);
    MemcpyCUDA(Q_CUDA, Q, N * K * sizeof(REAL), true);

    MemsetCUDA(Ni_CUDA, 0, K * KT * sizeof(int));
    MemsetCUDA(Mi_CUDA, 0, N * K * sizeof(int64));

    if (g_fastsingle_val == 1)
        UpdateZAdmixCUDAKernel<REAL, true > <<<nblock, nthread, block_share_size, GPUstream>>> (bayes_CUDA, m, useshare);
    else
        UpdateZAdmixCUDAKernel<REAL, false> <<<nblock, nthread, block_share_size, GPUstream>>> (bayes_CUDA, m, useshare);

    MemcpyCUDA(Ni, Ni_CUDA, K * KT * sizeof(int), false);
    MemcpyCUDA(Mi, Mi_CUDA, N * K * sizeof(int64), false);

    WaitCUDA();
}

//////////////////////////////////////////////////////////////////////////////////

template<typename REAL>
__global__ void UpdateQNoAdmixCUDAKernel(BAYESIAN<REAL>* bayes)
{
    extern __shared__ byte sdata[];
    double* buf2s = (double*)sdata;
    double& buf2t = buf2s[threadIdx.x];
    buf2t = 1;

    uint64 l = blockDim.x * blockIdx.x + threadIdx.x;

    if (l >= cuda_mem_device.L) return;

    int KT = cuda_mem_device.KT;
    int N = cuda_mem_device.N;
    int K = bayes->K;
    int maxploidy = cuda_mem_device.maxploidy;
    SLOCUS_CUDA* slocus = cuda_mem_device.sloc;

    int64* buf1 = (int64*)bayes->bufNK1_CUDA;
    double* buf2 = bayes->bufNK2_CUDA;
    REAL* p = bayes->Freq_CUDA + cuda_mem_device.allele_freq_offset[l];
    GENO_READER_CUDA rt(cuda_mem_device.loc_addr[l]);
    GENOTYPE_CUDA* gtab = GetLoc(l).GetGtab();

    for (int i = 0; i < N; ++i, buf1 += K, buf2 += K)
    {
        GENOTYPE_CUDA& gt = gtab[rt.Read()];
        int ploidy = 0; ushort* als;
        if (gt.Nalleles() > 0)
        {
            ploidy = gt.Ploidy();
            als = gt.GetAlleleArray();
        }

        for (int a = 0; a < maxploidy; ++a)
        {
            bool typed = a < ploidy;
            REAL* p2; if (typed) p2 = p + als[a];

            for (int k = 0; k < K; ++k, p2 += KT)
            {
                buf2t = typed ? *p2 : 1;

                double v1 = ReduceProd(buf2s);

                if (threadIdx.x == 0)
                    ChargeLogAtomicDevice(buf1[k], buf2[k], v1);
            }
        }
    }
}

template<typename REAL>
TARGET void BAYESIAN<REAL>::UpdateQNoAdmixCUDA()
{
    SetZero(Q, N * K);
    OpenLog((int64*)bufNK1, bufNK2, N * K);

    //add priori probability
    double* buf1 = bufNK1, * buf2 = bufNK2;
    if (locpriori) for (int i = 0; i < N; ++i, buf1 += K, buf2 += K)
    {
        if (ainds<REAL>[i]->vt == 0) continue;
        ChargeLog((int64*)buf1, buf2, Gamma + ainds<REAL>[i]->popid * K, K);
    }

    int nthread = 64;
    int block_share_size = sizeof(double) * nthread;
    int nblock = (nloc + nthread - 1) / nthread;

    MemcpyCUDA(bufNK1_CUDA, bufNK1, N * K * sizeof(double), true);
    MemcpyCUDA(bufNK2_CUDA, bufNK2, N * K * sizeof(double), true);

    UpdateQNoAdmixCUDAKernel <<<nblock, nthread, block_share_size, GPUstream >>> (bayes_CUDA);

    MemcpyCUDA(bufNK1, bufNK1_CUDA, N * K * sizeof(double), false);
    MemcpyCUDA(bufNK2, bufNK2_CUDA, N * K * sizeof(double), false);

    WaitCUDA();

    CloseLog((int64*)bufNK1, bufNK2, N * K);
    buf1 = bufNK1;
    REAL* q = Q;
    RNG<double> rng(seed + m, RNG_SALT_UPDATEQ);//double
    for (int i = 0; i < N; ++i, buf1 += K, q += K)
    {
        if (ainds<REAL>[i]->vt == 0) continue;
        ushort k2 = (ushort)rng.PolyLog(buf1, K);
        q[k2] = 1;
        Z[i] = k2;
    }
}

//////////////////////////////////////////////////////////////////////////////////

template<typename REAL, bool fast_fp32>
__global__ void UpdateQMetroCUDAKernel(BAYESIAN<REAL>* bayes)
{
    extern __shared__ byte sdata[];
    double* buf2s = (double*)sdata;
    double& buf2t = buf2s[threadIdx.x];
    buf2t = 1;

    uint64 l = blockDim.x * blockIdx.x + threadIdx.x;
    if (l >= cuda_mem_device.L) return;

    int KT = cuda_mem_device.KT;
    int N = cuda_mem_device.N;
    int K = bayes->K;
    int maxploidy = cuda_mem_device.maxploidy;
    SLOCUS_CUDA* slocus = cuda_mem_device.sloc;

    int64* buf1 = (int64*)bayes->bufN1_CUDA;
    double* buf2 = (double*)bayes->bufN2_CUDA;
    REAL* p = bayes->Freq_CUDA + cuda_mem_device.allele_freq_offset[l];
    REAL* q = bayes->Q_CUDA;
    REAL* bufi = (REAL*)bayes->bufNK1_CUDA;

    GENO_READER_CUDA rt(cuda_mem_device.loc_addr[l]);
    GENOTYPE_CUDA* gtab = GetLoc(l).GetGtab();

    for (int i = 0; i < N; ++i, q += K, bufi += K, buf1++, buf2++)
    {
        GENOTYPE_CUDA& gt = gtab[rt.Read()];
        int ploidy = 0; ushort* als;
        if (gt.Nalleles() > 0)
        {
            ploidy = gt.Ploidy();
            als = gt.GetAlleleArray();
        }

        for (int a = 0; a < maxploidy; ++a)
        {
            if constexpr (std::is_same_v<REAL, double> || !fast_fp32)
                buf2t = a < ploidy ? SumProdDiv(bufi, q, p + als[a], KT, K) : 1;
            else
                buf2t = a < ploidy ? SumProdDivx(bufi, q, p + als[a], KT, K) : 1;

            double v1 = ReduceProd(buf2s);

            if (threadIdx.x == 0)
                ChargeLogAtomicDevice(*buf1, *buf2, v1);
        }
    }
}

template<typename REAL>
TARGET void BAYESIAN<REAL>::UpdateQMetroCUDA()
{
    RNG<double> rng(seed + m, RNG_SALT_UPDATEQ);//REAL
    REAL* bufi = (REAL*)bufNK1;
    REAL* q = NULL;

    for (int i = 0; i < N; ++i, bufi += K)
    {
        if (ainds<REAL>[i]->vt == 0) continue;
        if (locpriori) rng.Dirichlet(bufi, AlphaLocal + ainds<REAL>[i]->popid * K, K);
        else           rng.Dirichlet(bufi, Alpha, K);
    }
    OpenLog((int64*)bufN1, bufN2, (int64)N);

    int nthread = 64;
    int block_share_size = sizeof(double) * nthread;
    int nblock = (nloc + nthread - 1) / nthread;

    MemcpyCUDA(Q_CUDA, Q, N * K * sizeof(REAL), true);
    MemcpyCUDA(bufNK1_CUDA, bufNK1, N * K * sizeof(REAL), true);
    MemsetCUDA(bufN1_CUDA, 0, N * sizeof(int64));
    MemcpyCUDA(bufN2_CUDA, bufN2, N * sizeof(double), true);

    if (g_fastsingle_val == 1)
        UpdateQMetroCUDAKernel<REAL, true > <<<nblock, nthread, block_share_size, GPUstream >>> (bayes_CUDA);
    else
        UpdateQMetroCUDAKernel<REAL, false> <<<nblock, nthread, block_share_size, GPUstream >>> (bayes_CUDA);

    MemcpyCUDA(bufN1, bufN1_CUDA, N * sizeof(double), false);
    MemcpyCUDA(bufN2, bufN2_CUDA, N * sizeof(double), false);

    WaitCUDA();

    CloseLog((int64*)bufN1, bufN2, (int64)N);

    bufi = (REAL*)bufNK1; q = Q;
    for (int i = 0; i < N; ++i, q += K, bufi += K)
    {
        if (ainds<REAL>[i]->vt == 0) continue;
        if (bufN1[i] >= NZERO || rng.Uniform() < exp(bufN1[i]))
            SetVal(q, bufi, K);
    }
}

//////////////////////////////////////////////////////////////////////////////////

template<typename REAL, bool isadmix, bool fast_fp32>
__global__ void RecordCUDAKernel(BAYESIAN<REAL>* bayes)
{
    extern __shared__ byte sdata[];
    double* buf2s = (double*)(sdata + sizeof(int) * blockDim.x);
    double& buf2t = buf2s[threadIdx.x];
    buf2t = 1;

    uint64 l = blockDim.x * blockIdx.x + threadIdx.x;
    if (l >= cuda_mem_device.L) return;

    int KT = cuda_mem_device.KT;
    int N = cuda_mem_device.N;
    int K = bayes->K;
    int maxploidy = cuda_mem_device.maxploidy;

    SLOCUS_CUDA* slocus = cuda_mem_device.sloc;
    ushort* Z;
    if constexpr (!isadmix) Z = bayes->Z_CUDA;

    REAL* p = bayes->Freq_CUDA + cuda_mem_device.allele_freq_offset[l], *q = bayes->Q_CUDA;
    int64 slog; double prod;
    if (threadIdx.x == 0) OpenLogDevice(slog, prod);
    GENO_READER_CUDA rt(cuda_mem_device.loc_addr[l]);
    GENOTYPE_CUDA* gtab = GetLoc(l).GetGtab();

    for (int i = 0; i < N; ++i, q += K)
    {
        GENOTYPE_CUDA& gt = gtab[rt.Read()];
        int ploidy = 0; ushort* als;
        if (gt.Nalleles() > 0)
        {
            ploidy = gt.Ploidy();
            als = gt.GetAlleleArray();
        }

        REAL* pp;
        if constexpr (!isadmix) pp = p + KT * Z[i];

        for (int a = 0; a < maxploidy; ++a)
        {
            if constexpr (isadmix)
            {
                if constexpr (std::is_same_v<REAL, double> || !fast_fp32)
                    buf2t = a < ploidy ? SumProd(q, p + als[a], KT, K) : 1;
                else
                    buf2t = a < ploidy ? SumProdx(q, p + als[a], KT, K) : 1;
            }
            else
                buf2t = a < ploidy ? pp[als[a]] : 1;

            double v1 = ReduceProd(buf2s);

            if (threadIdx.x == 0)
                ChargeLogDevice(slog, prod, v1);
        }
    }

    if (threadIdx.x == 0)
    {
        CloseLogDevice(slog, prod);
        atomicAdd(&bayes->bufN1_CUDA[l % N], prod);
    }
}

template<typename REAL>
TARGET void BAYESIAN<REAL>::RecordCUDA()
{
    int nthread = 64;
    int block_share_size = sizeof(double) * nthread;
    int nblock = (nloc + nthread - 1) / nthread;

    MemsetCUDA(bufN1_CUDA, 0, N * sizeof(int64));

    if (binaryq)
    {
        MemcpyCUDA(Z_CUDA, Z, N * sizeof(ushort), true);
        g_fastsingle_val == 1 ?
            RecordCUDAKernel<REAL, false, true > <<<nblock, nthread, block_share_size, GPUstream >>> (bayes_CUDA) :
            RecordCUDAKernel<REAL, false, false> <<<nblock, nthread, block_share_size, GPUstream >>> (bayes_CUDA);
    }
    else
    {
        MemcpyCUDA(Q_CUDA, Q, N * K * sizeof(REAL), true);
        g_fastsingle_val == 1 ?
            RecordCUDAKernel<REAL, true, true > <<<nblock, nthread, block_share_size, GPUstream >>> (bayes_CUDA) :
            RecordCUDAKernel<REAL, true, false> <<<nblock, nthread, block_share_size, GPUstream >>> (bayes_CUDA);
    }

    MemcpyCUDA(bufN1, bufN1_CUDA, N * sizeof(double), false);

    WaitCUDA();

    bufNK1[0] = Sum(bufN1, N);
}

#else

TARGETCUDA void WaitCUDA() { }

TARGETCUDA void* MallocHostCUDA(uint64 size) { return Malloc(size); }

TARGETCUDA void FreeHostCUDA(void* addr) { free(addr); }

TARGETCUDA void* MallocDeviceCUDA(uint64 size) { return NULL; }

TARGETCUDA void FreeCUDA(void* addr) { }

TARGETCUDA void MemsetCUDA(void* addr, int val, uint64 size) { }

TARGETCUDA void MemcpyCUDA(void* dst, void* src, uint64 size, bool todevice) { }

TARGETCUDA void FreeStructureMemory(int devID) { }

TARGETCUDA void CopyStructureMemory(int devID) { }

TARGETCUDA int GetDeviceCountCUDA() { return 0; }

TARGETCUDA void AllocMemoryCUDA() { }

TARGETCUDA void FreeMemoryCUDA() { }

TARGETCUDA void CreateStreamCUDA(int devID) { }

TARGETCUDA void DestroyStreamCUDA() { }

TARGETCUDA void ShowDevicesCUDA() { }

TARGETCUDA void ResetDeviceCUDA() { }

template<typename REAL>
TARGET void BAYESIAN<REAL>::UpdateZNoAdmixCUDA() { }

template<typename REAL>
TARGET void BAYESIAN<REAL>::UpdateZAdmixCUDA() { }

template<typename REAL>
TARGET void BAYESIAN<REAL>::UpdateQNoAdmixCUDA() { }

template<typename REAL>
TARGET void BAYESIAN<REAL>::UpdateQMetroCUDA() { }

template<typename REAL>
TARGET void BAYESIAN<REAL>::RecordCUDA() { }

TARGETCUDA void Eig64CUDA(double* A, double* U, double* V, int64 n) { }

TARGETCUDA void Eig32CUDA(float* A, float* U, float* V, int64 n) { }

TARGETCUDA void Svd64CUDA(double* A, double* U, double* S, double* VT, int64 m, int64 n) { }

TARGETCUDA void Svd32CUDA(float* A, float* U, float* S, float* VT, int64 m, int64 n) { }

TARGETCUDA void MatrixMul64CUDA(double* A, double* B, double* res, int64 m, int64 n, int64 p, bool Atrans, bool Btrans) { }

TARGETCUDA void MatrixMul32CUDA(float* A, float* B, float* res, int64 m, int64 n, int64 p, bool Atrans, bool Btrans) { }

#endif
