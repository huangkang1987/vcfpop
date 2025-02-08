/* AMX Instruction Set Functions */

#include "vcfpop.h"

#ifdef __aarch64__

#define AMX_NOP_OP_IMM5(op, imm5) \
    __asm("nop\nnop\nnop\n.word (0x201000 + (%0 << 5) + %1)" : : "i"(op), "i"(imm5) : "memory")

#define AMX_OP_GPR(op, gpr) \
    __asm(".word (0x201000 + (%0 << 5) + 0%1 - ((0%1 >> 4) * 6))" : : "i"(op), "r"((uint64_t)(gpr)) : "memory")

#define AMX_LDX(gpr)    AMX_OP_GPR( 0, gpr)
#define AMX_LDY(gpr)    AMX_OP_GPR( 1, gpr)
#define AMX_STX(gpr)    AMX_OP_GPR( 2, gpr)
#define AMX_STY(gpr)    AMX_OP_GPR( 3, gpr)
#define AMX_LDZ(gpr)    AMX_OP_GPR( 4, gpr)
#define AMX_STZ(gpr)    AMX_OP_GPR( 5, gpr)
#define AMX_LDZI(gpr)   AMX_OP_GPR( 6, gpr)
#define AMX_STZI(gpr)   AMX_OP_GPR( 7, gpr)
#define AMX_EXTRX(gpr)  AMX_OP_GPR( 8, gpr)
#define AMX_EXTRY(gpr)  AMX_OP_GPR( 9, gpr)
#define AMX_FMA64(gpr)  AMX_OP_GPR(10, gpr)
#define AMX_FMS64(gpr)  AMX_OP_GPR(11, gpr)
#define AMX_FMA32(gpr)  AMX_OP_GPR(12, gpr)
#define AMX_FMS32(gpr)  AMX_OP_GPR(13, gpr)
#define AMX_MAC16(gpr)  AMX_OP_GPR(14, gpr)
#define AMX_FMA16(gpr)  AMX_OP_GPR(15, gpr)
#define AMX_FMS16(gpr)  AMX_OP_GPR(16, gpr)
#define AMX_SET()       AMX_NOP_OP_IMM5(17, 0)
#define AMX_CLR()       AMX_NOP_OP_IMM5(17, 1)
#define AMX_VECINT(gpr) AMX_OP_GPR(18, gpr)
#define AMX_VECFP(gpr)  AMX_OP_GPR(19, gpr)
#define AMX_MATINT(gpr) AMX_OP_GPR(20, gpr)
#define AMX_MATFP(gpr)  AMX_OP_GPR(21, gpr)
#define AMX_GENLUT(gpr) AMX_OP_GPR(22, gpr)
#define PTR_ROW_FLAGS(ptr, row, flags) (((uint64_t)&*(ptr)) + (((uint64_t)((row) + (flags) * 64)) << 56))


TARGETAMX uint detect_amx_hardware_version() {
    __attribute__((aligned(256))) uint8_t buf[256];
    buf[64] = 3;
    buf[65] = 2;
    buf[129] = 1;
    AMX_SET(); // Set all of x/y/z to zero
    AMX_LDX(PTR_ROW_FLAGS(buf, 48, 1)); // On M1: copy buf[0:128] to x[0,1], on M2: copy buf[0:256] to x[0,1,2,3], on M3/M4: copy buf[0:256] to x[0,2,4,6]
    AMX_VECINT((2ull << 47) + (1ull << 31) + (129ull << 10)); // On M4: z[0] += x[2] + y[0], z[32] += x[3] + y[1]. Before M4: z[0] += concat(x[2][1:64], x[3][0:1]) + y[0], and similar for z[32] on M2/M3.
    AMX_STZ(PTR_ROW_FLAGS(buf, 0, 0)); // Copy z[0] to buf[0:64]
    AMX_CLR();
    return 1 + buf[0];
}

void ShowAMX_FP64(int id)
{
    double d[64][8] = {0};
    
    switch(id)
    {
        case 0: for (uint64 ii = 0; ii <  8; ++ii) AMX_STX((uint64)&d[ii][0] | (((uint64)ii) << 56)); break;
        case 1: for (uint64 ii = 0; ii <  8; ++ii) AMX_STY((uint64)&d[ii][0] | (((uint64)ii) << 56)); break;
        case 2: for (uint64 ii = 0; ii < 64; ++ii) AMX_STZ((uint64)&d[ii][0] | (((uint64)ii) << 56)); break;
    }
    
    switch(id)
    {
        case 0:
            printf("X\n");
            for(int ii = 0; ii < 8; ++ii)
            {
                for (int jj = 0; jj < 8; ++jj)
                    printf("%0.2lf\t", d[ii][jj]);
                printf("\n");
            }
            break;
        case 1:
            printf("Y\n");
            for(int ii = 0; ii < 8; ++ii)
            {
                for (int jj = 0; jj < 8; ++jj)
                    printf("%0.2lf\t", d[ii][jj]);
                printf("\n");
            }
            break;
        case 2:
            printf("Z\n");
            for(int ii = 0; ii < 64; ++ii)
            {
                for (int jj = 0; jj < 8; ++jj)
                    printf("%0.2lf\t", d[ii][jj]);
                printf("\n");
            }
            break;
    }
}

void ShowAMX_FP32(int id)
{
    float d[32][32] = {0};
    
    switch(id)
    {
        case 0: for (uint64 ii = 0; ii <  4; ++ii) AMX_STX((uint64)&d[ii][0] | (((uint64)ii) << 57) | (1ull << 62) ); break;
        case 1: for (uint64 ii = 0; ii <  4; ++ii) AMX_STY((uint64)&d[ii][0] | (((uint64)ii) << 57) | (1ull << 62) ); break;
        case 2: for (uint64 ii = 0; ii < 32; ++ii) AMX_STZ((uint64)&d[ii][0] | (((uint64)ii) << 57) | (1ull << 62) ); break;
    }
    
    switch(id)
    {
        case 0:
            printf("X\n");
            for(int ii = 0; ii < 4; ++ii)
            {
                for (int jj = 0; jj < 32; ++jj)
                    printf("%0.2f\t", d[ii][jj]);
                printf("\n");
            }
            break;
        case 1:
            printf("Y\n");
            for(int ii = 0; ii < 4; ++ii)
            {
                for (int jj = 0; jj < 32; ++jj)
                    printf("%0.2f\t", d[ii][jj]);
                printf("\n");
            }
            break;
        case 2:
            printf("Z\n");
            for(int ii = 0; ii < 32; ++ii)
            {
                for (int jj = 0; jj < 32; ++jj)
                    printf("%0.2f\t", d[ii][jj]);
                printf("\n");
            }
            break;
    }
}
    
void CheckAMX_FP64()
{
    ShowAMX_FP64(0);
    ShowAMX_FP64(1);
    ShowAMX_FP64(2);
}

void CheckAMX_FP32()
{
    ShowAMX_FP32(0);
    ShowAMX_FP32(1);
    ShowAMX_FP32(2);
}

TARGETAMX void DiagQuadFormAMX(double* res, double* A, double* D, int64 m, int64 n)
{
#define NA 16
#define NB 16
#define LDA     AMX_LDX( (uint64)Ax | (0ull << 56) | (1ull << 62) );
    
#define LDB     Mul(Bt, Bx, D[k], NB); AMX_LDY( (uint64)Bt | (0ull << 56) | (1ull << 62) );
    
#define FMA     UNROLL(4) AMX_FMA64( (kk << 20) | ((kk & 1) << 16) | ((kk >> 1) << 6) ); \

#define STZ     UNROLL(8) AMX_STZ( (uint64)&Z[kk + 0][0] | ((0ull + kk * 8) << 56) | (1ull << 62) ); \
                UNROLL(8) AMX_STZ( (uint64)&Z[kk + 8][0] | ((2ull + kk * 8) << 56) | (1ull << 62) ); 

#define CPZ     for (uint64 jj = 0; jj < Nb; ++jj) \
                    for (uint64 ii = 0; ii < Na; ++ii) \
                        res[(i + ii) * m + j + jj] = res[(j + jj) * m + i + ii] = Z[jj][ii];
    double Bt[NB];
    double Z[NB][NA] = { 0 };
    for (uint64 i = 0; i < m; i += NA) {
        uint64 Na = std::min((uint64)NA, m - i);
        for (uint64 j = 0; j <= i; j += NB) {
            uint64 Nb = std::min((uint64)NB, m - j);
            double* Ax = A + i, *Bx = A + j; AMX_SET();
            for (uint64 k = 0; k < n; k++, Ax += m, Bx += m)
            { LDA;  LDB;  FMA; }
            STZ; CPZ; AMX_CLR();
        }
    }
}

TARGETAMX void DiagQuadFormAMX(float* res, float* A, float* D, int64 m, int64 n)
{
#define NA 32
#define NB 32
    
#define LDA     AMX_LDX( (uint64)Ax | (0ull << 56) | (1ull << 62) );
    
#define LDB     Mul(Bt, Bx, D[k], NB); AMX_LDY( (uint64)Bt | (0ull << 56) | (1ull << 62) );
    
#define FMA     UNROLL(4) AMX_FMA32( (kk << 20) | ((kk & 1) << 16) | ((kk >> 1) << 6) ); 

#define STZ     UNROLL(16) AMX_STZ( (uint64)&Z[kk +  0][0] | ((0ull + kk * 4) << 56) | (1ull << 62) ); \
                UNROLL(16) AMX_STZ( (uint64)&Z[kk + 16][0] | ((2ull + kk * 4) << 56) | (1ull << 62) );
    
#define CPZ     for (uint64 jj = 0; jj < Nb; ++jj) \
                    for (uint64 ii = 0; ii < Na; ++ii) \
                        res[(i + ii) * m + j + jj] = res[(j + jj) * m + i + ii] = Z[jj][ii];
    float Bt[NB];
    float Z[NB][NA] = { 0 };
    for (uint64 i = 0; i < m; i += NA) {
        uint64 Na = std::min((uint64)NA, m - i);
        for (uint64 j = 0; j <= i; j += NB) {
            uint64 Nb = std::min((uint64)NB, m - j);
            float* Ax = A + i, *Bx = A + j; AMX_SET();
            for (uint64 k = 0; k < n; k++, Ax += m, Bx += m)
            { LDA;  LDB;  FMA; }
            STZ; CPZ; AMX_CLR();
        }
    }
}

TARGETAMX void MatrixMulAMX(double* res, double* A, double* B, int64 m, int64 n, int64 p)
{
#define NA 32
#define NB 16
#define LDA             AMX_LDX( (uint64)&Ax[ 0] | (0ull << 56) | (1ull << 62) ); \
                        AMX_LDX( (uint64)&Ax[16] | (2ull << 56) | (1ull << 62) );
    
#define LDB             AMX_LDY( (uint64)&Bx[ 0] | (0ull << 56) | (1ull << 62) );
    
#define FMA             UNROLL(8) AMX_FMA64( (kk << 20) | ((kk & 3) << 16) | ((kk >> 2) << 6) );

#define STZ             UNROLL(8) { AMX_STZ( (uint64)&Z[kk + 0][ 0] | ((4*kk+0ull) << 57) | (1ull << 62) );   \
                                    AMX_STZ( (uint64)&Z[kk + 0][16] | ((4*kk+1ull) << 57) | (1ull << 62) ); } \
                        UNROLL(8) { AMX_STZ( (uint64)&Z[kk + 8][ 0] | ((4*kk+2ull) << 57) | (1ull << 62) );   \
                                    AMX_STZ( (uint64)&Z[kk + 8][16] | ((4*kk+3ull) << 57) | (1ull << 62) ); }

#define CPZ             for (uint64 jj = 0; jj < Nb; ++jj) \
                            for (uint64 ii = 0; ii < Na; ++ii) \
                                res[(j + jj) * m + i + ii] = Z[jj][ii];

    int64 m_stride = (m + NA - 1) / NA * NA;
    int64 p_stride = (p + NB - 1) / NB * NB;

    double Z[NB][NA] = { 0 };
    for (uint64 i = 0; i < m; i += NA) {
        uint64 Na = std::min((uint64)NA, m - i);
        for (uint64 j = 0; j < p; j += NB) {
            uint64 Nb = std::min((uint64)NB, p - j);
            double* Ax = A + i, *Bx = B + j;
            AMX_SET();
            for (uint64 k = 0; k < n; k ++, Ax += m_stride, Bx += p_stride)
            { LDA;  LDB;  FMA; }
            STZ; CPZ; AMX_CLR();
        }
    }
}

TARGETAMX void MatrixMulAMX(float* res, float* A, float* B, int64 m, int64 n, int64 p)
{
#define NA 32
#define NB 32
#define LDA     AMX_LDX( (uint64)&Ax[ 0] | (0ull << 56) | (1ull << 62) );
    
#define LDB     AMX_LDY( (uint64)&Bx[ 0] | (0ull << 56) | (1ull << 62) );
    
#define FMA     AMX_FMA32( (0ull << 20) | (0ull << 16) | (0ull << 6) ); \
                AMX_FMA32( (1ull << 20) | (1ull << 16) | (0ull << 6) ); \
                AMX_FMA32( (2ull << 20) | (0ull << 16) | (1ull << 6) ); \
                AMX_FMA32( (3ull << 20) | (1ull << 16) | (1ull << 6) );
    
#define STZ    switch (Nb) { \
                case 32: AMX_STZ( (uint64)&res[i + (j + 31) * m] | (31ull << 57) | (1ull << 62) ); \
                case 31: AMX_STZ( (uint64)&res[i + (j + 30) * m] | (29ull << 57) | (1ull << 62) ); \
                case 30: AMX_STZ( (uint64)&res[i + (j + 29) * m] | (27ull << 57) | (1ull << 62) ); \
                case 29: AMX_STZ( (uint64)&res[i + (j + 28) * m] | (25ull << 57) | (1ull << 62) ); \
                case 28: AMX_STZ( (uint64)&res[i + (j + 27) * m] | (23ull << 57) | (1ull << 62) ); \
                case 27: AMX_STZ( (uint64)&res[i + (j + 26) * m] | (21ull << 57) | (1ull << 62) ); \
                case 26: AMX_STZ( (uint64)&res[i + (j + 25) * m] | (19ull << 57) | (1ull << 62) ); \
                case 25: AMX_STZ( (uint64)&res[i + (j + 24) * m] | (17ull << 57) | (1ull << 62) ); \
                case 24: AMX_STZ( (uint64)&res[i + (j + 23) * m] | (15ull << 57) | (1ull << 62) ); \
                case 23: AMX_STZ( (uint64)&res[i + (j + 22) * m] | (13ull << 57) | (1ull << 62) ); \
                case 22: AMX_STZ( (uint64)&res[i + (j + 21) * m] | (11ull << 57) | (1ull << 62) ); \
                case 21: AMX_STZ( (uint64)&res[i + (j + 20) * m] | ( 9ull << 57) | (1ull << 62) ); \
                case 20: AMX_STZ( (uint64)&res[i + (j + 19) * m] | ( 7ull << 57) | (1ull << 62) ); \
                case 19: AMX_STZ( (uint64)&res[i + (j + 18) * m] | ( 5ull << 57) | (1ull << 62) ); \
                case 18: AMX_STZ( (uint64)&res[i + (j + 17) * m] | ( 3ull << 57) | (1ull << 62) ); \
                case 17: AMX_STZ( (uint64)&res[i + (j + 16) * m] | ( 1ull << 57) | (1ull << 62) ); \
                case 16: AMX_STZ( (uint64)&res[i + (j + 15) * m] | (30ull << 57) | (1ull << 62) ); \
                case 15: AMX_STZ( (uint64)&res[i + (j + 14) * m] | (28ull << 57) | (1ull << 62) ); \
                case 14: AMX_STZ( (uint64)&res[i + (j + 13) * m] | (26ull << 57) | (1ull << 62) ); \
                case 13: AMX_STZ( (uint64)&res[i + (j + 12) * m] | (24ull << 57) | (1ull << 62) ); \
                case 12: AMX_STZ( (uint64)&res[i + (j + 11) * m] | (22ull << 57) | (1ull << 62) ); \
                case 11: AMX_STZ( (uint64)&res[i + (j + 10) * m] | (20ull << 57) | (1ull << 62) ); \
                case 10: AMX_STZ( (uint64)&res[i + (j +  9) * m] | (18ull << 57) | (1ull << 62) ); \
                case  9: AMX_STZ( (uint64)&res[i + (j +  8) * m] | (16ull << 57) | (1ull << 62) ); \
                case  8: AMX_STZ( (uint64)&res[i + (j +  7) * m] | (14ull << 57) | (1ull << 62) ); \
                case  7: AMX_STZ( (uint64)&res[i + (j +  6) * m] | (12ull << 57) | (1ull << 62) ); \
                case  6: AMX_STZ( (uint64)&res[i + (j +  5) * m] | (10ull << 57) | (1ull << 62) ); \
                case  5: AMX_STZ( (uint64)&res[i + (j +  4) * m] | ( 8ull << 57) | (1ull << 62) ); \
                case  4: AMX_STZ( (uint64)&res[i + (j +  3) * m] | ( 6ull << 57) | (1ull << 62) ); \
                case  3: AMX_STZ( (uint64)&res[i + (j +  2) * m] | ( 4ull << 57) | (1ull << 62) ); \
                case  2: AMX_STZ( (uint64)&res[i + (j +  1) * m] | ( 2ull << 57) | (1ull << 62) ); \
                case  1: AMX_STZ( (uint64)&res[i + (j +  0) * m] | ( 0ull << 57) | (1ull << 62) ); }
    
#define STZ2    UNROLL(16) AMX_STZ( (uint64)&Z[kk +  0][ 0] | ((kk*4 + 0ull) << 56) | (1ull << 62) ); \
                UNROLL(16) AMX_STZ( (uint64)&Z[kk + 16][ 0] | ((kk*4 + 2ull) << 56) | (1ull << 62) ); \
                for (uint64 jj = 0; jj < Nb; ++jj) SetVal(&res[i + (j + jj) * m], &Z[jj][0], Na);
    
int64 m_stride = (m + NA - 1) / NA * NA;
int64 p_stride = (p + NB - 1) / NB * NB;
    float Z[NB][NA] = { 0 };
    for (uint64 i = 0; i < m; i += NA)
    {  
        uint64 Na = std::min((uint64)NA, m - i);
        for (uint64 j = 0; j < p; j += NB)
        {
            uint64 Nb = std::min((uint64)NB, p - j);
            float* Ax = A + i, *Bx = B + j;  AMX_SET();
            for (uint64 k = 0; k < n; k++, Ax += m_stride, Bx += p_stride)
            { LDA;  LDB;  FMA; }
            STZ2; AMX_CLR();
        }
    }
}

#else

/* Get AMX version */
TARGETAMX uint detect_amx_hardware_version()
{
	Exit("Error: enter AMX functions.\n");
    return 0;
}

/* Quadratic form A D A' with D being a diagonal matrix, A is m*n, D is n*n, RowMajor with a 16 stride */
TARGETAMX void DiagQuadFormAMX(double* res, double* A, double* D, int64 m, int64 n)
{
	Exit("Error: enter AMX functions.\n");
}

/* Quadratic form A D A' with D being a diagonal matrix, A is m*n, D is n*n, RowMajor with a 16 stride */
TARGETAMX void DiagQuadFormAMX(float* res, float* A, float* D, int64 m, int64 n)
{
	Exit("Error: enter AMX functions.\n");
}

/* Matrix Muplification for A'B, A is n*m, B is n*p, RowMajor with 32 and 16 strides */
TARGETAMX void MatrixMulAMX(double* res, double* A, double* B, int64 m, int64 n, int64 p)
{
	Exit("Error: enter AMX functions.\n");
}

/* Matrix Muplification for A'B, A is n*m, B is n*p, RowMajor with a 32 stride */
TARGETAMX void MatrixMulAMX(float* res, float* A, float* B, int64 m, int64 n, int64 p)
{
	Exit("Error: enter AMX functions.\n");
}

#endif
