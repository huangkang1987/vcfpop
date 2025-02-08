/* AMX Instruction Set Functions */

#pragma once
#include "vcfpop.h"

/* Get AMX version */
TARGETAMX uint detect_amx_hardware_version();

/* Quadratic form A D A' with D being a diagonal matrix, A is m*n, D is n*n, RowMajor with a 16 stride */
TARGETAMX void DiagQuadFormAMX(double* res, double* A, double* D, int64 m, int64 n);

/* Quadratic form A D A' with D being a diagonal matrix, A is m*n, D is n*n, RowMajor with a 16 stride */
TARGETAMX void DiagQuadFormAMX(float* res, float* A, float* D, int64 m, int64 n);

/* Matrix Muplification for A'B, A is n*m, B is n*p, RowMajor with 32 and 16 strides */
TARGETAMX void MatrixMulAMX(double* res, double* A, double* B, int64 m, int64 n, int64 p);

/* Matrix Muplification for A'B, A is n*m, B is n*p, RowMajor with a 32 stride */
TARGETAMX void MatrixMulAMX(float* res, float* A, float* B, int64 m, int64 n, int64 p);
