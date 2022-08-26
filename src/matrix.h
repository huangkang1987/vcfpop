/* Matrix Functions */

#pragma once
#include "vcfpop.h"

/* Maximum diagonal element */
TARGET double MaxDiag(double* a, int m, int n);

/* Matrix multiplication G = C * D * C */
TARGET void MatrixMul2(double* side, double* mid, int N, double* res);

/* Eigen value decomposition for PCoA */
TARGET void EigenValueDecomp(double* mat, int N, double*& U, double*& V, int& maxp, int* idx);

TARGET void MatrixSPA(double* f, double* x, int spa_tn, int spa_np, double* a,
	int64 l, double* am, double* h, double* g, double* a2);

/* Matrix multiplication */
TARGET int MatrixMul(double* l, int lr, int lc, double* r, int rr, int rc, double* res);

/* Matrix inverstion */
TARGET void MatrixInv2(double* M, int m);

/* Matrix inverstion */
TARGET int MatrixInv(double* M, int m);

/* SVD decomposition */
TARGET void S(double fg[2], double cs[2]);

/* SVD decomposition */
TARGET void D(double* a, double* b, int m, int n, int k, double* c);

/* SVD decomposition */
TARGET void P(double* a, double* e, double* s, double* v, int m, int n);

/* SVD decomposition */
TARGET int MatrixSVD(double* a, int m, int n, double* u, double* v, double eps = EPSILON_SVD);

/* L2 norm */
TARGET double MatrixNorm(double* a, int m, int n);

/* Condition number */
TARGET double MatrixCond(double* a, int m, int n);

/* Solve Ax = B */
TARGET bool SolveEquation(double* A, double* B, double* x, int n);