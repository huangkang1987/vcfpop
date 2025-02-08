/* Matrix Functions */

#pragma once
#include "vcfpop.h"

/* econ SVD decomposition, eigen values are in descending order */
template<typename REAL>
TARGET void Svd(rmat& A, rmat& U, rcol& S, rmat& V);

/* econ EVD decomposition, singular values are in descending order */
template<typename REAL>
TARGET void Evd(rmat& A, rmat& U, rcol& V);

/* Use CUDA to perform SVD decomposition */
template<typename REAL>
TARGET void SvdCUDA(rmat& A, rmat& U, rcol& S, rmat& V);

/* Use CUDA to perform Eigen value decomposition */
template<typename REAL>
TARGET void EigCUDA(rmat& A, rmat& U, rcol& V);

/* Use CUDA to perform matrix multiplication */
template<typename REAL>
TARGET void MatrixMulCUDA(rmat& A, rmat& B, rmat& res, bool Atrans, bool Btrans);
