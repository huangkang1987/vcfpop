/* Matrix Functions */

#include "vcfpop.h"

template TARGET void Evd<double>(Mat<double>& A, Mat<double>& U, Col<double>& V);
template TARGET void Evd<float >(Mat<float >& A, Mat<float >& U, Col<float >& V);

template TARGET void Svd<double>(Mat<double>& A, Mat<double>& U, Col<double>& S, Mat<double>& V);
template TARGET void Svd<float >(Mat<float >& A, Mat<float >& U, Col<float >& S, Mat<float >& V);

template TARGET void SvdCUDA<double>(Mat<double>& A, Mat<double>& U, Col<double>& S, Mat<double>& V);
template TARGET void SvdCUDA<float >(Mat<float >& A, Mat<float >& U, Col<float >& S, Mat<float >& V);

template TARGET void EigCUDA<double>(Mat<double>& A, Mat<double>& U, Col<double>& V);
template TARGET void EigCUDA<float >(Mat<float >& A, Mat<float >& U, Col<float >& V);

template TARGET void MatrixMulCUDA<double>(Mat<double>& A, Mat<double>& B, Mat<double>& res, bool Atrans, bool Btrans);
template TARGET void MatrixMulCUDA<float >(Mat<float >& A, Mat<float >& B, Mat<float >& res, bool Atrans, bool Btrans);

/* econ EVD decomposition, eigen values are in descending order */
template<typename REAL>
TARGET void Evd(rmat& A, rmat& U, rcol& V)
{
	if (g_gpu_val == 1)
	{
		eig_sym(V, U, A, "dc");
		V = arma::flipud(V);
		U = arma::fliplr(U);
	}
	else
		EigCUDA(A, U, V);
}

/* econ SVD decomposition, singular values are in descending order */
template<typename REAL>
TARGET void Svd(rmat& A, rmat& U, rcol& S, rmat& V)
{
	if (g_gpu_val == 1)
		svd_econ(U, S, V, A, "both", "dc");
	else
		SvdCUDA(A, U, S, V);
		
	//svd(U, S, V, Gi, "dc");//10.536
	//svd_econ(U, S, V, Gi, "both", "dc");//6.164
	//SvdCUDA(Gi, U, S, V);		//5.744
	//SvdJacobCUDA(Gi, U, S, V);//7.5
}

/* Use CUDA to perform SVD decomposition */
template<typename REAL>
TARGET void SvdCUDA(rmat& A, rmat& U, rcol& S, rmat& V)
{
	int64 m = A.n_rows, n = A.n_cols, mn = std::min(m, n);

	U = zeros<rmat>(m, mn);
	S = zeros<rcol>(mn);
	V = zeros<rmat>(mn, n);

	if (std::is_same_v<REAL, double>)
	{
		Svd64CUDA((double*)A.memptr(), (double*)U.memptr(), (double*)S.memptr(), (double*)V.memptr(), m, n);
		inplace_trans(V);
	}
	else
	{
		Svd32CUDA((float*)A.memptr(), (float*)U.memptr(), (float*)S.memptr(), (float*)V.memptr(), m, n);
		inplace_trans(V);
	}
}

/* Use CUDA to perform Eigen value decomposition */
template<typename REAL>
TARGET void EigCUDA(rmat& A, rmat& U, rcol& V)
{
	int64 n = A.n_cols;

	A = (A + A.t()) * 0.5;
	rmat& u = U; rcol& s = V;
	u = zeros<rmat>(n, n);
	s = zeros<rcol>(n);
	rmat v = zeros<rmat>(n, n);

	if (std::is_same_v<REAL, double>)
		Svd64CUDA((double*)A.memptr(), (double*)u.memptr(), (double*)s.memptr(), (double*)v.memptr(), n, n);
	else
		Svd32CUDA((float*)A.memptr(), (float*)u.memptr(), (float*)s.memptr(), (float*)v.memptr(), n, n);
}

/* Use CUDA to perform matrix multiplication */
template<typename REAL>
TARGET void MatrixMulCUDA(rmat& A, rmat& B, rmat& res, bool Atrans, bool Btrans)
{
	int64 m  = Atrans ? A.n_cols : A.n_rows, n = Atrans ? A.n_rows : A.n_cols;
	int64 n2 = Btrans ? B.n_cols : B.n_rows, p = Btrans ? B.n_rows : B.n_cols;

	if (n != n2)
		Exit("\nError: MatrixMul dimension mismatch.\n");

	if (res.n_rows != m && res.n_cols != p)
		res = zeros<rmat>(m, p);

	if (std::is_same_v<REAL, double>)
		MatrixMul64CUDA((double*)A.memptr(), (double*)B.memptr(), (double*)res.memptr(), m, n, p, Atrans, Btrans);
	else
		MatrixMul32CUDA((float*)A.memptr(), (float*)B.memptr(), (float*)res.memptr(), m, n, p, Atrans, Btrans);
}
