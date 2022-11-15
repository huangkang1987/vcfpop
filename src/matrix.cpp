/* Matrix Functions */

#include "vcfpop.h"

template TARGET void MatrixSPA<double>(double* f, double* x, int spa_tn, int spa_np, double* a, int64 l, double* am, double* h, double* g, double* a2);
template TARGET void MatrixSPA<float >(double* f, double* x, int spa_tn, int spa_np, double* a, int64 l, double* am, double* h, double* g, double* a2);

#undef min
#undef max
#undef __AVX512F__
#undef __AVX512BW__
#include "Eigen/Eigen"
#include "Eigen/Core"
#include "Eigen/LU"
#include "Eigen/QR"
#include "Eigen/SVD"
#include "Eigen/Core"
#include "Spectra/SymEigsSolver.h"
#define __AVX512F__
#define __AVX512BW__

using namespace Eigen;
using namespace Spectra;

/* Maximum diagonal element */
TARGET double MaxDiag(double* a, int m, int n)
{
	double re = a[0];
	int e = Min(m, n);
	for (int i = 1; i < e; ++i)
		if (re < a[i * n + i])
			re = a[i * n + i];
	return re;
}

/* Matrix multiplication G = C * D * C */
TARGET void MatrixMul2(double* side, double* mid, int N, double* res)
{
	Map<MatrixXd> C(side, N, N);//symmetric
	Map<MatrixXd> D(mid, N, N);//symmetric
	D = C * D * C;

	for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j)
			res[i * N + j] = D(i, j);
}

/* Matrix multiplication G = C * D * C */
TARGET void MatrixMul2(float* side, float* mid, int N, float* res)
{
	Map<MatrixXf> C(side, N, N);//symmetric
	Map<MatrixXf> D(mid, N, N);//symmetric
	D = C * D * C;

	for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j)
			res[i * N + j] = D(i, j);
}

/* SPA function */
template<typename REAL>
TARGET void MatrixSPA(double* f, double* x, int spa_tn, int spa_np, double* a,
	int64 l, double * am, double* h, double* g, double* a2)
{
	Map<MatrixXd> B(f, spa_tn, 1);
	//may be wrong, matrix saved in column format
	Map<MatrixXd> X(x, spa_tn, spa_np);

	MatrixXd A = X.householderQr().solve(B);
	for (int i = 0; i < spa_np; ++i)
		a[i] = A(i, 0);

	//Downhill simplex
	int dim = spa_np;
	void* Param[] = { (void*)&l };
	CPOINT xx0 = CPOINT::DownHillSimplex(dim, 0, false, 0.0001, 15, SPA_Likelihood<REAL>, Param);
	SetVal(am, xx0.image, spa_np);
	am[spa_np + 2] = xx0.li;

	//Newton's method
	Map<MatrixXd> H(h, spa_np, spa_np);// symmetric
	Map<MatrixXd> G(g, spa_np, 1);// column vector
	SetVal(a2, a, spa_np);
	a2[0] += 1.0;
	double eps = 1e-6, likelihood = 0;
	for (int m = 0; m < 100; ++m)
	{
		likelihood = SPA_Hessian<REAL>(l, g, h, a, a2);
		if (likelihood > am[spa_np + 2])
		{
			SetVal(am, a, spa_np);
			am[spa_np + 2] = likelihood;
		}

		bool flag = true;
		for (int j = 0; j < spa_np; ++j)
			if (abs(a[j] - a2[j]) > eps)
				flag = false;
		SetVal(a2, a, spa_np);
		MatrixXd D = H.lu().solve(G);
		for (int j = 0; j < spa_np; ++j)
			a[j] -= D(j, 0);
		if (flag) break;
	}
}

/* Eigen value decomposition for PCoA */
TARGET void EigenValueDecomp(double *mat, int N, double *&U, double *&V, int& maxp, int *idx)
{
	static bool initflag = false;
	if (!initflag) Eigen::initParallel();
	setNbThreads(g_nthread_val);

	Map<MatrixXd> G(mat, N, N);//symmetric
	DenseSymMatProd<double> op(G);
	SymEigsSolver<DenseSymMatProd<double>> eigs(op, maxp, Min(maxp * 2, N));

	eigs.init();
	eigs.compute(SortRule::LargestAlge);

	// Retrieve results
	auto u = eigs.eigenvectors();
	auto v = eigs.eigenvalues(); 

	//sort by eigen values in descending order
	int nz = (int)v.size();

	for (int i = 0; i < nz; ++i)
		idx[i] = i;

	for (int i = 0; i < nz; ++i)
		for (int j = i + 1; j < nz; ++j)
			if (v(idx[i], 0) < v(idx[j], 0))
				Swap(idx[i], idx[j]);

	// calculate number of axises
	maxp = 0;
	for (int i = 0; i < nz; ++i)
	{
		if (v(idx[i], 0) > 0) maxp++;
		else break;
	}

	U = new double[N * maxp];
	V = new double[N];

	for (int i = 0; i < maxp; ++i)
	{
		V[i] = v(idx[i], 0);
		double s = MySqrt(V[i]);
		for (int j = 0; j < N; ++j)
			U[j * maxp + i] = s * u(j, idx[i]);
	}

	setNbThreads(1);
}

/* Eigen value decomposition for PCoA */
TARGET void EigenValueDecomp(float* mat, int N, float*& U, float*& V, int& maxp, int* idx)
{
	static bool initflag = false;
	if (!initflag) Eigen::initParallel();
	setNbThreads(g_nthread_val);

	Map<MatrixXf> G(mat, N, N);//symmetric
	DenseSymMatProd<float> op(G);
	SymEigsSolver<DenseSymMatProd<float>> eigs(op, maxp, Min(maxp * 2, N));

	eigs.init();
	eigs.compute(SortRule::LargestAlge);

	// Retrieve results
	auto u = eigs.eigenvectors();
	auto v = eigs.eigenvalues();

	//sort by eigen values in descending order
	int nz = (int)v.size();

	for (int i = 0; i < nz; ++i)
		idx[i] = i;

	for (int i = 0; i < nz; ++i)
		for (int j = i + 1; j < nz; ++j)
			if (v(idx[i], 0) < v(idx[j], 0))
				Swap(idx[i], idx[j]);

	// calculate number of axises
	maxp = 0;
	for (int i = 0; i < nz; ++i)
	{
		if (v(idx[i], 0) > 0) maxp++;
		else break;
	}

	U = new float[N * maxp];
	V = new float[N];

	for (int i = 0; i < maxp; ++i)
	{
		V[i] = v(idx[i], 0);
		float s = MySqrt(V[i]);
		for (int j = 0; j < N; ++j)
			U[j * maxp + i] = s * u(j, idx[i]);
	}

	setNbThreads(1);
}

/* Matrix multiplication */
TARGET int MatrixMul(double* l, int lr, int lc, double* r, int rr, int rc, double* res)
{
	if (lc != rr) return -1;
	int i, j, k;
	for (i = 0; i < lr; ++i)
	{
		for (j = 0; j < rc; ++j)
		{
			res[i * rc + j] = 0;
			for (k = 0; k < lc; ++k)
				res[i * rc + j] += l[i * lc + k] * r[k * rc + j];
		}
	}
	return 0;
}

/* Matrix multiplication */
TARGET int MatrixMul(float* l, int lr, int lc, float* r, int rr, int rc, float* res)
{
	if (lc != rr) return -1;
	int i, j, k;
	for (i = 0; i < lr; ++i)
	{
		for (j = 0; j < rc; ++j)
		{
			res[i * rc + j] = 0;
			for (k = 0; k < lc; ++k)
				res[i * rc + j] += l[i * lc + k] * r[k * rc + j];
		}
	}
	return 0;
}

/* Matrix inverstion */
TARGET void MatrixInv2(double* M, int n)
{
	Map<MatrixXd> MC(M, n, n);//symmetric
	MC = MC.inverse();
}

/* Matrix inverstion */
TARGET void MatrixInv2(float* M, int n)
{
	Map<MatrixXf> MC(M, n, n);//symmetric
	MC = MC.inverse();
}

/* Matrix inverstion */
TARGET int MatrixInv(double* M, int n)
{
	VLA_NEW(a, double, n * n);
	int is[100], js[100], i, j, k, l, u, v;
	double d, p;
	for (i = 0; i < n; ++i)
		for (j = 0; j < n; ++j)
			a[i * n + j] = M[i * n + j];

	for (k = 0; k <= n - 1; ++k)
	{
		d = 0.0;
		for (i = k; i <= n - 1; ++i)
			for (j = k; j <= n - 1; ++j)
			{
				l = i * n + j;
				p = fabs(a[l]);
				if (p > d)
				{
					d = p;
					is[k] = i;
					js[k] = j;
				}
			}
		if (d < 2.220446049250313e-016)
		{
			VLA_DELETE(a);
			return -1;
		}
		if (is[k] != k)
			for (j = 0; j <= n - 1; ++j)
			{
				u = k * n + j;
				v = is[k] * n + j;
				p = a[u];
				a[u] = a[v];
				a[v] = p;
			}
		if (js[k] != k)
			for (i = 0; i <= n - 1; ++i)
			{
				u = i * n + k;
				v = i * n + js[k];
				p = a[u];
				a[u] = a[v];
				a[v] = p;
			}
		l = k * n + k;
		a[l] = 1.0 / a[l];
		for (j = 0; j <= n - 1; ++j)
			if (j != k)
			{
				u = k * n + j;
				a[u] = a[u] * a[l];
			}
		for (i = 0; i <= n - 1; ++i)
			if (i != k)
				for (j = 0; j <= n - 1; ++j)
					if (j != k)
					{
						u = i * n + j;
						a[u] = a[u] - a[i * n + k] * a[k * n + j];
					}
		for (i = 0; i <= n - 1; ++i)
			if (i != k)
			{
				u = i * n + k;
				a[u] = -a[u] * a[l];
			}
	}
	for (k = n - 1; k >= 0; k--)
	{
		if (js[k] != k)
			for (j = 0; j <= n - 1; ++j)
			{
				u = k * n + j;
				v = js[k] * n + j;
				p = a[u];
				a[u] = a[v];
				a[v] = p;
			}
		if (is[k] != k)
			for (i = 0; i <= n - 1; ++i)
			{
				u = i * n + k;
				v = i * n + is[k];
				p = a[u];
				a[u] = a[v];
				a[v] = p;
			}
	}
	for (i = 0; i < n; ++i)
		for (j = 0; j < n; ++j)
			M[i * n + j] = a[i * n + j];
	VLA_DELETE(a);
	return 0;
}

/* SVD decomposition */
TARGET void S(double fg[2], double cs[2])
{
	double r, d;
	if ((fabs(fg[0]) + fabs(fg[1])) < EPSILON_SVD)
	{
		cs[0] = 1.0;
		cs[1] = 0.0;
		d = 0.0;
	}
	else
	{
		d = sqrt(fg[0] * fg[0] + fg[1] * fg[1]);
		if (fabs(fg[0]) > fabs(fg[1]))
		{
			d = fabs(d);
			if (fg[0] < 0.0) d = -d;
		}
		if (fabs(fg[1]) >= fabs(fg[0]))
		{
			d = fabs(d);
			if (fg[1] < 0.0) d = -d;
		}
		cs[0] = fg[0] / d;
		cs[1] = fg[1] / d;
	}
	r = 1.0;
	if (fabs(fg[0]) > fabs(fg[1])) r = cs[1];
	else if (fabs(cs[0]) > EPSILON_SVD) r = 1.0 / cs[0];
	fg[0] = d;
	fg[1] = r;
	return;
}

/* SVD decomposition */
TARGET void D(double* a, double* b, int m, int n, int k, double* c)
{
	int i, j, l, u;
	for (i = 0; i <= m - 1; ++i)
		for (j = 0; j <= k - 1; ++j)
		{
			u = i * k + j;
			c[u] = 0;
			for (l = 0; l <= n - 1; ++l)
				c[u] = c[u] + a[i * n + l] * b[l * k + j];
		}
	return;
}

/* SVD decomposition */
TARGET void P(double* a, double* e, double* s, double* v, int m, int n)
{
	int i, j, p, q;
	double d;
	if (m >= n) i = n;
	else i = m;
	for (j = 1; j <= i - 1; ++j)
	{
		a[(j - 1) * n + j - 1] = s[j - 1];
		a[(j - 1) * n + j] = e[j - 1];
	}
	a[(i - 1) * n + i - 1] = s[i - 1];
	if (m < n) a[(i - 1) * n + i] = e[i - 1];
	for (i = 1; i <= n - 1; ++i)
		for (j = i + 1; j <= n; ++j)
		{
			p = (i - 1) * n + j - 1;
			q = (j - 1) * n + i - 1;
			d = v[p];
			v[p] = v[q];
			v[q] = d;
		}
	return;
}

/* SVD decomposition */
TARGET int MatrixSVD(double* a, int m, int n, double* u, double* v, double eps)
{
	int ka = Max(m, n) + 1;
	int i, j, k, l, it, ll, kk, ix, iy, mm, nn, iz, ml, ks;
	double d, dd, t, sm, sml, eml, sk, ek, b, c, shh, fg[2], cs[2];
	VLA_NEW(s, double, ka);
	VLA_NEW(e, double, ka);
	VLA_NEW(w, double, ka);
	for (i = 1; i <= m; ++i)
	{
		ix = (i - 1) * m + i - 1;
		u[ix] = 0;
	}
	for (i = 1; i <= n; ++i)
	{
		iy = (i - 1) * n + i - 1;
		v[iy] = 0;
	}
	it = MAX_ITER_SVD;
	k = n;
	if (m - 1 < n) k = m - 1;
	l = m;
	if (n - 2 < m) l = n - 2;
	if (l < 0) l = 0;
	ll = k;
	if (l > k) ll = l;
	if (ll >= 1)
	{
		for (kk = 1; kk <= ll; ++kk)
		{
			if (kk <= k)
			{
				d = 0.0;
				for (i = kk; i <= m; ++i)
				{
					ix = (i - 1) * n + kk - 1;
					d = d + a[ix] * a[ix];
				}
				s[kk - 1] = sqrt(d);
				if (fabs(s[kk - 1]) > EPSILON_SVD)
				{
					ix = (kk - 1) * n + kk - 1;
					if (fabs(a[ix]) > EPSILON_SVD)
					{
						s[kk - 1] = fabs(s[kk - 1]);
						if (a[ix] < 0.0) s[kk - 1] = -s[kk - 1];
					}
					for (i = kk; i <= m; ++i)
					{
						iy = (i - 1) * n + kk - 1;
						a[iy] = a[iy] / s[kk - 1];
					}
					a[ix] = 1.0 + a[ix];
				}
				s[kk - 1] = -s[kk - 1];
			}
			if (n >= kk + 1)
			{
				for (j = kk + 1; j <= n; ++j)
				{
					if ((kk <= k) && (fabs(s[kk - 1]) > EPSILON_SVD))
					{
						d = 0.0;
						for (i = kk; i <= m; ++i)
						{
							ix = (i - 1) * n + kk - 1;
							iy = (i - 1) * n + j - 1;
							d = d + a[ix] * a[iy];
						}
						d = -d / a[(kk - 1) * n + kk - 1];
						for (i = kk; i <= m; ++i)
						{
							ix = (i - 1) * n + j - 1;
							iy = (i - 1) * n + kk - 1;
							a[ix] = a[ix] + d * a[iy];
						}
					}
					e[j - 1] = a[(kk - 1) * n + j - 1];
				}
			}
			if (kk <= k)
			{
				for (i = kk; i <= m; ++i)
				{
					ix = (i - 1) * m + kk - 1;
					iy = (i - 1) * n + kk - 1;
					u[ix] = a[iy];
				}
			}
			if (kk <= l)
			{
				d = 0.0;
				for (i = kk + 1; i <= n; ++i)
					d = d + e[i - 1] * e[i - 1];
				e[kk - 1] = sqrt(d);
				if (fabs(e[kk - 1]) > EPSILON_SVD)
				{
					if (fabs(e[kk]) > EPSILON_SVD)
					{
						e[kk - 1] = fabs(e[kk - 1]);
						if (e[kk] < 0.0)
							e[kk - 1] = -e[kk - 1];
					}
					for (i = kk + 1; i <= n; ++i)
						e[i - 1] = e[i - 1] / e[kk - 1];
					e[kk] = 1.0 + e[kk];
				}
				e[kk - 1] = -e[kk - 1];
				//if ((kk+1<=m)&&(e[kk-1]!=0.0))
				if ((kk + 1 <= m) && (fabs(e[kk - 1]) > EPSILON_SVD))
				{
					for (i = kk + 1; i <= m; ++i) w[i - 1] = 0.0;
					for (j = kk + 1; j <= n; ++j)
						for (i = kk + 1; i <= m; ++i)
							w[i - 1] = w[i - 1] + e[j - 1] * a[(i - 1) * n + j - 1];
					for (j = kk + 1; j <= n; ++j)
						for (i = kk + 1; i <= m; ++i)
						{
							ix = (i - 1) * n + j - 1;
							a[ix] = a[ix] - w[i - 1] * e[j - 1] / e[kk];
						}
				}
				for (i = kk + 1; i <= n; ++i)
					v[(i - 1) * n + kk - 1] = e[i - 1];
			}
		}
	}
	mm = n;
	if (m + 1 < n) mm = m + 1;
	if (k < n) s[k] = a[k * n + k];
	if (m < mm) s[mm - 1] = 0.0;
	if (l + 1 < mm) e[l] = a[l * n + mm - 1];
	e[mm - 1] = 0.0;
	nn = m;
	if (m > n) nn = n;
	if (nn >= k + 1)
	{
		for (j = k + 1; j <= nn; ++j)
		{
			for (i = 1; i <= m; ++i)
				u[(i - 1) * m + j - 1] = 0.0;
			u[(j - 1) * m + j - 1] = 1.0;
		}
	}
	if (k >= 1)
	{
		for (ll = 1; ll <= k; ++ll)
		{
			kk = k - ll + 1;
			iz = (kk - 1) * m + kk - 1;
			if (fabs(s[kk - 1]) > EPSILON_SVD)
			{
				if (nn >= kk + 1)
					for (j = kk + 1; j <= nn; ++j)
					{
						d = 0.0;
						for (i = kk; i <= m; ++i)
						{
							ix = (i - 1) * m + kk - 1;
							iy = (i - 1) * m + j - 1;
							d = d + u[ix] * u[iy] / u[iz];
						}
						d = -d;
						for (i = kk; i <= m; ++i)
						{
							ix = (i - 1) * m + j - 1;
							iy = (i - 1) * m + kk - 1;
							u[ix] = u[ix] + d * u[iy];
						}
					}
				for (i = kk; i <= m; ++i)
				{
					ix = (i - 1) * m + kk - 1;
					u[ix] = -u[ix];
				}
				u[iz] = 1.0 + u[iz];
				if (kk - 1 >= 1)
					for (i = 1; i <= kk - 1; ++i)
						u[(i - 1) * m + kk - 1] = 0.0;
			}
			else
			{
				for (i = 1; i <= m; ++i)
					u[(i - 1) * m + kk - 1] = 0.0;
				u[(kk - 1) * m + kk - 1] = 1.0;
			}
		}
	}
	for (ll = 1; ll <= n; ++ll)
	{
		kk = n - ll + 1;
		iz = kk * n + kk - 1;
		if ((kk <= l) && (fabs(e[kk - 1]) > EPSILON_SVD))
		{
			for (j = kk + 1; j <= n; ++j)
			{
				d = 0.0;
				for (i = kk + 1; i <= n; ++i)
				{
					ix = (i - 1) * n + kk - 1;
					iy = (i - 1) * n + j - 1;
					d = d + v[ix] * v[iy] / v[iz];
				}
				d = -d;
				for (i = kk + 1; i <= n; ++i)
				{
					ix = (i - 1) * n + j - 1;
					iy = (i - 1) * n + kk - 1;
					v[ix] = v[ix] + d * v[iy];
				}
			}
		}
		for (i = 1; i <= n; ++i)
			v[(i - 1) * n + kk - 1] = 0.0;
		v[iz - n] = 1.0;
	}
	for (i = 1; i <= m; ++i)
		for (j = 1; j <= n; ++j)
			a[(i - 1) * n + j - 1] = 0.0;
	ml = mm;
	it = MAX_ITER_SVD;
	for (;;)
	{
		if (!mm)
		{
			P(a, e, s, v, m, n);
			VLA_DELETE(s);
			VLA_DELETE(e);
			VLA_DELETE(w);
			return l;
		}
		if (!it)
		{
			P(a, e, s, v, m, n);
			VLA_DELETE(s);
			VLA_DELETE(e);
			VLA_DELETE(w);
			return -1;
		}
		kk = mm - 1;
		while ((kk != 0) && (fabs(e[kk - 1]) > EPSILON_SVD))
		{
			d = fabs(s[kk - 1]) + fabs(s[kk]);
			dd = fabs(e[kk - 1]);
			if (dd > eps * d)
				kk = kk - 1;
			else
				e[kk - 1] = 0.0;
		}
		if (kk == mm - 1)
		{
			kk = kk + 1;
			if (s[kk - 1] < 0.0)
			{
				s[kk - 1] = -s[kk - 1];
				for (i = 1; i <= n; ++i)
				{
					ix = (i - 1) * n + kk - 1;
					v[ix] = -v[ix];
				}
			}
			while ((kk != ml) && (s[kk - 1] < s[kk]))
			{
				d = s[kk - 1];
				s[kk - 1] = s[kk];
				s[kk] = d;
				if (kk < n)
					for (i = 1; i <= n; ++i)
					{
						ix = (i - 1) * n + kk - 1;
						iy = (i - 1) * n + kk;
						d = v[ix];
						v[ix] = v[iy];
						v[iy] = d;
					}
				if (kk < m)
					for (i = 1; i <= m; ++i)
					{
						ix = (i - 1) * m + kk - 1;
						iy = (i - 1) * m + kk;
						d = u[ix];
						u[ix] = u[iy];
						u[iy] = d;
					}
				kk = kk + 1;
			}
			it = MAX_ITER_SVD;
			mm = mm - 1;
		}
		else
		{
			ks = mm;
			while ((ks > kk) && (fabs(s[ks - 1]) > EPSILON_SVD))
			{
				d = 0.0;
				if (ks != mm)
					d = d + fabs(e[ks - 1]);
				if (ks != kk + 1) d = d + fabs(e[ks - 2]);
				dd = fabs(s[ks - 1]);
				if (dd > eps * d)
					ks = ks - 1;
				else
					s[ks - 1] = 0.0;
			}
			if (ks == kk)
			{
				kk = kk + 1;
				d = fabs(s[mm - 1]);
				t = fabs(s[mm - 2]);
				if (t > d)
					d = t;
				t = fabs(e[mm - 2]);
				if (t > d)
					d = t;
				t = fabs(s[kk - 1]);
				if (t > d)
					d = t;
				t = fabs(e[kk - 1]);
				if (t > d)
					d = t;
				sm = s[mm - 1] / d;
				sml = s[mm - 2] / d;
				eml = e[mm - 2] / d;
				sk = s[kk - 1] / d;
				ek = e[kk - 1] / d;
				b = ((sml + sm) * (sml - sm) + eml * eml) / 2.0;
				c = sm * eml;
				c = c * c;
				shh = 0.0;
				if ((fabs(b) > EPSILON_SVD) || (fabs(c) > EPSILON_SVD))
				{
					shh = sqrt(b * b + c);
					if (b < 0.0)
						shh = -shh;
					shh = c / (b + shh);
				}
				fg[0] = (sk + sm) * (sk - sm) - shh;
				fg[1] = sk * ek;
				for (i = kk; i <= mm - 1; ++i)
				{
					S(fg, cs);
					if (i != kk)
						e[i - 2] = fg[0];
					fg[0] = cs[0] * s[i - 1] + cs[1] * e[i - 1];
					e[i - 1] = cs[0] * e[i - 1] - cs[1] * s[i - 1];
					fg[1] = cs[1] * s[i];
					s[i] = cs[0] * s[i];
					if ((fabs(cs[0] - 1.0) > EPSILON_SVD) || (fabs(cs[1]) > EPSILON_SVD))
						for (j = 1; j <= n; ++j)
						{
							ix = (j - 1) * n + i - 1;
							iy = (j - 1) * n + i;
							d = cs[0] * v[ix] + cs[1] * v[iy];
							v[iy] = -cs[1] * v[ix] + cs[0] * v[iy];
							v[ix] = d;
						}
					S(fg, cs);
					s[i - 1] = fg[0];
					fg[0] = cs[0] * e[i - 1] + cs[1] * s[i];
					s[i] = -cs[1] * e[i - 1] + cs[0] * s[i];
					fg[1] = cs[1] * e[i];
					e[i] = cs[0] * e[i];
					if (i < m)
						if ((fabs(cs[0] - 1.0) > EPSILON_SVD) || (fabs(cs[1]) > EPSILON_SVD))
							for (j = 1; j <= m; ++j)
							{
								ix = (j - 1) * m + i - 1;
								iy = (j - 1) * m + i;
								d = cs[0] * u[ix] + cs[1] * u[iy];
								u[iy] = -cs[1] * u[ix] + cs[0] * u[iy];
								u[ix] = d;
							}
				}
				e[mm - 2] = fg[0];
				it = it - 1;
			}
			else
			{
				if (ks == mm)
				{
					kk = kk + 1;
					fg[1] = e[mm - 2];
					e[mm - 2] = 0.0;
					for (ll = kk; ll <= mm - 1; ++ll)
					{
						i = mm + kk - ll - 1;
						fg[0] = s[i - 1];
						S(fg, cs);
						s[i - 1] = fg[0];
						if (i != kk)
						{
							fg[1] = -cs[1] * e[i - 2];
							e[i - 2] = cs[0] * e[i - 2];
						}
						if ((fabs(cs[0] - 1.0) > EPSILON_SVD) || (fabs(cs[1]) > EPSILON_SVD))
							for (j = 1; j <= n; ++j)
							{
								ix = (j - 1) * n + i - 1;
								iy = (j - 1) * n + mm - 1;
								d = cs[0] * v[ix] + cs[1] * v[iy];
								v[iy] = -cs[1] * v[ix] + cs[0] * v[iy];
								v[ix] = d;
							}
					}
				}
				else
				{
					kk = ks + 1;
					fg[1] = e[kk - 2];
					e[kk - 2] = 0.0;
					for (i = kk; i <= mm; ++i)
					{
						fg[0] = s[i - 1];
						S(fg, cs);
						s[i - 1] = fg[0];
						fg[1] = -cs[1] * e[i - 1];
						e[i - 1] = cs[0] * e[i - 1];
						if ((fabs(cs[0] - 1.0) > EPSILON_SVD) || (fabs(cs[1]) > EPSILON_SVD))
							for (j = 1; j <= m; ++j)
							{
								ix = (j - 1) * m + i - 1;
								iy = (j - 1) * m + kk - 2;
								d = cs[0] * u[ix] + cs[1] * u[iy];
								u[iy] = -cs[1] * u[ix] + cs[0] * u[iy];
								u[ix] = d;
							}
					}
				}
			}
		}
	}
}

/* L2 norm */
TARGET double MatrixNorm(double* a, int m, int n)
{
	VLA_NEW(b, double, m * n);
	VLA_NEW(u, double, m * m);
	VLA_NEW(v, double, n * n);
	memmove(b, a, m * n * sizeof(double));
	MatrixSVD(b, m, n, u, v);
	double re = MaxDiag(b, m, n);
	VLA_DELETE(b);
	VLA_DELETE(u);
	VLA_DELETE(v);
	return re;
}

/* Condition number */
TARGET double MatrixCond(double* a, int m, int n)
{
	VLA_NEW(b, double, m * n);
	memmove(b, a, m * n * sizeof(double));
	MatrixInv(b, m);
	double re = MatrixNorm(a, m, n) * MatrixNorm(b, m, n);
	VLA_DELETE(b);
	return re;
}

/* Solve Ax = B */
TARGET bool SolveEquation(double* A, double* B, double* x, int n)
{
	/*
	Map<MatrixXd> Am(A, n, n);
	Map<MatrixXd> Bm(B, n, 1);
	//MatrixXd Xm1 = (Am * Am.transpose()).ldlt().solve(Am * Bm);
	//MatrixXd Xm2 = Am.transpose().colPivHouseholderQr().solve(Bm);
	//MatrixXd Xm3 = Am.transpose().bdcSvd(ComputeThinU | ComputeThinV).solve(Bm);

	MatrixXd Xm = Am.transpose().lu().solve(Bm);
	for (int i = 0; i < n; ++i)
		x[i] = Xm4(i, 0);
	*/
	//Solve A x = B
	if (MatrixInv(A, n) == -1)
	{
		Map<MatrixXd> Am(A, n, n);
		Map<MatrixXd> Bm(B, n, 1);
		MatrixXd Xm2 = Am.transpose().colPivHouseholderQr().solve(Bm);
		for (int i = 0; i < n; ++i)
			x[i] = Xm2(i, 0);
	}
	else
		MatrixMul(A, n, n, B, n, 1, x);

	double t = 0;
	switch (n)
	{
	case 1:
		t = x[0];
		break;
	case 2:
		t = x[0] + 0.5 * x[1];
		break;
	case 3:
		t = x[0] + 0.6666666666666666667 * x[1] + 0.3333333333333333333 * x[2];
		break;
	case 4:
		t = x[0] + 0.75 * x[1] + 0.5 * x[2] + 0.25 * x[3];
		break;
	case 5:
		t = x[0] + 0.8 * x[1] + 0.6 * x[2] + 0.4 * x[3] + 0.2 * x[4];
		break;
	case 6:
		t = x[0] + 0.8333333333333333333 * x[1] + 0.6666666666666666667 * x[2] + 0.5 * x[3] + 0.3333333333333333333 * x[4] + 0.1666666666666666667 * x[5];
		break;
	case 7:
		t = x[0] + 0.857142857142857143 * x[1] + 0.714285714285714286 * x[2] + 0.571428571428571429 * x[3] + 0.428571428571428571 * x[4] + 0.285714285714285714 * x[5] + 0.142857142857142857 * x[6];
		break;
	case 8:
		t = x[0] + 0.875 * x[1] + 0.75 * x[2] + 0.625 * x[3] + 0.5 * x[4] + 0.375 * x[5] + 0.25 * x[6] + 0.125 * x[7];
		break;
	}
	if (t > 1.001 || t < -16 || IsError(t))
		return false;
	return true;
}
