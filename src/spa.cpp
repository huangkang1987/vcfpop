/* Analysis of Spatial Structure, TEST, DO NOT USE */

#include "vcfpop.h"

template TARGET void CalcSPA<double>();
template TARGET void CalcSPA<float >();
template TARGET double SPA_Fij<double>(double* a, int i);
template TARGET double SPA_Fij<float >(double* a, int i);
template TARGET double SPA_Likelihood<double>(void* Param, CPOINT& xx, Mat<double>& G, Mat<double>& H);
template TARGET double SPA_Likelihood<float >(void* Param, CPOINT& xx, Mat<float >& G, Mat<float >& H);
template TARGET double SPA_Likelihood<double>(uint64 l, double* a);
template TARGET double SPA_Likelihood<float >(uint64 l, double* a);
template TARGET double SPA_Hessian<double>(uint64 l, double* G, double* H, double* a, double* a2);
template TARGET double SPA_Hessian<float >(uint64 l, double* G, double* H, double* a, double* a2);
template TARGET void SPA_Locus<double>(uint64 l, double* x, double* f, double* a, double* a2, double* am, double* g, double* h);
template TARGET void SPA_Locus<float >(uint64 l, double* x, double* f, double* a, double* a2, double* am, double* g, double* h);
template TARGET void MatrixSPA<double>(double* f, double* x, int spa_tn, int spa_np, double* a, int64 l, double* am, double* h, double* g, double* a2);
template TARGET void MatrixSPA<float >(double* f, double* x, int spa_tn, int spa_np, double* a, int64 l, double* am, double* h, double* g, double* a2);

#define extern 
extern thread_local bool* spa_valid;
extern thread_local int spa_tn;
extern double* spa_x;								//Spatial pattern, TEST, don't use, n * (dim + 1)
extern int spa_n;									//Spatial pattern, TEST, don't use, n
extern int spa_np;									//Spatial pattern, TEST, don't use, n
#undef extern 

/* Calculate spatical pattern */
template<typename REAL>
TARGET void CalcSPA()
{
	//locus parameter
	if (!spa) return;

	EvaluationBegin();
	if (spa_level_val == 1 && ad)
		Exit("\nError: SPA analysis at individual level (-spa_level=ind) is incompatible with allelic depth (-ad) option.\n");
	if (spa_level_val == 1 && nind <= spa_dim_val)
		Exit("\nError: the number of individuals should be greater than the dimension of coordinate %d.\n", spa_dim_val);
	if (spa_level_val == 2 && npop <= spa_dim_val)
		Exit("\nError: the number of populations should be greater than the dimension of coordinate %d.\n", spa_dim_val);

	OpenResFile("-spa", "Analysis of spatial structure");
	OpenTempFiles(g_nthread_val, ".spa");

	spa_n = spa_level_val == 1 ? nind : npop;
	spa_np = spa_dim_val + 1;
	spa_x = new double[spa_n * spa_np];
	SetVal(spa_x, -123456789.0, spa_n * spa_np);

	//Load coordinates
	char* t1 = (char*)spa_coord_val.c_str(), * t2 = 0;
	while (*t1)
	{
		while (*t1 == '\n' || *t1 == '\r' || *t1 == ' ' || *t1 == '\t') t1++;
		t2 = StrNextIdx(t1, ":", 1);

		if (spa_level_val == 1)
		{
			bool flag = false;
			for (int i = 0; i < nind; ++i)
			{
				if (memcmp(ainds<REAL>[i]->name, t1, t2 - t1) == 0 && (int)strlen(ainds<REAL>[i]->name) == (int)(t2 - t1))
				{
					if (spa_x[ainds<REAL>[i]->indid * spa_np] != -123456789.0)
						Exit("\nError: the coordinate of individual %s appear twice.\n", t1);
					flag = true;
					t1 = t2 + 1;
					for (int j = 0; j < spa_dim_val; ++j)
					{
						while (*t1 != '-' && *t1 != '.' && !(*t1 >= '0' && *t1 <= '9')) t1++;
						spa_x[ainds<REAL>[i]->indid * spa_np + j] = ReadDouble(t1);
					}
					break;
				}
			}

			if (!flag)
			{
				*t2 = 0;
				Exit("\nError: cannot find individual %s.\n", t1);
			}
		}
		else
		{
			bool flag = false;
			for (int p = 0; p < npop; ++p)
			{
				if (memcmp(apops<REAL>[p]->name, t1, t2 - t1) == 0 && (int)strlen(apops<REAL>[p]->name) == (int)(t2 - t1))
				{
					if (spa_x[apops<REAL>[p]->id * spa_np] != -123456789.0)
						Exit("\nError: the coordinate of population %s appear twice.\n", t1);
					flag = true;
					t1 = t2 + 1;
					for (int j = 0; j < spa_dim_val; ++j)
					{
						while (*t1 != '-' && *t1 != '.' && !(*t1 >= '0' && *t1 <= '9')) t1++;
						spa_x[apops<REAL>[p]->id * spa_np + j] = ReadDouble(t1);
					}
					break;
				}
			}

			if (!flag)
			{
				*t2 = 0;
				Exit("\nError: cannot find population %s.\n", t1);
			}
		}
	}

	for (int i = 0; i < spa_n; ++i)
	{
		if (spa_x[i * spa_np] == -123456789.0)
			Exit("\nError: the coordinate of %s %s is abscent.\n", spa_level_val == 1 ? "individual" : "population", ainds<REAL>[i]->name);
		spa_x[i * spa_np + spa_dim_val] = 1.0;
	}

	RunThreads(&SPAThread<REAL>, NULL, NULL, nloc, nloc, "\nCalculating SPA:\n", g_nthread_val, true);
	JoinTempFiles(g_nthread_val);
	CloseResFile();
	DEL(spa_x);

	EvaluationEnd("SPA");
}

/* Test. spatial pattern, do not use */
template<typename REAL>
TARGET double SPA_Fij(double* a, int i)
{
	double re = 1.0 / (exp(-SumProd(a, spa_x + i * spa_np, spa_np)) + 1.0);
	if (re < spa_truncate_val) re = spa_truncate_val;
	if (re > spa_truncate_val2) re = spa_truncate_val2;
	return re;
}

template<typename REAL>
TARGET double SPA_Likelihood(void* Param, CPOINT& xx, rmat& G, rmat& H)
{
	uint64 l = *(uint64*)Param;
	SetVal(xx.real_space, xx.unc_space, 8);
	return SPA_Likelihood<REAL>(l, xx.unc_space);
}

template<typename REAL>
TARGET double SPA_Likelihood(uint64 l, double* a)
{
	// L = sum_i (ad_1 lnfij + ad_2 ln(1-fij))
	double re = 0;
	double f1 = 0, f2 = 0;
	if (spa_level_val == 1)
	{
		GENOTYPE* gtab = GetLoc(l).GetGtab();
		GENO_READER rt(0, l);

		for (int i = 0; i < nind; ++i)
		{
			GENOTYPE& gt = gtab[rt.Read()];
			if (!spa_valid[i]) continue;
			f1 = SPA_Fij<REAL>(a, i);
			f2 = log(1.0 - f1);
			f1 = log(f1);

			ushort* als = gt.GetAlleleArray();
			for (int j = 0, v = gt.Ploidy(); j < v; ++j)
				re += als[j] == 0 ? f1 : f2;
		}
	}
	else for (int p = 0; p < npop; ++p)
	{
		if (!spa_valid[p]) continue;
		f1 = SPA_Fij<REAL>(a, p);
		f2 = log(1.0 - f1);
		f1 = log(f1);

		int nh = apops<REAL>[p]->loc_stat1[l].nhaplo;
		int n1 = (int)(apops<REAL>[p]->GetFreq(l, 0) * nh + 0.5);
		re += f1 * n1 + f2 * (nh - n1);
	}
	return re;
}

template<typename REAL>
TARGET double SPA_Hessian(uint64 l, double* G, double* H, double* a, double* a2)
{
	double eps = 1e-6;
	double j00 = SPA_Likelihood<REAL>(l, a);
	for (int j = 0; j < spa_np; ++j)
	{
		SetVal(a2, a, spa_np);
		a2[j] += eps;
		double j01 = SPA_Likelihood<REAL>(l, a2);
		double d11 = (j01 - j00) / eps;
		G[j] = d11;
		for (int k = j; k < spa_np; ++k)
		{
			SetVal(a2, a, spa_np);
			a2[k] += eps;
			double j10 = SPA_Likelihood<REAL>(l, a2);

			a2[j] += eps;
			double j11 = SPA_Likelihood<REAL>(l, a2);
			double d12 = (j11 - j10) / eps;
			H[j * spa_np + k] = H[k * spa_np + j] = (d12 - d11) / eps;
		}
	}
	return j00;
}

template<typename REAL>
TARGET void SPA_Locus(uint64 l, double* x, double* f, double* a, double* a2, double* am, double* g, double* h)
{
	char name_buf[NAME_BUF_LEN];
	//At least dim + 1 = 3 individuals should be genotyped

	// L = sum_i (ad_1 lnfij + ad_2 ln(1-fij))
	/*
		x1 a1 + y1 a2 + 1 a3 = -ln(1 / p1 - 1)
		x2 a1 + y2 a2 + 1 a3 = -ln(1 / p2 - 1)
		x3 a1 + y3 a2 + 1 a3 = -ln(1 / p3 - 1)
		x A = f
		A = inv(x) * b
	*/

	//Initial Value
	spa_tn = 0;
	if (spa_level_val == 1)
	{
		GENOTYPE* gtab = GetLoc(l).GetGtab();
		GENO_READER rt(0, l);

		for (int i = 0; i < nind; ++i)
		{
			GENOTYPE& gt = gtab[rt.Read()];
			if (gt.Nalleles() == 0)
			{
				spa_valid[i] = false;
				continue;
			}
			spa_valid[i] = true;
			int ac = 0, v = gt.Ploidy();
			ushort* als = gt.GetAlleleArray();
			for (int j = 0; j < v; ++j)
				if (als[j] == 0)
					ac++;
			f[spa_tn++] = ac / (double)v;
		}
	}
	else for (int p = 0; p < npop; ++p)
	{
		if (apops<REAL>[p]->loc_stat1[l].nhaplo == 0)
		{
			spa_valid[p] = false;
			continue;
		}
		spa_valid[p] = true;
		f[spa_tn++] = apops<REAL>[p]->GetFreq(l, 0);
	}
	if (spa_tn < spa_np) return;

	//rearrange coordinate
	for (int i = 0, c = 0; i < spa_n; ++i)
	{
		if (!spa_valid[i]) continue;
		for (int j = 0; j < spa_np; ++j)
			x[j * spa_tn + c] = spa_x[i * spa_np + j];
		c++;
	}

	//truncate
	for (int i = 0; i < spa_tn; ++i)
	{
		if (f[i] < spa_truncate_val) f[i] = spa_truncate_val;
		if (f[i] > spa_truncate_val2) f[i] = spa_truncate_val2;
		f[i] = -log(1.0 / f[i] - 1.0);
	}

	MatrixSPA<REAL>(f, x, spa_tn, spa_np, a, l, am, h, g, a2);

	double ex2 = 0, ex = 0;
	for (int i = 0; i < spa_n; ++i)
	{
		if (!spa_valid[i]) continue;
		double p = SPA_Fij<REAL>(am, i);
		ex2 += p * p;
		ex += p;
	}
	ex2 /= spa_tn;
	ex /= spa_tn;
	am[spa_np] = sqrt((ex2 - ex * ex) * spa_tn);
	am[spa_np + 1] = MySqrt((ex2 - ex * ex) * spa_tn / (spa_tn - 1));

	FILE* fout = TEMP_FILES[threadid];
	fprintf(fout, "%s", GetLoc(l).GetNameStr(name_buf));
	if (spa_odepth_val == 1)
	{
		if (spa_level_val == 1)
		{
			GENOTYPE* gtab = GetLoc(l).GetGtab();
			GENO_READER rt(0, l);

			for (int i = 0; i < nind; ++i)
			{
				GENOTYPE& gt = gtab[rt.Read()];
				if (!spa_valid[i])
				{
					fprintf(fout, "%c", g_delimiter_val);
					continue;
				}
				ushort* als = gt.GetAlleleArray();
				int v = gt.Ploidy(), n1 = 0;
				for (int j = 0; j < v; ++j)
					if (als[j] == 0) n1++;
				fprintf(fout, "%c", g_delimiter_val);
				fprintf(fout, "%d/%d", n1, v - n1);
			}

		}
		else for (int p = 0; p < npop; ++p)
		{
			if (!spa_valid[p])
			{
				fprintf(fout, "%c", g_delimiter_val);
				continue;
			}
			int nh = apops<REAL>[p]->loc_stat1[l].nhaplo;
			int n1 = (int)(apops<REAL>[p]->GetFreq(l, 0) * nh + 0.5);
			fprintf(fout, "%c", g_delimiter_val);
			fprintf(fout, "%d/%d", n1, nh - n1);
		}
	}
	if (spa_ofreq_val == 1) for (int i = 0; i < spa_n; ++i)
	{
		if (!spa_valid[i])
		{
			fprintf(fout, "%c", g_delimiter_val);
			continue;
		}
		double p = SPA_Fij<REAL>(am, i);
		fprintf(fout, "%c", g_delimiter_val);
		WriteReal(fout, p);
		fprintf(fout, "/");
		WriteReal(fout, 1.0 - p);
	}
	for (int i = 0; i < spa_np + 3; ++i)
	{
		fprintf(fout, "%c", g_delimiter_val);
		WriteReal(fout, am[i]);
	}
	fprintf(fout, "%s", g_linebreak_val);
}

/* SPA function */
template<typename REAL>
TARGET void MatrixSPA(double* f, double* x, int spa_tn, int spa_np, double* a,
	int64 l, double * am, double* h, double* g, double* a2)
{
	Mat<double> B(f, spa_tn, 1, false, true);
	//may be wrong, matrix saved in column format
	Mat<double> X(x, spa_tn, spa_np, false, true);

	Mat<double> A = solve(X, B);
	SetVal(a, A.memptr(), spa_np);

	//Downhill simplex
	int dim = spa_np;
	CPOINT xx0 = CPOINT::DownHillSimplex((void*)l, SPA_Likelihood<REAL>, dim);
	SetVal(am, xx0.unc_space, spa_np);
	am[spa_np + 2] = xx0.lnL;

	//Newton's method
	Mat<double> H(h, spa_np, spa_np, false, true);// symmetric
	Mat<double> G(g, spa_np, 1, false, true);// column vector
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
		Mat<double> D = solve(H, G);
		for (int j = 0; j < spa_np; ++j)
			a[j] -= D(j, 0);
		if (flag) break;
	}
}

/* Calculate spatial pattern using multiple threads */
THREAD2(SPAThread)
{
	FILE* fout = TEMP_FILES[threadid];
	if (threadid == 0)
	{
		fprintf(fout, "%s%sLocus", g_linebreak_val, g_linebreak_val);
		if (spa_odepth_val == 1)
		{
			if (spa_level_val == 1) for (int i = 0; i < nind; ++i)
			{
				fprintf(fout, "%c", g_delimiter_val);
				fprintf(fout, "%s_depth", ainds<REAL>[i]->name);
			}
			else for (int p = 0; p < npop; ++p)
			{
				fprintf(fout, "%c", g_delimiter_val);
				fprintf(fout, "%s_depth", apops<REAL>[p]->name);
			}
		}
		if (spa_ofreq_val == 1)
		{
			if (spa_level_val == 1) for (int i = 0; i < nind; ++i)
			{
				fprintf(fout, "%c", g_delimiter_val);
				fprintf(fout, "%s_freq", ainds<REAL>[i]->name);
			}
			else for (int p = 0; p < npop; ++p)
			{
				fprintf(fout, "%c", g_delimiter_val);
				fprintf(fout, "%s_freq", apops<REAL>[p]->name);
			}
		}
		for (int i = 0; i < spa_dim_val; ++i)
			fprintf(fout, "%ca%d", g_delimiter_val, i + 1);
		fprintf(fout, "%cb%cSPA%cSD%cLikelihood%s", g_delimiter_val, g_delimiter_val, g_delimiter_val, g_delimiter_val, g_linebreak_val);
	}

	int nthread = g_nthread_val;
	double nsec = nloc / (double)nthread + 1e-8;
	uint64 st1 = (uint64)(threadid * nsec), ed1 = (uint64)((threadid + 1) * nsec);

	VLA_NEW(x, double, spa_n * spa_np);
	VLA_NEW(f, double, spa_n);
	VLA_NEW(g, double, spa_np);
	VLA_NEW(h, double, spa_np * spa_np);
	VLA_NEW(a, double, spa_dim_val + 3);
	VLA_NEW(a2, double, spa_np);
	VLA_NEW(am, double, spa_dim_val + 4);
	spa_valid = new bool[spa_n];

	for (uint64 l = st1; l < ed1; ++l)
	{
		SPA_Locus<REAL>(l, x, f, a, a2, am, g, h);

		PROGRESS_VALUE++;
	}

	VLA_DELETE(x);
	VLA_DELETE(f);
	VLA_DELETE(g);
	VLA_DELETE(h);
	VLA_DELETE(a);
	VLA_DELETE(a2);
	VLA_DELETE(am);
	DEL(spa_valid);
}
