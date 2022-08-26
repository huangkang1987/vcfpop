/* Analysis of Spatial Structure, TEST, DO NOT USE */

#pragma once
#include "vcfpop.h"

#define extern 
extern _thread bool* spa_valid;
extern _thread int spa_tn;
extern double* spa_x;								//Spatial pattern, TEST, don't use, n * (dim + 1)
extern int spa_n;									//Spatial pattern, TEST, don't use, n
extern int spa_np;									//Spatial pattern, TEST, don't use, n
#undef extern 

/* Calculate spatical pattern */
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
				if (memcmp(ainds[i]->name, t1, t2 - t1) == 0 && (int)strlen(ainds[i]->name) == (int)(t2 - t1))
				{
					if (spa_x[ainds[i]->indid * spa_np] != -123456789.0)
						Exit("\nError: the coordinate of individual %s appear twice.\n", t1);
					flag = true;
					t1 = t2 + 1;
					for (int j = 0; j < spa_dim_val; ++j)
					{
						while (*t1 != '-' && *t1 != '.' && !(*t1 >= '0' && *t1 <= '9')) t1++;
						spa_x[ainds[i]->indid * spa_np + j] = ReadDouble(t1);
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
				if (memcmp(apops[p]->name, t1, t2 - t1) == 0 && (int)strlen(apops[p]->name) == (int)(t2 - t1))
				{
					if (spa_x[apops[p]->id * spa_np] != -123456789.0)
						Exit("\nError: the coordinate of population %s appear twice.\n", t1);
					flag = true;
					t1 = t2 + 1;
					for (int j = 0; j < spa_dim_val; ++j)
					{
						while (*t1 != '-' && *t1 != '.' && !(*t1 >= '0' && *t1 <= '9')) t1++;
						spa_x[apops[p]->id * spa_np + j] = ReadDouble(t1);
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
			Exit("\nError: the coordinate of %s %s is abscent.\n", spa_level_val == 1 ? "individual" : "population", ainds[i]->name);
		spa_x[i * spa_np + spa_dim_val] = 1.0;
	}

	RunThreads(&SPAThread, NULL, NULL, nloc, nloc, "\nCalculating SPA:\n", g_nthread_val, true);
	JoinTempFiles(g_nthread_val);
	CloseResFile();
	delete[] spa_x;

	EvaluationEnd("SPA");
}

/* Test. spatial pattern, do not use */
TARGET double SPA_Fij(double* a, int i)
{
	double re = 1.0 / (exp(-SumProd(a, spa_x + i * spa_np, spa_np)) + 1.0);
	if (re < spa_truncate_val) re = spa_truncate_val;
	if (re > spa_truncate_val2) re = spa_truncate_val2;
	return re;
}

TARGET double SPA_Likelihood(CPOINT& xx, void** Param)
{
	double* a = xx.image;
	uint64 l = *(uint64*)Param[0];
	return SPA_Likelihood(l, a);
}

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
			f1 = SPA_Fij(a, i);
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
		f1 = SPA_Fij(a, p);
		f2 = log(1.0 - f1);
		f1 = log(f1);

		int nh = apops[p]->loc_stat1[l].nhaplo;
		int n1 = (int)(apops[p]->GetFreq(l, 0) * nh + 0.5);
		re += f1 * n1 + f2 * (nh - n1);
	}
	return re;
}

TARGET double SPA_Hessian(uint64 l, double* G, double* H, double* a, double* a2)
{
	double eps = 1e-6;
	double j00 = SPA_Likelihood(l, a);
	for (int j = 0; j < spa_np; ++j)
	{
		SetVal(a2, a, spa_np);
		a2[j] += eps;
		double j01 = SPA_Likelihood(l, a2);
		double d11 = (j01 - j00) / eps;
		G[j] = d11;
		for (int k = j; k < spa_np; ++k)
		{
			SetVal(a2, a, spa_np);
			a2[k] += eps;
			double j10 = SPA_Likelihood(l, a2);

			a2[j] += eps;
			double j11 = SPA_Likelihood(l, a2);
			double d12 = (j11 - j10) / eps;
			H[j * spa_np + k] = H[k * spa_np + j] = (d12 - d11) / eps;
		}
	}
	return j00;
}

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
		if (apops[p]->loc_stat1[l].nhaplo == 0)
		{
			spa_valid[p] = false;
			continue;
		}
		spa_valid[p] = true;
		f[spa_tn++] = apops[p]->GetFreq(l, 0);
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

	MatrixSPA(f, x, spa_tn, spa_np, a, l, am, h, g, a2);

	double ex2 = 0, ex = 0;
	for (int i = 0; i < spa_n; ++i)
	{
		if (!spa_valid[i]) continue;
		double p = SPA_Fij(am, i);
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
			int nh = apops[p]->loc_stat1[l].nhaplo;
			int n1 = (int)(apops[p]->GetFreq(l, 0) * nh + 0.5);
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
		double p = SPA_Fij(am, i);
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

/* Calculate spatial pattern using multiple threads */
THREAD(SPAThread)
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
				fprintf(fout, "%s_depth", ainds[i]->name);
			}
			else for (int p = 0; p < npop; ++p)
			{
				fprintf(fout, "%c", g_delimiter_val);
				fprintf(fout, "%s_depth", apops[p]->name);
			}
		}
		if (spa_ofreq_val == 1)
		{
			if (spa_level_val == 1) for (int i = 0; i < nind; ++i)
			{
				fprintf(fout, "%c", g_delimiter_val);
				fprintf(fout, "%s_freq", ainds[i]->name);
			}
			else for (int p = 0; p < npop; ++p)
			{
				fprintf(fout, "%c", g_delimiter_val);
				fprintf(fout, "%s_freq", apops[p]->name);
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
		SPA_Locus(l, x, f, a, a2, am, g, h);

		PROGRESS_VALUE++;
	}

	VLA_DELETE(x);
	VLA_DELETE(f);
	VLA_DELETE(g);
	VLA_DELETE(h);
	VLA_DELETE(a);
	VLA_DELETE(a2);
	VLA_DELETE(am);
	delete[] spa_valid;
}