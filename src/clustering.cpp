/* Hierarchical Clustering Functions */

#pragma once
#include "vcfpop.h"

#ifndef _HCLUSTERING
	/* Initialize for distance matrix between individuals */
	TARGET HCLUSTERING::HCLUSTERING(double* d, IND** obj, int n, int m, MEMORY* _memory)
	{
		memory = _memory;
		method = m;
		nori = n;
		ncur = n;
		dori = d;
		memory->Alloc(dcur, n * n);
		memory->Alloc(dnew, n * n);
		SetVal(dcur, d, n * n);
		SetVal(dcur, 1e300, n, n + 1);

		// add some value to avoid float point error or same values of distance
		RNG rng(1);
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < i; ++j)
				dcur[i * n + j] = dcur[j * n + i] = 
					dcur[i * n + j] + dcur[i * n + j] * rng.Uniform() * 1e-6;

		node.SetSize(n);
		for (int i = 0; i < n; ++i)
		{
			HCLUSTER* tc;
			memory->Alloc(tc, 1);
			tc->isend = true;
			tc->endname = obj[i]->name;
			tc->x = i;
			tc->y = 0;
			tc->idlen = 1;
			memory->Alloc(tc->id, 1);
			tc->id[0] = (ushort)i;
			tc->left = NULL;
			tc->right = NULL;
			node.Push(tc);
		}
		Cluster();
	}

	/* Initialize for distance matrix between populations or regions */
	TARGET HCLUSTERING::HCLUSTERING(double* d, POP** obj, int n, int m, MEMORY* _memory)
	{
		memory = _memory;
		method = m;
		nori = n;
		ncur = n;
		dori = d;
		memory->Alloc(dcur, n * n);
		memory->Alloc(dnew, n * n);
		SetVal(dcur, d, n * n);
		SetVal(dcur, 1e300, n, n + 1);

		// add some value to avoid float point error or same values of distance
		RNG rng(1);
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < i; ++j)
				dcur[i * n + j] = dcur[j * n + i] =
				dcur[i * n + j] + dcur[i * n + j] * rng.Uniform() * 1e-6;

		node.SetSize(n);
		for (int i = 0; i < n; ++i)
		{
			HCLUSTER* tc;
			memory->Alloc(tc, 1);
			tc->isend = true;
			tc->endname = obj[i]->name;
			tc->x = i;
			tc->y = 0;
			tc->idlen = 1;
			memory->Alloc(tc->id, 1);
			tc->id[0] = (ushort)i;
			tc->left = NULL;
			tc->right = NULL;
			node.Push(tc);
		}
		Cluster();
	}

	/* Uninitialize */
	TARGET HCLUSTERING::~HCLUSTERING()
	{

	}

	/* Perform clustering */
	TARGET void HCLUSTERING::Cluster()
	{
		while (ncur > 1)
		{
			HCLUSTER* tc;
			memory->Alloc(tc, 1);
			int a, b;
			tc->isend = false;
			tc->endname = NULL;
			double mind = FindMinIdx(a, b);
			tc->y = mind * 0.5;
			tc->x = (node[a]->x + node[b]->x) * 0.5;
			tc->idlen = node[a]->idlen + node[b]->idlen;
			memory->Alloc(tc->id, tc->idlen);
			SetVal(tc->id, node[a]->id, node[a]->idlen);
			SetVal(tc->id + node[a]->idlen, node[b]->id, node[b]->idlen);
			tc->left = node[a];
			tc->right = node[b];

			ReduceMatrix(a, b);

			Swap(dcur, dnew);
			ncur--;

			node[a] = tc;
			node.Erase(b);
		}
	}

	/* Find index for minimum distance in the distance matrix */
	TARGET double HCLUSTERING::FindMinIdx(int& a, int& b)
	{
		double minv = 0;
		int id = GetMinIdx(dcur, ncur * ncur, minv);
		int a1 = id / ncur;
		int b1 = id % ncur;
		a = Min(a1, b1);
		b = Max(a1, b1);
		return minv;
	}

	/* Reduce distance matrix from nxn to (n-1)x(n-1) */
	TARGET void HCLUSTERING::ReduceMatrix(int _a, int _b)
	{
		int nnew = ncur - 1;
		int a = Min(_a, _b), b = Max(_a, _b);

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int i = 0; i < ncur; ++i)
		{
			if (i == b) continue; //bug fixed on 20220816
			int i2 = i > b ? i - 1 : i;
			SetVal(dnew + i2 * nnew, dcur + i * ncur, b);
			SetVal(dnew + i2 * nnew + b, dcur + i * ncur + b + 1, (ncur - b - 1));
		}

		if (method == 1)
		{
			//NEAREST
			for (int c = 0; c < ncur; ++c)
			{
				if (c == a || c == b) continue;
				int cnew = c > b ? c - 1 : c;
				double dac = dcur[a * ncur + c], dbc = dcur[b * ncur + c];
				dnew[cnew * nnew + a] = dnew[a * nnew + cnew] = Min(dac, dbc);
			}
		}

		if (method == 2)
		{
			//FURTHEST
			for (int c = 0; c < ncur; ++c)
			{
				if (c == a || c == b) continue;
				int cnew = c > b ? c - 1 : c;
				double dac = dcur[a * ncur + c], dbc = dcur[b * ncur + c];
				dnew[cnew * nnew + a] = dnew[a * nnew + cnew] = Max(dac, dbc);
			}
		}

		if (method == 3)
		{
			//UPGMA
			for (int c = 0; c < ncur; ++c)
			{
				if (c == a || c == b) continue;
				int cnew = c > b ? c - 1 : c;
				double dac = dcur[a * ncur + c], dbc = dcur[b * ncur + c];
				int na = node[a]->idlen, nb = node[b]->idlen;
				dnew[cnew * nnew + a] = dnew[a * nnew + cnew] = (dac * na + dbc * nb) / (na + nb);
			}
		}

		if (method == 4)
		{
			//WPGMA
			for (int c = 0; c < ncur; ++c)
			{
				if (c == a || c == b) continue;
				int cnew = c > b ? c - 1 : c;
				double dac = dcur[a * ncur + c], dbc = dcur[b * ncur + c];
				dnew[cnew * nnew + a] = dnew[a * nnew + cnew] = (dac + dbc) * 0.5;
			}
		}

		if (method == 5)
		{
			//UPGMC
			for (int c = 0; c < ncur; ++c)
			{
				if (c == a || c == b) continue;
				int cnew = c > b ? c - 1 : c;
				double dac = dcur[a * ncur + c], dbc = dcur[b * ncur + c], dab = dcur[a * ncur + b];
				int na = node[a]->idlen, nb = node[b]->idlen;
				dnew[cnew * nnew + a] = dnew[a * nnew + cnew] = MySqrt((na * dac * dac + nb * dbc * dbc) / (na + nb) - (dab * dab * na * nb) / (na + nb) / (na + nb));
			}
		}

		if (method == 6)
		{
			//WPGMC
			for (int c = 0; c < ncur; ++c)
			{
				if (c == a || c == b) continue;
				int cnew = c > b ? c - 1 : c;
				double dac = dcur[a * ncur + c], dbc = dcur[b * ncur + c], dab = dcur[a * ncur + b];
				dnew[cnew * nnew + a] = dnew[a * nnew + cnew] = MySqrt(dac * dac * 0.5 + dbc * dbc * 0.5 - dab * dab * 0.25);
			}
		}

		if (method == 7)
		{
			//WARD
			for (int c = 0; c < ncur; ++c)
			{
				if (c == a || c == b) continue;
				int cnew = c > b ? c - 1 : c;
				double dac = dcur[a * ncur + c], dbc = dcur[b * ncur + c], dab = dcur[a * ncur + b];
				int na = node[a]->idlen, nb = node[b]->idlen, nc = node[c]->idlen;
				dnew[cnew * nnew + a] = dnew[a * nnew + cnew] = MySqrt(((na + nc) * dac * dac + (nb + nc) * dbc * dbc - nc * dab * dab) / (na + nb + nc));
			}
		}
	}

	/* Print clustering results */
	TARGET void HCLUSTERING::PrintClustering(FILE* fout, HCLUSTER* c, double cy)
	{
		if (c == NULL)
		{
			c = node[0];
			if (c->isend)
			{
				fprintf(fout, "%s:", c->endname);
				WriteReal(fout, cy - c->y);
			}
			else
			{
				fprintf(fout, "(");
				PrintClustering(fout, c->left, c->y);
				fprintf(fout, ",");
				PrintClustering(fout, c->right, c->y);
				fprintf(fout, ");%s", g_linebreak_val);
			}
		}
		else if (c->isend)
		{
			fprintf(fout, "%s:", c->endname);
			WriteReal(fout, cy - c->y);
		}
		else
		{
			fprintf(fout, "(");
			PrintClustering(fout, c->left, c->y);
			fprintf(fout, ",");
			PrintClustering(fout, c->right, c->y);
			fprintf(fout, "):");
			WriteReal(fout, cy - c->y);
		}
	}
#endif

#define extern 
extern MEMORY* clustering_memory;					//Genetic distance memory class
extern double* clustering_matrix;					//Genetic distance array of genetic distance to perform PCoA
#undef extern 

/* Calculate hierarchical clustering */
TARGET void CalcClustering()
{
	if (!cluster) return;

	EvaluationBegin();
	GDIST_METHOD = 3;
	OpenResFile("-cluster", "#Hierarchical clustering");

	bool isfirst = true;
	int ntot = 0;
	if (cluster_level_val[1]) ntot += nind * (nind - 1) / 2;
	if (cluster_level_val[2]) ntot += npop * (npop - 1) / 2;
	if (cluster_level_val[3])
		for (int rl = 0; rl < lreg; ++rl)
			ntot += nreg[rl] * (nreg[rl] - 1) / 2;

	const char* level[] = { "", "individual", "population", "region level" };

	for (int m = 1; m <= 3; ++m)
	{
		if (cluster_level_val[m] == 0) continue;
		for (int rl = 0; rl < (m == 3 ? lreg : 1); ++rl)
		{
			gdindex[0] = m + rl;
			int n = m == 1 ? nind : (m == 2 ? npop : nreg[rl]);

			//Calculate genetic distance table between any two genotypes
			clustering_memory = new MEMORY[g_nthread_val];
			if (m == 1) GDIST::CacheIndGD();

			int nestimator = 0;
			for (int k = 1; k <= (m == 1 ? N_GD_ESTIMATOR - 2 * N_FST_ESTIMATOR : N_GD_ESTIMATOR); ++k)
				if (cluster_estimator_val[k])
					gdindex[k] = nestimator++;

			clustering_matrix = new double[n * n * nestimator];
			SetZero(clustering_matrix, n * n * nestimator);

			RunThreads(&ClusteringThread, NULL, NULL, ntot, n * (n - 1) / 2,
				"\nPerforming hierarchical clustering:\n", g_nthread_val, isfirst);
			isfirst = false;

			for (int k = 1; k <= (m == 1 ? N_GD_ESTIMATOR - 2 * N_FST_ESTIMATOR : N_GD_ESTIMATOR); ++k)
			{
				if (cluster_estimator_val[k] == 0) continue;
				for (int64 l = 1; l <= N_CLUSTER_METHOD; ++l)
				{
					if (cluster_method_val[l] == 0) continue;
					switch (m)
					{
					case 1:
					{
						fprintf(FRES, "#Level: %s, genetic distance: %s, method: %s%s", level[m], GD_ESTIMATOR[k], CLUSTER_METHOD[l], g_linebreak_val);
						HCLUSTERING tcluster(clustering_matrix + gdindex[k] * n * n, ainds, n, l, clustering_memory);
						tcluster.PrintClustering(FRES);
						break;
					}
					case 2:
					{
						fprintf(FRES, "#Level: %s, genetic distance: %s, method: %s%s", level[m], GD_ESTIMATOR[k], CLUSTER_METHOD[l], g_linebreak_val);
						HCLUSTERING tcluster(clustering_matrix + gdindex[k] * n * n, apops, n, l, clustering_memory);
						tcluster.PrintClustering(FRES);
						break;
					}
					case 3:
					{
						for (int rl2 = 0; rl2 < lreg; ++rl2)
						{
							fprintf(FRES, "#Level: %s %d, genetic distance: %s, method: %s%s", level[m], rl2 + 1, GD_ESTIMATOR[k], CLUSTER_METHOD[l], g_linebreak_val);
							HCLUSTERING tcluster(clustering_matrix + gdindex[k] * n * n, aregs[rl2], n, l, clustering_memory);
							tcluster.PrintClustering(FRES);
						}
						break;
					}
					}
				}
			}
			if (m == 1)
			{
				delete[] gd_tab[0];
				delete[] gd_tab;
			}
			delete[] clustering_memory;
			delete[] clustering_matrix;
		}
	}

	GDIST_METHOD = 0;
	CloseResFile();

	EvaluationEnd("Hierarchical clustering");

	if (cluster_plot_val == 1)
		RunRscript("cluster_plot.R");
}

/* Calculate genetic distance for Hierarchical clustering using multiple threads */
THREAD(ClusteringThread)
{
	//load ind
	int nthread = g_nthread_val;
	int type = gdindex[0];
	int n = type == 1 ? nind : (type == 2 ? npop : nreg[type - 3]);
	int64 progress = 0;
	byte* estimator = cluster_estimator_val;
	GDIST tbuf;

	VLA_NEW(p1, double, maxK);
	VLA_NEW(p2, double, maxK);
	for (int i = 0; i < n; ++i)
	{
		for (int j = i + 1; j < n; ++j)
		{
			if (progress++ % nthread != threadid) continue;
			switch (type)
			{
			case 1:  tbuf.CalcGD(ainds[i], ainds[j], p1, p2); break;
			case 2:  tbuf.CalcGD(apops[i], apops[j], p1); break;
			case 3:
			default: tbuf.CalcGD(aregs[type - 3][i], aregs[type - 3][j], p1);  break;
			}

			for (int k = 1; k <= (type == 1 ? N_GD_ESTIMATOR - 2 * N_FST_ESTIMATOR : N_GD_ESTIMATOR); ++k)
				if (estimator[k])
					*(clustering_matrix + n * n * gdindex[k] + j * n + i) =
					*(clustering_matrix + n * n * gdindex[k] + i * n + j) =
					*(&tbuf.Nei1972 + k - 1);

			PROGRESS_VALUE++;
		}
	}
	VLA_DELETE(p1);
	VLA_DELETE(p2);
}