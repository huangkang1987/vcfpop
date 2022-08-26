/* Genetic Distance Functions */

#pragma once
#include "vcfpop.h"

#ifndef _GDDIST
	/* Write column format header row for genetic distance estimation */
	TARGET void GDIST::ColumnPrintHeader()
	{
		fprintf(FRES, "%s%sA", g_linebreak_val, g_linebreak_val);
		for (int rl = gdist_type - 2; rl < lreg; ++rl)
			if (rl >= 0)
				fprintf(FRES, "%cRegL%d", g_delimiter_val, rl + 1);
			else
				fprintf(FRES, "%cPop", g_delimiter_val);

		fprintf(FRES, "%cB", g_delimiter_val);

		for (int rl = gdist_type - 2; rl < lreg; ++rl)
			if (rl >= 0)
				fprintf(FRES, "%cRegL%d", g_delimiter_val, rl + 1);
			else
				fprintf(FRES, "%cPop", g_delimiter_val);

		for (int k = 1; k <= (gdist_type == 1 ? N_GD_ESTIMATOR - 2 * N_FST_ESTIMATOR : N_GD_ESTIMATOR); ++k)
			if (gdist_estimator_val[k])
				fprintf(FRES, "%c%s", g_delimiter_val, GD_ESTIMATOR[k]);
	}

	/* Write column format result row for genetic distance estimation */
	TARGET void GDIST::ColumnPrintLine(int i, int j)
	{
		POP* tp = NULL;
		switch (gdist_type)
		{
		case 1:
			fprintf(FRES, "%s%s", g_linebreak_val, ainds[i]->name);
			tp = apops[ainds[i]->popid];
			break;
		case 2:
			fprintf(FRES, "%s%s", g_linebreak_val, apops[i]->name);
			tp = lreg >= 0 ? aregs[0][apops[i]->rid] : NULL;
			break;
		case 3:
		default:
			fprintf(FRES, "%s%s", g_linebreak_val, aregs[gdist_type - 3][i]->name);
			tp = aregs[gdist_type - 2][aregs[gdist_type - 3][i]->rid];
			break;
		}

		for (int rl = gdist_type - 2; rl < lreg; ++rl)
		{
			fprintf(FRES, "%c%s", g_delimiter_val, tp->name);
			tp = aregs[rl + 1][tp->rid];
		}

		switch (gdist_type)
		{
		case 1:
			fprintf(FRES, "%c%s", g_delimiter_val, ainds[j]->name);
			tp = apops[ainds[j]->popid];
			break;
		case 2:
			fprintf(FRES, "%c%s", g_delimiter_val, apops[j]->name);
			tp = lreg >= 0 ? aregs[0][apops[j]->rid] : NULL;
			break;
		case 3:
		default:
			fprintf(FRES, "%c%s", g_delimiter_val, aregs[gdist_type - 3][j]->name);
			tp = aregs[gdist_type - 2][aregs[gdist_type - 3][j]->rid];
			break;
		}

		for (int rl = gdist_type - 2; rl < lreg; ++rl)
		{
			fprintf(FRES, "%c%s", g_delimiter_val, tp->name);
			tp = aregs[rl + 1][tp->rid];
		}

		for (int k = 1; k <= (gdist_type == 1 ? N_GD_ESTIMATOR - 2 * N_FST_ESTIMATOR : N_GD_ESTIMATOR); ++k)
			if (gdist_estimator_val[k])
			{
				fprintf(FRES, "%c", g_delimiter_val);
				WriteReal(FRES, *((&Nei1972) + k - 1));
			}
	}

	/* Write matrix format header for genetic distance estimation */
	TARGET void GDIST::MatrixPrintMatrixHeader(int k, int n)
	{
		if (gdist_estimator_val[k] == 0) return;
		fprintf(TEMP_FILES[k], "%s%s%s", g_linebreak_val, g_linebreak_val, GD_ESTIMATOR[k]);
		for (int i = 0; i < n; ++i)
		{
			switch (gdist_type)
			{
			case 1: fprintf(TEMP_FILES[k], "%c%s", g_delimiter_val, ainds[i]->name); break;
			case 2: fprintf(TEMP_FILES[k], "%c%s", g_delimiter_val, apops[i]->name);  break;
			case 3:
			default:
				fprintf(TEMP_FILES[k], "%c%s", g_delimiter_val, aregs[gdist_type - 3][i]->name);  break;
			}
		}
	}

	/* Write matrix format row header for genetic distance estimation */
	TARGET void GDIST::MatrixPrintRowHeader(int k, int i)
	{
		if (gdist_estimator_val[k] == 0) return;
		switch (gdist_type)
		{
		case 1: fprintf(TEMP_FILES[k], "%s%s", g_linebreak_val, ainds[i]->name); break;
		case 2: fprintf(TEMP_FILES[k], "%s%s", g_linebreak_val, apops[i]->name);  break;
		case 3:
		default:
			fprintf(TEMP_FILES[k], "%s%s", g_linebreak_val, aregs[gdist_type - 3][i]->name);  break;
		}
	}

	/* Write matrix format grid for genetic distance estimation */
	TARGET void GDIST::MatrixPrintCell(int k)
	{
		if (gdist_estimator_val[k] == 0) return;
		fprintf(TEMP_FILES[k], "%c", g_delimiter_val);
		WriteReal(TEMP_FILES[k], *((&Nei1972) + k - 1));
	}

	/* Use total allele frequency as the missing data */
	TARGET void GDIST::GetMissingFreq(GENOTYPE& gt, int64 l, double* p, int k)
	{
		if (gt.Nalleles())
		{
			gt.GetFreq(p, k);
			return;
		}

		if (total_pop->loc_stat1[l].nhaplo)
			SetVal(p, total_pop->GetFreq(l), k);
	}

	/* Calculate genetic distance between genotypes and save in gdtab */
	TARGET void GDIST::CacheIndGD()
	{
		byte* estimator = NULL;
		switch (GDIST_METHOD)
		{
		default: break;
		case 1: estimator = gdist_estimator_val; break;
		case 2: estimator = pcoa_estimator_val; break;
		case 3: estimator = cluster_estimator_val; break;
		}

		if ((estimator[12] || estimator[20]) && abs(g_format_val) <= BCF)
			Exit("\nError: Reynolds_Slatkin1995 or Slatkin_Slatkin1995 genetic distance estimator uses Stepwise mutation model (smm), which can only be applied for non-vcf input file, and should use the allele size as the identifier. \n");
		if ((estimator[6]) && abs(g_format_val) <= BCF)
			Exit("\nError: Goldstein1995 genetic distance estimator uses Stepwise mutation model (smm), which can only be applied for non-vcf input file, and should use the allele size as the identifier. \n");

		if (ad) Exit("\nError: individual genetic distance (-gdist_level=ind -pcoa_level=ind -cluster_level=ind) is incompatible with allelic depth (-ad) option.\n");

		double* P1 = new double[maxK * g_nthread_val], * P2 = new double[maxK * g_nthread_val];
		gd_tab = new INDGD * [nloc];
		SetZero(gd_tab, nloc);

		int64 tabsize = 0;
		for (int64 l = 0; l < nloc; ++l)
		{
			int ngeno = GetLoc(l).ngeno;
			if (total_pop->loc_stat1[l].nhaplo == 0 || ngeno >= N_MAX_GDTAB) continue;
			tabsize += ((ngeno * (ngeno + 1)) >> 1);
		}

		INDGD* tabp = gd_tab[0] = new INDGD[tabsize];
		SetZero(tabp, tabsize);
		for (int64 l = 0; l < nloc; ++l)
		{
			int ngeno = GetLoc(l).ngeno;
			if (total_pop->loc_stat1[l].nhaplo == 0 || ngeno >= N_MAX_GDTAB) continue;
			gd_tab[l] = tabp;
			tabp += ((ngeno * (ngeno + 1)) >> 1);
		}

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int64 l = 0; l < nloc; ++l)
		{
			int ngeno = GetLoc(l).ngeno;
			if (total_pop->loc_stat1[l].nhaplo == 0 || ngeno >= N_MAX_GDTAB) continue;

			threadid = omp_get_thread_num();
			double* p1 = P1 + threadid * maxK, * p2 = P2 + threadid * maxK;

			int k = GetLoc(l).k;
			GENOTYPE* gtab = GetLoc(l).GetGtab();
			INDGD* tab = gd_tab[l];
			for (int gid1 = 0, idx = 0; gid1 < ngeno; ++gid1)
			{
				GENOTYPE& gt1 = gtab[gid1];
				for (int gid2 = gid1; gid2 < ngeno; ++gid2, ++idx)
				{
					GENOTYPE& gt2 = gtab[gid2];

					INDGD& tgd = tab[idx];
					tgd.ABtype = gt1.Nalleles() && gt2.Nalleles() ? 1 : 0;

					if (tgd.ABtype == 1 || gdist_weightmissing_val == 1)
					{
						GetMissingFreq(gt1, l, p1, k);
						GetMissingFreq(gt2, l, p2, k);

						double Sx1 = 0, Sx2 = 0;
						ushort* alen = GetLoc(l).GetAlenArray();

						for (int i = 0; i < k; ++i)
						{
							if (estimator[1] || estimator[7])
							{
								//Nei1972, Nei1973
								tgd.Jx1 += p1[i] * p1[i];
								tgd.Jx2 += p2[i] * p2[i];
								tgd.Jxy += p1[i] * p2[i];
							}

							if (estimator[2])
								tgd.Cavalli1967 += MySqrt(Max(DOUBLE_UNDERFLOW, p1[i] * p2[i]));

							if (estimator[3])
							{
								tgd.t1 += (p1[i] - p2[i]) * (p1[i] - p2[i]);
								tgd.t2 += p1[i] * p2[i];
							}

							if (estimator[4])
								tgd.Nei1983 += MySqrt(Max(DOUBLE_UNDERFLOW, p1[i] * p2[i]));

							if (estimator[5])
								tgd.Euclidean += (p1[i] - p2[i]) * (p1[i] - p2[i]);

							if (abs(g_format_val) > BCF && estimator[6])
							{
								Sx1 += p1[i] * alen[i];
								Sx2 += p2[i] * alen[i];
							}

							if (estimator[8])
								tgd.Roger1972 += (p1[i] - p2[i]) * (p1[i] - p2[i]) * 0.5;
						}

						if (estimator[5])
							tgd.Euclidean = Max(DOUBLE_UNDERFLOW, tgd.Euclidean);

						if (estimator[6])
							tgd.Goldstein1995 = (Sx1 - Sx2) * (Sx1 - Sx2);

						if (estimator[8])
							tgd.Roger1972 = MySqrt(Max(DOUBLE_UNDERFLOW, tgd.Roger1972));
					}
				}
			}
		}

		delete[] P1;
		delete[] P2;
	}

	/* Calculate genetic distance between two individuals */
	TARGET void GDIST::CalcGD(IND* a, IND* b, double* p1, double* p2)
	{
		byte* estimator = NULL;
		switch (GDIST_METHOD)
		{
		default: break;
		case 1: estimator = gdist_estimator_val; break;
		case 2: estimator = pcoa_estimator_val; break;
		case 3: estimator = cluster_estimator_val; break;
		}

		if ((estimator[12] || estimator[20]) && abs(g_format_val) <= BCF)
			Exit("\nError: Reynolds_Slatkin1995 or Slatkin_Slatkin1995 genetic distance estimator uses Stepwise mutation model (smm), which can only be applied for non-vcf input file, and should use the allele size as the identifier. \n");
		if ((estimator[6]) && abs(g_format_val) <= BCF)
			Exit("\nError: Goldstein1995 genetic distance estimator uses Stepwise mutation model (smm), which can only be applied for non-vcf input file, and should use the allele size as the identifier. \n");

		if (ad) Exit("\nError: individual genetic distance (-gdist_level=ind -pcoa_level=ind -cluster_level=ind) is incompatible with allelic depth (-ad) option.\n");

		SetZero(&Nei1972, N_GD_ESTIMATOR);

		if (a == b) return;
		INDGD sumgd, tgd;
		SetZero(&sumgd, 1);
		int aid = a->indid, bid = b->indid;
		int64 eL = 0;

		for (int64 l = 0; l < nloc; ++l)
		{
			if (total_pop->loc_stat1[l].nhaplo == 0) continue;

			int k = GetLoc(l).k;
			int gid1 = aid, gid2 = bid, ngeno = GetLoc(l).ngeno;
			IND::GetDyadGenotypeIdx(gid1, gid2, l);
			bool cis = gid1 <= gid2;
			INDGD* tab2 = gd_tab[l];

			if (ngeno < N_MAX_GDTAB)
			{
				INDGD& gd = tab2[GetLowerTriangularId(gid1, gid2, ngeno)];
				Add((double*)&sumgd, (double*)&gd, N_INDGD);
				if (!cis)
				{
					sumgd.Jx1 += gd.Jx2 - gd.Jx1;
					sumgd.Jx2 += gd.Jx1 - gd.Jx2;
				}
				if (gd.ABtype == 1 || gdist_weightmissing_val == 1)
					eL++;
				continue;
			}

			SetZero(&tgd, 1);

			GENOTYPE* gtab = GetLoc(l).GetGtab();
			GENOTYPE& gt1 = gtab[gid1], & gt2 = gtab[gid2];

			if (gt1.Nalleles() && gt2.Nalleles())
				tgd.ABtype = 1;

			if (tgd.ABtype == 1 || gdist_weightmissing_val == 1)
			{
				eL++;
				GetMissingFreq(gt1, l, p1, k);
				GetMissingFreq(gt2, l, p2, k);

				double Sx1 = 0, Sx2 = 0;
				ushort* alen = GetLoc(l).GetAlenArray();

				for (int i = 0; i < k; ++i)
				{
					if (estimator[1] || estimator[7])
					{
						//Nei1972, Nei1973
						tgd.Jx1 += p1[i] * p1[i];
						tgd.Jx2 += p2[i] * p2[i];
						tgd.Jxy += p1[i] * p2[i];
					}

					if (estimator[2])
						tgd.Cavalli1967 += MySqrt(Max(DOUBLE_UNDERFLOW, p1[i] * p2[i]));

					if (estimator[3])
					{
						tgd.t1 += (p1[i] - p2[i]) * (p1[i] - p2[i]);
						tgd.t2 += p1[i] * p2[i];
					}

					if (estimator[4])
						tgd.Nei1983 += MySqrt(Max(DOUBLE_UNDERFLOW, p1[i] * p2[i]));

					if (estimator[5])
						tgd.Euclidean += (p1[i] - p2[i]) * (p1[i] - p2[i]);

					if (abs(g_format_val) > BCF && estimator[6])
					{
						Sx1 += p1[i] * alen[i];
						Sx2 += p2[i] * alen[i];
					}

					if (estimator[8])
						tgd.Roger1972 += (p1[i] - p2[i]) * (p1[i] - p2[i]) * 0.5;
				}
				if (estimator[5])
					tgd.Euclidean = Max(DOUBLE_UNDERFLOW, tgd.Euclidean);

				if (estimator[6])
					tgd.Goldstein1995 = (Sx1 - Sx2) * (Sx1 - Sx2);

				if (estimator[8])
					tgd.Roger1972 = MySqrt(Max(DOUBLE_UNDERFLOW, tgd.Roger1972));

				Add((double*)&sumgd, (double*)&tgd, N_INDGD);
			}
		}

		if (eL == 0)
		{
			SetZero(&Nei1972, N_GD_ESTIMATOR);
			return;
		}

		double inveL = 1.0 / eL;
		if (estimator[1] || estimator[7])
		{
			sumgd.Jxy *= inveL;
			sumgd.Jx1 *= inveL;
			sumgd.Jx2 *= inveL;

			Nei1972 = -log(sumgd.Jxy / MySqrt(sumgd.Jx1 * sumgd.Jx2));
			Nei1974 = (sumgd.Jx1 + sumgd.Jx2) * 0.5 - sumgd.Jxy;
		}

		if (estimator[2])
			Cavalli1967 = TWODIVPISQ2 * MySqrt(1 - sumgd.Cavalli1967 * inveL);

		if (estimator[3])
		{
			sumgd.t2 = 2 * (eL - sumgd.t2);
			Reynolds1983 = sumgd.t2 > 0 ? MySqrt(sumgd.t1 / sumgd.t2) : 0;
		}

		if (estimator[4])
			Nei1983 = 1 - sumgd.Nei1983 * inveL;

		if (estimator[5])
			Euclidean = sumgd.Euclidean = MySqrt(sumgd.Euclidean);

		if (estimator[6])
			Goldstein1995 = sumgd.Goldstein1995 * inveL;

		if (estimator[8])
			Roger1972 = sumgd.Roger1972 * inveL;
	}

	/* Calculate genetic distance between two populations/regions */
	TARGET void GDIST::CalcGD(POP* a, POP* b, double* buf)
{
	byte* estimator = NULL;
	switch (GDIST_METHOD)
	{
	default: break;
	case 1: estimator = gdist_estimator_val; break;
	case 2: estimator = pcoa_estimator_val; break;
	case 3: estimator = cluster_estimator_val; break;
	}
	if ((estimator[12] || estimator[20]) && abs(g_format_val) <= BCF)
		Exit("\nError: Reynolds_Slatkin1995 or Slatkin_Slatkin1995 genetic distance estimator uses Stepwise mutation model (smm), which can only be applied for non-vcf input file, and should use the allele size as the identifier. \n");
	if ((estimator[6]) && abs(g_format_val) <= BCF)
		Exit("\nError: Goldstein1995 genetic distance estimator uses Stepwise mutation model (smm), which can only be applied for non-vcf input file, and should use the allele size as the identifier. \n");

	double Jx1 = 0, Jx2 = 0, Jxy = 0, t1 = 0, t2 = 0;
	SetZero(&Nei1972, N_GD_ESTIMATOR);

	if (a == b) return;
	double ABtype = 0;

	for (int64 l = 0; l < nloc; ++l)
	{
		LOCSTAT1* f1 = a->loc_stat1 + l, * f2 = b->loc_stat1 + l;
		if (f1->nhaplo && f2->nhaplo) ABtype++;
		else continue;

		double* p1 = a->GetFreq(l), * p2 = b->GetFreq(l);
		int k2 = GetLoc(l).k;
		ushort* alen = GetLoc(l).GetAlenArray();

		double Roger1972t = 0, Euclideant = 0, Sx1 = 0, Sx2 = 0;
		for (int i = 0; i < k2; ++i)
		{
			if (estimator[1] || estimator[7])
			{
				//Nei1972, Nei1973
				Jx1 += p1[i] * p1[i];
				Jx2 += p2[i] * p2[i];
				Jxy += p1[i] * p2[i];
			}

			if (estimator[2])
				Cavalli1967 += MySqrt(Max(DOUBLE_UNDERFLOW, p1[i] * p2[i]));

			if (estimator[3])
			{
				t1 += (p1[i] - p2[i]) * (p1[i] - p2[i]);
				t2 += p1[i] * p2[i];
			}

			if (estimator[4])
				Nei1983 += MySqrt(Max(DOUBLE_UNDERFLOW, p1[i] * p2[i]));

			if (estimator[5])
				Euclideant += (p1[i] - p2[i]) * (p1[i] - p2[i]);

			if (abs(g_format_val) > BCF && estimator[6])
			{
				Sx1 += p1[i] * alen[i];
				Sx2 += p2[i] * alen[i];
			}

			if (estimator[8])
				Roger1972t += (p1[i] - p2[i]) * (p1[i] - p2[i]) * 0.5;
		}
		Roger1972 += MySqrt(Max(DOUBLE_UNDERFLOW, Roger1972t));
		Euclidean += Max(DOUBLE_UNDERFLOW, Euclideant);

		if (estimator[6])
			Goldstein1995 += (Sx1 - Sx2) * (Sx1 - Sx2);
	}

	if (ABtype == 0)
	{
		SetZero(&Nei1972, N_GD_ESTIMATOR);
		return;
	}

	double inveL = 1.0 / ABtype;
	if (estimator[1] || estimator[7])
	{
		Jxy *= inveL;
		Jx1 *= inveL;
		Jx2 *= inveL;

		Nei1972 = -log(Jxy / MySqrt(Jx1 * Jx2));
		Nei1974 = (Jx1 + Jx2) * 0.5 - Jxy;
	}

	if (estimator[2])
		Cavalli1967 = TWODIVPISQ2 * MySqrt(1 - Cavalli1967 * inveL);

	if (estimator[3])
	{
		t2 = 2 * (ABtype - t2);
		Reynolds1983 = t2 > 0 ? MySqrt(t1 / t2) : 0;
	}

	if (estimator[4])
		Nei1983 = 1 - Nei1983 * inveL;

	if (estimator[5])
		Euclidean = MySqrt(Euclidean);

	if (estimator[6])
		Goldstein1995 *= inveL;

	if (estimator[8])
		Roger1972 *= inveL;

	POP* gs[] = { a, b };
	for (int i = 1; i <= N_FST_ESTIMATOR; ++i)
	{
		if (estimator[8 + i] || estimator[8 + N_FST_ESTIMATOR + i])
		{
			double Fst = FST::FstEstimator(gs, 2, i, NULL, buf);
			if (estimator[8 + i])				    *(double*)(&Slatkin_Nei1973 + i - 1) = Fst / (1 - Fst);
			if (estimator[8 + N_FST_ESTIMATOR + i]) *(double*)(&Reynolds_Nei1973 + i - 1) = -log(1 - Fst);
		}
	}
}
#endif

#define extern 
extern GDIST* gdist_buf;							//Circle buffer for genetic distance estimation, NBUF
extern int gdist_type;								//1 between inds, 2 between pops, 3 + between regions
extern int gdindex[N_GD_ESTIMATOR + 1];				//Index of ith used estimator
extern INDGD** gd_tab;								//Hash table saves the genetic distance between genotypes
#undef extern 

/* Calculate genetic distance */
TARGET void CalcDist()
{
	if (!gdist) return;
	EvaluationBegin();

	GDIST_METHOD = 1;
	OpenResFile("-gdist", "Genetic distance");

	bool isfirst = true;
	gdist_buf = new GDIST[NBUF];

	int64 ntot = 0;
	int64 n[] = { 0, nind * nind * 2, npop * npop * 2 * 100 };

	for (int m = 1; m <= 2; ++m)
		if (gdist_level_val[m])
			ntot += n[m];

	if (gdist_level_val[3])
		for (int rl = 0; rl < lreg; ++rl)
			ntot += nreg[rl] * nreg[rl] * 2 * 100;

	for (int m = 1; m <= 3; ++m)
	{
		if (gdist_level_val[m] == 0)
			continue;

		//Calculate genetic distance table between any two genotypes
		if (m == 1) GDIST::CacheIndGD();

		for (int rl = 0; rl < (m <= 2 ? 1 : lreg); ++rl)
		{
			int64 nthis = m <= 2 ? n[m] : nreg[rl] * nreg[rl] * 2 * 100;
			gdist_estimator_val[0] = true;
			OpenTempFiles(N_GD_ESTIMATOR + 1, ".gdist", gdist_estimator_val);
			gdist_type = m + rl;

			SetZero(gdist_buf, NBUF);

			RunThreads(&GeneticDistanceThread, &GeneticDistanceGuard2, &GeneticDistanceGuard1, ntot, nthis,
				"\nCalculating genetic distance:\n", g_nthread_val, isfirst);

			isfirst = false;
			JoinTempFiles(N_GD_ESTIMATOR + 1, gdist_estimator_val);
			gdist_estimator_val[0] = false;
		}

		if (m == 1)
		{
			delete[] gd_tab[0];
			delete[] gd_tab;
		}
	}
	delete[] gdist_buf;
	GDIST_METHOD = 0;
	CloseResFile();

	EvaluationEnd("Genetic distance estimation");

	if (gdist_plot_val == 1)
		RunRscript("gdist_plot.R");
}

/* Write column format genetic distance results in a guard thread */
THREAD(GeneticDistanceGuard1)
{
	int64& ii = progress1 = 0;
	int n = gdist_type == 1 ? nind : (gdist_type == 2 ? npop : nreg[gdist_type - 3]);
	int64 nadd2 = gdist_type == 1 ? 1 : 100;
	if (gdist_fmt_val[2])
		GDIST::ColumnPrintHeader();

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j, ++ii)
		{
			GUARD_BEGIN2

				if (j >= i && gdist_fmt_val[2])
					gdist_buf[ii % NBUF].ColumnPrintLine(i, j);

			PROGRESS_VALUE += nadd2;

			GUARD_END2
		}
	}
}

/* Write matrix format genetic distance results in a guard thread */
THREAD(GeneticDistanceGuard2)
{
	int64& ii = progress2 = 0;
	uint nk = gdist_type == 1 ? N_GD_ESTIMATOR - 2 * N_FST_ESTIMATOR : N_GD_ESTIMATOR;
	uint n = gdist_type == 1 ? nind : (gdist_type == 2 ? npop : nreg[gdist_type - 3]);
	int64 nadd2 = gdist_type == 1 ? 1 : 100;

	if (gdist_fmt_val[1])
		for (uint k = 1; k <= nk; ++k)
			GDIST::MatrixPrintMatrixHeader(k, n);

	for (uint i = 0; i < n; ++i)
	{
		if (gdist_fmt_val[1])
			for (uint k = 1; k <= nk; ++k)
				GDIST::MatrixPrintRowHeader(k, i);

		for (uint j = 0; j < n; ++j, ++ii)
		{
			GUARD_BEGIN2

				GDIST& gd = gdist_buf[ii % NBUF];
			if (gdist_fmt_val[1])
				for (uint k = 1; k <= nk; ++k)
					gd.MatrixPrintCell(k);

			PROGRESS_VALUE += nadd2;

			GUARD_END2
		}
	}
}

/* Calculate genetic distance using multiple threads */
THREAD(GeneticDistanceThread)
{
	int64 ii = 0;
	if (gdist_type == 1)
	{
		VLA_NEW(p1, double, maxK);
		VLA_NEW(p2, double, maxK);
		for (int i = 0; i < nind; ++i)
		{
			for (int j = 0; j < nind; ++j, ++ii)
			{
				THREAD_BEGIN2

					gdist_buf[ii % NBUF].CalcGD(ainds[i], ainds[j], p1, p2);

				THREAD_END2
			}
		}
		VLA_DELETE(p1);
		VLA_DELETE(p2);
	}
	else if (gdist_type >= 2)
	{
		int rl = gdist_type - 3;
		int n = gdist_type == 2 ? npop : nreg[rl];
		POP** tpop = gdist_type == 2 ? apops : aregs[rl];
		VLA_NEW(p1, double, maxK);
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j, ++ii)
			{
				THREAD_BEGIN2

					gdist_buf[ii % NBUF].CalcGD(tpop[i], tpop[j], p1);

				THREAD_END2
			}
		}
		VLA_DELETE(p1);
	}
}