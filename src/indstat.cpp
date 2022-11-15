/* Individual Statistics Functions */

#pragma once
#include "vcfpop.h"

template struct IND<double>;
template struct IND<float >;

template TARGET void IND<double>::IndividualStatisticsHeader(FILE* fout);
template TARGET void IND<float >::IndividualStatisticsHeader(FILE* fout);
template TARGET void IND<double>::PrintIndividualStatistics(FILE* fout);
template TARGET void IND<float >::PrintIndividualStatistics(FILE* fout);

template TARGET void CalcIndstat<double>();
template TARGET void CalcIndstat<float >();

#define extern 

#undef extern 

/* Calculate individual statistics */
template<typename REAL>
TARGET void CalcIndstat()
{
	if (!indstat) return;
	if (ad) Exit("\nError: individual statistics (-indstat) is incompatible with allelic depth (-ad) option.\n");
	EvaluationBegin();

	OpenResFile("-indstat", "Individual statistics");
	OpenTempFiles(g_nthread_val, ".indstat");

	IND<REAL>::IndividualStatisticsHeader(FRES);
	RunThreads(&IndividualStatisticsThread<REAL>, NULL, NULL, nind, nind,
		"\nCalculating individual statistics:\n", g_nthread_val, true);

	JoinTempFiles(g_nthread_val);
	CloseResFile();

	EvaluationEnd("Individual statistics estimation");
}

/* Calculate individual statistics using multiple threads */
THREAD2(IndividualStatisticsThread)
{
	uint nthread = g_nthread_val;
	double nsec = nind / (double)nthread + 1e-8;
	uint st = (uint64)(threadid * nsec), ed = (uint64)((threadid + 1) * nsec);

	for (uint i = st; i < ed; ++i)
	{
		ainds[i]->PrintIndividualStatistics(TEMP_FILES[threadid]);
		PROGRESS_VALUE++;
	}
}

/* Write header row for individual statistics */
template<typename REAL>
TARGET void IND<REAL>::IndividualStatisticsHeader(FILE* fout)
{
	char name_buf[NAME_BUF_LEN];
	fprintf(fout, "%s%s%c%c%c%c%c", g_linebreak_val, g_linebreak_val, g_delimiter_val, g_delimiter_val, g_delimiter_val, g_delimiter_val, g_delimiter_val);
	int tn = (indstat_type_val[1]) +
		CountNonZero(indstat_ref_val, N_MAX_OPTION) *
		(indstat_type_val[2] * CountNonZero(indstat_model_val, N_MAX_OPTION) +
			(indstat_type_val[3] + indstat_type_val[4]) * CountNonZero(indstat_estimator_val, N_MAX_OPTION)) - 1;

	int64 tloc = indstat_locus_val[1] + indstat_locus_val[2] * nloc;

	if (indstat_locus_val[1])
	{
		for (int rl = 0; rl < lreg; ++rl)
			fprintf(fout, "%c", g_delimiter_val);
		fprintf(fout, "%cAll loci", g_delimiter_val);
		for (int i = 0; i < tn; ++i)
			fprintf(fout, "%c", g_delimiter_val);
	}

	if (indstat_locus_val[2]) for (int64 l = 0; l < nloc; ++l)
	{
		fprintf(fout, "%c%s", g_delimiter_val, GetLoc(l).GetNameStr(name_buf));
		for (int i = 0; i < tn; ++i)
			fprintf(fout, "%c", g_delimiter_val);
	}

	fprintf(fout, "%sInd%cPop", g_linebreak_val, g_delimiter_val);
	for (int rl = 0; rl < lreg; ++rl)
		fprintf(fout, "%cRegL%d", g_delimiter_val, rl);
	fprintf(fout, "%c#typed%c#miss%cPloidy%c#Hap", g_delimiter_val, g_delimiter_val, g_delimiter_val, g_delimiter_val);

	for (int64 l = 0; l < tloc; ++l)
	{
		if (indstat_type_val[1])
			fprintf(fout, "%cH-idx", g_delimiter_val);

		for (int i = 1; i <= N_DRE_MODEL; ++i)
		{
			if (indstat_model_val[i] == 0) continue;
			if (indstat_type_val[2])
			{
				if (indstat_ref_val[1]) fprintf(fout, "%clnpg_pop_%s", g_delimiter_val, DRE_MODEL[i]);
				if (indstat_ref_val[2]) fprintf(fout, "%clnpg_reg_%s", g_delimiter_val, DRE_MODEL[i]);
				if (indstat_ref_val[3]) fprintf(fout, "%clnpg_tot_%s", g_delimiter_val, DRE_MODEL[i]);
			}
		}

		if (indstat_type_val[3])
		{
			if (indstat_ref_val[1])
			{
				if (indstat_estimator_val[1]) fprintf(fout, "%cF_pop_RI", g_delimiter_val);
				if (indstat_estimator_val[2]) fprintf(fout, "%cF_pop_LO", g_delimiter_val);
				if (indstat_estimator_val[3]) fprintf(fout, "%cF_pop_WE", g_delimiter_val);
			}
			if (indstat_ref_val[2])
			{
				if (indstat_estimator_val[1]) fprintf(fout, "%cF_reg_RI", g_delimiter_val);
				if (indstat_estimator_val[2]) fprintf(fout, "%cF_reg_LO", g_delimiter_val);
				if (indstat_estimator_val[3]) fprintf(fout, "%cF_reg_WE", g_delimiter_val);
			}
			if (indstat_ref_val[3])
			{
				if (indstat_estimator_val[1]) fprintf(fout, "%cF_tot_RI", g_delimiter_val);
				if (indstat_estimator_val[2]) fprintf(fout, "%cF_tot_LO", g_delimiter_val);
				if (indstat_estimator_val[3]) fprintf(fout, "%cF_tot_WE", g_delimiter_val);
			}
		}

		if (indstat_type_val[4])
		{
			if (indstat_ref_val[1])
			{
				if (indstat_estimator_val[1]) fprintf(fout, "%cTheta_pop_RI", g_delimiter_val);
				if (indstat_estimator_val[2]) fprintf(fout, "%cTheta_pop_LO", g_delimiter_val);
				if (indstat_estimator_val[3]) fprintf(fout, "%cTheta_pop_WE", g_delimiter_val);
			}
			if (indstat_ref_val[2])
			{
				if (indstat_estimator_val[1]) fprintf(fout, "%cTheta_reg_RI", g_delimiter_val);
				if (indstat_estimator_val[2]) fprintf(fout, "%cTheta_reg_LO", g_delimiter_val);
				if (indstat_estimator_val[3]) fprintf(fout, "%cTheta_reg_WE", g_delimiter_val);
			}
			if (indstat_ref_val[3])
			{
				if (indstat_estimator_val[1]) fprintf(fout, "%cTheta_tot_RI", g_delimiter_val);
				if (indstat_estimator_val[2]) fprintf(fout, "%cTheta_tot_LO", g_delimiter_val);
				if (indstat_estimator_val[3]) fprintf(fout, "%cTheta_tot_WE", g_delimiter_val);
			}
		}
	}
}

/* Write result row for individual statistics */
template<typename REAL>
TARGET void IND<REAL>::PrintIndividualStatistics(FILE* fout)
{
	int ntype = 0, minv = 99999, maxv = 0, nhaplo = 0, nmiss = 0;
	double hidx = 0;
	double f_tot_Ritland = 0, f_pop_Ritland = 0, f_reg_Ritland[N_MAX_REG];
	double f_tot_Loiselle = 0, f_pop_Loiselle = 0, f_reg_Loiselle[N_MAX_REG];
	double f_tot_Weir = 0, f_pop_Weir = 0, f_reg_Weir[N_MAX_REG];
	double t_tot_Ritland = 0, t_pop_Ritland = 0, t_reg_Ritland[N_MAX_REG];
	double t_tot_Loiselle = 0, t_pop_Loiselle = 0, t_reg_Loiselle[N_MAX_REG];
	double t_tot_Weir = 0, t_pop_Weir = 0, t_reg_Weir[N_MAX_REG];

	POP<REAL>* ttot = total_pop;
	POP<REAL>* tpop = apops[popid];
	POP<REAL>* treg[N_MAX_REG]; //20200505
	if (lreg >= 1) treg[0] = aregs[0][tpop->rid]; //20200505

	for (int rl = 1; rl < lreg; ++rl)
		treg[rl] = aregs[rl][treg[rl - 1]->rid];

	for (int64 l = 0; l < nloc; ++l)
	{
		GENOTYPE& gt = GetGenotype(l);//fine
		if (gt.Nalleles() == 0) { nmiss++; continue; }
		ntype++;

		int v = gt.Ploidy();
		nhaplo += v;
		if (v > maxv) maxv = v;
		if (v < minv) minv = v;
		hidx += gt.HIndex();
	}
	hidx /= ntype;

	fprintf(fout, "%s%s%c%s%c",
		g_linebreak_val, name, g_delimiter_val,
		tpop->name, g_delimiter_val);

	for (int rl = 0; rl < lreg; ++rl)
		fprintf(fout, "%s%c", treg[rl]->name, g_delimiter_val);

	fprintf(fout, "%d%c%d%c%d-%d%c%d",
		ntype, g_delimiter_val,
		nmiss, g_delimiter_val,
		minv, maxv, g_delimiter_val,
		nhaplo);

	if (indstat_locus_val[1])
	{
		if (indstat_type_val[1]) {
			fprintf(fout, "%c", g_delimiter_val); WriteReal(fout, hidx);
		}

		for (int i = 1; i <= N_DRE_MODEL; ++i)
		{
			if (indstat_model_val[i] == 0) continue;
			if (indstat_type_val[2]) {
				if (indstat_ref_val[1]) { fprintf(fout, "%c", g_delimiter_val); WriteReal(fout, GenoFreq(tpop, i, -1, 0)); }
				if (indstat_ref_val[2]) for (int rl = 0; rl < lreg; ++rl) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, GenoFreq(treg[rl], i, -1, 0)); }
				if (indstat_ref_val[3]) { fprintf(fout, "%c", g_delimiter_val); WriteReal(fout, GenoFreq(ttot, i, -1, 0)); }
			}
		}

		if (indstat_type_val[3] || indstat_type_val[4]) {
			if (indstat_ref_val[1])
				Theta(tpop, f_pop_Ritland, f_pop_Loiselle, f_pop_Weir, t_pop_Ritland, t_pop_Loiselle, t_pop_Weir);
			if (indstat_ref_val[2]) for (int rl = 0; rl < lreg; ++rl)
				Theta(treg[rl], f_reg_Ritland[rl], f_reg_Loiselle[rl], f_reg_Weir[rl], t_reg_Ritland[rl], t_reg_Loiselle[rl], t_reg_Weir[rl]);
			if (indstat_ref_val[3])
				Theta(ttot, f_tot_Ritland, f_tot_Loiselle, f_tot_Weir, t_tot_Ritland, t_tot_Loiselle, t_tot_Weir);
		}

		if (indstat_type_val[3]) {
			if (indstat_ref_val[1]) {
				if (indstat_estimator_val[1]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, f_pop_Ritland); }
				if (indstat_estimator_val[2]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, f_pop_Loiselle); }
				if (indstat_estimator_val[3]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, f_pop_Weir); }
			}
			if (indstat_ref_val[2]) for (int rl = 0; rl < lreg; ++rl) {
				if (indstat_estimator_val[1]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, f_reg_Ritland[rl]); }
				if (indstat_estimator_val[2]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, f_reg_Loiselle[rl]); }
				if (indstat_estimator_val[3]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, f_reg_Weir[rl]); }
			}
			if (indstat_ref_val[3]) {
				if (indstat_estimator_val[1]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, f_tot_Ritland); }
				if (indstat_estimator_val[2]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, f_tot_Loiselle); }
				if (indstat_estimator_val[3]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, f_tot_Weir); }
			}
		}

		if (indstat_type_val[4]) {
			if (indstat_ref_val[1]) {
				if (indstat_estimator_val[1]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, t_pop_Ritland); }
				if (indstat_estimator_val[2]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, t_pop_Loiselle); }
				if (indstat_estimator_val[3]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, t_pop_Weir); }
			}
			if (indstat_ref_val[2]) for (int rl = 0; rl < lreg; ++rl) {
				if (indstat_estimator_val[1]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, t_reg_Ritland[rl]); }
				if (indstat_estimator_val[2]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, t_reg_Loiselle[rl]); }
				if (indstat_estimator_val[3]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, t_reg_Weir[rl]); }
			}
			if (indstat_ref_val[3]) {
				if (indstat_estimator_val[1]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, t_tot_Ritland); }
				if (indstat_estimator_val[2]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, t_tot_Loiselle); }
				if (indstat_estimator_val[3]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, t_tot_Weir); }
			}
		}
	}

	/*******************************************************************/

	if (indstat_locus_val[2]) for (int64 l = 0; l < nloc; ++l)
	{
		if (indstat_type_val[1]) {
			fprintf(fout, "%c", g_delimiter_val); WriteReal(fout, GetGenotype(l).HIndex());//fine
		}

		for (int i = 1; i <= N_DRE_MODEL; ++i)
		{
			if (indstat_model_val[i] == 0) continue;
			if (indstat_type_val[2]) {
				if (indstat_ref_val[1]) { fprintf(fout, "%c", g_delimiter_val); WriteReal(fout, GenoFreq(tpop, i, l, 0)); }
				if (indstat_ref_val[2]) for (int rl = 0; rl < lreg; ++rl) { fprintf(fout, "%c", g_delimiter_val); WriteReal(fout, GenoFreq(treg[rl], i, l, 0)); }
				if (indstat_ref_val[3]) { fprintf(fout, "%c", g_delimiter_val); WriteReal(fout, GenoFreq(ttot, i, l, 0)); }
			}
		}

		if (indstat_type_val[3] || indstat_type_val[4]) {
			if (indstat_ref_val[1]) Theta(tpop, f_pop_Ritland, f_pop_Loiselle, f_pop_Weir, t_pop_Ritland, t_pop_Loiselle, t_pop_Weir, l);
			if (indstat_ref_val[2]) for (int rl = 0; rl < lreg; ++rl) Theta(treg[rl], f_reg_Ritland[rl], f_reg_Loiselle[rl], f_reg_Weir[rl], t_reg_Ritland[rl], t_reg_Loiselle[rl], t_reg_Weir[rl], l);
			if (indstat_ref_val[3]) Theta(ttot, f_tot_Ritland, f_tot_Loiselle, f_tot_Weir, t_tot_Ritland, t_tot_Loiselle, t_tot_Weir, l);
		}

		if (indstat_type_val[3]) {
			if (indstat_ref_val[1]) {
				if (indstat_estimator_val[1]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, f_pop_Ritland); }
				if (indstat_estimator_val[2]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, f_pop_Loiselle); }
				if (indstat_estimator_val[3]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, f_pop_Weir); }
			}
			if (indstat_ref_val[2]) for (int rl = 0; rl < lreg; ++rl) {
				if (indstat_estimator_val[1]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, f_reg_Ritland[rl]); }
				if (indstat_estimator_val[2]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, f_reg_Loiselle[rl]); }
				if (indstat_estimator_val[3]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, f_reg_Weir[rl]); }
			}
			if (indstat_ref_val[3]) {
				if (indstat_estimator_val[1]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, f_tot_Ritland); }
				if (indstat_estimator_val[2]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, f_tot_Loiselle); }
				if (indstat_estimator_val[3]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, f_tot_Weir); }
			}
		}

		if (indstat_type_val[4]) {
			if (indstat_ref_val[1]) {
				if (indstat_estimator_val[1]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, t_pop_Ritland); }
				if (indstat_estimator_val[2]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, t_pop_Loiselle); }
				if (indstat_estimator_val[3]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, t_pop_Weir); }
			}
			if (indstat_ref_val[2]) for (int rl = 0; rl < lreg; ++rl) {
				if (indstat_estimator_val[1]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, t_reg_Ritland[rl]); }
				if (indstat_estimator_val[2]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, t_reg_Loiselle[rl]); }
				if (indstat_estimator_val[3]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, t_reg_Weir[rl]); }
			}
			if (indstat_ref_val[3]) {
				if (indstat_estimator_val[1]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, t_tot_Ritland); }
				if (indstat_estimator_val[2]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, t_tot_Loiselle); }
				if (indstat_estimator_val[3]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, t_tot_Weir); }
			}
		}
	}
}

