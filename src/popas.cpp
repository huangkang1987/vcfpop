/* Population Assignment Functions */

#pragma once
#include "vcfpop.h"

template TARGET void CalcAssignment<double>();
template TARGET void CalcAssignment<float >();
template TARGET void IND<double>::AssignmentHeader(FILE* fout);
template TARGET void IND<float >::AssignmentHeader(FILE* fout);
template TARGET void IND<double>::WriteAssignment(FILE* fout);
template TARGET void IND<float >::WriteAssignment(FILE* fout);

#define extern 

#undef extern 

/* Calculate population assignment */
template<typename REAL>
TARGET void CalcAssignment()
{
	if (!popas) return;
	if (ad) Exit("\nError: population assignment (-popas) is incompatible with allelic depth (-ad) option.\n");

	EvaluationBegin();
	OpenResFile("-popas", "Population assignment");
	OpenTempFiles(g_nthread_val, ".popas");

	IND<REAL>::AssignmentHeader(FRES);
	RunThreads(&PopulationAssignmentThread<REAL>, NULL, NULL, nind, nind,
		"\nPerforming population assignment:\n", g_nthread_val, true);

	JoinTempFiles(g_nthread_val);
	CloseResFile();

	EvaluationEnd("Population assignment");

	if (popas_plot_val == 1)
		RunRscript("popas_plot.R");
}

/* Calculate population assignment using multiple threads */
THREAD2(PopulationAssignmentThread)
{
	//load ind
	double nsec = nind / (double)g_nthread_val + 1e-8;
	int st = (int)(threadid * nsec), ed = (int)((threadid + 1) * nsec);
	for (int i = st; i < ed; ++i)
	{
		ainds<REAL>[i]->WriteAssignment(TEMP_FILES[threadid]);

		PROGRESS_VALUE++;
	}
}

/* Write header row for population assignment */
template<typename REAL>
TARGET void IND<REAL>::AssignmentHeader(FILE* fout)
{
	fprintf(fout, "%s%sInd%cPop", g_linebreak_val, g_linebreak_val, g_delimiter_val);
	for (int rl = 0; rl < lreg; ++rl)
		fprintf(fout, "%cRegL%d", g_delimiter_val, rl + 1);
	fprintf(fout, "%c#typed%c#miss%cPloidy%c#Hap",
		g_delimiter_val, g_delimiter_val, g_delimiter_val, g_delimiter_val);
	for (int rl = -1; rl <= lreg; ++rl)
	{
		if (popas_level_val[rl == -1 ? 1 : 2] == 0) continue;
		POP<REAL>** grps = rl == -1 ? apops<REAL> : aregs<REAL>[rl];
		int ngrp = rl == -1 ? npop : nreg[rl];

		//level: pop reg
		for (int l2 = 1; l2 <= N_DRE_MODEL; ++l2)
		{
			if (popas_model_val[l2] == 0) continue;
			fprintf(fout, "%cassign_%s", g_delimiter_val, rl == -1 ? "pop" : "reg");
			if (rl >= 0) fprintf(fout, "L%d", rl + 1);
			fprintf(fout, "_%s", DRE_MODEL[l2]);
			for (int id = 0; id < ngrp; ++id)
				fprintf(fout, "%clnpg_%s_%s", g_delimiter_val, grps[id]->name, DRE_MODEL[l2]);
		}
	}
}

/* Write result row for population assignment */
template<typename REAL>
TARGET void IND<REAL>::WriteAssignment(FILE* fout)
{
	int64 ntype = 0, nhaplo = 0, nmiss = 0;
	int ms = npop;
	byte minv = 100, maxv = 0;

	VLA_NEW(lnPg, REAL, ms);
	SetZero(lnPg, ms);

	for (int64 l = 0; l < nloc; ++l)
	{
		GENOTYPE& gt = GetGenotype(l);//fine
		if (gt.Nalleles() == 0)
		{
			nmiss++;
			continue;
		}

		int v = gt.Ploidy();
		ntype++;
		nhaplo += v;

		if (v > maxv) maxv = (byte)v;
		if (v < minv) minv = (byte)v;
	}

	fprintf(fout, "%s%s%c%s",
		g_linebreak_val, name,
		g_delimiter_val, apops<REAL>[popid]->name);

	POP<REAL>* tr = lreg >= 0 ? aregs<REAL>[0][apops<REAL>[popid]->rid] : NULL;
	for (int rl = 0; rl < lreg; ++rl)
	{
		fprintf(fout, "%c%s", g_delimiter_val, tr->name);
		tr = aregs<REAL>[rl + 1][tr->rid];
	}

	fprintf(fout, "%c%lld%c%lld%c%d - %d%c%lld",
		g_delimiter_val, ntype,
		g_delimiter_val, nmiss,
		g_delimiter_val, minv, maxv,
		g_delimiter_val, nhaplo);

	for (int rl = -1; rl <= lreg; ++rl)
	{
		if (popas_level_val[rl == -1 ? 1 : 2] == 0) continue;
		POP<REAL>** grps = rl == -1 ? apops<REAL> : aregs<REAL>[rl];
		int ngrp = rl == -1 ? npop : nreg[rl];

		for (int m2 = 1; m2 <= N_DRE_MODEL; ++m2)
		{
			if (popas_model_val[m2] == 0) continue;

			SetZero(lnPg, ms);

			for (int j = 0; j < ngrp; ++j)
				lnPg[j] = GenoFreq(grps[j], m2, -1, popas_error_val);

			int64 mid = GetMaxID(lnPg, ngrp);
			fprintf(fout, "%c%s", g_delimiter_val, grps[mid]->name);
			for (int j = 0; j < ngrp; ++j)
			{
				fprintf(fout, "%c", g_delimiter_val);
				WriteReal(fout, lnPg[j]);
			}
		}
	}

	VLA_DELETE(lnPg);
}

