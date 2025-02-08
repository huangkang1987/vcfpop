/* Kinship Functions */

#include "vcfpop.h"

template struct KINSHIP<double>;
template struct KINSHIP<float >;

template TARGET void CalcKinship<double>();
template TARGET void CalcKinship<float >();
template TARGET void IND<double>::Theta(POP<double>* grp, double& f_ritland, double& f_loiselle, double& f_weir, double& f_vanraden, double& t_ritland, double& t_loiselle, double& t_weir, double& t_vanraden, int64 loc);
template TARGET void IND<float >::Theta(POP<float >* grp, double& f_ritland, double& f_loiselle, double& f_weir, double& f_vanraden, double& t_ritland, double& t_loiselle, double& t_weir, double& t_vanraden, int64 loc);

template<> KINSHIP<double>* kinship_buf<double>;
template<> KINSHIP<float >* kinship_buf<float >;

#ifndef _KINSHIP
/* Write header row for kinship estimation */
template<typename REAL>
TARGET void KINSHIP<REAL>::ColumnFormatHeader()
{
	fprintf(FRES, "%s%s%s%sA%cpop<REAL>",
		g_linebreak_val, g_linebreak_val,
		cpop<REAL>->name, g_linebreak_val,
		g_delimiter_val);
	for (int rl = 0; rl < lreg; ++rl)
		fprintf(FRES, "%cregL%d", g_delimiter_val, rl + 1);
	fprintf(FRES, "%cB%cpop<REAL>", g_delimiter_val, g_delimiter_val);
	for (int rl = 0; rl < lreg; ++rl)
		fprintf(FRES, "%cregL%d", g_delimiter_val, rl + 1);
	fprintf(FRES, "%cAB_typed%cA_typed%cB_typed", g_delimiter_val, g_delimiter_val, g_delimiter_val);

	for (int k = 1; k <= N_KINSHIP_ESTIMATOR; ++k)
		if (kinship_estimator_val[k])
			fprintf(FRES, "%c%s", g_delimiter_val, KINSHIP_ESTIMATOR[k]);
}

/* Write result row for kinship estimation */
template<typename REAL>
TARGET void KINSHIP<REAL>::ColumnFormatLine(int i, int j)
{
	fprintf(FRES, "%s%s%c%s%c",
		g_linebreak_val,
		ainds<REAL>[i]->name, g_delimiter_val,
		apops<REAL>[ainds<REAL>[i]->popid]->name, g_delimiter_val);

	POP<REAL>* tr = lreg >= 0 ? aregs<REAL>[0][apops<REAL>[ainds<REAL>[i]->popid]->rid] : NULL;
	for (int rl = 0; rl < lreg; ++rl)
	{
		fprintf(FRES, "%s%c", tr->name, g_delimiter_val);
		tr = aregs<REAL>[rl + 1][tr->rid];
	}

	fprintf(FRES, "%s%c%s%c",
		ainds<REAL>[j]->name, g_delimiter_val,
		apops<REAL>[ainds<REAL>[j]->popid]->name, g_delimiter_val);

	tr = lreg >= 0 ? aregs<REAL>[0][apops<REAL>[ainds<REAL>[j]->popid]->rid] : NULL;
	for (int rl = 0; rl < lreg; ++rl)
	{
		fprintf(FRES, "%s%c", tr->name, g_delimiter_val);
		tr = aregs<REAL>[rl + 1][tr->rid];
	}

	fprintf(FRES, "%d%c%d%c%d",
		ABtype, g_delimiter_val,
		Atype, g_delimiter_val,
		Btype);

	for (int k = 1; k <= N_KINSHIP_ESTIMATOR; ++k)
		if (kinship_estimator_val[k])
		{
			fprintf(FRES, "%c", g_delimiter_val);
			WriteReal(FRES, *((&Ritland1996) + k - 1));
		}
}

/* Write matrix format header for kinship estimation */
template<typename REAL>
TARGET void KINSHIP<REAL>::MatrixFormatHeader(int k, int n)
{
	if (kinship_estimator_val[k] == 0) return;
	fprintf(TEMP_FILES[k], "%s%s%s%s%s", g_linebreak_val, g_linebreak_val, cpop<REAL>->name, g_linebreak_val, KINSHIP_ESTIMATOR[k]);

	for (int i = 0; i < n; ++i)
		fprintf(TEMP_FILES[k], "%c%s", g_delimiter_val, cpop<REAL>->inds[i]->name);
}

/* Write matrix format row header for kinship estimation */
template<typename REAL>
TARGET void KINSHIP<REAL>::MatrixFormatRowHeader(int k, int i)
{
	if (kinship_estimator_val[k] == 0) return;
	fprintf(TEMP_FILES[k], "%s%s", g_linebreak_val, cpop<REAL>->inds[i]->name);
}

/* Write matrix format grid for kinship estimation */
template<typename REAL>
TARGET void KINSHIP<REAL>::MatrixFormatCell(int k)
{
	if (kinship_estimator_val[k] == 0) return;
	fprintf(TEMP_FILES[k], "%c", g_delimiter_val);
	WriteReal(TEMP_FILES[k], *((&Ritland1996) + k - 1));
}

/* Calculate relatedness coefficient */
template<typename REAL>
TARGET void KINSHIP<REAL>::CalcKinship(IND<REAL>* a, IND<REAL>* b)
{
	Ritland1996 = Loiselle1995 = Weir1996 = VanRaden2008 = 0;

	ABtype = Atype = Btype = 0;

	byte* estimator = kinship_estimator_val;

	for (int64 l = 0; l < nloc; ++l)
	{
		GENOTYPE& gt1 = a->GetGenotype(l), & gt2 = b->GetGenotype(l);//fine
		if (gt1.Nalleles()) Atype++;
		if (gt2.Nalleles()) Btype++;
		if (gt1.Nalleles() && gt2.Nalleles()) ABtype++;
	}

	if (ABtype == 0)
	{
		Ritland1996 = Loiselle1995 = Weir1996 = VanRaden2008 = NAN;
		return;
	}

	if (estimator[1]) Ritland1996 = RELATEDNESS<REAL>::R_Ritland1996(a, b, false, false);
	if (estimator[2]) Loiselle1995 = RELATEDNESS<REAL>::R_Loiselle1995(a, b, false, false);
	if (estimator[3]) Weir1996 = RELATEDNESS<REAL>::R_Weir1996(a, b, false);
	if (estimator[4]) VanRaden2008 = RELATEDNESS<REAL>::R_VanRaden2008(a, b, false);
}
#endif

#define extern 
template<typename REAL>
extern KINSHIP<REAL>* kinship_buf;			//Circle buffer for kinship estimation, NBUF
#undef extern 

/* Calculate kinship coefficient */
template<typename REAL>
TARGET void CalcKinship()
{
	if (!kinship) return;
	if (ad) Exit("\nError: kinship estimation (-kinship) is incompatible with allelic depth (-ad) option.\n");

	EvaluationBegin();
	OpenResFile("-kinship", "Kinship coefficient");

	bool isfirst = true;
	kinship_buf<REAL> = new KINSHIP<REAL>[NBUF];
	int64 ntot = 0;
	if (kinship_range_val[1]) for (int i = 0; i < npop; ++i)
		ntot += apops<REAL>[i]->nind * apops<REAL>[i]->nind;
	if (kinship_range_val[2])
		for (int rl = 0; rl < lreg; ++rl)
			for (int i = 0; i < nreg[rl]; ++i)
				ntot += aregs<REAL>[rl][i]->nind * aregs<REAL>[rl][i]->nind;
	if (kinship_range_val[3])
		ntot += total_pop<REAL>->nind * total_pop<REAL>->nind;
	ntot <<= 1;

	for (int k = 1; k <= 3; ++k)
	{
		if (kinship_range_val[k] == 0) continue;
		for (int rl = 0; rl < (k == 2 ? lreg : 1); ++rl)
		{
			int n = k == 1 ? npop : (k == 2 ? nreg[rl] : 1);
			for (int i = 0; i < n; ++i)
			{
				OpenTempFiles(N_KINSHIP_ESTIMATOR + 1, ".kinship");
				cpop<REAL> = (k == 1 ? apops<REAL>[i] : (k == 2 ? aregs<REAL>[rl][i] : total_pop<REAL>));
				SetZero(kinship_buf<REAL>, NBUF);

				RunThreads(&KinshipThread<REAL>, &KinshipGuard1<REAL>, &KinshipGuard2<REAL>, ntot, cpop<REAL>->nind * cpop<REAL>->nind * 2,
					"\nEstimating kinship coefficient between individuals:\n", g_nthread_val, isfirst);
				isfirst = false;
				JoinTempFiles(N_KINSHIP_ESTIMATOR + 1);
			}
		}
	}

	DEL(kinship_buf<REAL>);
	CloseResFile();

	EvaluationEnd("Kinship estimation");

	if (kinship_plot_val == 1)
		RunRscript("kinship_plot.R");
}

/* Write column format kinship coefficient results in a guard thread */
THREAD2(KinshipGuard1)
{
	int64& ii = progress1 = 0;
	int n = cpop<REAL>->nind;
	if (kinship_fmt_val[2])
		KINSHIP<REAL>::ColumnFormatHeader();

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j, ++ii)
		{
			GUARD_BEGIN2

			KINSHIP<REAL>& ks = kinship_buf<REAL>[progress1 % NBUF];
			if (j >= i && kinship_fmt_val[2])
				ks.ColumnFormatLine(i, j);

			PROGRESS_VALUE++;

			GUARD_END2
		}
	}
}

/* Write matrix format kinship coefficient results in a guard thread */
THREAD2(KinshipGuard2)
{
	int64& ii = progress2 = 0;
	int n = cpop<REAL>->nind;

	if (kinship_fmt_val[1])
		for (int k = 1; k <= N_KINSHIP_ESTIMATOR; ++k)
			KINSHIP<REAL>::MatrixFormatHeader(k, n);

	for (int i = 0; i < n; ++i)
	{
		if (kinship_fmt_val[1])
			for (int k = 1; k <= N_KINSHIP_ESTIMATOR; ++k)
				KINSHIP<REAL>::MatrixFormatRowHeader(k, i);

		for (int j = 0; j < n; ++j, ++progress2)
		{
			GUARD_BEGIN2

			KINSHIP<REAL>& ks = kinship_buf<REAL>[progress2 % NBUF];
			if (kinship_fmt_val[1])
				for (int k = 1; k <= N_KINSHIP_ESTIMATOR; ++k)
					ks.MatrixFormatCell(k);

			PROGRESS_VALUE++;

			GUARD_END2
		}
	}
}

/* Calculate kinship coefficient using multiple threads */
THREAD2(KinshipThread)
{
	int64 ii = 0;
	int ni = cpop<REAL>->nind;

	for (int i = 0; i < ni; ++i)
	{
		for (int j = 0; j < ni; ++j, ++ii)
		{
			THREAD_BEGIN2

			kinship_buf<REAL>[ii % NBUF].CalcKinship(ainds<REAL>[i], ainds<REAL>[j]);

			THREAD_END2
		}
	}
}

/* Calculate the individual kinship coefficient */
template<typename REAL>
TARGET void IND<REAL>::Theta(POP<REAL>* grp, double& f_ritland, double& f_loiselle, double& f_weir, double& f_vanraden, double& t_ritland, double& t_loiselle, double& t_weir, double& t_vanraden, int64 loc)
{
	double sw = 0, sr = 0, sr2 = 0, sw2 = 0, sr3 = 0, sw3 = 0, sr4 = 0, sw4 = 0;
	int ploidy = -1, N = 0;
	bool varploidy = false;
	int64 st = loc == (int64)-1 ? 0 : loc;
	int64 ed = loc == (int64)-1 ? nloc : loc + 1;

	for (int64 l = st; l < ed; ++l)
	{
		GENOTYPE& gt = GetGenotype(l);//fine
		if (gt.Nalleles() == 0) continue;
		REAL* p = grp->GetFreq(l);
		LOCSTAT1& stat1 = grp->loc_stat1[l];
		int k2 = stat1.k;
		if (k2 <= 1) continue;

		int v = gt.Ploidy();
		if (ploidy == -1) ploidy = v;
		if (ploidy != v) varploidy = true;

		int k = GetLoc(l).k;
		double tsw2 = 0, tsw3 = 1, t1 = 0, t2 = 0, t3 = 0;

		for (int i = 0; i < k; ++i)
		{
			double pr = p[i];
			if (p[i] * stat1.nhaplo <= 1e-5) continue;

			double af = gt.GetFreq<REAL>(i);

			t1 += af * af / pr;

			t2 += (af - pr) * (af - pr);
			tsw2 += pr * (1 - pr);

			t3 += af * af - pr * pr;
			tsw3 -= pr * pr;
		}

		//Ritland
		sr += t1 - 1;
		sw += k2 - 1;

		//Loiselle
		sr2 += t2;
		sw2 += tsw2;

		//Weir
		sr3 += t3;
		sw3 += tsw3;

		int i = std::max_element(p, p + k) - p;
		double pp = p[i];
		sr4 += (gt.GetFreq<REAL>(i) - pp) * (gt.GetFreq<REAL>(i) - pp);
		sw4 += pp - pp * pp;

		double t_ritland1 = (t1 - 1) / (k2 - 1);
		double t_loiselle1 = t2 / tsw2;
		double t_weir1 = t3 / tsw3;

		t_ritland1 = (v * t_ritland1 - 1) / (v - 1);
		t_loiselle1 = (v * t_loiselle1 - 1) / (v - 1);
		t_weir1 = (v * t_weir1 - 1) / (v - 1);

		f_ritland += t_ritland1;
		f_loiselle += t_loiselle1;
		f_weir += t_weir1;

		N++;
	}
	t_ritland = sr / sw;
	t_loiselle = sr2 / sw2;
	t_weir = sr3 / sw3;
	t_vanraden = sr4 / sw4;

	if (!varploidy)
	{
		f_ritland = ploidy > 1 ? (ploidy * t_ritland - 1) / (ploidy - 1) : NAN;
		f_loiselle = ploidy > 1 ? (ploidy * t_loiselle - 1) / (ploidy - 1) : NAN;
		f_weir = ploidy > 1 ? (ploidy * t_weir - 1) / (ploidy - 1) : NAN;
		f_vanraden = ploidy > 1 ? (ploidy * t_vanraden - 1) / (ploidy - 1) : NAN;
	}
	else
	{
		f_ritland = N > 0 ? f_ritland / N : NAN;
		f_loiselle = N > 0 ? f_loiselle / N : NAN;
		f_weir = N > 0 ? f_weir / N : NAN;
		f_vanraden = N > 0 ? f_vanraden / N : NAN;
	}
}
