/* Sliding Window Functions */

#pragma once
#include "vcfpop.h"
#define slide_table

template struct WINDOW<double>;
template struct WINDOW<float >;
template struct SWINDOW<double>;
template struct SWINDOW<float >;

template TARGET void WINDOW<double>::InitWindow();
template TARGET void WINDOW<float >::InitWindow();
template TARGET void WINDOW<double>::UnInitWindow();
template TARGET void WINDOW<float >::UnInitWindow();
template TARGET void WINDOW<double>::Write();
template TARGET void WINDOW<float >::Write();
template TARGET void WINDOW<double>::GetWindowId(int64 l, int& st, int& ed);
template TARGET void WINDOW<float >::GetWindowId(int64 l, int& st, int& ed);

template TARGET void SWINDOW<double>::CalcLocus(int64 l, double* buf);
template TARGET void SWINDOW<float >::CalcLocus(int64 l, double* buf);
//template TARGET void SWINDOW<double>::CalcLociPair(int64 l1, int nhaplo1, double* buf);
//template TARGET void SWINDOW<float >::CalcLociPair(int64 l1, int nhaplo1, double* buf);
template TARGET void SWINDOW<double>::InitSWindow();
template TARGET void SWINDOW<float >::InitSWindow();
template TARGET void SWINDOW<double>::UnInitSWindow();
template TARGET void SWINDOW<float >::UnInitSWindow();
template TARGET void SWINDOW<double>::Settle1();
template TARGET void SWINDOW<float >::Settle1();
//template TARGET void SWINDOW<double>::Settle2(int pst2, int st2);
//template TARGET void SWINDOW<float >::Settle2(int pst2, int st2);
template TARGET double SWINDOW<double>::SettleTajimaD(TABLE<HASH, int>& freq1, map<int, int>& freq2);
template TARGET double SWINDOW<float >::SettleTajimaD(TABLE<HASH, int>& freq1, map<int, int>& freq2);
template TARGET bool SWINDOW<double>::Settle(int i, int newid);
template TARGET bool SWINDOW<float >::Settle(int i, int newid);
template TARGET void SWINDOW<double>::GetWindowId(int64 l, int& st, int& ed);
template TARGET void SWINDOW<float >::GetWindowId(int64 l, int& st, int& ed);

template TARGET void CalcSlide<double>();
template TARGET void CalcSlide<float >();

#define extern 
template<> extern WINDOW<double> window<double>;
template<> extern WINDOW<float > window<float >;
extern vector<char*> chroms;
extern map<HASH, CHROM_PROP> chrom_sted;
#undef extern 

#ifndef _WINDOW
/* Initialize all windows */
template<typename REAL>
TARGET void WINDOW<REAL>::InitWindow()
{
	// assign cpop
	if (slide_pop_b && "total" == slide_pop_val || !slide_pop_b)
		cpop = total_pop;
	else if (slide_pop_b)
	{
		bool find = false;
		for (int i = 0; !find && i < npop; ++i)
			if (pop<REAL>[i].name == slide_pop_val)
			{
				find = true;
				cpop = &pop<REAL>[i];
			}

		if (!find)
		{
			for (uint rl = 0; !find && rl < reg<REAL>.size - 1; ++rl)
				for (uint i = 0; !find && i < reg<REAL>[rl].size; ++i)
					if (reg<REAL>[rl][i].name == slide_pop_val)
					{
						find = true;
						cpop = &reg<REAL>[rl][i];
					}
		}

		if (!find) Exit("\nError: Cannot find target population %d, check parameter -slide_pop.\n", slide_pop_val.c_str());
	}

	// assign grps
	{
		ngrps = 0;
		for (int i = 0; i < npop; ++i)
			if (apops[i]->nind > 0 && cpop->IsSubpop(apops[i]))
				ngrps++;

		grps = new POP<REAL>*[ngrps];
		int nc = 0;
		for (int i = 0; i < npop; ++i)
			if (apops[i]->nind > 0 && cpop->IsSubpop(apops[i]))
				grps[nc++] = apops[i];
	}

	//////////////////////////////////////////////////////////////////////

	//Calculate freq for cpop and diversity estimation
	RunThreads(&SlideFreqThread<REAL>, NULL, NULL, nloc * (int64)nind * (ngrps ? 2 : 1), nloc * (int64)nind * (ngrps ? 2 : 1),
		"\nPreparing allele frequency:\n", 1, true);

	RunThreads(&SlidePrepare<REAL>, NULL, NULL, nloc, nloc,
		"\nPreparing individual ploidy:\n", 1, true);

	//////////////////////////////////////////////////////////////////////

	int c = 1;

	cNei1973 = slide_estimator_val[c++] == 1;
	cWeir1984 = slide_estimator_val[c++] == 1;
	cHudson1992 = slide_estimator_val[c++] == 1;
	cHedrick2005 = slide_estimator_val[c++] == 1;
	cJost2008 = slide_estimator_val[c++] == 1;
	cHuang2021_aneu = slide_estimator_val[c++] == 1;
	cdxy = slide_estimator_val[c++] == 1;
	cpi = slide_estimator_val[c++] == 1;
	cthetaw = slide_estimator_val[c++] == 1;
	ctajimad = slide_estimator_val[c++] == 1;
	cr2 = slide_estimator_val[c++] == 1;
	cdprime = slide_estimator_val[c++] == 1;
	cr2delta = slide_estimator_val[c++] == 1;
	cdeltaprime = slide_estimator_val[c++] == 1;
	cfis = slide_estimator_val[c++] == 1;
	cho = slide_estimator_val[c++] == 1;
	che = slide_estimator_val[c++] == 1;
	cpic = slide_estimator_val[c++] == 1;
	cae = slide_estimator_val[c++] == 1;
	cI = slide_estimator_val[c++] == 1;
	clocipair = ctajimad || cr2 || cdprime || cr2delta || cdeltaprime;

	window_size = slide_windowsize_val;
	window_step = slide_windowstep_val;

	int nhaplo = cpop->nhaplotypes;
	if (nhaplo == 0) Exit("\nError: no individuals are genotyped.\n");

	slide_a1 = new double[nhaplo + 1];
	slide_a2 = new double[nhaplo + 1];
	slide_a1[0] = slide_a1[1] = slide_a2[0] = slide_a2[1] = 0;

	for (int i = 2; i <= nhaplo; ++i)
	{
		slide_a1[i] = slide_a1[i - 1] + 1.0 / i;
		slide_a2[i] = slide_a2[i - 1] + 1.0 / (i * i);
	}

	CHROM_PROP def_prop { 0, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF };

	// add chromosomes
	for (int64 l = 0; l < nloc; ++l)
	{
		char* chr = GetLoc(l).GetChrom();
		HASH ha = HashString(chr, (int)strlen(chr));
		
		if (chrom_sted.find(ha) == chrom_sted.end())
		{
			chrom_sted[ha] = def_prop;
			chroms.push_back(chr);
		}

		CHROM_PROP& prop = chrom_sted[ha];

		prop.min = Min(GetLocPos(l), prop.min);
		prop.max = Max(GetLocPos(l), prop.max);
	}

	// calculate st and ed window for each chromosome
	window_count = 0;
	for (int i = 0; i < chroms.size(); ++i)
	{
		char* chr = chroms[i];
		HASH ha = HashString(chr, (int)strlen(chr));
		CHROM_PROP& prop = chrom_sted[ha];

		int nwindow = Max(0, 1 + (prop.max - window_size + window_step) / window_step);
		prop.st = window_count;
		window_count += nwindow;
		prop.ed = window_count - 1;
	}

	swindow_count = (window_size + window_step - 1) / window_step;

#define NEW(x)  (x) = new REAL[window_count]; SetZero(x, window_count)

	// Allocate Memory
	K = new int[window_count]; SetZero(K, window_count);
	C = new int[window_count]; SetZero(C, window_count);
	lock = new LOCK[window_count]; SetZero(C, window_count);
	for (int i = 0; i < window_count; ++i) InitLock(lock[i]);

	// 1 Fst
	if (cNei1973)
	{
		NEW(Fst_Nei1973_sum1);
		NEW(Fst_Nei1973_sum2);
	}

	if (cWeir1984)
	{
		NEW(Fst_Weir1984_sum1);
		NEW(Fst_Weir1984_sum2);
	}

	if (cHudson1992)
	{
		NEW(Fst_Hudson1992_sum1);
		NEW(Fst_Hudson1992_sum2);
	}

	if (cHedrick2005)
	{
		NEW(Fst_Hedrick2005_sum1);
		NEW(Fst_Hedrick2005_sum2);
	}

	if (cJost2008)
	{
		NEW(Fst_Jost2008_sum1);
		NEW(Fst_Jost2008_sum2);
	}

	if (cHuang2021_aneu)
	{
		NEW(Fst_Huang2021_aneu_sum1);
		NEW(Fst_Huang2021_aneu_sum2);
	}

	// 2 Absolute divergence
	if (cdxy)
	{
		NEW(dxy_sum);
	}

	// 3 Nucleotide diversity
	if (cpi)
	{
		NEW(pi_sum);
	}

	// 4 Watterson¡¯s thetaw
	if (cthetaw)
	{
		NEW(thetaw_sum);
	}

	// 5 Tajima¡¯s D
	if (ctajimad)
	{
		NEW(d_sum);
		NEW(vd_sum);

		N = new int[window_count];	SetZero(N, window_count);

		int count = 0, pst = -1, ped = -1;
		for (int64 l = 0; l < nloc; ++l)
		{
			int st = -1, ed = -1;
			GetWindowId(l, st, ed);
			if (st != pst || ed != ped)
			{
				if (pst != -1)
					Add(N + pst, count, ped - pst + 1);
				count = 0;
				pst = st;
				ped = ed;
			}
			count++;
		}
		if (pst != -1)
			Add(N + pst, count, ped - pst + 1);

#ifdef slide_table
		nhaplo_freq1 = new TABLE<HASH, int>[window_count];
		SetZero(nhaplo_freq1, window_count);
#else
		nhaplo_freq2 = new map  <int,  int>[window_count];
		SetZero(nhaplo_freq2, window_count);
#endif
	}
	
	// 6 r2
	if (cr2)
	{
		NEW(D2_sum);
		NEW(Q_sum);
	}

	// 7 D'
	if (cdprime)
	{
		NEW(D_abs_sum);
		NEW(Dmax_abs_sum);
	}

	// 8 r2Delta
	if (cr2delta)
	{
		NEW(Delta2_sum);
		NEW(R_sum);
	}

	// 9 Delta'
	if (cdeltaprime)
	{
		NEW(Delta_abs_sum);
		NEW(Deltamax_abs_sum);
	}

	// 10 Fis
	if (cfis)
	{
		NEW(Fis_Nei1973_sum1);
		NEW(Fis_Nei1973_sum2);
	}

	// 11 Ho
	if (cho)
	{
		NEW(ho_sum);
	}

	// 12 He
	if (cho)
	{
		NEW(he_sum);
	}

	// 13 pic
	if (cpic)
	{
		NEW(pic_sum);
	}

	// 14 Ae
	if (cae)
	{
		NEW(ae_sum);
	}

	// 15 I
	if (cI)
	{
		NEW(I_sum);
	}
#undef NEW
}

/* Initialize all windows */
template<typename REAL>
TARGET void WINDOW<REAL>::UnInitWindow()
{
#define DEL(x) delete[] (x); (x) = NULL

	vector<char*>().swap(chroms);
	map<HASH, CHROM_PROP>().swap(chrom_sted);

	cpop->UnAllocFreq();
	for (int i = 0; i < ngrps; ++i)
		grps[i]->UnAllocFreq();

	delete[] allele_freq_offset;
	delete[] genotype_count_offset;

	allele_freq_offset = genotype_count_offset = NULL;

	delete[] slide_a1;
	delete[] slide_a2;
	delete[] grps;
	ngrps = 0;

	for (int i = 0; i < window_count; ++i) UnInitLock(lock[i]);
	delete[] K;
	delete[] C;
	delete[] lock;

	// 1 Fst
	if (cNei1973)
	{
		DEL(Fst_Nei1973_sum1);
		DEL(Fst_Nei1973_sum2);
	}

	if (cWeir1984)
	{
		DEL(Fst_Weir1984_sum1);
		DEL(Fst_Weir1984_sum2);
	}

	if (cHudson1992)
	{
		DEL(Fst_Hudson1992_sum1);
		DEL(Fst_Hudson1992_sum2);
	}

	if (cHedrick2005)
	{
		DEL(Fst_Hedrick2005_sum1);
		DEL(Fst_Hedrick2005_sum2);
	}

	if (cJost2008)
	{
		DEL(Fst_Jost2008_sum1);
		DEL(Fst_Jost2008_sum2);
	}

	if (cHuang2021_aneu)
	{
		DEL(Fst_Huang2021_aneu_sum1);
		DEL(Fst_Huang2021_aneu_sum2);
	}

	// 2 Absolute divergence
	if (cdxy)
	{
		DEL(dxy_sum);
	}

	// 3 Nucleotide diversity
	if (cpi)
	{
		DEL(pi_sum);
	}

	// 4 Watterson¡¯s thetaw
	if (cthetaw)
	{
		DEL(thetaw_sum);
	}

	// 5 Tajima¡¯s D
	if (ctajimad)
	{
		DEL(d_sum);
		DEL(vd_sum);
		DEL(N);
#ifdef slide_table
		delete[] nhaplo_freq1;
#else
		delete[] nhaplo_freq2;
#endif
	}

	// 6 r2
	if (cr2)
	{
		DEL(D2_sum);
		DEL(Q_sum);
	}

	// 7 D'
	if (cdprime)
	{
		DEL(D_abs_sum);
		DEL(Dmax_abs_sum);
	}

	// 8 r2Delta
	if (cr2delta)
	{
		DEL(Delta2_sum);
		DEL(R_sum);
	}

	// 9 Delta'
	if (cdeltaprime)
	{
		DEL(Delta_abs_sum);
		DEL(Deltamax_abs_sum);
	}

	// 10 Fis
	if (cfis)
	{
		DEL(Fis_Nei1973_sum1);
		DEL(Fis_Nei1973_sum2);
	}

	// 11 Ho
	if (cho)
	{
		DEL(ho_sum);
	}

	// 12 He
	if (che)
	{
		DEL(he_sum);
	}

	// 13 pic
	if (cpic)
	{
		DEL(pic_sum);
	}

	// 14 Ae
	if (cae)
	{
		DEL(ae_sum);
	}

	// 15 I
	if (cI)
	{
		DEL(I_sum);
	}

#undef DEL
}

/* Write result file */
template<typename REAL>
TARGET void WINDOW<REAL>::Write()
{
	//print header

	fprintf(FRES, "%s%sChrom", g_linebreak_val, g_linebreak_val);
	fprintf(FRES, "%cst", g_delimiter_val);
	fprintf(FRES, "%ced", g_delimiter_val);
	fprintf(FRES, "%c#variants", g_delimiter_val);

	// 1 Fst
	if (cNei1973)
		fprintf(FRES, "%cNei1973", g_delimiter_val);

	if (cWeir1984)
		fprintf(FRES, "%cWeir1984", g_delimiter_val);

	if (cHudson1992)
		fprintf(FRES, "%cHudson1992", g_delimiter_val);

	if (cHedrick2005)
		fprintf(FRES, "%cHedrick2005", g_delimiter_val);

	if (cJost2008)
		fprintf(FRES, "%cJost2008", g_delimiter_val);

	if (cHuang2021_aneu)
		fprintf(FRES, "%cHuang2021_aneu", g_delimiter_val);

	// 2 Absolute divergence
	if (cdxy)
		fprintf(FRES, "%cdxy", g_delimiter_val);

	// 3 Nucleotide diversity
	if (cpi)
		fprintf(FRES, "%cpi", g_delimiter_val);

	// 4 Watterson¡¯s thetaw
	if (cthetaw)
		fprintf(FRES, "%cthetaw", g_delimiter_val);

	// 5 Tajima¡¯s D
	if (ctajimad)
		fprintf(FRES, "%cTajimaD", g_delimiter_val);

	// 6 r2
	if (cr2)
		fprintf(FRES, "%cr2", g_delimiter_val);

	// 7 D'
	if (cdprime)
		fprintf(FRES, "%cD'", g_delimiter_val);

	// 8 r2Delta
	if (cr2delta)
		fprintf(FRES, "%cr2Delta", g_delimiter_val);

	// 9 Delta'
	if (cdeltaprime)
		fprintf(FRES, "%cDelta'", g_delimiter_val);

	// 10 Fis
	if (cfis)
		fprintf(FRES, "%cFis", g_delimiter_val);

	// 11 Ho
	if (cho)
		fprintf(FRES, "%cHo", g_delimiter_val);

	// 12 He
	if (che)
		fprintf(FRES, "%cHe", g_delimiter_val);

	// 13 pic
	if (cpic)
		fprintf(FRES, "%cPIC", g_delimiter_val);

	// 14 Ae
	if (cae)
		fprintf(FRES, "%cAe", g_delimiter_val);

	// 15 I
	if (cI)
		fprintf(FRES, "%cI", g_delimiter_val);

	fprintf(FRES, "%s", g_linebreak_val);

	for (char* chr : chroms)
	{
		HASH ha = HashString(chr, (int)strlen(chr));
		CHROM_PROP& prop = chrom_sted[ha];

		for (int i = prop.st; i <= prop.ed; ++i)
		{
			if (K[i] < slide_minvariants_val) continue;

			int64 pos_st = (i - prop.st) * (int64)window_step;
			int64 pos_ed = pos_st + window_size - 1;


			fprintf(FRES, "%s", chr);
			fprintf(FRES, "%c%lld", g_delimiter_val, pos_st + 1);
			fprintf(FRES, "%c%lld", g_delimiter_val, pos_ed + 1);
			fprintf(FRES, "%c%d",   g_delimiter_val, K[i]);

			double invK = 1.0 / K[i];

			// 1 Fst
			if (cNei1973)
			{
				fprintf(FRES, "%c", g_delimiter_val);
				WriteReal(FRES, Fst_Nei1973_sum1[i] / Fst_Nei1973_sum2[i]);
			}

			if (cWeir1984)
			{
				fprintf(FRES, "%c", g_delimiter_val);
				WriteReal(FRES, Fst_Weir1984_sum1[i] / Fst_Weir1984_sum2[i]);
			}

			if (cHudson1992)
			{
				fprintf(FRES, "%c", g_delimiter_val);
				WriteReal(FRES, Fst_Hudson1992_sum1[i] / Fst_Hudson1992_sum2[i]);
			}

			if (cHedrick2005)
			{
				fprintf(FRES, "%c", g_delimiter_val);
				WriteReal(FRES, Fst_Hedrick2005_sum1[i] / Fst_Hedrick2005_sum2[i]);
			}

			if (cJost2008)
			{
				fprintf(FRES, "%c", g_delimiter_val);
				WriteReal(FRES, Fst_Jost2008_sum1[i] / Fst_Jost2008_sum2[i]);
			}

			if (cHuang2021_aneu)
			{
				fprintf(FRES, "%c", g_delimiter_val);
				WriteReal(FRES, Fst_Huang2021_aneu_sum1[i] / Fst_Huang2021_aneu_sum2[i]);
			}

			// 2 Absolute divergence
			if (cdxy)
			{
				fprintf(FRES, "%c", g_delimiter_val);
				WriteReal(FRES, dxy_sum[i]);
			}

			// 3 Nucleotide diversity
			if (cpi)
			{
				fprintf(FRES, "%c", g_delimiter_val);
				WriteReal(FRES, pi_sum[i]);
			}

			// 4 Watterson¡¯s thetaw
			if (cthetaw)
			{
				fprintf(FRES, "%c", g_delimiter_val);
				WriteReal(FRES, thetaw_sum[i]);
			}

			// 5 Tajima¡¯s D
			if (ctajimad)
			{
				fprintf(FRES, "%c", g_delimiter_val);
				WriteReal(FRES, d_sum[i] / sqrt(vd_sum[i]));
			}

			// 6 r2
			if (cr2)
			{
				fprintf(FRES, "%c", g_delimiter_val);
				WriteReal(FRES, D2_sum[i] / Q_sum[i]);
			}

			// 7 D'
			if (cdprime)
			{
				fprintf(FRES, "%c", g_delimiter_val);
				WriteReal(FRES, D_abs_sum[i] / Dmax_abs_sum[i]);
			}

			// 8 r2Delta
			if (cr2delta)
			{
				fprintf(FRES, "%c", g_delimiter_val);
				WriteReal(FRES, Delta2_sum[i] / R_sum[i]);
			}

			// 9 Delta'
			if (cdeltaprime)
			{
				fprintf(FRES, "%c", g_delimiter_val);
				WriteReal(FRES, Delta_abs_sum[i] / Deltamax_abs_sum[i]);
			}

			// 10 Fis
			if (cfis)
			{
				fprintf(FRES, "%c", g_delimiter_val);
				WriteReal(FRES, Fis_Nei1973_sum1[i] / Fis_Nei1973_sum2[i]);
			}

			// 11 Ho
			if (cho)
			{
				fprintf(FRES, "%c", g_delimiter_val);
				WriteReal(FRES, ho_sum[i] * invK);
			}

			// 12 He
			if (che)
			{
				fprintf(FRES, "%c", g_delimiter_val);
				WriteReal(FRES, he_sum[i] * invK);
			}

			// 13 pic
			if (cpic)
			{
				fprintf(FRES, "%c", g_delimiter_val);
				WriteReal(FRES, pic_sum[i] * invK);
			}

			// 14 Ae
			if (cae)
			{
				fprintf(FRES, "%c", g_delimiter_val);
				WriteReal(FRES, ae_sum[i] * invK);
			}

			// 15 I
			if (cI)
			{
				fprintf(FRES, "%c", g_delimiter_val);
				WriteReal(FRES, I_sum[i] * invK);
			}

			fprintf(FRES, "%s", g_linebreak_val);
		}
	}
}

/* Get sliding window index range*/
template<typename REAL>
TARGET void WINDOW<REAL>::GetWindowId(int64 l, int& st, int& ed)
{
	char* chr = GetLoc(l).GetChrom();
	HASH ha = HashString(chr, (int)strlen(chr));
	CHROM_PROP& prop = chrom_sted[ha];

	if (prop.st == -1)
	{
		st = ed = -1;
		return;
	}

	st = prop.st + Max((int)GetLocPos(l) - (int)window_size + (int)window_step - 1, 0) / window_step;
	ed = Min(prop.st + ((int)GetLocPos(l) - 1) / (int)window_step, prop.ed);
}
#endif

#ifndef _SWINDOW

/* Calculate diveristy indices */
template<typename REAL>
TARGET void SWINDOW<REAL>::CalcLocus(int64 l, double* buf)
{
	int st = -1, ed = -1;
	GetWindowId(l, st, ed);

	if (st == -1) return;

	if (st != pst || ed != ped)
	{
		if (_C > 0) Settle1();

		for (int i = st; i <= ed; ++i)
			if (!Settle(i % swindow_count, i))
				break;

		for (int i = ed; i >= st; --i)
			if (!Settle(i % swindow_count, i))
				break;

		pst = st;
		ped = ed;
	}

//////////////////////////////////////////////////////////////////////////////////

	_C++;
	REAL* fre = cpop->GetFreq(l);
	ushort* gcount = cpop->GetGenoCount(l);

	double ho = 0, how = 0;
	double he = 0, ae = 0, I = 0;
	double a2 = 0, pic = 0;
	int nhaplo = 0, k1 = GetLoc(l).k, n = 0;

	GENOTYPE* gtab = GetLoc(l).GetGtab();
	int ngeno = GetLoc(l).ngeno;

	for (int gi = 0; gi < ngeno; ++gi)
	{
		GENOTYPE& gt = gtab[gi];
		if (gt.Nalleles() == 0 || gcount[gi] == 0) continue;
		uint c = gcount[gi];
		int v = gt.Ploidy();
		nhaplo += v * c;

		ho += gt.HIndex() * c * v * (v - 1);
		how +=              c * v * (v - 1);

		n += c;
	}

	int k = n == 0 ? 0 : (ushort)CountK(fre, k1);

	//set statistics to nan if k <= 2 to avoid involved in further average
	
	if (k < 2) return;

	ho = ho / how;
	he = 1;
	pic = 0;
	ae = 0;
	I = 0;

	for (int a = 0; a < k1; ++a)
	{
		REAL af = fre[a];
		if (af * nhaplo < 1e-5) continue;
		if (af > 1e-5) I += -af * log(af);
		a2 += af * af;

		for (int b = a + 1; b < k1; ++b)
			if (fre[b] > 0)
				pic += 2 * af * fre[b] * (1 - af * fre[b]);
	}

	he = he - a2;
	ae = 1 / a2;

	//////////////////////////////////////////////////////////////////////////////////
	
	_K++;

	// 1 Fst
	if (ngrps >= 2)
	{
		if (cNei1973)			FST<REAL>::Fst_Nei1973(grps, ngrps, NULL, buf, l, &_Fst_Nei1973_sum1, &_Fst_Nei1973_sum2); 
		if (cWeir1984)			FST<REAL>::Fst_Weir1984(grps, ngrps, NULL, l, &_Fst_Weir1984_sum1, &_Fst_Weir1984_sum2); 
		if (cHudson1992)		FST<REAL>::Fst_Hudson1992(grps, ngrps, NULL, buf, l, &_Fst_Hudson1992_sum1, &_Fst_Hudson1992_sum2); 
		if (cHedrick2005)		FST<REAL>::Fst_Hedrick2005(grps, ngrps, NULL, buf, l, &_Fst_Hedrick2005_sum1, &_Fst_Hedrick2005_sum2); 
		if (cJost2008)			FST<REAL>::Fst_Jost2008(grps, ngrps, NULL, buf, l, &_Fst_Jost2008_sum1, &_Fst_Jost2008_sum2); 
		if (cHuang2021_aneu)	FST<REAL>::Fst_Huang2021_aneu(grps, ngrps, 2, true, true, NULL, buf, l, &_Fst_Huang2021_aneu_sum1, &_Fst_Huang2021_aneu_sum2); 
	}

	// 2 Absolute divergence
	if (cdxy)
		_dxy_sum += FST<REAL>::dxy(grps, ngrps, buf, l);

	// 3 Nucleotide diversity
	if (cpi)
		_pi_sum += nhaplo / (nhaplo - 1) * he;

	// 4 Watterson¡¯s thetaw
	if (cthetaw)
		_thetaw_sum += 1.0 / slide_a1[nhaplo];

	// 5 Tajima¡¯s D
	if (ctajimad)
	{
		_d_sum += nhaplo / (nhaplo - 1) * he - 1.0 / slide_a1[nhaplo];

#ifdef slide_table
		if (!_nhaplo_freq1.ContainsKey(nhaplo))
			_nhaplo_freq1[nhaplo] = 1;
		else
			_nhaplo_freq1[nhaplo] ++;
#else
		if (_nhaplo_freq2.find(nhaplo) == _nhaplo_freq2.end())
			_nhaplo_freq2[nhaplo] = 1;
		else
			_nhaplo_freq2[nhaplo] ++;
#endif
		//_vard1_sum += ((nhaplo + 1) / (3 * (nhaplo - 1)) - 1.0 / slide_a1[nhaplo]) / slide_a1[nhaplo];
	}

	// 5 Tajima¡¯s D  6 r2  7 D'  8 r2Delta  9 Delta'
	//if (clocipair)
		//CalcLociPair(l, nhaplo, buf);

	// 10 Fis
	if (cfis)
	{
		_Fis_Nei1973_sum1 += he - ho;
		_Fis_Nei1973_sum2 += he;
	}

	// 11 Ho
	if (cho)
		_ho_sum += ho;

	// 12 He
	if (che)
		_he_sum += he;

	// 13 pic
	if (cpic)
		_pic_sum += pic;

	// 14 Ae
	if (cae)
		_ae_sum += ae;

	// 15 I
	if (cI)
		_I_sum += I;
}

/* Calculate Tajima's D, r2, D', r2Delta and Delta'
template<typename REAL>
TARGET void SWINDOW<REAL>::CalcLociPair(int64 l1, int nhaplo1, double* buf)
{
	double a1A = slide_a1[nhaplo1];
	double b1A = (nhaplo1 + 1) / (3 * (nhaplo1 - 1));
	double c1A = b1A - 1 / a1A;
	double a2A = slide_a2[nhaplo1];
	double b2A = (2 * (nhaplo1 * nhaplo1 + nhaplo1 + 3)) / (9 * nhaplo1 * (nhaplo1 - 1));
	double c2A = b2A - (nhaplo1 + 2) / (a1A * nhaplo1) + a2A / (a1A * a1A);

	// Linkage Disequilibrium
	int* mA  = (int*)buf;
	int* mB  = mA + maxK;
	int* mAB = mB + maxK;

	int k1 = GetLoc(l1).k;

	REAL* fre1 = cpop->GetFreq(l1);
	ushort* gcount1 = cpop->GetGenoCount(l1);
	LOCSTAT1& stat11 = cpop->loc_stat1[l1];
	GENOTYPE* gtab1 = GetLoc(l1).GetGtab();
	int ngeno1 = GetLoc(l1).ngeno;

	// Linkage Disequilibrium
	int* nAA = (int*)buf;	//maxK
	int* nBB = nAA + maxK;	//maxK
	int* nAX = nBB + maxK;	//maxK
	int* nBX = nAX + maxK;	//maxK
	int* nA  = nBX + maxK;	//maxK
	int* nB  = nA  + maxK;	//maxK
	int* nAB = nB  + maxK;	//maxK * maxK

	int pst2 = -1, ped2 = -1;

	for (int64 l2 = l1 + 1; l2 < nloc; ++l2)
	{
		int st2 = -1, ed2 = -1;
		GetWindowId(l2, st2, ed2);
		if (st2 == -1 || st2 > ped) break;
		if (pst2 == -1) { pst2 = st2; ped2 = ped2; }
		if (pst2 != st2) Settle2(pst2, st2);
		pst2 = st2;

		int k2 = GetLoc(l2).k;
		REAL* fre2 = cpop->GetFreq(l2);
		ushort* gcount2 = cpop->GetGenoCount(l2);
		LOCSTAT1& stat12 = cpop->loc_stat1[l2];
		GENOTYPE* gtab2 = GetLoc(l2).GetGtab();
		int ngeno2 = GetLoc(l2).ngeno;

		////////////////////////////////////////////////////////////
		if (ctajimad)
		{
			int nhaplo2 = 0;
			for (int gi2 = 0; gi2 < ngeno2; ++gi2)
			{
				GENOTYPE& gt2 = gtab2[gi2];
				if (gt2.Nalleles() == 0 || gcount2[gi2] == 0) continue;
				nhaplo2 += gt2.Ploidy() * gcount2[gi2];
			}

			double a1B = slide_a1[nhaplo2];
			double a2B = slide_a2[nhaplo2];
			double b2B = (2 * (nhaplo2 * nhaplo2 + nhaplo2 + 3)) / (9 * nhaplo2 * (nhaplo2 - 1));
			double c2B = b2B - (nhaplo2 + 2) / (a1B * nhaplo2) + a2B / (a1B * a1B);

			_vard2_sum += sqrt(c2A * c2B) / (sqrt(a2A * a2B) + a1A * a1B);
		}

		////////////////////////////////////////////////////////////
		// Linkage Disequilibrium
		if (cr2 || cdprime)
		{
			SetZero(mAB, maxK * maxK + maxK * 2);
			int nhaplo = 0;

			GENO_READER rt1(cpop->ind0id, l1), rt2(cpop->ind0id, l2);
			for (int i = 0; i < cpop->nind; ++i)
			{
				GENOTYPE& g1 = gtab1[rt1.Read()], &g2 = gtab2[rt2.Read()];
				if (g1.Nalleles() == 0 || g2.Nalleles() == 0) continue;

				int v1 = g1.Ploidy(), v2 = g2.Ploidy();
				if (v1 != v2) continue;

				nhaplo += v1;
				ushort* als1 = g1.GetAlleleArray(), * als2 = g2.GetAlleleArray();
				for (int a = 0; a < v1; ++a)
				{
					mA[als1[a]] ++;
					mB[als2[a]] ++;
					mAB[als1[a] * maxK + als2[a]] ++;
				}
			}

			double nhaploinv = 1.0 / nhaplo;

			for (int A = 0; A < k1; ++A)
			{
				if (mA[A] == 0) continue;
				for (int B = 0; B < k2; ++B)
				{
					if (mB[B] == 0) continue;
					double pA = mA[A] * nhaploinv, pB = mB[B] * nhaploinv;
					double DAB = mAB[A * maxK + B] * nhaploinv - pA * pB;
					_D2_sum += DAB * DAB;
					_Q_sum += pA * pB * (1 - pA) * (1 - pB);
					_D_abs_sum += fabs(DAB);
					_Dmax_abs_sum += (DAB >= 0 ?
						(pA < pB ? pA * (1 - pB) : pB * (1 - pA)) :
						(pA + pB <= 1 ? pA * pB : (1 - pA) * (1 - pB)));
				}
			}
		}

		////////////////////////////////////////////////////////

		if (cr2delta || cdeltaprime)
		{
			SetZero(nAB, maxK * maxK + maxK * 6);

			int sv2 = 0, nhaplo = 0, ntAABB = 0;

			GENO_READER rt1(cpop->ind0id, l1), rt2(cpop->ind0id, l2);
			for (int i = 0; i < cpop->nind; ++i)
			{
				GENOTYPE& g1 = gtab1[rt1.Read()], & g2 = gtab2[rt2.Read()];
				if (g1.Nalleles() == 0 || g2.Nalleles() == 0) continue;

				int nalleles1 = g1.Nalleles(), nalleles2 = g2.Nalleles();
				if (nalleles1 == 0 || nalleles2 == 0) continue;

				int v1 = g1.Ploidy(), v2 = g2.Ploidy();
				if (v1 != v2) continue;

				sv2 += v1 * v1;
				nhaplo += v1;
				ushort* als1 = g1.GetAlleleArray() + v1, * als2 = g2.GetAlleleArray() + v2;
				uint64 pattern1 = g1.GetPattern(), pattern2 = g2.GetPattern();

				uint64 p1 = pattern1;
				for (int a1 = 0; a1 < nalleles1; ++a1)
				{
					int ncopy1 = pattern1 & 0xF; pattern1 >>= 4;

					nA[als1[a1]] += ncopy1;
					nAA[als1[a1]] += ncopy1 * (ncopy1 - 1);

					uint64 p2 = pattern2;
					for (int a2 = 0; a2 < nalleles2; ++a2)
					{
						int ncopy2 = pattern2 & 0xF; pattern2 >>= 4;
						nAB[als1[a1] * maxK + als2[a2]] += ncopy1 * ncopy2;
					}
				}

				uint64 p2 = pattern2;
				for (int a2 = 0; a2 < nalleles2; ++a2)
				{
					int ncopy2 = pattern2 & 0xF; pattern2 >>= 4;
					nB[als2[a2]] += ncopy2;
					nBB[als2[a2]] += ncopy2 * (ncopy2 - 1);
				}
			}

			double invndpair = 1.0 / (sv2 - nhaplo);
			double invnhaplo = 1.0 / nhaplo;
			double vtilde = sv2 / (double)nhaplo;

			for (int A = 0; A < k1; ++A)
			{
				if (nA[A] == 0) continue;
				double pA = nA[A] * invnhaplo, pAA = nAA[A] * invndpair, pA2 = pA * pA, pAX = pA - pAA;

				for (int B = 0; B < k2; ++B)
				{
					if (nB[B] == 0) continue;
					double pB = nB[B] * invnhaplo, pBB = nBB[B] * invndpair, pB2 = pB * pB;
					double DeltaAB = nAB[A * maxK + B] * invnhaplo - vtilde * pA * pB;
					double lambda = 1 - pA - pB;

					_Delta2_sum += DeltaAB * DeltaAB;
					_R_sum += (pA + (vtilde - 1) * pAA - vtilde * pA2) * (pB + (vtilde - 1) * pBB - vtilde * pB2);
					_Delta_abs_sum += fabs(DeltaAB);
					_Deltamax_abs_sum += (DeltaAB >= 0 ?
						(pA < pB ?
							pA + (vtilde - 1) * (pAA + pAX * (pB - pA) / (1 - pA)) - vtilde * pA * pB :
							pB + (vtilde - 1) * (pAA * pB / pA) - vtilde * pA * pB) :
						(pA + pB <= 1 ?
							(vtilde - 1) * (pAX * pB / (1 - pA)) - vtilde * pA * pB :
							lambda + (vtilde - 1) * (pAA * lambda / pA + pAX * (pB - lambda) / (1 - pA)) - vtilde * pA * pB));
				}
			}
		}
	}

	Settle2(pst2, ped + 1);

	_vard2_sum = 0;
	_D2_sum = 0;
	_Q_sum = 0;
	_D_abs_sum = 0;
	_Dmax_abs_sum = 0;
	_Delta2_sum = 0;
	_R_sum = 0;
	_Delta_abs_sum = 0;
	_Deltamax_abs_sum = 0;
}
*/

/* Initialize a SWINDOW */
template<typename REAL>
TARGET void SWINDOW<REAL>::InitSWindow()
{
#define NEW(x) (x) = new double[swindow_count]; SetZero((x), swindow_count)

	//SetZero(this, 1);

	pst = ped = -1;
	_K = _C = 0;
	_Fst_Nei1973_sum1 = _Fst_Nei1973_sum2 = 0;
	_Fst_Weir1984_sum1 = _Fst_Weir1984_sum2 = 0;
	_Fst_Hudson1992_sum1 = _Fst_Hudson1992_sum2 = 0;
	_Fst_Hedrick2005_sum1 = _Fst_Hedrick2005_sum2 = 0;
	_Fst_Jost2008_sum1 = _Fst_Jost2008_sum2 = 0;
	_Fst_Huang2021_aneu_sum1 = _Fst_Huang2021_aneu_sum2 = 0;
	_dxy_sum = 0;
	_pi_sum = 0;
	_thetaw_sum = 0;
	_d_sum = 0; 
#ifdef slide_table
	new (&_nhaplo_freq1) TABLE<HASH, int>(true, NULL);
#endif

	_D2_sum = 0;
	_Q_sum = 0;
	_D_abs_sum = 0;
	_Dmax_abs_sum = 0;
	_Delta2_sum = 0;
	_R_sum = 0;
	_Delta_abs_sum = 0;
	_Deltamax_abs_sum = 0;
	_Fis_Nei1973_sum1 = 0;
	_Fis_Nei1973_sum2 = 0;
	_ho_sum = 0;
	_he_sum = 0;
	_pic_sum = 0;
	_ae_sum = 0;
	_I_sum = 0;

	int c = 1;

	cNei1973 = slide_estimator_val[c++] == 1;
	cWeir1984 = slide_estimator_val[c++] == 1;
	cHudson1992 = slide_estimator_val[c++] == 1;
	cHedrick2005 = slide_estimator_val[c++] == 1;
	cJost2008 = slide_estimator_val[c++] == 1;
	cHuang2021_aneu = slide_estimator_val[c++] == 1;
	cdxy = slide_estimator_val[c++] == 1;
	cpi = slide_estimator_val[c++] == 1;
	cthetaw = slide_estimator_val[c++] == 1;
	ctajimad = slide_estimator_val[c++] == 1;
	cr2 = slide_estimator_val[c++] == 1;
	cdprime = slide_estimator_val[c++] == 1;
	cr2delta = slide_estimator_val[c++] == 1;
	cdeltaprime = slide_estimator_val[c++] == 1;
	cfis = slide_estimator_val[c++] == 1;
	cho = slide_estimator_val[c++] == 1;
	che = slide_estimator_val[c++] == 1;
	cpic = slide_estimator_val[c++] == 1;
	cae = slide_estimator_val[c++] == 1;
	cI = slide_estimator_val[c++] == 1;
	clocipair = ctajimad || cr2 || cdprime || cr2delta || cdeltaprime;

	window_count = window<REAL>.window_count;
	swindow_count = window<REAL>.swindow_count;
	window_size = window<REAL>.window_size;
	window_step = window<REAL>.window_step;

	slide_a1 = window<REAL>.slide_a1;
	slide_a2 = window<REAL>.slide_a2;
	cpop_ = window<REAL>.cpop_;
	grps = window<REAL>.grps;
	ngrps = window<REAL>.ngrps;

	window_id	= new int[swindow_count];
	K			= new int[swindow_count];
	C			= new int[swindow_count];
	SetVal(window_id, -1, swindow_count);
	SetZero(K, swindow_count);
	SetZero(C, swindow_count);


	// Allocate Memory
	// 1 Fst
	if (cNei1973)
	{
		NEW(Fst_Nei1973_sum1);
		NEW(Fst_Nei1973_sum2);
	}

	if (cWeir1984)
	{
		NEW(Fst_Weir1984_sum1);
		NEW(Fst_Weir1984_sum2);
	}

	if (cHudson1992)
	{
		NEW(Fst_Hudson1992_sum1);
		NEW(Fst_Hudson1992_sum2);
	}

	if (cHedrick2005)
	{
		NEW(Fst_Hedrick2005_sum1);
		NEW(Fst_Hedrick2005_sum2);
	}

	if (cJost2008)
	{
		NEW(Fst_Jost2008_sum1);
		NEW(Fst_Jost2008_sum2);
	}

	if (cHuang2021_aneu)
	{
		NEW(Fst_Huang2021_aneu_sum1);
		NEW(Fst_Huang2021_aneu_sum2);
	}

	// 2 Absolute divergence
	if (cdxy)
	{
		NEW(dxy_sum);
	}

	// 3 Nucleotide diversity
	if (cpi)
	{
		NEW(pi_sum);
	}

	// 4 Watterson¡¯s thetaw
	if (cthetaw)
	{
		NEW(thetaw_sum);
	}

	// 5 Tajima¡¯s D
	if (ctajimad)
	{
		NEW(d_sum);
#ifdef slide_table
		nhaplo_freq1 = new TABLE<HASH, int>[swindow_count];
		for (int i = 0; i < swindow_count; ++i)
			new (&nhaplo_freq1[i]) TABLE<HASH, int>(true, NULL);
#else
		nhaplo_freq2 = new map<int, int>[swindow_count];
		for (int i = 0; i < swindow_count; ++i)
			new (&nhaplo_freq2[i]) map<int, int>();
#endif
	}

	// 6 r2
	if (cr2)
	{
		NEW(D2_sum);
		NEW(Q_sum);
	}

	// 7 D'
	if (cdprime)
	{
		NEW(D_abs_sum);
		NEW(Dmax_abs_sum);
	}

	// 8 r2Delta
	if (cr2delta)
	{
		NEW(Delta2_sum);
		NEW(R_sum);
	}

	// 9 Delta'
	if (cdeltaprime)
	{
		NEW(Delta_abs_sum);
		NEW(Deltamax_abs_sum);
	}

	// 10 Fis
	if (cfis)
	{
		NEW(Fis_Nei1973_sum1);
		NEW(Fis_Nei1973_sum2);
	}

	// 11 Ho
	if (cho)
	{
		NEW(ho_sum);
	}

	// 12 He
	if (che)
	{
		NEW(he_sum);
	}

	// 13 pic
	if (cpic)
	{
		NEW(pic_sum);
	}

	// 14 Ae
	if (cae)
	{
		NEW(ae_sum);
	}

	// 15 I
	if (cI)
	{
		NEW(I_sum);
	}
#undef NEW
}

/* Uninitialize a SWINDOW */
template<typename REAL>
TARGET void SWINDOW<REAL>::UnInitSWindow()
{
#define DEL(x) delete[] (x); (x) = NULL

	delete[] window_id;
	delete[] K;
	delete[] C;

	// 1 Fst
	if (cNei1973)
	{
		DEL(Fst_Nei1973_sum1);
		DEL(Fst_Nei1973_sum2);
	}

	if (cWeir1984)
	{
		DEL(Fst_Weir1984_sum1);
		DEL(Fst_Weir1984_sum2);
	}

	if (cHudson1992)
	{
		DEL(Fst_Hudson1992_sum1);
		DEL(Fst_Hudson1992_sum2);
	}

	if (cHedrick2005)
	{
		DEL(Fst_Hedrick2005_sum1);
		DEL(Fst_Hedrick2005_sum2);
	}

	if (cJost2008)
	{
		DEL(Fst_Jost2008_sum1);
		DEL(Fst_Jost2008_sum2);
	}

	if (cHuang2021_aneu)
	{
		DEL(Fst_Huang2021_aneu_sum1);
		DEL(Fst_Huang2021_aneu_sum2);
	}

	// 2 Absolute divergence
	if (cdxy)
	{
		DEL(dxy_sum);
	}

	// 3 Nucleotide diversity
	if (cpi)
	{
		DEL(pi_sum);
	}

	// 4 Watterson¡¯s thetaw
	if (cthetaw)
	{
		DEL(thetaw_sum);
	}

	// 5 Tajima¡¯s D
	if (ctajimad)
	{
		DEL(d_sum);
#ifdef slide_table
		for (int i = 0; i < swindow_count; ++i)
			nhaplo_freq1[i].~TABLE();
		delete[] nhaplo_freq1;
#else
		for (int i = 0; i < swindow_count; ++i)
			map<int, int>().swap(nhaplo_freq2[i]);
		delete[] nhaplo_freq2;
#endif
	}

	// 6 r2
	if (cr2)
	{
		DEL(D2_sum);
		DEL(Q_sum);
	}

	// 7 D'
	if (cdprime)
	{
		DEL(D_abs_sum);
		DEL(Dmax_abs_sum);
	}

	// 8 r2Delta
	if (cr2delta)
	{
		DEL(Delta2_sum);
		DEL(R_sum);
	}

	// 9 Delta'
	if (cdeltaprime)
	{
		DEL(Delta_abs_sum);
		DEL(Deltamax_abs_sum);
	}

	// 10 Fis
	if (cfis)
	{
		DEL(Fis_Nei1973_sum1);
		DEL(Fis_Nei1973_sum2);
	}

	// 11 Ho
	if (cho)
	{
		DEL(ho_sum);
	}

	// 12 He
	if (che)
	{
		DEL(he_sum);
	}

	// 13 pic
	if (cpic)
	{
		DEL(pic_sum);
	}

	// 14 Ae
	if (cae)
	{
		DEL(ae_sum);
	}

	// 15 I
	if (cI)
	{
		DEL(I_sum);
	}

	//SetZero(this, 1);
#undef DEL
}

/* When a new locus has a different range, distributed current results to SWINDOW */
template<typename REAL>
TARGET void SWINDOW<REAL>::Settle1()
{
	int stid = pst % swindow_count;
	int edid = ped % swindow_count;

	if (stid <= edid || ped - pst == swindow_count - 1)
	{
		// add to SWINDOW stid to edid
#define ADD(x) Add((x) + _st, _ ## x, _len); _ ## x = 0
#define ADD2(x,y) Add((x) + _st, y, _len); y = 0

		int _st = ped - pst == swindow_count - 1 ? 0 : stid;
		int _len = ped - pst == swindow_count - 1 ? swindow_count : edid + 1 - stid;

		ADD(C);
		if (_K == 0) return;
		ADD(K);

		// 1 Fst
		if (cNei1973)
		{
			ADD(Fst_Nei1973_sum1);
			ADD(Fst_Nei1973_sum2);
		}

		if (cWeir1984)
		{
			ADD(Fst_Weir1984_sum1);
			ADD(Fst_Weir1984_sum2);
		}

		if (cHudson1992)
		{
			ADD(Fst_Hudson1992_sum1);
			ADD(Fst_Hudson1992_sum2);
		}

		if (cHedrick2005)
		{
			ADD(Fst_Hedrick2005_sum1);
			ADD(Fst_Hedrick2005_sum2);
		}

		if (cJost2008)
		{
			ADD(Fst_Jost2008_sum1);
			ADD(Fst_Jost2008_sum2);
		}

		if (cHuang2021_aneu)
		{
			ADD(Fst_Huang2021_aneu_sum1);
			ADD(Fst_Huang2021_aneu_sum2);
		}

		// 2 Absolute divergence
		if (cdxy)
		{
			ADD(dxy_sum);
		}

		// 3 Nucleotide diversity
		if (cpi)
		{
			ADD(pi_sum);
		}

		// 4 Watterson¡¯s thetaw
		if (cthetaw)
		{
			ADD(thetaw_sum);
		}

		// 5 Tajima¡¯s D
		if (ctajimad)
		{
			ADD(d_sum);
#ifdef slide_table
			for (int i = 0; i < _nhaplo_freq1.size; ++i)
			{
				TABLE_ENTRY<HASH, int>& entry = _nhaplo_freq1.GetEntry(i);
				HASH nhaplo = entry.key;
				int count = entry.val;

				for (int jj = pst; jj <= ped; ++jj)
				{
					int j = jj % swindow_count;
					if (!nhaplo_freq1[j].ContainsKey(nhaplo))
						nhaplo_freq1[j][nhaplo] = count;
					else
						nhaplo_freq1[j][nhaplo] += count;
				}
			}
			_nhaplo_freq1.Clear();
#else
			for (auto entry : _nhaplo_freq2)
			{
				int nhaplo = entry.first;
				int count  = entry.second;

				for (int jj = pst; jj <= ped; ++jj)
				{
					int j = jj % swindow_count;
					if (nhaplo_freq2[j].find(nhaplo) == nhaplo_freq2[j].end())
						nhaplo_freq2[j][nhaplo] = count;
					else
						nhaplo_freq2[j][nhaplo] += count;
				}
			}
			_nhaplo_freq2.clear();
#endif
		}

		/*
		// 6 r2
		if (cr2)
		{
			ADD(D2_sum);
			ADD(Q_sum);
		}

		// 7 D'
		if (cdprime)
		{
			ADD(D_abs_sum);
			ADD(Dmax_abs_sum);
		}

		// 8 r2Delta
		if (cr2delta)
		{
			ADD(Delta2_sum);
			ADD(R_sum);
		}

		// 9 Delta'
		if (cdeltaprime)
		{
			ADD(Delta_abs_sum);
			ADD(Deltamax_abs_sum);
		}
		*/

		// 10 Fis
		if (cfis)
		{
			ADD(Fis_Nei1973_sum1);
			ADD(Fis_Nei1973_sum2);
		}

		// 11 Ho
		if (cho)
		{
			ADD(ho_sum);
		}

		// 12 He
		if (che)
		{
			ADD(he_sum);
		}

		// 13 pic
		if (cpic)
		{
			ADD(pic_sum);
		}

		// 14 Ae
		if (cae)
		{
			ADD(ae_sum);
		}

		// 15 I
		if (cI)
		{
			ADD(I_sum);
		}
#undef ADD
#undef ADD2
	}
	else
	{
	// add to SWINDOW (0 to edid) and (stid to swindow_count - 1)
#define ADD(x) Add((x), _ ## x, edid + 1); Add((x) + stid, _ ## x, _len); _ ## x = 0
#define ADD2(x,y) Add((x), y, edid + 1); Add((x) + stid, y, _len); y = 0

		int _len = swindow_count - stid;
		ADD(C);
		if (_K == 0) return;
		ADD(K);

		// 1 Fst
		if (cNei1973)
		{
			ADD(Fst_Nei1973_sum1);
			ADD(Fst_Nei1973_sum2);
		}

		if (cWeir1984)
		{
			ADD(Fst_Weir1984_sum1);
			ADD(Fst_Weir1984_sum2);
		}

		if (cHudson1992)
		{
			ADD(Fst_Hudson1992_sum1);
			ADD(Fst_Hudson1992_sum2);
		}

		if (cHedrick2005)
		{
			ADD(Fst_Hedrick2005_sum1);
			ADD(Fst_Hedrick2005_sum2);
		}

		if (cJost2008)
		{
			ADD(Fst_Jost2008_sum1);
			ADD(Fst_Jost2008_sum2);
		}

		if (cHuang2021_aneu)
		{
			ADD(Fst_Huang2021_aneu_sum1);
			ADD(Fst_Huang2021_aneu_sum2);
		}

		// 2 Absolute divergence
		if (cdxy)
		{
			ADD(dxy_sum);
		}

		// 3 Nucleotide diversity
		if (cpi)
		{
			ADD(pi_sum);
		}

		// 4 Watterson¡¯s thetaw
		if (cthetaw)
		{
			ADD(thetaw_sum);
		}

		// 5 Tajima¡¯s D
		if (ctajimad)
		{
			ADD(d_sum);
#ifdef slide_table
			for (int i = 0; i < _nhaplo_freq1.size; ++i)
			{
				TABLE_ENTRY<HASH, int>& entry = _nhaplo_freq1.GetEntry(i);
				HASH nhaplo = entry.key;
				int  count  = entry.val;

				for (int jj = pst; jj <= ped; ++jj)
				{
					int j = jj % swindow_count;
					if (!nhaplo_freq1[j].ContainsKey(nhaplo))
						nhaplo_freq1[j][nhaplo] = count;
					else
						nhaplo_freq1[j][nhaplo] += count;
				}
			}
			_nhaplo_freq1.Clear();
#else
			for (auto entry : _nhaplo_freq2)
			{
				int nhaplo = entry.first;
				int count  = entry.second;

				for (int jj = pst; jj <= ped; ++jj)
				{
					int j = jj % swindow_count;
					if (nhaplo_freq2[j].find(nhaplo) == nhaplo_freq2[j].end())
						nhaplo_freq2[j][nhaplo] = count;
					else
						nhaplo_freq2[j][nhaplo] += count;
				}
			}
			_nhaplo_freq2.clear();
#endif
		}

		/*
		// 6 r2
		if (cr2)
		{
			ADD(D2_sum);
			ADD(Q_sum);
		}

		// 7 D'
		if (cdprime)
		{
			ADD(D_abs_sum);
			ADD(Dmax_abs_sum);
		}

		// 8 r2Delta
		if (cr2delta)
		{
			ADD(Delta2_sum);
			ADD(R_sum);
		}

		// 9 Delta'
		if (cdeltaprime)
		{
			ADD(Delta_abs_sum);
			ADD(Deltamax_abs_sum);
		}
		*/

		// 10 Fis
		if (cfis)
		{
			ADD(Fis_Nei1973_sum1);
			ADD(Fis_Nei1973_sum2);
		}

		// 11 Ho
		if (cho)
		{
			ADD(ho_sum);
		}

		// 12 He
		if (che)
		{
			ADD(he_sum);
		}

		// 13 pic
		if (cpic)
		{
			ADD(pic_sum);
		}

		// 14 Ae
		if (cae)
		{
			ADD(ae_sum);
		}

		// 15 I
		if (cI)
		{
			ADD(I_sum);
		}
#undef ADD
#undef ADD2
	}
}

/* When the new loci pair has a different range, distributed current results to SWINDOW 
template<typename REAL>
TARGET void SWINDOW<REAL>::Settle2(int pst2, int st2)
{
	//from stid to edid, do not clear current results
	int stid = pst2 % swindow_count;
	int edid = (st2 - 1) % swindow_count;

	if (stid <= edid || (edid + 1) % swindow_count == stid)
	{
#define ADD(x) Add((x) + _st, _ ## x, _len)
#define ADD2(x, y) Add((x) + _st, y, _len)

		int _st =  (edid + 1) % swindow_count == stid ? 0 : stid;
		int _len = (edid + 1) % swindow_count == stid ? swindow_count : edid + 1 - stid;

		// 5 Tajima¡¯s D
		if (ctajimad)
		{
			ADD2(vard_sum, _vard2_sum);
		}

		// 6 r2
		if (cr2)
		{
			ADD(D2_sum);
			ADD(Q_sum);
		}

		// 7 D'
		if (cdprime)
		{
			ADD(D_abs_sum);
			ADD(Dmax_abs_sum);
		}

		// 8 r2Delta
		if (cr2delta)
		{
			ADD(Delta2_sum);
			ADD(R_sum);
		}

		// 9 Delta'
		if (cdeltaprime)
		{
			ADD(Delta_abs_sum);
			ADD(Deltamax_abs_sum);
		}
#undef ADD
	}
	else
	{
#define ADD(x) Add((x), _ ## x, edid + 1); Add((x) + stid, _ ## x, _len)
#define ADD2(x, y) Add((x), y, edid + 1); Add((x) + stid, y, _len)

		int _len = swindow_count - stid;

		// 5 Tajima¡¯s D
		if (ctajimad)
		{
			ADD2(vard_sum, _vard2_sum);
		}

		// 6 r2
		if (cr2)
		{
			ADD(D2_sum);
			ADD(Q_sum);
		}

		// 7 D'
		if (cdprime)
		{
			ADD(D_abs_sum);
			ADD(Dmax_abs_sum);
		}

		// 8 r2delta
		if (cr2delta)
		{
			ADD(Delta2_sum);
			ADD(R_sum);
		}

		// 9 delta'
		if (cdeltaprime)
		{
			ADD(Delta_abs_sum);
			ADD(Deltamax_abs_sum);
		}

#undef ADD
	}
}
*/

/* Settle Tajima D denominator */
template<typename REAL>
TARGET double SWINDOW<REAL>::SettleTajimaD(TABLE<HASH, int>& freq1, map<int, int>& freq2)
{
	double re = 0;

#ifdef slide_table
	for (int i = 0; i < freq1.size; ++i)
	{
		TABLE_ENTRY<HASH, int>& entry1 = freq1.GetEntry(i);
		HASH nhaplo1 = entry1.key;
		int  count1  = entry1.val;
		double a1A = slide_a1[nhaplo1];
		double b1A = (nhaplo1 + 1.0) / (3.0 * (nhaplo1 - 1.0));
		double c1A = b1A - 1.0 / a1A;
		double a2A = slide_a2[nhaplo1];
		double b2A = (2.0 * (nhaplo1 * nhaplo1 + nhaplo1 + 3.0)) / (9.0 * nhaplo1 * (nhaplo1 - 1.0));
		double c2A = b2A - (nhaplo1 + 2.0) / (a1A * nhaplo1) + a2A / (a1A * a1A);

		for (int j = i; j < freq1.size; ++j)
		{
			TABLE_ENTRY<HASH, int>& entry2 = freq1.GetEntry(j);
			HASH nhaplo2 = entry2.key;

			if (nhaplo1 == nhaplo2)
			{
				re += count1 * (double)(count1 - 1) * 0.5 * c2A / (a2A + a1A * a1A);
				re += count1 * c1A / a1A;
			}
			else
			{
				int count2 = entry2.val;
				double a1B = slide_a1[nhaplo2];
				double a2B = slide_a2[nhaplo2];
				double b2B = (2 * (nhaplo2 * nhaplo2 + nhaplo2 + 3)) / (9 * nhaplo2 * (nhaplo2 - 1));
				double c2B = b2B - (nhaplo2 + 2) / (a1B * nhaplo2) + a2B / (a1B * a1B);

				re += 2 * count1 * (double)count2 * sqrt(c2A * c2B) / (sqrt(a2A * a2B) + a1A * a1B);
			}
		}
	}
	freq1.Clear();
#else
	for (auto entry1 = freq2.begin(); entry1 != freq2.end(); ++entry1)
	{
		int nhaplo1 = entry1->first;
		int count1  = entry1->second;
		double a1A = slide_a1[nhaplo1];
		double b1A = (nhaplo1 + 1.0) / (3.0 * (nhaplo1 - 1.0));
		double c1A = b1A - 1.0 / a1A;
		double a2A = slide_a2[nhaplo1];
		double b2A = (2.0 * (nhaplo1 * nhaplo1 + nhaplo1 + 3.0)) / (9.0 * nhaplo1 * (nhaplo1 - 1.0));
		double c2A = b2A - (nhaplo1 + 2.0) / (a1A * nhaplo1) + a2A / (a1A * a1A);

		for (auto entry2 = entry1; entry2 != freq2.end(); ++entry2)
		{
			int nhaplo2 = entry2->first;
			if (nhaplo1 == nhaplo2)
			{
				re += count1 * (double)(count1 - 1) * 0.5 * c2A / (a2A + a1A * a1A);
				re += count1 * c1A / a1A;
			}
			else
			{
				int count2  = entry2->second;
				double a1B = slide_a1[nhaplo2];
				double a2B = slide_a2[nhaplo2];
				double b2B = (2 * (nhaplo2 * nhaplo2 + nhaplo2 + 3)) / (9 * nhaplo2 * (nhaplo2 - 1));
				double c2B = b2B - (nhaplo2 + 2) / (a1B * nhaplo2) + a2B / (a1B * a1B);

				re += 2 * count1 * (double)count2 * sqrt(c2A * c2B) / (sqrt(a2A * a2B) + a1A * a1B);
			}
		}
	}
	freq2.clear();
#endif

	return re;
}

/* Settle a SWINDOW fo WINDOW */
template<typename REAL>
TARGET bool SWINDOW<REAL>::Settle(int i, int newid)
{
	#define ADD(x) if (isadd) AtomicAddFloat(window<REAL>.x[owid], (REAL)x[i]); x[i] = 0

	if (window_id[i] == newid) return false;

	int owid = window_id[i];  window_id[i] = newid;
	bool isadd = owid != -1 && K[i] >= 1;

	atomic<int>& wC = *(atomic<int>*)&window<REAL>.C[owid];
	atomic<int>& wK = *(atomic<int>*)&window<REAL>.K[owid];

	if (isadd) wK += K[i]; K[i] = 0;

	// 1 Fst
	if (cNei1973)
	{
		ADD(Fst_Nei1973_sum1);
		ADD(Fst_Nei1973_sum2);
	}

	if (cWeir1984)
	{
		ADD(Fst_Weir1984_sum1);
		ADD(Fst_Weir1984_sum2);
	}

	if (cHudson1992)
	{
		ADD(Fst_Hudson1992_sum1);
		ADD(Fst_Hudson1992_sum2);
	}

	if (cHedrick2005)
	{
		ADD(Fst_Hedrick2005_sum1);
		ADD(Fst_Hedrick2005_sum2);
	}

	if (cJost2008)
	{
		ADD(Fst_Jost2008_sum1);
		ADD(Fst_Jost2008_sum2);
	}

	if (cHuang2021_aneu)
	{
		ADD(Fst_Huang2021_aneu_sum1);
		ADD(Fst_Huang2021_aneu_sum2);
	}

	// 2 Absolute divergence
	if (cdxy)
	{
		ADD(dxy_sum);
	}

	// 3 Nucleotide diversity
	if (cpi)
	{
		ADD(pi_sum);
	}

	// 4 Watterson¡¯s thetaw
	if (cthetaw)
	{
		ADD(thetaw_sum);
	}

	// 5 Tajima¡¯s D
	if (ctajimad)
	{
		ADD(d_sum);

		if (isadd)
		{
			Lock(window<REAL>.lock[owid]);
			if (owid != -1) wC += C[i]; C[i] = 0;
			bool issettle = window<REAL>.C[owid] == window<REAL>.N[owid];

#ifdef slide_table
			if (issettle && window<REAL>.nhaplo_freq1[owid].bucket == NULL)
#else	
			if (issettle && *(int*)&window<REAL>.nhaplo_freq2[owid] == NULL)
#endif
				window<REAL>.vd_sum[owid] = SettleTajimaD(nhaplo_freq1[i], nhaplo_freq2[i]);
			else
			{
#ifdef slide_table
				auto& freq1 = nhaplo_freq1[i];
				auto& freq2 = window<REAL>.nhaplo_freq1[owid];
				auto& freqo = window<REAL>.nhaplo_freq2[owid];

				if (freq2.bucket == NULL)
					new (&freq2) TABLE<HASH, int>(true, NULL);

				for (int j = 0; j < freq1.size; ++j)
				{
					TABLE_ENTRY<HASH, int>& entry = freq1.GetEntry(j);
					HASH nhaplo = entry.key;
					int  count  = entry.val;

					if (!freq2.ContainsKey(nhaplo))
						freq2[nhaplo] = count;
					else
						freq2[nhaplo] += count;
				}

				freq1.Clear();

				if (issettle)
				{
					window<REAL>.vd_sum[owid] = SettleTajimaD(freq2, freqo);
					freq2.~TABLE();
				}
#else
				auto& freq1 = nhaplo_freq2[i];
				auto& freq2 = window<REAL>.nhaplo_freq2[owid];
				auto& freqo = window<REAL>.nhaplo_freq1[owid];

				if (*(int*)&freq2 == NULL)
					new (&freq2) map<int, int>();

				for (auto entry : freq1)
				{
					int nhaplo = entry.first;
					int count = entry.second;

					if (freq2.find(nhaplo) == freq2.end())
						freq2[nhaplo] = count;
					else
						freq2[nhaplo] += count;
				}

				freq1.clear();

				if (issettle)
				{
					window<REAL>.vd_sum[owid] = SettleTajimaD(freqo, freq2);
					map<int, int>().swap(freq2);
					*(int*)&freq2 = NULL;
				}
#endif
			}
			UnLock(window<REAL>.lock[owid]);
		}
	}

	// 6 r2
	if (cr2)
	{
		ADD(D2_sum);
		ADD(Q_sum);
	}

	// 7 D'
	if (cdprime)
	{
		ADD(D_abs_sum);
		ADD(Dmax_abs_sum);
	}

	// 8 r2delta
	if (cr2delta)
	{
		ADD(Delta2_sum);
		ADD(R_sum);
	}

	// 9 delta'
	if (cdeltaprime)
	{
		ADD(Delta_abs_sum);
		ADD(Deltamax_abs_sum);
	}

	// 10 Fis
	if (cfis)
	{
		ADD(Fis_Nei1973_sum1);
		ADD(Fis_Nei1973_sum2);
	}

	// 11 Ho
	if (cho)
	{
		ADD(ho_sum);
	}

	// 12 He
	if (che)
	{
		ADD(he_sum);
	}

	// 13 pic
	if (cpic)
	{
		ADD(pic_sum);
	}

	// 14 Ae
	if (cae)
	{
		ADD(ae_sum);
	}

	// 15 I
	if (cI)
	{
		ADD(I_sum);
	}

	return true;
#undef ADD
}

/* Get sliding window index range*/
template<typename REAL>
TARGET void SWINDOW<REAL>::GetWindowId(int64 l, int& st, int& ed)
{
	char* chr = GetLoc(l).GetChrom();
	HASH ha = HashString(chr, (int)strlen(chr));
	CHROM_PROP& prop = chrom_sted[ha];

	if (prop.st == -1)
	{
		st = ed = -1;
		return;
	}

	st = prop.st + Max((int)GetLocPos(l) - (int)window_size + (int)window_step - 1, 0) / window_step;
	ed = Min(prop.st + ((int)GetLocPos(l) - 1) / (int)window_step, prop.ed);
}
#endif

/* Convert phased genotype into unphased genotype */
THREAD2(UpdateUnphaseGenotypes)
{
	// allocate memory
	MEMORY* omemory = locus_memory;
	MEMORY* nmemory = new MEMORY[g_nthread_val];
	MEMORY* tmemory = new MEMORY[g_nthread_val];
	SLOCUS* nslocus = new SLOCUS[nloc];

	// update genotypes
	TABLE<HASH, uint>* Gitab = new TABLE<HASH, uint>[g_nthread_val];
	for (int i = 0; i < g_nthread_val; ++i)
		new (&Gitab[i]) TABLE<HASH, uint>(false, &tmemory[i]);

	ushort* Gtmap = new ushort[maxG * g_nthread_val];
	SetZero(Gtmap, maxG * g_nthread_val);

	// calculate ngeno for each locus
#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
	for (int64 l = 0; l < nloc; ++l)
	{
		threadid = omp_get_thread_num();
		TABLE<HASH, uint>& gitab = Gitab[threadid];

		gitab.Clear();
		GENOTYPE* gtab = GetLoc(l).GetGtab();
		int ngeno = GetLoc(l).ngeno;
		ushort nals[N_MAX_PLOIDY];

		for (int gi = 0; gi < ngeno; ++gi)
		{
			GENOTYPE& gt = gtab[gi];
			int ploidy = gt.Ploidy();
			ushort* oals = gt.GetAlleleArray();
			HASH nha = 0;

			if (gt.Nalleles())
			{
				SetVal(nals, oals, ploidy);
				Sort(nals, ploidy); //unphase
				nha = HashGenotype(nals, ploidy);
			}
			else
				nha = missing_hash[ploidy];

			if (!gitab.ContainsKey(nha))
				gitab[nha] = gitab.size;
		}
		nslocus[l].ngeno = gitab.size;

		PROGRESS_VALUE++;
	}

	BUCKET nbucket;
	nbucket.CreateBucketGT(nslocus);

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
	for (int64 l = 0; l < nloc; ++l)
	{
		threadid = omp_get_thread_num();
		ushort* gtmap = Gtmap + threadid * maxG;
		TABLE<HASH, uint>& gitab = Gitab[threadid];
		new (&nslocus[l]) SLOCUS(nmemory[threadid], slocus[l], gitab, gtmap);

		GENO_READER rg(0, l);
		GENO_WRITER wg(l, &nbucket);

		for (int i = 0; i < nind; ++i)
			wg.Write(gtmap[rg.Read()]);
		wg.FinishWrite();

		PROGRESS_VALUE += 2;
	}

	delete[] Gitab;
	delete[] Gtmap;
	delete[] slocus; slocus = nslocus;
	delete[] omemory; locus_memory = nmemory;
	delete[] tmemory;

	GT = 0; maxG = 0;
	for (int64 l = 0; l < nloc; ++l)
	{
		GT += slocus[l].ngeno;
		maxG = Max(maxG, slocus[l].ngeno);
	}
	geno_bucket.Replace(nbucket);
}

/* Calculate sliding window */
template<typename REAL>
TARGET void CalcSlide()
{
	if (!slide) return;
	if (ad) Exit("\nError: sliding window (-slide) is incompatible with allelic depth (-ad) option.\n");
	if (abs(g_format_val) > BCF) Exit("\nError: sliding window (-slide) function can only be used for VCF/BCF input files.\n");

	EvaluationBegin();

	OpenResFile("-slide", "Sliding window");

	window<REAL>.InitWindow();

	RunThreads(&SlidingWindowThread<REAL>, NULL, NULL, nloc, nloc,"\nCalculating sliding window:\n", g_nthread_val, true);

	window<REAL>.Write();

	window<REAL>.UnInitWindow();

	CloseResFile();

	if (usephase) RunThreads(&UpdateUnphaseGenotypes<REAL>, NULL, NULL, nloc * 3, nloc * 3, "\nUnphasing genotypes:\n", 1, true);

	EvaluationEnd("Sliding window estimation");

	if (slide_plot_val == 1)
		RunRscript("slide_plot.R");
}

/* Calculate sliding window using multiple threads */
THREAD2(SlidingWindowThread)
{
	uint nthread = g_nthread_val;
	double nsec = nloc / (double)nthread + 1e-8;
	uint st = (uint64)(threadid * nsec), ed = (uint64)((threadid + 1) * nsec);

	int swindow_count = (slide_windowsize_val + slide_windowstep_val - 1) / slide_windowstep_val;
	SWINDOW<REAL> swins;
	swins.InitSWindow();

	double* buf = new double[maxK * maxK + 6 * maxK];
	swins.pst = swins.ped = 0;

	for (int64 l = st; l < ed; ++l)
	{
		swins.CalcLocus(l, buf);
		PROGRESS_VALUE++;
	}

	if (swins.pst != -1) swins.Settle1();
	for (int i = 0; i < swindow_count; ++i)
		swins.Settle(i, -1);

	swins.UnInitSWindow();
	delete[] buf;
}

/* Prepare allele frequency and individual ploidy */
THREAD2(SlidePrepare)
{
	int* Ploidy = new int[maxG * g_nthread_val];
	int* Nalleles = new int[maxG * g_nthread_val];
	POP<REAL>* spop = (*(POP<REAL>**)&window<REAL>.cpop_);
	IND<REAL>** inds = spop->inds;

	int n = spop->nind;
	int64* Vt = new int64[n * g_nthread_val];
	byte* Vmin = new byte[n * g_nthread_val];
	byte* Vmax = new byte[n * g_nthread_val];

	SetZero(Vt, n * g_nthread_val);
	SetVal(Vmin, (byte)100, n * g_nthread_val);
	SetZero(Vmax, n * g_nthread_val);

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
	for (int64 l = 0; l < nloc; ++l)
	{
		threadid = omp_get_thread_num();

		int* nalleles = Nalleles + maxG * threadid;
		int* ploidy = Ploidy + maxG * threadid;
		int64* vt = Vt + n * threadid;
		byte* vmin = Vmin + n * threadid;
		byte* vmax = Vmax + n * threadid;

		GENOTYPE* gtab = GetLoc(l).GetGtab();
		int ngeno = GetLoc(l).ngeno;

		for (int gi = 0; gi < ngeno; ++gi)
		{
			GENOTYPE& gt = gtab[gi];
			nalleles[gi] = gt.Nalleles();
			ploidy[gi] = gt.Ploidy();
		}

		//calculate allele frequency
		REAL* fre = spop->GetFreq(l);
		ushort* gcount = spop->GetGenoCount(l);
		GENO_READER rt(spop->ind0id, l);

		for (int j = 0; j < n; ++j)
		{
			int gid = rt.Read();
			if (nalleles[gid])
			{
				byte v = (byte)ploidy[gid];
				vt[j] += v;
				vmin[j] = Min(vmin[j], v);
				vmax[j] = Max(vmax[j], v);
			}
		}

		/////////////////////////////////////////////////////////////////

		for (int gi = 0; gi < ngeno; ++gi)
		{
			GENOTYPE& gt = gtab[gi];
			if (gt.Nalleles() == 0 || gcount[gi] == 0) continue;
			int c = (int)gcount[gi], v = gt.Ploidy();
			ushort* als = gt.GetAlleleArray();
			for (int i = 0; i < v; ++i)
				fre[als[i]] += c;
		}

		Unify(fre, GetLoc(l).k);
		PROGRESS_VALUE++;
	}

	//gather
	for (int tid = 1; tid < g_nthread_val; ++tid)
	{
		int64* vt = Vt + n * tid;
		byte* vmin = Vmin + n * tid;
		byte* vmax = Vmax + n * tid;

		for (int j = 0; j < n; ++j)
		{
			Vt[j] += vt[j];
			Vmin[j] = Min(Vmin[j], vmin[j]);
			Vmax[j] = Max(Vmax[j], vmax[j]);
		}
	}

	minploidy = 100; maxploidy = 0; maxvt = 0;
	for (int i = 0; i < n; ++i)
	{
		inds[i]->vt = Vt[i];
		inds[i]->vmin = Vmin[i] == 100 ? 0 : (byte)Vmin[i];
		inds[i]->vmax = Vmax[i];

		minploidy = Min(minploidy, inds[i]->vmin);
		maxploidy = Max(maxploidy, inds[i]->vmax);
		maxvt = Max(maxvt, inds[i]->vt);
		sumvt += inds[i]->vt;
	}

	delete[] Vt;
	delete[] Vmin;
	delete[] Vmax;
	delete[] Ploidy;
	delete[] Nalleles;

	POP<REAL>** grps = window<REAL>.grps;
	int ngrps = window<REAL>.ngrps;
	int tnhaplo = 0;
	for (int s = 0; s < ngrps; ++s)
	{
		int nhaplo = 0;
		IND<REAL>** inds2 = grps[s]->inds;
		for (int i = 0; i < grps[s]->nind; ++i)
			nhaplo += inds2[i]->vmax;
		grps[s]->nhaplotypes = nhaplo;

		tnhaplo += nhaplo;
	}
	spop->nhaplotypes = tnhaplo;
}

/* Calculate allele frequencies for each population and region */
THREAD2(SlideFreqThread)
{
	if (allele_freq_offset)      delete[] allele_freq_offset;
	if (genotype_count_offset)   delete[] genotype_count_offset;

	//Allocate new offset
	allele_freq_offset = new uint64[nloc];
	genotype_count_offset = new uint64[nloc];
	SetFF(allele_freq_offset, nloc);
	SetFF(genotype_count_offset, nloc);

	//Calculate new offset and total number of alleles and genotypes
	KT = GT = maxK = maxG = 0;
	for (int64 l = 0; l < nloc; ++l)
	{
		allele_freq_offset[l] = KT;
		genotype_count_offset[l] = GT;
		KT += GetLoc(l).k;
		GT += GetLoc(l).ngeno;
		maxK = Max((int)GetLoc(l).k, maxK);
		maxG = Max((int)GetLoc(l).ngeno, maxG);
	}

	/////////////////////////////////////////////////////////////
	
	POP<REAL>** grps = window<REAL>.grps;
	int ngrps = window<REAL>.ngrps;
	POP<REAL>* spop = (*(POP<REAL>**) & window<REAL>.cpop_);

	//Allocate memory for allele frequency and genotype count for each population and region
	spop->AllocFreq();
	spop->CalcFreqGcount();

	for (int s = 0; s < ngrps; ++s)
	{
		grps[s]->AllocFreq();
		grps[s]->CalcFreqGcount();
	}
}