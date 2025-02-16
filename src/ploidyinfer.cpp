/* Ploidy Inference Functions */

#include "vcfpop.h"

template struct SPF<double>;
template struct SPF<float >;

template TARGET void CalcPloidyInference<double>();
template TARGET void CalcPloidyInference<float >();

template TARGET void IND<double>::PloidyInferenceHeader(FILE* fout);
template TARGET void IND<float >::PloidyInferenceHeader(FILE* fout);
template TARGET void IND<double>::WritePloidyInference(FILE* fout);
template TARGET void IND<float >::WritePloidyInference(FILE* fout);
template TARGET double IND<double>::AvgG(int v, int m, double f);
template TARGET double IND<float >::AvgG(int v, int m, double f);
template TARGET void IND<double>::PloidyInferlnL3(umap<int64, SPF<double>>& depth, int v, double f0, double f1, double f2, double& l0, double& l1, double& l2);
template TARGET void IND<float >::PloidyInferlnL3(umap<int64, SPF<float >>& depth, int v, double f0, double f1, double f2, double& l0, double& l1, double& l2);
template TARGET double IND<double>::PloidyInferlnL(umap<int64, SPF<double>>& depth, int v, double f);
template TARGET double IND<float >::PloidyInferlnL(umap<int64, SPF<float >>& depth, int v, double f);
template TARGET void IND<double>::PloidyInference(int v, double& lnL, double& f, umap<int64, SPF<double>>& depth);
template TARGET void IND<float >::PloidyInference(int v, double& lnL, double& f, umap<int64, SPF<float >>& depth);

#define extern 

#undef extern 

/* Calculate ploidy inference using multiple threads */
THREAD2(PloidyInferenceThread)
{
	int nthread = g_nthread_val;
	double nsec = nind / (double)nthread + 1e-8;
	int st = (int)(threadid * nsec), ed = (int)((threadid + 1) * nsec);
	for (int i = st; i < ed; ++i)
	{
		ainds<REAL>[i]->WritePloidyInference(TEMP_FILES[threadid]);
		PROGRESS_VALUE++;
	}
}

/* Calculate individual ploidy inference */
template<typename REAL>
TARGET void CalcPloidyInference()
{
	if (!ploidyinfer) return;
	if (ad) Exit("\nError: ploidy inference (-ploidyinfer) is incompatible with allelic depth (-ad) option.\n");
	if (haplotype) Exit("\nError: ploidy inference (-ploidyinfer) is incompatible with haplotype extraction (-haplotype) option.\n");

	EvaluationBegin();
	OpenResFile("-ploidyinfer", "Ploidy inference");
	OpenTempFiles(g_nthread_val, ".ploidyinfer");

	IND<REAL>::PloidyInferenceHeader(FRES);
	RunThreads(&PloidyInferenceThread<REAL>, NULL, NULL, nind, nind,
		"\nCalculating individual ploiy level:\n", g_nthread_val, true);

	JoinTempFiles(g_nthread_val);
	CloseResFile();

	EvaluationEnd("Ploidy inference");
}

/* Write header row for ploidy inference, TEST */
template<typename REAL>
TARGET void IND<REAL>::PloidyInferenceHeader(FILE* fout)
{

	fprintf(fout, "%s%s%c%c%c%c%c%c%c", g_linebreak_val, g_linebreak_val, g_delimiter_val, g_delimiter_val, g_delimiter_val, g_delimiter_val, g_delimiter_val, g_delimiter_val, g_delimiter_val);
	for (int rl = 0; rl < lreg; ++rl)
		fprintf(fout, "%c", g_delimiter_val);
	for (int i = 1; i <= N_MAX_PLOIDY; ++i)
		if (ploidyinfer_type_val[i])
			fprintf(fout, "%s%c%c%c", PLOIDY_NAME[i], g_delimiter_val, g_delimiter_val, g_delimiter_val);
	if (ploidyinfer_histogram_val == 1)
		fprintf(fout, "Histogram Data");

	fprintf(fout, "%sInd%cPop", g_linebreak_val, g_delimiter_val);
	for (int rl = 0; rl < lreg; ++rl)
		fprintf(fout, "%cRegL%d", g_delimiter_val, rl);
	fprintf(fout, "%c#used%c#unused%c#miss%cEstPloidy%cPosterProb", g_delimiter_val, g_delimiter_val, g_delimiter_val, g_delimiter_val, g_delimiter_val);
	for (int i = 1; i <= N_MAX_PLOIDY; ++i)
		if (ploidyinfer_type_val[i])
			fprintf(fout, "%clnL%cf%cPosterProb", g_delimiter_val, g_delimiter_val, g_delimiter_val);

	if (ploidyinfer_histogram_val == 1)
	{
		double sep = 1.0 / ploidyinfer_nbins_val;
		for (int i = 0; i < ploidyinfer_nbins_val; ++i)
		{
			fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, i * sep);
			fprintf(fout, "~");
			WriteReal(fout, (i + 1) * sep);
		}
	}
}

/* Write result row for ploidy inference, TEST */
template<typename REAL>
TARGET void IND<REAL>::WritePloidyInference(FILE* fout)
{
	POP<REAL>* tpop = apops<REAL>[popid];
	POP<REAL>* treg[N_MAX_REG] = { 0 }; //20200505
	if (lreg >= 1) treg[0] = aregs<REAL>[0][tpop->rid]; //20200505

	for (int rl = 1; rl < lreg; ++rl)
		treg[rl] = aregs<REAL>[rl][treg[rl - 1]->rid];

	int estploidy = 0;
	double sumprob = 0, max_lnL = -1e30;
	double lnL[N_MAX_PLOIDY + 1] = { 0 }, f[N_MAX_PLOIDY + 1] = { 0 }, poster_prob[N_MAX_PLOIDY + 1] = { 0 };
	int64 used = 0, unused = 0, miss = 0;

	//extract ad
	umap<int64, SPF<REAL>> depth;
	SPF<REAL> tspf;
	tspf.count = 0;
	uint adb[16];
	uint* ad2 = adb;
	uint bins[101] = { 0 };
	char filename[PATH_LEN];
	sprintf(filename, "used2-%s.txt", name);
	FILE* f1 = fopen(filename, "wb");
	fprintf(f1, "dA\tdB\r\n");

	for (int64 l = 0; l < nloc; ++l)
	{
		GENOTYPE& gt = GetGenotype(l);//fine

		if (GetLoc(l).k != 2)
		{
			unused++;
			continue;
		}

		if (gt.Nalleles() == 0)
		{
			miss++;
			continue;
		}

		GetAlleleDepth(l, ad2);
		if (ad2[0] + ad2[1] == 0 ||
			(f_dp_b && (ad2[0] + ad2[1] < f_dp_min || ad2[0] + ad2[1] > f_dp_max)))
		{
			miss++;
			continue;
		}

		if (ad2[0] == 0 || ad2[1] == 0 || ad2[1] + ad2[0] < 10) { unused++; continue; }
		int64 ha = (int64)ad2[0] << 32 | ad2[1];
		if (depth.find(ha) == depth.end()) depth[ha] = tspf;
		depth[ha].count++;
		used++;
		fprintf(f1, "%d%c%d\r\n", ad2[0], g_delimiter_val, ad2[1]);
	}

	fclose(f1);

	for (int i = 1; i <= N_MAX_PLOIDY; ++i)
		if (ploidyinfer_type_val[i])
		{
			PloidyInference(i, lnL[i], f[i], depth);
			if (max_lnL < lnL[i])
			{
				max_lnL = lnL[i];
				estploidy = i;
			}
		}

	for (int i = 0; i <= N_MAX_PLOIDY; ++i)
		if (ploidyinfer_type_val[i])
			sumprob += exp(lnL[i] - max_lnL);

	for (int i = 0; i <= N_MAX_PLOIDY; ++i)
		if (ploidyinfer_type_val[i])
			poster_prob[i] = exp(lnL[i] - max_lnL) / sumprob;

	fprintf(fout, "%s%s%c%s", g_linebreak_val, name, g_delimiter_val, tpop->name);
	for (int rl = 0; rl < lreg; ++rl)
		fprintf(fout, "%c%s", g_delimiter_val, treg[rl]->name);

	fprintf(fout, "%c%lld%c%lld%c%lld%c%d%c", g_delimiter_val, used, g_delimiter_val, unused, g_delimiter_val, miss, g_delimiter_val, estploidy, g_delimiter_val);
	WriteReal(fout, poster_prob[estploidy]);

	for (int i = 0; i <= N_MAX_PLOIDY; ++i)
		if (ploidyinfer_type_val[i])
		{
			fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, lnL[i]);
			fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, f[i]);
			fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, poster_prob[i]);
		}

	if (ploidyinfer_histogram_val == 1)
	{
		for (auto d = depth.begin(); d != depth.end(); ++d)
		{
			int dA = d->first >> 32, dB = d->first & 0xFFFFFFFF;
			double fr = dA / (double)(dA + dB);
			int bid = (int)(fr * ploidyinfer_nbins_val);
			bins[bid] += d->second.count;
			bins[ploidyinfer_nbins_val - 1 - bid] += d->second.count;
		}

		for (int i = 0; i < ploidyinfer_nbins_val; ++i)
		{
			fprintf(fout, "%c", g_delimiter_val);
			fprintf(fout, "%d", bins[i]);
		}
	}

}

/* average genotypic frequency at a diallelic locus given m copies of A */
template<typename REAL>
TARGET double IND<REAL>::AvgG(int v, int m, double f)
{
	m = m > v - m ? v - m : m;

	switch (v)
	{
	case 1: return 0.5;
	case 2:
	{
		double div = 6;
		switch (m)
		{
		case 0: return (2 + f) / div;
		case 1: return (2 - 2 * f) / div;
		default: return -1;
		}
	}
	case 3:
	{
		double div = 4;
		switch (m)
		{
		case 0: return (1 + f) / div;
		case 1: return (1 - f) / div;
		default: return -1;
		}
	}
	case 4:
	{
		double div = 30 * (1 + f) * (1 + 2 * f);
		double f2 = f * f, f3 = f2 * f;
		switch (m)
		{
		case 0: return (6 + 27 * f + 38 * f2 + 19 * f3) / div;
		case 1: return (6 + 12 * f - 2 * f2 - 16 * f3) / div;
		case 2: return (6 + 12 * f - 12 * f2 - 6 * f3) / div;
		default: return -1;
		}
	}
	case 5:
	{
		double div = 12 * (1 + f) * (1 + 2 * f);
		double f2 = f * f, f3 = f2 * f;
		switch (m)
		{
		case 0: return (2 + 10 * f + 15 * f2 + 9 * f3) / div;
		case 1: return (2 + 4 * f + f2 - 7 * f3) / div;
		case 2: return (2 + 4 * f - 4 * f2 - 2 * f3) / div;
		default: return -1;
		}
	}
	case 6:
	{
		double div = 84 * (1 + f) * (1 + 2 * f) * (1 + 3 * f) * (1 + 4 * f);
		double f2 = f * f, f3 = f2 * f, f4 = f3 * f, f5 = f4 * f;
		switch (m)
		{
		case 0: return (12 + 150 * f + 708 * f2 + 1581 * f3 + 1726 * f4 + 863 * f5) / div;
		case 1: return (12 + 108 * f + 330 * f2 + 342 * f3 - 150 * f4 - 642 * f5) / div;
		case 2: return (12 + 108 * f + 288 * f2 + 153 * f3 - 402 * f4 - 159 * f5) / div;
		case 3: return (12 + 108 * f + 288 * f2 + 48 * f3 - 332 * f4 - 124 * f5) / div;
		default: return -1;
		}
	}
	case 7:
	{
		double div = 24 * (1 + f) * (1 + 2 * f) * (1 + 3 * f) * (1 + 4 * f);
		double f2 = f * f, f3 = f2 * f, f4 = f3 * f, f5 = f4 * f;
		switch (m)
		{
		case 0: return (3 + 39 * f + 190 * f2 + 438 * f3 + 495 * f4 + 275 * f5) / div;
		case 1: return (3 + 27 * f + 86 * f2 + 96 * f3 - 13 * f4 - 199 * f5) / div;
		case 2: return (3 + 27 * f + 72 * f2 + 54 * f3 - 111 * f4 - 45 * f5) / div;
		case 3: return (3 + 27 * f + 72 * f2 + 12 * f3 - 83 * f4 - 31 * f5) / div;
		default: return -1;
		}
	}
	case 8:
	{
		double div = 90 * (1 + f) * (1 + 2 * f) * (1 + 3 * f) * (1 + 4 * f) * (1 + 5 * f) * (1 + 6 * f);
		double f2 = f * f, f3 = f2 * f, f4 = f3 * f, f5 = f4 * f, f6 = f5 * f, f7 = f6 * f;
		switch (m)
		{
		case 0: return (10 + 245 * f + 2460 * f2 + 13075 * f3 + 39692 * f4 + 69459 * f5 + 67906 * f6 + 33953 * f7) / div;
		case 1: return (10 + 200 * f + 1590 * f2 + 6340 * f3 + 12854 * f4 + 10128 * f5 - 6998 * f6 - 24124 * f7) / div;
		case 2: return (10 + 200 * f + 1530 * f2 + 5590 * f3 + 9356 * f4 + 2622 * f5 - 14192 * f6 - 5116 * f7) / div;
		case 3: return (10 + 200 * f + 1530 * f2 + 5380 * f3 + 7718 * f4 - 1704 * f5 - 9866 * f6 - 3268 * f7) / div;
		case 4: return (10 + 200 * f + 1530 * f2 + 5380 * f3 + 6920 * f4 - 2250 * f5 - 8900 * f6 - 2890 * f7) / div;
		default: return -1;
		}
	}
	case 9:
	{
		double div = 20 * (1 + f) * (1 + 2 * f) * (1 + 3 * f) * (1 + 4 * f) * (1 + 5 * f) * (1 + 6 * f);
		double f2 = f * f, f3 = f2 * f, f4 = f3 * f, f5 = f4 * f, f6 = f5 * f, f7 = f6 * f;
		switch (m)
		{
		case 0: return (2 + 50 * f + 511 * f2 + 2761 * f3 + 8518 * f4 + 15178 * f5 + 15197 * f6 + 8183 * f7) / div;
		case 1: return (2 + 40 * f + 321 * f2 + 1301 * f3 + 2722 * f4 + 2316 * f5 - 961 * f6 - 5741 * f7) / div;
		case 2: return (2 + 40 * f + 306 * f2 + 1136 * f3 + 1966 * f4 + 864 * f5 - 3154 * f6 - 1160 * f7) / div;
		case 3: return (2 + 40 * f + 306 * f2 + 1076 * f3 + 1650 * f4 - 268 * f5 - 2102 * f6 - 704 * f7) / div;
		case 4: return (2 + 40 * f + 306 * f2 + 1076 * f3 + 1384 * f4 - 450 * f5 - 1780 * f6 - 578 * f7) / div;
		default: return -1;
		}
	}
	case 10:
	{
		double div = 132 * (1 + f) * (1 + 2 * f) * (1 + 3 * f) * (1 + 4 * f) * (1 + 5 * f) * (1 + 6 * f) * (1 + 7 * f) * (1 + 8 * f);
		double f2 = f * f, f3 = f2 * f, f4 = f3 * f, f5 = f4 * f, f6 = f5 * f, f7 = f6 * f, f8 = f7 * f, f9 = f8 * f;
		switch (m)
		{
		case 0: return (12 + 486 * f + 8440 * f2 + 82229 * f3 + 493806 * f4 + 1891753 * f5 + 4631876 * f6 + 7090179 * f7 + 6500866 * f8 + 3250433 * f9) / div;
		case 1: return (12 + 420 * f + 6218 * f2 + 50626 * f3 + 246174 * f4 + 721694 * f5 + 1192990 * f6 + 781206 * f7 - 739378 * f8 - 2259962 * f9) / div;
		case 2: return (12 + 420 * f + 6108 * f2 + 47931 * f3 + 219246 * f4 + 580443 * f5 + 776508 * f6 + 99393 * f7 - 1290522 * f8 - 439539 * f9) / div;
		case 3: return (12 + 420 * f + 6108 * f2 + 47436 * f3 + 210468 * f4 + 519492 * f5 + 567156 * f6 - 266940 * f7 - 827136 * f8 - 257016 * f9) / div;
		case 4: return (12 + 420 * f + 6108 * f2 + 47436 * f3 + 207960 * f4 + 489066 * f5 + 431064 * f6 - 312348 * f7 - 669000 * f8 - 200718 * f9) / div;
		case 5: return (12 + 420 * f + 6108 * f2 + 47436 * f3 + 207960 * f4 + 476592 * f5 + 393180 * f6 - 317892 * f7 - 627420 * f8 - 186396 * f9) / div;
		default: return -1;
		}
	}
	default: return -1;
	}
}

/* Calculate the likelihood for ploidy inference */
template<typename REAL>
TARGET void IND<REAL>::PloidyInferlnL3(umap<int64, SPF<REAL>>& depth, int v, double f0, double f1, double f2, double& l0, double& l1, double& l2)
{
	double re0 = 0, re1 = 0, re2 = 0;
	for (auto d = depth.begin(); d != depth.end(); ++d)
	{
		REAL* val = &(d->second.val1[0]);
		REAL* val2 = &(d->second.val2[0]);
		int dA = d->first >> 32, dB = d->first & 0xFFFFFFFF;
		double bicoef = LogBinomial(dA + dB, dA);
		double pr0A = 0, pr1A = 0, pr2A = 0;
		double pr0B = 1, pr1B = 1, pr2B = 1;
		for (int m = 0; m <= v; ++m)
		{
			double G = 0;
			G = AvgG(v, m, f0);
			pr0A += val[m] * G;
			pr0B -= val2[m] * G;
			G = AvgG(v, m, f1);
			pr1A += val[m] * G;
			pr1B -= val2[m] * G;
			G = AvgG(v, m, f2);
			pr2A += val[m] * G;
			pr2B -= val2[m] * G;
		}
		re0 += (log(pr0A / pr0B) + bicoef) * d->second.count;
		re1 += (log(pr1A / pr1B) + bicoef) * d->second.count;
		re2 += (log(pr2A / pr2B) + bicoef) * d->second.count;
	}
	l0 = re0;
	l1 = re1;
	l2 = re2;
}

/* Calculate the likelihood for ploidy inference */
template<typename REAL>
TARGET double IND<REAL>::PloidyInferlnL(umap<int64, SPF<REAL>>& depth, int v, double f)
{
	double re = 0;
	for (auto d = depth.begin(); d != depth.end(); ++d)
	{
		REAL* val = &(d->second.val1[0]);
		REAL* val2 = &(d->second.val2[0]);
		double prA = 0, prB = 1;
		int dA = d->first >> 32, dB = d->first & 0xFFFFFFFF;
		for (int m = 0; m <= v; ++m)
		{
			double G = AvgG(v, m, f);
			prA += val[m] * G;
			prB -= val2[m] * G;
		}
		re += (log(prA / prB) + LogBinomial(dA + dB, dA)) * d->second.count;
	}
	return re;
}

/* Infer the ploidy level from allelic depth distribution */
template<typename REAL>
TARGET void IND<REAL>::PloidyInference(int v, double& lnL, double& f, umap<int64, SPF<REAL>>& depth)
{
	static double LOGFRAC[N_MAX_PLOIDY + 1][N_MAX_PLOIDY + 1]		//Logarithm of fractions
		= { {0.0000000000000000E+00, -1.0000000000000000e+10, -1.0000000000000000e+10, -1.0000000000000000e+10, -1.0000000000000000e+10, -1.0000000000000000e+10, -1.0000000000000000e+10, -1.0000000000000000e+10, -1.0000000000000000e+10, -1.0000000000000000e+10, -1.0000000000000000e+10}, {1.0000000000000000e+10, 0.0000000000000000E+00, -6.9314718055994500E-01, -1.0986122886681100E+00, -1.3862943611198900E+00, -1.6094379124341000E+00, -1.7917594692280500E+00, -1.9459101490553100E+00, -2.0794415416798400E+00, -2.1972245773362200E+00, -2.3025850929940500E+00}, {1.0000000000000000e+10, 6.9314718055994500E-01, 0.0000000000000000E+00, -4.0546510810816400E-01, -6.9314718055994500E-01, -9.1629073187415500E-01, -1.0986122886681100E+00, -1.2527629684953700E+00, -1.3862943611198900E+00, -1.5040773967762700E+00, -1.6094379124341000E+00}, {1.0000000000000000e+10, 1.0986122886681100E+00, 4.0546510810816400E-01, 0.0000000000000000E+00, -2.8768207245178100E-01, -5.1082562376599100E-01, -6.9314718055994500E-01, -8.4729786038720400E-01, -9.8082925301172600E-01, -1.0986122886681100E+00, -1.2039728043259400E+00}, {1.0000000000000000e+10, 1.3862943611198900E+00, 6.9314718055994500E-01, 2.8768207245178100E-01, 0.0000000000000000E+00, -2.2314355131421000E-01, -4.0546510810816400E-01, -5.5961578793542300E-01, -6.9314718055994500E-01, -8.1093021621632900E-01, -9.1629073187415500E-01}, {1.0000000000000000e+10, 1.6094379124341000E+00, 9.1629073187415500E-01, 5.1082562376599100E-01, 2.2314355131421000E-01, 0.0000000000000000E+00, -1.8232155679395500E-01, -3.3647223662121300E-01, -4.7000362924573600E-01, -5.8778666490211900E-01, -6.9314718055994500E-01}, {1.0000000000000000e+10, 1.7917594692280500E+00, 1.0986122886681100E+00, 6.9314718055994500E-01, 4.0546510810816400E-01, 1.8232155679395500E-01, 0.0000000000000000E+00, -1.5415067982725800E-01, -2.8768207245178100E-01, -4.0546510810816400E-01, -5.1082562376599100E-01}, {1.0000000000000000e+10, 1.9459101490553100E+00, 1.2527629684953700E+00, 8.4729786038720400E-01, 5.5961578793542300E-01, 3.3647223662121300E-01, 1.5415067982725800E-01, 0.0000000000000000E+00, -1.3353139262452300E-01, -2.5131442828090600E-01, -3.5667494393873200E-01}, {1.0000000000000000e+10, 2.0794415416798400E+00, 1.3862943611198900E+00, 9.8082925301172600E-01, 6.9314718055994500E-01, 4.7000362924573600E-01, 2.8768207245178100E-01, 1.3353139262452300E-01, 0.0000000000000000E+00, -1.1778303565638400E-01, -2.2314355131421000E-01}, {1.0000000000000000e+10, 2.1972245773362200E+00, 1.5040773967762700E+00, 1.0986122886681100E+00, 8.1093021621632900E-01, 5.8778666490211900E-01, 4.0546510810816400E-01, 2.5131442828090600E-01, 1.1778303565638300E-01, 0.0000000000000000E+00, -1.0536051565782600E-01}, {1.0000000000000000e+10, 2.3025850929940500E+00, 1.6094379124341000E+00, 1.2039728043259400E+00, 9.1629073187415500E-01, 6.9314718055994500E-01, 5.1082562376599100E-01, 3.5667494393873200E-01, 2.2314355131421000E-01, 1.0536051565782600E-01, 0.0000000000000000E+00} };

	for (auto d = depth.begin(); d != depth.end(); ++d)
	{
		int dA = d->first >> 32, dB = d->first & 0xFFFFFFFF, dAB = dA + dB;
		REAL* val = &(d->second.val1[0]);
		REAL* val2 = &(d->second.val2[0]);
		for (int m = 0; m <= v; ++m)
		{
			val[m] = ((m == 0 && dA) || (m == v && dB)) ? 0 : exp(LOGFRAC[m][v] * dA + LOGFRAC[v - m][v] * dB);
			val2[m] = (m == 0 ? 0 : exp(LOGFRAC[m][v] * dAB)) + (m == v ? 0 : exp(LOGFRAC[v - m][v] * dAB));
		}
	}

	//Newton's method
	double dx = 1e-5, x0 = 0.1, maxy0 = -1e30, maxx0 = 0, lastx0 = -1;

	for (int m = 0; m < MAX_ITER_SVD; ++m)
	{
		double y0 = 0, y1 = 0, y2 = 0;
		PloidyInferlnL3(depth, v, x0, x0 + dx, x0 + dx + dx, y0, y1, y2);

		if (y0 > maxy0)
		{
			maxx0 = x0;
			maxy0 = y0;
		}

		double yd0 = (y1 - y0) / dx;
		double yd1 = (y2 - y1) / dx;
		double ydd0 = (yd1 - yd0) / dx;

		lastx0 = x0;
		x0 -= yd0 / ydd0;
		if (x0 <= 0) x0 = (double)1e-10;
		if (x0 >= 1) x0 = (double)0.9999999;

		if (fabs(lastx0 - x0) < 1e-6) break;
	}

	lnL = maxy0;
	f = maxx0;
}
