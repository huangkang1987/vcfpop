/* Linkage Disequilibrium Decay Functions */

#include "vcfpop.h"

template struct DECAY_INTERVAL<double>;
template struct DECAY_INTERVAL<float >;
template struct DECAY<double>;
template struct DECAY<float >;

template TARGET DECAY_INTERVAL<double>::DECAY_INTERVAL();
template TARGET DECAY_INTERVAL<float >::DECAY_INTERVAL();
template TARGET void DECAY_INTERVAL<double>::AddInterval(DECAY_INTERVAL<double>& interval);
template TARGET void DECAY_INTERVAL<float >::AddInterval(DECAY_INTERVAL<float >& interval);
template TARGET void DECAY_INTERVAL<double>::AddDecay(DECAY<double>& ld);
template TARGET void DECAY_INTERVAL<float >::AddDecay(DECAY<float >& ld);
template TARGET void DECAY_INTERVAL<double>::WriteHeader(FILE* fout);
template TARGET void DECAY_INTERVAL<float >::WriteHeader(FILE* fout);
template TARGET void DECAY_INTERVAL<double>::Write(FILE* fout, char* chrom, int id);
template TARGET void DECAY_INTERVAL<float >::Write(FILE* fout, char* chrom, int id);

template TARGET void DECAY<double>::CalcLD(int l1, int* l2, int* buf, int n, TABLE<uint64, int>* genopair);
template TARGET void DECAY<float >::CalcLD(int l1, int* l2, int* buf, int n, TABLE<uint64, int>* genopair);
template TARGET void DECAY<double>::WriteHeader(FILE* fout);
template TARGET void DECAY<float >::WriteHeader(FILE* fout);
template TARGET void DECAY<double>::Write(int id);
template TARGET void DECAY<float >::Write(int id);

template TARGET void CalcDecay<double>();
template TARGET void CalcDecay<float >();

template<> DECAY_INTERVAL<double>*		decay_global_intervals<double>;
template<> DECAY_INTERVAL<float >*		decay_global_intervals<float >;
template<> DECAY_INTERVAL<double>**		decay_chrom_intervals<double>;
template<> DECAY_INTERVAL<float >**		decay_chrom_intervals<float >;

umap<HASH, int>							decay_chrom_id;
vector<char*>							decay_chomas;
double									decay_interval_width;
constexpr int							decay_maxnloc2 = 256;
atomic<int>								decay_current_split = 0;

#ifndef _DECAY_INTERVAL
/* Set zero */
template<typename REAL>
TARGET DECAY_INTERVAL<REAL>::DECAY_INTERVAL()
{
	SetZero(this, 1);
}

/* Add thread specific interval to global specific interval */
template<typename REAL>
TARGET void DECAY_INTERVAL<REAL>::AddInterval(DECAY_INTERVAL<REAL>& interval)
{
	npair += interval.npair;
	Add(&r2A, &interval.r2A, 16);
}

/* Add LD measures to thread specific interval */
template<typename REAL>
TARGET void DECAY_INTERVAL<REAL>::AddDecay(DECAY<REAL>& ld)
{
	npair++;

	r2A			+= ld.r2A;
	r2A2		+= ld.r2A * ld.r2A / ld.r2B;
	r2B			+= ld.r2B;
	r2B2		+= ld.r2B * ld.r2B;

	r2DeltaA	+= ld.r2DeltaA;
	r2DeltaA2	+= ld.r2DeltaA * ld.r2DeltaA / ld.r2DeltaB;
	r2DeltaB	+= ld.r2DeltaB;
	r2DeltaB2	+= ld.r2DeltaB * ld.r2DeltaB;

	DpA			+= ld.DpA;
	DpA2		+= ld.DpA * ld.DpA / ld.DpB;
	DpB			+= ld.DpB;
	DpB2		+= ld.DpB * ld.DpB;

	DeltapA		+= ld.DeltapA;
	DeltapA2	+= ld.DeltapA * ld.DeltapA / ld.DeltapB;
	DeltapB		+= ld.DeltapB;
	DeltapB2	+= ld.DeltapB * ld.DeltapB;
}

/* Write table header */
template<typename REAL>
TARGET void DECAY_INTERVAL<REAL>::WriteHeader(FILE* fout)
{
	fprintf(fout, "%s%sChrom%cMinDist%cMaxDist%cnpair", g_linebreak_val, g_linebreak_val, g_delimiter_val, g_delimiter_val, g_delimiter_val);
	if (decay_estimator_val[1]) fprintf(fout, "%cr2%cSE",		g_delimiter_val, g_delimiter_val);
	if (decay_estimator_val[2]) fprintf(fout, "%cr2Delta%cSE",	g_delimiter_val, g_delimiter_val);
	if (decay_estimator_val[3]) fprintf(fout, "%cD'%cSE",		g_delimiter_val, g_delimiter_val);
	if (decay_estimator_val[4]) fprintf(fout, "%cDelta'%cSE",	g_delimiter_val, g_delimiter_val);

	fprintf(fout, "%s", g_linebreak_val);
}

/* Write table entry */
template<typename REAL>
TARGET void DECAY_INTERVAL<REAL>::Write(FILE* fout, char* chrom, int id)
{
	int mindist = (int)(id * decay_interval_width) + 1;
	int maxdist = (int)((id + 1) * decay_interval_width);

	if (npair)
	{
		if (block_estimator_val[1])
		{
			r2A /= r2B;
			r2A2 = sqrt((r2A2 / r2B - r2A * r2A) / (1 - r2B2 / (r2B * r2B)) * r2B2 / (r2B * r2B));
			//r2B = r2A / r2A2;
			//r2B2 = MinusLogPNormal(r2B);
		}

		if (block_estimator_val[2])
		{
			r2DeltaA /= r2DeltaB;
			r2DeltaA2 = sqrt((r2DeltaA2 / r2DeltaB - r2DeltaA * r2DeltaA) / (1 - r2DeltaB2 / (r2DeltaB * r2DeltaB)) * r2DeltaB2 / (r2DeltaB * r2DeltaB));
			//r2DeltaB = r2DeltaA / r2DeltaA2;
			//r2DeltaB2 = MinusLogPNormal(r2DeltaB);
		}

		if (block_estimator_val[3])
		{
			DpA /= DpB;
			DpA2 = sqrt((DpA2 / DpB - DpA * DpA) / (1 - DpB2 / (DpB * DpB)) * DpB2 / (DpB * DpB));
			//DpB = DpA / DpA2;
			//DpB2 = MinusLogPNormal(DpB);
		}

		if (block_estimator_val[4])
		{
			DeltapA /= DeltapB;
			DeltapA2 = sqrt((DeltapA2 / DeltapB - DeltapA * DeltapA) / (1 - DeltapB2 / (DeltapB * DeltapB)) * DeltapB2 / (DeltapB * DeltapB));
			//DeltapB = DeltapA / DeltapA2;
			//DeltapB2 = MinusLogPNormal(DeltapB);
		}
	}

	fprintf(fout, "%s%c", chrom, g_delimiter_val);
	fprintf(fout, "%d%c", mindist, g_delimiter_val);
	fprintf(fout, "%d%c", maxdist, g_delimiter_val);
	fprintf(fout, "%d", npair);

	//Weighted average
	if (decay_estimator_val[1])
	{
		if (npair == 0)
			fprintf(fout, "%cnan%cnan", g_delimiter_val, g_delimiter_val);
		else
		{
			fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, r2A);
			fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, r2A2);
			//fprintf(fout, "%c", g_delimiter_val);
			//WriteReal(fout, r2B2);
		}
	}

	if (decay_estimator_val[2])
	{
		if (npair == 0)
			fprintf(fout, "%cnan%cnan", g_delimiter_val, g_delimiter_val);
		else
		{
			fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, r2DeltaA);
			fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, r2DeltaA2);
			//fprintf(fout, "%c", g_delimiter_val);
			//WriteReal(fout, r2DeltaB2);
		}
	}

	if (decay_estimator_val[3])
	{
		if (npair == 0)
			fprintf(fout, "%cnan%cnan", g_delimiter_val, g_delimiter_val);
		else
		{
			fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, DpA);
			fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, DpA2);
			//fprintf(fout, "%c", g_delimiter_val);
			//WriteReal(fout, DpB2);
		}
	}

	if (decay_estimator_val[4])
	{
		if (npair == 0)
			fprintf(fout, "%cnan%cnan", g_delimiter_val, g_delimiter_val);
		else
		{
			fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, DeltapA);
			fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, DeltapA2);
			//fprintf(fout, "%c", g_delimiter_val);
			//WriteReal(fout, DeltapB2);
		}
	}
	fprintf(fout, "%s", g_linebreak_val);
}
#endif

#ifndef _DECAY
/* Calculate LD measures between loc1 and several loc2 */
template<typename REAL>
TARGET void DECAY<REAL>::CalcLD(int l1, int* l2, int* buf, int n, TABLE<uint64, int>* genopair)
{
	int k1 = GetLoc(l1).k;
	GENOTYPE* gtab1 = GetLoc(l1).GetGtab(), *gtab2[decay_maxnloc2];

	int* mA  = buf + maxK * 0;	//maxK
	int* mB  = buf + maxK * 1;	//maxK
	int* mAB = buf + maxK * 2;	//maxK * maxK

	int* nAA = buf + maxK * 0;	//maxK
	int* nBB = buf + maxK * 1;	//maxK
	int* nA  = buf + maxK * 2;	//maxK
	int* nB  = buf + maxK * 3;	//maxK
	int* nAB = buf + maxK * 4;	//maxK * maxK

	for (int kk = 0; kk < n; ++kk)
	{
		this[kk].loc1 = l1;
		this[kk].loc2 = l2[kk];
		this[kk].nhaplo = 0;
		gtab2[kk] = GetLoc(l2[kk]).GetGtab();
	}

	////////////////////////////////////////////////////////////
	//TABLE<uint64, int> genopair[decay_maxnloc2];
	//REP(n) new (&genopair[kk]) TABLE<uint64, int>(false, NULL, 64);

	REP(n) genopair[kk].Clear();
	GENO_READER rt1(cpop<REAL>->ind0id, l1), rt2[decay_maxnloc2];
	REP(n) new (&rt2[kk]) GENO_READER(cpop<REAL>->ind0id, l2[kk]);

	for (int i = 0; i < cpop<REAL>->nind; ++i)
	{
		uint64 id1 = rt1.Read();
		id1 = (id1 << 32) - id1;
		REP(n) genopair[kk].AddCount(id1 + rt2[kk].Read());
	}

	////////////////////////////////////////////////////////////
	if (decay_estimator_val[1] || decay_estimator_val[3])
	{
		for (int kk = 0; kk < n; ++kk)
		{
			SetZero(buf, maxK * maxK + maxK * 2);

			int cnhaplo = 0;
			int k2 = GetLoc(l2[kk]).k;

			for (int i = 0; i <= genopair[kk].mask; ++i)
			{
				TABLE_ENTRY<uint64, int>& kv = genopair[kk].bucket[i];
				if (kv.val == (uint64)-1) continue;

				uint64 key = kv.key;
				int count = kv.val;
				uint id1 = key / 0xFFFFFFFFull, id2 = key - id1 * 0xFFFFFFFFull;
				GENOTYPE& g1 = gtab1[id1], & g2 = gtab2[kk][id2];

				if (g1.Nalleles() == 0 || g2.Nalleles() == 0) continue;

				int v1 = g1.Ploidy(), v2 = g2.Ploidy();
				if (v1 != v2) continue;

				cnhaplo += v1 * count;
				ushort* als1 = g1.GetAlleleArray(), * als2 = g2.GetAlleleArray();
				for (int a = 0; a < v1; ++a)
				{
					mA[als1[a]] += count;
					mB[als2[a]] += count;
					mAB[als1[a] * maxK + als2[a]] += count;
				}
			}

			double nhaploinv = 1.0 / cnhaplo;
			double cr2A = 0, cr2B = 0, cDpA = 0, cDpB = 0;
			for (int A = 0; A < k1; ++A)
			{
				if (mA[A] == 0) continue;
				for (int B = 0; B < k2; ++B)
				{
					if (mB[B] == 0) continue;
					double pA = mA[A] * nhaploinv, pB = mB[B] * nhaploinv;
					double DAB = mAB[A * maxK + B] * nhaploinv - pA * pB;

					cr2A += DAB * DAB;
					cr2B += pA * pB * (1 - pA) * (1 - pB);

					cDpA += fabs(DAB);
					cDpB += (DAB >= 0 ?
						(pA < pB ? pA * (1 - pB) : pB * (1 - pA)) :
						(pA + pB <= 1 ? pA * pB : (1 - pA) * (1 - pB)));
				}
			}
			this[kk].nhaplo = cnhaplo;
			this[kk].r2A    = cr2A;
			this[kk].r2B    = cr2B;
			this[kk].DpA    = cDpA;
			this[kk].DpB    = cDpB;
		}
	}

	////////////////////////////////////////////////////////
	if (decay_estimator_val[2] || decay_estimator_val[4])
	{
		for (int kk = 0; kk < n; ++kk)
		{
			SetZero(buf, maxK * maxK + maxK * 4);

			int cnhaplo = 0;
			int sv2 = 0, ntAABB = 0;
			int k2 = GetLoc(l2[kk]).k;

			for (int i = 0; i <= genopair[kk].mask; ++i)
			{
				TABLE_ENTRY<uint64, int>& kv = genopair[kk].bucket[i];
				if (kv.val == (uint64)-1) continue;

				uint64 key = kv.key;
				int count = kv.val;
				uint id1 = key / 0xFFFFFFFFull, id2 = key - id1 * 0xFFFFFFFFull;
				GENOTYPE& g1 = gtab1[id1], & g2 = gtab2[kk][id2];
				if (g1.Nalleles() == 0 || g2.Nalleles() == 0) continue;

				int nalleles1 = g1.Nalleles(), nalleles2 = g2.Nalleles();
				if (nalleles1 == 0 || nalleles2 == 0) continue;

				int v1 = g1.Ploidy(), v2 = g2.Ploidy();
				if (v1 != v2) continue;

				sv2 += v1 * v1 * count;
				cnhaplo += v1 * count;
				ushort* als1 = g1.GetAlleleArray() + v1, * als2 = g2.GetAlleleArray() + v2;
				uint64 pattern1 = g1.GetPattern(), pattern2 = g2.GetPattern();

				uint64 p1 = pattern1;
				for (int a1 = 0; a1 < nalleles1; ++a1)
				{
					int ncopy1 = p1 & 0xF; p1 >>= 4;

					nA[als1[a1]] += ncopy1 * count;
					nAA[als1[a1]] += ncopy1 * (ncopy1 - 1) * count;

					uint64 p2 = pattern2;
					for (int a2 = 0; a2 < nalleles2; ++a2)
					{
						int ncopy2 = p2 & 0xF; p2 >>= 4;
						nAB[als1[a1] * maxK + als2[a2]] += ncopy1 * ncopy2 * count;
					}
				}

				uint64 p2 = pattern2;
				for (int a2 = 0; a2 < nalleles2; ++a2)
				{
					int ncopy2 = p2 & 0xF; p2 >>= 4;
					nB[als2[a2]] += ncopy2 * count;
					nBB[als2[a2]] += ncopy2 * (ncopy2 - 1) * count;
				}
			}

			double invndpair = 1.0 / (sv2 - cnhaplo);
			double invnhaplo = 1.0 / cnhaplo;
			double vtilde = sv2 / (double)cnhaplo;
			double cr2DeltaA = 0, cr2DeltaB = 0, cDeltapA = 0, cDeltapB = 0;

			for (int A = 0; A < k1; ++A)
			{
				if (nA[A] == 0) continue;
				double pA = nA[A] * invnhaplo, pAA = nAA[A] * invndpair, pA2 = pA * pA, pAX = pA - pAA;

				for (int B = 0; B < k2; ++B)
				{
					if (nB[B] == 0) continue;
					double pB = nB[B] * invnhaplo, pBB = nBB[B] * invndpair, pB2 = pB * pB;
					double DeltaAB = nAB[A * maxK + B] * invnhaplo - vtilde * pA * pB;
					double lambda = pA + pB - 1;

					cr2DeltaA = DeltaAB * DeltaAB;
					cr2DeltaB = (pA + (vtilde - 1) * pAA - vtilde * pA2) * (pB + (vtilde - 1) * pBB - vtilde * pB2);

					cDeltapA = fabs(DeltaAB);
					cDeltapB = fabs(DeltaAB >= 0 ?
						(pA < pB ?
							pA + (vtilde - 1) * (pAA + pAX * (pB - pA) / (1 - pA)) - vtilde * pA * pB :
							pB + (vtilde - 1) * (pAA * pB / pA) - vtilde * pA * pB) :
						(lambda <= 0 ?
							(vtilde - 1) * (pAX * pB / (1 - pA)) - vtilde * pA * pB :
							lambda + (vtilde - 1) * (pAA * lambda / pA + pAX) - vtilde * pA * pB));
				}
			}
			this[kk].nhaplo   = cnhaplo;
			this[kk].r2DeltaA = cr2DeltaA;
			this[kk].r2DeltaB = cr2DeltaB;
			this[kk].DeltapA  = cDeltapA;
			this[kk].DeltapB  = cDeltapB;
		}
	}
}

/* Write LD measures header */
template<typename REAL>
TARGET void DECAY<REAL>::WriteHeader(FILE* fout)
{
	fprintf(fout, "%s%sChrom%cLoc1%cLoc2%cnHap%cDist", g_linebreak_val, g_linebreak_val, g_delimiter_val, g_delimiter_val, g_delimiter_val, g_delimiter_val);

	if (decay_estimator_val[1]) fprintf(fout, "%cr2A%cr2B", g_delimiter_val, g_delimiter_val);
	if (decay_estimator_val[2]) fprintf(fout, "%cr2DeltaA%cr2DeltaB", g_delimiter_val, g_delimiter_val);
	if (decay_estimator_val[3]) fprintf(fout, "%cD'A%cD'B", g_delimiter_val, g_delimiter_val);
	if (decay_estimator_val[4]) fprintf(fout, "%cDelta'A%cDelta'B", g_delimiter_val, g_delimiter_val);

	fprintf(fout, "%s", g_linebreak_val);
}

/* Write LD measures to a temp file */
template<typename REAL>
TARGET void DECAY<REAL>::Write(int id)
{
	char name_buf[NAME_BUF_LEN];
	FILE* fout = TEMP_FILES[id];

	fprintf(fout, "%s%c", GetLoc(loc1).GetChrom(), g_delimiter_val);
	fprintf(fout, "%s%c", GetLoc(loc1).GetNameStr(name_buf), g_delimiter_val);
	fprintf(fout, "%s%c", GetLoc(loc2).GetNameStr(name_buf), g_delimiter_val);
	fprintf(fout, "%d%c", nhaplo, g_delimiter_val);
	fprintf(fout, "%lld", GetLocPos(loc2) - GetLocPos(loc1));

	if (decay_estimator_val[1]) 
	{ 
		fprintf(fout, "%c", g_delimiter_val); 
		WriteScientific(fout, r2A); 
		fprintf(fout, "%c", g_delimiter_val); 
		WriteScientific(fout, r2B);
	}

	if (decay_estimator_val[2]) 
	{ 
		fprintf(fout, "%c", g_delimiter_val); 
		WriteScientific(fout, r2DeltaA);
		fprintf(fout, "%c", g_delimiter_val); 
		WriteScientific(fout, r2DeltaB);
	}
	
	if (decay_estimator_val[3]) 
	{ 
		fprintf(fout, "%c", g_delimiter_val); 
		WriteScientific(fout, DpA);
		fprintf(fout, "%c", g_delimiter_val); 
		WriteScientific(fout, DpB);
	}
	
	if (decay_estimator_val[4]) 
	{ 
		fprintf(fout, "%c", g_delimiter_val); 
		WriteScientific(fout, DeltapA);
		fprintf(fout, "%c", g_delimiter_val);
		WriteScientific(fout, DeltapB);
	}

	fprintf(fout, "%s", g_linebreak_val);
}
#endif

/* Calculate LD decay */
template<typename REAL>
TARGET void CalcDecay()
{
	if (!decay) return;
	if (ad) Exit("\nError: LD decay (-decay) is incompatible with allelic depth (-ad) option.\n");
	if (abs(g_format_val) > BCF) Exit("\nError: LD decay (-decay) function can only be used for VCF/BCF input files.\n");

	EvaluationBegin();

	decay_interval_width = decay_maxdist_val / (double)decay_nintervals_val;

	OpenResFile("-decay", "LD decay");
	if (decay_pair_val == 1)
		OpenTempFiles(g_nthread_val * 4, ".decay");
	decay_current_split = 0;

	// assign cpop
	AssignPop<REAL>(decay_pop_b, decay_pop_val, "-decay_pop");

	RunThreads(&DecayFreqThread<REAL>, NULL, NULL, nloc * cpop<REAL>->nind, nloc * cpop<REAL>->nind,
		"\nPreparing allele frequency:\n", 1, true);

	// add chromosomes
	for (int64 l = 0; l < nloc; ++l)
	{
		char* chr = GetLoc(l).GetChrom();
		HASH ha = HashString(chr, (int)strlen(chr));
		if (decay_chrom_id.find(ha) == decay_chrom_id.end())
		{
			decay_chrom_id[ha] = decay_chomas.size();
			decay_chomas.push_back(chr);
		}
	}

	//allocate intervals
	bool multi_chrom = decay_chrom_id.size() > 1 && decay_chromosome_val == 1;
	decay_global_intervals<REAL> = new DECAY_INTERVAL<REAL>[decay_nintervals_val];

	if (multi_chrom)
	{
		decay_chrom_intervals<REAL> = new DECAY_INTERVAL<REAL> *[decay_chrom_id.size()];
		for (int i = 0; i < decay_chrom_id.size(); ++i)
			decay_chrom_intervals<REAL>[i] = new DECAY_INTERVAL<REAL>[decay_nintervals_val];
	}

	RunThreads(&DecayThread<REAL>, NULL, NULL, nloc, nloc, "\nCalculating LD decay:\n", g_nthread_val, true);

	DECAY_INTERVAL<REAL>::WriteHeader(FRES);

	for (int j = 0; j < decay_nintervals_val; ++j)
		decay_global_intervals<REAL>[j].Write(FRES, (char*)"All", j);

	if (multi_chrom)
	{
		for (int i = 0; i < decay_chrom_id.size(); ++i)
			for (int j = 0; j < decay_nintervals_val; ++j)
				decay_chrom_intervals<REAL>[i][j].Write(FRES, decay_chomas[i], j);
	}

	if (decay_pair_val == 1)
	{
		DECAY<REAL>::WriteHeader(FRES);
		JoinTempFiles(g_nthread_val * 4);
	}
	CloseResFile();

	if (multi_chrom)
	{
		for (int i = 0; i < decay_chrom_id.size(); ++i)
			DEL(decay_chrom_intervals<REAL>[i]);
		DEL(decay_chrom_intervals<REAL>);
	}

	DEL(decay_global_intervals<REAL>);

	umap<HASH, int>().swap(decay_chrom_id);
	vector<char*>().swap(decay_chomas);

	cpop<REAL>->UnAllocFreq();
	DEL(allele_freq_offset);
	DEL(genotype_count_offset);

	EvaluationEnd("LD decay analysis");

	if (decay_plot_val == 1)
		RunRscript("decay_plot.R");
}

/* Calculate LD decay using multiple threads */
THREAD2(DecayThread)
{
	//allocate intervals
	bool multi_chrom = decay_chrom_id.size() > 1 && decay_chromosome_val == 1;

	DECAY_INTERVAL<REAL>* global_intervals = NULL, ** chrom_intervals = NULL;
	global_intervals = new DECAY_INTERVAL<REAL>[decay_nintervals_val];

	if (multi_chrom)
	{
		chrom_intervals = new DECAY_INTERVAL<REAL>*[decay_chomas.size()];
		for (int i = 0; i < decay_chomas.size(); ++i)
			chrom_intervals[i] = new DECAY_INTERVAL<REAL>[decay_nintervals_val];
	}

	int loc2[decay_maxnloc2], * buf = new int[maxK * maxK + maxK * 4];
	int max_split = g_nthread_val * 4;
	TABLE<uint64, int> genopair[decay_maxnloc2];
	REP(decay_maxnloc2) new (&genopair[kk]) TABLE<uint64, int>(false, NULL, 64);

	for (int csplit = decay_current_split.fetch_add(1); csplit < max_split; csplit = decay_current_split.fetch_add(1))
	{
		double nsec = nloc / (double)max_split +1e-8;
		uint st = (uint64)(csplit * nsec), ed = (uint64)((csplit + 1) * nsec);

		for (uint l1 = st; l1 < ed; ++l1)
		{
			DECAY<REAL> decay[decay_maxnloc2];
			RNG<double> rng(g_seed_val + l1, RNG_SALT_LDDECAY);
			char* chr = GetLoc(l1).GetChrom();
			HASH ha = HashString(chr, (int)strlen(chr));
			int chrid = decay_chrom_id[ha];

			for (uint l2 = l1 + 1, nloc2 = 0; l2 < nloc; ++l2)
			{
				bool breakflag = false;
				if (strcmp(GetLoc(l1).GetChrom(), GetLoc(l2).GetChrom()) || GetLocPos(l2) - GetLocPos(l1) > decay_maxdist_val)
					breakflag = true;
				else if (decay_ratio_val == 1 || rng.Uniform() < decay_ratio_val)
					loc2[nloc2++] = l2;

				if ((breakflag || l2 == nloc - 1) && nloc2 || nloc2 == decay_maxnloc2 )
				{
					decay[0].CalcLD(l1, loc2, buf, nloc2, genopair);
					for (int j = 0; j < nloc2; ++j)
					{
						int dist2 = GetLocPos(loc2[j]) - GetLocPos(l1);
						int intervalid = dist2 == 0 ? 0 : (int)((dist2 - 1) / decay_interval_width);
						global_intervals[intervalid].AddDecay(decay[j]);
						if (multi_chrom) chrom_intervals[chrid][intervalid].AddDecay(decay[j]);
						if (decay_pair_val == 1) decay[j].Write(threadid);
					}
					nloc2 = 0;
					if (breakflag) break;
				}
			}

			PROGRESS_VALUE++;
		}
	}

	Lock(GLOCK1);
	if (multi_chrom)
	{
		for (int i = 0; i < decay_chomas.size(); ++i)
			for (int j = 0; j < decay_nintervals_val; ++j)
				decay_chrom_intervals<REAL>[i][j].AddInterval(chrom_intervals[i][j]);
	}

	for (int j = 0; j < decay_nintervals_val; ++j)
		decay_global_intervals<REAL>[j].AddInterval(global_intervals[j]);
	UnLock(GLOCK1);

	DEL(buf);

	if (multi_chrom)
	{
		for (int i = 0; i < decay_chomas.size(); ++i)
			DEL(chrom_intervals[i]);
		DEL(chrom_intervals);
	}

	DEL(global_intervals);
}

/* Calculate allele frequencies for cpop<REAL> */
THREAD2(DecayFreqThread)
{
	DEL(allele_freq_offset);
	DEL(genotype_count_offset);

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
		maxK = std::max((int)GetLoc(l).k, maxK);
		maxG = std::max((int)GetLoc(l).ngeno, maxG);
	}

	//Allocate memory for allele frequency and genotype count for each population and region
	cpop<REAL>->AllocFreq();
	cpop<REAL>->CalcFreqGcount();
}
