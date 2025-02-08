/* Genetic Diversity Functions */

#include "vcfpop.h"

template struct DIVERSITY<double>;
template struct DIVERSITY<float >;
template struct DIVSUM<double>;
template struct DIVSUM<float >;

template TARGET DIVERSITY<double>::DIVERSITY();
template TARGET DIVERSITY<float >::DIVERSITY();

template TARGET void CalcDiversity<double>();
template TARGET void CalcDiversity<float >();

template<> DIVERSITY<double>* diversity_buf<double>;
template<> DIVERSITY<float >* diversity_buf<float >;

template<> DIVSUM<double> diversity_sum<double>;
template<> DIVSUM<float > diversity_sum<float >;

#ifndef _DIVSUM
/* Initialize sum */
template<typename REAL>
TARGET void DIVSUM<REAL>::Init()
{
    SetZero(this, 1);
    NE1P = 1;
    NE2P = 1;
    NEPP = 1;
    NEID = 1;
    NESI = 1;
	minploidy = 100;
	maxploidy = 0;
    InitLock(lock);
}
/* Uninitialize sum */
template<typename REAL>
TARGET void DIVSUM<REAL>::UnInit()
{
    SetZero(this, 1);
    NE1P = 1;
    NE2P = 1;
    NEPP = 1;
    NEID = 1;
    NESI = 1;
	minploidy = 100;
	maxploidy = 0;
    UnInitLock(lock);
}

/* Add loc to the sum */
template<typename REAL>
TARGET void DIVSUM<REAL>::Add(DIVERSITY<REAL>& loc)
{
	Lock(lock);

	if (IsNormal(loc.NE1P)) NE1P *= loc.NE1P;
	if (IsNormal(loc.NE2P)) NE2P *= loc.NE2P;
	if (IsNormal(loc.NEPP)) NEPP *= loc.NEPP;
	if (IsNormal(loc.NEID)) NEID *= loc.NEID;
	if (IsNormal(loc.NESI)) NESI *= loc.NESI;

	if (loc.k >= 2)
	{
		ChargeSum(loc.k, k, kc);
		ChargeSum(loc.n, n, nc);
		ChargeSum(loc.nhaplo, nhaplo, nhaploc);
		ChargeSum(loc.ptype, ptype, ptypec);
		ChargeSum(loc.pval, pval, pvalc);

		ChargeWeight(loc.ho,  loc.how, ho,  how);
		ChargeWeight(loc.he,  loc.hew, he,  hew);
		ChargeWeight(loc.pi,  loc.piw, pi,  piw);
		ChargeWeight(loc.pic, loc.hew, pic, picw);
		ChargeWeight(loc.ae,  loc.hew, ae,  aew);
		ChargeWeight(loc.I,   loc.hew, I,   Iw);
	}

	if (loc.minploidy != 0 && loc.minploidy != 100)
		minploidy = std::min((int)loc.minploidy, minploidy);
	if (loc.maxploidy != 0 && loc.maxploidy != 100)
		maxploidy = std::max((int)loc.maxploidy, maxploidy);

	UnLock(lock);
}

/* Output diversity indices */
template<typename REAL>
TARGET void DIVSUM<REAL>::Write(FILE* fout, const char* name)
{
	k /= kc;
	n /= nc;
	nhaplo /= nhaploc;
	ptype /= ptypec;
	pval /= pvalc;

	ho /= how;
	he /= hew;
	pi /= piw;
	pic /= picw;
	ae /= aew;
	I /= Iw;

	fprintf(fout, "%s%s%c(%llu)%c%d-%d%c",
		g_linebreak_val, name,
		g_delimiter_val, (uint64)nloc,
		g_delimiter_val, minploidy,
		maxploidy, g_delimiter_val);

	WriteReal(fout, k);		fprintf(fout, "%c", g_delimiter_val);
	WriteReal(fout, n);		fprintf(fout, "%c", g_delimiter_val);
	WriteReal(fout, nhaplo);fprintf(fout, "%c", g_delimiter_val);

	WriteReal(fout, ho);	fprintf(fout, "%c", g_delimiter_val);
	WriteReal(fout, he);	fprintf(fout, "%c", g_delimiter_val);
	WriteReal(fout, pi);	fprintf(fout, "%c", g_delimiter_val);
	WriteReal(fout, pic);	fprintf(fout, "%c", g_delimiter_val);
	WriteReal(fout, ae);	fprintf(fout, "%c", g_delimiter_val);
	WriteReal(fout, I);		fprintf(fout, "%c", g_delimiter_val);

	WriteReal(fout, NE1P);	fprintf(fout, "%c", g_delimiter_val);
	WriteReal(fout, NE2P);	fprintf(fout, "%c", g_delimiter_val);
	WriteReal(fout, NEPP);	fprintf(fout, "%c", g_delimiter_val);
	WriteReal(fout, NEID);	fprintf(fout, "%c", g_delimiter_val);
	WriteReal(fout, NESI);	fprintf(fout, "%c", g_delimiter_val);

	WriteReal(fout, 1 - ho / he); fprintf(fout, "%c", g_delimiter_val);
	fprintf(fout, "-%c-%c-", g_delimiter_val, g_delimiter_val);
}

#endif

#ifndef _DIVERSITY
/* Initialize */
template<typename REAL>
TARGET DIVERSITY<REAL>::DIVERSITY()
{
	SetZero(this, 1);
}

/* Calculate diveristy indices */
template<typename REAL>
TARGET void DIVERSITY<REAL>::CalcDiversity(int64 _l)
{
	l = _l;
	minploidy = 100;
	maxploidy = 0;

	int k2 = GetLoc(l).k, ploidy = 0, ni = cpop<REAL>->nind;

	varploidy = false;
	n = 0;

	if (ni == 0)
	{
		nhaplo = 0;
		minploidy = 0;
		maxploidy = 0;

		bmaf = NAN;
		k = 0;
		ptype = NAN;
		pval = NAN;

		how = NAN;
		hew = NAN;
		piw = NAN;

		ho = NAN;
		he = NAN;
		pi = NAN;
		pic = NAN;
		ae = NAN;
		I = NAN;

		NE1P = NAN;
		NE2P = NAN;
		NEPP = NAN;
		NEID = NAN;
		NESI = NAN;

		return;
	}

	REAL* fre = cpop<REAL>->GetFreq(l);
	ho = 0; how = 0; nhaplo = 0;

	ushort* gcount = cpop<REAL>->GetGenoCount(l);
	LOCSTAT1& stat1 = cpop<REAL>->loc_stat1[l];

	//during diversity filter
	if (!reassigned)
	{
		GENO_READER rt(cpop<REAL>->ind0id, l);
		for (int i = 0; i < nind; ++i)
			gcount[rt.Read()]++;
	}

	GENOTYPE* gtab = GetLoc(l).GetGtab();
	int ngeno = GetLoc(l).ngeno;

	v2i = 0;
	if (ad == 0) for (int gi = 0; gi < ngeno; ++gi)
	{
		GENOTYPE& gt = gtab[gi];
		if (gt.Nalleles() == 0 || gcount[gi] == 0) continue;
		uint c = gcount[gi];
		int v = gt.Ploidy();
		nhaplo += v * c;

		if (v > maxploidy) maxploidy = (byte)v;
		if (v < minploidy) minploidy = (byte)v;
		if (ploidy == 0) ploidy = v;
		if (maxploidy != minploidy) varploidy = true;

		ho += gt.HIndex() * c * v * (v - 1);
		how += c * v * (v - 1);

		n += c;
		v2i += c * v * (v - 1) / 2;

		if (!reassigned)
		{
			ushort* als = gt.GetAlleleArray();
			for (int i = 0; i < v; ++i)
				fre[als[i]] += c;
		}
	}
	else for (int gi = 0; gi < ngeno; ++gi)
	{
		GENOTYPE& gt = gtab[gi];
		if (gt.Nalleles() == 0 || gcount[gi] == 0) continue;
		uint c = gcount[gi];
		int v = gt.Ploidy();

		if (v > maxploidy) maxploidy = (byte)v;
		if (v < minploidy) minploidy = (byte)v;
		if (ploidy == 0) ploidy = v;
		if (maxploidy != minploidy) varploidy = true;

		ho += gt.HIndex() * c * v * (v - 1);
		how += c * v * (v - 1);

		n += c;
		v2i += c * v * (v - 1) / 2;

		if (!reassigned)
		{
			ushort* als = gt.GetAlleleArray();
			for (int i = 0; i < v; ++i)
				fre[als[i]] += c;
		}
	}

	if (maxploidy == 0) minploidy = maxploidy = 0;

	ptype = n / (REAL)ni;

	if (!reassigned)
		Unify(fre, k2);

	k = stat1.k = n == 0 ? 0 : (ushort)CountK(fre, k2);

	//set statistics to nan if k <= 2 to avoid involved in further average
	if (k < 2)
	{
		how = NAN;
		hew = NAN;
		piw = NAN;

		ho = NAN;
		he = NAN;
		pi = NAN;
		pic = NAN;
		ae = NAN;
		I = NAN;
		NE1P = NAN;
		NE2P = NAN;
		NEPP = NAN;
		NESI = NAN;
		NEID = NAN;
		g = NAN;
		df = NAN;
		pval = NAN;
		return;
	}

	ho = ho / how;
	bmaf = 1;
	he = 1;
	pi = 1;
	pic = 0;
	ae = 0;
	I = 0;

	if (ad == 0) stat1.nhaplo = nhaplo;
	else nhaplo = stat1.nhaplo;
	int nh = nhaplo + 0.5;

	REAL a2 = 0, a3 = 0, a4 = 0, a5 = 0, a6 = 0;
	for (int a = 0; a < k2; ++a)
	{
		REAL af = fre[a];
		if (af * nh < 1e-5) continue;
		if (af > 1e-5) I += -af * log(af);
		if (bmaf > af) bmaf = af;
		a2 += af * af;
		a3 += af * af * af;
		a4 += af * af * af * af;
		a5 += af * af * af * af * af;
		a6 += af * af * af * af * af * af;

		for (int b = a + 1; b < k2; ++b)
			if (fre[b] > 0)
				pic += 2 * af * fre[b] * (1 - af * fre[b]);
	}

	if (bmaf > 0.5) bmaf = 1 - bmaf;

	hew = nhaplo * nhaplo;
	piw = nhaplo * (nhaplo - 1);
	he = 1 - a2;
	pi = nhaplo * he / (nhaplo - 1);

	ae = 1 / a2;
	NE1P = 1 - (1 - 4 * a2 + 2 * a2 * a2 + 4 * a3 - 3 * a4);
	NE2P = 1 - (1 - 2 * a2 * a2 - 2 * a2 + a3 + 3 * a2 * a3 - 3 * a5 + 2 * a4);
	NEPP = 1 - (1 + 4 * a4 - 4 * a5 - 3 * a6 - 8 * a2 * a2 + 8 * a2 * a3 + 2 * a3 * a3);
	NESI = 1 - (0.75 - 0.5 * a2 - 0.5 * a2 * a2 + 0.25 * a4);
	NEID = 1 - (1 + 3 * a4 - 4 * a2 * a2);

	//chi test
	g = NAN;
	df = NAN;
	pval = NAN;

	// perform HWE test for locus diversity filter
	if (ad == 0 && !varploidy) // disable !haplotype on 2023-4-4
	{
		// initial mapping
		ploidy = maxploidy;
		MEMORY tlocus_memory;
		TABLE<HASH, uint> gitab(false, &tlocus_memory);

		//original
		LOCUS* locus1 = NULL;
		if (useslocus)
			locus1 = new (tlocus_memory.Alloc(sizeof(LOCUS))) LOCUS(tlocus_memory, slocus[l]);
		else
			locus1 = &locus[l];

		VLA_NEW(fre1, REAL, k2);					SetVal(fre1, fre, k2);
		VLA_NEW(gcount1, ushort, ngeno);			SetVal(gcount1, gcount, ngeno);
		VLA_NEW(fre2, REAL, k2);
		VLA_NEW(gcount2, ushort, ngeno);
		VLA_NEW(amap, int, k2);

		ushort alleles[N_MAX_PLOIDY] = { 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF };

		bool flag = false;
		for (int kk = k2; kk >= 1 && !flag; kk = (kk <= k ? kk - 1 : k))
		{
			if (kk == 1)
			{
				g = NAN;
				df = NAN;
				pval = NAN;
				break;
			}

			TABLE<HASH, GENOTYPE*>& gftab1 = locus1->gftab;
			int ngeno1 = gftab1.size;
			SetZero(gcount2, ngeno1);
			SetZero(fre2, kk);

			if (kk <= k)
			{
				// test locus1
				flag = true;
				g = 0;

				for (int gi = 0; gi < ngeno1; ++gi)
				{
					GENOTYPE& gt = *gftab1(gi);
					if (gt.Nalleles() == 0 || gcount1[gi] == 0) continue;

					REAL obs = gcount1[gi];
					REAL exp = gt.GFZ(f_model_val, fre1) * n;
					if (exp < 5)
					{
						flag = false;
						break;
					}
					g += 2 * obs * log(obs / exp);
				}
			}

			if (!flag)
			{
				// fail in test, collapse
				if (kk <= k)
				{
					// secondary map
					// find two minor alleles
					int mid1 = (int)GetMinID(fre1, kk);
					fre1[mid1]++;
					int mid2 = (int)GetMinID(fre1, kk);
					fre1[mid1]--;
					if (mid1 > mid2) Swap(mid1, mid2);

					// create allele map
					for (int j = 0, p = 0; j < kk; ++j)
					{
						if (j == mid2)
						{
							amap[j] = mid1;
							fre2[mid1] = fre1[mid1] + fre1[mid2];
						}
						else
						{
							amap[j] = p;
							fre2[p++] = fre1[j];
						}
					}
				}
				else
				{
					// initial map, exclude abscent alleles
					int mid1 = GetMaxID(fre1, kk);
					for (int j = 0, p = 0; j < kk; ++j)
						amap[j] = fre1[j] < MIN_FREQ ? -1 : p++;
					for (int j = 0; j < kk; ++j)
						amap[j] = amap[j] == -1 ? amap[mid1] : amap[j];
				}

				// create dummy locus
				TABLE<HASH, TEMP_GENOTYPE> tab2(true, &tlocus_memory);
				int gasize = 0;

				// create new alleles
				for (int gi = 0; gi < ngeno1; ++gi)
				{
					GENOTYPE* gt1 = gftab1(gi);
					if (gt1->Nalleles() == 0 || gcount1[gi] == 0) continue;

					ushort* als = gt1->GetAlleleArray();
					for (int ai = 0; ai < ploidy; ++ai)
						alleles[ai] = (ushort)amap[als[ai]];

					Sort(alleles, ploidy);//unphase

					for (int j = 0; j < ploidy; ++j)
						fre2[alleles[j]] += gcount1[gi];

					HASH hash2 = HashGenotype(alleles, ploidy);
					if (!tab2.ContainsKey(hash2))
					{
						TEMP_GENOTYPE& gt2 = tab2[hash2];
						SetVal(gt2.alleles, alleles, N_MAX_PLOIDY);
						gt2.hash = hash2;
						gt2.ploidy = ploidy;
						gt2.gid = tab2.size - 1;
						gcount2[gt2.gid] += gcount1[gi];
						gasize += ploidy + GetNalleles(alleles, ploidy);
					}
					else
					{
						TEMP_GENOTYPE& gt2 = tab2[hash2];
						gcount2[gt2.gid] += gcount1[gi];
					}
				}


				// create dummy locus
				int ngeno2 = tab2.size;
				LOCUS* locus2 = new(tlocus_memory.Alloc(sizeof(LOCUS)))
					LOCUS(tlocus_memory, *locus1, 1, ngeno2, gasize, tab2);//chi2 test

				// swap
				Unify(fre2, kk);
				locus1 = locus2;
				Swap(gcount1, gcount2);
				Swap(fre1, fre2);
			}
			else
			{
				// pass test
				df = (int)(Binomial(ploidy + kk - 1, ploidy) - kk + 0.5);
				pval = 1 - ChiSquareDistCDF(g, df);
				break;
			}
		}


		VLA_DELETE(fre1);
		VLA_DELETE(fre2);
		VLA_DELETE(gcount1);
		VLA_DELETE(gcount2);
		VLA_DELETE(amap);
	}
}

/* Write header to the result file */
template<typename REAL>
TARGET void DIVERSITY<REAL>::WriteHeader(FILE* fout)
{
	fprintf(fout, "%s%sPop", g_linebreak_val, g_linebreak_val);
	fprintf(fout, "%cLocus", g_delimiter_val);
	fprintf(fout, "%cPloidy", g_delimiter_val);
	fprintf(fout, "%ck", g_delimiter_val);
	fprintf(fout, "%cn", g_delimiter_val);
	fprintf(fout, "%c#Hap", g_delimiter_val);
	fprintf(fout, "%cHo", g_delimiter_val);
	fprintf(fout, "%cHe", g_delimiter_val);
	fprintf(fout, "%cpi", g_delimiter_val);
	fprintf(fout, "%cPIC", g_delimiter_val);
	fprintf(fout, "%cAe", g_delimiter_val);
	fprintf(fout, "%cI", g_delimiter_val);
	fprintf(fout, "%cNE1P", g_delimiter_val);
	fprintf(fout, "%cNE2P", g_delimiter_val);
	fprintf(fout, "%cNEPP", g_delimiter_val);
	fprintf(fout, "%cNEID", g_delimiter_val);
	fprintf(fout, "%cNESID", g_delimiter_val);
	fprintf(fout, "%cFis", g_delimiter_val);
	fprintf(fout, "%cG", g_delimiter_val);
	fprintf(fout, "%cd.f.", g_delimiter_val);
	fprintf(fout, "%cP-val", g_delimiter_val);
}

/* Write diversity of a locus to the result file */
template<typename REAL>
TARGET void DIVERSITY<REAL>::WriteLocus(FILE* fout, const char* name)
{
	//sum pop, sum locus
	char name_buf[NAME_BUF_LEN];
	fprintf(fout, "%s%s%c%s", g_linebreak_val, name, g_delimiter_val, GetLoc(l).GetNameStr(name_buf));
	fprintf(fout, "%c%d-%d",  g_delimiter_val, minploidy, maxploidy);
	fprintf(fout, "%c%d",     g_delimiter_val, k);
	fprintf(fout, "%c%d",     g_delimiter_val, n);
	fprintf(fout, "%c%d%c",   g_delimiter_val, nhaplo, g_delimiter_val);

	WriteReal(fout, ho); fprintf(fout, "%c", g_delimiter_val);
	WriteReal(fout, he); fprintf(fout, "%c", g_delimiter_val);
	WriteReal(fout, pi); fprintf(fout, "%c", g_delimiter_val);
	WriteReal(fout, pic); fprintf(fout, "%c", g_delimiter_val);
	WriteReal(fout, ae);  fprintf(fout, "%c", g_delimiter_val);
	WriteReal(fout, I);   fprintf(fout, "%c", g_delimiter_val);

	WriteReal(fout, NE1P); fprintf(fout, "%c", g_delimiter_val);
	WriteReal(fout, NE2P); fprintf(fout, "%c", g_delimiter_val);
	WriteReal(fout, NEPP); fprintf(fout, "%c", g_delimiter_val);
	WriteReal(fout, NEID); fprintf(fout, "%c", g_delimiter_val);
	WriteReal(fout, NESI); fprintf(fout, "%c", g_delimiter_val);

	WriteReal(fout, 1 - ho / he); fprintf(fout, "%c", g_delimiter_val);
	WriteReal(fout, g);   fprintf(fout, "%c", g_delimiter_val);
	if (IsError(df))   fprintf(fout, "0%c", g_delimiter_val);
	else               fprintf(fout, "%0.0lf%c", df, g_delimiter_val);
	WriteReal(fout, pval);
}

#endif

#define extern 
template<typename REAL>
extern DIVERSITY<REAL>* diversity_buf;				//Circle buffer for diversity estimation, NBUF
template<typename REAL>
extern DIVSUM<REAL> diversity_sum; 					//Diversity sum
extern int diversity_stage;							//Diversity level, 3 total, 2 pop, 1 reg
#undef extern 

/* Calculate genetic diveristy indices */
template<typename REAL>
TARGET void CalcDiversity()
{
	if (!diversity) return;
	EvaluationBegin();

	OpenResFile("-diversity", "Individual statistics");
	OpenTempFiles(2, ".diversity");
	diversity_buf<REAL> = new DIVERSITY<REAL>[NBUF];

	bool isfirst = true;
	int64 ntot = 0;
	if (diversity_level_val[1] || diversity_level_val[4])
		for (int i = 0; i < npop; ++i)
			ntot += apops<REAL>[i]->nind * nloc;

	if (diversity_level_val[2] || diversity_level_val[5])
		for (int rl = lreg - 1; rl >= 0; --rl)
			for (int i = 0; i < nreg[rl]; ++i)
				ntot += aregs<REAL>[rl][i]->nind * nloc;

	if (diversity_level_val[3] || diversity_level_val[6])
		ntot += nind * nloc;

	//Total population
	if (diversity_level_val[3] || diversity_level_val[6])
	{
		diversity_stage = 3;
		cpop<REAL> = total_pop<REAL>;

		if (diversity_level_val[3]) DIVERSITY<REAL>::WriteHeader(TEMP_FILES[0]);
		if (diversity_level_val[6]) DIVERSITY<REAL>::WriteHeader(TEMP_FILES[1]);

		RunThreads(&DiversityThread<REAL>, &DiversityGuard<REAL>, NULL, ntot, cpop<REAL>->nind * nloc,
			"\nCalculating genetic diversity:\n", g_nthread_val, isfirst);
		isfirst = false;

		if (diversity_level_val[3])
			diversity_sum<REAL>.Write(TEMP_FILES[0], "Total");
        
        diversity_sum<REAL>.UnInit();
	}

	//Regions
	if (diversity_level_val[2] || diversity_level_val[5])
	{
		diversity_stage = 2;
		if (diversity_level_val[2]) DIVERSITY<REAL>::WriteHeader(TEMP_FILES[0]);
		if (diversity_level_val[5]) DIVERSITY<REAL>::WriteHeader(TEMP_FILES[1]);

		for (int rl = lreg - 1; rl >= 0; --rl)
			for (int i = 0; i < nreg[rl]; ++i)
			{
				cpop<REAL> = aregs<REAL>[rl][i];
				RunThreads(&DiversityThread<REAL>, &DiversityGuard<REAL>, NULL, ntot, cpop<REAL>->nind * nloc,
					"\nCalculating genetic diversity:\n", g_nthread_val, isfirst);
				isfirst = false;

				if (diversity_level_val[2])
					diversity_sum<REAL>.Write(TEMP_FILES[0], cpop<REAL>->name);
                
                diversity_sum<REAL>.UnInit();
			}
	}

	//Populations
	if (diversity_level_val[1] || diversity_level_val[4])
	{
		diversity_stage = 1;
		if (diversity_level_val[1]) DIVERSITY<REAL>::WriteHeader(TEMP_FILES[0]);
		if (diversity_level_val[4]) DIVERSITY<REAL>::WriteHeader(TEMP_FILES[1]);

		for (int i = 0; i < npop; ++i)
		{
			cpop<REAL> = apops<REAL>[i];
			RunThreads(&DiversityThread<REAL>, &DiversityGuard<REAL>, NULL, ntot, cpop<REAL>->nind * nloc,
				"\nCalculating genetic diversity:\n", g_nthread_val, isfirst);
			isfirst = false;

			if (diversity_level_val[1])
				diversity_sum<REAL>.Write(TEMP_FILES[0], cpop<REAL>->name);
            
            diversity_sum<REAL>.UnInit();
		}
	}

	JoinTempFiles(2);
	CloseResFile();
	CheckGenotypeId<REAL>();

	DEL(diversity_buf<REAL>);
	EvaluationEnd("Genetic diversity estimation");
}

/* Add and write genetic diversity */
THREAD2(DiversityGuard)
{
	int ni = cpop<REAL>->nind;
	diversity_sum<REAL>.Init();

	bool addsum = diversity_level_val[diversity_stage];
	bool writelocus = diversity_level_val[diversity_stage + 3];
	char* popname = cpop<REAL> == total_pop<REAL> ? (char*)"Total" : cpop<REAL>->name;

	for (int64& ii = progress1 = 0; ii < nloc; ii++)
	{
		GUARD_BEGIN

		PROGRESS_VALUE += ni;

		if (addsum)
			diversity_sum<REAL>.Add(diversity_buf<REAL>[ii % NBUF]);

		if (writelocus)
			diversity_buf<REAL>[ii % NBUF].WriteLocus(TEMP_FILES[1], popname);

		GUARD_END
	}
}

/* Calculate genetic diversity using multiple threads */
THREAD2(DiversityThread)
{
	for (int64 ii = 0; ii < nloc; ii++)
	{
		THREAD_BEGIN

		diversity_buf<REAL>[ii % NBUF].CalcDiversity(ii);

		THREAD_END
	}
}
