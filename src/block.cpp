/* Linkage Disequilibrium Decay Functions */

#include "vcfpop.h"

template TARGET void BLOCK_PAIR<double>::WriteHeader(FILE* fout);
template TARGET void BLOCK_PAIR<float >::WriteHeader(FILE* fout);
template TARGET void BLOCK_PAIR<double>::Write(FILE* fout);
template TARGET void BLOCK_PAIR<float >::Write(FILE* fout);

template TARGET void CalcBlock<double>();
template TARGET void CalcBlock<float >();

#ifndef _BLOCK_PAIR

/* Write table header */
template<typename REAL>
TARGET void BLOCK_PAIR<REAL>::WriteHeader(FILE* fout)
{
	fprintf(fout, "%s%sChrom%cID1%cPos1%cID2%cPos2%cnpair", g_linebreak_val, g_linebreak_val, g_delimiter_val, g_delimiter_val, g_delimiter_val, g_delimiter_val, g_delimiter_val);

	if (block_estimator_val[1]) fprintf(fout, "%cr2%cSE",		g_delimiter_val, g_delimiter_val);
	if (block_estimator_val[2]) fprintf(fout, "%cr2Delta%cSE",	g_delimiter_val, g_delimiter_val);
	if (block_estimator_val[3]) fprintf(fout, "%cD'%cSE",		g_delimiter_val, g_delimiter_val);
	if (block_estimator_val[4]) fprintf(fout, "%cDelta'%cSE",	g_delimiter_val, g_delimiter_val);

	fprintf(fout, "%s", g_linebreak_val);
}

/* Write table entry */
template<typename REAL>
TARGET void BLOCK_PAIR<REAL>::Write(FILE* fout)
{
	if (this->npair)
	{
		if (block_estimator_val[1])
		{
			this->r2A /= this->r2B;
			this->r2A2 = sqrt((this->r2A2 / this->r2B - this->r2A * this->r2A) / (1 - this->r2B2 / (this->r2B * this->r2B)) * this->r2B2 / (this->r2B * this->r2B));
			//this->r2B = this->r2A / this->r2A2;
			//this->r2B2 = MinusLogPNormal(this->r2B);
		}

		if (block_estimator_val[2])
		{
			this->r2DeltaA /= this->r2DeltaB;
			this->r2DeltaA2 = sqrt((this->r2DeltaA2 / this->r2DeltaB - this->r2DeltaA * this->r2DeltaA) / (1 - this->r2DeltaB2 / (this->r2DeltaB * this->r2DeltaB)) * this->r2DeltaB2 / (this->r2DeltaB * this->r2DeltaB));
			//this->r2DeltaB = this->r2DeltaA / this->r2DeltaA2; 
			//this->r2DeltaB2 = MinusLogPNormal(this->r2DeltaB);
		}

		if (block_estimator_val[3])
		{
			this->DpA /= this->DpB;
			this->DpA2 = sqrt((this->DpA2 / this->DpB - this->DpA * this->DpA) / (1 - this->DpB2 / (this->DpB * this->DpB)) * this->DpB2 / (this->DpB * this->DpB));
			//this->DpB = this->DpA / this->DpA2;
			//this->DpB2 = MinusLogPNormal(this->DpB);
		}

		if (block_estimator_val[4])
		{
			this->DeltapA /= this->DeltapB;
			this->DeltapA2 = sqrt((this->DeltapA2 / this->DeltapB - this->DeltapA * this->DeltapA) / (1 - this->DeltapB2 / (this->DeltapB * this->DeltapB)) * this->DeltapB2 / (this->DeltapB * this->DeltapB));
			//this->DeltapB = this->DeltapA / this->DeltapA2;
			//this->DeltapB2 = MinusLogPNormal(this->DeltapB);
		}
	}

	fprintf(fout, "%s%c%d%c%d%c", block1->chr, g_delimiter_val, block1->id + 1, g_delimiter_val, block1->pos, g_delimiter_val);
	fprintf(fout, "%d%c%d%c",                                   block2->id + 1, g_delimiter_val, block2->pos, g_delimiter_val);
	fprintf(fout, "%d", this->npair);

	//Weighted average
	if (block_estimator_val[1])
	{
		if (this->npair == 0)
			fprintf(fout, "%cnan%cnan", g_delimiter_val, g_delimiter_val);
		else
		{
			fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, this->r2A);
			fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, this->r2A2);
			//fprintf(fout, "%c", g_delimiter_val);
			//WriteReal(fout, this->r2B2);
		}
	}

	if (block_estimator_val[2])
	{
		if (this->npair == 0)
			fprintf(fout, "%cnan%cnan", g_delimiter_val, g_delimiter_val);
		else
		{
			fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, this->r2DeltaA);
			fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, this->r2DeltaA2);
			//fprintf(fout, "%c", g_delimiter_val);
			//WriteReal(fout, this->r2DeltaB2);
		}
	}

	if (block_estimator_val[3])
	{
		if (this->npair == 0)
			fprintf(fout, "%cnan%cnan", g_delimiter_val, g_delimiter_val);
		else
		{
			fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, this->DpA);
			fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, this->DpA2);
			//fprintf(fout, "%c", g_delimiter_val);
			//WriteReal(fout, this->DpB2);
		}
	}

	if (block_estimator_val[4])
	{
		if (this->npair == 0)
			fprintf(fout, "%cnan%cnan", g_delimiter_val, g_delimiter_val);
		else
		{
			fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, this->DeltapA);
			fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, this->DeltapA2);
			//fprintf(fout, "%c", g_delimiter_val);
			//WriteReal(fout, this->DeltapB2);
		}
	}
	fprintf(fout, "%s", g_linebreak_val);
}
#endif

template<> BLOCK_SINGLE<double>*	block_singles<double>;
template<> BLOCK_SINGLE<float >*	block_singles<float >;
template<> BLOCK_PAIR<double>*		block_pairs<double>;
template<> BLOCK_PAIR<float >*		block_pairs<float >;
vector<char*>						block_chromas;
umap<HASH, CHROM_PROP>				block_chrom_sted;
int									block_npair;
atomic<int>							block_cpair;
constexpr int						block_maxnloc2 = 256;

/* Calculate LD block */
template<typename REAL>
TARGET void CalcBlock()
{
	if (!block) return;
	if (ad) Exit("\nError: LD block (-block) is incompatible with allelic depth (-ad) option.\n");
	if (abs(g_format_val) > BCF) Exit("\nError: LD block (-block) function can only be used for VCF/BCF input files.\n");

	EvaluationBegin();

	OpenResFile("-block", "LD block");

	// assign cpop
	AssignPop<REAL>(block_pop_b, block_pop_val, "-block_pop");

	RunThreads(&DecayFreqThread<REAL>, NULL, NULL, nloc * cpop<REAL>->nind, nloc * cpop<REAL>->nind,
		"\nPreparing allele frequency:\n", 1, true);

	// add chromosomes
	CHROM_PROP def_prop{ 0, 0xFFFFFFFF, 0, 0xFFFFFFFF };
	for (int64 l = 0; l < nloc; ++l)
	{
		char* chr = GetLoc(l).GetChrom();
		HASH ha = HashString(chr, (int)strlen(chr));

		if (block_chrom_sted.find(ha) == block_chrom_sted.end())
		{
			block_chrom_sted[ha] = def_prop;
			block_chromas.push_back(chr);
		}

		CHROM_PROP& prop = block_chrom_sted[ha];

		prop.min = std::min(GetLocPos(l), prop.min);
		prop.max = std::max(GetLocPos(l), prop.max);
		prop.st =  std::min((uint64)l, prop.st);
		prop.ed =  std::max((uint64)l, prop.ed);
	}

	//add singles
	int nsingles = 0;
	{
		for (char* chr : block_chromas)
		{
			HASH ha = HashString(chr, (int)strlen(chr));
			CHROM_PROP& prop = block_chrom_sted[ha];
			nsingles += (int)ceil((prop.max - prop.min + 1) / (double)block_size_val);
		}

		block_singles<REAL> = new BLOCK_SINGLE<REAL>[nsingles];
		for (int chrid = 0, psingle = 0; chrid < block_chromas.size(); ++chrid)
		{
			char* chr = block_chromas[chrid];
			HASH ha = HashString(chr, (int)strlen(chr));
			CHROM_PROP& prop = block_chrom_sted[ha];
			int nsingle_chr = (int)ceil((prop.max - prop.min + 1) / (double)block_size_val);

			uint pos_st = prop.min, pos_ed = prop.min + block_size_val - 1;

			for (int i = 0, l = prop.st; i < nsingle_chr; ++i)
			{
				BLOCK_SINGLE<REAL> tsingle{ chr, -1, -1, -1, -1 };

				uint pos1 = GetLocPos(l);
				if (pos_st <= pos1 && pos1 <= pos_ed)
				{
					tsingle.loc_st = l;
					for (; l <= prop.ed && GetLocPos(l) <= pos_ed; ++l);
					tsingle.loc_ed = l - 1;
				}

				tsingle.id = psingle;
				tsingle.pos = pos_st;
				block_singles<REAL>[psingle++] = tsingle;
				pos_st += block_size_val;
				pos_ed += block_size_val;
			}
		}
	}

	//count block pairs
	block_npair = 0;  block_cpair = 0;
	for (int chrid = 0, psingle = 0; chrid < block_chromas.size(); ++chrid)
	{
		char* chr = block_chromas[chrid];
		HASH ha = HashString(chr, (int)strlen(chr));
		CHROM_PROP& prop = block_chrom_sted[ha];
		int nsingle_chr = (int)ceil((prop.max - prop.min + 1) / (double)block_size_val);

		for (int s1 = 0; s1 < nsingle_chr; ++s1)
			for (int s2 = s1 + 1; s2 < nsingle_chr && s2 - s1 <= block_maxd_val; ++s2)
				block_npair++;

		psingle += nsingle_chr;
	}

	//add block pairs
	block_pairs<REAL> = new BLOCK_PAIR<REAL>[block_npair];
	for (int chrid = 0, psingle = 0, ppair = 0; chrid < block_chromas.size(); ++chrid)
	{
		char* chr = block_chromas[chrid];
		HASH ha = HashString(chr, (int)strlen(chr));
		CHROM_PROP& prop = block_chrom_sted[ha];
		int nsingle_chr = (int)ceil((prop.max - prop.min + 1) / (double)block_size_val);

		for (int s1 = 0; s1 < nsingle_chr; ++s1)
		{
			for (int s2 = s1 + 1; s2 < nsingle_chr && s2 - s1 <= block_maxd_val; ++s2)
			{
				BLOCK_PAIR<REAL>& pair = block_pairs<REAL>[ppair++];
				pair.block1 = &block_singles<REAL>[psingle + s1];
				pair.block2 = &block_singles<REAL>[psingle + s2];
			}
		}
		psingle += nsingle_chr;
	}

	RunThreads(&BlockThread<REAL>, NULL, NULL, block_npair, block_npair, "\nCalculating LD block:\n", g_nthread_val, true);

	BLOCK_PAIR<REAL>::WriteHeader(FRES);

	for (int i = 0; i < block_npair; ++i)
		block_pairs<REAL>[i].Write(FRES);

	CloseResFile();

	DEL(block_singles<REAL>);
	DEL(block_pairs<REAL>);

	umap<HASH, CHROM_PROP>().swap(block_chrom_sted);
	vector<char*>().swap(block_chromas);

	cpop<REAL>->UnAllocFreq();
	DEL(allele_freq_offset);
	DEL(genotype_count_offset);

	EvaluationEnd("LD block analysis");

	if (block_plot_val == 1)
		RunRscript("block_plot.R");
}

/* Calculate LD block using multiple threads */
THREAD2(BlockThread)
{
	//allocate intervals
	TABLE<uint64, int> genopair[block_maxnloc2];
	REP(block_maxnloc2) new (&genopair[kk]) TABLE<uint64, int>(false, NULL, 64);
	int loc2[block_maxnloc2], * buf = new int[maxK * maxK + maxK * 4];
	
	for (int pairid = block_cpair.fetch_add(1); pairid < block_npair; pairid = block_cpair.fetch_add(1))
	{
		BLOCK_PAIR<REAL>& pair = block_pairs<REAL>[pairid];
		BLOCK_SINGLE<REAL>& block1 = *block_pairs<REAL>[pairid].block1;
		BLOCK_SINGLE<REAL>& block2 = *block_pairs<REAL>[pairid].block2;

		for (uint l1 = block1.loc_st; l1 <= block1.loc_ed; ++l1)
		{
			DECAY<REAL> decay[block_maxnloc2];
			RNG<double> rng(g_seed_val + l1, RNG_SALT_LDBLOCK);

			for (uint l2 = block2.loc_st, nloc2 = 0; l2 <= block2.loc_ed; ++l2)
			{
				if (block_ratio_val == 1 || rng.Uniform() < block_ratio_val)
					loc2[nloc2++] = l2;

				if (nloc2 == block_maxnloc2 || l2 == block2.loc_ed && nloc2)
				{
					decay[0].CalcLD(l1, loc2, buf, nloc2, genopair);
					for (int j = 0; j < nloc2; ++j)
						pair.AddDecay(decay[j]);
					nloc2 = 0;
				}
			}
		}

		PROGRESS_VALUE++;
	}

	DEL(buf);
}

