/* Filter Functions */

#pragma once
#include "vcfpop.h"

#ifndef _PCOA
#endif

#define extern 
extern bool diversity_filter;
extern bool genotype_filter;
extern bool individual_filter;
extern bool info_filter;
extern TABLE<HASH, TABLE<HASH, pair<LOCN, double>>> contig_diversity; //Max diversity indices in a shift-window
extern LOCN nfilter;								//Number of filtered loci
extern atomic<int> nfilterind;						//Number of filtered individuals
#undef extern 

/* Applying individual and diversity filters */
TARGET void ApplyFilter()
{
	EvaluationBegin();
	CheckGenotypeId();

	//1. Applying individual filter
	if (individual_filter)
	{
		RunThreads(&MarkerIndividual, NULL, NULL, nind + nloc * 2, nind + nloc,
			"\nApplying individual filter:\n", 1, true, (int)((nind + nloc) / (nind + nloc * 2) * g_progress_val));

		RunThreads(&RemoveIndividual, NULL, NULL, nind + nloc * 2, nloc,
			NULL, 1, false, g_progress_val - (int)((nind + nloc) / (nind + nloc * 2) * g_progress_val));
	}
	CheckGenotypeId();

	//Assign individuals into their vpop
	AssignVInds();
	CheckGenotypeId();

	//2. Apply diversity filter
	if (diversity_filter)
	{
		if (f_pop_b)
		{
			if ("total" == f_pop_val)
				cpop = total_pop;
			else
			{
				bool find = false;
				for (int i = 0; i < npop; ++i)
					if (pop[i].name == f_pop_val)
					{
						find = true;
						cpop = &pop[i];
					}

				if (!find) Exit("\nError: Cannot find target population %d, check parameter -f_pop.\n", f_pop_val.c_str());
			}
		}
		else if (f_region_b)
		{
			bool find = false;
			for (uint rl = 0; rl < reg.size - 1; ++rl)
				for (uint i = 0; i < reg[rl].size; ++i)
					if (reg[rl][i].name == f_region_val)
					{
						find = true;
						cpop = &reg[rl][i];
					}

			if (!find) Exit("\nError: Cannot find target region %s, check parameter -f_region.\n", f_region_val.c_str());
		}
		else
			cpop = total_pop;

		//Temp use allele frequency and genotype count for diversity calculation
		allele_freq_offset = new LOCN[nloc];
		genotype_count_offset = new LOCN[nloc];
		SetFF(allele_freq_offset, nloc);
		SetFF(genotype_count_offset, nloc);

		//Calculate offset
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

		//Calculate freq for cpop and diversity estimation
		cpop->AllocFreq();

		if (f_windowsize_b && abs(g_format_val) <= BCF)
			new(&contig_diversity) TABLE<HASH, TABLE<HASH, pair<LOCN, double>>>(true, NULL);

		RunThreads(&MarkerDiversity, NULL, NULL, nloc, nloc, "\nMarking locus diversity filter:\n", 1, true);

		if (f_windowsize_b && abs(g_format_val) <= BCF)
		{
			for (int64 l = 0; l < nloc; ++l)
				GetLoc(l).flag_pass = false;

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
			for (int c = 0; c < contig_diversity.size; ++c)
			{
				TABLE<HASH, pair<LOCN, double>>& contig = contig_diversity(c);
				for (uint w = 0; w < contig.size; ++w)
					GetLoc(contig(w).first).flag_pass = true;
				contig.~TABLE();
			}
			contig_diversity.~TABLE();//ok
			new (&contig_diversity) TABLE<HASH, TABLE<HASH, pair<LOCN, double>>>;
		}
		cpop->UnAllocFreq();

		//Release temp allele frequency and genotype count
		delete[] allele_freq_offset;
		delete[] genotype_count_offset;
		allele_freq_offset = genotype_count_offset = NULL;
	}

	CheckGenotypeId();

	//Mark monomorphic locus after previous filtering
	TABLE<int, int>* allele_table = new TABLE<int, int>[g_nthread_val];
	for (int i = 0; i < g_nthread_val; ++i)
		new(&allele_table[i]) TABLE<int, int>(true, NULL);

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
	for (int64 l = 0; l < nloc; ++l)
	{
		if (GetLoc(l).k < 2)
			GetLoc(l).flag_pass = false;
		else if (GetLoc(l).flag_pass == true)
		{
			threadid = omp_get_thread_num();
			TABLE<int, int>& atab = allele_table[threadid];
			atab.Clear();

			GENOTYPE* gtab = GetLoc(l).GetGtab();
			int ngeno = GetLoc(l).ngeno;
			bool pass_flag = false;

			for (int gi = 0; gi < ngeno; ++gi)
			{
				GENOTYPE& gt = gtab[gi];
				int nalleles = gt.Nalleles();
				if (nalleles == 0) continue;
				if (nalleles >= 2) { pass_flag = true; break; }
				atab[gt.GetAlleleArray()[0]] = 1;
				if (atab.size > 1) { pass_flag = true; break; }
			}

			GetLoc(l).flag_pass = pass_flag;
		}
	}
	delete[] allele_table;
	CheckGenotypeId();

	//3. Apply locus filter
	nfilter = 0;
	KT = GT = maxK = maxG = 0;
	for (int64 l = 0; l < nloc; ++l)
	{
		if (GetLoc(l).flag_pass)
		{
			nfilter++;
			KT += GetLoc(l).k;
			GT += GetLoc(l).ngeno;
			maxK = Max((int)GetLoc(l).k, maxK);
			maxG = Max((int)GetLoc(l).ngeno, maxG);
		}
	}

	CheckGenotypeId();

	if (nfilter == 0)
		Exit("\nError: all %d loci are excluded from analysis, try less strict filters.\n", nloc);

	if (nfilter != nloc)
		RunThreads(&RemoveLocus, NULL, NULL, nloc * (ad_bucket.base_addr ? 2 : 1), nloc * (ad_bucket.base_addr ? 2 : 1), "\nApplying locus diversity filter:\n", 1, true);

	CheckGenotypeId();

	if (locus_pos && !haplotype)
	{
		delete[] locus_pos;
		locus_pos = NULL;
	}

	EvaluationEnd("Filter genotypes, individuals and loci");
}

/* Recursive set vid, vpop, ind for each region */
TARGET void SetVReg(int rl, int i)
{
	if (rl >= 0)
	{
		if (reg[rl][i].npop == 0 || reg[rl][i].nind == 0) return;
		reg[rl][i].id = nreg[rl];
		reg[rl][i].vpop = rl >= 1 ? aregs[rl - 1] + nreg[rl - 1] : apops + npop;
		reg[rl][i].inds = rinds + nind;
		reg[rl][i].ind0id = nind;
		aregs[rl][nreg[rl]++] = &reg[rl][i];
	}

	if (rl > 0) for (uint j = 0; j < reg[rl - 1].size; ++j)
	{
		if (reg[rl - 1][j].rid != i) continue;
		reg[rl - 1][j].rid = reg[rl][i].id;
		SetVReg(rl - 1, j);
	}
	else for (uint j = 0; j < pop.size; ++j)
	{
		if (pop[j].rid != i || pop[j].nind == 0) continue;

		apops[npop] = &pop[j];
		pop[j].id = npop++;
		pop[j].rid = reg[rl][i].id;
		pop[j].vpop = NULL;
		pop[j].inds = rinds + nind;
		pop[j].ind0id = nind;
		nind += pop[j].nind;
	}
}

/* Sort individuals by population index to rinds array */
TARGET void AssignVInds()
{
	/* Original order is in ind and should not be changed to correctly obtain genotype index */

	//Assign pop/reg
	for (uint i = 0; i < pop.size; ++i)
		pop[i].nind = 0;

	for (uint rl = 0; rl < reg.size; ++rl)
		for (uint i = 0; i < reg[rl].size; ++i)
			reg[rl][i].npop = reg[rl][i].nind = 0;

	//First scan, count inds
	for (int j = 0; j < nind; ++j)
	{
		//default pop or input pop
		int i = (ainds[j]->popid < 1 || ainds[j]->popid >= pop.size) ? 0 : ainds[j]->popid;
		pop[i].nind++;
	}

	//Allocate vreg
	nregt = nregt2 = 0;
	lreg = reg.size;
	aregs = new POP * *[reg.size];
	SetZero(aregs, reg.size);
	SetZero(nreg, lreg);

	//Add up nind and npop to regions at lay 0
	for (uint i = 0; i < pop.size; ++i)
	{
		if (pop[i].nind == 0) continue;
		if (reg.size > 0 && reg[0].size > 0)
		{
			reg[0][pop[i].rid].npop++;
			reg[0][pop[i].rid].nind += pop[i].nind;
		}
	}

	//Add up nind and npop to regions at lay 1+
	for (uint rl = 0; rl < reg.size; ++rl)
	{
		for (uint i = 0; i < reg[rl].size; ++i)
		{
			if (reg[rl][i].npop > 0)
			{
				nreg[rl]++;  nregt++;
				if (rl < reg.size - 1)
				{
					reg[rl + 1][reg[rl][i].rid].npop++;
					reg[rl + 1][reg[rl][i].rid].nind += reg[rl][i].nind;
				}
			}
		}
		nregt2 += nreg[rl] * nreg[rl];
	}

	//Remove top lay while there is only one region in this lay
	while (lreg && nreg[lreg - 1] == 1)
	{
		nregt--;
		nregt2--;
		lreg--;
		if (lreg == 0) break;
	}

	//Allocate vpop
	npop = 0;
	for (uint i = 0; i < pop.size; ++i)
		if (pop[i].nind > 0)
			npop++;
	apops = new POP * [npop];

	//Contruct total pop
	if (lreg == 0 && npop == 1)
	{
		lreg = -1;
		for (uint i = 0; i < pop.size; ++i)
			if (pop[i].nind > 0)
				total_pop = &pop[i];
		nregt2 = nregt = 0;
	}
	else
	{
		for (uint i = 0; i < reg[lreg].size; ++i)
			if (reg[lreg][i].nind > 0)
				total_pop = &reg[lreg][i];
	}

	delete[] total_pop->name;
	total_pop->name = new char[6];
	sprintf(total_pop->name, "Total");
	total_pop->rid = -1;

	for (int rl = 0; rl <= lreg; ++rl)
	{
		aregs[rl] = new POP * [nreg[rl]];
		nreg[rl] = 0;
	}

	//Allocate rinds
	rinds = new IND * [nind];

	//Recursive arrange vpop, vreg and vind
	if (lreg == -1)
	{
		apops[0] = total_pop;
		total_pop->id = 0;
		total_pop->rid = -1;
		total_pop->vpop = NULL;
		total_pop->inds = rinds;
		total_pop->ind0id = 0;
	}
	else
	{
		nind = 0; npop = 0;
		SetVReg(lreg, total_pop->id);
	}

	//Set nind = 0
	for (uint i = 0; i < pop.size; ++i)
		pop[i].nind = 0;

	//Move ind to pop ind array and assign vpopid
	bool swap = false;
	for (int j = 0; j < nind; ++j)
	{
		//Default pop or input pop
		IND& i = *ainds[j];
		POP& p = pop[(i.popid < 1 || i.popid >= pop.size) ? 0 : i.popid];
		if (i.indid != p.ind0id + p.nind)
			swap = true;
		i.indid = p.ind0id + p.nind;
		p.inds[p.nind++] = ainds[j];
		i.popid = p.id;
	}

	if (!swap) return;

	//swap genotype
	ushort* gtbuf = new ushort[g_nthread_val * nind];

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
	for (int64 l = 0; l < nloc; ++l)
	{
		threadid = omp_get_thread_num();
		ushort* gtb = gtbuf + threadid * nind;

		GENO_READER rg(0, l);
		GENO_WRITER wg(l);

		for (int i = 0; i < nind; ++i)
			gtb[ainds[i]->indid] = (ushort)rg.Read();

		for (int i = 0; i < nind; ++i)
			wg.Write(gtb[i]);

		wg.FinishWrite();
	}

	delete[] gtbuf;

	if (ploidyinfer)
	{
		maxK = 0;
		for (int64 l = 0; l < nloc; ++l)
			maxK = Max((int)GetLoc(l).k, maxK);

		//swap allelic depth
		uint* dpbuf = new uint[g_nthread_val * nind * maxK];

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int64 l = 0; l < nloc; ++l)
		{
			threadid = omp_get_thread_num();
			uint* dpb = dpbuf + threadid * nind * maxK;

			GENO_READER rd(0, l, &ad_bucket);
			GENO_WRITER wd(l, &ad_bucket);
			int k2 = GetLoc(l).k;

			for (int i = 0; i < nind; ++i)
			{
				uint* dpb2 = dpb + ainds[i]->indid * k2;

				for (int k = 0; k < k2; ++k)
					dpb2[k] = rd.Read();
			}

			for (int i = 0, nend = k2 * nind; i < nend; ++i)
				wd.Write(*dpb++);
		}

		delete[] dpbuf;
	}

	delete[] ainds; ainds = rinds;
}

/* Marker individual filtered or not */
THREAD(MarkerIndividual)
{
	//Calculate offset
	maxG = 0;
	for (int64 l = 0; l < nloc; ++l)
		maxG = Max((int)GetLoc(l).ngeno, maxG);

	VLA_NEW(ploidytab, int, maxG * g_nthread_val);
	VLA_NEW(nallelestab, int, maxG * g_nthread_val);

	atomic<int64>* ntype = new atomic<int64>[nind];
	atomic<int64>* vt = new atomic<int64>[nind];
	atomic<byte>* vmin = new atomic<byte>[nind];
	atomic<byte>* vmax = new atomic<byte>[nind];

	SetZero(vt, nind);
	SetZero(ntype, nind);
	SetVal((byte*)vmin, (byte)0, nind);
	SetZero(vmax, nind);
	CheckGenotypeId();

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
	for (int64 l = 0; l < nloc; ++l)
	{
		threadid = omp_get_thread_num();

		int* ploidy = ploidytab + maxG * threadid;
		int* nalleles = nallelestab + maxG * threadid;

		GENOTYPE* gtab = GetLoc(l).GetGtab();
		int ngeno = GetLoc(l).ngeno;

		for (int gi = 0; gi < ngeno; ++gi)
		{
			GENOTYPE& gt = gtab[gi];
			nalleles[gi] = gt.Nalleles();
			ploidy[gi] = gt.Ploidy();
		}

		GENO_READER rt(0, l);
		for (int j = 0; j < nind; ++j)
		{
			int gid = rt.Read();
			if (nalleles[gid])
			{
				byte v = (byte)ploidy[gid];
				vt[j] += v;
				ntype[j]++;
				if (v < vmin[j]) AtomicMin(vmin[j], v);
				if (v > vmax[j]) AtomicMax(vmax[j], v);
			}
		}

		PROGRESS_VALUE++;
	}

	int64 f_ntype_min = (int64)(f_itype_min * nloc + 0.5);
	int64 f_ntype_max = (int64)(f_itype_max * nloc + 0.5);

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
	for (int i = 0; i < nind; ++i)
	{
		if (ainds[i] && f_itype_b && (f_ntype_min > ntype[i] || ntype[i] > f_ntype_max))
			ainds[i] = NULL;

		if (ainds[i] && f_iploidy_b && (f_iploidy_min > minploidy || minploidy > f_iploidy_max))
			ainds[i] = NULL;

		if (ainds[i] && f_iploidy_b && (f_iploidy_min > maxploidy || maxploidy > f_iploidy_max))
			ainds[i] = NULL;

		PROGRESS_VALUE++;
	}

	delete[] ntype;
	delete[] vt;
	delete[] vmin;
	delete[] vmax;

	VLA_DELETE(ploidytab);
	VLA_DELETE(nallelestab);
}

/* Marker locus filtered or not */
THREAD(MarkerDiversity)
{
	bool flag_window = f_windowsize_b && abs(g_format_val) <= BCF;
#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
	for (int64 l = 0; l < nloc; ++l)
	{
		DIVERSITY d;
		if (GetLoc(l).flag_pass)
		{
			d.CalcDiversity(l);

			if ((d.k == 0) ||
				(f_bmaf_b && d.k == 2 && (d.bmaf < f_bmaf_min || d.bmaf > f_bmaf_max)) ||
				(f_k_b && (d.k < f_k_min || d.k > f_k_max)) ||
				(f_n_b && (d.n < f_n_min || d.n > f_n_max)) ||
				(f_ptype_b && (d.ptype < f_ptype_min || d.ptype > f_ptype_max)) ||
				(f_pval_b && d.pval != -1 && (d.pval < f_pval_min || d.pval > f_pval_max)) ||
				(f_he_b && (d.he < f_he_min || d.he > f_he_max)) ||
				(f_ho_b && (d.ho < f_ho_min || d.ho > f_ho_max)) ||
				(f_pic_b && (d.pic < f_pic_min || d.pic > f_pic_max)) ||
				(f_ae_b && (d.ae < f_ae_min || d.ae > f_ae_max)) ||
				(f_I_b && (d.I < f_I_min || d.I > f_I_max)))
				GetLoc(l).flag_pass = false;

			//pass diversity filter, write its diversity stat to contig_diversity
			else if (flag_window)
			{
				double div = 0;
				switch (f_windowstat_val)
				{
				case 1: div = d.bmaf; break;
				case 2: div = d.k; break;
				case 3: div = d.n; break;
				case 4: div = d.ptype; break;
				case 5: div = d.he; break;
				case 6: div = d.ho; break;
				case 7: div = d.pic; break;
				case 8: div = d.ae; break;
				case 9: div = d.I; break;
				default: div = d.he; break;
				}

				HASH ha = HashString(GetLoc(l).GetChrom());
				if (!contig_diversity.ContainsKey(ha))
				{
					Lock(GLOCK2);
					if (!contig_diversity.ContainsKey(ha))
						new(&contig_diversity[ha]) TABLE<HASH, pair<LOCN, double>>(true, NULL);
					UnLock(GLOCK2);
				}

				TABLE<HASH, pair<LOCN, double>>& contig = contig_diversity[ha];
				while ((int64)contig.bucket == -1);

				LOCN w = (GetLocPos(l) - 1) / f_windowsize_val;
				if (!contig.ContainsKey(w))
				{
					Lock(GLOCK3[ha % 256]);
					if (!contig.ContainsKey(w))
						contig[w] = pair<LOCN, double>{ l, div };
					UnLock(GLOCK3[ha % 256]);
				}

				pair<LOCN, double>& window = contig[w];
				if (window.second < div)
				{
					window.first = l;
					window.second = div;
				}
			}

		}
		else
			GetLoc(l).flag_pass = false;

		PROGRESS_VALUE++;
	}
}

/* Remove individual fail to pass filter */
THREAD(RemoveIndividual)
{
	int nfind = 0;
	for (int i = 0; i < nind; ++i)
		if (ainds[i]) nfind++;

	if (nfind == nind)
	{
		PROGRESS_VALUE = PROGRESS_CEND;
		return;
	}

	if (nfind == 0)
		Exit("\nError: all %d individuals are excluded from analysis due to the aneuploid in the same contig, please check data.\n", nind);

	MEMORY* oindividual_memory = individual_memory;
	individual_memory = new MEMORY[1];

	IND** newind = new IND*[nfind];
	for (int i = 0, nid = 0; i < nind; ++i)
		if (ainds[i])
		{
			newind[nid] = new(individual_memory->Alloc(sizeof(IND))) IND(*ainds[i]);
			newind[nid]->indid = nid;
			nid++;
		}

	//move genotype
	BUCKET ngeno_bucket;
	ngeno_bucket.FilterIndividualGT(&geno_bucket, nfind, true, g_nthread_val);
	geno_bucket.Replace(ngeno_bucket);

	//allelic depth
	if (ad_bucket.base_addr)
	{
		BUCKET nad_bucket;
		nad_bucket.FilterIndividualAD(&ad_bucket, nfind, true, g_nthread_val);
		ad_bucket.Replace(nad_bucket);
	}

	delete[] ainds; ainds = newind; nind = nfind;
	delete[] oindividual_memory; 
}

/* Remove locus fail to pass filter */
THREAD(RemoveLocus)
{
	MEMORY* nlocus_memory = new MEMORY[g_nthread_val];
	SLOCUS* nslocus = new SLOCUS[nfilter];
	uint64* nlocus_pos = haplotype ? new uint64[nfilter] : NULL;

	if (ad_bucket.base_addr)
	{
		BUCKET nad_bucket;
		nad_bucket.FilterLocusAD(&ad_bucket, true, g_nthread_val);
		ad_bucket.Replace(nad_bucket);
	}

	{
		BUCKET ngeno_bucket;
	}

	{
		BUCKET ngeno_bucket;
		ngeno_bucket.FilterLocusGT(&geno_bucket, true, g_nthread_val, nlocus_memory, nslocus, nlocus_pos);
		geno_bucket.Replace(ngeno_bucket);
	}

	delete[] locus_memory;  locus_memory = nlocus_memory;
	delete[] slocus;        slocus = nslocus;

	if (haplotype)
	{
		delete[] locus_pos; locus_pos = nlocus_pos;
	}

	nloc = nfilter;
	//CheckGenotypeId();
}