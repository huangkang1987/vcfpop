/* Filter Functions */

#include "vcfpop.h"

template TARGET void ApplyFilter<double>();
template TARGET void ApplyFilter<float >();
template TARGET void SetVReg<double>(int rl, int i);
template TARGET void SetVReg<float >(int rl, int i);
template TARGET void AssignVInds<double>();
template TARGET void AssignVInds<float >();

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
template<typename REAL>
TARGET void ApplyFilter()
{
	EvaluationBegin();
	CheckGenotypeId<REAL>();

	//1. Applying individual filter
	if (individual_filter)
	{
		RunThreads(&MarkerIndividual<REAL>, NULL, NULL, nind + nloc * 2, nind + nloc,
			"\nApplying individual filter:\n", 1, true, (int)((nind + nloc) / (nind + nloc * 2) * g_progress_val));

		RunThreads(&RemoveIndividual<REAL>, NULL, NULL, nind + nloc * 2, nloc,
			NULL, 1, false, g_progress_val - (int)((nind + nloc) / (nind + nloc * 2) * g_progress_val));
	}
	CheckGenotypeId<REAL>();

	//Assign individuals into their vpop
	AssignVInds<REAL>();
	CheckGenotypeId<REAL>();

	//2. Apply diversity filter
	if (diversity_filter)
	{
		AssignPop<REAL>(f_pop_b, f_pop_val, "-f_pop");

		//Temp use allele frequency and genotype count for diversity calculation
		allele_freq_offset = new uint64[nloc];
		genotype_count_offset = new uint64[nloc];
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
			maxK = std::max((int)GetLoc(l).k, maxK);
			maxG = std::max((int)GetLoc(l).ngeno, maxG);
		}

		//Calculate freq for cpop and diversity estimation
		cpop<REAL>->AllocFreq();

		if (f_windowsize_b && abs(g_format_val) <= BCF)
			new (&contig_diversity) TABLE<HASH, TABLE<HASH, pair<LOCN, double>>>(true, NULL);

		RunThreads(&MarkerDiversity<REAL>, NULL, NULL, nloc, nloc, "\nMarking locus diversity filter:\n", 1, true);

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
				contig.bucket = NULL;
			}
			contig_diversity.~TABLE();//ok
			contig_diversity.bucket = NULL;
			new (&contig_diversity) TABLE<HASH, TABLE<HASH, pair<LOCN, double>>>;
		}
		cpop<REAL>->UnAllocFreq();

		//Release temp allele frequency and genotype count
		DEL(allele_freq_offset);
		DEL(genotype_count_offset);
	}

	CheckGenotypeId<REAL>();

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
	DEL(allele_table);
	CheckGenotypeId<REAL>();

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
			maxK = std::max((int)GetLoc(l).k, maxK);
			maxG = std::max((int)GetLoc(l).ngeno, maxG);
		}
	}

	CheckGenotypeId<REAL>();

	if (nfilter == 0)
		Exit("\nError: all %d loci are excluded from analysis, try less strict filters.\n", nloc);

	if (nfilter != nloc)
		RunThreads(&RemoveLocus<REAL>, NULL, NULL, nloc * (ad_bucket.base_addr ? 2 : 1), nloc * (ad_bucket.base_addr ? 2 : 1), "\nApplying locus diversity filter:\n", 1, true);

	CheckGenotypeId<REAL>();

	EvaluationEnd("Filter genotypes, individuals and loci");
}

/* Recursive set vid, vpop, ind for each region */
template<typename REAL>
TARGET void SetVReg(int rl, int i)
{
	if (rl >= 0)
	{
		if (reg<REAL>[rl][i].npop == 0 || reg<REAL>[rl][i].nind == 0) return;
		reg<REAL>[rl][i].id = nreg[rl];
		reg<REAL>[rl][i].vpop = rl >= 1 ? aregs<REAL>[rl - 1] + nreg[rl - 1] : apops<REAL> + npop;
		reg<REAL>[rl][i].inds = rinds<REAL> + nind;
		reg<REAL>[rl][i].ind0id = nind;
		aregs<REAL>[rl][nreg[rl]++] = &reg<REAL>[rl][i];
	}

	if (rl > 0) for (uint j = 0; j < reg<REAL>[rl - 1].size; ++j)
	{
		if (reg<REAL>[rl - 1][j].rid != i) continue;
		reg<REAL>[rl - 1][j].rid = reg<REAL>[rl][i].id;
		SetVReg<REAL>(rl - 1, j);
	}
	else for (uint j = 0; j < pop<REAL>.size; ++j)
	{
		if (pop<REAL>[j].rid != i || pop<REAL>[j].nind == 0) continue;

		apops<REAL>[npop] = &pop<REAL>[j];
		pop<REAL>[j].id = npop++;
		pop<REAL>[j].rid = reg<REAL>[rl][i].id;
		pop<REAL>[j].vpop = NULL;
		pop<REAL>[j].inds = rinds<REAL> + nind;
		pop<REAL>[j].ind0id = nind;
		nind += pop<REAL>[j].nind;
	}
}

/* Sort individuals by population index to rinds<REAL> array */
template<typename REAL>
TARGET void AssignVInds()
{
	/* Original order is in ind and should not be changed to correctly obtain genotype index */

	//Assign pop/reg
	for (uint i = 0; i < pop<REAL>.size; ++i)
		pop<REAL>[i].nind = 0;

	for (uint rl = 0; rl < reg<REAL>.size; ++rl)
		for (uint i = 0; i < reg<REAL>[rl].size; ++i)
			reg<REAL>[rl][i].npop = reg<REAL>[rl][i].nind = 0;

	//First scan, count inds
	for (int j = 0; j < nind; ++j)
	{
		//default pop or input pop
		int i = (ainds<REAL>[j]->popid < 1 || ainds<REAL>[j]->popid >= pop<REAL>.size) ? 0 : ainds<REAL>[j]->popid;
		pop<REAL>[i].nind++;
	}

	//Allocate vreg
	nregt = nregt2 = 0;
	lreg = reg<REAL>.size;
	aregs<REAL> = new POP<REAL>**[reg<REAL>.size];
	SetZero(aregs<REAL>, reg<REAL>.size);
	SetZero(nreg, lreg);

	//Add up nind and npop to regions at lay 0
	for (uint i = 0; i < pop<REAL>.size; ++i)
	{
		if (pop<REAL>[i].nind == 0) continue;
		if (reg<REAL>.size > 0 && reg<REAL>[0].size > 0)
		{
			reg<REAL>[0][pop<REAL>[i].rid].npop++;
			reg<REAL>[0][pop<REAL>[i].rid].nind += pop<REAL>[i].nind;
		}
	}

	//Add up nind and npop to regions at lay 1+
	for (uint rl = 0; rl < reg<REAL>.size; ++rl)
	{
		for (uint i = 0; i < reg<REAL>[rl].size; ++i)
		{
			if (reg<REAL>[rl][i].npop > 0)
			{
				nreg[rl]++;  nregt++;
				if (rl < reg<REAL>.size - 1)
				{
					reg<REAL>[rl + 1][reg<REAL>[rl][i].rid].npop++;
					reg<REAL>[rl + 1][reg<REAL>[rl][i].rid].nind += reg<REAL>[rl][i].nind;
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
	for (uint i = 0; i < pop<REAL>.size; ++i)
		if (pop<REAL>[i].nind > 0)
			npop++;
	apops<REAL> = new POP<REAL>* [npop];

	//Contruct total pop
	if (lreg == 0 && npop == 1)
	{
		lreg = -1;
		for (uint i = 0; i < pop<REAL>.size; ++i)
			if (pop<REAL>[i].nind > 0)
				total_pop<REAL> = &pop<REAL>[i];
		nregt2 = nregt = 0;
	}
	else
	{
		for (uint i = 0; i < reg<REAL>[lreg].size; ++i)
			if (reg<REAL>[lreg][i].nind > 0)
				total_pop<REAL> = &reg<REAL>[lreg][i];
	}

	DEL(total_pop<REAL>->name);
	total_pop<REAL>->name = new char[6];
	sprintf(total_pop<REAL>->name, "Total");
	total_pop<REAL>->rid = -1;

	for (int rl = 0; rl <= lreg; ++rl)
	{
		aregs<REAL>[rl] = new POP<REAL>* [nreg[rl]];
		nreg[rl] = 0;
	}

	//Allocate rinds<REAL>
	rinds<REAL> = new IND<REAL>* [nind];

	//Recursive arrange vpop, vreg and vind
	if (lreg == -1)
	{
		apops<REAL>[0] = total_pop<REAL>;
		total_pop<REAL>->id = 0;
		total_pop<REAL>->rid = -1;
		total_pop<REAL>->vpop = NULL;
		total_pop<REAL>->inds = rinds<REAL>;
		total_pop<REAL>->ind0id = 0;
	}
	else
	{
		nind = 0; npop = 0;
		SetVReg<REAL>(lreg, total_pop<REAL>->id);
	}

	//Set nind = 0
	for (uint i = 0; i < pop<REAL>.size; ++i)
		pop<REAL>[i].nind = 0;

	//Move ind to pop ind array and assign vpopid
	bool swap = false;
	for (int j = 0; j < nind; ++j)
	{
		//Default pop or input pop
		IND<REAL>& i = *ainds<REAL>[j];
		POP<REAL>& p = pop<REAL>[(i.popid < 1 || i.popid >= pop<REAL>.size) ? 0 : i.popid];
		if (i.indid != p.ind0id + p.nind)
			swap = true;
		i.indid = p.ind0id + p.nind;
		p.inds[p.nind++] = ainds<REAL>[j];
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
			gtb[ainds<REAL>[i]->indid] = (ushort)rg.Read();

		for (int i = 0; i < nind; ++i)
			wg.Write(gtb[i]);

		wg.FinishWrite();
	}

	DEL(gtbuf);

	if (ploidyinfer)
	{
		maxK = 0;
		for (int64 l = 0; l < nloc; ++l)
			maxK = std::max((int)GetLoc(l).k, maxK);

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
				uint* dpb2 = dpb + ainds<REAL>[i]->indid * k2;

				for (int k = 0; k < k2; ++k)
					dpb2[k] = rd.Read();
			}

			for (int i = 0, nend = k2 * nind; i < nend; ++i)
				wd.Write(*dpb++);
		}

		DEL(dpbuf);
	}

	DEL(ainds<REAL>); 
	ainds<REAL> = rinds<REAL>;
}

/* Marker individual filtered or not */
THREAD2(MarkerIndividual)
{
	//Calculate offset
	maxG = 0;
	for (int64 l = 0; l < nloc; ++l)
		maxG = std::max((int)GetLoc(l).ngeno, maxG);

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
	CheckGenotypeId<REAL>();

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
		if (ainds<REAL>[i] && f_itype_b && (f_ntype_min > ntype[i] || ntype[i] > f_ntype_max))
			ainds<REAL>[i] = NULL;

		if (ainds<REAL>[i] && f_iploidy_b && (f_iploidy_min > minploidy || minploidy > f_iploidy_max))
			ainds<REAL>[i] = NULL;

		if (ainds<REAL>[i] && f_iploidy_b && (f_iploidy_min > maxploidy || maxploidy > f_iploidy_max))
			ainds<REAL>[i] = NULL;

		PROGRESS_VALUE++;
	}

	DEL(ntype);
	DEL(vt);
	DEL(vmin);
	DEL(vmax);

	VLA_DELETE(ploidytab);
	VLA_DELETE(nallelestab);
}

/* Marker locus filtered or not */
THREAD2(MarkerDiversity)
{
	bool flag_window = f_windowsize_b && abs(g_format_val) <= BCF;
#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
	for (int64 l = 0; l < nloc; ++l)
	{
		DIVERSITY<REAL> d;
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

				pair<LOCN, double>& ww = contig[w];
				if (ww.second < div)
				{
					ww.first = l;
					ww.second = div;
				}
			}

		}
		else
			GetLoc(l).flag_pass = false;

		PROGRESS_VALUE++;
	}
}

/* Remove individual fail to pass filter */
THREAD2(RemoveIndividual)
{
	int nfind = 0;
	for (int i = 0; i < nind; ++i)
		if (ainds<REAL>[i]) nfind++;

	if (nfind == nind)
	{
		PROGRESS_VALUE = PROGRESS_CEND;
		return;
	}

	if (nfind == 0)
		Exit("\nError: all %d individuals are excluded from analysis due to the aneuploid in the same contig, please check data.\n", nind);

	MEMORY* oindividual_memory = individual_memory;
	individual_memory = new MEMORY[1];

	IND<REAL>** newind = new IND<REAL>*[nfind];
	for (int i = 0, nid = 0; i < nind; ++i)
		if (ainds<REAL>[i])
		{
			newind[nid] = new(individual_memory->Alloc(sizeof(IND<REAL>))) IND<REAL>(*ainds<REAL>[i]);
			newind[nid]->indid = nid;
			nid++;
		}

	//move genotype
    {
        BUCKET ngeno_bucket;
        ngeno_bucket.FilterIndividualGT<REAL>(&geno_bucket, nfind, true, g_nthread_val);
        geno_bucket.Replace(ngeno_bucket);
    }
    
	//allelic depth
	if (ad_bucket.base_addr)
	{
		BUCKET nad_bucket;
		nad_bucket.FilterIndividualAD<REAL>(&ad_bucket, nfind, true, g_nthread_val);
		ad_bucket.Replace(nad_bucket);
	}

	DEL(ainds<REAL>); ainds<REAL> = newind; nind = nfind;
	DEL(oindividual_memory);
}

/* Remove locus fail to pass filter */
THREAD2(RemoveLocus)
{
	MEMORY* nlocus_memory = new MEMORY[g_nthread_val];
	SLOCUS* nslocus = new SLOCUS[nfilter];
	uint64* nlocus_pos = uselocpos ? new uint64[nfilter] : NULL;

	if (ad_bucket.base_addr)
	{
		BUCKET nad_bucket;
		nad_bucket.FilterLocusAD(&ad_bucket, true, g_nthread_val);
		ad_bucket.Replace(nad_bucket);
	}

	{
		BUCKET ngeno_bucket;
		ngeno_bucket.FilterLocusGT(&geno_bucket, true, g_nthread_val, nlocus_memory, nslocus, nlocus_pos);
		geno_bucket.Replace(ngeno_bucket);
	}

	DEL(locus_memory);  locus_memory = nlocus_memory;
	DEL(slocus);        slocus = nslocus;

	if (uselocpos)
	{
		DEL(locus_pos);
		locus_pos = nlocus_pos;
	}

	nloc = nfilter;
	//CheckGenotypeId<REAL>();
}
