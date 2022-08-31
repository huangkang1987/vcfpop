/* Haplotype Extraction Functions */

#pragma once
#include "vcfpop.h"

#ifndef _EXHAPLO
	/* Do nothing */
	TARGET HAPLO_DUMMY_HAPLOTYPE::HAPLO_DUMMY_HAPLOTYPE()
	{

	}

	/* Extract the ith haplotype from an individual */
	TARGET void HAPLO_DUMMY_HAPLOTYPE::ExtractHaplotype(int vi, IND* ti, int64 st, int64 ed, int nvar, ushort aid, MEMORY& haplo_memory)
	{
		alleleid = aid;
		haplo_memory.Alloc(alleles, nvar);
		int sc = 0;

		for (int64 l = st; l <= ed; ++l)
		{
			int ta = ti->GetGenotype(GetLocId(l), GetLoc(l).GetGtab()).GetAlleleCopy(vi);//fine

			//encounter an missing allele, set all haplotype to missing
			if (ta == 0xFFFF)
			{
				SetFF(alleles, nvar);
				return;
			}
			alleles[sc++] = (ushort)ta;
		}
	}

	/* Print information for an extracted locus */
	TARGET void HAPLO_DUMMY_HAPLOTYPE::PrintHaplotype(FILE* f1, int64 st, int64 ed)
{
	fprintf(f1, "%s%d", g_linebreak_val, alleleid);
	for (int64 l = st, sc = 0; l <= ed; ++l)
		fprintf(f1, "%c%s", g_delimiter_val, GetLoc(l).GetAlleleName(alleles[sc++]));
}
#endif


#define extern 
extern atomic<int64> haplotype_contig;				//1st locus in the current contig 
extern LIST<QUICKSORT_PARAMETER> qslstack;			//Sort loci in haplotype extraction
extern LIST<HAPLO_DUMMY_LOCUS> haplotype_locus;	//Dummy locus information in haplotype extraction
extern SLOCUS* haplotype_nslocus;					//Extracted locus
extern uint64* locus_pos;							//Locus pos for small locus used in haplotype extraction
extern LOCN* locus_id;								//Locus id for small locus used in haplotype extraction
#undef extern 

/* Check aneuploid within contigs in haplotype extraction */
THREAD(CheckAneuploid)
{
	AssignPloidyThreadIn();
	nfilterind = 0;

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
	for (int i = 0; i < nind; ++i)
	{
		if (ainds[i]->vmin == ainds[i]->vmax)
		{
			nfilterind++;
			PROGRESS_VALUE++;
			continue;
		}

		byte vst = 0;
		for (int64 st = 0, ed = 0; ed < nloc; ++ed)
		{
			//if ed is in a different contig, reset st
			if (strcmp(GetLoc(st).GetChrom(), GetLoc(ed).GetChrom()))
			{
				st = ed;
				vst = 0;
				continue;
			}

			GENOTYPE& gt = ainds[i]->GetGenotype(GetLocId(ed), GetLoc(ed).GetGtab());//fine
			if (gt.Nalleles() == 0) continue;

			byte ved = (byte)gt.Ploidy();
			if (vst == 0) vst = ved;
			if (vst != ved)
			{
				ainds[i] = NULL;
				break;
			}
		}

		if (ainds[i]) nfilterind++;
		PROGRESS_VALUE++;
	}
}

/* Perform haplotype extraction */
TARGET void CalcHaplotype()
{
	if (abs(g_format_val) > BCF || !haplotype)
	{
		if (locus_pos)
		{
			delete[] locus_pos;
			locus_pos = NULL;
		}
		haplotype = false;
		return;
	}

	if (ad) Exit("\nError: haplotype extraction (-haplotype) is incompatible with allelic depth (-ad) option.\n");
	if (ploidyinfer) Exit("\nError: haplotype extraction (-haplotype) is incompatible with ploidy inference (-ploidyinfer) option.\n");

	EvaluationBegin();

	//Backup original locus id for sorting locus
	if (useslocus)
	{
		locus_id = new LOCN[nloc];
		for (int64 l = 0; l < nloc; ++l)
			GetLocId(l) = l;
	}

	//Sort locus according to contig and pos, so as to iterator locus
	QUICKSORT_PARAMETER par = { 0, nloc - 1 };
	qslstack.Push(par);
	RunThreads(&QSWorker, NULL, NULL, nloc, nloc, "\nSorting loci according to contig/chromosome and position:\n", g_nthread_val, true);

	//Marker aneuploids in the same contig
	RunThreads(&CheckAneuploid, NULL, NULL, nind + nloc, nind + nloc, "\nChecking ansioploids in the same contig:\n", 1, true);

	if (nfilterind == 0)
		Exit("\nError: all %d individuals are excluded from analysis due to the aneuploid in the same contig, please check data.\n", nind);

	//Remove aneuploids in the same contig
	if (nfilterind != nind)
		RunThreads(&RemoveIndividual, NULL, NULL, nloc, nloc, "\nFiltering ansioploids in the same contig to perform haplotype extraction:\n", 1, true);
	CheckGenotypeId();

	//CreateHaplotypeLocus
	haplotype_contig = -1;
	RunThreads(&CreateHaplotypeLocus, NULL, NULL, nloc, nloc, "\nCreating haplotype from VCF file:\n", g_nthread_val, true);
	CheckGenotypeId();
	
	//Sort extracted locus
	int64 nl = haplotype_locus.size;
	if (nl > 1)
	{
		QUICKSORT_PARAMETER par2 = { 0, nl - 1 };
		qslstack.Push(par2);
		RunThreads(&QSHapWorker, NULL, NULL, nl, nl, "\nSorting extracted loci according to contig/chromosome and position:\n", g_nthread_val, true);
	}

	if (nl == 0)
		Exit("\nError: Cannot extract haplotype, because all variants are not genotyped in all individuals or the conditions are too strict.\n");

	//Calculate new genotype index table offset
	haplo_bucket.CreateBucketGT(haplotype_locus);
	CheckGenotypeId();

	//Create new locus
	MEMORY* omemory = locus_memory;
	locus_memory = new MEMORY[g_nthread_val];
	haplotype_nslocus = new SLOCUS[nl];

	//Write haplotypes
	OpenResFile("-haplotype", "Extracted haplotype information");
	OpenTempFiles(g_nthread_val, ".haplotype");
	PROGRESS_VALUE = 0; PROGRESS_TOTAL = nl; PROGRESS_NOUTPUTED2 = PROGRESS_NOUTPUTED; PROGRESS_NOUTPUTED = 0;
	RunThreads(&WriteHaplotypeLocus, NULL, NULL, nl, nl, "\nWriting extracted loci:\n", g_nthread_val, true);
	CheckGenotypeId();

	ALLELE_IDENTIFIER = false;

	//Finish, release memory
	nloc = nl;

	delete[] slocus;
	slocus = haplotype_nslocus;

	delete[] locus_id;
	locus_id = NULL;

	delete[] locus_pos;
	locus_pos = NULL;

	delete[] omemory;

	haplotype_locus.~LIST();//ok
	new (&haplotype_locus) LIST<HAPLO_DUMMY_LOCUS>();

	geno_bucket.Replace(haplo_bucket);

	JoinTempFiles(g_nthread_val);
	CloseResFile();

	//Using unphased genotype and no allele name saved in the locus
	haplotype = false;

	//Calculate total number of alleles and genotypes
	KT = GT = maxK = maxG = 0;
	for (int64 l = 0; l < nloc; ++l)
	{
		KT += GetLoc(l).k;
		GT += GetLoc(l).ngeno;
		maxK = Max((int)GetLoc(l).k, maxK);
		maxG = Max((int)GetLoc(l).ngeno, maxG);
	}
	EvaluationEnd("Haplotype extraction");
}

/* Quick sort locus by contig and position */
TARGET void QSLocus(int64 left, int64 right)
{
	int64 i = left, j = right;

	if (right - left < 10)
	{
		for (int64 ii = left; ii <= right; ++ii)
			for (int64 jj = ii + 1; jj <= right; ++jj)
				if (strcmp(slocus[ii].GetChrom(), slocus[jj].GetChrom()) > 0 ||
					(strcmp(slocus[ii].GetChrom(), slocus[jj].GetChrom()) == 0 && GetLocPos(ii) > GetLocPos(jj)))
				{
					Swap(GetLocPos(ii), GetLocPos(jj));
					Swap(GetLocId(ii), GetLocId(jj));
					Swap(GetLoc(ii), GetLoc(jj));
				}

		PROGRESS_VALUE += right - left + 1;
		return;
	}

	int64 mid = (left + right) >> 1;
	SLOCUS& spivot = slocus[mid];

	while (left < j || i < right)
	{
		while (strcmp(slocus[i].GetChrom(), spivot.GetChrom()) < 0 ||
			(strcmp(slocus[i].GetChrom(), spivot.GetChrom()) == 0 && GetLocPos(i) < GetLocPos(mid))) i++;
		while (strcmp(slocus[j].GetChrom(), spivot.GetChrom()) > 0 ||
			(strcmp(slocus[j].GetChrom(), spivot.GetChrom()) == 0 && GetLocPos(j) > GetLocPos(mid))) j--;

		if (i <= j)
		{
			Swap(GetLocPos(i), GetLocPos(j));
			Swap(GetLocId(i), GetLocId(j));
			Swap(GetLoc(i), GetLoc(j));
			i++; j--;
		}

		if (i > j)
		{
			if (i == j + 2)
				PROGRESS_VALUE++;

			if (left < j)
			{
				QUICKSORT_PARAMETER par = { left, j };
				Lock(GLOCK2);
				qslstack.Push(par);
				UnLock(GLOCK2);
			}
			else if (left == j)
				PROGRESS_VALUE++;

			if (i < right)
			{
				QUICKSORT_PARAMETER par = { i, right };
				Lock(GLOCK2);
				qslstack.Push(par);
				UnLock(GLOCK2);
			}
			else if (i == right)
				PROGRESS_VALUE++;

			return;
		}
	}
}

/* Quick sort locus in a contig */
THREAD(QSWorker)
{
	QUICKSORT_PARAMETER par;
	while (PROGRESS_VALUE != PROGRESS_CEND)
	{
		bool hastask = false;
		Lock(GLOCK2);
		if (qslstack.size)
		{
			hastask = true;
			par = qslstack.Pop();
		}
		UnLock(GLOCK2);
		if (hastask)
			QSLocus(par.left, par.right);
	}
}

/* Quick sort extracted locus by contig and position */
TARGET void QSHapLocus(int64 left, int64 right)
{
	int64 i = left, j = right;

	int64 mid = (left + right) >> 1;
	HAPLO_DUMMY_LOCUS& pivot = haplotype_locus[mid];

	while (left < j || i < right)
	{
		while (haplotype_locus[i].chrid < pivot.chrid || (haplotype_locus[i].chrid == pivot.chrid && haplotype_locus[i].stpos < pivot.stpos)) i++;
		while (haplotype_locus[j].chrid > pivot.chrid || (haplotype_locus[j].chrid == pivot.chrid && haplotype_locus[j].stpos > pivot.stpos)) j--;

		if (i <= j)
		{
			Swap(haplotype_locus[i], haplotype_locus[j]);
			i++; j--;
		}

		if (i > j)
		{
			if (i == j + 2) PROGRESS_VALUE++;
			if (left < j)
			{
				QUICKSORT_PARAMETER par = { left, j };
				Lock(GLOCK2);
				qslstack.Push(par);
				UnLock(GLOCK2);
			}
			else if (left == j) PROGRESS_VALUE++;
			if (i < right)
			{
				QUICKSORT_PARAMETER par = { i, right };
				Lock(GLOCK2);
				qslstack.Push(par);
				UnLock(GLOCK2);
			}
			else if (i == right) PROGRESS_VALUE++;
			return;
		}
	}
}

/* Quick sort extracted locus in a contig */
THREAD(QSHapWorker)
{
	QUICKSORT_PARAMETER par;
	while (PROGRESS_VALUE != PROGRESS_CEND)
	{
		bool hastask = false;
		Lock(GLOCK2);
		if (qslstack.size)
		{
			hastask = true;
			par = qslstack.Pop();
		}
		UnLock(GLOCK2);

		if (hastask)
			QSHapLocus(par.left, par.right);
	}
}

/* Get number of alleles and genotypes at a dummy locus */
TARGET double GetDummyK(int64 st, int64 ed, TABLE<HASH, ushort>& hfidx, TABLE<HASH, ushort>& gfidx)
{
	hfidx.Clear(); gfidx.Clear();
	HASH hash[N_MAX_PLOIDY];
	ushort alleles[N_MAX_PLOIDY];
	bool usedploidy[N_MAX_PLOIDY + 1] = { 0 };
	int ntyped = 0;

	for (int i = 0; i < nind; ++i)
	{
		int v = 0;
		HashHaplotype(ainds[i], st, ed, hash, v);

		if (!usedploidy[v])
		{
			usedploidy[v] = true;
			gfidx.PushIndex(missing_hash[v]);
		}

		if (hash[0] != (HASH)-1)
		{
			for (int vi = 0; vi < v; ++vi)
				alleles[vi] = (ushort)hfidx.PushIndex(hash[vi]);
			Sort(alleles, v); //unphase
			gfidx.PushIndex(HashGenotype(alleles, v));
			ntyped++;
		}
	}

	return ntyped / (double)nind;
}

/* Create locus for haplotype extraction */
THREAD(CreateHaplotypeLocus)
{
	TABLE<HASH, ushort> hfidx(false, NULL), gfidx(false, NULL);
	LIST<HAPLO_DUMMY_LOCUS> dlocus(NULL);    //dummy locus

	//For each contig, find is begin and end variant
	for (int64 contig_st = -1, chrid = 0, contig_ed = -1; contig_st < nloc; )
	{
		int64 test_val = contig_st, next_val = contig_ed + 1;
		while (!haplotype_contig.compare_exchange_strong(test_val, next_val))
		{
			//jump a contig
			contig_st = contig_ed = contig_ed + 1;
			if (contig_st >= nloc) return;
			while (contig_ed < nloc && strcmp(GetLoc(contig_st).GetChrom(), GetLoc(contig_ed).GetChrom()) == 0) contig_ed++; contig_ed--;
			chrid++;
			test_val = contig_st; next_val = contig_ed + 1;
		}

		//1. Set i = 1 and j = 1 in the beginning.
		contig_st = contig_ed = contig_ed + 1;
		while (contig_ed < nloc && strcmp(GetLoc(contig_st).GetChrom(), GetLoc(contig_ed).GetChrom()) == 0) contig_ed++; contig_ed--;
		chrid++;
		int64 st = contig_st, ed = contig_st;

		//2. If the value of j exceeds the number of variants in this chromosome or contig, then terminate this algorithm. 
	step2:
		if (st > ed)
		{
			PROGRESS_VALUE += st - ed;
			ed = st;
		}
		if (ed > contig_ed) continue;

		//3. Calculate the values of several parameters (the length of each haplotype, the number of variants, the number of alleles and the number of genotypes). 
		int64 clen = GetLocPos(ed) - GetLocPos(st) + 1;
		int64 nvariants = ed - st + 1;
		double gtrate = GetDummyK(st, ed, hfidx, gfidx);

		//4. If any parameter exceeds the upper bound then increase i and go to step 2. 
		if (clen > haplotype_length_max ||
			nvariants > haplotype_variants_max ||
			gtrate < haplotype_ptype_min ||
			hfidx.size > haplotype_alleles_max ||
			gfidx.size > haplotype_genotypes_max)
		{
			st++;
			goto step2;
		}

		//5. If any parameter exceeds the lower bound then increase j and go to step 2. 
		if (clen < haplotype_length_min ||
			nvariants < haplotype_variants_min ||
			gtrate > haplotype_ptype_max ||
			hfidx.size < haplotype_alleles_min ||
			gfidx.size < haplotype_genotypes_min)
		{
			ed++;
			PROGRESS_VALUE++;
			goto step2;
		}

		//6. Combine all variants between the variants i and j, and extract the haplotypes.
		HAPLO_DUMMY_LOCUS tentry{ chrid, (int64)GetLocPos(st), st, ed, (int)nvariants, (int)hfidx.size, (int)gfidx.size };
		Lock(GLOCK2);
		haplotype_locus.Push(tentry);
		UnLock(GLOCK2);

		//7. Set i and j to the next applicable variant according to -haplotype_interval, and do to Step 2.
		int64 nextpos = (int64)GetLocPos(ed) + haplotype_interval_val;
		while ((int64)GetLocPos(st) <= nextpos) st++;
		goto step2;
	}
}

/* Output locus for haplotype extraction */
THREAD(WriteHaplotypeLocus)
{
	char name_buf[NAME_BUF_LEN];
	int nthread = g_nthread_val;
	int64 newloc = haplotype_locus.size;

	double nsec = newloc / (double)nthread + 1e-8;
	int64 st1 = (int64)(threadid * nsec), ed1 = (int64)((threadid + 1) * nsec);
	FILE* fout = TEMP_FILES[threadid];
	MEMORY haplo_memory;

	TABLE<HASH, HAPLO_DUMMY_HAPLOTYPE*> hftab(true, NULL);
	TABLE<HASH, TEMP_GENOTYPE> temptab(true, NULL);

	for (int64 l = st1; l < ed1; ++l)
	{
		hftab.Clear();    temptab.Clear();

		HAPLO_DUMMY_LOCUS& dlocus = haplotype_locus[l];
		HAPLO_DUMMY_LOCUS* dnext = l + 1 < newloc ? &haplotype_locus[l + 1] : NULL;

		int64 st = dlocus.st, ed = dlocus.ed, nvar = dlocus.nvar;
		int gasize = 0;
		GENO_WRITER wt(l, &haplo_bucket);

		for (int i = 0; i < nind; ++i)
		{
			HASH haplohash[N_MAX_PLOIDY];
			ushort alleles[N_MAX_PLOIDY] = { 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF };
			int ploidy = 0;

			//get haplotype hash
			HashHaplotype(ainds[i], st, ed, haplohash, ploidy);

			//missing data
			if (haplohash[0] == (HASH)-1)
				SetFF(alleles, ploidy);

			//map hash into allele index
			else for (int ai = 0; ai < ploidy; ++ai)
			{
				HASH ha = haplohash[ai];
				HAPLO_DUMMY_HAPLOTYPE* ht = NULL;

				if ((ht = hftab.Try(ha)) == NULL)
				{
					haplo_memory.Alloc(ht, 1);
					hftab[ha] = ht;
					ht->ExtractHaplotype(ai, ainds[i], st, ed, nvar, (ushort)hftab.size - 1, haplo_memory);
				}
				alleles[ai] = ht->alleleid;
			}

			//unphase
			Sort(alleles, ploidy);

			//calculate genotype hash
			HASH hash = HashGenotype(alleles, ploidy);

			//add missing genotype
			if (!temptab.ContainsKey(missing_hash[ploidy]))
			{
				TEMP_GENOTYPE& tgt = temptab[missing_hash[ploidy]];
				SetVal(tgt.alleles, missing_array, N_MAX_PLOIDY);
				tgt.hash = missing_hash[ploidy];
				tgt.ploidy = ploidy;
				tgt.gid = temptab.size - 1;
				gasize += 0;
			}

			//add genotype
			if (!temptab.ContainsKey(hash))
			{
				TEMP_GENOTYPE& tgt = temptab[hash];
				SetVal(tgt.alleles, alleles, N_MAX_PLOIDY);
				tgt.hash = hash;
				tgt.ploidy = ploidy;
				tgt.gid = temptab.size - 1;
				gasize += ploidy + GetNalleles(alleles, ploidy);
			}

			//write genotype id
			wt.Write(temptab[hash].gid);
		}
		wt.FinishWrite();

		//create locus
		int ngeno = temptab.size;

		if (ngeno != haplotype_locus[l].gsize)
			Exit("\nError: haplotype extraction find more genotypes than expect.");

		SLOCUS* sloc = new(&haplotype_nslocus[l]) SLOCUS(locus_memory[threadid], slocus[st], l, ngeno, gasize, temptab);

		//write results
		int hsize = hftab.size;
		sloc->k = (ushort)hsize;
		fprintf(fout, "%s%sLocus:%s", g_linebreak_val, g_linebreak_val, sloc->GetNameStr(name_buf));
		fprintf(fout, "%s#CHROM:%s", g_linebreak_val, sloc->GetChrom());
		fprintf(fout, "%s#Variants:%lld", g_linebreak_val, nvar);
		fprintf(fout, "%s#Haplotypes:%d", g_linebreak_val, hftab.size);
		fprintf(fout, "%sRange:%lld-%lld", g_linebreak_val, GetLocPos(dlocus.st), GetLocPos(dlocus.ed));
		fprintf(fout, "%sLength:%lld", g_linebreak_val, GetLocPos(dlocus.ed) - GetLocPos(dlocus.st) + 1);

		if (dnext && strcmp(GetLoc(dnext->st).GetChrom(), sloc->GetChrom()) == 0)
			fprintf(fout, "%sDistance to next extracted locus:%lld", g_linebreak_val, GetLocPos(dnext->st) - GetLocPos(dlocus.ed) - 1);
		else
			fprintf(fout, "%sNext extracted locus is in the different contig or chromosome", g_linebreak_val);
		fprintf(fout, "%sHapId", g_linebreak_val);

		for (int64 l2 = st; l2 <= ed; ++l2)
			fprintf(fout, "%c%s", g_delimiter_val, GetLoc(l2).GetNameStr(name_buf));

		for (int hi = 0; hi < hsize; ++hi)
			hftab(hi)->PrintHaplotype(fout, st, ed);
		haplo_memory.ClearMemory();

		PROGRESS_VALUE++;
	}
}
