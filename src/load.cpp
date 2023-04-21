/* Load File Functions */

#pragma once
#include "vcfpop.h"

template TARGET void CheckGenotypeId<double>();
template TARGET void CheckGenotypeId<float >();
template TARGET void AssignIndSub<double>(int indid, int popid, int& namelen, int& nind2);
template TARGET void AssignIndSub<float >(int indid, int popid, int& namelen, int& nind2);
template TARGET void AssignInds<double>();
template TARGET void AssignInds<float >();
template TARGET void LoadFile<double>();
template TARGET void LoadFile<float >();
template struct IND<double>;
template struct IND<float >;

template TARGET void IND<double>::genepop(char* t, bool iscount, GENOTYPE** gtab, ushort** gatab, GENO_WRITER* wt);
template TARGET void IND<float >::genepop(char* t, bool iscount, GENOTYPE** gtab, ushort** gatab, GENO_WRITER* wt);
template TARGET void IND<double>::spagedi(char* t, bool iscount, GENOTYPE** gtab, ushort** gatab, GENO_WRITER* wt);
template TARGET void IND<float >::spagedi(char* t, bool iscount, GENOTYPE** gtab, ushort** gatab, GENO_WRITER* wt);
template TARGET void IND<double>::cervus(char* t, bool iscount, GENOTYPE** gtab, ushort** gatab, GENO_WRITER* wt);
template TARGET void IND<float >::cervus(char* t, bool iscount, GENOTYPE** gtab, ushort** gatab, GENO_WRITER* wt);
template TARGET void IND<double>::arlequin(char* t, bool iscount, GENOTYPE** gtab, ushort** gatab, GENO_WRITER* wt);
template TARGET void IND<float >::arlequin(char* t, bool iscount, GENOTYPE** gtab, ushort** gatab, GENO_WRITER* wt);
template TARGET void IND<double>::structure(char* t, bool iscount, GENOTYPE** gtab, ushort** gatab, GENO_WRITER* wt);
template TARGET void IND<float >::structure(char* t, bool iscount, GENOTYPE** gtab, ushort** gatab, GENO_WRITER* wt);
template TARGET void IND<double>::polygene(char* t, bool iscount, GENOTYPE** gtab, ushort** gatab, GENO_WRITER* wt);
template TARGET void IND<float >::polygene(char* t, bool iscount, GENOTYPE** gtab, ushort** gatab, GENO_WRITER* wt);
template TARGET void IND<double>::polyrelatedness(char* t, bool iscount, GENOTYPE** gtab, ushort** gatab, GENO_WRITER* wt);
template TARGET void IND<float >::polyrelatedness(char* t, bool iscount, GENOTYPE** gtab, ushort** gatab, GENO_WRITER* wt);
template TARGET void IND<double>::genodive(char* t, bool iscount, GENOTYPE** gtab, ushort** gatab, GENO_WRITER* wt);
template TARGET void IND<float >::genodive(char* t, bool iscount, GENOTYPE** gtab, ushort** gatab, GENO_WRITER* wt);
template TARGET void IND<double>::plink(char* t, bool iscount, GENOTYPE** gtab, ushort** gatab, GENO_WRITER* wt);
template TARGET void IND<float >::plink(char* t, bool iscount, GENOTYPE** gtab, ushort** gatab, GENO_WRITER* wt);
template TARGET void IND<double>::AddBCFGenotype(int64 l, char*& gtstr, char*& gqstr, char*& dpstr, char*& adstr, int vlen, int asize, int gqlen, int dplen, int adlen, uint*& depth, TABLE<HASH, uint>& gfid, GENOTYPE*& gtab, ushort*& gatab, GENO_WRITER& wt);
template TARGET void IND<float >::AddBCFGenotype(int64 l, char*& gtstr, char*& gqstr, char*& dpstr, char*& adstr, int vlen, int asize, int gqlen, int dplen, int adlen, uint*& depth, TABLE<HASH, uint>& gfid, GENOTYPE*& gtab, ushort*& gatab, GENO_WRITER& wt);
template TARGET void IND<double>::AddVCFGenotype(char*& line, int64 l, uint*& depth, TABLE<HASH, uint>& gfid, GENOTYPE*& gtab, ushort*& gatab, GENO_WRITER& wt);
template TARGET void IND<float >::AddVCFGenotype(char*& line, int64 l, uint*& depth, TABLE<HASH, uint>& gfid, GENOTYPE*& gtab, ushort*& gatab, GENO_WRITER& wt);
template TARGET char* IND<double>::GetTagValue(char* re, int tagid);
template TARGET char* IND<float >::GetTagValue(char* re, int tagid);

#define extern 
extern int genotype_digit;
extern int genotype_extracol;
extern int genotype_missing;
extern int genotype_ambiguous;
extern bool load_complete;
extern bool usephase;
extern bool uselocpos;

extern INCBUFFER* load_buf;							//Circle buffer for loading vcf/bcf files, NBUF
extern char* vcf_header;							//Vcf header row
#undef extern 

/* Check Genotype index is correct */
template<typename REAL>
TARGET void CheckGenotypeId()
{
	return;
	for (int64 l = 0; l < nloc; ++l)
	{
		int ngeno = (useslocus ? slocus[l] : locus[l]).ngeno;
		GENO_READER rt(0, l);
		for (int i = 0; i < nind; ++i)
		{
			int gid1 = rt.Read();
			int gid2 = ainds[i]->GetGenotypeId(l);
			if (gid1 != gid2 || gid1 >= ngeno || gid1 < 0)
				Exit("\nError, genotype id mismatch.");
		}
	}
}

/* Assign individual indid to population popid */
template<typename REAL>
TARGET void AssignIndSub(int indid, int popid, int& namelen, int& nind2)
{
	if (ainds[indid]->popid != 0)
		Exit("\nError: populations %s and %s have the same individual %s.", pop<REAL>[ainds[indid]->popid].name, pop<REAL>[popid].name, ainds[indid]->name);

	ainds[indid]->popid = (ushort)popid;
	namelen += (int)strlen(ainds[indid]->name) + 1;
	nind2++;
}

/* Assign individuals to their populations */
template<typename REAL>
TARGET void AssignInds()
{
	//Assign all individuals to defpop initially and prepare name to indid hash table
	TABLE<HASH, int> name2id(false, NULL);
	for (int i = 0; i < nind; ++i)
	{
		ainds[i]->popid = 0;
		HASH ha = HashString(ainds[i]->name);
		if (name2id.ContainsKey(ha))
			Exit("\nError: individiuals #%d and #%d have the same name %s.", name2id[ha], i, ainds[i]->name);
		name2id[ha] = i;
	}

	//Parse indtext
	for (int pp = pop<REAL>.size - 1; pp >= 0; --pp)
	{
		int nind2 = 0;
		int namelen = 0;
		for (int j2 = 0; j2 < pop<REAL>[pp].nind; ++j2)
		{
			int v1 = 0, v2 = 0;
			ParseTwoNumber(pop<REAL>[pp].names[j2], v1, v2);

			//id format
			if (v1 == -1 && v2 == -1)
			{
				HASH ha = HashString(pop<REAL>[pp].names[j2]);
				if (!name2id.ContainsKey(ha))
					Exit("\nError: individiual %s in population %s is not in the genotype data.", pop<REAL>[pp].names[j2], pop<REAL>[pp].name);

				AssignIndSub<REAL>(name2id[ha], pp, namelen, nind2);
			}

			//#1 or id format
			else if (v2 == -1)
			{
				HASH ha = HashString(pop<REAL>[pp].names[j2]);
				int indid = -1;
				if (name2id.ContainsKey(ha))
					indid = name2id[ha];
				else
					indid = v1 - 1;

				if (indid >= 0 && indid < nind)
					AssignIndSub<REAL>(indid, pp, namelen, nind2);
				else
					Exit("\nError: can not found individual %s for population %s.", pop<REAL>[pp].names[j2], pop<REAL>[pp].name);
			}

			//#1-#2 format
			else
			{
				if (v1 <= 0 || v1 > nind)
					Exit("\nError: the first number of indtext %s for population %s exceeds range.", pop<REAL>[pp].names[j2], pop<REAL>[pp].name);
				if (v2 <= 0 || v2 > nind)
					Exit("\nError: the second number of indtext %s for population %s exceeds range.", pop<REAL>[pp].names[j2], pop<REAL>[pp].name);

				int vmin = Min(v1, v2), vmax = Max(v1, v2);
				for (int i = vmin - 1; i <= vmax - 1; ++i)
					AssignIndSub<REAL>(i, pp, namelen, nind2);
			}
		}

		if (pop<REAL>[pp].nind == 0) for (int i = 0; i < nind; ++i)
			if (ainds[i]->popid == pp)
			{
				ainds[i]->popid = (ushort)pp;
				namelen += (int)strlen(ainds[i]->name) + 1;
				nind2++;
			}

		//Allocate individuals names
		delete[] pop<REAL>[pp].names;
		pop<REAL>[pp].names = NULL;
		pop<REAL>[pp].nind = 0;
	}

	//pop<REAL>[i].inds will be allocated and copy individuals after filtered in function AssignVInds
}

/* Load data from input files */
template<typename REAL>
TARGET void LoadFile()
{
	EvaluationBegin();

	usephase = haplotype || (slide && (slide_estimator_val[12] || slide_estimator_val[13]));
	uselocpos = haplotype || slide || (diversity_filter && f_windowsize_b && abs(g_format_val) <= BCF) || (convert && convert_format_val[9]);

	//vcf/bcf format
	if (abs(g_format_val) <= BCF)
	{
		int64 buflen = 1024 * 64;

		INCBUFFER readbuf;
		vcf_header = new char[buflen]; vcf_header[0] = '\0';
		int64 title_buflen = buflen, title_len = 0;
		char*** header = new char** [g_input_row];

		for (int i = 0; i < g_input_row; ++i)
		{
			header[i] = new char* [g_input_col];

			for (int j = 0; j < g_input_col; ++j)
			{
				if (abs(g_format_val) == BCF)
				{
					FRead(readbuf.data, 1, 9, FILE_INFO[i][j].handle);
					if (*(uint*)readbuf.data != 0x02464342) Exit("\nError: Unsupported BCF format.\n");
					int headerlen = *(int*)(readbuf.data + 5);
					VLA_NEW(headerbuf, char, headerlen);
					FRead(headerbuf, 1, headerlen, FILE_INFO[i][j].handle);
					FILE_INFO[i][j].bcfheader = new BCFHEADER(headerbuf);
					VLA_DELETE(headerbuf);
					FSeek(FILE_INFO[i][j].handle, 9, SEEK_SET);
				}

				bool findhead = false;
				int64 readlen = 0;
				while (!FEof(FILE_INFO[i][j].handle))
				{
					readlen = readbuf.Gets(FILE_INFO[i][j].handle);
					if (LwrLineCmp("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	", readbuf.data) == 0)
					{
						if (readbuf.data[readlen - 1] == '\n') readbuf.data[--readlen] = '\0';
						if (readbuf.data[readlen - 1] == '\r') readbuf.data[--readlen] = '\0';
						header[i][j] = new char[readbuf.len + 1];
						strcpy(header[i][j], readbuf.data);
						findhead = true;
						if (abs(g_format_val) == BCF)
							FILE_INFO[i][j].bcfheader->nsample = CountChar(header[i][j], '\t') - 8;
						break;
					}
				}

				if (i == 0)
				{
					if (title_len + readlen + 1 >= title_buflen)
					{
						title_buflen <<= 1;
						char* tbuf2 = new char[title_buflen];
						strcpy(tbuf2, vcf_header);
						delete[] vcf_header;
					}
					char* tbuf = StrNextIdx(header[i][j], '\t', 9);
					int64 tlen = readlen - (tbuf - header[i][j]);
					memmove(vcf_header + title_len, tbuf, tlen); //concatenate individual names
					title_len += tlen;
					vcf_header[title_len] = '\0';
				}

				if ((abs(g_format_val) == BCF) && FILE_INFO[i][j].bcfheader->format_gqid != FILE_INFO[i][0].bcfheader->format_gqid)
					Exit("\nError: BCF format error, FORMAT GQ tag indices are different among BCFs. \n", FILE_INFO[i][j].name.c_str(), FILE_INFO[i][0].name.c_str());

				if ((abs(g_format_val) == BCF) && FILE_INFO[i][j].bcfheader->format_dpid != FILE_INFO[i][0].bcfheader->format_dpid)
					Exit("\nError: BCF format error, FORMAT DP tag indices are different among BCFs. \n", FILE_INFO[i][j].name.c_str(), FILE_INFO[i][0].name.c_str());

				if ((abs(g_format_val) == BCF) && FILE_INFO[i][j].bcfheader->format_adid != FILE_INFO[i][0].bcfheader->format_adid)
					Exit("\nError: BCF format error, FORMAT AD tag indices are different among BCFs. \n", FILE_INFO[i][j].name.c_str(), FILE_INFO[i][0].name.c_str());

				if (i && strcmp(header[i][j], header[0][j]))
					Exit("\nError: VCF/BCF format error, individuals in %s and %s are different. \n", FILE_INFO[i][j].name.c_str(), FILE_INFO[i][0].name.c_str());

				if (!findhead)
					Exit("\nError: VCF/BCF format error, cannot find genotype table header in %s.\n", FILE_INFO[i][j].name.c_str());

				if (abs(g_format_val) == BCF)
					FSeek(FILE_INFO[i][j].handle, 1, SEEK_CUR);
			}
		}

		for (int i = 0; i < g_input_row; ++i)
		{
			for (int j = 0; j < g_input_col; ++j)
				delete[] header[i][j];
			delete[] header[i];
		}
		delete[] header;

		//create new locus list
		new (&locus_list) VMEMORY(1);
		locus = (LOCUS*)locus_list.base_addr;

		//create individuals
		vcf_header[title_len + 1] = '\0';
		nind = CountChar(vcf_header, '\t');
		char* title = vcf_header + 1;
		ainds = new IND<REAL>*[nind];

		//construct geno_bucket
		geno_bucket.CreateBucket();
		if (ploidyinfer)
			ad_bucket.CreateBucket();

		int indc = 0;
		do
		{
			individual_memory->Alloc(ainds[indc], 1);
			new(ainds[indc]) IND<REAL>(title, indc);
			indc++;
		} while (*title);
		delete[] vcf_header; vcf_header = (char*)"\0";
		load_buf = new INCBUFFER[NBUF];

		AssignInds<REAL>();

		progress1 = 0;

		PROGRESS_VALUE = 0; PROGRESS_TOTAL = nloc; PROGRESS_NOUTPUTED2 = PROGRESS_NOUTPUTED; PROGRESS_NOUTPUTED = 0;

		load_complete = false;
		if (abs(g_format_val) == VCF)
			RunThreads(&LoadVCF<REAL>, &LoadVCFGuard<REAL>, NULL, TOTLEN_COMPRESS, TOTLEN_COMPRESS, "\nLoading VCF:\n", g_nthread_val, true, g_progress_val);
		else
			RunThreads(&LoadBCF<REAL>, &LoadBCFGuard<REAL>, NULL, TOTLEN_COMPRESS, TOTLEN_COMPRESS, "\nLoading BCF:\n", g_nthread_val, true, g_progress_val);

		locus = new LOCUS[nloc];
		SetVal(locus, (LOCUS*)locus_list.base_addr, nloc);

		if (g_missingploidy_b)
			ReplaceMissingGenotypes();

		//use allele depth as allele frequency, each ind represents a population
		if (ad)
		{
			ad++;
			progress1 = progress2 = 0;
			allele_freq_offset = new uint64[nloc];
			genotype_count_offset = new uint64[nloc];
			SetFF(allele_freq_offset, nloc);
			SetFF(genotype_count_offset, nloc);
			KT = GT = maxK = maxG = 0;

			for (int64 l = 0; l < nloc; ++l)
			{
				allele_freq_offset[l] = KT;
				genotype_count_offset[l] = GT;
				KT += locus[l].k;
				GT += locus[l].ngeno;
				maxK = Max((int)locus[l].k, maxK);
				maxG = Max((int)locus[l].ngeno, maxG);
			}

			// use allele frequency only 
			for (uint i = 0; i < pop<REAL>.size; ++i)
			{
				pop<REAL>[i].AllocFreq();
				delete[] pop<REAL>[i].loc_stat1;
				delete[] pop<REAL>[i].genocount;
				pop<REAL>[i].genocount = NULL;
				pop<REAL>[i].loc_stat1 = NULL;
			}

			load_complete = false;
			if (abs(g_format_val) == VCF)
				RunThreads(&LoadVCF<REAL>, &LoadVCFGuard<REAL>, NULL, TOTLEN_COMPRESS, TOTLEN_COMPRESS, "\nLoading VCF:\n", g_nthread_val, true, g_progress_val);
			else
				RunThreads(&LoadBCF<REAL>, &LoadBCFGuard<REAL>, NULL, TOTLEN_COMPRESS, TOTLEN_COMPRESS, "\nLoading BCF:\n", g_nthread_val, true, g_progress_val);

			delete[] allele_freq_offset;
			delete[] genotype_count_offset;
			allele_freq_offset = genotype_count_offset = NULL;
		}

		delete[] load_buf;
	}
	//non-vcf format: genepop|spagedi|cervus|arlequin|structure|polygene|polyrelatedness|genodive|plink
	else
	{
		if (g_format_val < 0)
		{
			TOTLEN_DECOMPRESS = 0;
			byte* buf = new byte[1024 * 1024];
			for (int i = 0; i < g_input_row; ++i)
			{
				for (int j = 0; j < g_input_col; ++j)
				{
					int64 offset = FTell(FILE_INFO[i][j].handle);
					while (!FEof(FILE_INFO[i][j].handle))
						FRead(buf, 1, 1024 * 1024, FILE_INFO[i][j].handle);
					FILE_INFO[i][j].decompressed_len = FTell(FILE_INFO[i][j].handle);
					TOTLEN_DECOMPRESS += FILE_INFO[i][j].decompressed_len;
					FSeek(FILE_INFO[i][j].handle, offset, SEEK_SET);
				}
			}
			delete[] buf;
		}

		if (ad) Exit("\nError: non-vcf and non-bcf format are incompatible with allelic depth (-ad) option.\n");
		if (usephase) Exit("\nError: non-vcf and non-bcf formats do not contains phased genotype that will be used in haplotype extraction or calculate r2/D' for sliding windows.\n");

		switch (abs(g_format_val))
		{
		case GENEPOP:
			RunThreads(&LoadGenepop<REAL>, NULL, NULL, TOTLEN_DECOMPRESS << 1, TOTLEN_DECOMPRESS << 1, "\nLoading genepop file:\n", 1, true);
			break;
		case SPAGEDI:
			RunThreads(&LoadSpagedi<REAL>, NULL, NULL, TOTLEN_DECOMPRESS << 1, TOTLEN_DECOMPRESS << 1, "\nLoading spagedi file:\n", 1, true);
			break;
		case CERVUS:
			RunThreads(&LoadCervus<REAL>, NULL, NULL, TOTLEN_DECOMPRESS << 1, TOTLEN_DECOMPRESS << 1, "\nLoading cervus file:\n", 1, true);
			break;
		case ARLEQUIN:
			RunThreads(&LoadArlequin<REAL>, NULL, NULL, TOTLEN_DECOMPRESS << 1, TOTLEN_DECOMPRESS << 1, "\nLoading structure file:\n", 1, true);
			break;
		case STRUCTURE:
			RunThreads(&LoadStructure<REAL>, NULL, NULL, TOTLEN_DECOMPRESS << 1, TOTLEN_DECOMPRESS << 1, "\nLoading structure file:\n", 1, true);
			break;
		case POLYGENE:
			RunThreads(&LoadPolyGene<REAL>, NULL, NULL, TOTLEN_DECOMPRESS << 1, TOTLEN_DECOMPRESS << 1, "\nLoading polygene file:\n", 1, true);
			break;
		case POLYRELATEDNESS:
			RunThreads(&LoadPolyRelatedness<REAL>, NULL, NULL, TOTLEN_DECOMPRESS << 1, TOTLEN_DECOMPRESS << 1, "\nLoading polyrelatedness file:\n", 1, true);
			break;
		case GENODIVE:
			RunThreads(&LoadGenoDive<REAL>, NULL, NULL, TOTLEN_DECOMPRESS << 1, TOTLEN_DECOMPRESS << 1, "\nLoading genodive file:\n", 1, true);
			break;
		case PLINK:
			RunThreads(&LoadPlink<REAL>, NULL, NULL, TOTLEN_DECOMPRESS << 1, TOTLEN_DECOMPRESS << 1, "\nLoading plink file:\n", 1, true);
			break;
		}
		AssignInds<REAL>();
	}

	// Close input files
	for (int i = 0; i < g_input_row; ++i)
		for (int j = 0; j < g_input_col; ++j)
		{
			FClose(FILE_INFO[i][j].handle);
			if (abs(g_format_val) == BCF && FILE_INFO[i][j].bcfheader)
			{
				delete FILE_INFO[i][j].bcfheader;
				FILE_INFO[i][j].bcfheader = NULL;
			}
		}

	CheckGenotypeId<REAL>();

	//compress locus into slocus
	MEMORY* nlocus_memory = new MEMORY[g_nthread_val];
	slocus = new SLOCUS[nloc];
	if (uselocpos) locus_pos = new uint64[nloc];

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
	for (int64 l = 0; l < nloc; ++l)
	{
		threadid = omp_get_thread_num();
		if (uselocpos) GetLocPos(l) = locus[l].pos;
		new(&slocus[l]) SLOCUS(nlocus_memory[threadid], locus[l]);
	}

	delete[] locus;
	delete[] locus_memory;
	locus_memory = nlocus_memory;
	useslocus = true;
	CheckGenotypeId<REAL>();
	EvaluationEnd("Load input file");
}

/* Create genotype table for non-vcf/bcf input files */
TARGET void CreateGenoIndexTable(GENO_WRITER* wt)
{
	geno_bucket.CreateBucketGT(locus);

	for (int64 l = 0; l < nloc; ++l)
		new(&wt[l]) GENO_WRITER(l);
}

/* load from Genepop input format */
THREAD2(LoadGenepop)
{
	//load genepop
	INCBUFFER buf;
	buf.Gets(FILE_INFO[0][0].handle); //ignore firstline
	bool findhead = false;

	//count locus
	int64 loc_pos = FTell(FILE_INFO[0][0].handle); //backup position
	while (!FEof(FILE_INFO[0][0].handle))
	{
		buf.Gets(FILE_INFO[0][0].handle);
		if (LwrLineCmp("pop", buf.data) == 0)
		{
			findhead = true;
			break;
		}
		nloc++;
	}

	if (!findhead)
		Exit("\nError: Cannot parse input file.\n");

	if (nloc == 0)
		Exit("\nError: there are no loci in this file.\n");

	int64 ind_pos = FTell(FILE_INFO[0][0].handle); //backup position
	while (!FEof(FILE_INFO[0][0].handle))
	{
		int64 rlen = buf.Gets(FILE_INFO[0][0].handle);
		if (LwrLineCmp("pop", buf.data) && rlen > 4ll)
			nind++;
	}

	if (nind == 0)
		Exit("\nError: there are no individuals in this file.\n");

	ainds = new IND<REAL>* [nind];
	/*   1   */
	ushort** gatab = new ushort * [nloc];
	GENOTYPE** gtab = new GENOTYPE * [nloc];
	GENO_WRITER* wt = new GENO_WRITER[nloc];

	locus = new LOCUS[nloc];
	nvcf_memory = new MEMORY[g_nthread_val];
	nvcf_gfid = new TABLE<HASH, uint>[nloc];

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
	for (int64 l = 0; l < nloc; ++l)
	{
		threadid = omp_get_thread_num();
		new(&nvcf_gfid[l]) TABLE<HASH, GENOTYPE*>(true, &nvcf_memory[threadid]);
	}

	for (int64 j = 1; j >= 0; --j)
	{
		FSeek(FILE_INFO[0][0].handle, ind_pos, SEEK_SET);
		for (int i = 0; i < nind; )
		{
			buf.Gets(FILE_INFO[0][0].handle);
			if (LwrLineCmp("pop", buf.data) == 0) continue;
			if (j != 0) individual_memory->Alloc(ainds[i], 1);
			/*   2   */
			new(ainds[i]) IND<REAL>(buf.data, j, i, gtab, gatab, wt);
			i++;
			PROGRESS_VALUE = FTell(FILE_INFO[0][0].handle) + (j ? 0 : FILE_INFO[0][0].decompressed_len);
		}

		if (j != 0)
		{
			/*   3   */
			FSeek(FILE_INFO[0][0].handle, loc_pos, SEEK_SET);
			for (int64 l = 0; l < nloc; ++l)
			{
				buf.Gets(FILE_INFO[0][0].handle);
				if (LwrLineCmp("pop", buf.data) == 0) continue;
				new(&locus[l]) LOCUS(locus_memory[threadid], buf.data, l, nvcf_gfid[l].size, gtab[l], gatab[l]);
			}
			CreateGenoIndexTable(wt);
		}
	}

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
	for (int64 l = 0; l < nloc; ++l)
		wt[l].FinishWrite();

	PROGRESS_VALUE = PROGRESS_CEND;
	IndexAlleleLength();

	/*   4   */
	delete[] wt;
	delete[] gatab;
	delete[] gtab;
	delete[] nvcf_gfid;
	delete[] nvcf_memory;
}

/* load from Spagedi input format */
THREAD2(LoadSpagedi)
{
	//load spagedi
	INCBUFFER buf;

	int64 tpos = 0;

	do
	{
		tpos = FTell(FILE_INFO[0][0].handle);
		buf.Gets(FILE_INFO[0][0].handle); //ignore firstline
	} while (buf.data[0] == '/' && buf.data[1] == '/');

	char* bufp = buf.data;
	for (int i = 0; i < 5; ++i)
		genotype_digit = ReadInteger(bufp);
	buf.Gets(FILE_INFO[0][0].handle); //-3
	int extracol = ReadIntegerKeep(buf.data) + 3;
	genotype_extracol = extracol;

	//count locus
	int64 loc_pos = FTell(FILE_INFO[0][0].handle); //backup position
	buf.Gets(FILE_INFO[0][0].handle);
	bufp = buf.data + strlen(buf.data);
	while (*bufp == '\n' || *bufp == '\r' || *bufp == '\t' || *bufp == ' ')
		*bufp-- = '\0';
	nloc = CountChar(buf.data, '\t') + 1 - 2 - extracol;

	if (nloc == 0)
		Exit("\nError: there are no loci in this file.\n");

	//count individual
	int64 ind_pos = FTell(FILE_INFO[0][0].handle); //backup position
	while (!FEof(FILE_INFO[0][0].handle))
	{
		buf.Gets(FILE_INFO[0][0].handle);
		if (LwrLineCmp("end", buf.data) == 0) break;
		nind++;
	}

	if (nind == 0)
		Exit("\nError: there are no individuals in this file.\n");

	FSeek(FILE_INFO[0][0].handle, ind_pos, SEEK_SET);

	ainds = new IND<REAL>* [nind];
	/*   1   */
	ushort** gatab = new ushort * [nloc];
	GENOTYPE** gtab = new GENOTYPE * [nloc];
	GENO_WRITER* wt = new GENO_WRITER[nloc];

	locus = new LOCUS[nloc];
	nvcf_memory = new MEMORY[g_nthread_val];
	nvcf_gfid = new TABLE<HASH, uint>[nloc];
#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
	for (int64 l = 0; l < nloc; ++l)
	{
		threadid = omp_get_thread_num();
		new(&nvcf_gfid[l]) TABLE<HASH, GENOTYPE*>(true, &nvcf_memory[threadid]);
	}

	for (int64 j = 1; j >= 0; --j)
	{
		FSeek(FILE_INFO[0][0].handle, ind_pos, SEEK_SET);
		for (int i = 0; i < nind; )
		{
			buf.Gets(FILE_INFO[0][0].handle);
			if (LwrLineCmp("end", buf.data) == 0) break;
			if (j != 0) individual_memory->Alloc(ainds[i], 1);
			/*   2   */
			new(ainds[i]) IND<REAL>(buf.data, j, i, gtab, gatab, wt);
			i++;
			PROGRESS_VALUE = FTell(FILE_INFO[0][0].handle) + (j ? 0 : TOTLEN_DECOMPRESS);
		}

		if (j != 0)
		{
			/*   3   */
			FSeek(FILE_INFO[0][0].handle, loc_pos, SEEK_SET);
			buf.Gets(FILE_INFO[0][0].handle);
			char* locus_name = StrNextIdx(buf.data, '\t', 2 + extracol) + 1;
			for (int64 l = 0; l < nloc; ++l)
			{
				new(&locus[l]) LOCUS(locus_memory[threadid], locus_name, l, nvcf_gfid[l].size, gtab[l], gatab[l]);
				locus_name = StrNextIdx(locus_name, '\t', 1) + 1;
			}
			CreateGenoIndexTable(wt);
		}
	}

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
	for (int64 l = 0; l < nloc; ++l)
		wt[l].FinishWrite();

	PROGRESS_VALUE = PROGRESS_CEND;
	IndexAlleleLength();

	/*   4   */
	delete[] wt;
	delete[] gatab;
	delete[] gtab;
	delete[] nvcf_gfid;
	delete[] nvcf_memory;
}

/* load from Cervus input format */
THREAD2(LoadCervus)
{
	//load cervus
	INCBUFFER buf;
	char* bufp = buf.data;

	//count locus
	int64 loc_pos = FTell(FILE_INFO[0][0].handle); //backup position
	buf.Gets(FILE_INFO[0][0].handle);

	bufp = buf.data + strlen(buf.data);
	while (*bufp == '\n' || *bufp == '\r' || *bufp == ',' || *bufp == ' ')
		*bufp-- = '\0';
	nloc = CountChar(buf.data, ',') + 1;
	int extracol = (nloc - ((nloc - 1) / 2) * 2) - 1;
	nloc = (nloc - 1) / 2;
	genotype_extracol = extracol;

	if (nloc == 0)
		Exit("\nError: there are no loci in this file.\n");

	//count individual
	int64 ind_pos = FTell(FILE_INFO[0][0].handle); //backup position
	while (!FEof(FILE_INFO[0][0].handle))
		if (buf.Gets(FILE_INFO[0][0].handle) > 3) nind++;

	if (nind == 0)
		Exit("\nError: there are no individuals in this file.\n");

	FSeek(FILE_INFO[0][0].handle, ind_pos, SEEK_SET);

	ainds = new IND<REAL>* [nind];
	/*   1   */
	ushort** gatab = new ushort * [nloc];
	GENOTYPE** gtab = new GENOTYPE * [nloc];
	GENO_WRITER* wt = new GENO_WRITER[nloc];

	locus = new LOCUS[nloc];
	nvcf_memory = new MEMORY[g_nthread_val];
	nvcf_gfid = new TABLE<HASH, uint>[nloc];
#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
	for (int64 l = 0; l < nloc; ++l)
	{
		threadid = omp_get_thread_num();
		new(&nvcf_gfid[l]) TABLE<HASH, GENOTYPE*>(true, &nvcf_memory[threadid]);
	}

	for (int64 j = 1; j >= 0; --j)
	{
		FSeek(FILE_INFO[0][0].handle, ind_pos, SEEK_SET);
		for (int i = 0; i < nind; )
		{
			if (buf.Gets(FILE_INFO[0][0].handle) <= 3)  continue;
			if (j != 0) individual_memory->Alloc(ainds[i], 1);
			/*   2   */
			new(ainds[i]) IND<REAL>(buf.data, j, i, gtab, gatab, wt);
			i++;
			PROGRESS_VALUE = FTell(FILE_INFO[0][0].handle) + (j ? 0 : TOTLEN_DECOMPRESS);
		}
		if (j != 0)
		{
			/*   3   */
			FSeek(FILE_INFO[0][0].handle, loc_pos, SEEK_SET);
			buf.Gets(FILE_INFO[0][0].handle);
			char* locus_name = StrNextIdx(buf.data, ',', extracol + 1) + 1;

			for (int64 l = 0; l < nloc; ++l)
			{
				char* locus_name2 = StrNextIdx(locus_name, ',', 1);
				*(locus_name2 - 1) = '\0';
				new(&locus[l]) LOCUS(locus_memory[threadid], locus_name, l, nvcf_gfid[l].size, gtab[l], gatab[l]);
				*(locus_name2 - 1) = 'A';
				locus_name = StrNextIdx(locus_name, ',', 2) + 1;
			}
			CreateGenoIndexTable(wt);
		}
	}

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
	for (int64 l = 0; l < nloc; ++l)
		wt[l].FinishWrite();

	PROGRESS_VALUE = PROGRESS_CEND;
	IndexAlleleLength();

	/*   4   */
	delete[] wt;
	delete[] gatab;
	delete[] gtab;
	delete[] nvcf_gfid;
	delete[] nvcf_memory;
}

/* load from Arlequin input format */
THREAD2(LoadArlequin)
{
	//load structure
	INCBUFFER buf;
	char*& bufp = buf.derived, * bufp2;

	//count locus
	int64 loc_pos = 0; //backup position
	bool header = false;
	for (;;)
	{
		buf.Gets(FILE_INFO[0][0].handle);
		bufp = buf.data;
		while (*bufp == '\t' || *bufp == ' ') bufp++;
		if (LwrLineCmp("SampleData", bufp) == 0)
		{
			header = true;
			loc_pos = FTell(FILE_INFO[0][0].handle);
			break;
		}
	}

	buf.Gets(FILE_INFO[0][0].handle);
	ReplaceChar(buf.data, '\t', ' ');
	bufp = buf.data;
	while (*bufp == ' ') bufp++; //name
	bufp = StrNextIdx(bufp, ' ', 1); while (*bufp == ' ') bufp++; //1
	bufp = StrNextIdx(bufp, ' ', 1); while (*bufp == ' ') bufp++; //genotype
	bufp2 = bufp + strlen(bufp) - 1;
	while (*bufp2 == ' ' || *bufp2 == '\n' || *bufp2 == '\r') *bufp2-- = '\0';
	nloc = CountChar(bufp, ' ') + 1;
	if (nloc == 0) Exit("\nError: there are no loci in this file.\n");

	int64 ind_pos = loc_pos;
	FSeek(FILE_INFO[0][0].handle, ind_pos, SEEK_SET);
	while (!FEof(FILE_INFO[0][0].handle))
	{
		//read name
		if (buf.Gets(FILE_INFO[0][0].handle) == 0) break;
		bufp = buf.data;
		ReplaceChar(buf.data, '\t', ' ');
		while (*bufp == ' ') bufp++;
		if (*bufp == '}') for (;;)
		{
			if (buf.Gets(FILE_INFO[0][0].handle) == 0) break;
			ReplaceChar(buf.data, '\t', ' ');
			bufp = buf.data;
			while (*bufp == ' ') bufp++;
			if (LwrLineCmp("SampleData", bufp) == 0) break;
		}
		else
		{
			if (buf.Gets(FILE_INFO[0][0].handle, bufp) == 0) break;
			nind++;
		}
	}

	if (nind == 0)
		Exit("\nError: there are no individuals in this file.\n");

	ainds = new IND<REAL>* [nind];
	/*   1   */
	ushort** gatab = new ushort * [nloc];
	GENOTYPE** gtab = new GENOTYPE * [nloc];
	GENO_WRITER* wt = new GENO_WRITER[nloc];

	locus = new LOCUS[nloc];
	nvcf_memory = new MEMORY[g_nthread_val];
	nvcf_gfid = new TABLE<HASH, uint>[nloc];
#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
	for (int64 l = 0; l < nloc; ++l)
	{
		threadid = omp_get_thread_num();
		new(&nvcf_gfid[l]) TABLE<HASH, GENOTYPE*>(true, &nvcf_memory[threadid]);
	}

	for (int64 j = 1; j >= 0; --j)
	{
		FSeek(FILE_INFO[0][0].handle, ind_pos, SEEK_SET);
		for (int i = 0; i < nind; )
		{
			//read name
			if (buf.Gets(FILE_INFO[0][0].handle) == 0) break;
			bufp = buf.data;
			ReplaceChar(buf.data, '\t', ' ');
			while (*bufp == ' ') bufp++;
			if (*bufp == '}') for (;;)
			{
				if (buf.Gets(FILE_INFO[0][0].handle) == 0) break;
				ReplaceChar(buf.data, '\t', ' ');
				bufp = buf.data;
				while (*bufp == ' ') bufp++;
				if (LwrLineCmp("SampleData", bufp) == 0) break;
			}
			else
			{
				bufp = buf.data + strlen(buf.data);
				if (buf.Gets(FILE_INFO[0][0].handle, bufp) == 0) break;
				if (j != 0) individual_memory->Alloc(ainds[i], 1);
				/*   2   */
				new(ainds[i]) IND<REAL>(buf.data, j, i, gtab, gatab, wt);
				i++;
				PROGRESS_VALUE = FTell(FILE_INFO[0][0].handle) + (j ? 0 : TOTLEN_DECOMPRESS);
			}
		}
		if (j != 0)
		{
			/*   3   */
			FSeek(FILE_INFO[0][0].handle, ind_pos, SEEK_SET);
			for (int64 l = 0; l < nloc; ++l)
				new(&locus[l]) LOCUS(locus_memory[threadid], (char*)NULL, l, nvcf_gfid[l].size, gtab[l], gatab[l]);
			CreateGenoIndexTable(wt);
		}
	}

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
	for (int64 l = 0; l < nloc; ++l)
		wt[l].FinishWrite();

	PROGRESS_VALUE = PROGRESS_CEND;
	IndexAlleleLength();

	/*   4   */
	delete[] wt;
	delete[] gatab;
	delete[] gtab;
	delete[] nvcf_gfid;
	delete[] nvcf_memory;
}

/* load from Structure input format */
THREAD2(LoadStructure)
{
	//load structure
	INCBUFFER buf;
	char* bufp = buf.data;

	//count locus
	int64 loc_pos = 0; //backup position
	buf.Gets(FILE_INFO[0][0].handle);
	int buflen2 = (int)strlen(buf.data);
	while (buf.data[buflen2 - 1] == ' ' || buf.data[buflen2 - 1] == '\r' || buf.data[buflen2 - 1] == '\n')
		buf.data[buflen2-- - 1] = '\0';
	ReplaceChar(buf.data, '\t', ' ');
	nloc = CountChar(buf.data, ' ') + 1;

	buf.Gets(FILE_INFO[0][0].handle);
	ReplaceChar(buf.data, '\t', ' ');
	int64 ncol = CountChar(buf.data, ' ') + 1;

	bool headerrow = false;
	if (nloc != ncol)
		headerrow = true;
	else
	{
		if (!g_extracol_b)
			Exit("\nError: no header row detected for structure input, please specify -g_extracol.\n");
		nloc = ncol - 1 - g_extracol_val;
	}

	g_extracol_val = ncol - nloc - 1;
	if (nloc == 0) Exit("\nError: there are no loci in this file.\n");

	FSeek(FILE_INFO[0][0].handle, 0, SEEK_SET);

	//jump 1
	buf.Gets(FILE_INFO[0][0].handle);
	ReplaceChar(buf.data, '\t', ' ');
	bufp = buf.data;

	//count individual
	if (headerrow == false) FSeek(FILE_INFO[0][0].handle, 0, SEEK_SET);
	int64 ind_pos = FTell(FILE_INFO[0][0].handle); //backup position

	while (!FEof(FILE_INFO[0][0].handle))
	{
		//read name
		int64 tpos1 = FTell(FILE_INFO[0][0].handle);
		if (buf.Gets(FILE_INFO[0][0].handle, 1024, buf.data) == 0) break;

		ReplaceChar(buf.data, '\t', ' ');
		bufp = buf.data;
		while (*bufp != ' ') bufp++; //skip name
		*bufp++ = '\0'; *bufp = '\0';
		FSeek(FILE_INFO[0][0].handle, tpos1, SEEK_SET);

		//read lines
		while (!FEof(FILE_INFO[0][0].handle))
		{
			//read one more section
			int64 tpos2 = FTell(FILE_INFO[0][0].handle);
			buf.Gets(FILE_INFO[0][0].handle, 1024, bufp);
			ReplaceChar(bufp, '\t', ' ');
			FSeek(FILE_INFO[0][0].handle, tpos2, SEEK_SET);
			if (LineCmpABterm(buf.data, bufp, ' ')) break; //not the same individual, break

			//the same individual, add length
			buf.Gets(FILE_INFO[0][0].handle, bufp);
			ReplaceChar(bufp, '\t', ' ');
		}

		//end line
		nind++;
	}

	if (nind == 0)
		Exit("\nError: there are no individuals in this file.\n");

	ainds = new IND<REAL>* [nind];
	/*   1   */
	ushort** gatab = new ushort * [nloc];
	GENOTYPE** gtab = new GENOTYPE * [nloc];
	GENO_WRITER* wt = new GENO_WRITER[nloc];

	locus = new LOCUS[nloc];
	nvcf_memory = new MEMORY[g_nthread_val];
	nvcf_gfid = new TABLE<HASH, uint>[nloc];
#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
	for (int64 l = 0; l < nloc; ++l)
	{
		threadid = omp_get_thread_num();
		new(&nvcf_gfid[l]) TABLE<HASH, GENOTYPE*>(true, &nvcf_memory[threadid]);
	}

	for (int64 j = 1; j >= 0; --j)
	{
		FSeek(FILE_INFO[0][0].handle, ind_pos, SEEK_SET);
		for (int i = 0; i < nind; )
		{
			//read name
			int64 tpos1 = FTell(FILE_INFO[0][0].handle);
			if (buf.Gets(FILE_INFO[0][0].handle, 1024, buf.data) == 0) break;
			ReplaceChar(buf.data, '\t', ' ');

			bufp = buf.data;
			while (*bufp != ' ') bufp++; //skip name
			*bufp++ = '\0'; *bufp = '\0';
			FSeek(FILE_INFO[0][0].handle, tpos1, SEEK_SET);

			//read lines
			char* bufp2 = bufp;
			while (!FEof(FILE_INFO[0][0].handle))
			{
				//read one more section
				int64 tpos2 = FTell(FILE_INFO[0][0].handle);
				buf.Gets(FILE_INFO[0][0].handle, 1024, bufp2);
				ReplaceChar(bufp2, '\t', ' ');
				FSeek(FILE_INFO[0][0].handle, tpos2, SEEK_SET);

				if (LineCmpABterm(buf.data, bufp2, ' '))
				{
					*bufp2 = '\0';
					break; //not the same individual, break
				}

				//the same individual, add length
				while (!FEof(FILE_INFO[0][0].handle))
				{
					buf.Gets(FILE_INFO[0][0].handle, bufp2);
					ReplaceChar(bufp2, '\t', ' ');
					while (*bufp2) bufp2++;
					if (*(bufp2 - 1) == '\n')
						break;
				}
			}

			//end line
			if (j != 0) individual_memory->Alloc(ainds[i], 1);
			/*   2   */
			new(ainds[i]) IND<REAL>(bufp, j, i, gtab, gatab, wt);
			i++;
			PROGRESS_VALUE = FTell(FILE_INFO[0][0].handle) + (j ? 0 : TOTLEN_DECOMPRESS);
		}
		if (j != 0)
		{
			/*   3   */
			FSeek(FILE_INFO[0][0].handle, loc_pos, SEEK_SET);
			buf.Gets(FILE_INFO[0][0].handle);
			ReplaceChar(buf.data, '\t', ' ');
			char* locus_name = buf.data;

			if (headerrow) for (int64 l = 0; l < nloc; ++l)
			{
				char* locus_name2 = StrNextIdx(locus_name, ' ', 1);
				if (locus_name2) *locus_name2 = '\0';
				new(&locus[l]) LOCUS(locus_memory[threadid], locus_name, l, nvcf_gfid[l].size, gtab[l], gatab[l]);
				if (locus_name2) *locus_name2 = ' ';
				locus_name = locus_name2 + 1;
			}
			else for (int64 l = 0; l < nloc; ++l)
				new(&locus[l]) LOCUS(locus_memory[threadid], (char*)NULL, l, nvcf_gfid[l].size, gtab[l], gatab[l]);
			CreateGenoIndexTable(wt);
		}
	}

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
	for (int64 l = 0; l < nloc; ++l)
		wt[l].FinishWrite();

	PROGRESS_VALUE = PROGRESS_CEND;
	IndexAlleleLength();

	/*   4   */
	delete[] wt;
	delete[] gatab;
	delete[] gtab;
	delete[] nvcf_gfid;
	delete[] nvcf_memory;
}

/* load from PolyGene input format */
THREAD2(LoadPolyGene)
{
	INCBUFFER buf;

	int64 tlen = buf.Gets(FILE_INFO[0][0].handle); //ignore firstline
	while (buf.data[tlen] == '\t') buf.data[tlen--] = '\0';

	nloc = CountChar(buf.data, '\t') - 2;  //ind pop ploidy
	if (nloc == '\0')
		Exit("\nError: there are no loci in this file.\n");

	//Read locus
	char* locus_name = StrNextIdx(buf.data, '\t', 3) + 1;
	char* locus_name_b = new char[strlen(locus_name) + 1];
	strcpy(locus_name_b, locus_name);
	locus_name = locus_name_b;

	//count individual
	int64 ind_pos = FTell(FILE_INFO[0][0].handle); //backup position
	while (!FEof(FILE_INFO[0][0].handle))
		if (buf.Gets(FILE_INFO[0][0].handle) > 5)
			nind++;

	if (nind == 0)
		Exit("\nError: there are no individuals in this file.\n");

	FSeek(FILE_INFO[0][0].handle, ind_pos, SEEK_SET);

	ainds = new IND<REAL>* [nind];
	/*   1   */
	ushort** gatab = new ushort * [nloc];
	GENOTYPE** gtab = new GENOTYPE * [nloc];
	GENO_WRITER* wt = new GENO_WRITER[nloc];

	locus = new LOCUS[nloc];
	nvcf_memory = new MEMORY[g_nthread_val];
	nvcf_gfid = new TABLE<HASH, uint>[nloc];

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
	for (int64 l = 0; l < nloc; ++l)
	{
		threadid = omp_get_thread_num();
		new(&nvcf_gfid[l]) TABLE<HASH, GENOTYPE*>(true, &nvcf_memory[threadid]);
	}

	for (int64 j = 1; j >= 0; --j)
	{
		ind_pos = FTell(FILE_INFO[0][0].handle);
		for (int i = 0; i < nind; )
		{
			buf.Gets(FILE_INFO[0][0].handle);
			if (j != 0) individual_memory->Alloc(ainds[i], 1);
			/*   2   */
			new(ainds[i]) IND<REAL>(buf.data, j, i, gtab, gatab, wt);
			i++;
			PROGRESS_VALUE = FTell(FILE_INFO[0][0].handle) + (j ? 0 : TOTLEN_DECOMPRESS);
		}
		FSeek(FILE_INFO[0][0].handle, ind_pos, SEEK_SET);
		if (j != 0)
		{
			/*   3   */
			for (int64 l = 0; l < nloc; ++l)
			{
				new(&locus[l]) LOCUS(locus_memory[threadid], locus_name, l, nvcf_gfid[l].size, gtab[l], gatab[l]);
				locus_name = StrNextIdx(locus_name, '\t', 1) + 1;
			}
			CreateGenoIndexTable(wt);
		}
	}

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
	for (int64 l = 0; l < nloc; ++l)
		wt[l].FinishWrite();

	PROGRESS_VALUE = PROGRESS_CEND;
	IndexAlleleLength();

	/*   4   */
	delete[] locus_name_b;
	delete[] wt;
	delete[] gatab;
	delete[] gtab;
	delete[] nvcf_gfid;
	delete[] nvcf_memory;
}

/* load from PolyRelatedness input format */
THREAD2(LoadPolyRelatedness)
{
	//load polyrelatedness
	INCBUFFER buf;

	int64 tpos = 0;
	do
	{
		tpos = FTell(FILE_INFO[0][0].handle);
		buf.Gets(FILE_INFO[0][0].handle); //ignore firstline
	} while (buf.data[0] == '/' && buf.data[1] == '/');

	char* bufp = buf.data;
	genotype_digit = ReadInteger(bufp);
	genotype_missing = ReadInteger(bufp);
	genotype_missing = ReadInteger(bufp);
	genotype_ambiguous = ReadInteger(bufp);

	//count locus
	bool findhead = false;
	while (!FEof(FILE_INFO[0][0].handle))
	{
		buf.Gets(FILE_INFO[0][0].handle);
		if (LwrLineCmp("//genotype", buf.data) == 0)
		{
			findhead = true;
			break;
		}
	}
	int64 loc_pos = FTell(FILE_INFO[0][0].handle); //backup position
	buf.Gets(FILE_INFO[0][0].handle);

	bufp = buf.data + strlen(buf.data);
	while (*bufp == '\n' || *bufp == '\r' || *bufp == '\t' || *bufp == ' ')
		*bufp-- = '\0';
	nloc = CountChar(buf.data, '\t') + 1 - 2;

	if (nloc == '\0')
		Exit("\nError: there are no loci in this file.\n");

	// alloc memory
	FSeek(FILE_INFO[0][0].handle, loc_pos, SEEK_SET);

	//jump 2
	buf.Gets(FILE_INFO[0][0].handle);
	char* locus_name = StrNextIdx(buf.data, '\t', 2) + 1, * locus_name_b = NULL;
	if (locus_name != (char*)1)
	{
		locus_name_b = new char[strlen(locus_name) + 1];
		strcpy(locus_name_b, locus_name);
		locus_name = locus_name_b;
	}

	//count individual
	int64 ind_pos = FTell(FILE_INFO[0][0].handle); //backup position
	while (!FEof(FILE_INFO[0][0].handle))
	{
		buf.Gets(FILE_INFO[0][0].handle);
		if (LwrLineCmp("//end of file", buf.data) == 0) break;
		nind++;
	}

	if (nind == 0)
		Exit("\nError: there are no individuals in this file.\n");

	ainds = new IND<REAL>* [nind];
	/*   1   */
	ushort** gatab = new ushort * [nloc];
	GENOTYPE** gtab = new GENOTYPE * [nloc];
	GENO_WRITER* wt = new GENO_WRITER[nloc];

	locus = new LOCUS[nloc];
	nvcf_memory = new MEMORY[g_nthread_val];
	nvcf_gfid = new TABLE<HASH, uint>[nloc];

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
	for (int64 l = 0; l < nloc; ++l)
	{
		threadid = omp_get_thread_num();
		new(&nvcf_gfid[l]) TABLE<HASH, GENOTYPE*>(true, &nvcf_memory[threadid]);
	}

	for (int64 j = 1; j >= 0; --j)
	{
		FSeek(FILE_INFO[0][0].handle, ind_pos, SEEK_SET);
		for (int i = 0; i < nind; )
		{
			buf.Gets(FILE_INFO[0][0].handle);
			if (LwrLineCmp("//end of file", buf.data) == 0) break;
			if (j != 0) individual_memory->Alloc(ainds[i], 1);
			/*   2   */
			new(ainds[i]) IND<REAL>(buf.data, j, i, gtab, gatab, wt);
			i++;
			PROGRESS_VALUE = FTell(FILE_INFO[0][0].handle) + (j ? 0 : TOTLEN_DECOMPRESS);
		}

		if (j != 0)
		{
			/*   3   */
			FSeek(FILE_INFO[0][0].handle, loc_pos, SEEK_SET);
			for (int64 l = 0; l < nloc; ++l)
			{
				new(&locus[l]) LOCUS(locus_memory[threadid], locus_name == (char*)1 ? NULL : locus_name, l, nvcf_gfid[l].size, gtab[l], gatab[l]);
				if (locus_name != (char*)1) locus_name = StrNextIdx(locus_name, '\t', 1) + 1;
			}
			CreateGenoIndexTable(wt);
		}
	}

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
	for (int64 l = 0; l < nloc; ++l)
		wt[l].FinishWrite();

	PROGRESS_VALUE = PROGRESS_CEND;
	IndexAlleleLength();

	/*   4   */
	if (locus_name_b) delete[] locus_name_b;
	delete[] wt;
	delete[] gatab;
	delete[] gtab;
	delete[] nvcf_gfid;
	delete[] nvcf_memory;
}

/* load from Genodive input format */
THREAD2(LoadGenoDive)
{
	INCBUFFER buf;

	//ignore firstline
	buf.Gets(FILE_INFO[0][0].handle);

	//second line
	buf.Gets(FILE_INFO[0][0].handle);
	char* bufp = buf.data;

	nind = ReadInteger(bufp);
	int _npop = ReadInteger(bufp) + 0;
	nloc = ReadInteger(bufp);
	maxploidy = ReadInteger(bufp) + 0;
	genotype_digit = ReadInteger(bufp);

	if (nloc == '\0')
		Exit("\nError: there are no loci in this file.\n");
	if (nind == '\0')
		Exit("\nError: there are no individuals in this file.\n");

	//skip p lines
	for (int p = 0; p < _npop; ++p)
		buf.Gets(FILE_INFO[0][0].handle);

	//header row
	int64 loc_pos = FTell(FILE_INFO[0][0].handle); //backup position
	buf.Gets(FILE_INFO[0][0].handle);

	//calculate extra columns
	buf.Gets(FILE_INFO[0][0].handle);
	char* name0 = StrNextChar(StrNextSpace(StrNextChar(buf.data)));

	//whether exist extra columns
	g_extracol_val = name0[0] >= '0' && name0[0] <= '9' ? 1 : 0;

	// alloc memory
	FSeek(FILE_INFO[0][0].handle, loc_pos, SEEK_SET);

	//jump 2 + g_extracol_val columns
	buf.Gets(FILE_INFO[0][0].handle);
	char* locus_name = StrNextIdx(buf.data, '\t', 2 + g_extracol_val) + 1, * locus_name_b = NULL;
	if (locus_name != (char*)1)
	{
		locus_name_b = new char[strlen(locus_name) + 1];
		strcpy(locus_name_b, locus_name);
		locus_name = locus_name_b;
	}

	int64 ind_pos = FTell(FILE_INFO[0][0].handle); //backup position
	ainds = new IND<REAL>* [nind];

	/*   1   */
	ushort** gatab = new ushort * [nloc];
	GENOTYPE** gtab = new GENOTYPE * [nloc];
	GENO_WRITER* wt = new GENO_WRITER[nloc];

	locus = new LOCUS[nloc];
	nvcf_memory = new MEMORY[g_nthread_val];
	nvcf_gfid = new TABLE<HASH, uint>[nloc];

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
	for (int64 l = 0; l < nloc; ++l)
	{
		threadid = omp_get_thread_num();
		new(&nvcf_gfid[l]) TABLE<HASH, GENOTYPE*>(true, &nvcf_memory[threadid]);
	}

	for (int64 j = 1; j >= 0; --j)
	{
		FSeek(FILE_INFO[0][0].handle, ind_pos, SEEK_SET);
		for (int i = 0; i < nind; )
		{
			buf.Gets(FILE_INFO[0][0].handle);
			if (j != 0) individual_memory->Alloc(ainds[i], 1);
			/*   2   */
			new(ainds[i]) IND<REAL>(buf.data, j, i, gtab, gatab, wt);
			i++;
			PROGRESS_VALUE = FTell(FILE_INFO[0][0].handle) + (j ? 0 : TOTLEN_DECOMPRESS);
		}

		if (j != 0)
		{
			/*   3   */
			FSeek(FILE_INFO[0][0].handle, loc_pos, SEEK_SET);
			for (int64 l = 0; l < nloc; ++l)
			{
				new(&locus[l]) LOCUS(locus_memory[threadid], locus_name == (char*)1 ? NULL : locus_name, l, nvcf_gfid[l].size, gtab[l], gatab[l]);
				if (locus_name != (char*)1) locus_name = StrNextIdx(locus_name, '\t', 1) + 1;
			}
			CreateGenoIndexTable(wt);
		}
	}

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
	for (int64 l = 0; l < nloc; ++l)
		wt[l].FinishWrite();

	PROGRESS_VALUE = PROGRESS_CEND;
	IndexAlleleLength();

	/*   4   */
	if (locus_name_b) delete[] locus_name_b;
	delete[] wt;
	delete[] gatab;
	delete[] gtab;
	delete[] nvcf_gfid;
	delete[] nvcf_memory;
}

/* load from Plink input format */
THREAD2(LoadPlink)
{
	//load plink
	INCBUFFER buf, locbuf;
	char* bufp = buf.data;
	char mapname[4096];
	strcpy(mapname, g_input_val.c_str());
	FILE* fmap = NULL;
	int64 loc_pos = 0; //backup position

	char* map_ext = mapname + strlen(mapname);
	while (*map_ext != '.') map_ext--; 
	if ((map_ext[0] == 'g' || map_ext[0] == 'G') && (map_ext[1] == 'z' || map_ext[1] == 'Z') && map_ext[2] == '\0')
	{
		map_ext--; 
		while (*map_ext != '.') 
			map_ext--;
	}

	map_ext++;
	const char* extensions[] = {"map", "map.gz", "txt", "loc"};
	for (int i = 0; i < 4; ++i)
	{
		sprintf(map_ext, "%s", extensions[i]);
		if (FileExists(mapname))
		{
			fmap = FOpen(mapname, "rb");
			for (;;)
			{
				loc_pos = FTell(fmap);
				locbuf.Gets(fmap);
				if (locbuf.data[0] != '#') break;
			}
			break;
		}
	}

	//count locus
	int64 ind_pos = 0;
	for (;;)
	{
		ind_pos = FTell(FILE_INFO[0][0].handle); //backup position
		buf.Gets(FILE_INFO[0][0].handle);
		if (buf.data[0] != '#') break;
	}
	ReplaceChar(buf.data, '\t', ' ');

	bufp = buf.data + strlen(buf.data);
	while (*bufp == '\n' || *bufp == '\r' || *bufp == ',' || *bufp == ' ')
		*bufp-- = '\0';
	nloc = (CountChar(buf.data, ' ') + 1 - 6) / 2;
	genotype_extracol = 0;

	if (nloc == 0)
		Exit("\nError: there are no loci in this file.\n");

	//count individual
	FSeek(FILE_INFO[0][0].handle, ind_pos, SEEK_SET);
	while (!FEof(FILE_INFO[0][0].handle))
		if (buf.Gets(FILE_INFO[0][0].handle) > 3) nind++;

	if (nind == 0)
		Exit("\nError: there are no individuals in this file.\n");

	FSeek(FILE_INFO[0][0].handle, ind_pos, SEEK_SET);

	ainds = new IND<REAL>*[nind];
	/*   1   */
	ushort** gatab = new ushort * [nloc];
	GENOTYPE** gtab = new GENOTYPE * [nloc];
	GENO_WRITER* wt = new GENO_WRITER[nloc];

	locus = new LOCUS[nloc];
	nvcf_memory = new MEMORY[g_nthread_val];
	nvcf_gfid = new TABLE<HASH, uint>[nloc];
#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
	for (int64 l = 0; l < nloc; ++l)
	{
		threadid = omp_get_thread_num();
		new(&nvcf_gfid[l]) TABLE<HASH, GENOTYPE*>(true, &nvcf_memory[threadid]);
	}

	for (int64 j = 1; j >= 0; --j)
	{
		FSeek(FILE_INFO[0][0].handle, ind_pos, SEEK_SET);
		for (int i = 0; i < nind; )
		{
			if (buf.Gets(FILE_INFO[0][0].handle) <= 3)  continue;
			if (j != 0) individual_memory->Alloc(ainds[i], 1);
			/*   2   */
			new(ainds[i]) IND<REAL>(buf.data, j, i, gtab, gatab, wt);
			i++;
			PROGRESS_VALUE = FTell(FILE_INFO[0][0].handle) + (j ? 0 : TOTLEN_DECOMPRESS);
		}
		if (j != 0)
		{
			/*   3   */

			if (fmap)
			{
				FSeek(fmap, loc_pos, SEEK_SET);

				for (int64 l = 0; l < nloc; ++l)
				{
					locbuf.Gets(fmap);
					ReplaceChar(locbuf.data, '\t', ' ');
					char* locus_name = StrNextIdx(locbuf.data, ' ', 1) + 1;
					char* locus_name2 = StrNextIdx(locus_name, ' ', 1);
					*locus_name2 = '\0';
					new(&locus[l]) LOCUS(locus_memory[threadid], locus_name, l, nvcf_gfid[l].size, gtab[l], gatab[l]);
				}
			}
			else
			{
				char tlocname[100];
				for (int64 l = 0; l < nloc; ++l)
				{
					sprintf(tlocname, "loc%lld", l + 1);
					new(&locus[l]) LOCUS(locus_memory[threadid], tlocname, l, nvcf_gfid[l].size, gtab[l], gatab[l]);
				}
			}

			CreateGenoIndexTable(wt);

		}
	}

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
	for (int64 l = 0; l < nloc; ++l)
		wt[l].FinishWrite();

	PROGRESS_VALUE = PROGRESS_CEND;
	IndexAlleleLength();

	/*   4   */
	if (fmap) FClose(fmap);
	delete[] wt;
	delete[] gatab;
	delete[] gtab;
	delete[] nvcf_gfid;
	delete[] nvcf_memory;
}

/* Indexing integer alleles to zero based for non-vcf input, with allele identifier being the size */
TARGET void IndexAlleleLength()
{
	//paralleled

	TABLE<ushort, ushort>* aftab = new TABLE<ushort, ushort>[g_nthread_val];
	for (int i = 0; i < g_nthread_val; ++i)
		new(&aftab[i]) TABLE<ushort, ushort>(true, NULL);

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
	for (int64 l = 0; l < nloc; ++l)
	{
		threadid = omp_get_thread_num();
		TABLE<ushort, ushort>& atab = aftab[threadid];
		atab.Clear();

		GENOTYPE* gtab = locus[l].GetGtab();
		int ngenotype = locus[l].ngeno;

		for (int gi = 0; gi < ngenotype; ++gi)
		{
			GENOTYPE& gt = gtab[gi];
			if (gt.Nalleles() == 0) continue;
			int v = gt.Ploidy();
			ushort* als = gt.GetAlleleArray();
			for (int vi = 0; vi < v; ++vi)
				atab.PushIndex(als[vi]);
		}

		ushort*& alen = locus[l]._alen;
		auto bucket = atab.bucket;
		auto index = atab.index;
		locus_memory[threadid].Alloc(locus[l]._alen, atab.size);
		for (uint i = 0; i < atab.size; ++i)
			alen[i] = bucket[index[i]].key;

		locus[l].k = (ushort)atab.size;
		locus[l].flag_alen = 1;
		if (locus[l].k < 2) locus[l].flag_pass = false;

		for (int gi = 0; gi < ngenotype; ++gi)
		{
			GENOTYPE& gt = gtab[gi];
			if (gt.Nalleles() == 0) continue;
			int v = gt.Ploidy();
			ushort alleles[N_MAX_PLOIDY] = { 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF };
			ushort* als = gt.GetAlleleArray();

			//mapping alleles
			for (int vi = 0; vi < v; ++vi)
				alleles[vi] = atab[als[vi]];

			//non-vcf genotypes are unphased, can be sort
			Sort(alleles, v);

			ushort* nullgatab = NULL;

			//sorted, reconstruct alleles, hash is invalid
			new(&gt) GENOTYPE(nullgatab, alleles, v);
		}

	}

	delete[] aftab;
}

/* Process lines from memory */
THREAD2(LoadBCF)
{
	TABLE<HASH, uint> gfid(false, NULL);
	int64& pval = *(int64*)&PROGRESS_VALUE;
	int64& cend = *(int64*)&PROGRESS_CEND;

	for (int64 ii = 0; !load_complete || pval != cend; ii++)
	{
		while (ii >= progress1)
		{
			Sleep(SLEEP_TIME_TINY);
			if (pval == cend && load_complete) return;
		}

		int64 avail_val = ii * 4 + 2;

		if (state_lock[ii % NBUF].compare_exchange_strong(avail_val, avail_val + 1)) {

			gfid.Clear();

			LOADLINE& line = *(LOADLINE*)load_buf[ii % NBUF].data;
			byte byteflag = line.flag;

			char* str = line.data;
			char* dst = StrNextIdx(str, '\t', 9) + 1;//GT
			byte type = *dst++;

			//higher 4 bits, number of elements, i.e., ploidy level
			int maxv = type >> 4;
			if (maxv == 0xF) maxv = ReadTypedInt(dst);
			if (maxv >= 10) Exit("\nError: do not support a ploidy level greater than 10.\n");

			//lower 4 bits, number of bytes or data type
			type &= 0xF;
			int asize = 0, nploidy = 0;
			switch (type)
			{
			case 1: asize = 1; break;
			case 2: asize = 2; break;
			case 3: asize = 4; Exit("\nError: do not support 4 byte allele identifier in BCF format, the maximum number of alleles allowed is 65535.\n"); break; break;
			default: Exit("\nError: GT format is should encode as integers.\n"); break;
			}

			bool usedploidy[N_MAX_PLOIDY + 1] = { 0 };
			char* genostr = dst;
			for (int i = 0; i < nind; ++i)
			{
				ushort alleles[N_MAX_PLOIDY] = { 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF };
				bool phased = true;
				int v = 0;
				ReadBCFGenoString(genostr, alleles, phased, v, asize, maxv, ii, NULL);
				genostr += asize * maxv;

				//will use phased genotype, do not need unphased genotype Warning
				//if (usephase && !phased)
					//locus[ii].flag_pass = false;

				//add missing at ploidy v
				if (!usedploidy[v])
				{
					int tid = gfid.size;
					gfid[missing_hash[v]] = tid;
					usedploidy[v] = true;
					if (g_missingploidy_b) nploidy++;
				}

				//add gfid if not missing
				if (alleles[0] != 0xFFFF && !(usephase && !phased))
				{
					//phased, do not sort alleles
					if (!usephase) Sort(alleles, v);

					HASH ha = HashGenotype(alleles, v);
					if (!gfid.ContainsKey(ha))
					{
						//add to hash and alleles array size
						int tid = gfid.size;
						gfid[ha] = tid;
						locus[ii].gasize += v + GetNalleles(alleles, v);
					}
				}
			}

			//add missing at ploidy 1-10
			byte pad_missing[N_MAX_PLOIDY + 1] = { 0 };
			if (g_missingploidy_b)
				for (int v = 1; v <= N_MAX_PLOIDY; ++v)
				{
					if (g_missingploidy_val[v] && !gfid.ContainsKey(missing_hash[v])) //bug fixed
					{
						int tid = gfid.size;
						gfid[missing_hash[v]] = tid;
						usedploidy[v] = true;
						pad_missing[v] = true;
						nploidy++;
					}
				}

			OFFSET off = geno_bucket.AddOffsetGT(CeilLog2((int64)gfid.size), ii);

			/*******************************************************************/

			//alloc in one piece
			GENOTYPE* gtab = NULL;  ushort* gatab;
			if (ad != 2)
				new(&locus[ii]) LOCUS(locus_memory[threadid], str, (uint64)byteflag, gfid.size, gtab, gatab);//LoadBCF
			else
				str = StrNextIdx(str, '\t', 9) + 1;
			locus[ii].nploidy = nploidy;

			char* gt = str + 1, * buf1 = str + maxv * asize * nind + 1, * gq = NULL, * dp = NULL, * adval = NULL;
			int gqlen = 0, dplen = 0, adlen = 0;
			if (locus[ii].gqid != 0xFFFF)
			{
				gq = buf1;
				switch (*gq++ & 0xF)
				{
				case 1: gqlen = 1; break;
				case 2: gqlen = 2; break;
				case 3: gqlen = 4; break;
				}
				buf1 += gqlen * nind + 1;
			}
			if (locus[ii].dpid != 0xFFFF)
			{
				dp = buf1;
				switch (*dp++ & 0xF)
				{
				case 1: dplen = 1; break;
				case 2: dplen = 2; break;
				case 3: dplen = 4; break;
				}
				buf1 += dplen * nind + 1;
			}
			if (locus[ii].adid != 0xFFFF)
			{
				adval = buf1;
				switch (*adval++ & 0xF)
				{
				case 1: adlen = 1; break;
				case 2: adlen = 2; break;
				case 3: adlen = 4; break;
				}
				buf1 += adlen * nind + 1;
			}

			uint* depth = ploidyinfer ? new uint[locus[ii].k * nind] : NULL, * depth2 = depth;

			//add missing for used ploidy
			TABLE<HASH, GENOTYPE*>& gftab = locus[ii].gftab;
			for (int v = 1; v <= N_MAX_PLOIDY; ++v)
				if (usedploidy[v])
					gftab[missing_hash[v]] = new(gtab++) GENOTYPE(gatab, missing_genotype[v]);

			//geno_bucket.offset[ii] = off;
			GENO_WRITER wt(ii);
			for (int j = 0; j < nind; ++j)
				ainds[j]->AddBCFGenotype(ii, gt, gq, dp, adval, maxv, asize, gqlen, dplen, adlen, depth, gfid, gtab, gatab, wt);
			wt.FinishWrite();

			//add missing at ploidy 1-10
			if (g_missingploidy_b)
				for (int v = 1; v <= N_MAX_PLOIDY; ++v)
					if (pad_missing[v])
						gftab[missing_hash[v]] = new(gtab++) GENOTYPE(gatab, missing_genotype[v]); //bug fixed


			if (ploidyinfer)
			{
				//add allele depth
				ad_bucket.AddOffsetAD((uint64)CeilLog2((int)locus[ii].maxdepth + 1), ii, locus[ii].k);

				IND<REAL>::SetAlleleDepth(ii, depth2, nind * locus[ii].k, 0u);
				delete[] depth2;
			}

			//geno_bucket.offset[ii] = off;
			PROGRESS_VALUE++;

			state_lock[ii % NBUF] = (ii + NBUF) * 4;
		}
	}
}

/* Read lines from BCF file */
 THREAD2(LoadBCFGuard)
{
	int64& ii = progress1 = 0;
	VLA_NEW(next_line, int64, g_input_col);
	VLA_NEW(original_position, int64, g_input_row * g_input_col);
	VLA_NEW(previous_offset, int64, g_input_row * g_input_col);
	VLA_NEW(jloc, int64, g_input_col);

	INCBUFFER buf, bufgt, bufgq, bufdp, bufad, locinfo;

	for (int i = 0; i < g_input_row; ++i)
	{
		for (int j = 0; j < g_input_col; ++j)
		{
			next_line[j] = original_position[i * g_input_col + j] = FTell(FILE_INFO[i][j].handle);
			previous_offset[i * g_input_col + j] = FOffset(FILE_INFO[i][j].handle);
			PROGRESS_VALUE += previous_offset[i * g_input_col + j];
			jloc[j] = 0;
		}

		//read a whole text line across bcfs in the same row
		for (;;)
		{
			while ((state_lock[ii % NBUF] >> 2) != ii)
				Sleep(SLEEP_TIME_TINY);

			state_lock[ii % NBUF] = ii * 4 + 1;

			float qual = FLT_MAX;
			bool pass = true;
			char *&bufj = buf.derived, *& bufjgt = bufgt.derived, *& bufjgq = bufgq.derived, *& bufjdp = bufdp.derived, *& bufjad = bufad.derived;
			bufj = buf.data; bufjgt = bufgt.data; bufjgq = bufgq.data; bufjdp = bufdp.data; bufjad = bufad.data;

			uint64 byteflag = 0;
			bool breakflag = false;

			for (int j = 0; j < g_input_col; ++j)
			{
				if (FEof(FILE_INFO[i][j].handle))
				{
					breakflag = true;
					break;
				}

				uint64 coffset = FOffset(FILE_INFO[i][j].handle);
				PROGRESS_VALUE += coffset - previous_offset[i * g_input_col + j];
				previous_offset[i * g_input_col + j] = coffset;

				int infolen = FGetUint(FILE_INFO[i][j].handle);

				if (FEof(FILE_INFO[i][j].handle))
				{
					breakflag = true;
					break;
				}

				int64 fmtpos = next_line[j] + 8 + infolen;
				next_line[j] = fmtpos + FGetUint(FILE_INFO[i][j].handle);

				//extend buf
				buf.Expand(bufj - buf.data + infolen + 4096);

				//CHROM
				int chromid = FGetUint(FILE_INFO[i][j].handle);
				if (j == 0) { AppendString(bufj, FILE_INFO[i][j].bcfheader->contig_name[chromid]); *bufj++ = '\t'; }

				//POS
				int pos = FGetUint(FILE_INFO[i][j].handle);
				if (j == 0) { AppendInt(bufj, pos); *bufj++ = '\t'; }

				//RLEN
				int rlen = FGetUint(FILE_INFO[i][j].handle);

				//QUAL
				float f1 = FGetFloat(FILE_INFO[i][j].handle);
				qual = Min(f1, qual);

				//#INFO
				int ninfo = (int)FGetUshort(FILE_INFO[i][j].handle);

				//#ALLELE
				int nallele = (int)FGetUshort(FILE_INFO[i][j].handle);

				//#SAMPLE
				int nsamp = (int)FGetUshort(FILE_INFO[i][j].handle) | ((int)FGetByte(FILE_INFO[i][j].handle) << 16);

				//#FORMAT
				int nformat = (int)FGetByte(FILE_INFO[i][j].handle);

				//ID
				int typelen = (int)FGetByte(FILE_INFO[i][j].handle) >> 4;

				if (typelen == 0xF)
					typelen = FGetTypedInt(FILE_INFO[i][j].handle);
				if (j == 0)
				{
					FGet(FILE_INFO[i][j].handle, bufj, typelen);
					bufj += typelen;
					*bufj++ = '\t';
				}
				else
					FSeek(FILE_INFO[i][j].handle, typelen, SEEK_CUR);

				//REF ALT
				for (int a = 0; a < nallele; ++a)
				{
					typelen = FGetByte(FILE_INFO[i][j].handle) >> 4;
					if (typelen == 0xF) typelen = FGetTypedInt(FILE_INFO[i][j].handle);
					if (j == 0)
					{
						FGet(FILE_INFO[i][j].handle, bufj, typelen);
						bufj += typelen;
						*bufj++ = a == 0 ? '\t' : ',';
					}
					else
						FSeek(FILE_INFO[i][j].handle, typelen, SEEK_CUR);
				}
				if (j == 0) *(bufj - 1) = '\t';

				if (g_input_col > 1)
				{
					if (j != 0)
					{
						jloc[j]++;

						int tinfolen = bufj - buf.data + 1;
						if (memcmp(locinfo.data, buf.data, tinfolen))
							Exit("\nError: variant information in %s and %s mismatch, at lines %d:\n%s\nvs\n%s\n",
								FILE_INFO[i][0].name.c_str(), FILE_INFO[i][j].name.c_str(), jloc[j], locinfo.data, bufj);
					}
					else if (g_input_col > 1)
					{
						jloc[0]++;

						locinfo.Move(buf.data, bufj - buf.data + 1);
					}
				}

				if (f_qual_b)
				{
					//qual filter
					if (f_qual_b && IsNormal(qual) && qual < f_qual_min && qual > f_qual_max)
						byteflag |= 0x4;
				}

				//FILTER
				int filtertype = (int)FGetByte(FILE_INFO[i][j].handle);
				int numfilter = filtertype >> 4; filtertype &= 0xF;

				if (f_original_b && f_original_val == 1)
				{
					for (int k = 0; k < numfilter; ++k)
					{
						//original filter
						int filterid = 0;
						switch (filtertype)
						{
						case 1: filterid = (int)FGetByte(FILE_INFO[i][j].handle); break;
						case 2: filterid = (int)FGetUshort(FILE_INFO[i][j].handle); break;
						case 3: filterid = (int)FGetUint(FILE_INFO[i][j].handle); break;
						}

						if (filterid != FILE_INFO[i][j].bcfheader->filter_passidx)
						{
							byteflag |= 0x2;
							pass = false;
						}
					}
				}

				//SKIP INFO
				FSeek(FILE_INFO[i][j].handle, fmtpos, SEEK_SET);
				int gtid = FILE_INFO[i][j].bcfheader->format_gtid, gqid = FILE_INFO[i][j].bcfheader->format_gqid, dpid = FILE_INFO[i][j].bcfheader->format_dpid, adid = FILE_INFO[i][j].bcfheader->format_adid;

				/*******************************************************************/

				//FORMAT
				for (int a = 0; a < nformat; ++a)
				{
					int fkey = (int)FGetTypedInt(FILE_INFO[i][j].handle);
					int ftype = (int)FGetByte(FILE_INFO[i][j].handle);
					int flen = (ftype >> 4 == 0xF) ? FGetTypedInt(FILE_INFO[i][j].handle) : (ftype >> 4);

					switch (ftype & 0xF)
					{
					case 1: break;
					case 2: flen <<= 1; break;
					case 3: flen <<= 2; break;
					case 5: flen <<= 2; break;
					case 7: break;
					}
					int datalen = FILE_INFO[i][j].bcfheader->nsample * flen;

					if (fkey == gtid)
					{
						bufgt.Expand(bufjgt - bufgt.data + datalen + 4096);
						if (j == 0)
						{
							*(byte*)bufjgt++ = (byte)ftype;
							if ((ftype >> 4) == 0xF) AppendTypedInt(bufjgt, flen);
						}
						FGet(FILE_INFO[i][j].handle, bufjgt, datalen);
						bufjgt += datalen;
					}
					else if (fkey == gqid)
					{
						bufgq.Expand(bufjgq - bufgq.data + datalen + 4096);
						if (j == 0)
						{
							*(byte*)bufjgq++ = (byte)ftype;
							if ((ftype >> 4) == 0xF) AppendTypedInt(bufjgq, flen);
						}
						FGet(FILE_INFO[i][j].handle, bufjgq, datalen);
						bufjgq += datalen;
					}
					else if (fkey == dpid)
					{
						bufdp.Expand(bufjdp - bufdp.data + datalen + 4096);
						if (j == 0)
						{
							*(byte*)bufjdp++ = (byte)ftype;
							if ((ftype >> 4) == 0xF) AppendTypedInt(bufjdp, flen);
						}
						FGet(FILE_INFO[i][j].handle, bufjdp, datalen);
						bufjdp += datalen;
					}
					else if (fkey == adid)
					{
						bufad.Expand(bufjad - bufad.data + datalen + 4096);
						if (j == 0)
						{
							*(byte*)bufjad++ = (byte)ftype;
							if ((ftype >> 4) == 0xF) AppendTypedInt(bufjad, flen);
						}
						FGet(FILE_INFO[i][j].handle, bufjad, datalen);
						bufjad += datalen;
					}
					else
						FSeek(FILE_INFO[i][j].handle, datalen, SEEK_CUR);
				}

				FSeek(FILE_INFO[i][j].handle, next_line[j], SEEK_SET);

			}

			if (breakflag) 
				break;

			sprintf(bufj, "%lf", qual); while (*bufj) bufj++;
			if (pass) AppendString(bufj, "\tPASS\t\tGT");
			else AppendString(bufj, "\t.\t\tGT");

			if (bufjgq > bufgq.data) AppendString(bufj, ":GQ");
			if (bufjdp > bufdp.data) AppendString(bufj, ":DP");
			if (bufjad > bufad.data) AppendString(bufj, ":AD");

			*bufj++ = '\t';
			int64 writelen = sizeof(LOADLINE) +
				(bufj - buf.data) + 
				(bufjgt - bufgt.data) +
				(bufjgq - bufgq.data) +
				(bufjdp - bufdp.data) +
				(bufjad - bufad.data);

			if (bufjgq > bufgq.data) byteflag |= 0x08; //hasgq
			if (bufjdp > bufdp.data) byteflag |= 0x10; //hasdp
			if (bufjad > bufad.data) byteflag |= 0x20; //hasad

			if (load_buf[ii % NBUF].len < writelen	)
				load_buf[ii % NBUF].Expand(writelen + writelen / 2);

			LOADLINE& line = *(LOADLINE*)load_buf[ii % NBUF].data;
			line.size = writelen;
			line.flag = byteflag;

			char* writep = line.data;
			writep[writelen - 1] = '\0';
			SetVal(writep, buf.data,   bufj   - buf.data);     writep += bufj   - buf.data;
			SetVal(writep, bufgt.data, bufjgt - bufgt.data);   writep += bufjgt - bufgt.data;
			SetVal(writep, bufgq.data, bufjgq - bufgq.data);   writep += bufjgq - bufgq.data;
			SetVal(writep, bufdp.data, bufjdp - bufdp.data);   writep += bufjdp - bufdp.data;
			SetVal(writep, bufad.data, bufjad - bufad.data);   writep += bufjad - bufad.data;

			nloc++;
			if (&locus[nloc] > (LOCUS*)locus_list.tail_addr)
				locus_list.Alloc((byte*)&locus[nloc]);

			PROGRESS_TOTAL++;
			PROGRESS_CEND++;

			state_lock[ii % NBUF] = ii * 4 + 2;

			ii++;
		}

		for (int j = 0; j < g_input_col; ++j)
		{
			FGetUint(FILE_INFO[i][j].handle);
			if (!FEof(FILE_INFO[i][j].handle))
				Exit("\nError: number of variants in %s mismatch that in %s.\n", FILE_INFO[i][j].name.c_str(), FILE_INFO[i][0].name.c_str());
			if (jloc[j] != jloc[0])
				Exit("\nError: number of variants in %s mismatch that in %s.\n", FILE_INFO[i][j].name.c_str(), FILE_INFO[i][0].name.c_str());

			uint64 coffset = FOffset(FILE_INFO[i][j].handle);
			PROGRESS_VALUE += coffset - previous_offset[i * g_input_col + j];
			previous_offset[i * g_input_col + j] = coffset;
		}
	}

	TOTLEN_DECOMPRESS = 0;
	for (int i = 0; i < g_input_row; ++i)
		for (int j = 0; j < g_input_col; ++j)
		{
			FILE_INFO[i][j].decompressed_len = FTell(FILE_INFO[i][j].handle);
			TOTLEN_DECOMPRESS += FILE_INFO[i][j].decompressed_len;
			FSeek(FILE_INFO[i][j].handle, original_position[i * g_input_col + j], SEEK_SET);
		}

	VLA_DELETE(next_line);
	VLA_DELETE(original_position);
	VLA_DELETE(previous_offset);
	VLA_DELETE(jloc);
	load_complete = true;
}

/* Process lines from memory */
THREAD2(LoadVCF)
{
	TABLE<HASH, uint> gfid(false, NULL);
	bool fmt_filter[1024] = { 0 };
	int64& pval = *(int64*)&PROGRESS_VALUE;
	int64& cend = *(int64*)&PROGRESS_CEND;
	MEMORY& memory = locus_memory[threadid];

	for (int64 ii = 0; !load_complete || pval != cend; ii++)
	{
		while (ii >= progress1)
		{
			Sleep(SLEEP_TIME_TINY);
			if (pval == cend && load_complete) return;
		}

		int64 avail_val = ii * 4 + 2;
		if (state_lock[ii % NBUF].compare_exchange_strong(avail_val, avail_val + 1)) {

			gfid.Clear();
			int fmt_count = 0;

			LOADLINE& line = *(LOADLINE*)load_buf[ii % NBUF].data;
			int64 res = line.size;
			int64 byteflag = line.flag;

			char* str = line.data;
			char* src1 = str;
			char* src2 = StrNextIdx(src1, '\t', 7, res - (src1 - str)) + 1;	    //INFO
			char* dst = src2;
			src1 = StrNextIdx(src2, '\t', 1, res - (src2 - str)) + 1;			//FORMAT
			src2 = StrNextIdx(src1, '\t', 1, res - (src1 - str));				//GENOTYPES

			*dst++ = '\t';
			fmt_count = 0;
			bool writed = false;
			int gttag = -1;

			while (src1 > (char*)1 && src1 < src2)
			{
				ushort tagid = *(ushort*)src1;
				char tagnext = src1[2];

				if ((tagnext == ':' || tagnext == '\t') &&
					(tagid == 0x5447 || //GT
						tagid == 0x5044 || //DP
						tagid == 0x5147 || //GQ
					(tagid == 0x4441 && (ad || ploidyinfer)))) //AD
				{
					if (writed) *dst++ = ':';
					*(ushort*)dst = *(ushort*)src1;
					dst += 2; src1 += 2;
					if (tagid == 0x5447) gttag = fmt_count;
					fmt_filter[fmt_count++] = writed = true;
				}
				else
					fmt_filter[fmt_count++] = false;

				while (*src1 != ':' && *src1 != '\t')
					src1++;
				src1++;
			}

			*dst++ = '\t';
			src1 = src2 + 1;
			bool usedploidy[N_MAX_PLOIDY + 1] = { 0 };
			int count = 0, nploidy = 0;

			while (*src1)
			{
				count++;
				writed = false;
				for (int i = 0; i < fmt_count; ++i)
				{
					if (fmt_filter[i])
					{
						char* genostr = src1;
						if (writed) *dst++ = ':';

						while (*src1 != ':' && *src1 != '\0' && *src1 != '\t' && *src1 != '\n')
							*dst++ = *src1++;

						if (gttag == i)
						{
							// genotype
							int len = (int)(src1 - genostr);
							int v = CountPloidy(genostr, len);
							if (v <= 0 || v > N_MAX_PLOIDY)
								Exit("\nError: ploidy level greater than %d in individual %s at %d-th locus.\n", N_MAX_PLOIDY, ainds[count]->name, ii + 1);

							//add missing at ploidy v
							if (!usedploidy[v])
							{
								int tid = gfid.size;
								gfid[missing_hash[v]] = tid;
								usedploidy[v] = true;
								if (g_missingploidy_b) nploidy++;
							}

							//add gfid if not missing
							if (genostr[0] != '.' && !(usephase && ContainsChar(genostr, '/', len)))
							{
								ushort alleles[N_MAX_PLOIDY] = { 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF };
								ReadVCFGenoString(alleles, genostr, v, ii, NULL);

								//phased, do not sort alleles
								if (!usephase) Sort(alleles, v);

								HASH ha = HashGenotype(alleles, v);
								if (!gfid.ContainsKey(ha))
								{
									//add to hash and alleles array size
									int tid = gfid.size;
									gfid[ha] = tid;
									locus[ii].gasize += v + GetNalleles(alleles, v);
								}
							}
						}

						if (*src1 == '\t' || *src1 == '\n')
						{
							src1++;
							break;
						}
						else if (*src1 == '\0') break;

						src1++;
						writed = true;
					}
					else
					{
						//for (register char c = *src1; c != ':' && c != '\0' && c != '\t' && c != '\n'; c = *++src1);
						while (*src1 != ':' && *src1 != '\0' && *src1 != '\t' && *src1 != '\n') src1++;
						if (*src1 == '\t' || *src1 == '\n')
						{
							src1++;
							break;
						}
						else if (*src1 == '\0') break;
						src1++;
					}
				}
				if (count == nind) break;
				*dst++ = '\t';
			}

			*dst++ = '\n';
			*dst++ = '\0';

			//add missing at ploidy 1-10
			byte pad_missing[N_MAX_PLOIDY + 1] = { 0 };
			if (g_missingploidy_b)
				for (int v = 1; v <= N_MAX_PLOIDY; ++v)
				{
					if (g_missingploidy_val[v] && !gfid.ContainsKey(missing_hash[v])) //bug fixed
					{
						int tid = gfid.size;
						gfid[missing_hash[v]] = tid;
						usedploidy[v] = true;
						pad_missing[v] = true;
						nploidy++;
					}
				}

			OFFSET off = geno_bucket.AddOffsetGT(CeilLog2((int64)gfid.size), ii);
            
			//alloc in one piece
			GENOTYPE* gtab = NULL;  ushort* gatab;
			if (ad != 2)
				new(&locus[ii]) LOCUS(memory, str, (uint64)byteflag, gfid.size, gtab, gatab);
			else
				str = StrNextIdx(str, '\t', 9) + 1;

			locus[ii].nploidy = nploidy;
			uint* depth = ploidyinfer ? new uint[locus[ii].k * nind] : NULL, *depth2 = depth;

			//geno_bucket.offset[ii] = off;
			TABLE<HASH, GENOTYPE*>& gftab = locus[ii].gftab;
			GENO_WRITER wt(ii);
			for (int j = 0; j < nind; ++j)
				ainds[j]->AddVCFGenotype(str, ii, depth, gfid, gtab, gatab, wt);
			wt.FinishWrite();

			//add missing at ploidy 1-10
			if (g_missingploidy_b)
				for (int v = 1; v <= N_MAX_PLOIDY; ++v)
					if (pad_missing[v])
						gftab[missing_hash[v]] = new(gtab++) GENOTYPE(gatab, missing_genotype[v]); //bug fixed

			if (ploidyinfer)
			{
				ad_bucket.AddOffsetAD((uint64)CeilLog2((int)locus[ii].maxdepth + 1), ii, locus[ii].k);
				IND<REAL>::SetAlleleDepth(ii, depth2, nind * locus[ii].k, 0u);
				delete[] depth2;
			}

			//geno_bucket.offset[ii] = off;
			PROGRESS_VALUE++;

			state_lock[ii % NBUF] = (ii + NBUF) * 4;
		}
	}
}

/* Read lines from VCF file */
THREAD2(LoadVCFGuard)
{
	int64& ii = progress1 = 0;

	VLA_NEW(original_position, int64, g_input_row * g_input_col);
	VLA_NEW(vbuf, VCFBUFFER, g_input_col);

	INCBUFFER locinfo;
	for (int j = 0; j < g_input_col; ++j)
		vbuf[j].Expand(LINE_BUFFER + 1);

	for (int i = 0; i < g_input_row; ++i)
	{
		for (int j = 0; j < g_input_col; ++j)
		{
			original_position[i * g_input_col + j] = FTell(FILE_INFO[i][j].handle);
			vbuf[j].Read1(FILE_INFO[i][j].handle);
		}

		for (;;)
		{
			while ((state_lock[ii % NBUF] >> 2) != ii)
				Sleep(SLEEP_TIME_TINY);

			state_lock[ii % NBUF] = ii * 4 + 1;

			char* first_line = NULL, * first_fmt = NULL;
			int64 byte_flag = 0, writelen = sizeof(LOADLINE);

			for (int j = 0; j < g_input_col; ++j)
			{
				char*& this_line = vbuf[j].line_pos;
				this_line += vbuf[j].line_len;
				int64 remain_size = vbuf[j].data_size - (vbuf[j].line_pos - vbuf[j].data);
				char* next_line = StrNextIdx(this_line, '\n', 1, remain_size) + 1;

				//do not have next line in buffer, move remaining data to the beginning of buffer and read a new batch of lines
				if (next_line == (char*)1)
				{
					vbuf[j].Read2(FILE_INFO[i][j].handle);
					remain_size = vbuf[j].data_size - (vbuf[j].line_pos - vbuf[j].data);
					next_line = StrNextIdx(this_line, '\n', 1, remain_size) + 1;

					if (next_line == (char*)1)
						next_line = this_line + vbuf[j].data_size;
				}

				vbuf[j].line_len = next_line - this_line;

				//nothing to read
				if (vbuf[j].line_len == 0 || this_line[0] == '\0')
					break;

				char* qualstr = StrNextIdx(this_line, '\t', 5, remain_size) + 1;

				if (g_input_col > 1)
				{
					vbuf[j].jloc++;

					if (j)
					{
						//compare CHROM POS ID REF ALT
						if (LineCmp(locinfo.data, this_line))
						{
							Exit("\nError: variant information in %s and %s mismatch, at lines %d:\n%s\nvs\n%s\n",
								FILE_INFO[i][0].name.c_str(), FILE_INFO[i][j].name.c_str(), vbuf[j].jloc, locinfo.data, this_line);
						}

					}
					else
					{
						//compare CHROM POS ID REF ALT
						locinfo.Move(this_line, qualstr - this_line - 1);
					}
				}

				if (j == 0)
				{
					//save first format
					int tinfolen = qualstr - this_line - 1;
					first_line = this_line;
					first_fmt = StrNextIdx(this_line + tinfolen + 1, '\t', 3, remain_size - tinfolen) + 1;
				}

				int nlinebreak = 0;
				//set line break to \0
				if (this_line[vbuf[j].line_len - 1] == '\n')
				{
					this_line[vbuf[j].line_len - 1] = '\0';
					nlinebreak = 1;
				}

				if (this_line[vbuf[j].line_len - 2] == '\r')
				{
					this_line[vbuf[j].line_len - 2] = '\0';
					nlinebreak = 2;
				}

				//test qual filter
				if (f_qual_b)
				{
					double qual = -1;
					if (*qualstr != '.')
						qual = ReadDoubleKeep(qualstr);
					if (f_qual_b && qual != -1 && qual < f_qual_min && qual > f_qual_max)
						byte_flag |= 0x4;	//fail to pass qual filter
					//avg_qual += qual;
				}

				//test original filter
				if (f_original_b && f_original_val == 1)
				{
					char* filtstr = StrNextIdx(qualstr, '\t', 1, remain_size - (qualstr - this_line)) + 1;
					if (LwrLineCmp("pass", filtstr))
						byte_flag |= 0x2;	//fail to pass original filter
				}

				char* fmt = StrNextIdx(qualstr, '\t', 3, remain_size) + 1;

				vbuf[j].geno_pos = vbuf[j].line_pos;
				vbuf[j].geno_len = vbuf[j].line_len - nlinebreak;

				//compare FORMAT
				if (g_input_col > 1 && j)
				{
					vbuf[j].geno_pos = StrNextIdx(fmt, '\t', 1, remain_size);
					vbuf[j].geno_len = vbuf[j].line_len - nlinebreak - (vbuf[j].geno_pos - vbuf[j].line_pos);

					if (LineCmpAterm(first_fmt, fmt, '\t'))
					{
						*StrNextIdx(fmt, '\t', 1) = '\0';
						*StrNextIdx(first_fmt, '\t', 1) = '\0';
						*StrNextIdx(first_line, '\t', 3) = '\0';
						Exit("\nError: FORMAT field in %s and %s are different, at line %d, variant %s: \n%s\nvs\n%s\n",
							FILE_INFO[i][0].name.c_str(), FILE_INFO[i][j].name.c_str(), ii + 1, first_line, first_fmt, fmt);
					}
				}

				writelen += vbuf[j].geno_len;
			}

			//read complete, break
			if (vbuf[0].line_len == 0 || first_line == NULL || first_line[0] == '\0')
				break;

			if (load_buf[ii % NBUF].len < writelen)
				load_buf[ii % NBUF].Expand(writelen + writelen / 2);

			LOADLINE& line = *(LOADLINE*)load_buf[ii % NBUF].data;
			line.size = writelen;
			line.flag = byte_flag;
			char* writebuf = line.data;
			for (int j = 0; j < g_input_col; ++j)
			{
				memcpy(writebuf, vbuf[j].geno_pos, vbuf[j].geno_len);
				writebuf += vbuf[j].geno_len;
			}
			*writebuf++ = '\0';

			nloc++;
			if (&locus[nloc] > (LOCUS*)locus_list.tail_addr)
				locus_list.Alloc((byte*)&locus[nloc]);

			PROGRESS_TOTAL++;
			PROGRESS_CEND++;

			state_lock[ii % NBUF] = ii * 4 + 2;
			ii++;
		}
	}

	TOTLEN_DECOMPRESS = 0;
	for (int i = 0; i < g_input_row; ++i)
		for (int j = 0; j < g_input_col; ++j)
		{
			FILE_INFO[i][j].decompressed_len = FTell(FILE_INFO[i][j].handle);
			TOTLEN_DECOMPRESS += FILE_INFO[i][j].decompressed_len;
			FSeek(FILE_INFO[i][j].handle, original_position[i * g_input_col + j], SEEK_SET);
		}

	VLA_DELETE(vbuf);
	VLA_DELETE(original_position);
	load_complete = true;
}

/* Replace missing genotype with vmin */
TARGET void ReplaceMissingGenotypes()
{
	atomic<uint> maxg{ 0 };
#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
	for (int64 l = 0; l < nloc; ++l)
		if (locus[l].ngeno > maxg)
			AtomicMax(maxg, locus[l].ngeno);

	VLA_NEW(ploidytab, int, maxg * g_nthread_val);
	VLA_NEW(nallelestab, int, maxg * g_nthread_val);

	atomic<byte>* vmin_typed = new atomic<byte>[nind];
	atomic<byte>* vmin_all = new atomic<byte>[nind];
	SetVal((byte*)vmin_typed, (byte)100, nind);
	SetVal((byte*)vmin_all, (byte)100, nind);

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
	for (int64 l = 0; l < nloc; ++l)
	{
		threadid = omp_get_thread_num();

		int* nalleles = nallelestab + maxg * threadid;
		int* ploidy = ploidytab + maxg * threadid;

		GENOTYPE* gtab = locus[l].GetGtab();
		int ngeno = locus[l].ngeno;

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
			byte v = (byte)ploidy[gid];

			if (nalleles[gid] && v < vmin_typed[j])
				AtomicMin(vmin_typed[j], v);
			if (v < vmin_all[j])
				AtomicMin(vmin_all[j], v);
		}
	}

	//replace missing genotypes
#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
	for (int64 l = 0; l < nloc; ++l)
	{
		threadid = omp_get_thread_num();

		GENOTYPE* gtab = locus[l].GetGtab();
		int ngeno = locus[l].ngeno;
		int missingid[N_MAX_PLOIDY + 1] = { 0 };
		for (int gid = 0; gid < ngeno; ++gid)
		{
			GENOTYPE& gt = gtab[gid];
			if (gt.Nalleles() == 0)
				missingid[gt.Ploidy()] = gid;
		}

		GENO_READER rt(0, l);
		GENO_WRITER wt(l);

		for (int i = 0; i < nind; ++i)
		{
			int gid = rt.Read();
			GENOTYPE& gt = gtab[gid];
			wt.Write(gt.Nalleles() == 0 ? 
				missingid[vmin_typed[i] == 100 ? 
					((byte)vmin_all[i] == 100 ? (byte)2 : (byte)vmin_all[i]) : (byte)vmin_typed[i]] : gid);
		}
		wt.FinishWrite();
	}

	delete[] vmin_typed;
	delete[] vmin_all;
	VLA_DELETE(ploidytab);
	VLA_DELETE(nallelestab);
}

/* Create individual from genepop */
template<typename REAL>
TARGET void IND<REAL>::genepop(char* t, bool iscount, GENOTYPE** gtab, ushort** gatab, GENO_WRITER* wt)
{
	//two rounds, first round obtain the number of genotypes
	name = t;
	char* genstr = NULL;

	genstr = StrNextIdx(t, ',', 1);
	char* tname = genstr - 1;
	*genstr++ = '\0';

	if (!iscount)
	{
		while (*tname == ' ') *tname-- = '\0';
		individual_memory->Alloc(tname, (int)strlen(name) + 1);
		strcpy(tname, name);
		name = tname;
	}
	else name = NULL;

	int ploidy = 2;
	for (int64 l = 0; l < nloc; ++l)
	{
		TABLE<HASH, GENOTYPE*>& gftab = locus[l].gftab;
		TABLE<HASH, uint>& gfid = nvcf_gfid[l];
		bool ismissing = false;

		while (*genstr == ' ' || *genstr == '\t') genstr++;
		char* gend = genstr;
		while (*gend >= '0' && *gend <= '9') gend++;
		*gend = '\0';

		int64 glen = strlen(genstr);
		ushort alleles[N_MAX_PLOIDY] = { 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF };
		if (glen == 4)
		{
			char b = genstr[2];
			genstr[2] = '\0';
			alleles[0] = (ushort)ReadIntegerKeep(genstr);
			genstr[2] = b;
			alleles[1] = (ushort)ReadIntegerKeep(genstr + 2);
		}
		else if (glen == 6)
		{
			char b = genstr[3];
			genstr[3] = '\0';
			alleles[0] = (ushort)ReadIntegerKeep(genstr);
			genstr[3] = b;
			alleles[1] = (ushort)ReadIntegerKeep(genstr + 3);
		}
		else
			Exit("\nError: Format error in individual %s at locus %s.\n", name, locus[l].GetName());

		if (alleles[0] == '\0' && alleles[1] == '\0')
		{
			ismissing = true;
			SetFF(alleles, ploidy);
		}
		Sort(alleles, ploidy);//unphase

		if (iscount)
		{
			HASH mha = missing_hash[ploidy];
			if (!gfid.ContainsKey(mha))
			{
				int tid = gfid.size;
				gfid[mha] = tid;
				locus[l].gasize += 0;
			}
			if (!ismissing)
			{
				HASH hash = HashGenotype(alleles, ploidy);
				if (!gfid.ContainsKey(hash))
				{
					int tid = gfid.size;
					gfid[hash] = tid;
					locus[l].gasize += ploidy + GetNalleles(alleles, ploidy);
				}
			}
		}
		else
		{
			HASH mha = missing_hash[ploidy];
			uint mid = gfid[mha];
			if (!gftab.ContainsKey(mha))
				gftab[mha] = new(gtab[l]++) GENOTYPE(gatab[l], missing_genotype[ploidy]);

			HASH hash = HashGenotype(alleles, ploidy);
			uint gid = gfid[hash];
			if (!gftab.ContainsKey(hash))
				gftab[hash] = new(gtab[l]++) GENOTYPE(gatab[l], alleles, ploidy);

			wt[l].Write(genotype_filter && f_ploidy_b && (f_ploidy_min > ploidy || ploidy > f_ploidy_max) ? mid : gid);
		}

		genstr = gend + 1;
	}
}

/* Create individual from spagedi */
template<typename REAL>
TARGET void IND<REAL>::spagedi(char* t, bool iscount, GENOTYPE** gtab, ushort** gatab, GENO_WRITER* wt)
{
	name = t;
	char* genstr = NULL;

	genstr = StrNextIdx(t, '\t', 1);
	char* tname = genstr - 1;
	*genstr++ = '\0';

	if (!iscount)
	{
		while (*tname == ' ') *tname-- = '\0';
		individual_memory->Alloc(tname, (int)strlen(name) + 1);
		//tname = new char[strlen(name) + 1];
		strcpy(tname, name);
		name = tname;
	}
	else name = NULL;

	genstr = StrNextIdx(genstr, '\t', genotype_extracol + 1) + 1; //1 (pop) + extracols

	int ndigit = genotype_digit;
	int maxv = 0;
	for (int64 l = 0; l < nloc; ++l)
	{
		TABLE<HASH, GENOTYPE*>& gftab = locus[l].gftab;
		TABLE<HASH, uint>& gfid = nvcf_gfid[l];

		while (*genstr == '\t') genstr++;
		ushort alleles[N_MAX_PLOIDY] = { 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF };
		int v = 0;
		bool ismissing = false;
		while (*genstr != '\t' && *genstr != '\r' && *genstr != '\n' && *genstr != '\0')
		{
			alleles[v] = (ushort)ReadIntegerSpagedi(genstr, ndigit);
			if (alleles[v] == 0)
				ismissing = true;
			if (++v > N_MAX_PLOIDY)
				Exit("\nError: ploidy level greater than %d in individual %s at locus %s.\n", N_MAX_PLOIDY, name, locus[l].GetName());
		}
		maxv = Max(v, maxv);

		if (!iscount && ismissing)
		{
			v = maxploidy;
			SetFF(alleles, v);
		}
		Sort(alleles, v);//unphase

		if (iscount)
		{
			HASH mha = missing_hash[v];
			if (!gfid.ContainsKey(mha))
			{
				int tid = gfid.size;
				gfid[mha] = tid;
				locus[l].gasize += 0;
			}
			if (!ismissing)
			{
				HASH hash = HashGenotype(alleles, v);
				if (!gfid.ContainsKey(hash))
				{
					int tid = gfid.size;
					gfid[hash] = tid;
					locus[l].gasize += v + GetNalleles(alleles, v);
				}
			}
		}
		else
		{
			HASH mha = missing_hash[v];
			uint mid = gfid[mha];
			if (!gftab.ContainsKey(mha))
				gftab[mha] = new(gtab[l]++) GENOTYPE(gatab[l], missing_genotype[v]);

			HASH hash = HashGenotype(alleles, v);
			uint gid = gfid[hash];
			if (!gftab.ContainsKey(hash))
				gftab[hash] = new(gtab[l]++) GENOTYPE(gatab[l], alleles, v);

			wt[l].Write(genotype_filter && f_ploidy_b && (f_ploidy_min > v || v > f_ploidy_max) ? mid : gid);
		}
	}
	maxploidy = (byte)maxv;
}

/* Create individual from cervus */
template<typename REAL>
TARGET void IND<REAL>::cervus(char* t, bool iscount, GENOTYPE** gtab, ushort** gatab, GENO_WRITER* wt)
{
	int extracol = genotype_extracol;
	name = t;
	char* genstr = NULL;

	genstr = StrNextIdx(t, ',', 1);
	char* tname = genstr - 1;
	*genstr++ = '\0';

	if (!iscount)
	{
		while (*tname == ' ') *tname-- = '\0';
		individual_memory->Alloc(tname, (int)strlen(name) + 1);
		strcpy(tname, name);
		name = tname;
	}
	else name = NULL;

	genstr = StrNextIdx(genstr, ',', extracol) + 1;

	int v = 2;
	for (int64 l = 0; l < nloc; ++l)
	{
		TABLE<HASH, GENOTYPE*>& gftab = locus[l].gftab;
		TABLE<HASH, uint>& gfid = nvcf_gfid[l];
		bool ismissing = false;

		ushort alleles[N_MAX_PLOIDY] = { 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF };
		if (*genstr != ',' && *genstr != '\r' && *genstr != '\n' && *genstr != '\0')
			alleles[0] = (ushort)ReadInteger(genstr);
		genstr++;
		if (*genstr != ',' && *genstr != '\r' && *genstr != '\n' && *genstr != '\0')
			alleles[1] = (ushort)ReadInteger(genstr);
		genstr++;

		if (alleles[0] == 0xFFFF && alleles[1] == 0xFFFF)
		{
			ismissing = true;
			SetFF(alleles, v);
		}
		Sort(alleles, v);//unphase

		if (iscount)
		{
			HASH mha = missing_hash[v];
			if (!gfid.ContainsKey(mha))
			{
				int tid = gfid.size;
				gfid[mha] = tid;
				locus[l].gasize += 0;
			}
			if (!ismissing)
			{
				HASH hash = HashGenotype(alleles, v);
				if (!gfid.ContainsKey(hash))
				{
					int tid = gfid.size;
					gfid[hash] = tid;
					locus[l].gasize += v + GetNalleles(alleles, v);
				}
			}
		}
		else
		{
			HASH mha = missing_hash[v];
			uint mid = gfid[mha];
			if (!gftab.ContainsKey(mha))
				gftab[mha] = new(gtab[l]++) GENOTYPE(gatab[l], missing_genotype[v]);

			HASH hash = HashGenotype(alleles, v);
			uint gid = gfid[hash];
			if (!gftab.ContainsKey(hash))
				gftab[hash] = new(gtab[l]++) GENOTYPE(gatab[l], alleles, v);

			wt[l].Write(genotype_filter && f_ploidy_b && (f_ploidy_min > v || v > f_ploidy_max) ? mid : gid);
		}
	}
}

/* Create individual from arlequin */
template<typename REAL>
TARGET void IND<REAL>::arlequin(char* t, bool iscount, GENOTYPE** gtab, ushort** gatab, GENO_WRITER* wt)
{
	ReplaceChar(t, '\t', ' ');
	while (*t == ' ') t++;
	name = t;
	char* genstr = t + strlen(t) - 1;
	while (*genstr == ' ' || *genstr == '\r' || *genstr == '\n') *genstr-- = '\0';

	int64 nline = CountChar(t, '\n') + 1;
	if (t[strlen(t) - 1] == '\n')  nline--;

	char* lines[N_MAX_PLOIDY];
	genstr = t;
	for (int i = 0; i < nline; ++i)
	{
		while (*genstr == ' ') genstr++;
		if (i == 0) genstr = StrNextIdx(genstr, ' ', 2) + 1; //name 1 geno
		lines[i] = genstr;
		genstr = StrNextIdx(genstr, '\n', 1) + 1;
	}

	genstr = StrNextIdx(t, ' ', 1);
	char* tname = genstr - 1;
	*genstr++ = '\0';

	if (!iscount)
	{
		while (*tname == ' ') *tname-- = '\0';
		individual_memory->Alloc(tname, (int)strlen(name) + 1);
		//tname = new char[strlen(name) + 1];
		strcpy(tname, name);
		name = tname;
	}
	else name = NULL;

	int v = (int)nline;
	for (int64 l = 0; l < nloc; ++l)
	{
		TABLE<HASH, GENOTYPE*>& gftab = locus[l].gftab;
		TABLE<HASH, uint>& gfid = nvcf_gfid[l];

		if (v > N_MAX_PLOIDY)
			Exit("\nError: ploidy level greater than %d in individual %s at locus %s.\n", N_MAX_PLOIDY, name, locus[l].GetName());

		ushort alleles[N_MAX_PLOIDY] = { 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF };
		bool ismissing = false;

		for (int j = 0; j < v; ++j)
		{
			alleles[j] = (ushort)ReadInteger(lines[j]);
			if (alleles[j] == 0xFFFF)
				ismissing = true;
		}

		if (ismissing) SetFF(alleles, v);
		Sort(alleles, v);//unphase

		if (iscount)
		{
			HASH mha = missing_hash[v];
			if (!gfid.ContainsKey(mha))
			{
				int tid = gfid.size;
				gfid[mha] = tid;
				locus[l].gasize += 0;
			}
			if (!ismissing)
			{
				HASH hash = HashGenotype(alleles, v);
				if (!gfid.ContainsKey(hash))
				{
					int tid = gfid.size;
					gfid[hash] = tid;
					locus[l].gasize += v + GetNalleles(alleles, v);
				}
			}
		}
		else
		{
			HASH mha = missing_hash[v];
			uint mid = gfid[mha];
			if (!gftab.ContainsKey(mha))
				gftab[mha] = new(gtab[l]++) GENOTYPE(gatab[l], missing_genotype[v]);

			HASH hash = HashGenotype(alleles, v);
			uint gid = gfid[hash];
			if (!gftab.ContainsKey(hash))
				gftab[hash] = new(gtab[l]++) GENOTYPE(gatab[l], alleles, v);

			wt[l].Write(genotype_filter && f_ploidy_b && (f_ploidy_min > v || v > f_ploidy_max) ? mid : gid);
		}
	}
}

/* Create individual from structure */
template<typename REAL>
TARGET void IND<REAL>::structure(char* t, bool iscount, GENOTYPE** gtab, ushort** gatab, GENO_WRITER* wt)
{
	name = t;

	int64 nline = CountChar(t, '\n') + 1;
	if (t[strlen(t) - 1] == '\n')  nline--;

	char* lines[N_MAX_PLOIDY];
	char* genstr = t;
	for (int i = 0; i < nline; ++i)
	{
		genstr = StrNextIdx(genstr, ' ', 1 + g_extracol_val) + 1; //skip name, pop
		lines[i] = genstr;
		genstr = StrNextIdx(genstr, '\n', 1) + 1;
	}

	genstr = StrNextIdx(t, ' ', 1);
	char* tname = genstr - 1;
	*genstr++ = '\0';

	if (!iscount)
	{
		while (*tname == ' ') *tname-- = '\0';
		individual_memory->Alloc(tname, (int)strlen(name) + 1);
		//tname = new char[strlen(name) + 1];
		strcpy(tname, name);
		name = tname;
	}
	else name = NULL;

	int v = (int)nline;
	for (int64 l = 0; l < nloc; ++l)
	{
		TABLE<HASH, GENOTYPE*>& gftab = locus[l].gftab;
		TABLE<HASH, uint>& gfid = nvcf_gfid[l];

		if (v > N_MAX_PLOIDY)
			Exit("\nError: ploidy level greater than %d in individual %s at locus %s.\n", N_MAX_PLOIDY, name, locus[l].GetName());

		ushort alleles[N_MAX_PLOIDY] = { 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF };
		bool ismissing = false;

		for (int j = 0; j < v; ++j)
		{
			alleles[j] = (ushort)ReadInteger(lines[j]);
			if (alleles[j] == 0xFFF7)//-9 missing
				ismissing = true;
		}

		if (ismissing) SetFF(alleles, v);
		Sort(alleles, v);//unphase

		if (iscount)
		{
			HASH mha = missing_hash[v];
			if (!gfid.ContainsKey(mha))
			{
				int tid = gfid.size;
				gfid[mha] = tid;
				locus[l].gasize += 0;
			}
			if (!ismissing)
			{
				HASH hash = HashGenotype(alleles, v);
				if (!gfid.ContainsKey(hash))
				{
					int tid = gfid.size;
					gfid[hash] = tid;
					locus[l].gasize += v + GetNalleles(alleles, v);
				}
			}
		}
		else
		{
			HASH mha = missing_hash[v];
			uint mid = gfid[mha];
			if (!gftab.ContainsKey(mha))
				gftab[mha] = new(gtab[l]++) GENOTYPE(gatab[l], missing_genotype[v]);

			HASH hash = HashGenotype(alleles, v);
			uint gid = gfid[hash];
			if (!gftab.ContainsKey(hash))
				gftab[hash] = new(gtab[l]++) GENOTYPE(gatab[l], alleles, v);

			wt[l].Write(genotype_filter && f_ploidy_b && (f_ploidy_min > v || v > f_ploidy_max) ? mid : gid);
		}
	}
}

/* Create individual from polygene */
template<typename REAL>
TARGET void IND<REAL>::polygene(char* t, bool iscount, GENOTYPE** gtab, ushort** gatab, GENO_WRITER* wt)
{
	name = t;
	char* genstr = NULL;

	genstr = StrNextIdx(t, '\t', 1);
	char* tname = genstr - 1;
	*genstr++ = '\0';

	if (!iscount)
	{
		while (*tname == ' ') *tname-- = '\0';
		individual_memory->Alloc(tname, (int)strlen(name) + 1);
		//tname = new char[strlen(name) + 1];
		strcpy(tname, name);
		name = tname;
	}
	else name = NULL;

	genstr = StrNextIdx(genstr, '\t', 1) + 1;
	int ploidy = ReadInteger(genstr);
	if (ploidy <= 0 || ploidy > N_MAX_PLOIDY)
		Exit("\nError: ploidy level greater than %d in individual %s.\n", N_MAX_PLOIDY, name);

	for (int64 l = 0; l < nloc; ++l)
	{
		TABLE<HASH, GENOTYPE*>& gftab = locus[l].gftab;
		TABLE<HASH, uint>& gfid = nvcf_gfid[l];
		if (*genstr == '\t') genstr++;
		ushort alleles[N_MAX_PLOIDY] = { 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF };

		int64 len = StrNextIdx(genstr, '\t', 1) ? StrNextIdx(genstr, '\t', 1) - genstr : strlen(genstr);
		int v = (int)CountChar(genstr, ',', len);

		bool ismissing = false;
		if (len == 0)
			ismissing = true;
		else if (v && v != ploidy - 1)
			Exit("\nError: allelic phenotype in polygene format is not supported.\n");

		v = 0;
		while (*genstr != '\t' && *genstr != '\r' && *genstr != '\n' && *genstr != ' ' && *genstr != '\0')
			alleles[v++] = (ushort)ReadInteger(genstr);

		v = ploidy;
		if (ismissing) SetFF(alleles, v);
		Sort(alleles, v);//unphase

		if (iscount)
		{
			HASH mha = missing_hash[v];
			if (!gfid.ContainsKey(mha))
			{
				int tid = gfid.size;
				gfid[mha] = tid;
				locus[l].gasize += 0;
			}

			if (!ismissing) 
			{
				HASH hash = HashGenotype(alleles, v);
				if (!gfid.ContainsKey(hash))
				{
					int tid = gfid.size;
					gfid[hash] = tid;
					locus[l].gasize += v + GetNalleles(alleles, v);
				}
			}
		}
		else
		{
			HASH mha = missing_hash[v];
			uint mid = gfid[mha];
			if (!gftab.ContainsKey(mha))
				gftab[mha] = new(gtab[l]++) GENOTYPE(gatab[l], missing_genotype[v]);

			HASH hash = HashGenotype(alleles, v);
			uint gid = gfid[hash];
			if (!gftab.ContainsKey(hash))
				gftab[hash] = new(gtab[l]++) GENOTYPE(gatab[l], alleles, v);

			wt[l].Write(genotype_filter && f_ploidy_b && (f_ploidy_min > v || v > f_ploidy_max) ? mid : gid);
		}
	}
}

/* Create individual from polyrelatedness */
template<typename REAL>
TARGET void IND<REAL>::polyrelatedness(char* t, bool iscount, GENOTYPE** gtab, ushort** gatab, GENO_WRITER* wt)
{
	name = t;
	char* genstr = NULL;

	genstr = StrNextIdx(t, '\t', 1);
	char* tname = genstr - 1;
	*genstr++ = '\0';

	if (!iscount)
	{
		while (*tname == ' ') *tname-- = '\0';
		individual_memory->Alloc(tname, (int)strlen(name) + 1);
		//tname = new char[strlen(name) + 1];
		strcpy(tname, name);
		name = tname;
	}
	else name = NULL;

	genstr = StrNextIdx(genstr, '\t', 1) + 1;

	int ndigit = genotype_digit;
	for (int64 l = 0; l < nloc; ++l)
	{
		TABLE<HASH, GENOTYPE*>& gftab = locus[l].gftab;
		TABLE<HASH, uint>& gfid = nvcf_gfid[l];

		while (*genstr == '\t') genstr++;
		ushort alleles[N_MAX_PLOIDY] = { 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF };
		int v = 0;
		bool ismissing = false;
		while (*genstr != '\t' && *genstr != '\r' && *genstr != '\n' && *genstr != ' ' && *genstr != '\0')
		{
			alleles[v] = (ushort)ReadIntegerSpagedi(genstr, ndigit);
			if (alleles[v] == genotype_missing)
				ismissing = true;
			if (alleles[v] == genotype_ambiguous)
				Exit("\nError: do not support ambiguous alleles in individual %s at locus %s\n.", name, locus[l].GetName());
			if (++v > N_MAX_PLOIDY)
				Exit("\nError: ploidy level greater than %d in individual %s at locus %s.\n", N_MAX_PLOIDY, name, locus[l].GetName());
		}

		if (ismissing) SetFF(alleles, v);
		Sort(alleles, v);//unphase

		if (iscount)
		{
			HASH mha = missing_hash[v];
			if (!gfid.ContainsKey(mha))
			{
				int tid = gfid.size;
				gfid[mha] = tid;
				locus[l].gasize += 0;
			}

			if (!ismissing) 
			{
				HASH hash = HashGenotype(alleles, v);
				if (!gfid.ContainsKey(hash))
				{
					int tid = gfid.size;
					gfid[hash] = tid;
					locus[l].gasize += v + GetNalleles(alleles, v);
				}
			}
		}
		else
		{
			HASH mha = missing_hash[v];
			uint mid = gfid[mha];
			if (!gftab.ContainsKey(mha))
				gftab[mha] = new(gtab[l]++) GENOTYPE(gatab[l], missing_genotype[v]);

			HASH hash = HashGenotype(alleles, v);
			uint gid = gfid[hash];
			if (!gftab.ContainsKey(hash))
				gftab[hash] = new(gtab[l]++) GENOTYPE(gatab[l], alleles, v);

			wt[l].Write(genotype_filter && f_ploidy_b && (f_ploidy_min > v || v > f_ploidy_max) ? mid : gid);
		}
	}
}

/* Create individual from genodive */
template<typename REAL>
TARGET void IND<REAL>::genodive(char* t, bool iscount, GENOTYPE** gtab, ushort** gatab, GENO_WRITER* wt)
{
	//2nd column
	name = StrNextChar(StrNextSpace(StrNextChar(t)));

	//if exist extra column
	if (name[0] >= '0' && name[0] <= '9')
		name = StrNextChar(StrNextSpace(name));

	//1st char of genotype
	char* genstr = StrNextChar(StrNextSpace(name));

	if (!iscount)
	{
		char* tname = StrNextSpace(name);
		*tname = '\0';
		individual_memory->Alloc(tname, tname - name + 1);
		//tname = new char[strlen(name) + 1];
		strcpy(tname, name);
		name = tname;
	}
	else name = NULL;

	//set end of name
	genstr[-1] = '\0';

	int ndigit = genotype_digit;
	int maxv = 0;
	for (int64 l = 0; l < nloc; ++l)
	{
		TABLE<HASH, GENOTYPE*>& gftab = locus[l].gftab;
		TABLE<HASH, uint>& gfid = nvcf_gfid[l];

		while (*genstr == '\t') genstr++;
		ushort alleles[N_MAX_PLOIDY] = { 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF };
		int v = 0;
		bool ismissing = false;
		while (*genstr != '\t' && *genstr != '\r' && *genstr != '\n' && *genstr != '\0')
		{
			alleles[v] = (ushort)ReadIntegerSpagedi(genstr, ndigit);
			if (alleles[v] == 0)
				ismissing = true;
			if (++v > N_MAX_PLOIDY)
				Exit("\nError: ploidy level greater than %d in individual %s at locus %s.\n", N_MAX_PLOIDY, name, locus[l].GetName());
		}
		maxv = Max(v, maxv);

		if (!iscount && ismissing)
		{
			v = maxploidy;
			SetFF(alleles, v);
		}
		Sort(alleles, v);//unphase

		if (iscount)
		{
			HASH mha = missing_hash[v];
			if (!gfid.ContainsKey(mha))
			{
				int tid = gfid.size;
				gfid[mha] = tid;
				locus[l].gasize += 0;
			}

			if (!ismissing)
			{
				HASH hash = HashGenotype(alleles, v);
				if (!gfid.ContainsKey(hash))
				{
					int tid = gfid.size;
					gfid[hash] = tid;
					locus[l].gasize += v + GetNalleles(alleles, v);
				}
			}
		}
		else
		{
			HASH mha = missing_hash[v];
			uint mid = gfid[mha];
			if (!gftab.ContainsKey(mha))
				gftab[mha] = new(gtab[l]++) GENOTYPE(gatab[l], missing_genotype[v]);

			HASH hash = HashGenotype(alleles, v);
			uint gid = gfid[hash];
			if (!gftab.ContainsKey(hash))
				gftab[hash] = new(gtab[l]++) GENOTYPE(gatab[l], alleles, v);
			wt[l].Write(genotype_filter && f_ploidy_b && (f_ploidy_min > v || v > f_ploidy_max) ? mid : gid);
		}
	}
	maxploidy = (byte)maxv;
}

/* Create individual from plink */
template<typename REAL>
TARGET void IND<REAL>::plink(char* t, bool iscount, GENOTYPE** gtab, ushort** gatab, GENO_WRITER* wt)
{
	ReplaceChar(t, '\t', ' ');
	int extracol = genotype_extracol;
	name = StrNextIdx(t, ' ', 1) + 1;
	char* genstr = NULL;

	genstr = StrNextIdx(name, ' ', 1);
	char* tname = genstr - 1;
	*genstr++ = '\0';

	if (!iscount)
	{
		while (*tname == ' ') *tname-- = '\0';
		individual_memory->Alloc(tname, (int)strlen(name) + 1);
		strcpy(tname, name);
		name = tname;
	}
	else name = NULL;

	genstr = StrNextIdx(genstr, ' ', 4) + 1;

	int v = 2;
	for (int64 l = 0; l < nloc; ++l)
	{
		TABLE<HASH, GENOTYPE*>& gftab = locus[l].gftab;
		TABLE<HASH, uint>& gfid = nvcf_gfid[l];
		bool ismissing = false;

		ushort alleles[N_MAX_PLOIDY] = { 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF };
		if (*genstr != ',' && *genstr != '\r' && *genstr != '\n' && *genstr != '\0')
		{
			switch (*genstr)
			{
			case 'A': alleles[0] = 1; break;
			case 'T': alleles[0] = 2; break;
			case 'C': alleles[0] = 3; break;
			case 'G': alleles[0] = 4; break;
			default:  alleles[0] = (ushort)ReadInteger(genstr); break;
			}
		}
		while (*genstr == ' ') genstr++;

		if (*genstr != ',' && *genstr != '\r' && *genstr != '\n' && *genstr != '\0')
		{
			switch (*genstr)
			{
			case 'A': alleles[1] = 1; break;
			case 'T': alleles[1] = 2; break;
			case 'C': alleles[1] = 3; break;
			case 'G': alleles[1] = 4; break;
			default:  alleles[1] = (ushort)ReadInteger(genstr); break;
			}
		}
		while (*genstr == ' ') genstr++;

		if (alleles[0] == 0 || alleles[0] == 0xFFFF || alleles[1] == 0 || alleles[1] == 0xFFFF)
		{
			ismissing = true;
			SetFF(alleles, v);
		}
		Sort(alleles, v);//unphase

		if (iscount)
		{
			HASH mha = missing_hash[v];
			if (!gfid.ContainsKey(mha))
			{
				int tid = gfid.size;
				gfid[mha] = tid;
				locus[l].gasize += 0;
			}
			if (!ismissing)
			{
				HASH hash = HashGenotype(alleles, v);
				if (!gfid.ContainsKey(hash))
				{
					int tid = gfid.size;
					gfid[hash] = tid;
					locus[l].gasize += v + GetNalleles(alleles, v);
				}
			}
		}
		else
		{
			HASH mha = missing_hash[v];
			uint mid = gfid[mha];
			if (!gftab.ContainsKey(mha))
				gftab[mha] = new(gtab[l]++) GENOTYPE(gatab[l], missing_genotype[v]);

			HASH hash = HashGenotype(alleles, v);
			uint gid = gfid[hash];
			if (!gftab.ContainsKey(hash))
				gftab[hash] = new(gtab[l]++) GENOTYPE(gatab[l], alleles, v);

			wt[l].Write(genotype_filter && f_ploidy_b && (f_ploidy_min > v || v > f_ploidy_max) ? mid : gid);
		}
	}
}

/* Read and set genotype from bcf input */
template<typename REAL>
TARGET void IND<REAL>::AddBCFGenotype(int64 l, char*& gtstr, char*& gqstr, char*& dpstr, char*& adstr, int vlen, int asize, int gqlen, int dplen, int adlen, uint*& depth, TABLE<HASH, uint>& gfid, GENOTYPE*& gtab, ushort*& gatab, GENO_WRITER& wt)
{
	LOCUS& loc = locus[l];
	TABLE<HASH, GENOTYPE*>& gftab = loc.gftab;
	ushort alleles[N_MAX_PLOIDY] = { 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF };

	//convert hash
	char* gqnext = gqstr ? gqstr + gqlen : NULL, * dpnext = dpstr ? dpstr + dplen : NULL, * adnext = adstr ? adstr + adlen : NULL;
	bool phase = true;
	int v = 0;
	ReadBCFGenoString(gtstr, alleles, phase, v, asize, vlen, l, name);

	//if use phased genotype, exclude unphased genotype as missing
	if (usephase && !phase)
	{
		wt.Write(gfid[missing_hash[v]]);
		gtstr += asize * vlen; gqstr = gqnext; dpstr = dpnext; adstr = adnext;
		return;
	}

	if (!usephase) Sort(alleles, v); //phaseed, do not sort alleles

	/*******************************************************************/

	HASH mha = missing_hash[v], hash = HashGenotype(alleles, v);
	uint mid = gfid[mha], gid = gfid[hash];

	if (!gftab.ContainsKey(mha))
		gftab[mha] = new(gtab++) GENOTYPE(gatab, missing_genotype[v]);

	if (!gftab.ContainsKey(hash))
		gftab[hash] = new(gtab++) GENOTYPE(gatab, alleles, v);

	/*******************************************************************/

	if (genotype_filter)
	{
		if (gid != mid && f_dp_b && dpstr && locus[l].dpid != 0xFFFF)
		{
			uint dp = ReadBinInteger(dpstr, dplen);
			if (dp != -1 && (dp < f_dp_min || dp > f_dp_max))
				gid = mid;
		}

		if (gid != mid && f_gq_b && gqstr && locus[l].gqid != 0xFFFF)
		{
			int gq2 = ReadBinInteger(gqstr, gqlen);
			if (gq2 != -1 && (gq2 < f_gq_min || gq2 > f_gq_max))
				gid = mid;
		}

		if (gid != mid && f_ploidy_b && (v < f_ploidy_min || v > f_ploidy_max))
			gid = mid;
	}

	if (ad == 2 && gid != mid && locus[l].adid != 0xFFFF)
	{
		REAL* fre = pop<REAL>[popid].GetFreq(l);
		for (int j = 0; j < locus[l].k; ++j)
			fre[j] += ReadBinInteger(adstr, adlen);
	}
	else if (ploidyinfer && locus[l].adid != 0xFFFF)
	{
		uint maxdepth = 0;
		for (int j = 0; j < locus[l].k; ++j)
		{
			*depth = ReadBinInteger(adstr, adlen);
			maxdepth = Max(*depth, maxdepth);
			depth++;
		}
		locus[l].maxdepth = maxdepth;
	}

	wt.Write(gid);
	gtstr += asize * vlen; gqstr = gqnext; dpstr = dpnext; adstr = adnext;
}

/* Read and set alleles from vcf input */
template<typename REAL>
TARGET void IND<REAL>::AddVCFGenotype(char*& line, int64 l, uint*& depth, TABLE<HASH, uint>& gfid, GENOTYPE*& gtab, ushort*& gatab, GENO_WRITER& wt)
{
	LOCUS& loc = locus[l];
	TABLE<HASH, GENOTYPE*>& gftab = loc.gftab;
	int fmtc = 0, gtid = loc.gtid;
	char* linebegin = line;
	for (;;)
	{
		fmtc++;
		while (*line != ':' && *line != '\t' && *line != '\n' && *line != '\0') line++;
		if (*line == ':')
			*line++ = '\0';
		else
		{
			*line++ = '\0';
			break;
		}
	}

	if (fmtc == 0)
		Exit("\nError: Format error in individual %s at locus %s_%s.\n", name, loc.GetChrom(), loc.GetName());

	//convert hash
	char* genostr = GetTagValue(linebegin, gtid);
	int v = CountPloidy(genostr);
	if (v > N_MAX_PLOIDY) Exit("\nError: ploidy level greater than %d in individual %s at locus %s_%s.\n", N_MAX_PLOIDY, name, loc.GetChrom(), loc.GetName());
	
	//if use phased genotype, exclude unphased genotype as missing
	if (usephase && ContainsChar(genostr, '/'))
	{
		wt.Write(gfid[missing_hash[v]]);
		return;
	}

	ushort alleles[N_MAX_PLOIDY] = { 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF };
	ReadVCFGenoString(alleles, genostr, v, l, name);

	if (!usephase) Sort(alleles, v); //phased, do not sort alleles

	/*******************************************************************/

	HASH mha = missing_hash[v], hash = HashGenotype(alleles, v);
	uint mid = gfid[mha], gid = gfid[hash];

	if (!gftab.ContainsKey(mha))
		gftab[mha] = new(gtab++) GENOTYPE(gatab, missing_genotype[v]);

	if (!gftab.ContainsKey(hash))
		gftab[hash] = new(gtab++) GENOTYPE(gatab, alleles, v);

	/*******************************************************************/

	if (genotype_filter)
	{
		if (gid != mid && f_dp_b && locus[l].dpid != 0xFFFF && locus[l].dpid < fmtc)
		{
			uint dp = (uint)ReadIntegerKeep(GetTagValue(linebegin, locus[l].dpid));
			if (dp != -1 && (dp < f_dp_min || dp > f_dp_max))
				gid = mid;
		}

		if (gid != mid && f_gq_b && locus[l].gqid != 0xFFFF && locus[l].gqid < fmtc)
		{
			int gq = ReadIntegerKeep(GetTagValue(linebegin, locus[l].gqid));
			if (gq != -1 && (gq < f_gq_min || gq > f_gq_max))
				gid = mid;
		}

		if (gid != mid && f_ploidy_b && (v < f_ploidy_min || v > f_ploidy_max))
			gid = mid;
	}

	if (ad == 2 && gid != mid && locus[l].adid != 0xFFFF && locus[l].adid < fmtc)
	{
		REAL* fre = pop<REAL>[popid].GetFreq(l);
		char* adval = GetTagValue(linebegin, locus[l].adid) - 1;
		for (int j = 0; j < locus[l].k; ++j)
			fre[j] += ReadInteger(++adval);
	}
	else if (ploidyinfer && locus[l].adid != 0xFFFF)
	{
		if (locus[l].adid < fmtc)
		{
			uint maxdepth = 0;
			char* adstr = GetTagValue(linebegin, locus[l].adid) - 1;
			for (int j = 0; j < locus[l].k; ++j)
			{
				*depth = ReadInteger(++adstr);
				maxdepth = Max(*depth, maxdepth);
				depth++;
			}
			locus[l].maxdepth = maxdepth;
		}
		else SetZero(depth, locus[l].k);
	}

	wt.Write(gid);
}

/* Get tag value */
template<typename REAL>
TARGET char* IND<REAL>::GetTagValue(char* re, int tagid)
{
	if (tagid == -1) return 0;
	for (int i = 0; i < tagid; ++i)
		while (*re++);
	return re;
}

/* Count ploidy level from VCF genotype string */
__forceinline TARGET int CountPloidy(char* str)
{
	int count = 1;
	for (int i = 1; str[i + 1]; ++i)
		count += (str[i] == '|') | (str[i] == '/');
	return count;
}

/* Count ploidy level from VCF genotype string */
__forceinline TARGET int CountPloidy(char* str, int len)
{
	int re = 1;
	switch (len)
	{
	case  0: return 0;
	default :
	{
		for (int i = 1; i < len - 1; ++i)
			re += (str[i] == '|') | (str[i] == '/');
		return re;
	}
	case 20: re += (str[18] == '|') | (str[18] == '/');
	case 19: re += (str[17] == '|') | (str[17] == '/');
	case 18: re += (str[16] == '|') | (str[16] == '/');
	case 17: re += (str[15] == '|') | (str[15] == '/');
	case 16: re += (str[14] == '|') | (str[14] == '/');
	case 15: re += (str[13] == '|') | (str[13] == '/');
	case 14: re += (str[12] == '|') | (str[12] == '/');
	case 13: re += (str[11] == '|') | (str[11] == '/');
	case 12: re += (str[10] == '|') | (str[10] == '/');
	case 11: re += (str[ 9] == '|') | (str[ 9] == '/');
	case 10: re += (str[ 8] == '|') | (str[ 8] == '/');
	case  9: re += (str[ 7] == '|') | (str[ 7] == '/');
	case  8: re += (str[ 6] == '|') | (str[ 6] == '/');
	case  7: re += (str[ 5] == '|') | (str[ 5] == '/');
	case  6: re += (str[ 4] == '|') | (str[ 4] == '/');
	case  5: re += (str[ 3] == '|') | (str[ 3] == '/');
	case  4: re += (str[ 2] == '|') | (str[ 2] == '/');
	case  3: re += (str[ 1] == '|') | (str[ 1] == '/');
	case  2: ;
	case  1: ;
	}
	return re;
}
