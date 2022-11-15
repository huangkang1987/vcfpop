/* File Conversion Functions */

#pragma once
#include "vcfpop.h"

template TARGET void ConvertFile<double>();
template TARGET void ConvertFile<float >();
template TARGET void ConvertGenepop<double>(int ntot, bool& isfirst);
template TARGET void ConvertGenepop<float >(int ntot, bool& isfirst);
template TARGET void ConvertSpagedi<double>(int ntot, bool& isfirst);
template TARGET void ConvertSpagedi<float >(int ntot, bool& isfirst);
template TARGET void ConvertCervus<double>(int ntot, bool& isfirst);
template TARGET void ConvertCervus<float >(int ntot, bool& isfirst);
template TARGET void ConvertArlequin<double>(int ntot, bool& isfirst);
template TARGET void ConvertArlequin<float >(int ntot, bool& isfirst);
template TARGET void ConvertStructure<double>(int ntot, bool& isfirst);
template TARGET void ConvertStructure<float >(int ntot, bool& isfirst);
template TARGET void ConvertPolygene<double>(int ntot, bool& isfirst);
template TARGET void ConvertPolygene<float >(int ntot, bool& isfirst);
template TARGET void ConvertPolyRelatedness<double>(int ntot, bool& isfirst);
template TARGET void ConvertPolyRelatedness<float >(int ntot, bool& isfirst);
template TARGET void ConvertGenoDive<double>(int ntot, bool& isfirst);
template TARGET void ConvertGenoDive<float >(int ntot, bool& isfirst);

#define extern 
extern MEMORY* conversion_memory;					//Memory class for conversion_string
extern MEMORY* conversion_memory2;					//Memory class for genotype_string and convert_buf
extern LIST<char*>* conversion_string;				//Genotype string for each genotype for file conversion
extern FILE* convert_file;							//Convert file pointer
extern char** convert_buf;							//Circle buffer for file conversion, NBUF
extern int64 convert_linesize;						//Max length of each line in converted file
#undef extern 

/* Convert into other genotype formats */
template<typename REAL>
TARGET void ConvertFile()
{
	if (!convert) return;
	if (ad) Exit("\nError: file convertion (-convert) is incompatible with allelic depth (-ad) option.\n");

	NBUF = BIG_FILE ? g_nthread_val * 2 : CALC_THREAD_BUFFER;

	EvaluationBegin();

	bool isfirst = true;
	int ntot = 0;
	for (int i = 1; i <= N_CONVERTER; ++i)
		if (convert_format_val[i])
			ntot += nind;

	for (int i = 1; i <= N_CONVERTER; ++i)
	{
		if (convert_format_val[i] == 0) continue;

		//Alloc memory to buffer genotype string
		convert_buf = new char* [NBUF];
		SetZero(convert_buf, NBUF);
		conversion_memory = new MEMORY[g_nthread_val];
		conversion_memory2 = new MEMORY[g_nthread_val];
		conversion_string = (LIST<char*>*)malloc(nloc * sizeof(LIST<char*>));

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int64 l = 0; l < nloc; ++l)
		{
			threadid = omp_get_thread_num();
			new(&conversion_string[l]) LIST<char*>(&conversion_memory[threadid]);
		}

		switch (i)
		{
		case 1: ConvertGenepop<REAL>(ntot, isfirst); break;
		case 2: ConvertSpagedi<REAL>(ntot, isfirst); break;
		case 3: ConvertCervus<REAL>(ntot, isfirst); break;
		case 4: ConvertArlequin<REAL>(ntot, isfirst); break;
		case 5: ConvertStructure<REAL>(ntot, isfirst); break;
		case 6: ConvertPolygene<REAL>(ntot, isfirst); break;
		case 7: ConvertPolyRelatedness<REAL>(ntot, isfirst); break;
		case 8: ConvertGenoDive<REAL>(ntot, isfirst); break;
		}

		free(conversion_string);
		delete[] conversion_memory2;
		delete[] conversion_memory;
		delete[] convert_buf;
	}

	EvaluationEnd("File conversion");
	NBUF = CALC_THREAD_BUFFER * g_nthread_val;
}

/* Write convert genepop genotypes in a guard thread */
THREAD2(ConvertGenepopGuard)
{
	for (int64& ii = progress1 = 0; ii < nind; ++ii, ++PROGRESS_VALUE)
	{
		GUARD_BEGIN

		if (ii == 0 || rinds[ii]->popid != rinds[ii - 1]->popid)
			fprintf(convert_file, "pop\r\n");

		fwrite(convert_buf[ii % NBUF], (int)strlen(convert_buf[ii % NBUF]), 1, convert_file);

		GUARD_END
	}
}

/* Write convert arlequin genotypes in a guard thread */
THREAD2(ConvertArlequinGuard)
{
	for (int64& ii = progress1 = 0; ii < nind; ++ii, ++PROGRESS_VALUE)
	{
		GUARD_BEGIN

		if (ii == 0 || rinds[ii]->popid != rinds[ii - 1]->popid)
		{
			ushort popid = (ushort)rinds[ii]->popid;
			if (ii) fprintf(convert_file, "\t\t}\r\n\r\n");
			fprintf(convert_file, "\t\tSampleName=\"%s\"\r\n\t\tSampleSize=%d\r\n\t\tSampleData={\r\n",
				apops[popid]->name, apops[popid]->nind);
		}

		fwrite(convert_buf[ii % NBUF], (int)strlen(convert_buf[ii % NBUF]), 1, convert_file);

		GUARD_END
	}
}

/* Write convert genotypes in a guard thread */
THREAD2(ConvertGuard)
{
	for (int64& ii = progress1 = 0; ii < nind; ++ii, ++PROGRESS_VALUE)
	{
		GUARD_BEGIN

		fwrite(convert_buf[ii % NBUF], (int)strlen(convert_buf[ii % NBUF]), 1, convert_file);

		GUARD_END
	}
}

/* Convert genotype string */
TARGET void PrepareGenotypeString(int format)
{
#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
	for (int64 l = 0; l < nloc; ++l)
	{
		threadid = omp_get_thread_num();
		GENOTYPE* gtab = GetLoc(l).GetGtab();
		int ngeno = GetLoc(l).ngeno;

		conversion_string[l].Clear();
		for (int gi = 0; gi < ngeno; ++gi)
		{
			char* str = NULL;

			switch (format)
			{
			case 1: str = gtab[gi].GetGenepopStr(); break;
			case 2: str = gtab[gi].GetSpagediStr(); break;
			case 3: str = gtab[gi].GetCervusStr(); break;
			case 4: str = gtab[gi].GetArlequinStr(); break;
			case 5: str = gtab[gi].GetStructureStr(); break;
			case 6: str = gtab[gi].GetPolygeneStr(); break;
			case 7: str = gtab[gi].GetPolyRelatednessStr(); break;
			case 8: str = gtab[gi].GetGenoDiveStr(); break;
			}
			conversion_string[l].Push(str);
		}
	}
}

/* Convert into genepop format */
template<typename REAL>
TARGET void ConvertGenepop(int ntot, bool& isfirst)
{
	char name_buf[NAME_BUF_LEN];
	struct tm* t1;
	time_t tt1;
	time(&tt1);
	t1 = localtime(&tt1);
	char filename[PATH_LEN];
	convert_file = FOpen(filename, "wb", "%s%s", g_output_val.c_str(), ".convert.genepop.txt");
	fprintf(convert_file, "Genepop data file created by vcfpop %s on %04d-%02d-%02d %02d:%02d:%02d\r\n",
		VERSION, t1->tm_year + 1900, t1->tm_mon + 1, t1->tm_mday, t1->tm_hour, t1->tm_min, t1->tm_sec);

	for (int64 l = 0; l < nloc; ++l)
		fprintf(convert_file, "%s\r\n", GetLoc(l).GetNameStr(name_buf));

	convert_linesize = IND_NAME_LEN + 7 * nloc;
	for (int j = 0; j < NBUF; ++j)
		conversion_memory2->Alloc(convert_buf[j], convert_linesize);

	PrepareGenotypeString(1);

	RunThreads(&ConvertGenepopInd<REAL>, &ConvertGenepopGuard<REAL>, NULL, ntot, nind,
		"\nConverting population genetics software format:\n", g_nthread_val, isfirst);

	isfirst = false;
	fclose(convert_file);
}

/* Convert individual genotypes into genepop format in multiple threads */
THREAD2(ConvertGenepopInd)
{
	for (int64 ii = 0; ii < nind; ++ii)
	{
		THREAD_BEGIN

		char* str = convert_buf[ii % NBUF];
		IND<REAL>& ind = *rinds[ii];

		AppendString(str, ind.name);
		AppendString(str, ",");

		for (int64 l = 0; l < nloc; ++l)
			AppendString(str, conversion_string[l][ind.GetGenotypeId(l)]);
		AppendString(str, "\r\n");

		*str++ = '\0';

		THREAD_END
	}
}

/* Convert into spagedi format */
template<typename REAL>
TARGET void ConvertSpagedi(int ntot, bool& isfirst)
{
	char name_buf[NAME_BUF_LEN];
	struct tm* t1;
	time_t tt1;
	time(&tt1);
	t1 = localtime(&tt1);
	char filename[PATH_LEN];
	convert_file = FOpen(filename, "wb", "%s%s", g_output_val.c_str(), ".convert.spagedi.txt");
	fprintf(convert_file, "//spagedi data file created by vcfpop %s on %04d-%02d-%02d %02d:%02d:%02d\r\n",
		VERSION, t1->tm_year + 1900, t1->tm_mon + 1, t1->tm_mday, t1->tm_hour, t1->tm_min, t1->tm_sec);
	fprintf(convert_file, "//#individuals	#categories	#coordinates	#loci	#digits/allele	max ploidy\r\n");
	fprintf(convert_file, "%d	%d	0	%llu	3	%d\r\n-3\r\nind	pop", nind, npop, (uint64)nloc, maxploidy);

	for (int64 l = 0; l < nloc; ++l)
		fprintf(convert_file, "\t%s", GetLoc(l).GetNameStr(name_buf));
	fprintf(convert_file, "\r\n");

	convert_linesize = IND_NAME_LEN + (maxploidy * 3 + 1) * nloc;
	for (int j = 0; j < NBUF; ++j)
		conversion_memory2->Alloc(convert_buf[j], convert_linesize);

	PrepareGenotypeString(2);

	RunThreads(&ConvertSpagediInd<REAL>, &ConvertGuard<REAL>, NULL, ntot, nind,
		"\nConverting population genetics software format:\n", g_nthread_val, isfirst);

	isfirst = false;
	fprintf(convert_file, "END");
	fclose(convert_file);
}

/* Convert individual genotypes into spagedi format in multiple threads */
THREAD2(ConvertSpagediInd)
{
	for (int64 ii = 0; ii < nind; ++ii)
	{
		THREAD_BEGIN

		char* str = convert_buf[ii % NBUF];
		IND<REAL>& ind = *rinds[ii];

		AppendInt(str, ind.popid + 1);
		AppendString(str, "\t");
		AppendString(str, ind.name);

		for (int64 l = 0; l < nloc; ++l)
			AppendString(str, conversion_string[l][ind.GetGenotypeId(l)]);
		AppendString(str, "\r\n");

		*str++ = '\0';

		THREAD_END
	}
}

/* Convert into cervus format */
template<typename REAL>
TARGET void ConvertCervus(int ntot, bool& isfirst)
{
	char filename[PATH_LEN];
	char name_buf[NAME_BUF_LEN];
	convert_file = FOpen(filename, "wb", "%s%s", g_output_val.c_str(), ".convert.cervus.csv");
	fprintf(convert_file, "ind,pop");

	for (int64 l = 0; l < nloc; ++l)
		fprintf(convert_file, ",%sA,%sB", GetLoc(l).GetNameStr(name_buf), GetLoc(l).GetNameStr(name_buf));
	fprintf(convert_file, "\r\n");

	convert_linesize = IND_NAME_LEN + 8 * nloc;
	for (int j = 0; j < NBUF; ++j)
		conversion_memory2->Alloc(convert_buf[j], convert_linesize);

	PrepareGenotypeString(3);

	RunThreads(&ConvertCervusInd<REAL>, &ConvertGuard<REAL>, NULL, ntot, nind,
		"\nConverting population genetics software format:\n", g_nthread_val, isfirst);
	isfirst = false;
	fclose(convert_file);
}

/* Convert individual genotypes into cervus format in multiple threads */
THREAD2(ConvertCervusInd)
{
	for (int64 ii = 0; ii < nind; ++ii)
	{
		THREAD_BEGIN

		char* str = convert_buf[ii % NBUF];
		IND<REAL>& ind = *rinds[ii];

		AppendString(str, ind.name);
		AppendString(str, ",");
		AppendString(str, apops[ind.popid]->name);

		for (int64 l = 0; l < nloc; ++l)
			AppendString(str, conversion_string[l][ind.GetGenotypeId(l)]);
		AppendString(str, "\r\n");

		*str++ = '\0';

		THREAD_END
	}
}

/* Convert into arlequin format */
template<typename REAL>
TARGET void ConvertArlequin(int ntot, bool& isfirst)
{
	char name_buf[NAME_BUF_LEN];
	struct tm* t1;
	time_t tt1;
	time(&tt1);
	t1 = localtime(&tt1);
	char filename[PATH_LEN];
	convert_file = FOpen(filename, "wb", "%s%s", g_output_val.c_str(), ".convert.arlequin.arp");
	fprintf(convert_file, "#arlequin data file created by vcfpop %s on %04d-%02d-%02d %02d:%02d:%02d\r\n",
		VERSION, t1->tm_year + 1900, t1->tm_mon + 1, t1->tm_mday, t1->tm_hour, t1->tm_min, t1->tm_sec);
	fprintf(convert_file, "[Profile]\r\n\tTitle=\"vcfpop\"\r\n\tNbSamples=%d\r\n\tGenotypicData=1\r\n\tGameticPhase=0\r\n\tDataType=MICROSAT\r\n\tLocusSeparator=WHITESPACE\r\n\tMissingData=\"?\"\r\n\r\n# Locus Name", npop);

	for (int64 l = 0; l < nloc; ++l)
		fprintf(convert_file, "\r\n# %s", GetLoc(l).GetNameStr(name_buf));
	fprintf(convert_file, "\r\n\r\n[Data]\r\n\t[[Samples]]\r\n");

	convert_linesize = IND_NAME_LEN + 9 * nloc;
	for (int j = 0; j < NBUF; ++j)
		conversion_memory2->Alloc(convert_buf[j], convert_linesize);

	PrepareGenotypeString(4);

	RunThreads(&ConvertArlequinInd<REAL>, &ConvertArlequinGuard<REAL>, NULL, ntot, nind,
		"\nConverting population genetics software format:\n", g_nthread_val, isfirst);
	isfirst = false;

	fprintf(convert_file, "\t\t}\r\n\r\n");
	fprintf(convert_file, "[[Structure]]\r\n\r\n\t\tStructureName=\"Region\"\r\n\t\tNbGroups=%d\r\n\r\n", nreg[0]);
	for (int i = 0; i < nreg[0]; ++i)
	{
		POP<REAL>* r = lreg >= 0 ? aregs[0][i] : NULL;
		POP<REAL>** pops = r->vpop;
		fprintf(convert_file, "\t\tGroup={\r\n");
		for (int j = 0; j < r->npop; ++j)
			fprintf(convert_file, "\t\t\t\"%s\"\r\n", pops[j]->name);
		fprintf(convert_file, "\t\t}\r\n\r\n");
	}
	fclose(convert_file);
}

/* Convert individual genotypes into arlequin format in multiple threads */
THREAD2(ConvertArlequinInd)
{
	for (int64 ii = 0; ii < nind; ++ii)
	{
		THREAD_BEGIN

		char* str = convert_buf[ii % NBUF];
		IND<REAL>& ind = *rinds[ii];

		AppendString(str, "\t\t\t");
		AppendString(str, ind.name);
		AppendString(str, " 1");

		for (int64 l = 0; l < nloc; ++l)
			AppendString(str, conversion_string[l][ind.GetGenotypeId(l)]);
		AppendString(str, "\r\n\t\t\t");

		for (int64 l = 0; l < nloc; ++l)
			AppendString(str, conversion_string[l][ind.GetGenotypeId(l)] + 5);
		AppendString(str, "\r\n");

		*str++ = '\0';
		THREAD_END
	}
}

/* Convert into structure format */
template<typename REAL>
TARGET void ConvertStructure(int ntot, bool& isfirst)
{
	char name_buf[NAME_BUF_LEN];
	struct tm* t1;
	time_t tt1;
	time(&tt1);
	t1 = localtime(&tt1);
	char filename[PATH_LEN];
	convert_file = FOpen(filename, "wb", "%s%s", g_output_val.c_str(), ".convert.structure.txt");

	for (int64 l = 0; l < nloc; ++l)
		fprintf(convert_file, "%s ", GetLoc(l).GetNameStr(name_buf));
	FSeek(convert_file, -1, SEEK_CUR);
	fprintf(convert_file, "\r\n");

	convert_linesize = IND_NAME_LEN + (maxploidy * 4) * nloc;
	for (int j = 0; j < NBUF; ++j)
		conversion_memory2->Alloc(convert_buf[j], convert_linesize);

	PrepareGenotypeString(5);

	RunThreads(&ConvertStructureInd<REAL>, &ConvertGuard<REAL>, NULL, ntot, nind,
		"\nConverting population genetics software format:\n", g_nthread_val, isfirst);
	isfirst = false;
	fclose(convert_file);
}

/* Convert individual genotypes into structure format in multiple threads */
THREAD2(ConvertStructureInd)
{
	for (int64 ii = 0; ii < nind; ++ii)
	{
		THREAD_BEGIN

		char* str = convert_buf[ii % NBUF];
		IND<REAL>& ind = *rinds[ii];

		if (ind.vmin != ind.vmax)
			Exit("\nError: Cannot convert structure format due to it is a aneuploid.\n", ind.name);

		for (int h = 0; h < ind.vmin; ++h)
		{
			AppendString(str, ind.name);
			AppendString(str, " ");
			AppendString(str, apops[ind.popid]->name);

			for (int64 l = 0; l < nloc; ++l)
			{
				AppendString(str, " ");
				AppendString(str, conversion_string[l][ind.GetGenotypeId(l)] + h * 5);
			}
			AppendString(str, "\r\n");
		}

		*str++ = '\0';

		THREAD_END
	}
}

/* Convert into polygene format */
template<typename REAL>
TARGET void ConvertPolygene(int ntot, bool& isfirst)
{
	char name_buf[NAME_BUF_LEN];
	char filename[PATH_LEN];
	convert_file = FOpen(filename, "wb", "%s%s", g_output_val.c_str(), ".convert.polygene.txt");
	fprintf(convert_file, "ID\tPop\tPloidy");

	for (int64 l = 0; l < nloc; ++l)
		fprintf(convert_file, "\t%s", GetLoc(l).GetNameStr(name_buf));
	fprintf(convert_file, "\r\n");

	convert_linesize = IND_NAME_LEN + maxploidy * 4 * nloc;
	for (int j = 0; j < NBUF; ++j)
		conversion_memory2->Alloc(convert_buf[j], convert_linesize);

	PrepareGenotypeString(6);

	RunThreads(&ConvertPolygeneInd<REAL>, &ConvertGuard<REAL>, NULL, ntot, nind,
		"\nConverting population genetics software format:\n", g_nthread_val, isfirst);
	isfirst = false;
	fclose(convert_file);
}

/* Convert individual genotypes into polygene format in multiple threads */
THREAD2(ConvertPolygeneInd)
{
	for (int64 ii = 0; ii < nind; ++ii)
	{
		THREAD_BEGIN

		char* str = convert_buf[ii % NBUF];
		IND<REAL>& ind = *rinds[ii];

		if (ind.vmin != ind.vmax)
			Exit("\nError: Cannot convert structure format due to it is a aneuploid.\n", ind.name);

		AppendString(str, ind.name);
		AppendString(str, "\t");
		AppendString(str, apops[ind.popid]->name);
		AppendString(str, "\t");
		AppendInt(str, ind.vmin);

		for (int64 l = 0; l < nloc; ++l)
			AppendString(str, conversion_string[l][ind.GetGenotypeId(l)]);
		AppendString(str, "\r\n");

		*str++ = '\0';
		THREAD_END
	}
}

/* Convert into polyrelatedness format */
template<typename REAL>
TARGET void ConvertPolyRelatedness(int ntot, bool& isfirst)
{
	char name_buf[NAME_BUF_LEN];
	struct tm* t1;
	time_t tt1;
	time(&tt1);
	t1 = localtime(&tt1);
	char filename[PATH_LEN];
	convert_file = FOpen(filename, "wb", "%s%s", g_output_val.c_str(), ".convert.polyrelatedness.txt");
	fprintf(convert_file, "//configuration\r\n//#alleledigits(1~4)\t#outputdigits(0~10)\t#missingallele\t#ambiguousallele\t#nthreads(1~64)\r\n3\t8\t000\t999\t2\r\n//genotype\r\nInd\tPop");
	for (int64 l = 0; l < nloc; ++l)
		fprintf(convert_file, "\t%s", GetLoc(l).GetNameStr(name_buf));
	fprintf(convert_file, "\r\n");

	convert_linesize = IND_NAME_LEN + (maxploidy * 3 + 1) * nloc;
	for (int j = 0; j < NBUF; ++j)
		conversion_memory2->Alloc(convert_buf[j], convert_linesize);

	PrepareGenotypeString(7);

	RunThreads(&ConvertPolyRelatednessInd<REAL>, &ConvertGuard<REAL>, NULL, ntot, nind,
		"\nConverting population genetics software format:\n", g_nthread_val, isfirst);
	isfirst = false;

	fprintf(convert_file, "//end of file");
	fclose(convert_file);
}

/* Convert individual genotypes into polyrelatedness format in multiple threads */
THREAD2(ConvertPolyRelatednessInd)
{
	for (int64 ii = 0; ii < nind; ++ii)
	{
		THREAD_BEGIN

		char* str = convert_buf[ii % NBUF];
		IND<REAL>& ind = *rinds[ii];

		AppendString(str, ind.name);
		AppendString(str, "\t");
		AppendString(str, apops[ind.popid]->name);

		for (int64 l = 0; l < nloc; ++l)
			AppendString(str, conversion_string[l][ind.GetGenotypeId(l)]);
		AppendString(str, "\r\n");

		*str++ = '\0';
		THREAD_END
	}
}

/* Convert into genodive format */
template<typename REAL>
TARGET void ConvertGenoDive(int ntot, bool& isfirst)
{
	char name_buf[NAME_BUF_LEN];
	struct tm* t1;
	time_t tt1;
	time(&tt1);
	t1 = localtime(&tt1);
	char filename[PATH_LEN];
	convert_file = FOpen(filename, "wb", "%s%s", g_output_val.c_str(), ".convert.genodive.txt");
	fprintf(convert_file, "//genodive data file created by vcfpop %s on %04d-%02d-%02d %02d:%02d:%02d\r\n",
		VERSION, t1->tm_year + 1900, t1->tm_mon + 1, t1->tm_mday, t1->tm_hour, t1->tm_min, t1->tm_sec);
	fprintf(convert_file, "%d	%d	%d	%u	3\r\n", nind, npop, nloc, maxploidy);

	for (int i = 0; i < npop; ++i)
	{
		POP<REAL>* tp = apops[i];
		fprintf(convert_file, "%s", tp->name);
		for (int rl = 0; rl <= lreg; ++rl)
		{
			tp = aregs[rl][tp->rid];
			fprintf(convert_file, "\t%s", tp->name);
		}
		fprintf(convert_file, "\r\n");
	}

	fprintf(convert_file, "Population	Individual");
	for (int64 l = 0; l < nloc; ++l)
		fprintf(convert_file, "\t%s", GetLoc(l).GetNameStr(name_buf));
	fprintf(convert_file, "\r\n");

	convert_linesize = IND_NAME_LEN + (maxploidy * 3 + 1) * nloc;
	for (int j = 0; j < NBUF; ++j)
		conversion_memory2->Alloc(convert_buf[j], convert_linesize);

	PrepareGenotypeString(8);

	RunThreads(&ConvertSpagediInd<REAL>, &ConvertGuard<REAL>, NULL, ntot, nind,
		"\nConverting population genetics software format:\n", g_nthread_val, isfirst);

	isfirst = false;
	fprintf(convert_file, "END");
	fclose(convert_file);
}

/* Convert individual genotypes into genodive format in multiple threads */
THREAD2(ConvertGenoDiveInd)
{
	for (int64 ii = 0; ii < nind; ++ii)
	{
		THREAD_BEGIN

		char* str = convert_buf[ii % NBUF];
		IND<REAL>& ind = *rinds[ii];

		AppendString(str, ind.name);
		AppendString(str, "\t");
		AppendString(str, apops[ind.popid]->name);

		for (int64 l = 0; l < nloc; ++l)
			AppendString(str, conversion_string[l][ind.GetGenotypeId(l)]);
		AppendString(str, "\r\n");

		*str++ = '\0';

		THREAD_END
	}
}
