/* File Conversion Functions */

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
template TARGET void ConvertPlink<double>(int ntot, bool& isfirst);
template TARGET void ConvertPlink<float >(int ntot, bool& isfirst);

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
		case 9: ConvertPlink<REAL>(ntot, isfirst); break;
		}

		free(conversion_string);
		DEL(conversion_memory2);
		DEL(conversion_memory);
		DEL(convert_buf);
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

		if (ii == 0 || rinds<REAL>[ii]->popid != rinds<REAL>[ii - 1]->popid)
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

		if (ii == 0 || rinds<REAL>[ii]->popid != rinds<REAL>[ii - 1]->popid)
		{
			ushort popid = (ushort)rinds<REAL>[ii]->popid;
			if (ii) fprintf(convert_file, "\t\t}\r\n\r\n");

			int noutput = apops<REAL>[popid]->nind;
			if (convert_mode_val >= 4)
			{
				noutput = 0;
				for (int j = 0; j < apops<REAL>[popid]->nind; ++j)
				{
					IND<REAL>* ind = apops<REAL>[popid]->inds[j];
					noutput += ind->vmax <= 2 ? 1 : ind->vmax / 2;
				}
			}
			fprintf(convert_file, "\t\tSampleName=\"%s\"\r\n\t\tSampleSize=%d\r\n\t\tSampleData={\r\n",
				apops<REAL>[popid]->name, noutput);
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
			case 9: str = gtab[gi].GetPlinkStr(); break;
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

	convert_linesize = convert_mode_val <= 3 ? 
		IND_NAME_LEN + (2 * MAX_ALLELE_BYTE + 1) * nloc :
		(IND_NAME_LEN + (2 * MAX_ALLELE_BYTE + 1) * nloc) * (maxploidy / 2);

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
		IND<REAL>& ind = *rinds<REAL>[ii];

		switch (convert_mode_val)
		{
		default:
		case 1:
		case 2://disable,truncate
		{
			AppendString(str, ind.name);
			AppendString(str, ",");

			for (int64 l = 0; l < nloc; ++l)
				AppendString(str, conversion_string[l][ind.GetGenotypeId(l)]);

			AppendString(str, "\r\n");
			break;
		}
		case 3://choose
		{
			RNG<double> rng(g_seed_val + ii, RNG_SALT_CONVERT);
			AppendString(str, ind.name);
			AppendString(str, ",");

			for (int64 l = 0; l < nloc; ++l)
			{
				GENOTYPE& gt = ind.GetGenotype(l);
				ushort* als = gt.GetAlleleArray();
				int ploidy = gt.Ploidy(), nalleles = gt.Nalleles();
				if (nalleles == 0 || ploidy == 1) 
					AppendString(str, " 000000");
				else
				{
					int a1 = rng.Next(ploidy), a2 = rng.NextAvoid(ploidy, a1);
					sprintf(str, " %03d%03d", als[a1] + 1, als[a2] + 1);
					while (*str) str++; 
				}
			}

			AppendString(str, "\r\n");
			break;
		}
		case 4://split
		{
			int vmax = ind.vmax, nsplit = vmax / 2, nsplit2 = nsplit * 2, permals[N_MAX_PLOIDY];
			ushort* permals16 = (ushort*)permals;

			for (int j = 0; j < nsplit2; j += 2)
			{
				AppendString(str, ind.name);
				AppendString(str, "A,"); str[-2] += j / 2;

				for (int64 l = 0; l < nloc; ++l)
				{
					GENOTYPE& gt = ind.GetGenotype(l);
					ushort* als = gt.GetAlleleArray();
					int ploidy = gt.Ploidy(), nalleles = gt.Nalleles();
					if (nalleles == 0 || ploidy == 1 || j + 1 >= ploidy) 
						AppendString(str, " 000000");
					else
					{
						sprintf(str, " %03d%03d", als[j] + 1, als[j + 1] + 1);
						while (*str) str++;
					}
				}

				AppendString(str, "\r\n");
			}
			break;
		}
		case 5://shuffle
		{
			int vmax = ind.vmax, nsplit = vmax / 2, nsplit2 = nsplit * 2, permals[N_MAX_PLOIDY];
			ushort* permals16 = (ushort*)permals;

			for (int j = 0; j < nsplit2; j += 2)
			{
				RNG<double> rng(g_seed_val + ii, RNG_SALT_CONVERT);
				AppendString(str, ind.name);
				AppendString(str, "A,"); str[-2] += j / 2;

				for (int64 l = 0; l < nloc; ++l)
				{
					GENOTYPE& gt = ind.GetGenotype(l);
					ushort* als = gt.GetAlleleArray();
					int ploidy = gt.Ploidy(), nalleles = gt.Nalleles();

					if (nalleles == 0 || ploidy == 1) 
						AppendString(str, " 000000");
					else
					{
						SetVal(permals16, als, ploidy);
						rng.Permute(permals16, ploidy);
						SetFF(permals16 + ploidy - ploidy % 2, nsplit2 - ploidy + ploidy % 2);
						rng.Permute(permals, nsplit);
						sprintf(str, " %03d%03d", (short)permals16[j] + 1, (short)permals16[j + 1] + 1);
						while (*str) str++;
					}
				}

				AppendString(str, "\r\n");
			}
			break;
		}
		}

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

	int noutput = nind;
	if (convert_mode_val >= 4)
	{
		noutput = 0;
		for (int j = 0; j < nind; ++j)
		{
			IND<REAL>* ind = ainds<REAL>[j];
			noutput += ind->vmax <= 2 ? 1 : ind->vmax / 2;
		}
	}
	
	fprintf(convert_file, "%d	%d	0	%llu	3	%d\r\n-3\r\nind	pop", noutput, npop, (uint64)nloc, maxploidy);

	for (int64 l = 0; l < nloc; ++l)
		fprintf(convert_file, "\t%s", GetLoc(l).GetNameStr(name_buf));
	fprintf(convert_file, "\r\n");

	convert_linesize =
		convert_mode_val == 1 ? IND_NAME_LEN + (maxploidy * MAX_ALLELE_BYTE + 1) * nloc :
		(convert_mode_val <= 3 ? IND_NAME_LEN + (2 * MAX_ALLELE_BYTE + 1) * nloc :
		(IND_NAME_LEN + (2 * MAX_ALLELE_BYTE + 1) * nloc) * (maxploidy / 2));

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
		IND<REAL>& ind = *rinds<REAL>[ii];

		switch (convert_mode_val)
		{
		default:
		case 1:
		case 2://disable,truncate
		{
			AppendString(str, ind.name);
			AppendString(str, "\t");
			AppendString(str, apops<REAL>[ind.popid]->name);

			for (int64 l = 0; l < nloc; ++l)
				AppendString(str, conversion_string[l][ind.GetGenotypeId(l)]);

			AppendString(str, "\r\n");
			break;
		}
		case 3://choose
		{
			RNG<double> rng(g_seed_val + ii, RNG_SALT_CONVERT);
			AppendString(str, ind.name);
			AppendString(str, "\t");
			AppendString(str, apops<REAL>[ind.popid]->name);

			for (int64 l = 0; l < nloc; ++l)
			{
				GENOTYPE& gt = ind.GetGenotype(l);
				ushort* als = gt.GetAlleleArray();
				int ploidy = gt.Ploidy(), nalleles = gt.Nalleles();
				if (nalleles == 0 || ploidy == 1)
					AppendString(str, "\t000000");
				else
				{
					int a1 = rng.Next(ploidy), a2 = rng.NextAvoid(ploidy, a1);
					sprintf(str, "\t%03d%03d", als[a1] + 1, als[a2] + 1);
					while (*str) str++;
				}
			}

			AppendString(str, "\r\n");
			break;
		}
		case 4://split
		{
			int vmax = ind.vmax, nsplit = vmax / 2, nsplit2 = nsplit * 2, permals[N_MAX_PLOIDY];
			ushort* permals16 = (ushort*)permals;

			for (int j = 0; j < nsplit2; j += 2)
			{
				AppendString(str, ind.name);
				AppendString(str, "A\t"); str[-2] += j / 2;
				AppendString(str, apops<REAL>[ind.popid]->name);

				for (int64 l = 0; l < nloc; ++l)
				{
					GENOTYPE& gt = ind.GetGenotype(l);
					ushort* als = gt.GetAlleleArray();
					int ploidy = gt.Ploidy(), nalleles = gt.Nalleles();
					if (nalleles == 0 || ploidy == 1 || j + 1 >= ploidy)
						AppendString(str, "\t000000");
					else
					{
						sprintf(str, "\t%03d%03d", als[j] + 1, als[j + 1] + 1);
						while (*str) str++;
					}
				}

				AppendString(str, "\r\n");
			}
			break;
		}
		case 5://shuffle
		{
			int vmax = ind.vmax, nsplit = vmax / 2, nsplit2 = nsplit * 2, permals[N_MAX_PLOIDY];
			ushort* permals16 = (ushort*)permals;

			for (int j = 0; j < nsplit2; j += 2)
			{
				RNG<double> rng(g_seed_val + ii, RNG_SALT_CONVERT);
				AppendString(str, ind.name);
				AppendString(str, "A\t"); str[-2] += j / 2;
				AppendString(str, apops<REAL>[ind.popid]->name);

				for (int64 l = 0; l < nloc; ++l)
				{
					GENOTYPE& gt = ind.GetGenotype(l);
					ushort* als = gt.GetAlleleArray();
					int ploidy = gt.Ploidy(), nalleles = gt.Nalleles();

					if (nalleles == 0 || ploidy == 1)
						AppendString(str, "\t000000");
					else
					{
						SetVal(permals16, als, ploidy);
						rng.Permute(permals16, ploidy);
						SetFF(permals16 + ploidy - ploidy % 2, nsplit2 - ploidy + ploidy % 2);
						rng.Permute(permals, nsplit);
						if (permals16[j] == 0xFFFF || permals16[j + 1] == 0xFFFF)
							AppendString(str, "\t000000");
						else
						{
							sprintf(str, "\t%03d%03d", (short)permals16[j] + 1, (short)permals16[j + 1] + 1);
							while (*str) str++;
						}
					}
				}

				AppendString(str, "\r\n");
			}
			break;
		}
		}

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

	convert_linesize = convert_mode_val <= 3 ?
		IND_NAME_LEN + (2 * MAX_ALLELE_BYTE + 2) * nloc :
		(IND_NAME_LEN + (2 * MAX_ALLELE_BYTE + 2) * nloc) * (maxploidy / 2);

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
		IND<REAL>& ind = *rinds<REAL>[ii];

		switch (convert_mode_val)
		{
		default:
		case 1:
		case 2://disable,truncate
		{
			AppendString(str, ind.name);
			AppendString(str, ",");
			AppendString(str, apops<REAL>[ind.popid]->name);

			for (int64 l = 0; l < nloc; ++l)
				AppendString(str, conversion_string[l][ind.GetGenotypeId(l)]);

			AppendString(str, "\r\n");
			break;
		}
		case 3://choose
		{
			RNG<double> rng(g_seed_val + ii, RNG_SALT_CONVERT);
			AppendString(str, ind.name);
			AppendString(str, ",");
			AppendString(str, apops<REAL>[ind.popid]->name);

			for (int64 l = 0; l < nloc; ++l)
			{
				GENOTYPE& gt = ind.GetGenotype(l);
				ushort* als = gt.GetAlleleArray();
				int ploidy = gt.Ploidy(), nalleles = gt.Nalleles();
				if (nalleles == 0 || ploidy == 1) 
					AppendString(str, ",,");
				else
				{
					int a1 = rng.Next(ploidy), a2 = rng.NextAvoid(ploidy, a1);
					sprintf(str, ",%d,%d", als[a1] + 1, als[a2] + 1);
					while (*str) str++;
				}
			}

			AppendString(str, "\r\n");
			break;
		}
		case 4://split
		{
			int vmax = ind.vmax, nsplit = vmax / 2, nsplit2 = nsplit * 2, permals[N_MAX_PLOIDY];
			ushort* permals16 = (ushort*)permals;

			for (int j = 0; j < nsplit2; j += 2)
			{
				AppendString(str, ind.name);
				AppendString(str, "A,"); str[-2] += j / 2;
				AppendString(str, apops<REAL>[ind.popid]->name);

				for (int64 l = 0; l < nloc; ++l)
				{
					GENOTYPE& gt = ind.GetGenotype(l);
					ushort* als = gt.GetAlleleArray();
					int ploidy = gt.Ploidy(), nalleles = gt.Nalleles();
					if (nalleles == 0 || ploidy == 1 || j + 1 >= ploidy) 
						AppendString(str, ",,");
					else
					{
						sprintf(str, ",%d,%d", als[j] + 1, als[j + 1] + 1);
						while (*str) str++;
					}
				}

				AppendString(str, "\r\n");
			}
			break;
		}
		case 5://shuffle
		{
			int vmax = ind.vmax, nsplit = vmax / 2, nsplit2 = nsplit * 2, permals[N_MAX_PLOIDY];
			ushort* permals16 = (ushort*)permals;

			for (int j = 0; j < nsplit2; j += 2)
			{
				RNG<double> rng(g_seed_val + ii, RNG_SALT_CONVERT);
				AppendString(str, ind.name);
				AppendString(str, "A,"); str[-2] += j / 2;
				AppendString(str, apops<REAL>[ind.popid]->name);

				for (int64 l = 0; l < nloc; ++l)
				{
					GENOTYPE& gt = ind.GetGenotype(l);
					ushort* als = gt.GetAlleleArray();
					int ploidy = gt.Ploidy(), nalleles = gt.Nalleles();

					if (nalleles == 0 || ploidy == 1) 
						AppendString(str, ",,");
					else
					{
						SetVal(permals16, als, ploidy);
						rng.Permute(permals16, ploidy);
						SetFF(permals16 + ploidy - ploidy % 2, nsplit2 - ploidy + ploidy % 2);
						rng.Permute(permals, nsplit);
						if (permals16[j] == 0xFFFF || permals16[j + 1] == 0xFFFF)
							AppendString(str, ",,");
						else
						{
							sprintf(str, ",%d,%d", (short)permals16[j] + 1, (short)permals16[j + 1] + 1);
							while (*str) str++;
						}
					}
				}

				AppendString(str, "\r\n");
			}
			break;
		}
		}

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

	convert_linesize = convert_mode_val <= 3 ?
		IND_NAME_LEN + 12 + (2 * MAX_ALLELE_BYTE + 2) * nloc :
		(IND_NAME_LEN + 12 + (2 * MAX_ALLELE_BYTE + 2) * nloc) * (maxploidy / 2);

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
		POP<REAL>* r = lreg >= 0 ? aregs<REAL>[0][i] : NULL;
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
		IND<REAL>& ind = *rinds<REAL>[ii];

		switch (convert_mode_val)
		{
		default:
		case 1:
		case 2://disable,truncate
		{
			AppendString(str, "\t\t\t");
			AppendString(str, ind.name);
			AppendString(str, " 1");

			for (int o = 0; o < 2; ++o)
			{
				for (int64 l = 0; l < nloc; ++l)
					AppendString(str, conversion_string[l][ind.GetGenotypeId(l)] + o * (MAX_ALLELE_BYTE + 2));

				AppendString(str, o == 0 ? "\r\n\t\t\t       " : "\r\n");
			}

			break;
		}
		case 3://choose
		{
			AppendString(str, "\t\t\t");
			AppendString(str, ind.name);
			AppendString(str, " 1");

			for (int o = 0; o < 2; ++o)
			{
				RNG<double> rng(g_seed_val + ii, RNG_SALT_CONVERT);
				for (int64 l = 0; l < nloc; ++l)
				{
					GENOTYPE& gt = ind.GetGenotype(l);
					ushort* als = gt.GetAlleleArray();
					int ploidy = gt.Ploidy(), nalleles = gt.Nalleles();
					if (nalleles == 0 || ploidy == 1)
						AppendString(str, " ?");
					else
					{
						int a[2];
						a[0] = rng.Next(ploidy);
						a[1] = rng.NextAvoid(ploidy, a[0]);
						sprintf(str, " %d", als[a[o]] + 1);
						while (*str) str++;
					}
				}

				AppendString(str, o == 0 ? "\r\n\t\t\t       " : "\r\n");
			}
			break;
		}
		case 4://split
		{
			int vmax = ind.vmax, nsplit = vmax / 2, nsplit2 = nsplit * 2, permals[N_MAX_PLOIDY];
			ushort* permals16 = (ushort*)permals;

			for (int j = 0; j < nsplit2; j += 2)
			{
				AppendString(str, "\t\t\t");
				AppendString(str, ind.name);
				AppendString(str, "A"); str[-1] += j / 2;
				AppendString(str, " 1");

				for (int o = 0; o < 2; ++o)
				{
					for (int64 l = 0; l < nloc; ++l)
					{
						GENOTYPE& gt = ind.GetGenotype(l);
						ushort* als = gt.GetAlleleArray();
						int ploidy = gt.Ploidy(), nalleles = gt.Nalleles();
						if (nalleles == 0 || ploidy == 1 || j + 1 >= ploidy)
							AppendString(str, " ?");
						else
						{
							sprintf(str, " %d", als[j + o] + 1);
							while (*str) str++;
						}
					}

					AppendString(str, o == 0 ? "\r\n\t\t\t       " : "\r\n");
				}
			}
			break;
		}
		case 5://shuffle
		{
			int vmax = ind.vmax, nsplit = vmax / 2, nsplit2 = nsplit * 2, permals[N_MAX_PLOIDY];
			ushort* permals16 = (ushort*)permals;

			for (int j = 0; j < nsplit2; j += 2)
			{
				AppendString(str, "\t\t\t");
				AppendString(str, ind.name);
				AppendString(str, "A"); str[-1] += j / 2;
				AppendString(str, " 1");

				for (int o = 0; o < 2; ++o)
				{
					RNG<double> rng(g_seed_val + ii, RNG_SALT_CONVERT);
					for (int64 l = 0; l < nloc; ++l)
					{
						GENOTYPE& gt = ind.GetGenotype(l);
						ushort* als = gt.GetAlleleArray();
						int ploidy = gt.Ploidy(), nalleles = gt.Nalleles();

						if (nalleles == 0 || ploidy == 1)
							AppendString(str, " ?");
						else
						{
							SetVal(permals16, als, ploidy);
							rng.Permute(permals16, ploidy);
							SetFF(permals16 + ploidy - ploidy % 2, nsplit2 - ploidy + ploidy % 2);
							rng.Permute(permals, nsplit);
							if (permals16[j] == 0xFFFF || permals16[j + 1] == 0xFFFF)
								AppendString(str, " ?");
							else
							{
								sprintf(str, " %d", (short)permals16[j + o] + 1);
								while (*str) str++;
							}
						}
					}

					AppendString(str, o == 0 ? "\r\n\t\t\t       " : "\r\n");
				}
			}
			break;
		}
		}

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

	convert_linesize = 
		convert_mode_val == 1 ? (IND_NAME_LEN + (MAX_ALLELE_BYTE + 1) * nloc) * maxploidy :
		(convert_mode_val <= 3 ? (IND_NAME_LEN + (MAX_ALLELE_BYTE + 1) * nloc) * 2 :
		(IND_NAME_LEN + (MAX_ALLELE_BYTE + 1) * nloc) * 2 * (maxploidy / 2));

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
		IND<REAL>& ind = *rinds<REAL>[ii];

		switch (convert_mode_val)
		{
		default:
		case 1://disable
			if (ind.vmin != ind.vmax)
				Exit("\nError: Cannot convert structure format due to it is a aneuploid.\n", ind.name);
		case 2://truncate
		{	
			for (int h = 0; h < (convert_mode_val == 1 ? ind.vmin : 2); ++h)
			{
				AppendString(str, ind.name);
				AppendString(str, " ");
				AppendString(str, apops<REAL>[ind.popid]->name);

				for (int64 l = 0; l < nloc; ++l)
				{
					AppendString(str, " ");
					AppendString(str, conversion_string[l][ind.GetGenotypeId(l)] + h * (MAX_ALLELE_BYTE + 2));
				}

				AppendString(str, "\r\n");
			}
			break;
		}
		case 3://choose
		{
			for (int h = 0; h < 2; ++h)
			{
				AppendString(str, ind.name);
				AppendString(str, " ");
				AppendString(str, apops<REAL>[ind.popid]->name);
				RNG<double> rng(g_seed_val + ii, RNG_SALT_CONVERT);

				for (int64 l = 0; l < nloc; ++l)
				{
					GENOTYPE& gt = ind.GetGenotype(l);
					ushort* als = gt.GetAlleleArray();
					int ploidy = gt.Ploidy(), nalleles = gt.Nalleles();
					if (nalleles == 0 || ploidy == 1)
						AppendString(str, " -9");
					else
					{
						int a[2];
						a[0] = rng.Next(ploidy);
						a[1] = rng.NextAvoid(ploidy, a[0]);
						sprintf(str, " %d", als[a[h]] + 1);
						while (*str) str++;
					}
				}

				AppendString(str, "\r\n");
			}
			break;
		}
		case 4://split
		{
			int vmax = ind.vmax, nsplit = vmax / 2, nsplit2 = nsplit * 2, permals[N_MAX_PLOIDY];
			ushort* permals16 = (ushort*)permals;

			for (int j = 0; j < nsplit2; j += 2)
			{
				for (int h = 0; h < 2; ++h)
				{
					AppendString(str, ind.name);
					AppendString(str, "A "); str[-2] += j / 2;
					AppendString(str, apops<REAL>[ind.popid]->name);

					for (int64 l = 0; l < nloc; ++l)
					{
						GENOTYPE& gt = ind.GetGenotype(l);
						ushort* als = gt.GetAlleleArray();
						int ploidy = gt.Ploidy(), nalleles = gt.Nalleles();
						if (nalleles == 0 || ploidy == 1 || j + 1 >= ploidy)
							AppendString(str, " -9");
						else
						{
							sprintf(str, " %d", als[j + h] + 1);
							while (*str) str++;
						}
					}

					AppendString(str, "\r\n");
				}
			}
			break;
		}
		case 5://shuffle
		{
			int vmax = ind.vmax, nsplit = vmax / 2, nsplit2 = nsplit * 2, permals[N_MAX_PLOIDY];
			ushort* permals16 = (ushort*)permals;

			for (int j = 0; j < nsplit2; j += 2)
			{
				for (int h = 0; h < 2; ++h)
				{
					RNG<double> rng(g_seed_val + ii, RNG_SALT_CONVERT);
					AppendString(str, ind.name);
					AppendString(str, "A "); str[-2] += j / 2;
					AppendString(str, apops<REAL>[ind.popid]->name);

					for (int64 l = 0; l < nloc; ++l)
					{
						GENOTYPE& gt = ind.GetGenotype(l);
						ushort* als = gt.GetAlleleArray();
						int ploidy = gt.Ploidy(), nalleles = gt.Nalleles();

						if (nalleles == 0 || ploidy == 1)
							AppendString(str, " -9");
						else
						{
							SetVal(permals16, als, ploidy);
							rng.Permute(permals16, ploidy);
							SetFF(permals16 + ploidy - ploidy % 2, nsplit2 - ploidy + ploidy % 2);
							rng.Permute(permals, nsplit);

							if (permals16[j] == 0xFFFF || permals16[j + 1] == 0xFFFF)
								AppendString(str, " -9");
							else
							{
								sprintf(str, " %d", (short)permals16[j + h] + 1);
								while (*str) str++;
							}
						}
					}

					AppendString(str, "\r\n");
				}
			}
			break;
		}
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

	convert_linesize =
		convert_mode_val == 1 ? IND_NAME_LEN + maxploidy * (MAX_ALLELE_BYTE + 1) * nloc :
		(convert_mode_val <= 3 ? IND_NAME_LEN + 2 * (MAX_ALLELE_BYTE + 1) * nloc :
		(IND_NAME_LEN + 2 * (MAX_ALLELE_BYTE + 1) * nloc) * (maxploidy / 2));

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
		IND<REAL>& ind = *rinds<REAL>[ii];

		switch (convert_mode_val)
		{
		default:
		case 1://disable,
			if (ind.vmin != ind.vmax)
				Exit("\nError: Cannot convert structure format due to it is a aneuploid.\n", ind.name);
		case 2://truncate
		{
			AppendString(str, ind.name);
			AppendString(str, "\t");
			AppendString(str, apops<REAL>[ind.popid]->name);
			AppendString(str, "\t");
			AppendInt(str, convert_mode_val == 1 ? ind.vmin : 2);

			for (int64 l = 0; l < nloc; ++l)
				AppendString(str, conversion_string[l][ind.GetGenotypeId(l)]);

			AppendString(str, "\r\n");
			break;
		}
		case 3://choose
		{
			RNG<double> rng(g_seed_val + ii, RNG_SALT_CONVERT);
			AppendString(str, ind.name);
			AppendString(str, "\t");
			AppendString(str, apops<REAL>[ind.popid]->name);
			AppendString(str, "\t2");

			for (int64 l = 0; l < nloc; ++l)
			{
				GENOTYPE& gt = ind.GetGenotype(l);
				ushort* als = gt.GetAlleleArray();
				int ploidy = gt.Ploidy(), nalleles = gt.Nalleles();
				if (nalleles == 0 || ploidy == 1)
					AppendString(str, "\t");
				else
				{
					int a1 = rng.Next(ploidy), a2 = rng.NextAvoid(ploidy, a1);
					sprintf(str, "\t%d,%d", als[a1] + 1, als[a2] + 1);
					while (*str) str++;
				}
			}

			AppendString(str, "\r\n");
			break;
		}
		case 4://split
		{
			int vmax = ind.vmax, nsplit = vmax / 2, nsplit2 = nsplit * 2, permals[N_MAX_PLOIDY];
			ushort* permals16 = (ushort*)permals;

			for (int j = 0; j < nsplit2; j += 2)
			{
				AppendString(str, ind.name);
				AppendString(str, "A\t"); str[-2] += j / 2;
				AppendString(str, apops<REAL>[ind.popid]->name);
				AppendString(str, "\t2");

				for (int64 l = 0; l < nloc; ++l)
				{
					GENOTYPE& gt = ind.GetGenotype(l);
					ushort* als = gt.GetAlleleArray();
					int ploidy = gt.Ploidy(), nalleles = gt.Nalleles();
					if (nalleles == 0 || ploidy == 1 || j + 1 >= ploidy)
						AppendString(str, "\t");
					else
					{
						sprintf(str, "\t%d,%d", als[j] + 1, als[j + 1] + 1);
						while (*str) str++;
					}
				}

				AppendString(str, "\r\n");
			}
			break;
		}
		case 5://shuffle
		{
			int vmax = ind.vmax, nsplit = vmax / 2, nsplit2 = nsplit * 2, permals[N_MAX_PLOIDY];
			ushort* permals16 = (ushort*)permals;

			for (int j = 0; j < nsplit2; j += 2)
			{
				RNG<double> rng(g_seed_val + ii, RNG_SALT_CONVERT);
				AppendString(str, ind.name);
				AppendString(str, "A\t"); str[-2] += j / 2;
				AppendString(str, apops<REAL>[ind.popid]->name);
				AppendString(str, "\t2");

				for (int64 l = 0; l < nloc; ++l)
				{
					GENOTYPE& gt = ind.GetGenotype(l);
					ushort* als = gt.GetAlleleArray();
					int ploidy = gt.Ploidy(), nalleles = gt.Nalleles();

					if (nalleles == 0 || ploidy == 1)
						AppendString(str, "\t");
					else
					{
						SetVal(permals16, als, ploidy);
						rng.Permute(permals16, ploidy);
						SetFF(permals16 + ploidy - ploidy % 2, nsplit2 - ploidy + ploidy % 2);
						rng.Permute(permals, nsplit);
						if (permals16[j] == 0xFFFF || permals16[j + 1] == 0xFFFF)
							AppendString(str, "\t");
						else
						{
							sprintf(str, "\t%d,%d", (short)permals16[j] + 1, (short)permals16[j + 1] + 1);
							while (*str) str++;
						}
					}
				}

				AppendString(str, "\r\n");
			}
			break;
		}
		}

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

	convert_linesize =
		convert_mode_val == 1 ? IND_NAME_LEN + (MAX_ALLELE_BYTE * maxploidy + 1) * nloc :
		(convert_mode_val <= 3 ? IND_NAME_LEN + (MAX_ALLELE_BYTE * 2 + 1) * nloc :
		(IND_NAME_LEN + (MAX_ALLELE_BYTE * 2 + 1) * nloc) * (maxploidy / 2));

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
		IND<REAL>& ind = *rinds<REAL>[ii];

		switch (convert_mode_val)
		{
		default:
		case 1:
		case 2://disable,truncate
		{
			AppendString(str, ind.name);
			AppendString(str, "\t");
			AppendString(str, apops<REAL>[ind.popid]->name);

			for (int64 l = 0; l < nloc; ++l)
				AppendString(str, conversion_string[l][ind.GetGenotypeId(l)]);

			AppendString(str, "\r\n");
			break;
		}
		case 3://choose
		{
			RNG<double> rng(g_seed_val + ii, RNG_SALT_CONVERT);
			AppendString(str, ind.name);
			AppendString(str, "\t");
			AppendString(str, apops<REAL>[ind.popid]->name);

			for (int64 l = 0; l < nloc; ++l)
			{
				GENOTYPE& gt = ind.GetGenotype(l);
				ushort* als = gt.GetAlleleArray();
				int ploidy = gt.Ploidy(), nalleles = gt.Nalleles();
				if (nalleles == 0 || ploidy == 1)
					AppendString(str, "\t000000");
				else
				{
					int a1 = rng.Next(ploidy), a2 = rng.NextAvoid(ploidy, a1);
					sprintf(str, "\t%03d%03d", als[a1] + 1, als[a2] + 1);
					while (*str) str++;
				}
			}

			AppendString(str, "\r\n");
			break;
		}
		case 4://split
		{
			int vmax = ind.vmax, nsplit = vmax / 2, nsplit2 = nsplit * 2, permals[N_MAX_PLOIDY];
			ushort* permals16 = (ushort*)permals;

			for (int j = 0; j < nsplit2; j += 2)
			{
				AppendString(str, ind.name);
				AppendString(str, "A\t"); str[-2] += j / 2;
				AppendString(str, apops<REAL>[ind.popid]->name);

				for (int64 l = 0; l < nloc; ++l)
				{
					GENOTYPE& gt = ind.GetGenotype(l);
					ushort* als = gt.GetAlleleArray();
					int ploidy = gt.Ploidy(), nalleles = gt.Nalleles();
					if (nalleles == 0 || ploidy == 1 || j + 1 >= ploidy)
						AppendString(str, "\t000000");
					else
					{
						sprintf(str, "\t%03d%03d", als[j] + 1, als[j + 1] + 1);
						while (*str) str++;
					}
				}

				AppendString(str, "\r\n");
			}
			break;
		}
		case 5://shuffle
		{
			int vmax = ind.vmax, nsplit = vmax / 2, nsplit2 = nsplit * 2, permals[N_MAX_PLOIDY];
			ushort* permals16 = (ushort*)permals;

			for (int j = 0; j < nsplit2; j += 2)
			{
				RNG<double> rng(g_seed_val + ii, RNG_SALT_CONVERT);
				AppendString(str, ind.name);
				AppendString(str, "A\t"); str[-2] += j / 2;
				AppendString(str, apops<REAL>[ind.popid]->name);

				for (int64 l = 0; l < nloc; ++l)
				{
					GENOTYPE& gt = ind.GetGenotype(l);
					ushort* als = gt.GetAlleleArray();
					int ploidy = gt.Ploidy(), nalleles = gt.Nalleles();

					if (nalleles == 0 || ploidy == 1)
						AppendString(str, "\t000000");
					else
					{
						SetVal(permals16, als, ploidy);
						rng.Permute(permals16, ploidy);
						SetFF(permals16 + ploidy - ploidy % 2, nsplit2 - ploidy + ploidy % 2);
						rng.Permute(permals, nsplit);
						sprintf(str, "\t%03d%03d", (short)permals16[j] + 1, (short)permals16[j + 1] + 1);
						while (*str) str++;
					}
				}

				AppendString(str, "\r\n");
			}
			break;
		}
		}

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

	int noutput = nind;
	if (convert_mode_val >= 4)
	{
		noutput = 0;
		for (int j = 0; j < nind; ++j)
		{
			IND<REAL>* ind = ainds<REAL>[j];
			noutput += ind->vmax <= 2 ? 1 : ind->vmax / 2;
		}
	}

	fprintf(convert_file, "%d	%d	%d	%u	3\r\n", noutput, npop, nloc, maxploidy);

	for (int i = 0; i < npop; ++i)
	{
		POP<REAL>* tp = apops<REAL>[i];
		fprintf(convert_file, "%s", tp->name);
		for (int rl = 0; rl <= lreg; ++rl)
		{
			tp = aregs<REAL>[rl][tp->rid];
			fprintf(convert_file, "\t%s", tp->name);
		}
		fprintf(convert_file, "\r\n");
	}

	fprintf(convert_file, "Population	Individual");
	for (int64 l = 0; l < nloc; ++l)
		fprintf(convert_file, "\t%s", GetLoc(l).GetNameStr(name_buf));
	fprintf(convert_file, "\r\n");

	convert_linesize =
		convert_mode_val == 1 ? IND_NAME_LEN + (maxploidy * MAX_ALLELE_BYTE + 1) * nloc :
		(convert_mode_val <= 3 ? IND_NAME_LEN + (MAX_ALLELE_BYTE * 2 + 1) * nloc :
		(IND_NAME_LEN + (MAX_ALLELE_BYTE * 2 + 1) * nloc) * (maxploidy / 2));

	for (int j = 0; j < NBUF; ++j)
		conversion_memory2->Alloc(convert_buf[j], convert_linesize);

	PrepareGenotypeString(8);

	RunThreads(&ConvertGenoDiveInd<REAL>, &ConvertGuard<REAL>, NULL, ntot, nind,
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
		IND<REAL>& ind = *rinds<REAL>[ii];

		switch (convert_mode_val)
		{
		default:
		case 1:
		case 2://disable,truncate
		{
			AppendInt(str, ind.popid + 1);
			AppendString(str, "\t");
			AppendString(str, ind.name);

			for (int64 l = 0; l < nloc; ++l)
				AppendString(str, conversion_string[l][ind.GetGenotypeId(l)]);

			AppendString(str, "\r\n");
			break;
		}
		case 3://choose
		{
			RNG<double> rng(g_seed_val + ii, RNG_SALT_CONVERT);
			AppendInt(str, ind.popid + 1);
			AppendString(str, "\t");
			AppendString(str, ind.name);

			for (int64 l = 0; l < nloc; ++l)
			{
				GENOTYPE& gt = ind.GetGenotype(l);
				ushort* als = gt.GetAlleleArray();
				int ploidy = gt.Ploidy(), nalleles = gt.Nalleles();
				if (nalleles == 0 || ploidy == 1)
					AppendString(str, "\t000000");
				else
				{
					int a1 = rng.Next(ploidy), a2 = rng.NextAvoid(ploidy, a1);
					sprintf(str, "\t%03d%03d", als[a1] + 1, als[a2] + 1);
					while (*str) str++;
				}
			}

			AppendString(str, "\r\n");
			break;
		}
		case 4://split
		{
			int vmax = ind.vmax, nsplit = vmax / 2, nsplit2 = nsplit * 2, permals[N_MAX_PLOIDY];
			ushort* permals16 = (ushort*)permals;

			for (int j = 0; j < nsplit2; j += 2)
			{
				AppendInt(str, ind.popid + 1);
				AppendString(str, "\t");
				AppendString(str, ind.name);
				AppendString(str, "A"); str[-1] += j / 2;

				for (int64 l = 0; l < nloc; ++l)
				{
					GENOTYPE& gt = ind.GetGenotype(l);
					ushort* als = gt.GetAlleleArray();
					int ploidy = gt.Ploidy(), nalleles = gt.Nalleles();
					if (nalleles == 0 || ploidy == 1 || j + 1 >= ploidy)
						AppendString(str, "\t000000");
					else
					{
						sprintf(str, "\t%03d%03d", als[j] + 1, als[j + 1] + 1);
						while (*str) str++;
					}
				}

				AppendString(str, "\r\n");
			}
			break;
		}
		case 5://shuffle
		{
			int vmax = ind.vmax, nsplit = vmax / 2, nsplit2 = nsplit * 2, permals[N_MAX_PLOIDY];
			ushort* permals16 = (ushort*)permals;

			for (int j = 0; j < nsplit2; j += 2)
			{
				RNG<double> rng(g_seed_val + ii, RNG_SALT_CONVERT);
				AppendInt(str, ind.popid + 1);
				AppendString(str, "\t");
				AppendString(str, ind.name);
				AppendString(str, "A"); str[-1] += j / 2;

				for (int64 l = 0; l < nloc; ++l)
				{
					GENOTYPE& gt = ind.GetGenotype(l);
					ushort* als = gt.GetAlleleArray();
					int ploidy = gt.Ploidy(), nalleles = gt.Nalleles();

					if (nalleles == 0 || ploidy == 1)
						AppendString(str, "\t000000");
					else
					{
						SetVal(permals16, als, ploidy);
						rng.Permute(permals16, ploidy);
						SetFF(permals16 + ploidy - ploidy % 2, nsplit2 - ploidy + ploidy % 2);
						rng.Permute(permals, nsplit);
						sprintf(str, "\t%03d%03d", (short)permals16[j] + 1, (short)permals16[j + 1] + 1);
						while (*str) str++;
					}
				}

				AppendString(str, "\r\n");
			}
			break;
		}
		}

		*str++ = '\0';

		THREAD_END
	}
}

/* Convert into plink format */
template<typename REAL>
TARGET void ConvertPlink(int ntot, bool& isfirst)
{
	char filename[PATH_LEN];
	char name_buf[NAME_BUF_LEN];

	time_t time1;
	time(&time1);
	FRES_TIME = localtime(&time1);
	FRES = FOpen(FRES_NAME, "wb", "%s.%s", g_output_val.c_str(), "convert.plink.map");
	OpenTempFiles(g_nthread_val, ".convert.plink");
	convert_file = FOpen(filename, "wb", "%s.%s", g_output_val.c_str(), "convert.plink.ped");

	convert_linesize = convert_mode_val <= 3 ?
		IND_NAME_LEN + (MAX_ALLELE_BYTE * 2 + 2) * nloc + 12 :
		(IND_NAME_LEN + (MAX_ALLELE_BYTE * 2 + 2) * nloc + 12) * (maxploidy / 2);

	for (int j = 0; j < NBUF; ++j)
		conversion_memory2->Alloc(convert_buf[j], convert_linesize);

	PrepareGenotypeString(9);

	RunThreads(&ConvertPlinkInd<REAL>, &ConvertGuard<REAL>, NULL, ntot, nind,
		"\nConverting population genetics software format:\n", g_nthread_val, isfirst);
	
	isfirst = false;
	fclose(convert_file);
	JoinTempFiles(g_nthread_val);
	CloseResFile();
}

/* Convert individual genotypes into plink format in multiple threads */
THREAD2(ConvertPlinkInd)
{
	// map file
	int64 lst = threadid * nloc / g_nthread_val, led = (threadid + 1) * nloc / g_nthread_val;
	FILE* ftmp = TEMP_FILES[threadid];	
	char name_buf[NAME_BUF_LEN];

	for (int64 l = lst; l < led; ++l)
	{
		//chrom	locus	genetic pos	physical pos
		fprintf(ftmp, "%s\t%s\t0\t%lld\r\n",
			GetLoc(l).GetChrom(),
			GetLoc(l).GetNameStr(name_buf),
			uselocpos ? GetLocPos(l) : 0ll);
	}

	// individual genotype file
	for (int64 ii = 0; ii < nind; ++ii)
	{
		THREAD_BEGIN

		char* str = convert_buf[ii % NBUF];
		IND<REAL>& ind = *rinds<REAL>[ii];

		switch (convert_mode_val)
		{
		default:
		case 1:
		case 2://disable,truncate
		{
			AppendString(str, apops<REAL>[ind.popid]->name);
			AppendString(str, "\t");
			AppendString(str, ind.name);
			AppendString(str, "\t0\t0\t0\t0");

			for (int64 l = 0; l < nloc; ++l)
				AppendString(str, conversion_string[l][ind.GetGenotypeId(l)]);

			AppendString(str, "\r\n");

			break;
		}
		case 3://choose
		{
			RNG<double> rng(g_seed_val + ii, RNG_SALT_CONVERT);
			AppendString(str, apops<REAL>[ind.popid]->name);
			AppendString(str, "\t");
			AppendString(str, ind.name);
			AppendString(str, "\t0\t0\t0\t0");

			for (int64 l = 0; l < nloc; ++l)
			{
				GENOTYPE& gt = ind.GetGenotype(l);
				ushort* als = gt.GetAlleleArray();
				int ploidy = gt.Ploidy(), nalleles = gt.Nalleles();
				if (nalleles == 0 || ploidy == 1)
					AppendString(str, "\t0 0");
				else
				{
					int a1 = rng.Next(ploidy), a2 = rng.NextAvoid(ploidy, a1);
					sprintf(str, "\t%d %d", als[a1] + 1, als[a2] + 1);
					while (*str) str++;
				}
			}

			AppendString(str, "\r\n");
			break;
		}
		case 4://split
		{
			int vmax = ind.vmax, nsplit = vmax / 2, nsplit2 = nsplit * 2, permals[N_MAX_PLOIDY];
			ushort* permals16 = (ushort*)permals;

			for (int j = 0; j < nsplit2; j += 2)
			{
				AppendString(str, apops<REAL>[ind.popid]->name);
				AppendString(str, "\t");
				AppendString(str, ind.name);
				AppendString(str, "A\t0\t0\t0\t0"); str[-9] += j / 2;

				for (int64 l = 0; l < nloc; ++l)
				{
					GENOTYPE& gt = ind.GetGenotype(l);
					ushort* als = gt.GetAlleleArray();
					int ploidy = gt.Ploidy(), nalleles = gt.Nalleles();
					if (nalleles == 0 || ploidy == 1 || j + 1 >= ploidy)
						AppendString(str, "\t0 0");
					else
					{
						sprintf(str, "\t%d %d", als[j] + 1, als[j + 1] + 1);
						while (*str) str++;
					}
				}

				AppendString(str, "\r\n");
			}
			break;
		}
		case 5://shuffle
		{
			int vmax = ind.vmax, nsplit = vmax / 2, nsplit2 = nsplit * 2, permals[N_MAX_PLOIDY];
			ushort* permals16 = (ushort*)permals;

			for (int j = 0; j < nsplit2; j += 2)
			{
				RNG<double> rng(g_seed_val + ii, RNG_SALT_CONVERT);
				AppendString(str, apops<REAL>[ind.popid]->name);
				AppendString(str, "\t");
				AppendString(str, ind.name);
				AppendString(str, "A\t0\t0\t0\t0"); str[-9] += j / 2;

				for (int64 l = 0; l < nloc; ++l)
				{
					GENOTYPE& gt = ind.GetGenotype(l);
					ushort* als = gt.GetAlleleArray();
					int ploidy = gt.Ploidy(), nalleles = gt.Nalleles();

					if (nalleles == 0 || ploidy == 1)
						AppendString(str, "\t0 0");
					else
					{
						SetVal(permals16, als, ploidy);
						rng.Permute(permals16, ploidy);
						SetFF(permals16 + ploidy - ploidy % 2, nsplit2 - ploidy + ploidy % 2);
						rng.Permute(permals, nsplit);
						sprintf(str, "\t%d %d", (short)permals16[j] + 1, (short)permals16[j + 1] + 1);
						while (*str) str++;
					}
				}

				AppendString(str, "\r\n");
			}
			break;
		}
		}

		*str++ = '\0';

		THREAD_END
	}
}