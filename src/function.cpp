/* Functions */

#include "vcfpop.h"

template struct IND<double>;
template struct IND<float >;
template struct LOCSTAT2<double>;
template struct LOCSTAT2<float >;
template struct POP<double>;
template struct POP<float >;

template TARGET GENOTYPE& IND<double>::GetGenotype(int64 l);
template TARGET GENOTYPE& IND<float >::GetGenotype(int64 l);
template TARGET GENOTYPE& IND<double>::GetGenotype(int64 l, GENOTYPE* gtab);
template TARGET GENOTYPE& IND<float >::GetGenotype(int64 l, GENOTYPE* gtab);
template TARGET int IND<double>::GetGenotypeId(int64 l);
template TARGET int IND<float >::GetGenotypeId(int64 l);
template TARGET int IND<double>::GetGenotypeId(int64 l, byte* bucket, OFFSET* offset);
template TARGET int IND<float >::GetGenotypeId(int64 l, byte* bucket, OFFSET* offset);
template TARGET IND<double>::IND();
template TARGET IND<float >::IND();
template TARGET IND<double>::IND(char* t, bool iscount, int id, GENOTYPE** gtab, ushort** gatab, GENO_WRITER* wt);
template TARGET IND<float >::IND(char* t, bool iscount, int id, GENOTYPE** gtab, ushort** gatab, GENO_WRITER* wt);
template TARGET IND<double>::IND(char*& title, int id);
template TARGET IND<float >::IND(char*& title, int id);
template TARGET IND<double>::IND(IND<double>& ref);
template TARGET IND<float >::IND(IND<float >& ref);
template TARGET double IND<double>::GenoFreq(POP<double>* grp, int model, int64 loc, double e);
template TARGET double IND<float >::GenoFreq(POP<float >* grp, int model, int64 loc, double e);
template TARGET void IND<double>::GetDyadGenotypeIdx(int& id1, int& id2, int64 l);
template TARGET void IND<float >::GetDyadGenotypeIdx(int& id1, int& id2, int64 l);

template TARGET double* POP<double>::GetFreq(int64 l);
template TARGET float * POP<float >::GetFreq(int64 l);
template TARGET double POP<double>::GetFreq(int64 l, int a);
template TARGET float  POP<float >::GetFreq(int64 l, int a);
template TARGET ushort* POP<double>::GetGenoCount(int64 l);
template TARGET ushort* POP<float >::GetGenoCount(int64 l);
template TARGET ushort POP<double>::GetGenoCount(int64 l, int id);
template TARGET ushort POP<float >::GetGenoCount(int64 l, int id);
template TARGET void POP<double>::AllocFreq();
template TARGET void POP<float >::AllocFreq();
template TARGET void POP<double>::UnAllocFreq();
template TARGET void POP<float >::UnAllocFreq();
template TARGET void POP<double>::GetLocStat2(LOCSTAT2<double>* loc_stat2);
template TARGET void POP<float >::GetLocStat2(LOCSTAT2<float >* loc_stat2);
template TARGET bool POP<double>::IsSubpop(POP<double>* pop);
template TARGET bool POP<float >::IsSubpop(POP<float >* pop);
template TARGET void POP<double>::CalcFreqGcount();
template TARGET void POP<float >::CalcFreqGcount();

template TARGET void Initialize<double>();
template TARGET void Initialize<float >();
template TARGET void UnInitialize1<double>();
template TARGET void UnInitialize1<float >();
template TARGET void UnInitialize2<double>();
template TARGET void UnInitialize2<float >();
template TARGET void AssignPloidy<double>();
template TARGET void AssignPloidy<float >();
template TARGET void CalcFreq<double>();
template TARGET void CalcFreq<float >();
template TARGET void EstimatePES<double>();
template TARGET void EstimatePES<float >();
template TARGET void Calculate<double>();
template TARGET void Calculate<float >();
template TARGET double GENOTYPE::GetFreq<double>(int k2);
template TARGET float  GENOTYPE::GetFreq<float >(int k2);
template TARGET void GENOTYPE::GetFreq<double>(double* p, int k2);
template TARGET void GENOTYPE::GetFreq<float >(float * p, int k2);

#define extern 

/* Misc */
extern void* cpop_;									//Current population
#define cpop (*(POP<REAL>**)&cpop_)
extern ushort missing_array[N_MAX_PLOIDY];			//Allele array of the missing genotypes
extern GENOTYPE missing_genotype[N_MAX_PLOIDY + 1];	//Missing genotype at different ploidy level
extern HASH missing_hash[N_MAX_PLOIDY + 1];			//Hash of missing genotype

/* Global */
extern SLOCUS* slocus;								//SLocus information
extern bool useslocus;								//Use small locus
extern TABLE<HASH, GENOTYPE*> emptygftab;			//Empty genotype table

/* Load File */
template<> extern LIST<POP<double>> pop<double>;						//Original input population
template<> extern LIST<POP<float >> pop<float >;						//Original input population
template<> extern LIST<LIST<POP<double>>> reg<double>; 				//Original input regions
template<> extern LIST<LIST<POP<float >>> reg<float >; 				//Original input regions

extern LOCUS* locus;								//Locus information
extern MEMORY* individual_memory;					//Individual memory class
extern MEMORY* locus_memory;						//Locus memory class
extern MEMORY* nvcf_memory;							//Locus memory for first round counting ngeno
extern TABLE<HASH, uint>* nvcf_gfid;				//Hash table for counting ngeno 

/* Allele frequency and genotype count offset */
extern uint64* allele_freq_offset;					//Allele frequency array at locus l
extern int maxK;									//Max number of alleles
extern int64 KT;									//Total number of alleles

extern uint64* genotype_count_offset;					//Genotype count array at locus l
extern int maxG;									//Max number of genotypes at all loci
extern int64 GT;									//Total number of genotypes


/* Genotype bucket */
extern VMEMORY locus_list;
extern BUCKET geno_bucket;
extern BUCKET haplo_bucket;
extern BUCKET ad_bucket;

/* Reassign individuals and populations */
extern bool reassigned;								//Is ploidy assigned, to distiguish diversity filter and diversity estimation
extern void* rinds_;						//Rearranged individuals according to population source
#define rinds (*(IND<REAL>***)&rinds_)
#undef extern

#ifndef _BCFHEADER
	
	/* Uninitialize BCF header */
	TARGET BCFHEADER::~BCFHEADER()
	{
		if (contig_name == NULL) return;
		for (int64 i = 0; i < contig_size; ++i)
		{
			delete[] contig_name[i];
			contig_name[i] = NULL;
		}
		delete[] contig_name;
		contig_name = NULL;
	}

	/* Initialize BCF header */
	TARGET BCFHEADER::BCFHEADER(char* header)
	{
		int64 line = 0;
		contig_size = 0;
		contig_name = NULL;
		filter_passidx = format_gtid = format_gqid = format_dpid = format_adid = 0xFFFFFFFF;
		char* bak = header;

		int filteridx = 0, fmtidx = 0, contigidx = 0;
		while (*header)
		{
			line++;
			char* hbak = header;
			while (*header == '#') header++;

			if (LwrLineCmp("filter=<", header) == 0)
			{
				header = LineNextIdx(header, "ID=", 1);
				if (header == NULL)
				{
					header = StrNextIdx(header, '\n', 1); if (header) *header = '\0';
					Exit("\nError: bcf header format error, at line %d:\n%s\n", line, hbak);
				}
				header += 3;
				if (LwrLineCmp("pass", header) == 0)
				{
					header = LineNextIdx(header, "IDX=", 1) + 4;
					if (header == (char*)4)
						filter_passidx = ReadInteger(header);
					else
						filter_passidx = filteridx;
				}
				filteridx++;
			}
			else if (LwrLineCmp("contig=<", header) == 0)
			{
				header = LineNextIdx(header, "IDX=", 1) + 4;
				int64 idx = header == (char*)4 ? contigidx : ReadInteger(header);
				contig_size = Max(idx + 1, contig_size);
				contigidx++;
			}
			else if (LwrLineCmp("format=<", header) == 0)
			{
				header = LineNextIdx(header, "ID=", 1);
				if (header == NULL)
				{
					header = StrNextIdx(header, '\n', 1); if (header) *header = '\0';
					Exit("\nError: bcf header format error, at line %d:\n%s\n", line, hbak);
				}
				header += 3;

				if ((LwrLineCmp("gt", header) && LwrLineCmp("dp", header) && LwrLineCmp("gq", header) && LwrLineCmp("ad", header)) == 0)
				{
					char ftype = header[1];
					header = LineNextIdx(header, "IDX=", 1) + 4;
					switch (ftype)
					{
					case 'T':
					case 't': format_gtid = header == (char*)4 ? fmtidx : ReadInteger(header); break;
					case 'Q':
					case 'q': format_gqid = header == (char*)4 ? fmtidx : ReadInteger(header); break;
					case 'P':
					case 'p': format_dpid = header == (char*)4 ? fmtidx : ReadInteger(header); break;
					case 'D':
					case 'd': format_adid = header == (char*)4 ? fmtidx : ReadInteger(header); break;
					}
				}
				fmtidx++;
			}
			header = StrNextIdx(header, '\n', 1) + 1;
		}

		if (format_gtid == -1)
			Exit("\nError: bcf header format error, no GT format field found.", line);

		if (contig_size == 0)
			Exit("\nError: bcf header format error, no contig/chromosome found.", line);

		contig_name = new char*[contig_size];

		line = 0; contigidx = 0; header = bak;
		while (*header)
		{
			line++;
			while (*header == '#') header++;

			if (LwrLineCmp("contig=<", header) == 0)
			{
				header = LineNextIdx(header, "ID=", 1) + 3;
				int64 tlen = StrNextIdx(header, ',', 1) - header;
				char* tstr = new char[tlen + 1];
				memmove(tstr, header, tlen);
				tstr[tlen] = '\0';

				header = LineNextIdx(header, "IDX=", 1) + 4;
				int64 idx = 0;
				if (header == (char*)4)
					idx = contigidx;
				else
					idx = ReadInteger(header);
				contig_name[idx] = tstr;
				contigidx++;
			}
			header = StrNextIdx(header, '\n', 1) + 1;
		}
	}
#endif

#ifndef _GENOTYPE
	/* Do nothing */
	TARGET GENOTYPE::GENOTYPE()
	{

	}

	/* Create genotype from alleles and hash */
	TARGET GENOTYPE::GENOTYPE(ushort*& gatab, ushort* alleles, int ploidy)
	{
		//memory == NULL for reconstruct: index alleles for non-vcf input
		//alleles should not be sorted for phased genoytpe
		
		patternid = 0;
		
		//zero ploidy, not set
		if (ploidy == 0)
		{
			SetAlleleArray(0);
			return;
		}

		if (alleles[0] == 0xFFFF)
		{
			//missing but have ploidy info
			patternid = N_PATTERN_END + (uint)ploidy;
			SetAlleleArray((ushort*)-1);
			return;
		}

		byte acount[N_MAX_PLOIDY] = { 0 };
		ushort akeys[N_MAX_PLOIDY] = { 0 };
		int nalleles = 0;

		//counting number of alleles
		for (int i = 0; i < ploidy; ++i)
		{
			bool isnew = true;
			ushort asi = alleles[i];
			for (int j = i - 1; j >= 0; --j)
				if (asi == alleles[j])
				{
					isnew = false;
					break;
				}
			if (!isnew) continue;

			//asi is a new allele, count dosage
			nalleles++;

			int ac = 1;
			akeys[nalleles - 1] = asi;
			for (int j = i + 1; j < ploidy; ++j)
				if (asi == alleles[j]) ac++;
			acount[nalleles - 1] = (byte)ac;
		}

		//sorting alleles in descending order of dosage
		for (int i = 0; i < nalleles; ++i)
			for (int j = i + 1; j < nalleles; ++j)
				if (acount[i] < acount[j] ||
					(acount[i] == acount[j] && akeys[i] < akeys[j]))
				{
					Swap(acount[i], acount[j]);
					Swap(akeys[i], akeys[j]);
				}

		ushort* als;
		if (gatab)
		{
			als = gatab;
			gatab += ploidy + nalleles;
			SetAlleleArray(als);
		}
		else //reconstruct, replace in indexing alleles, do not allloc new memory
			als = GetAlleleArray();

		SetVal(als, alleles, ploidy);
		SetVal(als + ploidy, akeys, nalleles);

		//allele pattern
		uint64 pattern = 0;
		for (int i = 0; i < nalleles; ++i)
			pattern = (pattern << 4) | acount[i];

		switch (pattern)
		{
		default: patternid = 0; break; //missing but have ploidy info
		case 0x1ull: patternid = 1; break;
		case 0x2ull: patternid = 2; break;
		case 0x11ull: patternid = 3; break;
		case 0x3ull: patternid = 4; break;
		case 0x21ull: patternid = 5; break;
		case 0x111ull: patternid = 6; break;
		case 0x4ull: patternid = 7; break;
		case 0x31ull: patternid = 8; break;
		case 0x22ull: patternid = 9; break;
		case 0x211ull: patternid = 10; break;
		case 0x1111ull: patternid = 11; break;
		case 0x5ull: patternid = 12; break;
		case 0x41ull: patternid = 13; break;
		case 0x32ull: patternid = 14; break;
		case 0x311ull: patternid = 15; break;
		case 0x221ull: patternid = 16; break;
		case 0x2111ull: patternid = 17; break;
		case 0x11111ull: patternid = 18; break;
		case 0x6ull: patternid = 19; break;
		case 0x51ull: patternid = 20; break;
		case 0x42ull: patternid = 21; break;
		case 0x411ull: patternid = 22; break;
		case 0x33ull: patternid = 23; break;
		case 0x321ull: patternid = 24; break;
		case 0x3111ull: patternid = 25; break;
		case 0x222ull: patternid = 26; break;
		case 0x2211ull: patternid = 27; break;
		case 0x21111ull: patternid = 28; break;
		case 0x111111ull: patternid = 29; break;
		case 0x7ull: patternid = 30; break;
		case 0x61ull: patternid = 31; break;
		case 0x52ull: patternid = 32; break;
		case 0x511ull: patternid = 33; break;
		case 0x43ull: patternid = 34; break;
		case 0x421ull: patternid = 35; break;
		case 0x4111ull: patternid = 36; break;
		case 0x331ull: patternid = 37; break;
		case 0x322ull: patternid = 38; break;
		case 0x3211ull: patternid = 39; break;
		case 0x31111ull: patternid = 40; break;
		case 0x2221ull: patternid = 41; break;
		case 0x22111ull: patternid = 42; break;
		case 0x211111ull: patternid = 43; break;
		case 0x1111111ull: patternid = 44; break;
		case 0x8ull: patternid = 45; break;
		case 0x71ull: patternid = 46; break;
		case 0x62ull: patternid = 47; break;
		case 0x611ull: patternid = 48; break;
		case 0x53ull: patternid = 49; break;
		case 0x521ull: patternid = 50; break;
		case 0x5111ull: patternid = 51; break;
		case 0x44ull: patternid = 52; break;
		case 0x431ull: patternid = 53; break;
		case 0x422ull: patternid = 54; break;
		case 0x4211ull: patternid = 55; break;
		case 0x41111ull: patternid = 56; break;
		case 0x332ull: patternid = 57; break;
		case 0x3311ull: patternid = 58; break;
		case 0x3221ull: patternid = 59; break;
		case 0x32111ull: patternid = 60; break;
		case 0x311111ull: patternid = 61; break;
		case 0x2222ull: patternid = 62; break;
		case 0x22211ull: patternid = 63; break;
		case 0x221111ull: patternid = 64; break;
		case 0x2111111ull: patternid = 65; break;
		case 0x11111111ull: patternid = 66; break;
		case 0x9ull: patternid = 67; break;
		case 0x81ull: patternid = 68; break;
		case 0x72ull: patternid = 69; break;
		case 0x711ull: patternid = 70; break;
		case 0x63ull: patternid = 71; break;
		case 0x621ull: patternid = 72; break;
		case 0x6111ull: patternid = 73; break;
		case 0x54ull: patternid = 74; break;
		case 0x531ull: patternid = 75; break;
		case 0x522ull: patternid = 76; break;
		case 0x5211ull: patternid = 77; break;
		case 0x51111ull: patternid = 78; break;
		case 0x441ull: patternid = 79; break;
		case 0x432ull: patternid = 80; break;
		case 0x4311ull: patternid = 81; break;
		case 0x4221ull: patternid = 82; break;
		case 0x42111ull: patternid = 83; break;
		case 0x411111ull: patternid = 84; break;
		case 0x333ull: patternid = 85; break;
		case 0x3321ull: patternid = 86; break;
		case 0x33111ull: patternid = 87; break;
		case 0x3222ull: patternid = 88; break;
		case 0x32211ull: patternid = 89; break;
		case 0x321111ull: patternid = 90; break;
		case 0x3111111ull: patternid = 91; break;
		case 0x22221ull: patternid = 92; break;
		case 0x222111ull: patternid = 93; break;
		case 0x2211111ull: patternid = 94; break;
		case 0x21111111ull: patternid = 95; break;
		case 0x111111111ull: patternid = 96; break;
		case 0xAull: patternid = 97; break;
		case 0x91ull: patternid = 98; break;
		case 0x82ull: patternid = 99; break;
		case 0x811ull: patternid = 100; break;
		case 0x73ull: patternid = 101; break;
		case 0x721ull: patternid = 102; break;
		case 0x7111ull: patternid = 103; break;
		case 0x64ull: patternid = 104; break;
		case 0x631ull: patternid = 105; break;
		case 0x622ull: patternid = 106; break;
		case 0x6211ull: patternid = 107; break;
		case 0x61111ull: patternid = 108; break;
		case 0x55ull: patternid = 109; break;
		case 0x541ull: patternid = 110; break;
		case 0x532ull: patternid = 111; break;
		case 0x5311ull: patternid = 112; break;
		case 0x5221ull: patternid = 113; break;
		case 0x52111ull: patternid = 114; break;
		case 0x511111ull: patternid = 115; break;
		case 0x442ull: patternid = 116; break;
		case 0x4411ull: patternid = 117; break;
		case 0x433ull: patternid = 118; break;
		case 0x4321ull: patternid = 119; break;
		case 0x43111ull: patternid = 120; break;
		case 0x4222ull: patternid = 121; break;
		case 0x42211ull: patternid = 122; break;
		case 0x421111ull: patternid = 123; break;
		case 0x4111111ull: patternid = 124; break;
		case 0x3331ull: patternid = 125; break;
		case 0x3322ull: patternid = 126; break;
		case 0x33211ull: patternid = 127; break;
		case 0x331111ull: patternid = 128; break;
		case 0x32221ull: patternid = 129; break;
		case 0x322111ull: patternid = 130; break;
		case 0x3211111ull: patternid = 131; break;
		case 0x31111111ull: patternid = 132; break;
		case 0x22222ull: patternid = 133; break;
		case 0x222211ull: patternid = 134; break;
		case 0x2221111ull: patternid = 135; break;
		case 0x22111111ull: patternid = 136; break;
		case 0x211111111ull: patternid = 137; break;
		case 0x1111111111ull: patternid = 138; break;
		}
	}

	/* Copy from a reference genotype */
	TARGET GENOTYPE::GENOTYPE(ushort*& gatab, GENOTYPE& ref)
	{
		patternid = ref.patternid;
		int nalleles = Nalleles();

		if (nalleles)
		{
			int ploidy = Ploidy();
			ushort* als = gatab;
			gatab += ploidy + nalleles;
			SetVal(als, ref.GetAlleleArray(), ploidy + nalleles);
			SetAlleleArray(als);
		}
		else
			SetAlleleArray((ushort*)-1);
	}

	/* Get allele copy at ith haplotype */
	TARGET ushort GENOTYPE::GetAlleleCopy(int i)
	{
		switch (offset)
		{
		case 0:
			return NULL;
		case 0xFFFFFF:
			return 0xFFFF;
		default:
			return *((ushort*)((byte*)this + offset) + i);
		}
	}

	/* Get allele array */
	TARGET ushort* GENOTYPE::GetAlleleArray()
	{
		switch (offset)
		{
		case 0:
			return NULL;
		case 0xFFFFFF:
			return missing_array;
		default:
			return (ushort*)((byte*)this + offset);
		}
	}

	/* Set allele array */
	TARGET void GENOTYPE::SetAlleleArray(ushort* alleles)
	{
		if (alleles == NULL) //not set
			offset = 0;
		else if (alleles == (ushort*)-1)//empty
			offset = 0xFFFFFF;
		else
			offset = (uint)((byte*)alleles - (byte*)this);
	}

	/* Get pattern code */
	TARGET uint64 GENOTYPE::GetPattern()
	{
		static uint64 PT_PATTERN[150] =									//Pattern index to pattern
		{ 0x0ull, 0x1ull, 0x2ull, 0x11ull, 0x3ull, 0x21ull, 0x111ull, 0x4ull, 0x31ull, 0x22ull, 0x211ull, 0x1111ull, 0x5ull, 0x41ull, 0x32ull, 0x311ull, 0x221ull, 0x2111ull, 0x11111ull, 0x6ull, 0x51ull, 0x42ull, 0x411ull, 0x33ull, 0x321ull, 0x3111ull, 0x222ull, 0x2211ull, 0x21111ull, 0x111111ull, 0x7ull, 0x61ull, 0x52ull, 0x511ull, 0x43ull, 0x421ull, 0x4111ull, 0x331ull, 0x322ull, 0x3211ull, 0x31111ull, 0x2221ull, 0x22111ull, 0x211111ull, 0x1111111ull, 0x8ull, 0x71ull, 0x62ull, 0x611ull, 0x53ull, 0x521ull, 0x5111ull, 0x44ull, 0x431ull, 0x422ull, 0x4211ull, 0x41111ull, 0x332ull, 0x3311ull, 0x3221ull, 0x32111ull, 0x311111ull, 0x2222ull, 0x22211ull, 0x221111ull, 0x2111111ull, 0x11111111ull, 0x9ull, 0x81ull, 0x72ull, 0x711ull, 0x63ull, 0x621ull, 0x6111ull, 0x54ull, 0x531ull, 0x522ull, 0x5211ull, 0x51111ull, 0x441ull, 0x432ull, 0x4311ull, 0x4221ull, 0x42111ull, 0x411111ull, 0x333ull, 0x3321ull, 0x33111ull, 0x3222ull, 0x32211ull, 0x321111ull, 0x3111111ull, 0x22221ull, 0x222111ull, 0x2211111ull, 0x21111111ull, 0x111111111ull, 0xAull, 0x91ull, 0x82ull, 0x811ull, 0x73ull, 0x721ull, 0x7111ull, 0x64ull, 0x631ull, 0x622ull, 0x6211ull, 0x61111ull, 0x55ull, 0x541ull, 0x532ull, 0x5311ull, 0x5221ull, 0x52111ull, 0x511111ull, 0x442ull, 0x4411ull, 0x433ull, 0x4321ull, 0x43111ull, 0x4222ull, 0x42211ull, 0x421111ull, 0x4111111ull, 0x3331ull, 0x3322ull, 0x33211ull, 0x331111ull, 0x32221ull, 0x322111ull, 0x3211111ull, 0x31111111ull, 0x22222ull, 0x222211ull, 0x2221111ull, 0x22111111ull, 0x211111111ull, 0x1111111111ull, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

		return PT_PATTERN[patternid];
		/*
		switch (patternid)
		{
		default: 	return 0ull;
		case	1: 	return 0x1ull;
		case	2: 	return 0x2ull;
		case	3: 	return 0x11ull;
		case	4: 	return 0x3ull;
		case	5: 	return 0x21ull;
		case	6: 	return 0x111ull;
		case	7: 	return 0x4ull;
		case	8: 	return 0x31ull;
		case	9: 	return 0x22ull;
		case	10: 	return 0x211ull;
		case	11: 	return 0x1111ull;
		case	12: 	return 0x5ull;
		case	13: 	return 0x41ull;
		case	14: 	return 0x32ull;
		case	15: 	return 0x311ull;
		case	16: 	return 0x221ull;
		case	17: 	return 0x2111ull;
		case	18: 	return 0x11111ull;
		case	19: 	return 0x6ull;
		case	20: 	return 0x51ull;
		case	21: 	return 0x42ull;
		case	22: 	return 0x411ull;
		case	23: 	return 0x33ull;
		case	24: 	return 0x321ull;
		case	25: 	return 0x3111ull;
		case	26: 	return 0x222ull;
		case	27: 	return 0x2211ull;
		case	28: 	return 0x21111ull;
		case	29: 	return 0x111111ull;
		case	30: 	return 0x7ull;
		case	31: 	return 0x61ull;
		case	32: 	return 0x52ull;
		case	33: 	return 0x511ull;
		case	34: 	return 0x43ull;
		case	35: 	return 0x421ull;
		case	36: 	return 0x4111ull;
		case	37: 	return 0x331ull;
		case	38: 	return 0x322ull;
		case	39: 	return 0x3211ull;
		case	40: 	return 0x31111ull;
		case	41: 	return 0x2221ull;
		case	42: 	return 0x22111ull;
		case	43: 	return 0x211111ull;
		case	44: 	return 0x1111111ull;
		case	45: 	return 0x8ull;
		case	46: 	return 0x71ull;
		case	47: 	return 0x62ull;
		case	48: 	return 0x611ull;
		case	49: 	return 0x53ull;
		case	50: 	return 0x521ull;
		case	51: 	return 0x5111ull;
		case	52: 	return 0x44ull;
		case	53: 	return 0x431ull;
		case	54: 	return 0x422ull;
		case	55: 	return 0x4211ull;
		case	56: 	return 0x41111ull;
		case	57: 	return 0x332ull;
		case	58: 	return 0x3311ull;
		case	59: 	return 0x3221ull;
		case	60: 	return 0x32111ull;
		case	61: 	return 0x311111ull;
		case	62: 	return 0x2222ull;
		case	63: 	return 0x22211ull;
		case	64: 	return 0x221111ull;
		case	65: 	return 0x2111111ull;
		case	66: 	return 0x11111111ull;
		case	67: 	return 0x9ull;
		case	68: 	return 0x81ull;
		case	69: 	return 0x72ull;
		case	70: 	return 0x711ull;
		case	71: 	return 0x63ull;
		case	72: 	return 0x621ull;
		case	73: 	return 0x6111ull;
		case	74: 	return 0x54ull;
		case	75: 	return 0x531ull;
		case	76: 	return 0x522ull;
		case	77: 	return 0x5211ull;
		case	78: 	return 0x51111ull;
		case	79: 	return 0x441ull;
		case	80: 	return 0x432ull;
		case	81: 	return 0x4311ull;
		case	82: 	return 0x4221ull;
		case	83: 	return 0x42111ull;
		case	84: 	return 0x411111ull;
		case	85: 	return 0x333ull;
		case	86: 	return 0x3321ull;
		case	87: 	return 0x33111ull;
		case	88: 	return 0x3222ull;
		case	89: 	return 0x32211ull;
		case	90: 	return 0x321111ull;
		case	91: 	return 0x3111111ull;
		case	92: 	return 0x22221ull;
		case	93: 	return 0x222111ull;
		case	94: 	return 0x2211111ull;
		case	95: 	return 0x21111111ull;
		case	96: 	return 0x111111111ull;
		case	97: 	return 0xAull;
		case	98: 	return 0x91ull;
		case	99: 	return 0x82ull;
		case	100: 	return 0x811ull;
		case	101: 	return 0x73ull;
		case	102: 	return 0x721ull;
		case	103: 	return 0x7111ull;
		case	104: 	return 0x64ull;
		case	105: 	return 0x631ull;
		case	106: 	return 0x622ull;
		case	107: 	return 0x6211ull;
		case	108: 	return 0x61111ull;
		case	109: 	return 0x55ull;
		case	110: 	return 0x541ull;
		case	111: 	return 0x532ull;
		case	112: 	return 0x5311ull;
		case	113: 	return 0x5221ull;
		case	114: 	return 0x52111ull;
		case	115: 	return 0x511111ull;
		case	116: 	return 0x442ull;
		case	117: 	return 0x4411ull;
		case	118: 	return 0x433ull;
		case	119: 	return 0x4321ull;
		case	120: 	return 0x43111ull;
		case	121: 	return 0x4222ull;
		case	122: 	return 0x42211ull;
		case	123: 	return 0x421111ull;
		case	124: 	return 0x4111111ull;
		case	125: 	return 0x3331ull;
		case	126: 	return 0x3322ull;
		case	127: 	return 0x33211ull;
		case	128: 	return 0x331111ull;
		case	129: 	return 0x32221ull;
		case	130: 	return 0x322111ull;
		case	131: 	return 0x3211111ull;
		case	132: 	return 0x31111111ull;
		case	133: 	return 0x22222ull;
		case	134: 	return 0x222211ull;
		case	135: 	return 0x2221111ull;
		case	136: 	return 0x22111111ull;
		case	137: 	return 0x211111111ull;
		case	138: 	return 0x1111111111ull;
		}
		*/
	}

	/* Number of alleles */
	TARGET int GENOTYPE::Nalleles()
	{
		static int PT_NALLELES[150] = 									//Pattern index to number of alleles
		{ 0, 1, 1, 2, 1, 2, 3, 1, 2, 2, 3, 4, 1, 2, 2, 3, 3, 4, 5, 1, 2, 2, 3, 2, 3, 4, 3, 4, 5, 6, 1, 2, 2, 3, 2, 3, 4, 3, 3, 4, 5, 4, 5, 6, 7, 1, 2, 2, 3, 2, 3, 4, 2, 3, 3, 4, 5, 3, 4, 4, 5, 6, 4, 5, 6, 7, 8, 1, 2, 2, 3, 2, 3, 4, 2, 3, 3, 4, 5, 3, 3, 4, 4, 5, 6, 3, 4, 5, 4, 5, 6, 7, 5, 6, 7, 8, 9, 1, 2, 2, 3, 2, 3, 4, 2, 3, 3, 4, 5, 2, 3, 3, 4, 4, 5, 6, 3, 4, 3, 4, 5, 4, 5, 6, 7, 4, 4, 5, 6, 5, 6, 7, 8, 5, 6, 7, 8, 9, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

		return PT_NALLELES[patternid];
		/*
		switch (pattern)
		{
			default:	return 0;
			case 0x1ull:	return 1;
			case 0x2ull:	return 1;
			case 0x11ull:	return 2;
			case 0x3ull:	return 1;
			case 0x21ull:	return 2;
			case 0x111ull:	return 3;
			case 0x4ull:	return 1;
			case 0x31ull:	return 2;
			case 0x22ull:	return 2;
			case 0x211ull:	return 3;
			case 0x1111ull:	return 4;
			case 0x5ull:	return 1;
			case 0x41ull:	return 2;
			case 0x32ull:	return 2;
			case 0x311ull:	return 3;
			case 0x221ull:	return 3;
			case 0x2111ull:	return 4;
			case 0x11111ull:	return 5;
			case 0x6ull:	return 1;
			case 0x51ull:	return 2;
			case 0x42ull:	return 2;
			case 0x411ull:	return 3;
			case 0x33ull:	return 2;
			case 0x321ull:	return 3;
			case 0x3111ull:	return 4;
			case 0x222ull:	return 3;
			case 0x2211ull:	return 4;
			case 0x21111ull:	return 5;
			case 0x111111ull:	return 6;
			case 0x7ull:	return 1;
			case 0x61ull:	return 2;
			case 0x52ull:	return 2;
			case 0x511ull:	return 3;
			case 0x43ull:	return 2;
			case 0x421ull:	return 3;
			case 0x4111ull:	return 4;
			case 0x331ull:	return 3;
			case 0x322ull:	return 3;
			case 0x3211ull:	return 4;
			case 0x31111ull:	return 5;
			case 0x2221ull:	return 4;
			case 0x22111ull:	return 5;
			case 0x211111ull:	return 6;
			case 0x1111111ull:	return 7;
			case 0x8ull:	return 1;
			case 0x71ull:	return 2;
			case 0x62ull:	return 2;
			case 0x611ull:	return 3;
			case 0x53ull:	return 2;
			case 0x521ull:	return 3;
			case 0x5111ull:	return 4;
			case 0x44ull:	return 2;
			case 0x431ull:	return 3;
			case 0x422ull:	return 3;
			case 0x4211ull:	return 4;
			case 0x41111ull:	return 5;
			case 0x332ull:	return 3;
			case 0x3311ull:	return 4;
			case 0x3221ull:	return 4;
			case 0x32111ull:	return 5;
			case 0x311111ull:	return 6;
			case 0x2222ull:	return 4;
			case 0x22211ull:	return 5;
			case 0x221111ull:	return 6;
			case 0x2111111ull:	return 7;
			case 0x11111111ull:	return 8;
			case 0x9ull:	return 1;
			case 0x81ull:	return 2;
			case 0x72ull:	return 2;
			case 0x711ull:	return 3;
			case 0x63ull:	return 2;
			case 0x621ull:	return 3;
			case 0x6111ull:	return 4;
			case 0x54ull:	return 2;
			case 0x531ull:	return 3;
			case 0x522ull:	return 3;
			case 0x5211ull:	return 4;
			case 0x51111ull:	return 5;
			case 0x441ull:	return 3;
			case 0x432ull:	return 3;
			case 0x4311ull:	return 4;
			case 0x4221ull:	return 4;
			case 0x42111ull:	return 5;
			case 0x411111ull:	return 6;
			case 0x333ull:	return 3;
			case 0x3321ull:	return 4;
			case 0x33111ull:	return 5;
			case 0x3222ull:	return 4;
			case 0x32211ull:	return 5;
			case 0x321111ull:	return 6;
			case 0x3111111ull:	return 7;
			case 0x22221ull:	return 5;
			case 0x222111ull:	return 6;
			case 0x2211111ull:	return 7;
			case 0x21111111ull:	return 8;
			case 0x111111111ull:	return 9;
			case 0xAull:	return 1;
			case 0x91ull:	return 2;
			case 0x82ull:	return 2;
			case 0x811ull:	return 3;
			case 0x73ull:	return 2;
			case 0x721ull:	return 3;
			case 0x7111ull:	return 4;
			case 0x64ull:	return 2;
			case 0x631ull:	return 3;
			case 0x622ull:	return 3;
			case 0x6211ull:	return 4;
			case 0x61111ull:	return 5;
			case 0x55ull:	return 2;
			case 0x541ull:	return 3;
			case 0x532ull:	return 3;
			case 0x5311ull:	return 4;
			case 0x5221ull:	return 4;
			case 0x52111ull:	return 5;
			case 0x511111ull:	return 6;
			case 0x442ull:	return 3;
			case 0x4411ull:	return 4;
			case 0x433ull:	return 3;
			case 0x4321ull:	return 4;
			case 0x43111ull:	return 5;
			case 0x4222ull:	return 4;
			case 0x42211ull:	return 5;
			case 0x421111ull:	return 6;
			case 0x4111111ull:	return 7;
			case 0x3331ull:	return 4;
			case 0x3322ull:	return 4;
			case 0x33211ull:	return 5;
			case 0x331111ull:	return 6;
			case 0x32221ull:	return 5;
			case 0x322111ull:	return 6;
			case 0x3211111ull:	return 7;
			case 0x31111111ull:	return 8;
			case 0x22222ull:	return 5;
			case 0x222211ull:	return 6;
			case 0x2221111ull:	return 7;
			case 0x22111111ull:	return 8;
			case 0x211111111ull:	return 9;
			case 0x1111111111ull:	return 10;
		}
		*/
		/*switch (patternid)
		{
			default: 	return 0;
			case	1: 	return 1;
			case	2: 	return 1;
			case	3: 	return 2;
			case	4: 	return 1;
			case	5: 	return 2;
			case	6: 	return 3;
			case	7: 	return 1;
			case	8: 	return 2;
			case	9: 	return 2;
			case	10: 	return 3;
			case	11: 	return 4;
			case	12: 	return 1;
			case	13: 	return 2;
			case	14: 	return 2;
			case	15: 	return 3;
			case	16: 	return 3;
			case	17: 	return 4;
			case	18: 	return 5;
			case	19: 	return 1;
			case	20: 	return 2;
			case	21: 	return 2;
			case	22: 	return 3;
			case	23: 	return 2;
			case	24: 	return 3;
			case	25: 	return 4;
			case	26: 	return 3;
			case	27: 	return 4;
			case	28: 	return 5;
			case	29: 	return 6;
			case	30: 	return 1;
			case	31: 	return 2;
			case	32: 	return 2;
			case	33: 	return 3;
			case	34: 	return 2;
			case	35: 	return 3;
			case	36: 	return 4;
			case	37: 	return 3;
			case	38: 	return 3;
			case	39: 	return 4;
			case	40: 	return 5;
			case	41: 	return 4;
			case	42: 	return 5;
			case	43: 	return 6;
			case	44: 	return 7;
			case	45: 	return 1;
			case	46: 	return 2;
			case	47: 	return 2;
			case	48: 	return 3;
			case	49: 	return 2;
			case	50: 	return 3;
			case	51: 	return 4;
			case	52: 	return 2;
			case	53: 	return 3;
			case	54: 	return 3;
			case	55: 	return 4;
			case	56: 	return 5;
			case	57: 	return 3;
			case	58: 	return 4;
			case	59: 	return 4;
			case	60: 	return 5;
			case	61: 	return 6;
			case	62: 	return 4;
			case	63: 	return 5;
			case	64: 	return 6;
			case	65: 	return 7;
			case	66: 	return 8;
			case	67: 	return 1;
			case	68: 	return 2;
			case	69: 	return 2;
			case	70: 	return 3;
			case	71: 	return 2;
			case	72: 	return 3;
			case	73: 	return 4;
			case	74: 	return 2;
			case	75: 	return 3;
			case	76: 	return 3;
			case	77: 	return 4;
			case	78: 	return 5;
			case	79: 	return 3;
			case	80: 	return 3;
			case	81: 	return 4;
			case	82: 	return 4;
			case	83: 	return 5;
			case	84: 	return 6;
			case	85: 	return 3;
			case	86: 	return 4;
			case	87: 	return 5;
			case	88: 	return 4;
			case	89: 	return 5;
			case	90: 	return 6;
			case	91: 	return 7;
			case	92: 	return 5;
			case	93: 	return 6;
			case	94: 	return 7;
			case	95: 	return 8;
			case	96: 	return 9;
			case	97: 	return 1;
			case	98: 	return 2;
			case	99: 	return 2;
			case	100: 	return 3;
			case	101: 	return 2;
			case	102: 	return 3;
			case	103: 	return 4;
			case	104: 	return 2;
			case	105: 	return 3;
			case	106: 	return 3;
			case	107: 	return 4;
			case	108: 	return 5;
			case	109: 	return 2;
			case	110: 	return 3;
			case	111: 	return 3;
			case	112: 	return 4;
			case	113: 	return 4;
			case	114: 	return 5;
			case	115: 	return 6;
			case	116: 	return 3;
			case	117: 	return 4;
			case	118: 	return 3;
			case	119: 	return 4;
			case	120: 	return 5;
			case	121: 	return 4;
			case	122: 	return 5;
			case	123: 	return 6;
			case	124: 	return 7;
			case	125: 	return 4;
			case	126: 	return 4;
			case	127: 	return 5;
			case	128: 	return 6;
			case	129: 	return 5;
			case	130: 	return 6;
			case	131: 	return 7;
			case	132: 	return 8;
			case	133: 	return 5;
			case	134: 	return 6;
			case	135: 	return 7;
			case	136: 	return 8;
			case	137: 	return 9;
			case	138: 	return 10;
		}*/
	}

	/* Ploidy level */
	TARGET int GENOTYPE::Ploidy()
	{
		static int PT_PLOIDY[150] = 									//Pattern index to ploidy level
		{ 0, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };

		return PT_PLOIDY[patternid];
		/*
		switch (pattern)
		{
		default: return 0;
		case 0x1ull: return 1;
		case 0x2ull: 
		case 0x11ull: return 2;
		case 0x3ull: 
		case 0x21ull: 
		case 0x111ull: return 3;
		case 0x4ull: 
		case 0x31ull: 
		case 0x22ull: 
		case 0x211ull: 
		case 0x1111ull: return 4;
		case 0x5ull: 
		case 0x41ull: 
		case 0x32ull: 
		case 0x311ull: 
		case 0x221ull: 
		case 0x2111ull: 
		case 0x11111ull: return 5;
		case 0x6ull: 
		case 0x51ull: 
		case 0x42ull: 
		case 0x411ull: 
		case 0x33ull: 
		case 0x321ull: 
		case 0x3111ull: 
		case 0x222ull:
		case 0x2211ull: 
		case 0x21111ull: 
		case 0x111111ull: return 6;
		case 0x7ull: 
		case 0x61ull: 
		case 0x52ull: 
		case 0x511ull: 
		case 0x43ull: 
		case 0x421ull: 
		case 0x4111ull: 
		case 0x331ull: 
		case 0x322ull:
		case 0x3211ull: 
		case 0x31111ull: 
		case 0x2221ull: 
		case 0x22111ull: 
		case 0x211111ull: 
		case 0x1111111ull: return 7;
		case 0x8ull: 
		case 0x71ull: 
		case 0x62ull: 
		case 0x611ull:
		case 0x53ull: 
		case 0x521ull: 
		case 0x5111ull: 
		case 0x44ull: 
		case 0x431ull: 
		case 0x422ull: 
		case 0x4211ull: 
		case 0x41111ull: 
		case 0x332ull: 
		case 0x3311ull: 
		case 0x3221ull: 
		case 0x32111ull: 
		case 0x311111ull: 
		case 0x2222ull: 
		case 0x22211ull: 
		case 0x221111ull: 
		case 0x2111111ull: 
		case 0x11111111ull: return 8;
		case 0x9ull: 
		case 0x81ull: 
		case 0x72ull: 
		case 0x711ull: 
		case 0x63ull: 
		case 0x621ull: 
		case 0x6111ull: 
		case 0x54ull: 
		case 0x531ull: 
		case 0x522ull: 
		case 0x5211ull: 
		case 0x51111ull: 
		case 0x441ull: 
		case 0x432ull: 
		case 0x4311ull: 
		case 0x4221ull: 
		case 0x42111ull: 
		case 0x411111ull: 
		case 0x333ull: 
		case 0x3321ull: 
		case 0x33111ull: 
		case 0x3222ull: 
		case 0x32211ull: 
		case 0x321111ull: 
		case 0x3111111ull: 
		case 0x22221ull: 
		case 0x222111ull: 
		case 0x2211111ull: 
		case 0x21111111ull: 
		case 0x111111111ull: return 9;
		case 0xAull: 
		case 0x91ull: 
		case 0x82ull: 
		case 0x811ull: 
		case 0x73ull: 
		case 0x721ull: 
		case 0x7111ull: 
		case 0x64ull: 
		case 0x631ull: 
		case 0x622ull: 
		case 0x6211ull: 
		case 0x61111ull: 
		case 0x55ull: 
		case 0x541ull: 
		case 0x532ull: 
		case 0x5311ull: 
		case 0x5221ull: 
		case 0x52111ull: 
		case 0x511111ull: 
		case 0x442ull: 
		case 0x4411ull: 
		case 0x433ull: 
		case 0x4321ull: 
		case 0x43111ull: 
		case 0x4222ull: 
		case 0x42211ull: 
		case 0x421111ull: 
		case 0x4111111ull: 
		case 0x3331ull: 
		case 0x3322ull: 
		case 0x33211ull: 
		case 0x331111ull: 
		case 0x32221ull: 
		case 0x322111ull: 
		case 0x3211111ull: 
		case 0x31111111ull:
		case 0x22222ull: 
		case 0x222211ull: 
		case 0x2221111ull: 
		case 0x22111111ull: 
		case 0x211111111ull: 
		case 0x1111111111ull: return 10;
		}
		*/
		/*switch (patternid)
		{
		default: return patternid >= N_PATTERN_END ? patternid - N_PATTERN_END : 0;
			case 1: return 1;
			case 2:
			case 3: return 2;
			case 4:
			case 5:
			case 6: return 3;
			case 7:
			case 8:
			case 9:
			case 10:
			case 11: return 4;
			case 12:
			case 13:
			case 14:
			case 15:
			case 16:
			case 17:
			case 18: return 5;
			case 19:
			case 20:
			case 21:
			case 22:
			case 23:
			case 24:
			case 25:
			case 26:
			case 27:
			case 28:
			case 29: return 6;
			case 30:
			case 31:
			case 32:
			case 33:
			case 34:
			case 35:
			case 36:
			case 37:
			case 38:
			case 39:
			case 40:
			case 41:
			case 42:
			case 43:
			case 44: return 7;
			case 45:
			case 46:
			case 47:
			case 48:
			case 49:
			case 50:
			case 51:
			case 52:
			case 53:
			case 54:
			case 55:
			case 56:
			case 57:
			case 58:
			case 59:
			case 60:
			case 61:
			case 62:
			case 63:
			case 64:
			case 65:
			case 66: return 8;
			case 67:
			case 68:
			case 69:
			case 70:
			case 71:
			case 72:
			case 73:
			case 74:
			case 75:
			case 76:
			case 77:
			case 78:
			case 79:
			case 80:
			case 81:
			case 82:
			case 83:
			case 84:
			case 85:
			case 86:
			case 87:
			case 88:
			case 89:
			case 90:
			case 91:
			case 92:
			case 93:
			case 94:
			case 95:
			case 96: return 9;
			case 97:
			case 98:
			case 99:
			case 100:
			case 101:
			case 102:
			case 103:
			case 104:
			case 105:
			case 106:
			case 107:
			case 108:
			case 109:
			case 110:
			case 111:
			case 112:
			case 113:
			case 114:
			case 115:
			case 116:
			case 117:
			case 118:
			case 119:
			case 120:
			case 121:
			case 122:
			case 123:
			case 124:
			case 125:
			case 126:
			case 127:
			case 128:
			case 129:
			case 130:
			case 131:
			case 132:
			case 133:
			case 134:
			case 135:
			case 136:
			case 137:
			case 138: return 10;
		}
		*/
	}

	/* Crc32 hash */
	TARGET HASH GENOTYPE::Hash()
	{
		return HashGenotype(GetAlleleArray(), Ploidy());
	}

	/* Heterozygosity in this genotype */
	TARGET double GENOTYPE::HIndex()
	{
		//array is 50% faster than switch
		static double PT_HINDEX[150] = 									//Pattern index to h-index
		{ NAN, 0.000000000000000, 0.000000000000000, 1.000000000000000, 0.000000000000000, 0.666666666666667, 1.000000000000000, 0.000000000000000, 0.500000000000000, 0.666666666666667, 0.833333333333333, 1.000000000000000, 0.000000000000000, 0.400000000000000, 0.600000000000000, 0.700000000000000, 0.800000000000000, 0.900000000000000, 1.000000000000000, 0.000000000000000, 0.333333333333333, 0.533333333333333, 0.600000000000000, 0.600000000000000, 0.733333333333333, 0.800000000000000, 0.800000000000000, 0.866666666666667, 0.933333333333333, 1.000000000000000, 0.000000000000000, 0.285714285714286, 0.476190476190476, 0.523809523809524, 0.571428571428572, 0.666666666666667, 0.714285714285714, 0.714285714285714, 0.761904761904762, 0.809523809523809, 0.857142857142857, 0.857142857142857, 0.904761904761905, 0.952380952380952, 1.000000000000000, 0.000000000000000, 0.250000000000000, 0.428571428571429, 0.464285714285714, 0.535714285714286, 0.607142857142857, 0.642857142857143, 0.571428571428571, 0.678571428571429, 0.714285714285714, 0.750000000000000, 0.785714285714286, 0.750000000000000, 0.785714285714286, 0.821428571428571, 0.857142857142857, 0.892857142857143, 0.857142857142857, 0.892857142857143, 0.928571428571429, 0.964285714285714, 1.000000000000000, 0.000000000000000, 0.222222222222222, 0.388888888888889, 0.416666666666667, 0.500000000000000, 0.555555555555556, 0.583333333333333, 0.555555555555556, 0.638888888888889, 0.666666666666667, 0.694444444444444, 0.722222222222222, 0.666666666666667, 0.722222222222222, 0.750000000000000, 0.777777777777778, 0.805555555555556, 0.833333333333333, 0.750000000000000, 0.805555555555556, 0.833333333333333, 0.833333333333333, 0.861111111111111, 0.888888888888889, 0.916666666666667, 0.888888888888889, 0.916666666666667, 0.944444444444444, 0.972222222222222, 1.000000000000000, 0.000000000000000, 0.200000000000000, 0.355555555555555, 0.377777777777778, 0.466666666666667, 0.511111111111111, 0.533333333333333, 0.533333333333333, 0.600000000000000, 0.622222222222222, 0.644444444444444, 0.666666666666667, 0.555555555555556, 0.644444444444444, 0.688888888888889, 0.711111111111111, 0.733333333333333, 0.755555555555555, 0.777777777777778, 0.711111111111111, 0.733333333333333, 0.733333333333333, 0.777777777777778, 0.800000000000000, 0.800000000000000, 0.822222222222222, 0.844444444444444, 0.866666666666667, 0.800000000000000, 0.822222222222222, 0.844444444444444, 0.866666666666667, 0.866666666666667, 0.888888888888889, 0.911111111111111, 0.933333333333333, 0.888888888888889, 0.911111111111111, 0.933333333333333, 0.955555555555555, 0.977777777777778, 1.000000000000000, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN };

		return PT_HINDEX[patternid];
		/*switch (patternid)
		{
			default: 	return NAN;
			case	1: 	return 0.000000000000000;
			case	2: 	return 0.000000000000000;
			case	3: 	return 1.000000000000000;
			case	4: 	return 0.000000000000000;
			case	5: 	return 0.666666666666667;
			case	6: 	return 1.000000000000000;
			case	7: 	return 0.000000000000000;
			case	8: 	return 0.500000000000000;
			case	9: 	return 0.666666666666667;
			case	10: 	return 0.833333333333333;
			case	11: 	return 1.000000000000000;
			case	12: 	return 0.000000000000000;
			case	13: 	return 0.400000000000000;
			case	14: 	return 0.600000000000000;
			case	15: 	return 0.700000000000000;
			case	16: 	return 0.800000000000000;
			case	17: 	return 0.900000000000000;
			case	18: 	return 1.000000000000000;
			case	19: 	return 0.000000000000000;
			case	20: 	return 0.333333333333333;
			case	21: 	return 0.533333333333333;
			case	22: 	return 0.600000000000000;
			case	23: 	return 0.600000000000000;
			case	24: 	return 0.733333333333333;
			case	25: 	return 0.800000000000000;
			case	26: 	return 0.800000000000000;
			case	27: 	return 0.866666666666667;
			case	28: 	return 0.933333333333333;
			case	29: 	return 1.000000000000000;
			case	30: 	return 0.000000000000000;
			case	31: 	return 0.285714285714286;
			case	32: 	return 0.476190476190476;
			case	33: 	return 0.523809523809524;
			case	34: 	return 0.571428571428572;
			case	35: 	return 0.666666666666667;
			case	36: 	return 0.714285714285714;
			case	37: 	return 0.714285714285714;
			case	38: 	return 0.761904761904762;
			case	39: 	return 0.809523809523809;
			case	40: 	return 0.857142857142857;
			case	41: 	return 0.857142857142857;
			case	42: 	return 0.904761904761905;
			case	43: 	return 0.952380952380952;
			case	44: 	return 1.000000000000000;
			case	45: 	return 0.000000000000000;
			case	46: 	return 0.250000000000000;
			case	47: 	return 0.428571428571429;
			case	48: 	return 0.464285714285714;
			case	49: 	return 0.535714285714286;
			case	50: 	return 0.607142857142857;
			case	51: 	return 0.642857142857143;
			case	52: 	return 0.571428571428571;
			case	53: 	return 0.678571428571429;
			case	54: 	return 0.714285714285714;
			case	55: 	return 0.750000000000000;
			case	56: 	return 0.785714285714286;
			case	57: 	return 0.750000000000000;
			case	58: 	return 0.785714285714286;
			case	59: 	return 0.821428571428571;
			case	60: 	return 0.857142857142857;
			case	61: 	return 0.892857142857143;
			case	62: 	return 0.857142857142857;
			case	63: 	return 0.892857142857143;
			case	64: 	return 0.928571428571429;
			case	65: 	return 0.964285714285714;
			case	66: 	return 1.000000000000000;
			case	67: 	return 0.000000000000000;
			case	68: 	return 0.222222222222222;
			case	69: 	return 0.388888888888889;
			case	70: 	return 0.416666666666667;
			case	71: 	return 0.500000000000000;
			case	72: 	return 0.555555555555556;
			case	73: 	return 0.583333333333333;
			case	74: 	return 0.555555555555556;
			case	75: 	return 0.638888888888889;
			case	76: 	return 0.666666666666667;
			case	77: 	return 0.694444444444444;
			case	78: 	return 0.722222222222222;
			case	79: 	return 0.666666666666667;
			case	80: 	return 0.722222222222222;
			case	81: 	return 0.750000000000000;
			case	82: 	return 0.777777777777778;
			case	83: 	return 0.805555555555556;
			case	84: 	return 0.833333333333333;
			case	85: 	return 0.750000000000000;
			case	86: 	return 0.805555555555556;
			case	87: 	return 0.833333333333333;
			case	88: 	return 0.833333333333333;
			case	89: 	return 0.861111111111111;
			case	90: 	return 0.888888888888889;
			case	91: 	return 0.916666666666667;
			case	92: 	return 0.888888888888889;
			case	93: 	return 0.916666666666667;
			case	94: 	return 0.944444444444444;
			case	95: 	return 0.972222222222222;
			case	96: 	return 1.000000000000000;
			case	97: 	return 0.000000000000000;
			case	98: 	return 0.200000000000000;
			case	99: 	return 0.355555555555555;
			case	100: 	return 0.377777777777778;
			case	101: 	return 0.466666666666667;
			case	102: 	return 0.511111111111111;
			case	103: 	return 0.533333333333333;
			case	104: 	return 0.533333333333333;
			case	105: 	return 0.600000000000000;
			case	106: 	return 0.622222222222222;
			case	107: 	return 0.644444444444444;
			case	108: 	return 0.666666666666667;
			case	109: 	return 0.555555555555556;
			case	110: 	return 0.644444444444444;
			case	111: 	return 0.688888888888889;
			case	112: 	return 0.711111111111111;
			case	113: 	return 0.733333333333333;
			case	114: 	return 0.755555555555555;
			case	115: 	return 0.777777777777778;
			case	116: 	return 0.711111111111111;
			case	117: 	return 0.733333333333333;
			case	118: 	return 0.733333333333333;
			case	119: 	return 0.777777777777778;
			case	120: 	return 0.800000000000000;
			case	121: 	return 0.800000000000000;
			case	122: 	return 0.822222222222222;
			case	123: 	return 0.844444444444444;
			case	124: 	return 0.866666666666667;
			case	125: 	return 0.800000000000000;
			case	126: 	return 0.822222222222222;
			case	127: 	return 0.844444444444444;
			case	128: 	return 0.866666666666667;
			case	129: 	return 0.866666666666667;
			case	130: 	return 0.888888888888889;
			case	131: 	return 0.911111111111111;
			case	132: 	return 0.933333333333333;
			case	133: 	return 0.888888888888889;
			case	134: 	return 0.911111111111111;
			case	135: 	return 0.933333333333333;
			case	136: 	return 0.955555555555555;
			case	137: 	return 0.977777777777778;
			case	138: 	return 1.000000000000000;
		}*/
	}

	/* SS within genotype under IAM model */
	TARGET double GENOTYPE::SS_IAM()
	{
		static double PT_SSIAM[150] = 									//Pattern index to ss smm
		{ NAN, NAN, 0.000000000000000, 0.500000000000000, 0.000000000000000, 0.666666666666667, 1.000000000000000, 0.000000000000000, 0.750000000000000, 1.000000000000000, 1.250000000000000, 1.500000000000000, 0.000000000000000, 0.800000000000000, 1.200000000000000, 1.400000000000000, 1.600000000000000, 1.800000000000000, 2.000000000000000, 0.000000000000000, 0.833333333333332, 1.333333333333330, 1.500000000000000, 1.500000000000000, 1.833333333333330, 2.000000000000000, 2.000000000000000, 2.166666666666670, 2.333333333333330, 2.500000000000000, 0.000000000000000, 0.857142857142858, 1.428571428571430, 1.571428571428570, 1.714285714285720, 2.000000000000000, 2.142857142857140, 2.142857142857140, 2.285714285714290, 2.428571428571430, 2.571428571428570, 2.571428571428570, 2.714285714285710, 2.857142857142860, 3.000000000000000, 0.000000000000000, 0.875000000000000, 1.500000000000000, 1.625000000000000, 1.875000000000000, 2.125000000000000, 2.250000000000000, 2.000000000000000, 2.375000000000000, 2.500000000000000, 2.625000000000000, 2.750000000000000, 2.625000000000000, 2.750000000000000, 2.875000000000000, 3.000000000000000, 3.125000000000000, 3.000000000000000, 3.125000000000000, 3.250000000000000, 3.375000000000000, 3.500000000000000, 0.000000000000000, 0.888888888888888, 1.555555555555560, 1.666666666666670, 2.000000000000000, 2.222222222222220, 2.333333333333330, 2.222222222222220, 2.555555555555560, 2.666666666666670, 2.777777777777780, 2.888888888888890, 2.666666666666670, 2.888888888888890, 3.000000000000000, 3.111111111111110, 3.222222222222220, 3.333333333333330, 3.000000000000000, 3.222222222222220, 3.333333333333330, 3.333333333333330, 3.444444444444440, 3.555555555555560, 3.666666666666670, 3.555555555555560, 3.666666666666670, 3.777777777777780, 3.888888888888890, 4.000000000000000, 0.000000000000000, 0.900000000000000, 1.600000000000000, 1.700000000000000, 2.100000000000000, 2.300000000000000, 2.400000000000000, 2.400000000000000, 2.700000000000000, 2.800000000000000, 2.900000000000000, 3.000000000000000, 2.500000000000000, 2.900000000000000, 3.100000000000000, 3.200000000000000, 3.300000000000000, 3.400000000000000, 3.500000000000000, 3.200000000000000, 3.300000000000000, 3.300000000000000, 3.500000000000000, 3.600000000000000, 3.600000000000000, 3.700000000000000, 3.800000000000000, 3.900000000000000, 3.600000000000000, 3.700000000000000, 3.800000000000000, 3.900000000000000, 3.900000000000000, 4.000000000000000, 4.100000000000000, 4.200000000000000, 4.000000000000000, 4.100000000000000, 4.200000000000000, 4.300000000000000, 4.400000000000000, 4.500000000000000, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN };

		return PT_SSIAM[patternid];
		/*switch (patternid)
		{
		default: return NAN;
		case 2: return 0.0000000000000000;
		case 3: return 0.5000000000000000;
		case 4: return 0.0000000000000000;
		case 5: return 0.6666666666666670;
		case 6: return 1.0000000000000000;
		case 7: return 0.0000000000000000;
		case 8: return 0.7500000000000000;
		case 9: return 1.0000000000000000;
		case 10: return 1.2500000000000000;
		case 11: return 1.5000000000000000;
		case 12: return 0.0000000000000000;
		case 13: return 0.8000000000000000;
		case 14: return 1.2000000000000000;
		case 15: return 1.4000000000000000;
		case 16: return 1.6000000000000000;
		case 17: return 1.8000000000000000;
		case 18: return 2.0000000000000000;
		case 19: return 0.0000000000000000;
		case 20: return 0.8333333333333320;
		case 21: return 1.3333333333333300;
		case 22: return 1.5000000000000000;
		case 23: return 1.5000000000000000;
		case 24: return 1.8333333333333300;
		case 25: return 2.0000000000000000;
		case 26: return 2.0000000000000000;
		case 27: return 2.1666666666666700;
		case 28: return 2.3333333333333300;
		case 29: return 2.5000000000000000;
		case 30: return 0.0000000000000000;
		case 31: return 0.8571428571428580;
		case 32: return 1.4285714285714300;
		case 33: return 1.5714285714285700;
		case 34: return 1.7142857142857200;
		case 35: return 2.0000000000000000;
		case 36: return 2.1428571428571400;
		case 37: return 2.1428571428571400;
		case 38: return 2.2857142857142900;
		case 39: return 2.4285714285714300;
		case 40: return 2.5714285714285700;
		case 41: return 2.5714285714285700;
		case 42: return 2.7142857142857100;
		case 43: return 2.8571428571428600;
		case 44: return 3.0000000000000000;
		case 45: return 0.0000000000000000;
		case 46: return 0.8750000000000000;
		case 47: return 1.5000000000000000;
		case 48: return 1.6250000000000000;
		case 49: return 1.8750000000000000;
		case 50: return 2.1250000000000000;
		case 51: return 2.2500000000000000;
		case 52: return 2.0000000000000000;
		case 53: return 2.3750000000000000;
		case 54: return 2.5000000000000000;
		case 55: return 2.6250000000000000;
		case 56: return 2.7500000000000000;
		case 57: return 2.6250000000000000;
		case 58: return 2.7500000000000000;
		case 59: return 2.8750000000000000;
		case 60: return 3.0000000000000000;
		case 61: return 3.1250000000000000;
		case 62: return 3.0000000000000000;
		case 63: return 3.1250000000000000;
		case 64: return 3.2500000000000000;
		case 65: return 3.3750000000000000;
		case 66: return 3.5000000000000000;
		case 67: return 0.0000000000000000;
		case 68: return 0.8888888888888880;
		case 69: return 1.5555555555555600;
		case 70: return 1.6666666666666700;
		case 71: return 2.0000000000000000;
		case 72: return 2.2222222222222200;
		case 73: return 2.3333333333333300;
		case 74: return 2.2222222222222200;
		case 75: return 2.5555555555555600;
		case 76: return 2.6666666666666700;
		case 77: return 2.7777777777777800;
		case 78: return 2.8888888888888900;
		case 79: return 2.6666666666666700;
		case 80: return 2.8888888888888900;
		case 81: return 3.0000000000000000;
		case 82: return 3.1111111111111100;
		case 83: return 3.2222222222222200;
		case 84: return 3.3333333333333300;
		case 85: return 3.0000000000000000;
		case 86: return 3.2222222222222200;
		case 87: return 3.3333333333333300;
		case 88: return 3.3333333333333300;
		case 89: return 3.4444444444444400;
		case 90: return 3.5555555555555600;
		case 91: return 3.6666666666666700;
		case 92: return 3.5555555555555600;
		case 93: return 3.6666666666666700;
		case 94: return 3.7777777777777800;
		case 95: return 3.8888888888888900;
		case 96: return 4.0000000000000000;
		case 97: return 0.0000000000000000;
		case 98: return 0.9000000000000000;
		case 99: return 1.6000000000000000;
		case 100: return 1.7000000000000000;
		case 101: return 2.1000000000000000;
		case 102: return 2.3000000000000000;
		case 103: return 2.4000000000000000;
		case 104: return 2.4000000000000000;
		case 105: return 2.7000000000000000;
		case 106: return 2.8000000000000000;
		case 107: return 2.9000000000000000;
		case 108: return 3.0000000000000000;
		case 109: return 2.5000000000000000;
		case 110: return 2.9000000000000000;
		case 111: return 3.1000000000000000;
		case 112: return 3.2000000000000000;
		case 113: return 3.3000000000000000;
		case 114: return 3.4000000000000000;
		case 115: return 3.5000000000000000;
		case 116: return 3.2000000000000000;
		case 117: return 3.3000000000000000;
		case 118: return 3.3000000000000000;
		case 119: return 3.5000000000000000;
		case 120: return 3.6000000000000000;
		case 121: return 3.6000000000000000;
		case 122: return 3.7000000000000000;
		case 123: return 3.8000000000000000;
		case 124: return 3.9000000000000000;
		case 125: return 3.6000000000000000;
		case 126: return 3.7000000000000000;
		case 127: return 3.8000000000000000;
		case 128: return 3.9000000000000000;
		case 129: return 3.9000000000000000;
		case 130: return 4.0000000000000000;
		case 131: return 4.1000000000000000;
		case 132: return 4.2000000000000000;
		case 133: return 4.0000000000000000;
		case 134: return 4.1000000000000000;
		case 135: return 4.2000000000000000;
		case 136: return 4.3000000000000000;
		case 137: return 4.4000000000000000;
		case 138: return 4.5000000000000000;
		}*/
		//return HIndex() * ((Ploidy() - 1) * 0.5);
	}

	/* SS within genotype under SMM model */
	TARGET double GENOTYPE::SS_SMM(ushort* alen)
	{
		double smm = 0;
		ushort* als = GetAlleleArray();
		int ploidy = Ploidy();
		for (int i1 = 0; i1 < ploidy; ++i1)
			for (int i2 = 0; i2 < i1; ++i2)
				smm += ((int)alen[als[i1]] - (int)alen[als[i2]]) * ((int)alen[als[i1]] - (int)alen[als[i2]]);
		return smm / ploidy;
	}

	/* The multinomial coefficient for HWE/RCS genotype frequency */
	TARGET double GENOTYPE::HWECoef()
	{
		static double PT_HWECOEF[150] = 								//Pattern index to hwe coef
		{ NAN, 1, 1, 2, 1, 3, 6, 1, 4, 6, 12, 24, 1, 5, 10, 20, 30, 60, 120, 1, 6, 15, 30, 20, 60, 120, 90, 180, 360, 720, 1, 7, 21, 42, 35, 105, 210, 140, 210, 420, 840, 630, 1260, 2520, 5040, 1, 8, 28, 56, 56, 168, 336, 70, 280, 420, 840, 1680, 560, 1120, 1680, 3360, 6720, 2520, 5040, 10080, 20160, 40320, 1, 9, 36, 72, 84, 252, 504, 126, 504, 756, 1512, 3024, 630, 1260, 2520, 3780, 7560, 15120, 1680, 5040, 10080, 7560, 15120, 30240, 60480, 22680, 45360, 90720, 181440, 362880, 1, 10, 45, 90, 120, 360, 720, 210, 840, 1260, 2520, 5040, 252, 1260, 2520, 5040, 7560, 15120, 30240, 3150, 6300, 4200, 12600, 25200, 18900, 37800, 75600, 151200, 16800, 25200, 50400, 100800, 75600, 151200, 302400, 604800, 113400, 226800, 453600, 907200, 1814400, 3628800, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN };

		return PT_HWECOEF[patternid];
		/*
		switch (pattern)
		{
		default: return NAN;
		case 0x1ull: return 1;
		case 0x2ull: return 1;
		case 0x11ull: return 2;
		case 0x3ull: return 1;
		case 0x21ull: return 3;
		case 0x111ull: return 6;
		case 0x4ull: return 1;
		case 0x31ull: return 4;
		case 0x22ull: return 6;
		case 0x211ull: return 12;
		case 0x1111ull: return 24;
		case 0x5ull: return 1;
		case 0x41ull: return 5;
		case 0x32ull: return 10;
		case 0x311ull: return 20;
		case 0x221ull: return 30;
		case 0x2111ull: return 60;
		case 0x11111ull: return 120;
		case 0x6ull: return 1;
		case 0x51ull: return 6;
		case 0x42ull: return 15;
		case 0x411ull: return 30;
		case 0x33ull: return 20;
		case 0x321ull: return 60;
		case 0x3111ull: return 120;
		case 0x222ull: return 90;
		case 0x2211ull: return 180;
		case 0x21111ull: return 360;
		case 0x111111ull: return 720;
		case 0x7ull: return 1;
		case 0x61ull: return 7;
		case 0x52ull: return 21;
		case 0x511ull: return 42;
		case 0x43ull: return 35;
		case 0x421ull: return 105;
		case 0x4111ull: return 210;
		case 0x331ull: return 140;
		case 0x322ull: return 210;
		case 0x3211ull: return 420;
		case 0x31111ull: return 840;
		case 0x2221ull: return 630;
		case 0x22111ull: return 1260;
		case 0x211111ull: return 2520;
		case 0x1111111ull: return 5040;
		case 0x8ull: return 1;
		case 0x71ull: return 8;
		case 0x62ull: return 28;
		case 0x611ull: return 56;
		case 0x53ull: return 56;
		case 0x521ull: return 168;
		case 0x5111ull: return 336;
		case 0x44ull: return 70;
		case 0x431ull: return 280;
		case 0x422ull: return 420;
		case 0x4211ull: return 840;
		case 0x41111ull: return 1680;
		case 0x332ull: return 560;
		case 0x3311ull: return 1120;
		case 0x3221ull: return 1680;
		case 0x32111ull: return 3360;
		case 0x311111ull: return 6720;
		case 0x2222ull: return 2520;
		case 0x22211ull: return 5040;
		case 0x221111ull: return 10080;
		case 0x2111111ull: return 20160;
		case 0x11111111ull: return 40320;
		case 0x9ull: return 1;
		case 0x81ull: return 9;
		case 0x72ull: return 36;
		case 0x711ull: return 72;
		case 0x63ull: return 84;
		case 0x621ull: return 252;
		case 0x6111ull: return 504;
		case 0x54ull: return 126;
		case 0x531ull: return 504;
		case 0x522ull: return 756;
		case 0x5211ull: return 1512;
		case 0x51111ull: return 3024;
		case 0x441ull: return 630;
		case 0x432ull: return 1260;
		case 0x4311ull: return 2520;
		case 0x4221ull: return 3780;
		case 0x42111ull: return 7560;
		case 0x411111ull: return 15120;
		case 0x333ull: return 1680;
		case 0x3321ull: return 5040;
		case 0x33111ull: return 10080;
		case 0x3222ull: return 7560;
		case 0x32211ull: return 15120;
		case 0x321111ull: return 30240;
		case 0x3111111ull: return 60480;
		case 0x22221ull: return 22680;
		case 0x222111ull: return 45360;
		case 0x2211111ull: return 90720;
		case 0x21111111ull: return 181440;
		case 0x111111111ull: return 362880;
		case 0xAull: return 1;
		case 0x91ull: return 10;
		case 0x82ull: return 45;
		case 0x811ull: return 90;
		case 0x73ull: return 120;
		case 0x721ull: return 360;
		case 0x7111ull: return 720;
		case 0x64ull: return 210;
		case 0x631ull: return 840;
		case 0x622ull: return 1260;
		case 0x6211ull: return 2520;
		case 0x61111ull: return 5040;
		case 0x55ull: return 252;
		case 0x541ull: return 1260;
		case 0x532ull: return 2520;
		case 0x5311ull: return 5040;
		case 0x5221ull: return 7560;
		case 0x52111ull: return 15120;
		case 0x511111ull: return 30240;
		case 0x442ull: return 3150;
		case 0x4411ull: return 6300;
		case 0x433ull: return 4200;
		case 0x4321ull: return 12600;
		case 0x43111ull: return 25200;
		case 0x4222ull: return 18900;
		case 0x42211ull: return 37800;
		case 0x421111ull: return 75600;
		case 0x4111111ull: return 151200;
		case 0x3331ull: return 16800;
		case 0x3322ull: return 25200;
		case 0x33211ull: return 50400;
		case 0x331111ull: return 100800;
		case 0x32221ull: return 75600;
		case 0x322111ull: return 151200;
		case 0x3211111ull: return 302400;
		case 0x31111111ull: return 604800;
		case 0x22222ull: return 113400;
		case 0x222211ull: return 226800;
		case 0x2221111ull: return 453600;
		case 0x22111111ull: return 907200;
		case 0x211111111ull: return 1814400;
		case 0x1111111111ull: return 3628800;
		}
		*/
		/*switch (patternid)
		{
			default: 	return NAN;
			case	1: 	return 1;
			case	2: 	return 1;
			case	3: 	return 2;
			case	4: 	return 1;
			case	5: 	return 3;
			case	6: 	return 6;
			case	7: 	return 1;
			case	8: 	return 4;
			case	9: 	return 6;
			case	10: 	return 12;
			case	11: 	return 24;
			case	12: 	return 1;
			case	13: 	return 5;
			case	14: 	return 10;
			case	15: 	return 20;
			case	16: 	return 30;
			case	17: 	return 60;
			case	18: 	return 120;
			case	19: 	return 1;
			case	20: 	return 6;
			case	21: 	return 15;
			case	22: 	return 30;
			case	23: 	return 20;
			case	24: 	return 60;
			case	25: 	return 120;
			case	26: 	return 90;
			case	27: 	return 180;
			case	28: 	return 360;
			case	29: 	return 720;
			case	30: 	return 1;
			case	31: 	return 7;
			case	32: 	return 21;
			case	33: 	return 42;
			case	34: 	return 35;
			case	35: 	return 105;
			case	36: 	return 210;
			case	37: 	return 140;
			case	38: 	return 210;
			case	39: 	return 420;
			case	40: 	return 840;
			case	41: 	return 630;
			case	42: 	return 1260;
			case	43: 	return 2520;
			case	44: 	return 5040;
			case	45: 	return 1;
			case	46: 	return 8;
			case	47: 	return 28;
			case	48: 	return 56;
			case	49: 	return 56;
			case	50: 	return 168;
			case	51: 	return 336;
			case	52: 	return 70;
			case	53: 	return 280;
			case	54: 	return 420;
			case	55: 	return 840;
			case	56: 	return 1680;
			case	57: 	return 560;
			case	58: 	return 1120;
			case	59: 	return 1680;
			case	60: 	return 3360;
			case	61: 	return 6720;
			case	62: 	return 2520;
			case	63: 	return 5040;
			case	64: 	return 10080;
			case	65: 	return 20160;
			case	66: 	return 40320;
			case	67: 	return 1;
			case	68: 	return 9;
			case	69: 	return 36;
			case	70: 	return 72;
			case	71: 	return 84;
			case	72: 	return 252;
			case	73: 	return 504;
			case	74: 	return 126;
			case	75: 	return 504;
			case	76: 	return 756;
			case	77: 	return 1512;
			case	78: 	return 3024;
			case	79: 	return 630;
			case	80: 	return 1260;
			case	81: 	return 2520;
			case	82: 	return 3780;
			case	83: 	return 7560;
			case	84: 	return 15120;
			case	85: 	return 1680;
			case	86: 	return 5040;
			case	87: 	return 10080;
			case	88: 	return 7560;
			case	89: 	return 15120;
			case	90: 	return 30240;
			case	91: 	return 60480;
			case	92: 	return 22680;
			case	93: 	return 45360;
			case	94: 	return 90720;
			case	95: 	return 181440;
			case	96: 	return 362880;
			case	97: 	return 1;
			case	98: 	return 10;
			case	99: 	return 45;
			case	100: 	return 90;
			case	101: 	return 120;
			case	102: 	return 360;
			case	103: 	return 720;
			case	104: 	return 210;
			case	105: 	return 840;
			case	106: 	return 1260;
			case	107: 	return 2520;
			case	108: 	return 5040;
			case	109: 	return 252;
			case	110: 	return 1260;
			case	111: 	return 2520;
			case	112: 	return 5040;
			case	113: 	return 7560;
			case	114: 	return 15120;
			case	115: 	return 30240;
			case	116: 	return 3150;
			case	117: 	return 6300;
			case	118: 	return 4200;
			case	119: 	return 12600;
			case	120: 	return 25200;
			case	121: 	return 18900;
			case	122: 	return 37800;
			case	123: 	return 75600;
			case	124: 	return 151200;
			case	125: 	return 16800;
			case	126: 	return 25200;
			case	127: 	return 50400;
			case	128: 	return 100800;
			case	129: 	return 75600;
			case	130: 	return 151200;
			case	131: 	return 302400;
			case	132: 	return 604800;
			case	133: 	return 113400;
			case	134: 	return 226800;
			case	135: 	return 453600;
			case	136: 	return 907200;
			case	137: 	return 1814400;
			case	138: 	return 3628800;
		}*/
	}

	/* Genotypic frequency for zygotes under inbreeding */
	TARGET double GENOTYPE::GFZ(int* allele_count, int sum, double f)
	{
		// allele_count[i] / sum is the allele frequency

		int ploidy = Ploidy(), nalleles = Nalleles();
		if (nalleles == 0) return 1;

		double re = HWECoef(), re2 = 1, F = (1.0 / f - 1);
		if (F < 1e-20) F = 1e-20;  if (F > 1e20) F = 1e20;
		double FP = F / sum;
		uint64 ap = GetPattern();
		ushort* als = GetAlleleArray();

		for (int a = nalleles - 1; a >= 0; --a)
		{
			for (int j = 0; j < (ap & 0xF); ++j)
				re *= allele_count[als[ploidy + a]] * FP + j;
			ap >>= 4;
		}
		for (int j = 0; j < ploidy; ++j)
			re2 *= F + j;
		return re / re2;
	}

	/* Genotypic frequency for zygotes under specific double-reduction model */
	TARGET double GENOTYPE::GFZ(int DR_MODE, double* f)
	{
		int ploidy = Ploidy(), nalleles = Nalleles();
		if (nalleles == 0) return 1;

		ushort* a = GetAlleleArray() + ploidy;
		if (DR_MODE == 0 || ploidy <= 2 || ploidy % 2 == 1 || ploidy > N_MAX_PLOIDY)
		{
			double re = HWECoef();
			uint64 ap = GetPattern();
			for (int ai = nalleles - 1; ai >= 0; --ai)
			{
				re *= IntegerPower(f[a[ai]], (int)(ap & 0xF));
				ap >>= 4;
			}
			return re;
		}

		double* alpha = &ALPHA[DR_MODE][ploidy][0];
		switch (patternid)
		{
			case	7: 		return GFZ4_iiii(alpha[1], f[a[0]]);
			case	8: 		return GFZ4_iiij(alpha[1], f[a[0]], f[a[1]]);
			case	9: 		return GFZ4_iijj(alpha[1], f[a[0]], f[a[1]]);
			case	10: 	return GFZ4_iijk(alpha[1], f[a[0]], f[a[1]], f[a[2]]);
			case	11: 	return GFZ4_ijkl(alpha[1], f[a[0]], f[a[1]], f[a[2]], f[a[3]]);

			case	19: 	return GFZ6_iiiiii(alpha[1], f[a[0]]);
			case	20: 	return GFZ6_iiiiij(alpha[1], f[a[0]], f[a[1]]);
			case	21: 	return GFZ6_iiiijj(alpha[1], f[a[0]], f[a[1]]);
			case	22: 	return GFZ6_iiiijk(alpha[1], f[a[0]], f[a[1]], f[a[2]]);
			case	23: 	return GFZ6_iiijjj(alpha[1], f[a[0]], f[a[1]]);
			case	24: 	return GFZ6_iiijjk(alpha[1], f[a[0]], f[a[1]], f[a[2]]);
			case	25: 	return GFZ6_iiijkl(alpha[1], f[a[0]], f[a[1]], f[a[2]], f[a[3]]);
			case	26: 	return GFZ6_iijjkk(alpha[1], f[a[0]], f[a[1]], f[a[2]]);
			case	27: 	return GFZ6_iijjkl(alpha[1], f[a[0]], f[a[1]], f[a[2]], f[a[3]]);
			case	28: 	return GFZ6_iijklm(alpha[1], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]]);
			case	29: 	return GFZ6_ijklmn(alpha[1], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]], f[a[5]]);

			case	45: 	return GFZ8_iiiiiiii(alpha[1], alpha[2], f[a[0]]);
			case	46: 	return GFZ8_iiiiiiij(alpha[1], alpha[2], f[a[0]], f[a[1]]);
			case	47: 	return GFZ8_iiiiiijj(alpha[1], alpha[2], f[a[0]], f[a[1]]);
			case	48: 	return GFZ8_iiiiiijk(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]]);
			case	49: 	return GFZ8_iiiiijjj(alpha[1], alpha[2], f[a[0]], f[a[1]]);
			case	50: 	return GFZ8_iiiiijjk(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]]);
			case	51: 	return GFZ8_iiiiijkl(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]]);
			case	52: 	return GFZ8_iiiijjjj(alpha[1], alpha[2], f[a[0]], f[a[1]]);
			case	53: 	return GFZ8_iiiijjjk(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]]);
			case	54: 	return GFZ8_iiiijjkk(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]]);
			case	55: 	return GFZ8_iiiijjkl(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]]);
			case	56: 	return GFZ8_iiiijklm(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]]);
			case	57: 	return GFZ8_iiijjjkk(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]]);
			case	58: 	return GFZ8_iiijjjkl(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]]);
			case	59: 	return GFZ8_iiijjkkl(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]]);
			case	60: 	return GFZ8_iiijjklm(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]]);
			case	61: 	return GFZ8_iiijklmn(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]], f[a[5]]);
			case	62: 	return GFZ8_iijjkkll(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]]);
			case	63: 	return GFZ8_iijjkklm(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]]);
			case	64: 	return GFZ8_iijjklmn(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]], f[a[5]]);
			case	65: 	return GFZ8_iijklmno(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]], f[a[5]], f[a[6]]);
			case	66: 	return GFZ8_ijklmnop(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]], f[a[5]], f[a[6]], f[a[7]]);

			case	97: 	return GFZ10_iiiiiiiiii(alpha[1], alpha[2], f[a[0]]);
			case	98: 	return GFZ10_iiiiiiiiij(alpha[1], alpha[2], f[a[0]], f[a[1]]);
			case	99: 	return GFZ10_iiiiiiiijj(alpha[1], alpha[2], f[a[0]], f[a[1]]);
			case	100: 	return GFZ10_iiiiiiiijk(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]]);
			case	101: 	return GFZ10_iiiiiiijjj(alpha[1], alpha[2], f[a[0]], f[a[1]]);
			case	102: 	return GFZ10_iiiiiiijjk(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]]);
			case	103: 	return GFZ10_iiiiiiijkl(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]]);
			case	104: 	return GFZ10_iiiiiijjjj(alpha[1], alpha[2], f[a[0]], f[a[1]]);
			case	105: 	return GFZ10_iiiiiijjjk(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]]);
			case	106: 	return GFZ10_iiiiiijjkk(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]]);
			case	107: 	return GFZ10_iiiiiijjkl(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]]);
			case	108: 	return GFZ10_iiiiiijklm(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]]);
			case	109: 	return GFZ10_iiiiijjjjj(alpha[1], alpha[2], f[a[0]], f[a[1]]);
			case	110: 	return GFZ10_iiiiijjjjk(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]]);
			case	111: 	return GFZ10_iiiiijjjkk(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]]);
			case	112: 	return GFZ10_iiiiijjjkl(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]]);
			case	113: 	return GFZ10_iiiiijjkkl(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]]);
			case	114: 	return GFZ10_iiiiijjklm(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]]);
			case	115: 	return GFZ10_iiiiijklmn(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]], f[a[5]]);
			case	116: 	return GFZ10_iiiijjjjkk(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]]);
			case	117: 	return GFZ10_iiiijjjjkl(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]]);
			case	118: 	return GFZ10_iiiijjjkkk(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]]);
			case	119: 	return GFZ10_iiiijjjkkl(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]]);
			case	120: 	return GFZ10_iiiijjjklm(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]]);
			case	121: 	return GFZ10_iiiijjkkll(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]]);
			case	122: 	return GFZ10_iiiijjkklm(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]]);
			case	123: 	return GFZ10_iiiijjklmn(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]], f[a[5]]);
			case	124: 	return GFZ10_iiiijklmno(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]], f[a[5]], f[a[6]]);
			case	125: 	return GFZ10_iiijjjkkkl(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]]);
			case	126: 	return GFZ10_iiijjjkkll(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]]);
			case	127: 	return GFZ10_iiijjjkklm(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]]);
			case	128: 	return GFZ10_iiijjjklmn(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]], f[a[5]]);
			case	129: 	return GFZ10_iiijjkkllm(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]]);
			case	130: 	return GFZ10_iiijjkklmn(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]], f[a[5]]);
			case	131: 	return GFZ10_iiijjklmno(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]], f[a[5]], f[a[6]]);
			case	132: 	return GFZ10_iiijklmnop(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]], f[a[5]], f[a[6]], f[a[7]]);
			case	133: 	return GFZ10_iijjkkllmm(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]]);
			case	134: 	return GFZ10_iijjkkllmn(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]], f[a[5]]);
			case	135: 	return GFZ10_iijjkklmno(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]], f[a[5]], f[a[6]]);
			case	136: 	return GFZ10_iijjklmnop(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]], f[a[5]], f[a[6]], f[a[7]]);
			case	137: 	return GFZ10_iijklmnopq(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]], f[a[5]], f[a[6]], f[a[7]], f[a[8]]);
			case	138: 	return GFZ10_ijklmnopqr(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]], f[a[5]], f[a[6]], f[a[7]], f[a[8]], f[a[9]]);

			default	:	return -1;
		}
	}

	/* Genotypic frequency for zygotes under specific double-reduction model */
	TARGET double GENOTYPE::GFZ(int DR_MODE, float* f)
	{
		int ploidy = Ploidy(), nalleles = Nalleles();
		if (nalleles == 0) return 1;

		ushort* a = GetAlleleArray() + ploidy;
		if (DR_MODE == 0 || ploidy <= 2 || ploidy % 2 == 1 || ploidy > N_MAX_PLOIDY)
		{
			double re = HWECoef();
			uint64 ap = GetPattern();
			for (int ai = nalleles - 1; ai >= 0; --ai)
			{
				re *= IntegerPower((double)f[a[ai]], (int)(ap & 0xF));
				ap >>= 4;
			}
			return re;
		}

		double* alpha = &ALPHA[DR_MODE][ploidy][0];
		switch (patternid)
		{
			case	7: 		return GFZ4_iiii(alpha[1], (double)f[a[0]]);
			case	8: 		return GFZ4_iiij(alpha[1], (double)f[a[0]], (double)f[a[1]]);
			case	9: 		return GFZ4_iijj(alpha[1], (double)f[a[0]], (double)f[a[1]]);
			case	10: 	return GFZ4_iijk(alpha[1], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]]);
			case	11: 	return GFZ4_ijkl(alpha[1], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]]);

			case	19: 	return GFZ6_iiiiii(alpha[1], (double)f[a[0]]);
			case	20: 	return GFZ6_iiiiij(alpha[1], (double)f[a[0]], (double)f[a[1]]);
			case	21: 	return GFZ6_iiiijj(alpha[1], (double)f[a[0]], (double)f[a[1]]);
			case	22: 	return GFZ6_iiiijk(alpha[1], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]]);
			case	23: 	return GFZ6_iiijjj(alpha[1], (double)f[a[0]], (double)f[a[1]]);
			case	24: 	return GFZ6_iiijjk(alpha[1], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]]);
			case	25: 	return GFZ6_iiijkl(alpha[1], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]]);
			case	26: 	return GFZ6_iijjkk(alpha[1], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]]);
			case	27: 	return GFZ6_iijjkl(alpha[1], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]]);
			case	28: 	return GFZ6_iijklm(alpha[1], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]], (double)f[a[4]]);
			case	29: 	return GFZ6_ijklmn(alpha[1], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]], (double)f[a[4]], (double)f[a[5]]);

			case	45: 	return GFZ8_iiiiiiii(alpha[1], alpha[2], (double)f[a[0]]);
			case	46: 	return GFZ8_iiiiiiij(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]]);
			case	47: 	return GFZ8_iiiiiijj(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]]);
			case	48: 	return GFZ8_iiiiiijk(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]]);
			case	49: 	return GFZ8_iiiiijjj(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]]);
			case	50: 	return GFZ8_iiiiijjk(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]]);
			case	51: 	return GFZ8_iiiiijkl(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]]);
			case	52: 	return GFZ8_iiiijjjj(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]]);
			case	53: 	return GFZ8_iiiijjjk(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]]);
			case	54: 	return GFZ8_iiiijjkk(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]]);
			case	55: 	return GFZ8_iiiijjkl(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]]);
			case	56: 	return GFZ8_iiiijklm(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]], (double)f[a[4]]);
			case	57: 	return GFZ8_iiijjjkk(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]]);
			case	58: 	return GFZ8_iiijjjkl(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]]);
			case	59: 	return GFZ8_iiijjkkl(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]]);
			case	60: 	return GFZ8_iiijjklm(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]], (double)f[a[4]]);
			case	61: 	return GFZ8_iiijklmn(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]], (double)f[a[4]], (double)f[a[5]]);
			case	62: 	return GFZ8_iijjkkll(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]]);
			case	63: 	return GFZ8_iijjkklm(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]], (double)f[a[4]]);
			case	64: 	return GFZ8_iijjklmn(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]], (double)f[a[4]], (double)f[a[5]]);
			case	65: 	return GFZ8_iijklmno(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]], (double)f[a[4]], (double)f[a[5]], (double)f[a[6]]);
			case	66: 	return GFZ8_ijklmnop(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]], (double)f[a[4]], (double)f[a[5]], (double)f[a[6]], (double)f[a[7]]);

			case	97: 	return GFZ10_iiiiiiiiii(alpha[1], alpha[2], (double)f[a[0]]);
			case	98: 	return GFZ10_iiiiiiiiij(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]]);
			case	99: 	return GFZ10_iiiiiiiijj(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]]);
			case	100: 	return GFZ10_iiiiiiiijk(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]]);
			case	101: 	return GFZ10_iiiiiiijjj(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]]);
			case	102: 	return GFZ10_iiiiiiijjk(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]]);
			case	103: 	return GFZ10_iiiiiiijkl(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]]);
			case	104: 	return GFZ10_iiiiiijjjj(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]]);
			case	105: 	return GFZ10_iiiiiijjjk(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]]);
			case	106: 	return GFZ10_iiiiiijjkk(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]]);
			case	107: 	return GFZ10_iiiiiijjkl(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]]);
			case	108: 	return GFZ10_iiiiiijklm(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]], (double)f[a[4]]);
			case	109: 	return GFZ10_iiiiijjjjj(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]]);
			case	110: 	return GFZ10_iiiiijjjjk(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]]);
			case	111: 	return GFZ10_iiiiijjjkk(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]]);
			case	112: 	return GFZ10_iiiiijjjkl(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]]);
			case	113: 	return GFZ10_iiiiijjkkl(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]]);
			case	114: 	return GFZ10_iiiiijjklm(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]], (double)f[a[4]]);
			case	115: 	return GFZ10_iiiiijklmn(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]], (double)f[a[4]], (double)f[a[5]]);
			case	116: 	return GFZ10_iiiijjjjkk(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]]);
			case	117: 	return GFZ10_iiiijjjjkl(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]]);
			case	118: 	return GFZ10_iiiijjjkkk(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]]);
			case	119: 	return GFZ10_iiiijjjkkl(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]]);
			case	120: 	return GFZ10_iiiijjjklm(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]], (double)f[a[4]]);
			case	121: 	return GFZ10_iiiijjkkll(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]]);
			case	122: 	return GFZ10_iiiijjkklm(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]], (double)f[a[4]]);
			case	123: 	return GFZ10_iiiijjklmn(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]], (double)f[a[4]], (double)f[a[5]]);
			case	124: 	return GFZ10_iiiijklmno(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]], (double)f[a[4]], (double)f[a[5]], (double)f[a[6]]);
			case	125: 	return GFZ10_iiijjjkkkl(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]]);
			case	126: 	return GFZ10_iiijjjkkll(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]]);
			case	127: 	return GFZ10_iiijjjkklm(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]], (double)f[a[4]]);
			case	128: 	return GFZ10_iiijjjklmn(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]], (double)f[a[4]], (double)f[a[5]]);
			case	129: 	return GFZ10_iiijjkkllm(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]], (double)f[a[4]]);
			case	130: 	return GFZ10_iiijjkklmn(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]], (double)f[a[4]], (double)f[a[5]]);
			case	131: 	return GFZ10_iiijjklmno(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]], (double)f[a[4]], (double)f[a[5]], (double)f[a[6]]);
			case	132: 	return GFZ10_iiijklmnop(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]], (double)f[a[4]], (double)f[a[5]], (double)f[a[6]], (double)f[a[7]]);
			case	133: 	return GFZ10_iijjkkllmm(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]], (double)f[a[4]]);
			case	134: 	return GFZ10_iijjkkllmn(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]], (double)f[a[4]], (double)f[a[5]]);
			case	135: 	return GFZ10_iijjkklmno(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]], (double)f[a[4]], (double)f[a[5]], (double)f[a[6]]);
			case	136: 	return GFZ10_iijjklmnop(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]], (double)f[a[4]], (double)f[a[5]], (double)f[a[6]], (double)f[a[7]]);
			case	137: 	return GFZ10_iijklmnopq(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]], (double)f[a[4]], (double)f[a[5]], (double)f[a[6]], (double)f[a[7]], (double)f[a[8]]);
			case	138: 	return GFZ10_ijklmnopqr(alpha[1], alpha[2], (double)f[a[0]], (double)f[a[1]], (double)f[a[2]], (double)f[a[3]], (double)f[a[4]], (double)f[a[5]], (double)f[a[6]], (double)f[a[7]], (double)f[a[8]], (double)f[a[9]]);

			default	:	return -1;
		}
	}

	/* Number of copies of target allele */
	TARGET int GENOTYPE::GetAlleleCount(int a)
	{
		int ploidy = Ploidy(), nalleles = Nalleles();
		if (nalleles == 0)
			return 0;

		int count = 0;
		ushort* als = GetAlleleArray(), usa = (ushort)a;
		for (int i = 0; i < ploidy; ++i)
			count += (als[i] == usa);
		return count;
	}

	/* Frequency of target allele in this genotype */
	template<typename REAL>
	TARGET REAL GENOTYPE::GetFreq(int a)
	{
		static REAL PT_FREQ[171] =
		{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, -1, 1.0, 0.50, 0.333333333333333, 0.250, 0.20, 0.166666666666667, 0.142857142857143, 0.1250, 0.111111111111111, 0.10, -1, -1, -1, -1, -1, -1, -1, 1.0, 0.666666666666667, 0.50, 0.40, 0.333333333333333, 0.285714285714286, 0.250, 0.222222222222222, 0.20, -1, -1, -1, -1, -1, -1, -1, -1, 1.0, 0.750, 0.60, 0.50, 0.428571428571429, 0.3750, 0.333333333333333, 0.30, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1.0, 0.80, 0.666666666666667, 0.571428571428571, 0.50, 0.444444444444444, 0.40, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1.0, 0.833333333333333, 0.714285714285714, 0.6250, 0.555555555555556, 0.50, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1.0, 0.857142857142857, 0.750, 0.666666666666667, 0.60, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1.0, 0.8750, 0.777777777777778, 0.70, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1.0, 0.888888888888889, 0.80, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1.0, 0.90, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1.0 };

		int ploidy = Ploidy(), nalleles = Nalleles();
		if (nalleles == 0) return NAN;

		int count = 0;
		ushort* als = GetAlleleArray(), usa = (ushort)a;;
		for (int i = 0; i < ploidy; ++i)
			count += (als[i] == usa);

		return PT_FREQ[(count << 4) | ploidy];
		/*switch ((count << 4) | ploidy)
		{
		default: return NAN;
		case 0x00: return 0.000000000000000;
		case 0x01: return 0.000000000000000;
		case 0x02: return 0.000000000000000;
		case 0x03: return 0.000000000000000;
		case 0x04: return 0.000000000000000;
		case 0x05: return 0.000000000000000;
		case 0x06: return 0.000000000000000;
		case 0x07: return 0.000000000000000;
		case 0x08: return 0.000000000000000;
		case 0x09: return 0.000000000000000;
		case 0x0A: return 0.000000000000000;
		case 0x11: return 1.000000000000000;
		case 0x12: return 0.500000000000000;
		case 0x13: return 0.333333333333333;
		case 0x14: return 0.250000000000000;
		case 0x15: return 0.200000000000000;
		case 0x16: return 0.166666666666667;
		case 0x17: return 0.142857142857143;
		case 0x18: return 0.125000000000000;
		case 0x19: return 0.111111111111111;
		case 0x1A: return 0.100000000000000;
		case 0x22: return 1.000000000000000;
		case 0x23: return 0.666666666666667;
		case 0x24: return 0.500000000000000;
		case 0x25: return 0.400000000000000;
		case 0x26: return 0.333333333333333;
		case 0x27: return 0.285714285714286;
		case 0x28: return 0.250000000000000;
		case 0x29: return 0.222222222222222;
		case 0x2A: return 0.200000000000000;
		case 0x33: return 1.000000000000000;
		case 0x34: return 0.750000000000000;
		case 0x35: return 0.600000000000000;
		case 0x36: return 0.500000000000000;
		case 0x37: return 0.428571428571429;
		case 0x38: return 0.375000000000000;
		case 0x39: return 0.333333333333333;
		case 0x3A: return 0.300000000000000;
		case 0x44: return 1.000000000000000;
		case 0x45: return 0.800000000000000;
		case 0x46: return 0.666666666666667;
		case 0x47: return 0.571428571428571;
		case 0x48: return 0.500000000000000;
		case 0x49: return 0.444444444444444;
		case 0x4A: return 0.400000000000000;
		case 0x55: return 1.000000000000000;
		case 0x56: return 0.833333333333333;
		case 0x57: return 0.714285714285714;
		case 0x58: return 0.625000000000000;
		case 0x59: return 0.555555555555556;
		case 0x5A: return 0.500000000000000;
		case 0x66: return 1.000000000000000;
		case 0x67: return 0.857142857142857;
		case 0x68: return 0.750000000000000;
		case 0x69: return 0.666666666666667;
		case 0x6A: return 0.600000000000000;
		case 0x77: return 1.000000000000000;
		case 0x78: return 0.875000000000000;
		case 0x79: return 0.777777777777778;
		case 0x7A: return 0.700000000000000;
		case 0x88: return 1.000000000000000;
		case 0x89: return 0.888888888888889;
		case 0x8A: return 0.800000000000000;
		case 0x99: return 1.000000000000000;
		case 0x9A: return 0.900000000000000;
		case 0xAA: return 1.000000000000000;
		}*/
		//return count / (double)ploidy;
	}

	/* Frequencies of all alleles in this genotype */
	template<typename REAL>
	TARGET void GENOTYPE::GetFreq(REAL* p, int k2)
	{
		int ploidy = Ploidy(), nalleles = Nalleles();
		SetZero(p, k2);
		if (nalleles == 0) return;

		REAL f = 1.0 / ploidy;
		ushort* als = GetAlleleArray();
		for (int i = 0; i < ploidy; ++i)
			p[als[i]] += f;
	}

	/* Obtain Genepop genotype string */
	TARGET char* GENOTYPE::GetGenepopStr()
	{
		int ploidy = Ploidy(), nalleles = Nalleles();
		if (convert_mode_val == 1 && ploidy != 2) Exit("\nError: Genepop supports diploids only.\n");
		char* str = NULL;

		if (convert_mode_val == 1 || convert_mode_val == 2) //disable or truncate
		{
			conversion_memory2[threadid].Alloc(str, MAX_ALLELE_BYTE * 2 + 2);

			ushort* als = GetAlleleArray();
			if (nalleles == 0 || convert_mode_val != 1 && ploidy == 1)
				sprintf(str, " 000000");
			else
				sprintf(str, " %03d%03d", als[0] + 1, als[1] + 1);
		}

		return str;
	}

	/* Obtain Spagedi genotype string */
	TARGET char* GENOTYPE::GetSpagediStr()
	{
		int ploidy = Ploidy(), nalleles = Nalleles(); 
		char* str = NULL;

		if (convert_mode_val == 1 || convert_mode_val == 2) //disable or truncate
		{
			ushort* als = GetAlleleArray();
			if (convert_mode_val == 1)
			{
				conversion_memory2[threadid].Alloc(str, MAX_ALLELE_BYTE * ploidy + 2);
				char* re2 = str;
				sprintf(re2++, "\t");
				for (int i = 0; i < ploidy; ++i)
				{
					if (nalleles == 0)
						sprintf(re2, "000");
					else
						sprintf(re2, "%03d", (short)als[i] + 1);
					while (*re2) re2++;
				}
			}
			else if (convert_mode_val == 2)
			{
				conversion_memory2[threadid].Alloc(str, MAX_ALLELE_BYTE * 2 + 2);
				if (ploidy == 1)
					sprintf(str, "\t000000");
				else
					sprintf(str, "\t%03d%03d", (short)als[0] + 1, (short)als[1] + 1);
			}
		}
		return str;
	}

	/* Obtain Cervus genotype string */
	TARGET char* GENOTYPE::GetCervusStr()
	{
		int ploidy = Ploidy(), nalleles = Nalleles();
		if (convert_mode_val == 1 && ploidy != 2) Exit("\nError: Cervus supports diploids only.\n");
		char* str = NULL;

		if (convert_mode_val == 1 || convert_mode_val == 2) //disable or truncate
		{
			ushort* als = GetAlleleArray();
			if (nalleles == 0 || convert_mode_val != 1 && ploidy == 1)
			{
				conversion_memory2[threadid].Alloc(str, 3);
				sprintf(str, ",,");
			}
			else
			{
				conversion_memory2[threadid].Alloc(str, MAX_ALLELE_BYTE * 2 + 3);
				sprintf(str, ",%d,%d", als[0] + 1, als[1] + 1);
			}
		}
		return str;
	}

	/* Obtain Arlequin genotype string */
	TARGET char* GENOTYPE::GetArlequinStr()
	{
		int ploidy = Ploidy(), nalleles = Nalleles();
		if (convert_mode_val == 1 && ploidy != 2) Exit("\nError: Arlequin supports diploids only.\n");
		char* str = NULL;

		if (convert_mode_val == 1 || convert_mode_val == 2) //disable or truncate
		{
			conversion_memory2[threadid].Alloc(str, (MAX_ALLELE_BYTE + 2) * 2);
			ushort* als = GetAlleleArray();

			if (nalleles == 0 || convert_mode_val != 1 && ploidy == 1)
			{
				sprintf(str,                       " ?");
				sprintf(str + MAX_ALLELE_BYTE + 2, " ?");
			}
			else
			{
				sprintf(str,                       " %d", als[0] + 1);
				sprintf(str + MAX_ALLELE_BYTE + 2, " %d", als[1] + 1);
			}
		}

		return str;
	}

	/* Obtain Structure genotype string */
	TARGET char* GENOTYPE::GetStructureStr()
	{
		int ploidy = Ploidy(), nalleles = Nalleles();
		char* str = NULL;

		if (convert_mode_val == 1)
		{
			ushort* als = GetAlleleArray();
			conversion_memory2[threadid].Alloc(str, (MAX_ALLELE_BYTE + 2) * ploidy);

			if (nalleles == 0)
				for (int i = 0; i < ploidy; ++i)
					sprintf(str + i * (MAX_ALLELE_BYTE + 2), " -9");
			else
				for (int i = 0; i < ploidy; ++i)
					sprintf(str + i * (MAX_ALLELE_BYTE + 2), " %d", als[i] + 1);
		}
		else if (convert_mode_val == 2)
		{
			ushort* als = GetAlleleArray();
			conversion_memory2[threadid].Alloc(str, (MAX_ALLELE_BYTE + 2) * 2);

			if (nalleles == 0 || ploidy == 1)
			{
				sprintf(str,                         " -9");
				sprintf(str + (MAX_ALLELE_BYTE + 2), " -9");
			}
			else
			{
				sprintf(str,                         " %d", als[0] + 1);
				sprintf(str + (MAX_ALLELE_BYTE + 2), " %d", als[1] + 1);
			}
		}

		return str;
	}

	/* Obtain PolyGene genotype string */
	TARGET char* GENOTYPE::GetPolygeneStr()
	{
		int ploidy = Ploidy(), nalleles = Nalleles();
		char* str = NULL;

		if (convert_mode_val == 1 || convert_mode_val == 2) //disable or truncate
		{
			ushort* als = GetAlleleArray();

			if (nalleles == 0 || convert_mode_val != 1 && ploidy == 1)
			{
				conversion_memory2[threadid].Alloc(str, 2);
				sprintf(str, "\t");
			}
			else if (convert_mode_val == 1)
			{
				conversion_memory2[threadid].Alloc(str, (MAX_ALLELE_BYTE + 1) * ploidy + 1);
				char* re2 = str;
				sprintf(re2++, "\t");
				for (int i = 0; i < ploidy; ++i)
				{
					if (i != 0) sprintf(re2, ",%d", als[i] + 1);
					else sprintf(re2, "%d", als[i] + 1);
					while (*re2) re2++;
				}
			}
			else if (convert_mode_val == 2)
			{
				conversion_memory2[threadid].Alloc(str, (MAX_ALLELE_BYTE + 1) * 2 + 1);
				sprintf(str, "\t%d,%d", als[0] + 1, als[1] + 1);
			}
		}
		return str;
	}

	/* Obtain PolyRelatedness genotype string */
	TARGET char* GENOTYPE::GetPolyRelatednessStr()
	{
		int ploidy = Ploidy(), nalleles = Nalleles();
		char* str = NULL;

		if (convert_mode_val == 1)
		{
			ushort* als = GetAlleleArray();
			conversion_memory2[threadid].Alloc(str, MAX_ALLELE_BYTE * ploidy + 2);

			char* re2 = str;
			sprintf(re2++, "\t"); 
			if (nalleles == 0)
				for (int i = 0; i < ploidy; ++i)
				{ 
					sprintf(re2, "000"); 
					re2 += 3;
				}
			else
				for (int i = 0; i < ploidy; ++i)
				{ 
					sprintf(re2, "%03d", als[i] + 1); 
					while (*re2) re2++;
				}
		}
		else if (convert_mode_val == 2)
		{
			ushort* als = GetAlleleArray();
			conversion_memory2[threadid].Alloc(str, MAX_ALLELE_BYTE * 2 + 2);

			if (nalleles == 0 || convert_mode_val != 1 && ploidy == 1)
				sprintf(str, "\t000000");
			else
				sprintf(str, "\t%03d%03d", als[0] + 1, als[1] + 1);
		}

		return str;
	}

	/* Obtain GenoDive genotype string */
	TARGET char* GENOTYPE::GetGenoDiveStr()
	{
		int ploidy = Ploidy(), nalleles = Nalleles();
		char* str = NULL;

		if (convert_mode_val == 1)
		{
			conversion_memory2[threadid].Alloc(str, MAX_ALLELE_BYTE * ploidy + 2);
			ushort* als = GetAlleleArray();

			char* re2 = str;
			sprintf(re2++, "\t"); 
			if (nalleles == 0) 
				for (int i = 0; i < ploidy; ++i)
				{ 
					sprintf(re2, "000"); 
					re2 += 3; 
				}
			else
				for (int i = 0; i < ploidy; ++i)
				{ 
					sprintf(re2, "%03d", als[i] + 1);
					while (*re2) re2++;
				}
		}
		else if (convert_mode_val == 2)
		{
			conversion_memory2[threadid].Alloc(str, MAX_ALLELE_BYTE * 2 + 2);
			ushort* als = GetAlleleArray();

			if (nalleles == 0 || convert_mode_val != 1 && ploidy == 1)
				sprintf(str, "\t000000");
			else
				sprintf(str, "\t%03d%03d", als[0] + 1, als[1] + 1);
		}
		return str;
	}

	/* Obtain Plink genotype string */
	TARGET char* GENOTYPE::GetPlinkStr()
	{
		int ploidy = Ploidy(), nalleles = Nalleles();
		if (convert_mode_val == 1 && ploidy != 2) Exit("\nError: Plink supports diploids only.\n");
		char* str = NULL;

		if (convert_mode_val == 1 || convert_mode_val == 2)
		{
			ushort* als = GetAlleleArray();
			if (nalleles == 0 || convert_mode_val != 1 && ploidy == 1)
			{
				conversion_memory2[threadid].Alloc(str, 5);
				sprintf(str, "\t0 0");
			}
			else
			{
				conversion_memory2[threadid].Alloc(str, (MAX_ALLELE_BYTE + 1) * 2 + 1);
				sprintf(str, "\t%d %d", als[0] + 1, als[1] + 1);
			}
		}

		return str;
	}
#endif

#ifndef _SLOCUS
	TARGET SLOCUS::SLOCUS()
	{
		SetZero(this, 1);
	}

	/* Convert from LOCUS */
	TARGET SLOCUS::SLOCUS(MEMORY& memory, LOCUS& ref)
	{
		SetVal(this, (SLOCUS*)&ref, 1);

		int gsize = (int)ngeno;
		int asize = flag_alen ? (int)k : 0;
		int ssize = ref.GetEnd() - ref.GetChrom();
		int gasize = (int)ref.GetGenoAlleleSize();//genotype allele array
		int tsize = gsize * sizeof(GENOTYPE) + asize * sizeof(ushort) + ssize + gasize * sizeof(ushort);

		byte* bucket = NULL;
		memory.Alloc(bucket, tsize);
		bits1 = (uint64)bucket;

		//copy genotype array
		SetVal(GetGtab(), ref.GetGtab(), gsize);

		//copy allele length array
		if (flag_alen)
			SetVal(GetAlenArray(), ref.GetAlenArray(), asize);

		//copy chrom ...
		SetVal(GetChrom(), ref.GetChrom(), ssize);

		//copy genotype allele array
		SetVal(GetGenoAlleleArray(), ref.GetGenoAlleleArray(), gasize);

		//set genotype allele offset
		GENOTYPE* gtab = GetGtab();
		ushort* gatab = GetGenoAlleleArray();
		for (uint gi = 0; gi < ngeno; ++gi)
		{
			GENOTYPE& gt = gtab[gi];
			if (gt.Nalleles())
			{
				gt.SetAlleleArray(gatab);
				gatab += gt.Nalleles() + gt.Ploidy();
			}
			else  //missing genotype
				gt.SetAlleleArray((ushort*)-1);
		}
	}

	/* Create unphase locus */
	TARGET SLOCUS::SLOCUS(MEMORY& memory, SLOCUS& ref, TABLE<HASH, uint>& gitab, ushort* gtmap)
	{
		SetVal(this, &ref, 1);
		ushort alleles[N_MAX_PLOIDY] = { 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF };
		gitab.Clear();
		GENOTYPE* refgtab = ref.GetGtab();

		int gasize = 0, gsize = 0;
		for (uint i = 0; i < ref.ngeno; ++i)
		{
			GENOTYPE& gt = refgtab[i];
			int ploidy = gt.Ploidy(), nalleles = gt.Nalleles();
			HASH ha;
			if (nalleles)
			{
				ushort* als = gt.GetAlleleArray();
				SetVal(alleles, als, ploidy);
				Sort(alleles, ploidy);
				ha = HashGenotype(alleles, ploidy);
			}
			else
				ha = missing_hash[ploidy];

			if (!gitab.ContainsKey(ha))
			{
				gitab[ha] = gtmap[i] = (ushort)gitab.size;
				gasize += ploidy + nalleles;
				gsize++;
			}
		}

		int asize = flag_alen ? (int)k : 0;
		int ssize = ref.GetEnd() - ref.GetChrom();
		int tsize = gsize * sizeof(GENOTYPE) + asize * sizeof(ushort) + ssize + gasize * sizeof(ushort);

		byte* bucket = NULL;
		memory.Alloc(bucket, tsize);
		bits1 = (uint64)bucket;
		ngeno = (uint)gitab.size;

		//copy genotype array
		SetVal(GetGtab(), ref.GetGtab(), gsize);

		//copy allele length array
		if (flag_alen)
			SetVal(GetAlenArray(), ref.GetAlenArray(), asize);

		//copy chrom ...
		SetVal(GetChrom(), ref.GetChrom(), ssize);

		//copy genotype allele array
		SetVal(GetGenoAlleleArray(), ref.GetGenoAlleleArray(), gasize);

		//construct gtab
		GENOTYPE* gtab = GetGtab();
		ushort* gatab = GetGenoAlleleArray();
		for (uint i = 0, ng = 0; i < ref.ngeno; ++i)
		{
			GENOTYPE& gt = refgtab[i];
			int ploidy = gt.Ploidy(), nalleles = gt.Nalleles();

			HASH ha;
			if (nalleles)
			{
				ushort* als = gt.GetAlleleArray();
				SetVal(alleles, als, ploidy);
				Sort(alleles, ploidy);
				ha = HashGenotype(alleles, ploidy);
			}
			else
			{
				ha = missing_hash[ploidy];
				SetFF(alleles, ploidy);
			}

			if (gtmap[i] == ng)
				new(&gtab[ng++]) GENOTYPE(gatab, alleles, ploidy);
		}
	}

	/* Create locus for haplotype extraction and Chi-square test */
	TARGET SLOCUS::SLOCUS(MEMORY& memory, SLOCUS& ref, int64 _id, int _ngeno, int _gasize, TABLE<HASH, TEMP_GENOTYPE>& temptab)
	{
		//do not need copy alen
		SetZero(this, 1);

		ngeno = (int)_ngeno;
		flag_pass = true;
		flag_alen = false;

		int gsize = ngeno = (int)_ngeno;
		int asize = 0;
		int ssize = (int)strlen(ref.GetChrom()) + 1 + 3 + CeilLog10((uint64)_id + 1) + 1 + 1;
		int gasize = _gasize;
		int tsize = gsize * sizeof(GENOTYPE) + asize * sizeof(ushort) + ssize + gasize * sizeof(ushort);

		byte* bucket = memory.Alloc(tsize);
		bits1 = (uint64)bucket;

		//no alen
		
		//chrom
		SetZero(GetChrom(), ssize);
		strcpy(GetChrom(), ref.GetChrom());
		sprintf(GetName(), "Loc%llu", _id + 1);//20220717
		
		//no allele identifier!!!
		GENOTYPE* gtab = GetGtab();
		ushort* gatab = (ushort*)(bucket + gsize * sizeof(GENOTYPE) + asize * sizeof(ushort) + ssize);
		for (uint gi = 0; gi < ngeno; ++gi)
		{
			TEMP_GENOTYPE& tgt = temptab(gi);
			new(gtab++) GENOTYPE(gatab, tgt.alleles, tgt.ploidy);
		}
	}

	/* Deep copy from SLOCUS */
	TARGET SLOCUS::SLOCUS(MEMORY& memory, SLOCUS& ref)
	{
		SetVal(this, &ref, 1);
		int tsize = 0;

		//genotype alleles array
		GENOTYPE* refgtab = ref.GetGtab();
		for (uint gi = 0; gi < ngeno; ++gi)
		{
			int nalleles = refgtab[gi].Nalleles();
			tsize += nalleles ? nalleles + refgtab[gi].Ploidy() : 0;
		}
		tsize *= sizeof(ushort);
		tsize += (int)((uint64)ref.GetEnd() - ref.bits1);

		byte* bucket = NULL;
		memory.Alloc(bucket, tsize);
		bits1 = (uint64)bucket;
		SetVal(bucket, (byte*)ref.GetGtab(), tsize);
	}

	/* Get Genotype array */
	TARGET GENOTYPE* SLOCUS::GetGtab()
	{
		return (GENOTYPE*)bits1;
	}

	/* Get end of chrom \0 name \0 {(allele identifiers \0)[k] */
	TARGET char* SLOCUS::GetEnd()
	{
		return StrNextIdx0(GetChrom(), 2 + (ALLELE_IDENTIFIER ? (int64)k : 0)) + 1;
	}

	/* Get chrom string */
	TARGET char* SLOCUS::GetChrom()
	{
		return (char*)bits1 + ngeno * sizeof(GENOTYPE) + (flag_alen ? k * sizeof(ushort) : 0);
	}

	/* Get locus identifier */
	TARGET char* SLOCUS::GetName()
	{
		return StrNextIdx0(GetChrom(), 1) + 1;
	}

	/* Get locus identifier */
	TARGET char* SLOCUS::GetNameStr(char* buf)
	{
		char* chr = GetChrom();

		if (g_locusname_val == 1)
		{
			//chr
			strcpy(buf, chr);
		}
		else if (g_locusname_val == 2)
		{
			//pos
			char* pos = StrNextIdx0(chr, 1) + 1;
			char* ref = StrNextIdx(pos, '_', 1);
			if (ref) *ref = 0;
			strcpy(buf, pos);
			if (ref) *ref = '_';
		}
		else if (g_locusname_val == 3)
		{
			//chr_pos
			char* pos = StrNextIdx0(chr, 1) + 1;
			char* ref = StrNextIdx(pos, '_', 1);
			if (ref) *ref = 0;
			int chrlen = (int)strlen(chr);
			if (chrlen == 0)
				strcpy(buf, pos);
			else
			{
				strcpy(buf, chr);
				buf[chrlen] = '_';
				strcpy(buf + chrlen + 1, pos);
			}
			if (ref) *ref = '_';
		}
		else if (g_locusname_val == 4)
		{
			//pos_ref_alt
			char* pos = StrNextIdx0(chr, 1) + 1;
			char* ref = StrNextIdx(pos, '_', 1);
			int chrlen = (int)strlen(chr);
			if (chrlen == 0)
				strcpy(buf, pos);
			else
			{
				strcpy(buf, chr);
				buf[chrlen] = '_';
				strcpy(buf + chrlen + 1, ref + 1);
			}
		}
		else if (g_locusname_val == 5)
		{
			//chr_pos_ref_alt
			char* pos = StrNextIdx0(chr, 1) + 1;
			int chrlen = (int)strlen(chr);
			if (chrlen == 0)
				strcpy(buf, pos);
			else
			{
				strcpy(buf, chr);
				buf[chrlen] = '_';
				strcpy(buf + chrlen + 1, pos);
			}
		}
		return buf;
	}

	/* Get allele name for vcf/bcf */
	TARGET char* SLOCUS::GetAlleleName(int a)
	{
		return ALLELE_IDENTIFIER ? StrNextIdx0(GetChrom(), 2 + a) + 1 : NULL;
	}

	/* Get alen array */
	TARGET ushort* SLOCUS::GetAlenArray()
	{
		return flag_alen ? (ushort*)(bits1 + ngeno * sizeof(GENOTYPE)) : NULL;
	}

	/* Get genotype allele array */
	TARGET ushort* SLOCUS::GetGenoAlleleArray()
	{
		return (ushort*)GetEnd();
	}

	/* Get Genotype alleles array size Sum(ploidy+nalleles)*/
	TARGET uint SLOCUS::GetGenoAlleleSize()
	{
		uint gasize = 0;
		GENOTYPE* gtab = GetGtab();
		for (uint gi = 0; gi < ngeno; ++gi)
		{
			GENOTYPE gt = gtab[gi];
			if (gt.Nalleles() != 0)
				gasize += gt.Nalleles() + gt.Ploidy();
		}
		return gasize;
	}

	/* Get SMM distance */
	TARGET double SLOCUS::GetSMMDist(int a, int b)
	{
		ushort* alen = GetAlenArray();
		return ((int)alen[a] - (int)alen[b]) * ((int)alen[a] - (int)alen[b]);
	}
#endif

#ifndef _LOCUS

	/* Initialize */
	TARGET LOCUS::LOCUS()
	{
		SetZero(this, 1);
	}

	/* Deep Copy Locus */
	TARGET LOCUS::LOCUS(MEMORY& memory, int64 _id, LOCUS& ref)
	{
		SetVal(this, &ref, 1);
		id = (LOCN)_id;

		//copy SLOCUS
		new(this) SLOCUS(memory, ref);

		GENOTYPE* refgtab = ref.GetGtab();
		int gsize = (int)ref.ngeno; gasize = 0;
		for (uint gi = 0; gi < ngeno; ++gi)
		{
			GENOTYPE& gt = refgtab[gi];
			int ploidy = gt.Ploidy(), nalleles = gt.Nalleles();
			if (nalleles)
				gasize += ploidy + nalleles;
		}
		int asize = flag_alen ? (int)k : 0;
		int ssize = ref.GetEnd() - ref.GetChrom();
		int tsize = gsize * sizeof(GENOTYPE) + asize * sizeof(ushort) + ssize + gasize * sizeof(ushort);

		//alloc memory
		byte* bucket = NULL;
		memory.Alloc(bucket, tsize);
		bits1 = (uint64)bucket;
		ngeno = (uint)gsize;

		//set interal variables
		_chrom = (char*)bits1 + ngeno * sizeof(GENOTYPE) + (flag_alen ? k * sizeof(ushort) : 0);
		_alen = flag_alen ? (ushort*)(bits1 + ngeno * sizeof(GENOTYPE)) : 0;
		_genoallele = (ushort*)(StrNextIdx0(GetChrom(), 2 + (ALLELE_IDENTIFIER ? (int64)k : 0)) + 1);

		//copy genotype array
		SetVal(GetGtab(), ref.GetGtab(), gsize);

		//copy allele length array
		if (flag_alen)
			SetVal(GetAlenArray(), ref.GetAlenArray(), asize);

		//copy chrom ...
		SetVal(GetChrom(), ref.GetChrom(), ssize);

		//copy genotype allele array
		SetVal(GetGenoAlleleArray(), ref.GetGenoAlleleArray(), gasize);

		//set genotype allele offset
		GENOTYPE* gtab = GetGtab();
		ushort* gatab = GetGenoAlleleArray();
		for (uint gi = 0; gi < ngeno; ++gi)
		{
			GENOTYPE& gt = gtab[gi];
			if (gt.Nalleles())
			{
				gt.SetAlleleArray(gatab);
				gatab += gt.Nalleles() + gt.Ploidy();
			}
			else  //missing genotype
				gt.SetAlleleArray((ushort*)-1);
		}

		//construct gftab
		new(&gftab) TABLE<HASH, GENOTYPE*>(true, &memory, (int)ngeno);
		auto bucket2 = ref.gftab.bucket;
		auto index2  = ref.gftab.index;
		for (uint gi = 0; gi < ngeno; ++gi)
			gftab[bucket2[index2[gi]].key] = &gtab[gi];
	}

	/* Create unphase locus */
	TARGET LOCUS::LOCUS(MEMORY& memory, LOCUS& ref, TABLE<HASH, uint>& gitab, ushort* gtmap)
	{
		SetVal(this, &ref, 1);

		ushort alleles[N_MAX_PLOIDY] = { 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF };
		gitab.Clear();

		GENOTYPE* refgtab = ref.GetGtab();
		int gsize = 0; gasize = 0;
		for (uint gi = 0; gi < ngeno; ++gi)
		{
			GENOTYPE& gt = refgtab[gi];
			int ploidy = gt.Ploidy(), nalleles = gt.Nalleles();
			HASH ha;
			if (nalleles)
			{
				ushort* als = gt.GetAlleleArray();
				SetVal(alleles, als, ploidy);
				Sort(alleles, ploidy);
				ha = HashGenotype(alleles, ploidy);
			}
			else
				ha = missing_hash[ploidy];

			if (!gitab.ContainsKey(ha))
			{
				gitab[ha] = gtmap[gi] = (ushort)gitab.size;
				gasize += ploidy + nalleles;
				gsize++;
			}
		}

		int asize = flag_alen ? (int)k : 0;
		int ssize = ref.GetEnd() - ref.GetChrom();
		int tsize = gsize * sizeof(GENOTYPE) + asize * sizeof(ushort) + ssize + gasize * sizeof(ushort);

		//alloc memory
		byte* bucket = NULL;
		memory.Alloc(bucket, tsize);
		bits1 = (uint64)bucket;
		ngeno = (uint)gsize;

		//set internal variables
		_chrom = (char*)bits1 + ngeno * sizeof(GENOTYPE) + (flag_alen ? k * sizeof(ushort) : 0);
		_alen = flag_alen ? (ushort*)(bits1 + ngeno * sizeof(GENOTYPE)) : 0;
		_genoallele = (ushort*)(StrNextIdx0(GetChrom(), 2 + (ALLELE_IDENTIFIER ? (int64)k : 0)) + 1);

		//copy genotype array
		SetVal(GetGtab(), ref.GetGtab(), gsize);

		//copy allele length array
		if (flag_alen)
			SetVal(GetAlenArray(), ref.GetAlenArray(), asize);

		//copy chrom ...
		SetVal(GetChrom(), ref.GetChrom(), ssize);

		//copy genotype allele array
		SetVal(GetGenoAlleleArray(), ref.GetGenoAlleleArray(), gasize);

		//construct gftab
		new(&gftab) TABLE<HASH, GENOTYPE*>(true, &memory, (int)ngeno);

		//construct genotypes
		GENOTYPE* gtab = GetGtab();
		ushort* gatab = GetGenoAlleleArray();
		for (uint i = 0, ng = 0; i < ref.ngeno; ++i)
		{
			GENOTYPE& gt = refgtab[i];
			int ploidy = gt.Ploidy(), nalleles = gt.Nalleles();

			HASH ha;
			if (nalleles)
			{
				ushort* als = gt.GetAlleleArray();
				SetVal(alleles, als, ploidy);
				Sort(alleles, ploidy);
				ha = HashGenotype(alleles, ploidy);
			}
			else
			{
				ha = missing_hash[ploidy];
				SetFF(alleles, ploidy);
			}

			if (gtmap[i] == ng)
			{
				new(&gtab[ng]) GENOTYPE(gatab, alleles, ploidy);
				gftab[ha] = &gtab[ng++];
			}
		}
	}

	/* Create locus for haplotype extraction and Chi-square test */
	TARGET LOCUS::LOCUS(MEMORY& memory, LOCUS& ref, int64 _id, int _ngeno, int _gasize, TABLE<HASH, TEMP_GENOTYPE>& temptab)
	{
		//do not need copy alen
		SetZero(this, 1);
		ngeno = (uint)_ngeno;
		flag_pass = true;
		flag_alen = false;

		id = (LOCN)_id;
		pos = ref.pos;

		//int gsize = ngeno;
		//int asize = 0;
		int ssize = (int)strlen(ref.GetChrom()) + 1 + 3 + CeilLog10(id + 1) + 1 + 1;

		gasize = _gasize;

		byte* bucket = memory.Alloc(ngeno * sizeof(GENOTYPE) + gasize * sizeof(ushort));
		bits1 = (uint64)bucket;
		GENOTYPE* gtab = (GENOTYPE*)bits1;

		_genoallele = (ushort*)(gtab + ngeno);
		_alen = NULL;
		 memory.Alloc(_chrom, ssize);

		//no alen

		//chrom
		SetZero(GetChrom(), ssize);
		strcpy(GetChrom(), ref.GetChrom());
		sprintf(GetName(), "Loc%llu", (uint64)id + 1); //20220717

		//no allele identifier

		//gftab
		new(&gftab) TABLE<HASH, GENOTYPE*>(true, &memory, ngeno);
		ushort* gatab = GetGenoAlleleArray();
		for (uint gi = 0; gi < ngeno; ++gi)
		{
			TEMP_GENOTYPE& tgt = temptab(gi);
			gftab[tgt.hash] = new(gtab++) GENOTYPE(gatab, tgt.alleles, tgt.ploidy);
		}
	}

	/* For dummy locus for collapse alleles during testing genotype distributions */
	TARGET LOCUS::LOCUS(MEMORY& memory, SLOCUS& ref)
	{
		//do not need copy alen
		SetZero(this, 1);
		ngeno = ref.ngeno;
		id = (LOCN)-1;
		flag_pass = true;
		flag_original = true;
		flag_qual = true;

		//use original genotype and genotype allele memory
		bits1 = (uint64)ref.GetGtab();
		_genoallele = ref.GetGenoAlleleArray();

		//set chrom, alen, genoallelearray
		_chrom = ref.GetChrom();
		_alen = ref.GetAlenArray();
		_genoallele = ref.GetGenoAlleleArray();
		gasize = 0;

		//create gftab
		GENOTYPE* gtab = GetGtab();
		new(&gftab) TABLE<HASH, GENOTYPE*>(true, &memory, ngeno);
		for (uint gi = 0; gi < ngeno; ++gi)
		{
			GENOTYPE& gt = gtab[gi];
			HASH hash = gt.Hash();
			if (gt.Nalleles())
				gasize += gt.Nalleles() + gt.Ploidy();
			if (!gftab.ContainsKey(hash))
				gftab[hash] = &gt;
		}
	}

	/* For non-vcf input, set locus name and id, SLOCUS */
	TARGET LOCUS::LOCUS(MEMORY& memory, char* line, int64 _id, int _ngenotype, GENOTYPE*& gtab, ushort*& gatab)
	{
		//_alen is set at IndexAlleleLength
		int _gasize = gasize;
		SetZero(this, 1);
		ngeno = _ngenotype;
		gasize = _gasize;
		id = _id;
		flag_pass = true;
		flag_original = true;
		flag_qual = true;
		flag_indel = false;

		gtab = (GENOTYPE*)memory.Alloc(ngeno * sizeof(GENOTYPE) + gasize * sizeof(ushort));
		bits1 = (uint64)gtab;
		_genoallele = gatab = (ushort*)(gtab + ngeno);

		new(&gftab) TABLE<HASH, GENOTYPE*>(true, &memory, ngeno);

		if (line)
		{
			//no chrom, set locus name
			int namelen = 0;
			while (line[namelen] &&
				line[namelen] != ' '  && line[namelen] != '\t' &&
				line[namelen] != '\r' && line[namelen] != '\n')
				namelen++;

			int tsize = 1 + namelen + 1 + 1;

			memory.Alloc(_chrom, tsize);
			SetZero(_chrom, namelen + 3);
			SetVal(_chrom + 1, line, namelen);
		}
		else
		{

			int tsize = 1 + 3 + CeilLog10(id + 1) + 1 + 1;

			//no chrom, no locus name, name it
			memory.Alloc(_chrom, tsize);
			SetZero(_chrom, CeilLog10(id + 1) + 5);
			sprintf(_chrom + 1, "Loc%llu", (uint64)id + 1);
		}
	}
	
	/* For vcf input, set locus name and id, SLOCUS */
	TARGET LOCUS::LOCUS(MEMORY& memory, char*& line, uint64 _mask, int _ngenotype, GENOTYPE*& gtab, ushort*& gatab)
	{
		//do not need copy alen
		int _gasize = gasize;
		SetZero(this, 1);
		ngeno = _ngenotype;
		gasize = _gasize;
		pos = 0xFFFFFFFFFFFFFFFF;
		flag_pass = true;
		gtab = (GENOTYPE*)memory.Alloc(ngeno * sizeof(GENOTYPE) + gasize * sizeof(ushort));
		bits1 = (uint64)gtab;
		_genoallele = gatab = (ushort*)(gtab + ngeno);

		new(&gftab) TABLE<HASH, GENOTYPE*>(true, &memory, ngeno);

		//locate address of each string
		char* base = line, *pstr = NULL;
		while (*line != '\t') line++; *line++ = '\0';

		if (line[0] != '.') pos = ReadLong(line);
		while (*line != '\t') line++; *line++ = '\0';

		char* name = line;
		while (*line != '\t') line++; *line++ = '\0';

		char* ref = line;
		while (*line != '\t') line++; *line++ = '\0';

		char* alt = line;
		while (*line != '\t') line++; *line++ = '\0';

		//write chrom string
		if ((name[0] == '.' && name[1] == '\0') || name[0] == '\0')
		{
			//no loc name, name it
			int slen = (int)(strlen(base) + 1 + strlen(ref) + 1 + strlen(alt) + 1) * 2 + CeilLog10(pos) + 1;

			//alloc
			memory.Alloc(_chrom, slen);
			pstr = _chrom;
			SetZero(pstr, slen);

			//copy chrom
			strcpy(pstr, base);
			pstr += strlen(base) + 1;

			//set name
			sprintf(pstr, "%lld_%s_%s", pos, ref, alt);
		}
		else
		{
			//has loc name
			int slen = (int)(strlen(base) + 1 + strlen(name) + 1 + strlen(ref) + 1 + strlen(alt) + 1);

			//alloc
			memory.Alloc(_chrom, slen);
			pstr = _chrom;
			SetZero(pstr, slen);

			//copy chrom
			strcpy(pstr, base);
			pstr += strlen(base) + 1;

			//copy name
			strcpy(pstr, name);
		}

		//relocate name
		name = pstr;

		//copy allele identifiers
		pstr += strlen(name) + 1;
		sprintf(pstr, "%s,%s", ref, alt);

		//obtain number of alleles
		k = (ushort)(CountChar(alt, ',') + 2);
		ReplaceChar(pstr, ',', '\0');
		ReplaceChar(alt, ',', '\0');

		//SNP or indel
		int maxal = 0, minal = 65536000;
		for (int i = 0; i < k; ++i)
		{
			int tl = (int)strlen(pstr);
			if (tl > maxal) maxal = tl;
			if (tl < minal) minal = tl;
			pstr += tl + 1;
		}
		flag_indel = minal != maxal;

		//QUAL to FILTER
		while (*line != '\t') line++; *line++ = '\0';

		if (info_filter)								//MASK INFO Filter
		{
			if (f_qual_b && (_mask & 0x4))				//fail in qual filter
				flag_pass = false;
			if (f_type_b)
			{
				if (f_type_val == 1 && flag_indel)		//snp
					flag_pass = false;
				else if (f_type_val == 2 && !flag_indel) //indel
					flag_pass = false;
			}
			if (f_original_b && f_original_val == 1 && (_mask & 0x2)) //fail in original filter
				flag_pass = false;
		}

		//FILTER to INFO
		while (*line != '\t') line++; *line++ = '\0';

		//INFO to FORMAT
		for (;;)
		{
			while (*line != ';' && *line != '\t') line++;
			if (*line == ';')
				*line++ = '\0';
			else
			{
				*line++ = '\0';
				break;
			}
		}

		char* end = StrNextIdx(line + 1, '\t', 1);
		*end = '\0';
		format_size = (ushort)CountChar(line, ':') + 1;
		*end = '\t';

		//find offset of each format string (from base)
		ushort format_offset[1024];
		//VLA_NEW(format_offset, ushort, format_size);
		for (int i = 0; line < end; ++i)
		{
			format_offset[i] = (ushort)(line - base);
			while (*line != ':' && *line != '\t') line++;
			if (*line == ':')
				*line++ = '\0';
			else
			{
				*line++ = '\0';
				break;
			}
		}

		dpid = gqid = adid = 0xFFFF;
		if (abs( g_format_val ) <= BCF)
		{
			gtid = GetFormatId(base, (char*)"GT", format_offset);
			if (gtid == 0xFFFF) //must have gt
				Exit("\nExit: variant %s does not have GT format_offset tag. \n", name);

			gqid = GetFormatId(base, (char*)"GQ", format_offset);
			flag_hasgq = gqid != 0xFFFF;

			dpid = GetFormatId(base, (char*)"DP", format_offset);
			flag_hasdp = dpid != 0xFFFF;

			adid = GetFormatId(base, (char*)"AD", format_offset);
			flag_hasad = adid != 0xFFFF;
		}
		//VLA_DELETE(format_offset);
	}

	/* Get end of chrom \0 name \0 {(allele identifiers \0)[k] */
	TARGET char* LOCUS::GetEnd()
	{
		return StrNextIdx0(GetChrom(), 2 + (ALLELE_IDENTIFIER ? k : 0)) + 1;
	}

	/* Get chrom string */
	TARGET char* LOCUS::GetChrom()
	{
		return _chrom;
	}

	/* Get locus identifier */
	TARGET char* LOCUS::GetName()
	{
		return StrNextIdx0(GetChrom(), 1) + 1;
	}

	/* Get locus identifier */
	TARGET char* LOCUS::GetNameStr(char* buf)
	{
		char* chr = GetChrom();

		if (g_locusname_val == 1)
		{
			//chr
			strcpy(buf, chr);
		}
		else if (g_locusname_val == 2)
		{
			//pos
			char* pos2 = StrNextIdx0(chr, 1) + 1;
			char* ref = StrNextIdx(chr, '_', 1);
			if (ref) *ref = 0;
			strcpy(buf, pos2);
			if (ref) *ref = '_';
		}
		else if (g_locusname_val == 3)
		{
			//chr_pos
			char* pos2 = StrNextIdx0(chr, 1) + 1;
			char* ref = StrNextIdx(chr, '_', 1);
			if (ref) *ref = 0;
			int chrlen = (int)strlen(chr);
			if (chrlen == 0)
				strcpy(buf, pos2);
			else
			{
				strcpy(buf, chr);
				buf[chrlen] = '_';
				strcpy(buf + chrlen + 1, pos2);
			}
			if (ref) *ref = '_';
		}
		else if (g_locusname_val == 4)
		{
			//chr_ref_alt
			//char* pos2 = StrNextIdx0(chr, 1) + 1;
			char* ref = StrNextIdx(chr, '_', 1);
			int chrlen = (int)strlen(chr);
			if (chrlen == 0)
				strcpy(buf, ref + 1);
			else
			{
				strcpy(buf, chr);
				buf[chrlen] = '_';
				strcpy(buf + chrlen + 1, ref + 1);
			}
		}
		else if (g_locusname_val == 5)
		{
			//pos_ref_alt
			char* pos2 = StrNextIdx0(chr, 1) + 1;
			strcpy(buf, pos2);
		}
		else if (g_locusname_val == 6)
		{
			//chr_pos_ref_alt
			char* pos2 = StrNextIdx0(chr, 1) + 1;
			int chrlen = (int)strlen(chr);
			if (chrlen == 0)
				strcpy(buf, pos2);
			else
			{
				strcpy(buf, chr);
				buf[chrlen] = '_';
				strcpy(buf + chrlen + 1, pos2);
			}
		}
		return buf;
	}

	/* Get allele name for vcf/bcf */
	TARGET char* LOCUS::GetAlleleName(int a)
	{
		return ALLELE_IDENTIFIER ? StrNextIdx0(GetChrom(), 2 + a) + 1 : NULL;
	}

	/* Get alen array */
	TARGET ushort* LOCUS::GetAlenArray()
	{
		return _alen;
	}

	/* Get genotype allele array */
	TARGET ushort* LOCUS::GetGenoAlleleArray()
	{
		return _genoallele;
	}

	/* Get Genotype alleles array size Sum(ploidy+nalleles)*/
	TARGET uint LOCUS::GetGenoAlleleSize()
	{
		if (gasize == 0)
		{
			GENOTYPE* gtab = GetGtab();
			for (uint gi = 0; gi < ngeno; ++gi)
			{
				GENOTYPE gt = gtab[gi];
				if (gt.Nalleles() != 0)
					gasize += gt.Nalleles() + gt.Ploidy();
			}
		}
		return gasize;
	}

	/* Get index of a target format */
	TARGET ushort LOCUS::GetFormatId(char* base, char* fmt_name, ushort* fmt_offset)
	{
		for (int i = 0; i < format_size; ++i)
			if (LwrParCmp(fmt_name, base + fmt_offset[i]) == 0)
				return (ushort)i;
		return 0xFFFF;
	}
#endif

#ifndef _GENO_READER

	/* Do nothing */
	TARGET GENO_READER::GENO_READER()
	{

	}

	/* Initialize reader */
	TARGET GENO_READER::GENO_READER(int indid, int64 l, BUCKET* bucket)
	{
		//set pos and size
		SetZero(this, 1);

		//set bucket from default bucket or assigned bucket
		if (bucket == NULL) bucket = &geno_bucket;
		pos = (uint64*)(bucket->base_addr + bucket->offset[l].offset);
		size = bucket->offset[l].size;

		if (indid)
		{
			//skip indid * size bits
			int nreadbits = indid * size;

			//pos move nreadbits / 64
			pos += nreadbits >> 6;

			//remain bits to read
			nreadbits &= 63;

			//read 64 bits to data and read remain bits from data
			data = *pos++ >> nreadbits;

			//set nbits
			nbits = 64 - nreadbits;
		}
		else
		{
			//read 64 bits
			data = *pos++;

			//set nbits
			nbits = 64;
		}
	}

	/* Get id of next ind (order by indid) */
	TARGET int GENO_READER::Read()
	{
		// if data is empty
		if (nbits < size) [[unlikely]]
		{
			// move nbits data to gid
			int gid = (int)data;
			
			// remain number of bits to read
			int rbits = size - nbits;

			// read 64 bits to data
			data = *pos++;

			// read rbits from data and concate to higher bits in gid
			gid |= ((uint)data & ((1u << rbits) - 1u)) << nbits;

			//shift right
			data >>= rbits;
			nbits = 64 - rbits;
			return gid;
		}
		else [[likely]]
		{
			//read size bits
			int gid = (int)data & ((1u << size) - 1u);

			//shift right
			data >>= size;
			nbits -= size;
			return gid;
		}
	}
#endif

#ifndef _GENO_WRITER

	TARGET GENO_WRITER::GENO_WRITER()
	{

	}

	/* Initialize writer */
	TARGET GENO_WRITER::GENO_WRITER(int64 l, BUCKET* bucket)
	{
		SetZero(this, 1);
		if (bucket == NULL) bucket = &geno_bucket;
		pos = (uint64*)(bucket->base_addr + bucket->offset[l].offset);
		size = (uint)bucket->offset[l].size;
	}

	/* Write id of next ind to buffer */
	TARGET void GENO_WRITER::Write(uint gid)
	{
		gid &= (1u << size) - 1u;

		//if data is full
		if (nbits + size > 64) [[unlikely]]
		{
			if (nbits == 64)
			{
				*pos++ = data;
				data = (uint64)gid;
				nbits = size;
			}
			else
			{
				// number of bits can write to data
				int wbits = 64 - nbits;

				// write wbits to data
				data |= (uint64)gid << nbits;

				// write 64 bit data to pos
				*pos++ = data;

				// write remain bits to data
				data = gid >> wbits;
				nbits = size - wbits;
			}
		}
		else [[likely]]
		{
			data |= (uint64)gid << nbits;
			nbits += size;
		}
	}

	/* Write all remaining bits to buffer */
	TARGET void GENO_WRITER::FinishWrite()
	{
		byte* pos2 = (byte*)pos, *data2 = (byte*)&data;
		for (int i = 0; i * 8 < nbits; ++i)
			pos2[i] = data2[i];
	}

#endif

#ifndef _IND
	/* Initialize */
	template<typename REAL>
	TARGET IND<REAL>::IND()
	{
		indid = 0xFFFF;
		name = NULL;
	}

	/* Create individual for non-vcf input */
	template<typename REAL>
	TARGET IND<REAL>::IND(char* t, bool iscount, int id, GENOTYPE** gtab, ushort** gatab, GENO_WRITER* wt)
	{
		indid = id;
		vmin = vmax = 0;
		switch (abs( g_format_val ))//genepop|spagedi|cervus|arlequin|structure|polygene|polyrelatedness|genodive|plink
		{
		case GENEPOP: genepop(t, iscount, gtab, gatab, wt); break;
		case SPAGEDI: spagedi(t, iscount, gtab, gatab, wt); break;
		case CERVUS: cervus(t, iscount, gtab, gatab, wt); break;
		case ARLEQUIN: arlequin(t, iscount, gtab, gatab, wt); break;
		case STRUCTURE: structure(t, iscount, gtab, gatab, wt); break;
		case POLYGENE: polygene(t, iscount, gtab, gatab, wt); break;
		case POLYRELATEDNESS: polyrelatedness(t, iscount, gtab, gatab, wt); break;
		case GENODIVE: genodive(t, iscount, gtab, gatab, wt); break;
		case PLINK: plink(t, iscount, gtab, gatab, wt); break;
		}
	}

	/* Create individual for vcf/bcf input */
	template<typename REAL>
	TARGET IND<REAL>::IND(char*& title, int id)
	{
		//from VCF line
		name = title;

		while (*title != '\t' && *title != '\n' && *title != '\0') title++;
		*title++ = '\0';

		char* tname;
		individual_memory->Alloc(tname, (int)strlen(name) + 1);
		//tname = new char[strlen(name) + 1];
		strcpy(tname, name);
		name = tname;
		indid = id;
		vmin = vmax = 0;
	}

	/* Create individual from a reference individual */
	template<typename REAL>
	TARGET IND<REAL>::IND(IND& ref)
	{
		SetVal(this, &ref, 1);
		individual_memory->Alloc(name, (int)strlen(ref.name) + 1);
		strcpy(name, ref.name);
	}

	/* Unnitialize */
	template<typename REAL>
	TARGET IND<REAL>::~IND()
	{

	}

	/* Set individual genotype with default bucket */
	template<typename REAL>
	TARGET void IND<REAL>::SetGenotype(int64 l, uint gid)
	{
		uint size = (uint)geno_bucket.offset[l].size;
		uint64 offset = size * indid;
		uint* pos = (uint*)(geno_bucket.base_addr + geno_bucket.offset[l].offset + (offset >> 3));
		offset &= 7;
		*pos = (*pos & (~(((1u << size) - 1u) << offset))) | (gid << offset);
	}

	/* Set individual genotype with local bucket */
	template<typename REAL>
	TARGET void IND<REAL>::SetGenotype(int64 l, uint gid, OFFSET* _offset, byte* bucket)
	{
		uint size = (uint)_offset[l].size;
		uint64 offset = size * indid;
		uint* pos = (uint*)(bucket + _offset[l].offset + (offset >> 3));
		offset &= 7;
		*pos = (*pos & (~(((1u << size) - 1u) << offset))) | (gid << offset);
	}

	/* Get index for a pair of genotype */
	template<typename REAL>
	TARGET void IND<REAL>::GetDyadGenotypeIdx(int& id1, int& id2, int64 l)
	{
		OFFSET ot = geno_bucket.offset[l];
		uint size = (uint)ot.size;
		byte* bucket = geno_bucket.base_addr + ot.offset;
		uint64 o1 = size * id1, o2 = size * id2;
		uint* pos1 = (uint*)(bucket + (o1 >> 3));
		uint* pos2 = (uint*)(bucket + (o2 >> 3));
		id1 = (*pos1 >> (o1 & 7)) & ((1u << size) - 1u);
		id2 = (*pos2 >> (o2 & 7)) & ((1u << size) - 1u);
	}

	/* Get individual genotype index from default table */
	template<typename REAL>
	TARGET int IND<REAL>::GetGenotypeId(int64 l, byte* _bucket, OFFSET* _offset)
	{
		uint size = (uint)_offset[l].size;
		uint64 offset = size * indid;
		uint* pos = (uint*)(_bucket + _offset[l].offset + (offset >> 3));
		return (*pos >> (offset & 7)) & ((1u << size) - 1u);
	}

	/* Get individual genotype index from default table */
	template<typename REAL>
	TARGET int IND<REAL>::GetGenotypeId(int64 l)
	{
		uint size = (uint)geno_bucket.offset[l].size;
		uint64 offset = size * indid;
		uint* pos = (uint*)(geno_bucket.base_addr + geno_bucket.offset[l].offset + (offset >> 3));
		return (*pos >> (offset & 7)) & ((1u << size) - 1u);   //5% faster than GENO_MASK32[size]
	}

	/* Get individual genotype from default table */
	template<typename REAL>
	TARGET GENOTYPE& IND<REAL>::GetGenotype(int64 l)
	{
		return GetLoc(l).GetGtab()[GetGenotypeId(l)];
	}

	/* Get individual genotype from local table */
	template<typename REAL>
	TARGET GENOTYPE& IND<REAL>::GetGenotype(int64 l, TABLE<HASH, GENOTYPE*>& gftab)
	{
		return *gftab(GetGenotypeId(l));
	}

	/* Get individual genotype from local table */
	template<typename REAL>
	TARGET GENOTYPE& IND<REAL>::GetGenotype(int64 l, GENOTYPE* gtab)
	{
		return gtab[GetGenotypeId(l)];
	}

	/* Calculate the genotypic frequency */
	template<typename REAL>
	TARGET double IND<REAL>::GenoFreq(POP<REAL>* grp, int model, int64 loc, double e)
	{
		//Logarithm
		int64 slog = 0; double prod = 1;
		int64 st = loc == (int64)-1 ? 0 : loc;
		int64 ed = loc == (int64)-1 ? nloc : loc + 1;
		double e2 = e;
		double e1 = 1 - e2;

		OpenLog(slog, prod);//slog,prod
		for (int64 l = st; l < ed; ++l)
		{
			GENOTYPE& gt = GetGenotype(l);//fine
			if (gt.Nalleles() == 0) continue;
			int m = model == 4 ? GetLoc(l).pes_model : model;
			double pr = gt.GFZ(m, grp->GetFreq(l));
			ChargeLog(slog, prod, e1 * pr + e2 * gt.GFZ(m, total_pop->GetFreq(l)) * (1 - pr));
		}
		CloseLog(slog, prod);

		if (loc >= 0 && prod == 0) prod = NAN;
		return prod;
	}
#endif

#ifndef _POP
	/* Initialize */
	template<typename REAL>
	TARGET POP<REAL>::POP()
	{
		SetZero(this, 1);
	}

	/* Create a pop */
	template<typename REAL>
	TARGET POP<REAL>::POP(char* _name, char** _names, int _nind, int _regid, int _npop, int _id, bool _ispop)
	{
		ispop = _ispop;
		id = _id;
		name = _name;
		names = _names;
		nind = _nind;
		rid = _regid;
		inds = NULL;
		loc_stat1 = NULL;
		allelefreq = NULL;
		genocount = NULL;

		npop = _npop;
		vpop = NULL;
	}

	/* Get genotype count array */
	template<typename REAL>
	TARGET ushort* POP<REAL>::GetGenoCount(int64 l)
	{
		return genocount + genotype_count_offset[l];
	}

	/* Get genotype count of a genotype */
	template<typename REAL>
	TARGET ushort POP<REAL>::GetGenoCount(int64 l, int gid)
	{
		return *(genocount + genotype_count_offset[l] + gid);
	}

	/* Get allele frequencies array */
	template<typename REAL>
	TARGET REAL* POP<REAL>::GetFreq(int64 l)
	{
		return allelefreq + allele_freq_offset[l];
	}

	/* Get allele frequency of an allele */
	template<typename REAL>
	TARGET REAL POP<REAL>::GetFreq(int64 l, int a)
	{
		return *(allelefreq + allele_freq_offset[l] + a);
	}

	/* Uninitialize a pop, unalloc memory except population name */
	template<typename REAL>
	TARGET void POP<REAL>::Uninitialize1()
	{
		if (allelefreq && allelefreq != total_pop->allelefreq) 
		{
			delete[] allelefreq; 
			allelefreq = NULL;
		}

		if (genocount)
		{
			delete[] genocount;
			genocount = NULL;
		}

		if (loc_stat1 != total_pop->loc_stat1)
		{
			delete[] loc_stat1; 
			loc_stat1 = NULL;
		}
		
		if (names) 
		{ 
			delete[] names;
			names = NULL;
		}
	}

	/* Uninitialize a pop, unalloc population name */
	template<typename REAL>
	TARGET void POP<REAL>::Uninitialize2()
	{
		if (allelefreq)
		{
			delete[] allelefreq;
			allelefreq = NULL;
		}

		if (genocount)
		{
			delete[] genocount;
			genocount = NULL;
		}

		if (loc_stat1)
		{
			delete[] loc_stat1;
			loc_stat1 = NULL;
		}

		if (names)
		{
			delete[] names;
			names = NULL;
		}

		if (name)
		{
			delete[] name;
			name = NULL;
		}
	}

	/* Uncllocate memory for locstat, allele frequency, genotype count */
	template<typename REAL>
	TARGET void POP<REAL>::UnAllocFreq()
	{
		if (loc_stat1) 
		{ 
			delete[] loc_stat1; 
			loc_stat1 = NULL;
		}

		if (allelefreq) 
		{ 
			delete[] allelefreq; 
			allelefreq = NULL; 
		}

		if (genocount) 
		{ 
			delete[] genocount; 
			genocount = NULL; 
		}
	}

	/* Allocate memory for locstat, allele frequency, genotype count */
	template<typename REAL>
	TARGET void POP<REAL>::AllocFreq()
	{
		//to apply diversity filter and test genotype distributions
		loc_stat1 = new LOCSTAT1[nloc];		SetZero(loc_stat1, nloc);
		allelefreq = new REAL[KT];			SetZero(allelefreq, KT);
		genocount = new ushort[GT];			SetZero(genocount, GT);
	}

	/* Move after filter locus */
	template<typename REAL>
	TARGET void POP<REAL>::MoveFreq(LOCN* nafoffset, LOCN* ngcoffset)
	{
		LOCSTAT1 *nloc_stat = NULL;
		if (loc_stat1)
		{
			nloc_stat = new LOCSTAT1[nfilter];
			SetZero(nloc_stat, nfilter);
		}

		REAL* nallelefreq = NULL;
		if (allelefreq)
		{
			nallelefreq = new REAL[KT];
			SetZero(nallelefreq, KT);
		}

		ushort* ngenocount = NULL;
		if (genocount)
		{
			ngenocount = new ushort[GT];
			SetZero(ngenocount, GT);
		}

		//move freq and count to another compat memory
		for (int64 l = 0, nl = 0; l < nloc; ++l)
		{
			if (GetLoc(l).flag_pass)
			{
				if (loc_stat1)
					nloc_stat[nl] = loc_stat1[l];
				if (allelefreq)
					SetVal(nallelefreq + nafoffset[nl], GetFreq(l), GetLoc(l).k);
				if (genocount)
					SetVal(ngenocount + ngcoffset[nl], GetGenoCount(l), GetLoc(l).ngeno);
				nl++;
			}
		}

		if (loc_stat1)
		{
			delete[] loc_stat1;
			loc_stat1 = nloc_stat;
		}

		if (allelefreq)
		{
			delete[] allelefreq;
			allelefreq = nallelefreq;
		}

		if (genocount)
		{
			delete[] genocount;
			genocount = ngenocount;
		}
	}

	/* Clear memory for locstat, allele frequency, genotype count */
	template<typename REAL>
	TARGET void POP<REAL>::ClearFreqGcount()
	{
		SetZero(loc_stat1, nloc);
		SetZero(allelefreq, KT);
		SetZero(genocount, GT);
	}

	/* Calculate loc stat, allele frequency, genotype count */
	template<typename REAL>
	TARGET void POP<REAL>::CalcFreqGcount()
	{
		if (nind == 0) return;

		if (ad == 0)
		{
#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
			for (int64 l = 0; l < nloc; ++l)
			{
				GENOTYPE* gtab = GetLoc(l).GetGtab();
				int ngeno = GetLoc(l).ngeno;

				REAL* p = GetFreq(l);
				ushort* gcount = GetGenoCount(l);
				LOCSTAT1 &stat1 = loc_stat1[l];

				GENO_READER rt(ind0id, l);
				for (int i = 0; i < nind; ++i)
					gcount[rt.Read()]++;

				int nhaplo = 0;
				for (int gi = 0; gi < ngeno; ++gi)
				{
					int gc = gcount[gi];
					if (gc == 0) continue;

					GENOTYPE& gt = gtab[gi]; 
					if (gt.Nalleles() == 0) continue;

					ushort* als = gt.GetAlleleArray();
					int v = gt.Ploidy();
					nhaplo += v * gc;

					for (int j = 0; j < v; ++j)
						p[als[j]] += gc;
				}

				int k2 = GetLoc(l).k;
				stat1.k = (ushort)CountK(p, k2);
				stat1.nhaplo = nhaplo;
				Unify(p, k2);

				PROGRESS_VALUE += nind;
			}
		}
		else
		{
			loc_stat1 = new LOCSTAT1[nloc];
			genocount = new ushort[GT];
			SetZero(genocount, GT);

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
			for (int64 l = 0; l < nloc; ++l)
			{
				REAL* p = GetFreq(l);
				ushort* gcount = GetGenoCount(l);
				LOCSTAT1 &stat1 = loc_stat1[l];

				GENO_READER rt(ind0id, l);
				for (int i = 0; i < nind; ++i)
					gcount[rt.Read()]++; 

				int k2 = GetLoc(l).k;
				stat1.k = (ushort)CountK(p, k2);
				stat1.nhaplo = (int)(Sum(p, k2) + 0.5);
				Unify(p, k2);
				PROGRESS_VALUE += nind;
			}
		}
	}

	/* Calculate loc_stat2 in a pre-allocated buffer, used in relatedness estimation */
	template<typename REAL>
	TARGET void POP<REAL>::GetLocStat2(LOCSTAT2<REAL>* loc_stat2)
	{
		SetZero(relatedness_loc_stat, nloc);

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int64 l = 0; l < nloc; ++l)
		{
			REAL* p = GetFreq(l);
			LOCSTAT2<REAL>& stat2 = loc_stat2[l];
			int k2 = GetLoc(l).k, nhaplo = loc_stat1[l].nhaplo;

			stat2.s2 = stat2.s3 = stat2.s4 = 0;
			for (int a = 0; a < k2; ++a)
			{
				REAL af = p[a];
				if (af * nhaplo <= 1e-5) continue;
				stat2.s2 += af * af;
				stat2.s3 += af * af * af;
				stat2.s4 += af * af * af * af;
			}
		}
	}

	/* Is a argument a subpop of this */
	template<typename REAL>
	TARGET bool POP<REAL>::IsSubpop(POP<REAL>* subpop)
	{
		if (ispop) return false;
		for (int i = 0; i < npop; ++i)
		{
			if (vpop[i]->IsSubpop(subpop) || vpop[i] == subpop)
				return true;
		}
		return false;
	}

#endif

/* Initialize */
template<typename REAL>
TARGET void Initialize()
{
	if (g_eval_val == 1 && FileExists((g_output_val + ".eval.txt").c_str()))
		FileDelete((g_output_val + ".eval.txt").c_str());

	if (g_eval_val == 1)
	{
		FILE* f1 = fopen((g_output_val + ".eval.txt").c_str(), "wb");
		fputs("Time expense (s)\tAnalysis", f1);
		fputs(g_linebreak_val, f1);
		fclose(f1);
	}

	// Initialize variables
	individual_memory = NULL;
	locus_memory = NULL;
	conversion_memory = NULL;
	conversion_memory2 = NULL;
	conversion_string = NULL;

	clustering_memory = NULL;
	gd_tab = NULL;

	allele_freq_offset = NULL;
	maxK = 0;
	KT = 0;

	genotype_count_offset = NULL;
	maxG = 0;
	GT = 0;

	geno_bucket.offset.~LIST();
	ad_bucket.offset.~LIST();
	haplo_bucket.offset.~LIST();
	new(&geno_bucket)	BUCKET();
	new(&ad_bucket)		BUCKET();
	new(&haplo_bucket)	BUCKET();

	load_buf = NULL;
	vcf_header = NULL;

	pop<REAL>.Clear();
	for (uint i = 0; i < reg<REAL>.size; ++i)
		reg<REAL>[i].Clear();
	reg<REAL>.Clear();

	slocus = NULL;
	useslocus = false;

	ainds = NULL;
	apops = NULL;
	aregs = NULL;
	npop = 0;
	lreg = 0;
	SetZero(nreg, N_MAX_REG);
	nregt = 0;
	nregt2 = 0;

	total_pop = NULL;

	nloc = 0;
	nfilter = 0;
	nind = 0;
	nfilterind = 0;
	reassigned = 0;
	rinds = NULL;
	cpop = NULL;
	maxploidy = 0;
	minploidy = 0;

	progress1 = 0;
	progress2 = 0;
	state_lock = NULL;

	qslstack.Clear();
	haplotype_locus.Clear();

	convert_file = NULL;
	convert_buf = NULL;
	convert_linesize = 0;

	diversity_buf = NULL;
	SetZero(&diversity_sum<REAL>, 1);
	diversity_stage = 0;

	SetZero(fst_buf, 6);
	fst_type = 0;

	gdist_buf = NULL;
	gdist_type = 0;

	amova_buf = NULL;

	relatedness_buf = NULL;
	relatedness_loc_stat = NULL;

	kinship_buf = NULL;
	SetZero(gdindex, N_GD_ESTIMATOR + 1);

	clustering_matrix = NULL;
	structure_par = NULL;
	structure_totalruns = 0;

	spa_x = NULL;
	spa_n = 0;
	spa_np = 0;

	if (g_input_row * g_input_col == 0)
		Exit("\nError: no input files detected.\n");

	TOTLEN_DECOMPRESS = 0;
	TOTLEN_COMPRESS = 0;
	for (int i = 0; i < g_input_row; ++i)
	{
		for (int j = 0; j < g_input_col; ++j)
		{
			if ((FILE_INFO[i][j].handle = fopen(FILE_INFO[i][j].path.c_str(), "rb")) == NULL)
				Exit("\nError: input file %s does not exist.\n", FILE_INFO[i][j].path.c_str());

			FILEINFO& tf = FILE_INFO[i][j];
			path tpath(tf.path);
			tf.name = canonical(tpath).string();
			tf.compressed_len = GetFileLen((char*)FILE_INFO[i][j].path.c_str());
			TOTLEN_COMPRESS += tf.compressed_len;

			uint magic = FGetUshort(tf.handle);
			fclose(tf.handle);

			if (magic == 0x8B1F)//*.gz
			{
				g_format_val = -abs( g_format_val );
				if ((tf.handle = fopen(tf.path.c_str(), "rb")) == NULL)
					Exit("\nError: input file %s does not exist.\n", tf.path.c_str());
				fseeko64(tf.handle, 0, SEEK_END);
				tf.compressed_len = ftello64(tf.handle);
				fclose(tf.handle);

				tf.handle = FOpen(tf.path.c_str(), "rb");
			}
			else
			{
				if ((tf.handle = FOpen(tf.path.c_str(), "rb")) == NULL)
					Exit("\nError: input file %s does not exist.\n", tf.path.c_str());
				tf.compressed_len = GetFileLen((char*)tf.path.c_str());
				TOTLEN_DECOMPRESS += tf.compressed_len;
			}

			if (tf.handle == NULL)
				Exit("\nError: cannot open input file %s.\n", tf.path.c_str());
		}
	}

	/*
	if (g_format_val > BCF)
	{
		g_filebuf = new char[LINE_BUFFER];
		setvbuf(FILE_INFO[0][0].handle, g_filebuf, _IOFBF, LINE_BUFFER);
	}
	else
		g_filebuf = NULL;
	*/

	ALLELE_IDENTIFIER = abs( g_format_val ) <= BCF;

	FRES_BUF = new char[LINE_BUFFER];
	FRES_NAME = new char[PATH_LEN];
	omp_set_num_threads(g_nthread_val);

	InitBinomial();
	InitAlpha();
	InitCryptTable();

	if (g_input_row * g_input_col != 1 && abs( g_format_val ) > BCF)
		Exit("\nError: multiple input files are only supported for vcf format.\n");

	individual_memory = new MEMORY[1];
	locus_memory = new MEMORY[g_nthread_val];
	for (int i = 0; i < g_nthread_val; ++i)
		locus_memory[i].Alloc(1);
	for (int i = 0; i < g_nthread_val; ++i)
		locus_memory[i].Alloc(1);
	
	//initial missing genotypes
	SetFF(missing_array, N_MAX_PLOIDY);
	ushort* nullgatab = NULL;
	for (int i = 0; i <= N_MAX_PLOIDY; ++i)
	{
		missing_hash[i] = HashGenotype(missing_array, i);
		new(&missing_genotype[i]) GENOTYPE(nullgatab, missing_array, i);//sorted
	}

	POP<REAL> tp;
	LIST<POP<REAL>> tlist(NULL);
	pop<REAL>.Push(tp);
	char* defpopname = new char[7];
	strcpy(defpopname, "DefPop");
	new(&pop<REAL>[0]) POP<REAL>(defpopname, NULL, 0, 0, 0, 0, true);

	//assign pop 1
	if (g_indtext_b)
	{
		int indtextlen = (int)g_indtext_val.size();
		char* t1 = (char*)g_indtext_val.c_str(), *t2 = t1 - 1;
		char* indtextend = t1 + indtextlen;
		ReplaceChar(t1, '\r', '\n');
		ReplaceChar(t1, ' ', '\n');
		int rl = -2;
		while (*t1 && t1 < indtextend)
		{
			t1 = t2 + 1;
			while (*t1 == '\n' || *t1 == '\r' || *t1 == ' ' || *t1 == '\t') t1++;
			if (t1 >= indtextend) break;

			if (reg<REAL>.size == 0 || LwrLineCmp("#REG", t1) == 0)
			{
				//add reg
				rl++;
				if (reg<REAL>.size)
				{
					t1 = StrNextIdx(t1, "\n", 1) + 1;
					while (*t1 == '\n' || *t1 == '\r' || *t1 == ' ' || *t1 == '\t') t1++;
				}
				POP<REAL> tr2;
				reg<REAL>.Push(tlist);
				reg<REAL>[rl + 1].Push(tr2);
				char* defregname = new char[15];
				sprintf(defregname, "DefRegL%d", rl + 2);
				new(&reg<REAL>[rl + 1][0]) POP<REAL>(defregname, NULL, 0, 0, 0, 0, false);
			}

			if (rl == -1)
			{
				t2 = StrNextIdx(t1, ":", 1);
				if (t2 == NULL) break;
				*t2 = '\0';

				tp.name = new char[strlen(t1) + 1];
				strcpy(tp.name, t1);

				for (uint i = 0; i < pop<REAL>.size; ++i)
					if (strcmp(pop<REAL>[i].name, t1) == 0)
						Exit("\nError: Two populations have the same name: %s\n", t1);

				t1 = t2 + 1;
				t2 = StrNextIdx(t1, "\n", 1);
				if (t2) *t2 = '\0';
				else t2 = StrNextIdx0(t1, 1);

				tp.ispop = true;
				tp.id = pop<REAL>.size;
				tp.names = SplitStr(t1, ',', tp.nind);//deleted
				tp.inds = NULL;
				tp.allelefreq = NULL;
				tp.genocount = NULL;
				tp.loc_stat1 = NULL;
				tp.rid = 0;
				tp.vpop = NULL;
				tp.npop = 0;
				pop<REAL>.Push(tp);
			}
			else
			{
				t2 = StrNextIdx(t1, ":", 1);
				if (t2 == NULL) break;
				*t2 = '\0';

				POP<REAL> tr;
				tr.ispop = false;
				tr.id = reg<REAL>[rl].size;
				tr.name = new char[strlen(t1) + 1];
				strcpy(tr.name, t1);

				t1 = t2 + 1;
				t2 = StrNextIdx(t1, "\n", 1);
				if (t2) *t2 = '\0';
				else t2 = StrNextIdx0(t1, 1);

				tr.names = SplitStr(t1, ',', tr.npop);//deleted
				tr.nind = 0;
				tr.inds = NULL;
				tr.allelefreq = NULL;
				tr.genocount = NULL;
				tr.loc_stat1 = NULL;
				tr.vpop = NULL;
				tr.rid = 0;

				//check subreg/pop name
				for (int rl2 = 0; rl2 <= rl; ++rl2)
					for (uint i = 0; i < reg<REAL>[rl2].size; ++i)
						if (strcmp(reg<REAL>[rl2][i].name, tr.name) == 0)
							Exit("\nError: repeat population/region name: %s.\n", tr.name);

				for (uint i = 0; i < pop<REAL>.size; ++i)
					if (strcmp(pop<REAL>[i].name, tr.name) == 0)
						Exit("\nError: repeat population/region name: %s.\n", tr.name);

				//assign pop
				LIST<POP<REAL>>& popv = rl == 0 ? pop<REAL> : reg<REAL>[rl - 1];
				int npop2 = 0, namelen = 0;

				for (int j2 = 0; j2 < tr.npop; ++j2)
				{
					int v1 = 0, v2 = 0;
					ParseTwoNumber(tr.names[j2], v1, v2);

					// #1-#2 format
					if (v1 != -1 && v2 != -1)
					{
						if (v1 <= 0 || v1 >= popv.size)
							Exit("\nError: the first number of indtext %s for region %s exceeds range.", tr.names[j2], tr.name);
						if (v2 <= 0 || v2 >= popv.size)
							Exit("\nError: the second number of indtext %s for region %s exceeds range.", tr.names[j2], tr.name);

						int vmin = Min(v1, v2), vmax = Max(v1, v2);
						for (int i = vmin; i <= vmax; ++i)
						{
							if (popv[i].rid == tr.id)
								Exit("\nError: Population/subregion appear twice in Region %s.\n", popv[i].name, tr.name);
							if (popv[i].rid != 0)
								Exit("\nError: Two regions %s and %s have the same population/subregion: %s.\n", reg<REAL>[rl][popv[i].rid].name, tr.name, popv[i].name);
							namelen += (int)strlen(popv[i].name) + 1;
							popv[i].rid = (ushort)tr.id;
							npop2++;
						}
					}
					// #1 or id format
					else 
					{
						bool found = false;
						for (uint i = 0; i < popv.size; ++i)
							if (v1 == i || strcmp(popv[i].name, tr.names[j2]) == 0)
							{
								if (popv[i].rid == tr.id)
									Exit("\nError: Population/subregion appear twice in Region %s.\n", popv[i].name, tr.name);
								if (popv[i].rid != 0)
									Exit("\nError: Two regions %s and %s have the same population/subregion: %s.\n", reg<REAL>[rl][popv[i].rid].name, tr.name, popv[i].name);
								namelen += (int)strlen(popv[i].name) + 1;
								popv[i].rid = (ushort)tr.id;
								npop2++;
								found = true;
								break;
							}

						if (!found)
							Exit("\nError: Can not found population/subregion %s in region %s.\n", tr.names[j2], tr.name);
					}
				}

				delete[] tr.names;
				tr.names = NULL;
				tr.npop = npop2;
				reg<REAL>[rl].Push(tr);
			}
		}
	}

	NBUF = CALC_THREAD_BUFFER * g_nthread_val;
	state_lock = new atomic<int64>[NBUF];
	GLOCK3 = new LOCK[256];
	for (int i = 0; i < 256; ++i)
		InitLock(GLOCK3[i]);


	// check f_pop and slide_pop existence
	if (diversity_filter && f_pop_b)
	{
		if ("total" != f_pop_val)
		{
			bool find = false;
			for (int i = 0; i < npop; ++i)
				if (pop<REAL>[i].name == f_pop_val)
					find = true;

			if (!find)
			{
				for (uint rl = 0; rl < reg<REAL>.size - 1; ++rl)
					for (uint i = 0; i < reg<REAL>[rl].size; ++i)
						if (reg<REAL>[rl][i].name == f_pop_val)
							find = true;
				if (!find) Exit("\nError: Cannot find target population %d, check parameter -f_pop.\n", f_pop_val.c_str());
			}
		}
	}

	if (slide && slide_pop_b)
	{
		if ("total" != slide_pop_val)
		{
			bool find = false;
			for (int i = 0; i < npop; ++i)
				if (pop<REAL>[i].name == slide_pop_val)
					find = true;

			if (!find)
			{
				for (uint rl = 0; rl < reg<REAL>.size - 1; ++rl)
					for (uint i = 0; i < reg<REAL>[rl].size; ++i)
						if (reg<REAL>[rl][i].name == slide_pop_val)
							find = true;
				if (!find) Exit("\nError: Cannot find target population %d, check parameter -f_pop.\n", slide_pop_val.c_str());
			}
		}
	}
}

/* UnInitialize before Bayesian clustering*/
template<typename REAL>
TARGET void UnInitialize1()
{
	//Exit("\nCalculation complete!\n");

	vector<vector<FILEINFO>>().swap(FILE_INFO);
	g_input_row = g_input_col = 0;
	g_input_val = "";
	TOTLEN_DECOMPRESS = 0;

	if (lreg == -1 && npop == 1)
	{
		delete[] total_pop->name;
		total_pop->name = NULL;
	}

	for (uint rl = 0; rl < reg<REAL>.size; ++rl)
		for (uint i = 0; i < reg<REAL>[rl].size; ++i)
			reg<REAL>[rl][i].Uninitialize1();

	for (uint i = 0; i < pop<REAL>.size; ++i)
		pop<REAL>[i].Uninitialize1();

	qslstack.~LIST();           new (&qslstack) LIST<QUICKSORT_PARAMETER>();
	haplo_bucket.~BUCKET();		new (&haplo_bucket) BUCKET();
	ad_bucket.~BUCKET();		new (&ad_bucket) BUCKET();

	if (genotype_count_offset)  delete[] genotype_count_offset;    genotype_count_offset = NULL;
	if (cryptTable)				delete[] cryptTable;				cryptTable = NULL;
    
    if (GLOCK3)
    {
        for (int i = 0; i < 256; ++i)
            UnInitLock(GLOCK3[i]);
                                delete[] GLOCK3;					GLOCK3 = NULL;
    }
}

/* Final unInitialize */
template<typename REAL>
TARGET void UnInitialize2()
{
	if (ainds || rinds)
	{
		if (ainds) delete[] ainds;
		if (rinds && ainds != rinds) delete[] rinds;
		ainds = rinds = NULL;
	}

	for (uint rl = 0; rl < reg<REAL>.size; ++rl)
		for (uint i = 0; i < reg<REAL>[rl].size; ++i)
			reg<REAL>[rl][i].Uninitialize2();

	if (aregs)
	{
		for (uint rl = 0; rl < reg<REAL>.size; ++rl)
			if (aregs[rl])
			{
				delete[] aregs[rl];
				aregs[rl] = NULL;
			}
		delete[] aregs; aregs = NULL;
	}

	for (uint i = 0; i < pop<REAL>.size; ++i)
		pop<REAL>[i].Uninitialize2();

	if (apops) delete[] apops; apops = NULL;

	for (uint rl = 0; rl < reg<REAL>.size; ++rl)
		reg<REAL>[rl].~LIST();
	reg<REAL>.~LIST();			new (&reg<REAL>) LIST<LIST<POP<REAL>>>();
	pop<REAL>.~LIST();			new (&pop<REAL>) LIST<POP<REAL>>();
	geno_bucket.~BUCKET();		new (&geno_bucket) BUCKET();

	if (slocus)					delete[] slocus;					slocus = NULL;
	if (allele_freq_offset)     delete[] allele_freq_offset;		allele_freq_offset = NULL;
	if (state_lock)				delete[] state_lock;				state_lock = NULL;
	if (locus_memory)			delete[] locus_memory;				locus_memory = NULL;
	if (individual_memory)		delete[] individual_memory;			individual_memory = NULL;
	if (FRES_BUF)				delete[] FRES_BUF;			FRES_BUF = NULL;
	if (FRES_NAME)				delete[] FRES_NAME;			FRES_NAME = NULL;
}

/* Calculate individual mininum and maximum ploidy, and sum ploidy levels */
template<typename REAL>
TARGET void AssignPloidy()
{
	EvaluationBegin();
	RunThreads(&AssignPloidyThread<REAL>, NULL, NULL, nloc, nloc,
		"\nAssigning individual ploidy:\n", 1, true);

	//Calculate the total number of haplotypes in each population
	for (int i = 0; i < npop; ++i)
		apops[i]->nhaplotypes = 0;

	for (int i = 0; i < nind; ++i)
		apops[ainds[i]->popid]->nhaplotypes += ainds[i]->vmax;

	for (int rl = 0; rl <= lreg; ++rl)
		for (int i = 0; i < nreg[rl]; ++i)
		{
			POP<REAL>* tp = aregs[rl][i];
			tp->nhaplotypes = 0;
			for (int j = 0; j < tp->npop; ++j)
				tp->nhaplotypes += tp->vpop[j]->nhaplotypes;
		}

	reassigned = true;
	EvaluationEnd("Assign individual ploidy");
}

/* Calculate allele frequencies for each population and region for further use */
template<typename REAL>
TARGET void CalcFreq()
{
	EvaluationBegin();
	//Calculate allele frequency and genotype count
	RunThreads(&CalcFreqThread<REAL>, NULL, NULL, nloc * (int64)nind * (lreg + 2), nloc * (int64)nind * (lreg + 2),
		"\nCalculating allele frequencies:\n", 1, true);
	EvaluationEnd("Allele frequency estimation");
	CheckGenotypeId<REAL>();
}

/* Estimate PES model */
template<typename REAL>
TARGET void EstimatePES()
{
	if ((diversity && diversity_model_val[4] == 1) || (indstat && indstat_model_val[4] == 1) || (popas && popas_model_val[4] == 1))
	{
		EvaluationBegin();
		RunThreads(&GetLocusPESModel<REAL>, NULL, NULL, nloc, nloc, "\nEstimating double-reduction rate for PES model:\n", 1, true);
		EvaluationEnd("PES model estimation");
	}
}

/* Calculate various analyses */
template<typename REAL>
TARGET void Calculate()
{
	Initialize<REAL>();

	// 1. Load
	LoadFile<REAL>();

	// 2. Filter
	ApplyFilter<REAL>();//allow ad, locpos

	EvaluationStat();

	// 3. Sliding window
	CalcSlide<REAL>(); //forbid ad, locpos

	// 4. Haplotype
	CalcHaplotype<REAL>();//forbid ad, locpos

	AssignPloidy<REAL>();//allow ad

	if (diversity || indstat || fst || gdist || 
		amova || popas || relatedness || kinship || 
		pcoa || cluster || spa || ploidyinfer ||
		(structure && structure_f_val == 1))
	CalcFreq<REAL>();//allow ad

	//    Estimate pes model
	EstimatePES<REAL>();//allow ad

	// 5. Convert
	ConvertFile<REAL>(); //forbid ad, circle buf, locpos

	if (locus_pos)
	{
		delete[] locus_pos;
		locus_pos = NULL;
	}

	// 6. Genetic diversity
	CalcDiversity<REAL>(); //allow ad, circle buf

	// 7. Individual statistics
	CalcIndstat<REAL>(); //forbid ad, circle buf

	// 8. Genetic differentiation
	CalcDiff<REAL>(); //allow ad

	// 9. Genetic distance
	CalcDist<REAL>(); //allow ad, circle buf

	// 10. AMOVA
	CalcAMOVA<REAL>(); //forbid ad

	// 11. Population assignment
	CalcAssignment<REAL>(); //forbid ad

	// 12. Relatedness coefficient
	CalcRelatedness<REAL>(); //forbid ad, circle buf

	// 13. Kinship coefficient
	CalcKinship<REAL>(); //forbid ad, circle buf

	// 14. Principal coordinate analysis
	CalcPCOA<REAL>(); //allow ad

	// 15. Hierarchical clustering
	CalcClustering<REAL>(); //allow ad

	//TEST SPA
	CalcSPA<REAL>();//allow ad

	// 16. Ploidy Inference
	CalcPloidyInference<REAL>(); //forbid ad

	// Free allele freq
	UnInitialize1<REAL>();

	// 17. Bayesian clustering
	if (structure_eval_val == 1)
	{
		structure_nadmburnin_val = 0;
		structure_nburnin_val = 0;
		structure_nreps_val = 1;
		structure_nthinning_val = 1;
		structure_nruns_val = 1;
		structure_krange_min = 5;
		structure_krange_max = 5;
		g_nthread_val = 4;
		structure_nstream_val = 1;

#ifdef CUDA
			g_gpu_val = 2;
			nGPU = GetDeviceCountCUDA();
#else
			g_gpu_val = 1;
			nGPU = 0;
#endif

		for (structure_admix_val = 1; structure_admix_val <= 2; ++structure_admix_val)
		{
			g_fastsingle_val = 1;
			CalcBayesian<float >(); //forbid ad

			g_fastsingle_val = 2;
			CalcBayesian<float >(); //forbid ad

			CalcBayesian<double>(); //forbid ad
		}
	}
	else
		CalcBayesian<REAL>(); //forbid ad

	// Free allele freq
	UnInitialize2<REAL>();
}

/* Calculate allele frequencies for each population and region */
THREAD2(CalcFreqThread)
{
	if (allele_freq_offset)      delete[] allele_freq_offset;
	if (genotype_count_offset)   delete[] genotype_count_offset;

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
		maxK = Max((int)GetLoc(l).k, maxK);
		maxG = Max((int)GetLoc(l).ngeno, maxG);
	}

	//Allocate memory for allele frequency and genotype count for each population and region
	for (int i = 0; i < npop; ++i)
		if (apops[i]->nind && ad == 0)
			apops[i]->AllocFreq();

	for(int rl = 0; rl < lreg; ++rl)
		for (int i = 0; i < nreg[rl]; ++i)
				aregs[rl][i]->AllocFreq();

	if (lreg > -1) total_pop->AllocFreq();

	if (ad)
	{
		for (int rl = 0; rl < lreg; ++rl)
			for (int i = 0; i < nreg[rl]; ++i)
			{
				POP<REAL>* cp = aregs[rl][i];
				POP<REAL>* *vp = cp->vpop;
				for (int p = 0; p < cp->npop; ++p)
					Add(cp->allelefreq, vp[p]->allelefreq, KT);
			}

		if (lreg > -1)
		{
			POP<REAL>* cp = total_pop;
			POP<REAL>* *vp = cp->vpop;
			for (int p = 0; p < cp->npop; ++p)
				Add(cp->allelefreq, vp[p]->allelefreq, KT);
		}
	}

	for (int i = 0; i < npop; ++i)
		apops[i]->CalcFreqGcount();

	for (int rl = 0; rl < lreg; ++rl)
		for (int i = 0; i < nreg[rl]; ++i)
			aregs[rl][i]->CalcFreqGcount();

	if (lreg > -1) total_pop->CalcFreqGcount();
}

/* Calculate individual minimum and maximum ploidy, and sum ploidy levels */
THREAD2(AssignPloidyThread)
{
	int* Ploidy = new int[maxG * g_nthread_val];
	int* Nalleles = new int[maxG * g_nthread_val];
	int64* Vt = new int64[nind * g_nthread_val];
	byte* Vmin = new byte[nind * g_nthread_val];
	byte* Vmax = new byte[nind * g_nthread_val];

	SetZero(Vt, nind * g_nthread_val);
	SetVal(Vmin, (byte)100, nind * g_nthread_val);
	SetZero(Vmax, nind * g_nthread_val);

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
	for (int64 l = 0; l < nloc; ++l)
	{
		threadid = omp_get_thread_num();

		int* nalleles = Nalleles + maxG * threadid;
		int* ploidy = Ploidy + maxG * threadid;
		int64* vt = Vt + nind * threadid;
		byte* vmin = Vmin + nind * threadid;
		byte* vmax = Vmax + nind * threadid;

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
				vmin[j] = Min(vmin[j], v);
				vmax[j] = Max(vmax[j], v);
			}
		}
		PROGRESS_VALUE++;
	}

	//gather
	for (int tid = 1; tid < g_nthread_val; ++tid)
	{
		int64* vt = Vt + nind * tid;
		byte* vmin = Vmin + nind * tid;
		byte* vmax = Vmax + nind * tid;

		for (int j = 0; j < nind; ++j)
		{
			Vt[j] += vt[j];
			Vmin[j] = Min(Vmin[j], vmin[j]);
			Vmax[j] = Max(Vmax[j], vmax[j]);
		}
	}

	minploidy = 100; maxploidy = 0; maxvt = 0;
	for (int i = 0; i < nind; ++i)
	{
		ainds[i]->vt = Vt[i];
		ainds[i]->vmin = Vmin[i] == 100 ? 0 : (byte)Vmin[i];
		ainds[i]->vmax = Vmax[i];

		minploidy = Min(minploidy, ainds[i]->vmin);
		maxploidy = Max(maxploidy, ainds[i]->vmax);
		maxvt = Max(maxvt, ainds[i]->vt);
		sumvt += ainds[i]->vt;
	}

	delete[] Vt;
	delete[] Vmin;
	delete[] Vmax;
	delete[] Ploidy;
	delete[] Nalleles;
}

/* Find the optimal Partial Equational Segregation model for each locus */
THREAD2(GetLocusPESModel)
{
#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
	for (int64 l = 0; l < nloc; ++l)
	{
		GENOTYPE* gtab = GetLoc(l).GetGtab();
		int ngeno = GetLoc(l).ngeno;

		int model = 0;
		double maxli = -DBL_MAX;
		for (int m = 4; m <= N_DRE_MODELT; ++m)
		{
			double slogd = 0, prod = 1; 
			int64& slog = *(int64*) & slogd;

			OpenLog(slog, prod);
			for (int p = 0; p < npop; ++p)
			{
				ushort* gcount = cpop->GetGenoCount(l);
				REAL* freq = apops[p]->GetFreq(l);

				for (int gi = 0; gi < ngeno; ++gi)
				{
					GENOTYPE& gt = gtab[gi];
					if (gt.Nalleles() == 0 || gcount[gi] == 0) continue;

					double gfz = gt.GFZ(m, freq);
					gfz = gfz > 0 ? gfz : 1;

					if (gcount[gi] < 10)
						for (int ci = 0; ci < gcount[gi]; ++ci)
							ChargeLog(slog, prod, gfz);
					else
						slog += MyLog(gcount[gi]) * gcount[gi];
				}
			}
			CloseLog(slog, prod);

			if (slogd > maxli)
			{
				maxli = slog;
				model = m;
			}
		}

		GetLoc(l).pes_model = (byte)model;
		VLA_DELETE(gtab);

		PROGRESS_VALUE++;
	}
}
