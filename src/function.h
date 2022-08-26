/* Functions */

#pragma once
#include "vcfpop.h"

#pragma pack(push, 1)

/* Temporatory genotype, will be compressed into GENOTYPE */
struct TEMP_GENOTYPE
{
	HASH hash;								//Genotype hash
	int ploidy;								//Ploidy level
	int gid;								//Genotype id
	ushort alleles[N_MAX_PLOIDY];			//Indexed alleles
};

/* Header class of BCF file */
class BCFHEADER
{
public:
	int format_gtid;						//Index of GT field
	int format_gqid;						//Index of GQ field
	int format_dpid;						//Index of DP field
	int format_adid;						//Index of AD field
	int filter_passidx;						//Pass filter index
	char** contig_name;						//Config names
	int64 contig_size;						//Contig size
	int64 nsample;							//Number of sample

	//##FILTER=<ID=PASS,Description="All filters passed",IDX=0>
	//##contig=<ID=1,assembly=b37,length=249250621,IDX=0>
	//##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype",IDX=1>
	//##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants",IDX=2>

	/* Uninitialize BCF header */
	TARGET ~BCFHEADER();

	/* Initialize BCF header */
	TARGET BCFHEADER(char* header);
};

struct GENOTYPE
{
			//4 bytes
	uint offset : 24;						//this pointer + offset * 2 is address of alleles array
			//ushort* alleles;				//Alleles copy in ascending order [ploidy]
											//Unique alleles [nalleles] order by dosage descending
	uint patternid : 8;						//Pattern index

	/* Do nothing */
	TARGET GENOTYPE();

	/* Create genotype from alleles and hash */
	TARGET GENOTYPE(ushort*& gatab, ushort* alleles, int ploidy);

	/* Copy from a reference genotype */
	TARGET GENOTYPE(ushort*& gatab, GENOTYPE& ref);

	/* Get allele copy at ith haplotype */
	TARGET ushort GetAlleleCopy(int i);

	/* Get allele array */
	TARGET ushort* GetAlleleArray();

	/* Set allele array */
	TARGET void SetAlleleArray(ushort* alleles);

	/* Get pattern code */
	TARGET uint64 GetPattern();

	/* Number of alleles */
	TARGET int Nalleles();

	/* Ploidy level */
	TARGET int Ploidy();

	/* Crc32 hash */
	TARGET HASH Hash();

	/* Heterozygosity in this genotype */
	TARGET double HIndex();

	/* SS within genotype under IAM model */
	TARGET double SS_IAM();

	/* SS within genotype under IAM model */
	TARGET double SS_SMM(ushort* alen);

	/* The multinomial coefficient for HWE/RCS genotype frequency */
	TARGET double HWECoef();

	/* Genotypic frequency for zygotes under inbreeding */
	TARGET double GFZ(int* allele_count, int sum, double f);

	/* Genotypic frequency for zygotes under specific double-reduction model */
	TARGET double GFZ(int DR_MODE, double* f);

	/* Number of copies of target allele */
	TARGET int GetAlleleCount(int a);

	/* Frequency of target allele in this genotype */
	TARGET double GetFreq(int a);

	/* Frequencies of all alleles in this genotype */
	TARGET void GetFreq(double* p, int k2);

	/* Obtain Genepop genotype string */
	TARGET char* GetGenepopStr();

	/* Obtain Spagedi genotype string */
	TARGET char* GetSpagediStr();

	/* Obtain Cervus genotype string */
	TARGET char* GetCervusStr();

	/* Obtain Arlequin genotype string */
	TARGET char* GetArlequinStr();

	/* Obtain Arlequin genotype string */
	TARGET char* GetStructureStr();

	/* Obtain Polygene genotype string */
	TARGET char* GetPolygeneStr();

	/* Obtain PolyRelatedness genotype string */
	TARGET char* GetPolyRelatednessStr();

	/* Obtain GenoDive genotype string */
	TARGET char* GetGenoDiveStr();
};

class SLOCUS
{
	//small locus, 12 bytes

public:
	uint64 bits1 : 48;						//Genotype table[ngeno]
											//alen[k] {for non-vcf SMM distance}
											//chrom \0 name \0 {(allele identifiers \0)[k] for vcf/bcf)}
											//genotype alleles[gasize]
	uint64 k : 16;							//Number of alleles

	uint   ngeno : 22;						//Number of genotypes
	uint   flag_pass : 1;					//0 pass filter
	uint   flag_alen : 1;					//Has alen table?
	uint   pes_model : 8;					//1 for RCS, 2 for PRCS, 3 for CES, 4+ for PES

	/* Initialize */
	TARGET SLOCUS();

	/* Convert from LOCUS */
	TARGET SLOCUS(MEMORY& memory, LOCUS& ref);

	/* Create unphase locus */
	TARGET SLOCUS(MEMORY& memory, SLOCUS& ref, TABLE<HASH, uint>& gitab, ushort* gtmap);

	/* Create locus for haplotype extraction and Chi-square test */
	TARGET SLOCUS(MEMORY& memory, SLOCUS& ref, int64 _id, int _ngeno, int _gasize, TABLE<HASH, TEMP_GENOTYPE>& temptab);

	/* Deep copy from SLOCUS */
	TARGET SLOCUS(MEMORY& memory, SLOCUS& ref);

	/* Get Genotype array */
	TARGET GENOTYPE* GetGtab();

	/* Get end of chrom \0 name \0 {(allele identifiers \0)[k] */
	TARGET char* GetEnd();

	/* Get chrom string */
	TARGET char* GetChrom();

	/* Get locus identifier */
	TARGET char* GetName();

	/* Get locus identifier */
	TARGET char* GetNameStr(char* buf);

	/* Get allele name for vcf/bcf */
	TARGET char* GetAlleleName(int a);

	/* Get alen array */
	TARGET ushort* GetAlenArray();

	/* Get genotype allele array */
	TARGET ushort* GetGenoAlleleArray();

	/* Get Genotype alleles array size Sum(ploidy+nalleles)*/
	TARGET uint GetGenoAlleleSize();

	/* Get SMM distance */
	TARGET double GetSMMDist(int a, int b);
};

class LOCUS : public SLOCUS
{
	//+ 32 + 35 bytes
public:
	TABLE<HASH, GENOTYPE*> gftab;			//Genotype table, do not need release
	char* _chrom;							//Chrom \0 name \0 {(allele identifiers \0)[k] for vcf/bcf)}
	ushort* _alen;							//alen[k] {for non-vcf SMM distance}
	ushort* _genoallele;
	
	uint64 pos : 48;						//Position in chromosome	
	uint64 nploidy : 16;					//Number of ploidy levels
	LOCN id;								//Index
	uint maxdepth;							//Max allele depth
	uint gasize;							//Size of genotype alleles sum(ploidy + nalleles)

	ushort gtid;							//Index of GT tag
	ushort gqid;							//Index of GQ tag
	ushort dpid;							//Index of DP tag
	ushort adid;							//Index of AD tag
	ushort format_size;						//Format offset

	byte flag_original : 1;					//1 pass original filter
	byte flag_qual : 1;						//2 qual
	byte flag_indel : 1;					//3 isindel
	byte flag_hasgq : 1;					//4 hasgq
	byte flag_hasdp : 1;					//5 hasdp
	byte flag_hasad : 1;					//7 hasad

	/* Initialize */
	TARGET LOCUS();

	/* Deep Copy Locus, SLOCUS */
	TARGET LOCUS(MEMORY& memory, int64 _id, LOCUS& ref);

	/* Create unphase locus */
	TARGET LOCUS(MEMORY& memory, LOCUS& ref, TABLE<HASH, uint>& gitab, ushort* gtmap);

	/* Create locus for haplotype extraction and Chi-square test, SLOCUS */
	TARGET LOCUS(MEMORY& memory, LOCUS& ref, int64 _id, int _ngeno, int _gasize, TABLE<HASH, TEMP_GENOTYPE>& alstab);

	/* For dummy locus for collapse alleles during testing genotype distributions */
	TARGET LOCUS(MEMORY& memory, SLOCUS& ref);

	/* For non-vcf input, set locus name and id, SLOCUS */
	TARGET LOCUS(MEMORY& memory, char* line, int64 _id, int _ngenotype, GENOTYPE*& gtab, ushort*& gatab);

	/* For vcf input, set locus name and id, SLOCUS */
	TARGET LOCUS(MEMORY& memory, char*& line, uint64 _mask, int _ngenotype, GENOTYPE*& gtab, ushort*& gatab);

	/* Get end of chrom \0 name \0 {(allele identifiers \0)[k] */
	TARGET char* GetEnd();

	/* Get chrom string */
	TARGET char* GetChrom();

	/* Get locus identifier */
	TARGET char* GetName();

	/* Get locus identifier */
	TARGET char* GetNameStr(char* buf);

	/* Get allele name for vcf/bcf */
	TARGET char* GetAlleleName(int a);

	/* Get alen array */
	TARGET ushort* GetAlenArray();

	/* Get genotype allele array */
	TARGET ushort* GetGenoAlleleArray();

	/* Get Genotype alleles array size Sum(ploidy+nalleles)*/
	TARGET uint GetGenoAlleleSize();

	/* Get index of a target format */
	TARGET ushort GetFormatId(char* base, char* fmt_name, ushort* fmt_offset);
};

class GENO_READER
{
public:
	uint64* pos;							//Current read pointer
	uint64 data;							//Readed bits
	byte size;								//Number of bits a genotype id used
	byte nbits;								//Number of bits remaining in data

	/* Do nothing */
	TARGET GENO_READER();

	/* Initialize reader */
	TARGET GENO_READER(int indid, int64 l, BUCKET* bucket = NULL);

	/* Get id of next ind (order by indid) */
	TARGET int Read();
};

class GENO_WRITER
{
public:
	uint64* pos;							//Current read pointer
	uint64 data;							//Readed bits
	byte size;								//Number of bits a genotype id used
	byte nbits;								//Number of bits in data

	/* Do nothing */
	TARGET GENO_WRITER();

	/* Initialize writer */
	TARGET GENO_WRITER(int64 l, BUCKET* bucket = NULL);

	/* Write id of next ind to buffer */
	TARGET void Write(uint gtid);

	/* Write all remaining bits to buffer */
	TARGET void FinishWrite();
};

class IND
{
public:
	int indid;								//Individual id
	int popid;								//Population id
	char* name;								//Individual identifier
	int64 vt;								//Sum of allele copies across loci, missing genotype do not account
	byte vmax;								//Maximum ploidy level
	byte vmin;								//Minimum ploidy level

	/* Initialize */
	TARGET IND();

	/* Create individual for non-vcf input */
	TARGET IND(char* t, bool iscount, int id, GENOTYPE** gtab, ushort** gatab, GENO_WRITER* wt);

	/* Create individual for vcf input */
	TARGET IND(char*& title, int id);

	/* Create individual from a reference */
	TARGET IND(IND& ref);

	/* Unnitialize */
	TARGET ~IND();

	/* Set allele sequencing depth, for ad, test */
	TARGET static void SetAlleleDepth(int64 l, uint* depth, int K, int indid);

	/* Set allele sequencing depth, for ad, test */
	TARGET void SetAlleleDepth(int64 l, uint* depth, int K);

	/* Set allele sequencing depth, for ad, test */
	TARGET void SetAlleleDepth(int64 l, uint* depth, int K, OFFSET* _offset, byte* bucket);

	/* Set allele sequencing depth, for ad, test */
	TARGET void GetAlleleDepth(int64 l, uint* depth);

	/* Set individual genotype with default bucket */
	TARGET void SetGenotype(int64 l, uint gid);

	/* Set individual genotype with local bucket */
	TARGET void SetGenotype(int64 l, uint gid, OFFSET* _offset, byte* bucket);

	/* Get index for a pair of genotype */
	TARGET static void GetDyadGenotypeIdx(int& id1, int& id2, int64 l);

	/* Get individual genotype index from default table */
	TARGET int GetGenotypeId(int64 l, byte* bucket, OFFSET* offset);

	/* Get individual genotype index from default table */
	TARGET int GetGenotypeId(int64 l);

	/* Get individual genotype from default table */
	TARGET GENOTYPE& GetGenotype(int64 l);

	/* Get individual genotype from local table */
	TARGET GENOTYPE& GetGenotype(int64 l, TABLE<HASH, GENOTYPE*>& gftab);

	/* Get individual genotype from local table */
	TARGET GENOTYPE& GetGenotype(int64 l, GENOTYPE* gtab);

	/* Create individual from genepop */
	TARGET void genepop(char* t, bool iscount, GENOTYPE** gtab, ushort** gatab, GENO_WRITER* wt);

	/* Create individual from spagedi */
	TARGET void spagedi(char* t, bool iscount, GENOTYPE** gtab, ushort** gatab, GENO_WRITER* wt);

	/* Create individual from cervus */
	TARGET void cervus(char* t, bool iscount, GENOTYPE** gtab, ushort** gatab, GENO_WRITER* wt);

	/* Create individual from arlequin */
	TARGET void arlequin(char* t, bool iscount, GENOTYPE** gtab, ushort** gatab, GENO_WRITER* wt);

	/* Create individual from structure */
	TARGET void structure(char* t, bool iscount, GENOTYPE** gtab, ushort** gatab, GENO_WRITER* wt);

	/* Create individual from polyrelatedness */
	TARGET void polyrelatedness(char* t, bool iscount, GENOTYPE** gtab, ushort** gatab, GENO_WRITER* wt);

	/* Create individual from polygene */
	TARGET void polygene(char* t, bool iscount, GENOTYPE** gtab, ushort** gatab, GENO_WRITER* wt);

	/* Create individual from genodive */
	TARGET void genodive(char* t, bool iscount, GENOTYPE** gtab, ushort** gatab, GENO_WRITER* wt);

	/* Calculate the likelihood of genotype data */
	TARGET double GenoFreq(POP* grp, int model, int64 loc, double e);

	/* Calculate the individual kinship coefficient */
	TARGET void Theta(POP* grp, double& f_ritland, double& f_loiselle, double& f_weir, double& t_ritland, double& t_loiselle, double& t_weir, int64 loc = -1);

	/* Write header row for individual statistics */
	TARGET static void IndividualStatisticsHeader(FILE* fout);

	/* Write result row for individual statistics */
	TARGET void PrintIndividualStatistics(FILE* fout);

	/* Write header row for ploidy inference, test */
	TARGET static void PloidyInferenceHeader(FILE* fout);

	/* Write result row for ploidy inference, test */
	TARGET void PrintPloidyInference(FILE* fout);

	/* Calculate the likelihood for ploidy inference */
	TARGET void PloidyInferlnL3(map<int64, SPF>& depth, int v, double f0, double f1, double f2, double& l0, double& l1, double& l2);

	/* Calculate the likelihood for ploidy inference */
	TARGET double PloidyInferlnL(map<int64, SPF>& depth, int v, double f);

	/* Infer the ploidy level from allelic depth distribution */
	TARGET void PloidyInference(int v, double& lnL, double& f, map<int64, SPF>& depth);

	/* average genotypic frequency at a diallelic locus given m copies of A */
	TARGET double AvgG(int v, int m, double f);

	/* Write header row for population assignment */
	TARGET static void AssignmentHeader(FILE* fout);

	/* Write result row for population assignment */
	TARGET void PrintAssignment(FILE* fout);

	/* Read and set genotype from bcf input */
	TARGET void AddBCFGenotype(int64 l, char*& gtstr, char*& gqstr, char*& dpstr, char*& adstr, int vlen, int asize, int gqlen, int dplen, int adlen, uint*& depth, TABLE<HASH, uint>& gfid, GENOTYPE*& gtab, ushort*& gatab, GENO_WRITER& wt);

	/* Read and set genotype from vcf input */
	TARGET void AddVCFGenotype(char*& line, int64 l, uint*& depth, TABLE<HASH, uint>& gfid, GENOTYPE*& gtab, ushort*& gatab, GENO_WRITER& wt);

	/* Get tag value */
	TARGET char* GetTagValue(char* re, int tagid);
};

/* Statistics of locus */
struct LOCSTAT1
{
	ushort k;								//Number of alleles in target pop
	int nhaplo;								//Number of allele copies
};

/* Statistics of locus */
struct LOCSTAT2
{
	double s2;								//Sum(pi^2)
	double s3;								//Sum(pi^3)
	double s4;								//Sum(pi^4)
};

struct POP
{
	int id;									//Population id, will be rearrange using tree strucutre of total population
	char* name;								//Population name
	char** names;							//Individual names, used during parsing indtext
	int nind;								//Number of individuals
	int ind0id;								//Indid of the first individual, used to fast  genotype id iteration in this pop
	IND* *inds;								//Individuals

	LOCSTAT1* loc_stat1;					//Locus diversity and statistics
	double* allelefreq;						//+ allele_freq_offset[l] is the allele frequency array at locus l
	ushort* genocount;						//+ genotype_count_offset[l] is the genotype count array at locus l, consider missing genotype

	int rid;								//Index of the region it belongs

	POP* *vpop;								//Subpopulations
	int npop;								//Number of subpopulations
	int nhaplotypes;						//Number of haplotypes, used in amova homoploid method

	bool ispop;								//Is a population or a region

	/* Initialize */
	TARGET POP();

	/* Create a pop */
	TARGET POP(char* _name, char** _names, int _nind, int _regid, int _npop, int _id, bool _ispop);

	/* Get genotype count array */
	TARGET ushort* GetGenoCount(int64 l);

	/* Get genotype count of a genotype */
	TARGET ushort GetGenoCount(int64 l, int id);

	/* Get allele frequencies array */
	TARGET double* GetFreq(int64 l);

	/* Get allele frequency of an allele */
	TARGET double GetFreq(int64 l, int a);

	/* Uninitialize a pop */
	TARGET void Uninitialize();

	/* Uncllocate memory for locstat, allele frequency, genotype count */
	TARGET void UnAllocFreq();

	/* Allocate memory for locstat, allele frequency, genotype count */
	TARGET void AllocFreq();

	/* Move after filter monomorphic locus in bayesian clustering */
	TARGET void MoveFreq(LOCN* nafoffset, LOCN* ngcoffset);

	/* Clear memory for locstat, allele frequency, genotype count */
	TARGET void ClearFreqGcount();

	/* Calculate loc stat, allele frequency, genotype count */
	TARGET void CalcFreqGcount();

	/* Calculate loc_stat2 in a pre-allocated buffer, used in relatedness estimation */
	TARGET void GetLocStat2(LOCSTAT2* loc_stat2);
};

/* Initialize */
TARGET void Initialize();

/* UnInitialize */
TARGET void UnInitialize();

/* Calculate individual minimum and maximum ploidy, and sum ploidy levels */
TARGET void AssignPloidy();

/* Calculate allele frequencies for each population and region for further use */
TARGET void CalcFreq();

/* Estimate PES model */
TARGET void EstimatePES();

/* Calculate genetic diveristy indices */
TARGET void CalcDiversity();

/* Calculate various analyses */
TARGET void Calculate();

/* Calculate allele frequencies for each population and region */
THREADH(CalcAlleleFreq);

/* Calculate individual minimum and maximum ploidy, and sum ploidy levels */
THREADH(AssignPloidyThread);

/* Find the optimal Partial Equational Segregation model for each locus */
THREADH(GetLocusPESModel);

/* Functions */

/* Misc */
extern POP* cpop;									//Current population
extern ushort missing_array[N_MAX_PLOIDY];			//Allele array of the missing genotypes
extern GENOTYPE missing_genotype[N_MAX_PLOIDY + 1];	//Missing genotype at different ploidy level
extern HASH missing_hash[N_MAX_PLOIDY + 1];			//Hash of missing genotype

/* Global */
extern SLOCUS* slocus;								//SLocus information
extern bool useslocus;								//Use small locus
extern TABLE<HASH, GENOTYPE*> emptygftab;			//Empty genotype table


/* Load File */
extern LIST<POP> pop;								//Original input population
extern LIST<LIST<POP>> reg;							//Original input regions
extern LOCUS* locus;								//Locus information

extern MEMORY* individual_memory;					//Individual memory class
extern MEMORY* locus_memory;						//Locus memory class
extern MEMORY* nvcf_memory;							//Locus memory for first round counting ngeno
extern TABLE<HASH, uint>* nvcf_gfid;				//Hash table for counting ngeno 

/* Allele frequency and genotype count offset */
extern LOCN* allele_freq_offset;					//Allele frequency array at locus l
extern int maxK;									//Max number of alleles
extern int64 KT;									//Total number of alleles

extern LOCN* genotype_count_offset;					//Genotype count array at locus l
extern int maxG;									//Max number of genotypes at all loci
extern int64 GT;									//Total number of genotypes


/* Genotype bucket */

extern VMEMORY locus_list;
extern BUCKET geno_bucket, haplo_bucket, ad_bucket;

/* Reassign individuals and populations */
extern bool reassigned;								//Is ploidy assigned, to distiguish diversity filter and diversity estimation
extern IND** rinds;									//Rearranged individuals according to population source

#pragma pack(pop)