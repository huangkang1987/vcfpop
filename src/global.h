/* Constants and Global variables */

#pragma once
#include "vcfpop.h"

struct OFFSET;
struct LOADLINE;
struct FILEINFO;

#ifdef _WIN64
#define _thread __declspec(thread)
#else 
#define _thread __thread
#endif

#define VERSION						("1.07b")
#define DATE						("2023-5-15")

#define N_FUNC						(18)						//Number of main options
#define N_MAX_OPTION				(35)						//Maximum options
#define N_MAX_SLIDEPLOT				(5)							//Maximum number of circles for sliding window plot
#define N_MAX_REG					(10)						//Maximum region levels
#define N_MAX_GDTAB					(20)						//Maximum number of genotypes a locus have to buffer individual genetic distance
#define N_CONVERTER					(9)							//Number of file converters
#define N_GD_ESTIMATOR				(24)						//Number of genetic distance estimators
#define N_INDGD						(11)						//Number of genetic distances in individual gd structure
#define N_CLUSTER_METHOD			(7)							//Number of hierarchical clustering estimators
#define N_FST_ESTIMATOR				(8)							//Number of Fst estimators
#define N_RELATEDNESS_ESTIMATOR 	(16)						//Number of relatedness estimators
#define N_KINSHIP_ESTIMATOR			(3)							//Number of kinship estimators
#define N_DRE_MODEL					(4)							//Number of double-reductio models
#define N_DRE_MODELT				(104)						//Number of double-reductio models with additional 1% gridents for PES models
#define N_MAX_PLOIDY				(10)						//Maximum ploidy level supported
#define N_MAXGAMMALN				(1024)						//Number of gammaln bufferred 
#define N_PATTERN_END				(139)						//Number of patterns, for pid>=139, pid-139 is ploidy level

#define M_PI						(3.14159265358979323846)
#define TWODIVPISQ2					(0.900316316157106)
#define NZERO						(-1e-10)
#define LIKELIHOOD_TERM				(1e-8)						//Threshold of likelihood difference to terminate iteration
#define MAX_ITER_DOWNHILL			(1200)						//Max iterations for down-hill simplex
#define MAX_ITER_SVD				(60)						//Max iterations for SVD decompsition
#define MAX_ALLELE_BYTE				(3)							//Max number of bytes an allele used in converting files
#define EPSILON_SVD					(1e-30)						//Minimum positive value for SVD decompsition
#define DOUBLE_UNDERFLOW			(1e-100)					//Lower bound of real number
#define DOUBLE_OVERFLOW				(1e100)						//Upper bound of real number
#define SINGLE_UNDERFLOW			(1e-10)						//Lower bound of real number
#define SINGLE_OVERFLOW				(1e10)						//Upper bound of real number
#define MIN_FREQ					(1e-10)						//Minimum allele freq to avoid being zero in unify

#define PATH_LEN					(8192)						//Buffer size for file name
#define NAME_BUF_LEN				(2048)						//Buffer size for locus name
#define IND_NAME_LEN				(256)						//Buffer size for individual identifier
#define CALC_THREAD_BUFFER			(512)						//Buffer size for kinship, relatedness, distance threads
#define LINE_BUFFER					(1048576)					//Buffer size for a line
#define INCREASE_BUFFER				(4096)						//Buffer size can be increased
#define SLEEP_TIME_TINY				(1)							//A small sleep time to wait other threads

#define RNG_SALT_SHUFFLE			0x8E3B0DAEFB191075
#define RNG_SALT_INITADM			0xDA310D821043E1D3
#define RNG_SALT_INITFIX			0x2C5E436525D873B5
#define RNG_SALT_UPDATEQ			0xE7853EB09C92935A
#define RNG_SALT_UPDATEZ			0x9B83516C63E249C3
#define RNG_SALT_UPDATEP			0xD4EA1DB746DF14E3
#define RNG_SALT_UPDATEA			0x49913AD9EDBD2EB8
#define RNG_SALT_UPDATELoc			0x95EE7991226E2D14
#define RNG_SALT_UPDATEAlpha		0x8A26E2C4D1BC5E54
#define RNG_SALT_UPDATELambda		0x67099E2FE2510092
#define RNG_SALT_UPDATEF			0xF82104A5ADDAAF9F
#define RNG_SALT_AMOVAHOMO			0xD1028FE10E8C369B
#define RNG_SALT_AMOVAANEU			0xC409B844103CAFE8
#define RNG_SALT_AMOVAML			0x03ED109FD0917108
#define RNG_SALT_CLUSTERING			0xC1BE4D107D3A8DA9
#define RNG_SALT_SIMDBENCHMARK		0x2EC7FCD12D417E90
#define RNG_SALT_CONVERT			0xFB768FFE7CB91E61


enum FILE_FORMAT
{
	VCF = 1, BCF, GENEPOP, SPAGEDI, CERVUS, ARLEQUIN, STRUCTURE, POLYGENE, POLYRELATEDNESS, GENODIVE, PLINK
};


/* Global Variables */
extern timepoint EVAL_BEGIN;
extern LOCK GLOCK1, GLOCK2, *GLOCK3;
extern bool ALLELE_IDENTIFIER;									//No allele name string for non-vcf/bcf input and vcf/bcf that performed haplotype extraction
extern int SIMD_TYPE;											//Single-instruction-multiple-data instruction used
extern int GDIST_METHOD;										//Current GD method, 1 genetic distance, 2 pcoa or 3 hierarchical clustering

extern string EXEDIR;											//Executable directory
extern string CURDIR;											//Current directory
extern string OUTDIR;											//Output directory
extern string OUTFILE;											//Output file

extern double ALPHA[N_DRE_MODELT + 1][N_MAX_PLOIDY + 1][3];		//Double reduction rates
extern double BINOMIAL[N_MAX_PLOIDY + 1][N_MAX_PLOIDY + 1];		//Binomial coefficients

extern atomic<int64> PROGRESS_VALUE;							//Progress value 
extern int64 PROGRESS_TOTAL;									//Total tasks in this function
extern int64 PROGRESS_CSTART;									//Start value of this batch
extern int64 PROGRESS_CEND;										//End value of this batch
extern int64 PROGRESS_NOUTPUTED;								//Number of outputed characters
extern int64 PROGRESS_NOUTPUTED2;								//Number of previously outputed characters

extern const char* GD_ESTIMATOR[];								//Genetic distance estimator names
extern const char* CLUSTER_METHOD[];							//Hierarchical method names
extern const char* FST_ESTIMATOR[];								//Fst estimator names
extern const char* RELATEDNESS_ESTIMATOR[];						//Relatedness estimator names
extern const char* KINSHIP_ESTIMATOR[];							//Kinship estimator names
extern const char* DRE_MODEL[];									//Double-reduction model names


/* Virtual Memory */

extern bool BIG_FILE;											//use > 10Gib genotype id memory
extern atomic<int64> VMEM_SIZE;									//number of bytes allocated
extern atomic<int64> VMEM_NBLOCK;								//number of virtual memory blocks allocated


/* File information and length */

extern vector<vector<FILEINFO>> FILE_INFO;
extern int64 TOTLEN_DECOMPRESS;
extern int64 TOTLEN_COMPRESS;


/* Temple and results file */

extern char* FRES_NAME;											//Result file name
extern char* FRES_BUF;											//Result buffer
extern tm* FRES_TIME;											//Result time
extern FILE* FRES;												//Result file pointers
extern FILE** TEMP_FILES;										//Temporatory file pointers
extern char** TEMP_BUFS;										//Temporatory buffer
extern char** TEMP_NAMES;										//Temporatory file names


/* Misc */

extern uint* cryptTable;										//crypt table for calculate hash for genotypes
extern int64 progress1, progress2;								//Progress value for read / write in a circle buffer
extern atomic<int64>* state_lock;								//State of circle buffer
extern int NBUF;												//CALC_THREAD_BUFFER * g_nthread_val;
extern _thread int threadid;									//Thread index


/* Population info */

extern void* ainds_;											//Individuals
#define ainds (*(IND<REAL>***)&ainds_)

extern void* apops_;											//Rearranged populations
#define apops (*(POP<REAL>***)&apops_)

extern void* aregs_;											//Rearranged regions
#define aregs (*(POP<REAL>****)&aregs_)

extern void* total_pop_;										//Total population
#define total_pop (*(POP<REAL>**)&total_pop_)

extern int npop;												//Number of populations
extern int lreg;												//Level of regions
extern int nreg[N_MAX_REG];										//Number of regions in each level
extern int nregt;												//Total number of regions
extern int nregt2;												//Total number of region pairs across levels

extern LOCN nloc;												//Number of loci
extern int nind;												//Number of individuals
extern byte maxploidy;											//Max ploidy in all genotypes
extern byte minploidy;											//Min ploidy in all genotypes
extern int64 maxvt;												//Max number of allele copies in an indivdiaul
extern int64 sumvt;												//Total number of allele copies in all indivdiauls


/* GPU */

extern int nGPU;												//Number of GPU devices

/* Memory offset for compat bit-wise storage */

struct OFFSET
{
	uint64 offset : 48;											//Offset, some base pointer add this offset can obtain data at a locus
	uint64 size : 16;											//Size in bits of each piece of data, use piece wise storage
};

struct LOADLINE
{
	uint64 size : 56;
	uint64 flag : 8;
	char data[1];
};

struct FILEINFO
{
	string path;
	string name;
	int64 decompressed_len;
	int64 compressed_len;
	FILE* handle;
	BCFHEADER* bcfheader;
	int64 decompressed_offset;
	int64 compressed_offset;
};