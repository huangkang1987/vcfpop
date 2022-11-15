/* Constants and Global variables */

#include "vcfpop.h"

#define extern 

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

extern const char* GD_ESTIMATOR[] =								//Genetic distance estimator names
{ "", "Nei1972", "Cavalli-Sforza1967", "Reynolds1983", "Nei1983", "Euclidean", "Goldstein1995", "Nei1974", "Roger1972", "Slatkin_Nei1973", "Slatkin_Weir1984", "Slatkin_Hudson1992", "Slatkin_Slatkin1995", "Slatkin_Hedrick2005", "Slatkin_Jost2008", "Slatkin_Huang2021_homo", "Slatkin_Huang2021_aneu", "Reynolds_Nei1973", "Reynolds_Weir1984", "Reynolds_Hudson1992", "Reynolds_Slatkin1995", "Reynolds_Hedrick2005", "Reynolds_Jost2008", "Reynolds_Huang2021_homo", "Reynolds_Huang2021_aneu" };
extern const char* CLUSTER_METHOD[] = 							//Hierarchical method names
{ "", "NEAREST", "FURTHEST", "UPGMA", "WPGMA", "UPGMC", "WPGMC", "WARD" };
extern const char* FST_ESTIMATOR[] = 							//Fst estimator names
{ "", "Nei1973", "Weir1984", "Hudson1992", "Slatkin1995", "Hedrick2005", "Jost2008", "Huang2021_homo", "Huang2021_aneu" };
extern const char* RELATEDNESS_ESTIMATOR[] = 					//Relatedness estimator names
{ "", "Lynch1999", "Wang2002", "Thomas2010", "Li1993", "Queller1989", "Huang2016A", "Huang2016B", "Milligan2003", "Anderson2007", "Huang2014", "Huang2015", "Ritland1996_modified", "Loiselle1995_modified", "Ritland1996", "Loiselle1995", "Weir1996" };
extern const char* KINSHIP_ESTIMATOR[] = 						//Kinship estimator names
{ "", "Ritland1996", "Loiselle1995", "Weir1996" };
extern const char* DRE_MODEL[] = 								//Double-reduction model names
{ "", "rcs", "prcs", "ces", "pes" };


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

#undef extern 