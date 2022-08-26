/* Bayesian clustering functions */

#pragma once
#include "vcfpop.h"

#pragma pack(push, 1)

struct STRUCTURE_RUNINFO
{
	atomic_flag flag;						//Is task been allocated to a thread
	int k;									//Number of ancestral clusters
	int id;									//Run index
	int rep;								//Replicate index
	double MeanlnL;							//Mean lnLSetAlleleDepth
	double VarlnL;							//Var lnL
	double lnPD;							//Ln Prob of Data
};

class BAYESIAN_READER
{
public:
	union
	{
		uint64 data64[8];					//Readed bits of eight individuals
#ifdef __aarch64__
		uint64x2_t data128[4];				//Readed bits of eight individuals
#else 
		union
		{
			__m128i data128[4];				//Readed bits of eight individuals
			__m256i data256[2];				//Readed bits of eight individuals
			__m512i data512;				//Readed bits of eight individuals
		};
#endif
	};
	uint64* pos;							//Current read pointer
	int nbits;								//Number of bits remaining in data

	/* Do nothing */
	TARGET BAYESIAN_READER();

	/* Initialize reader */
	TARGET BAYESIAN_READER(int indid);
};

class BAYESIAN_WRITER
{
public:
	union
	{
		uint64 data64[8];					//Writed bits of eight individuals
#ifdef __aarch64__
		uint64x2_t data128[4];				//Writed bits of two individuals
#else
		union
		{
			__m128i data128[4];				//Writed bits of eight individuals
			__m256i data256[2];				//Writed bits of eight individuals
			__m512i data512;				//Writed bits of eight individuals
		};
#endif
	};
	uint64* pos;							//Current read pointer
	int nbits;								//Number of bits in data

	/* Do nothing */
	TARGET BAYESIAN_WRITER();

	/* Initialize writer */
	TARGET BAYESIAN_WRITER(int indid);

	/* Write id of next ind to buffer */
	TARGET void Write(int size, int aid);

	/* Write id of next ind to buffer */
	TARGET void Write8(int size, uint64 aid[STRUCTURE_NPACK]);

	/* Write all remaining bits to buffer */
	TARGET void FinishWrite();

	/* Write all remaining bits to buffer */
	TARGET void FinishWrite8();
};

struct SCLUSTER
{
	double* bucket;							//Allele frequency bucket, KT elements

	/* Get allele frequency array */
	TARGET double* GetFreq(int64 l);

	/* Get allele frequency */
	TARGET double GetFreq(int64 l, int allele);

	/* Set allele frequency pointer */
	TARGET void SetFreq(double* _bucket);
};

class BAYESIAN
{
public:
	/* Parameter */
	int S;									//Number of sampling locations
	int K;									//Number of clusters
	int N;									//Number of individuals
	int64 L;								//Number of loci

	/* Model */
	bool admix;								//Admix model
	bool locpriori;							//Locpriori model
	bool fmodel;							//F model

	/* MCMC */
	int nburnin;							//Number of burnin iterations
	int nreps;								//Number of iteration after burnin
	int nthinning;							//Sampling interval
	int nruns;								//Number of independent runs for each K
	int nadmburnin;							//Number of admburnin iterations

	/* Misc */
	double lambda;							//Dirichlet parameter to update the allele frequencies 
	double stdlambda;						//Standard deviation of new lambda
	double maxlambda;						//Maximum of new lambda
	bool inferlambda;						//Updated lambda in each iteration
	bool difflambda;						//Use separate lambda for each cluster
	int diversity;							//Output diversity

	/* Admix */
	double alpha;							//The initial alpha, the priori Dirichlet parameter of admixture proportions Q
	bool inferalpha;						//Update alpha in ADMIX model
	bool diffalpha;							//Use separate alpha for each cluster
	bool uniformalpha;						//Priori distribution for alpha, uniform or gamma alpha
	double stdalpha;						//Standard deviation of uniform priori distribution of alpha
	double maxalpha;						//Maximum of uniform priori distribution of alpha
	double alphapriora;						//One gamma priori distribution parameter
	double alphapriorb;						//The other gamma priori distribution parameter
	int metrofreq;							//Frequency of Metropolis-Hastings update of admixture proportions Q

	/* Locpriori model */
	double r;								//Initial value of r, where r evaluates the informativeness of data for the sampling location
	double maxr;							//Maximum of new r, where r evaluates the informativeness of data for the sampling location
	double epsr;							//Max step value of new r
	double epseta;							//Max step value of new eta for the non-ADMIXTURE model
	double epsgamma;						//Max step value of new gamma for the non-ADMIXTURE model

	/* F model */
	double pmeanf;							//Priori mean F
	double pstdf;							//Priori standard deviation of F
	double stdf;							//Standard deviation of new F
	bool fsame;								//Use the same F in all clusters

	double** bufK;							//Buffer with length K
	double* bufNK1;							//Buffer with length N*K, align 64 bytes
	double* bufNK2;							//Buffer with length N*K, align 64 bytes
	byte* bufNK1o;						//Buffer with length N*K
	byte* bufNK2o;						//Buffer with length N*K
	double* bufN1;							//Buffer with length N
	double* bufN2;							//Buffer with length N


	/* Fixed */
	SCLUSTER* cluster;						//Clusters, K elements
	SCLUSTER* clusterb;						//Additioanl clusters to write temp allele freq, K elements
	double* Base;							//(2 * K + 3) * KT elements, cluster[K*KT], clusterb[K*KT], PA[KT], PA1[KT], PA2[KT] 

	double* Lambda;							//K, Priori dirichlet parameter for each cluster
	int64* MiSum;							//N*K, MiSum added by Mi in each iteration
	ushort* Z;								//N, current culster of individual i
	int64* Mi;								//N*K, Mi[i,k] is the number of allele copies of individuali assigned to cluster k
	double* Q;								//N*K, Q[i,k] is the priori proportion of indidivual i's gene from cluster k
	int* Ni;								//K*KT, Ni[k, k2] is the number of allele copies in each cluster
	double* Alpha;							//K, the priori Dirichlet parameter of admixture proportions Q for each cluster

	/* F model, epsilon, ancestral population */
	double* f;								//K,   (1-F)/F
	double* F;								//K,    Fst
	SCLUSTER PA, PA1, PA2;

	/* LocPriori model */
	double* Eta;							//K, eta for each cluster
	double* Gamma;							//S*K
	double* Di;								//S*K, number of individual from s into k

	double* AlphaLocal;						//S*K
	double* SumAlpha;						//S

	/* Misc */
	RNG rng;								//Random number generator
#ifndef __aarch64__
	RNGSSE rngSSE;							//Random number generator
	RNGAVX rngAVX;							//Random number generator
	RNG512 rng512;							//Random number generator
#else
	RNGNEO rngNEO;							//Random number generator
#endif

	int m;									//Current itation id
	bool binaryq;							//Is q only takes from 0 and 1
	bool singlez;							//noadmix, noadmixburnin, nolocipriori
	STRUCTURE_RUNINFO* par2;				//Structure parameter
	int nr;									//Number of recorded records
	double* rout;							//Output buffer
	int rlen;								//Rout size
	int* kdis;								//kdis[k] = number of loci has k alleles


	/* Set all bits to 0 */
	TARGET BAYESIAN();

	/* Write results for a run */
	TARGET void PrintStructure();

	/* Write results summary for all runs */
	TARGET static void PrintSummary(STRUCTURE_RUNINFO* sp, int len);

	/* Copy parameters */
	TARGET void ReadPar(STRUCTURE_RUNINFO* _par2);

	/* Initialize MCMC */
	TARGET void InitFix();

	/* Initialize MCMC for admix model */
	TARGET void InitAdmix();

	/* Initialize MCMC for locpriori model */
	TARGET void InitLocPriori();

	/* Initialize MCMC for F model */
	TARGET void InitFmodel();

	/* Update allele frequency for all clusters */
	TARGET void UpdateP();

	/* Update a priori ancetral proportion for non-admix model */
	TARGET void UpdateQNoAdmix();

	/* Update a priori ancetral proportion for non-admix model */
	TARGET void UpdateQNoAdmixDecompress();

	/* Update a priori ancetral proportion for admix model */
	TARGET void UpdateQAdmix();

	/* Update a priori ancetral proportion by Metropolis-Hastings for admix model*/
	TARGET void UpdateQMetro();

	/* Update a priori ancetral proportion */
	TARGET void UpdateQ();

	/* Update locpriori parameters */
	TARGET void UpdateLocPriori();

	/* Update ancestral proportion for each allele or for each individual */
	TARGET void UpdateZ();

	/* Update Dirichlet parameter alpha (to draw admixture proportion Q) in the admix model */
	TARGET void UpdateAlpha();

	/* Update Dirichlet parameter lambda (to draw allele frequency) */
	TARGET void UpdateLambda();

	/* Update allele frequency of ancestral population for the F model */
	TARGET void UpdatePA();

	/* Update population-specific Fst for the F model */
	TARGET void UpdateF();

	/* Finalize records */
	TARGET void Arrange();

	/* Record updated MCMC parameters */
	TARGET void Record();

	/* Free memory */
	TARGET void Uninit();

	/* Perform MCMC */
	TARGET void MCMC();
};

#pragma pack(pop)

extern STRUCTURE_RUNINFO* structure_par;						//Genetic distance for pcoa and hierarchical clustering
extern int structure_totalruns;									//Total number of runs in Bayesian clustering
extern uint64** structure_allele;								//compressed allele array for each individual
extern byte* structure_size;									//allele data size (in bits) for each locus
extern int64 structure_indnbytes;								//number of bytes used in decompoess model of structure

/* Calculate bayesian clustering */
TARGET void CalcBayesian();

/* Calculate Bayesian clustering using multiple threads */
THREADH(BayesianThread);
