/* Bayesian clustering functions */

#pragma once
#include "vcfpop.h"

#define ClusterFreq(k)              (Freq   + (k) * KT)
#define ClusterMeanFreq(k)          (FreqM  + (k) * KT)
#define ClusterLocusFreq(k,l)       (Freq   + (k) * KT + allele_freq_offset[l])
#define ClusterAlleleFreq(k,l,a)    (*(Freq + (k) * KT + allele_freq_offset[l] + (a)))

#define AncestralFreqA              (FreqA)
#define AncestralFreqB              (FreqA +     KT)
#define AncestralFreqC              (FreqA + 2 * KT)
#define AncestralLocusFreqA(l)      (FreqA          + allele_freq_offset[l])
#define AncestralLocusFreqB(l)      (FreqA +     KT + allele_freq_offset[l])
#define AncestralLocusFreqC(l)      (FreqA + 2 * KT + allele_freq_offset[l])

struct STRUCTURE_RUNINFO;
template<typename REAL> struct SCLUSTER;
template<typename REAL> struct BAYESIAN;

#pragma pack(push, 8)

struct STRUCTURE_RUNINFO
{
	double MeanlnL;							//Mean lnLSetAlleleDepth
	double VarlnL;							//Var lnL
	double lnPD;							//Ln Prob of Data
    int k;									//Number of ancestral clusters
    int id;									//Taskid
    int rep;								//Replicate index
    atomic_flag flag;						//Is task been allocated to a thread
};

template<typename REAL>
struct SCLUSTER
{
    REAL* bucket;

    /* Get allele frequency array */
    TARGET REAL* GetFreq(int64 l);

    /* Get allele frequency */
    TARGET REAL GetFreq(int64 l, int allele);

    /* Set allele frequency pointer */
    TARGET void SetFreq(REAL* _bucket);
};

template<typename REAL>
struct BAYESIAN
{
    /* Parameter */
    int id;									//Taskid 
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
    int diversity;							//Output diversity
    bool inferlambda;						//Updated lambda in each iteration
    bool difflambda;						//Use separate lambda for each cluster

    /* Admix */
    double alpha;							//The initial alpha, the priori Dirichlet parameter of admixture proportions Q
    double stdalpha;						//Standard deviation of uniform priori distribution of alpha
    double maxalpha;						//Maximum of uniform priori distribution of alpha
    double alphapriora;						//One gamma priori distribution parameter
    double alphapriorb;						//The other gamma priori distribution parameter
    int metrofreq;							//Frequency of Metropolis-Hastings update of admixture proportions Q
    bool inferalpha;						//Update alpha in ADMIX model
    bool diffalpha;							//Use separate alpha for each cluster
    bool uniformalpha;						//Priori distribution for alpha, uniform or gamma alpha

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
    REAL* fbuf1;
    REAL* fbuf2;

    double* bufKthread;                     //Buffer with length K*#subthreads
    double* bufNK1;							//Buffer with length N*K*#subthreads, align 64 bytes
    double* bufNK2;							//Buffer with length N*K*#subthreads, align 64 bytes
    byte* bufNK1o;							//Buffer with length N*K*#subthreads
    byte* bufNK2o;							//Buffer with length N*K*#subthreads
    double* bufN1;							//Buffer with length N*#subthreads
    double* bufN2;							//Buffer with length N*#subthreads

    /* Fixed */
    REAL* Freq;							    //K*KT, current allele frequency for each cluster 
    REAL* FreqM;							//K*KT, mean allele frequency in sampled iterations
    REAL* FreqA;							//3*KT, ancestral allele frequency

    double* Lambda;							//K, Priori dirichlet parameter for each cluster
    int64* MiSum;							//N*K, MiSum added by Mi in each iteration
    ushort* Z;								//N, current culster of individual i
    int64* Mi;								//N*K*#subthreads, Mi[i,k] is the number of allele copies of individuali assigned to cluster k
    REAL* Q;								//N*K, Q[i,k] is the priori proportion of indidivual i's gene from cluster k
    int* Ni;								//K*KT, Ni[k, k2] is the number of allele copies in each cluster
    double* Alpha;							//K, the priori Dirichlet parameter of admixture proportions Q for each cluster

    /* F model, epsilon, ancestral population */
    double* f;								//K,   (1-F)/F
    double* F;								//K,    Fst

    /* LocPriori model */
    double* Eta;							//K, eta for each cluster
    double* Gamma;							//S*K
    REAL* Di;								//S*K, number of individual from s into k

    double* AlphaLocal;						//S*K
    double* SumAlpha;						//S

    /* CUDA */
    BAYESIAN<REAL>* bayes_CUDA;				//this in GPU memory
    REAL* Freq_CUDA;
    ushort* Z_CUDA;							//N, current culster of individual i
    int* Ni_CUDA;							//K*KT, Ni[k, k2] is the number of allele copies in each cluster
    int64* Mi_CUDA;							//N*K, Mi[i,k] is the number of allele copies of individuali assigned to cluster k
    REAL* Q_CUDA;							//N*K, Q[i,k] is the priori proportion of indidivual i's gene from cluster k
    double* bufNK1_CUDA;					//Buffer with length N*K
    double* bufNK2_CUDA;					//Buffer with length N*K
    double* bufN1_CUDA;						//Buffer with length N
    double* bufN2_CUDA;						//Buffer with length N
    REAL* FreqA_CUDA;
    double* f_CUDA;							//K,   (1-F)/F
    double* Lambda_CUDA;					//K, Priori dirichlet parameter for each cluster

    uint64 seed;

    /* Misc */
    STRUCTURE_RUNINFO* par2;				//Structure parameter
    int m;									//Current itation id
    int nr;									//Number of recorded records
    int nlnL;                               //Number of likelihoods recorded
    double* lnL;                            //Likelihoods recorded
    double* rout;							//Output buffer
    int64* kdis;							//kdis[k] = number of loci has k alleles
    int rlen;								//Rout size
    bool iscudastream;                      //Use GPU in this interation
    bool usecuda;                           //Use GPU in this interation
    bool binaryq;							//Is q only takes from 0 and 1
    bool singlez;							//noadmix, noadmixburnin, nolocipriori
    atomic<int64> l_atomic[32];             //loop variable avoid thread conflict

    /* Set all bits to 0 */
    TARGET BAYESIAN();

    /* Write results for a run */
    TARGET void PrintStructure();

    /* Write results summary for all runs */
    TARGET static void PrintSummary(STRUCTURE_RUNINFO* sp, int len);

    /* Copy parameters */
    TARGET void ReadPar(STRUCTURE_RUNINFO* _par2);

    /* Initialize MCMC */
    TARGET void InitFix(int tid = -1);

    /* Initialize MCMC for admix model */
    TARGET void InitAdmix(int tid = -1);

    /* Initialize MCMC for locpriori model */
    TARGET void InitLocPriori();

    /* Initialize MCMC for F model */
    TARGET void InitFmodel();

    /* Allocate and copy GPU memory */
    TARGET void ManageCUDA(bool isrelease);

    /* Allocate and copy GPU memory */
    TARGET void InitCUDA();

    /* Update allele frequency for all clusters */
    TARGET void UpdateP(int tid = -1);
    TARGET void UpdatePCUDA();

    /* Update a priori ancetral proportion for non-admix model */
    TARGET void UpdateQNoAdmix();

    /* Update a priori ancetral proportion for admix model */
    TARGET void UpdateQAdmix();

    /* Update a priori ancetral proportion by Metropolis-Hastings for admix model*/
    TARGET void UpdateQMetro();

    /* Update a priori ancetral proportion */
    TARGET void UpdateQ();

    /* Update locpriori parameters */
    TARGET void UpdateLocPriori();

    /* Update ancestral proportion for non-admix model */
    TARGET void UpdateZNoAdmix();

    /* Update ancestral proportion for admix model */
    TARGET void UpdateZAdmix();

    /* Update ancestral proportion for each allele or for each individual */
    TARGET void UpdateZ();

    /* Update Dirichlet parameter alpha (to draw admixture proportion Q) in the admix model */
    TARGET void UpdateAlpha();

    /* Update Dirichlet parameter lambda (to draw allele frequency) */
    TARGET void UpdateLambda();

    /* Update allele frequency of ancestral population for the F model */
    TARGET void UpdateA();
    TARGET void UpdateA1(int tid = -1);
    TARGET void UpdateA2(int tid = -1);

    /* Update population-specific Fst for the F model */
    TARGET void UpdateF();

    /* Finalize records */
    TARGET void Arrange();

    /* Record updated MCMC parameters */
    TARGET    void Record();

    template<bool isadmix = true, bool fast_fp32 = true>
    TARGET    void RecordCPU(int tid = -1);
    template<bool isadmix = true, bool fast_fp32 = true>
    TARGET512 void Record512(int tid = -1);
    template<bool isadmix = true, bool fast_fp32 = true>
    TARGETAVX void RecordAVX(int tid = -1);
    template<bool isadmix = true, bool fast_fp32 = true>
    TARGETSSE void RecordSSE(int tid = -1);
    template<bool isadmix = true, bool fast_fp32 = true>
    TARGETNEO void RecordNEO(int tid = -1);
    TARGET    void RecordCUDA();

    /* Free memory */
    TARGET void Uninit();

    /* Free GPU memory */
    TARGET void UninitCUDA();

    /* Perform MCMC */
    TARGET void MCMC();

    TARGET512 void UpdateQNoAdmix512(int tid = -1);
    TARGETAVX void UpdateQNoAdmixAVX(int tid = -1);
    TARGETSSE void UpdateQNoAdmixSSE(int tid = -1);
    TARGETNEO void UpdateQNoAdmixNEO(int tid = -1);
    TARGET    void UpdateQNoAdmixCPU(int tid = -1);
    TARGET    void UpdateQNoAdmixCUDA();

    template<bool fast_fp32 = true>
    TARGET512 void UpdateQMetro512(int tid = -1);
    template<bool fast_fp32 = true>
    TARGETAVX void UpdateQMetroAVX(int tid = -1);
    template<bool fast_fp32 = true>
    TARGETSSE void UpdateQMetroSSE(int tid = -1);
    template<bool fast_fp32 = true>
    TARGETNEO void UpdateQMetroNEO(int tid = -1);
    template<bool fast_fp32 = true>
    TARGET    void UpdateQMetroCPU(int tid = -1);
    TARGET    void UpdateQMetroCUDA();

    template<bool fast_fp32 = true>
    TARGET512 void UpdateZAdmix512(int tid = -1);
    template<bool fast_fp32 = true>
    TARGETAVX void UpdateZAdmixAVX(int tid = -1);
    template<bool fast_fp32 = true>
    TARGETSSE void UpdateZAdmixSSE(int tid = -1);
    template<bool fast_fp32 = true>
    TARGETNEO void UpdateZAdmixNEO(int tid = -1);
    template<bool fast_fp32 = true>
    TARGET    void UpdateZAdmixCPU(int tid = -1);
    TARGET    void UpdateZAdmixCUDA();

    TARGET512 void UpdateZNoAdmix512(int tid = -1);
    TARGETAVX void UpdateZNoAdmixAVX(int tid = -1);
    TARGETSSE void UpdateZNoAdmixSSE(int tid = -1);
    TARGETNEO void UpdateZNoAdmixNEO(int tid = -1);
    TARGET    void UpdateZNoAdmixCPU(int tid = -1);
    TARGET    void UpdateZNoAdmixCUDA();
};


#pragma pack(pop)

extern STRUCTURE_RUNINFO* structure_par;						//Genetic distance for pcoa and hierarchical clustering
extern int structure_totalruns;									//Total number of runs in Bayesian clustering
extern int structure_nsubthread;								//Total number of runs in Bayesian clustering
extern int structure_loc_size_min, structure_loc_size_max;      
extern int64* structure_loc_lend;
extern int64 structure_loc_coprime;
extern int64 structure_loc_coprime64[32];
extern int64* structure_loc_original_idx;
extern int structure_navailable_cuda;                           //number of free cuda devices
extern atomic<int>* structure_cuda_taskid;                      //current taskids using cuda
extern int structure_nsimult;                                   //number of task simultaneous run on a single cuda device

/* Quick sort locus by genotype bitsize */
TARGET void QSLocus2(int64 left, int64 right);

/* Quick sort locus by genotype bitsize */
THREADH(QSWorker2);

/* Calculate bayesian clustering */
template<typename REAL>
TARGET void CalcBayesian();

/* Calculate Bayesian clustering using multiple threads */
THREAD2H(BayesianThread);

/* Calculate Bayesian clustering using multiple threads */
THREAD2H(ArrangeLocus);
