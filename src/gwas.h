/* Genome-Wide Association Studies Functions */

#pragma once
#include "vcfpop.h"

template<typename REAL> struct GWAS;

#define GWAS_GG(i,j) Gb[(i) * m + (j)]
#define REF_MAT(x) new (&x)	   rmat (ref->x.memptr(), ref->x.n_rows, ref->x.n_cols, false, true);
#define REF_COL(x) new (&x)	   rcol (ref->x.memptr(), ref->x.n_elem, false, true);
#define REF_ROW(x) new (&x)	   rrow (ref->x.memptr(), ref->x.n_elem, false, true);
#define nan GWAS_NAN<REAL>

//#pragma pack(push, 1)

template<typename REAL>
struct RESULT_WALD
{
	REAL lnL_REML;
	REAL W;
	REAL mlogP;
	REAL P;

	TARGET void Write(FILE* fout);
};

template<typename REAL>
struct RESULT_LRT
{
	REAL lnL_ML;
	REAL LR;
	REAL mlogP;
	REAL P;

	TARGET void Write(FILE* fout);
};

template<typename REAL>
struct RESULT_SCORE
{
	REAL lnL0_ML;
	REAL LM;
	REAL mlogP;
	REAL P;

	TARGET void Write(FILE* fout);
};

template<typename REAL>
struct GWAS_NULL_MODEL
{
	double r;							//optimized parameter
	double sig;							//optimized parameter
	CPOINT xx_ml0;						//optimized parameter
	rmat ma0;							//solved slope
};

template<typename REAL>
struct GWAS
{
	int64 tid;							//Thread id
	int64 n;               				//Number of individuals
	int64 m;                			//Number of locus
	int64 cl;                			//Current locus
	
	int k;               				// = mQt.n_cols;
	int nk;               				// = n - k;
	int kX;               				//d.f. of environmental variables without intercept
	int kXI;							//d.f. of environmental variables and intercepts 
	int kXIG;							//d.f. of environmental variables, intercepts and slopes

	int nnumeric;          				//Number of numeric environmental variables
	int nnominal;           			//Number of nominal environmental variables
	int nresponse;          			//Number of response variable

	int nploidy;            			//Number of ploidies
	int nallele;            			//Number of alleles
	int64* noffset;						//Allele offset

	int ploidy[N_MAX_PLOIDY + 1];		//Ploidy levels
	int ploidyidx[N_MAX_PLOIDY+1];		//Index of ploidy levels

	byte* Gb;							//n*m		Share, Genotype array, GWAS_DEBUG
	rmat oG;							//n*m		Share, Genotype matrix in a batch
	int64 oG_offset0;					//          Allele offset that the first col of oG

	rmat oY;							//n*ky		Share, Response variable
	rmat oX;							//n*kx		Share, Environmental covariates
	rmat oI;							//n*ki		Specific, Intercepts
	umap<HASH, rmat> Isbase;
	vector<HASH> Ihashbase;				

	umap<HASH, rmat>* Is;				//			Specific, Intercepts in a batch
	vector<HASH>* Ihash;				//			Specific, Intercept hashs in a batch
	int64 oI_lst;						//          Locus id that the first elements of Ihash

	rmat oR;							//n*n		Share, Relatedness matrix
	rmat* tR;							//n*n		Share, temporatory relatedness matrix
	rmat dR;							//n*1		Share, Diagonal elements of R

	rmat mY;							//n*ky		Specific, Response variable
	rmat mX;							//n*kx		Specific, Environmental covariates
	rmat mI;							//n*ki		Specific, Intercepts
	rcol cR;							//n*1		Specific, Diagonal elements for R
	rmat mG;							//n*kg		Specific, Genotype matrix, expanding multiple alleles into multiple columns
	rmat mXI;							//n*kxi		Specific, [X I]
	rmat mXIG;							//n*kxig	Specific, [X I G]

	rmat mU;							//n*n		Specific, Eigen vectors for R
	rmat mU_stride;						//n*n		Specific, Eigen vectors for R
	rcol cV;							//n*1		Specific, Eigen values for R
	rmat mPt;							//n*1		Specific, P = Y' * U
	rmat mQt;							//n*kt		Specific, Q = U' * Xpg
	
	rcol cVr;							//n*1		Specific, Vr = r2 * V + I
	rcol cW01;							//n*1		Specific, Wij = V ^ i / (Vr ^ j);
	rcol cr;							//n*1		Specific, Vr = r2 * V + I
	rcol cW12;							//n*1		Specific, Wij = V ^ i / (Vr ^ j);

	rmat mE01;							//k*k		Specific, Eij = Q*Wij*Qt
	rmat mF01;							//k*1		Specific, Fij = Q*Wij*Pt
	rmat mJ01;							//1*1		Specific, Jij = P*Wij*Pt
	rmat mE12;							//k*k		Specific, Eij = Q*Wij*Qt
	rmat mF12;							//k*k		Specific, Eij = Q*Wij*Qt
	rmat mJ12;							//k*k		Specific, Eij = Q*Wij*Qt
	rmat mE23;							//k*k		Specific, Eij = Q*Wij*Qt
	rmat mF23;							//k*k		Specific, Eij = Q*Wij*Qt
	rmat mJ23;							//k*k		Specific, Eij = Q*Wij*Qt
	rmat iE01;							//k*k		Specific, inv(E01), also the covariance matrix of betas
	rmat iEF01;							//k*1		Specific, E01 \ F01 also the slope (betas)

	// Parameters
	double sig2_from_r;					//Profiled sig2 from r

	// Results
	FILE* fout;							//Output file handle
	umap<int, GWAS_NULL_MODEL<REAL>> lnL0_Res;//Result point for the null model inconsider genotypes

	int64* sample_size;
	int64* degrees_of_freedom;
	RESULT_WALD <REAL>* wald;
	RESULT_LRT  <REAL>* lrt;
	RESULT_SCORE<REAL>* score;

	/* Write result header */
	TARGET void WriteGWASHeader();

	/* Write row header */
	TARGET void WriteGWASRowHeader();

	/* Write a cell */
	TARGET void WriteGWASCell(int p, int pid, int id);

	/* Write results for a locus */
	TARGET void WriteLocus(int64 _l);

	/* Do nothing */
	TARGET GWAS();

	/* Destructor */
	TARGET ~GWAS();

	/* Clone */
	TARGET GWAS(GWAS<REAL>* ref, int _tid);

	/* Perform Score Test */
	TARGET void ScoreTest(CPOINT& xx_ml0, int resid, GWAS_NULL_MODEL<REAL>* null_model);

	/* Perform LRT Test */
	TARGET void LRTTest(CPOINT& xx_ml, CPOINT& xx_ml0, int resid);

	/* Perform Wald Test */
	TARGET void WaldTest(CPOINT& xx_reml, int resid);

	/* Prepare GWAS */
	TARGET void Prepare();

	/* Copy Is Ihash oG oI_lst */
	TARGET void CopyRef(GWAS<REAL>* ref);

	/* Perform GWAS for locus l */
	TARGET void CalcLocus(int64 _l);

	/* Read a batch of genotypes */
	TARGET void ReadBatch(int64 lst, int64 led);

	/* Set independent variables */
	TARGET void SetIndVar(rmat& _mX);

	/* Set Vr W01 E01 F01 J01 */
	TARGET void SetVar01(double r2);

	/* Set Vr W01 E01 F01 J01 W12 E12 F12 J12 */
	TARGET void SetVar12(double r2);

	/* Optimizer Warpper */
	TARGET static double lnL_reml_profile1(void* Param, CPOINT& xx, rmat& G, rmat& H);

	/* Optimizer Warpper */
	TARGET static double lnL_ml_profile1(void* Param, CPOINT& xx, rmat& G, rmat& H);

	/* Optimizer Warpper */
	TARGET static double lnL_reml_profile2(void* Param, CPOINT& xx, rmat& G, rmat& H);

	/* Optimizer Warpper */
	TARGET static double lnL_ml_profile2(void* Param, CPOINT& xx, rmat& G, rmat& H);

	/* LogLikelihood using REML criterion */
	TARGET double lnL_reml_profile1(double r, rmat* G, rmat* H, bool issquare);

	/* LogLikelihood using ML criterion */
	TARGET double lnL_ml_profile1(double r, rmat* G, rmat* H, bool issquare);

	/* LogLikelihood using REML criterion */
	TARGET double lnL_reml_profile2(double r, double sig, rmat* G, rmat* H, bool issquare);

	/* LogLikelihood using REML criterion */
	TARGET double lnL_ml_profile2(double r, double sig, rmat* G, rmat* H, bool issquare);

	/* LogLikelihood without profile */
	TARGET double lnL2(double r, double sig, rmat& ma, rmat* G, rmat* H);

	/* Variance-covariance matrix by OPG estimator */
	TARGET rmat GetV_OPG2(double sig, double tau, rmat& beta, rmat& X);

	/* Variance-covariance matrix by Hessian estimator */
	TARGET rmat GetV_Hessian2(double sig, double tau, rmat& X);
};

template<typename REAL>
struct GWAS_MEM
{
	// Used to share eigen-values and vectors of relatedness matrix
	// when missing individuals are different among loci or response variables

	HASH hash;						// hashit
	REAL* buf;						// buffer
	atomic<int> nref;				// number of references
	atomic<double> last_access;		// last access time
	atomic_flag prepared;			// calculating
	atomic_flag del;				// not in the tab or list

	/* Find the entry from the dictionary */
	TARGET static GWAS_MEM<REAL>& FindEntry(uvec& valid, uvec& misid, double rnd);
	
	/* Do nothing */
	TARGET GWAS_MEM();

	/* Allocate memory and perform Eigen-value decomposition */
	TARGET void Calc(uvec& valid, HASH _hash, double rnd);

	/* Destructor */
	TARGET ~GWAS_MEM();
};

//#pragma pack(pop)

template<typename REAL> 
extern REAL GWAS_NAN;
template<typename REAL> 
extern GWAS<REAL> GWAS_Root;
template<typename REAL> 
extern GWAS<REAL>* GWAS_Threads;
template<typename REAL> 
extern atomic<REAL>* GWAS_Rmem;
template<typename REAL> 
extern umap<HASH, GWAS_MEM<REAL>> GWAS_umap;
template<typename REAL> 
extern map<double, GWAS_MEM<REAL>*> GWAS_map;
template<typename REAL> 
extern rmat* GWAS_tR;

extern timepoint GWAS_begin;
extern vector<string> GWAS_colname;
extern vector<string> GWAS_coltype;
extern umap<string, int> GWAS_indid;
extern atomic<int64> GWAS_batch_index;
extern atomic<int> GWAS_ploidy_presence[N_MAX_PLOIDY + 1];

/* Is a NA or NAN string */
TARGET bool IsNaNStr(char* ptr);

/* Read a line from csv file */
TARGET bool ReadCsvLine(vector<string>& row, std::ifstream& file);

/* Remove duplicated columns in a design matrix */
template<typename REAL>
TARGET void RemoveDupCol(rmat& mX);

/* mean vector discarding nan values */
template<typename REAL>
TARGET rrow NanMean(rmat& x);

/* var vector discarding nan values */
template<typename REAL>
TARGET rrow NanVar(rmat& x, bool ispop = false);

/* covariance matrix discarding nan values */
template<typename REAL>
TARGET rmat NanCov(rmat& x, bool ispop = false);

/* Generate multivariate normal distributed random vectors */
template<typename REAL>
TARGET rmat MVNRand(RNG<double>& rng, rmat& Mu, rmat& Sigma, int64 n = 1);

/* Fill missing data */
template<typename REAL>
TARGET void Imputation(int stage, int method, rmat& G, int64 lst, int64 led, int64* offset, bool shrinked, rmat& R);

/* Calculate GWAS */
template<typename REAL>
TARGET void CalcGWAS();

/* Calculate GWAS using multiple threads */
THREAD2H(GWASRelatednessThread);

/* Calculate GWAS using multiple threads */
THREAD2H(GWASMatrixRelatednessThread);

/* Calculate GWAS using multiple threads */
THREAD2H(GWASThread);