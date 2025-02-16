/* Sliding Window Functions */

#pragma once
#include "vcfpop.h"

template<typename REAL> struct WINDOW;

#pragma pack(push, 1)

struct CHROM_PROP
{
	uint64 min;
	uint64 max; 
	uint64 st;
	uint64 ed;
};

template<typename REAL>
struct SWINDOW
{
	uint window_count;
	uint swindow_count;
	uint window_size;
	uint window_step;

	double* slide_a1;
	double* slide_a2;
	POP<REAL>* cpop_;
	POP<REAL>** grps;
	int ngrps;

	bool cNei1973;
	bool cWeir1984;
	bool cHudson1992;
	bool cHedrick2005;
	bool cJost2008;
	bool cHuang2021_aneu;
	bool cdxy;
	bool cpi;
	bool cthetaw;
	bool ctajimad;
	bool cr2;
	bool cdprime;
	bool cr2delta;
	bool cdeltaprime;
	bool cfis;
	bool cho;
	bool che;
	bool cpic;
	bool cae;
	bool cI;
	bool clocipair;

	int pst, ped;
	int _K;					//number of polymorphic variants
	int _C;					//number of variants

	/* Fst */
	double _Fst_Nei1973_sum1;
	double _Fst_Nei1973_sum2;
	double _Fst_Weir1984_sum1;
	double _Fst_Weir1984_sum2;
	double _Fst_Hudson1992_sum1;
	double _Fst_Hudson1992_sum2;
	double _Fst_Hedrick2005_sum1;
	double _Fst_Hedrick2005_sum2;
	double _Fst_Jost2008_sum1;
	double _Fst_Jost2008_sum2;
	double _Fst_Huang2021_aneu_sum1;
	double _Fst_Huang2021_aneu_sum2;

	/* Absolute divergence */
	double _dxy_sum;

	/* 1 Nucleotide diversity */
	double _pi_sum;

	/* 2 Watterson's thetaw */
	double _thetaw_sum;

	/* 3 Tajima's D */
	double _d_sum;
	TABLE<HASH, int> _nhaplo_freq1;
	umap  <int,  int> _nhaplo_freq2;

	/* 4 r2 */
	double _D2_sum;
	double _Q_sum;

	/* 5 D' */
	double _D_abs_sum;
	double _Dmax_abs_sum;

	/* 6 r2Delta */
	double _Delta2_sum;
	double _R_sum;

	/* 7 Delta' */
	double _Delta_abs_sum;
	double _Deltamax_abs_sum;

	/* 8 Fis */
	double _Fis_Nei1973_sum1;
	double _Fis_Nei1973_sum2;

	/* 9 Diversity */
	double _ho_sum;
	double _he_sum;
	double _pic_sum;
	double _ae_sum;
	double _I_sum;

	int* window_id; 
	int* K;					//number of polymorphic variants
	int* C;					//number of variants
	LOCK* lock;

	/* Fst */
	double* Fst_Nei1973_sum1;
	double* Fst_Nei1973_sum2;
	double* Fst_Weir1984_sum1;
	double* Fst_Weir1984_sum2;
	double* Fst_Hudson1992_sum1;
	double* Fst_Hudson1992_sum2;
	double* Fst_Hedrick2005_sum1;
	double* Fst_Hedrick2005_sum2;
	double* Fst_Jost2008_sum1;
	double* Fst_Jost2008_sum2;
	double* Fst_Huang2021_aneu_sum1;
	double* Fst_Huang2021_aneu_sum2;

	/* Absolute divergence */
	double* dxy_sum;

	/* 1 Nucleotide diversity */
	double* pi_sum;

	/* 2 Watterson's thetaw */
	double* thetaw_sum;

	/* 3 Tajima's D */
	double* d_sum;
	//double* vd_sum;
	TABLE<HASH, int>* nhaplo_freq1;
	umap  <int,  int>* nhaplo_freq2;

	/* 4 r2 */
	double* D2_sum;
	double* Q_sum;

	/* 5 D' */
	double* D_abs_sum;
	double* Dmax_abs_sum;

	/* 6 r2Delta */
	double* Delta2_sum;
	double* R_sum;

	/* 7 Delta' */
	double* Delta_abs_sum;
	double* Deltamax_abs_sum;

	/* 8 Fis */
	double* Fis_Nei1973_sum1;
	double* Fis_Nei1973_sum2;

	/* 9 Diversity */
	double* ho_sum;
	double* he_sum;
	double* pic_sum;
	double* ae_sum;
	double* I_sum;

	/* Calculate diveristy indices */
	TARGET void CalcLocus(int64 l, double* buf);

	/* Calculate Tajima's D */
	//TARGET void CalcLociPair(int64 l1, int nhaplo1, double* buf);

	/* Initialize a SWINDOW */
	TARGET void InitSWindow();

	/* Uninitialize a SWINDOW */
	TARGET void UnInitSWindow();

	/* When the new locus has a different range, distributed current results to SWINDOW */
	TARGET void Settle1();

	/* When the new loci pair has a different range, distributed current results to SWINDOW */
	//TARGET void Settle2(int pst2, int st2);

	/* Settle Tajima D denominator */
	TARGET double SettleTajimaD(TABLE<HASH, int>& freq1, umap<int, int>& freq2);

	/* Settle a SWINDOW */
	TARGET bool Settle(int i, int newid);

	/* Get sliding window index range*/
	TARGET void GetWindowId(int64 l, int& st, int& ed);
};

template<typename REAL>
struct WINDOW
{
	uint window_count;
	uint swindow_count;
	uint window_size;
	uint window_step;

	double* slide_a1;
	double* slide_a2;
	POP<REAL>* cpop_;
	POP<REAL>** grps;
	int ngrps;

	bool cNei1973;
	bool cWeir1984;
	bool cHudson1992;
	bool cHedrick2005;
	bool cJost2008;
	bool cHuang2021_aneu;
	bool cdxy;
	bool cpi;
	bool cthetaw;
	bool ctajimad;
	bool cr2;
	bool cdprime;
	bool cr2delta;
	bool cdeltaprime;
	bool cfis;
	bool cho;
	bool che;
	bool cpic;
	bool cae;
	bool cI;
	bool clocipair;

	char** chrom;
	
	int* K;					//number of polymorphic variants
	int* C;					//number of variants
	LOCK* lock;				

	/* Fst */
	REAL* Fst_Nei1973_sum1;
	REAL* Fst_Nei1973_sum2;
	REAL* Fst_Weir1984_sum1;
	REAL* Fst_Weir1984_sum2;
	REAL* Fst_Hudson1992_sum1;
	REAL* Fst_Hudson1992_sum2;
	REAL* Fst_Hedrick2005_sum1;
	REAL* Fst_Hedrick2005_sum2;
	REAL* Fst_Jost2008_sum1;
	REAL* Fst_Jost2008_sum2;
	REAL* Fst_Huang2021_aneu_sum1;
	REAL* Fst_Huang2021_aneu_sum2;

	/* Absolute divergence */
	REAL* dxy_sum;

	/* 1 Nucleotide diversity */
	REAL* pi_sum;

	/* 2 Watterson's thetaw */
	REAL* thetaw_sum;

	/* 3 Tajima's D */
	REAL* d_sum;
	REAL* vd_sum;
	int* N;
	TABLE<HASH, int>* nhaplo_freq1;
	umap  <int,  int>* nhaplo_freq2;

	/* 4 r2 */
	REAL* D2_sum;
	REAL* Q_sum;

	/* 5 D' */
	REAL* D_abs_sum;
	REAL* Dmax_abs_sum;

	/* 6 r2Delta */
	REAL* Delta2_sum;
	REAL* R_sum;

	/* 7 Delta' */
	REAL* Delta_abs_sum;
	REAL* Deltamax_abs_sum;

	/* 8 Fis */
	REAL* Fis_Nei1973_sum1;
	REAL* Fis_Nei1973_sum2;

	/* 9 Diversity */
	REAL* ho_sum;
	REAL* he_sum;
	REAL* pic_sum;
	REAL* ae_sum;
	REAL* I_sum;

	/* Initialize all windows */
	TARGET void InitWindow();

	/* UnInitialize all windows */
	TARGET void UnInitWindow();

	/* Write result file */
	TARGET void Write();

	/* Get sliding window index range*/
	TARGET void GetWindowId(int64 l, int& st, int& ed);
};

#pragma pack(pop)

template<typename REAL> extern WINDOW<REAL> window;
extern vector<char*> slide_chromas;
extern umap<HASH, CHROM_PROP> slide_chrom_sted;

/* Convert phased genotype into unphased genotype */
THREAD2H(UpdateUnphaseGenotypes);

/* Calculate sliding window */
template<typename REAL>
TARGET void CalcSlide();

/* Calculate sliding window using multiple threads */
THREAD2H(SlidingWindowThread);

/* Prepare allele frequency and individual ploidy */
THREAD2H(SlidePrepare);

/* Calculate allele frequencies for each population and region */
THREAD2H(SlideFreqThread);