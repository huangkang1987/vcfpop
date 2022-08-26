/* Analysis of Molecular Variances Functions */


#pragma once
#include "vcfpop.h"

#pragma pack(push, 1)

/* Used to permute vessels in AMOVA */
class VESSEL_ITERATOR
{
public:
	int relative_id[N_MAX_REG + 3];			//Index in its nesting vessles (from 0~nsubunits)
	int universal_id[N_MAX_REG + 3];		//Index among all vessles at this lay (from 0~nhaplo)
	VESSEL* trace[N_MAX_REG + 3];			//Current and its higher level vessels
	int lay;								//Current lay

	/* Go to start */
	TARGET void Rewind(int nlay);

	/* Copy from a reference*/
	TARGET void Copy(VESSEL_ITERATOR& ref, int nlay);

	/* Initialize */
	TARGET VESSEL_ITERATOR();

	/* Initialize */
	TARGET VESSEL_ITERATOR(int _lay, VESSEL& root, int nlay);

	/* Uninitialize */
	TARGET ~VESSEL_ITERATOR();

	/* Go to next vessel */
	TARGET void Next(int nlay);

	/* Get haplotype index to calculate genetic distance */
	TARGET int GetHapId();

	/* Get allele to fetch calculate distance */
	TARGET ushort GetAllele();

	/* Get subpopulation in print SS */
	TARGET POP* GetSubpop(int nlay, int tlay);

	/* Get individual in print SS */
	TARGET IND* GetInd(int nlay, int tlay);
};

/* Vessel of genes in AMOVA */
class VESSEL
{
public:
	VESSEL** subunits;						//Vessels nested within this vessel
	int* nhaplos;							//Number of haplotypes at locus l
	int* allelecount;						//KT elements, + allele_freq_offset[l] is the count of alleles at locus l in this vessel
	int nsubunits;							//Number of subunits
	int nhaplo;								//Number of haplotypes, homoploid model
	int hid;								//Haplotype id, -1 for non-allele vessel
	short lay;								//Level
	ushort allele;							//Allele for aneu model

	/* Uninitialize */
	TARGET ~VESSEL();

	/* Initialize */
	TARGET VESSEL();

	/* Deep copy a vessel */
	TARGET VESSEL(VESSEL& r);

	/* Create vessel from population */
	TARGET VESSEL(POP* s, int _lay, int& _hid, int64 loc, int method);

	/* Create vessel from individual */
	TARGET VESSEL(IND* s, int _lay, int& _hid, int64 loc, int method);

	/* Create vessel from haplotype */
	TARGET VESSEL(int _lay, int& _hid, ushort _allele);

	/* Get the allele count array at locus l*/
	TARGET int* GetAlleleCount(int l);

	/* Save all vellels in level fa into an array */
	TARGET void GetVessels(VESSEL** vs, int& nvessels, int fa);

	/* Replace with shuffled vessels */
	TARGET int Replace(VESSEL** vs, int& nvessels, int fa, int method);

	/* Shuffle fa level vessels among fb level vessels */
	TARGET void Shuffle(RNG& rng, int fa, int fb, int method, VESSEL** buf);

	/* Calculate matrix C for maximum-likelihood method */
	TARGET void GetCML(double* C, int64 l, int* tid, double* tw, int Nh, int nlay, double** W);

	/* Initialize matrix C, the coefficient matrix of S = CV */
	TARGET void InitC(double* C, int* tid, int Nh, int nlay);

	/* Calculate matrix C */
	TARGET void GetC(double* tw, int* tid, double* C, int nlay, double** W, int64 l);

	/* Count number of vessels in each hierarchy */
	TARGET void CountVessels(int* count);

	/* Initialize W, W[lay][tid[lay]] = 1 / nhaplo */
	TARGET void InitW(MEMORY& mem, double**& W, int nlay);

	/* Calculate SS for homoploid method */
	TARGET void GetSSHomo(double* SS, double* gd, int Nh, double** W, int nlay, VESSEL_ITERATOR& ve1, VESSEL_ITERATOR& ve2);

	/* Calculate SS for aneuploid method */
	TARGET void GetSSAneu(ushort* hap_bucket, double* SS, bool isiam, int nh, int64 l, int k, double* missing0, double** W, int nlay, VESSEL_ITERATOR& ve1, VESSEL_ITERATOR& ve2);

	/* Calculate SS for aneuploid method */
	TARGET void GetSSAneu(double* SS, bool isiam, int nh, int k, double* missing0, double** W, int nlay, VESSEL_ITERATOR& ve1, VESSEL_ITERATOR& ve2);

	/* Calculate variance component matrix V */
	TARGET static void GetV(double* C, double* SS, double*& V, int nlay);

	/* Calculate F-statistics */
	TARGET static void GetF(double* V, double* F, double* vs, int nlay);

	/* Calculate F-statistics */
	TARGET static void GetF(double* Fi, double* F, int nlay);
};

/* Use to avoid thread conflict when calculate distance matrix */
struct AMOVA_DISTBUF
{
	double* pointer;
	double val;
};

struct AMOVA
{
	//One time
	double* V;								//Variance components at each level, nlay elements
	double* F;								//F-statistics, nlay*nlay elements
	int* G;									//Number of permuations with F'>F, nlay*nlay elements
	int* E;									//Number of permuations with F'=F, nlay*nlay elements

	double* EF2;							//Mean permuated squared F', nlay*nlay elements
	double* EF;								//Mean permuated F', nlay*nlay elements

	double* SS;								//SS within each level, nlay elements
	double* DF;								//Degrees-of-freedom for each level, nlay elements
	double** SSW;							//SS within each vessel
	int* nSSW;								//Size of each SSW, nlay elements

	int nlay;								//Number of levels
	int method;								//Method: 1 homoploid, 2 aneuploid or 4 likelihood
	int npermed;							//Number of permuated iterations
	int Lind;								//Consider individual level?

	/* Initialize */
	TARGET AMOVA();

	/* Extract dummy haplotype for homoploid method */
	TARGET void GetHaplotype(ushort* bucket);

	/* Perform AMOVA using homoploid method */
	TARGET void CalcAMOVA_homo();

	/* Perform AMOVA using aneuploid method */
	TARGET void CalcAMOVA_aneu();

	/* Calculate likelihood for permuated data */
	TARGET static double Likelihood(CPOINT& xx, void** Param);

	/* Perform AMOVA using maximum-likelihood method */
	TARGET void CalcAMOVA_ml();

	/* Destructor */
	TARGET ~AMOVA();

	/* Write results */
	TARGET void PrintAMOVA(FILE* fout);
};

#pragma pack(pop)

extern _thread MEMORY* amova_memory;				//memory class for amova vessels
extern AMOVA* amova_buf;							//Read/Write buffer, nthread
extern double* amova_matrix;						//Genetic distance matrix used in AMOVA

/* Calculate analysis of molecular variance */
TARGET void CalcAMOVA();

/* Calculate AMOVA using multiple threads */
THREADH(AMOVAThread);