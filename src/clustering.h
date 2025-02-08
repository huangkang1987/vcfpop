/* Hierarchical Clustering and PCoA Functions */

#pragma once
#include "vcfpop.h"

template<typename REAL> struct HCLUSTER;
template<typename REAL> struct HCLUSTERING;

#pragma pack(push, 1)

template<typename REAL>
struct HCLUSTER
{
	bool isend;								//Is a leaf node
	char* endname;							//Object name
	double x;								//Node coordinate x
	double y;								//Node coordinate y
	uint* id;								//Objects ids
	uint  idlen;							//Objects size
	HCLUSTER<REAL>* left;					//Left node
	HCLUSTER<REAL>* right;					//Right node
};

template<typename REAL>
struct HCLUSTERING
{
	int method;								//Clustering method
	REAL* dori; 							//Original distance matrix
	REAL* dcur;								//New distance matrix
	REAL* dnew;								//Current distance matrix
	int nori;								//Dimension of dori
	int ncur;								//Dimension of dcur
	LIST<HCLUSTER<REAL>*> node;				//Nodes
	MEMORY* memory;							//Memory class

	/* Initialize for distance matrix between individuals */
	TARGET HCLUSTERING(REAL* d, IND<REAL>** obj, int n, int m, MEMORY* _memory);

	/* Initialize for distance matrix between populations or regions */
	TARGET HCLUSTERING(REAL* d, POP<REAL>** obj, int n, int m, MEMORY* _memory);

	/* Uninitialize */
	TARGET ~HCLUSTERING();

	/* Perform clustering */
	TARGET void Cluster(); 

	/* Find index for minimum distance in the distance matrix */
	TARGET REAL FindMinIdx(int& a, int& b);

	/* Reduce distance matrix from nxn to (n-1)x(n-1) */
	TARGET void ReduceMatrix(int _a, int _b);

	/* Write clustering results */
	TARGET void WriteClustering(FILE* fout, HCLUSTER<REAL>* c = NULL, REAL cy = 0);
};

#pragma pack(pop)

extern MEMORY* clustering_memory;						//Genetic distance memory class
template<typename REAL>
extern REAL* clustering_matrix;							//Genetic distance array of genetic distance to perform PCoA

/* Calculate hierarchical clustering */
template<typename REAL>
TARGET void CalcClustering();

/* Calculate genetic distance for Hierarchical clustering using multiple threads */
THREAD2H(ClusteringThread);