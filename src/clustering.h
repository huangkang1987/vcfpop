/* Hierarchical Clustering and PCoA Functions */

#pragma once
#include "vcfpop.h"

#pragma pack(push, 1)

struct HCLUSTER
{
	bool isend;								//Is a leaf node
	char* endname;							//Object name
	double x;								//Node coordinate x
	double y;								//Node coordinate y
	uint* id;								//Objects ids
	uint  idlen;							//Objects size
	HCLUSTER* left;							//Left node
	HCLUSTER* right;						//Right node
};

class HCLUSTERING
{
public:
	int method;								//Clustering method
	double* dori; 							//Original distance matrix
	double* dcur;							//New distance matrix
	double* dnew;							//Current distance matrix
	int nori;								//Dimension of dori
	int ncur;								//Dimension of dcur
	LIST<HCLUSTER*> node;					//Nodes
	MEMORY* memory;							//Memory class

	/* Initialize for distance matrix between individuals */
	TARGET HCLUSTERING(double* d, IND** obj, int n, int m, MEMORY* _memory);

	/* Initialize for distance matrix between populations or regions */
	TARGET HCLUSTERING(double* d, POP** obj, int n, int m, MEMORY* _memory);

	/* Uninitialize */
	TARGET ~HCLUSTERING();

	/* Perform clustering */
	TARGET void Cluster(); 

	/* Find index for minimum distance in the distance matrix */
	TARGET double FindMinIdx(int& a, int& b);

	/* Reduce distance matrix from nxn to (n-1)x(n-1) */
	TARGET void ReduceMatrix(int _a, int _b);

	/* Print clustering results */
	TARGET void PrintClustering(FILE* fout, HCLUSTER* c = NULL, double cy = 0);
};

#pragma pack(pop)

extern MEMORY* clustering_memory;					//Genetic distance memory class
extern double* clustering_matrix;					//Genetic distance array of genetic distance to perform PCoA

/* Calculate hierarchical clustering */
TARGET void CalcClustering();

/* Calculate genetic distance for Hierarchical clustering using multiple threads */
THREADH(ClusteringThread);