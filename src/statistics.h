/* Statistics Functions */

#pragma once
#include "vcfpop.h"
#pragma pack(push, 1)

/* Gamma function */
TARGET double Gamma1(double x);

/* Natural logarithm of Gamma function */
TARGET double LogGamma1(double x);

/* Incomplete Gamma function */
TARGET double Gamma2(double a, double x);

/* Natural logarithm of incomplete Gamma function */
TARGET double LogGamma2(double a, double x);

/* Right tail probability of Chi-square distribution */
TARGET double ChiSquareProb(double x2, double df);

/* Random number generator */

struct RNG
{
public:
	uint64 x;
	uint64 y;
	uint64 seed;

	double U1, U2;
	bool state;

	/* Initialize rng */
	TARGET RNG();

	/* Initialize rng */
	TARGET RNG(uint64 s);

	/* Draw a uniform distriubted interger */
	TARGET uint64 XorShift128p();

	/* Get a random sequence from 0 ~ n-1 */
	TARGET void GetRandSeq(int* td, int n);

	/* Shuffle an array */
	template<typename T>
	TARGET void Permute(T* val, int n)
	{
		//https://lemire.me/blog/2016/06/30/fast-random-shuffling/
		for (int i = n; i > 1; --i)
			Swap(val[i - 1], val[Next(i)]);
	}

	/* Draw a uniform distriubted real number */
	TARGET double Uniform();

	/* Draw a uniform distriubted real number */
	TARGET double Uniform(double min, double max);

	/* Draw a uniform distriubted real number */
	TARGET double Uniform(double max);

	/* Draw a uniform distriubted interger */
	TARGET int64 Next(int64 min, int64 max);

	/* Draw a uniform distriubted interger */
	TARGET int64 Next(int64 max);

	/* Draw a uniform distriubted interger and avoid sample av */
	TARGET int64 NextAvoid(int64 max, int64 av);

	/* Draw a normally distriubted real number */
	TARGET double Normal();

	/* Draw a normally distriubted real number */
	TARGET double Normal(double mean, double std);

	/* Draw a polynormial distriubted integer */
	TARGET int Poly(double* a, int n);

	/* Draw a polynormial distriubted integer with propoirtions in natural logarithm */
	TARGET int PolyLog(double* a, int n);
	
	/* Draw a polynormial distriubted integer with propoirtions in natural logarithm */
	TARGET int PolyLog(double* a, int n, int sep);

	/* Draw a real number from gamma distribution */
	TARGET double Gamma(double alpha, double beta = 1);

	/* Draw a real number from beta distribution */
	TARGET double Beta(double a, double b);

	/* Draw a vector from Dirichlet distribution D(a1 f, a2 f, ...) */
	TARGET void Dirichlet(double* res, double* a, int n, double f);

	/* Draw a vector from Dirichlet distribution D(a1, a2, ...) */
	TARGET void Dirichlet(double* res, double* a, int n);

	/* Draw a vector from Dirichlet distribution D(a1 + b1, a2 + b2, ...) */
	TARGET void Dirichlet(double* res, double* a, int64* b, int n);

	/* Draw a vector from Dirichlet distribution D(a1 + b1, a2 + b2, ...) */
	TARGET void Dirichlet(double* res, double* a, int* b, int n);
};

/* Chi-square test, combine tables with any grid with an expectition < 5 */

/* Input a vector, return proportion of grids with exp < threshold and two mininum indices */
TARGET double FindMinIndex(double* exp, int m, int& i1, int& i2, double threshold);

/* Input a column, return proportion of grids with exp < threshold and two mininum indices */
TARGET double FindMinIndex(double* exp, int m, int n, int& i1, int& j1, int& i2, int& j2, double threshold);

/* strcmp two rows */
TARGET int CompareRow(double* obs, int n, int r1, int r2);

/* strcmp two columns */
TARGET int CompareCol(double* obs, int m, int n, int c1, int c2);

/* Sort columns */
TARGET void SortCol(double* obs, int m, int n);

/* Chi-square test, combine tables with any grid with an expectition < 5 */
TARGET void CombineTable(double* obs, int m, int n, double& g, int& df, double& p, bool test, double* obs2, double* exp, double* exp2, double* rowsum, double* colsum);

/* Factorial */
TARGET double Factorial(int n);

/* Binomial coefficient */
TARGET double Binomial(int n, int r);

/* Natural logarithm of binomial coefficient */
TARGET double LogBinomial(int n, int r);

/* Initialize BINOMIAL global variable */
TARGET void InitBinomial();

#pragma pack(pop)