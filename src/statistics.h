/* Statistics Functions */

#pragma once
#include "vcfpop.h"

#pragma pack(push, 1)

/* Random number generator */
template<typename REAL>
struct RNG
{

};

template<>
struct RNG<double>
{
	uint64 x;
	uint64 y;

	double U1;
	double U2;

	/* Initialize rng */
	TARGET RNG();

	/* Initialize rng */
	TARGET RNG(uint64 s, uint64 salt);

	/* Draw a uniform distriubted interger */
	TARGET uint64 XorShift();

	/* Draw a polynormial distriubted integer */
	TARGET int Poly(double* a, int n, int sep = 1);

	/* Draw a polynormial distriubted integer with propoirtions in natural logarithm */
	TARGET int PolyLog(double* a, int n, int sep = 1);

	/* Get a random sequence from 0 ~ n-1 */
	template<typename T>
	TARGET void GetRandSeq(T* td, T n);

	/* Draw a uniform distriubted real number */
	TARGET double Uniform(double min = 0, double max = 1);

	/* Draw uniform distriubted real numbers */
	template<typename T1>
	TARGET void Uniform(T1* a, int n, double min = 0, double max = 1);

	/* Draw a normal distriubted real number */
	TARGET double Normal(double mean = 0, double std = 1);

	/* Draw normal distriubted real numbers */
	template<typename T1>
	TARGET void Normal(T1* a, int n, double mean = 0, double std = 1);

	/* Draw a real number from gamma distribution */
	TARGET double Gamma(double alpha, double beta = 1);

	/* Draw a real number from beta distribution */
	TARGET double Beta(double a, double b);

	/* Draw a vector from Dirichlet distribution D(a1 f, a2 f, ...) */
	template<typename T1, typename T2>
	TARGET void Dirichlet(T1* res, T2* a, int n, double f = 1);

	/* Draw a vector from Dirichlet distribution D(a1 + b1, a2 + b2, ...) */
	template<typename T1, typename T2, typename T3>
	TARGET void Dirichlet(T1* res, T2* a, T3* b, int n);

	/* Shuffle an array */
	template<typename T, typename T2>
	TARGET void Permute(T* val, T2 n);

	/* Draw uniform distriubted intergers */
	template<typename T1>
	TARGET void Integer(T1* arr, int n, uint64 min = 0, uint64 max = (uint64)-1);

	/* Draw a uniform distriubted interger */
	TARGET uint64 Next(uint64 min, uint64 max);

	/* Draw a uniform distriubted interger */
	TARGET uint64 Next(uint64 max);

	/* Draw a uniform distriubted interger and avoid sample av */
	TARGET uint64 NextAvoid(uint64 max, uint64 av);
};

template<>
struct RNG<float>
{
	uint x;
	uint y;
	uint z;

	/* Initialize rng */
	TARGET RNG();

	/* Initialize rng */
	TARGET RNG(uint64 s, uint64 salt);

	/* Draw a uniform distriubted interger */
	TARGET uint XorShift();

	/* Draw a polynormial distriubted integer */
	TARGET int Poly(float* a, int n, int sep = 1);

	/* Draw a uniform distriubted real number */
	TARGET float Uniform(float min = 0, float max = 1);

#ifdef _xxx

	float U1;
	float U2;
	/* Get a random sequence from 0 ~ n-1 */
	template<typename T>
	TARGET void GetRandSeq(T* td, T n);

	/* Draw uniform distriubted real numbers */
	TARGET void Uniform(float* arr, int n, float min, float max = 1);

	/* Draw a normal distriubted real number */
	TARGET float Normal(float mean = 0, float std = 1);

	/* Draw normal distriubted real numbers */
	TARGET void Normal(float* arr, int n, float mean = 0, float std = 1);

	/* Draw a polynormial distriubted integer with propoirtions in natural logarithm */
	TARGET int PolyLog(float* a, int n, int sep = 1);

	/* Draw a real number from gamma distribution */
	TARGET float Gamma(float alpha, float beta = 1);

	/* Draw a real number from beta distribution */
	TARGET float Beta(float a, float b);

	/* Draw a vector from Dirichlet distribution D(a1 f, a2 f, ...) */
	template<typename T1, typename T2>
	TARGET void Dirichlet(T1* res, T2* a, int n, double f);

	/* Draw a vector from Dirichlet distribution D(a1, a2, ...) */
	template<typename T1, typename T2>
	TARGET void Dirichlet(T1* res, T2* a, int n);

	/* Draw a vector from Dirichlet distribution D(a1 + b1, a2 + b2, ...) */
	template<typename T1, typename T2, typename T3>
	TARGET void Dirichlet(T1* res, T2* a, T3* b, int n);

	/* Shuffle an array */
	template<typename T, typename T2>
	TARGET void Permute(T* val, T2 n);

	/* Draw a uniform distriubted interger */
	TARGET uint Next(uint min, uint max);

	/* Draw a uniform distriubted interger */
	TARGET uint Next(uint max);

	/* Draw a uniform distriubted interger and avoid sample av */
	TARGET uint NextAvoid(uint max, uint av);
#endif
};

#pragma pack(pop)


/* Gamma function */
TARGET double Gamma(double x);

/* Natural logarithm of Gamma function */
TARGET double LogGamma(double x);

/* Regularized incomplete Gamma function */
TARGET double GammaIncRegularized(double a, double x);

/* Incomplete Gamma function */
TARGET double GammaInc(double a, double x);

/* Natural logarithm of regularized incomplete Gamma function */
TARGET double LogGammaIncRegularized(double a, double x);

/* Natural logarithm of incomplete Gamma function */
TARGET double LogGammaInc(double a, double x);

/* Beta function */
TARGET double Beta(double a, double b);

/* Natural logarithm of Beta function */
TARGET double LogBeta(double a, double b);

/* Regularized incomplete Beta function */
TARGET double BetaIncRegularized(double a, double b, double z);

/* Incomplete Beta function */
TARGET double BetaInc(double a, double b, double z);

/* Natural logarithm of regularized incomplete Beta function */
TARGET double LogBetaIncRegularized(double a, double b, double z);

/* Natural logarithm of incomplete Beta function */
TARGET double LogBetaInc(double a, double b, double z);

/* Two tailled probabiliy of normal distribution */
TARGET double MinusLogPNormal(double x);

/* Minus Log(P) for large Chi-square static */
TARGET double MinusLogPChi2(double x2, double df);

/* Minus Log(P) for large F-static */
TARGET double MinusLogPF(double f, double df1, double df2);

/* Minus Log(P) for large t-static */
TARGET double MinusLogPT(double t, double df);

/* Cumulative distribution function for Chi-square distribution */
TARGET double ChiSquareDistCDF(double x2, double df);

// Cumulative distribution function for standard normal distribution
TARGET double NormalDistCDF(double x);

// Cumulative distribution function for T-distribution
TARGET double TDistCDF(double t, double df);

// Cumulative distribution function for F-distribution
TARGET double FDistCDF(double f, double df1, double df2);

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
