/* Genome-Wide Association Studies Functions */

#include "vcfpop.h"

template struct PCOA<double>;
template struct PCOA<float >;

template TARGET void CalcPCOA<double>();
template TARGET void CalcPCOA<float >();

template<> double* pcoa_matrix<double>;
template<> float * pcoa_matrix<float >;

#ifndef _PCOA
/* Do nothing */
template<typename REAL>
TARGET PCOA<REAL>::PCOA()
{

}

/* Destructor */
template<typename REAL>
TARGET PCOA<REAL>::~PCOA()
{
}

/* Perform PCoA */
template<typename REAL>
TARGET int PCOA<REAL>::CalcPCoA(int _maxp)
{
	// performing pcoa
	//http://www.esapubs.org/archive/ecol/E084/011/CAP_UserNotes.pdf
	maxp = _maxp = N - 1 > _maxp ? _maxp : N - 1;

	REAL* D1 = new REAL[N * N];
	REAL* C = new REAL[N * N];
	REAL* G = new REAL[N * N];

	SetVal(D1, D, N * N);

	// Maximum normal distance
	REAL ma = (REAL)-1e20, ex = 0;
	int npair = 0;
	for (int i = 0; i < N; ++i)
		for (int j = i; j < N; ++j)
		{
			REAL val = D1[i * N + j];
			if (IsError(val)) continue;

			if (ma < val && val < 1e20)
				ma = val;
			ex += val;
			npair++;
		}
	ex /= npair;

	// Check distances
	for (int i = 0; i < N; ++i)
		for (int j = i; j < N; ++j)
		{
			if (D1[i * N + j] < 0 || IsError(D1[i * N + j]))
				D1[i * N + j] = D1[j * N + i] = ex;
			if (D1[i * N + j] > 1e20)
				D1[i * N + j] = D1[j * N + i] = ma * 1.2;
		}

	// Total variance
	Vt = 0;
	for (int i = 0; i < N; ++i)
		for (int j = 0; j < i; ++j)
			Vt += D1[i * N + j] * D1[i * N + j];
	Vt /= N * (N - 1);

	// Centering
	for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j)
			D1[i * N + j] = -0.5 * D1[i * N + j] * D1[i * N + j];

	SetVal(C, (REAL)(-1.0 / N), N * N);
	for (int i = 0; i < N; ++i)
		C[i * N + i] = 1 - 1.0 / N;

	//MatrixMul2(C, D1, N, G);
	REAL* tU = new REAL[N * N];
	REAL* tV = new REAL[N];

	rmat mC(C, N, N, false, true);
	rmat mD(D1,N, N, false, true);
	rmat mG(G, N, N, false, true);
	rmat mU(tU, N, N, false, true);
	rcol cV(tV, N, false, true);
	mG = mC * mD * mC;

	if (N > 2)
	{
		int* idx = (int*)C;
		Evd(mG, mU, cV);
		int nz = sum(cV > 0), ncols = mU.n_cols;
		maxp = std::min(maxp, nz);

		U = new REAL[N * maxp];
		V = new REAL[maxp];
		for (int i = 0; i < maxp; ++i)
		{
			V[i] = cV(i, 0);
			SetVal(U + i * N, mU.memptr() + i * N, N);
			Mul(U + i * N, MySqrt(V[i]) * (U[i * N] >= 0 ? 1 : -1), N);
		}
	}
	else if (N == 2)
	{
		REAL b = abs(G[0]);
		maxp = 1;
		U = new REAL[N * maxp];
		V = new REAL[N];
		V[0] = b * 2;
		U[0 * maxp + 0] = b * MySqrt(2);
		U[1 * maxp + 0] = b * MySqrt(2);
	}
	else
	{
		maxp = 0;
		U = new REAL[1];
		V = new REAL[1];
		V[0] = U[0] = 0;
	}

	DEL(D1);
	DEL(C);
	DEL(G);
	DEL(tU);
	DEL(tV);
	return maxp;
}

/* Write PCoA */
template<typename REAL>
TARGET void PCOA<REAL>::WritePCoA(FILE* fout, REAL* d, int n, int _est, int _type)
{
	D = d;
	N = n;
	estimator = _est;
	type = _type;
	U = V = 0;

	CalcPCoA(pcoa_dim_val);

	fprintf(fout, "%s%s%s%sTotal variance%c", g_linebreak_val, g_linebreak_val, GD_ESTIMATOR[estimator], g_linebreak_val, g_delimiter_val);
	WriteReal(fout, Vt);

	fprintf(fout, "%sVariance", g_linebreak_val);
	for (int i = 0; i < maxp; ++i)
	{
		fprintf(fout, "%c", g_delimiter_val);
		WriteReal(fout, V[i] / (N - 1));
	}

	switch (type)
	{
	case 1: fprintf(fout, "%sInd%cPop", g_linebreak_val, g_delimiter_val); break;
	case 2: fprintf(fout, "%sPop%cRegL1", g_linebreak_val, g_delimiter_val); break;
	case 3:
	default:
		fprintf(fout, "%sRegL%d%cRegL%d", g_linebreak_val, type - 2, g_delimiter_val, type - 1); break;
	}

	for (int i = 0; i < maxp; ++i)
		fprintf(fout, "%cPC%d", g_delimiter_val, i + 1);

	if (type == 1)
	{
		int nn = nind;
		for (int i = 0; i < nn; ++i)
		{
			fprintf(fout, "%s%s%c%s", g_linebreak_val, ainds<REAL>[i]->name, g_delimiter_val, apops<REAL>[ainds<REAL>[i]->popid]->name);
			for (int j = 0; j < maxp; ++j)
			{
				fprintf(fout, "%c", g_delimiter_val);
				WriteReal(fout, U[i + j * nn]);
			}
		}
	}

	if (type == 2)
	{
		int nn = npop;
		for (int i = 0; i < nn; ++i)
		{
			fprintf(fout, "%s%s%c%s", g_linebreak_val, apops<REAL>[i]->name, g_delimiter_val,
				apops<REAL>[i]->rid == -1 ? "Total" : aregs<REAL>[0][apops<REAL>[i]->rid]->name);
			for (int j = 0; j < maxp; ++j)
			{
				fprintf(fout, "%c", g_delimiter_val);
				WriteReal(fout, U[i + j * nn]);
			}
		}
	}

	if (type >= 3)
	{
		int nn = nreg[type - 3];
		for (int i = 0; i < nn; ++i)
		{
			fprintf(fout, "%s%s%c%s", g_linebreak_val, aregs<REAL>[type - 3][i]->name, g_delimiter_val,
				aregs<REAL>[type - 3][i]->rid == -1 ? "Total" : aregs<REAL>[type - 2][apops<REAL>[i]->rid]->name);
			for (int j = 0; j < maxp; ++j)
			{
				fprintf(fout, "%c", g_delimiter_val);
				WriteReal(fout, U[i + j * nn]);
			}
		}
	}

	DEL(U);
	DEL(V);
}
#endif

#define extern 
template<typename REAL>
extern REAL* pcoa_matrix;					//Genetic distance array of genetic distance to perform PCoA
#undef extern 

/* Calculate principal coordinate analysis */
template<typename REAL>
TARGET void CalcPCOA()
{
	if (!pcoa) return;

	EvaluationBegin();
	GDIST_METHOD = 2;
	OpenResFile("-pcoa", "Principal coordinate analysis");

	bool isfirst = true;
	int ntot = 0;
	if (pcoa_level_val[1]) ntot += nind * (nind - 1) / 2;
	if (pcoa_level_val[2]) ntot += npop * (npop - 1) / 2;
	if (pcoa_level_val[3])
		for (int rl = 0; rl < lreg; ++rl)
			ntot += nreg[rl] * (nreg[rl] - 1) / 2;

	for (int m = 1; m <= 3; ++m)
	{
		if (pcoa_level_val[m] == 0) continue;
		for (int rl = 0; rl < (m < 3 ? 1 : lreg); ++rl)
		{
			PCOA<REAL> tpcoa;
			gdindex[0] = m + rl;
			int n = m == 1 ? nind : (m == 2 ? npop : nreg[rl]);

			int nestimator = 0;
			for (int k = 1; k <= (m == 1 ? N_GD_ESTIMATOR - 2 * N_FST_ESTIMATOR : N_GD_ESTIMATOR); ++k)
				if (pcoa_estimator_val[k])
					gdindex[k] = nestimator++;

			pcoa_matrix<REAL> = new REAL[n * n * nestimator];
			SetZero(pcoa_matrix<REAL>, n * n * nestimator);

			//Calculate genetic distance table between any two genotypes
			if (m == 1) GDIST<REAL>::CacheIndGD();

			RunThreads(&PCoAThread<REAL>, NULL, NULL, ntot, n * (n - 1) / 2,
				"\nPerforming principal coordinate analysis:\n", g_nthread_val, isfirst);
			isfirst = false;

			for (int k = 1; k <= (m == 1 ? N_GD_ESTIMATOR - 2 * N_FST_ESTIMATOR : N_GD_ESTIMATOR); ++k)
				if (pcoa_estimator_val[k])
					tpcoa.WritePCoA(FRES, pcoa_matrix<REAL> + gdindex[k] * n * n, n, k, m + rl);

			if (m == 1)
			{
				DEL(gd_tab<REAL>[0]);
				DEL(gd_tab<REAL>);
			}
			DEL(pcoa_matrix<REAL>);
		}
	}

	CloseResFile();
	GDIST_METHOD = 0;

	EvaluationEnd("PCoA");

	if (pcoa_plot_val == 1)
		RunRscript("pcoa_plot.R");
}

/* Calculate genetic distance for PCoA using multiple threads */
THREAD2(PCoAThread)
{
	//load ind
	int nthread = g_nthread_val;
	int type = gdindex[0];
	int n = type == 1 ? nind : (type == 2 ? npop : nreg[type - 3]);
	int64 progress = 0;
	byte* estimator = pcoa_estimator_val;
	GDIST<REAL> tbuf;

	VLA_NEW(p1, REAL, maxK * 2);
	VLA_NEW(p2, REAL, maxK * 2);
	for (int i = 0; i < n; ++i)
	{
		for (int j = i + 1; j < n; ++j)
		{
			if (progress++ % nthread != threadid) continue;
			switch (type)
			{
			case 1:  tbuf.CalcGD(ainds<REAL>[i], ainds<REAL>[j], p1, p2); break;
			case 2:  tbuf.CalcGD(apops<REAL>[i], apops<REAL>[j], (double*)p1); break;
			case 3:
			default: tbuf.CalcGD(aregs<REAL>[type - 3][i], aregs<REAL>[type - 3][j], (double*)p1);  break;
			}

			for (int k = 1; k <= (type == 1 ? N_GD_ESTIMATOR - 2 * N_FST_ESTIMATOR : N_GD_ESTIMATOR); ++k)
				if (estimator[k])
					*(pcoa_matrix<REAL> + n * n * gdindex[k] + j * n + i) =
					*(pcoa_matrix<REAL> + n * n * gdindex[k] + i * n + j) =
					*(&tbuf.Nei1972 + k - 1);

			PROGRESS_VALUE++;
		}
	}
	VLA_DELETE(p1);
	VLA_DELETE(p2);
}
