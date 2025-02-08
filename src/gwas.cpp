/* Genome-Wide Association Studies Functions */

#include "vcfpop.h"

using ifstream = std::ifstream;
using stringstream = std::stringstream;

template struct GWAS<double>;
template struct GWAS<float >;

template<> double GWAS_NAN<double> = (double)-9.87654321012345e-300;
template<> float  GWAS_NAN<float > = (float )-9.8765432e-30;

template<> GWAS<double> GWAS_Root<double>;
template<> GWAS<float > GWAS_Root<float >;

template<> GWAS<double>* GWAS_Threads<double>;
template<> GWAS<float >* GWAS_Threads<float >;

template<> umap<HASH, GWAS_MEM<double>> GWAS_umap<double>;
template<> umap<HASH, GWAS_MEM<float >> GWAS_umap<float >;

template<> map<double, GWAS_MEM<double>*> GWAS_map<double>;
template<> map<double, GWAS_MEM<float >*> GWAS_map<float >;

template TARGET void CalcGWAS<double>();
template TARGET void CalcGWAS<float >();

constexpr bool GWAS_DEBUG = false;

TARGET CPOINT BestPoint(CPOINT& a, CPOINT&& b)
{
	return a.lnL > b.lnL ? a : b;
}

#ifndef _GWAS_RESULTS

template<typename REAL>
TARGET void RESULT_WALD<REAL>::Write(FILE* fout)
{
	fprintf(fout, "%c", g_delimiter_val);	WriteReal(fout, lnL_REML);
	fprintf(fout, "%c", g_delimiter_val);	WriteReal(fout, W);
	//fprintf(fout, "%c", g_delimiter_val);	WriteReal(fout, df);
	fprintf(fout, "%c", g_delimiter_val);	WriteReal(fout, mlogP);
	fprintf(fout, "%c", g_delimiter_val);	WriteReal(fout, P);
}

template<typename REAL>
TARGET void RESULT_LRT<REAL>::Write(FILE* fout)
{
	fprintf(fout, "%c", g_delimiter_val);	WriteReal(fout, lnL_ML);
	fprintf(fout, "%c", g_delimiter_val);	WriteReal(fout, LR);
	//fprintf(fout, "%c", g_delimiter_val);	WriteReal(fout, df);
	fprintf(fout, "%c", g_delimiter_val);	WriteReal(fout, mlogP);
	fprintf(fout, "%c", g_delimiter_val);	WriteReal(fout, P);
}

template<typename REAL>
TARGET void RESULT_SCORE<REAL>::Write(FILE* fout)
{
	fprintf(fout, "%c", g_delimiter_val);	WriteReal(fout, lnL0_ML);
	fprintf(fout, "%c", g_delimiter_val);	WriteReal(fout, LM);
	//fprintf(fout, "%c", g_delimiter_val);	WriteReal(fout, df);
	fprintf(fout, "%c", g_delimiter_val);	WriteReal(fout, mlogP);
	fprintf(fout, "%c", g_delimiter_val);	WriteReal(fout, P);
}

#endif

#ifndef _GWAS

/* Write result header */
template<typename REAL>
TARGET void GWAS<REAL>::WriteGWASHeader()
{
	fprintf(fout, "Chrom%c", g_delimiter_val);
	fprintf(fout, "Position%c", g_delimiter_val);
	fprintf(fout, "Locus%c", g_delimiter_val);
	fprintf(fout, "#Alleles%c", g_delimiter_val);
	fprintf(fout, "#Ploidy");
	for (int i = 0; i < GWAS_coltype.size(); ++i)
	{
		if (GWAS_coltype[i] == "Normal")
		{
			fprintf(fout, "%c%s:n", g_delimiter_val, GWAS_colname[i].c_str());
			fprintf(fout, "%cdf", g_delimiter_val);

			if (gwas_test_val[1])
			{
				//Wald
				fprintf(fout, "%cWald:lnL_REML", g_delimiter_val);
				fprintf(fout, "%cW", g_delimiter_val);
				fprintf(fout, "%c-lnP", g_delimiter_val);
				fprintf(fout, "%cP", g_delimiter_val);
			}
			if (gwas_test_val[2])
			{
				//LRT
				fprintf(fout, "%cLRT:lnL_ML", g_delimiter_val);
				fprintf(fout, "%cLR", g_delimiter_val);
				fprintf(fout, "%c-lnP", g_delimiter_val);
				fprintf(fout, "%cP", g_delimiter_val);
			}
			if (gwas_test_val[3])
			{
				//Score
				fprintf(fout, "%cScore:lnL0_ML", g_delimiter_val);
				fprintf(fout, "%cLM", g_delimiter_val);
				fprintf(fout, "%c-lnP", g_delimiter_val);
				fprintf(fout, "%cP", g_delimiter_val);
			}
		}
	}
}

/* Write row header */
template<typename REAL>
TARGET void GWAS<REAL>::WriteGWASRowHeader()
{
	if (GWAS_DEBUG)
	{
		fprintf(fout, "%s", g_linebreak_val);
		fprintf(fout, "Chr1%c", g_delimiter_val);
		fprintf(fout, "%lld%c", cl + 1, g_delimiter_val);
		fprintf(fout, "%s%c", "Loc", g_delimiter_val);
		fprintf(fout, "%d%c", nallele, g_delimiter_val);
		fprintf(fout, "%d", nploidy);
	}
	else
	{
		auto loc = GetLoc(cl);
		fprintf(fout, "%s", g_linebreak_val);
		fprintf(fout, "%s%c", loc.GetChrom(), g_delimiter_val);
		fprintf(fout, "%lld%c", GetLocPos(cl), g_delimiter_val);
		fprintf(fout, "%s%c", loc.GetName(), g_delimiter_val);
		fprintf(fout, "%d%c", nallele, g_delimiter_val);
		fprintf(fout, "%d", nploidy);
	}
}

/* Write a cell */
template<typename REAL>
TARGET void GWAS<REAL>::WriteGWASCell(int p, int pid, int id)
{
	if (GWAS_coltype[p] == "Normal")
	{
		int resid = GWAS_Root<REAL>.nresponse * id + pid;

		fprintf(fout, "%c%lld%c%lld", g_delimiter_val, sample_size[resid], g_delimiter_val, degrees_of_freedom[resid]);
		if (gwas_test_val[1]) wald[resid].Write(fout);
		if (gwas_test_val[2]) lrt[resid].Write(fout);
		if (gwas_test_val[3]) score[resid].Write(fout);
	}
}

/* Write results for a locus */
template<typename REAL>
TARGET void GWAS<REAL>::WriteLocus(int64 _l)
{
	cl = _l;
	int id = cl % gwas_batch_val;

	nallele = noffset[cl + 1] - noffset[cl];
	n = GWAS_Root<REAL>.n;

	WriteGWASRowHeader();
	for (int yi = 1, yp = 0; yi < GWAS_coltype.size(); ++yi)
		if (GWAS_coltype[yi] == "Normal")
			WriteGWASCell(yi, yp++, id);
}

/* Do nothing */
template<typename REAL>
TARGET GWAS<REAL>::GWAS()
{

}

/* Destructor */
template<typename REAL>
TARGET GWAS<REAL>::~GWAS()
{
	if (tid == -1)
	{
		DEL(Gb);
		DEL(noffset);

		DEL(sample_size);
		DEL(degrees_of_freedom);
		DEL(wald);
		DEL(lrt);
		DEL(score);
	}
	umap<int, GWAS_NULL_MODEL<REAL>>().swap(lnL0_Res);
}

/* Clone */
template<typename REAL>
TARGET GWAS<REAL>::GWAS(GWAS<REAL>* ref, int _tid)
{
	tid = _tid;
	n = ref->n;
	m = ref->m;
	cl = -1;

	kX = ref->kX;
	kXI = ref->kXI;
	kXIG = ref->kXIG;

	nnumeric = ref->nnumeric;
	nnominal = ref->nnominal;
	nresponse = ref->nresponse;

	nploidy = ref->nploidy;
	nallele = ref->nallele;
	noffset = ref->noffset;

	SetVal(ploidy, ref->ploidy, N_MAX_PLOIDY);
	SetVal(ploidyidx, ref->ploidyidx, N_MAX_PLOIDY);

	Gb = ref->Gb;
	if (ref->oG.n_elem) REF_MAT(oG);

	REF_MAT(oY);
	REF_MAT(oX);
	REF_MAT(oI);
	REF_MAT(oR);
	REF_MAT(dR);

	sample_size = ref->sample_size;
	degrees_of_freedom = ref->degrees_of_freedom;
	wald  = ref->wald;
	lrt   = ref->lrt;
	score = ref->score;
}

/* Perform Score Test */
template<typename REAL>
TARGET void GWAS<REAL>::ScoreTest(CPOINT& xx_ml0, int resid, GWAS_NULL_MODEL<REAL>* null_model)
{
	rmat ma0; 
	double r = 0, sig = 0;
	if (null_model)
	{
		r = null_model->r;
		sig = null_model->sig;
		ma0 = null_model->ma0;
	}
	else
	{
		lnL_ml_profile1(fabs(xx_ml0.real_space[0]), NULL, NULL, false);
		r = fabs(xx_ml0.real_space[0]);
		sig = sqrt(sig2_from_r);
		ma0 = iEF01; 
	}

	//Get Grad vector
	ma0.reshape(kXIG, 1);
	SetIndVar(mXIG);

	// (r, sig) scheme
	rmat G, H;
	lnL2(r, sig, ma0, &G, &H);

	//The distribution of P under H0 is uniform so I use Hessian estimator for V
	rmat V = GetV_Hessian2(r, sig, mXIG);
	score[resid].lnL0_ML = xx_ml0.lnL;
	score[resid].LM = trace(G.t() * V * G);
	//score[resid].df = kXIG - kXI;
	score[resid].mlogP = MinusLogPChi2(score[resid].LM, degrees_of_freedom[resid]);
	score[resid].P = exp(-score[resid].mlogP);
}

/* Perform LRT Test */
template<typename REAL>
TARGET void GWAS<REAL>::LRTTest(CPOINT& xx_ml, CPOINT& xx_ml0, int resid)
{
	lrt[resid].lnL_ML = xx_ml.lnL;
	lrt[resid].LR = 2 * std::max(xx_ml.lnL - xx_ml0.lnL, 0.0);
	//lrt[resid].df = kXIG - kXI;
	lrt[resid].mlogP = MinusLogPChi2(lrt[resid].LR, degrees_of_freedom[resid]);
	lrt[resid].P = exp(-lrt[resid].mlogP);
}

/* Perform ANOVA
template<typename REAL>
TARGET void GWAS<REAL>::ANOVA(CPOINT& xx_reml, int resid)
{
	rmat G, H;
	lnL_reml_profile2(fabs(xx_reml.real_space[0]), NAN, &G, &H, false);
	double r = fabs(xx_reml.real_space[0]), sig2 = sig2_from_r, sig = sqrt(sig2);

	// t test
	if (0)
	{
		rmat var_con = diagvec(iE01) * sig2;
		rmat& slope = iEF01;
		rmat t = slope / sqrt(var_con);
		// Satterthwaite correction for d.f.
		rmat grad = join_horiz(diagvec(iE01 * mE12 * iE01) * (sig2 * 2 * r),
						            diagvec(iE01) * (2 * sig)).t();
		rmat denom = diagvec(grad.t() * solve(-H, grad, solve_opts::force_sym));
		rmat df = 2 * (var_con % var_con / denom);
		rmat mlogPt = t * 0;
		for (int i = 0; i < t.n_elem; ++i)
			mlogPt(i) = MinusLogPT(t(i), df(i));
		rmat Pt = exp(-mlogPt);
	}

	// Satterthwaite correction for d.f.
	int np = kXIG - kXI;
	rmat Jac1 = iE01 * mE12 * iE01; 
	Jac1 = Jac1.submat(kXI, kXI, kXIG - 1, kXIG - 1) * (sig2 * 2 * r);
	rmat Jac2 = iE01.submat(kXI, kXI, kXIG - 1, kXIG - 1) * (2 * sig);
	rmat mP; rcol cD;
	eig_sym(cD, mP, iE01.submat(kXI, kXI, kXIG - 1, kXIG - 1) * sig2, "dc");
	rmat JJ =  join_horiz(diagvec(mP.t() * Jac1 * mP), diagvec(mP.t() * Jac2 * mP));
	rcol nu = 2 * (cD % cD) / (diagvec(JJ * solve(-H, JJ.t(), solve_opts::force_sym)));
	double E = (double)sum(nu / (nu - 2));

	// ANOVA
	wald[resid].lnL_REML = xx_reml.lnL;
	wald[resid].df1 = np;
	wald[resid].df2 = 2 * E / (E - np);
	wald[resid].F = sum(pow(mP.t() * iEF01.rows(kXI, kXIG - 1), 2) / cD) / np;
	wald[resid].mlogP = MinusLogPF(wald[resid].F, wald[resid].df1, wald[resid].df2);
	wald[resid].P = exp(-wald[resid].mlogP);
}
*/

/* Perform Wald Test */
template<typename REAL>
TARGET void GWAS<REAL>::WaldTest(CPOINT& xx_reml, int resid)
{
	lnL_reml_profile1(fabs(xx_reml.real_space[0]), NULL, NULL, false);
	iE01 = inv_sympd(mE01);
	rmat betaG = iEF01.rows(kXI, kXIG - 1);

	// Wald test
	wald[resid].lnL_REML = xx_reml.lnL;
	//wald[resid].df = kXIG - kXI;
	wald[resid].W = trace(betaG.t() * inv_sympd(iE01.submat(kXI, kXI, kXIG - 1, kXIG - 1) * sig2_from_r) * betaG);  
	wald[resid].mlogP = MinusLogPChi2(wald[resid].W, degrees_of_freedom[resid]);
	wald[resid].P = exp(-wald[resid].mlogP);
}

/* Perform GWAS */
template<typename REAL>
TARGET void GWAS<REAL>::Prepare()
{
	tid = -1;

	if (GWAS_DEBUG)
	{
		//Read from ped file, gwas_debug
		FILE* f1 = fopen("C:\\Desktop\\gwas\\gemma_clone\\clone.ped", "rb");
		FILE* f2 = fopen("C:\\Desktop\\gwas\\gemma_clone\\clone.env", "rb");

		// Allocate memory, gwas_debug
		nploidy = 1;
		ploidy[0] = 2;
		ploidyidx[2] = 0;

		n = 100;
		nloc = m = 1000;

		kX = 3;
		kXI = 4;
		kXIG = 4;

		oY = zeros<rmat>(n, 1);  //
		oX = zeros<rmat>(n, kX);  //

		// Read data, gwas_debug
		noffset = new int64[m + 1];
		Gb = new byte[n * m * 2];
		char buf[100];

		for (int i = 0; i < n; ++i)
			for (int j = 0; j < kX; ++j)
				fscanf(f2, sizeof(REAL) == 8 ? "%lf" : "%f", &oX(i, j));

		for (int i = 0; i < n; ++i)
		{
			fscanf(f1, "%s", buf); // Family ID
			fscanf(f1, "%s", buf); // Individual ID
			fscanf(f1, "%s", buf); // Paternal ID
			fscanf(f1, "%s", buf); // Maternal ID
			fscanf(f1, "%s", buf); // Sex
			fscanf(f1, sizeof(REAL) == 8 ? "%lf" : "%f", &oY(i)); // Phenotype

			// Read Genotype, gwas_debug
			for (int64 l = 0; l < m; ++l)
			{
				char c1;
				for (;;)
				{
					c1 = fgetc(f1);
					if (isgraph(c1))
					{
						fgetc(f1);
						break;
					}
				}
				Gb[(l * n + i) * 2 + 0] = c1;

				for (;;)
				{
					c1 = fgetc(f1);
					if (isgraph(c1))
					{
						fgetc(f1);
						break;
					}
				}
				Gb[(l * n + i) * 2 + 1] = c1;
			}
		}

		// count alleles, gwas_debug
		int coffset = 0;
		vector<byte> valleles;
		for (int64 l = 0; l < m; ++l)
		{
			int acount[256] = { 0 }, alleles[256] = { 0 };
			byte* gb = &Gb[l * n * 2];
			for (int i = 0; i < 2 * n; ++i)
				acount[gb[i]]++;

			nallele = 0;
			for (byte i = 0; i < 255; ++i)
				if (acount[i] && i != '0')
				{
					alleles[nallele++] = i;
					valleles.push_back(i);
				}

			noffset[l] = coffset;
			coffset += nallele;
		}
		noffset[m] = coffset;

		// create genotype matrix, gwas_debug
		oG = zeros<rmat>(n, noffset[m]);

		// read genotypes
		for (int64 l = 0; l < m; ++l)
		{
			int nals = noffset[l + 1] - noffset[l];
			int64 offset = noffset[l];
			int aid[256] = { 0 }, ntype = 0;
			byte* gb = &Gb[l * n * 2];

			for (int a = 0; a < nals; ++a)
				aid[valleles[offset + a]] = a;

			for (int i = 0; i < n; ++i)
			{
				byte a1 = gb[i * 2 + 0], a2 = gb[i * 2 + 1];
				if (a1 == '0' || a2 == '0')
					for (int a = 0; a < nals; ++a)
						oG(i, offset + a) = nan;
				else
				{
					ntype++;
					if (a1 == a2)
						oG(i, offset + aid[a1]) = 1;
					else
					{
						oG(i, offset + aid[a1]) = 0.5;
						oG(i, offset + aid[a2]) = 0.5;
					}
				}
			}
		}
		fclose(f1);
		fclose(f2);

		// fill missing, gwas_debug
		int64 nbatch = (m + gwas_batch_val - 1) / gwas_batch_val;
#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int64 batch = 0; batch < nbatch; ++batch)
		{
			int64 lst = batch * gwas_batch_val;
			int64 led = std::min((batch + 1) * gwas_batch_val, m);
			int64 colst = noffset[lst];
			int64 coled = noffset[led] - 1;

			rmat tG(oG.memptr() + colst * n, n, coled - colst + 1, false, true);
			Imputation<REAL>(0, gwas_imputeG_val, tG, lst, led, noffset, false, *(rmat*)NULL);
		}

		// estimate relatedness, gwas_debug
		if (gwas_restimator_val == 1 || gwas_restimator_val == 2)
		{
			rrow mW = zeros<rrow>(1, noffset[m]);
			for (int64 l = 0; l < m; ++l)
			{
				int na = noffset[l + 1] - noffset[l];
				for (int64 j = noffset[l]; j < noffset[l + 1]; ++j)
					mW(j) = na;
			}
			mW = sqrt(1.0 / (m * mW));

			rrow oGmean = mean(oG);
			rrow oGsd = stddev(oG, 1);
			rmat mGz = gwas_restimator_val == 1 ?
				(oG.each_row() - oGmean) :
				(oG.each_row() - oGmean).each_row() / oGsd;

			mGz.each_row() %= mW;

			oR = mGz * mGz.t();
		}

		GWAS_coltype.push_back("Ind");		GWAS_coltype.push_back("Normal");
		GWAS_colname.push_back("Ind");		GWAS_colname.push_back("Phenotype");

		// Convert proportion into dosage, debug
		if (gwas_dosage_val == 1)
			oG = oG * 2;
	}
	else
	{
		nploidy = 0;
		m = nloc;
		ifstream file(gwas_input_val);

		if (!file.is_open())
			Exit("Error: cannot open GWAS input file.");

		if (!ReadCsvLine(GWAS_colname, file))
			Exit("Error: cannot read first row of GWAS input file.");

		if (!ReadCsvLine(GWAS_coltype, file))
			Exit("Error: cannot read second row of GWAS input file.");

		if (GWAS_coltype.size() != GWAS_colname.size())
			Exit("Error: #cols are wrong in second row of GWAS input file.");

		if (GWAS_coltype[0] != "Ind" || GWAS_coltype.size() <= 2)
			Exit("Error: format error in GWAS input file.");
		
		for (int j = 1; j < GWAS_coltype.size(); ++j)
		{
			if (GWAS_coltype[j] == "Normal")		nresponse++;
			if (GWAS_coltype[j] == "Bernoulli")		nresponse++;
			if (GWAS_coltype[j] == "Poisson")		nresponse++;
			if (GWAS_coltype[j] == "Gamma")			nresponse++;
			if (GWAS_coltype[j] == "Numeric")		nnumeric++;
			if (GWAS_coltype[j] == "Nominal")		nnominal++;
		}

		n = cpop<REAL>->nind;
		for (int i = 0; i < n; ++i)
			GWAS_indid[cpop<REAL>->inds[i]->name] = i;

		vector<vector<string>> cells;
		for (int nr = 3; ; nr++)
		{
			vector<string> row;
			if (!ReadCsvLine(row, file)) break;

			if (row.size() != GWAS_coltype.size())
				Exit("Error: #cols are wrong in row %d of GWAS input file.", nr);

			cells.push_back(row);
		}
		file.close();

		if (cells.size() < (uint64)n)
			Exit("Error: #individual in GWAS input file mismatches that in target population.");

		// check individuals, gwas_release
		vector<int> csvids;//id of ith individual in the csv file
		byte* nread = new byte[n];
		SetZero(nread, n);
		for (int i = 0; i < n; ++i)
		{
			string& name = cells[i][0];

			if (GWAS_indid.find(name) == GWAS_indid.end())
				Exit("Error: individual %s in GWAS input file is not find in target population.", name.c_str());

			int rid = GWAS_indid[name];
			csvids.push_back(rid);

			if (nread[rid] == 0)
				nread[rid] = 1;
			else
				Exit("Error: individual %s in GWAS input file appears more than once at line %d.", name.c_str(), 3 + i);
		}
		DEL(nread);

		// count #level in each factor, gwas_release
		oY = zeros<rmat>(n, nresponse);
		oX = zeros<rmat>(n, 0);
		for (int j = 1, cy = 0, cx = 0; j < GWAS_coltype.size(); ++j)
		{
			if (GWAS_coltype[j] == "Normal" || GWAS_coltype[j] == "Bernoulli" || GWAS_coltype[j] == "Poisson" || GWAS_coltype[j] == "Gamma")
			{
				for (int i = 0; i < n; ++i)
				{
					int id = csvids[i];
					char* ptr = (char*)cells[i][j].c_str();
					oY(id, cy) = IsNaNStr(ptr) ? nan : ReadDoubleKeep(ptr);
				}
				cy++;
			}
			if (GWAS_coltype[j] == "Numeric")
			{
				rmat tX = zeros<rmat>(n, 1);
				for (int i = 0; i < n; ++i)
				{
					int id = csvids[i];
					char* ptr = (char*)cells[i][j].c_str();
					tX(id, 0) = IsNaNStr(ptr) ? nan : ReadDoubleKeep(ptr);
				}
				oX = join_horiz(oX, tX);
				cx++;
			}
			if (GWAS_coltype[j] == "Nominal")
			{
				umap<string, int> factorid;
				for (int i = 0; i < n; ++i)
				{
					if (!IsNaNStr((char*)cells[i][j].c_str()) && factorid.find(cells[i][j]) == factorid.end())
						factorid[cells[i][j]] = (int)factorid.size();
				}

				//this factor has more than one levels, gwas_release
				if (factorid.size() > 1)
				{
					rmat tX = zeros<rmat>(n, factorid.size());
					for (int i = 0; i < n; ++i)
					{
						int id = csvids[i];
						char* ptr = (char*)cells[i][j].c_str();
						if (IsNaNStr(ptr))
							tX.row(id) = nan;
						else
							tX(id, factorid[cells[i][j]]) = 1;
					}
					oX = join_horiz(oX, tX.cols(0, tX.n_cols - 2));
					cx += tX.n_cols - 1;
				}
			}
		}

		kX = oX.n_cols;
		kXI = kX;
		kXIG = kXI;

		noffset = new int64[m + 1];
		int64 coffset = 0;
		for (int64 l = 0; l < m; ++l)
		{
			noffset[l] = coffset;
			coffset += GetLoc(l).k;
		}
		noffset[m] = coffset;

		// create genotype matrix, gwas_release
		SetZero(GWAS_ploidy_presence, N_MAX_PLOIDY + 1);

		if (gwas_restimator_val <= 2)
		{
			//estimate relatedness by matrix operations, gwas_release
			tR = new rmat[g_nthread_val];
			for (int i = 0; i < g_nthread_val; ++i)
				tR[i] = zeros<rmat>(n, n);
			GWAS_batch_index = 0;

			RunThreads(&GWASMatrixRelatednessThread<REAL>, NULL, NULL,
				m + 10 * ((m + gwas_batch_val - 1) / gwas_batch_val) + g_nthread_val,
				m + 10 * ((m + gwas_batch_val - 1) / gwas_batch_val) + g_nthread_val,
				"\nEstimating relatedness coefficient between individuals:\n", g_nthread_val, true);

			oR = tR[0];
			tR[0].clear();
			for (int i = 1; i < g_nthread_val; ++i)
			{
				oR += tR[i];
				tR[i].clear();
			}
			DEL(tR);
		}
		else
		{
			// estimate relatedness by other estimators, gwas_release
			relatedness_loc_stat<REAL> = new LOCSTAT2<REAL>[m];
			SetZero(relatedness_loc_stat<REAL>, m);
			cpop<REAL>->GetLocStat2(relatedness_loc_stat<REAL>);
			oR = zeros<rmat>(n, n);

			RunThreads(&GWASRelatednessThread<REAL>, NULL, NULL, 
				n * (n + 1) / 2 + 10 * g_nthread_val,
				n * (n + 1) / 2 + 10 * g_nthread_val,
				"\nEstimating relatedness coefficient between individuals:\n", g_nthread_val, true);

			DEL(relatedness_loc_stat<REAL>);
		}

		//assign ploidy, gwas_release
		nploidy = 0;
		for (int v = 1; v <= N_MAX_PLOIDY; ++v)
		{
			if (GWAS_ploidy_presence[v])
			{
				ploidy[nploidy] = v;
				ploidyidx[v] = nploidy;
				nploidy++;
			}
		}
	}

	// Extract diagonal elements
	dR = oR.diag();

	// remove duplicated columns in oX
	RemoveDupCol(oX);
	kX = oX.n_cols;

	// add intercept columns
	if (gwas_intercept_val == 1)
	{
		oI = zeros<rmat>(n, nploidy);
		GENO_READER rt(0, 0);
		GENOTYPE* gtab = GetLoc(0).GetGtab();
		for (int i = 0; i < n; ++i)
			oI(i, ploidyidx[gtab[rt.Read()].Ploidy()]) = 1;
	}
	else
		oI = ones<rmat>(n, 1);

	// imputate X and Y
	ucol invalid = find(sum(oY == nan, 1) + sum(oX == nan, 1) + sum(oI == nan, 1)) + 1;
	if (invalid.n_elem > 0)
	{
		rmat D = oI.n_cols > 1 ? join_horiz(oY, oX, oI) : join_horiz(oY, oX);
		Imputation(0, gwas_imputeXY_val, D, 0, D.n_cols, NULL, false, oR);
		oY = D.cols(0, oY.n_cols - 1);
		oX = D.cols(oY.n_cols, oY.n_cols + oX.n_cols - 1);
		if (oI.n_cols > 1)
			oI = D.cols(oY.n_cols + oX.n_cols, oY.n_cols + oX.n_cols + oI.n_cols - 1);
	}

	// allocate result 
	wald = NULL; 
	lrt = NULL; 
	score = NULL;
	sample_size = new int64[nresponse * gwas_batch_val];
	degrees_of_freedom = new int64[nresponse * gwas_batch_val];
	if (gwas_test_val[1]) wald  = new RESULT_WALD <REAL>[nresponse * gwas_batch_val];
	if (gwas_test_val[2]) lrt   = new RESULT_LRT  <REAL>[nresponse * gwas_batch_val];
	if (gwas_test_val[3]) score = new RESULT_SCORE<REAL>[nresponse * gwas_batch_val];
}

/* Copy Is Ihash oG oI_lst */
template<typename REAL>
TARGET void GWAS<REAL>::CopyRef(GWAS<REAL>* ref)
{
	oG_offset0 = ref->oG_offset0;
	Is = ref->Is;
	Ihash = ref->Ihash;
	oI_lst = ref->oI_lst;
	REF_MAT(oG);
}

/* Perform GWAS for locus cl */
template<typename REAL>
TARGET void GWAS<REAL>::CalcLocus(int64 _l)
{
	cl = _l;
	int id = cl % gwas_batch_val;
	void* Param = (void*)this;

	nallele = noffset[cl + 1] - noffset[cl];
	n = GWAS_Root<REAL>.n;

	umap<HASH, rmat>& Isbase = Is[0];
	vector<HASH>& Ihashbase = Ihash[0];

	//read genotypes
	SetZero(GWAS_ploidy_presence, N_MAX_PLOIDY + 1);
	if (GWAS_DEBUG)
	{
		mG = oG.cols(noffset[cl], noffset[cl] + nallele - 2);
		mI = oI;
		GWAS_ploidy_presence[2] = 1;
		nploidy = 1;
		ploidy[0] = 2;
		ploidyidx[2] = 0;
	}
	else
	{
		//load imputated genotype matrix at this locus
		mG = oG.cols(noffset[cl] - oG_offset0, noffset[cl] + nallele - 2 - oG_offset0);

		//set ploidies
		nploidy = 0;
		SetZero(ploidy, N_MAX_PLOIDY + 1);
		SetZero(ploidyidx, N_MAX_PLOIDY + 1);

		//using different intercepts among ploidies
		if (GWAS_Root<REAL>.nploidy > 1 && gwas_intercept_val == 1)
		{
			/*
			rmat mI = zeros<rmat>(n, GWAS_Root<REAL>.nploidy);
			GENO_READER rt(cpop<REAL>->ind0id, cl);
			GENOTYPE* gtab = GetLoc(cl).GetGtab();

			for (int i = 0; i < n; ++i)
				mI(i, GWAS_Root<REAL>.ploidyidx[gtab[rt.Read()].Ploidy()]) = 1;
			*/

			// load saved mI
			mI = Isbase[Ihashbase[cl - oI_lst]];
			
			// ploidy and ploidyidx
			rrow sI = sum(mI);
			nploidy = sum(sI > 0);

			for (int i = 0, ip = 0; i < GWAS_Root<REAL>.nploidy; ++i)
			{
				if (sI(i) > 0)
				{
					int v = GWAS_Root<REAL>.ploidy[i];
					ploidy[ip] = v;
					ploidyidx[v] = ip;
					ip++;
				}
			}

			//remove empty col for mI
			mI = mI.cols(find(sI > 0));
		}
		//using the same intercept among ploidies
		else
		{
			mI = oI;
			nploidy = 1;
			SetVal(ploidy, GWAS_Root<REAL>.ploidy, N_MAX_PLOIDY + 1);
			SetVal(ploidyidx, GWAS_Root<REAL>.ploidyidx, N_MAX_PLOIDY + 1);
		}
	}

	// mark missing x and g
	ucol valid_xg = all(mG != nan, 1) && all(oX != nan, 1);
	
	// different slopes among ploidies
	if (gwas_slope_val == 1)
	{
		int np = mI.n_cols;
		rmat mGt = zeros<rmat>(n, np * mG.n_cols);
		for (int i = 0; i < np; ++i)
			mGt.cols(i * np, (i + 1) * np - 1) = mG.each_col() % mI.col(i);
		mG = mGt;
		RemoveDupCol(mG);
	}

	// enumerate each Y
	for (int yi = 1, yp = 0; yi < GWAS_coltype.size(); ++yi)
	{
		// normal distributed Y
		if (GWAS_coltype[yi] == "Normal")
		{
			int resid = nresponse * id + yp;
			RNG<double> rng(g_seed_val + cl * nresponse + yp, RNG_SALT_GWAS);
			ucol valid = find(valid_xg      && oY.col(yp) != nan);
			ucol misid = find(valid_xg == 0 || oY.col(yp) == nan);
			ucol yp2 = { (uint64) yp++ };
            
			sample_size[resid] = n = valid.n_elem;
			degrees_of_freedom[resid] = 0;

			// prepare matrices U V
			GWAS_MEM<REAL>& tmem = GWAS_MEM<REAL>::FindEntry(valid, misid, rng.Uniform());
			new (&cV)  rcol(tmem.buf, n, false, true);
            new (&mU)  rmat(tmem.buf + n, n, n, false, true);
            
			cVr		= zeros<rcol>(n, 1);
			cW01	= zeros<rcol>(n, 1);
			cr		= zeros<rcol>(n, 1);
			cW12	= zeros<rcol>(n, 1);

			cR = dR(valid);
			mY = oY.submat(valid, yp2);
            mPt = mU.t() * mY;
            
			mXI = join_horiz(oX.rows(valid), mI.rows(valid));
			RemoveDupCol(mXI);
			kXI = mXI.n_cols;

			// insufficient information, terminate
			if (n <= kXI) continue;

			mXIG = join_horiz(mXI, mG.rows(valid));
			RemoveDupCol(mXIG);
			kXIG = mXIG.n_cols;
			degrees_of_freedom[resid] = kXIG - kXI;

			// genotype provide no additional information, terminate
			if (kXI == kXIG || n <= kXIG) continue;

			CPOINT xx_reml, xx_ml0, xx_ml;
			GWAS_NULL_MODEL<REAL>* null_model = NULL;
            
			if (gwas_imputeG_val != 10 && lnL0_Res.find(yp) != lnL0_Res.end())
				null_model = &lnL0_Res[yp];

			if (gwas_test_val[2] || gwas_test_val[3])
			{
				SetIndVar(mXI);

				if (null_model)
					xx_ml0 = null_model->xx_ml0;
				else
				{
					xx_ml0 = CPOINT::GradientDescent(Param, GWAS::lnL_ml_profile2, 2);
					xx_ml0 = BestPoint(xx_ml0, CPOINT::GradientDescent(Param, GWAS::lnL_ml_profile1, 1));
					
					if (gwas_imputeG_val != 10)
					{
						lnL_ml_profile1(fabs(xx_ml0.real_space[0]), NULL, NULL, false);
						null_model = &lnL0_Res[yp];
						null_model->r = fabs(xx_ml0.real_space[0]);
						null_model->sig = fabs(sig2_from_r);
						null_model->xx_ml0 = xx_ml0;
						null_model->ma0 = iEF01;
					}
				}
			}
            
			// Score test
			if (gwas_test_val[3])
				ScoreTest(xx_ml0, resid, null_model);
			else if (gwas_test_val[1] || gwas_test_val[2])
				SetIndVar(mXIG);

			// LRT test
			if (gwas_test_val[2])
			{
				xx_ml = CPOINT::GradientDescent(Param, GWAS::lnL_ml_profile2, 2);
				xx_ml = BestPoint(xx_ml, CPOINT::GradientDescent(Param, GWAS::lnL_ml_profile1, 1, false, xx_ml.unc_space));
				
				if (xx_ml.lnL < xx_ml0.lnL)
				{	
					xx_ml = BestPoint(xx_ml, CPOINT::GradientDescent(Param, GWAS::lnL_ml_profile1, 1));
					LRTTest(xx_ml, xx_ml0, resid);
				}
				else
				{
					LRTTest(xx_ml, xx_ml0, resid);
					if (lrt[resid].P < 0.01)
					{
						xx_ml = BestPoint(xx_ml, CPOINT::GradientDescent(Param, GWAS::lnL_ml_profile1, 1));
						if (lrt[resid].lnL_ML != xx_ml.lnL)
							LRTTest(xx_ml, xx_ml0, resid);
					}
				}
			}
            
			// Wald test
			if (gwas_test_val[1])
			{
				xx_reml = CPOINT::GradientDescent(Param, GWAS::lnL_reml_profile2, 2);
				xx_reml = BestPoint(xx_reml, CPOINT::GradientDescent(Param, GWAS::lnL_reml_profile1, 1, false, xx_reml.unc_space));
				WaldTest(xx_reml, resid);

				if (wald[resid].P < 0.01)
				{
					xx_reml = BestPoint(xx_reml, CPOINT::GradientDescent(Param, GWAS::lnL_reml_profile1, 1));
					if (wald[resid].lnL_REML != xx_reml.lnL)
						WaldTest(xx_reml, resid);
				}
			}			
            
			tmem.nref--;
			tmem.~GWAS_MEM();
		}
	}
}

/* Read a batch of genotypes */
template<typename REAL>
TARGET void GWAS<REAL>::ReadBatch(int64 lst, int64 led)
{
	if (GWAS_DEBUG) return;

	if (gwas_test_val[1]) SetVal((REAL*)wald,  (REAL)NAN, gwas_batch_val * sizeof(RESULT_WALD <REAL>) / sizeof(REAL));
	if (gwas_test_val[2]) SetVal((REAL*)lrt,   (REAL)NAN, gwas_batch_val * sizeof(RESULT_LRT  <REAL>) / sizeof(REAL));
	if (gwas_test_val[3]) SetVal((REAL*)score, (REAL)NAN, gwas_batch_val * sizeof(RESULT_SCORE<REAL>) / sizeof(REAL));

	openblas_set_num_threads(g_nthread_val);
	oG_offset0 = noffset[lst]; oI_lst = lst; Isbase.clear(); Ihashbase.clear();
	Is = &Isbase; Ihash = &Ihashbase;
	int64 tnallele = noffset[led] - noffset[lst];

	oG = zeros<rmat>(n, tnallele);

	for (int64 l = lst; l < led; ++l)
	{
		GENO_READER rt(cpop<REAL>->ind0id, l);
		GENOTYPE* gtab = GetLoc(l).GetGtab();
		int64 nals = noffset[l + 1] - noffset[l];
		int64 pals = noffset[l] - oG_offset0;
		rmat tI = zeros<rmat>(n, GWAS_Root<REAL>.nploidy);

		for (int i = 0; i < n; ++i)
		{
			GENOTYPE& gt = gtab[rt.Read()];
			int vi = gt.Ploidy();
			if (vi > 2)
				vi = vi;
			REAL iv = 1.0 / vi;

			if (gt.Nalleles() == 0)
			{
				for (int a = 0; a < nals; ++a)
					oG(i, pals + a) = nan;
			}
			else
			{
				ushort* als = gt.GetAlleleArray();
				for (int a = 0; a < vi; ++a)
					oG(i, pals + als[a]) += iv;
			}

			tI(i, GWAS_Root<REAL>.ploidyidx[vi]) = 1;
		}

		if (GWAS_Root<REAL>.nploidy > 1 && gwas_intercept_val == 1)
		{
			HASH ha = HashString((char*)tI.memptr(), tI.n_elem * sizeof(REAL));
			if (Isbase.find(ha) == Isbase.end())
				Isbase[ha] = tI;
			Ihashbase.push_back(ha);
		}
	}

	Imputation<REAL>(1, gwas_imputeG_val, oG, lst, led, noffset, false, oR);
	openblas_set_num_threads(1);
}

/* Set independent variables */
template<typename REAL>
TARGET void GWAS<REAL>::SetIndVar(rmat& _mX)
{
	k = _mX.n_cols;
	nk = n - k;
	if (mQt.n_rows != n || mQt.n_cols != k)
		mQt = zeros<rmat>(n, k);

    MatrixMul(mQt.memptr(), mU.memptr(), _mX.memptr(), n, n, k);

	if (k != mE01.n_cols)
	{
		mE01 = zeros<rmat>(k, k);
		mF01 = zeros<rmat>(k, 1);
		mJ01 = zeros<rmat>(1, 1);

		mE12 = zeros<rmat>(k, k);
		mF12 = zeros<rmat>(k, 1);
		mJ12 = zeros<rmat>(1, 1);

		mE23 = zeros<rmat>(k, k);
		mF23 = zeros<rmat>(k, 1);
		mJ23 = zeros<rmat>(1, 1);
	}
}

/* Set Vr W01 E01 F01 J01 */
template<typename REAL>
TARGET void GWAS<REAL>::SetVar01(double r2)
{
	//cVr = 1.0 + cV * r2;
	SetVal(cVr.memptr(), (REAL)1, n);
	AddProd(cVr.memptr(), cV.memptr(), r2, n);

	//cW01 = 1.0 / cVr;
	Div(cW01.memptr(), 1, cVr.memptr(), n);

	DiagQuadForm(mE01.memptr(), mQt.memptr(), cW01.memptr(), k, n);
	DiagQuadForm(mF01.memptr(), mQt.memptr(), cW01.memptr(), mPt.memptr(), k, n);
	DiagQuadForm(mJ01.memptr(), mPt.memptr(), cW01.memptr(), n);

	/*
	mE01 = mQt.t() * diagmat(cW01) * mQt;
	mF01 = mQt.t() * diagmat(cW01) * mPt;
	mJ01 = mPt.t() * diagmat(cW01) * mPt;
	*/
}

/* Set Vr W01 E01 F01 J01 W12 E12 F12 J12 */
template<typename REAL>
TARGET void GWAS<REAL>::SetVar12(double r2)
{
	//cVr = 1.0 + cV * r2;
	SetVal(cVr.memptr(), (REAL)1, n);
	AddProd(cVr.memptr(), cV.memptr(), r2, n);

	//cW01 = 1.0 / cVr;
	Div(cW01.memptr(), 1, cVr.memptr(), n);

	Mul(cr.memptr(), cV.memptr(), cW01.memptr(), n);
	Mul(cW12.memptr(), cr.memptr(), cW01.memptr(), n);

	DiagQuadForm(mE01.memptr(), mQt.memptr(), cW01.memptr(), k, n);
	DiagQuadForm(mF01.memptr(), mQt.memptr(), cW01.memptr(), mPt.memptr(), k, n);
	DiagQuadForm(mJ01.memptr(), mPt.memptr(), cW01.memptr(), n);

	DiagQuadForm(mE12.memptr(), mQt.memptr(), cW12.memptr(), k, n);
	DiagQuadForm(mF12.memptr(), mQt.memptr(), cW12.memptr(), mPt.memptr(), k, n);
	DiagQuadForm(mJ12.memptr(), mPt.memptr(), cW12.memptr(), n);

	return;
}

/* Optimizer Warpper */
template<typename REAL>
TARGET double GWAS<REAL>::lnL_reml_profile1(void* Param, CPOINT& xx, rmat& G, rmat& H)
{
	xx.real_space[0] = xx.unc_space[0];
	return ((GWAS<REAL>*)Param)->lnL_reml_profile1(xx.real_space[0], &G, &H, false);
}

/* Optimizer Warpper */
template<typename REAL>
TARGET double GWAS<REAL>::lnL_ml_profile1(void* Param, CPOINT& xx, rmat& G, rmat& H)
{
	xx.real_space[0] = xx.unc_space[0];
	return ((GWAS<REAL>*)Param)->lnL_ml_profile1(xx.real_space[0], &G, &H, false);
}

/* Optimizer Warpper */
template<typename REAL>
TARGET double GWAS<REAL>::lnL_reml_profile2(void* Param, CPOINT& xx, rmat& G, rmat& H)
{
	xx.real_space[0] = xx.unc_space[0];
	xx.real_space[1] = xx.unc_space[1];
	return ((GWAS<REAL>*)Param)->lnL_reml_profile2(xx.real_space[0], xx.real_space[1], &G, &H, false);
}

/* Optimizer Warpper */
template<typename REAL>
TARGET double GWAS<REAL>::lnL_ml_profile2(void* Param, CPOINT& xx, rmat& G, rmat& H)
{
	xx.real_space[0] = xx.unc_space[0];
	xx.real_space[1] = xx.unc_space[1];
	return ((GWAS<REAL>*)Param)->lnL_ml_profile2(xx.real_space[0], xx.real_space[1], &G, &H, false);
}

/* LogLikelihood using REML criterion */
template<typename REAL>
TARGET double GWAS<REAL>::lnL_reml_profile1(double r, rmat* G, rmat* H, bool issquare)
{
	double r2 = r * r;

	if (G)
	{
		SetVar12(r2);
		iE01 = inv_sympd(mE01);
		iEF01 = iE01 * mF01;
	}
	else
	{
		SetVar01(r2);
		iEF01 = solve(mE01, mF01, solve_opts::force_sym);
	}

	double ldVr = LogProd(cVr.memptr(), n);
	double sig2 = sig2_from_r = trace(mJ01 - mF01.t() * iEF01) / nk;
	double f = -0.5 * (nk * 1.83787706640934548 + nk * log(sig2) + ldVr + real(log_det(mE01)) + nk);

	if (G)
	{
		rmat iJ = inv_sympd(mJ01 - mF01.t() * iEF01);
		rmat iS = iJ * (-mJ12 + 2 * iEF01.t() * mF12 - iEF01.t() * mE12 * iEF01);
		G[0] = -0.5 * (sum(cr) - trace(iE01 * mE12) + nk * iS);

		if (H)
		{
			rcol cW23 = cW12 % cr;
			DiagQuadForm(mE23.memptr(), mQt.memptr(), cW23.memptr(), k, n);
			DiagQuadForm(mF23.memptr(), mQt.memptr(), cW23.memptr(), mPt.memptr(), k, n);
			DiagQuadForm(mJ23.memptr(), mPt.memptr(), cW23.memptr(), n);

			rmat iEF12 = iE01 * mF12, iEF23 = iE01 * mF23;
			H[0] = -0.5 * (-nk * iS * iS + 2 * nk * iJ * (mJ23 - mF12.t() * iEF12 + iEF01.t() * (mE12 * (2 * iEF12 - iE01 * mE12 * iEF01) + mE23 * iEF01 - 2 * mF23)) - sum(cr % cr) + 2 * trace(iE01 * mE23) - trace(iE01 * mE12 * iE01 * mE12));
		}

		if (!issquare)
		{
			if (H) H[0] = 4 * r2 * H[0] + 2 * G[0];
			G[0] = 2 * r * G[0];
		}
	}

	return f;
}

/* LogLikelihood using ML criterion */
template<typename REAL>
TARGET double GWAS<REAL>::lnL_ml_profile1(double r, rmat* G, rmat* H, bool issquare)
{
	double r2 = r * r;

	if (G)
	{
		SetVar12(r2);
		iE01 = inv_sympd(mE01);
		iEF01 = iE01 * mF01;
	}
	else
	{
		SetVar01(r2);
		iEF01 = solve(mE01, mF01, solve_opts::force_sym);
	}

	double ldVr = LogProd(cVr.memptr(), n);
	double sig2 = sig2_from_r = trace(mJ01 - mF01.t() * iEF01) / n;
	double f = -0.5 * (n * 1.83787706640934548 + n * log(sig2) + ldVr + n);

	if (G)
	{
		rmat iJ = inv_sympd(mJ01 - mF01.t() * iEF01);
		rmat iS = iJ * (-mJ12 + 2 * iEF01.t() * mF12 - iEF01.t() * mE12 * iEF01);
		G[0] = -0.5 * (sum(cr) + n * iS);

		if (H)
		{
			rcol cW23 = cW12 % cr;
			DiagQuadForm(mE23.memptr(), mQt.memptr(), cW23.memptr(), k, n);
			DiagQuadForm(mF23.memptr(), mQt.memptr(), cW23.memptr(), mPt.memptr(), k, n);
			DiagQuadForm(mJ23.memptr(), mPt.memptr(), cW23.memptr(), n);

			rmat iEF12 = iE01 * mF12, iEF23 = iE01 * mF23;
			H[0] = -0.5 * (-n * iS * iS + 2 * n * iJ * (mJ23 - mF12.t() * iEF12 + iEF01.t() * (mE12 * (2 * iEF12 - iE01 * mE12 * iEF01) + mE23 * iEF01 - 2 * mF23)) - sum(cr % cr));
		}

		if (!issquare)
		{
			if (H) H[0] = 4 * r2 * H[0] + 2 * G[0];
			G[0] = 2 * r * G[0];
		}
	}

	return f;
}

/* LogLikelihood using REML criterion */
template<typename REAL>
TARGET double GWAS<REAL>::lnL_reml_profile2(double r, double sig, rmat* G, rmat* H, bool issquare)
{
	double r2 = r * r, sig2 = sig * sig;

	if (G)
	{
		SetVar12(r2);
		iE01 = inv_sympd(mE01);
		iEF01 = iE01 * mF01;
	}
	else
	{
		SetVar01(r2);
		iEF01 = solve(mE01, mF01, solve_opts::force_sym);
	}

	if (!IsNormal(sig))
	{
		sig2 = sig2_from_r = trace(mJ01 - mF01.t() * iEF01) / nk;
		sig = sqrt(sig2);
	}

	double ldVr = LogProd(cVr.memptr(), n);
	double f = -0.5 * (nk * 1.83787706640934548 + nk * log(sig2) + ldVr + real(log_det(mE01)) + trace(mJ01) / sig2 - trace(mF01.t() * iEF01) / sig2);

	if (G)
	{
		double isig2 = 1 / sig2, isig4 = 1 / (sig2 * sig2), isig6 = 1 / (sig2 * sig2 * sig2);

		rmat iEF12 = iE01 * mF12;

		G[0] = zeros<rmat>(2, 1);
		G[0](0, 0) = -0.5 * trace(sum(cr) - trace(iE01 * mE12) - isig2 * mJ12 + (2 * isig2) * (mF01.t() * iEF12) - isig2 * (iEF01.t() * mE12 * iEF01));
		G[0](1, 0) = -0.5 * trace(isig2 * nk - isig4 * mJ01 + isig4 * (mF01.t() * iEF01));

		if (H)
		{
			rcol cW23 = cW12 % cr;
			DiagQuadForm(mE23.memptr(), mQt.memptr(), cW23.memptr(), k, n);
			DiagQuadForm(mF23.memptr(), mQt.memptr(), cW23.memptr(), mPt.memptr(), k, n);
			DiagQuadForm(mJ23.memptr(), mPt.memptr(), cW23.memptr(), n);

			H[0] = zeros<rmat>(2, 2);
			H[0](0, 0) = -0.5 * trace(-sum(cr % cr) - trace(iE01 * mE12 * iE01 * mE12) + 2 * trace(iE01 * mE23) + (2 * isig2) * (mJ23 - mF12.t() * iEF12 + iEF01.t() * (mE12 * (2 * iEF12 - iE01 * mE12 * iEF01) + mE23 * iEF01 - 2 * mF23)));
			H[0](0, 1) = -0.5 * trace(isig4 * mJ12 - (2 * isig4) * (mF01.t() * iEF12) + isig4 * (iEF01.t() * mE12 * iEF01));
			H[0](1, 0) = H[0](0, 1);
			H[0](1, 1) = -0.5 * trace(-isig4 * nk + (2 * isig6) * mJ01 - (2 * isig6) * (mF01.t() * iEF01));
		}

		if (!issquare)
		{
			rcol mtheta = { {(REAL)r}, {(REAL)sig} };
			if (H) H[0] = 4 * H[0] % (mtheta * mtheta.t()) + 2 * diagmat(G[0]);
			G[0] = 2 * G[0] % mtheta;
		}
	}

	return f;
}

/* LogLikelihood using ML criterion */
template<typename REAL>
TARGET double GWAS<REAL>::lnL_ml_profile2(double r, double sig, rmat* G, rmat* H, bool issquare)
{
	double r2 = r * r, sig2 = sig * sig;

	if (G)
	{
		SetVar12(r2);
		iE01 = inv_sympd(mE01);
		iEF01 = iE01 * mF01;
	}
	else
	{
		SetVar01(r2);
		iEF01 = solve(mE01, mF01, solve_opts::force_sym);
	}

	if (!IsNormal(sig))
	{
		sig2 = sig2_from_r = trace(mJ01 - mF01.t() * iEF01) / n;
		sig = sqrt(sig2);
	}

	double ldVr = LogProd(cVr.memptr(), n);
	double f = -0.5 * (n * 1.83787706640934548 + n * log(sig2) + ldVr + trace(mJ01) / sig2 - trace(mF01.t() * iEF01) / sig2);

	if (G)
	{
		double isig2 = 1 / sig2, isig4 = 1 / (sig2 * sig2), isig6 = 1 / (sig2 * sig2 * sig2);

		rmat iEF12 = iE01 * mF12;

		G[0] = zeros<rmat>(2, 1);
		G[0](0, 0) = -0.5 * trace(sum(cr) - isig2 * mJ12 + (2 * isig2) * (mF01.t() * iEF12) - isig2 * (iEF01.t() * mE12 * iEF01));
		G[0](1, 0) = -0.5 * trace(isig2 * n - isig4 * mJ01 + isig4 * (mF01.t() * iEF01));

		if (H)
		{
			rcol cW23 = cW12 % cr;
			DiagQuadForm(mE23.memptr(), mQt.memptr(), cW23.memptr(), k, n);
			DiagQuadForm(mF23.memptr(), mQt.memptr(), cW23.memptr(), mPt.memptr(), k, n);
			DiagQuadForm(mJ23.memptr(), mPt.memptr(), cW23.memptr(), n);

			H[0] = zeros<rmat>(2, 2);
			H[0](0, 0) = -0.5 * trace(-sum(cr % cr) + (2 * isig2) * (mJ23 - mF12.t() * iEF12 + iEF01.t() * (mE12 * (2 * iEF12 - iE01 * mE12 * iEF01) + mE23 * iEF01 - 2 * mF23)));
			H[0](0, 1) = -0.5 * trace(isig4 * mJ12 - (2 * isig4) * (mF01.t() * iEF12) + isig4 * (iEF01.t() * mE12 * iEF01));
			H[0](1, 0) = H[0](0, 1);
			H[0](1, 1) = -0.5 * trace(-isig4 * n + (2 * isig6) * mJ01 - (2 * isig6) * (mF01.t() * iEF01));
		}

		if (!issquare)
		{
			rcol mtheta = { {(REAL)r}, {(REAL)sig} };
			if (H) H[0] = 4 * H[0] % (mtheta * mtheta.t()) + 2 * diagmat(G[0]);
			G[0] = 2 * G[0] % mtheta;
		}
	} 

	return f;
}

/* LogLikelihood without profile */
template<typename REAL>
TARGET double GWAS<REAL>::lnL2(double r, double sig, rmat& ma, rmat* G, rmat* H)
{
	double r2 = r * r, sig2 = sig * sig;

	SetVar01(r2);

	double isig2 = 1 / sig2, isig4 = 1 / (sig2 * sig2), isig6 = 1 / (sig2 * sig2 * sig2);
	double ldVr = LogProd(cVr.memptr(), n);
	double f = -0.5 *
		(
			n * 1.83787706640934548 + ldVr + isig2 * trace(mJ01) - 2 * isig2 * trace(mF01.t() * ma) + isig2 * trace(ma.t() * mE01 * ma)
			);

	G[0] = isig2 * (mF01 - mE01 * ma);
	H[0] = - isig2 * mE01;

	return f;
}

/* Variance-covariance matrix by OPG estimator */
template<typename REAL>
TARGET rmat GWAS<REAL>::GetV_OPG2(double r, double sig, rmat& beta, rmat& X)
{
	double r2 = r * r;
	rmat G = (X.each_col() % ((mY - X * beta) / (sig + (sig * r2) * cR)));
	return inv_sympd(G.t() * G);
}

/* Variance-covariance matrix by Hessian estimator */
template<typename REAL>
TARGET rmat GWAS<REAL>::GetV_Hessian2(double r, double sig, rmat& X)
{
	double sig2 = sig * sig, r2 = r * r;
	rmat H = X.t() * (X.each_col() % (-1 / (sig2 + sig2 * r2 * cR)));
	return inv_sympd(-H);
}

#ifdef _xxx
/* LogLikelihood using REML criterion */
template<typename REAL>
TARGET double GWAS<REAL>::lnL_reml_profile2(double r, double sig)
{
	double r2 = r * r, sig2 = sig * sig;

	SetVar01(r2);

	double ldVr = LogProd(cVr.memptr(), n);
	iEF01 = solve(mE01, mF01, solve_opts::force_sym);

	if (!IsNormal(sig))
	{
		sig2 = sig2_from_r = trace(mJ01 - mF01.t() * iEF01) / nk;
		sig = sqrt(sig2);
	}

	double f = -0.5 * (nk * 1.83787706640934548 + nk * log(sig2) + ldVr + real(log_det(mE01)) + trace(mJ01) / sig2 - trace(mF01.t() * iEF01) / sig2);

	return f;
}

/* LogLikelihood using REML criterion */
template<typename REAL>
TARGET double GWAS<REAL>::lnL_reml_profile3(double sig, double tau)
{
	double sig2 = sig * sig, tau2 = tau * tau;

	SetVal(cVr.memptr(), (REAL)sig2, n);
	AddProd(cVr.memptr(), cV.memptr(), tau2, n);
	Div(1, cVr.memptr(), cW01.memptr(), n);
	double ldVr = LogProd(cVr.memptr(), n);

	mE01 = mQt.t() * diagmat(cW01) * mQt;
	mF01 = mQt.t() * diagmat(cW01) * mPt;
	mJ01 = mPt.t() * diagmat(cW01) * mPt;

	rmat iEF01 = solve(mE01, mF01, solve_opts::force_sym);

	double f = -0.5 * (nk * 1.83787706640934548 + ldVr + real(log_det(mE01)) + trace(mJ01) - SumProd(iEF01.memptr(), mF01.memptr(), k));

	return f;
}

/* LogLikelihood using REML criterion */
template<typename REAL>
TARGET double GWAS<REAL>::lnL_reml_profile3(double sig, double tau, rmat& G, rmat& H, bool issquare)
{
	double sig2 = sig * sig, tau2 = tau * tau;
	cVr = tau2 * cV + sig2;

	rmat mQ = mQt.t(), mP = mPt.t();
	rcol iVr1 = 1 / cVr, cV2 = cV % cV, iVr2 = iVr1 % iVr1, iVr3 = iVr2 % iVr1;
	rmat mW01 = diagmat(iVr1), mW11 = diagmat(cV % iVr1);
	rmat mW02 = diagmat(iVr2), mW12 = diagmat(cV % iVr2), mW22 = diagmat(cV2 % iVr2);
	rmat mW03 = diagmat(iVr3), mW13 = diagmat(cV % iVr3), mW23 = diagmat(cV2 % iVr3);

	// 5. Assign Eij, Fij, Jij, slow
	rmat mE01 = mQ * mW01 * mQt, mE02 = mQ * mW02 * mQt, mE12 = mQ * mW12 * mQt, mE03 = mQ * mW03 * mQt, mE13 = mQ * mW13 * mQt, mE23 = mQ * mW23 * mQt;
	rmat mJ01 = mP * mW01 * mPt, mJ02 = mP * mW02 * mPt, mJ12 = mP * mW12 * mPt, mJ03 = mP * mW03 * mPt, mJ13 = mP * mW13 * mPt, mJ23 = mP * mW23 * mPt;
	rmat mF01 = mQ * mW01 * mPt, mF02 = mQ * mW02 * mPt, mF12 = mQ * mW12 * mPt, mF03 = mQ * mW03 * mPt, mF13 = mQ * mW13 * mPt, mF23 = mQ * mW23 * mPt;

	// 6. Assign
	rmat mA = inv_sympd(mE01);
	rmat mF01tA = mF01.t() * mA, mF02tA = mF02.t() * mA, mF12tA = mF12.t() * mA;
	rmat mE02A = mE02 * mA, mE12A = mE12 * mA, mE03A = mE03 * mA, mE13A = mE13 * mA, mE23A = mE23 * mA;
	rmat mE02A02A = mE02A * mE02A, mE02A12A = mE02A * mE12A, mE12A12A = mE12A * mE12A;
	rmat mFm = mF12 - mE12A * mF01;

	// 7. Assign G and H
	G = zeros<rmat>(2, 1);
	H = zeros<rmat>(2, 2);
	G(0, 0) = -0.5 * (trace(mW01) - trace(mE02A) - trace(mJ02) - trace(mF01tA * (mE02A * mF01 - 2 * mF02)));
	G(1, 0) = -0.5 * (trace(mW11) - trace(mE12A) - trace(mJ12) - trace(mF01tA * (mE12A * mF01 - 2 * mF12)));

	H(0, 0) = -0.5 * (-trace(mW02) + 2 * trace(mE03A) - trace(mE02A02A) + trace(2 * mJ03 - 2 * mF01tA * ((mE02A02A - mE03A) * mF01 - mE02A * mF02 + 2 * mF03) - 2 * trace(mF02tA * (mF02 - mE02A * mF01))));
	H(0, 1) = -0.5 * (-trace(mW12) + 2 * trace(mE13A) - trace(mE02A12A) + trace(2 * mJ13 - 2 * mF01tA * ((mE02A12A - mE13A) * mF01 - mE02A * mF12 + 2 * mF13) - 2 * trace(mF02tA * mFm)));
	H(1, 1) = -0.5 * (-trace(mW22) + 2 * trace(mE23A) - trace(mE12A12A) + trace(2 * mJ23 - 2 * mF01tA * ((mE12A12A - mE23A) * mF01 - mE12A * mF12 + 2 * mF23) - 2 * trace(mF12tA * mFm)));
	H(1, 0) = H(0, 1);

	if (!issquare)
	{
		rcol mtheta = { {(REAL)sig}, {(REAL)tau} };
		H = 4 * H % (mtheta * mtheta.t()) + 2 * diagmat(G);
		G = 2 * G % mtheta;
	}

	double ldVr = LogProd(cVr.memptr(), n);
	double f = -0.5 * (nk * 1.83787706640934548 + ldVr + real(log_det(mE01)) - trace(mF01tA * mF01 - mJ01));

	return f;
}

/* LogLikelihood using ML criterion */
template<typename REAL>
TARGET double GWAS<REAL>::lnL_ml_profile1(double r)
{
	double r2 = r * r;

	SetVar01(r2);

	iEF01 = solve(mE01, mF01, solve_opts::force_sym);

	double ldVr = LogProd(cVr.memptr(), n);
	double sig2 = sig2_from_r = trace(mJ01 - mF01.t() * iEF01) / n;
	double f = -0.5 * (n * 1.83787706640934548 + n * log(sig2) + ldVr + n);

	return f;
}

/* LogLikelihood using ML criterion */
template<typename REAL>
TARGET double GWAS<REAL>::lnL_ml_profile2(double r, double sig)
{
	double r2 = r * r, sig2 = sig * sig;

	SetVar01(r2);

	double ldVr = LogProd(cVr.memptr(), n);
	iEF01 = solve(mE01, mF01, solve_opts::force_sym);

	if (!IsNormal(sig))
	{
		sig2 = sig2_from_r = trace(mJ01 - mF01.t() * iEF01) / n;
		sig = sqrt(sig2);
	}

	double f = -0.5 * (n * 1.83787706640934548 + n * log(sig2) + ldVr + trace(mJ01) / sig2 - trace(mF01.t() * iEF01) / sig2);

	return f;
}

/* LogLikelihood using ML criterion */
template<typename REAL>
TARGET double GWAS<REAL>::lnL_ml_profile3(double sig, double tau)
{
	double sig2 = sig * sig, tau2 = tau * tau;

	SetVal(cVr.memptr(), (REAL)sig2, n);
	AddProd(cVr.memptr(), cV.memptr(), tau2, n);
	Div(1, cVr.memptr(), cW01.memptr(), n);
	double ldVr = LogProd(cVr.memptr(), n);

	mE01 = mQt.t() * diagmat(cW01) * mQt;
	mF01 = mQt.t() * diagmat(cW01) * mPt;
	mJ01 = mPt.t() * diagmat(cW01) * mPt;

	rmat iEF01 = solve(mE01, mF01, solve_opts::force_sym);

	double f = -0.5 * (n * 1.83787706640934548 + ldVr + trace(mJ01) - SumProd(iEF01.memptr(), mF01.memptr(), k));

	return f;
}

/* LogLikelihood using ML criterion */
template<typename REAL>
TARGET double GWAS<REAL>::lnL_ml_profile3(double sig, double tau, rmat& G, rmat& H, bool issquare)
{
	double sig2 = sig * sig, tau2 = tau * tau;
	cVr = tau2 * cV + sig2;

	rcol iVr1 = 1 / cVr, cV2 = cV % cV, iVr2 = iVr1 % iVr1, iVr3 = iVr2 % iVr1;
	rmat mW01 = diagmat(iVr1), mW11 = diagmat(cV % iVr1);
	rmat mW02 = diagmat(iVr2), mW12 = diagmat(cV % iVr2), mW22 = diagmat(cV2 % iVr2);
	rmat mW03 = diagmat(iVr3), mW13 = diagmat(cV % iVr3), mW23 = diagmat(cV2 % iVr3);

	// 5. Assign Eij, Fij, Jij, slow
	rmat mE01 = mQt.t() * mW01 * mQt, mE02 = mQt.t() * mW02 * mQt, mE12 = mQt.t() * mW12 * mQt, mE03 = mQt.t() * mW03 * mQt, mE13 = mQt.t() * mW13 * mQt, mE23 = mQt.t() * mW23 * mQt;
	rmat mJ01 = mPt.t() * mW01 * mPt, mJ02 = mPt.t() * mW02 * mPt, mJ12 = mPt.t() * mW12 * mPt, mJ03 = mPt.t() * mW03 * mPt, mJ13 = mPt.t() * mW13 * mPt, mJ23 = mPt.t() * mW23 * mPt;
	rmat mF01 = mQt.t() * mW01 * mPt, mF02 = mQt.t() * mW02 * mPt, mF12 = mQt.t() * mW12 * mPt, mF03 = mQt.t() * mW03 * mPt, mF13 = mQt.t() * mW13 * mPt, mF23 = mQt.t() * mW23 * mPt;

	// 6. Assign
	rmat mA = inv_sympd(mE01);
	rmat mF01tA = mF01.t() * mA, mF02tA = mF02.t() * mA, mF12tA = mF12.t() * mA;
	rmat mE02A = mE02 * mA, mE12A = mE12 * mA, mE03A = mE03 * mA, mE13A = mE13 * mA, mE23A = mE23 * mA;
	rmat mE02A02A = mE02A * mE02A, mE02A12A = mE02A * mE12A, mE12A12A = mE12A * mE12A;
	rmat mFm = mF12 - mE12A * mF01;

	// 7. Assign G and H
	G = zeros<rmat>(2, 1);
	H = zeros<rmat>(2, 2);
	G(0, 0) = -0.5 * (trace(mW01) - trace(mJ02) - trace(mF01tA * (mE02A * mF01 - 2 * mF02)));
	G(1, 0) = -0.5 * (trace(mW11) - trace(mJ12) - trace(mF01tA * (mE12A * mF01 - 2 * mF12)));

	H(0, 0) = -0.5 * (-trace(mW02) + 2 * trace(mJ03) - 2 * trace(mF01tA * ((mE02A02A - mE03A) * mF01 - mE02A * mF02 + 2 * mF03)) - 2 * trace(mF02tA * (mF02 - mE02A * mF01)));
	H(0, 1) = -0.5 * (-trace(mW12) + 2 * trace(mJ13) - 2 * trace(mF01tA * ((mE02A12A - mE13A) * mF01 - mE02A * mF12 + 2 * mF13)) - 2 * trace(mF02tA * mFm));
	H(1, 1) = -0.5 * (-trace(mW22) + 2 * trace(mJ23) - 2 * trace(mF01tA * ((mE12A12A - mE23A) * mF01 - mE12A * mF12 + 2 * mF23)) - 2 * trace(mF12tA * mFm));
	H(1, 0) = H(0, 1);

	if (!issquare)
	{
		rmat mtheta = { {(REAL)sig}, {(REAL)tau} };
		H = 4 * H % (mtheta * mtheta.t()) + 2 * diagmat(G);
		G = 2 * G % mtheta;
	}

	double ldVr = LogProd(cVr.memptr(), n);
	double f = -0.5 * (n * 1.83787706640934548 + ldVr + trace(mJ01) - trace(mF01tA * mF01));

	return f;
}

/* LogLikelihood without profile */
template<typename REAL>
TARGET double GWAS<REAL>::lnL3(double sig, double tau, rmat& ma, rmat& G, rmat& H)
{
	double sig2 = sig * sig, tau2 = tau * tau;

	cVr = tau2 * cV + sig2;
	rcol ciVr1 = 1 / cVr;
	rmat mW01 = diagmat(ciVr1), mE01 = mQt.t() * mW01 * mQt, mF01 = mQt.t() * mW01 * mPt, mJ01 = mPt.t() * mW01 * mPt;

	double ldVr = LogProd(cVr.memptr(), n);
	double f = -0.5 *
		(
			n * 1.83787706640934548 + ldVr + trace(mJ01) - 2 * trace(mF01.t() * ma) + trace(ma.t() * mE01 * ma)
			);

	G = mF01 - mE01 * ma;
	H = -mE01;

	return f;
}

/* Variance-covariance matrix by OPG estimator */
template<typename REAL>
TARGET rmat GWAS<REAL>::GetV_OPG3(double sig, double tau, rmat& ma, rmat& mX)
{
	double sig2 = sig * sig, tau2 = tau * tau;
	rmat G = mX.each_col() % ((mY - mX * ma) / (tau2 * cR + sig2));
	return inv_sympd(G.t() * G);
}

/* Variance-covariance matrix by Hessian estimator */
template<typename REAL>
TARGET rmat GWAS<REAL>::GetV_Hessian3(double sig, double tau, rmat& ma, rmat& mX)
{
	double sig2 = sig * sig, tau2 = tau * tau;
	rmat H = -(mX.t() * (mX.each_col() % (1 / (tau2 * cR + sig2))));
	return inv_sympd(-H);
}
#endif

#endif

#ifndef _MEM

/* Find the entry from the dictionary */
template<typename REAL>
TARGET GWAS_MEM<REAL>& GWAS_MEM<REAL>::FindEntry(ucol& valid, ucol& misid, double rnd)
{
	umap<HASH, GWAS_MEM<REAL>>& tab = GWAS_umap<REAL>;
	map<double, GWAS_MEM<REAL>*>& list = GWAS_map<REAL>;

	HASH hash = HashString((char*)misid.memptr(), misid.n_elem * sizeof(misid(0)));

	//find entry with the same hash
	if (tab.find(hash) != tab.end())
	{
		// is there is the entry, add the ref count, update access time and move it to the end of lists
		GWAS_MEM<REAL>& tmem = tab[hash];
		while(!tmem.prepared.test()) Sleep(SLEEP_TIME_TINY);

		Lock(GLOCK1);
		double last = tmem.last_access;
		tmem.nref++; 
		AtomicMax(tmem.last_access, floor(GetElapse(GWAS_begin) * 1e3) + rnd);
		list.erase(last);
		list[tmem.last_access] = &tmem;
		UnLock(GLOCK1);

		return tmem;
	}
	else
	{
		//otherwise create a new entry
		Lock(GLOCK1);
		auto [it, inserted] = tab.try_emplace(hash);
		UnLock(GLOCK1);

		GWAS_MEM<REAL>& tmem = it->second;
		if (inserted)
		{
			AtomicMax(tmem.last_access, floor(GetElapse(GWAS_begin) * 1e3) + rnd);
			list[tmem.last_access] = &tmem;
			tmem.Calc(valid, hash, rnd);

			Lock(GLOCK1);
			double last = tmem.last_access;
			tmem.nref++;
			AtomicMax(tmem.last_access, floor(GetElapse(GWAS_begin) * 1e3) + rnd);
			list.erase(last);
			list[tmem.last_access] = &tmem;
			if (tab.size() >= 100)
			{
				for (int i = 0; i < 10; i++)
				{
					for (auto& ii : list)
					{
						GWAS_MEM<REAL>& tmem2 = *ii.second;
						atomic_flag del;
						del.test_and_set();

						tab.erase(tmem2.hash);
						list.erase(ii.first);
						tmem2.~GWAS_MEM();
						break;
					}
				}
			}
			UnLock(GLOCK1);
		}
		else
		{
			while (!tmem.prepared.test()) Sleep(SLEEP_TIME_TINY);

			Lock(GLOCK1);
			double last = tmem.last_access;
			tmem.nref++;
			AtomicMax(tmem.last_access, floor(GetElapse(GWAS_begin) * 1e3) + rnd);
			list.erase(last);
			list[tmem.last_access] = &tmem;
			UnLock(GLOCK1);
		}

		return tmem;
	}
}

/* Do nothing */
template<typename REAL>
TARGET GWAS_MEM<REAL>::GWAS_MEM()
{

}

/* Allocate memory and perform Eigen-value decomposition */
template<typename REAL>
TARGET void GWAS_MEM<REAL>::Calc(ucol& valid, HASH _hash, double rnd)
{
	prepared.clear();
	int64 n = valid.n_elem;

	hash = _hash;
	buf = new REAL[n + n * n];
	nref = 0;
	last_access = 0;
	prepared.clear();
	del.clear();

	rmat mR = GWAS_Root<REAL>.oR.submat(valid, valid);
	rcol cV(buf    , n, false, true);
	rmat mU(buf + n, n, n, false, true);
	
	if (gwas_imputeG_val != 10) openblas_set_num_threads(g_nthread_val);
	Evd(mR, mU, cV);
	if (gwas_imputeG_val != 10) openblas_set_num_threads(1);

	cV.elem(find(cV < (REAL)0.0001)).fill((REAL)0.0001);

	AtomicMax(last_access, floor(GetElapse(GWAS_begin) * 1e3) + rnd);

	prepared.test_and_set();
}

/* Destructor */
template<typename REAL>
TARGET GWAS_MEM<REAL>::~GWAS_MEM()
{
	if (nref == 0 && del.test())
	{
		hash = 0;
		if (buf) DEL(buf);
		nref = 0;
		last_access = 0;
		prepared.clear();
		del.clear();
	}
}
#endif

#define extern 
template<typename REAL>
extern REAL GWAS_NAN;
template<typename REAL>
extern GWAS<REAL> GWAS_Root;
template<typename REAL>
extern GWAS<REAL>* GWAS_Threads;
template<typename REAL>
extern umap<HASH, GWAS_MEM<REAL>> GWAS_umap;
template<typename REAL>
extern map<double, GWAS_MEM<REAL>*> GWAS_map;

extern timepoint GWAS_begin;
extern vector<string> GWAS_colname;
extern vector<string> GWAS_coltype;
extern umap<string, int> GWAS_indid;
extern atomic<int64> GWAS_batch_index;
extern atomic<int> GWAS_ploidy_presence[N_MAX_PLOIDY + 1];

#undef extern 

/* Is a NA or NAN string */
TARGET bool IsNaNStr(char* ptr)
{
	if (LwrParCmp(ptr, "nan") == 0 || LwrParCmp(ptr, "na") == 0)
		return true;
	return false;
}

/* Read a line from csv file */
TARGET bool ReadCsvLine(vector<string>& row, ifstream& file)
{
	row.clear();
	string line;
	if (!getline(file, line)) return false;

    line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
    
	stringstream ss(line);
	string cell;

	while (getline(ss, cell, ','))
		row.push_back(cell);

	return true;
}

/* Remove duplicated columns in a design matrix */
template<typename REAL>
TARGET void RemoveDupCol(rmat &mX)
{
	int r = rank(mX);
	if (r < mX.n_cols)
	{
		rmat mQ, mR;
		qr(mQ, mR, mX);
		mX = mX.cols(find(abs(mR.diag()) > (REAL)1e-10));
	}
}

/* mean vector discarding nan values */
template<typename REAL>
TARGET rrow NanMean(rmat& x)
{
	uint64 n = x.n_rows;
	urow ntype = sum(x != nan);
	rrow ex = (sum(x) - (n - ntype) * nan) / ntype;

	ex(find(ntype < 1)).fill(0);
	return ex;
}

/* var vector discarding nan values */
template<typename REAL>
TARGET rrow NanVar(rmat& x, bool ispop)
{
	urow ntype = sum(x != nan);
	urow nmiss = sum(x == nan);
	rrow ex = (sum(x) - nmiss * nan) / ntype;
	rrow ex2 = (sum(x % x) - nmiss * nan * nan) / ntype;
	ex2 = ex2 - ex % ex;
	
	if (!ispop)
	{
		//sample variance
		ex2(ntype <= 1).fill(1);
	}
	else
	{
		//population variance
		ex2 = ex2 % (ntype / (ntype - 1));
		ex2(ntype < 1).fill(1);
	}

	return ex2;
}

/* covariance matrix discarding nan values */
template<typename REAL>
TARGET rmat NanCov(rmat& x, bool ispop)
{
	int64 m = x.n_cols;
	rmat type = conv_to<rmat>::from(x != nan);
	rmat ntype2 = type.t() * type;
	rmat x2  = x % type;
	rmat exy = (x2.t() * x2) / ntype2;
	rmat ex  = (x2.t() * type) / ntype2;
	rmat C   = exy - ex.t() % ex;

	if (!ispop)
	{
		//sample variance
		C = (ntype2 / (ntype2 - 1)) % C;
		C(find(ntype2 <= 1)).fill(0);
		for (int64 i = 0; i < m; ++i)
			if (ntype2(i, i) <= 1)
				C(i, i) = 1;
	}
	else
	{
		//population variance
		C(find(ntype2 < 1)).fill(0);
		for (int64 i = 0; i < m; ++i)
			if (ntype2(i, i) < 1)
				C(i, i) = 1;
	}

	return C;
}

/* Generate multivariate normal distributed random vectors */
template<typename REAL>
TARGET rmat MVNRand(RNG<double>& rng, rmat& Mu, rmat& Sigma, int64 n)
{
	rmat re, Sigmasq;
	n = Mu.n_rows == 1 ? n : Mu.n_rows;

	rmat U; rcol V;
	eig_sym(V, U, Sigma, "dc");
	ucol vid = find(V > (REAL)1e-10);
	Sigmasq = U.cols(vid);
	Sigmasq.each_row() %= sqrt(V.rows(vid).t());
	Sigmasq = Sigmasq.t();

	rmat xn(n, Sigmasq.n_rows);
	rng.Normal(xn.memptr(), xn.n_elem);
	re = xn * Sigmasq;

	if (Mu.n_rows == 1)
		re.each_row() += Mu;
	else
		re += Mu;

	return re;
}

/* Fill missing data */
template<typename REAL>
TARGET void Imputation(int stage, int method, rmat& G, int64 lst, int64 led, int64* offset, bool shrinked, rmat& R)
{
#define GetCol(x) (offset == NULL ? x : offset[x] - (shrinked ? x : 0) - col0)

	if (method == 10 && stage == 0) method = 1;
	if (method == 10) return;

	// mean|median|sample|norm|mvn|cmean|cmvn|knn|svd|discard
	int64 n = G.n_rows, m = G.n_cols, col0 = (offset == NULL ? 0 : offset[lst]);

	// 1 mean
 	if (method == 1)
	{
		urow nType = sum(G != nan);
		rrow Gmean = NanMean(G);
		
		for (int64 a = 0; a < m; ++a)
		{
			if (nType(a) == (uint64)n) continue;
			G.col(a).replace(nan, Gmean(a));
		}
	}

	// 2 median
	if (method == 2)
	{
		urow nType = sum(G != nan);

		for (int64 a = 0; a < m; ++a)
		{
			if (nType(a) == (uint64)n) continue;
			rmat Gi(G.memptr() + a * n, n, 1, false, true);
			REAL med = nType(a) >= 1 ? median(Gi.elem(find(Gi != nan))) : 0;
			Gi.replace(nan, med);
		}
	}

	// 3 sample
	if (method == 3)
	{
		urow nType = sum(G != nan);
		rrow Gmean = NanMean(G);

		for (int64 l = lst; l < led; l++)
		{
			int64 colst = GetCol(l), coled = GetCol(l + 1);
			if (nType(colst) == (uint64)n) continue;

			int64 ncol = coled - colst;
			rmat Gi(G.memptr() + colst * n, n, ncol, false, true);
			rcol Gi0(G.memptr() + colst * n, n, false, true);

			RNG<double> rng(g_seed_val + l, RNG_SALT_GWAS);
			ucol misid = find(Gi0 == nan);
			rmat Gvalid = Gi.rows(find(Gi0 != nan));
			ucol sampid = misid;
			rng.Integer(sampid.memptr(), sampid.n_elem, 0LL, Gvalid.n_rows);
			Gi.rows(misid) = Gvalid.rows(sampid);
		}
	}

	// 4 norm
	if (method == 4)
	{
		urow nType = sum(G != nan);
		rrow Gmean = NanMean(G), Gstd = sqrt(NanVar(G));

		for (int64 l = lst; l < led; l++)
		{
			int64 colst = GetCol(l), coled = GetCol(l + 1);
			if (nType(colst) == (uint64)n) continue;

			int64 ncol = coled - colst;
			rmat Gi(G.memptr() + colst * n, n, ncol, false, true);
			rcol Gi0(G.memptr() + colst * n, n, false, true);
			rrow Gm(Gmean.memptr() + colst, ncol, false, true);
			rrow Gs(Gstd.memptr() + colst, ncol, false, true);

			RNG<double> rng(g_seed_val + l, RNG_SALT_GWAS);
			ucol misid = find(Gi0 == nan);
			rmat Gf = zeros<rmat>(misid.n_elem, ncol);
			rng.Normal(Gf.memptr(), Gf.n_elem);
			Gi.rows(misid) = (Gf.each_row() % Gs).each_row() + Gm;
		}
	}

	// 5 mvn
	if (method == 5)
	{
		urow nType = sum(G != nan);
		rrow Gmean = NanMean(G);

		for (int64 l = lst; l < led; l++)
		{
			int64 colst = GetCol(l), coled = GetCol(l + 1);
			if (nType(colst) == (uint64)n) continue;

			int64 ncol = coled - colst;
			rmat Gi(G.memptr() + colst * n, n, ncol, false, true);
			rcol Gi0(G.memptr() + colst * n, n, false, true);
			rrow Gm(Gmean.memptr() + colst, ncol, false, true);

			RNG<double> rng(g_seed_val + l, RNG_SALT_GWAS);
			ucol misid = find(Gi0 == nan);
			rmat Gcov = NanCov(Gi);
			Gi.rows(misid) = MVNRand(rng, Gm, Gcov, misid.n_elem);
		}
	}

	// 6 cmean
	// 7 cmvn
	if (method == 6 || method == 7)
	{
		urow nType = sum(G != nan);
		rrow tMu = NanMean(G);
		rmat tSigma = NanCov(G);

		ucol ColidO, Colid;
		
		if (offset == NULL)
			Colid = regspace<ucol>(0, m - 1);
		else
		{
			new(&ColidO) ucol((uint64*)offset + lst, led - lst, false, true);
			Colid = ColidO - col0;
			if (shrinked) Colid = Colid - regspace<ucol>(lst, led - 1);
		}

		//extract first column for each locus
		rmat Gfirst = G.cols(Colid);
		nType = nType.cols(Colid);

		//whether typed
		umat Type1 = Gfirst != nan;
		umat Miss1 = Gfirst == nan;
		imat Type((int64*)Type1.memptr(), Type1.n_rows, Type1.n_cols, false, true);
		imat Miss((int64*)Miss1.memptr(), Miss1.n_rows, Miss1.n_cols, false, true);

		imat LeftCol = Type.each_row() % regspace<irow>(0, led - 1 - lst) - Miss; //Miss as -1
		imat RightCol = Type.each_row() % regspace<irow>(0, led - 1 - lst) + led * Miss;//Miss as led

		for (int64 l = lst; l < led; l++)
		{
			if (nType(l - lst) == (uint64)n) continue;
			int64 colst = GetCol(l), coled = GetCol(l + 1), ncol = coled - colst;
			rmat Gi(G.memptr() + colst * n, n, ncol, false, true);
			rcol Gi0(G.memptr() + colst * n, n, false, true);
			ucol lrange = regspace<ucol>(colst, coled - 1);

			icol LeftId = {}, RightId = {};

			if (l == lst)
				LeftId = -ones<icol>(n, 1);
			else
				LeftId = max(LeftCol.head_cols(l - lst), 1);

			if (l == led - 1)
				RightId = -ones<icol>(n, 1);
			else
				RightId = min(RightCol.tail_cols(led - 1 - l), 1);

			icol RefId = LeftId * 0x100000000 + RightId;
			icol UniqId = unique(RefId);

			RNG<double> rng(g_seed_val + l, RNG_SALT_GWAS);

			for (int64 u : UniqId)
			{
				ucol inds = find(RefId == u && Gi0 == nan), leftrange = {}, rightrange = {};
				uint* uu = (uint*)&u, left = uu[1], right = uu[0];

				if (left != 0xFFFFFFFF)
					leftrange = regspace<ucol>(GetCol(left), GetCol(left + 1) - 2 + shrinked);
				if (right != 0xFFFFFFFF)
					rightrange = regspace<ucol>(GetCol(right), GetCol(right + 1) - 2 + shrinked);

				//conditional mean and covariance
				rmat Muc, Sigmac;
				if (left != 0xFFFFFFFF || right != 0xFFFFFFFF)
				{
					//has reference loci
					ucol rrange = join_vert(leftrange, rightrange);
					rrow Mu1 = tMu.cols(lrange);
					rrow Mu2 = tMu.cols(rrange);
					rmat ref = G(inds, rrange);
					rmat S11 = tSigma.submat(lrange, lrange);
					rmat S21 = tSigma.submat(rrange, lrange);
					rmat S22 = tSigma.submat(rrange, rrange);
					rmat iS22S21 = solve(S22, S21, solve_opts::force_sym);
					// solve is faster than pinv and solve_opts::fast

					Muc = (ref.each_row() - Mu2) * iS22S21;
					Muc.each_row() += Mu1;
					if (method == 7) 
						Sigmac = S11 - S21.t() * iS22S21;
				}
				else
				{
					//no reference loci
					Muc = tMu.cols(lrange);
					if (method == 7) 
						Sigmac = tSigma.submat(lrange, lrange);
				}

				//cmean
				if (method == 6)
					G(inds, lrange) = Muc;
				//cmvn
				else
					G(inds, lrange) = MVNRand(rng, Muc, Sigmac);
			}
		}
	}

	// 8 knn
	if (method == 8)
	{
		ucol iType = sum(G != nan, 1);
		urow nType = sum(G != nan);
		rrow Gmean = NanMean(G);
		rmat G_imp = G;

		for (int64 i = 0; i < n; ++i)
		{
			if (iType(i) == (uint64)m) continue;

			rcol Rcol = R.col(i);

			Rcol(i) = (REAL)-1e30; //exclude i itself
			ucol idx = sort_index(Rcol, "descend");

			for (int64 l = lst; l < led; l++) 
			{
				int64 colst = GetCol(l), coled = GetCol(l + 1);
				if (G(i, colst) != nan) continue;

				// this locus need imputation
				int ncol = coled - colst, n_neigh = 0;
				span l1span(colst, coled - 1);
				rrow sum_neigh = zeros<rrow>(1, ncol);

				for (auto j : idx)
				{
					if (G(j, colst) != nan)
					{
						sum_neigh += G(j, l1span);
						n_neigh++;
						if (n_neigh >= gwas_nneighbor_val) break;
					}
				}
				if (n_neigh == 0)
					G_imp(i, l1span) = 0;
				else
					G_imp(i, l1span) = sum_neigh / n_neigh;
			}
		}

		G = G_imp;
	}

	// 9 svd
	if (method == 9)
	{
		rmat Gi = G;
		urow nType = sum(G != nan);
		rmat Gmean = (sum(G) - nType * nan) / nType;

		for (int64 a = 0; a < m; ++a)
		{
			if (nType(a) == (uint64)n) continue;
			Gi.col(a).replace(nan, Gmean(a));
		}

		rmat U, V; rcol S;
		Svd(Gi, U, S, V);

		int nsvd = std::min(gwas_nsvd_val, (int)S.n_elem);
		Gi = U.cols(0, nsvd - 1) * diagmat(S.rows(0, nsvd - 1)) * V.cols(0, nsvd - 1).t();
		G = G % (G != nan) + Gi % (G == nan);
	}
}

/* Calculate GWAS */
template<typename REAL>
TARGET void CalcGWAS()
{
	if (!gwas) return;

	GWAS_begin = GetNow();

	EvaluationBegin();
	OpenResFile("-gwas", "Genome-Wide Association Studies");
	fprintf(FRES, "%s%s", g_linebreak_val, g_linebreak_val);
	OpenTempFiles(g_nthread_val, ".gwas");

	AssignPop<REAL>(gwas_pop_b, gwas_pop_val, "-gwas_pop");
	GWAS_Root<REAL>.Prepare();

	GWAS_Threads<REAL> = new GWAS<REAL>[g_nthread_val];
	for (int i = 0; i < g_nthread_val; ++i)
		new (&GWAS_Threads<REAL>[i]) GWAS<REAL>(&GWAS_Root<REAL>, i);

	GWAS_batch_index = 0;
	RunThreads(&GWASThread<REAL>, NULL, NULL, GWAS_Root<REAL>.m, GWAS_Root<REAL>.m,
		"\nCalculating Genome-Wide Association Studies:\n", 1, true);

	JoinTempFiles(g_nthread_val);
	CloseResFile();
	EvaluationEnd("GWAS");

	if (gwas_plot_val == 1)
		RunRscript("gwas_plot.R");

	for (auto& it : GWAS_map<REAL>)
	{
		GWAS_MEM<REAL>& tmem2 = *it.second;
		tmem2.del.test_and_set();
		tmem2.~GWAS_MEM();
	}

	umap<HASH, GWAS_MEM<REAL>>().swap(GWAS_umap<REAL>);
	map<double, GWAS_MEM<REAL>*>().swap(GWAS_map<REAL>);
	vector<string>().swap(GWAS_colname);
	vector<string>().swap(GWAS_coltype);
	umap<string, int>().swap(GWAS_indid);
}

/* Calculate Relatedness using multiple threads */
THREAD2(GWASRelatednessThread)
{
	//load ind
	int64 nc = 0;
	int n = cpop<REAL>->nind, estimator = gwas_restimator_val - 2;
	IND<REAL>** inds = cpop<REAL>->inds;
	RELATEDNESS<REAL> Rc;
	atomic<REAL>* Rmem = (atomic<REAL>*)GWAS_Root<REAL>.oR.memptr();
	int ploidy_presence2[N_MAX_PLOIDY + 1] = { 0 };
	int64 m = GWAS_Root<REAL>.m;

	if (relatedness_estimator_val[8] || relatedness_estimator_val[9])  Anderson2007_Coef = new double[m * 3];
	if (relatedness_estimator_val[11]) Huang2015_Coef = new double[m * 9];

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j <= n; ++j, ++nc)
		{
			if (nc % g_nthread_val == threadid)
			{
				REAL t = Rc.RelatednessEstimator(estimator, inds[i], inds[j]);
				Rmem[i * n + j] = Rmem[j * n + i] = t;

				PROGRESS_VALUE++;
			}
		}
	}

	if (relatedness_estimator_val[8] || relatedness_estimator_val[9])  DEL(Anderson2007_Coef);
	if (relatedness_estimator_val[11]) DEL(Huang2015_Coef);

	//calculate ploidy presence
	if (maxploidy == minploidy)
		ploidy_presence2[maxploidy] = 1;
	else
	{
		for (int64 l = threadid; l < m; l += g_nthread_val)
		{
			GENO_READER rt(cpop<REAL>->ind0id, l);
			GENOTYPE* gtab = GetLoc(l).GetGtab();
			umap<int, int> gtid;
			for (int i = 0; i < n; ++i)
				gtid[rt.Read()] = 1;

			for (auto kv : gtid)
			{
				GENOTYPE& gt = gtab[kv.first];
				ploidy_presence2[gt.Ploidy()] = 1;
			}
		}
	}

	for (int v = 1; v <= N_MAX_PLOIDY; ++v)
		GWAS_ploidy_presence[v] |= ploidy_presence2[v];

	PROGRESS_VALUE += 10;
}

/* Calculate Relatedness using multiple threads */
THREAD2(GWASMatrixRelatednessThread)
{
	int64 m = GWAS_Root<REAL>.m;
	int64 n = GWAS_Root<REAL>.n;
	int64* noffset = GWAS_Root<REAL>.noffset;
	int ploidy_presence2[N_MAX_PLOIDY + 1] = { 0 };
	rmat& mr = GWAS_Root<REAL>.tR[threadid];
	SetZero(mr.memptr(), n * n);

	for (int64 lst = GWAS_batch_index.fetch_add(gwas_batch_val); lst < m; lst = GWAS_batch_index.fetch_add(gwas_batch_val))
	{
		int64 col0 = noffset[lst], led = std::min(lst + gwas_batch_val, m);
		int64 tnallele = noffset[led] - noffset[lst];
		rmat tG(n, tnallele);
		rrow tNallele(tnallele);

		for (int64 l = lst; l < led; ++l)
		{
			int64 colst = noffset[l], coled = noffset[l + 1] - 1;
			GENO_READER rt(cpop<REAL>->ind0id, l);
			GENOTYPE* gtab = GetLoc(l).GetGtab();
			int nallele = coled - colst + 1;
			int pallele = colst - col0;
			tNallele.cols(pallele, pallele + nallele - 1).fill(nallele);

			for (int i = 0; i < n; ++i)
			{
				GENOTYPE& gt = gtab[rt.Read()];
				int vi = gt.Ploidy();
				REAL iv = 1.0 / vi;
				ploidy_presence2[vi] = 1;

				if (gt.Nalleles() == 0)
				{
					for (int a = 0; a < nallele; ++a)
						tG(i, pallele + a) = nan;
				}
				else
				{
					ushort* als = gt.GetAlleleArray();
					for (int a = 0; a < vi; ++a)
						tG(i, pallele + als[a]) += iv;
				}
			}

			PROGRESS_VALUE++;
		}

		// fill missing, gwas_release
		Imputation<REAL>(0, gwas_imputeG_val, tG, lst, led, noffset, false, *(rmat*)NULL);

		rrow tMean = mean(tG), tStd  = stddev(tG);
		tG.each_row() -= tMean;
		tG.each_row() /= sqrt(m * tNallele);
		if (gwas_restimator_val == 2)
			tG.each_row() /= tStd;
		rmat tG2 = tG.cols(tStd > (REAL)1e-9);
		mr += tG2 * tG2.t();

		PROGRESS_VALUE += 10;
	}

	for (int v = 1; v <= N_MAX_PLOIDY; ++v)
		GWAS_ploidy_presence[v] |= ploidy_presence2[v];

	PROGRESS_VALUE += 1;
}

/* Calculate GWAS using multiple threads */
THREAD2(GWASThread)
{
	int64 m = GWAS_Root<REAL>.m;
	GWAS<REAL>* threads = GWAS_Threads<REAL>;
	threads[0].fout = TEMP_FILES[0];
	threads[0].WriteGWASHeader();
	byte* state = new byte[gwas_batch_val];
	
	for (int64 lst = GWAS_batch_index.fetch_add(gwas_batch_val); lst < m; lst = GWAS_batch_index.fetch_add(gwas_batch_val))
	{
		int64 led = std::min(lst + gwas_batch_val, m);

		threads[0].ReadBatch(lst, led);

		for (int i = 1; i < g_nthread_val; ++i)
			threads[i].CopyRef(&threads[0]);

		int64 lwrite = lst;
		SetZero(state, sizeof(byte) * gwas_batch_val);

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int64 l = lst; l < led; ++l)
		{
			int tid = omp_get_thread_num();
			threads[tid].CalcLocus(l);
			state[l - lst] = 1;

			if (tid == 0)
			{
				for (; lwrite < led && state[lwrite - lst]; ++lwrite)
					threads[0].WriteLocus(lwrite);
			}

			PROGRESS_VALUE++;
		}

		for (; lwrite < led; ++lwrite)
			threads[0].WriteLocus(lwrite);
	}
}
