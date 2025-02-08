/* Genetic Differentiation Functions */

#include "vcfpop.h"

template struct FST<double>;
template struct FST<float >;

template TARGET void CalcDiff<double>();
template TARGET void CalcDiff<float >();
template TARGET double FST<double>::dxy(POP<double>** grps, int n, double* buf, int64 l);
template TARGET double FST<float >::dxy(POP<float >** grps, int n, double* buf, int64 l);

template<> FST<double>* fst_buf<double>[6];
template<> FST<float >* fst_buf<float >[6];

#ifndef _FST
/* Fst estimator warpper */
template<typename REAL>
TARGET double FST<REAL>::FstEstimator(POP<REAL>** grps, int n, int e, double* each, double* buf)
{
	switch (e)
	{
	case 1: return Fst_Nei1973(grps, n, each, buf);
	case 2: return Fst_Weir1984(grps, n, each);
	case 3: return Fst_Hudson1992(grps, n, each, buf);
	case 4: return Fst_Slatkin1995(grps, n, each, buf);
	case 5: return Fst_Hedrick2005(grps, n, each, buf);
	case 6: return Fst_Jost2008(grps, n, each, buf);
	case 7: return Fst_Huang2021_homo(grps, n, 2, true, each);//homoploid
	case 8: return Fst_Huang2021_aneu(grps, n, 2, true, true, each, buf);//aneuploid
	}
	return 0;
}

/* Estimate Fst and test differentiation for two pops */
template<typename REAL>
TARGET void FST<REAL>::CalcFst(POP<REAL>* a, POP<REAL>* b)
{
	POP<REAL>* c[] = { a, b };
	CalcFst(c, 2);
}

/* Estimate Fst and test differentiation for multiple pops */
template<typename REAL>
TARGET void FST<REAL>::CalcFst(POP<REAL>** grps, int n)
{
	Genotype_G = Genotype_P = Allele_G = Allele_P = 0; Genotype_DF = Allele_DF = 0;//Difference
	Nei1973 = Weir1984 = Hudson1992 = Slatkin1995 = Hedrick2005 = Jost2008 = Huang2021_homo = Huang2021_aneu = 0;

	VLA_NEW(buf, double, maxK);
	byte* estimator = fst_estimator_val;
	SetZero(&Nei1973, N_FST_ESTIMATOR);

	for (int i = 1; i <= N_FST_ESTIMATOR; ++i)
	{
		if (estimator[i] == 0) continue;
		*(&Nei1973 + i - 1) = fst_locus_val[2] ? new double[nloc] : NULL;
		*(&Nei1973T + i - 1) = FstEstimator(grps, n, i, *(&Nei1973 + i - 1), buf);
	}

	//test Fst
	if (fst_locus_val[1] == 0 && fst_locus_val[2] == 0)
		Genotype_PT = Allele_PT = NAN;

	if ((fst_locus_val[1] || fst_locus_val[2]) && fst_test_val[1])
	{
		Genotype_DFT = 0; Genotype_GT = 0;
		Genotype_G = new double[nloc];
		Genotype_DF = new int[nloc];
		Genotype_P = new double[nloc];

		int64 gcsize = maxG + 2;
		VLA_NEW(obs, double, gcsize * (4 * n + 1));
		double* obs2 = obs + gcsize * n;
		double* exp = obs2 + gcsize * n;
		double* exp2 = exp + gcsize * n;
		double* colsum = exp2 + gcsize * n;
		VLA_NEW(rowsum, double, n);

		for (int64 l = 0; l < nloc; ++l)
		{
			bool sameploidy = true;
			int ploidy = -1;
			GENOTYPE* gtab = GetLoc(l).GetGtab();
			int ngeno = GetLoc(l).ngeno;

			SetZero(obs, ngeno * n);
			if (ngeno >= 2 && n >= 2)
			{
				for (int gi = 0; gi < ngeno; ++gi)
				{
					GENOTYPE& gt = gtab[gi];
					if (gt.Nalleles() == 0) continue;

					int v = gt.Ploidy();
					if (ploidy == -1) ploidy = v;
					if (ploidy != v)
					{
						sameploidy = false;
						break;
					}

					for (int p = 0; p < n; ++p)
						obs[p * ngeno + gi] = grps[p]->GetGenoCount(l, gi);
				}
				if (sameploidy)
				{
					CombineTable(obs, n, ngeno, Genotype_G[l], Genotype_DF[l], Genotype_P[l], (bool)fst_locus_val[2], obs2, exp, exp2, rowsum, colsum);
					Genotype_DFT += Genotype_DF[l];
					Genotype_GT += Genotype_G[l];
				}
			}
			if (ngeno < 2 || n < 2 || !sameploidy)
			{
				Genotype_DF[l] = 0;
				Genotype_G[l] = 0;
				Genotype_P[l] = NAN;
			}
		}
		VLA_DELETE(obs);
		VLA_DELETE(rowsum);
		Genotype_PT = Genotype_DFT > 0 ? 1 - ChiSquareDistCDF(Genotype_GT, Genotype_DFT) : NAN;
	}

	if ((fst_locus_val[1] || fst_locus_val[2]) && fst_test_val[2])
	{
		Allele_DFT = 0; Allele_GT = 0;
		Allele_G = new double[nloc];
		Allele_DF = new int[nloc];
		Allele_P = new double[nloc];
		VLA_NEW(obs, double, maxK * n * 4 + maxK + n);
		double* obs2 = obs + maxK * n;
		double* exp = obs2 + maxK * n;
		double* exp2 = exp + maxK * n;
		double* colsum = exp2 + maxK * n;
		double* rowsum = colsum + maxK;

		for (int64 l = 0; l < nloc; ++l)
		{
			int m = grps[0]->loc_stat1[l].k;

			if (m >= 2 && n >= 2)
			{
				SetZero(obs, m * n);
				for (int p = 0; p < n; ++p)
				{
					REAL* freq = grps[p]->GetFreq(l);
					int nhaplo = grps[p]->loc_stat1[l].nhaplo;
					for (int i = 0; i < m; ++i)
						obs[p * m + i] = nhaplo * freq[i];
				}
				CombineTable(obs, n, m, Allele_G[l], Allele_DF[l], Allele_P[l], (bool)fst_locus_val[2], obs2, exp, exp2, rowsum, colsum);
				Allele_DFT += Allele_DF[l];
				Allele_GT += Allele_G[l];
			}
			else
			{
				Allele_DF[l] = 0; Allele_G[l] = 0;
				Allele_P[l] = NAN;
			}
		}
		Allele_PT = Allele_DFT > 0 ? 1 - ChiSquareDistCDF(Allele_GT, Allele_DFT) : NAN;
		VLA_DELETE(obs);
	}

	PROGRESS_VALUE++;

	VLA_DELETE(buf);
	return;
}

/* Uninitialize */
template<typename REAL>
TARGET void FST<REAL>::Uninitialize()
{
	DEL(Nei1973);
	DEL(Weir1984);
	DEL(Hudson1992);
	DEL(Slatkin1995);
	DEL(Hedrick2005);
	DEL(Jost2008);
	DEL(Huang2021_homo);
	DEL(Huang2021_aneu);

	DEL(Genotype_G);
	DEL(Genotype_DF);
	DEL(Genotype_P);
	DEL(Allele_G);
	DEL(Allele_DF);
	DEL(Allele_P);
}

/* Nei 1973 Gst estimator based on heterozgysotiy */
template<typename REAL>
TARGET double FST<REAL>::Fst_Nei1973(POP<REAL>** grps, int n, double* each, double* buf, int64 _l, double* frac1, double* frac2)
{
	//weight by pop size chakraborty 1974, allow ad
	double Dst = 0, Ht = 0;
	double* pi = buf;
	for (int64 l = (_l == -1 ? 0 : _l); l < (_l == -1 ? nloc : _l + 1); ++l)
	{
		if (each) each[l] = NAN;
		int k2 = GetLoc(l).k;
		if (k2 < 2) continue;

		double ht = 0, hs = 0, sw = 0;
		int s = 0;
		SetZero(pi, k2);
		for (int i = 0; i < n; ++i)
		{
			int nhaplo = grps[i]->loc_stat1[l].nhaplo;
			REAL* p = grps[i]->GetFreq(l);
			if (nhaplo == 0) continue;
			s++;
			AddProd(pi, p, nhaplo, k2);
			sw += nhaplo;
			hs += nhaplo * SumSquare(p, k2);
		}

		if (s < 2) continue;
		Unify(pi, k2);
		hs = 1 - hs / sw;
		ht = 1 - SumSquare(pi, k2);
		Dst += ht - hs;
		Ht += ht;
		if (each) each[l] = IsError((ht - hs) / ht) ? NAN : (ht - hs) / ht;
	}
	if (_l != -1) { *frac1 = Dst; *frac2 = Ht; }
	return Dst / Ht;
}

/* Weir 1984 Fst estimator based on variance components */
template<typename REAL>
TARGET double FST<REAL>::Fst_Weir1984(POP<REAL>** grps, int n, double* each, int64 _l, double* frac1, double* frac2)
{
	//for diploid only, allow ad
	if (minploidy != 2 || maxploidy != 2)
		Exit("\nError: Weir1984 Fst estimator can only be applied for diploids\n");

	double theta_W1 = 0, theta_W2 = 0;
	for (int64 l = (_l == -1 ? 0 : _l); l < (_l == -1 ? nloc : _l + 1); ++l)
	{
		if (each) each[l] = NAN;
		int k2 = GetLoc(l).k;
		if (k2 < 2) continue;

		double nt = 0, nb = 0, nc = 0;
		double r = 0, ni2 = 0;
		for (int i = 0; i < n; ++i)
		{
			int nhaplo = grps[i]->loc_stat1[l].nhaplo;
			if (nhaplo == 0) continue;
			r++;
			int ni = nhaplo >> 1;
			nt += ni;
			ni2 += ni * ni;
		}
		if (r < 2 || nt == 0) continue;

		nb = nt / r;
		nc = (r * nb - ni2 / nt) / (r - 1);

		double w1 = 0, w2 = 0;
		for (int a = 0; a < k2; ++a)
		{
			double p1 = 0, p2 = 0, s2 = 0, hb = 0;
			for (int i = 0; i < n; ++i)
			{
				int nhaplo = grps[i]->loc_stat1[l].nhaplo;
				if (nhaplo == 0) continue;
				double p = grps[i]->GetFreq(l, a);
				int ni = nhaplo >> 1;
				p1 += p * ni;
				p2 += p * p * ni;
			}
			p1 /= nt;
			p2 /= nt;
			s2 = (p2 - p1 * p1) * r / (r - 1);
			hb = 2 * (p1 - p2);

			if (nb <= 1 || hb <= 0) break;

			double va = nb / nc * (s2 - 1 / (nb - 1) * (p1 * (1 - p1) - (r - 1) / r * s2 - 0.25 * hb));
			double vb = nb / (nb - 1) * (p1 * (1 - p1) - (r - 1) / r * s2 - (2 * nb - 1) / (4 * nb) * hb);
			double vc = 0.5 * hb;

			if (IsError(va) || IsError(vb) || IsError(vc)) break;

			w1 += va;
			w2 += va + vb + vc;
		}
		theta_W1 += w1;
		theta_W2 += w2;
		if (each) each[l] = IsError(w1 / w2) ? NAN : w1 / w2;
	}

	if (_l != -1) { *frac1 = theta_W1; *frac2 = theta_W2; }
	return theta_W1 / theta_W2;
}

/* Hudson 1992 Fst estimator based on mean allele difference */
template<typename REAL>
TARGET double FST<REAL>::Fst_Hudson1992(POP<REAL>** grps, int n, double* each, double* buf, int64 _l, double* frac1, double* frac2)
{
	//do not need to weight, allow ad
	double* Ni = buf;
	double f1 = 0, f2 = 0;
	for (int64 l = (_l == -1 ? 0 : _l); l < (_l == -1 ? nloc : _l + 1); ++l)
	{
		if (each) each[l] = NAN;
		int k2 = GetLoc(l).k;
		if (k2 < 2) continue;

		int Nh = 0, Np = 0;
		double Dw = 0, Nw = 0;
		SetZero(Ni, k2);
		for (int i = 0; i < n; ++i)
		{
			int nhaplo = grps[i]->loc_stat1[l].nhaplo;
			if (nhaplo == 0) continue;
			REAL* p = grps[i]->GetFreq(l);
			AddProd(Ni, p, nhaplo, k2);

			Nh += nhaplo;
			Np++;

			int nh = nhaplo;
			Dw += nh * (nh - 1); //2x
			Nw += nh * (nh - 1); //2x
			for (int a = 0; a < k2; ++a)
			{
				double na = nh * p[a];
				if (na > 1) Dw -= na * (na - 1); //substrate same
			}
		}

		if (Np < 2) continue;

		double hw = Dw / Nw; //mean diff
		double Db = (double)(Nh * (Nh - 1)), Nb = Db;
		for (int a = 0; a < k2; ++a)
			Db -= Ni[a] * (Ni[a] - 1);
		double hb = (Db - Dw) / (Nb - Nw);

		if (each) each[l] = IsError(1 - hw / hb) ? NAN : 1 - hw / hb;
		if (IsNormal(hw) && IsNormal(hb))
		{
			f1 += hb - hw;
			f2 += hb;
		}
	}
	if (_l != -1) { *frac1 = f1; *frac2 = f2; }
	return f1 / f2;
}

/* Slatkin 1995 Fst estimator based on allele size */
template<typename REAL>
TARGET double FST<REAL>::Fst_Slatkin1995(POP<REAL>** grps, int n, double* each, double* buf, int64 _l, double* frac1, double* frac2)
{
	if (abs(g_format_val) <= BCF)
		Exit("\nError: Slatkin1995 Fst estimator can only be applied for non-vcf input file, and should use size as allele identifier. \n");

	//do not need to weight, allow ad
	double* nt = buf;
	double f1 = 0, f2 = 0;
	for (int64 l = (_l == -1 ? 0 : _l); l < (_l == -1 ? nloc : _l + 1); ++l)
	{
		double Sw1 = 0, Sw2 = 0, St1 = 0, St2 = 0;
		if (each) each[l] = NAN;
		int k2 = GetLoc(l).k;
		if (k2 < 2) continue;

		int Nh = 0, Np = 0;
		SetZero(nt, k2);
		ushort* alen = GetLoc(l).GetAlenArray();

		for (int i = 0; i < n; ++i)
		{
			int nhaplo = grps[i]->loc_stat1[l].nhaplo;
			if (nhaplo == 0) continue;
			REAL* p = grps[i]->GetFreq(l);
			AddProd(nt, p, nhaplo, k2);

			Nh += nhaplo;
			Np++;
			Sw2 += nhaplo * (nhaplo - 1) * 0.5;

			for (int ai = 0; ai < k2; ++ai)
			{
				double ni = nhaplo * p[ai];
				for (int aj = 0; aj < ai; ++aj)
				{
					double nj = nhaplo * p[aj];
					Sw1 += ni * nj * ((int)alen[ai] - (int)alen[aj]) * ((int)alen[ai] - (int)alen[aj]);
				}
			}
		}

		if (Np < 2) continue;
		St2 = Nh * (Nh - 1) * 0.5;
		for (int i = 0; i < k2; ++i)
			for (int j = 0; j < i; ++j)
				St1 += nt[i] * nt[j] * ((int)alen[i] - (int)alen[j]) * ((int)alen[i] - (int)alen[j]);

		double Sw = Sw1 / Sw2, St = St1 / St2;
		if (IsNormal(Sw) && IsNormal(St))
		{
			if (each) each[l] = IsError((St - Sw) / St) ? NAN : (St - Sw) / St;
			f1 += St - Sw;
			f2 += St;
		}
	}

	if (_l != -1) { *frac1 = f1; *frac2 = f2; }
	return f1 / f2;
}

/* Hedrick 2005 G'st */
template<typename REAL>
TARGET double FST<REAL>::Fst_Hedrick2005(POP<REAL>** grps, int n, double* each, double* buf, int64 _l, double* frac1, double* frac2)
{
	//weight by pop size chakraborty 1974, allow ad
	double DA = 0, DB = 0;
	double* pi = buf;
	for (int64 l = (_l == -1 ? 0 : _l); l < (_l == -1 ? nloc : _l + 1); ++l)
	{
		if (each) each[l] = NAN;
		int k2 = GetLoc(l).k;
		if (k2 < 2) continue;

		double ht = 0, hs = 0, sw = 0;
		int s = 0;
		SetZero(pi, k2);
		for (int i = 0; i < n; ++i)
		{
			int nhaplo = grps[i]->loc_stat1[l].nhaplo;
			if (nhaplo == 0) continue;
			REAL* p = grps[i]->GetFreq(l);
			s++;
			sw += nhaplo;
			AddProd(pi, p, nhaplo, k2);
			hs += nhaplo * SumSquare(p, k2);
		}
		if (s < 2) continue;
		Unify(pi, k2);
		hs = 1 - hs / sw;
		ht = 1 - SumSquare(pi, k2);
		double da = (ht - hs) * (s - 1 + hs);
		double db = ht * (s - 1) * (1 - hs);
		DA += da;
		DB += db;
		if (each) each[l] = IsError(da / db) ? NAN : da / db;
	}

	if (_l != -1) { *frac1 = DA; *frac2 = DB; }
	return DA / DB;
}

/* Jost 2008 D */
template<typename REAL>
TARGET double FST<REAL>::Fst_Jost2008(POP<REAL>** grps, int n, double* each, double* buf, int64 _l, double* frac1, double* frac2)
{
	//weight by pop size chakraborty 1974, allow ad
	double DA = 0, DB = 0;
	double* pi = buf;
	for (int64 l = (_l == -1 ? 0 : _l); l < (_l == -1 ? nloc : _l + 1); ++l)
	{
		if (each) each[l] = NAN;
		int k2 = GetLoc(l).k;
		if (k2 < 2) continue;

		double ht = 0, hs = 0, sw = 0;
		int s = 0;
		SetZero(pi, k2);
		for (int i = 0; i < n; ++i)
		{
			int nhaplo = grps[i]->loc_stat1[l].nhaplo;
			if (nhaplo == 0) continue;
			REAL* p = grps[i]->GetFreq(l);
			s++;
			sw += nhaplo;
			AddProd(pi, p, nhaplo, k2);
			hs += nhaplo * SumSquare(p, k2);
		}
		if (s < 2) continue;
		Unify(pi, k2);
		hs = 1 - hs / sw;
		ht = 1 - SumSquare(pi, k2);
		double da = (ht - hs) * s;
		double db = ht * (s - 1);
		DA += da;
		DB += db;
		if (each) each[l] = IsError(da / db) ? NAN : da / db;
	}

	if (_l != -1) { *frac1 = DA; *frac2 = DB; }
	return DA / DB;
}

/* Huang 2021 Fst estimator based on multi-level amova */
template<typename REAL>
TARGET double FST<REAL>::Fst_Huang2021_homo(POP<REAL>** grps, int n, int layer, bool isiam, double* each)
{
	//2 or 3 Level, homoploids, forbid ad
	if (ad && layer == 3) Exit("\nError: Huang2021 Fst estimator (-fst_estimator=Huang2021_homo) is incompatible with allelic depth (-ad) option.\n");
	double svi2dvp = 0, svi2dvt = 0, svp2dvt = 0;
	int Ni = 0, Np = n, Nh = 0;
	VLA_NEW(pop_nhaplo, int, n);
	SetZero(pop_nhaplo, n);

	for (int p = 0; p < n; ++p)
	{
		Ni += grps[p]->nind;
		IND<REAL>** vind = grps[p]->inds;
		for (int i = 0; i < grps[p]->nind; ++i)
		{
			int v1 = vind[i]->vmin, v2 = vind[i]->vmax;
			if (v1 != v2)
				Exit("\nError: Huang2021 Fst estimator (-fst_estimator=Huang2021_homo) do not support aneuploids, in individual %s.\n", vind[i]->name);
			pop_nhaplo[p] += v1;
			Nh += v1;
		}
	}

	struct HAP
	{
		double invi;
		double invp;
		int indid;
		ushort popid;
	};

	VLA_NEW(hap, HAP, Nh);
	double invt = 1.0 / Nh;
	ushort* hap_bucket = new ushort[Nh * nloc];
	SetFF(hap_bucket, Nh * nloc);

	//place alleles
	for (int p = 0, ph = 0; p < n; ++p)
	{
		IND<REAL>** vind = grps[p]->inds;
		int iend = grps[p]->nind;

		for (int64 l = 0; l < nloc; ++l)
		{
			GENO_READER rt(grps[p]->ind0id, l);
			GENOTYPE* gtab = GetLoc(l).GetGtab();

			for (int i = 0, ph2 = ph; i < iend; ++i)
			{
				ushort* als = gtab[rt.Read()].GetAlleleArray();
				for (int j = 0, vi = vind[i]->vmin; j < vi; ++j, ++ph2)
					hap_bucket[ph2 * nloc + l] = als[j];
			}
		}


		double invp = 1.0 / pop_nhaplo[p];
		for (int i = 0; i < iend; ++i)
		{
			int vi = vind[i]->vmin;

			hap[ph].indid = i;
			hap[ph].popid = (ushort)p;
			hap[ph].invi = 1.0 / vi;
			hap[ph].invp = invp;

			for (int j = 1; j < vi; ++j)
				hap[ph + j] = hap[ph];

			ph += vi;
		}
	}

	VLA_DELETE(pop_nhaplo);

	//calc ns
	for (int p = 0; p < n; ++p)
	{
		int vp = 0;
		double psvi2 = 0;
		IND<REAL>** vind2 = grps[p]->inds;
		for (int i = 0, iend = grps[p]->nind; i < iend; ++i)
		{
			int vi = vind2[i]->vmin;
			psvi2 += vi * vi;
			vp += vi;
		}
		svi2dvt += psvi2;
		svi2dvp += psvi2 / vp;
		svp2dvt += vp * vp;
	}
	svi2dvt /= Nh;
	svp2dvt /= Nh;

	double n11 = 0;
	double n21 = 0, n22 = 0;
	double n31 = 0, n32 = 0, n33 = 0;

	switch (layer)
	{
	case 2:
		n11 = Nh - Np;
		n21 = Nh - svp2dvt;
		n22 = Nh - 1;
		break;
	case 3:
		n11 = Nh - Ni;
		n21 = Nh - svi2dvp; n22 = Nh - Np;
		n31 = Nh - svp2dvt; n32 = Nh - svi2dvt; n33 = Nh - 1;
		break;
	}

	VLA_NEW(lsswi, double, each ? nloc : 0);
	VLA_NEW(lsswp, double, each ? nloc : 0);
	VLA_NEW(lsstot, double, each ? nloc : 0);
	double sswi = 0, sswp = 0, sstot = 0;
	double vwi = 0, vai = 0, vwp = 0, vap = 0, vtot = 0;

	for (int i = 0; i < Nh; ++i)
	{
		HAP& hi = hap[i];
		ushort* ai = hap_bucket + i * nloc;
		for (int j = 0; j < i; ++j)
		{
			HAP& hj = hap[j];
			ushort* aj = hap_bucket + j * nloc;
			double gd = 0;
			for (int64 l = 0; l < nloc; ++l)
			{
				SLOCUS loc = GetLoc(l);
				int k2 = loc.k;

				double tgd = 0;
				if (isiam)
				{
					if (ai[l] == 0xFFFF && aj[l] == 0xFFFF)
						tgd = 1 - SumProd(grps[hap[i].popid]->GetFreq(l), grps[hap[j].popid]->GetFreq(l), k2);
					else if (ai[l] == 0xFFFF)
						tgd = 1 - grps[hap[i].popid]->GetFreq(l, aj[l]);
					else if (aj[l] == 0xFFFF)
						tgd = 1 - grps[hap[j].popid]->GetFreq(l, ai[l]);
					else if (ai[l] != aj[l])
						tgd = 1;
				}
				else
				{
					if (ai[l] == 0xFFFF && aj[l] == 0xFFFF)
						tgd = SumProdSMM(loc.GetAlenArray(), grps[hap[i].popid]->GetFreq(l), grps[hap[j].popid]->GetFreq(l), k2);
					else if (ai[l] == 0xFFFF)
						tgd = SumProdSMM(loc.GetAlenArray(), grps[hap[i].popid]->GetFreq(l), aj[l], k2);
					else if (aj[l] == 0xFFFF)
						tgd = SumProdSMM(loc.GetAlenArray(), grps[hap[j].popid]->GetFreq(l), ai[l], k2);
					else if (ai[l] != aj[l])
						tgd = loc.GetSMMDist(ai[l], aj[l]);
				}

				if (each)
				{
					if (hi.indid == hj.indid) lsswi[l] += tgd * hi.invi;
					if (hi.popid == hj.popid) lsswp[l] += tgd * hi.invp;
					lsstot[l] += tgd * invt;
				}
				gd += tgd;
			}
			if (hi.indid == hj.indid) sswi += gd * hi.invi;
			if (hi.popid == hj.popid) sswp += gd * hi.invp;
			sstot += gd * invt;
		}
	}
	DEL(hap_bucket);
	VLA_DELETE(hap);

	if (each) for (int64 l = 0; l < nloc; ++l)
	{
		switch (layer)
		{
		case 2:
			vwp = lsswp[l] / n11;
			vap = (n11 * lsstot[l] - n22 * lsswp[l]) / (n11 * n21);
			vtot = vwp + vap;
			break;
		case 3:
			vwi = lsswi[l] / n11;
			vai = (n11 * lsswp[l] - n22 * lsswi[l]) / (n11 * n21);
			vap = (n22 * n32 * lsswi[l] - n21 * n33 * lsswi[l] - n11 * n32 * lsswp[l] + n11 * n21 * lsstot[l]) / (n11 * n21 * n31);
			vtot = vwi + vai + vap;
			break;
		}
		each[l] = vap / vtot;
	}
	VLA_DELETE(lsswi);
	VLA_DELETE(lsswp);
	VLA_DELETE(lsstot);

	switch (layer)
	{
	case 2:
		vwp = sswp / n11;
		vap = (n11 * sstot - n22 * sswp) / (n11 * n21);
		vtot = vwp + vap;
		break;
	case 3:
		vwi = sswi / n11;
		vai = (n11 * sswp - n22 * sswi) / (n11 * n21);
		vap = (n22 * n32 * sswi - n21 * n33 * sswi - n11 * n32 * sswp + n11 * n21 * sstot) / (n11 * n21 * n31);
		vtot = vwi + vai + vap;
		break;
	}
	return vap / vtot;
}

/* Huang 2021 Fst estimator based on multi-level amova */
template<typename REAL>
TARGET double FST<REAL>::Fst_Huang2021_aneu(POP<REAL>** grps, int n, int layer, bool isiam, bool sumss, double* each, double* buf, int64 _l, double* frac1, double* frac2)
{
	//2 or 3 Level, ansioploids
	if (ad && layer == 3) Exit("\nError: Huang2021 Fst estimator (-fst_estimator=Huang2021_aneu) is incompatible with allelic depth (-ad) option.\n");
	double Vap = 0, Vtot = 0;
	double* Na = buf;
	double SStot = 0, SSwp = 0, SSwi = 0;
	double N11 = 0, N21 = 0, N22 = 0, N31 = 0, N32 = 0, N33 = 0;
	for (int64 l = (_l == -1 ? 0 : _l); l < (_l == -1 ? nloc : _l + 1); ++l)
	{
		if (each) each[l] = NAN;
		int k2 = GetLoc(l).k;
		ushort* alen = GetLoc(l).GetAlenArray();

		if (k2 < 2) continue;

		//assign allele
		int Nh = 0, Np = 0, Ni = 0;
		double n11 = 0, n21 = 0, n22 = 0, n31 = 0, n32 = 0, n33 = 0;
		double svi2dvp = 0, svi2dvt = 0, svp2dvt = 0;
		double sstot = 0, sswi = 0, sswp = 0;

		SetZero(Na, k2);
		for (int p = 0; p < n; ++p)
		{
			int vp = grps[p]->loc_stat1[l].nhaplo;
			if (vp == 0) continue;

			REAL* fp = grps[p]->GetFreq(l);
			double vi2 = 0;

			if (layer == 3)
			{
				ushort* gcount = grps[p]->GetGenoCount(l);
				GENOTYPE* gtab = GetLoc(l).GetGtab();
				int ngeno = GetLoc(l).ngeno;

				for (int gi = 0; gi < ngeno; ++gi)
				{
					int gc = gcount[gi];
					if (gc == 0) continue;

					GENOTYPE& gt = gtab[gi];
					if (gt.Nalleles() == 0) continue;

					sswi += gc * (isiam ? gt.SS_IAM() : gt.SS_SMM(alen));
					Ni += gc;
					vi2 += gc * gt.Ploidy() * gt.Ploidy();
				}
			}

			Np++;
			Nh += vp;
			svi2dvp += vi2 / vp;
			svi2dvt += vi2;
			svp2dvt += vp * vp;
			AddProd(Na, fp, vp, k2);
			sswp += SSP(fp, k2, vp, true, NULL);
		}

		if (Np < 2) continue;

		svi2dvt /= Nh;
		svp2dvt /= Nh;
		sstot = SSC(Na, k2, true, NULL);

		//ns
		if (layer == 3)
		{
			n11 = Nh - Ni;
			n21 = Nh - svi2dvp; n22 = Nh - Np;
			n31 = Nh - svp2dvt; n32 = Nh - svi2dvt; n33 = Nh - 1;

			SSwi += sswp; SSwp += sswp; SStot += sstot;

			N11 += n11;
			N21 += n21; N22 += n22;
			N31 += n31; N32 += n32; N33 += n33;

			double vwi = sswi / n11;
			double vai = (n11 * sswp - n22 * sswi) / (n11 * n21);
			double vap = (n22 * n32 * sswi - n21 * n33 * sswi - n11 * n32 * sswp + n11 * n21 * sstot) / (n11 * n21 * n31);

			if (IsError(vwi) || IsError(vai) || IsError(vap)) continue;

			Vap += vap;
			Vtot += vwi + vai + vap;
			if (each) each[l] = IsError(vap / (vwi + vai + vap)) ? NAN : vap / (vwi + vai + vap);
		}
		else
		{
			n11 = Nh - Np; n21 = Nh - svp2dvt; n22 = Nh - 1;

			SSwp += sswp; SStot += sstot;

			N11 += n11; N21 += n21; N22 += n22;

			double vwp = sswp / n11;
			double vap = (n11 * sstot - n22 * sswp) / (n11 * n21);

			if (IsError(vwp) || IsError(vap)) continue;

			Vap += vap;
			Vtot += vwp + vap;
			if (each) each[l] = IsError(vap / (vwp + vap)) ? NAN : vap / (vwp + vap);
		}

	}

	if (layer == 3 && sumss)
	{
		double Vwi = SSwi / N11;
		double Vai = (N11 * SSwp - N22 * SSwi) / (N11 * N21);
		Vap = (N22 * N32 * SSwi - N21 * N33 * SSwi - N11 * N32 * SSwp + N11 * N21 * SStot) / (N11 * N21 * N31);
		Vtot = Vwi + Vai + Vap;
	}
	else if (sumss)
	{
		double Vwp = SSwp / N11;
		Vap = (N11 * SStot - N22 * SSwp) / (N11 * N21);
		Vtot = Vwp + Vap;
	}

	if (_l != -1) { *frac1 = Vap; *frac2 = Vtot; }
	return Vap / Vtot;
}

/* Absolute divergence */
template<typename REAL>
TARGET double FST<REAL>::dxy(POP<REAL>** grps, int n, double* buf, int64 l)
{
	//do not need to weight, allow ad
	double* Ni = buf;
	double f1 = 0, f2 = 0;
	int k2 = GetLoc(l).k;
	if (k2 < 2) return -1;

	int Nh = 0, Np = 0;
	double Dw = 0, Nw = 0;
	SetZero(Ni, k2);
	for (int i = 0; i < n; ++i)
	{
		int nhaplo = grps[i]->loc_stat1[l].nhaplo;
		if (nhaplo == 0) continue;

		REAL* p = grps[i]->GetFreq(l);
		AddProd(Ni, p, nhaplo, k2);

		Nh += nhaplo;
		Np++;

		int nh = nhaplo;
		Dw += nh * (nh - 1); //2x
		Nw += nh * (nh - 1); //2x
		for (int a = 0; a < k2; ++a)
		{
			double na = nh * p[a];
			if (na > 1) Dw -= na * (na - 1); //substrate same
		}
	}

	if (Np < 2) return -1;

	double hw = Dw / Nw; //mean diff
	double Db = (double)(Nh * (Nh - 1)), Nb = Db;
	for (int a = 0; a < k2; ++a)
		Db -= Ni[a] * (Ni[a] - 1);
	double hb = (Db - Dw) / (Nb - Nw);

	if (IsNormal(hw) && IsNormal(hb))
	{
		f1 += hb - hw;
		f2 += hb;
		return f1 / f2;
	}
	
	return -1;
	
}

/* Write results file in column format */
template<typename REAL>
TARGET void FST<REAL>::ColumnFormat(FILE* fout)
{
	byte* estimator = fst_estimator_val;
	char name_buf[NAME_BUF_LEN];
#define FORBEGIN \
		for (int type = 1; type <= 5; ++type)\
		{\
			if (fst_level_val[type] == 0) continue;\
			if ((type == 1 || type == 4 || type == 3) && lreg == 0) continue;\
			FST<REAL>* Fst = fst_buf<REAL>[type];\
			int n0 = 0, n1 = 0, n2 = 0;\
			switch (type)\
			{\
			case 1: n0 = lreg; n1 = 1; n2 = 1; break;\
			case 2: n0 = 1; n1 = 1; n2 = 1; break;\
			case 3: n0 = lreg; n1 = 1; n2 = 1; break;\
			case 4: n0 = lreg; n1 = nregt2; n2 = nregt2; break;\
			case 5: n0 = 1; n1 = npop; n2 = npop; break;\
			}\
			for (int rl = (int)n0 - 1; rl >= 0; --rl)\
				{\
					if (type == 3) n2 = nreg[rl];\
					if (type == 4) n1 = n2 = nreg[rl];\
					for (int i = 0; i < n1; ++i)\
					{\
						for (int j = type <= 3 ? 0 : i + 1; j < n2; ++j)\
						{\
							FST* g = Fst + i * n1 + j;
#define FOREND }}}}

	//type 1 among reg in tot, 2 among pop in tot, 3 among pop/reg in reg, 4 between reg, 5 between pop

	fprintf(fout, "%s%s", g_linebreak_val, g_linebreak_val);

	//Line 1
	fprintf(fout, "Locus%cA", g_delimiter_val);
	FORBEGIN
		switch (type)
		{
		case 1: fprintf(fout, "%cAmong all regL%d", g_delimiter_val, rl + 1); break;
		case 2: fprintf(fout, "%cAmong all pops", g_delimiter_val); break;
		case 3:
			if (rl == 0) fprintf(fout, "%cAmong pops in %s", g_delimiter_val, aregs<REAL>[rl][j]->name);
			else		 fprintf(fout, "%cAmong regsL%d in %s", g_delimiter_val, rl, aregs<REAL>[rl][j]->name);
			break;
		case 4:  fprintf(fout, "%c%s", g_delimiter_val, aregs<REAL>[rl][i]->name); break;
		case 5:  fprintf(fout, "%c%s", g_delimiter_val, apops<REAL>[i]->name); break;
		}
	FOREND
		fprintf(fout, "%s", g_linebreak_val);

	//Line 2
	fprintf(fout, "%cB", g_delimiter_val);
	FORBEGIN
		switch (type)
		{
		case 1: fprintf(fout, "%c", g_delimiter_val); break;
		case 2: fprintf(fout, "%c", g_delimiter_val); break;
		case 3:
			if (rl == 0) fprintf(fout, "%c", g_delimiter_val);
			else		 fprintf(fout, "%c", g_delimiter_val);
			break;
		case 4:  fprintf(fout, "%c%s", g_delimiter_val, aregs<REAL>[rl][j]->name); break;
		case 5:  fprintf(fout, "%c%s", g_delimiter_val, apops<REAL>[j]->name); break;
		}
	FOREND
		fprintf(fout, "%s", g_linebreak_val);

	//Line 3
	if (fst_locus_val[1])
	{
		fprintf(fout, "All loci");

		for (int k = 1; k <= N_FST_ESTIMATOR; ++k)
		{
			if (estimator[k] == 0) continue;
			fprintf(fout, "%c%s", g_delimiter_val, FST_ESTIMATOR[k]);

			FORBEGIN
				fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, *((&g->Nei1973T) + k - 1));
			FOREND
				fprintf(fout, "%s", g_linebreak_val);
		}

		if (fst_test_val[1])
		{
			fprintf(fout, "%cGenotype G", g_delimiter_val);
			FORBEGIN
				fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, g->Genotype_GT);
			FOREND
				fprintf(fout, "%s", g_linebreak_val);

			fprintf(fout, "%cd.f.", g_delimiter_val);
			FORBEGIN
				fprintf(fout, "%c", g_delimiter_val);
			fprintf(fout, "%d", g->Genotype_DFT);
			FOREND
				fprintf(fout, "%s", g_linebreak_val);

			fprintf(fout, "%cP", g_delimiter_val);
			FORBEGIN
				fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, g->Genotype_PT);
			FOREND
				fprintf(fout, "%s", g_linebreak_val);
		}

		if (fst_test_val[2])
		{
			fprintf(fout, "%cAllele G", g_delimiter_val);
			FORBEGIN
				fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, g->Allele_GT);
			FOREND
				fprintf(fout, "%s", g_linebreak_val);

			fprintf(fout, "%cd.f.", g_delimiter_val);
			FORBEGIN
				fprintf(fout, "%c", g_delimiter_val);
			fprintf(fout, "%d", g->Allele_DFT);
			FOREND
				fprintf(fout, "%s", g_linebreak_val);

			fprintf(fout, "%cP", g_delimiter_val);
			FORBEGIN
				fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, g->Allele_PT);
			FOREND
				fprintf(fout, "%s", g_linebreak_val);
		}
	}


	if (fst_locus_val[2]) for (int64 l = 0; l < nloc; ++l)
	{
		fprintf(fout, "%s", GetLoc(l).GetNameStr(name_buf));

		for (int k = 1; k <= N_FST_ESTIMATOR; ++k)
		{
			if (estimator[k] == 0) continue;
			fprintf(fout, "%c%s", g_delimiter_val, FST_ESTIMATOR[k]);

			FORBEGIN
				fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, *(*((&g->Nei1973) + k - 1) + l));
			FOREND
				fprintf(fout, "%s", g_linebreak_val);
		}


		if (fst_test_val[1])
		{
			fprintf(fout, "%cGenotype G", g_delimiter_val);
			FORBEGIN
				fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, g->Genotype_G[l]);
			FOREND
				fprintf(fout, "%s", g_linebreak_val);

			fprintf(fout, "%cd.f.", g_delimiter_val);
			FORBEGIN
				fprintf(fout, "%c", g_delimiter_val);
			fprintf(fout, "%d", g->Genotype_DF[l]);
			FOREND
				fprintf(fout, "%s", g_linebreak_val);

			fprintf(fout, "%cP", g_delimiter_val);
			FORBEGIN
				fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, g->Genotype_P[l]);
			FOREND
				fprintf(fout, "%s", g_linebreak_val);
		}

		if (fst_test_val[2])
		{
			fprintf(fout, "%cAllele G", g_delimiter_val);
			FORBEGIN
				fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, g->Allele_G[l]);
			FOREND
				fprintf(fout, "%s", g_linebreak_val);

			fprintf(fout, "%cd.f.", g_delimiter_val);
			FORBEGIN
				fprintf(fout, "%c", g_delimiter_val);
			fprintf(fout, "%d", g->Allele_DF[l]);
			FOREND
				fprintf(fout, "%s", g_linebreak_val);

			fprintf(fout, "%cP", g_delimiter_val);
			FORBEGIN
				fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, g->Allele_P[l]);
			FOREND
				fprintf(fout, "%s", g_linebreak_val);
		}
	}
}

/* Write results file in matrix format */
template<typename REAL>
TARGET void FST<REAL>::MatrixFormat(FILE* fout, FST* Fst, int n, int type)
{
	byte* estimator = fst_estimator_val;
	for (int rl = type == 4 ? lreg - 1 : 0; rl >= 0; --rl)
	{
		int n1 = type == 4 ? nreg[rl] : n;
		for (int k = 1; k <= N_FST_ESTIMATOR; ++k)
		{
			if (k <= N_FST_ESTIMATOR && estimator[k] == 0) continue;

			fprintf(fout, "%s%s%s", g_linebreak_val, g_linebreak_val, FST_ESTIMATOR[k]);

			for (int i = 0; i < n1; ++i)
			{
				switch (type)
				{
				case 4: fprintf(fout, "%c%s", g_delimiter_val, aregs<REAL>[rl][i]->name);  break;
				case 5: fprintf(fout, "%c%s", g_delimiter_val, apops<REAL>[i]->name);  break;
				}
			}

			FST* g = Fst;
			for (int i = 0; i < n1; ++i)
			{
				switch (type)
				{
				case 4: fprintf(fout, "%s%s", g_linebreak_val, aregs<REAL>[rl][i]->name);  break;
				case 5: fprintf(fout, "%s%s", g_linebreak_val, apops<REAL>[i]->name);  break;
				}

				for (int j = 0; j < n1; ++j)
				{
					fprintf(fout, "%c", g_delimiter_val);
					WriteReal(fout, *((&g->Nei1973T) + k - 1));
					g++;
				}
			}
		}
		Fst += n1 * n1;
	}
}
#endif

#define extern 
template<typename REAL>
extern FST<REAL>* fst_buf[6];							//Read/Write buffer for each type of fst
extern int fst_type;									//1 Among regions, 2 among pops, 3 among pops/regs in region, 4 between regions, 5 between pops 
#undef extern 


/* Calculate genetic differentiation */
template<typename REAL>
TARGET void CalcDiff()
{
	if (!fst) return;

	EvaluationBegin();
	OpenResFile("-fst", "Genetic differentiation");
	OpenTempFiles(1, ".fst");

	bool isfirst = true; int64 ntot = 0;
	if (fst_level_val[1] && lreg > 0) ntot += lreg;
	if (fst_level_val[2]) ntot += 1;
	if (fst_level_val[3] && lreg > 0) ntot += nregt;
	if (fst_level_val[4] && lreg > 0) ntot += (nregt2 - nregt) / 2;
	if (fst_level_val[5]) ntot += npop * (npop - 1) / 2;

	for (fst_type = 1; fst_type <= 5; ++fst_type)
	{
		if (fst_level_val[fst_type] == 0) continue;
		if ((fst_type == 1 || fst_type == 4 || fst_type == 3) && lreg == 0) continue;

		int n1 = 0, n2 = 0;
		switch (fst_type)
		{
		case 1: n1 = lreg >= 0 ? lreg : 0;  n2 = lreg >= 0 ? lreg : 0; break;
		case 2: n1 = 1;					    n2 = 1;    break;
		case 3: n1 = nregt;					n2 = nregt; break;
		case 4: n1 = nregt2;				n2 = (nregt2 - nregt) / 2;  break;
		case 5: n1 = npop * npop;			n2 = npop * (npop - 1) / 2; break;
		}

		if (n1 == 0 || n2 == 0) continue;

		fst_buf<REAL>[fst_type] = new FST<REAL>[n1];
		SetZero(fst_buf<REAL>[fst_type], n1);
		RunThreads(&GeneticDifferentiationThread<REAL>, NULL, NULL, ntot, n2,
			"\nCalculating genetic differentiation:\n", g_nthread_val, isfirst);
		isfirst = false;

		if ((fst_type == 4 || fst_type == 5) && fst_fmt_val[1])
			FST<REAL>::MatrixFormat(TEMP_FILES[0], fst_buf<REAL>[fst_type], npop, fst_type);
	}

	FST<REAL>::ColumnFormat(FRES);

	JoinTempFiles(1);
	CloseResFile();

	for (fst_type = 1; fst_type <= 5; ++fst_type)
	{
		if (fst_level_val[fst_type] == 0) continue;
		if ((fst_type == 1 || fst_type == 4) && lreg == 0) continue;

		int n1 = 0, n2 = 0;
		switch (fst_type)
		{
		case 1: n1 = lreg;  n2 = lreg; break;
		case 2: n1 = 1;     n2 = 1;    break;
		case 3: n1 = nregt; n2 = nregt; break;
		case 4: n1 = nregt2; n2 = (nregt2 - nregt) / 2;  break;
		case 5: n1 = npop * npop; n2 = npop * (npop - 1) / 2; break;
		}

		if (fst_type <= 3)
		{
			for (int ni = 0; ni < n1; ++ni)
				fst_buf<REAL>[fst_type][ni].Uninitialize();
		}
		else if (fst_type == 4)
		{
			auto fst_buf2 = fst_buf<REAL>[fst_type];
			for (int rl = lreg - 1; rl >= 0; --rl)
			{
				for (int i = 0; i < nreg[rl]; ++i)
					for (int j = i + 1; j < nreg[rl]; ++j)
						fst_buf2[i * nreg[rl] + j].Uninitialize();
				fst_buf2 += nreg[rl] * nreg[rl];
			}
		}
		else if (fst_type == 5)
		{
			auto fst_buf2 = fst_buf<REAL>[fst_type];
			for (int i = 0; i < npop; ++i)
				for (int j = i + 1; j < npop; ++j)
					fst_buf2[i * npop + j].Uninitialize();
		}
		DEL(fst_buf<REAL>[fst_type]);
	}

	EvaluationEnd("Genetic differentiation estimation");

	if (fst_plot_val == 1)
		RunRscript("fst_plot.R");
}

/* Calculate genetic differentiation using multiple threads */
THREAD2(GeneticDifferentiationThread)
{
	//load ind
	int nthread = g_nthread_val;
	int64 progress = 0;

	FST<REAL>* fst_buf2 = fst_buf<REAL>[fst_type];
	if (fst_type == 1)
	{
		//among regions
		for (int rl = 0; rl < lreg; ++rl)
			if (progress++ % nthread == threadid)
				fst_buf2[lreg - 1 - rl].CalcFst(aregs<REAL>[rl], nreg[rl]);
	}
	if (fst_type == 2)
	{
		//among pops
		if (progress++ % nthread == threadid)
			fst_buf2[0].CalcFst(apops<REAL>, npop);
	}
	if (fst_type == 3)
	{
		//among pops/regs in region
		for (int rl = lreg - 1, p = 0; rl >= 0; --rl)
			for (int i = 0; i < nreg[rl]; ++i, ++p)
				if (progress++ % nthread == threadid)
					fst_buf2[p].CalcFst(aregs<REAL>[rl][i]->vpop, aregs<REAL>[rl][i]->npop);
	}
	if (fst_type == 4)
	{
		//between regions
		for (int rl = lreg - 1; rl >= 0; --rl)
		{
			int n = nreg[rl];
			for (int i = 0; i < n; ++i)
				for (int j = i + 1; j < n; ++j)
					if (progress++ % nthread == threadid)
					{
						fst_buf2[i * n + j].CalcFst(aregs<REAL>[rl][i], aregs<REAL>[rl][j]);
						fst_buf2[j * n + i] = fst_buf2[i * n + j];
					}
			fst_buf2 += n * n;
		}
	}
	if (fst_type == 5)
	{
		//between pops
		int n = npop;
		for (int i = 0; i < n; ++i)
			for (int j = i + 1; j < n; ++j)
				if (progress++ % nthread == threadid)
				{
					fst_buf2[i * n + j].CalcFst(apops<REAL>[i], apops<REAL>[j]);
					fst_buf2[j * n + i] = fst_buf2[i * n + j];
				}
	}
}