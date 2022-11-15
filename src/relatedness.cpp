/* Relatedness Functions */

#pragma once
#include "vcfpop.h"

struct Huang2015ENTRY;
template struct RELATEDNESS<double>;
template struct RELATEDNESS<float >;
template TARGET void CalcRelatedness<double>();
template TARGET void CalcRelatedness<float >();

#ifndef _RELATEDNESS
/* Write header row for relatedness estimation */
template<typename REAL>
TARGET void RELATEDNESS<REAL>::ColumnPrintHeader()
{
	fprintf(FRES, "%s%s%s%sA%cpop",
		g_linebreak_val, g_linebreak_val,
		cpop->name, g_linebreak_val,
		g_delimiter_val);
	for (int rl = 0; rl < lreg; ++rl)
		fprintf(FRES, "%cregL%d", g_delimiter_val, rl + 1);
	fprintf(FRES, "%cB%cpop", g_delimiter_val, g_delimiter_val);
	for (int rl = 0; rl < lreg; ++rl)
		fprintf(FRES, "%cregL%d", g_delimiter_val, rl + 1);
	fprintf(FRES, "%cAB_typed%cA_typed%cB_typed", g_delimiter_val, g_delimiter_val, g_delimiter_val);

	for (int k = 1; k <= N_RELATEDNESS_ESTIMATOR; ++k)
		if (relatedness_estimator_val[k])
			fprintf(FRES, "%c%s", g_delimiter_val, RELATEDNESS_ESTIMATOR[k]);
}

/* Write result row for relatedness estimation */
template<typename REAL>
TARGET void RELATEDNESS<REAL>::ColumnPrintLine(int i, int j)
{
	fprintf(FRES, "%s%s%c%s%c",
		g_linebreak_val,
		ainds[i]->name, g_delimiter_val,
		apops[ainds[i]->popid]->name, g_delimiter_val);

	POP<REAL>* tr = lreg >= 0 ? aregs[0][apops[ainds[i]->popid]->rid] : NULL;
	for (int rl = 0; rl < lreg; ++rl)
	{
		fprintf(FRES, "%s%c", tr->name, g_delimiter_val);
		tr = aregs[rl + 1][tr->rid];
	}

	fprintf(FRES, "%s%c%s%c",
		ainds[j]->name, g_delimiter_val,
		apops[ainds[j]->popid]->name, g_delimiter_val);

	tr = lreg >= 0 ? aregs[0][apops[ainds[j]->popid]->rid] : NULL;
	for (int rl = 0; rl < lreg; ++rl)
	{
		fprintf(FRES, "%s%c", tr->name, g_delimiter_val);
		tr = aregs[rl + 1][tr->rid];
	}

	fprintf(FRES, "%d%c%d%c%d",
		ABtype, g_delimiter_val,
		Atype, g_delimiter_val,
		Btype);

	for (int k = 1; k <= N_RELATEDNESS_ESTIMATOR; ++k)
		if (relatedness_estimator_val[k])
		{
			fprintf(FRES, "%c", g_delimiter_val);
			WriteReal(FRES, *((&Lynch1999) + k - 1));
		}
}

/* Write matrix format header for relatedness estimation */
template<typename REAL>
TARGET void RELATEDNESS<REAL>::MatrixPrintMatrixHeader(int k, int n)
{
	if (relatedness_estimator_val[k] == 0) return;
	fprintf(TEMP_FILES[k], "%s%s%s%s%s", g_linebreak_val, g_linebreak_val, cpop->name, g_linebreak_val, RELATEDNESS_ESTIMATOR[k]);

	for (int i = 0; i < n; ++i)
		fprintf(TEMP_FILES[k], "%c%s", g_delimiter_val, cpop->inds[i]->name);
}

/* Write matrix format row header for relatedness estimation */
template<typename REAL>
TARGET void RELATEDNESS<REAL>::MatrixPrintRowHeader(int k, int i)
{
	if (relatedness_estimator_val[k] == 0) return;
	fprintf(TEMP_FILES[k], "%s%s", g_linebreak_val, cpop->inds[i]->name);
}

/* Write matrix format grid for relatedness estimation */
template<typename REAL>
TARGET void RELATEDNESS<REAL>::MatrixPrintCell(int k)
{
	if (relatedness_estimator_val[k] == 0) return;
	fprintf(TEMP_FILES[k], "%c", g_delimiter_val);
	WriteReal(TEMP_FILES[k], *((&Lynch1999) + k - 1));
}

/* Calculate relatedness coefficient */
template<typename REAL>
TARGET void RELATEDNESS<REAL>::CalcRelatedness(IND<REAL>* x, IND<REAL>* y)
{
	SetZero(&Lynch1999, N_RELATEDNESS_ESTIMATOR);
	ABtype = Atype = Btype = 0;
	byte* estimator = relatedness_estimator_val;

	for (int64 l = 0; l < nloc; ++l)
	{
		GENOTYPE& gt1 = x->GetGenotype(l), & gt2 = y->GetGenotype(l);//fine
		if (gt1.Nalleles()) Atype++;
		if (gt2.Nalleles()) Btype++;
		if (gt1.Nalleles() && gt2.Nalleles()) ABtype++;
	}

	if (ABtype == 0)
	{
		for (int k = 1; k <= N_RELATEDNESS_ESTIMATOR; ++k)
			*(&Lynch1999 + k - 1) = NAN;
		return;
	}

	for (int k = 1; k <= N_RELATEDNESS_ESTIMATOR; ++k)
		if (estimator[k])
			*(&Lynch1999 + k - 1) = RelatednessEstimator(k, x, y);
}

/* Relatedness estimator Warpper */
template<typename REAL>
TARGET double RELATEDNESS<REAL>::RelatednessEstimator(int k, IND<REAL>* x, IND<REAL>* y)
{
	switch (k)
	{
	case 1: return R_Lynch1999(x, y);
	case 2: return R_Wang2002(x, y);
	case 3: return R_Thomas2010(x, y);
	case 4: return R_Li1993(x, y);
	case 5: return R_Queller1989(x, y);
	case 6: return R_Huang2016A(x, y);
	case 7: return R_Huang2016B(x, y);
	case 8: return R_Milligan2003(x, y);
	case 9: return R_Anderson2007(x, y, true);
	case 10: return R_Huang2014(x, y);
	case 11: return R_Huang2015(x, y);
	case 12: return R_Ritland1996(x, y, true, true);
	case 13: return R_Loiselle1995(x, y, true, true);
	case 14: return R_Ritland1996(x, y, false, true);
	case 15: return R_Loiselle1995(x, y, false, true);
	case 16: return R_Weir1996(x, y, true);
	}
	return NAN;
}

/* Lynch 1999 relatedness estimator */
template<typename REAL>
TARGET double RELATEDNESS<REAL>::R_Lynch1999(IND<REAL>* x, IND<REAL>* y)
{
	if (x->vmax > 2 || x->vmin < 2)
		Exit("\nError: Lynch1999 relatedness estimators only supports diploids, in individual %s.\n", x->name);
	if (y->vmax > 2 || y->vmin < 2)
		Exit("\nError: Lynch1999 relatedness estimators only supports diploids, in individual %s.\n", y->name);

	double srx = 0, swx = 0, sry = 0, swy = 0;

	for (int64 l = 0; l < nloc; ++l)
	{
		GENOTYPE& gx = x->GetGenotype(l), & gy = y->GetGenotype(l);//fine
		if (gx.Nalleles() == 0 || gy.Nalleles() == 0) continue;
		if (cpop->loc_stat1[l].k <= 1) continue;

		REAL* p = cpop->GetFreq(l);
		int a = gx.GetAlleleCopy(0), b = gx.GetAlleleCopy(1), c = gy.GetAlleleCopy(0), d = gy.GetAlleleCopy(1);
		int dab = a == b, dcd = c == d, dbc = b == c, dbd = b == d, dac = a == c, dad = a == d;
		double wx = ((1 + dab) * (p[a] + p[b]) - 4 * p[a] * p[b]) / (2 * p[a] * p[b]);
		double rx = (p[a] * (dbc + dbd) + p[b] * (dac + dad) - 4 * p[a] * p[b]) / ((1 + dab) * (p[a] + p[b]) - 4 * p[a] * p[b]);
		if (IsNormal(rx) && IsNormal(wx))
		{
			srx += wx * rx;
			swx += wx;
		}
		double wy = ((1 + dcd) * (p[c] + p[d]) - 4 * p[c] * p[d]) / (2 * p[c] * p[d]);
		double ry = (p[c] * (dad + dbd) + p[d] * (dac + dbc) - 4 * p[c] * p[d]) / ((1 + dcd) * (p[c] + p[d]) - 4 * p[c] * p[d]);
		if (IsNormal(ry) && IsNormal(wy))
		{
			sry += wy * ry;
			swy += wy;
		}
	}

	return (srx / swx + sry / swy) * 0.5;
}

/* Wang 2002 relatedness estimator */
template<typename REAL>
TARGET double RELATEDNESS<REAL>::R_Wang2002(IND<REAL>* x, IND<REAL>* y)
{
	if (x->vmax > 2 || x->vmin < 2)
		Exit("\nError: Wang2002 relatedness estimators only supports diploids, in individual %s.\n", x->name);
	if (y->vmax > 2 || y->vmin < 2)
		Exit("\nError: Wang2002 relatedness estimators only supports diploids, in individual %s.\n", y->name);

	double sw = 0, P1 = 0, P2 = 0, P3 = 0, P4 = 0;
	double a2 = 0, a3 = 0, a4 = 0, a22 = 0;

	for (int64 l = 0; l < nloc; ++l)
	{
		GENOTYPE& gx = x->GetGenotype(l), & gy = y->GetGenotype(l);//fine
		if (gx.Nalleles() == 0 || gy.Nalleles() == 0) continue;
		if (cpop->loc_stat1[l].k <= 1) continue;
		LOCSTAT2<REAL>& stat2 = relatedness_loc_stat[l];

		double w = 1.0 / (2 * stat2.s2 - stat2.s3);
		sw += w;

		a2 += w * stat2.s2;
		a3 += w * stat2.s3;
		a4 += w * stat2.s4;
		a22 += w * stat2.s2 * stat2.s2;

		int a = gx.GetAlleleCopy(0), b = gx.GetAlleleCopy(1), c = gy.GetAlleleCopy(0), d = gy.GetAlleleCopy(1);

		if ((a == c && b == d) || (a == d && b == c)) P1 += w;
		else if ((a == b || c == d) && (a == c || b == d || a == d || b == c)) P2 += w;
		else if (a == c || a == d || b == c || b == d) P3 += w;
		else P4 += w;
	}

	a2 /= sw; a3 /= sw; a4 /= sw; a22 /= sw;
	P1 /= sw; P2 /= sw; P3 /= sw; P4 /= sw;

	double b = 2 * a22 - a4;//e1
	double c = a2 - 2 * a22 + a4;//e5
	double d = 4 * (a3 - a4);//e2
	double e = 2 * (a2 - 3 * a3 + 2 * a4);//e4
	double f = 8 * (a4 - a3) + 4 * (a2 - a22);//e3
	double g = 1 - 3 * a2 + 2 * a3 - f;//e6
	double V = (1 - b) * (1 - b) * (e * e * f + d * g * g) - (1 - b) * (e * f - d * g) * (e * f - d * g) + 2 * c * d * f * (1 - b) * (g + e) + c * c * d * f * (d + f);
	double phi = ((d * f * ((e + g) * (1 - b) + c * (d + f)) * (P1 - 1)
		+ d * (1 - b) * (g * (1 - b - d) + f * (c + e)) * P3
		+ f * (1 - b) * (e * (1 - b - f) + d * (c + g)) * P2) / V);
	double delta = ((c * d * f * (e + g) * (P1 + 1 - 2 * b)
		+ ((1 - b) * (f * e * e + d * g * g) - (e * f - d * g) * (e * f - d * g)) * (P1 - b)
		+ c * (d * g - e * f) * (d * P3 - f * P2) - c * c * d * f * (P3 + P2 - d - f) - c * (1 - b) * (d * g * P3 + e * f * P2)) / V);

	return delta + phi / 2;
}

/* Thomas 2010 relatedness estimator */
template<typename REAL>
TARGET double RELATEDNESS<REAL>::R_Thomas2010(IND<REAL>* x, IND<REAL>* y)
{
	if (x->vmax > 2 || x->vmin < 2)
		Exit("\nError: Thomas2010 relatedness estimators only supports diploids, in individual %s.\n", x->name);
	if (y->vmax > 2 || y->vmin < 2)
		Exit("\nError: Thomas2010 relatedness estimators only supports diploids, in individual %s.\n", y->name);

	int t1, t2;
	double sr = 0, sw = 0;//, D = 0, w4 = 0;
	double P[3] = { 0 };

	for (int64 l = 0; l < nloc; ++l)
	{
		GENOTYPE& gx = x->GetGenotype(l), & gy = y->GetGenotype(l);//fine
		if (gx.Nalleles() == 0 || gy.Nalleles() == 0) continue;
		if (cpop->loc_stat1[l].k <= 1) continue;
		LOCSTAT2<REAL>& stat2 = relatedness_loc_stat[l];

		int va = gx.GetAlleleCopy(0), vb = gx.GetAlleleCopy(1), vc = gy.GetAlleleCopy(0), vd = gy.GetAlleleCopy(1);
		t1 = (va == vc) + (vb == vd);
		t2 = (va == vd) + (vb == vc);
		P[0] = P[1] = P[2] = 0;
		P[t1 >= t2 ? 2 - t1 : 2 - t2] = 1;

		double b = 2 * stat2.s2 * stat2.s2 - stat2.s4;
		double c = stat2.s2 - 2 * stat2.s2 * stat2.s2 + stat2.s4;
		double d = 4 * (stat2.s2 - stat2.s2 * stat2.s2 - stat2.s3 + stat2.s4);
		double e = 1 - 5 * stat2.s2 + 4 * stat2.s2 * stat2.s2 + 4 * stat2.s3 - 4 * stat2.s4;
		double div = c * d + (1 - b) * e;
		double vr = (4 * b * e * e - 4 * (c * d - b * e) * (c * d - b * e) - 4 * c * d * e + d * (1 - b) * (1 - b) - d * d * (1 - b)) / (4 * div * div);

		double rr = ((P[0] - b) * (d * 0.5 + e) + (P[1] - d) * (0.5 - 0.5 * b - c)) / vr / div;
		if (IsNormal(rr) && IsNormal(vr))
		{
			sw += 1 / vr;
			sr += rr;
		}
		//double v4 = (d * c * c + b * e * e - (c * d - e * b) * (c * d - e * b)) / (div * div);
		//D += (e * (P[0] - b) - c * (P[1] - d)) / v4 / div;
		//w4 += 1 / v4;
	}
	//D /= w4;
	return sr / sw;
}

/* Li 1993 relatedness estimator */
template<typename REAL>
TARGET double RELATEDNESS<REAL>::R_Li1993(IND<REAL>* x, IND<REAL>* y)
{
	if (x->vmax > 2 || x->vmin < 2)
		Exit("\nError: Li1993 relatedness estimators only supports diploids, in individual %s.\n", x->name);
	if (y->vmax > 2 || y->vmin < 2)
		Exit("\nError: Li1993 relatedness estimators only supports diploids, in individual %s.\n", y->name);

	double sr = 0, sw = 0;

	for (int64 l = 0; l < nloc; ++l)
	{
		GENOTYPE& gx = x->GetGenotype(l), & gy = y->GetGenotype(l);//fine
		if (gx.Nalleles() == 0 || gy.Nalleles() == 0) continue;
		if (cpop->loc_stat1[l].k <= 1) continue;
		LOCSTAT2<REAL>& stat2 = relatedness_loc_stat[l];

		int va = gx.GetAlleleCopy(0), vb = gx.GetAlleleCopy(1), vc = gy.GetAlleleCopy(0), vd = gy.GetAlleleCopy(1);

		double a2 = stat2.s2, a3 = stat2.s3;
		double S0 = 2 * a2 - a3;
		double Sxy = 0;

		if ((va == vc && vb == vd) || (va == vd && vb == vc))								Sxy = 1;
		else if ((va == vb || vc == vd) && (va == vc || vb == vd || va == vd || vb == vc))	Sxy = 0.75;
		else if (va == vc || va == vd || vb == vc || vb == vd)								Sxy = 0.5;
		else																				Sxy = 0;

		double rr = (Sxy - S0) / (1 - S0);
		if (IsNormal(rr)) { sr += (Sxy - S0) / (1 - S0); sw++; }
	}
	return sr / sw;
}

/* Queller 1989 relatedness estimator */
template<typename REAL>
TARGET double RELATEDNESS<REAL>::R_Queller1989(IND<REAL>* x, IND<REAL>* y)
{
	if (x->vmax > 2 || x->vmin < 2)
		Exit("\nError: Queller1989 relatedness estimators only supports diploids, in individual %s.\n", x->name);
	if (y->vmax > 2 || y->vmin < 2)
		Exit("\nError: Queller1989 relatedness estimators only supports diploids, in individual %s.\n", y->name);

	double srx = 0, swx = 0, sry = 0, swy = 0;

	for (int64 l = 0; l < nloc; ++l)
	{
		GENOTYPE& gx = x->GetGenotype(l), & gy = y->GetGenotype(l);//fine
		if (gx.Nalleles() == 0 || gy.Nalleles() == 0) continue;
		if (cpop->loc_stat1[l].k <= 2) continue;//cannot be used for diallelic loci

		int va = gx.GetAlleleCopy(0), vb = gx.GetAlleleCopy(1), vc = gy.GetAlleleCopy(0), vd = gy.GetAlleleCopy(1);
		int Sac = va == vc, Sad = va == vd, Sbc = vb == vc, Sbd = vb == vd, Sab = va == vb, Scd = vc == vd;
		REAL* p = cpop->GetFreq(l);

		double rx = (0.5 * (Sac + Sad + Sbc + Sbd) - p[va] - p[vb]) / (1 + Sab - p[va] - p[vb]);
		if (IsNormal(rx)) { srx += rx; swx++; }
		double ry = (0.5 * (Sac + Sad + Sbc + Sbd) - p[vc] - p[vd]) / (1 + Scd - p[vc] - p[vd]);
		if (IsNormal(ry)) { sry += ry; swy++; }
	}
	return (srx / swx + sry / swy) * 0.5;
}

/* Huang 2016 relatedness estimator A */
template<typename REAL>
TARGET double RELATEDNESS<REAL>::R_Huang2016A(IND<REAL>* x, IND<REAL>* y)
{
	if (x->vmax > 2 || x->vmin < 2)
		Exit("\nError: Huang2016A relatedness estimators only supports diploids, in individual %s.\n", x->name);
	if (y->vmax > 2 || y->vmin < 2)
		Exit("\nError: Huang2016A relatedness estimators only supports diploids, in individual %s.\n", y->name);

	double c1 = 0, c2 = 0, c3 = 0, c4 = 0, c5 = 0, c6 = 0;
	double S = 0, S2 = 0, sw = 0;

	for (int64 l = 0; l < nloc; ++l)
	{
		GENOTYPE& gx = x->GetGenotype(l), & gy = y->GetGenotype(l);//fine
		if (gx.Nalleles() == 0 || gy.Nalleles() == 0) continue;
		if (cpop->loc_stat1[l].k <= 1) continue;
		LOCSTAT2<REAL>& stat2 = relatedness_loc_stat[l];

		int va = gx.GetAlleleCopy(0), vb = gx.GetAlleleCopy(1), vc = gy.GetAlleleCopy(0), vd = gy.GetAlleleCopy(1);
		double w = 1.0 / (2 * stat2.s2 - stat2.s3);
		if (IsError(w)) continue;

		sw += w;
		double s = 0;
		if ((va == vc && vb == vd) || (va == vd && vb == vc))								s = 1;
		else if ((va == vb || vc == vd) && (va == vc || vb == vd || va == vd || vb == vc))	s = 0.75;
		else if (va == vc || va == vd || vb == vc || vb == vd)								s = 0.5;
		else																				s = 0;

		S += w * s;
		S2 += w * s * s;

		c1 += w * (1 - 2 * stat2.s2 + stat2.s3);
		c2 += w * (0.5 * (1 - 2 * stat2.s2 + stat2.s3));
		c3 += w * (2 * stat2.s2 - stat2.s3);
		c4 += w * (0.25 * (4 - 4 * stat2.s2 - 4 * stat2.s2 * stat2.s2 - stat2.s3 + 5 * stat2.s4));
		c5 += w * (0.125 * (2 + 3 * stat2.s2 - 8 * stat2.s2 * stat2.s2 - 7 * stat2.s3 + 10 * stat2.s4));
		c6 += w * (stat2.s2 + stat2.s2 * stat2.s2 + 0.25 * (stat2.s3 - 5 * stat2.s4));
	}
	S /= sw; S2 /= sw;
	c1 /= sw; c2 /= sw; c3 /= sw; c4 /= sw; c5 /= sw; c6 /= sw;
	double delta = (c5 * S - c2 * S2 + c2 * c6 - c3 * c5) / (c1 * c5 - c2 * c4);
	double phi = (-c4 * S + c1 * S2 - c1 * c6 + c3 * c4) / (c1 * c5 - c2 * c4);
	return phi / 2 + delta;
}

/* Huang 2016 relatedness estimator B */
template<typename REAL>
TARGET double RELATEDNESS<REAL>::R_Huang2016B(IND<REAL>* x, IND<REAL>* y)
{
	if (x->vmax > 2 || x->vmin < 2)
		Exit("\nError: Huang2016B relatedness estimators only supports diploids, in individual %s.\n", x->name);
	if (y->vmax > 2 || y->vmin < 2)
		Exit("\nError: Huang2016B relatedness estimators only supports diploids, in individual %s.\n", y->name);

#define RELAT(X) (((2 * c5 - c4) * ((X) - c3) - (2 * c2 - c1) * ((X)*(X) - c6)) / (c1 * c5 - c2 * c4) / 2)

	double c1 = 0, c2 = 0, c3 = 0, c4 = 0, c5 = 0, c6 = 0;
	double srx = 0, swx = 0, sry = 0, swy = 0;

	for (int64 l = 0; l < nloc; ++l)
	{
		GENOTYPE& gx = x->GetGenotype(l), & gy = y->GetGenotype(l);//fine
		if (gx.Nalleles() == 0 || gy.Nalleles() == 0) continue;
		if (cpop->loc_stat1[l].k <= 1) continue;

		int va = gx.GetAlleleCopy(0), vb = gx.GetAlleleCopy(1), vc = gy.GetAlleleCopy(0), vd = gy.GetAlleleCopy(1);
		REAL* p = cpop->GetFreq(l);
		double S = 0;
		if ((va == vc && vb == vd) || (va == vd && vb == vc)) S = 1;
		else if ((va == vb || vc == vd) && (va == vc || vb == vd || va == vd || vb == vc)) S = 0.75;
		else if (va == vc || va == vd || vb == vc || vb == vd) S = 0.5;
		else S = 0;

		double pi = 0, pj = 0, px = 0, X = 0, X2 = 0, w = 0;

		if (va != vb)
		{
			pi = p[va]; pj = p[vb]; px = 1 - pi - pj;
			double b = 2 * pi * pj, d = pj * pj + pi * pi, f = 2 * px * (pi + pj), c = (pi + pj) / 2, e = c, g = px;

			c1 = 1 - b - 0.75 * d - 0.5 * f; c4 = 1 - b - 0.5625 * d - 0.25 * f; c2 = c - b + 0.75 * (e - d) + 0.5 * (g - f);
			c5 = c - b + 0.5625 * (e - d) + 0.25 * (g - f); c3 = b + 0.75 * d + 0.5 * f; c6 = b + 0.5625 * d + 0.25 * f;

			double x1 = RELAT(1), x2 = RELAT(0.75), x3 = RELAT(0.5), x4 = RELAT(0);
			X = x1 * b + x2 * d + x3 * f + x4 * (1 - b - d - f);
			X2 = x1 * x1 * b + x2 * x2 * d + x3 * x3 * f + x4 * x4 * (1 - b - d - f);
			w = 1.0 / (X2 - X * X);
		}
		else
		{
			pi = p[va]; px = 1 - pi;
			double b = pi * pi, d = 2 * pi * px, c = (2 * pi * pi) / (2 * pi), e = (2 * pi * px) / (2 * pi);

			c1 = 1 - b - 0.75 * d; c4 = 1 - b - 0.5625 * d; c2 = c - b + 0.75 * (e - d);
			c5 = c - b + 0.5625 * (e - d); c3 = b + 0.75 * d; c6 = b + 0.5625 * d;

			double x1 = RELAT(1), x2 = RELAT(0.75), x3 = RELAT(0);
			X = x1 * b + x2 * d + x3 * (1 - b - d);
			X2 = x1 * x1 * b + x2 * x2 * d + x3 * x3 * (1 - b - d);
			w = 1.0 / (X2 - X * X);
		}

		double rx = w * RELAT(S);
		if (abs(X) < 1e-10 && IsNormal(w) && IsNormal(rx))
		{
			srx += rx;
			swx += w;
		}

		if (vc != vd)
		{
			pi = p[vc]; pj = p[vd]; px = 1 - pi - pj;
			double b = 2 * pi * pj, d = pj * pj + pi * pi, f = 2 * px * (pi + pj), c = (pi + pj) / 2, e = c, g = px;

			c1 = 1 - b - 0.75 * d - 0.5 * f; c4 = 1 - b - 0.5625 * d - 0.25 * f; c2 = c - b + 0.75 * (e - d) + 0.5 * (g - f);
			c5 = c - b + 0.5625 * (e - d) + 0.25 * (g - f); c3 = b + 0.75 * d + 0.5 * f; c6 = b + 0.5625 * d + 0.25 * f;

			double x1 = RELAT(1), x2 = RELAT(0.75), x3 = RELAT(0.5), x4 = RELAT(0);
			X = x1 * b + x2 * d + x3 * f + x4 * (1 - b - d - f);
			X2 = x1 * x1 * b + x2 * x2 * d + x3 * x3 * f + x4 * x4 * (1 - b - d - f);
			w = 1.0 / (X2 - X * X);
		}
		else
		{
			pi = p[vc]; px = 1 - pi;
			double b = pi * pi, d = 2 * pi * px, c = (2 * pi * pi) / (2 * pi), e = (2 * pi * px) / (2 * pi);

			c1 = 1 - b - 0.75 * d; c4 = 1 - b - 0.5625 * d; c2 = c - b + 0.75 * (e - d);
			c5 = c - b + 0.5625 * (e - d); c3 = b + 0.75 * d; c6 = b + 0.5625 * d;

			double x1 = RELAT(1), x2 = RELAT(0.75), x3 = RELAT(0);
			X = x1 * b + x2 * d + x3 * (1 - b - d);
			X2 = x1 * x1 * b + x2 * x2 * d + x3 * x3 * (1 - b - d);
			w = 1.0 / (X2 - X * X);
		}
		double ry = w * RELAT(S);
		if (abs(X) < 1e-10 && IsNormal(w) && IsNormal(ry))
		{
			sry += ry;
			swy += w;
		}
	}
	return 0.5 * (srx / swx + sry / swy);
#undef RELAT
}

/* Initialize Anderson 2007 relatedness estimator */
template<typename REAL>
TARGET void RELATEDNESS<REAL>::R_AndersonInitialize(IND<REAL>* x, IND<REAL>* y)
{
	if (x->vmax > 2 || x->vmin < 2)
		Exit("\nError: Anderson2007 and Milligan2003 relatedness estimators only supports diploids, in individual %s.\n", x->name);
	if (y->vmax > 2 || y->vmin < 2)
		Exit("\nError: Anderson2007 and Milligan2003 relatedness estimators only supports diploids, in individual %s.\n", y->name);

	SetFF(Anderson2007_Coef, nloc * 3);
	for (int64 l = 0; l < nloc; ++l)
	{
		GENOTYPE& gx = x->GetGenotype(l), & gy = y->GetGenotype(l);//fine
		if (gx.Nalleles() == 0 || gy.Nalleles() == 0) continue;
		if (cpop->loc_stat1[l].k <= 1) continue;

		REAL* p = cpop->GetFreq(l);
		int va = gx.GetAlleleCopy(0), vb = gx.GetAlleleCopy(1), vc = gy.GetAlleleCopy(0), vd = gy.GetAlleleCopy(1);
		int ibs = 0;
		double pi = 0, pj = 0, pk = 0, pl = 0;

		if (va == vb) {
			if (vc == vd) {
				if (vc == va) { ibs = 1; pi = p[va]; }
				else { ibs = 2; pi = p[va]; pj = p[vc]; }
			}
			else {
				if (vc == va) { ibs = 3; pi = p[va]; pj = p[vd]; }
				else if (vd == va) { ibs = 3; pi = p[va]; pj = p[vc]; }
				else { ibs = 4; pi = p[va]; pj = p[vc]; pk = p[vd]; }
			}
		}
		else if (vc == vd) {
			if (vc == va) { ibs = 5; pi = p[vc]; pj = p[vb]; }
			else if (vc == vb) { ibs = 5; pi = p[vc]; pj = p[va]; }
			else { ibs = 6; pi = p[vc]; pj = p[va]; pk = p[vb]; }
		}
		else {
			if ((va == vc && vb == vd) || (va == vd && vb == vc))
			{
				ibs = 7; pi = p[va]; pj = p[vb];
			}
			else if (va == vc) { ibs = 8; pi = p[va]; pj = p[vb]; pk = p[vd]; }
			else if (vb == vd) { ibs = 8; pi = p[vb]; pj = p[va]; pk = p[vc]; }
			else if (va == vd) { ibs = 8; pi = p[va]; pj = p[vb]; pk = p[vc]; }
			else if (vb == vc) { ibs = 8; pi = p[vb]; pj = p[va]; pk = p[vd]; }
			else { ibs = 9; pi = p[va]; pj = p[vb]; pk = p[vc]; pl = p[vd]; }
		}

		switch (ibs)
		{
		case 1:  Anderson2007_Coef[l * 3 + 0] = pi * pi;	 Anderson2007_Coef[l * 3 + 1] = pi * pi * pi;			Anderson2007_Coef[l * 3 + 2] = pi * pi * pi * pi; break;
		case 2:  Anderson2007_Coef[l * 3 + 0] = 0;			 Anderson2007_Coef[l * 3 + 1] = 0;						Anderson2007_Coef[l * 3 + 2] = pi * pj * pi * pj; break;
		case 3:  Anderson2007_Coef[l * 3 + 0] = 0;			 Anderson2007_Coef[l * 3 + 1] = pi * pi * pj;			Anderson2007_Coef[l * 3 + 2] = 2 * pi * pi * pi * pj; break;
		case 4:  Anderson2007_Coef[l * 3 + 0] = 0;			 Anderson2007_Coef[l * 3 + 1] = 0;						Anderson2007_Coef[l * 3 + 2] = 2 * pi * pi * pj * pk; break;
		case 5:  Anderson2007_Coef[l * 3 + 0] = 0;			 Anderson2007_Coef[l * 3 + 1] = pi * pi * pj;			Anderson2007_Coef[l * 3 + 2] = 2 * pi * pi * pi * pj; break;
		case 6:  Anderson2007_Coef[l * 3 + 0] = 0;			 Anderson2007_Coef[l * 3 + 1] = 0;						Anderson2007_Coef[l * 3 + 2] = 2 * pi * pi * pj * pk; break;
		case 7:  Anderson2007_Coef[l * 3 + 0] = 2 * pi * pj; Anderson2007_Coef[l * 3 + 1] = pi * pj * (pi + pj);	Anderson2007_Coef[l * 3 + 2] = 4 * pi * pi * pj * pj; break;
		case 8:  Anderson2007_Coef[l * 3 + 0] = 0;			 Anderson2007_Coef[l * 3 + 1] = pi * pj * pk;			Anderson2007_Coef[l * 3 + 2] = 4 * pi * pi * pj * pk; break;
		case 9:  Anderson2007_Coef[l * 3 + 0] = 0;			 Anderson2007_Coef[l * 3 + 1] = 0;						Anderson2007_Coef[l * 3 + 2] = 4 * pi * pj * pk * pl; break;
		default: break;
		}
	}
}

/* Milligan 2003 relatedness estimator */
template<typename REAL>
TARGET double RELATEDNESS<REAL>::R_Milligan2003(IND<REAL>* x, IND<REAL>* y)
{
	return R_Anderson2007(x, y, false);
}

/* Anderson 2007 relatedness estimator */
template<typename REAL>
TARGET double RELATEDNESS<REAL>::R_Anderson2007(IND<REAL>* x, IND<REAL>* y, bool confine)
{
	R_AndersonInitialize(x, y);
	int dim = 2;
	CPOINT xx0 = CPOINT::DownHillSimplex(dim, 0, confine, 0.0001, 15, L_Anderson, NULL);
	xx0.Image2Real();
	return xx0.real[0] + xx0.real[1] / 2;
}

/* Calculate Anderson 2007 likelihood */
template<typename REAL>
TARGET double RELATEDNESS<REAL>::L_Anderson(CPOINT& x, void** unusued)
{
	x.Image2Real();
	double* S = &x.real[0];
	int64 slog = 0; double prod = 1;

	OpenLog(slog, prod);//slog,prod
	for (int64 l = 0; l < nloc; ++l)
	{
		if (*(uint*)(Anderson2007_Coef + l * 3) == 0xFFFFFFFF) continue;
		ChargeLog(slog, prod, SumProd(Anderson2007_Coef + l * 3, S, 3));
	}
	CloseLog(slog, prod);

	x.li = prod < -1e100 ? -1e100 : prod;
	return x.li;
}

/* Ritland 1996 kinship estimator, convert into relatedness */
template<typename REAL>
TARGET double RELATEDNESS<REAL>::R_Ritland1996(IND<REAL>* x, IND<REAL>* y, bool ismodified, bool isrelatedness)
{
	double sr = 0, sw = 0;

	for (int64 l = 0; l < nloc; ++l)
	{
		GENOTYPE& gx = x->GetGenotype(l), & gy = y->GetGenotype(l);//fine
		if (gx.Nalleles() == 0 || gy.Nalleles() == 0) continue;
		LOCSTAT1& stat1 = cpop->loc_stat1[l];
		if (stat1.k <= 1) continue;
		LOCSTAT2<REAL>& stat2 = relatedness_loc_stat[l];

		int vx = gx.Ploidy(), vy = gy.Ploidy();
		int k = GetLoc(l).k;
		REAL* p = cpop->GetFreq(l);
		double tx = -1, ty = -1, txy = -1;

		for (int i = 0; i < k; ++i)
		{
			if (p[i] * stat1.nhaplo <= 1e-5) continue;
			double ax = gx.GetFreq<REAL>(i);
			double ay = gy.GetFreq<REAL>(i);
			double invpi = 1.0 / p[i];
			tx += ax * ax * invpi;
			ty += ay * ay * invpi;
			txy += ax * ay * invpi;
		}

		int minv = Min(vx, vy), maxv = Min(vx, vy);
		double r = 0, w = 0;

		if (isrelatedness)
		{
			if (ismodified)
			{
				//w = 1.0 / stat2.s2;
				//r = minv / (double)(minv + maxv) * (txy / tx + txy / ty) / stat2.s2;

				w = 1.0 / stat2.s2;
				r = (w * minv * txy * (tx + ty)) / ((maxv + minv) * tx * ty);
			}
			else
			{
				r = minv * txy;
				w = stat1.k - 1;
			}
		}
		else
		{
			r = txy;
			w = stat1.k - 1;
		}

		if (IsNormal(r) && IsNormal(w))
		{
			sr += r;
			sw += w;
		}
	}
	return sr / sw;
}

/* Loiselle 1995 kinship estimator, convert into relatedness */
template<typename REAL>
TARGET double RELATEDNESS<REAL>::R_Loiselle1995(IND<REAL>* x, IND<REAL>* y, bool ismodified, bool isrelatedness)
{
	double sr = 0, sw = 0;
	for (int64 l = 0; l < nloc; ++l)
	{
		GENOTYPE& gx = x->GetGenotype(l), & gy = y->GetGenotype(l);//fine
		if (gx.Nalleles() == 0 || gy.Nalleles() == 0) continue;
		LOCSTAT1& stat1 = cpop->loc_stat1[l];
		if (stat1.k <= 1) continue;
		LOCSTAT2<REAL>& stat2 = relatedness_loc_stat[l];

		int k = GetLoc(l).k;
		REAL* p = cpop->GetFreq(l);
		double txy = 0, tb = 0, tx = 0, ty = 0;
		int vx = gx.Ploidy(), vy = gy.Ploidy();

		for (int i = 0; i < k; ++i)
		{
			if (p[i] * stat1.nhaplo <= 1e-5) continue;
			double ax = gx.GetFreq<REAL>(i);
			double ay = gy.GetFreq<REAL>(i);
			txy += (ax - p[i]) * (ay - p[i]);
			tx += (ax - p[i]) * (ax - p[i]);
			ty += (ay - p[i]) * (ay - p[i]);
			tb += p[i] * (1 - p[i]);
		}

		int minv = Min(vx, vy), maxv = Min(vx, vy);
		double r = 0, w = 0;

		if (isrelatedness)
		{
			if (ismodified)
			{
				//w = 1.0 / stat2.s2;
				//r = minv / (double)(minv + maxv) * (txy / tx + txy / ty) / stat2.s2;

				w = 1.0 / stat2.s2;
				r = (w * minv * txy * (tx + ty)) / ((maxv + minv) * tx * ty);
			}
			else
			{
				r = minv * txy;
				w = tb;
			}
		}
		else
		{
			r = txy;
			w = tb;
		}

		if (IsNormal(r) && IsNormal(w))
		{
			sr += r;
			sw += w;
		}
	}
	return sr / sw;
}

/* Weir 1996 kinship estimator, convert into relatedness */
template<typename REAL>
TARGET double RELATEDNESS<REAL>::R_Weir1996(IND<REAL>* x, IND<REAL>* y, bool isrelatedness)
{
	double sr = 0, sw = 0;

	for (int64 l = 0; l < nloc; ++l)
	{
		GENOTYPE& gx = x->GetGenotype(l), & gy = y->GetGenotype(l);//fine
		if (gx.Nalleles() == 0 || gy.Nalleles() == 0) continue;
		LOCSTAT1& stat1 = cpop->loc_stat1[l];
		if (stat1.k <= 1) continue;

		int vx = gx.Ploidy(), vy = gy.Ploidy();
		int k = GetLoc(l).k;
		int minv = Min(vx, vy);
		REAL* p = cpop->GetFreq(l);

		double r = 0, w = 0;
		for (int i = 0; i < k; ++i)
		{
			if (p[i] * stat1.nhaplo <= 1e-5) continue;
			r += (isrelatedness ? minv : 1.0) * (gx.GetFreq<REAL>(i) * gy.GetFreq<REAL>(i) - p[i] * p[i]);
			w += p[i] * p[i];
		}

		if (IsNormal(r) && IsNormal(w))
		{
			sr += r;
			sw += 1 - w;
		}
	}
	return sr / sw;
}

/* Huang 2014 relatedness estimator */
template<typename REAL>
TARGET double RELATEDNESS<REAL>::R_Huang2014(IND<REAL>* x, IND<REAL>* y)
{
	double srx = 0, swx = 0, sry = 0, swy = 0;
	if (x->vmax > 8) Exit("\nError: Huang2014 estimator do not support ploidy level > 8, in individual %s.\n", x->name);
	if (y->vmax > 8) Exit("\nError: Huang2014 estimator do not support ploidy level > 8, in individual %s.\n", y->name);

	for (int64 l = 0; l < nloc; ++l)
	{
		GENOTYPE& gx = x->GetGenotype(l), & gy = y->GetGenotype(l);//fine
		if (gx.Nalleles() == 0 || gy.Nalleles() == 0) continue;
		if (cpop->loc_stat1[l].k <= 1) continue;

		double r = 0, w = 0;
		r = HuangMoment(gx, gy, l, w);
		if (IsNormal(r))
		{
			srx += w * r;
			swx += w;
		}

		r = HuangMoment(gy, gx, l, w);
		if (IsNormal(r))
		{
			sry += w * r;
			swy += w;
		}
	}
	return (srx / swx + sry / swy) * 0.5;
}

/* Huang 2014 relatedness estimator : similarity index */
template<typename REAL>
TARGET double RELATEDNESS<REAL>::S_Index(int* c, int* d, int ploidyx, int ploidyy)
{
	int a[8] = { 0 };
	int b[8] = { 0 };
	SetVal(a, c, ploidyx);
	SetVal(b, d, ploidyy);
	int S = 0;
	for (int i = 0; i < ploidyx; ++i)
	{
		if (a[i] >= 0)
		{
			for (int j = 0; j < ploidyy; ++j)
			{
				if (a[i] == b[j])
				{
					a[i] = b[j] = 0xFFFFFFFF;
					S++;
					break;
				}
			}
		}
	}
	return S * 1.0 / Max(ploidyx, ploidyy);
}

/* Huang 2014 relatedness estimator : get genotype pattern for reference individual */
template<typename REAL>
TARGET int RELATEDNESS<REAL>::GetRefMode(int* a, int ploidy)
{
	int score = 0;
	for (int i = 0; i < ploidy; ++i)
		for (int j = i + 1; j < ploidy; ++j)
			if (a[i] == a[j])
				score++;

	int kind = 1;
	int previous = 0;
	int tkind[8], tl = 0;
	int tkind2[8], tl2 = 0;

	for (int i = 1; i < ploidy; ++i)
		if (a[i] != a[i - 1])
		{
			kind++;
			tkind[tl++] = i - previous;
			tkind2[tl2++] = a[previous];
			previous = i;
		}
	tkind[tl++] = ploidy - previous;
	tkind2[tl2++] = a[previous];

	for (int i = 0; i < kind; ++i)
	{
		for (int j = i + 1; j < kind; ++j)
		{
			if (tkind[i] < tkind[j])
			{
				Swap(tkind[i], tkind[j]);
				Swap(tkind2[i], tkind2[j]);
			}
		}
	}

	int tc = 0;
	for (int i = 0; i < kind; ++i)
		for (int j = 0; j < tkind[i]; ++j)
			a[tc++] = tkind2[i];

	switch (ploidy)
	{
	case 1: switch (score)
	{
	case 0: return 1;
	default: return 0;
	}
	case 2: switch (score)
	{
	case 0: return 1;
	case 1: return 2;
	default: return 0;
	}
	case 3: switch (score)
	{
	case 0: return 1;
	case 1: return 2;
	case 3: return 3;
	default: return 0;
	}
	case 4: switch (score)
	{
	case 0: return 1;
	case 1: return 2;
	case 2: return 3;
	case 3: return 4;
	case 6: return 5;
	default: return 0;
	}
	case 5: switch (score)
	{
	case 0: return 1;
	case 1: return 2;
	case 2: return 3;
	case 3: return 4;
	case 4: return 5;
	case 6: return 6;
	case 10: return 7;
	default: return 0;
	}
	case 6: switch (score)
	{
	case 0: return 1;
	case 1: return 2;
	case 2: return 3;
	case 3: return kind == 3 ? 4 : 5;
	case 4: return 6;
	case 6: return kind == 2 ? 7 : 8;
	case 7: return 9;
	case 10: return 10;
	case 15: return 11;
	default: return 0;
	}
	case 7: switch (score)
	{
	case 0: return 1;
	case 1: return 2;
	case 2: return 3;
	case 3: return kind == 4 ? 4 : 5;
	case 4: return 6;
	case 5: return 7;
	case 6: return kind == 3 ? 8 : 9;
	case 7: return 10;
	case 9: return 11;
	case 10: return 12;
	case 11: return 13;
	case 15: return 14;
	case 21: return 15;
	default: return 0;
	}
	case 8: switch (score)
	{
	case 0: return 1;
	case 1: return 2;
	case 2: return 3;
	case 3: return kind == 5 ? 4 : 6;
	case 4: return kind == 4 ? 5 : 7;
	case 5: return 8;
	case 6: return kind == 4 ? 9 : 11;
	case 7: return kind == 3 ? 10 : 12;
	case 8: return 13;
	case 9: return 14;
	case 10: return 16;
	case 11: return 17;
	case 12: return 15;
	case 13: return 18;
	case 15: return 19;
	case 16: return 20;
	case 21: return 21;
	case 28: return 22;
	default: return 0;
	}
	default: return 0;
	}
}

/* Huang 2014 relatedness estimator : calculate relatedness */
template<typename REAL>
TARGET double RELATEDNESS<REAL>::HuangMoment(GENOTYPE& gx, GENOTYPE& gy, int64 l, double& weight)
{
	weight = 0;
	int vx = gx.Ploidy(), vy = gy.Ploidy();
	int minv = Min(vx, vy), maxv = Max(vx, vy);
	int cp = maxv, cpp = maxv + 1;
	int start = cp - minv;

	int xx[8], yy[8];

	double M[64] = { 0 };
	double A[8] = { 0 };
	double E[8] = { 0 };

	double e[81] = { 0 };
	double P[9] = { 0 };
	double Delta[9] = { 0 };
	double Sol[9] = { 0 };
	double Diff[9] = { 0 };

	double E1 = 0, E2 = 0;
	double tmb[8];

	ushort* gxals = gx.GetAlleleArray(), * gyals = gy.GetAlleleArray();

	for (int j = 0; j < cp; ++j)
	{
		tmb[j] = 1.0;
		xx[j] = gxals[j];
		yy[j] = gyals[j];
	}

	int refmode = 10000 * vx + 100 * vx + GetRefMode(xx, vx);

	E[0] = S_Index(xx, yy, vx, vy);
	for (int j = 1; j < cp; ++j)
		E[j] = E[j - 1] * E[0];

	MOMRelatednessAssign(cp, refmode, e, cpop->GetFreq(l), (int*)xx);

	for (int j1 = 0; j1 < cpp; ++j1)
		for (int j2 = 0; j2 < cp; ++j2)
			e[j1 * cpp + j2] -= e[j1 * cpp + cp];

	for (int j1 = 0; j1 < cp; ++j1)
	{
		for (int j2 = 0; j2 < cp; ++j2)
			tmb[j2] *= (cp - j2) / (double)cp;

		double t1 = 0;
		for (int j2 = 0; j2 < cp; ++j2)
		{
			double t = 0;
			for (int j3 = 0; j3 < cp; ++j3)
				t += e[j3 * cpp + j2] * tmb[j3];
			M[j1 * cp + j2] += t;
			t1 += e[j2 * cpp + cp] * tmb[j2];
		}
		A[j1] += t1;
	}

	for (int j1 = 0; j1 <= cp; ++j1)
		P[j1] += e[j1 * cpp + cp];

	for (int j1 = 0; j1 < minv; ++j1)
		Diff[j1] = E[j1] - A[j1];

	for (int j1 = minv; j1 < cp; ++j1)
	{
		Diff[j1] = 0;
		for (int j2 = 0; j2 < cp; ++j2)
			M[j1 * cp + j2] = M[j2 * cp + cp - 1 - j1] = 0;
	}

	for (int j1 = minv; j1 < cp; ++j1)
		M[j1 * cp + cp - 1 - j1] = 1;

	if (!SolveEquation((double*)M, Diff, Delta, cp))
		return 0;

	for (int j = 0; j < cp; ++j)
		Delta[cp] += Delta[j] * (cp - j) / (double)cp;

	//calc variance
	for (int j1 = start; j1 < cpp; ++j1)
	{
		double s = (cp - j1) / (double)cp;

		for (int j2 = 0; j2 < minv; ++j2)
			Diff[j2] = mp(s, j2 + 1) - A[j2];
		for (int j2 = minv; j2 < cp; ++j2)
			Diff[j2] = 0;

		MatrixMul((double*)M, cp, cp, (double*)Diff, cp, 1, (double*)Sol);

		Sol[cp] = 0;
		for (int j = 0; j < cp; ++j)
			Sol[cp] += Sol[j] * (cp - j) / (double)cp;

		E1 += P[j1] * Sol[cp];
		E2 += P[j1] * Sol[cp] * Sol[cp];
	}

	weight = 1 / (E2 - E1 * E1);
	if (weight > DOUBLE_OVERFLOW) weight = DOUBLE_OVERFLOW;

	return Delta[cp];
}

/* Initialize Huang 2015 relatedness estimator */
template<typename REAL>
TARGET void RELATEDNESS<REAL>::Huang2015_Initialize()
{
#pragma pack(push, 1)
	Huang2015_maps = new TABLE<int, Huang2015ENTRY>[9];
	struct ts {
		uint64 pattern;
		uint hash;
	};
	int len[] = { 0, 2, 9, 31, 109, 339, 1043, 2998, 8405 };
	ts* tdata = (ts*)&mlbin_data[0];//151Kib

	for (int p = 1, count = 0; p <= 8; ++p)
	{
		new(&Huang2015_maps[p]) TABLE<int, Huang2015ENTRY>(false, NULL, len[p]);
		TABLE<int, Huang2015ENTRY>& m = Huang2015_maps[p];
		for (int i = 1; i <= len[p]; ++i)
		{
			Huang2015ENTRY tentry = { i, tdata[count].pattern };
			m[tdata[count++].hash] = tentry;
		}
	}
#pragma pack(pop)
}

/* Uninitialize Huang 2015 relatedness estimator */
template<typename REAL>
TARGET void RELATEDNESS<REAL>::Huang2015_Uninitialize()
{
	delete[] Huang2015_maps;
}

/* Huang 2015 relatedness estimator */
template<typename REAL>
TARGET double RELATEDNESS<REAL>::R_Huang2015(IND<REAL>* x, IND<REAL>* y)
{
	int vx = x->vmax, vy = y->vmax;
	int minv = vx > vy ? vy : vx, maxv = vy > vx ? vy : vx;
	if (vx > 8)
		Exit("\nError: Huang2015 relatedness estimator do not support ploidy level > 8, in individual %s at locus %s.\n", x->name);
	if (vy > 8)
		Exit("\nError: Huang2015 relatedness estimator do not support ploidy level > 8, in individual %s at locus %s.\n", y->name);
	if (x->vmin != x->vmax)
		Exit("\nError: Huang2015 relatedness estimator do not support aneuploids, in individual %s.\n", x->name);
	if (y->vmin != y->vmax)
		Exit("\nError: Huang2015 relatedness estimator do not support aneuploids, in individual %s.\n", y->name);

	for (int64 l = 0; l < nloc; ++l)
	{
		Huang2015_Coef[l * 9] = -999;

		GENOTYPE& gx = x->GetGenotype(l), & gy = y->GetGenotype(l);//fine
		if (gx.Nalleles() == 0 || gy.Nalleles() == 0) continue;
		if (cpop->loc_stat1[l].k <= 1) continue; //monomorphic in current population

		int xx[8] = { 0 }, yy[8] = { 0 }, alleles[16] = { 0 };
		ushort* xals = gx.GetAlleleArray(), * yals = gy.GetAlleleArray();
		for (int i = 0; i < vx; ++i) xx[i] = xals[i];
		for (int i = 0; i < vy; ++i) yy[i] = yals[i];

		Huang2015ENTRY& entry = Huang2015_maps[maxv][GetHuang2015Hash(xx, yy, maxv)];
		if (entry.ibs == 0) continue;

		Huang2015_MatchAllele(entry.pattern, xx, yy, alleles, maxv);
		MLRelatednessAssign(maxv, cpop->GetFreq(l), alleles, Huang2015_Coef + l * 9, entry.ibs);
	}

	int dim = minv, diff = abs(vx - vy);
	bool confine = vx == vy && vx % 2 == 0;
	CPOINT xx0 = CPOINT::DownHillSimplex(dim, diff, confine, 0.0001, 15, RELATEDNESS<REAL>::L_Huang2015, NULL);
	xx0.Image2Real();

	double re = 0;
	for (int i = 0; i < minv; ++i)
		re += xx0.real[i] * (minv - i) / maxv;

	return re;
}

/* Calculate Huang 2015 likelihood */
template<typename REAL>
TARGET double RELATEDNESS<REAL>::L_Huang2015(CPOINT& xx, void** unusued)
{
	xx.Image2Real();

	int64 slog = 0; double prod = 1;
	OpenLog(slog, prod);//slog,prod
	for (int64 l = 0; l < nloc; ++l)
	{
		if (Huang2015_Coef[l * 9] == -999) continue;
		double lt = 0;
		for (int k = 0; k <= xx.dim; ++k)
			lt += Huang2015_Coef[l * 9 + xx.diff + k] * xx.real[k] / BINOMIAL[xx.diff + k][xx.diff];
		ChargeLog(slog, prod, lt);
	}
	CloseLog(slog, prod);

	xx.li = prod < -1e100 ? -1e100 : prod;
	return xx.li;
}

/* Huang 2015 likelihood estimator: Match genotype-pair pattern and assign alleles */
template<typename REAL>
TARGET void RELATEDNESS<REAL>::Huang2015_MatchAllele(int64 pattern, int* gx, int* gy, int* alleles, int p)
{
	for (int i = 0; i < 16; ++i)
		alleles[i] = 0;

	int n = Max((int)pattern & 0xF, (int)(pattern >> (p * 4)) & 0xF) + 1;
	int a1[16], a2[16], a22[16], a22n = 0;
	int cx[8] = { 0 }, cy[8] = { 0 };
	for (int i = p - 1; i >= 0; --i)
	{
		cy[i] = pattern & 0xF;
		pattern >>= 4;
	}
	for (int i = p - 1; i >= 0; --i)
	{
		cx[i] = pattern & 0xF;
		pattern >>= 4;
	}
	for (int i = 0; i < n; ++i)
	{
		int a = 0, b = 0;
		for (int j = 0; j < p; ++j)
		{
			if (cx[j] == i) a++;
			if (cy[j] == i) b++;
		}
		a1[i] = (a << 20) | (b << 16) | i;
	}

	for (int i = 0; i < p + p; ++i)
	{
		int ta = i >= p ? gy[i - p] : gx[i];
		//find is ta in a22
		bool flag = false;
		for (int j = 0; j < a22n; ++j)
			if (a22[j] == ta) flag = true;
		if (flag) continue;
		a22[a22n] = ta;
		int a = 0, b = 0;
		for (int j = 0; j < p; ++j)
		{
			if (gx[j] == ta) a++;
			if (gy[j] == ta) b++;
		}
		a2[a22n++] = (a << 20) | (b << 16) | ta;
	}

	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
		{
			if (a1[j] < a1[i]) Swap(a1[i], a1[j]);
			if (a2[j] < a2[i]) Swap(a2[i], a2[j]);
		}

	for (int i = 0; i < n; ++i)
		alleles[a1[i] & 0xF] = a2[i] & 0xFFFF;
}
#endif

#define extern 
extern TABLE<int, Huang2015ENTRY>* Huang2015_maps;				//Huang2015_maps[ploidylevel][hash] is a entry saves the ibs modex index and genotype pair pattern
extern _thread double* Anderson2007_Coef;						//Anderson 2007 relatedness estimator coefficients
extern _thread double* Huang2015_Coef;							//Huang 2015 relatedness estimator coefficients
extern void* relatedness_buf_;									//Circle buffer for relatedness estimation, NBUF
#define relatedness_buf (*(RELATEDNESS<REAL>**)&relatedness_buf_)
extern void* relatedness_loc_stat_;							//Locus information temporatorily used for cpop
#define relatedness_loc_stat (*(LOCSTAT2<REAL>**)&relatedness_loc_stat_)
#undef extern 

/* Calculate relatedness coefficient */
template<typename REAL>
TARGET void CalcRelatedness()
{
	if (!relatedness) return;
	if (ad) Exit("\nError: relatedness estimation (-relatedness) is incompatible with allelic depth (-ad) option.\n");

	EvaluationBegin();
	if (relatedness_estimator_val[11]) RELATEDNESS<REAL>::Huang2015_Initialize();
	OpenResFile("-relatedness", "Relatedness coefficient");

	bool isfirst = true;
	relatedness_buf = new RELATEDNESS<REAL>[NBUF];
	relatedness_loc_stat = new LOCSTAT2<REAL>[nloc];

	int64 ntot = 0;
	if (relatedness_range_val[1]) for (int i = 0; i < npop; ++i)
		ntot += apops[i]->nind * apops[i]->nind;
	if (relatedness_range_val[2])
		for (int rl = 0; rl < lreg; ++rl)
			for (int i = 0; i < nreg[rl]; ++i)
				ntot += aregs[rl][i]->nind * aregs[rl][i]->nind;
	if (relatedness_range_val[3])
		ntot += total_pop->nind * total_pop->nind;
	ntot <<= 1;

	for (int k = 1; k <= 3; ++k)
	{
		if (relatedness_range_val[k] == 0) continue;
		for (int rl = 0; rl < (k == 2 ? lreg : 1); ++rl)
		{
			int n = k == 1 ? npop : (k == 2 ? nreg[rl] : 1);
			for (int i = 0; i < n; ++i)
			{
				OpenTempFiles(N_RELATEDNESS_ESTIMATOR + 1, ".relatedness");
				cpop = (k == 1 ? apops[i] : (k == 2 ? aregs[rl][i] : total_pop));
				SetZero(relatedness_buf, NBUF);
				cpop->GetLocStat2(relatedness_loc_stat);

				RunThreads(&RelatednessThread<REAL>, &RelatednessGuard1<REAL>, &RelatednessGuard2<REAL>, ntot, cpop->nind * cpop->nind * 2,
					"\nEstimating relatedness coefficient between individuals:\n", g_nthread_val, isfirst);
				isfirst = false;
				JoinTempFiles(N_RELATEDNESS_ESTIMATOR + 1);
			}
		}
	}

	delete[] relatedness_buf;
	delete[] relatedness_loc_stat;
	CloseResFile();
	if (relatedness_estimator_val[11]) RELATEDNESS<REAL>::Huang2015_Uninitialize();

	EvaluationEnd("Relatedness estimation");

	if (relatedness_plot_val == 1)
		RunRscript("relatedness_plot.R");
}

/* Write column format relatedness coefficient results in a guard thread */
THREAD2(RelatednessGuard1)
{
	int64& ii = progress1 = 0;
	int n = cpop->nind;
	if (relatedness_fmt_val[2])
		RELATEDNESS<REAL>::ColumnPrintHeader();

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j, ++ii)
		{
			GUARD_BEGIN2

			RELATEDNESS<REAL>& re = relatedness_buf[ii % NBUF];
			if (j >= i && relatedness_fmt_val[2])
				re.ColumnPrintLine(i, j);

			PROGRESS_VALUE++;

			GUARD_END2
		}
	}
}

/* Write matrix format relatedness coefficient results in a guard thread */
THREAD2(RelatednessGuard2)
{
	int64& ii = progress2 = 0;
	int n = cpop->nind;

	if (relatedness_fmt_val[1])
		for (int k = 1; k <= N_RELATEDNESS_ESTIMATOR; ++k)
			RELATEDNESS<REAL>::MatrixPrintMatrixHeader(k, n);

	for (int i = 0; i < n; ++i)
	{
		if (relatedness_fmt_val[1])
			for (int k = 1; k <= N_RELATEDNESS_ESTIMATOR; ++k)
				RELATEDNESS<REAL>::MatrixPrintRowHeader(k, i);

		for (int j = 0; j < n; ++j, ++ii)
		{
			GUARD_BEGIN2

			RELATEDNESS<REAL>& re = relatedness_buf[ii % NBUF];

			if (relatedness_fmt_val[1])
				for (int k = 1; k <= N_RELATEDNESS_ESTIMATOR; ++k)
					re.MatrixPrintCell(k);

			PROGRESS_VALUE++;

			GUARD_END2
		}
	}
}

/* Calculate relatedness coefficient using multiple threads */
THREAD2(RelatednessThread)
{
	//load ind
	int64 ii = 0;
	int ni = cpop->nind;

	if (relatedness_estimator_val[8] || relatedness_estimator_val[9])  Anderson2007_Coef = new double[nloc * 3];
	if (relatedness_estimator_val[11]) Huang2015_Coef = new double[nloc * 9];

	for (int i = 0; i < ni; ++i)
	{
		for (int j = 0; j < ni; ++j, ++ii)
		{
			THREAD_BEGIN2

			relatedness_buf[ii % NBUF].CalcRelatedness(ainds[i], ainds[j]);

			THREAD_END2
		}
	}

	if (relatedness_estimator_val[8] || relatedness_estimator_val[9])  delete[] Anderson2007_Coef;
	if (relatedness_estimator_val[11]) delete[] Huang2015_Coef;
}
