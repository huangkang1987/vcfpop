/* Analysis of Molecular Variances Functions */

#pragma once
#include "vcfpop.h"

#ifndef _VESSEL_ITERATOR
	/* Go to start */
	TARGET void VESSEL_ITERATOR::Rewind(int nlay)
	{
		SetZero(relative_id, nlay + 1);
		SetZero(universal_id, nlay + 1);
		for (int clay = nlay - 1; clay >= lay; --clay)
			trace[clay] = trace[clay + 1]->subunits[0];
	}

	/* Copy from a reference*/
	TARGET void VESSEL_ITERATOR::Copy(VESSEL_ITERATOR& ref, int nlay)
	{
		SetVal(relative_id, ref.relative_id, nlay + 1);
		SetVal(universal_id, ref.universal_id, nlay + 1);
		SetVal(trace, ref.trace, nlay + 1);
		lay = ref.lay;
	}

	/* Initialize */
	TARGET VESSEL_ITERATOR::VESSEL_ITERATOR()
	{

	}

	/* Initialize */
	TARGET VESSEL_ITERATOR::VESSEL_ITERATOR(int _lay, VESSEL& root, int nlay)
	{
		lay = _lay;
		SetZero(relative_id, nlay + 1);
		SetZero(universal_id, nlay + 1);
		trace[nlay] = &root;
		trace[nlay + 1] = NULL;
		for (int clay = nlay - 1; clay >= lay; --clay)
			trace[clay] = trace[clay + 1]->subunits[0];
	}

	/* Uninitialize */
	TARGET VESSEL_ITERATOR::~VESSEL_ITERATOR()
	{

	}

	/* Go to next vessel */
	TARGET void VESSEL_ITERATOR::Next(int nlay)
	{
		//Higher level vessel have remaining vessels, use rest vessels
		if (trace[lay + 1] && trace[lay + 1]->nsubunits > relative_id[lay] + 1)
		{
			relative_id[lay]++;
			universal_id[lay]++;
			trace[lay] = trace[lay + 1]->subunits[relative_id[lay]];
		}
		else
		{
			//Higher level vessel do not have remaining vessels
			for (int clay = lay; clay < nlay; ++clay)
			{
				if (trace[clay + 1]->nsubunits > relative_id[clay] + 1)
				{
					universal_id[clay]++;
					relative_id[clay]++;
					trace[clay] = trace[clay + 1]->subunits[relative_id[clay]];
					for (int clay2 = clay - 1; clay2 >= lay; --clay2)
						trace[clay2] = trace[clay2 + 1]->subunits[relative_id[clay2]];
					break;
				}
				else
				{
					universal_id[clay]++;
					relative_id[clay] = 0;
					trace[lay] = NULL;
				}
			}
		}
		//return trace[lay];
	}

	/* Get haplotype index to calculate genetic distance */
	TARGET int VESSEL_ITERATOR::GetHapId()
	{
		return trace[lay]->hid;
	}

	/* Get allele to fetch calculate distance */
	TARGET ushort VESSEL_ITERATOR::GetAllele()
	{
		return trace[lay]->allele;
	}

	/* Get subpopulation in print SS */
	TARGET POP* VESSEL_ITERATOR::GetSubpop(int nlay, int tlay)
	{
		POP* tpop = total_pop;
		for (int clay = nlay - 1; clay >= tlay; --clay)
			tpop = tpop->vpop[relative_id[clay]];
		return tpop;
	}

	/* Get individual in print SS */
	TARGET IND* VESSEL_ITERATOR::GetInd(int nlay, int tlay)
{
	POP* tpop = total_pop;
	for (int clay = nlay - 1; clay >= tlay; --clay)
		if (clay > tlay)
			tpop = tpop->vpop[relative_id[clay]];
		else
			return tpop->inds[relative_id[clay]];
	return NULL;
}
#endif

#ifndef _VESSEL
	/* Uninitialize */
	TARGET VESSEL::~VESSEL()
	{

	}

	/* Initialize */
	TARGET VESSEL::VESSEL()
	{
		subunits = NULL;
		nhaplos = NULL;
		allelecount = NULL;
		nsubunits = nhaplo = 0;
		lay = -1;
		hid = -1;
		allele = 0xFFFF;
	}

	/* Deep copy a vessel */
	TARGET VESSEL::VESSEL(VESSEL& r)
	{
		subunits = NULL;
		nhaplos = NULL;
		allelecount = NULL;
		nhaplo = r.nhaplo;
		lay = r.lay;
		hid = r.hid;
		allele = r.allele;
		nsubunits = r.nsubunits;

		if (r.subunits != NULL)
		{
			amova_memory->Alloc(subunits, nsubunits);
			for (int i = 0; i < nsubunits; ++i)
				subunits[i] = new(amova_memory->Alloc(sizeof(VESSEL))) VESSEL(*r.subunits[i]);
		}

		if (r.nhaplos != NULL)
		{
			amova_memory->Alloc(nhaplos, nloc);
			SetVal(nhaplos, r.nhaplos, nloc);
		}

		if (r.allelecount != NULL)
		{
			amova_memory->Alloc(allelecount, KT);
			SetVal(allelecount, r.allelecount, KT);
		}
	}

	/* Create vessel from population */
	TARGET VESSEL::VESSEL(POP* s, int _lay, int& _hid, int64 loc, int method)
	{
		subunits = NULL;
		nhaplos = NULL;
		allelecount = NULL;
		lay = (ushort)_lay;
		nsubunits = 0;
		nhaplo = 0;
		hid = -1;
		allele = 0xFFFF;

		if (!s->ispop)
		{
			if (loc == -1)
			{
				//region homo, ml
				nsubunits = s->npop;
				amova_memory->Alloc(subunits, nsubunits);
				nhaplo = s->nhaplotypes;

				if (method == 4)
				{
					amova_memory->Alloc(nhaplos, nloc);
					SetZero(nhaplos, nloc);
					amova_memory->Alloc(allelecount, KT);
					SetZero(allelecount, KT);
				}

				for (int i = 0; i < nsubunits; ++i)
				{
					subunits[i] = new(amova_memory->Alloc(sizeof(VESSEL))) VESSEL(s->vpop[i], lay - 1, _hid, loc, method);
					if (method == 4)
					{
						Add(nhaplos, subunits[i]->nhaplos, nloc);
						Add(allelecount, subunits[i]->allelecount, KT);
					}
				}
			}
			else
			{
				//region aneu
				for (int i = 0; i < s->npop; ++i)
					if (s->vpop[i]->loc_stat1[loc].nhaplo > 0)
						nsubunits++;

				amova_memory->Alloc(subunits, nsubunits);
				nhaplo = s->loc_stat1[loc].nhaplo;

				for (int i = 0, ic = 0; i < s->npop; ++i)
				{
					if (s->vpop[i]->loc_stat1[loc].nhaplo > 0)
						subunits[ic++] = new(amova_memory->Alloc(sizeof(VESSEL))) VESSEL(s->vpop[i], lay - 1, _hid, loc, method);
					else
						_hid += s->vpop[i]->nhaplotypes;
				}
			}
		}
		else if (amova_cind_val == 1)
		{
			//pop, ind
			if (loc == -1)
			{
				//pop, ind, homo/ml
				nsubunits = s->nind;
				amova_memory->Alloc(subunits, nsubunits);
				nhaplo = s->nhaplotypes;
				if (method < 3) //homo, aneu
					for (int i = 0; i < nsubunits; ++i)
						subunits[i] = new(amova_memory->Alloc(sizeof(VESSEL))) VESSEL(s->inds[i], lay - 1, _hid, loc, method);
				else
				{
					//ml
					amova_memory->Alloc(nhaplos, nloc);
					SetZero(nhaplos, nloc);
					amova_memory->Alloc(allelecount, KT);
					SetZero(allelecount, KT);

					for (int i = 0; i < nsubunits; ++i)
					{
						subunits[i] = new(amova_memory->Alloc(sizeof(VESSEL))) VESSEL(s->inds[i], lay - 1, _hid, loc, method);
						Add(nhaplos, subunits[i]->nhaplos, nloc);
					}

					for (int64 l = 0; l < nloc; ++l)
					{
						int* acount = GetAlleleCount(l);
						ushort* gcount = s->GetGenoCount(l);
						GENOTYPE* gtab = GetLoc(l).GetGtab();
						int ngeno = GetLoc(l).ngeno;

						for (int gi = 0; gi < ngeno; ++gi)
						{
							int gc = gcount[gi];
							if (gc == 0) continue;

							GENOTYPE& gt = gtab[gi];
							if (gt.Nalleles() == 0) continue;

							ushort* als = gt.GetAlleleArray();
							for (int j = 0, vi = gt.Ploidy(); j < vi; ++j)
								acount[als[j]] += gc;
						}
					}
				}
			}
			else
			{
				//pop, ind, aneu
				ushort* gcount = s->GetGenoCount(loc);
				GENOTYPE* gtab = GetLoc(loc).GetGtab();
				int ngeno = GetLoc(loc).ngeno;

				for (int gi = 0; gi < ngeno; ++gi)
				{
					int gc = gcount[gi];
					if (gc == 0) continue;

					GENOTYPE& gt = gtab[gi];
					if (gt.Nalleles() == 0) continue;
					nsubunits += gc;
				}

				amova_memory->Alloc(subunits, nsubunits);
				nhaplo = s->loc_stat1[loc].nhaplo;

				GENO_READER rt(s->ind0id, loc);
				for (int i = 0, ic = 0, iend = s->nind; i < iend; ++i)
				{
					GENOTYPE& gt = gtab[rt.Read()];
					if (gt.Nalleles())
						subunits[ic++] = new(amova_memory->Alloc(sizeof(VESSEL))) VESSEL(s->inds[i], lay - 1, _hid, loc, method);
					else
						_hid += gt.Ploidy();
				}
			}
		}
		else
		{
			//pop, no ind
			if (loc == -1)
			{
				//pop, no ind, homo/ml
				if (method < 3)
				{
					//homo
					nsubunits = nhaplo = s->nhaplotypes;
					amova_memory->Alloc(subunits, nsubunits);
					for (int i = 0, ic = 0, iend = s->nind; i < iend; ++i)
						for (int j = 0, vi = s->inds[i]->vmin; j < vi; ++j)
							subunits[ic++] = new(amova_memory->Alloc(sizeof(VESSEL))) VESSEL(lay - 1, _hid, 0xFFFF);
				}
				else
				{
					//ml
					nhaplo = s->nhaplotypes;
					nsubunits = s->nind;
					amova_memory->Alloc(subunits, nsubunits);

					amova_memory->Alloc(nhaplos, nloc);
					SetZero(nhaplos, nloc);
					amova_memory->Alloc(allelecount, KT);
					SetZero(allelecount, KT);

					for (int i = 0; i < nsubunits; ++i)
					{
						subunits[i] = new(amova_memory->Alloc(sizeof(VESSEL))) VESSEL(s->inds[i], lay - 1, _hid, loc, method);
						Add(nhaplos, subunits[i]->nhaplos, nloc);
					}

					for (int64 l = 0; l < nloc; ++l)
					{
						int* acount = GetAlleleCount(l);
						ushort* gcount = s->GetGenoCount(l);
						GENOTYPE* gtab = GetLoc(l).GetGtab();
						int ngeno = GetLoc(l).ngeno;

						for (int gi = 0; gi < ngeno; ++gi)
						{
							int gc = gcount[gi];
							if (gc == 0) continue;

							GENOTYPE& gt = gtab[gi];
							if (gt.Nalleles() == 0) continue;

							ushort* als = gt.GetAlleleArray();
							for (int j = 0, vi = gt.Ploidy(); j < vi; ++j)
								acount[als[j]] += gc;
						}
					}
				}
			}
			else
			{
				//pop, no ind, aneu
				ushort* gcount = s->GetGenoCount(loc);
				GENOTYPE* gtab = GetLoc(loc).GetGtab();
				int ngeno = GetLoc(loc).ngeno;

				for (int gi = 0; gi < ngeno; ++gi)
				{
					int gc = gcount[gi];
					if (gc == 0) continue;

					GENOTYPE& gt = gtab[gi];
					if (gt.Nalleles() == 0) continue;
					nsubunits += gc * gt.Ploidy();
				}

				amova_memory->Alloc(subunits, nsubunits);
				nhaplo = s->loc_stat1[loc].nhaplo;

				GENO_READER rt(s->ind0id, loc);
				for (int i = 0, ic = 0, iend = s->nind; i < iend; ++i)
				{
					GENOTYPE& gt = gtab[rt.Read()];
					ushort* alleles = gt.GetAlleleArray();
					if (gt.Nalleles()) for (int j = 0, vi = gt.Ploidy(); j < vi; ++j)
						subunits[ic++] = new(amova_memory->Alloc(sizeof(VESSEL))) VESSEL(lay - 1, _hid, alleles[j]);
					else
						_hid += gt.Ploidy();
				}
			}
		}
	}

	/* Create vessel from individual */
	TARGET VESSEL::VESSEL(IND* s, int _lay, int& _hid, int64 loc, int method)
	{
		subunits = NULL;
		nhaplos = NULL;
		allelecount = NULL;
		lay = (ushort)_lay;
		nhaplo = nsubunits = 0;
		hid = -1;
		allele = 0xFFFF;

		if (method < 3)
		{
			if (loc == -1)
			{
				//homo
				nsubunits = nhaplo = s->vmin;
				amova_memory->Alloc(subunits, nhaplo);
				for (int i = 0; i < nhaplo; ++i)
					subunits[i] = new(amova_memory->Alloc(sizeof(VESSEL))) VESSEL(lay - 1, _hid, 0xFFFF);
			}
			else
			{
				//aneu
				GENOTYPE& gt = s->GetGenotype(loc);
				ushort* alleles = gt.GetAlleleArray();
				nsubunits = nhaplo = gt.Ploidy();
				amova_memory->Alloc(subunits, nhaplo);
				for (int i = 0; i < nhaplo; ++i)
					subunits[i] = new(amova_memory->Alloc(sizeof(VESSEL))) VESSEL(lay - 1, _hid, alleles[i]);
			}
		}
		else
		{
			//ml
			amova_memory->Alloc(nhaplos, nloc);
			amova_memory->Alloc(allelecount, KT);
			SetZero(allelecount, KT);

			hid = _hid++;
			nhaplo = s->vmin;
			for (int64 l = 0; l < nloc; ++l)
			{
				int* acount = GetAlleleCount(l);
				GENOTYPE& gt = s->GetGenotype(l);//fine
				nhaplos[l] = gt.Nalleles() ? gt.Ploidy() : 0;
				if (gt.Nalleles())
				{
					ushort* als = gt.GetAlleleArray();
					for (int j = 0, vi = gt.Ploidy(); j < vi; ++j)
						acount[als[j]]++;
				}
			}
		}
	}

	/* Create vessel from haplotype */
	TARGET VESSEL::VESSEL(int _lay, int& _hid, ushort _allele)
	{
		subunits = NULL;
		nhaplos = NULL;
		allelecount = NULL;
		nsubunits = 0;
		nhaplo = 1;
		lay = (ushort)_lay;
		hid = _hid++;
		allele = _allele;
	}

	/* Get the allele count array at locus l*/
	TARGET int* VESSEL::GetAlleleCount(int l)
	{
		return allelecount + allele_freq_offset[l];
	}

	/* Save all vellels in level fa into an array */
	TARGET void VESSEL::GetVessels(VESSEL** vs, int& nvessels, int fa)
	{
		if (fa == lay) for (int i = 0; i < nsubunits; ++i)
			vs[nvessels++] = subunits[i];
		else for (int i = 0; i < nsubunits; ++i)
			subunits[i]->GetVessels(vs, nvessels, fa);
	}

	/* Replace with shuffled vessels */
	TARGET int VESSEL::Replace(VESSEL** vs, int& nvessels, int fa, int method)
	{
		if (fa == lay)
		{
			nhaplo = 0;
			if (method == 4)
			{
				SetZero(nhaplos, nloc);
				SetZero(allelecount, KT);
			}

			for (int i = 0; i < nsubunits; ++i)
			{
				subunits[i] = vs[nvessels++];
				nhaplo += subunits[i]->nhaplo;

				if (method == 4)
				{
					Add(nhaplos, subunits[i]->nhaplos, nloc);
					Add(allelecount, subunits[i]->allelecount, KT);
				}
			}
		}
		else
		{
			nhaplo = 0;
			if (method == 4)
			{
				SetZero(nhaplos, nloc);
				SetZero(allelecount, KT);
			}

			for (int i = 0; i < nsubunits; ++i)
			{
				nhaplo += subunits[i]->Replace(vs, nvessels, fa, method);

				if (method == 4)
				{
					Add(nhaplos, subunits[i]->nhaplos, nloc);
					Add(allelecount, subunits[i]->allelecount, KT);
				}
			}
		}
		return nhaplo;
	}

	/* Shuffle fa level vessels among fb level vessels */
	TARGET void VESSEL::Shuffle(RNG& rng, int fa, int fb, int method, VESSEL** buf)
	{
		if (lay > fb) for (int i = 0; i < nsubunits; ++i)
			subunits[i]->Shuffle(rng, fa, fb, method, buf);
		else if (fb == lay)
		{
			int nvessels = 0;
			//VLA_NEW(vs, VESSEL*, nhaplo);
			GetVessels(buf, nvessels, fa);
			rng.Permute(buf, nvessels);
			nvessels = 0;
			Replace(buf, nvessels, fa, method);
			//VLA_DELETE(vs);
		}
	}

	/* Calculate matrix C for maximum-likelihood method */
	TARGET void VESSEL::GetCML(double* C, int64 l, int* tid, double* tw, int Nh, int nlay, double** W)
	{
		InitC(C, tid, Nh, nlay);
		GetC(tw, tid, C, nlay, W, l);
	}

	/* Initialize matrix C, the coefficient matrix of S = CV */
	TARGET void VESSEL::InitC(double* C, int* tid, int Nh, int nlay)
	{
		SetZero(C, nlay * nlay);
		for (int i = 0; i < nlay; ++i)
			for (int j = 0; j <= i; ++j)
				C[i * nlay + j] = Nh;
		SetZero(tid, nlay + 1);
		return;
	}

	/* Calculate matrix C */
	TARGET void VESSEL::GetC(double* tw, int* tid, double* C, int nlay, double** W, int64 l)
	{
		if (l == -1)
		{
			if (nhaplo == 0) return;
			tw[lay] = 1.0 / nhaplo;
			for (int i = lay; i < nlay; ++i)
				C[i * nlay + lay] -= nhaplo * nhaplo * tw[i + 1];
			for (int i = 0; i < nsubunits; ++i)
				subunits[i]->GetC(tw, tid, C, nlay, W, -1);
			W[lay][tid[lay]++] = tw[lay];
		}
		else
		{
			if (nhaplos[l] == 0) { tid[lay]++; return; }
			tw[lay] = 1.0 / nhaplos[l];
			for (int i = lay; i < nlay; ++i)
				C[i * nlay + lay] -= nhaplos[l] * nhaplos[l] * tw[i + 1];

			if (amova_cind_val == 2)
			{
				if (lay > 1) for (int i = 0; i < nsubunits; ++i)
					subunits[i]->GetC(tw, tid, C, nlay, W, l);
				else for (int i = lay - 1; i < nlay; ++i)
					C[i * nlay + lay - 1] -= nhaplos[l] * tw[i + 1];
			}
			else
			{
				if (subunits != NULL) for (int i = 0; i < nsubunits; ++i)
					subunits[i]->GetC(tw, tid, C, nlay, W, l);
				else for (int i = lay - 1; i < nlay; ++i)
					C[i * nlay + lay - 1] -= nhaplos[l] * tw[i + 1];
			}
			W[lay][tid[lay]++] = tw[lay];
		}
	}

	/* Count number of vessels in each level */
	TARGET void VESSEL::CountVessels(int* count)
	{
		count[lay]++;
		for (int i = 0; i < nsubunits; ++i)
			subunits[i]->CountVessels(count);
	}

	/* Initialize W, W[lay][tid[lay]] = 1 / nhaplo */
	TARGET void VESSEL::InitW(MEMORY& mem, double**& W, int nlay)
	{
		VLA_NEW(count, int, nlay + 1);
		SetZero(count, nlay + 1);
		CountVessels(count);
		mem.Alloc(W, nlay + 1);
		for (int clay = 0; clay <= nlay; ++clay)
			mem.Alloc(W[clay], count[clay]);
		VLA_DELETE(count);
	}

	/* Calculate SS for homoploid method */
	TARGET void VESSEL::GetSSHomo(double* SS, double* gd, int Nh, double** W, int nlay, VESSEL_ITERATOR& ve1, VESSEL_ITERATOR& ve2)
	{
		SetZero(SS, nlay);
		ve1.Rewind(nlay);
		for (int i = 0; i < Nh; ++i)
		{
			double* gd1 = gd + ve1.GetHapId() * Nh;
			ve2.Copy(ve1, nlay); ve2.Next(nlay);
			for (int j = i + 1; j < Nh; ++j)
			{
				double GDt = gd1[ve2.GetHapId()];
				for (int clay = 1; clay <= nlay; ++clay)
					if (ve1.universal_id[clay] == ve2.universal_id[clay])
						SS[clay - 1] += GDt * W[clay][ve1.universal_id[clay]];
				ve2.Next(nlay);
			}
			ve1.Next(nlay);
		}
	}

	/* Calculate SS for aneuploid method */
	TARGET void VESSEL::GetSSAneu(ushort* hap_bucket, double* SS, bool isiam, int nh, int64 l, int k, double* missing0, double** W, int nlay, VESSEL_ITERATOR& ve1, VESSEL_ITERATOR& ve2)
	{
		SetZero(SS, nlay);
		ve1.Rewind(nlay);
		for (int i = 0; i < nh; ++i)
		{
			ve2.Copy(ve1, nlay); ve2.Next(nlay);
			ushort* hi = hap_bucket + ve1.GetHapId() * nloc;
			for (int j = i + 1; j < nh; ++j)
			{
				ushort* hj = hap_bucket + ve2.GetHapId() * nloc;
				double tdist = isiam ? hi[l] != hj[l] : missing0[hi[l] * k + hj[l]];
				if (tdist != 0)
					for (int clay = 1; clay <= nlay; ++clay)
						if (ve1.universal_id[clay] == ve2.universal_id[clay])
							SS[clay - 1] += tdist * W[clay][ve1.universal_id[clay]];
				ve2.Next(nlay);
			}
			ve1.Next(nlay);
		}
	}

	/* Calculate SS for aneuploid method */
	TARGET void VESSEL::GetSSAneu(double* SS, bool isiam, int nh, int k, double* missing0, double** W, int nlay, VESSEL_ITERATOR& ve1, VESSEL_ITERATOR& ve2)
	{
		SetZero(SS, nlay);
		ve1.Rewind(nlay);
		for (int i = 0; i < nh; ++i)
		{
			ve2.Copy(ve1, nlay); ve2.Next(nlay);
			ushort a1 = ve1.GetAllele();
			for (int j = i + 1; j < nh; ++j)
			{
				ushort a2 = ve2.GetAllele();
				double tdist = isiam ? a1 != a2 : missing0[a1 * k + a2];
				if (tdist != 0)
					for (int clay = 1; clay <= nlay; ++clay)
						if (ve1.universal_id[clay] == ve2.universal_id[clay])
							SS[clay - 1] += tdist * W[clay][ve1.universal_id[clay]];
				ve2.Next(nlay);
			}
			ve1.Next(nlay);
		}
	}

	/* Calculate variance component matrix V */
	TARGET void VESSEL::GetV(double* C, double* SS, double*& V, int nlay)
	{
		MatrixInv2(C, nlay);
		MatrixMul(C, nlay, nlay, SS, nlay, 1, V);
	}

	/* Calculate F-statistics */
	TARGET void VESSEL::GetF(double* V, double* F, double* vs, int nlay)
	{
		//Force non-negative variance component
		if (amova_trunc_val == 1)
			Truncate(V, nlay);

		vs[0] = V[0];
		for (int i = 1; i < nlay; ++i)
			vs[i] = vs[i - 1] + V[i];
		for (int i = 0; i < nlay; ++i)
			for (int j = i + 1; j < nlay; ++j)
				F[i * nlay + j] = 1.0 - vs[i] / vs[j];
	}

	/* Calculate F-statistics */
	TARGET void VESSEL::GetF(double* Fi, double* F, int nlay)
{
	if (amova_cind_val == 2)
	{
		for (int i = 0; i < nlay; ++i)
			for (int j = i + 1; j < nlay; ++j)
				F[i * nlay + j] = 1 - (1 - Fi[j]) / (1 - Fi[i]);
	}
	else
	{
		for (int i = 1; i < nlay; ++i)
			F[0 * nlay + i] = Fi[i];
		for (int i = 1; i < nlay; ++i)
			for (int j = i + 1; j < nlay; ++j)
				F[i * nlay + j] = 1 - (1 - Fi[j]) / (1 - Fi[i]);
	}

}
#endif

#ifndef _AMOVA

	/* Initialize */
	TARGET AMOVA::AMOVA()
	{
		Lind = amova_cind_val == 1;
		nlay = Lind + 2 + lreg;
		SSW = new double* [nlay];									nSSW = new int[nlay];
		if (Lind)
		{
			SSW[0] = new double[nind];								SetZero(SSW[0], nSSW[0] = nind);
			SSW[1] = new double[npop];								SetZero(SSW[1], nSSW[1] = npop);
			for (int clay = 2; clay < nlay; ++clay)
			{
				SSW[clay] = new double[nreg[clay - 2]];				SetZero(SSW[clay], nSSW[clay] = nreg[clay - 2]);
			}
		}
		else
		{
			SSW[0] = new double[npop];								SetZero(SSW[0], nSSW[0] = npop);
			for (int clay = 1; clay < nlay; ++clay)
			{
				SSW[clay] = new double[nreg[clay - 1]];				SetZero(SSW[clay], nSSW[clay] = nreg[clay - 1]);
			}
		}

		DF = new double[nlay + 1];									SetZero(DF, nlay + 1);
		V = new double[nlay];										SetZero(V, nlay);
		SS = new double[nlay];										SetZero(SS, nlay);
		F = new double[nlay * nlay]; 								SetZero(F, nlay * nlay);
		EF = new double[nlay * nlay];								SetZero(EF, nlay * nlay);
		EF2 = new double[nlay * nlay];								SetZero(EF2, nlay * nlay);
		G = new int[nlay * nlay];									SetZero(G, nlay * nlay);
		E = new int[nlay * nlay];									SetZero(E, nlay * nlay);
	}

	/* Extract dummy haplotype for homoploid method */
	TARGET void AMOVA::GetHaplotype(ushort* bucket)
	{
		for (int i = 0; i < nind; ++i)
		{
			IND* tind = ainds[i];
			int vi = tind->vmin, vi2 = tind->vmax;
			if (vi != vi2)
				Exit("\nError: Homoploid amova estimator can only be used for homoploids. Error in individual %s.\n", tind->name);
		}

		for (int64 l = 0; l < nloc; ++l)
		{
			GENOTYPE* gtab = GetLoc(l).GetGtab();
			GENO_READER rt(0, l);

			for (int i = 0, ph = 0; i < nind; ++i)
			{
				GENOTYPE& gt = gtab[rt.Read()];
				ushort* als = gt.GetAlleleArray();
				for (int j = 0, vi = gt.Ploidy(); j < vi; ++j, ++ph)
					bucket[ph * nloc + l] = als[j];
			}
		}
	}

	/* Perform AMOVA using homoploid method */
	TARGET void AMOVA::CalcAMOVA_homo()
	{
		//single thread
		MEMORY mem1;		VLA_NEW(Gmem2, MEMORY, g_nthread_val);
		int64 nperm = amova_nperm_val * (amova_test_val == 1);

		method = amova_cmethod_val;
		Lind = amova_cind_val == 1;
		nlay = Lind + 2 + lreg;

		//Dummy haplotype, homoploids
		bool isiam = amova_cmutation_val == 1;
		if (!isiam && abs(g_format_val) <= BCF)
			Exit("\nError: Stepwise mutation model (smm) in AMOVA can only be applied for non-vcf input file, and should use size as allele identifier. \n");

		for (int i = 0; i < nind; ++i)
			if (ainds[i]->vmin != ainds[i]->vmax)
				Exit("\nError: Homoploid AMOVA method do not support aneuploids, in individual %s.\n", ainds[i]->name);

		int Nh = total_pop->nhaplotypes;
		VESSEL** Gpermbuf = new VESSEL * [Nh * g_nthread_val];		SetZero(Gpermbuf, Nh * g_nthread_val);

#define distW(x,y) distw[(int64)(x)*Nh+(y)]
		double* distw = new double[Nh * Nh];						SetZero(distw, Nh * Nh);

		VLA_NEW(Gss, double, nlay * g_nthread_val);					SetZero(Gss, nlay * g_nthread_val);
		VLA_NEW(Gv, double, nlay * g_nthread_val);					SetZero(Gv, nlay * g_nthread_val);
		VLA_NEW(Gvs, double, nlay * g_nthread_val);					SetZero(Gvs, nlay * g_nthread_val);
		VLA_NEW(Gtid, int, (nlay + 1) * g_nthread_val);				SetZero(Gtid, (nlay + 1) * g_nthread_val);
		VLA_NEW(Gtw, double, (nlay + 1) * g_nthread_val);			SetZero(Gtw, (nlay + 1) * g_nthread_val);
		VLA_NEW(Gf, double, nlay * nlay * g_nthread_val);			SetZero(Gf, nlay * nlay * g_nthread_val);
		VLA_NEW(GC, double, nlay * nlay * g_nthread_val);			SetZero(GC, nlay * nlay * g_nthread_val);
		VLA_NEW(Gg, int, nlay * nlay * g_nthread_val);				SetZero(Gg, nlay * nlay * g_nthread_val);
		VLA_NEW(Ge, int, nlay * nlay * g_nthread_val);				SetZero(Ge, nlay * nlay * g_nthread_val);
		VLA_NEW(Gef, double, nlay * nlay * g_nthread_val);			SetZero(Gef, nlay * nlay * g_nthread_val);
		VLA_NEW(Gef2, double, nlay * nlay * g_nthread_val);			SetZero(Gef2, nlay * nlay * g_nthread_val);

		int tref = 0;
		amova_memory = &mem1; amova_memory->ClearMemory();
		VESSEL _vs(total_pop, nlay, tref, -1, 1);

		double** W;  _vs.InitW(mem1, W, nlay);
		_vs.InitC(GC, Gtid, Nh, nlay);
		_vs.GetC(Gtw, Gtid, GC, nlay, W, -1);

		for (int i = 0; i < nlay; ++i)
			DF[i] = round(i > 0 ? (GC[i * nlay + 0] - GC[(i - 1) * nlay + 0]) : (GC[i * nlay + 0]));
		DF[nlay] = GC[(nlay - 1) * nlay + 0];

		double* MISSING2 = NULL, * MISSING1 = NULL, * MISSING0 = NULL;
		MISSING2 = new double[g_nthread_val * npop * npop];
		MISSING1 = new double[g_nthread_val * npop * maxK];
		MISSING0 = isiam ? NULL : new double[g_nthread_val * maxK * maxK];

		//Calculate genetic distance
		int64 progress_sep = (int64)(Binomial(total_pop->nhaplotypes, 2) * 0.1 + 0.5);
#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int64 l = 0; l < nloc; ++l)//
		{
			threadid = omp_get_thread_num();
			SLOCUS loc = GetLoc(l);

			int k = loc.k;
			double* missing2 = MISSING2 + threadid * npop * npop;
			double* missing1 = MISSING1 + threadid * npop * maxK;
			double* missing0 = isiam ? NULL : MISSING0 + threadid * maxK * maxK;
			ushort* alen = isiam ? NULL : loc.GetAlenArray();

			//prepare SMM distance matrix between alleles
			if (!isiam) for (int i = 0; i < k; ++i)
			{
				missing0[i * k + i] = 0;
				for (int j = 0; j < i; ++j)
					missing0[i * k + j] = missing0[j * k + i] = ((int)alen[i] - (int)alen[j]) * ((int)alen[i] - (int)alen[j]);
			}

			SetVal(missing2, 0.0, npop * npop);
			SetVal(missing1, 0.0, npop * maxK);

			GENOTYPE* gtab = loc.GetGtab();
			GENO_READER r1(0, l), r1b;
			for (int i1 = 0, h1 = 0; i1 < nind; ++i1)
			{
				r1b = r1;
				int p1id = ainds[i1]->popid;
				POP* p1 = apops[p1id];
				ushort* als1 = gtab[r1.Read()].GetAlleleArray();

				for (int j1 = 0, v1 = ainds[i1]->vmin; j1 < v1; ++j1, ++h1)
				{
					ushort a1 = als1[j1];
					GENO_READER r2 = r1b;

					for (int i2 = i1, h2 = h1 + 1; i2 < nind; ++i2)
					{
						int p2id = ainds[i2]->popid;
						POP* p2 = apops[p2id];
						ushort* als2 = gtab[r2.Read()].GetAlleleArray();

						for (int j2 = i2 == i1 ? j1 + 1 : 0, v2 = ainds[i2]->vmin; j2 < v2; ++j2, ++h2)
						{
							ushort a2 = als2[j2];

							double tdist = 0;
							if (isiam)
							{
								if (a1 == 0xFFFF && a2 == 0xFFFF)	tdist = (missing2[p1id * npop + p2id] == 0 ? (missing2[p1id * npop + p2id] = 1 - SumProd(p1->GetFreq(l), p2->GetFreq(l), k)) : missing2[p1id * npop + p2id]);
								else if (a1 == 0xFFFF)				tdist = (missing1[a2 * npop + p1id] == 0 ? (missing1[a2 * npop + p1id] = 1 - p1->GetFreq(l, a2)) : missing1[a2 * npop + p1id]);
								else if (a2 == 0xFFFF)				tdist = (missing1[a1 * npop + p2id] == 0 ? (missing1[a1 * npop + p2id] = 1 - p2->GetFreq(l, a1)) : missing1[a1 * npop + p2id]);
								else if (a1 != a2)					tdist = 1;
							}
							else
							{
								if (a1 == 0xFFFF && a2 == 0xFFFF)	tdist = (missing2[p1id * npop + p2id] == 0 ? (missing2[p1id * npop + p2id] = SumProdSMM(alen, p1->GetFreq(l), p2->GetFreq(l), k)) : missing2[p1id * npop + p2id]);
								else if (a1 == 0xFFFF)				tdist = (missing1[a2 * npop + p1id] == 0 ? (missing1[a2 * npop + p1id] = SumProdSMM(alen, p1->GetFreq(l), a2, k)) : missing1[a2 * npop + p1id]);
								else if (a2 == 0xFFFF)				tdist = (missing1[a1 * npop + p2id] == 0 ? (missing1[a1 * npop + p2id] = SumProdSMM(alen, p2->GetFreq(l), a1, k)) : missing1[a1 * npop + p2id]);
								else if (a1 != a2)					tdist = missing0[a1 * k + a2];
							}

							if (tdist != 0)  AtomicAddD(distW(h1, h2), tdist);
						}
					}
				}
			}
			PROGRESS_VALUE += progress_sep;
		}

		delete[] MISSING2;
		delete[] MISSING1;
		if (MISSING0) delete[] MISSING0;

		//Calculate SS within each vessel
#pragma omp parallel  num_threads(g_nthread_val)
		{
			threadid = omp_get_thread_num();
			VESSEL_ITERATOR ve1(0, _vs, nlay), ve2(0, _vs, nlay);
			for (int i = 0; i < Nh; ++i)
			{
				if (i % g_nthread_val != threadid) { ve1.Next(nlay); continue; }

				int hid1 = ve1.GetHapId();
				ve2.Copy(ve1, nlay);  ve2.Next(nlay);

				for (int j = i + 1; j < Nh; ++j)
				{
					int hid2 = ve2.GetHapId();
					double tdist = distW(hid2, hid1) = distW(hid1, hid2);
					if (tdist != 0)
						for (int clay = 1; clay <= nlay; ++clay)
							if (ve1.universal_id[clay] == ve2.universal_id[clay])
								AtomicAddD(SSW[clay - 1][ve1.universal_id[clay]], tdist * W[clay][ve1.universal_id[clay]]);

					ve2.Next(nlay);
				}
				ve1.Next(nlay);

				PROGRESS_VALUE += Nh - i - 1;
			}
		}

		//Calculate initial variance components and F-statistics
		{
			VESSEL_ITERATOR ve1(0, _vs, nlay), ve2(0, _vs, nlay);
			_vs.InitC(GC, Gtid, Nh, nlay);
			_vs.GetC(Gtw, Gtid, GC, nlay, W, -1);
			_vs.GetSSHomo(Gss, distw, Nh, W, nlay, ve1, ve2);

			VESSEL::GetV(GC, Gss, Gv, nlay);
			VESSEL::GetF(Gv, Gf, Gvs, nlay);

			SetVal(SS, Gss, nlay);
			SetVal(V, Gv, nlay);
			SetVal(F, Gf, nlay * nlay);
		}

		for (int fa = 1; fa <= nlay; ++fa)
			for (int fb = fa + 1; fb <= nlay; ++fb)
			{
				int pairid = (fa - 1) * nlay + (fb - 1);
				VLA_NEW(_vs2, VESSEL, g_nthread_val);

				for (int i = 0; i < g_nthread_val; ++i)
				{
					amova_memory = &Gmem2[i];	amova_memory->ClearMemory();
					new(&_vs2[i]) VESSEL(_vs);
				}

				//Begin permutations
#pragma omp parallel  for num_threads(g_nthread_val)  schedule(static, 1)
				for (int64 m = 0; m < nperm; ++m)
				{
					threadid = omp_get_thread_num();
					double* ss = Gss + nlay * threadid;
					double* v = Gv + nlay * threadid;
					double* vs = Gvs + nlay * threadid;
					int* tid = Gtid + (nlay + 1) * threadid;
					double* tw = Gtw + (nlay + 1) * threadid;
					double* f = Gf + nlay * nlay * threadid;
					double* C = GC + nlay * nlay * threadid;
					int* g = Gg + nlay * nlay * threadid;
					int* e = Ge + nlay * nlay * threadid;
					double* ef = Gef + nlay * nlay * threadid;
					double* ef2 = Gef2 + nlay * nlay * threadid;
					VESSEL** permbuf = Gpermbuf + Nh * threadid;

					RNG rng(g_seed_val + pairid * nperm + m);
					amova_memory = &Gmem2[threadid];	amova_memory->ClearMemory();
					VESSEL& vs2 = _vs2[threadid];
					VESSEL_ITERATOR ve1(0, vs2, nlay), ve2(0, vs2, nlay);

					vs2.Shuffle(rng, fa, fb, 1, permbuf);
					vs2.InitC(C, tid, Nh, nlay);
					vs2.GetC(tw, tid, C, nlay, W, -1);
					vs2.GetSSHomo(ss, distw, Nh, W, nlay, ve1, ve2);

					VESSEL::GetV(C, ss, v, nlay);
					VESSEL::GetF(v, f, vs, nlay);

					if (f[pairid] > F[pairid] + 1e-7) g[pairid]++;
					else if (f[pairid] > F[pairid] - 1e-7) e[pairid]++;
					ef[pairid] += f[pairid];
					ef2[pairid] += f[pairid] * f[pairid];

					PROGRESS_VALUE += 1000;
				}

				//Sum results of multiple threads
				for (int i = 0; i < g_nthread_val; ++i)
				{
					int* g = Gg + nlay * nlay * i;
					int* e = Ge + nlay * nlay * i;
					double* ef = Gef + nlay * nlay * i;
					double* ef2 = Gef2 + nlay * nlay * i;

					G[pairid] += g[pairid];
					E[pairid] += e[pairid];
					EF[pairid] += ef[pairid];
					EF2[pairid] += ef2[pairid];
				}

				VLA_DELETE(_vs2);
			}

		VLA_DELETE(Gss);
		VLA_DELETE(Gv);
		VLA_DELETE(Gvs);
		VLA_DELETE(Gtid);
		VLA_DELETE(Gtw);
		VLA_DELETE(Gf);
		VLA_DELETE(GC);
		VLA_DELETE(Gg);
		VLA_DELETE(Ge);
		VLA_DELETE(Gef);
		VLA_DELETE(Gef2);
		VLA_DELETE(Gmem2);
		delete[] Gpermbuf;
		delete[] distw;
	}

	/* Perform AMOVA using aneuploid method */
	TARGET void AMOVA::CalcAMOVA_aneu()
	{
		//Aneuploids, sum SS and C across locus
		method = amova_cmethod_val;
		Lind = amova_cind_val == 1;
		nlay = Lind + 2 + lreg;

		bool isiam = amova_cmutation_val == 1;
		if (!isiam && abs(g_format_val) <= BCF)
			Exit("\nError: Stepwise mutation model (smm) in AMOVA can only be applied for non-vcf input file, and should use size as allele identifier. \n");

		//Number of pseudo permuations at each locus
		int64 M = amova_pseudo_val * (amova_test_val == 1);
		int64 nperm = amova_nperm_val * (amova_test_val == 1);
		bool PseudoPerm = amova_pseudo_val > 0;
		if (!PseudoPerm) M = nperm;

		int npair = BINOMIAL[nlay][2];
		double* SSL, * CL;

		//Choose the approach with the least memory

		//Approach 1: CacheLocus, cache results (SS,C) for each locus
		//and each real permute will randomly select one result out of M 
		//results at each locus and take their sum
#define SSLWLocus(lid,pair,mid) SSL[((lid) * npair * M + (pair) * M + (mid)) * nlay]
#define  CLWLocus(lid,pair,mid)  CL[((lid) * npair * M + (pair) * M + (mid)) * nlay * nlay]
		bool CacheLocus = PseudoPerm && (nloc * M < nperm);
		if (CacheLocus)
		{
			//[L][npair][M]
			//l * npair * M + pairid * M
			int npt = nloc * npair * M * nlay;
			SSL = new double[npt];										SetZero(SSL, npt);
			CL = new double[npt * nlay];								SetZero(CL, npt * nlay);
		}

		//Approach 2: CachePerm, cache results (SS,C) for each permute
		//each locus have M results, randomly distriubte to amova_nperm_val permutes
#define SSLWPerm(pair,perm) SSL[((pair) * nperm + (perm)) * nlay]
#define  CLWPerm(pair,perm)  CL[((pair) * nperm + (perm)) * nlay * nlay]
		bool CachePerm = !CacheLocus;
		if (CachePerm)
		{
			//[npair][nperm]
			//pairid * nperm
			int npt = npair * nperm * nlay;
			SSL = new double[npt];										SetZero(SSL, npt);
			CL = new double[npt * nlay];								SetZero(CL, npt * nlay);
		}

		//Allocate single thread memory
		VLA_NEW(C, double, nlay * nlay);								SetZero(C, nlay * nlay);

		//Allocate memory for each thread
		int Nht = total_pop->nhaplotypes;
		VESSEL** Gpermbuf = new VESSEL * [Nht * g_nthread_val];			SetZero(Gpermbuf, Nht * g_nthread_val);

		VLA_NEW(Gss, double, nlay * g_nthread_val);						SetZero(Gss, nlay * g_nthread_val);
		VLA_NEW(Gv, double, nlay * g_nthread_val);						SetZero(Gv, nlay * g_nthread_val);
		VLA_NEW(Gvs, double, nlay * g_nthread_val);						SetZero(Gvs, nlay * g_nthread_val);
		VLA_NEW(Gf, double, nlay * nlay * g_nthread_val);				SetZero(Gf, nlay * nlay * g_nthread_val);
		VLA_NEW(Gc, double, nlay * nlay * g_nthread_val);				SetZero(Gc, nlay * nlay * g_nthread_val);
		VLA_NEW(Gtid, int, (nlay + 1) * g_nthread_val);					SetZero(Gtid, (nlay + 1) * g_nthread_val);
		VLA_NEW(Gtw, double, (nlay + 1) * g_nthread_val);				SetZero(Gtw, (nlay + 1) * g_nthread_val);
		VLA_NEW(Gg, int, nlay * nlay * g_nthread_val);					SetZero(Gg, nlay * nlay * g_nthread_val);
		VLA_NEW(Ge, int, nlay * nlay * g_nthread_val);					SetZero(Ge, nlay * nlay * g_nthread_val);
		VLA_NEW(Gef, double, nlay * nlay * g_nthread_val);				SetZero(Gef, nlay * nlay * g_nthread_val);
		VLA_NEW(Gef2, double, nlay * nlay * g_nthread_val);				SetZero(Gef2, nlay * nlay * g_nthread_val);
		VLA_NEW(Gmem1, MEMORY, g_nthread_val);
		double* Gmss = NULL, * Gmc = NULL;

		if (CachePerm && PseudoPerm)
		{
			Gmss = new double[nlay * M * g_nthread_val];				SetZero(Gmss, nlay * M * g_nthread_val);
			Gmc = new double[nlay * nlay * M * g_nthread_val];			SetZero(Gmc, nlay * nlay * M * g_nthread_val);
		}

		double* MISSING0 = isiam ? NULL : new double[g_nthread_val * maxK * maxK];

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int64 l = 0; l < nloc; ++l)
		{
			threadid = omp_get_thread_num();
			SLOCUS loc = GetLoc(l);
			int k = loc.k;
			double* missing0 = isiam ? NULL : MISSING0 + threadid * maxK * maxK;
			ushort* alen = isiam ? NULL : loc.GetAlenArray();

			//prepare SMM distance matrix between alleles
			if (!isiam) for (int i = 0; i < k; ++i)
			{
				missing0[i * k + i] = 0;
				for (int j = 0; j < i; ++j)
					missing0[i * k + j] = missing0[j * k + i] = ((int)alen[i] - (int)alen[j]) * ((int)alen[i] - (int)alen[j]);
			}

			int nh = total_pop->loc_stat1[l].nhaplo;
			MEMORY& mem1 = Gmem1[threadid];
			amova_memory = &mem1; amova_memory->ClearMemory();

			double* ss = Gss + nlay * threadid;
			//double* v   = Gv   + nlay * threadid;
			//double* vs  = Gvs  + nlay * threadid;
			//double* f   = Gf   + nlay * nlay * threadid;
			double* c = Gc + nlay * nlay * threadid;
			int* tid = Gtid + (nlay + 1) * threadid;
			double* tw = Gtw + (nlay + 1) * threadid;
			double* mss = Gmss + nlay * M * threadid;
			double* mc = Gmc + nlay * nlay * M * threadid;
			VESSEL** permbuf = Gpermbuf + Nht * threadid;

			int tref = 0;
			VESSEL _vs(total_pop, nlay, tref, l, 1);
			double** W;  _vs.InitW(mem1, W, nlay);

			_vs.InitC(c, tid, nh, nlay);
			_vs.GetC(tw, tid, c, nlay, W, -1);

			for (int i = 0; i < nlay; ++i)
				AtomicAddD(DF[i], round(i > 0 ? (c[i * nlay + 0] - c[(i - 1) * nlay + 0]) : (c[i * nlay + 0])));
			AtomicAddD(DF[nlay], c[(nlay - 1) * nlay + 0]);

			//Calculate SS within each vessel
			{
				VESSEL_ITERATOR ve1(0, _vs, nlay), ve2(0, _vs, nlay);
				GENOTYPE* gtab = GetLoc(l).GetGtab();
				GENO_READER r1(0, l), r1b;

				for (int i1 = 0, h1 = 0; i1 < nind; ++i1)
				{
					r1b = r1;
					//POP* p1 = apops[ainds[i1]->popid];
					GENOTYPE& g1 = gtab[r1.Read()];
					ushort* als1 = g1.GetAlleleArray();

					for (int j1 = 0, v1 = g1.Ploidy(), na1 = g1.Nalleles(); j1 < v1; ++j1, ++h1)
					{
						if (na1 == 0) continue;

						ushort a1 = als1[j1];
						GENO_READER r2 = r1b;
						ve2.Copy(ve1, nlay); ve2.Next(nlay);

						for (int i2 = i1, h2 = h1 + 1; i2 < nind; ++i2)
						{
							//POP* p2 = apops[ainds[i2]->popid];
							GENOTYPE& g2 = gtab[r2.Read()];
							ushort* als2 = g2.GetAlleleArray();

							for (int j2 = i2 == i1 ? j1 + 1 : 0, v2 = g2.Ploidy(), na2 = g2.Nalleles(); j2 < v2; ++j2, ++h2)
							{
								if (na2 == 0) continue;
								ushort a2 = als2[j2];
								double tdist = isiam ? a1 != a2 : missing0[a1 * k + a2];
								if (tdist != 0)
									for (int clay = 1; clay <= nlay; ++clay)
										if (ve1.universal_id[clay] == ve2.universal_id[clay])
											AtomicAddD(SSW[clay - 1][ve1.universal_id[clay]], tdist * W[clay][ve1.universal_id[clay]]);

								ve2.Next(nlay);
							}
						}
						ve1.Next(nlay);
					}
				}
			}

			//Calculate variance components and F-statistics
			{
				VESSEL& vs1 = _vs;
				VESSEL_ITERATOR ve1(0, vs1, nlay), ve2(0, vs1, nlay);

				vs1.InitC(c, tid, nh, nlay);
				vs1.GetC(tw, tid, c, nlay, W, -1);
				vs1.GetSSAneu(ss, isiam, nh, k, missing0, W, nlay, ve1, ve2);

				AtomicAddD(SS, ss, nlay);
				AtomicAddD(C, c, nlay * nlay);
			}

			//Calculate SS and C of each pseudo permutation at this locus
			for (int fa = 1, pairid = 0; fa <= nlay; ++fa)
			{
				for (int fb = fa + 1; fb <= nlay; ++fb, ++pairid)
				{
					VESSEL vs2(_vs);
					VESSEL_ITERATOR ve1(0, vs2, nlay), ve2(0, vs2, nlay);

					double* ssl, * cl;

					if (CacheLocus)
					{
						ssl = &SSLWLocus(l, pairid, 0);	 //SSL + (l * npair * M + pairid * M) * nlay;
						cl = &CLWLocus(l, pairid, 0);    // CL + (l * npair * M + pairid * M) * nlay * nlay;
					}

					if (CachePerm)
					{
						ssl = &SSLWPerm(pairid, 0);		 //SSL + (pairid * nperm) * nlay;
						cl = &CLWPerm(pairid, 0);		 // CL + (pairid * nperm) * nlay * nlay;
					}

					for (int64 m2 = 0; m2 < M; ++m2)
					{
						RNG rng(g_seed_val + pairid * nloc * M + l * M + m2);

						vs2.Shuffle(rng, fa, fb, 2, permbuf);
						vs2.InitC(c, tid, nh, nlay);
						vs2.GetC(tw, tid, c, nlay, W, -1);
						vs2.GetSSAneu(ss, isiam, nh, k, missing0, W, nlay, ve1, ve2);

						//Save results for M pseudo permuations and will be add to real permuations after all loci are finished
						if (CacheLocus)
						{
							//ss matrix of (m2+1)-th permutation at locus l
							SetVal(ssl + nlay * m2, ss, nlay);//locus specific, do not need atomic operations
							SetVal(cl  + nlay * nlay * m2, c, nlay * nlay);
						}

						//Save results for M pseudo permuations and will be distributed to nperm real permuations
						if (CachePerm && PseudoPerm)
						{
							//ss matrix of (m2+1)-th pseudo-permutation at locus l
							SetVal(mss + nlay * m2, ss, nlay);//thread specific, do not need atomic operations
							SetVal(mc  + nlay * nlay * m2, c, nlay * nlay);
						}

						//Save results for nperm real permuations
						if (CachePerm && !PseudoPerm)
						{
							//the ss matrix of a real permuation is the sum of ss and c matrices across L locui
							AtomicAddD(ssl + nlay * m2, ss, nlay);
							AtomicAddD(cl  + nlay * nlay * m2, c, nlay * nlay);
						}

						PROGRESS_VALUE += 500;
					}

					//Distribute M pseudo perms into nperm real perms
					if (CachePerm && PseudoPerm)
					{
						for (int64 m = 0; m < nperm; ++m)
						{
							//the ss matrix of a real permuation is the sum of ss and c matrices across L locui
							uint64 m2 = (Hash64ULong(g_seed_val) ^ Hash64ULong(pairid * nloc * nperm + l * nperm + m)) % M;
							AtomicAddD(ssl + nlay * m, mss + nlay * m2, nlay);
							AtomicAddD(cl  + nlay * nlay * m, mc + nlay * nlay * m2, nlay * nlay);
						}
					}
				}
			}
		}

		if (MISSING0) delete[] MISSING0;
		
		/*******************************************************************/

		//Sum SS and C across loci for each permuation
		VESSEL::GetV(C, SS, V, nlay);
		VESSEL::GetF(V, F, Gvs, nlay);

		VLA_NEW(GpermSS, double, nlay * g_nthread_val);
		VLA_NEW(GpermC, double, nlay * nlay * g_nthread_val);

		for (int fa = 1, pairid = 0; fa <= nlay; ++fa)
		{
			for (int fb = fa + 1; fb <= nlay; ++fb, ++pairid)
			{
				int idx = (fa - 1) * nlay + (fb - 1);
				double* clPerm = &CLWPerm(pairid, 0);		// CL + pairid * (nperm * nlay * nlay);
				double* sslPerm = &SSLWPerm(pairid, 0);		//SSL + pairid * (nperm * nlay);

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
				for (int64 m = 0; m < nperm; ++m)
				{
					threadid = omp_get_thread_num();

					double* v = Gv + nlay * threadid;
					double* f = Gf + nlay * nlay * threadid;
					double* vs = Gvs + nlay * threadid;
					int* g = Gg + nlay * nlay * threadid;
					int* e = Ge + nlay * nlay * threadid;
					double* ef = Gef + nlay * nlay * threadid;
					double* ef2 = Gef2 + nlay * nlay * threadid;

					if (CacheLocus)
					{
						double* permSS = GpermSS + nlay * threadid;
						double* permC = GpermC + nlay * nlay * threadid;
						SetZero(permSS, nlay);
						SetZero(permC, nlay * nlay);

						//Random sample one SS and C out of M results at each locus and sum them across loci
						for (int64 l = 0; l < nloc; ++l)
						{
							uint64 m2 = (Hash64ULong(g_seed_val) ^ Hash64ULong(pairid * nloc * nperm + l * nperm + m)) % M;
							double* ssl = &SSLWLocus(l, pairid, 0);//SSL + (l * npair * M + pairid * M) * nlay;
							double* cl = &CLWLocus(l, pairid, 0);// CL + (l * npair * M + pairid * M) * nlay * nlay;
							Add(permSS, ssl + nlay * m2, nlay);    //thread specific, do not need atomic operations
							Add(permC, cl + nlay * nlay * m2, nlay * nlay);
						}

						VESSEL::GetV(permC, permSS, v, nlay);
					}

					if (CachePerm)
						VESSEL::GetV(clPerm + nlay * nlay * m, sslPerm + nlay * m, v, nlay);

					VESSEL::GetF(v, f, vs, nlay);

					if (f[idx] > F[idx] + 1e-7)  g[idx]++;
					else if (f[idx] > F[idx] - 1e-7)  e[idx]++;

					ef[idx] += f[idx];
					ef2[idx] += f[idx] * f[idx];

					PROGRESS_VALUE += 10;
				}

				//Sum results of multiple threads
				for (int i = 0; i < g_nthread_val; ++i)
				{
					int* g = Gg + nlay * nlay * i;
					int* e = Ge + nlay * nlay * i;
					double* ef = Gef + nlay * nlay * i;
					double* ef2 = Gef2 + nlay * nlay * i;

					G[idx] += g[idx];
					E[idx] += e[idx];
					EF[idx] += ef[idx];
					EF2[idx] += ef2[idx];
				}
			}
		}

		VLA_DELETE(GpermSS);
		VLA_DELETE(GpermC);
		VLA_DELETE(Gss);
		VLA_DELETE(Gv);
		VLA_DELETE(Gvs);
		VLA_DELETE(Gf);
		VLA_DELETE(Gc);
		VLA_DELETE(Gtid);
		VLA_DELETE(Gtw);
		VLA_DELETE(Gg);
		VLA_DELETE(Ge);
		VLA_DELETE(Gef);
		VLA_DELETE(Gef2);
		VLA_DELETE(C);
		VLA_DELETE(Gmem1);
		delete[] Gpermbuf;
		delete[] SSL;
		delete[] CL;

		if (CachePerm && PseudoPerm)
		{
			delete[] Gmss;
			delete[] Gmc;
		}
	}

	/* Calculate likelihood for permuated data */
	TARGET double AMOVA::Likelihood(CPOINT& xx, void** Param)
	{
		int clay = *(int*)Param[0];
		int nlay = *(int*)Param[1];
		VESSEL_ITERATOR& ve = *(VESSEL_ITERATOR*)Param[2];

		xx.Image2RealSelfing();
		ve.Rewind(nlay);
		int64 slog = 0; double prod = 1;
		double f = xx.real[0];
		OpenLog(slog, prod);//slog,prod

		for (int i = 0; i < nind; ++i)
		{
			IND* ind = ainds[ve.GetHapId()];
			for (int64 l = 0; l < nloc; ++l)
			{
				GENOTYPE& gt = ind->GetGenotype(l);//fine
				if (gt.Nalleles())
				{
					double gfz = gt.GFZ(ve.trace[clay + 1]->GetAlleleCount(l), ve.trace[clay + 1]->nhaplos[l], f);
					ChargeLog(slog, prod, gfz);
				}
			}
			ve.Next(nlay);
		}

		CloseLog(slog, prod);
		return prod;
	}

	/* Perform AMOVA using maximum-likelihood method */
	TARGET void AMOVA::CalcAMOVA_ml()
	{
		MEMORY mem1, mem2;

		method = amova_cmethod_val;
		Lind = amova_cind_val == 1;
		nlay = Lind + 2 + lreg;

		//Dummy haplotype, homoploids
		bool isiam = amova_cmutation_val == 1;
		if (!isiam && abs(g_format_val) <= BCF)
			Exit("\nError: Stepwise mutation model (smm) in AMOVA can only be applied for non-vcf input file, and should use size as allele identifier. \n");

		for (int i = 0; i < nind; ++i)
			if (ainds[i]->vmin != ainds[i]->vmax)
				Exit("\nError: Likelihood AMOVA method do not support aneuploids, in individual %s.\n", ainds[i]->name);

		int tref = 0;
		amova_memory = &mem1; amova_memory->ClearMemory();
		VESSEL _vs(total_pop, nlay, tref, -1, 4);
		double** W;  _vs.InitW(mem1, W, nlay);

		int Nht = total_pop->nhaplotypes;
		VESSEL** Gpermbuf = new VESSEL * [Nht * g_nthread_val];		SetZero(Gpermbuf, Nht * g_nthread_val);

		VLA_NEW(Fi, double, nlay);									SetZero(Fi, nlay);
		VLA_NEW(C, double, nlay * nlay);							SetZero(C, nlay * nlay);
		VLA_NEW(Gfi, double, nlay * g_nthread_val);					SetZero(Gfi, nlay * g_nthread_val);
		VLA_NEW(Gf, double, nlay * nlay * g_nthread_val);			SetZero(Gf, nlay * nlay * g_nthread_val);
		VLA_NEW(Gc, double, nlay * nlay * g_nthread_val);			SetZero(Gc, nlay * nlay * g_nthread_val);
		VLA_NEW(Gc2, double, nlay * nlay * g_nthread_val);			SetZero(Gc2, nlay * nlay * g_nthread_val);
		VLA_NEW(Gtid, int, (nlay + 1) * g_nthread_val);				SetZero(Gtid, (nlay + 1) * g_nthread_val);
		VLA_NEW(Gtw, double, (nlay + 1) * g_nthread_val);			SetZero(Gtw, (nlay + 1) * g_nthread_val);
		VLA_NEW(Gg, int, nlay * nlay * g_nthread_val);				SetZero(Gg, nlay * nlay * g_nthread_val);
		VLA_NEW(Ge, int, nlay * nlay * g_nthread_val);				SetZero(Ge, nlay * nlay * g_nthread_val);
		VLA_NEW(Gef, double, nlay * nlay * g_nthread_val);			SetZero(Gef, nlay * nlay * g_nthread_val);
		VLA_NEW(Gef2, double, nlay * nlay * g_nthread_val);			SetZero(Gef2, nlay * nlay * g_nthread_val);

		//Calculate original F-statistics
		for (int clay = Lind; clay < nlay; ++clay)
		{
			VESSEL_ITERATOR ve(amova_cind_val == 2 ? 0 : 1, _vs, nlay);
			void* Param[] = { (void*)&clay, (void*)&nlay, (void*)&ve };
			CPOINT xx = CPOINT::DownHillSimplex(1, 0, false, 0.1, 10, AMOVA::Likelihood, Param);
			xx.Image2RealSelfing();
			Fi[clay] = xx.real[0];
		}
		VESSEL::GetF(Fi, F, nlay);

		//int nthread = amova_test_val == 1 ? g_nthread_val : 1;
		int64 nperm = amova_nperm_val * (amova_test_val == 1);
		VLA_NEW(_vs2, VESSEL, g_nthread_val);

		for (int fa = 1 + Lind; fa <= nlay; ++fa)
			for (int fb = fa + 1; fb <= nlay; ++fb)
			{
				int pairid = (fa - 1) * nlay + (fb - 1);

				for (int i = 0; i < g_nthread_val; ++i)
					new(&_vs2[i]) VESSEL(_vs);

				//Permute
#pragma omp parallel  for num_threads(g_nthread_val)  schedule(static, 1)
				for (int64 m = 0; m < nperm; ++m)
				{
					RNG rng(g_seed_val + pairid * nperm + m);

					threadid = omp_get_thread_num();
					double* fi = Gfi + nlay * threadid;
					double* f = Gf + nlay * nlay * threadid;
					int* g = Gg + nlay * nlay * threadid;
					int* e = Ge + nlay * nlay * threadid;
					double* ef = Gef + nlay * nlay * threadid;
					double* ef2 = Gef2 + nlay * nlay * threadid;
					VESSEL** permbuf = Gpermbuf + Nht * threadid;

					VESSEL& vs2 = _vs2[threadid];
					VESSEL_ITERATOR ve(amova_cind_val == 2 ? 0 : 1, vs2, nlay);

					vs2.Shuffle(rng, fa, fb, 4, permbuf);
					for (int clay = Lind; clay < nlay; ++clay)
					{
						void* Param[] = { (void*)&clay, (void*)&nlay, (void*)&ve };
						CPOINT xx = CPOINT::DownHillSimplex(1, 0, false, 0.1, 10, AMOVA::Likelihood, Param);
						xx.Image2RealSelfing();
						fi[clay] = xx.real[0];
					}

					VESSEL::GetF(fi, f, nlay);

					if (f[pairid] > F[pairid] + 1e-7) g[pairid]++;
					else if (f[pairid] > F[pairid] - 1e-7) e[pairid]++;
					ef[pairid] += f[pairid];
					ef2[pairid] += f[pairid] * f[pairid];

					PROGRESS_VALUE += 5000;
				}

			}
		VLA_DELETE(_vs2);

		if (Lind)
		{
			V[0] = 1 - Fi[nlay - 1];
			V[1] = nlay == 2 ? 1 - V[0] : Fi[1] * (1.0 - Fi[nlay - 1]) / (1.0 - Fi[1]);
			if (nlay > 2) V[nlay - 1] = (Fi[nlay - 1] - Fi[nlay - 2]) / (1.0 - Fi[nlay - 2]);
			for (int clay = 2; clay < nlay - 1; ++clay)
				V[clay] = (Fi[clay] - Fi[clay - 1]) * (1 - Fi[nlay - 1]) / ((1 - Fi[clay - 1]) * (1 - Fi[clay]));
		}
		else
		{
			double VV[100];
			VV[0] = 1 - Fi[nlay - 1];
			VV[1] = nlay == 2 ? 1 - VV[0] : Fi[0] * (1.0 - Fi[nlay - 1]) / (1.0 - Fi[0]);
			if (nlay > 1) VV[nlay] = (Fi[nlay - 1] - Fi[nlay - 2]) / (1.0 - Fi[nlay - 2]);
			for (int clay = 1; clay < nlay - 1; ++clay)
				VV[clay + 1] = (Fi[clay] - Fi[clay - 1]) * (1 - Fi[nlay - 1]) / ((1 - Fi[clay - 1]) * (1 - Fi[clay]));
			for (int i = 0; i < nlay; ++i)
				V[i] = VV[i + 1];
			V[0] += VV[0];
		}

		double SSTOT = 0;
#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int64 l = 0; l < nloc; ++l)
		{
			double nh = total_pop->loc_stat1[l].nhaplo;
			if (nh == 0) continue;

			double* c = Gc + nlay * nlay * threadid;
			double* c2 = Gc2 + nlay * nlay * threadid;
			int* tid = Gtid + (nlay + 1) * threadid;
			double* tw = Gtw + (nlay + 1) * threadid;

			AtomicAddD(SSTOT, (1.0 - SumSquare(total_pop->GetFreq(l), GetLoc(l).k) /* loc_stat1[l].a2 */) * 0.5 * nh);
			_vs.InitC(c, tid, total_pop->loc_stat1[l].nhaplo, nlay);
			_vs.GetC(tw, tid, c, nlay, W, l);
			Add(c2, c, nlay * nlay);
		}

		for (int i = 0; i < g_nthread_val; ++i)
		{
			Add(C, Gc2 + nlay * nlay * i, nlay * nlay);
			Add(G, Gg + nlay * nlay * i, nlay * nlay);
			Add(E, Ge + nlay * nlay * i, nlay * nlay);
			Add(EF, Gef + nlay * nlay * i, nlay * nlay);
			Add(EF2, Gef2 + nlay * nlay * i, nlay * nlay);
		}

		for (int i = 0; i < nlay; ++i)
			DF[i] = round(i > 0 ? (C[i * nlay + 0] - C[(i - 1) * nlay + 0]) : (C[i * nlay + 0]));
		DF[nlay] = C[(nlay - 1) * nlay + 0];

		double sstot = 0;
		for (int i = 0; i < nlay; ++i)
			sstot += C[(nlay - 1) * nlay + i] * V[i];
		double VTOT = SSTOT / sstot;
		Mul(V, VTOT, nlay);

		MatrixMul(C, nlay, nlay, V, nlay, 1, SS);

		VLA_DELETE(Fi);
		VLA_DELETE(C);
		VLA_DELETE(Gfi);
		VLA_DELETE(Gf);
		VLA_DELETE(Gc);
		VLA_DELETE(Gc2);
		VLA_DELETE(Gtid);
		VLA_DELETE(Gtw);
		VLA_DELETE(Gg);
		VLA_DELETE(Ge);
		VLA_DELETE(Gef);
		VLA_DELETE(Gef2);
		delete[] Gpermbuf;
	}

	/* Destructor */
	TARGET AMOVA::~AMOVA()
	{
		for (int clay = 0; clay < nlay; ++clay)
			delete[] SSW[clay];
		delete[] nSSW;
		delete[] SSW;
		delete[] EF;
		delete[] EF2;
		delete[] G;
		delete[] E;
		delete[] DF;
		delete[] SS;
		delete[] V;
		delete[] F;
	}

	/* Write results */
	TARGET void AMOVA::PrintAMOVA(FILE* fout)
	{
		MEMORY mem1;
		amova_memory = &mem1;
		Lind = amova_cind_val == 1;
		double Vtot = Sum(V, nlay);
		VLA_NEW(SSB, double, nlay + 1);
		VLA_NEW(MS, double, nlay + 1);
		for (int i = 0; i < nlay; ++i)
		{
			SSB[i] = i > 0 ? SS[i] - SS[i - 1] : SS[i];
			MS[i] = SSB[i] / DF[i];
		}
		SSB[nlay] = Sum(SSB, nlay);
		MS[nlay] = SSB[nlay] / DF[nlay];

		switch (method)
		{
		case 1:
			fprintf(fout, "%s%sAMOVA Summary, method: homoploid, mutation model: %s, ind-level=%s",
				g_linebreak_val, g_linebreak_val,
				amova_cmutation_val == 1 ? "IAM" : "SMM",
				Lind == 1 ? "yes" : "no");
			break;
		case 2:
			fprintf(fout, "%s%sAMOVA Summary, method: aneuploid, mutation model: %s, ind-level=%s",
				g_linebreak_val, g_linebreak_val,
				amova_cmutation_val == 1 ? "IAM" : "SMM",
				Lind == 1 ? "yes" : "no");
			break;
		case 3:
			fprintf(fout, "%s%sAMOVA Summary, method: likelihood, mutation model: IAM, ind-level=%s",
				g_linebreak_val, g_linebreak_val,
				Lind == 1 ? "yes" : "no");
			break;
		}

		fprintf(fout, "%sSource%cd.f.%cSS%cMS%cVar%cPercentage", g_linebreak_val, g_delimiter_val, g_delimiter_val, g_delimiter_val, g_delimiter_val, g_delimiter_val);

		double inv_npermed = 1.0 / (amova_nperm_val * (amova_test_val == 1) + 1);

		VLA_NEW(level_name, char*, nlay);
		VLA_NEW(level_short, char*, nlay);

		if (Lind)
		{
			level_name[0] = (char*)"Individuals";
			level_name[1] = (char*)"Populations";
			level_name[nlay - 1] = (char*)"Total";
			level_short[0] = (char*)"I";
			level_short[1] = (char*)"S";
			level_short[nlay - 1] = (char*)"T";
			for (int rl = 0; rl < lreg; ++rl)
			{
				level_name[rl + 2] = new char[20];
				level_short[rl + 2] = new char[6];
				sprintf(level_name[rl + 2], "Regions Level %d", rl + 1);
				sprintf(level_short[rl + 2], "C%d", rl + 1);
			}
		}
		else
		{
			level_name[0] = (char*)"Populations";
			level_name[nlay - 1] = (char*)"Total";
			level_short[0] = (char*)"S";
			level_short[nlay - 1] = (char*)"T";
			for (int rl = 0; rl < lreg; ++rl)
			{
				level_name[rl + 1] = new char[20];
				level_short[rl + 1] = new char[6];
				sprintf(level_name[rl + 1], "Regions Level %d", rl + 1);
				sprintf(level_short[rl + 1], "C%d", rl + 1);
			}
		}

		fprintf(fout, "%sWithin %s%c", g_linebreak_val, level_name[0], g_delimiter_val);
		fprintf(fout, "%0.0lf", DF[0]);										fprintf(fout, "%c", g_delimiter_val);
		WriteReal(fout, SSB[0]);											fprintf(fout, "%c", g_delimiter_val);
		WriteReal(fout, MS[0]);												fprintf(fout, "%c", g_delimiter_val);
		WriteReal(fout, V[0]);												fprintf(fout, "%c", g_delimiter_val);
		WriteReal(fout, V[0] * 100.0 / Vtot);								fprintf(fout, "%c", g_delimiter_val);

		for (int i = 1; i < nlay; ++i)
		{
			fprintf(fout, "%sAmong %s%c", g_linebreak_val, level_name[i - 1], g_delimiter_val);
			fprintf(fout, "%0.0lf", DF[i]);									fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, SSB[i]);										fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, MS[i]);											fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, V[i]);											fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, V[i] * 100.0 / Vtot);							fprintf(fout, "%c", g_delimiter_val);
		}
		fprintf(fout, "%sTotal%c", g_linebreak_val, g_delimiter_val);
		fprintf(fout, "%0.0lf", DF[nlay]);									fprintf(fout, "%c", g_delimiter_val);
		WriteReal(fout, SSB[nlay]);											fprintf(fout, "%c", g_delimiter_val);
		WriteReal(fout, MS[nlay]);											fprintf(fout, "%c", g_delimiter_val);
		WriteReal(fout, Vtot);												fprintf(fout, "%c", g_delimiter_val);
		WriteReal(fout, 100.0);												fprintf(fout, "%c", g_delimiter_val);

		fprintf(fout, "%s%sF-statistics", g_linebreak_val, g_linebreak_val);
		fprintf(fout, "%sStatistics%cValue%cPermute Mean%cPermute Var%cPr(rand>obs)%cPr(rand=obs)", g_linebreak_val, g_delimiter_val, g_delimiter_val, g_delimiter_val, g_delimiter_val, g_delimiter_val);

		for (int i = 0; i < nlay; ++i)
			for (int j = i + 1; j < nlay; ++j)
			{
				fprintf(fout, "%sF%s%s%c", g_linebreak_val, level_short[i], level_short[j], g_delimiter_val);
				WriteReal(fout, F[i * nlay + j]);							fprintf(fout, "%c", g_delimiter_val);
				if (method != 3 || (method == 3 && Lind != 0 && i > 0) || (method == 3 && Lind == 0))
				{
					WriteReal(fout, EF[i * nlay + j] * inv_npermed);		fprintf(fout, "%c", g_delimiter_val);
					WriteReal(fout, EF2[i * nlay + j] * inv_npermed - EF[i * nlay + j] * inv_npermed * EF[i * nlay + j] * inv_npermed);
					fprintf(fout, "%c", g_delimiter_val);
					WriteReal(fout, G[i * nlay + j] * inv_npermed);			fprintf(fout, "%c", g_delimiter_val);
					WriteReal(fout, E[i * nlay + j] * inv_npermed);			fprintf(fout, "%c", g_delimiter_val);
				}
				else
					fprintf(fout, "-%c-%c-%c-%c", g_delimiter_val, g_delimiter_val, g_delimiter_val, g_delimiter_val);
			}

		if (method != 3 && amova_printss_val == 1)
		{
			for (int rl = lreg - 1; rl >= 0; --rl)
			{
				int tref = 0;
				VESSEL vs(total_pop, nlay, tref, -1, 1);
				VESSEL_ITERATOR ve1(2 + Lind + rl, vs, nlay);

				fprintf(fout, "%s%sSS within regions level %d %s", g_linebreak_val, g_linebreak_val, rl + 1, g_linebreak_val);
				for (int rl2 = rl; rl2 < lreg; ++rl2)
					fprintf(fout, "RegL%d%c", rl2 + 1, g_delimiter_val);
				fprintf(fout, "SS");
				for (int i = 0; i < nreg[rl]; ++i)
				{
					fprintf(fout, "%s", g_linebreak_val);
					POP* tr = ve1.GetSubpop(nlay, 2 + Lind + rl);
					for (int rl2 = rl; rl2 < lreg; ++rl2)
					{
						fprintf(fout, "%s%c", tr->name, g_delimiter_val);
						tr = aregs[rl2 + 1][tr->rid];
					}
					WriteReal(fout, SSW[Lind + rl + 1][i]);
					ve1.Next(nlay);
				}
			}


			{
				int tref = 0;
				VESSEL vs(total_pop, nlay, tref, -1, 1);
				VESSEL_ITERATOR ve1(1 + Lind, vs, nlay);

				fprintf(fout, "%s%sSS within populations%sPop", g_linebreak_val, g_linebreak_val, g_linebreak_val);
				for (int rl = 0; rl < lreg; ++rl)
					fprintf(fout, "%cRegL%d", g_delimiter_val, rl + 1);
				fprintf(fout, "%cSS", g_delimiter_val);
				for (int i = 0; i < npop; ++i)
				{
					POP* tr = ve1.GetSubpop(nlay, 1 + Lind);
					fprintf(fout, "%s%s%c", g_linebreak_val, tr->name, g_delimiter_val);
					for (int rl = 0; rl < lreg; ++rl)
					{
						tr = aregs[rl][tr->rid];
						fprintf(fout, "%s%c", tr->name, g_delimiter_val);
					}
					WriteReal(fout, SSW[Lind][i]);
					ve1.Next(nlay);
				}
			}

			if (Lind)
			{
				int tref = 0;
				VESSEL vs(total_pop, nlay, tref, -1, 1);
				VESSEL_ITERATOR ve1(1, vs, nlay);

				fprintf(fout, "%s%sSS within individuals%sInd%cPop", g_linebreak_val, g_linebreak_val, g_linebreak_val, g_delimiter_val);
				for (int rl = 0; rl < lreg; ++rl)
					fprintf(fout, "%cRegL%d", g_delimiter_val, rl + 1);
				fprintf(fout, "%cSS", g_delimiter_val);
				for (int i = 0; i < nind; ++i)
				{
					POP* tr = ve1.GetSubpop(nlay, 2);
					IND* ti = ve1.GetInd(nlay, 1);
					fprintf(fout, "%s%s%c%s%c", g_linebreak_val, ti->name, g_delimiter_val, tr->name, g_delimiter_val);
					for (int rl = 0; rl < lreg; ++rl)
					{
						tr = aregs[rl][tr->rid];
						fprintf(fout, "%s%c", tr->name, g_delimiter_val);
					}
					WriteReal(fout, SSW[0][i]);
					ve1.Next(nlay);
				}
			}
		}

		VLA_DELETE(SSB);
		VLA_DELETE(MS);
		for (int rl = 0; rl < lreg; ++rl)
		{
			delete[] level_name[rl + 1 + Lind];
			delete[] level_short[rl + 1 + Lind];
		}
		VLA_DELETE(level_name);
		VLA_DELETE(level_short);
	}
#endif

#define extern 
extern _thread MEMORY* amova_memory;				//memory class for amova vessels
extern AMOVA* amova_buf;							//Read/Write buffer, nthread
extern double* amova_matrix;						//Genetic distance matrix used in AMOVA
#undef extern 

/* Calculate analysis of molecular variance */
TARGET void CalcAMOVA()
{
	if (!amova) return;
	if (ad) Exit("\nError: AMOVA (-amova) is incompatible with allelic depth (-ad) option.\n");

	EvaluationBegin();
	OpenResFile("-amova", "Analysis of molecular variances");

	bool isfirst = true;
	int64 nt = 0;
	for (int i1 = 1; i1 <= 3; ++i1) if (amova_method_val[i1])
		for (int i2 = 1; i2 <= (i1 == 3 ? 1 : 2); ++i2) if (amova_mutation_val[i2])
			for (int i3 = 1; i3 <= 2; ++i3) if (amova_ind_val[i3])
			{
				int64 nperm = amova_nperm_val * (amova_test_val == 1);

				switch (i1)
				{
				case 1: nt += nloc * (int64)(Binomial(total_pop->nhaplotypes, 2) * 0.1 + 0.5) +
					(int64)(Binomial(total_pop->nhaplotypes, 2) + 0.5) +
					nperm * (int64)(BINOMIAL[(i3 == 1) + 2 + lreg][2] + 0.5) * 1000; break;
				case 2:
				{
					if (amova_pseudo_val > 0 && amova_pseudo_val < 10)
						Exit("\nError: In aneuploid AMOVA method, the number of pseudo-permutations should be greater than 10. \n");

					int64 M = amova_pseudo_val * (amova_test_val == 1);
					bool PseudoPerm = amova_pseudo_val > 0;
					if (!PseudoPerm) M = nperm;
					int npair = BINOMIAL[(i3 == 1) + 2 + lreg][2] + 0.5;

					nt += nloc * M * (int64)npair * 500 + nperm * npair * 10;
					break;
				}
				case 3: nt += (int64)(nperm * BINOMIAL[2 + lreg][2] * 5000 + 0.5); break;
				}
			}

	PROGRESS_VALUE2 = PROGRESS_VALUE3 = 0;

	for (int i1 = 1; i1 <= 3; ++i1) if (amova_method_val[i1])
		for (int i2 = 1; i2 <= (i1 == 3 ? 1 : 2); ++i2) if (amova_mutation_val[i2])
			for (int i3 = 1; i3 <= 2; ++i3) if (amova_ind_val[i3])
			{
				int nthread = 1;
				int64 nperm = amova_nperm_val * (amova_test_val == 1);

				amova_cmethod_val = i1;
				amova_cmutation_val = i2;
				amova_cind_val = i3;
				amova_buf = new AMOVA();

				int64 nced = 0;
				switch (i1)
				{
				case 1: nced += nloc * (int64)(Binomial(total_pop->nhaplotypes, 2) * 0.1 + 0.5) +
					(int64)(Binomial(total_pop->nhaplotypes, 2) + 0.5) +
					nperm * (int64)(BINOMIAL[(i3 == 1) + 2 + lreg][2] + 0.5) * 1000; break;
				case 2:
				{
					int64 M = amova_pseudo_val * (amova_test_val == 1);
					bool PseudoPerm = amova_pseudo_val > 0;
					if (!PseudoPerm) M = nperm;
					int npair = BINOMIAL[(i3 == 1) + 2 + lreg][2] + 0.5;
					nced += nloc * M * (int64)npair * 500 + nperm * npair * 10;
					break;
				}
				case 3: nced += (int64)(nperm * BINOMIAL[2 + lreg][2] * 5000 + 0.5); break;
				}

				RunThreads(&AMOVAThread, NULL, NULL, nt, nced,
					"\nPerforming analysis of molecular variance:\n", nthread, isfirst);

				isfirst = false;
				int nlay = amova_buf[0].nlay;
				for (int i = 1; i < nthread; ++i)
				{
					Add(amova_buf[0].G, amova_buf[i].G, nlay * nlay);
					Add(amova_buf[0].E, amova_buf[i].E, nlay * nlay);
					Add(amova_buf[0].EF2, amova_buf[i].EF2, nlay * nlay);
					if (amova_buf[i].nSSW) for (int j = 0; j < nlay; ++j)
						Add(amova_buf[0].SSW[j], amova_buf[i].SSW[j], amova_buf[i].nSSW[j]);
				}
				amova_buf[0].PrintAMOVA(FRES);

				delete amova_buf;
			}

	CloseResFile();
	EvaluationEnd("AMOVA");
}

/* Calculate AMOVA using multiple threads */
THREAD(AMOVAThread)
{
	switch (amova_cmethod_val)
	{
	case 1: amova_buf[threadid].CalcAMOVA_homo(); break;
	case 2: amova_buf[threadid].CalcAMOVA_aneu(); break;
	case 3: amova_buf[threadid].CalcAMOVA_ml(); break;
	}
}
