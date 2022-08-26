/* Bayesian clustering functions */

#include "vcfpop.h"

#ifndef _BAYESIAN_READER

    /* Do nothing */
    TARGET BAYESIAN_READER::BAYESIAN_READER()
    {

    }

    /* Initialize reader */
    TARGET BAYESIAN_READER::BAYESIAN_READER(int indid)
    {
        SetZero(this, 1);
        pos = (uint64*)structure_allele[indid];
    }

#endif

#ifndef _BAYESIAN_WRITER

    TARGET BAYESIAN_WRITER::BAYESIAN_WRITER()
    {

    }

    /* Initialize writer */
    TARGET BAYESIAN_WRITER::BAYESIAN_WRITER(int indid)
    {
        SetZero(this, 1);
        pos = (uint64*)structure_allele[indid];
    }

    /* Write id of next ind to buffer */
    TARGET void BAYESIAN_WRITER::Write(int size, int aid)
    {
        aid &= (1u << size) - 1u;

        //if data is full
        if (nbits + size > 64) [[unlikely]]
        {
            if (nbits == 64)
            {
                *pos++ = data64[0];
                data64[0] = (uint64)aid;
                nbits = size;
            }
            else
            {
                // number of bits can write to data
                int wbits = 64 - nbits;

                // write wbits to data
                data64[0] |= (uint64)aid << nbits;

                // write 64 bit data to pos
                *pos++ = data64[0];

                // write remain bits to data
                data64[0] = aid >> wbits;
                nbits = size - wbits;
            }
        }
        else [[likely]]
        {
            data64[0] |= (uint64)aid << nbits;
            nbits += size;
        }
    }

    /* Write id of next ind to buffer */
#ifdef __aarch64__
    TARGET void BAYESIAN_WRITER::Write8(int size, uint64 aid[8])
    {
        uint64x2_t* aid128 = (uint64x2_t*)aid;

        uint64x2_t mask = vdupq_n_u64((1ull << size) - 1ull);
        int64x2_t nlbits;

        aid128[0] = vandq_u64(aid128[0], mask);
        aid128[1] = vandq_u64(aid128[1], mask);
        aid128[2] = vandq_u64(aid128[2], mask);
        aid128[3] = vandq_u64(aid128[3], mask);

        //if data is full
        if (nbits + size > 64) [[unlikely]]
        {
            if (nbits == 64)
            {
                uint64* tpos = pos;

                *tpos = data64[0]; tpos += structure_indnbytes;
                *tpos = data64[1]; tpos += structure_indnbytes;
                *tpos = data64[2]; tpos += structure_indnbytes;
                *tpos = data64[3]; tpos += structure_indnbytes;
                *tpos = data64[4]; tpos += structure_indnbytes;
                *tpos = data64[5]; tpos += structure_indnbytes;
                *tpos = data64[6]; tpos += structure_indnbytes;
                *tpos = data64[7]; tpos += structure_indnbytes;

                pos++;

                data128[0] = aid128[0];
                data128[1] = aid128[1];
                data128[2] = aid128[2];
                data128[3] = aid128[3];

                nbits = size;
            }
            else
            {
                // number of bits can write to data
                const int wbits = 64 - nbits;

                // write wbits to data
                data64[0] |= (uint64)aid[0] << nbits;
                data64[1] |= (uint64)aid[1] << nbits;
                data64[2] |= (uint64)aid[2] << nbits;
                data64[3] |= (uint64)aid[3] << nbits;
                data64[4] |= (uint64)aid[4] << nbits;
                data64[5] |= (uint64)aid[5] << nbits;
                data64[6] |= (uint64)aid[6] << nbits;
                data64[7] |= (uint64)aid[7] << nbits;

                // write 64 bit data to pos
                uint64* tpos = pos;

                *tpos = data64[0]; tpos += structure_indnbytes;
                *tpos = data64[1]; tpos += structure_indnbytes;
                *tpos = data64[2]; tpos += structure_indnbytes;
                *tpos = data64[3]; tpos += structure_indnbytes;
                *tpos = data64[4]; tpos += structure_indnbytes;
                *tpos = data64[5]; tpos += structure_indnbytes;
                *tpos = data64[6]; tpos += structure_indnbytes;
                *tpos = data64[7]; tpos += structure_indnbytes;

                pos++;

                // write remain bits to data
                nlbits = vdupq_n_s64(-wbits);  /* right move wbits (negative) */
                data128[0] = vshlq_u64(aid128[0], nlbits);
                data128[1] = vshlq_u64(aid128[1], nlbits);
                data128[2] = vshlq_u64(aid128[2], nlbits);
                data128[3] = vshlq_u64(aid128[3], nlbits);

                nbits = size - wbits;
            }
        }
        else [[likely]]
        {
            nlbits = vdupq_n_s64(nbits); /* left move nbits (positive) */
            data128[0] = vorrq_u64(data128[0], vshlq_u64(aid128[0], nlbits));
            data128[1] = vorrq_u64(data128[1], vshlq_u64(aid128[1], nlbits));
            data128[2] = vorrq_u64(data128[2], vshlq_u64(aid128[2], nlbits));
            data128[3] = vorrq_u64(data128[3], vshlq_u64(aid128[3], nlbits));

            nbits += size;
        }
    }
#else
    TARGET void BAYESIAN_WRITER::Write8(int size, uint64 aid[STRUCTURE_NPACK])
    {
        __m128i* aid128 = (__m128i*)aid;

        __m128i mask = _mm_set1_epi64x((1ull << size) - 1ull);
        aid128[0] = _mm_and_si128(aid128[0], mask);
        aid128[1] = _mm_and_si128(aid128[1], mask);
        aid128[2] = _mm_and_si128(aid128[2], mask);
        aid128[3] = _mm_and_si128(aid128[3], mask);

        //if data is full
        if (nbits + size > 64) [[unlikely]]
        {
            if (nbits == 64)
            {
                uint64* tpos = pos;

                *tpos = data64[0]; tpos += structure_indnbytes;
                *tpos = data64[1]; tpos += structure_indnbytes;
                *tpos = data64[2]; tpos += structure_indnbytes;
                *tpos = data64[3]; tpos += structure_indnbytes;
                *tpos = data64[4]; tpos += structure_indnbytes;
                *tpos = data64[5]; tpos += structure_indnbytes;
                *tpos = data64[6]; tpos += structure_indnbytes;
                *tpos = data64[7]; tpos += structure_indnbytes;

                pos++;

                data128[0] = aid128[0];
                data128[1] = aid128[1];
                data128[2] = aid128[2];
                data128[3] = aid128[3];

                nbits = size;
            }
            else
            {
                // number of bits can write to data
                int wbits = 64 - nbits;

                // write wbits to data
                data64[0] |= (uint64)aid[0] << nbits;
                data64[1] |= (uint64)aid[1] << nbits;
                data64[2] |= (uint64)aid[2] << nbits;
                data64[3] |= (uint64)aid[3] << nbits;
                data64[4] |= (uint64)aid[4] << nbits;
                data64[5] |= (uint64)aid[5] << nbits;
                data64[6] |= (uint64)aid[6] << nbits;
                data64[7] |= (uint64)aid[7] << nbits;

                // write 64 bit data to pos
                uint64* tpos = pos;

                *tpos = data64[0]; tpos += structure_indnbytes;
                *tpos = data64[1]; tpos += structure_indnbytes;
                *tpos = data64[2]; tpos += structure_indnbytes;
                *tpos = data64[3]; tpos += structure_indnbytes;
                *tpos = data64[4]; tpos += structure_indnbytes;
                *tpos = data64[5]; tpos += structure_indnbytes;
                *tpos = data64[6]; tpos += structure_indnbytes;
                *tpos = data64[7]; tpos += structure_indnbytes;

                pos++;

                // write remain bits to data
                data128[0] = _mm_srli_epi64(aid128[0], wbits);
                data128[1] = _mm_srli_epi64(aid128[1], wbits);
                data128[2] = _mm_srli_epi64(aid128[2], wbits);
                data128[3] = _mm_srli_epi64(aid128[3], wbits);

                nbits = size - wbits;
            }
        }
        else [[likely]]
        {
            data128[0] = _mm_or_si128(data128[0], _mm_slli_epi64(aid128[0], nbits));
            data128[1] = _mm_or_si128(data128[1], _mm_slli_epi64(aid128[1], nbits));
            data128[2] = _mm_or_si128(data128[2], _mm_slli_epi64(aid128[2], nbits));
            data128[3] = _mm_or_si128(data128[3], _mm_slli_epi64(aid128[3], nbits));

            nbits += size;
        }
    }
#endif
    /* Write all remaining bits to buffer */
    TARGET void BAYESIAN_WRITER::FinishWrite()
    {
        byte* pos2 = (byte*)(pos);
        byte* data2 = (byte*)&data64[0];
        for (int i = 0; i * 8 < nbits; ++i)
            pos2[i] = data2[i];
    }

    /* Write all remaining bits to buffer */
    TARGET void BAYESIAN_WRITER::FinishWrite8()
    {
        byte* pos2, *data2;
        for (int k = 0; k < 8; ++k)
        {
            pos2 = (byte*)(pos + structure_indnbytes * k);
            data2 = (byte*)&data64[k];
            for (int i = 0; i * 8 < nbits; ++i)
                pos2[i] = data2[i];
        }
    }

#endif

#ifndef _BAYESIAN

    /* Set all bits to 0 */
    TARGET BAYESIAN::BAYESIAN()
    {
        SetZero(this, 1);
    }

    /* Write results for a run */
    TARGET void BAYESIAN::PrintStructure()
    {
        char filename[PATH_LEN];
        char name_buf[NAME_BUF_LEN];
        sprintf(filename, "%s.structure.k=%d_rep=%d_id=%d.txt", g_output_val.c_str(), K, par2->rep, par2->id);
        FILE* fout = fopen(filename, "wb");

        fprintf(fout, "%s%sParameters:%s", g_linebreak_val, g_linebreak_val, g_linebreak_val);
        fprintf(fout, "Seed=%lld%s", rng.seed, g_linebreak_val);
        fprintf(fout, "Model=%s,%s,%s%s",
            admix ? "ADM" : "NOADM",
            locpriori ? "LOC" : "NOLOC",
            fmodel ? "F" : "NOF", g_linebreak_val);

        /*
        fprintf(fout, "#Inds=%d%s", N, g_linebreak_val);
        fprintf(fout, "#Loci=%d%s", L, g_linebreak_val);
        fprintf(fout, "#Location=%d%s", S, g_linebreak_val);
        fprintf(fout, "#K=%d%s", K, g_linebreak_val);
        fprintf(fout, "#Burnin=%d%s", nburnin, g_linebreak_val);
        fprintf(fout, "#Reps=%d%s", nreps, g_linebreak_val);
        fprintf(fout, "#Thinning=%d%s", nthinning, g_linebreak_val);
        fprintf(fout, "#Runs=%d%s", nruns, g_linebreak_val);
        fprintf(fout, "#Admburnin=%d%s", nadmburnin, g_linebreak_val);

        fprintf(fout, "lambda="); WriteReal(fout, lambda); fprintf(fout, "%s", g_linebreak_val);
        fprintf(fout, "stdlambda="); WriteReal(fout, stdlambda); fprintf(fout, "%s", g_linebreak_val);
        fprintf(fout, "maxlambda="); WriteReal(fout, maxlambda); fprintf(fout, "%s", g_linebreak_val);
        fprintf(fout, "inferlambda=%s%s", inferlambda ? "yes" : "no", g_linebreak_val);
        fprintf(fout, "difflambda=%s%s", difflambda ? "yes" : "no", g_linebreak_val);
        fprintf(fout, "diversity=%s%s", diversity == 1 ? "yes" : "no", g_linebreak_val);

        fprintf(fout, "alpha="); WriteReal(fout, alpha); fprintf(fout, "%s", g_linebreak_val);
        fprintf(fout, "inferalpha=%s%s", inferalpha ? "yes" : "no", g_linebreak_val);
        fprintf(fout, "diffalpha=%s%s", diffalpha ? "yes" : "no", g_linebreak_val);
        fprintf(fout, "uniformalpha=%s%s", uniformalpha ? "yes" : "no", g_linebreak_val);
        fprintf(fout, "stdalpha="); WriteReal(fout, stdalpha); fprintf(fout, "%s", g_linebreak_val);
        fprintf(fout, "maxalpha="); WriteReal(fout, maxalpha); fprintf(fout, "%s", g_linebreak_val);
        fprintf(fout, "alphapriora="); WriteReal(fout, alphapriora); fprintf(fout, "%s", g_linebreak_val);
        fprintf(fout, "alphapriorb="); WriteReal(fout, alphapriorb); fprintf(fout, "%s", g_linebreak_val);
        fprintf(fout, "metrofreq=%d%s", metrofreq, g_linebreak_val);

        fprintf(fout, "r="); WriteReal(fout, r); fprintf(fout, "%s", g_linebreak_val);
        fprintf(fout, "maxr="); WriteReal(fout, maxr); fprintf(fout, "%s", g_linebreak_val);
        fprintf(fout, "epsr="); WriteReal(fout, epsr); fprintf(fout, "%s", g_linebreak_val);
        fprintf(fout, "epseta="); WriteReal(fout, epseta); fprintf(fout, "%s", g_linebreak_val);
        fprintf(fout, "epsgamma="); WriteReal(fout, epsgamma); fprintf(fout, "%s", g_linebreak_val);

        fprintf(fout, "pmeanf="); WriteReal(fout, pmeanf); fprintf(fout, "%s", g_linebreak_val);
        fprintf(fout, "pstdf="); WriteReal(fout, pstdf); fprintf(fout, "%s", g_linebreak_val);
        fprintf(fout, "stdf="); WriteReal(fout, stdf); fprintf(fout, "%s", g_linebreak_val);
        fprintf(fout, "fsame=%s%s", inferalpha ? "yes" : "no", g_linebreak_val);
        */


        par2->MeanlnL = rout[0];
        par2->VarlnL = rout[2];
        par2->lnPD = rout[3];

        fprintf(fout, "Number of cluster=%d%s", K, g_linebreak_val);
        fprintf(fout, "Replicate=%d%s", par2->rep, g_linebreak_val);
        fprintf(fout, "Run id=%d%s", par2->id, g_linebreak_val);
        fprintf(fout, "Mean value of ln likelihood="); WriteReal(fout, par2->MeanlnL); fprintf(fout, "%s", g_linebreak_val);
        fprintf(fout, "Variance of ln likelihood="); WriteReal(fout, par2->VarlnL); fprintf(fout, "%s", g_linebreak_val);
        fprintf(fout, "Estimated Ln Prob of Data="); WriteReal(fout, par2->lnPD); fprintf(fout, "%s%s", g_linebreak_val, g_linebreak_val);


        int nf = fmodel ? (fsame ? 1 : K) : 0;
        int nl = difflambda ? K : 1;

        if (inferlambda)
        {
            if (!difflambda)
            {
                fprintf(fout, "%sMean value of lambda=", g_linebreak_val); WriteReal(fout, rout[4]);
            }
            else for (int j = 0; j < K; ++j)
            {
                fprintf(fout, "%sMean value of lambda_%d=", g_linebreak_val, j + 1); WriteReal(fout, rout[4 + j]);
            }
        }

        if (locpriori)
        {
            fprintf(fout, "%sMean value of r=", g_linebreak_val); WriteReal(fout, rout[4 + nl]);
            fprintf(fout, "%sMean value of global %s", g_linebreak_val, admix ? "alpha" : "eta");
            fprintf(fout, "%s%cCluster%s", g_linebreak_val, g_delimiter_val, g_linebreak_val);
            for (int j = 0; j < K; ++j)
                fprintf(fout, "%c%d", g_delimiter_val, j + 1);
            fprintf(fout, "%s", g_linebreak_val);
            for (int j = 0; j < K; ++j)
            {
                fprintf(fout, "%c", g_delimiter_val); WriteReal(fout, rout[4 + nl + 1 + j]);
            }
            fprintf(fout, "%sMean value of local %s for each location", g_linebreak_val, admix ? "alpha" : "gamma");
            fprintf(fout, "%s%cCluster%sPop", g_linebreak_val, g_delimiter_val, g_linebreak_val);
            for (int j = 0; j < K; ++j)
                fprintf(fout, "%c%d", g_delimiter_val, j + 1);
            for (int i = 0; i < S; ++i)
            {
                fprintf(fout, "%s%s", g_linebreak_val, apops[i]->name);
                for (int j = 0; j < K; ++j)
                {
                    fprintf(fout, "%c", g_delimiter_val); WriteReal(fout, rout[4 + nl + 1 + K + i * K + j]);
                }
            }
        }
        else if (admix && inferalpha)
        {
            if (!diffalpha)
            {
                fprintf(fout, "%sMean value of alpha=", g_linebreak_val); WriteReal(fout, rout[4 + nl]);
            }
            else for (int j = 0; j < K; ++j)
            {
                fprintf(fout, "%sMean value of alpha_%d=", g_linebreak_val, j + 1); WriteReal(fout, rout[4 + nl + j]);
            }
        }

        if (fmodel)
        {
            if (fsame)
            {
                fprintf(fout, "%sMean value of Fst=", g_linebreak_val); WriteReal(fout, rout[rlen - 1]);
            }
            else
            {
                fprintf(fout, "%sMean value of Fst", g_linebreak_val);
                fprintf(fout, "%s%cCluster%s", g_linebreak_val, g_delimiter_val, g_linebreak_val);
                for (int j = 0; j < K; ++j)
                    fprintf(fout, "%c%d", g_delimiter_val, j + 1);
                fprintf(fout, "%s", g_linebreak_val);
                for (int i = 0; i < nf; ++i)
                {
                    fprintf(fout, "%c", g_delimiter_val); WriteReal(fout, rout[rlen - nf + i]);
                }
            }
        }

        VLA_NEW(O, double, S * K);
        SetZero(O, S * K);
        fprintf(fout, "%s%sProportion of membership of each pre-defined population in each of the %d clusters", g_linebreak_val, g_linebreak_val, K);
        fprintf(fout, "%s%cCluster%sPop", g_linebreak_val, g_delimiter_val, g_linebreak_val);
        for (int i = 0; i < N; ++i)
            Add(O + ainds[i]->popid * K, Q + i * K, K);
        Unify(O, S, K);

        for (int k = 0; k < K; ++k)
            fprintf(fout, "%c%d", g_delimiter_val, k + 1);

        for (int s = 0; s < S; ++s)
        {
            fprintf(fout, "%s%s", g_linebreak_val, apops[s]->name);
            for (int k = 0; k < K; ++k)
            {
                fprintf(fout, "%c", g_delimiter_val);
                WriteReal(fout, O[s * K + k]);
            }
        }

        fprintf(fout, "%s%sInferred ancestry of individuals", g_linebreak_val, g_linebreak_val);
        fprintf(fout, "%s%c%cCluster", g_linebreak_val, g_delimiter_val, g_delimiter_val);
        fprintf(fout, "%sInd%cPop%cRegL1", g_linebreak_val, g_delimiter_val, g_delimiter_val);
        for (int j = 0; j < K; ++j)
            fprintf(fout, "%c%d", g_delimiter_val, j + 1);

        double* MiSumDouble = (double*)MiSum;
        for (int i = 0; i < N; ++i)
        {
            fprintf(fout, "%s%s%c%s%c%s", g_linebreak_val, ainds[i]->name,
                g_delimiter_val, apops[ainds[i]->popid]->name,
                g_delimiter_val, apops[ainds[i]->popid]->rid == -1 ? "Total" : aregs[0][apops[ainds[i]->popid]->rid]->name
            );
            for (int k = 0; k < K; ++k)
            {
                fprintf(fout, "%c", g_delimiter_val);
                WriteReal(fout, MiSumDouble[i * K + k]);
            }
        }
        VLA_DELETE(O);

        VLA_NEW(H, double, K);
        VLA_NEW(D, double, K * K);
        double* H2 = new double[K * L];

        SetZero(H, K);
        for (int k = 0; k < K; ++k)
        {
            double* p = clusterb[k].bucket;
            for (int64 l = 0; l < L; ++l)
            {
                int k2 = GetLoc(l).k;
                double h = 1 - SumSquare(p, k2);
                H2[k * L + l] = h;
                H[k] += h;
                p += k2;
            }
            H[k] /= L;

            for (int j = 0; j <= k; ++j)
            {
                double* p1 = clusterb[k].bucket;
                double* p2 = clusterb[j].bucket;
                double dt = 0;
                for (int64 l = 0; l < L; ++l)
                {
                    int k2 = GetLoc(l).k;
                    dt += 1 - SumProd(p1, p2, k2);
                    p1 += k2;
                    p2 += k2;
                }
                D[j * K + k] = D[k * K + j] = dt / L - (H[k] + H[j]) / 2;
            }
        }

        fprintf(fout, "%s%sAllele-frequency divergence among pops (Net nucleotide distance)", g_linebreak_val, g_linebreak_val);
        fprintf(fout, "%s%cCluster%sCluster", g_linebreak_val, g_delimiter_val, g_linebreak_val);
        for (int j = 0; j < K; ++j)
            fprintf(fout, "%c%d", g_delimiter_val, j + 1);

        for (int i = 0; i < K; ++i)
        {
            fprintf(fout, "%s%d", g_linebreak_val, i + 1);
            for (int j = 0; j < K; ++j)
            {
                fprintf(fout, "%c", g_delimiter_val);
                WriteReal(fout, D[i * K + j]);
            }
        }

        if (structure_diversity_val == 1)
        {
            fprintf(fout, "%s%sEstimated genetic diversity in each cluster", g_linebreak_val, g_linebreak_val);
            fprintf(fout, "%sLocus%cCluster%cHeterozygosity", g_linebreak_val, g_delimiter_val, g_delimiter_val);
            for (int a = 0; a < maxK; ++a)
                fprintf(fout, "%cAllele%d", g_delimiter_val, a);

            for (int k = 0; k < K; ++k)
            {
                fprintf(fout, "%s%s%c%d", g_linebreak_val, k ? "" : "Average", g_delimiter_val, k + 1);
                WriteReal(fout, H[k]);
                bufK[k] = clusterb[k].bucket;
            }

            for (int64 l = 0; l < L; ++l)
            {
                for (int k = 0; k < K; ++k)
                {
                    fprintf(fout, "%s%s%c%d%c", g_linebreak_val, k ? "" : GetLoc(l).GetNameStr(name_buf), g_delimiter_val, k + 1, g_delimiter_val);
                    WriteReal(fout, H2[k * L + l]);
                    for (int a = 0; a < GetLoc(l).k; ++a)
                    {
                        fprintf(fout, "%c", g_delimiter_val);
                        WriteReal(fout, *bufK[k]++);
                    }
                }
            }
        }

        VLA_DELETE(H);
        VLA_DELETE(D);
        delete[] H2;

        fclose(fout);
        Uninit();
    }

    /* Write results summary for all runs */
    TARGET void BAYESIAN::PrintSummary(STRUCTURE_RUNINFO* sp, int len)
    {
        OpenResFile("-structure", "#Bayesian clustering");

        fprintf(FRES, "%s%sID%ck%crep%cMean lnL%cVar lnL%clnP(D)",
            g_linebreak_val, g_linebreak_val, g_delimiter_val, g_delimiter_val, g_delimiter_val, g_delimiter_val, g_delimiter_val);
        for (int i = 0; i < len; ++i)
        {
            fprintf(FRES, "%s%d%c%d%c%d%c",
                g_linebreak_val, sp[i].id,
                g_delimiter_val, sp[i].k,
                g_delimiter_val, sp[i].rep, g_delimiter_val);
            WriteReal(FRES, sp[i].MeanlnL);
            fprintf(FRES, "%c", g_delimiter_val);
            WriteReal(FRES, sp[i].VarlnL);
            fprintf(FRES, "%c", g_delimiter_val);
            WriteReal(FRES, sp[i].lnPD);
        }

        CloseResFile();
    }

    /* Copy parameters */
    TARGET void BAYESIAN::ReadPar(STRUCTURE_RUNINFO* _par2)
    {
        par2 = _par2;
        new(&rng) RNG(g_seed_val + par2->id);

        S = npop;

        K = par2->k;
        N = nind;
        L = nloc;

        admix = structure_admix_val == 1;
        locpriori = structure_locpriori_val == 1;
        fmodel = structure_f_val == 1;

        nburnin = structure_nburnin_val;
        nreps = structure_nreps_val;
        nthinning = structure_nthinning_val;
        nruns = structure_nruns_val;
        nadmburnin = structure_nadmburnin_val;

        lambda = structure_lambda_val;
        stdlambda = structure_stdlambda_val;
        maxlambda = structure_maxlambda_val;
        inferlambda = structure_inferlambda_val == 0;
        difflambda = structure_difflambda_val == 1;
        diversity = structure_diversity_val;

        alpha = structure_alpha_val;
        inferalpha = structure_inferalpha_val == 1;
        diffalpha = structure_diffalpha_val == 1;
        uniformalpha = structure_uniformalpha_val == 1;
        stdalpha = structure_stdalpha_val;
        maxalpha = structure_maxalpha_val;
        alphapriora = structure_alphapriora_val;
        alphapriorb = structure_alphapriorb_val;
        metrofreq = structure_metrofreq_val;

        r = structure_r_val;
        maxr = structure_maxr_val;
        epsr = structure_stdr_val;
        epseta = structure_epseta_val;
        epsgamma = structure_epsgamma_val;

        pmeanf = structure_pmeanf_val;
        pstdf = structure_pstdf_val;
        stdf = structure_stdf_val;
        fsame = structure_singlef_val == 1;

        kdis = new int[maxK + 1];
        SetZero(kdis, maxK + 1);

        for (int64 l = 0; l < L; ++l)
            kdis[GetLoc(l).k]++;
    }

    /* Initialize MCMC */
    TARGET void BAYESIAN::InitFix()
    {
        nr = 0;
        r = 1;
        bufK = new double* [K];

        //bufNK1 = (double*) operator new[](sizeof(double) * N * K, (std::align_val_t)(64));
        //bufNK2 = (double*) operator new[](sizeof(double) * N * K, (std::align_val_t)(64));
        bufNK1o = (byte*)malloc(sizeof(double) * N * K + 64); //alignment for SIMD 
        bufNK2o = (byte*)malloc(sizeof(double) * N * K + 64); //alignment for SIMD 
        bufNK1 = (double*)((uint64)bufNK1o + (64 - ((uint64)bufNK1o % 64)) % 64);
        bufNK2 = (double*)((uint64)bufNK2o + (64 - ((uint64)bufNK2o % 64)) % 64);

        bufN1 = new double[N];
        bufN2 = new double[N];

        MiSum = new int64[N * K];
        SetZero(MiSum, N * K);

        Base = new double[(K * 2 + 3) * KT];
        SetZero(Base, (K * 2 + 3) * KT);

        cluster = new SCLUSTER[K];
        clusterb = new SCLUSTER[K];
        for (int k = 0; k < K; ++k)
        {
            cluster[k].SetFreq(Base + k * KT);
            clusterb[k].SetFreq(Base + k * KT + K * KT);
            double* p = cluster[k].bucket;
            for (int64 l = 0; l < L; ++l)
            {
                int k2 = GetLoc(l).k;
                SetVal(p, 1.0 / k2, k2);
                p += k2;
            }
        }

        Ni = new int[K * KT];
        SetZero(Ni, K * KT);

        Mi = new int64[N * K];
        SetZero(Mi, N * K);

        Q = new double[N * K];
        SetVal(Q, 1.0 / K, N * K);

        Alpha = new double[K];
        SetVal(Alpha, alpha, K);

        Lambda = new double[K];
        SetVal(Lambda, lambda, K);

        int nf = fmodel ? (fsame ? 1 : K) : 0;
        int na = admix ? (diffalpha ? K : 1) : 0;
        int nl = difflambda ? K : 1;

        if (locpriori) rlen = 4 + nl + 1 + K + K * S + nf;
        else rlen = 4 + nl + na + nf;

        rout = new double[rlen];
        SetZero(rout, rlen);
    }

    /* Initialize MCMC for admix model */
    TARGET void BAYESIAN::InitAdmix()
    {
        // Init Z, Z, Mi and Ni

        //Adm
        Z = new ushort[N];
        SetZero(Z, N);

        if (nadmburnin == 0 && !locpriori && !admix)
        {
            singlez = true;
            for (int i = 0; i < N; ++i)
            {
                Z[i] = (ushort)rng.Next(K);
                Mi[i * K + Z[i]] = ainds[i]->vt;
            }

            for (int64 l = 0, o = 0; l < L; ++l, o += GetLoc(l).k)
            {
                GENO_READER rt(0, l);
                GENOTYPE* gtab = GetLoc(l).GetGtab();

                for (int i = 0; i < N; ++i)
                {
                    GENOTYPE& gt = gtab[rt.Read()];
                    if (gt.Nalleles() == 0) continue;

                    int* ni = Ni + Z[i] * KT + o;
                    ushort* als = gt.GetAlleleArray();
                    for (int a = 0, vi = gt.Ploidy(); a < vi; ++a)
                        ni[als[a]]++;
                }
            }
        }
        else
        {
            for (int64 l = 0, o = 0; l < L; ++l, o += GetLoc(l).k)
            {
                GENO_READER rt(0, l);
                GENOTYPE* gtab = GetLoc(l).GetGtab();
                int64* mi = Mi;

                for (int i = 0; i < N; ++i, mi += K)
                {
                    GENOTYPE& gt = gtab[rt.Read()];
                    if (gt.Nalleles() == 0) continue;

                    int* ni = Ni + o;
                    ushort* als = gt.GetAlleleArray();
                    for (int a = 0, vi = gt.Ploidy(); a < vi; ++a)
                    {
                        ushort k2 = (ushort)rng.Next(K);
                        mi[k2]++;
                        ni[k2 * KT + als[a]]++;
                    }
                }
            }
        }
    }

    /* Initialize MCMC for locpriori model */
    TARGET void BAYESIAN::InitLocPriori()
    {
        if (locpriori)
        {
            SumAlpha = new double[S];
            SetZero(SumAlpha, S);

            AlphaLocal = new double[S * K];
            SetVal(AlphaLocal, alpha, S * K);

            if (!admix)
            {
                Eta = new double[K];
                SetVal(Eta, 1.0 / K, K);

                Gamma = new double[S * K];
                SetVal(Gamma, 1.0 / K, S * K);

                Di = new double[S * K];
                SetZero(Di, S * K);
            }
        }
    }

    /* Initialize MCMC for F model */
    TARGET void BAYESIAN::InitFmodel()
    {
        if (fmodel)
        {
            f = new double[K];
            F = new double[K];

            PA.SetFreq(Base + 2 * K * KT);
            PA1.SetFreq(Base + 2 * K * KT + KT);
            PA2.SetFreq(Base + 2 * K * KT + 2 * KT);

            for (int k = 0; k < K; ++k)
            {
                F[k] = pmeanf;
                f[k] = (1 - F[k]) / F[k];
            }

            for (int k = 0; k < K; ++k)
                SetVal(cluster[k].bucket, Lambda[k], KT);

            double* p = PA.bucket;
            SetVal(p, lambda, KT);
            for (int64 l = 0; l < L; ++l)
            {
                int k2 = GetLoc(l).k;
                double* p2 = total_pop->GetFreq(l);
                double nhaplo = total_pop->loc_stat1[l].nhaplo;
                for (int a = 0; a < k2; ++a)
                    p[a] = (lambda + p2[a] * nhaplo) / (k2 * lambda + nhaplo);
                p += k2;
            }
        }
    }

    /* Update allele frequency for all clusters */
    TARGET void BAYESIAN::UpdateP()
    {
        //set lambda, dirichlet distribution priori parameter
        if (fmodel)
        {
            if (fsame)
            {
                Mul(cluster[0].bucket, PA.bucket, f[0], KT);
                for (int k = 1; k < K; ++k)
                    SetVal(cluster[k].bucket, cluster[0].bucket, KT);
            }
            else for (int k = 0; k < K; ++k)
                Mul(cluster[k].bucket, PA.bucket, f[k], KT);
        }
        else for (int k = 0; k < K; ++k)
            SetVal(cluster[k].bucket, Lambda[k], KT);

        int* ni = Ni;
        for (int k = 0; k < K; ++k)
        {
            double* p = cluster[k].bucket;
            for (int64 l = 0; l < L; ++l)
            {
                int k2 = GetLoc(l).k;
                rng.Dirichlet(p, p, ni, k2);
                p += k2;
                ni += k2;
            }
        }
    }

    /* Update a priori ancetral proportion for non-admix model */
    TARGET void BAYESIAN::UpdateQNoAdmix()//new 3x faster
    {
        if (STRUCTURE_DECOMP)
        {
            SetZero(Q, N * K);

            switch (SIMD_TYPE)
            {
#ifdef __aarch64__
            case 2: UpdateQNoAdmixNEO(Q, K, bufNK1, bufNK2, Gamma, rngNEO, locpriori, Base, Z); break;
#else
            case 4: UpdateQNoAdmix512(Q, K, bufNK1, bufNK2, Gamma, rng512, locpriori, Base, Z); break;
            case 3: UpdateQNoAdmixAVX(Q, K, bufNK1, bufNK2, Gamma, rngAVX, locpriori, Base, Z); break;
            case 2: UpdateQNoAdmixSSE(Q, K, bufNK1, bufNK2, Gamma, rngSSE, locpriori, Base, Z); break;
#endif
            }

            binaryq = true;
            return;
        }

        SetZero(Q, N * K);
        OpenLog((int64*)bufNK1, bufNK2, N * K);

        //add priori probability
        double* bufi = bufNK1, * bufbi = bufNK2;
        if (locpriori) for (int i = 0; i < N; ++i, bufi += K, bufbi += K)
        {
            if (ainds[i]->vt == 0) continue;
            ChargeLog((int64*)bufi, bufbi, Gamma + ainds[i]->popid * K, K);
        }

        double* p = Base;
        for (int64 l = 0; l < L; ++l, p += GetLoc(l).k)
        {
            GENO_READER rt(0, l);
            GENOTYPE* gtab = GetLoc(l).GetGtab();

            double* bufl = bufNK1, * bufbl = bufNK2;
            for (int i = 0; i < N; ++i, bufl += K, bufbl += K)
            {
                GENOTYPE& gt = gtab[rt.Read()];
                if (gt.Nalleles() == 0) continue;
                ushort* als = gt.GetAlleleArray();
                for (int a = 0, vi = gt.Ploidy(); a < vi; ++a)
                    ChargeLog((int64*)bufl, bufbl, p + als[a], K, KT);
            }
        }
        CloseLog((int64*)bufNK1, bufNK2, N * K);

        bufi = bufNK1;
        double* q = Q;
        for (int i = 0; i < N; ++i, bufi += K, q += K)
        {
            if (ainds[i]->vt == 0) continue;
            ushort k2 = (ushort)rng.PolyLog(bufi, K);
            q[k2] = 1;
            Z[i] = k2;
        }

        binaryq = true;
    }

    /* Update a priori ancetral proportion for admix model */
    TARGET void BAYESIAN::UpdateQAdmix()
    {
        //used by noadmix (in admburnin) and admix model
        if (locpriori)
        {
            for (int i = 0; i < N; ++i)
                if (ainds[i]->vt)
                    rng.Dirichlet(Q + i * K, AlphaLocal + ainds[i]->popid * K, Mi + i * K, K);
        }
        else
        {
            for (int i = 0; i < N; ++i)
                if (ainds[i]->vt)
                    rng.Dirichlet(Q + i * K, Alpha, Mi + i * K, K);
        }
        binaryq = false;
    }

    /* Update a priori ancetral proportion by Metropolis-Hastings for admix model */
    TARGET void BAYESIAN::UpdateQMetro()
    {
        if (STRUCTURE_DECOMP)
        {
            double* q0 = Q, * b0 = bufNK1;
            for (int i = 0; i < nind; ++i, b0 += K)
            {
                if (ainds[i]->vt == 0) continue;
                if (locpriori) rng.Dirichlet(b0, AlphaLocal + ainds[i]->popid * K, K);
                else           rng.Dirichlet(b0, Alpha, K);
            }

            switch (SIMD_TYPE)
            {
#ifdef __aarch64__
            case 2: UpdateQMetroNEO(Q, K, bufNK1, bufN1, bufN2, Base); break;
#else
            case 4: UpdateQMetro512(Q, K, bufNK1, bufN1, bufN2, Base); break;
            case 3: UpdateQMetroAVX(Q, K, bufNK1, bufN1, bufN2, Base); break;
            case 2: UpdateQMetroSSE(Q, K, bufNK1, bufN1, bufN2, Base); break;
#endif
            }

            b0 = bufNK1; q0 = Q;
            for (int i = 0; i < N; ++i, q0 += K, b0 += K)
            {
                if (ainds[i]->vt == 0) continue;
                if (bufN1[i] >= NZERO || rng.Uniform() < exp(bufN1[i]))
                    SetVal(q0, b0, K);
            }

            binaryq = false;
            return;
        }

        //add priori probability
        double* bufi = bufNK1;
        for (int i = 0; i < N; ++i, bufi += K)
        {
            if (ainds[i]->vt == 0) continue;
            if (locpriori) rng.Dirichlet(bufi, AlphaLocal + ainds[i]->popid * K, K);
            else           rng.Dirichlet(bufi, Alpha, K);
        }

        OpenLog((int64*)bufN1, bufN2, N);
        double* p = Base, * q = NULL;
    
        for (int64 l = 0; l < L; ++l, p += GetLoc(l).k)
        {
            GENO_READER rt(0, l);
            GENOTYPE* gtab = GetLoc(l).GetGtab();
            bufi = bufNK1; q = Q;
        
            for (int i = 0; i < N; ++i, q += K, bufi += K)
            {
                GENOTYPE& gt = gtab[rt.Read()];
                if (gt.Nalleles() == 0) continue;

                ushort* als = gt.GetAlleleArray();
                for (int a = 0, vi = gt.Ploidy(); a < vi; ++a)
                    ChargeLog(*(int64*)&bufN1[i], bufN2[i], SumProdDiv(bufi, q, p + als[a], KT, K));
            }
        }
        CloseLog((int64*)bufN1, bufN2, N);

        bufi = bufNK1; q = Q;
        for (int i = 0; i < N; ++i, q += K, bufi += K)
        {
            if (ainds[i]->vt == 0) continue;
            if (bufN1[i] >= NZERO || rng.Uniform() < exp(bufN1[i]))
                SetVal(q, bufi, K);
        }

        binaryq = false;
    }

    /* Update a priori ancetral proportion */
    TARGET void BAYESIAN::UpdateQ()
    {
        //draw gene proportion for each individual
        if (!admix)
        {
            if (m >= nadmburnin || locpriori)
                UpdateQNoAdmix();
            else
                UpdateQAdmix();

            if (locpriori)
            {
                SetZero(Di, S * K);
                double* q = Q;
                for (int i = 0; i < N; ++i, q += K)
                    if (ainds[i]->vt) //bug fixed for individuals with no data
                        Add(Di + ainds[i]->popid * K, q, K);
            }
        }
        else
        {
            if (metrofreq > 0 && m % metrofreq == 0)
                UpdateQMetro();
            else
                UpdateQAdmix();
        }
    }

    /* Update locpriori parameters */
    TARGET void BAYESIAN::UpdateLocPriori()
    {
        if (!locpriori) return;
        if (admix)
        {
            double rm = rng.Uniform(r - epsr, r + epsr);
            if (rm > 0 && rm < maxr)
            {
                double dlnL = 0, d = rm - r, dt1 = rm * MyLog(rm) - r * MyLog(r);
                for (int k = 0; k < K; ++k)
                {
                    dlnL += S * (Alpha[k] * dt1 - LogGamma1(rm * Alpha[k]) + LogGamma1(r * Alpha[k]));
                    dlnL += Alpha[k] * d * LogProd(AlphaLocal + k, S, K) - Sum(AlphaLocal + k, S, K) * d; 
                }
                if (dlnL >= NZERO || rng.Uniform() < exp(dlnL))
                    r = rm;
            }
        }
        else
        {
            double rm = rng.Uniform(r - epsr, r + epsr);
            if (rm > 0 && rm < maxr)
            {
                double dlnL = LogGamma1(rm) - LogGamma1(r);
                for (int k = 0; k < K; ++k)
                    dlnL += LogGamma1(r * Eta[k]) - LogGamma1(rm * Eta[k]);
                dlnL *= S;

                for (int k = 0; k < K; ++k)
                    dlnL += (rm - r) * Eta[k] * LogProd(Gamma + k, S, K);

                if (dlnL >= NZERO || rng.Uniform() < exp(dlnL))
                    r = rm;
            }

            //update eta
            if (K >= 2)
            {
                int64 i1 = rng.Next(K), j1 = rng.NextAvoid(K, i1);
                double delta = rng.Uniform(epseta);
                double e1 = Eta[i1] + delta;
                double e2 = Eta[j1] - delta;
                if (e1 < 1 && e2 > 0)
                {
                    double dlnL = S * (LogGamma1(r * Eta[i1]) + LogGamma1(r * Eta[j1]) - LogGamma1(r * e1) - LogGamma1(r * e2));
                    dlnL += r * delta * LogProdDiv(Gamma + i1, Gamma + j1, S, K);

                    if (dlnL >= NZERO || rng.Uniform() < exp(dlnL))
                    {
                        Eta[i1] += delta;
                        Eta[j1] -= delta;
                    }
                }
            }

            //update gamma
            if (K >= 2)
            {
                for (int s = 0; s < S; ++s)
                {
                    int64 i1 = rng.Next(K), j1 = rng.NextAvoid(K, i1);
                    double delta = rng.Uniform(epsgamma);

                    if (Gamma[s * K + i1] + delta < 1 && Gamma[s * K + j1] - delta > 0)
                    {
                        double dlnL = (r * Eta[i1] - 1.0 + Di[s * K + i1]) * MyLog(1 + delta / Gamma[s * K + i1]) +
                            (r * Eta[j1] - 1.0 + Di[s * K + j1]) * MyLog(1 - delta / Gamma[s * K + j1]);

                        if (dlnL >= NZERO || rng.Uniform() < exp(dlnL))
                        {
                            Gamma[s * K + i1] += delta;
                            Gamma[s * K + j1] -= delta;
                        }
                    }
                }
            }

            if (m % 10 == 0)
            {
                Unify(Eta, K);
                Unify(Gamma, S, K);
            }
        }
    }

    /* Update ancestral proportion for each allele or for each individual */
    TARGET void BAYESIAN::UpdateZ()
    {
        //Update Z, Mi, Ni
        SetZero(Mi, N * K);
        SetZero(Ni, K * KT);

        if (binaryq)
        {
            // Z in already updated in update Q
            if (STRUCTURE_DECOMP)
            {
                //Update Mi, count number of allele copies in individual i and cluster k
                for (int i = 0; i < nind; ++i)
                    Mi[i * K + Z[i]] = ainds[i]->vt;

                //update Ni, count the number of different alleles in each cluster
                switch (SIMD_TYPE)
                {
#ifdef __aarch64__
                case 2: UpdateZBinaryNEO(Ni, Z); break;
#else
                case 4: UpdateZBinary512(Ni, Z); break;
                case 3: UpdateZBinaryAVX(Ni, Z); break;
                case 2: UpdateZBinarySSE(Ni, Z); break;
#endif
                }

                return;
            }

            //Count number of allele copies in individual i and cluster k
            for (int i = 0; i < N; ++i)
                Mi[i * K + Z[i]] = ainds[i]->vt;

            //Count number of allele k2 in cluster k
            for (int64 l = 0, o = 0; l < L; ++l, o += GetLoc(l).k)
            {
                GENO_READER rt(0, l);
                GENOTYPE* gtab = GetLoc(l).GetGtab();

                for (int i = 0; i < N; ++i)
                {
                    GENOTYPE& gt = gtab[rt.Read()];
                    if (gt.Nalleles() == 0) continue;

                    ushort* als = gt.GetAlleleArray();
                    int* ni = Ni + Z[i] * KT + o;
                    for (int a = 0, vi = gt.Ploidy(); a < vi; ++a)
                        ni[als[a]]++;
                }
            }
        }
        else
        {
            if (STRUCTURE_DECOMP)
            {
                //draw Z according to Q and freq, and update Mi and Ni
                switch (SIMD_TYPE)
                {
#ifdef __aarch64__
                case 2: UpdateZAdmixNEO(Q, K, Mi, Ni, bufNK1, bufNK2, Base, rngNEO); break;
#else
                case 4: UpdateZAdmix512(Q, K, Mi, Ni, bufNK1, bufNK2, Base, rng512); break;
                case 3: UpdateZAdmixAVX(Q, K, Mi, Ni, bufNK1, bufNK2, Base, rngAVX); break;
                case 2: UpdateZAdmixSSE(Q, K, Mi, Ni, bufNK1, bufNK2, Base, rngSSE); break;
#endif
                }

                return;
            }

            for (int64 l = 0, o = 0; l < L; ++l, o += GetLoc(l).k)
            {
                GENO_READER rt(0, l);
                GENOTYPE* gtab = GetLoc(l).GetGtab();

                for (int i = 0; i < N; ++i)
                {
                    GENOTYPE& gt = gtab[rt.Read()];
                    if (gt.Nalleles() == 0) continue;

                    double* q = Q + i * K;
                    ushort* als = gt.GetAlleleArray();
                    int64* mi = Mi + i * K;
                    int* ni = Ni + o;

                    for (int a = 0, vi = gt.Ploidy(); a < vi; ++a)
                    {
                        //The same to previous allele, do not update multinomial prob to save time
                        if (a == 0 || als[a] != als[a - 1])
                            for (int k = 0; k < K; ++k)
                                bufNK1[k] = q[k] * cluster[k].GetFreq(l, als[a]);

                        //draw cluster for each allele copy
                        ushort k2 = (ushort)rng.Poly(bufNK1, K);

                        //Update Mi, NI
                        mi[k2]++;
                        ni[k2 * KT + als[a]]++;
                    }
                }
            }
        }
    }

    /* Update Dirichlet parameter alpha (to draw admixture proportion Q) in the admix model */
    TARGET void BAYESIAN::UpdateAlpha()
    {
        if (locpriori && admix)
        {
            //UpdateAlphaLocPrior
            //update AlphaGlobal
            for (int k = 0; k < K; ++k)
            {
                double am = rng.Normal(Alpha[k], stdalpha);
                if (am > 0 && am < maxalpha)
                {
                    double d = am - Alpha[k];
                    double dlnL = S * (d * r * MyLog(r) + LogGamma1(r * Alpha[k]) - LogGamma1(r * am));

                    dlnL += r * d * LogProd(AlphaLocal + k, S, K);

                    if (dlnL >= NZERO || rng.Uniform() < exp(dlnL))
                        Alpha[k] = am;
                }
            }

            //update SumAlpha
            if (m % 10 == 0) 
                for (int s = 0; s < S; ++s)
                    SumAlpha[s] = Sum(&AlphaLocal[s * K], K);

            //update AlphaLocal
            for (int s = 0; s < S; ++s)
            {
                int snind = apops[s]->nind;
                double* sq = &Q[apops[s]->ind0id * K + 0];

                for (int k = 0; k < K; ++k)
                {
                    double ao = AlphaLocal[s * K + k], am = rng.Normal(ao, stdalpha);

                    if (am > 0 && am < maxalpha)
                    {
                        double d = am - ao;
                        double dlnL = (r * Alpha[k] - 1) * MyLog(am / ao) - d * r +
                                      apops[s]->nind * (LogGamma1(SumAlpha[s] + d) -
                                      LogGamma1(am) - LogGamma1(SumAlpha[s]) + LogGamma1(ao));

                        dlnL += d * LogProd(sq + k, snind, K);

                        if (dlnL >= NZERO || rng.Uniform() < exp(dlnL))
                        {
                            AlphaLocal[s * K + k] = am;
                            SumAlpha[s] += d;
                        }
                    }
                }
            }
        }
        else if (inferalpha)
        {
            if (!admix && m >= nadmburnin) return;
            //for admix and not locpriori

            if (diffalpha)
            {
                double sumalpha = 0;
                for (int k = 0; k < K; ++k)
                    sumalpha += Alpha[k];

                for (int k = 0; k < K; ++k)
                {
                    double am = rng.Normal(Alpha[k], stdalpha);
                    if (am <= 0 || (am >= maxalpha && uniformalpha)) continue;
                    double ao = Alpha[k];
                    double d = am - ao;
                    double dlnL = uniformalpha ? 0 : ((alphapriora - 1) * MyLog(am / ao) + (ao - am) / alphapriorb);

                    dlnL += (am - ao) * LogProd(Q + k, N, K) + (LogGamma1(sumalpha + d) - LogGamma1(sumalpha) - LogGamma1(am) + LogGamma1(ao)) * N;
                    if (dlnL >= NZERO || rng.Uniform() < exp(dlnL))
                        Alpha[k] = am;
                }
            }
            else
            {
                double am = rng.Normal(Alpha[0], stdalpha);
                if (am <= 0 || (am >= maxalpha && uniformalpha)) return;
                double ao = Alpha[0];
                double dlnL = uniformalpha ? 0 : ((alphapriora - 1) * MyLog(am / ao) + (ao - am) / alphapriorb);

                dlnL += (am - ao) * LogProd(Q, N * K) + (LogGamma1(K * am) - LogGamma1(K * ao) - K * LogGamma1(am) + K * LogGamma1(ao)) * N;
                if (dlnL >= NZERO || rng.Uniform() < exp(dlnL))
                    for (int k = 0; k < K; ++k)
                        Alpha[k] = am;
            }
        }
    }

    /* Update Dirichlet parameter lambda (to draw allele frequency) */
    TARGET void BAYESIAN::UpdateLambda()
    {
        if (!inferlambda) return;
        if (difflambda) for (int k = 0; k < K; ++k)
        {
            double lo = Lambda[k];
            double lm = rng.Normal(lo, stdlambda);
            if (lm <= 0 || lm >= maxlambda) continue;
            double dlnL = 0;
            double gd = LogGamma1(lo) - LogGamma1(lm);

            for (int i = 2; i <= maxK; ++i)
                dlnL += kdis[i] * (LogGamma1(i * lm) - LogGamma1(i * lo) + i * gd);

            dlnL += (lm - lo) * LogProd(cluster[k].bucket, KT);

            if (dlnL >= NZERO || rng.Uniform() < exp(dlnL))
                Lambda[k] = lm;
        }
        else
        {
            int np = fmodel ? 1 : K;
            double lo = Lambda[0];
            double lm = rng.Normal(lo, stdlambda);
            if (lm <= 0 || lm >= maxlambda) return;
            double dlnL = 0;
            double gd = LogGamma1(lo) - LogGamma1(lm);

            for (int i = 2; i <= maxK; ++i)
                dlnL += kdis[i] * np * (LogGamma1(i * lm) - LogGamma1(i * lo) + i * gd);

            dlnL += (lm - lo) * (fmodel ? LogProd(PA.bucket, KT) : LogProd(Base, K * KT));

            if (dlnL >= NZERO || rng.Uniform() < exp(dlnL))
                SetVal(Lambda, lm, K);
        }
    }

    /* Update allele frequency of ancestral population for the F model */
    TARGET void BAYESIAN::UpdatePA()
    {
        if (!fmodel) return;
        double lam = Lambda[0];//single lambda in F model?

        if (rng.Uniform() < 0.5)
        {
            //Independent update PA
            SetVal(PA1.bucket, lam, KT);

            for (int k = 0; k < K; ++k)
                AddProd(PA1.bucket, cluster[k].bucket, f[k], KT);

            double* pa = PA.bucket, * t1 = PA1.bucket, * t2 = PA2.bucket;

            //unable to optimize
            for (int k = 0; k < K; ++k)
                bufK[k] = cluster[k].bucket;

            for (int64 l = 0; l < L; ++l)
            {
                int k2 = GetLoc(l).k;
                rng.Dirichlet(t2, t1, k2);

                double dlnL = 0;
                for (int a = 0; a < k2; ++a)
                    dlnL += (t1[a] - lam) * (MyLog(pa[a] / t2[a]));

                for (int k = 0; k < K; ++k)
                {
                    double fr = f[k];
                    for (int a = 0; a < k2; ++a)
                        dlnL += LogGamma1(fr * pa[a]) - LogGamma1(fr * t2[a])
                        + fr * (t2[a] - pa[a]) * MyLog(*bufK[k]++);
                }

                if (dlnL >= NZERO || rng.Uniform() < exp(dlnL))
                    SetVal(pa, t2, k2);

                pa += k2; t1 += k2; t2 += k2;
            }
        }
        else
        {
            double* pa = PA.bucket, * t1 = PA1.bucket, * t2 = PA2.bucket;
            double* p = Base;

            for (int64 l = 0; l < L; ++l)
            {
                int k2 = GetLoc(l).k;
                int m1 = rng.Next(k2), n1 = rng.NextAvoid(k2, m1);

                double delta = rng.Uniform(pow(N, -0.5));
                double pm0 = pa[m1], pm1 = pm0 + delta, pn0 = pa[n1], pn1 = pn0 - delta;
                if (pm1 >= 1 || pn1 <= 0)
                {
                    pa += k2; t1 += k2; t2 += k2;
                    continue;
                }
                double dlnL = (lam - 1) * MyLog(pm1 * pn1 / pm0 / pn0);

                for (int k = 0; k < K; ++k)
                {
                    double fr = f[k];
                    dlnL += LogGamma1(fr * pm0) + LogGamma1(fr * pn0)
                        - LogGamma1(fr * pm1) - LogGamma1(fr * pn1) +
                        (fr * delta) * MyLog(p[k * KT + m1] / p[k * KT + n1]);
                    p += k2;
                }

                if (dlnL >= NZERO || rng.Uniform() < exp(dlnL))
                {
                    pa[m1] = pm1;
                    pa[n1] = pn1;
                }

                pa += k2; t1 += k2; t2 += k2;
            }
        }
    }

    /* Update population-specific Fst for the F model */
    TARGET void BAYESIAN::UpdateF()
    {
        //F model
        if (!fmodel) return;

        // Update F/Fst
        for (int k = 0; k < K; ++k)
        {
            double Fo = F[k];
            double Fm = rng.Normal(Fo, stdf);
            if (Fm < 0 || Fm > 1)
            {
                if (fsame) break;
                else continue;
            }
            double fo = f[k], fm = (1 - Fm) / Fm;
            double dlnL = pmeanf * (Fo - Fm) / (pstdf * pstdf) +
                (pmeanf * pmeanf / (pstdf * pstdf) - 1) * MyLog(Fm / Fo) +
                (fsame ? K : 1) * nloc * (LogGamma1(fm) - LogGamma1(fo));//eL add in 2018.09

            for (int k2 = k; k2 < K; ++k2)
            {
                double* p = cluster[k2].bucket;
                double* pa = PA.bucket;//2018.09
                for (int a = 0; a < KT; ++a)
                {
                    double pa2 = *pa++;
                    dlnL += pa2 * (fm - fo) * MyLog(*p++) + LogGamma1(fo * pa2) - LogGamma1(fm * pa2);
                }
                if (!fsame) break;
            }

            if (dlnL >= NZERO || rng.Uniform() < exp(dlnL))
            {
                if (fsame)
                {
                    SetVal(F, Fm, K);
                    SetVal(f, fm, K);
                }
                else
                {
                    F[k] = Fm;
                    f[k] = fm;
                }
            }

            if (fsame) break;
        }
    }

    /* Finalize records */
    TARGET void BAYESIAN::Arrange()
    {
        UnifyInt64ToDouble(MiSum, N, K);//cast to double
        Mul(rout, 1.0 / nr, rlen);
        rout[2] = rout[1] - rout[0] * rout[0];
        rout[3] = rout[0] - rout[2] / 2;
        Mul(Base + K * KT, 1.0 / nr, K * KT);
    }

    /* Record updated MCMC parameters */
    TARGET void BAYESIAN::Record()
    {
        if (m >= nburnin && (m - nburnin) % nthinning == 0)
        {
            nr++;
            Add(MiSum, Mi, N * K);
            int64 slog = 0; double prod = 1;
            double* p = Base;

            OpenLog(slog, prod);
            for (int64 l = 0; l < L; ++l)
            {
                GENO_READER rt(0, l);
                GENOTYPE* gtab = GetLoc(l).GetGtab();
                double* q = Q;

                for (int i = 0; i < N; ++i, q += K)
                {
                    GENOTYPE& gt = gtab[rt.Read()];
                    if (gt.Nalleles() == 0) continue;

                    ushort* als = gt.GetAlleleArray();
                    for (int a = 0, vi = gt.Ploidy(); a < vi; ++a)
                        ChargeLog(slog, prod, SumProd(q, p + als[a], KT, K));
                }
                p += GetLoc(l).k;
            }
            CloseLog(slog, prod);

            //add to result
            rout[0] += prod;
            rout[1] += prod * prod;

            int nf = fmodel ? (fsame ? 1 : K) : 0;
            int na = admix ? (diffalpha ? K : 1) : 0;
            int nl = difflambda ? K : 1;

            Add(rout + 4, Lambda, nl);

            if (locpriori)
            {
                rout[4 + nl] += r;
                Add(rout + 4 + nl + 1, admix ? Alpha : Eta, K);
                Add(rout + 4 + nl + 1 + K, admix ? AlphaLocal : Gamma, S * K);
            }
            else if (admix)
                Add(rout + 4 + nl, Alpha, na);
            Add(rout + rlen - nf, F, nf);
            Add(Base + K * KT, Base, K * KT);
        }
        PROGRESS_VALUE++;
    }

    /* Free memory */
    TARGET void BAYESIAN::Uninit()
    {
        //Normal
        TryDelete(MiSum);
        TryDelete(rout);
        TryDelete(cluster);
        TryDelete(clusterb);
        TryDelete(Base);
        TryDelete(Lambda);
        TryDelete(bufK);
        TryDelete(bufN1);
        TryDelete(bufN2);

        //if (bufNK1) operator delete[](bufNK1, sizeof(double) * N * K, (std::align_val_t)(64));
        //if (bufNK2) operator delete[](bufNK2, sizeof(double) * N * K, (std::align_val_t)(64));
        if (bufNK1) { free(bufNK1o); bufNK1 = NULL; bufNK1o = NULL; }
        if (bufNK2) { free(bufNK2o); bufNK2 = NULL; bufNK2o = NULL; }
        bufNK1 = NULL; bufNK1o = NULL; bufNK2 = NULL; bufNK2o = NULL;

        //Adm
        TryDelete(Z);
        TryDelete(Q);
        TryDelete(Mi);
        TryDelete(Ni);

        //Fmodel
        TryDelete(F);
        TryDelete(f);

        //Loc
        TryDelete(Alpha);
        TryDelete(SumAlpha);
        TryDelete(AlphaLocal);
        TryDelete(Di);
        TryDelete(Eta);
        TryDelete(Gamma);
        TryDelete(kdis);
    }

    /* Perform MCMC */
    TARGET void BAYESIAN::MCMC()
    {
        InitFix();
        InitAdmix();
        InitLocPriori();
        InitFmodel();

        nadmburnin = (!locpriori && !admix) ? nadmburnin : 0;
        nburnin += nadmburnin;
        int ntrep = nburnin + nreps;
        uint64 seed = (g_seed_val + par2->id);

        for (m = 0; m < ntrep; ++m)
        {
            if ((m & 0xFFFF) == 0)
            {
                //reset seeds
                uint s = HashULong(seed + m);
                new(&rng) RNG(s);
                switch (SIMD_TYPE)
                {
#ifdef __aarch64__
                case 2: new(&rngNEO) RNGNEO(s); break;
#else
                case 4: new(&rng512) RNG512(s); break;
                case 3: new(&rngAVX) RNGAVX(s); break;
                case 2: new(&rngSSE) RNGSSE(s); break;
#endif
                }
            }

            binaryq = false;
            UpdateP();//Allele frequency
            UpdateQ();//Individual priori gene proportion
            UpdateLocPriori();//LocPriori r
            UpdateZ();//Individual gene origin
            UpdateAlpha();//LocPriori+Admix
            UpdateLambda();//
            UpdatePA();//
            UpdateF();//
            Record();
        }
        Arrange();
    }

#endif

#ifndef _STCLUSTER
    /* Get allele frequency array */
    TARGET double* SCLUSTER::GetFreq(int64 l)
    {
        return bucket + allele_freq_offset[l];
    }

    /* Get allele frequency */
    TARGET double SCLUSTER::GetFreq(int64 l, int allele)
    {
        return bucket[allele_freq_offset[l] + allele];
    }

    /* Set allele frequency pointer */
    TARGET void SCLUSTER::SetFreq(double* _bucket)
    {
        bucket = _bucket;
    }
#endif

#define extern
extern STRUCTURE_RUNINFO* structure_par;                        //Genetic distance for pcoa and hierarchical clustering
extern int structure_totalruns;                                    //Total number of runs in Bayesian clustering
extern uint64** structure_allele;                                //compressed allele array for each individual
extern byte* structure_size;                                    //allele data size (in bits) for each locus
extern int64 structure_indnbytes;                                //number of bytes used in decompoess model of structure
#undef extern

/* Calculate bayesian clustering */
TARGET void CalcBayesian()
{
    if (!structure) return;
    if (ad) Exit("\nError: Bayesian clustering (-structure) is incompatible with allelic depth (-ad) option.\n");

    EvaluationBegin();

    //alloc tasks
    structure_totalruns = (structure_krange_max + 1 - structure_krange_min) * structure_nruns_val;
    structure_par = new STRUCTURE_RUNINFO[structure_totalruns];
    int lc = 0;

    STRUCTURE_DECOMP = structure_totalruns && structure_decompress_val == 1 && SIMD_TYPE >= 2 && nind > STRUCTURE_NPACK; //minploidy == maxploidy

    //decompressed model, homoploid only
    if (STRUCTURE_DECOMP)
    {
        structure_size = new byte[nloc];
        uint64 haploid_nbits = 0;
        for (int64 l = 0; l < nloc; ++l)
        {
            structure_size[l] = (byte)CeilLog2((int64)GetLoc(l).k + 1);
            haploid_nbits += structure_size[l];
        }

        //individual memory
        structure_indnbytes = (7 + haploid_nbits * maxploidy) >> 3;
        structure_indnbytes = ((structure_indnbytes + 7) >> 3) << 3; //align 8 bytes
        structure_allele = new uint64*[nind];
        //structure_allele[0] = (uint64*) operator new[](sizeof(uint64)* structure_indnbytes* nind, (std::align_val_t)(8));
        structure_allele[0] = new uint64[structure_indnbytes * nind];
        
        for (int i = 1; i < nind; ++i)
            structure_allele[i] = structure_allele[i - 1] + structure_indnbytes;

        //decompress genotype
        int npack = 8, nloop = (nind + 7) / npack;

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
        for (int j = 0; j < nloop; ++j)
        {
            int i = j * npack;

            if (i + 7 >= nind)
                i = nind - npack;

            BAYESIAN_WRITER wt(i);
            IND* ii[8] = { ainds[i + 0], ainds[i + 1], ainds[i + 2], ainds[i + 3],
                           ainds[i + 4], ainds[i + 5], ainds[i + 6], ainds[i + 7] };

            for (int64 l = 0; l < nloc; ++l)
            {
                int size = structure_size[l];

                GENOTYPE& gt1 = ii[0]->GetGenotype(l), &gt2 = ii[1]->GetGenotype(l);
                GENOTYPE& gt3 = ii[2]->GetGenotype(l), &gt4 = ii[3]->GetGenotype(l);
                GENOTYPE& gt5 = ii[4]->GetGenotype(l), &gt6 = ii[5]->GetGenotype(l);
                GENOTYPE& gt7 = ii[6]->GetGenotype(l), &gt8 = ii[7]->GetGenotype(l);

                int v1 = gt1.Ploidy(), v2 = gt2.Ploidy(), v3 = gt3.Ploidy(), v4 = gt4.Ploidy();
                int v5 = gt5.Ploidy(), v6 = gt6.Ploidy(), v7 = gt7.Ploidy(), v8 = gt8.Ploidy();

                ushort* als[8] = {
                    gt1.GetAlleleArray(), gt2.GetAlleleArray(), gt3.GetAlleleArray(), gt4.GetAlleleArray(),
                    gt5.GetAlleleArray(), gt6.GetAlleleArray(), gt7.GetAlleleArray(), gt8.GetAlleleArray(),
                };

                for (int a = 0; a < maxploidy; ++a)
                {
                    uint64 aid[8] = { 
                        (uint64)(a < v1 ? als[0][a] : 0xFFFFu), (uint64)(a < v2 ? als[1][a] : 0xFFFFu),
                        (uint64)(a < v3 ? als[2][a] : 0xFFFFu), (uint64)(a < v4 ? als[3][a] : 0xFFFFu),
                        (uint64)(a < v5 ? als[4][a] : 0xFFFFu), (uint64)(a < v6 ? als[5][a] : 0xFFFFu),
                        (uint64)(a < v7 ? als[6][a] : 0xFFFFu), (uint64)(a < v8 ? als[7][a] : 0xFFFFu) };

                    wt.Write8(size, aid);
                }
            }

            wt.FinishWrite8();
        }
    }

    //time test xxx
    //structure_nburnin_val = 10;
    //structure_nreps_val = 10;

    for (int i = structure_krange_min; i <= structure_krange_max; ++i)
        for (int j = 1; j <= structure_nruns_val; ++j)
        {
            structure_par[lc].flag.clear();
            structure_par[lc].k = i;
            structure_par[lc].id = lc + 1;
            structure_par[lc].rep = j;
            structure_par[lc].MeanlnL = 0;
            structure_par[lc].VarlnL = 0;
            structure_par[lc].lnPD = 0;
            lc++;
        }

    int nrun = ((structure_locpriori_val == 2 && structure_admix_val == 2) ? structure_nadmburnin_val : 0) + structure_nreps_val + structure_nburnin_val;
    RunThreads(&BayesianThread, NULL, NULL, structure_totalruns * nrun, structure_totalruns * nrun,
        "\nPerforming bayesian clustering:\n", g_nthread_val, true);
    BAYESIAN::PrintSummary(structure_par, structure_totalruns);
    delete[] structure_par;

    //decompressed model
    if (STRUCTURE_DECOMP)
    {
        //operator delete[](structure_allele[0], sizeof(uint64)* structure_indnbytes * nind, (std::align_val_t)(8));
        delete[] structure_allele[0];
        delete[] structure_allele;
        delete[] structure_size;
    }

    EvaluationEnd("Bayesian clustering");

    if (structure_plot_val == 1)
        RunRscript("structure_plot.R");
}

/* Calculate Bayesian clustering using multiple threads */
THREAD(BayesianThread)
{
    for (int i = 0; i < structure_totalruns; ++i)
    {
        if (structure_par[i].flag.test_and_set()) continue;

        BAYESIAN Structure;
        Structure.ReadPar(structure_par + i);
        Structure.MCMC();
        Structure.PrintStructure();
    }
}
