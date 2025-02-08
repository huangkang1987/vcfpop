/* Bayesian clustering functions */

#include "vcfpop.h"

#pragma pack(push, 1)
#define TIME_TEST

template struct SCLUSTER<double>;
template struct SCLUSTER<float >;
template struct BAYESIAN<double>;
template struct BAYESIAN<float >;

template TARGET void CalcBayesian<double>();
template TARGET void CalcBayesian<float >();

template TARGET void BAYESIAN<double>::UpdateQMetroCPU<true >(int tid);
template TARGET void BAYESIAN<double>::UpdateQMetroCPU<false>(int tid);
template TARGET void BAYESIAN<float >::UpdateQMetroCPU<true >(int tid);
template TARGET void BAYESIAN<float >::UpdateQMetroCPU<false>(int tid);

template TARGET void BAYESIAN<double>::UpdateQNoAdmixCPU(int tid);
template TARGET void BAYESIAN<float >::UpdateQNoAdmixCPU(int tid);

template TARGET void BAYESIAN<double>::UpdateZAdmixCPU<true >(int tid);
template TARGET void BAYESIAN<double>::UpdateZAdmixCPU<false>(int tid);
template TARGET void BAYESIAN<float >::UpdateZAdmixCPU<true >(int tid);
template TARGET void BAYESIAN<float >::UpdateZAdmixCPU<false>(int tid);

template TARGET void BAYESIAN<double>::UpdateZNoAdmixCPU(int tid);
template TARGET void BAYESIAN<float >::UpdateZNoAdmixCPU(int tid);

template TARGET void BAYESIAN<double>::InitAdmix(int tid);
template TARGET void BAYESIAN<float >::InitAdmix(int tid);

template TARGET void BAYESIAN<double>::RecordCPU<true , true >(int tid);
template TARGET void BAYESIAN<double>::RecordCPU<true , false>(int tid);
template TARGET void BAYESIAN<double>::RecordCPU<false, true >(int tid);
template TARGET void BAYESIAN<double>::RecordCPU<false, false>(int tid);
template TARGET void BAYESIAN<float >::RecordCPU<true , true >(int tid);
template TARGET void BAYESIAN<float >::RecordCPU<true , false>(int tid);
template TARGET void BAYESIAN<float >::RecordCPU<false, true >(int tid);
template TARGET void BAYESIAN<float >::RecordCPU<false, false>(int tid);

template TARGET void BAYESIAN<double>::UpdateQMetro();
template TARGET void BAYESIAN<float >::UpdateQMetro();

#ifndef _SCLUSTER

/* Get allele frequency array */
template<typename REAL>
TARGET REAL* SCLUSTER<REAL>::GetFreq(int64 l)
{
    return bucket + allele_freq_offset[l];
}

/* Get allele frequency */
template<typename REAL>
TARGET REAL SCLUSTER<REAL>::GetFreq(int64 l, int allele)
{
    return bucket[allele_freq_offset[l] + allele];
}

/* Set allele frequency pointer */
template<typename REAL>
TARGET void SCLUSTER<REAL>::SetFreq(REAL* _bucket)
{
    bucket = _bucket;
}
#endif

#ifndef _BAYESIAN

#ifdef __aarch64__
const char* simdstr[] = { "", "CPU" , "NEO",  "GPU"};
#else
const char* simdstr[] = { "", "CPU" , "SSE" , "AVX" , "512",  "GPU"};
#endif

/* Set all bits to 0 */
template<typename REAL>
TARGET BAYESIAN<REAL>::BAYESIAN()
{
    SetZero(this, 1);
}

/* Write results for a run */
template<typename REAL>
TARGET void BAYESIAN<REAL>::WriteStructure()
{
    char filename[PATH_LEN];
    char name_buf[NAME_BUF_LEN];
    FILE* fout;

    //write lnL
    if (structure_writelnl_val == 1)
    {
        sprintf(filename, "%s.structure.k=%d_rep=%d_id=%d.lnl.txt", g_output_val.c_str(), K, par2->rep, par2->id);
        fout = fopen(filename, "wb");
        fprintf(fout, "%d%s", nburnin, g_linebreak_val);
        fprintf(fout, "Iteration%clnL%s", g_delimiter_val, g_linebreak_val);
        double* lnl = lnL;
        for (int i = 0; i < nburnin; i += nthinning)
        {
            fprintf(fout, "%d%c", i, g_delimiter_val);
            WriteReal(fout, *lnl++);
            fprintf(fout, "%s", g_linebreak_val);
        }
        for (int i = 0; i < nreps; i += nthinning)
        {
            fprintf(fout, "%d%c", i + nburnin, g_delimiter_val);
            WriteReal(fout, *lnl++);
            fprintf(fout, "%s", g_linebreak_val);
        }
        fclose(fout);
    }

    sprintf(filename, "%s.structure.k=%d_rep=%d_id=%d.txt", g_output_val.c_str(), K, par2->rep, par2->id);
    fout = fopen(filename, "wb");

    fprintf(fout, "%s%sParameters:%s", g_linebreak_val, g_linebreak_val, g_linebreak_val);
    fprintf(fout, "Seed=%lld%s", seed, g_linebreak_val);
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
            fprintf(fout, "%s%s", g_linebreak_val, apops<REAL>[i]->name);
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

    VLA_NEW(O, REAL, S * K);
    SetZero(O, S * K);
    fprintf(fout, "%s%sProportion of membership of each pre-defined population in each of the %d clusters", g_linebreak_val, g_linebreak_val, K);
    fprintf(fout, "%s%cCluster%sPop", g_linebreak_val, g_delimiter_val, g_linebreak_val);
    for (int i = 0; i < N; ++i)
        Add(O + ainds<REAL>[i]->popid * K, Q + i * K, K);
    Unify(O, S, K);

    for (int k = 0; k < K; ++k)
        fprintf(fout, "%c%d", g_delimiter_val, k + 1);

    for (int s = 0; s < S; ++s)
    {
        fprintf(fout, "%s%s", g_linebreak_val, apops<REAL>[s]->name);
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
        fprintf(fout, "%s%s%c%s%c%s", g_linebreak_val, ainds<REAL>[i]->name,
            g_delimiter_val, apops<REAL>[ainds<REAL>[i]->popid]->name,
            g_delimiter_val, apops<REAL>[ainds<REAL>[i]->popid]->rid == -1 ? "Total" : aregs<REAL>[0][apops<REAL>[ainds<REAL>[i]->popid]->rid]->name
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
        REAL* p = ClusterMeanFreq(k);
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
            REAL* p1 = ClusterMeanFreq(k), *p2 = ClusterMeanFreq(j);
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

        REAL** bufK = (REAL**)bufNK1;
        for (int k = 0; k < K; ++k)
        {
            fprintf(fout, "%s%s%c%d", g_linebreak_val, k ? "" : "Average", g_delimiter_val, k + 1);
            WriteReal(fout, H[k]);
            bufK[k] = ClusterMeanFreq(k);
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
    DEL(H2);

    fclose(fout);
    Uninit();
}

/* Write results summary for all runs */
template<typename REAL>
TARGET void BAYESIAN<REAL>::WriteStructureSummary(STRUCTURE_RUNINFO* sp, int len)
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
template<typename REAL>
TARGET void BAYESIAN<REAL>::ReadPar(STRUCTURE_RUNINFO* _par2)
{
    par2 = _par2;

    iscudastream = false;
    usecuda = false;

    id = par2->id;
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

    kdis = new int64[maxK + 1]; 
    SetZero(kdis, maxK + 1);
    for (int64 l = 0; l < L; ++l)
        kdis[GetLoc(l).k]++;

    nadmburnin = (!locpriori && !admix) ? structure_nadmburnin_val : 0;

    nburnin += nadmburnin;
    seed = HashULong(g_seed_val + par2->id);
}

/* Initialize MCMC */
template<typename REAL>
TARGET void BAYESIAN<REAL>::InitFix(int tid)
{
    if (tid == -1)
    {
        nr = 0;
        r = 1;

        Freq = (REAL*)MallocHostCUDA(K * KT * sizeof(REAL));
        SetZero(Freq, K * KT);

        FreqM = (REAL*)Malloc(K * KT * sizeof(REAL));
        SetZero(FreqM, K * KT);

        FreqA = (REAL*)Malloc(3 * KT * sizeof(REAL));
        SetZero(FreqA, 3 * KT);

        bufKthread = new double[structure_nsubthread * K]; 

        bufNK1o = (byte*)MallocHostCUDA((std::max(N, 64) * K * sizeof(double) + 63) * structure_nsubthread);  //align for SIMD 
        bufNK2o = (byte*)MallocHostCUDA((std::max(N, 64) * K * sizeof(double) + 63) * structure_nsubthread);  //align for SIMD 
        bufNK1 = (double*)Align64(bufNK1o);
        bufNK2 = (double*)Align64(bufNK2o);
        bufN1 = (double*)MallocHostCUDA(N * structure_nsubthread * sizeof(double));
        bufN2 = (double*)MallocHostCUDA(N * structure_nsubthread * sizeof(double));

        MiSum = new int64[N * K]; 
        SetZero(MiSum, N * K);

        SetVal(Freq, (REAL)0.5, K * KT);

        Ni = (int*)MallocHostCUDA(K * KT * sizeof(int));
        SetZero(Ni, K * KT);

        Mi = (int64*)MallocHostCUDA(N * K * structure_nsubthread * sizeof(int64));
        SetZero(Mi, N * K * structure_nsubthread);

        Q = (REAL*)MallocHostCUDA(N * K * sizeof(REAL));
        SetVal(Q, (REAL)(1.0 / K), N * K);

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

        if (structure_writelnl_val == 1)
        {
            lnL = new double[(nburnin / nthinning) + (nreps / nthinning) + 1];
            nlnL = 0;
        }

        //////////////////////////////////////////////////////////

        l_atomic[0] = 0;

        InitFix(0);

        return;
    }

#pragma omp parallel num_threads(structure_nsubthread)
    {
        for (int64 ll = l_atomic[0].fetch_add(1); ll < L; ll = l_atomic[0].fetch_add(1))
        {
            int64 l = ll * structure_loc_coprime % L;
            int k2 = GetLoc(l).k;
            if (k2 != 2)
            {
                REAL* p = Freq + allele_freq_offset[l], val = (REAL)(1.0 / k2);
                for (int k = 0; k < K; ++k, p += KT)
                    SetVal(p, val, k2);
                p += k2;
            }
        }
    }
}

/* Initialize MCMC for admix model */
template<typename REAL>
TARGET void BAYESIAN<REAL>::InitAdmix(int tid)
{
    // Init Z, Z, Mi and Ni

    if (tid == -1)
    {
        //Adm
        Z = (ushort*)MallocHostCUDA((N + 32) * sizeof(ushort)); //align for SIMD
        SetZero(Z, N);

        singlez = nadmburnin == 0 && !admix;

        if (singlez)
        {
            RNG<double> rng(seed - 1, RNG_SALT_INITADM);//double
            for (int i = 0; i < N; ++i)
            {
                Z[i] = (ushort)rng.Next(K);
                Mi[i * K + Z[i]] = ainds<REAL>[i]->vt;
            }
        }

        //////////////////////////////////////////////////////////

        l_atomic[0] = 0;

        InitAdmix(0);

        //avoid thread-conflict
        if (!singlez)
            for (int i = 1; i < structure_nsubthread; ++i)
                Add(Mi, Mi + N * K * i, N * K);
        return;
    }

    atomic<int> thread_counter = 0;
#pragma omp parallel num_threads(structure_nsubthread)
    {
        int tid2 = thread_counter.fetch_add(1);

        if (singlez)
        {
            for (int64 ll = l_atomic[0].fetch_add(1); ll < L; ll = l_atomic[0].fetch_add(1))
            {
                int64 l = ll * structure_loc_coprime % L;
                int64 o = allele_freq_offset[l];
                GENO_READER rt(0, l);
                GENOTYPE* gtab = GetLoc(l).GetGtab();

                for (int i = 0; i < N; ++i)
                {
                    GENOTYPE& gt = gtab[rt.Read()];
                    if (gt.Nalleles() == 0) continue;

                    int* ni = Ni + Z[i] * KT + o;
                    ushort* als = gt.GetAlleleArray();
                    for (int a = 0; a < gt.Ploidy(); ++a)
                        ni[als[a]]++;
                }
            }
        }
        else
        {
            for (int64 ll = l_atomic[0].fetch_add(1); ll < L; ll = l_atomic[0].fetch_add(1))
            {
                int64 l = ll * structure_loc_coprime % L;
                int64 o = allele_freq_offset[l];
                GENO_READER rt(0, l);
                GENOTYPE* gtab = GetLoc(l).GetGtab();

                //avoid thread-conflict
                int64* mi = Mi + N * K * tid2;
                RNG<double> rng(seed + l, RNG_SALT_INITADM);//double

                for (int i = 0; i < N; ++i, mi += K)
                {
                    GENOTYPE& gt = gtab[rt.Read()];
                    if (gt.Nalleles() == 0) continue;

                    int* ni = Ni + o;
                    ushort* als = gt.GetAlleleArray();
                    for (int a = 0; a < gt.Ploidy(); ++a)
                    {
                        ushort k2 = (ushort)rng.Next(K);
                        //avoid thread-conflict
                        mi[k2]++;
                        ni[k2 * KT + als[a]]++;
                    }
                }
            }
        }
    }
}

/* Initialize MCMC for locpriori model */
template<typename REAL>
TARGET void BAYESIAN<REAL>::InitLocPriori()
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
            SetVal(Eta, (double)(1.0 / K), K);

            Gamma = new double[S * K];
            SetVal(Gamma, (double)(1.0 / K), S * K);

            Di = new REAL[S * K];
            SetZero(Di, S * K);
        }
    }
}

/* Initialize MCMC for F model */
template<typename REAL>
TARGET void BAYESIAN<REAL>::InitFmodel()
{
    if (fmodel)
    {
        f = new double[K];
        F = new double[K];
        fbuf1 = new REAL[(maxK + 16) * 32]; 
        fbuf2 = new REAL[(maxK + 16) * 32]; 

        for (int k = 0; k < K; ++k)
        {
            F[k] = pmeanf;
            f[k] = (1 - F[k]) / F[k];
        }
        for (int k = 0; k < K; ++k)
            SetVal(ClusterFreq(k), (REAL)Lambda[k], KT);

        REAL* p = AncestralFreqA;
        SetVal(p, (REAL)lambda, KT);
        for (int64 l = 0; l < L; ++l)
        {
            int k2 = GetLoc(l).k;
            REAL* p2 = total_pop<REAL>->GetFreq(l);
            int nhaplo = total_pop<REAL>->loc_stat1[l].nhaplo;
            for (int a = 0; a < k2; ++a)
                p[a] = (lambda + p2[a] * nhaplo) / (k2 * lambda + nhaplo);
            p += k2;
        }
    }
}

template<typename REAL>
TARGET void BAYESIAN<REAL>::ManageCUDA(bool isrelease)
{
    if (g_gpu_val == 1) return;

    //run complete, release CUDA resources
    if (isrelease)
    {
        if (usecuda)
        {
            UninitCUDA();
            DestroyStreamCUDA();

            for (int devID = 0; devID < nGPU * structure_nsimult; ++devID)
                if (structure_cuda_taskid[devID] == id)
                    structure_cuda_taskid[devID] = 0;

            if (!iscudastream)
            {
                atomic<int>& navail = *(atomic<int>*) & structure_navailable_cuda;
                navail++;
            }
        }
        return;
    }

    if (usecuda || (!iscudastream && structure_navailable_cuda == 0)) return;

    //find a new available cuda device
    for (int devID = 0; devID < nGPU * structure_nsimult; ++devID)
    {
        int zero = 0;
        if (structure_cuda_taskid[devID].compare_exchange_strong(zero, id))
        {
            usecuda = true;

            if (!iscudastream)
            {
                atomic<int>& navail = *(atomic<int>*) & structure_navailable_cuda;
                navail--;
            }
            CreateStreamCUDA(devID / structure_nsimult);
            InitCUDA();
            break;
        }
    }
}

/* Allocate and copy GPU memory */
template<typename REAL>
TARGET void BAYESIAN<REAL>::InitCUDA()
{
    if (g_gpu_val == 1) return;

    int64 tsize = 0; 

    /*Bayes*/    tsize += Align16(sizeof(BAYESIAN<REAL>));
    /*Freq*/     tsize += Align16(K * KT * sizeof(REAL));
    /*Ni*/       tsize += Align16(K * KT * sizeof(int));
    /*Mi*/       tsize += Align16(N * K * sizeof(int64));
    /*Q*/        tsize += Align16(N * K * sizeof(REAL));
    /*Z*/        tsize += Align16(N * sizeof(ushort));
    /*bufNK1*/   tsize += Align16(N * K * sizeof(double));
    /*bufNK2*/   tsize += Align16(N * K * sizeof(double));
    /*bufN1*/    tsize += Align16(N * sizeof(double));
    /*bufN2*/    tsize += Align16(N * sizeof(double));

    byte* p = (byte*)MallocDeviceCUDA(tsize);
    /*Bayes*/    bayes_CUDA = (BAYESIAN<REAL>*)p;       p += Align16(sizeof(BAYESIAN<REAL>));
    /*Freq*/     Freq_CUDA = (REAL*)p;                  p += Align16(K * KT * sizeof(REAL));
    /*Ni*/       Ni_CUDA = (int*)p;                     p += Align16(K * KT * sizeof(int));
    /*Mi*/       Mi_CUDA = (int64*)p;                   p += Align16(N * K * sizeof(int64));
    /*Q*/        Q_CUDA = (REAL*)p;                     p += Align16(N * K * sizeof(REAL));
    /*Z*/        Z_CUDA = (ushort*)p;                   p += Align16(N * sizeof(ushort));
    /*bufNK1*/   bufNK1_CUDA = (double*)p;              p += Align16(N * K * sizeof(double));
    /*bufNK2*/   bufNK2_CUDA = (double*)p;              p += Align16(N * K * sizeof(double));
    /*bufN1*/    bufN1_CUDA = (double*)p;               p += Align16(N * sizeof(double));
    /*bufN2*/    bufN2_CUDA = (double*)p;               p += Align16(N * sizeof(double));

    /*Bayes*/    MemcpyCUDA(bayes_CUDA, this, sizeof(BAYESIAN<REAL>), true);
    /*Freq*/     MemcpyCUDA(Freq_CUDA, Freq, K * KT * sizeof(REAL), true);
    /*Ni*/       MemcpyCUDA(Ni_CUDA, Ni, K * KT * sizeof(int), true);
    /*Mi*/       MemcpyCUDA(Mi_CUDA, Mi, N * K * sizeof(int64), true);
    /*Q*/        MemcpyCUDA(Q_CUDA, Q, N * K * sizeof(REAL), true);
    /*Z*/        MemcpyCUDA(Z_CUDA, Z, N * sizeof(ushort), true);
    /*bufNK1*/   MemsetCUDA(bufNK1_CUDA, 0, N * K * sizeof(double));
    /*bufNK2*/   MemsetCUDA(bufNK2_CUDA, 0, N * K * sizeof(double));
    /*bufN1*/    MemsetCUDA(bufN1_CUDA, 0, N * sizeof(double));
    /*bufN2*/    MemsetCUDA(bufN2_CUDA, 0, N * sizeof(double));
}

/* Free GPU memory */
template<typename REAL>
TARGET void BAYESIAN<REAL>::UninitCUDA()
{
    if (g_gpu_val == 1) return;

    FreeCUDA(bayes_CUDA);
}

/* Update allele frequency for all clusters */
template<typename REAL>
TARGET void BAYESIAN<REAL>::UpdateP(int tid)
{
    if (tid == -1)
    {
        //set lambda, dirichlet distribution priori parameter
        if (fmodel) //checked
        {
            if (fsame) //checked
            {
                Mul(ClusterFreq(0), AncestralFreqA, f[0], KT);
                for (int k = 1; k < K; ++k)
                    SetVal(ClusterFreq(k), ClusterFreq(0), KT);
            }
            else for (int k = 0; k < K; ++k)
                Mul(ClusterFreq(k), AncestralFreqA, f[k], KT);
        }
        else for (int k = 0; k < K; ++k) //checked
            SetVal(ClusterFreq(k), (REAL)Lambda[k], KT);

        //////////////////////////////////////////////////////////

        l_atomic[0] = 0;

        UpdateP(0);

        if (usecuda)
            MemcpyCUDA(Freq_CUDA, Freq, K * KT * sizeof(REAL), true);

        return;
    }

#pragma omp parallel num_threads(structure_nsubthread)
    {
        //split into 128 block, each use a single random number generator and a thread
        //so as to use multi-thread acceleration and keep the same results for differnet #threads
        for (int blockid = l_atomic[0].fetch_add(1); blockid < 128; blockid = l_atomic[0].fetch_add(1))
        {
            int64 kl_start = blockid * K * L / 128;
            int64 kl_end = (blockid + 1) * K * L / 128;
            RNG<double> rng(seed + m * 128 + blockid, RNG_SALT_UPDATEP);//REAL

            for (int64 kl = kl_start; kl < kl_end; ++kl)
            {
                int64 k = kl / L, l = kl % L;
                int* ni = Ni + k * KT + allele_freq_offset[l];
                REAL* p = ClusterLocusFreq(k, l);
                rng.Dirichlet(p, p, ni, (int)GetLoc(l).k);
            }
        }
    }
}

/* Update a priori ancetral proportion for non-admix model */
template<typename REAL>
TARGET void BAYESIAN<REAL>::UpdateQNoAdmixCPU(int tid)
{
    if (tid == -1)
    {
        SetZero(Q, N * K);
        OpenLog((int64*)bufNK1, bufNK2, N * K * structure_nsubthread);

        //add priori probability
        double* buf1 = bufNK1, * buf2 = bufNK2;
        if (locpriori) for (int i = 0; i < N; ++i, buf1 += K, buf2 += K)
        {
            if (ainds<REAL>[i]->vt == 0) continue;
            ChargeLog((int64*)buf1, buf2, Gamma + ainds<REAL>[i]->popid * K, K);
        }

        //////////////////////////////////////////////////////////

        l_atomic[0] = 0;

        UpdateQNoAdmixCPU(0);

        //avoid thread-conflict
        for (int i = 1; i < structure_nsubthread; ++i)
            Add(bufNK1, bufNK1 + N * K * i, N * K);

        buf1 = bufNK1;
        REAL* q = Q;
        RNG<double> rng(seed + m, RNG_SALT_UPDATEQ);//double
        for (int i = 0; i < N; ++i, buf1 += K, q += K)
        {
            if (ainds<REAL>[i]->vt == 0) continue;
            ushort k2 = (ushort)rng.PolyLog(buf1, K);
            q[k2] = 1;
            Z[i] = k2;
        }

        binaryq = true;
        return;
    }

    atomic<int> thread_counter = 0;
#pragma omp parallel num_threads(structure_nsubthread)
    {
        int tid2 = thread_counter.fetch_add(1);

        for (int64 ll = l_atomic[0].fetch_add(1); ll < L; ll = l_atomic[0].fetch_add(1))
        {
            int64 l = ll * structure_loc_coprime % L;
            REAL* p = Freq + allele_freq_offset[l];
            GENO_READER rt(0, l);
            GENOTYPE* gtab = GetLoc(l).GetGtab();

            //avoid thread-conflict
            double* buf1 = bufNK1 + N * K * tid2, * buf2 = bufNK2 + N * K * tid2;
            for (int i = 0; i < N; ++i, buf1 += K, buf2 += K)
            {
                GENOTYPE& gt = gtab[rt.Read()];
                if (gt.Nalleles() == 0) continue;
                ushort* als = gt.GetAlleleArray();
                for (int a = 0; a < gt.Ploidy(); ++a)
                    ChargeLog((int64*)buf1, buf2, p + als[a], K, KT);
            }
        }

        //avoid thread-conflict
        CloseLog((int64*)bufNK1 + N * K * tid2, bufNK2 + N * K * tid2, N * K);
    }
}

/* Update a priori ancetral proportion for non-admix model */
template<typename REAL>
TARGET void BAYESIAN<REAL>::UpdateQNoAdmix()
{
    binaryq = true;

    if (structure_eval_val == 2)
    {
        if (usecuda)
            UpdateQNoAdmixCUDA();
        else switch (SIMD_TYPE)
        {
#ifdef __aarch64__
        case 2: UpdateQNoAdmixNEO(); return;
#else
        case 4: UpdateQNoAdmix512(); return;
        case 3: UpdateQNoAdmixAVX(); return;
        case 2: UpdateQNoAdmixSSE(); return;
#endif
        default: UpdateQNoAdmixCPU(); return;
        }

        return;
    }

    //evaluate
#ifdef __aarch64__
    for (int i = 2; i >= 1; --i)
#else
    for (int i = 5; i >= 1; --i)
#endif
    {
        timepoint begin = GetNow();
        switch (i)
        {
#ifdef __aarch64__
        case 2: UpdateQNoAdmixNEO(); break;
#else
        case 5: UpdateQNoAdmixCUDA(); break;
        case 4: UpdateQNoAdmix512(); break;
        case 3: UpdateQNoAdmixAVX(); break;
        case 2: UpdateQNoAdmixSSE(); break;
#endif
        case 1: UpdateQNoAdmixCPU(); break;
        }
        printf("QNoAdmix%s %0.5f %x %0.15e\n", simdstr[i], GetElapse(begin), HashString((char*)Z, N * 2), LogProd(bufNK1, N * K));
    }
}

/* Update a priori ancetral proportion for admix model */
template<typename REAL>
TARGET void BAYESIAN<REAL>::UpdateQAdmix()
{
    //used by noadmix (in admburnin) and admix model
    binaryq = false;
    RNG<double> rng(seed + m, RNG_SALT_UPDATEQ);//REAL

    if (locpriori)
    {
        for (int i = 0; i < N; ++i)
            if (ainds<REAL>[i]->vt)
                rng.Dirichlet(Q + i * K, AlphaLocal + ainds<REAL>[i]->popid * K, Mi + i * K, K);
    }
    else
    {
        for (int i = 0; i < N; ++i)
            if (ainds<REAL>[i]->vt)
                rng.Dirichlet(Q + i * K, Alpha, Mi + i * K, K);
    }
}

template<typename REAL>
template<bool fast_fp32>
TARGET void BAYESIAN<REAL>::UpdateQMetroCPU(int tid)
{
    if (tid == -1)
    {
        RNG<double> rng(seed + m, RNG_SALT_UPDATEQ);//REAL
        REAL* bufi = (REAL*)bufNK1;
        REAL* q = NULL;

        for (int i = 0; i < N; ++i, bufi += K)
        {
            if (ainds<REAL>[i]->vt == 0) continue;
            if (locpriori) rng.Dirichlet(bufi, AlphaLocal + ainds<REAL>[i]->popid * K, K);
            else           rng.Dirichlet(bufi, Alpha, K);
        }

        OpenLog((int64*)bufN1, bufN2, N * structure_nsubthread);

        //////////////////////////////////////////////////////////

        l_atomic[0] = 0;

        g_fastsingle_val == 1 ? UpdateQMetroCPU<true >(0) : UpdateQMetroCPU<false>(0);

        //avoid thread-conflict
        for (int i = 1; i < structure_nsubthread; ++i)
            Add(bufN1, bufN1 + N * i, N);

        bufi = (REAL*)bufNK1; q = Q;
        for (int i = 0; i < N; ++i, q += K, bufi += K)
        {
            if (ainds<REAL>[i]->vt == 0) continue;
            if (bufN1[i] >= NZERO || rng.Uniform() < exp(bufN1[i]))
                SetVal(q, bufi, K);
        }

        return;
    }

    atomic<int> thread_counter = 0;
#pragma omp parallel num_threads(structure_nsubthread)
    {
        int tid2 = thread_counter.fetch_add(1);

#ifdef PARFOR_PARTITION
        for (int64 l = tid2 * L / structure_nsubthread; l < (tid2 + 1) * L / structure_nsubthread; ++l)
        {
#else
        for (int64 ll = l_atomic[0].fetch_add(1); ll < L; ll = l_atomic[0].fetch_add(1))
        {
            int64 l = ll * structure_loc_coprime % L;
#endif
            GENO_READER rt(0, l);
            GENOTYPE* gtab = GetLoc(l).GetGtab();
            REAL* p = Freq + allele_freq_offset[l];
            REAL* bufi = (REAL*)bufNK1;
            REAL* q = Q;

            double* buf1 = bufN1 + N * tid2, * buf2 = bufN2 + N * tid2;
            for (int i = 0; i < N; ++i, q += K, bufi += K, buf1++, buf2++)
            {
                GENOTYPE& gt = gtab[rt.Read()];
                if (gt.Nalleles() == 0)
                    continue;

                ushort* als = gt.GetAlleleArray();
                for (int a = 0; a < gt.Ploidy(); ++a)
                {
                    if constexpr (std::is_same_v<REAL, double> || !fast_fp32)
                        ChargeLog(*(int64*)buf1, *buf2, SumProdDiv(bufi, q, p + als[a], KT, K));
                    else
                        ChargeLog(*(int64*)buf1, *buf2, SumProdDivx(bufi, q, p + als[a], KT, K));
                }
            }
        }

        //avoid thread-conflict
        CloseLog((int64*)bufN1 + N * tid2, bufN2 + N * tid2, N);
    }
}

/* Update a priori ancetral proportion by Metropolis-Hastings for admix model */
template<typename REAL>
TARGET void BAYESIAN<REAL>::UpdateQMetro()
{
    binaryq = false;

    if (structure_eval_val == 2)
    {
        if (usecuda)
            UpdateQMetroCUDA();
        else switch (SIMD_TYPE)
        {
#ifdef __aarch64__
        case 2: UpdateQMetroNEO(); return;
#else
        case 4: UpdateQMetro512(); return;
        case 3: UpdateQMetroAVX(); return;
        case 2: UpdateQMetroSSE(); return;
#endif
        default: UpdateQMetroCPU(); return;
        }

        return;
    }

    //evaluate
#ifdef __aarch64__
    for (int i = 2; i >= 1; --i)
#else
    for (int i = 5; i >= 1; --i)
#endif
    {
        SetVal(Q, (REAL)(1.0 / K), N * K);
        timepoint begin = GetNow();
        switch (i)
        {
#ifdef __aarch64__
        case 2: UpdateQMetroNEO(); break;
#else
        case 5: UpdateQMetroCUDA(); break;
        case 4: UpdateQMetro512(); break;
        case 3: UpdateQMetroAVX(); break;
        case 2: UpdateQMetroSSE(); break;
#endif
        case 1: UpdateQMetroCPU(); break;
        }
        printf("QMetro%s %0.5f %0.15e\n", simdstr[i], GetElapse(begin), Sum(bufN1, N));
    }
}

/* Update a priori ancetral proportion */
template<typename REAL>
TARGET void BAYESIAN<REAL>::UpdateQ()
{
    //draw gene proportion for each individual
    if (!admix)
    {
        if (m >= nadmburnin) //checked
            UpdateQNoAdmix();
        else
            UpdateQAdmix();

        if (locpriori)
        {
            SetZero(Di, S * K);
            REAL* q = Q;
            for (int i = 0; i < N; ++i, q += K)
                if (ainds<REAL>[i]->vt) //bug fixed for individuals with no data
                    Add(Di + ainds<REAL>[i]->popid * K, q, K);
        }
    }
    else
    {
        if (metrofreq > 0 && m % metrofreq == 0) //checked
            UpdateQMetro();
        else
            UpdateQAdmix();
    }
}

/* Update locpriori parameters */
template<typename REAL>
TARGET void BAYESIAN<REAL>::UpdateLocPriori()
{
    if (!locpriori) return;

    RNG<double> rng(seed + m, RNG_SALT_UPDATELoc);//double
    if (admix)
    {
        double rm = rng.Uniform(r - epsr, r + epsr);
        if (rm > 0 && rm < maxr)
        {
            double dlnL = 0, d = rm - r, dt1 = rm * MyLog(rm) - r * MyLog(r);
            for (int k = 0; k < K; ++k)
            {
                dlnL += S * (Alpha[k] * dt1 - LogGamma(rm * Alpha[k]) + LogGamma(r * Alpha[k]));
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
            double dlnL = LogGamma(rm) - LogGamma(r);
            for (int k = 0; k < K; ++k)
                dlnL += LogGamma(r * Eta[k]) - LogGamma(rm * Eta[k]);
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
            double delta = rng.Uniform(0, epseta);
            double e1 = Eta[i1] + delta;
            double e2 = Eta[j1] - delta;
            if (e1 < 1 && e2 > 0)
            {
                double dlnL = S * (LogGamma(r * Eta[i1]) + LogGamma(r * Eta[j1]) - LogGamma(r * e1) - LogGamma(r * e2));
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
                double delta = rng.Uniform(0, epsgamma);

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
template<typename REAL>
TARGET void BAYESIAN<REAL>::UpdateZNoAdmixCPU(int tid)
{
    if (tid == -1)
    {
        SetZero(Mi, N * K);
        SetZero(Ni, K * KT);

        // Z in already updated in update Q
        //Count number of allele copies in individual i and cluster k
        for (int i = 0; i < N; ++i)
            Mi[i * K + Z[i]] = ainds<REAL>[i]->vt;

        //////////////////////////////////////////////////////////

        l_atomic[0] = 0;

        UpdateZNoAdmixCPU(0);

        return;
    }

#pragma omp parallel num_threads(structure_nsubthread)
    {
        for (int64 ll = l_atomic[0].fetch_add(1); ll < L; ll = l_atomic[0].fetch_add(1))
        {
            int64 l = ll * structure_loc_coprime % L;
            int64 o = allele_freq_offset[l];
            GENO_READER rt(0, l);
            GENOTYPE* gtab = GetLoc(l).GetGtab();

            for (int i = 0; i < N; ++i)
            {
                GENOTYPE& gt = gtab[rt.Read()];
                if (gt.Nalleles() == 0) continue;

                ushort* als = gt.GetAlleleArray();
                int* ni = Ni + Z[i] * KT + o;
                for (int a = 0; a < gt.Ploidy(); ++a)
                    ni[als[a]]++;
            }
        }
    }
}

/* Update ancestral proportion for each allele or for each individual */
template<typename REAL>
TARGET void BAYESIAN<REAL>::UpdateZNoAdmix()
{
    if (structure_eval_val == 2)
    {
        if (usecuda)
            UpdateZNoAdmixCUDA();
        else switch (SIMD_TYPE)
        {
#ifdef __aarch64__
        case 2: UpdateZNoAdmixNEO(); return;
#else
        case 4: UpdateZNoAdmix512(); return;
        case 3: UpdateZNoAdmixAVX(); return;
        case 2: UpdateZNoAdmixSSE(); return;
#endif
        default: UpdateZNoAdmixCPU(); return;
        }

        return;
    }

    //evaluate
#ifdef __aarch64__
    for (int i = 2; i >= 1; --i)
#else
    for (int i = 5; i >= 1; --i)
#endif
    {
        timepoint begin = GetNow();
        switch (i)
        {
#ifdef __aarch64__
        case 2: UpdateZNoAdmixNEO(); break;
#else
        case 5: UpdateZNoAdmixCUDA(); break;
        case 4: UpdateZNoAdmix512(); break;
        case 3: UpdateZNoAdmixAVX(); break;
        case 2: UpdateZNoAdmixSSE(); break;
#endif
        case 1: UpdateZNoAdmixCPU(); break;
        }
        printf("UpdateZNoAdmix%s %0.5f %x\n", simdstr[i], GetElapse(begin), HashString((char*)Ni, K * KT * 4));
    }
}

/* Update ancestral proportion for each allele or for each individual */
template<typename REAL>
template<bool fast_fp32>
TARGET void BAYESIAN<REAL>::UpdateZAdmixCPU(int tid)
{
    if (tid == -1)
    {
        SetZero(Mi, N * K * structure_nsubthread);
        SetZero(Ni, K * KT);

        //////////////////////////////////////////////////////////

        l_atomic[0] = 0;

        g_fastsingle_val == 1 ? UpdateZAdmixCPU<true >(0) : UpdateZAdmixCPU<false>(0);

        //avoid thread-conflict
        for (int i = 1; i < structure_nsubthread; ++i)
            Add(Mi, Mi + N * K * i, N * K);

        return;
    }

    atomic<int> thread_counter = 0;
#pragma omp parallel num_threads(structure_nsubthread)
    {
        int tid2 = thread_counter.fetch_add(1);

        for (int64 ll = l_atomic[0].fetch_add(1); ll < L; ll = l_atomic[0].fetch_add(1))
        {
            int64 l = ll * structure_loc_coprime % L;
            GENO_READER rt(0, l);
            GENOTYPE* gtab = GetLoc(l).GetGtab();
            int64 o = allele_freq_offset[l];
            RNG<double> rngd; RNG<float > rngs;

            if constexpr (std::is_same_v<REAL, double> || !fast_fp32)
                new(&rngd) RNG<double>(seed + m * L + l, RNG_SALT_UPDATEZ);
            else
                new(&rngs) RNG<float >(seed + m * L + l, RNG_SALT_UPDATEZ);

            //avoid thread-conflict
            double* bufkd = bufKthread + tid2 * K;
            float* bufks = (float*)bufkd;
            int64* mi = Mi + N * K * tid2;
            REAL* q = Q;

            for (int i = 0; i < N; ++i, mi += K, q += K)
            {
                GENOTYPE& gt = gtab[rt.Read()];
                int ploidy = gt.Nalleles() ? gt.Ploidy() : 0;
                ushort* als = gt.GetAlleleArray();
                int* ni = Ni + o;

                for (int a = 0; a < maxploidy; ++a)
                {
                    //alway draw to keey results consistent
                    if (a >= ploidy)
                    {
                        if constexpr (std::is_same_v<REAL, double> || !fast_fp32)
                            rngd.XorShift();
                        else
                            rngs.XorShift();
                        continue;
                    }

                    //The same to previous allele, do not update multinomial prob to save time
                    if (a == 0 || als[a] != als[a - 1])
                        for (int k = 0; k < K; ++k)
                            if constexpr (std::is_same_v<REAL, double> || !fast_fp32)
                                bufkd[k] = (double)q[k] * (double)ClusterAlleleFreq(k, l, als[a]);
                            else
                                bufks[k] = q[k] * ClusterAlleleFreq(k, l, als[a]);

                    ushort k2;
                    //draw cluster for each allele copy
                    if constexpr (std::is_same_v<REAL, double> || !fast_fp32)
                        k2 = (ushort)rngd.Poly(bufkd, K);
                    else
                        k2 = (ushort)rngs.Poly(bufks, K);

                    //Update Mi, NI
                    mi[k2]++;
                    ni[k2 * KT + als[a]]++;
                }
            }
        }
    }
}

/* Update ancestral proportion for each allele or for each individual */
template<typename REAL>
TARGET void BAYESIAN<REAL>::UpdateZAdmix()
{
    if (structure_eval_val == 2)
    {
        if (usecuda)
            UpdateZAdmixCUDA();
        else switch (SIMD_TYPE)
        {
#ifdef __aarch64__
        case 2: UpdateZAdmixNEO(); return;
#else
        case 4: UpdateZAdmix512(); return;
        case 3: UpdateZAdmixAVX(); return;
        case 2: UpdateZAdmixSSE(); return;
#endif
        default: UpdateZAdmixCPU(); return;
        }

        return;
    }

    //evaluate
#ifdef __aarch64__
    for (int i = 2; i >= 1; --i)
#else
    for (int i = 5; i >= 1; --i)
#endif
    {
        SetZero(bufNK1, N * K);
        timepoint begin = GetNow();
        switch (i)
        {
#ifdef __aarch64__
        case 2: UpdateZAdmixNEO(); break;
#else
        case 5: UpdateZAdmixCUDA(); break;
        case 4: UpdateZAdmix512(); break;
        case 3: UpdateZAdmixAVX(); break;
        case 2: UpdateZAdmixSSE(); break;
#endif
        case 1: UpdateZAdmixCPU(); break;
        }
        printf("ZAdmix%s %0.5f %x %x %0.15e\n", simdstr[i], GetElapse(begin), HashString((char*)Ni, K * KT * 4), HashString((char*)Mi, N * K * 8), Sum(bufNK1, structure_nsubthread));

    }
}

/* Update ancestral proportion for each allele or for each individual */
template<typename REAL>
TARGET void BAYESIAN<REAL>::UpdateZ()
{
    if (binaryq)//checked
        UpdateZNoAdmix();
    else
        UpdateZAdmix();
}

/* Update Dirichlet parameter alpha (to draw admixture proportion Q) in the admix model */
template<typename REAL>
TARGET void BAYESIAN<REAL>::UpdateAlpha()
{
    RNG<double> rng(seed + m, RNG_SALT_UPDATEAlpha);//double

    if (locpriori && admix)//checked
    {
        //UpdateAlphaLocPrior
        
        //update AlphaGlobal
        for (int k = 0; k < K; ++k)
        {
            double am = rng.Normal(Alpha[k], stdalpha);
            if (am > 0 && am < maxalpha)
            {
                double d = am - Alpha[k];
                double dlnL = S * (d * r * MyLog(r) + LogGamma(r * Alpha[k]) - LogGamma(r * am));

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
            int snind = apops<REAL>[s]->nind;
            REAL* sq = &Q[apops<REAL>[s]->ind0id * K + 0];

            for (int k = 0; k < K; ++k)
            {
                double ao = AlphaLocal[s * K + k], am = rng.Normal(ao, stdalpha);

                if (am > 0 && am < maxalpha)
                {
                    double d = am - ao;
                    double dlnL = (r * Alpha[k] - 1) * MyLog(am / ao) - d * r +
                        apops<REAL>[s]->nind * (LogGamma(SumAlpha[s] + d) -
                            LogGamma(am) - LogGamma(SumAlpha[s]) + LogGamma(ao));

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
    else if (inferalpha)//checked
    {
        if (!admix && m >= nadmburnin) return;//checked

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

                dlnL += (am - ao) * LogProd(Q + k, N, K) + (LogGamma(sumalpha + d) - LogGamma(sumalpha) - LogGamma(am) + LogGamma(ao)) * N;
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

            dlnL += (am - ao) * LogProd(Q, N * K) + (LogGamma(K * am) - LogGamma(K * ao) - K * LogGamma(am) + K * LogGamma(ao)) * N;
            if (dlnL >= NZERO || rng.Uniform() < exp(dlnL))
                for (int k = 0; k < K; ++k)
                    Alpha[k] = am;
        }
    }
}

/* Update Dirichlet parameter lambda (to draw allele frequency) */
template<typename REAL>
TARGET void BAYESIAN<REAL>::UpdateLambda()
{
    if (!inferlambda) return;

    RNG<double> rng(seed + m, RNG_SALT_UPDATELambda);//checked

    if (difflambda) for (int k = 0; k < K; ++k)
    {
        double lo = Lambda[k], lm = rng.Normal(lo, stdlambda);
        if (lm <= 0 || lm >= maxlambda) continue;

        double dlnL = 0, gd = LogGamma(lo) - LogGamma(lm);

        for (int i = 2; i <= maxK; ++i)
            dlnL += kdis[i] * (LogGamma(i * lm) - LogGamma(i * lo) + i * gd);

        dlnL += (lm - lo) * LogProd(ClusterFreq(k), KT);

        if (dlnL >= NZERO || rng.Uniform() < exp(dlnL))
            Lambda[k] = lm;
    }
    else
    {
        int np = fmodel ? 1 : K;
        double lo = Lambda[0], lm = rng.Normal(lo, stdlambda);
        if (lm <= 0 || lm >= maxlambda) return;

        double dlnL = 0, gd = LogGamma(lo) - LogGamma(lm);

        for (int i = 2; i <= maxK; ++i)
            dlnL += kdis[i] * np * (LogGamma(i * lm) - LogGamma(i * lo) + i * gd);

        dlnL += (lm - lo) * (fmodel ? LogProd(AncestralFreqA, KT) : LogProd(Freq, K * KT));

        if (dlnL >= NZERO || rng.Uniform() < exp(dlnL))
            SetVal(Lambda, lm, K);
    }
}

/* Update allele frequency of ancestral population for the F model */
template<typename REAL>
TARGET void BAYESIAN<REAL>::UpdateA()
{
    if (!fmodel) return;

    if (HashULong(seed + m) & 1)
        UpdateA1();
    else
        UpdateA2();
}

template<typename REAL>
TARGET void BAYESIAN<REAL>::UpdateA1(int tid)
{
    if (tid == -1)
    {
        l_atomic[0] = 0;

        UpdateA1(0);

        return;
    }

#pragma omp parallel num_threads(structure_nsubthread)
    {
        double lam = Lambda[0];//single lambda in F model?

        //split into 32 block, each use a single random number generator and a thread
        //so as to use multi-thread acceleration and keep the same results for differnet #threads
        for (int blockid = l_atomic[0].fetch_add(1); blockid < 32; blockid = l_atomic[0].fetch_add(1))
        {
            int64 l_start = blockid * L / 32;
            int64 l_end = (blockid + 1) * L / 32;
            RNG<double> rng(seed + m * 32 + blockid, RNG_SALT_UPDATEA);//REAL

            for (int64 l = l_start; l < l_end; ++l)
            {
                int k2 = GetLoc(l).k;
                REAL* pa = AncestralLocusFreqA(l);
                REAL* pb = fbuf1 + (maxK + 16) * blockid;
                REAL* pc = fbuf2 + (maxK + 16) * blockid;

                SetVal(pb, (REAL)lam, k2);
                for (int k = 0; k < K; ++k)
                    AddProd(pb, ClusterLocusFreq(k, l), f[k], k2);

                rng.Dirichlet(pc, pb, k2);

                double dlnL = 0;
                for (int a = 0; a < k2; ++a)
                    dlnL += (pb[a] - lam) * (MyLog(pa[a] / pc[a]));

                for (int k = 0; k < K; ++k)
                {
                    REAL* pk = ClusterLocusFreq(k, l);
                    double fr = f[k];
                    for (int a = 0; a < k2; ++a)
                        dlnL += LogGamma(fr * pa[a]) - LogGamma(fr * pc[a])
                        + fr * (pc[a] - pa[a]) * MyLog(pk[a]);
                }

                if (dlnL >= NZERO || rng.Uniform() < exp(dlnL))
                    SetVal(pa, pc, k2);
            }
        }

        /* //Original version
        {
            REAL** bufK = (REAL**)bufNK1;
            //Independent update ancestral allele frequency
            SetVal(AncestralFreqB, (REAL)lam, KT);

            for (int k = 0; k < K; ++k)
                AddProd(AncestralFreqB, ClusterFreq(k), f[k], KT);

            REAL* pa = AncestralFreqA, * pb = AncestralFreqB, * pc = AncestralFreqC;

            for (int k = 0; k < K; ++k)
                bufK[k] = ClusterFreq(k);

            for (int64 l = 0; l < L; ++l)
            {
                int k2 = GetLoc(l).k;
                rng.Dirichlet(pc, pb, k2);

                double dlnL = 0;
                for (int a = 0; a < k2; ++a)
                    dlnL += (pb[a] - lam) * (MyLog(pa[a] / pc[a]));

                for (int k = 0; k < K; ++k)
                {
                    double fr = f[k];
                    for (int a = 0; a < k2; ++a)
                        dlnL += LogGamma(fr * pa[a]) - LogGamma(fr * pc[a])
                        + fr * (pc[a] - pa[a]) * MyLog(*bufK[k]++);
                }

                if (dlnL >= NZERO || rng.Uniform() < exp(dlnL))
                    SetVal(pa, pc, k2);

                pa += k2; pb += k2; pc += k2;
            }
        }
        */
    }
}

template<typename REAL>
TARGET void BAYESIAN<REAL>::UpdateA2(int tid)
{
    if (tid == -1)
    {
        l_atomic[0] = 0;

        UpdateA2(0);

        return;
    }

#pragma omp parallel num_threads(structure_nsubthread)
    {
        double lam = Lambda[0];//single lambda in F model?

        for (int blockid = l_atomic[0].fetch_add(1); blockid < 32; blockid = l_atomic[0].fetch_add(1))
        {
            int64 l_start = blockid * L / 32;
            int64 l_end = (blockid + 1) * L / 32;
            RNG<double> rng(seed + m * 32 + blockid, RNG_SALT_UPDATEA);//double

            uint64 loffset = allele_freq_offset[l_start];
            REAL* pa = AncestralFreqA + loffset, * pb = AncestralFreqB + loffset, * pc = AncestralFreqC + loffset;
            REAL* p = Freq + loffset;

            for (int64 l = l_start; l < l_end; ++l)
            {
                int k2 = GetLoc(l).k;
                int m1 = rng.Next(k2), n1 = rng.NextAvoid(k2, m1);

                double delta = rng.Uniform(0, pow(N, -0.5));
                double pm0 = pa[m1], pm1 = pm0 + delta, pn0 = pa[n1], pn1 = pn0 - delta;
                if (pm1 >= 1 || pn1 <= 0)
                {
                    pa += k2; pb += k2; pc += k2;
                    continue;
                }

                double dlnL = (lam - 1) * MyLog(pm1 * pn1 / pm0 / pn0);

                for (int k = 0; k < K; ++k)
                {
                    double fr = f[k];
                    dlnL += LogGamma(fr * pm0) + LogGamma(fr * pn0)
                          - LogGamma(fr * pm1) - LogGamma(fr * pn1)
                          + (fr * delta) * MyLog(p[k * KT + m1] / p[k * KT + n1]);
                }

                if (dlnL >= NZERO || rng.Uniform() < exp(dlnL))
                {
                    pa[m1] = pm1;
                    pa[n1] = pn1;
                }

                p += k2; pa += k2; pb += k2; pc += k2;
            }
        }
    }
}

/* /Update cluster-specific Fst for the F model */
template<typename REAL>
TARGET void BAYESIAN<REAL>::UpdateF()
{
    //F model
    if (!fmodel) return;

    RNG<double> rng(seed + m, RNG_SALT_UPDATEF);//checked

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
            (fsame ? K : 1) * nloc * (LogGamma(fm) - LogGamma(fo));//eL add in 2018.09

        for (int k2 = k; k2 < K; ++k2)
        {
            REAL* p = ClusterFreq(k2);
            REAL* pa = AncestralFreqA;//2018.09
            for (int a = 0; a < KT; ++a)
            {
                REAL pa2 = *pa++;
                dlnL += pa2 * (fm - fo) * MyLog(*p++) + LogGamma(fo * pa2) - LogGamma(fm * pa2);
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
template<typename REAL>
TARGET void BAYESIAN<REAL>::Arrange()
{
    UnifyInt64ToDouble(MiSum, N, K);//cast to double
    Mul(rout, 1.0 / nr, rlen);
    rout[2] = rout[1] - rout[0] * rout[0]; 
    rout[3] = rout[0] - rout[2] / 2;
    Mul(FreqM, (REAL)(1.0 / nr), K * KT);
}

/* Record updated MCMC parameters */
template<typename REAL>
template<bool isadmix, bool fast_fp32>
TARGET void BAYESIAN<REAL>::RecordCPU(int tid)
{
    if (tid == -1)
    {
        Add(MiSum, Mi, N * K);

        //////////////////////////////////////////////////////////

        l_atomic[0] = 0;

        switch ((int)binaryq * 10 + g_fastsingle_val)
        {
        case 01: BAYESIAN<REAL>::RecordCPU<true , true >(0); break;
        case 02: BAYESIAN<REAL>::RecordCPU<true , false>(0); break;
        case 11: BAYESIAN<REAL>::RecordCPU<false, true >(0); break;
        case 12: BAYESIAN<REAL>::RecordCPU<false, false>(0); break;
        }

        bufNK1[0] = Sum(bufNK1, structure_nsubthread);
        return;
    }

    atomic<int> thread_counter = 0;
#pragma omp parallel num_threads(structure_nsubthread)
    {
        int tid2 = thread_counter.fetch_add(1);

        int64 slog = 0; double prod = 1;
        OpenLog(slog, prod);

        for (int64 ll = l_atomic[0].fetch_add(1); ll < L; ll = l_atomic[0].fetch_add(1))
        {
            int64 l = ll * structure_loc_coprime % L;
            GENO_READER rt(0, l);
            GENOTYPE* gtab = GetLoc(l).GetGtab();
            REAL* p = Freq + allele_freq_offset[l], * q = Q;

            for (int i = 0; i < N; ++i, q += K)
            {
                GENOTYPE& gt = gtab[rt.Read()];
                if (gt.Nalleles() == 0) continue;

                REAL* pp;
                if constexpr (!isadmix)
                    pp = p + KT * Z[i];

                ushort* als = gt.GetAlleleArray();
                for (int a = 0; a < gt.Ploidy(); ++a)
                {
                    if constexpr (isadmix)
                    {
                        if constexpr (std::is_same_v<REAL, double> || !fast_fp32)
                            ChargeLog(slog, prod, SumProd(q, p + als[a], KT, K));
                        else
                            ChargeLog(slog, prod, SumProdx(q, p + als[a], KT, K));
                    }
                    else
                        ChargeLog(slog, prod, pp[als[a]]);
                }
            }
        }

        CloseLog(slog, prod);
        bufNK1[tid2] = prod;
    }
}

/* Record updated MCMC parameters */
template<typename REAL>
TARGET void BAYESIAN<REAL>::Record()
{
    if ((structure_writelnl_val == 1 && m < nburnin && m % nthinning == 0) || (m >= nburnin && (m - nburnin) % nthinning == 0))
    {
        if (structure_eval_val == 2)
        {
            if (usecuda)
                RecordCUDA();
            else switch (SIMD_TYPE)
            {
#ifdef __aarch64__
            case 2: RecordNEO(); break;
#else
            case 4: Record512(); break;
            case 3: RecordAVX(); break;
            case 2: RecordSSE(); break;
#endif
            default: RecordCPU(); break;
            }
        }
        else
        {
            //evaluate
#ifdef __aarch64__
            for (int i = 2; i >= 1; --i)
#else
            for (int i = 5; i >= 1; --i)
#endif
            {
                SetZero(bufNK1, N * K);
                timepoint begin = GetNow();
                switch (i)
                {
#ifdef __aarch64__
                case 2: RecordNEO(); break;
#else
                case 5: RecordCUDA(); break;
                case 4: Record512(); break;
                case 3: RecordAVX(); break;
                case 2: RecordSSE(); break;
#endif
                case 1: RecordCPU(); break;
                }
                printf("Record%s %0.5f %0.15e\n", simdstr[i], GetElapse(begin), bufNK1[0]);
            }
        }

        /*
        Add(MiSum, Mi, N * K);
        int64 slog = 0; double prod = 1;
        REAL* p = Freq;
        OpenLog(slog, prod);
        for (int64 l = 0; l < L; p += GetLoc(l).k, ++l)
        {
            GENO_READER rt(0, l);
            GENOTYPE* gtab = GetLoc(l).GetGtab();
            REAL* q = Q;

            for (int i = 0; i < N; ++i, q += K)
            {
                GENOTYPE& gt = gtab[rt.Read()];
                if (gt.Nalleles() == 0) continue;

                ushort* als = gt.GetAlleleArray();
                for (int a = 0; a < gt.Ploidy(); ++a)
                    ChargeLog(slog, prod, SumProd(q, p + als[a], KT, K));
            }
        }
        CloseLog(slog, prod);
        */

        if (structure_writelnl_val == 1) 
            lnL[nlnL++] = bufNK1[0];

        if (m >= nburnin)
        {
            nr++;
            //add to result
            double likelihood = bufNK1[0];
            rout[0] += likelihood;
            rout[1] += likelihood * likelihood;

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
            Add(FreqM, Freq, K * KT);
        }
    }

    PROGRESS_VALUE += K;
}

/* Free memory */
template<typename REAL>
TARGET void BAYESIAN<REAL>::Uninit()
{
    //Normal
    DEL(kdis);
    DEL(bufKthread);
    DEL(MiSum);
    DEL(Alpha);
    DEL(Lambda);
    DEL(rout);
    FREE(FreqA);
    FREE(FreqM);

    FREECUDA(Freq);
    FREECUDA(Ni);
    FREECUDA(Mi);
    FREECUDA(Q);
    FREECUDA(bufNK1o);          bufNK1 = NULL; 
    FREECUDA(bufNK2o);          bufNK2 = NULL; 
    FREECUDA(bufN1);
    FREECUDA(bufN2);

    //Write lnL
    DEL(lnL);

    //Adm
    FREECUDA(Z);

    //Fmodel
    if (fmodel)
    {
        DEL(F);
        DEL(f);
        DEL(fbuf1);
        DEL(fbuf2);
    }

    //Loc
    if (locpriori)
    {
        DEL(SumAlpha);
        DEL(AlphaLocal);

        if (!admix)
        {
            DEL(Eta);
            DEL(Gamma);
            DEL(Di);
        }
    }
}

template<typename T>
TARGET int IdxCmp(T* a, T* b, int size)
{
    if (!memcmp(a, b, size * sizeof(T))) return -1;
    for (int i = 0; i < size; ++i)
        if (a[i] != b[i])
            return i;
    return -1;
}

/* Perform MCMC */
template<typename REAL>
TARGET void BAYESIAN<REAL>::MCMC()
{
    double singleiter = 0;
    timepoint begin1 = GetNow();

    InitFix();
    InitAdmix();
    InitLocPriori();
    InitFmodel();

    for (m = 0; m < nburnin + nreps; ++m)
    {
        timepoint begin2 = GetNow();

        //Manage CUDA devices
        ManageCUDA(false);

        //Update allele frequency for all clusters
        UpdateP();

        //Update a priori ancetral proportion
        UpdateQ();
        
        //Update locpriori parameters r
        UpdateLocPriori();

        //Update ancestral proportion for each allele or for each individual
        UpdateZ();
        
        //Update Dirichlet parameter alpha (to draw admixture proportion Q) in the admix model
        UpdateAlpha();

        //Update Dirichlet parameter lambda (to draw allele frequency)
        UpdateLambda();

        //Update allele frequency of ancestral population for the F model
        UpdateA();

        //Update cluster-specific Fst for the F model
        UpdateF();

        //Record likelihood and Q, Z
        Record();

        if (structure_eval_val == 1) 
            printf("StructureSingleIteration %0.5f s\n", singleiter = GetElapse(begin2));
    }

    Arrange();

    ManageCUDA(true);

    if (structure_eval_val == 1)
        printf("StructureTaskOverhead %0.5f s\n", GetElapse(begin1) - singleiter);
}

#endif

#pragma pack(pop)

#define extern
extern STRUCTURE_RUNINFO* structure_par;                        //Genetic distance for pcoa and hierarchical clustering
extern int structure_totalruns;                                 //Total number of runs in Bayesian clustering
extern int structure_nsubthread;								//Total number of runs in Bayesian clustering
extern int structure_loc_size_min, structure_loc_size_max;
extern int64* structure_loc_lend;
extern int64 structure_loc_coprime;
extern int64 structure_loc_coprime64[32];
extern int64* structure_loc_original_idx;
extern int structure_navailable_cuda;                           //number of free cuda devices
extern atomic<int>* structure_cuda_taskid;                      //current taskids using cuda
extern int structure_nsimult;                                   //number of task simultaneous run on a single cuda device
#undef extern

/* Is two numbers are coprime*/
bool IsCoprime(int64 x, int64 y)
{
    int64 z = 1;
    while (x % y != 0) 
    {
        z = x % y;
        x = y;
        y = z;
    }
    if (z == 1)
        return true;
    return false;
}

/* Quick sort locus by genotype bitsize */
TARGET void QSLocus2(int64 left, int64 right)
{
    int64 i = left, j = right;

    if (right - left < 10)
    {
        for (int64 ii = left; ii <= right; ++ii)
            for (int64 jj = ii + 1; jj <= right; ++jj)
                if (geno_bucket.offset[ii].size > geno_bucket.offset[jj].size ||
                    geno_bucket.offset[ii].size == geno_bucket.offset[jj].size &&
                    structure_loc_original_idx[ii] > structure_loc_original_idx[jj])
                {
                    Swap(GetLoc(ii), GetLoc(jj));
                    Swap(geno_bucket.offset[ii], geno_bucket.offset[jj]);
                    Swap(structure_loc_original_idx[ii], structure_loc_original_idx[jj]);
                }

        PROGRESS_VALUE += right - left + 1;
        return;
    }

    int64 mid = (left + right) >> 1;
    uint64 midsize = geno_bucket.offset[mid].size;
    uint64 mididx = structure_loc_original_idx[mid];

    while (left < j || i < right)
    {
        while (geno_bucket.offset[i].size < midsize ||
            geno_bucket.offset[i].size == midsize &&
            structure_loc_original_idx[i] < mididx) i++;

        while (geno_bucket.offset[j].size > midsize ||
            geno_bucket.offset[j].size == midsize &&
            structure_loc_original_idx[j] > mididx) j--;

        if (i <= j)
        {
            Swap(GetLoc(i), GetLoc(j));
            Swap(geno_bucket.offset[i], geno_bucket.offset[j]);
            Swap(structure_loc_original_idx[i], structure_loc_original_idx[j]);
            i++; j--;
        }

        if (i > j)
        {
            if (i == j + 2)
                PROGRESS_VALUE++;

            if (left < j)
            {
                QUICKSORT_PARAMETER par = { left, j };
                Lock(GLOCK2);
                qslstack.Push(par);
                UnLock(GLOCK2);
            }
            else if (left == j)
                PROGRESS_VALUE++;

            if (i < right)
            {
                QUICKSORT_PARAMETER par = { i, right };
                Lock(GLOCK2);
                qslstack.Push(par);
                UnLock(GLOCK2);
            }
            else if (i == right)
                PROGRESS_VALUE++;

            return;
        }
    }
}

/* Quick sort locus by genotype bitsize */
THREAD(QSWorker2)
{
    QUICKSORT_PARAMETER par;
    while (PROGRESS_VALUE != PROGRESS_CEND)
    {
        bool hastask = false;

        Lock(GLOCK2);
        if (qslstack.size)
        {
            hastask = true;
            par = qslstack.Pop();
        }
        UnLock(GLOCK2);

        if (hastask)
            QSLocus2(par.left, par.right);
    }
}

/* Rearrange locus */
THREAD2(ArrangeLocus)
{
    //force Genotype bucket in ascending order
    {
        MEMORY* nlocus_memory = new MEMORY[g_nthread_val];
        SLOCUS* nslocus = new SLOCUS[nfilter];
        {
            BUCKET ngeno_bucket;
            ngeno_bucket.FilterLocusGT(&geno_bucket, true, g_nthread_val, nlocus_memory, nslocus, NULL);
            geno_bucket.Replace(ngeno_bucket);
        }
        DEL(locus_memory);  locus_memory = nlocus_memory;
        DEL(slocus);        slocus = nslocus;
    }

    //force Gtab in ascending order
    {
        MEMORY* nlocus_memory = new MEMORY[g_nthread_val];
        SLOCUS* nslocus = new SLOCUS[nloc];
        
        uint64 tsize = 0;
        for (int tid = 0; tid < g_nthread_val; ++tid)
            for (int bid = 0; bid < locus_memory[tid].nblocks; ++bid)
                tsize += locus_memory[tid].blocks[bid].used;
        
        nlocus_memory->blocks[0].size = tsize;
        DEL(nlocus_memory->blocks[0].bucket);
        
        nlocus_memory->blocks[0].bucket = new bool[tsize];
        for (int64 l = 0; l < nloc; ++l)
            new(&nslocus[l]) SLOCUS(nlocus_memory[0], slocus[l]);
        
        DEL(locus_memory);  locus_memory = nlocus_memory;
        DEL(slocus);        slocus = nslocus;
    }

    nloc = nfilter;
    CheckGenotypeId<REAL>();
}

/* Calculate bayesian clustering */
template<typename REAL>
TARGET void CalcBayesian()
{
    //evaluate the time usage of Bayesian clustering
    timepoint begin = GetNow();

    if (!structure) return;
    if (ad) Exit("\nError: Bayesian clustering (-structure) is incompatible with allelic depth (-ad) option.\n");

    EvaluationBegin();

    //Sort locus according to size
    structure_loc_original_idx = new int64[nloc];
    for (int64 l = 0; l < nloc; ++l)
        structure_loc_original_idx[l] = l;
    QUICKSORT_PARAMETER par = { 0, nloc - 1 };
    qslstack.Push(par);
    RunThreads(&QSWorker2, NULL, NULL, nloc, nloc, "\nSorting loci according to genotype bitsize:\n", g_nthread_val, true);
    DEL(structure_loc_original_idx);

    structure_loc_size_min = geno_bucket.offset[0].size;
    structure_loc_size_max = geno_bucket.offset[nloc - 1].size;
    structure_loc_lend = new int64[structure_loc_size_max + 1]; 
    SetZero(structure_loc_lend, structure_loc_size_max + 1);
    uint64 csize = structure_loc_size_min;
    for (int64 l = 0; l < nloc; ++l)
        while (geno_bucket.offset[l].size != csize)
        {
            structure_loc_lend[csize] = l;
            csize++;
        }
    structure_loc_lend[structure_loc_size_max] = nloc;

    RunThreads(&ArrangeLocus<REAL>, NULL, NULL, nloc, nloc, "\nArrange genotype bucket:\n", 1, true);

    //update allele_freq_offset
    KT = 0;
    if (allele_freq_offset == NULL) 
        allele_freq_offset = new uint64[nloc];
    for (int64 l = 0; l < nloc; ++l)
    {
        allele_freq_offset[l] = KT;
        KT += GetLoc(l).k;
    }

    //shuffle allele to reduce conflict for cores in different ccx
    structure_loc_coprime = nloc <= 67 ? 1 : 67;
    while (!IsCoprime(nloc, structure_loc_coprime)) structure_loc_coprime++;
    if (structure_loc_coprime >= nloc) structure_loc_coprime = 1;

    for (int lsize = structure_loc_size_min; lsize <= structure_loc_size_max; ++lsize)
    {
        int64 lstart = structure_loc_lend[lsize - 1], lend = structure_loc_lend[lsize];
        int64 lend0 = (lend - lstart + 63) / 64;
        structure_loc_coprime64[lsize] = lend0 <= 11 ? 1 : 11;
        while (!IsCoprime(lend0, structure_loc_coprime64[lsize])) structure_loc_coprime64[lsize]++;
        if (structure_loc_coprime >= nloc) structure_loc_coprime64[lsize] = 1;
    }

    //alloc tasks
    structure_nsubthread = (g_nthread_val + structure_nstream_val - 1) / structure_nstream_val;
    structure_totalruns = (structure_krange_max + 1 - structure_krange_min) * structure_nruns_val;
    structure_par = new STRUCTURE_RUNINFO[structure_totalruns];

    for (int i = structure_krange_min, lc = 0; i <= structure_krange_max; ++i)
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

    if (g_gpu_val == 2)
    {
        ResetDeviceCUDA();
        AllocMemoryCUDA();
        for (int devID = 0; devID < nGPU; ++devID)
            CopyStructureMemory(0);
    }

    int64 niter = ((structure_locpriori_val == 2 && structure_admix_val == 2) ? structure_nadmburnin_val : 0) + structure_nreps_val + structure_nburnin_val;
    int64 niterK = (structure_krange_max + structure_krange_min) * (structure_krange_max + 1 - structure_krange_min) / 2 * structure_nruns_val * niter;

    structure_nsimult = 1; //> 1 can increase GPU occupancy
    structure_navailable_cuda = 0;
    structure_cuda_taskid = new atomic<int>[nGPU * structure_nsimult];
    SetZero(structure_cuda_taskid, nGPU * structure_nsimult);

    if (structure_eval_val == 1)
    {
        printf("\nStructureFunctionOverhead %0.5f s\n", GetElapse(begin));
        begin = GetNow();
    }

    RunThreads(&BayesianThread<REAL>, NULL, NULL, niterK, niterK,
        "\nPerforming bayesian clustering:\n", structure_nstream_val, true);

    if (structure_eval_val == 1)
    {
        printf("\nStructureFunctionTotal %0.5f s\n", GetElapse(begin));
        begin = GetNow();
    }

    DEL(structure_cuda_taskid);

    BAYESIAN<REAL>::WriteStructureSummary(structure_par, structure_totalruns);

    DEL(structure_par);
    DEL(structure_loc_lend);

    if (g_gpu_val == 2)
    {
        for (int devID = 0; devID < nGPU; ++devID)
            FreeStructureMemory(devID);
        FreeMemoryCUDA();
        ResetDeviceCUDA();
    }

    EvaluationEnd("Bayesian clustering");

    if (structure_plot_val == 1)
        RunRscript("structure_plot.R");
}

/* Calculate Bayesian clustering using multiple threads */
THREAD2(BayesianThread)
{
    if (threadid < nGPU * structure_nsimult)
    {
        //GPU stream
        for (int i = structure_totalruns - 1; i >= 0; --i)
        {
            if (structure_par[i].flag.test_and_set()) continue;

            BAYESIAN<REAL> Structure;
            Structure.ReadPar(&structure_par[i]);
            Structure.iscudastream = true;
            Structure.MCMC();
            Structure.WriteStructure();
        }

        atomic<int>& navail = *(atomic<int>*) & structure_navailable_cuda;
        navail++;
    }
    else
    {
        //CPU stream
        for (int i = 0; i < structure_totalruns; ++i)
        {
            if (structure_par[i].flag.test_and_set()) continue;

            BAYESIAN<REAL> Structure;
            Structure.iscudastream = false;
            Structure.ReadPar(&structure_par[i]);
            Structure.MCMC();
            Structure.WriteStructure();
        }
    }
}
  