/* Parameter Functions */

#include "vcfpop.h"

struct indtab_node
{
	int level; //ind = 0, pop = 1, reg >= 2;
	string name;
	map<string, indtab_node*> subunits;

	indtab_node(string _name, int _level)
	{
		name = _name;
		level = _level;
	}

	void Write(string* buf)
	{
		if (level == 0) return;
		buf[level].append(" ");
		buf[level].append(name);
		int count = 0;
		for (auto kv : subunits)
		{
			buf[level].append(count++ == 0 ? ":" : ",");
			buf[level].append(kv.first);
			kv.second->Write(buf);
		}
	}
};

#define extern

extern vector<string> argv;

/* Parameter file */
extern bool p_b;								extern string p_val;

/* Global settings */
extern bool g_decimal_b;						extern int g_decimal_val;							extern char g_decimal_str[16];
extern bool g_scientific_b;						extern int g_scientific_val;
extern bool g_nthread_b;						extern int g_nthread_val;
extern bool g_simd_b;							extern int g_simd_val;
extern bool g_gpu_b;							extern int g_gpu_val;
extern bool g_float_b;							extern int g_float_val;
extern bool g_fastsingle_b;						extern int g_fastsingle_val;
extern bool g_seed_b;							extern int g_seed_val;
extern bool g_tmpdir_b;							extern string g_tmpdir_val;
extern bool g_progress_b;						extern int g_progress_val;

extern bool g_input_b;							extern string g_input_val;
extern int g_input_row, g_input_col;

extern bool g_format_b;							extern int g_format_val;
extern bool g_locusname_b;						extern int g_locusname_val;
extern bool g_extracol_b;						extern int g_extracol_val;
extern bool g_output_b;							extern string g_output_val;
extern bool g_indtext_b;						extern string g_indtext_val;
extern bool g_indtab_b;							extern string g_indtab_val;
extern bool g_delimiter_b;						extern char g_delimiter_val;
extern bool g_linebreak_b;						extern char* g_linebreak_val;
extern bool g_rscript_b;						extern string g_rscript_val;

extern bool g_benchmark_b;						extern int g_benchmark_val;
extern bool g_replot_b;							extern int g_replot_val;
extern bool g_eval_b;							extern int g_eval_val;
extern bool g_missingploidy_b;					extern byte g_missingploidy_val[N_MAX_OPTION];
extern bool g_maxlen_b;							extern int64 g_maxlen_val;

/* Filters */
extern bool f_filter;
extern bool f_qual_b;							extern double f_qual_min, f_qual_max;
extern bool f_type_b;							extern int f_type_val;
extern bool f_original_b;						extern int f_original_val;
extern bool f_chrprefix_b;						extern vector<string> f_chrprefix_val;
extern bool f_chrname_b;						extern vector<string> f_chrname_val;
extern bool f_pop_b;							extern string f_pop_val;
extern bool f_bmaf_b;							extern double f_bmaf_min, f_bmaf_max;
extern bool f_k_b;								extern int f_k_min, f_k_max;
extern bool f_n_b;								extern int f_n_min, f_n_max;
extern bool f_ptype_b;							extern double f_ptype_min, f_ptype_max;
extern bool f_pval_b;							extern double f_pval_min, f_pval_max;
extern bool f_model_b;							extern int f_model_val;
extern bool f_he_b;								extern double f_he_min, f_he_max;
extern bool f_ho_b;								extern double f_ho_min, f_ho_max;
extern bool f_pic_b;							extern double f_pic_min, f_pic_max;
extern bool f_ae_b;								extern double f_ae_min, f_ae_max;
extern bool f_I_b;								extern double f_I_min, f_I_max;
extern bool f_dp_b;								extern uint f_dp_min, f_dp_max;
extern bool f_gq_b;								extern int f_gq_min, f_gq_max;
extern bool f_ploidy_b;							extern int f_ploidy_min, f_ploidy_max;
extern bool f_itype_b;							extern double f_itype_min, f_itype_max;
extern bool f_iploidy_b;						extern int f_iploidy_min, f_iploidy_max;
extern bool f_windowsize_b;						extern int f_windowsize_val;
extern bool f_windowstat_b;						extern int f_windowstat_val;

/* Haplotype extraction */
extern bool haplotype;
extern bool haplotype_ptype_b;					extern double haplotype_ptype_min, haplotype_ptype_max;
extern bool haplotype_length_b;					extern int64 haplotype_length_min, haplotype_length_max;
extern bool haplotype_variants_b;				extern int haplotype_variants_min, haplotype_variants_max;
extern bool haplotype_interval_b;				extern int64 haplotype_interval_val;
extern bool haplotype_alleles_b;				extern int haplotype_alleles_min, haplotype_alleles_max;
extern bool haplotype_genotypes_b;				extern int haplotype_genotypes_min, haplotype_genotypes_max;

/* File conversion */
extern bool convert;
extern bool convert_format_b;					extern byte convert_format_val[N_MAX_OPTION];
extern bool convert_mode_b;						extern int convert_mode_val;

/* Individual statistics */
extern bool indstat;
extern bool indstat_type_b;						extern byte indstat_type_val[N_MAX_OPTION];
extern bool indstat_model_b;					extern byte indstat_model_val[N_MAX_OPTION];
extern bool indstat_estimator_b;				extern byte indstat_estimator_val[N_MAX_OPTION];
extern bool indstat_ref_b;						extern byte indstat_ref_val[N_MAX_OPTION];
extern bool indstat_locus_b;					extern byte indstat_locus_val[N_MAX_OPTION];

/* Genetic Diversity */
extern bool diversity;
extern bool diversity_level_b;					extern byte diversity_level_val[N_MAX_OPTION];
extern bool diversity_model_b;					extern byte diversity_model_val[N_MAX_OPTION];

/* Genetic differentiation */
extern bool fst;
extern bool fst_plot_b;							extern int fst_plot_val;
extern bool fst_level_b;						extern byte fst_level_val[N_MAX_OPTION];
extern bool fst_estimator_b;					extern byte fst_estimator_val[N_MAX_OPTION];
extern bool fst_fmt_b;							extern byte fst_fmt_val[N_MAX_OPTION];
extern bool fst_locus_b;						extern byte fst_locus_val[N_MAX_OPTION];
extern bool fst_test_b;							extern byte fst_test_val[N_MAX_OPTION];

/* Genetic distance */
extern bool gdist;
extern bool gdist_plot_b;						extern int gdist_plot_val;
extern bool gdist_level_b;						extern byte gdist_level_val[N_MAX_OPTION];
extern bool gdist_weightmissing_b;				extern int gdist_weightmissing_val;
extern bool gdist_estimator_b;					extern byte gdist_estimator_val[N_MAX_OPTION];
extern bool gdist_fmt_b;						extern byte gdist_fmt_val[N_MAX_OPTION];

/* Analysis of molecular variances */
extern bool amova;
extern bool amova_method_b;						extern byte amova_method_val[N_MAX_OPTION];					extern int amova_cmethod_val;
extern bool amova_mutation_b;					extern byte amova_mutation_val[N_MAX_OPTION];				extern int amova_cmutation_val;
extern bool amova_ind_b;						extern byte amova_ind_val[N_MAX_OPTION];					extern int amova_cind_val;
extern bool amova_trunc_b;						extern int amova_trunc_val;
extern bool amova_test_b;						extern int amova_test_val;
extern bool amova_nperm_b;						extern int amova_nperm_val;
extern bool amova_pseudo_b;						extern int amova_pseudo_val;
extern bool amova_printss_b;					extern int amova_printss_val;

/* Sliding window */
extern bool slide;
extern bool slide_plot_b;						extern int slide_plot_val;
extern bool slide_plot_columns_b;				extern int slide_plot_columns_val[N_MAX_SLIDEPLOT]; //values begins from 1
extern bool slide_plot_styles_b;				extern int slide_plot_styles_val[N_MAX_SLIDEPLOT];
extern bool slide_windowsize_b;					extern int slide_windowsize_val;
extern bool slide_windowstep_b;					extern int slide_windowstep_val;
extern bool slide_minvariants_b;				extern int slide_minvariants_val;
extern bool slide_estimator_b;					extern byte slide_estimator_val[N_MAX_OPTION];
extern bool slide_pop_b;						extern string slide_pop_val;

/* Population assignment */
extern bool popas;
extern bool popas_plot_b;						extern int popas_plot_val;
extern bool popas_model_b;						extern byte popas_model_val[N_MAX_OPTION];
extern bool popas_level_b;						extern byte popas_level_val[N_MAX_OPTION];
extern bool popas_error_b;						extern double popas_error_val;

/* Relatedness coefficient estimation */
extern bool relatedness;
extern bool relatedness_plot_b;					extern int relatedness_plot_val;
extern bool relatedness_range_b;				extern byte relatedness_range_val[N_MAX_OPTION];
extern bool relatedness_fmt_b;					extern byte relatedness_fmt_val[N_MAX_OPTION];
extern bool relatedness_estimator_b;			extern byte relatedness_estimator_val[N_MAX_OPTION];

/* Kinship coefficient estimation */
extern bool kinship;
extern bool kinship_plot_b;						extern int kinship_plot_val;
extern bool kinship_range_b;					extern byte kinship_range_val[N_MAX_OPTION];
extern bool kinship_fmt_b;						extern byte kinship_fmt_val[N_MAX_OPTION];
extern bool kinship_estimator_b;				extern byte kinship_estimator_val[N_MAX_OPTION];

/* Principal coordinate analysis */
extern bool pcoa;
extern bool pcoa_plot_b;						extern int pcoa_plot_val;
extern bool pcoa_level_b;						extern byte pcoa_level_val[N_MAX_OPTION];
extern bool pcoa_dim_b;							extern int pcoa_dim_val;
extern bool pcoa_estimator_b;					extern byte pcoa_estimator_val[N_MAX_OPTION];

/* Hierarchical clustering */
extern bool cluster;
extern bool cluster_plot_b;						extern int cluster_plot_val;
extern bool cluster_level_b;					extern byte cluster_level_val[N_MAX_OPTION];
extern bool cluster_method_b;					extern byte cluster_method_val[N_MAX_OPTION];
extern bool cluster_estimator_b;				extern byte cluster_estimator_val[N_MAX_OPTION];

/* Bayesian clustering */
extern bool structure;
extern bool structure_nstream_b;				extern int structure_nstream_val;
extern bool structure_plot_b;					extern int structure_plot_val;
extern bool structure_eval_b;					extern int structure_eval_val;
extern bool structure_writelnl_b;				extern int structure_writelnl_val;

/* Bayesian clustering: Model */
extern bool structure_admix_b;					extern int structure_admix_val;
extern bool structure_locpriori_b;				extern int structure_locpriori_val;
extern bool structure_f_b;						extern int structure_f_val;

/* Bayesian clustering: MCMC */
extern bool structure_krange_b;					extern int structure_krange_min, structure_krange_max;
extern bool structure_nburnin_b;				extern int structure_nburnin_val;
extern bool structure_nreps_b;					extern int structure_nreps_val;
extern bool structure_nthinning_b;				extern int structure_nthinning_val;
extern bool structure_nruns_b;					extern int structure_nruns_val;
extern bool structure_nadmburnin_b;				extern int structure_nadmburnin_val;

/* Bayesian clustering: Misc */
extern bool structure_lambda_b;					extern double structure_lambda_val;
extern bool structure_stdlambda_b;				extern double structure_stdlambda_val;
extern bool structure_maxlambda_b;				extern double structure_maxlambda_val;
extern bool structure_inferlambda_b;			extern int structure_inferlambda_val;
extern bool structure_difflambda_b;				extern int structure_difflambda_val;
extern bool structure_diversity_b;				extern int structure_diversity_val;

/* Bayesian clustering: Admix */
extern bool structure_alpha_b;					extern double structure_alpha_val;
extern bool structure_stdalpha_b;				extern double structure_stdalpha_val;
extern bool structure_maxalpha_b;				extern double structure_maxalpha_val;
extern bool structure_inferalpha_b;				extern int structure_inferalpha_val;
extern bool structure_diffalpha_b;				extern int structure_diffalpha_val;
extern bool structure_uniformalpha_b;			extern int structure_uniformalpha_val;
extern bool structure_alphapriora_b;			extern double structure_alphapriora_val;
extern bool structure_alphapriorb_b;			extern double structure_alphapriorb_val;
extern bool structure_metrofreq_b;				extern int structure_metrofreq_val;

/* Bayesian clustering: locpriori */
extern bool structure_r_b;						extern double structure_r_val;
extern bool structure_maxr_b;					extern double structure_maxr_val;
extern bool structure_stdr_b;					extern double structure_stdr_val;
extern bool structure_epseta_b;					extern double structure_epseta_val;
extern bool structure_epsgamma_b;				extern double structure_epsgamma_val;

/* Bayesian clustering: fmodel */
extern bool structure_pmeanf_b;					extern double structure_pmeanf_val;
extern bool structure_pstdf_b;					extern double structure_pstdf_val;
extern bool structure_stdf_b;					extern double structure_stdf_val;
extern bool structure_singlef_b;				extern int structure_singlef_val;

/* Ploidy inference from allelic depth distribution, test, do not use */
extern bool ploidyinfer;
extern bool ploidyinfer_type_b;					extern byte ploidyinfer_type_val[N_MAX_OPTION];
extern bool ploidyinfer_histogram_b;			extern int ploidyinfer_histogram_val;
extern bool ploidyinfer_nbins_b;				extern int ploidyinfer_nbins_val;

/* Allelic depth analysis, test, do not use */
extern int ad;

/* Analysis of spatial structure, test, do not use */
extern bool spa;
extern bool spa_dim_b;							extern int spa_dim_val;
extern bool spa_level_b;						extern int spa_level_val;
extern bool spa_odepth_b;						extern int spa_odepth_val;
extern bool spa_ofreq_b;						extern int spa_ofreq_val;
extern bool spa_coord_b;						extern string spa_coord_val;
extern bool spa_truncate_b;						extern double spa_truncate_val, spa_truncate_val2;//n
#undef extern


/* Replace path deliminator*/
TARGET string ReplacePathDelim(string& str, bool ispath)
{
	static string t1{ PATH_DELIM_REVERSE }, t2{ PATH_DELIM };
	static string t3{ PATH_DOUBLE_DELIM_REVERSE }, t4{ PATH_DOUBLE_DELIM };
	string re;
	re = TrimParQuote(str);
	re = ReplaceStr(re, t1, t2);
	re = ReplaceStr(re, t3, t4);
	if (ispath && re[re.size() - 1] != PATH_DELIM)
		re.push_back(PATH_DELIM);
	return re;
}

/* Initialize parameters */
TARGET void SetDefaultParameters()
{
#ifdef __aarch64__
	SIMD_TYPE = 2;
#else
	SIMD_TYPE = 3;
#endif
	srand((unsigned int)time(NULL));

	InitLock(GLOCK1);
	InitLock(GLOCK2);

	CURDIR = GetCurDir();

	diversity_filter = 0;
	genotype_filter = 0;
	individual_filter = 0;
	info_filter = 0;

	ad = 0;

	spa = false;
	spa_dim_b = false;				spa_dim_val = 2;
	spa_level_b = false;			spa_level_val = 1;
	spa_odepth_b = false;			spa_odepth_val = 2;
	spa_ofreq_b = false;			spa_ofreq_val = 2;
	spa_coord_b = false;			spa_coord_val = "";
	spa_truncate_b = false;			spa_truncate_val = 0.001;		spa_truncate_val2 = 1.0 - spa_truncate_val;

	f_filter = false;
	f_qual_b = false;				f_qual_min = 0;					f_qual_max = 0;
	f_type_b = false;				f_type_val = 0;
	f_original_b = false;			f_original_val = 0;
	f_chrprefix_b = false;			
	f_chrname_b = false;			
	f_pop_b = false;				f_pop_val = "total";
	f_bmaf_b = false;				f_bmaf_min = 0;					f_bmaf_max = 0;
	f_k_b = false;					f_k_min = 0;					f_k_max = 0;
	f_n_b = false;					f_n_min = 0;					f_n_max = 0;
	f_ptype_b = false;				f_ptype_min = 0; 				f_ptype_max = 0;
	f_pval_b = false; 				f_pval_min = 0; 				f_pval_max = 0;
	f_model_b = false; 				f_model_val = 0;
	f_he_b = false; 				f_he_min = 0; 					f_he_max = 0;
	f_ho_b = false; 				f_ho_min = 0; 					f_ho_max = 0;
	f_pic_b = false; 				f_pic_min = 0; 					f_pic_max = 0;
	f_ae_b = false; 				f_ae_min = 0; 					f_ae_max = 0;
	f_I_b = false; 					f_I_min = 0; 					f_I_max = 0;
	f_dp_b = false; 				f_dp_min = 0; 					f_dp_max = 0;
	f_gq_b = false; 				f_gq_min = 0; 					f_gq_max = 0;
	f_ploidy_b = false; 			f_ploidy_min = 0; 				f_ploidy_max = 0;
	f_itype_b = false;				f_itype_min = 0;				f_itype_max = 0;
	f_iploidy_b = false;			f_iploidy_min = 0;				f_iploidy_max = 0;
	f_windowsize_b = false;			f_windowsize_val = 0;
	f_windowstat_b = false;			f_windowstat_val = 5;

	g_decimal_b = false;			g_decimal_val = 5;				sprintf(g_decimal_str, "%%0.%df", 5);
	g_scientific_b = false;			g_scientific_val = 2;
	g_nthread_b = false;			g_nthread_val = 2;
#ifdef __aarch64__ 
	g_simd_b = false;				g_simd_val = 1;
#else
	g_simd_b = false;				g_simd_val = 3;
#endif
	g_gpu_b = false;				g_gpu_val = 1;
	g_float_b = false;				g_float_val = 2;
	g_fastsingle_b = false;			g_fastsingle_val = 1;
	g_seed_b = false;				g_seed_val = rand();
	g_tmpdir_b = false;				g_tmpdir_val = CURDIR;
	g_progress_b = false;			g_progress_val = 80;

	g_input_b = false;				g_input_val = "";
	g_input_row = 0;				g_input_col = 0;
	
	g_format_b = false;				g_format_val = 1;
	g_locusname_b = false;			g_locusname_val = 3;
	g_extracol_b = false;			g_extracol_val = 0;
	g_output_b = false;				g_output_val = "vcfpop.out";
	g_indtext_b = false;			g_indtext_val = "";
	g_indtab_b = false;				g_indtab_val = "";
	g_delimiter_b = false;			g_delimiter_val = '\t';
	g_linebreak_b = false;			g_linebreak_val = (char*)"\n";

	g_rscript_b = false;			g_rscript_val = "";
	g_benchmark_b = false;			g_benchmark_val = 2;
	g_replot_b = false;				g_replot_val = 2;
	g_eval_b = false;				g_eval_val = 2;
	g_missingploidy_b = false;		SetZero(g_missingploidy_val, N_MAX_OPTION);
	g_maxlen_b = false;				g_maxlen_val = 99999999999999999;

	haplotype = false;
	haplotype_ptype_b = false;		haplotype_ptype_min = 0.8;						haplotype_ptype_max = 1;
	haplotype_length_b = false;		haplotype_length_min = 1ull;					haplotype_length_max = 1000000ull;
	haplotype_variants_b = false;	haplotype_variants_min = 5;						haplotype_variants_max = 20;
	haplotype_interval_b = false;	haplotype_interval_val = 0ll;
	haplotype_alleles_b = false;	haplotype_alleles_min = 2;						haplotype_alleles_max = 65535;
	haplotype_genotypes_b = false;	haplotype_genotypes_min = 2;					haplotype_genotypes_max = 65535;

	convert = false;
	convert_format_b = false;		SetZero(convert_format_val, N_MAX_OPTION);		convert_format_val[1] = 2;
	convert_mode_b = false;		    convert_mode_val = 1;

	indstat = false;
	indstat_type_b = false;			SetZero(indstat_type_val, N_MAX_OPTION);		indstat_type_val[1] = indstat_type_val[2] = 1;
	indstat_model_b = false;		SetZero(indstat_model_val, N_MAX_OPTION);		indstat_model_val[1] = 1;
	indstat_estimator_b = false;	SetZero(indstat_estimator_val, N_MAX_OPTION);	indstat_estimator_val[1] = 1;
	indstat_ref_b = false;			SetZero(indstat_ref_val, N_MAX_OPTION);			indstat_ref_val[3] = 1;
	indstat_locus_b = false;		SetZero(indstat_locus_val, N_MAX_OPTION);		indstat_locus_val[1] = 1;

	diversity = false;
	diversity_level_b = false;		SetZero(diversity_level_val, N_MAX_OPTION);		diversity_level_val[1] = diversity_level_val[4] = 1;
	diversity_model_b = false;		SetZero(diversity_model_val, N_MAX_OPTION);		diversity_model_val[1] = 1;

	fst = false;
	fst_level_b = false;			SetZero(fst_level_val, N_MAX_OPTION);			fst_level_val[5] = 1;
	fst_estimator_b = false;		SetZero(fst_estimator_val, N_MAX_OPTION);		fst_estimator_val[1] = 1;
	fst_fmt_b = false;				SetZero(fst_fmt_val, N_MAX_OPTION);				fst_fmt_val[1] = 1;
	fst_locus_b = false;			SetZero(fst_locus_val, N_MAX_OPTION);			fst_locus_val[1] = 1;
	fst_test_b = false;				SetZero(fst_test_val, N_MAX_OPTION);

	gdist = false;
	gdist_level_b = false;			SetZero(gdist_level_val, N_MAX_OPTION);			gdist_level_val[2] = 1;
	gdist_weightmissing_b = false;	gdist_weightmissing_val = 1;
	gdist_estimator_b = false;		SetZero(gdist_estimator_val, N_MAX_OPTION);		gdist_estimator_val[1] = 1;
	gdist_fmt_b = false;			SetZero(gdist_fmt_val, N_MAX_OPTION);			gdist_fmt_val[1] = 1;

	amova = false;
	amova_method_b = false;			amova_method_val[1] = 1;
	amova_mutation_b = false;		amova_mutation_val[1] = 1;
	amova_ind_b = false;			amova_ind_val[1] = 1;
	amova_trunc_b = false;			amova_trunc_val = 2;
	amova_test_b = false;			amova_test_val = 1;
	amova_nperm_b = false;			amova_nperm_val = 9999;
	amova_pseudo_b = false;			amova_pseudo_val = 50;
	amova_printss_b = false;		amova_printss_val = 2;

	slide = false;
	slide_plot_b = false;			slide_plot_val = 2;
	slide_plot_columns_b = false;	slide_plot_columns_val[0] = 1; slide_plot_columns_val[1] = 2; slide_plot_columns_val[2] = 3; slide_plot_columns_val[3] = 4; slide_plot_columns_val[4] = 5;
	slide_plot_styles_b = false;	slide_plot_styles_val[0] = 1;  slide_plot_styles_val[1] = 2;  slide_plot_styles_val[2] = 3;  slide_plot_styles_val[3] = 4;  slide_plot_styles_val[4] = 5;
	slide_windowsize_b = false;		slide_windowsize_val = 1000000;
	slide_windowstep_b = false;		slide_windowstep_val = 100000;
	slide_minvariants_b = false;	slide_minvariants_val = 10;
	slide_estimator_b = false;		SetZero(slide_estimator_val, N_MAX_OPTION);	    slide_estimator_val[1] = 1;
	slide_pop_b = false;			slide_pop_val = "total";

	popas = false;
	popas_model_b = false;			SetZero(popas_model_val, N_MAX_OPTION);			popas_model_val[1] = 1;
	popas_level_b = false;			SetZero(popas_level_val, N_MAX_OPTION);			popas_level_val[1] = 1;
	popas_error_b = false;			popas_error_val = 0.01;

	relatedness = false;
	relatedness_range_b = false;	SetZero(relatedness_range_val, N_MAX_OPTION);	relatedness_range_val[3] = 1;
	relatedness_fmt_b = false;		SetZero(relatedness_fmt_val, N_MAX_OPTION);		relatedness_fmt_val[1] = 1;
	relatedness_estimator_b = false; SetZero(relatedness_estimator_val, N_MAX_OPTION); relatedness_estimator_val[10] = 1;

	kinship = false;
	kinship_range_b = false;		SetZero(kinship_range_val, N_MAX_OPTION);		kinship_range_val[3] = 1;
	kinship_fmt_b = false;			SetZero(kinship_fmt_val, N_MAX_OPTION);			kinship_fmt_val[1] = 1;
	kinship_estimator_b = false;	SetZero(kinship_estimator_val, N_MAX_OPTION);	kinship_estimator_val[9] = 1;

	pcoa = false;
	pcoa_plot_b = false;			pcoa_plot_val = 2;
	pcoa_level_b = false;			SetZero(pcoa_level_val, N_MAX_OPTION);			pcoa_level_val[1] = 1;
	pcoa_dim_b = false;				pcoa_dim_val = 3;
	pcoa_estimator_b = false;		SetZero(pcoa_estimator_val, N_MAX_OPTION);		pcoa_estimator_val[1] = 1;

	cluster = false;
	cluster_plot_b = false;			cluster_plot_val = 2;
	cluster_level_b = false;		SetZero(cluster_level_val, N_MAX_OPTION);		cluster_level_val[1] = 1;
	cluster_method_b = false;		SetZero(cluster_method_val, N_MAX_OPTION);		cluster_method_val[3] = 1;
	cluster_estimator_b = false;	SetZero(cluster_estimator_val, N_MAX_OPTION);	cluster_estimator_val[1] = 1;

	structure = false;
	structure_nstream_b = false;    structure_nstream_val = 4;
	structure_plot_b = false;		structure_plot_val = 2;
	structure_eval_b = false;		structure_eval_val = 2;
	structure_writelnl_b = false;	structure_writelnl_val = 2;

	structure_admix_b = false;		structure_admix_val = 2;
	structure_locpriori_b = false;	structure_locpriori_val = 2;
	structure_f_b = false;			structure_f_val = 2;
	structure_krange_b = false;		structure_krange_min = 1;						structure_krange_max = 5;

	structure_nburnin_b = false;	structure_nburnin_val = 1000;
	structure_nreps_b = false;		structure_nreps_val = 10000;
	structure_nthinning_b = false;	structure_nthinning_val = 10;
	structure_nruns_b = false;		structure_nruns_val = 1;
	structure_nadmburnin_b = false; structure_nadmburnin_val = 500;

	structure_lambda_b = false;		structure_lambda_val = 1;
	structure_stdlambda_b = false;	structure_stdlambda_val = 0.3;
	structure_maxlambda_b = false;	structure_maxlambda_val = 10;
	structure_inferlambda_b = false;structure_inferlambda_val = 2;
	structure_difflambda_b = false; structure_difflambda_val = 2;
	structure_diversity_b = false;	structure_diversity_val = 2;

	structure_alpha_b = false;		structure_alpha_val = 1;
	structure_stdalpha_b = false;	structure_stdalpha_val = 0.025;
	structure_maxalpha_b = false;	structure_maxalpha_val = 10;
	structure_inferalpha_b = false;	structure_inferalpha_val = 1;
	structure_diffalpha_b = false;	structure_diffalpha_val = 0;
	structure_uniformalpha_b = false; structure_uniformalpha_val = 1;
	structure_alphapriora_b = false; structure_alphapriora_val = 0.05;
	structure_alphapriorb_b = false; structure_alphapriorb_val = 0.001;

	structure_metrofreq_b = false;	structure_metrofreq_val = 10;
	structure_r_b = false;			structure_r_val = 1;
	structure_maxr_b = false;		structure_maxr_val = 20;
	structure_stdr_b = false;		structure_stdr_val = 0.1;
	structure_epseta_b = false;		structure_epseta_val = 0.025;
	structure_epsgamma_b = false;	structure_epsgamma_val = 0.025;

	structure_pmeanf_b = false;		structure_pmeanf_val = 0.01;
	structure_pstdf_b = false;		structure_pstdf_val = 0.05;
	structure_stdf_b = false;		structure_stdf_val = 0.05;
	structure_singlef_b = false;	structure_singlef_val = 2;

	ploidyinfer = false;
	ploidyinfer_type_b = false;		SetZero(ploidyinfer_type_val, N_MAX_OPTION);	ploidyinfer_type_val[2] = ploidyinfer_type_val[4] = 1;
	ploidyinfer_histogram_b = false; ploidyinfer_histogram_val = 2;
	ploidyinfer_nbins_b = false;	ploidyinfer_nbins_val = 20;

	diversity_filter = false;
	genotype_filter = false;
	individual_filter = false;
	info_filter = false;

	nGPU = 0;
}

/* Set parameters */
TARGET void SetParameters(bool isparfile)
{
	if (!isparfile) SetDefaultParameters();

	for (int i = 0; i < argv.size(); ++i)
	{
		if (!isparfile && !LwrLineCmp("-p=", argv[i]))
		{
			if (p_b) Exit("\nError: parameter %s has been assigned twice.\n", argv[i].c_str());
			p_b = true;
			p_val = ReplacePathDelim(argv[i], false);
		}
		else if (!LwrLineCmp("-ad", argv[i]))
		{
			ad = 1;
		}
		else if (!LwrLineCmp("-spa", argv[i]))
		{
			if (!LwrStrCmp("-spa", argv[i]))
				GetParBool(argv[i], spa);
			else if (!LwrLineCmp("-spa_dim=", argv[i]))
				GetParInteger(argv[i], spa_dim_b, spa_dim_val, 1, 15);
			else if (!LwrLineCmp("-spa_level=", argv[i]))
				GetParString(argv[i], "ind|pop", spa_level_b, spa_level_val);
			else if (!LwrLineCmp("-spa_odepth=", argv[i]))
				GetParString(argv[i], "yes|no", spa_odepth_b, spa_odepth_val);
			else if (!LwrLineCmp("-spa_ofreq=", argv[i]))
				GetParString(argv[i], "yes|no", spa_ofreq_b, spa_ofreq_val);
			else if (!LwrLineCmp("-spa_truncate=", argv[i]))
			{
				GetParDouble(argv[i], spa_truncate_b, spa_truncate_val, 0.0000000001, 0.01);
				spa_truncate_val2 = 1.0 - spa_truncate_val;
			}
			else if (!LwrLineCmp("-spa_coord=", argv[i]))
			{
				if (spa_coord_b) Exit("\nError: parameter %s has been assigned twice.\n", argv[i].c_str());
				spa_coord_b = true;
				spa_coord_val = ReplaceStr(argv[i].substr(11), "#n", "\n");
			}
		}
		else if (!LwrLineCmp("-f_", argv[i]) || !LwrStrCmp("-f", argv[i]))
		{
			if (!LwrStrCmp("-f", argv[i]))
				GetParBool(argv[i], f_filter);
			else if (!LwrLineCmp("-f_qual=", argv[i]))
				GetRangeParDouble(argv[i], f_qual_b, f_qual_min, f_qual_max, 0, 255);
			else if (!LwrLineCmp("-f_type=", argv[i]))
				GetParString(argv[i], "snp|indel|both", f_type_b, f_type_val);
			else if (!LwrLineCmp("-f_original=", argv[i]))
				GetParString(argv[i], "yes|no", f_original_b, f_original_val);
			else if (!LwrLineCmp("-f_chrprefix=", argv[i]))
			{
				if (f_chrprefix_b) Exit("\nError: parameter %s has been assigned twice.\n", argv[i].c_str());
				f_chrprefix_b = true;
				f_chrprefix_val = SplitStr(TrimParQuote(argv[i]), ",");
			}
			else if (!LwrLineCmp("-f_chrname=", argv[i]))
			{
				if (f_chrname_b) Exit("\nError: parameter %s has been assigned twice.\n", argv[i].c_str());
				f_chrname_b = true;
				f_chrname_val = SplitStr(TrimParQuote(argv[i]), ",");
			}
			else if (!LwrLineCmp("-f_pop=", argv[i]))
			{
				if (f_pop_b) Exit("\nError: parameter %s has been assigned twice.\n", argv[i].c_str());
				f_pop_b = true;
				f_pop_val = TrimParQuote(argv[i]);
			}
			else if (!LwrLineCmp("-f_bmaf=", argv[i]))
				GetRangeParDouble(argv[i], f_bmaf_b, f_bmaf_min, f_bmaf_max, 0, 0.5);
			else if (!LwrLineCmp("-f_k=", argv[i]))
				GetRangeParInteger(argv[i], f_k_b, f_k_min, f_k_max, 1, 65535);
			else if (!LwrLineCmp("-f_n=", argv[i]))
				GetRangeParInteger(argv[i], f_n_b, f_n_min, f_n_max, 0, 65535);
			else if (!LwrLineCmp("-f_ptype=", argv[i]))
				GetRangeParDouble(argv[i], f_ptype_b, f_ptype_min, f_ptype_max, 0, 1);
			else if (!LwrLineCmp("-f_pval=", argv[i]))
				GetRangeParDouble(argv[i], f_pval_b, f_pval_min, f_pval_max, 0, 1);
			else if (!LwrLineCmp("-f_model=", argv[i]))
				GetParString(argv[i], "rcs|prcs|ces|pes", f_model_b, f_model_val);
			else if (!LwrLineCmp("-f_he=", argv[i]))
				GetRangeParDouble(argv[i], f_he_b, f_he_min, f_he_max, 0, 1);
			else if (!LwrLineCmp("-f_ho=", argv[i]))
				GetRangeParDouble(argv[i], f_ho_b, f_ho_min, f_ho_max, 0, 1);
			else if (!LwrLineCmp("-f_pic=", argv[i]))
				GetRangeParDouble(argv[i], f_pic_b, f_pic_min, f_pic_max, 0, 1);
			else if (!LwrLineCmp("-f_ae=", argv[i]))
				GetRangeParDouble(argv[i], f_ae_b, f_ae_min, f_ae_max, 0, 99999999);
			else if (!LwrLineCmp("-f_I=", argv[i]))
				GetRangeParDouble(argv[i], f_I_b, f_I_min, f_I_max, 0, 99999999);
			else if (!LwrLineCmp("-f_dp=", argv[i]))
				GetRangeParInteger(argv[i], f_dp_b, f_dp_min, f_dp_max, 0, 0x7FFFFFFF);
			else if (!LwrLineCmp("-f_gq=", argv[i]))
				GetRangeParInteger(argv[i], f_gq_b, f_gq_min, f_gq_max, 0, 100);
			else if (!LwrLineCmp("-f_ploidy=", argv[i]))
				GetRangeParInteger(argv[i], f_ploidy_b, f_ploidy_min, f_ploidy_max, 1, N_MAX_PLOIDY);
			else if (!LwrLineCmp("-f_itype=", argv[i]))
				GetRangeParDouble(argv[i], f_itype_b, f_itype_min, f_itype_max, 0, 1);
			else if (!LwrLineCmp("-f_iploidy=", argv[i]))
				GetRangeParInteger(argv[i], f_iploidy_b, f_iploidy_min, f_iploidy_max, 1, N_MAX_PLOIDY);
			else if (!LwrLineCmp("-f_windowsize=", argv[i]))
				GetParInteger(argv[i], f_windowsize_b, f_windowsize_val, 1, 1000000000);
			else if (!LwrLineCmp("-f_windowstat=", argv[i]))
				GetParString(argv[i], "bmaf|k|n|ptype|he|ho|pic|ae", f_windowstat_b, f_windowstat_val);
			else
				Exit("\nError: Unrecognized parameter: %s\n", argv[i].c_str());
		}
		else if (!LwrLineCmp("-g_", argv[i]))
		{
			if (!LwrLineCmp("-g_decimal=", argv[i]))
			{
				GetParInteger(argv[i], g_decimal_b, g_decimal_val, 1, 15);
				sprintf(g_decimal_str, g_scientific_val == 1 ? "%%0.%de" : "%%0.%df", g_decimal_val);
			}
			else if (!LwrLineCmp("-g_scientific=", argv[i]))
			{
				GetParString(argv[i], "yes|no", g_scientific_b, g_scientific_val);
				sprintf(g_decimal_str, g_scientific_val == 1 ? "%%0.%de" : "%%0.%df", g_decimal_val);
			}
			else if (!LwrLineCmp("-g_nthread=", argv[i]))
				GetParInteger(argv[i], g_nthread_b, g_nthread_val, 1, 4096);
			else if (!LwrLineCmp("-g_simd=", argv[i]))
			{
#ifdef __aarch64__
				GetParString(argv[i], "none|neon", g_simd_b, g_simd_val);
#else
				GetParString(argv[i], "none|sse|avx|avx512", g_simd_b, g_simd_val);
#endif
				SIMD_TYPE = g_simd_val;
			}
			else if (!LwrLineCmp("-g_gpu=", argv[i]))
			{
				GetParString(argv[i], "none|cuda", g_gpu_b, g_gpu_val);
				if (g_gpu_val == 2)
				{
					nGPU = GetDeviceCountCUDA();
					if (nGPU == 0) Exit("\nError: cannot found CUDA devices, check graphic card driver or use -g_gpu=none.\n");
				}
				else nGPU = 0;
			}
			else if (!LwrLineCmp("-g_float=", argv[i]))
				GetParString(argv[i], "single|double", g_float_b, g_float_val);
			else if (!LwrLineCmp("-g_fastsingle=", argv[i]))
				GetParString(argv[i], "yes|no", g_fastsingle_b, g_fastsingle_val);
			else if (!LwrLineCmp("-g_seed=", argv[i]))
			{
				GetParInteger(argv[i], g_seed_b, g_seed_val, 0, 0x7FFFFFFF);
				if (g_seed_val == 0) g_seed_val = (int)time(NULL);
			}
			else if (!LwrLineCmp("-g_tmpdir=", argv[i]))
			{
				if (g_tmpdir_b) Exit("\nError: parameter %s has been assigned twice.\n", argv[i].c_str());
				g_tmpdir_b = true;
				g_tmpdir_val = ReplacePathDelim(argv[i], true);

				if (!FileExists(g_tmpdir_val.c_str()))
					DirectoryCreate(g_tmpdir_val.c_str());
				if (!FileExists(g_tmpdir_val.c_str()))
					Exit("\nError: cannot create temp directory %s.\n", g_tmpdir_val.c_str());
			}
			else if (!LwrLineCmp("-g_progress=", argv[i]))
				GetParInteger(argv[i], g_progress_b, g_progress_val, 10, 1000000);
			else if (!LwrLineCmp("-g_input=", argv[i]))
			{
				if (g_input_b) Exit("\nError: parameter %s has been assigned twice.\n", argv[i].c_str());
				g_input_b = true;
				g_input_val = ReplacePathDelim(argv[i], false);

				vector<string> ts = SplitStr(g_input_val, "&");
				g_input_row = (int)ts.size();
				g_input_col = CountChar(ts[0], '|') + 1;

				TOTLEN_COMPRESS = 0;
				FILE_INFO.clear();
				for (int ii = 0; ii < g_input_row; ++ii)
				{
					vector<string> row_files = SplitStr(ts[ii], "|");
					FILE_INFO.push_back(vector<FILEINFO>());

					if (row_files.size() != g_input_col)
						Exit("\nError: multiple input vcf files at each row should have the same number of columns.\n");

					for (int j = 0; j < g_input_col; ++j)
					{
						FILEINFO tf;
						path tpath(row_files[j]);
						tf.path = absolute(tpath).string();
						FILE_INFO[ii].push_back(tf);
					}
				}
			}
			else if (!LwrLineCmp("-g_format=", argv[i]))
				GetParString(argv[i], "vcf|bcf|genepop|spagedi|cervus|arlequin|structure|polygene|polyrelatedness|genodive|plink", g_format_b, g_format_val);
			else if (!LwrLineCmp("-g_locusname=", argv[i]))
				GetParString(argv[i], "chr|pos|chr_pos|chr_ref_alt|pos_ref_alt|chr_pos_ref_alt", g_locusname_b, g_locusname_val);
			else if (!LwrLineCmp("-g_extracol=", argv[i]))
				GetParInteger(argv[i], g_extracol_b, g_extracol_val, 0, 4096);
			else if (!LwrLineCmp("-g_output=", argv[i]))
			{
				if (g_output_b) Exit("\nError: parameter %s has been assigned twice.\n", argv[i].c_str());
				g_output_b = true;
				g_output_val = ReplacePathDelim(argv[i], false);

				path tpath(g_output_val);
				OUTFILE = absolute(tpath).string();
				OUTDIR = absolute(tpath).parent_path().string();
				OUTDIR.push_back(PATH_DELIM);
				if (!FileExists(OUTDIR.c_str()))
					DirectoryCreate(OUTDIR.c_str());
			}
			else if (!LwrLineCmp("-g_indtext=", argv[i]))
			{
				if (g_indtext_b) Exit("\nError: parameter %s has been assigned twice.\n", argv[i].c_str());
				if (g_indtab_b) Exit("\nError: options -g_indtab and -g_indtext are exclusive.\n");
				g_indtext_b = true;
				g_indtext_val = ReplaceStr(argv[i].substr(11), "#n", "\n");
			}
			else if (!LwrLineCmp("-g_indtab=", argv[i]))
			{
				if (g_indtab_b) Exit("\nError: parameter %s has been assigned twice.\n", argv[i].c_str());
				if (g_indtext_b) Exit("\nError: options -g_indtab and -g_indtext are exclusive.\n");
				g_indtab_b = true;
				g_indtab_val = ReplaceStr(argv[i].substr(10), "#n", "\n");
				g_indtext_b = true;
				ConvertIndTable();
			}
			else if (!LwrLineCmp("-g_delimiter=", argv[i]))
			{
				int tval = 0;
				GetParString(argv[i], "comma|tab", g_delimiter_b, tval);
				if (tval == 1) g_delimiter_val = ',';
				else g_delimiter_val = '\t';
			}
			else if (!LwrLineCmp("-g_linebreak=", argv[i]))
			{
				int tval = 0;
				GetParString(argv[i], "unix|win", g_linebreak_b, tval);
				if (tval == 1) g_linebreak_val = (char*)"\n";
				else g_linebreak_val = (char*)"\r\n";
			}
			else if (!LwrLineCmp("-g_rscript=", argv[i]))
			{
				if (g_rscript_b) Exit("\nError: parameter %s has been assigned twice.\n", argv[i].c_str());
				g_rscript_b = true;
				g_rscript_val = TrimParQuote(argv[i]);
			}
			else if (!LwrLineCmp("-g_benchmark=", argv[i]))
				GetParString(argv[i], "yes|no", g_benchmark_b, g_benchmark_val);
			else if (!LwrLineCmp("-g_replot=", argv[i]))
				GetParString(argv[i], "yes|no", g_replot_b, g_replot_val);
			else if (!LwrLineCmp("-g_eval=", argv[i]))
				GetParString(argv[i], "yes|no", g_eval_b, g_eval_val);
			else if (!LwrLineCmp("-g_missingploidy=", argv[i]))
				GetParStringMultiSel(argv[i], "1|2|3|4|5|6|7|8|9|10", g_missingploidy_b, g_missingploidy_val);
			else if (!LwrLineCmp("-g_maxlen=", argv[i]))
				GetParLong(argv[i], g_maxlen_b, g_maxlen_val, 0, 99999999999999999);
			else
				Exit("\nError: Unrecognized parameter: %s\n", argv[i].c_str());
		}
		else if (!LwrLineCmp("-haplotype", argv[i]))
		{
			if (!LwrStrCmp("-haplotype", argv[i]))
				GetParBool(argv[i], haplotype);
			else if (!LwrLineCmp("-haplotype_ptype=", argv[i]))
				GetRangeParDouble(argv[i], haplotype_ptype_b, haplotype_ptype_min, haplotype_ptype_max, 0.1, 1);
			else if (!LwrLineCmp("-haplotype_length=", argv[i]))
				GetRangeParLong(argv[i], haplotype_length_b, haplotype_length_min, haplotype_length_max, 1ll, 10000000000ll);
			else if (!LwrLineCmp("-haplotype_variants=", argv[i]))
				GetRangeParInteger(argv[i], haplotype_variants_b, haplotype_variants_min, haplotype_variants_max, 1, 1000000000);
			else if (!LwrLineCmp("-haplotype_interval=", argv[i]))
				GetParLong(argv[i], haplotype_interval_b, haplotype_interval_val, 0ll, 1000000000ll);
			else if (!LwrLineCmp("-haplotype_alleles=", argv[i]))
				GetRangeParInteger(argv[i], haplotype_alleles_b, haplotype_alleles_min, haplotype_alleles_max, 2, 65535);
			else if (!LwrLineCmp("-haplotype_genotypes=", argv[i]))
				GetRangeParInteger(argv[i], haplotype_genotypes_b, haplotype_genotypes_min, haplotype_genotypes_max, 2, 65535);
			else
				Exit("\nError: Unrecognized parameter: %s\n", argv[i].c_str());
		}
		else if (!LwrLineCmp("-slide", argv[i]))
		{
			if (!LwrStrCmp("-slide", argv[i]))
				GetParBool(argv[i], slide);
			else if (!LwrLineCmp("-slide_plot=", argv[i]))
				GetParString(argv[i], "yes|no", slide_plot_b, slide_plot_val);
			else if (!LwrLineCmp("-slide_plot_columns=", argv[i]))
				GetParIntegerArray(argv[i], slide_plot_columns_b, slide_plot_columns_val, 1, 100, N_MAX_SLIDEPLOT);
			else if (!LwrLineCmp("-slide_plot_styles=", argv[i]))
				GetParStringArray(argv[i], "dot|bar|line|heat", slide_plot_styles_b, slide_plot_styles_val, N_MAX_SLIDEPLOT);
			else if (!LwrLineCmp("-slide_windowsize=", argv[i]))
				GetParInteger(argv[i], slide_windowsize_b, slide_windowsize_val, 1000, 1000000000);
			else if (!LwrLineCmp("-slide_windowstep=", argv[i]))
				GetParInteger(argv[i], slide_windowstep_b, slide_windowstep_val, 1000, 1000000000);
			else if (!LwrLineCmp("-slide_minvariants=", argv[i]))
				GetParInteger(argv[i], slide_minvariants_b, slide_minvariants_val, 10, 1000000000);
			else if (!LwrLineCmp("-slide_estimator=", argv[i]))
				GetParStringMultiSel(argv[i], "Nei1973|Weir1984|Hudson1992|Hedrick2005|Jost2008|Huang2021_aneu|dxy|pi|thetaw|TajimaD|r2|D'|r2D|Delta'|fis|ho|he|pic|ae|I", slide_estimator_b, slide_estimator_val);
			else if (!LwrLineCmp("-slide_pop=", argv[i]))
			{
				if (slide_pop_b) Exit("\nError: parameter %s has been assigned twice.\n", argv[i].c_str());
				slide_pop_b = true;
				slide_pop_val = TrimParQuote(argv[i]);
			}
			else
				Exit("\nError: Unrecognized parameter: %s\n", argv[i].c_str());
		}
		else if (!LwrLineCmp("-convert", argv[i]))
		{
			if (!LwrStrCmp("-convert", argv[i]))
				GetParBool(argv[i], convert);
			else if (!LwrLineCmp("-convert_format=", argv[i]))
				GetParStringMultiSel(argv[i], "genepop|spagedi|cervus|arlequin|structure|polygene|polyrelatedness|genodive|plink", convert_format_b, convert_format_val);
			else if (!LwrLineCmp("-convert_mode=", argv[i]))
				GetParString(argv[i], "disable|truncate|choose|split|shuffle", convert_mode_b, convert_mode_val);
			else
				Exit("\nError: Unrecognized parameter: %s\n", argv[i].c_str());
		}
		else if (!LwrLineCmp("-indstat", argv[i]))
		{
			if (!LwrStrCmp("-indstat", argv[i]))
				GetParBool(argv[i], indstat);
			else if (!LwrLineCmp("-indstat_type=", argv[i]))
				GetParStringMultiSel(argv[i], "hidx|lnpg|f|theta", indstat_type_b, indstat_type_val);
			else if (!LwrLineCmp("-indstat_model=", argv[i]))
				GetParStringMultiSel(argv[i], "rcs|prcs|ces|pes", indstat_model_b, indstat_model_val);
			else if (!LwrLineCmp("-indstat_estimator=", argv[i]))
				GetParStringMultiSel(argv[i], "Ritland1996|Loiselle1995|Weir1996", indstat_estimator_b, indstat_estimator_val);
			else if (!LwrLineCmp("-indstat_ref=", argv[i]))
				GetParStringMultiSel(argv[i], "pop|reg|total", indstat_ref_b, indstat_ref_val);
			else if (!LwrLineCmp("-indstat_locus=", argv[i]))
				GetParStringMultiSel(argv[i], "all|each", indstat_locus_b, indstat_locus_val);
			else
				Exit("\nError: Unrecognized parameter: %s\n", argv[i].c_str());
		}
		else if (!LwrLineCmp("-diversity", argv[i]))
		{
			if (!LwrStrCmp("-diversity", argv[i]))
				GetParBool(argv[i], diversity);
			else if (!LwrLineCmp("-diversity_level=", argv[i]))
				GetParStringMultiSel(argv[i], "pop|reg|tot|popXloc|regXloc|totXloc", diversity_level_b, diversity_level_val);
			else if (!LwrLineCmp("-diversity_model=", argv[i]))
				GetParStringMultiSel(argv[i], "rcs|prcs|ces|pes", diversity_model_b, diversity_model_val);
			else
				Exit("\nError: Unrecognized parameter: %s\n", argv[i].c_str());
		}
		else if (!LwrLineCmp("-fst", argv[i]))
		{
			if (!LwrStrCmp("-fst", argv[i]))
				GetParBool(argv[i], fst);
			else if (!LwrLineCmp("-fst_plot=", argv[i]))
				GetParString(argv[i], "yes|no", fst_plot_b, fst_plot_val);
			else if (!LwrLineCmp("-fst_level=", argv[i]))
				GetParStringMultiSel(argv[i], "regXtot|popXtot|popXreg|reg|pop", fst_level_b, fst_level_val);
			else if (!LwrLineCmp("-fst_estimator=", argv[i]))
				GetParStringMultiSel(argv[i], "Nei1973|Weir1984|Hudson1992|Slatkin1995|Hedrick2005|Jost2008|Huang2021_homo|Huang2021_aneu", fst_estimator_b, fst_estimator_val);
			else if (!LwrLineCmp("-fst_fmt=", argv[i]))
				GetParStringMultiSel(argv[i], "matrix|table", fst_fmt_b, fst_fmt_val);
			else if (!LwrLineCmp("-fst_locus=", argv[i]))
				GetParStringMultiSel(argv[i], "all|each", fst_locus_b, fst_locus_val);
			else if (!LwrLineCmp("-fst_test=", argv[i]))
				GetParStringMultiSel(argv[i], "genotype|allele", fst_test_b, fst_test_val);
			else
				Exit("\nError: Unrecognized parameter: %s\n", argv[i].c_str());
		}
		else if (!LwrLineCmp("-gdist", argv[i]))
		{
			if (!LwrStrCmp("-gdist", argv[i]))
				GetParBool(argv[i], gdist);
			else if (!LwrLineCmp("-gdist_plot=", argv[i]))
				GetParString(argv[i], "yes|no", gdist_plot_b, gdist_plot_val);
			else if (!LwrLineCmp("-gdist_level=", argv[i]))
				GetParStringMultiSel(argv[i], "ind|pop|reg", gdist_level_b, gdist_level_val);
			else if (!LwrLineCmp("-gdist_weightmissing=", argv[i]))
				GetParString(argv[i], "yes|no", gdist_weightmissing_b, gdist_weightmissing_val);
			else if (!LwrLineCmp("-gdist_estimator=", argv[i]))
				GetParStringMultiSel(argv[i], "Nei1972|Cavalli-Sforza1967|Reynolds1983|Nei1983|Euclidean|Goldstein1995|Nei1974|Roger1972|Slatkin_Nei1973|Slatkin_Weir1984|Slatkin_Hudson1992|Slatkin_Slatkin1995|Slatkin_Hedrick2005|Slatkin_Jost2008|Slatkin_Huang2021_homo|Slatkin_Huang2021_aneu|Reynolds_Nei1973|Reynolds_Weir1984|Reynolds_Hudson1992|Reynolds_Slatkin1995|Reynolds_Hedrick2005|Reynolds_Jost2008|Reynolds_Huang2021_homo|Reynolds_Huang2021_aneu", gdist_estimator_b, gdist_estimator_val);
			else if (!LwrLineCmp("-gdist_fmt=", argv[i]))
				GetParStringMultiSel(argv[i], "matrix|table", gdist_fmt_b, gdist_fmt_val);
			else
				Exit("\nError: Unrecognized parameter: %s\n", argv[i].c_str());
		}
		else if (!LwrLineCmp("-amova", argv[i]))
		{
			if (!LwrStrCmp("-amova", argv[i]))
				GetParBool(argv[i], amova);
			else if (!LwrLineCmp("-amova_method=", argv[i]))
				GetParStringMultiSel(argv[i], "homoploid|aneuploid|likelihood", amova_method_b, amova_method_val);
			else if (!LwrLineCmp("-amova_mutation=", argv[i]))
				GetParStringMultiSel(argv[i], "iam|smm", amova_mutation_b, amova_mutation_val);
			else if (!LwrLineCmp("-amova_ind=", argv[i]))
				GetParStringMultiSel(argv[i], "yes|no", amova_ind_b, amova_ind_val);
			else if (!LwrLineCmp("-amova_trunc=", argv[i]))
				GetParString(argv[i], "yes|no", amova_trunc_b, amova_trunc_val);
			else if (!LwrLineCmp("-amova_test=", argv[i]))
				GetParString(argv[i], "yes|no", amova_test_b, amova_test_val);
			else if (!LwrLineCmp("-amova_nperm=", argv[i]))
				GetParInteger(argv[i], amova_nperm_b, amova_nperm_val, 1, 99999999);
			else if (!LwrLineCmp("-amova_pseudo=", argv[i]))
				GetParInteger(argv[i], amova_pseudo_b, amova_pseudo_val, 0, 9999);
			else if (!LwrLineCmp("-amova_printss=", argv[i]))
				GetParString(argv[i], "yes|no", amova_printss_b, amova_printss_val);
			else
				Exit("\nError: Unrecognized parameter: %s\n", argv[i].c_str());
		}
		else if (!LwrLineCmp("-popas", argv[i]))
		{
			if (!LwrStrCmp("-popas", argv[i]))
				GetParBool(argv[i], popas);
			else if (!LwrLineCmp("-popas_plot=", argv[i]))
				GetParString(argv[i], "yes|no", popas_plot_b, popas_plot_val);
			else if (!LwrLineCmp("-popas_model=", argv[i]))
				GetParStringMultiSel(argv[i], "rcs|prcs|ces|pes", popas_model_b, popas_model_val);
			else if (!LwrLineCmp("-popas_level=", argv[i]))
				GetParStringMultiSel(argv[i], "pop|reg", popas_level_b, popas_level_val);
			else if (!LwrLineCmp("-popas_error=", argv[i]))
				GetParDouble(argv[i], popas_error_b, popas_error_val, 0, 0.2);
			else
				Exit("\nError: Unrecognized parameter: %s\n", argv[i].c_str());
		}
		else if (!LwrLineCmp("-relatedness", argv[i]))
		{
			if (!LwrStrCmp("-relatedness", argv[i]))
				GetParBool(argv[i], relatedness);
			else if (!LwrLineCmp("-relatedness_plot=", argv[i]))
				GetParString(argv[i], "yes|no", relatedness_plot_b, relatedness_plot_val);
			else if (!LwrLineCmp("-relatedness_range=", argv[i]))
				GetParStringMultiSel(argv[i], "pop|reg|total", relatedness_range_b, relatedness_range_val);
			else if (!LwrLineCmp("-relatedness_fmt=", argv[i]))
				GetParStringMultiSel(argv[i], "matrix|table", relatedness_fmt_b, relatedness_fmt_val);
			else if (!LwrLineCmp("-relatedness_estimator=", argv[i]))
				GetParStringMultiSel(argv[i], "Lynch1999|Wang2002|Thomas2010|Li1993|Queller1989|Huang2016A|Huang2016B|Milligan2003|Anderson2007|Huang2014|Huang2015|Ritland1996_modified|Loiselle1995_modified|Ritland1996|Loiselle1995|Weir1996", relatedness_estimator_b, relatedness_estimator_val);
			else
				Exit("\nError: Unrecognized parameter: %s\n", argv[i].c_str());
		}
		else if (!LwrLineCmp("-kinship", argv[i]))
		{
			if (!LwrStrCmp("-kinship", argv[i]))
				GetParBool(argv[i], kinship);
			else if (!LwrLineCmp("-kinship_plot=", argv[i]))
				GetParString(argv[i], "yes|no", kinship_plot_b, kinship_plot_val);
			else if (!LwrLineCmp("-kinship_range=", argv[i]))
				GetParStringMultiSel(argv[i], "pop|reg|total", kinship_range_b, kinship_range_val);
			else if (!LwrLineCmp("-kinship_fmt=", argv[i]))
				GetParStringMultiSel(argv[i], "matrix|table", kinship_fmt_b, kinship_fmt_val);
			else if (!LwrLineCmp("-kinship_estimator=", argv[i]))
				GetParStringMultiSel(argv[i], "Ritland1996|Loiselle1995|Weir1996", kinship_estimator_b, kinship_estimator_val);
			else
				Exit("\nError: Unrecognized parameter: %s\n", argv[i].c_str());
		}
		else if (!LwrLineCmp("-pcoa", argv[i]))
		{
			if (!LwrStrCmp("-pcoa", argv[i]))
				GetParBool(argv[i], pcoa);
			else if (!LwrLineCmp("-pcoa_plot=", argv[i]))
				GetParString(argv[i], "yes|no", pcoa_plot_b, pcoa_plot_val);
			else if (!LwrLineCmp("-pcoa_level=", argv[i]))
				GetParStringMultiSel(argv[i], "ind|pop|reg", pcoa_level_b, pcoa_level_val);
			else if (!LwrLineCmp("-pcoa_dim=", argv[i]))
				GetParInteger(argv[i], pcoa_dim_b, pcoa_dim_val, 1, 4096);
			else if (!LwrLineCmp("-pcoa_estimator=", argv[i]))
				GetParStringMultiSel(argv[i], "Nei1972|Cavalli-Sforza1967|Reynolds1983|Nei1983|Euclidean|Goldstein1995|Nei1974|Roger1972|Slatkin_Nei1973|Slatkin_Weir1984|Slatkin_Hudson1992|Slatkin_Slatkin1995|Slatkin_Hedrick2005|Slatkin_Jost2008|Slatkin_Huang2021_homo|Slatkin_Huang2021_aneu|Reynolds_Nei1973|Reynolds_Weir1984|Reynolds_Hudson1992|Reynolds_Slatkin1995|Reynolds_Hedrick2005|Reynolds_Jost2008|Reynolds_Huang2021_homo|Reynolds_Huang2021_aneu", pcoa_estimator_b, pcoa_estimator_val);
			else
				Exit("\nError: Unrecognized parameter: %s\n", argv[i].c_str());
		}
		else if (!LwrLineCmp("-cluster", argv[i]))
		{
			if (!LwrStrCmp("-cluster", argv[i]))
				GetParBool(argv[i], cluster);
			else if (!LwrLineCmp("-cluster_plot=", argv[i]))
				GetParString(argv[i], "yes|no", cluster_plot_b, cluster_plot_val);
			else if (!LwrLineCmp("-cluster_level=", argv[i]))
				GetParStringMultiSel(argv[i], "ind|pop|reg", cluster_level_b, cluster_level_val);
			else if (!LwrLineCmp("-cluster_method=", argv[i]))
				GetParStringMultiSel(argv[i], "NEAREST|FURTHEST|UPGMA|WPGMA|UPGMC|WPGMC|WARD", cluster_method_b, cluster_method_val);
			else if (!LwrLineCmp("-cluster_estimator=", argv[i]))
				GetParStringMultiSel(argv[i], "Nei1972|Cavalli-Sforza1967|Reynolds1983|Nei1983|Euclidean|Goldstein1995|Nei1974|Roger1972|Slatkin_Nei1973|Slatkin_Weir1984|Slatkin_Hudson1992|Slatkin_Slatkin1995|Slatkin_Hedrick2005|Slatkin_Jost2008|Slatkin_Huang2021_homo|Slatkin_Huang2021_aneu|Reynolds_Nei1973|Reynolds_Weir1984|Reynolds_Hudson1992|Reynolds_Slatkin1995|Reynolds_Hedrick2005|Reynolds_Jost2008|Reynolds_Huang2021_homo|Reynolds_Huang2021_aneu", cluster_estimator_b, cluster_estimator_val);
			else
				Exit("\nError: Unrecognized parameter: %s\n", argv[i].c_str());
		}
		else if (!LwrLineCmp("-structure", argv[i]))
		{
			if (!LwrStrCmp("-structure", argv[i]))
				GetParBool(argv[i], structure);
			else if (!LwrLineCmp("-structure_nstream=", argv[i]))
				GetParInteger(argv[i], structure_nstream_b, structure_nstream_val, 1, 32);
			else if (!LwrLineCmp("-structure_plot=", argv[i]))
				GetParString(argv[i], "yes|no", structure_plot_b, structure_plot_val);
			else if (!LwrLineCmp("-structure_writelnl=", argv[i]))
				GetParString(argv[i], "yes|no", structure_writelnl_b, structure_writelnl_val);
			else if (!LwrLineCmp("-structure_eval=", argv[i]))
				GetParString(argv[i], "yes|no", structure_eval_b, structure_eval_val);

			else if (!LwrLineCmp("-structure_admix=", argv[i]))
				GetParString(argv[i], "yes|no", structure_admix_b, structure_admix_val);
			else if (!LwrLineCmp("-structure_locpriori=", argv[i]))
				GetParString(argv[i], "yes|no", structure_locpriori_b, structure_locpriori_val);
			else if (!LwrLineCmp("-structure_f=", argv[i]))
				GetParString(argv[i], "yes|no", structure_f_b, structure_f_val);

			else if (!LwrLineCmp("-structure_krange=", argv[i]))
				GetRangeParInteger(argv[i], structure_krange_b, structure_krange_min, structure_krange_max, 1, 1000);
			else if (!LwrLineCmp("-structure_nburnin=", argv[i]))
				GetParInteger(argv[i], structure_nburnin_b, structure_nburnin_val, 0, 10000000);
			else if (!LwrLineCmp("-structure_nreps=", argv[i]))
				GetParInteger(argv[i], structure_nreps_b, structure_nreps_val, 0, 10000000);
			else if (!LwrLineCmp("-structure_nthinning=", argv[i]))
				GetParInteger(argv[i], structure_nthinning_b, structure_nthinning_val, 1, 10000);
			else if (!LwrLineCmp("-structure_nruns=", argv[i]))
				GetParInteger(argv[i], structure_nruns_b, structure_nruns_val, 1, 1000);

			else if (!LwrLineCmp("-structure_nadmburnin=", argv[i]))
				GetParInteger(argv[i], structure_nadmburnin_b, structure_nadmburnin_val, 0, 10000000);
			else if (!LwrLineCmp("-structure_lambda=", argv[i]))
				GetParDouble(argv[i], structure_lambda_b, structure_lambda_val, 0, 10000);
			else if (!LwrLineCmp("-structure_stdlambda=", argv[i]))
				GetParDouble(argv[i], structure_stdlambda_b, structure_stdlambda_val, 0, 10000);
			else if (!LwrLineCmp("-structure_maxlambda=", argv[i]))
				GetParDouble(argv[i], structure_maxlambda_b, structure_maxlambda_val, 0, 10000);
			else if (!LwrLineCmp("-structure_inferlambda=", argv[i]))
				GetParString(argv[i], "yes|no", structure_inferlambda_b, structure_inferlambda_val);
			else if (!LwrLineCmp("-structure_difflambda=", argv[i]))
				GetParString(argv[i], "yes|no", structure_difflambda_b, structure_difflambda_val);
			else if (!LwrLineCmp("-structure_diversity=", argv[i]))
				GetParString(argv[i], "yes|no", structure_diversity_b, structure_diversity_val);

			else if (!LwrLineCmp("-structure_alpha=", argv[i]))
				GetParDouble(argv[i], structure_alpha_b, structure_alpha_val, 0, 10000);
			else if (!LwrLineCmp("-structure_inferalpha=", argv[i]))
				GetParString(argv[i], "yes|no", structure_inferalpha_b, structure_inferalpha_val);
			else if (!LwrLineCmp("-structure_diffalpha=", argv[i]))
				GetParString(argv[i], "yes|no", structure_diffalpha_b, structure_diffalpha_val);
			else if (!LwrLineCmp("-structure_uniformalpha=", argv[i]))
				GetParString(argv[i], "yes|no", structure_uniformalpha_b, structure_uniformalpha_val);
			else if (!LwrLineCmp("-structure_stdalpha=", argv[i]))
				GetParDouble(argv[i], structure_stdalpha_b, structure_stdalpha_val, 0, 10000);
			else if (!LwrLineCmp("-structure_maxalpha=", argv[i]))
				GetParDouble(argv[i], structure_maxalpha_b, structure_maxalpha_val, 0, 10000);
			else if (!LwrLineCmp("-structure_alphapriora=", argv[i]))
				GetParDouble(argv[i], structure_alphapriora_b, structure_alphapriora_val, 0, 10000);
			else if (!LwrLineCmp("-structure_alphapriorb=", argv[i]))
				GetParDouble(argv[i], structure_alphapriorb_b, structure_alphapriorb_val, 0, 10000);
			else if (!LwrLineCmp("-structure_metrofreq=", argv[i]))
				GetParInteger(argv[i], structure_metrofreq_b, structure_metrofreq_val, 0, 1000000);

			else if (!LwrLineCmp("-structure_r=", argv[i]))
				GetParDouble(argv[i], structure_r_b, structure_r_val, 0, 10000);
			else if (!LwrLineCmp("-structure_maxr=", argv[i]))
				GetParDouble(argv[i], structure_maxr_b, structure_maxr_val, 0, 10000);
			else if (!LwrLineCmp("-structure_stdr=", argv[i]))
				GetParDouble(argv[i], structure_stdr_b, structure_stdr_val, 0, 10000);
			else if (!LwrLineCmp("-structure_epseta=", argv[i]))
				GetParDouble(argv[i], structure_epseta_b, structure_epseta_val, 0, 10000);
			else if (!LwrLineCmp("-structure_epsgamma=", argv[i]))
				GetParDouble(argv[i], structure_epsgamma_b, structure_epsgamma_val, 0, 10000);

			else if (!LwrLineCmp("-structure_pmeanf=", argv[i]))
				GetParDouble(argv[i], structure_pmeanf_b, structure_pmeanf_val, 0, 10000);
			else if (!LwrLineCmp("-structure_pstdf=", argv[i]))
				GetParDouble(argv[i], structure_pstdf_b, structure_pstdf_val, 0, 10000);
			else if (!LwrLineCmp("-structure_stdf=", argv[i]))
				GetParDouble(argv[i], structure_stdf_b, structure_stdf_val, 0, 10000);
			else if (!LwrLineCmp("-structure_singlef=", argv[i]))
				GetParString(argv[i], "yes|no", structure_singlef_b, structure_singlef_val);
			else
				Exit("\nError: Unrecognized parameter: %s\n", argv[i].c_str());
		}
		else if (!LwrLineCmp("-ploidyinfer", argv[i]))
		{
			if (!LwrStrCmp("-ploidyinfer", argv[i]))
				GetParBool(argv[i], ploidyinfer);
			else if (!LwrLineCmp("-ploidyinfer_type=", argv[i]))
				GetParStringMultiSel(argv[i], "1|2|3|4|5|6|7|8|9|10", ploidyinfer_type_b, ploidyinfer_type_val);
			else if (!LwrLineCmp("-ploidyinfer_histogram=", argv[i]))
				GetParString(argv[i], "yes|no", ploidyinfer_histogram_b, ploidyinfer_histogram_val);
			else if (!LwrLineCmp("-ploidyinfer_nbins=", argv[i]))
				GetParInteger(argv[i], ploidyinfer_nbins_b, ploidyinfer_nbins_val, 10, 100);
			else
				Exit("\nError: Unrecognized parameter: %s\n", argv[i].c_str());
		}
		else continue;

		if (LwrLineCmp("-g_indtab=", argv[i]))
			printf("%s\r\n", argv[i].c_str());
	}

	if (f_filter && (f_bmaf_b || f_k_b || f_n_b || f_ptype_b || f_pval_b || f_he_b || f_ho_b || f_pic_b || f_ae_b || f_I_b || f_windowsize_b))
		diversity_filter = true;
	if (f_filter && (f_dp_b || f_gq_b || f_ploidy_b))
		genotype_filter = true;
	if (f_filter && (f_itype_b || f_iploidy_b))
		individual_filter = true;
	if (f_filter && (f_qual_b || f_type_b || f_original_b || f_chrprefix_b || f_chrname_b))
		info_filter = true;

	if (!isparfile && p_b)
	{
		string partext = ReadAllText(p_val);
		partext = ReplaceStr(partext, "\r\n", "\n");
		partext = ReplaceStr(partext, "\n", "#DELIM#", '"');
		partext = ReplaceStr(partext, "\t", "#DELIM#", '"');
		partext = ReplaceStr(partext, " ", "#DELIM#", '"');
		partext = ReplaceStr(partext, "#DELIM##DELIM##DELIM##DELIM##DELIM#", "#DELIM#");
		partext = ReplaceStr(partext, "#DELIM##DELIM##DELIM##DELIM##DELIM#", "#DELIM#");
		partext = ReplaceStr(partext, "#DELIM##DELIM##DELIM##DELIM##DELIM#", "#DELIM#");
		partext = ReplaceStr(partext, "#DELIM##DELIM##DELIM##DELIM##DELIM#", "#DELIM#");
		partext = ReplaceStr(partext, "#DELIM##DELIM#", "#DELIM#");
		partext = ReplaceStr(partext, "#DELIM##DELIM#", "#DELIM#");
		argv = SplitStr(partext, "#DELIM#");

		for (int i = 0; i < argv.size(); ++i)
			argv[i] = TrimQuote(argv[i]);
		SetParameters(true);
		return;
	}

	if ((	structure_plot_val == 1 || 
			cluster_plot_val == 1 ||
			pcoa_plot_val == 1 ||
			kinship_plot_val == 1 ||
			relatedness_plot_val == 1 ||
			gdist_plot_val == 1 ||
			fst_plot_val == 1 ||
			popas_plot_val == 1 ||
			g_replot_val == 1) &&
		!FileExists(g_rscript_val.c_str()))
		Exit("\nError: Rscript binary executable is not found.\n");
}

/* Release parameter memory */
TARGET void ReleaseParameters()
{
	vector<vector<FILEINFO>>().swap(FILE_INFO);
	vector<string>().swap(argv);
	string().swap(p_val);
	string().swap(g_indtext_val);
	string().swap(g_indtab_val);
	string().swap(g_input_val);
	string().swap(g_rscript_val);
	string().swap(g_output_val);
	string().swap(g_tmpdir_val);
	string().swap(f_pop_val);
	string().swap(slide_pop_val);
	string().swap(spa_coord_val);
	string().swap(OUTDIR);
	string().swap(EXEDIR);
	string().swap(CURDIR);
	string().swap(OUTFILE);
    UnInitLock(GLOCK1);
    UnInitLock(GLOCK2);
}

/* Write parameters to result file */
TARGET void WriteParameters(FILE* f1, char* type, char* prefix)
{
	fprintf(f1, "%s  Parameter: %s%s", prefix, argv[0].c_str(), g_linebreak_val);
	for (int i = 1; i < argv.size(); ++i)
	{
		if (f_filter &&!LwrLineCmp("-f_", argv[i]) || 
			(LwrLineCmp("-g_", argv[i]) == 0 && LwrLineCmp("-g_indtext", argv[i]) != 0 && LwrLineCmp("-g_indtab", argv[i]) != 0) ||
			LwrLineCmp(type, argv[i]) == 0)
			fprintf(f1, "%s    %s%s", prefix, argv[i].c_str(), g_linebreak_val);
	}
}

/* Convert indtab into indtext*/
TARGET void ConvertIndTable()
{
	g_indtab_val = ReplaceStr(g_indtab_val, "\r\n", "\n");
	g_indtab_val = ReplaceStr(g_indtab_val, "\t", " ");

	vector<string> lines = SplitStr(g_indtab_val, "\n");
	int nlines = (int)lines.size(), ncols0 = 0;

	indtab_node* root = NULL;
	vector<indtab_node*> nodes;
	for (int i = 0; i < nlines; ++i)
	{
		vector<string> grids = SplitStr(lines[i], " ");

		if (ncols0 == 0)
		{
			ncols0 = (int)grids.size();
			root = new indtab_node("total", ncols0);
		}

		if (ncols0 != grids.size())
			Exit("\nError: -g_indtab have different number of columns in lines %d.\n", i + 1);

		indtab_node* pnode = root;
		for (int j = (int)grids.size() - 1; j >= 0; --j)
		{
			if (pnode->subunits.find(grids[j]) == pnode->subunits.end())
			{
				indtab_node* tnode = new indtab_node(grids[j], j);
				nodes.push_back(tnode);
				pnode->subunits[grids[j]] = tnode;
			}
			pnode = pnode->subunits[grids[j]];
		}
	}

	//check nodes with the same name
	for (int i = 0; i < nodes.size(); ++i)
		for (int j = i + 1; j < nodes.size(); ++j)
			if (nodes[i]->name == nodes[j]->name)
				Exit("\nError: two objects (individuals, populations or regions) have the same identifier: %s, %s",
					nodes[i]->name.c_str(), nodes[j]->name.c_str());

	//prepare -g_indtext
	string* buf = new string[ncols0 + 1];
	indtab_node* pnode = root;
	while (pnode->subunits.size() == 1 && pnode->level > 1)
		pnode = pnode->subunits.begin()->second;

	for (int i = 2; i <= pnode->level; ++i)
		buf[i].append(" #REG");

	pnode->Write(buf);

	//write -g_indtext
	g_indtext_val = "";
	for (int i = 1; i <= pnode->level; ++i)
		g_indtext_val.append(buf[i]);
	delete[] buf;

	for (int i = 0; i < nodes.size(); ++i)
		delete nodes[i];

	if (root) delete root;
}
