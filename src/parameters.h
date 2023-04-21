/* Parameter Functions */

#pragma once
#include "vcfpop.h"

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
extern bool g_fastsinglet_b;					extern int g_fastsingle_val;
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

/* Filters */
extern bool f_filter;
extern bool f_qual_b;							extern double f_qual_min, f_qual_max;
extern bool f_type_b;							extern int f_type_val;
extern bool f_original_b;						extern int f_original_val;
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
																								  

/* Initialize parameters */
TARGET void SetDefaultParameters();

/* Replace path deliminator*/
TARGET string ReplacePathDelim(string& str, bool ispath);

/* Set parameters */
TARGET void SetParameters(bool isparfile);

/* Release parameter memory */
TARGET void ReleaseParameters();

/* Write parameters to result file */
TARGET void WriteParameters(FILE* f1, char* type, char* prefix);

/* Convert indtab into indtext*/
TARGET void ConvertIndTable();
