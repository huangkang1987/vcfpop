/* Menu */

#include "vcfpop.h"

/* Print help information */
TARGET bool PrintHelp()
{
	bool helpcat[N_FUNC] = { 0 };
	uint helpcount = 0;
	bool gpar = false;
	for (uint i = 0; i < argv.size(); ++i)
	{
		if (!LwrLineCmp("-help", argv[i]) || !LwrLineCmp("-h", argv[i]) ||
			!LwrLineCmp("--help", argv[i]) || !LwrLineCmp("--h", argv[i]) ||
			!LwrLineCmp("/help", argv[i]) || !LwrLineCmp("/h", argv[i]))
			helpcat[0] = true;
		if (!LwrLineCmp("-p=", argv[i]))
			gpar = true;

		if ("-g" == argv[i] || !LwrLineCmp("-g_", argv[i]))
			helpcat[1] = true;
		if ("-f" == argv[i] || !LwrLineCmp("-f_", argv[i]))
			helpcat[2] = true;
		if (!LwrLineCmp("-haplotype", argv[i]))
			helpcat[3] = true;
		if (!LwrLineCmp("-convert", argv[i]))
			helpcat[4] = true;
		if (!LwrLineCmp("-diversity", argv[i]))
			helpcat[5] = true;
		if (!LwrLineCmp("-indstat", argv[i]))
			helpcat[6] = true;
		if (!LwrLineCmp("-fst", argv[i]))
			helpcat[7] = true;
		if (!LwrLineCmp("-gdist", argv[i]))
			helpcat[8] = true;
		if (!LwrLineCmp("-amova", argv[i]))
			helpcat[9] = true;
		if (!LwrLineCmp("-popas", argv[i]))
			helpcat[10] = true;
		if (!LwrLineCmp("-relatedness", argv[i]))
			helpcat[11] = true;
		if (!LwrLineCmp("-kinship", argv[i]))
			helpcat[12] = true;
		if (!LwrLineCmp("-pcoa", argv[i]))
			helpcat[13] = true;
		if (!LwrLineCmp("-cluster", argv[i]))
			helpcat[14] = true;
		if (!LwrLineCmp("-structure", argv[i]))
			helpcat[15] = true;
		if (!LwrLineCmp("-ploidyinfer", argv[i]))
			helpcat[16] = true;
	}

	for (uint i = 0; i < N_FUNC; ++i)
		if (helpcat[i])
			helpcount++;

	if ((!helpcat[0] && helpcount) || gpar)
		return false;

	//0. Help
	if (helpcat[0] || !helpcount)
	{
		printf("0. Help and Parameters set\n");
		printf("-h, -help, --h, --help, /h, /help\n");
		printi("Print help information for corresponding commands.\n");
		printi("Basic commands, type -h -function_name to view detail help information:\n");
		printi("  1. General settings: -g\n");
		printi("  2. Filter for individual, locus or genotype: -f\n");
		printi("  3. Haplotype extraction: -haplotype\n");
		printi("  4. File conversion: -convert\n");
		printi("  5. Genetic diversity indices: -diversity\n");
		printi("  6. Individual statistics: -indstat\n");
		printi("  7. Genetic differentiation: -fst\n");
		printi("  8. Genetic distance: -gdist\n");
		printi("  9. Analysis of molecular variance: -amova\n");
		printi("  10. Population assignment: -popas\n");
		printi("  11. Relatedness coefficient: -relatedness\n");
		printi("  12. Kinship coefficient: -kinship\n");
		printi("  13. Principal coordinate analysis: -pcoa\n");
		printi("  14. Hierarchical clustering: -cluster\n");
		printi("  15. Bayesian clustering: -structure\n");
		//printi("  16. Ploidy inference: -ploidyinfer\n");
		printf("  -p=pars.txt\n");
		printi("Load parameters set file, the formats are the same as command-line parameters, and linebreak can be used as the parameter sepreator in the parameters set file.\n");
		printf("\n");
	}

	//2. General settings
	if (helpcat[1])
	{
		printf("1. General settings\n");
		printf("-g_decimal=0~15, integer, default:5\n");
		printi("Decimal places of output real numbers.\n");
		printf("-g_scientific=yes|no, string, default:no\n");
		printi("Use scientific notation to output real numbers.\n");
		printf("-g_nthread=1~4096, integer, default:4\n");
		printi("Number of threads used in calculation.\n");
#ifdef __aarch64__
		printf("-g_simd=none|neon, integer, default:neon\n");
		printi("SIMD (Single Instruction Multiple Data) instruction sets to accelerate float point or vector operations, where neon can handle 128 bits simultaneously, respectively.\n");
#else
		printf("-g_simd=none|sse|avx|avx512, integer, default:avx\n");
		printi("SIMD (Single Instruction Multiple Data) instruction sets to accelerate float point or vector operations, where mmx, sse (using SSE1.0 to SSE4.2 instructions), avx (using AVX, AVX2 and FMA instructions) and avx512 (using AVX512 F & BW instructions) can handle 64, 128, 256, 512 bits simultaneously, respectively.\n");
#endif
		printf("-g_seed=0~2147483647, integer, default:0\n");
		printi("Random number generator seed, 0 denotes using system time as the seed.\n");
		printf("-g_tmpdir=path, string, default:current_directory\n");
		printi("Directory placing the temporary files.\n");
		printf("-g_progress=10~100000, integer, default:80\n");
		printi("The number of characters used for the progress bar.\n");
		printf("-g_input=file_path, string\n");
		printi("Input file. Multiple VCF/BCF files using '|' and '&' as the column and row separators, respectively, e.g., var1-4ind1-3.vcf|var1-4ind4-6.vcf&var5-9ind1-3.vcf|var5-9ind4-6.vcf. The 'FORMAT' field in each file should be equal.\n");
		printf("-g_format=vcf|bcf|genepop|spagedi|cervus|arlequin|structure|polygene|polyrelatedness, string, default:vcf\n");
		printi("Input file format. The population and region should be defined in g_indtext, and the within input file will not be used. For GENEPOP format, VCFPOP do not support extra information; for CERVUS, at most one extra column (population or sex) is allowed; for SPAGEDI, multiple extra columns (population or coordinate) are allowed; for STRUCTURE, the number of extra columns can be specified.\n");
		printf("-g_locusname=chr|pos|chr_pos|chr_ref_alt|pos_ref_alt|chr_pos_ref_alt, string, default:chr_pos\n");
		printi("Locus name definition, only appliable for VCF/BCF input files.\n");
		printf("-g_extracol=0~4096, integer, default:0\n");
		printi("Regarding the number of extra columns between individuals and genotypes in the STRUCTURE input file, if there is a header row for locus, the number of extra columns can be automatically detected.\n");
		printf("-g_output=file_path, string, default:vcfpop.out\n");
		printi("Prefix of output files.\n");
		printf("-g_indtext=text, string, optional\n");
		printi("This assigns the population of each individual and region of each population (where regions can be nested to perform multi-level AMOVA). If not specified, all individuals are assigned to a default population. A population (or a region) name begins with an identifier and a colon, the individuals (or populations) of this population (or region) are separated by commas. The individuals, populations, and regions not included are assigned to a default population or a default region. '#REG' is used to separate the regions at different levels. In a command line, #n can be used as a line break, whilst in parameter files, line breaks can be used. For example:\n");
		
		printi("-g_indtext=\"pop1:ind1,ind2,ind3\n");
		printi("pop2:#4-#6\n");
		printi("pop3:ind7,ind8,ind9\n");
		printi("pop4:ind10,ind11,ind12\n");
		printi("#REG\n");
		printi("reg1:#1-#2\n");
		printi("reg2:pop3,pop4\n");
		printi("#REG\n");
		printi("Total:reg1,reg2\"\n");
		printf("-g_indtab=text, string, optional\n");
		printi("Tabular format of -g_indtext. In command line, #n can be used as line break, while in parameter files, line breaks can be used:\n");
		printi("-g_indtab=\"ind1	pop1	reg1	Total\n");
		printi("ind2	pop1	reg1	Total\n");
		printi("ind3	pop1	reg1	Total\n");
		printi("ind4	pop2	reg1	Total\n");
		printi("ind5	pop2	reg1	Total\n");
		printi("ind6	pop2	reg1	Total\n");
		printi("ind7	pop3	reg2	Total\n");
		printi("ind8	pop3	reg2	Total\n");
		printi("ind9	pop3	reg2	Total\n");
		printi("ind10	pop4	reg2	Total\n");
		printi("ind11	pop4	reg2	Total\n");
		printi("ind12	pop4	reg2	Total\"\n");
		printf("-g_delimitator=comma|tab, string, default:tab\n");
		printi("Column delimiter style, where a comma is used in the CSV format, and a tab is used in the text editors.\n");
		printf("-g_linebreak=unix|win, string, default:unix\n");
		printi("Line break style, where '\\n' is used for unix, and '\\r\\n' is used for windows.\n");
		printf("-g_rscript=file_path, string, optional\n");
		printi("Path of Rscript binary for plotting results.\n");
		printf("-g_benchmark=yes|no, string, default:no\n");
		printi("Evaluates the performance of each SIMD instructions set from mmx to specified type before calculation.\n");
		printf("-g_replot=yes|no, string, default:no\n");
		printi("Plots figures for previously calculated results instead of analyzing data.\n");
		printf("-g_eval=yes|no, string, default:no\n");
		printi("Evaluates the time expense for each function.\n");
		printf("-g_missingploidy=1|2|3|4|5|6|7|8|9|10, string, multiple selections, optional\n");
		printi("If the VCF/BCF file is combined from various VCF/BCF files for mixed-ploidy populations, the missing data may be extended from the called genotypes at this variant and the ploidy levels are misassigned. This option assigns the ploidy levels of the missing genotypes to the minimum ploidy level of the individual and is applicable for VCF/BCF input files. The value of this option should be all possible ploidy levels in the samples, e.g., -g_missingploidy=2,4,6.\n");
		printf("\n");
	}

	//2. Filter
	if (helpcat[2])
	{
		printf("2. Filter for variants, genotypes, individuals, and loci\n");

		printf("-f\n");
		printi("Enable filters.\n");

		printf("Variant information filters (applied during file load):\n");
		printf("-f_qual=[min_val,max_val], real range, optional\n");
		printi("Range of variant quality. If multiple VCF/BCF files are used, the variant is filtered when at least one QUAL field is out of the range. Only applicable to VCF/BCF input files.\n");
		printf("-f_type=snp|indel|both, string, optional\n");
		printi("Type of variants used in calculations. Only applicable to VCF/BCF input files.\n");
		printf("-f_original=yes|no, string, optional\n");
		printi("Use original filter of VCF/BCF file. If multiple VCF/BCF files are used, the variant is filtered when at least one original filter is not a 'PASS'. Only applicable to VCF/BCF input files.\n");

		printf("Genotype filters (applied during file load):\n");
		printf("-f_dp=[min_val,max_val], integer range, optional\n");
		printi("Range of sequencing depth. Only applicable to VCF/BCF input files.\n");
		printf("-f_gq=[min_val,max_val], integer range, optional\n");
		printi("Range of genotype quality. Only applicable to VCF/BCF input files.\n");
		printf("-f_ploidy=[min_val,max_val], integer range, optional\n");
		printi("Range of ploidy level for genotypes.\n");
		printf("Individual filters (applied after file load):\n");
		printf("-f_itype=[min_val,max_val], real range, optional\n");
		printi("Range of typed ratio of an individual.\n");
		printf("-f_iploidy=[min_val,max_val], integer range, optional\n");
		printi("Range of ploidy level for individuals.\n");

		printf("Locus diversity filters (applied after individual filter):\n");
		printf("-f_pop=pop_identifier|total, string, default:total\n");
		printi("Target population used to calculate diversity and apply diversity filters.\n");
		printf("-f_region=region_identifier, string, optional\n");
		printi("Target region used to calculate diversity and apply diversity filters.\n");
		printf("-f_bmaf=[min_val,max_val], real range, optional\n");
		printi("Range of frequencies of minor alleles for biallelic loci.\n");
		printf("-f_k=[min_val,max_val], integer range, optional\n");
		printi("Range of number of alleles.\n");
		printf("-f_n=[min_val,max_val], integer range, optional\n");
		printi("Range of number of typed individuals.\n");
		printf("-f_ptype=[min_val,max_val], real range, optional\n");
		printi("Range of typed ratio of a locus.\n");
		printf("-f_pval=[min_val,max_val], real range, optional\n");
		printi("Range of P values in equilibrium tests.\n");
		printf("-f_model=rcs|prcs|ces|pes, string, default:rcs\n");
		printi("Double-reduction model to calculate genotypic frequencies for polyploids.\n");
		printf("-f_he=[min_val,max_val], real range, optional\n");
		printi("Range of expected heterozygosity.\n");
		printf("-f_ho=[min_val,max_val], real range, optional\n");
		printi("Range of observed heterozygosity.\n");
		printf("-f_pic=[min_val,max_val], real range, optional\n");
		printi("Range of polymorphic information content.\n");
		printf("-f_ae=[min_val,max_val], real range, optional\n");
		printi("Range of effective number of alleles.\n");
		printf("-f_I=[min_val,max_val], real range, optional\n");
		printi("Range of Shannon's information index.\n");
		printf("-f_windowsize=1~1000000000, integer, optional\n");
		printi("Uses the most polymorphic variant in a window; only applicable to VCF/BCF input files.\n");
		printf("-f_windowstat=bmaf|k|n|ptype|he|ho|pic|ae|I, string, default:he\n");
		printi("The diversity statistic used to select the most polymorphic variant in a window; only applicable to VCF/BCF input files.\n");
		printf("\n");
	}

	//3. Haplotype
	if (helpcat[3])
	{
		printf("3. Haplotype extraction\n");
		printf("-haplotype\n");
		printi("Extracts haplotypes from phased genotypes, then use the haplotypes as alleles for further analysis. Note that all genotypes must be phased and only the variants genotyped in all individuals are used. The haplotype definitions are saved in *.haplotype.txt.\n");
		printf("-haplotype_ptype=[min_val,max_val], real range, default:[0.8,1]\n");
		printi("Range of genotype rate at extract loci.\n");
		printf("-haplotype_length=[min_val,max_val], integer range, default:[1,1000000]\n");
		printi("Range of haplotype size (in bp).\n");
		printf("-haplotype_variants=[min_val,max_val], integer range, default:[5,20]\n");
		printi("Range of number of variants in the haplotype.\n");
		printf("-haplotype_interval=0~1000000000, integer, default:0\n");
		printi("Minimum interval between extracted loci (in bp).\n");
		printf("-haplotype_alleles=[min_val,max_val], integer range, default:[2,65535]\n");
		printi("Range of number of alleles at the extracted locus.\n");
		printf("-haplotype_genotypes=[min_val,max_val], integer range, default:[2,65535]\n");
		printi("Range of number of genotypes at the extracted locus.\n");
		printf("\n");
	}

	//4. Conversion
	if (helpcat[4])
	{
		printf("4. Conversion\n");
		printf("-convert\n");
		printi("Converts filtered data (and extracted haplotype) into the input format of other software. The result is saved in *.convert.genepop.txt.\n");
		printf("-convert_format=genepop|spagedi|cervus|arlequin|structure|polygene|polyrelatedness|genodive, string, multiple selections, default:spagedi\n");
		printi("Target format, where genepop, cervus and arlequin formats only support diploids.\n");
		printf("\n");
	}

	//5. Genetic diversity indices
	if (helpcat[5])
	{
		printf("5. Genetic diversity indices\n");
		printf("-diversity\n");
		printi("Estimates the genetic diversity indices. Results are saved in *.diversity.txt.\n");
		printf("-diversity_level=pop|reg|tot|popXloc|regXloc|totXloc, string, multiple selections, default:loc,pop\n");
		printi("Output mean diversity across all loci in each population, each region or in the total population, or output diversity for each locus in each population, in each region or in the total population.\n");
		printf("-diversity_model=rcs|prcs|ces|pes, string, multiple selections, default:rcs\n");
		printi("Double-reduction model to calculate genotypic frequencies for polyploids.\n");
		printf("\n");
	}

	//6. Individual statistics
	if (helpcat[6])
	{
		printf("6. Individual statistics\n");
		printf("-indstat\n");
		printi("Calculates individual statistics (e.g., inbreeding coefficient, heterozygosity, and etc). Results are saved in *.indstat.txt.\n");
		printf("-indstat_type=hidx|lnpg|f|theta, string, multiple selections, default:hidx,lnpg\n");
		printi("Output statistics: heterozygosity index, natural logarithm of genotypic frequency, inbreeding coefficient and kinship coefficient.\n");
		printf("-indstat_model=rcs|prcs|ces|pes, string, multiple selections, default:rcs\n");
		printi("Double-reduction model to calculate genotypic frequencies for polyploids.\n");
		printf("-indstat_estimator=Ritland1996|Loiselle1995|Weir1996, string, multiple selections, default:Ritland1996\n");
		printi("Inbreeding coefficient and kinship coefficient (within an individual itself) estimators.\n");
		printf("-indstat_ref=pop|reg|total, string, multiple selections, default:pop\n");
		printi("Reference population: in the population, the region or the total population.\n");
		printf("-indstat_locus=all|each, string, multiple selections, default:all\n");
		printi("Output individual statistics for all loci or for each locus.\n");
		printf("\n");
	}

	//7. Genetic differentiation
	if (helpcat[7])
	{
		printf("7. Genetic differentiation\n");
		printf("-fst\n");
		printi("Estimates the genetic differentiation statistics. Results are saved in *.fst.txt.\n");
		printf("-fst_plot=yes|no, string, default:no\n");
		printi("Draws a heatmap for the results of the genetic differentiation estimation. Results are saved in *.fst.pdf.\n");
		printf("-fst_level=regXtot|popXtot|popXreg|reg|pop, string, multiple selections, default:pop\n");
		printi("Estimates the genetic differentiation among all regions, among all populations, among populations in each region, between any two regions, and between any two populations.\n");
		printf("-fst_estimator=Nei1973|Weir1984|Hudson1992|Slatkin1995|Hedrick2005|Jost2008|Huang2021_homo|Huang2021_aneu, string, multiple selections, default:Nei1973\n");
		printi("Genetic differentiation estimator, Nei1973 (Gst; Nei 1973, PNAS), Weir1984 (variance decomposition method, Weir & Cockerham 1984, Evolution), Hudson1992 (mean difference method, Hudson et al. 1992, Genetics), Slatkin1995 (Rst, Slatkin 1995, Genetics, for non-VCF/BCF input file only), Hedrick2005 (G'st; Hedrick 2005, Evolution), Jost2008 (D; Jost 2008, Molecular Ecology), Huang2021 (variance decomposition method for polyploid or aneuploid, Integrative Zoology).\n");
		printf("-fst_fmt=matrix|table, string, multiple selections, default:matrix\n");
		printi("Output format.\n");
		printf("-fst_locus=all|each, string, multiple selections, default:all\n");
		printi("Calculates genetic differentiation and performs a test for all loci or for each locus.\n");
		printf("-fst_test=genotype|allele, string, multiple selections, optional\n");
		printi("Tests the significance of differentiation by Fisher's G-test based on genotype distributions or allele distributions.\n");
		printf("\n");
	}

	//8. Genetic distance
	if (helpcat[8])
	{
		printf("8. Genetic distance\n");
		printf("-gdist\n");
		printi("Estimates the genetic distance. Results are saved in *.gdist.txt.\n");
		printf("-gdist_plot=yes|no, string, default:no\n");
		printi("Draws a heatmap for the results of genetic distance estimation. Results are saved in *.gdist.pdf.\n");
		printf("-gdist_weightmissing=yes|no, string, default:yes\n");
		printi("Use population/region allele frequency for missing data.\n");
		printf("-gdist_level=ind|pop|reg, string, multiple selections, default:pop\n");
		printi("Estimates the genetic distance between individuals, populations or regions.\n");
		printf("-gdist_estimator=Nei1972|Cavalli-Sforza1967|Reynolds1983|Nei1983|Euclidean|Goldstein1995|Nei1974|Roger1972|Slatkin_Nei1973|Slatkin_Weir1984|Slatkin_Hudson1992|Slatkin_Slatkin1995|Slatkin_Hedrick2005|Slatkin_Jost2008|Slatkin_Huang2021_homo|Slatkin_Huang2021_aneu|Reynolds_Nei1973|Reynolds_Weir1984|Reynolds_Hudson1992|Reynolds_Slatkin1995|Reynolds_Hedrick2005|Reynolds_Jost2008|Reynolds_Huang2021_homo|Reynolds_Huang2021_aneu, string, multiple selections, default:Nei1972\n");
		printi("Genetic distance estimators: Nei1972 (Ds, Nei 1972, Am Nat), Cavalli-Sforza1967 (Cavalli-Sforza & Edwards 1967, Am J Human Genet), Reynolds1983 (thetaW, Reynolds et al. Genetics, 1983), Nei1983  (Da, Nei 1983, J Mol Evol), Euclidean, Goldstein1995 (dmu2, Goldstein 1995, PNAS), Nei1974 (Dm, Nei & Roychoudhury 1974, Am J Human Genet), Roger1972 (Rogers 1972, Studies in Genetics), the Slatkin's transform d = Fst/(1-Fst) converts the range of Fst from [0,1] to [0, infinity), and the Reynolds's transformation d = -ln(1 - Fst).\n");
		printf("-gdist_fmt=matrix|table, string, multiple selections, default:matrix\n");
		printi("Output format.\n");
		printf("\n");
	}

	//9. Analysis of molecular variance
	if (helpcat[9])
	{
		printf("9. Analysis of molecular variance\n");
		printf("-amova\n");
		printi("Performs analysis of molecular variance. Results are saved in *.amova.txt.\n");
		printf("-amova_method=homoploid|aneuploid|likelihood, string, multiple selections, default:homoploid\n");
		printi("The homoploid method requires that all individuals are homoploids, and performs AMOVA and tests by extracting and permuting the dummy haplotypes. The aneuploid method supports aneuploids and permutes the alleles at each locus. The likelihood method also supports aneuploids, and uses the maximum-likelihood estimator to estimate F-statistics (Fis, Fic, Fit).\n");
		printf("-amova_mutation=iam|smm, string, multiple selections, default:iam\n");
		printi("Allele mutation model, iam denotes infinity alleles model (Fst like, distance between alleles is a binary variable) and smm denotes stepwise mutation model (Rst like, distance between alleles is the absoulte value of their difference in sizes). The smm model can only be applied for non-vcf input file and should use size as the allele identifier.\n");
		printf("-amova_ind=yes|no, string, multiple selections, default:yes\n");
		printi("Includes the individual level during AMOVA.\n");
		printf("-amova_trunc=yes|no, string, default:no\n");
		printi("Truncates negative variance component estimates to zero.\n");
		printf("-amova_test=yes|no, string, default:yes\n");
		printi("Evaluates the significance of each variance component and F-statistics (Fis, Fic, Fit, Fsc, Fst, Fct).\n");
		printf("-amova_nperm=99~99999999, integer, default:9999\n");
		printi("Number of permutations.\n");
		printf("-amova_pseudo=0~9999, integer, default:50\n");
		printi("Number of pseudo-permutations for the aneuploid method. Zero-value disables the pseudo-permutation.\n");
		printf("-amova_printss=yes|no, string, default:no\n");
		printi("Prints SS within individuals, populations and regions.\n");
		printf("\n");
	}

	//10. Population assignment
	if (helpcat[10])
	{
		printf("10. Population assignment\n");
		printf("-popas\n");
		printi("Assigns individuals to their natal population according to their genotypic frequencies in each population. Results are saved in *.popas.txt.\n");
		printf("-popas_plot=yes|no, string, default:no\n");
		printi("Draws a barplot for the results of population assignment. Results are saved in *.popas.pdf.\n");
		printf("-popas_model=rcs|prcs|ces|pes, string, multiple selections, default:rcs\n");
		printi("Double-reduction model to calculate genotypic frequencies for polyploids.\n");
		printf("-popas_level=pop|reg, string, multiple selections, default:pop\n");
		printi("Assigns individuals to populations or regions.\n");
		printf("-popas_error=0~0.2, real, default:0.01\n");
		printi("Mistype rate, used to avoid the probability of being zero.\n");
		printf("\n");
	}

	//11. Relatedness coefficient
	if (helpcat[11])
	{
		printf("11. Relatedness coefficient\n");
		printf("-relatedness\n");
		printi("Estimates pairwise relatedness between individuals. Results are saved in *.relatedness.txt.\n");
		printf("-relatedness_plot=yes|no, string, default:no\n");
		printi("Draws a heatmap for the results of relatedness estimation. Results are saved in *.relatedness.pdf.\n");
		printf("-relatedness_range=pop|reg|total, string, multiple selections, default:total\n");
		printi("Estimates pairwise relatedness between members within the same population, the same region or the total population.\n");
		printf("-relatedness_fmt=matrix|table, string, multiple selections, default:matrix\n");
		printi("Output format.\n");
		printf("-relatedness_estimator=Lynch1999|Wang2002|Thomas2010|Li1993|Queller1989|Huang2016A|Huang2016B|Milligan2003|Anderson2007|Huang2014|Huang2015|Ritland1996_modified|Loiselle1995_modified|Ritland1996|Loiselle1995|Weir1996, string, multiple selections, default=Huang2014\n");
		printi("Relatedness estimators: Huang2014 and Huang2015 support ploidy level <= 8, Ritland1996, Loiselle1995 and Weir1996 estimators support ploidy level <= 10, and other estimators only support diploids. Milligan2003, Anderson2007 and Huang2015 are maximum-likelihood estimators, and other estimators are method-of-moment estimators. Unbiased Ritland1996 and Loiselle1995 relatedness estimates are converted from kinship coefficient by eqn (7) of Huang et al. (2015, Mol Ecol Resour).\n");
		printf("\n");
	}

	//12. Kinship coefficient
	if (helpcat[12])
	{
		printf("12. Kinship coefficient\n");
		printf("-kinship\n");
		printi("Estimates kinship coefficient between individuals. Results are saved in *.kinship.txt.\n");
		printf("-kinship_plot=yes|no, string, default:no\n");
		printi("Draws a heatmap for the results of kinship estimation. Results are saved in *.kinship.pdf.\n");
		printf("-kinship_range=pop|reg|total, string, multiple selections, default:total\n");
		printi("Estimates the kinship coefficient between members within the same population, the same region or the total population.\n");
		printf("-kinship_fmt=matrix|table, string, multiple selections, default:matrix\n");
		printi("Output format.\n");
		printf("-kinship_estimator=Ritland1996|Loiselle1995|Weir1996, string, multiple selections, default:Ritland1996\n");
		printi("Kinship estimators. Supports a maximum level of ploidy of 10.\n");
		printf("\n");
	}

	//13. Principal coordinate analysis
	if (helpcat[13])
	{
		printf("13. Principal coordinate analysis\n");
		printf("-pcoa\n");
		printi("Performs a principal coordinate analysis for individuals, populations or regions. Results are saved in *.pcoa.txt.\n");
		printf("-pcoa_plot=yes|no, string, default:no\n");
		printi("Draws a scatter plot for the results of principal coordinate analysis. Results are saved in *.pcoa.pdf.\n");
		printf("-pcoa_level=ind|pop|reg, string, multiple selections, default:ind\n");
		printi("Ordinate individuals, populations or regions.\n");
		printf("-pcoa_dim=1~4096, default:3\n");
		printi("Number of dimensions to output.\n");
		printf("-pcoa_estimator=Nei1972|Cavalli-Sforza1967|Reynolds1983|Nei1983|Euclidean|Goldstein1995|Nei1974|Roger1972|Slatkin_Nei1973|Slatkin_Weir1984|Slatkin_Hudson1992|Slatkin_Slatkin1995|Slatkin_Hedrick2005|Slatkin_Jost2008|Slatkin_Huang2021_homo|Slatkin_Huang2021_aneu|Reynolds_Nei1973|Reynolds_Weir1984|Reynolds_Hudson1992|Reynolds_Slatkin1995|Reynolds_Hedrick2005|Reynolds_Jost2008|Reynolds_Huang2021_homo|Reynolds_Huang2021_aneu, string, multiple selections, default:Nei1972\n");
		printi("Genetic distance estimators. Results of Euclidean distance are equivalent to PCA.\n");
		printf("\n");
	}

	//14. Hierarchical clustering
	if (helpcat[14])
	{
		printf("14. Hierarchical clustering\n");
		printf("-cluster\n");
		printi("Perform hierarchical clustering for individuals, populations or regions. Results are saved in *.cluster.txt in standard tree format.\n");
		printf("-cluster_plot=yes|no, string, default:no\n");
		printi("Draws a dendrogram for the results of hierarchical clustering. Results are saved in *.cluster.pdf.\n");
		printf("-cluster_level=ind|pop|reg, string, multiple selections, default:ind\n");
		printi("Level of object in clustering: individuals, populations or regions.\n");
		printf("-cluster_method=NEAREST|FURTHEST|UPGMA|WPGMA|UPGMC|WPGMC|WARD, string, multiple selections, default:UPGMA\n");
		printi("Clustering methods.\n");
		printf("-cluster_estimator=Nei1972|Cavalli-Sforza1967|Reynolds1983|Nei1983|Euclidean|Goldstein1995|Nei1974|Roger1972|Slatkin_Nei1973|Slatkin_Weir1984|Slatkin_Hudson1992|Slatkin_Slatkin1995|Slatkin_Hedrick2005|Slatkin_Jost2008|Slatkin_Huang2021_homo|Slatkin_Huang2021_aneu|Reynolds_Nei1973|Reynolds_Weir1984|Reynolds_Hudson1992|Reynolds_Slatkin1995|Reynolds_Hedrick2005|Reynolds_Jost2008|Reynolds_Huang2021_homo|Reynolds_Huang2021_aneu, string, multiple selections, default:Nei1972\n");
		printi("Genetic distance estimators.\n");
		printf("\n");
	}

	//15. Bayesian clustering
	if (helpcat[15])
	{
		printf("15. Bayesian clustering\n");
		printf("-structure\n");
		printi("Perform Bayesian clustering. Results are saved in *.structure.txt and *.structure.k=?.rep=?_id=?.txt. The former is the summary and the latter is the result of each run.\n");
		printf("-structure_plot=yes|no, string, default:no\n");
		printi("Draws barplots for the results of Bayesian clustering. Results are saved in *.structure.pdf.\n");
		printf("-structure_decompress=yes|no, string, default:no\n");
		printi("Decompresses allele copies to accelerate calculations.\n");

		printf("Model:\n");
		printf("-structure_admix=yes|no, string, default:no\n");
		printi("ADMIX model assumes each allele copy at each locus within the same individual can be drawn from the different clusters. Otherwise, all allele copies within the same individual are drawn from a cluster in each iteration.\n");
		printf("-structure_locpriori=yes|no, string, default:no\n");
		printi("LOCPRIORI model uses the sample population to cluster individuals.\n");
		printf("-structure_f=yes|no, string, default:no\n");
		printi("F model assumes the allele frequencies in each cluster are correlated with that in the ancestral population.\n");

		printf("MCMC:\n");
		printf("-structure_krange=[min_val,max_val], integer range, default:[1,5]\n");
		printi("Range of K (number of clusters).\n");
		printf("-structure_nburnin=0~10000000, integer, default:1000\n");
		printi("Number of burn-in cycles.\n");
		printf("-structure_nreps=0~10000000, integer, default:10000\n");
		printi("Number of iterations after burn-in.\n");
		printf("-structure_nthinning=1~10000, integer, default:10\n");
		printi("Sampling interval to dememorize.\n");
		printf("-structure_nruns=1~1000, integer, default:1\n");
		printi("Number of independent runs for each value of K.\n");
		printf("-structure_nadmburnin=0~10000000, integer, default:500\n");
		printi("Number of admixture burn-in cycles. This parameter is used for the non-ADMIXTURE and non-LOCPRIORI models, which generates a proper initial state to prevent the Markov chain to be blocked in the local maxima.\n");
		
		printf("Misc:\n");
		printf("-structure_lambda=0~10000, real, default:1\n");
		printi("The initial value of lambda.\n");
		printf("-structure_inferlambda=yes|no, string, default:no\n");
		printi("Updated lambda in each iteration.\n");
		printf("-structure_stdlambda=0~10000, real, default:0.3\n");
		printi("Standard deviation of new lambda.\n");
		printf("-structure_maxlambda=0~10000, real, default:10\n");
		printi("Maximum of new lambda.\n");
		printf("-structure_difflambda=yes|no, string, default:no\n");
		printi("Use separate lambda for each cluster.\n");
		printf("-structure_diversity=yes|no, default:no\n");
		printi("Output diversity parameters for each cluster.\n");
		
		printf("ADMIX:\n");
		printf("-structure_alpha=0~10000, real, default:1\n");
		printi("The initial alpha, the priori Dirichlet parameter of admixture proportions Q.\n");
		printf("-structure_inferalpha=yes|no, string, default:yes\n");
		printi("Update alpha in ADMIX model.\n");
		printf("-structure_diffalpha=yes|no, string, default:no\n");
		printi("Use separate alpha for each cluster.\n");
		printf("-structure_uniformalpha=yes|no, string, default:yes\n");
		printi("Priori distribution for alpha, yes for uniform distribution and no for gamma distribution.\n");
		printf("-structure_stdalpha=0~10000, real, default:0.025\n");
		printi("Standard deviation of uniform priori distribution of alpha.\n");
		printf("-structure_maxalpha=0~10000, real, default:10\n");
		printi("Maximum of uniform priori distribution of alpha.\n");
		printf("-structure_alphapriora=0~10000, real, default:0.05\n");
		printi("One gamma priori distribution parameter.\n");
		printf("-structure_alphapriorb=0~10000, real, default:0.001\n");
		printi("The other gamma priori distribution parameter.\n");
		printf("-structure_metrofreq=0~1000000, integer, default:10\n");
		printi("Frequency of Metropolis-Hastings update of admixture proportions Q, set 0 to disable Metropolis-Hastings update.\n");
		
		printf("LOCPRIORI:\n");
		printf("-structure_r=0~10000, real, default:1\n");
		printi("Initial value of r, where r evaluates the informativeness of data for the sampling location.\n");
		printf("-structure_maxr=0~10000, real, default:20\n");
		printi("Maximum of new r.\n");
		printf("-structure_epsr=0~10000, real, default:0.1\n");
		printi("Max step value of new r.\n");
		printf("-structure_epseta=0~10000, real, default:0.025\n");
		printi("Max step value of new eta for the non-ADMIXTURE model, where eta reflects the relative proportion of individuals assigned to a cluster.\n");
		printf("-structure_epsgamma=0~10000, real, default:0.025\n");
		printi("Max step value of new gamma for the non-ADMIXTURE model, where gamma reflects the relative proportion of individuals sampled from a location and assigned to a cluster.\n");
		
		printf("FMODEL:\n");
		printf("-structure_pmeanf=0~10000, real, default:0.01\n");
		printi("Priori mean F, where F is the amount of drift from the ancestral population to the cluster k in the F model.\n");
		printf("-structure_pstdf=0~10000, real, default:0.05\n");
		printi("Priori standard deviation of F.\n");
		printf("-structure_stdf=0~10000, real, default:0.05\n");
		printi("Standard deviation of new F.\n");
		printf("-structure_singlef=yes|no, string, default:no\n");
		printi("Use the same F in all clusters.\n");
		printf("\n");
	}

	//16. PloidyInfer
	if (helpcat[16])
	{
		printf("16. Ploidy Inference\n");
		printf("-ploidyinfer\n");
		printi("Inference on individual ploidy levels (assuming autopolyploids). Results are saved in '*.ploidyinfer.txt'.\n");
		printf("-ploidyinfer_histogram=yes|no, string, default:no\n");
		printi("Output histogram data.\n");
		printf("-ploidyinfer_nbins=10~100, integer, default:20\n");
		printi("Number of bins used to plot the histogram.\n");
		printf("-ploidyinfer_type=1|2|3|4|5|6|7|8|9|10, string, multiple selections, default:2,4\n");
		printi("Possible ploidy levels.\n");
		printf("\n");
	}
	return true;
}

/* Run benchmark for SIMD instruction set */
TARGET void SimdBenchmark()
{
	printf("Begin SIMD benchmark test...\n");
	int SIMDTYPEBAK = SIMD_TYPE;
	RNG rng(g_seed_val);

#ifdef __aarch64__
	const char* simdstr[] = { "", "compile" , "   neon" };
#else
	const char* simdstr[] = { "", "compile" , "    sse" , "    avx" , " avx512" };
#endif

	InitCryptTable();
	int len = 65535, sep = 1023, rep = 500;

	double* a = new double[len];
	double* b = new double[len * sep];
	double* c = new double[len * sep];
	byte* ploidy = new byte[len * sep];

	double* cc[] = {
		c + sep * 0, c + sep * 1, c + sep * 2, c + sep * 3,
		c + sep * 4, c + sep * 5, c + sep * 6, c + sep * 7,
		c + sep * 8, c + sep * 9, c + sep * 10, c + sep * 11,
		c + sep * 12, c + sep * 13, c + sep * 14, c + sep * 15,
		c + sep * 16, c + sep * 17, c + sep * 18, c + sep * 19,
		c + sep * 20, c + sep * 21, c + sep * 22, c + sep * 23,
		c + sep * 24, c + sep * 25, c + sep * 26, c + sep * 27,
		c + sep * 28, c + sep * 29, c + sep * 30, c + sep * 31
	};

	double* cc2[] = {
		c + sep * 0, c + sep * 1, c + sep * 2, c + sep * 3,
		c + sep * 4, c + sep * 5, c + sep * 6, c + sep * 7,
		c + sep * 8, c + sep * 9, c + sep * 10, c + sep * 11,
		c + sep * 12, c + sep * 13, c + sep * 14, c + sep * 15,
		c + sep * 16, c + sep * 17, c + sep * 18, c + sep * 19,
		c + sep * 20, c + sep * 21, c + sep * 22, c + sep * 23,
		c + sep * 24, c + sep * 25, c + sep * 26, c + sep * 27,
		c + sep * 28, c + sep * 29, c + sep * 30, c + sep * 31
	};

	SetVal(ploidy, (byte)10, len * sep);
	for (int i = 0; i < len; ++i)
	{
		a[i] = rng.Uniform();
		b[i] = rng.Uniform();
		c[i] = rng.Uniform();
		b[i * sep] = rng.Uniform();
		c[i * sep] = rng.Uniform();
		c[i + sep * 0] = rng.Uniform();
		c[i + sep * 1] = rng.Uniform();
		c[i + sep * 2] = rng.Uniform();
		c[i + sep * 3] = rng.Uniform();
		c[i + sep * 4] = rng.Uniform();
		c[i + sep * 5] = rng.Uniform();
		c[i + sep * 6] = rng.Uniform();
		c[i + sep * 7] = rng.Uniform();
		c[i + sep * 8] = rng.Uniform();
		c[i + sep * 9] = rng.Uniform();
		c[i + sep * 10] = rng.Uniform();
		c[i + sep * 11] = rng.Uniform();
		c[i + sep * 12] = rng.Uniform();
		c[i + sep * 13] = rng.Uniform();
		c[i + sep * 14] = rng.Uniform();
		c[i + sep * 15] = rng.Uniform();
		ploidy[i] = (byte)rng.Next(10);
	}

	timepoint begin;
	
	for (SIMD_TYPE = 1; SIMD_TYPE <= SIMDTYPEBAK; ++SIMD_TYPE)
	{
		begin = GetNow();

		uint* A = (uint*)a;
		ushort* C = (ushort*)b;
		double re = 0;
		SetZero(A, len);
		C[0] = 0;
		for (int j = 0; j < rep * 200; ++j)
		{
			C[0] ++;
			SetVal(A, C, len);
			re += A[0] + A[1] + A[2] + A[3] + A[len - 4] + A[len - 3] + A[len - 2] + A[len - 1];
		}
		re += Sum((byte*)A, len * 4);
		double duration = GetElapse(begin);
		printf("SetVal (ushort to uint) %s: %0.5lf s, res = %u\n", simdstr[SIMD_TYPE], duration, *(uint*)(a + 32));
	}
	printf("\n");

	for (SIMD_TYPE = 1; SIMD_TYPE <= SIMDTYPEBAK; ++SIMD_TYPE)
	{
		begin = GetNow();
		b[0] = 3;
		double val = 0;
		uint64 re = 0;
		for (int j = 0; j < rep * 50; ++j)
		{
			b[0]++;
			re += GetMinIdx(b, len, val);
		}
		double duration = GetElapse(begin);
		printf("GetMinIdx %s: %0.5lf s, res = %lld, val = %.16le\n", simdstr[SIMD_TYPE], duration, re, val);
	}
	printf("\n");
	
	for (SIMD_TYPE = 1; SIMD_TYPE <= SIMDTYPEBAK; ++SIMD_TYPE)
	{
		begin = GetNow();
		b[0] = 3;
		double val1 = 0, val2 = 0;
		double v1 = 0, v2 = 0;
		for (int j = 0; j < rep * 150; ++j)
		{
			b[0]++;
			GetMinMaxVal(b, len, v1, v2);
			val1 += v1;
			val2 += v2;
		}
		double duration = GetElapse(begin);
		printf("GetMinMaxVal %s: %0.5lf s, res = %.16le, %.16le\n", simdstr[SIMD_TYPE], duration, val1, val2);
	}
	printf("\n");
	
	for (SIMD_TYPE = 1; SIMD_TYPE <= SIMDTYPEBAK; ++SIMD_TYPE)
	{
		begin = GetNow();
		b[0] = 3;
		double val = 0;
		for (int j = 0; j < rep * 150; ++j)
		{
			b[0]++;
			val += GetMinVal(b, len);
		}
		double duration = GetElapse(begin);
		printf("GetMinVal %s: %0.5lf s, res = %.16le\n", simdstr[SIMD_TYPE], duration, val);
	}
	printf("\n");

	for (SIMD_TYPE = 1; SIMD_TYPE <= SIMDTYPEBAK; ++SIMD_TYPE)
	{
		begin = GetNow();
		b[0] = 3;
		int64 val = 0;
		for (int j = 0; j < rep * 150; ++j)
		{
			b[0]++;
			val += GetMinVal((int64*)b, len);
		}
		double duration = GetElapse(begin);
		printf("GetMinVal %s: %0.5lf s, res = %lld\n", simdstr[SIMD_TYPE], duration, val);
	}
	printf("\n");

	for (SIMD_TYPE = 1; SIMD_TYPE <= SIMDTYPEBAK; ++SIMD_TYPE)
	{
		begin = GetNow();
		b[0] = 3;
		double val = 0;
		for (int j = 0; j < rep * 150; ++j)
		{
			b[0]++;
			val += GetMaxVal(b, len);
		}
		double duration = GetElapse(begin);
		printf("GetMaxVal (step = 1) %s: %0.5lf s, val = %.16le\n", simdstr[SIMD_TYPE], duration, val);
	}
	printf("\n");
	
	for (SIMD_TYPE = 1; SIMD_TYPE <= SIMDTYPEBAK; ++SIMD_TYPE)
	{
		begin = GetNow();
		b[0] = 3;
		double val = 0;
		for (int j = 0; j < rep * 3; ++j)
		{
			b[0]++;
			val += GetMaxVal(b, len, sep);
		}
		double duration = GetElapse(begin);
		printf("GetMaxVal (step = %d) %s: %0.5lf s, val = %.16le\n", sep, simdstr[SIMD_TYPE], duration, val);
	}
	printf("\n");

	for (SIMD_TYPE = 1; SIMD_TYPE <= SIMDTYPEBAK; ++SIMD_TYPE)
	{
		begin = GetNow();
		int64 re = 0;
		ploidy[0] = 0;
		for (int j = 0; j < rep / 22.7; ++j)
		{
			ploidy[0] = (ploidy[0] % 10) + 1;
			re += CountNonZero(ploidy, len * sep);
		}
		double duration = GetElapse(begin);
		printf("CountNonZero %s: %0.5lf s, res = %lld\n", simdstr[SIMD_TYPE], duration, re);
	}
	printf("\n");

	for (SIMD_TYPE = 1; SIMD_TYPE <= SIMDTYPEBAK; ++SIMD_TYPE)
	{
		begin = GetNow();
		double re = 0;
		b[0] = 0;
		for (int j = 0; j < rep * 240; ++j)
		{
			b[0] ++;
			re += Sum(b, len);
		}
		double duration = GetElapse(begin);
		printf("Sum (double) %s: %0.5lf s, res = %.16le\n", simdstr[SIMD_TYPE], duration, re);
	}
	printf("\n");

	for (SIMD_TYPE = 1; SIMD_TYPE <= SIMDTYPEBAK; ++SIMD_TYPE)
	{
		begin = GetNow();
		double re = 0;
		c[0] = 0;
		for (int j = 0; j < rep * 3; ++j)
		{
			c[0] ++;
			re = Sum(c, len, sep);
		}
		double duration = GetElapse(begin);
		printf("Sum (double, step = %d) %s: %0.5lf s, res = %.16le\n", sep, simdstr[SIMD_TYPE], duration, re);
	}
	printf("\n");

	for (SIMD_TYPE = 1; SIMD_TYPE <= SIMDTYPEBAK; ++SIMD_TYPE)
	{
		begin = GetNow();
		double re = 0;
		ploidy[0] = 0;
		for (int j = 0; j < rep * 180; ++j)
		{
			ploidy[0] = 10;
			re = Sum(ploidy, len);
		}
		double duration = GetElapse(begin);
		printf("Sum (byte) %s: %0.5lf s, res = %.16le\n", simdstr[SIMD_TYPE], duration, re);
	}
	printf("\n");

	for (SIMD_TYPE = 1; SIMD_TYPE <= SIMDTYPEBAK; ++SIMD_TYPE)
	{
		c[0] = 0;
		SetZero(a, len);
		begin = GetNow();
		double re = 0;
		for (int j = 0; j < rep * 2; ++j)
		{
			c[0] ++;
			Sum(a, cc, 32, len);
			SetVal(cc, cc2, 32);
			re += a[0 * 11] + a[1 * 11] + a[2 * 11] + a[3 * 11] + a[4 * 11] + a[5 * 11] + a[6 * 11] + a[7 * 11];
		}
		re += Sum(a, len);
		double duration = GetElapse(begin);
		printf("Sum (double, 32 arrays) %s: %0.5lf s, res = %.16le\n", simdstr[SIMD_TYPE], duration, re);
	}
	printf("\n");

	SIMD_TYPE = 1;
	Add(b, 0.542, len);//avoid prod being zero
	for (SIMD_TYPE = 1; SIMD_TYPE <= SIMDTYPEBAK; ++SIMD_TYPE)
	{
		double re = 0;
		b[0] = 0;
		begin = GetNow();
		for (int j = 0; j < rep * 180; ++j)
		{
			b[0] ++;
			re += Prod(b, len);
		}
		double duration = GetElapse(begin);
		printf("Prod (double) %s: %0.5lf s, res = %.16le\n", simdstr[SIMD_TYPE], duration, re);
	}
	printf("\n");

	SIMD_TYPE = 1;
	Add(c, 0.542, len * sep);//avoid prod being zero
	for (SIMD_TYPE = 1; SIMD_TYPE <= SIMDTYPEBAK; ++SIMD_TYPE)
	{
		begin = GetNow();
		double re = 0;
		c[0] = 0;
		for (int j = 0; j < rep * 3; ++j)
		{
			c[0] ++;
			re = Prod(c, len, sep);
		}
		double duration = GetElapse(begin);
		printf("Prod (double, step = %d) %s: %0.5lf s, res = %0.16le\n", sep, simdstr[SIMD_TYPE], duration, re);
	}
	printf("\n");

	for (SIMD_TYPE = 1; SIMD_TYPE <= SIMDTYPEBAK; ++SIMD_TYPE)
	{
		begin = GetNow();
		double re = 0;
		ploidy[0] = 0;
		for (int j = 0; j < rep * 160; ++j)
		{
			ploidy[0] = 10;
			re = SumSquare(ploidy, len);
		}
		double duration = GetElapse(begin);
		printf("SumSquare (byte) %s: %0.5lf s, res = %.16le\n", simdstr[SIMD_TYPE], duration, re);
	}
	printf("\n");

	for (SIMD_TYPE = 1; SIMD_TYPE <= SIMDTYPEBAK; ++SIMD_TYPE)
	{
		begin = GetNow();
		double re = 0;
		b[0] = 0;
		for (int j = 0; j < rep * 300; ++j)
		{
			b[0] ++;
			re = SumSquare(b, len);
		}
		double duration = GetElapse(begin);
		printf("SumSquare (double) %s: %0.5lf s, res = %.16le\n", simdstr[SIMD_TYPE], duration, re);
	}
	printf("\n"); 

	for (SIMD_TYPE = 1; SIMD_TYPE <= SIMDTYPEBAK; ++SIMD_TYPE)
	{
		begin = GetNow();
		double s1 = 0, s2 = 0;
		b[0] = 0;
		for (int j = 0; j < rep * 35; ++j)
		{
			b[0] ++;
			SumSumSquare(b, len, s1, s2);
		}
		double duration = GetElapse(begin);
		printf("SumSumSquare (double) %s: %0.5lf s, res = %.16le, %.16le\n", simdstr[SIMD_TYPE], duration, s1, s2);
	}
	printf("\n");

	for (SIMD_TYPE = 1; SIMD_TYPE <= SIMDTYPEBAK; ++SIMD_TYPE)
	{
		begin = GetNow();
		double re = 0;
		c[0] = 0;
		for (int j = 0; j < rep * 3; ++j)
		{
			c[0]++;
			re = SumProdDiv(a, b, c, sep, len);
		}
		double duration = GetElapse(begin);
		printf("SumProdDiv (double, step = %d) %s: %0.5lf s, res = %.16le\n", sep, simdstr[SIMD_TYPE], duration, log(re));
	}
	printf("\n");

	for (SIMD_TYPE = 1; SIMD_TYPE <= SIMDTYPEBAK; ++SIMD_TYPE)
	{
		begin = GetNow();
		double re = 0;
		c[0] = 0;
		for (int j = 0; j < rep * 200; ++j)
		{
			c[0]++;
			re = SumProd(b, c, len);
		}
		double duration = GetElapse(begin);
		printf("SumProd (double, step = 1) %s: %0.5lf s, res = %.16le\n", simdstr[SIMD_TYPE], duration, re);
	}
	printf("\n");

	for (SIMD_TYPE = 1; SIMD_TYPE <= SIMDTYPEBAK; ++SIMD_TYPE)
	{
		begin = GetNow();
		double re = 0;
		c[0] = 0;
		for (int j = 0; j < rep * 3.5; ++j)
		{
			c[0]++;
			re = SumProd(b, c, sep, len);
		}
		double duration = GetElapse(begin);
		printf("SumProd (double, step = %d) %s: %0.5lf s, res = %.16le\n", sep, simdstr[SIMD_TYPE], duration, re);
	}
	printf("\n");

	for (SIMD_TYPE = 1; SIMD_TYPE <= SIMDTYPEBAK; ++SIMD_TYPE)
	{
		SetZero(a, len);
		begin = GetNow();
		c[0] = 0;
		double re = 0;
		for (int j = 0; j < rep * 150; ++j)
		{
			c[0] ++;
			SetZero(a, len);
			Add(a, c, len);
			re += a[0] + a[1] + a[2] + a[3] + a[len - 4] + a[len - 3] + a[len - 2] + a[len - 1];
		}
		re += Sum(a, len);
		double duration = GetElapse(begin);
		printf("Add (double, array) %s: %0.5lf s, res = %.16le\n", simdstr[SIMD_TYPE], duration, re);
	}
	printf("\n");

	for (SIMD_TYPE = 1; SIMD_TYPE <= SIMDTYPEBAK; ++SIMD_TYPE)
	{
		SetZero(a, len);
		begin = GetNow();
		c[0] = 0;
		double re = 0;
		for (int j = 0; j < rep * 150; ++j)
		{
			c[0] ++;
			SetZero(a, len);
			Add(a, c[0], len);
			re += a[0] + a[1] + a[2] + a[3] + a[len - 4] + a[len - 3] + a[len - 2] + a[len - 1];
		}
		re += Sum(a, len);
		double duration = GetElapse(begin);
		printf("Add (double, constant) %s: %0.5lf s, res = %.16le\n", simdstr[SIMD_TYPE], duration, re);
	}
	printf("\n");

	for (SIMD_TYPE = 1; SIMD_TYPE <= SIMDTYPEBAK; ++SIMD_TYPE)
	{
		SetZero(a, len);
		begin = GetNow();
		int64 re = 0;
		int* A = (int*)a, * C = (int*)c;
		C[0] = 0;
		for (int j = 0; j < rep * 150; ++j)
		{
			C[0] ++;
			SetZero(A, len);
			Add(A, C, len);
			re += A[0] + A[1] + A[2] + A[3] + A[len - 4] + A[len - 3] + A[len - 2] + A[len - 1];
		}
		re += HashString((char*)A, len * 8);
		double duration = GetElapse(begin);
		printf("Add (int, array) %s: %0.5lf s, res = %lld\n", simdstr[SIMD_TYPE], duration, re);
	}
	printf("\n");

	for (SIMD_TYPE = 1; SIMD_TYPE <= SIMDTYPEBAK; ++SIMD_TYPE)
	{
		SetZero(a, len);
		begin = GetNow();
		int64 re = 0;
		int64* A = (int64*)a, * C = (int64*)c;
		C[0] = 0;
		for (int j = 0; j < rep * 100; ++j)
		{
			C[0] ++;
			SetZero(A, len);
			Add(A, C, len);
			re += A[0] + A[1] + A[2] + A[3] + A[len - 4] + A[len - 3] + A[len - 2] + A[len - 1];
		}
		re += HashString((char*)A, len * 8);
		double duration = GetElapse(begin);
		printf("Add (int64, array) %s: %0.5lf s, res = %lld\n", simdstr[SIMD_TYPE], duration, re);
	}
	printf("\n");

	for (SIMD_TYPE = 1; SIMD_TYPE <= SIMDTYPEBAK; ++SIMD_TYPE)
	{
		SetZero(a, len);
		SetVal(a, b, len);
		begin = GetNow();
		double re = 0;
		for (int j = 0; j < rep * 120; ++j)
		{
			a[0] ++;
			Mul(a, 1.000001, len);
			re += a[0] + a[1] + a[2] + a[3] + a[len - 4] + a[len - 3] + a[len - 2] + a[len - 1];
		}
		re += Sum(a, len);
		double duration = GetElapse(begin);
		printf("Mul (double, constant) %s: %0.5lf s, res = %.16le\n", simdstr[SIMD_TYPE], duration, re);
	}
	printf("\n");

	for (SIMD_TYPE = 1; SIMD_TYPE <= SIMDTYPEBAK; ++SIMD_TYPE)
	{
		SetZero(a, len);
		begin = GetNow();
		b[0] = 0;
		double re = 0;
		for (int j = 0; j < rep * 100; ++j)
		{
			b[0] ++;
			Mul(a, b, 1.41421356237, len);
			re += a[0] + a[1] + a[2] + a[3] + a[len - 4] + a[len - 3] + a[len - 2] + a[len - 1];
		}
		re += Sum(a, len);
		double duration = GetElapse(begin);
		printf("Mul (double, 1 array) %s: %0.5lf s, res = %.16le\n", simdstr[SIMD_TYPE], duration, re);
	}
	printf("\n");

	for (SIMD_TYPE = 1; SIMD_TYPE <= SIMDTYPEBAK; ++SIMD_TYPE)
	{
		SetZero(a, len);
		begin = GetNow();
		b[0] = 0;
		double re = 0;
		for (int j = 0; j < rep * 70; ++j)
		{
			b[0] ++;
			Mul(a, b, c, len);
			re += a[0] + a[1] + a[2] + a[3] + a[len - 4] + a[len - 3] + a[len - 2] + a[len - 1];
		}
		re += Sum(a, len);
		double duration = GetElapse(begin);
		printf("Mul (double, 2 array) %s: %0.5lf s, res = %.16le\n", simdstr[SIMD_TYPE], duration, re);
	}
	printf("\n");

	for (SIMD_TYPE = 1; SIMD_TYPE <= SIMDTYPEBAK; ++SIMD_TYPE)
	{
		double re = 0;
		b[0] = 0.1;
		begin = GetNow();
		for (int j = 0; j < rep * 22; ++j)
		{
			b[0] += 0.1;
			re += LogProd(b, len);
		}
		double duration = GetElapse(begin);
		printf("LogProd (double) %s: %0.5lf s, res = %.16le\n", simdstr[SIMD_TYPE], duration, re);
	}
	printf("\n");

	for (SIMD_TYPE = 1; SIMD_TYPE <= SIMDTYPEBAK; ++SIMD_TYPE)
	{
		begin = GetNow();
		double re = 0;
		c[0] = 0.1;
		for (int j = 0; j < rep * 2; ++j)
		{
			c[0] += 0.1;
			re = LogProd(c, len, sep);
		}
		double duration = GetElapse(begin);
		printf("LogProd (double, step = %d) %s: %0.5lf s, res = %.16le\n", sep, simdstr[SIMD_TYPE], duration, re);
	}
	printf("\n");

	for (SIMD_TYPE = 1; SIMD_TYPE <= SIMDTYPEBAK; ++SIMD_TYPE)
	{
		begin = GetNow();
		double re = 0;
		c[0] = 0.1;
		for (int j = 0; j < rep * 1.2; ++j)
		{
			c[0] += 0.1;
			re = LogProdDiv(b, c, len, sep);
		}
		double duration = GetElapse(begin);
		printf("LogProdDiv (double, step %d) %s: %0.5lf s, res = %.16le\n", sep, simdstr[SIMD_TYPE], duration, re);
	}
	printf("\n");

	for (SIMD_TYPE = 1; SIMD_TYPE <= SIMDTYPEBAK; ++SIMD_TYPE)
	{
		SetZero(a, len);
		begin = GetNow();
		b[0] = 0;
		double re = 0;
		for (int j = 0; j < rep * 60; ++j)
		{
			b[0] ++;
			AddProd(a, b, c, len);
			re += a[0] + a[1] + a[2] + a[3] + a[len - 4] + a[len - 3] + a[len - 2] + a[len - 1];
		}
		re += Sum(a, len);
		double duration = GetElapse(begin);
		printf("AddProd (double, array) %s: %0.5lf s, res = %.16le\n", simdstr[SIMD_TYPE], duration, re);
	}
	printf("\n");

	for (SIMD_TYPE = 1; SIMD_TYPE <= SIMDTYPEBAK; ++SIMD_TYPE)
	{
		SetZero(a, len);
		begin = GetNow();
		b[0] = 0;
		double re = 0;
		for (int j = 0; j < rep * 90; ++j)
		{
			b[0] ++;
			AddProd(a, b, c[0], len);
			re += a[0] + a[1] + a[2] + a[3] + a[len - 4] + a[len - 3] + a[len - 2] + a[len - 1];
		}
		re += Sum(a, len);
		double duration = GetElapse(begin);
		printf("AddProd (double, constant) %s: %0.5lf s, res = %.16le\n", simdstr[SIMD_TYPE], duration, re);
	}
	printf("\n");

	for (SIMD_TYPE = 1; SIMD_TYPE <= SIMDTYPEBAK; ++SIMD_TYPE)
	{
		SetZero(a, len);
		begin = GetNow();
		b[0] = 0;
		double re = 0;
		for (int j = 0; j < rep * 55; ++j)
		{
			b[0] ++;
			SetVal(a, b, len);
			Unify(a, len);
			re += a[0] + a[1] + a[2] + a[3] + a[len - 4] + a[len - 3] + a[len - 2] + a[len - 1];
		}
		re += Sum(a, len);
		double duration = GetElapse(begin);
		printf("Unify %s: %0.5lf s, res = %.16le\n", simdstr[SIMD_TYPE], duration, re);
	}
	printf("\n");

	ploidy[len * sep - 1] = 0;
	for (SIMD_TYPE = 1; SIMD_TYPE <= SIMDTYPEBAK; ++SIMD_TYPE)
	{
		begin = GetNow();
		int64 re = 0;
		ploidy[0] = 0;
		for (int j = 0; j < rep * 30; ++j)
		{
			ploidy[0] = (ploidy[0] % 10) + 1;
			re += (int64)StrNextIdx((char*)ploidy, '\n', 10000, len * sep);
		}
		double duration = GetElapse(begin);
		printf("StrNextIdx %s: %0.5lf s, res = %lld\n", simdstr[SIMD_TYPE], duration, re);
	}
	printf("\n");

	ploidy[len * sep - 1] = 0;
	for (SIMD_TYPE = 1; SIMD_TYPE <= SIMDTYPEBAK; ++SIMD_TYPE)
	{
		begin = GetNow();
		int64 re = 0;
		ploidy[0] = 0;
		for (int j = 0; j < rep * 0.06; ++j)
		{
			ploidy[0] = (ploidy[0] % 10) + 1;
			re += CountChar((char*)ploidy, '\n', len * sep);
		}
		double duration = GetElapse(begin);
		printf("CountChar %s: %0.5lf s, res = %lld\n", simdstr[SIMD_TYPE], duration, re);
	}
	printf("\n");

	SIMD_TYPE = SIMDTYPEBAK;
	delete[] a;
	delete[] b;
	delete[] c;
	delete[] ploidy;
}
