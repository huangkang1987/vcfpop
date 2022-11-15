/*
 vcfpop v1.05
 -- Perform population genetics analyses based on NGS data for haploids, diploids and polyploids.

 Author      Huang Kang
 Affiliation Northwest University
 Email       huangkang@nwu.edu.cn
 Update      2022/11/14
 */

#include "vcfpop.h"

 /* Main function */
TARGET int main(int _argc, char** _argv)
{
	//set executable path
	path exe_path(_argv[0]);
	EXEDIR = exe_path.parent_path().string();
	EXEDIR.push_back(PATH_DELIM);

	//copy parameters
	for (int i = 0; i < _argc; ++i)
		argv.push_back(_argv[i]);
	
	//dev-c++ debug argv.push_back("-p=par_1G.txt");
	
	//print splash
	printf("\nvcfpop v %s   %s\n", VERSION, DATE);
	printf("    Perform population genetics analyses for haploids, diploids, polyploids and ansioploids based on NGS data.\n");
	printf("    Huang Kang, Ph.D., Associate Prof.\n");
	printf("    Northwest University\n");

	//parse parameters
	if (PrintHelp()) return 0;

	//load parameters
	SetParameters(false);

	//clear temp files
	ClearTempFiles(g_tmpdir_val);

	//show directories
	printf("\n");
	printf("Working directory: \n  %s\n", CURDIR.c_str());
	printf("Temp directory: \n  %s\n", g_tmpdir_val.c_str());
	if (g_gpu_val == 2) ShowDevicesCUDA();
	printf("\n");

	//simd benchmark test
	threadid = 9999999;
	if (g_benchmark_val == 1)
	{
		SimdBenchmark<double>();
		SimdBenchmark<float >();
	}

	timepoint begin = GetNow();

	if (g_replot_val == 1)
	{
		//replot previous results
		string scripts[] = { "fst_plot.R", "gdist_plot.R", "popas_plot.R", "relatedness_plot.R",
							 "kinship_plot.R", "pcoa_plot.R", "cluster_plot.R", "structure_plot.R" };

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int i = 0; i < 8; ++i)
			RunRscript(scripts[i]);
	}
	else
	{
		//calculate current results
		if (g_float_val == 1)
			Calculate<float >();
		else
			Calculate<double>();
	}

	printf("\n\nCalculations finished in %0.3lf seconds.\n", GetElapse(begin));

	ReleaseParameters();

	return 0;
	//_ASSERTE(_CrtCheckMemory());
}
