/*
 vcfpop v1.03
 -- Perform population genetics analyses based on NGS data for haploids, diploids and polyploids.

 Author      Huang Kang
 Affiliation Northwest University
 Email       huangkang@nwu.edu.cn
 Update      2022/08/25
 */

#include "vcfpop.h"

 /* Main function */
TARGET int main(int _argc, char** _argv)
{
	try 
	{
		{
			path exe_path(_argv[0]);
			EXEDIR = exe_path.parent_path().string();
			EXEDIR.push_back(PATH_DELIM);

			for (int i = 0; i < _argc; ++i)
				argv.push_back(_argv[i]);
		}

		//title
		printf("\nvcfpop v %s   %s\n", VERSION, DATE);
		printf("    Perform population genetics analyses for haploids, diploids, polyploids and ansioploids based on NGS data.\n");
		printf("    Huang Kang, Ph.D., Associate Prof.\n");
		printf("    Northwest University\n");

		//parse parameters
		if (PrintHelp()) return 0;

		//load parameters
		SetParameters(false);

		ClearTempFiles(g_tmpdir_val);

		//show directories
		printf("    Working directory: %s\n", CURDIR.c_str());
		printf("    Temp directory: %s\n\n", g_tmpdir_val.c_str());

		threadid = 9999999;
		if (g_benchmark_val == 1) 
			SimdBenchmark();

		timepoint begin = GetNow();

		if (g_replot_val == 1)
		{
			string scripts[] = { "fst_plot.R", "gdist_plot.R", "popas_plot.R", "relatedness_plot.R",
								 "kinship_plot.R", "pcoa_plot.R", "cluster_plot.R", "structure_plot.R" };

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
			for (int i = 0; i < 8; ++i)
				RunRscript(scripts[i]);
		}
		else 
			Calculate();

		printf("\n\nCalculations finished in %0.3lf seconds.\n", GetElapse(begin));

		ReleaseParameters();
	}
	catch (std::exception& e)
	{
		FILE* f1 = fopen((g_output_val + ".err.txt").c_str(), "wb");
		if (f1)
		{
			fprintf(f1, "\nError: %s\n", e.what());
			fclose(f1);
		}
		Exit("\nError: %s\n", e.what());
	}
	catch (...)
	{
		Exit("\nError: unknown program error. \n");
	}

	return 0;
	//_ASSERTE(_CrtCheckMemory());
}