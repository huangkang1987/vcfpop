/*
	 vcfpop v1.08
	 -- Perform population genetics analyses based on NGS data for haploids, diploids and polyploids.

	 Author      Huang Kang
	 Affiliation Northwest University
	 Email       huangkang@nwu.edu.cn
	 Update      2025/02/01
 */

#include "vcfpop.h"

 /* Main function */
TARGET int main(int _argc, char** _argv)
{
	//disable openblas warning
	fclose(stderr);

	//disable multiple-thread mode of openblas
	openblas_set_num_threads(1);

	char buf[PATH_LEN];
#ifdef _WIN64
	GetModuleFileNameA(NULL, buf, PATH_LEN - 1);
#elif defined(__APPLE__)
	uint len = (uint)PATH_LEN;
	_NSGetExecutablePath(buf, &len);
#else
	readlink("/proc/self/exe", buf, PATH_LEN - 1);
#endif

	path exe_path(buf);
	EXEDIR = exe_path.parent_path().string();
	EXEDIR.push_back(PATH_DELIM);

	//copy parameters
	for (int i = 0; i < _argc; ++i)
		argv.push_back(_argv[i]);
	
	//print splash
	printf("\nvcfpop v %s   %s\n", VERSION, DATE);
	printf("    Perform population genetics analyses for haploids, diploids, polyploids and ansioploids based on NGS data.\n");
	printf("    Huang Kang, Ph.D., Prof.\n");
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
		if (slide && !ad && abs(g_format_val) <= BCF && slide_plot_val == 1)
			RunRscript("slide_plot.R");
		if (decay && !ad && abs(g_format_val) <= BCF && decay_plot_val == 1)
			RunRscript("decay_plot.R");
		if (block && !ad && abs(g_format_val) <= BCF && block_plot_val == 1)
			RunRscript("block_plot.R");
		if (gwas && !ad && abs(g_format_val) <= BCF && gwas_plot_val == 1)
			RunRscript("gwas_plot.R");
		if (gdist && gdist_plot_val == 1)
			RunRscript("gdist_plot.R");
		if (pcoa && pcoa_plot_val == 1)
			RunRscript("pcoa_plot.R");
		if (fst && fst_plot_val == 1)
			RunRscript("fst_plot.R");
		if (cluster && cluster_plot_val == 1)
			RunRscript("cluster_plot.R");
		if (relatedness && !ad && relatedness_plot_val == 1)
			RunRscript("relatedness_plot.R");
		if (kinship && !ad && kinship_plot_val == 1)
			RunRscript("kinship_plot.R");
		if (popas && !ad && popas_plot_val == 1)
			RunRscript("popas_plot.R");
		if (structure && !ad && structure_plot_val == 1)
			RunRscript("structure_plot.R");
	}
	else
	{
		//calculate current results
		if (g_float_val == 1)
			Calculate<float >();
		else
			Calculate<double>();
	}

	ReleaseParameters();

	//Wait untile all figures are plotted
	if (PLOT_THREAD.size())
	{
		bool alive = false;
		for (thread& thread : PLOT_THREAD)
			if (thread.joinable())
				alive = true;

		if (alive)
			printf("\nWaiting for output figure plotted.\n");
		
		for (thread& thread : PLOT_THREAD)
			if (thread.joinable())
				thread.join();
	}

	printf("\n\nCalculations finished in %0.3lf seconds.\n", GetElapse(begin));
	
	return 0;
	//_ASSERTE(_CrtCheckMemory());
}
