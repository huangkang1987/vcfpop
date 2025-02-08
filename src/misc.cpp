/* Misc */

#include "vcfpop.h"

template TARGET void BUCKET::FilterIndividualGT<double>(BUCKET* obucket, int nfilter, bool progress, int nthread);
template TARGET void BUCKET::FilterIndividualGT<float >(BUCKET* obucket, int nfilter, bool progress, int nthread);
template TARGET void BUCKET::FilterIndividualAD<double>(BUCKET* obucket, int nfilter, bool progress, int nthread);
template TARGET void BUCKET::FilterIndividualAD<float >(BUCKET* obucket, int nfilter, bool progress, int nthread);
template TARGET CPOINT CPOINT::GradientDescent<double>(void* Param, double (*func)(void* Param, CPOINT&, Mat<double>&, Mat<double>&), int dim, bool IsMinimize, double* Init);
template TARGET CPOINT CPOINT::GradientDescent<float >(void* Param, double (*func)(void* Param, CPOINT&, Mat<float >&, Mat<float >&), int dim, bool IsMinimize, double* Init);
template TARGET CPOINT CPOINT::DownHillSimplex<double>(void* Param, double (*func)(void* Param, CPOINT&, Mat<double>&, Mat<double>&), int dim, double sep, int nrep, bool IsMinimize);
template TARGET CPOINT CPOINT::DownHillSimplex<float >(void* Param, double (*func)(void* Param, CPOINT&, Mat<float >&, Mat<float >&), int dim, double sep, int nrep, bool IsMinimize);

#pragma pack(push, 1)

/* Set a bit in a variable */
TARGET void SetBit(byte& b, int pos, bool val)
{
	if (val)
		b |= (1 << pos);
	else
		b &= ~(1 << pos);
}

/* Get a bit in a variable */
TARGET bool GetBit(byte b, int pos)
{
	return b & (1 << pos);
}

/* Get number of different elements */
TARGET int GetNalleles(ushort* alleles, int ploidy)
{
	int nalleles = 0;
	for (int i = 0; i < ploidy; ++i)
	{
		bool isnew = true;
		ushort asi = alleles[i];
		for (int j = i - 1; j >= 0; --j)
		{
			if (asi == alleles[j])
			{
				isnew = false;
				break;
			}
		}
		if (isnew) nalleles++;
	}
	return nalleles;
}

/* Get current directory */
TARGET string GetCurDir()
{
	string re = std::filesystem::current_path().string();
	re.push_back(PATH_DELIM);
	return re;
}

/* Set current directory */
TARGET void SetCurDir(string dir)
{
	current_path(dir);
}

/* Get absoulte path */
TARGET string GetAbsPath(string file)
{
	string re = weakly_canonical(path(file)).string();
	return re;
}

/* Get parent path */
TARGET string GetParentPath(string file)
{
	string re = weakly_canonical(path(file)).parent_path().string();
	re.push_back(PATH_DELIM);
	return re;
}

/* Clear temp files */
TARGET void ClearTempFiles(string& dir)
{
    for (const auto& entry : fs::directory_iterator(dir))
    {
        if (entry.is_regular_file() && entry.path().filename().string().find("vcfpop_") == 0 && entry.path().extension() == ".tmp")
        {
			try
			{
                fs::remove(entry.path());
                std::cout << "Deleted: " << entry.path() << std::endl;
            }
			catch (const fs::filesystem_error& e)
			{
				//std::cerr << "Error: " << e.what() << std::endl;
			}
        }
    }
}

/* Pause console */
TARGET void Pause(void)
{
#ifdef _WIN64
	system("pause");
#else
	printf("Press ENTER to continue...\n");
	getchar();
#endif
}

/* Start a timer */
thread_local timepoint tic_toc_now;
TARGET void tic()
{
	tic_toc_now = std::chrono::steady_clock::now();
}

/* Show elapse time */
TARGET void toc()
{
	printf("Elapsed time is %0.6lf seconds.\n", GetElapse(tic_toc_now));
}

/* Get time */
TARGET timepoint GetNow()
{
	return std::chrono::steady_clock::now();
}

/* Get elapse time since begin */
TARGET double GetElapse(timepoint& begin)
{
	timepoint end = GetNow();
	return std::chrono::duration<double>(end - begin).count();
}

/* Record time when begin evaluation */
TARGET void EvaluationBegin()
{
	if (g_eval_val != 1) return;
	EVAL_BEGIN = GetNow();
}

/* Record time when evaluation is finished and append results to the evaluation file */
TARGET void EvaluationEnd(const char* text)
{
	if (g_eval_val != 1) return;
	double duration = GetElapse(EVAL_BEGIN);
	FILE* f1 = fopen((g_output_val + ".eval.txt").c_str(), "ab");
	fprintf(f1, "%s%c%0.3lf s\n", text, g_delimiter_val, duration);
	fclose(f1);
}

/* Return memory usage */
TARGET int64 GetMemoryUsage()
{
#if defined(_WIN64)
		/* Windows -------------------------------------------------- */
		PROCESS_MEMORY_COUNTERS info;
		GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
		return (size_t)info.WorkingSetSize;

#elif defined(__APPLE__)
		/* OSX ------------------------------------------------------ */
		struct mach_task_basic_info info;
		mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
		if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO,
			(task_info_t)&info, &infoCount) != KERN_SUCCESS)
			return (size_t)0L;      /* Can't access? */
		return (size_t)info.resident_size;

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
		/* Linux ---------------------------------------------------- */
		long rss = 0L;
		FILE* fp = NULL;
		if ((fp = fopen("/proc/self/statm", "r")) == NULL)
			return (size_t)0L;      /* Can't open? */
		if (fscanf(fp, "%*s%ld", &rss) != 1)
		{
			fclose(fp);
			return (size_t)0L;      /* Can't read? */
		}
		fclose(fp);
		return (size_t)rss * (size_t)sysconf(_SC_PAGESIZE);

#else
		/* AIX, BSD, Solaris, and Unknown OS ------------------------ */
		return (size_t)0L;          /* Unsupported. */
#endif
}

/* Return memory usage string */
TARGET void GetUsageString(char* str, int64 val)
{
	double size = val;
	if      (size > 1099511627776.0)
		sprintf(str, "%0.2lf TiB", size / 1099511627776.0);
	else if (size > 1073741824.0)
		sprintf(str, "%0.2lf GiB", size / 1073741824.0);
	else if (size > 1048576.0)
		sprintf(str, "%0.2lf MiB", size / 1048576.0);
	else if (size > 1024.0)
		sprintf(str, "%0.2lf KiB", size / 1024.0);
	else
		sprintf(str, "%0.0lf B  ", size);
}

/* Write number of individuals, populations, loci, raw data size and memory usage */
TARGET void EvaluationStat()
{
	if (g_eval_val != 1) return;
	FILE* f1 = fopen((g_output_val + ".eval.txt").c_str(), "ab");

	char memstr[30], sizestr[30];
	GetUsageString(memstr, GetMemoryUsage());
	GetUsageString(sizestr, TOTLEN_DECOMPRESS);

	fprintf(f1, "nind%c%d\n", g_delimiter_val, nind);
	fprintf(f1, "npop%c%d\n", g_delimiter_val, npop);
	fprintf(f1, "nloc%c%d\n", g_delimiter_val, nloc);
	fprintf(f1, "raw%c%s\n",  g_delimiter_val, sizestr);
	fprintf(f1, "memory%c%s\n", g_delimiter_val, memstr);

	fclose(f1);
}

/* Call an R script */
TARGET void RunRscript(string script)
{
	if (!FileExists(g_rscript_val.c_str()))
		Exit("\nError: Rscript binary executable is not found.\n");

	if (!FileExists((EXEDIR + "rscripts" + PATH_DELIM + script).c_str()))
		Exit("\nError: R script file %s is not found.\n", (EXEDIR + "rscripts" + PATH_DELIM + script).c_str());

	printf("\nPloting figure using %s ... \n", script.c_str());

#ifdef _WIN64
	string cmd = "\"\"" + g_rscript_val + "\" \"" +
		EXEDIR + "rscripts" + PATH_DELIM + script + "\" \"" +
		OUTFILE + "\"\"";
#else
	string cmd = "\"" + g_rscript_val + "\" \"" +
		EXEDIR + "rscripts" + PATH_DELIM + script + "\" \"" +
		OUTFILE + "\"";
#endif

	PLOT_THREAD.emplace_back(thread([cmd]() {
		std::system(cmd.c_str());
	}));
}

/* Exit program with a message */
TARGET void Exit(const char* fmt, ...)
{
	static atomic_flag first = ATOMIC_FLAG_INIT;
	if (!first.test_and_set())
	{
		va_list ap;
		va_start(ap, fmt);
		vprintf(fmt, ap);
		va_end(ap);
		Pause();
		exit(0);
	}
	else
	{
		Sleep(200000);
	}
}

#ifndef _CPOINT
TARGET CPOINT::CPOINT(int _dim)
{
	SetZero(this, 1);
	dim = _dim;
}

TARGET bool CPOINT::operator >(CPOINT& a)
{
	return lnL > a.lnL;
}

TARGET bool CPOINT::operator >=(CPOINT& a)
{
	return lnL >= a.lnL;
}

TARGET bool CPOINT::operator <(CPOINT& a)
{
	return lnL < a.lnL;
}

TARGET bool CPOINT::operator <=(CPOINT& a)
{
	return lnL <= a.lnL;
}

TARGET bool CPOINT::operator ==(CPOINT& a)
{
	return lnL == a.lnL;
}

TARGET bool CPOINT::operator !=(CPOINT& a)
{
	return lnL != a.lnL;
}

TARGET CPOINT CPOINT::operator =(const CPOINT& a)
{
	SetVal(this, (CPOINT*)&a, 1);
	return *this;
}

TARGET CPOINT CPOINT::operator +(const CPOINT& b)
{
	CPOINT re(*this);
	for (int i = 0; i < dim; ++i)
		re.unc_space[i] += b.unc_space[i];
	return re;
}

TARGET CPOINT CPOINT::operator -(const CPOINT& b)
{
	CPOINT re(*this);
	for (int i = 0; i < dim; ++i)
		re.unc_space[i] -= b.unc_space[i];
	return re;
}

TARGET CPOINT CPOINT::operator *(const double b)
{
	CPOINT re(*this);
	for (int i = 0; i < dim; ++i)
		re.unc_space[i] *= b;
	return re;
}

TARGET CPOINT CPOINT::operator /(const double b)
{
	CPOINT re(*this);
	for (int i = 0; i < dim; ++i)
		re.unc_space[i] /= b;
	return re;
}

TARGET CPOINT& CPOINT::operator +=(CPOINT& a)
{
	for (int i = 0; i < dim; ++i)
		unc_space[i] += a.unc_space[i];
	return *this;
}

TARGET CPOINT& CPOINT::operator -=(CPOINT& a)
{
	for (int i = 0; i < dim; ++i)
		unc_space[i] -= a.unc_space[i];
	return *this;
}

TARGET CPOINT& CPOINT::operator *=(double a)
{
	for (int i = 0; i < dim; ++i)
		unc_space[i] *= a;
	return *this;
}

TARGET CPOINT& CPOINT::operator /=(double a)
{
	for (int i = 0; i < dim; ++i)
		unc_space[i] /= a;
	return *this;
}

/* Distance in real space */
TARGET double CPOINT::DistanceReal(CPOINT& a)
{
	double s = 0;
	for (int i = 0; i < dim; ++i)
		s += (real_space[i] - a.real_space[i]) * (real_space[i] - a.real_space[i]);
	return MySqrt(s);
}

/* Distance in unc_space space */
TARGET double CPOINT::DistanceUnc(CPOINT& a)
{
	double s = 0;
	for (int i = 0; i < dim; ++i)
		s += (unc_space[i] - a.unc_space[i]) * (unc_space[i] - a.unc_space[i]);
	return MySqrt(s);
}

/* Break iteration if distance between points are smaller than eps */
TARGET bool CPOINT::IsBreak(CPOINT* xx, int dim, double eps)
{
	return xx[0].DistanceReal(xx[dim]) < eps;
}

/* Sort points */
TARGET void CPOINT::Order(CPOINT* xx, int dim)
{
	for (int i = 0; i <= dim; ++i)
		for (int j = i + 1; j <= dim; ++j)
			if (xx[i] < xx[j])
				Swap(xx[i], xx[j]);
}

// Gradient descent algorithm, need Gradient matrix but not Hessian matrix
template<typename REAL>
TARGET CPOINT CPOINT::GradientDescent(void* Param, double (*func)(void* Param, CPOINT&, rmat&, rmat&), int dim, bool IsMinimize, double *Init)
{
	CPOINT re0 = { 0 }, re1 = { 0 };
	SetZero(&re0, 1); SetZero(&re1, 1);
	re0.dim = re1.dim = dim;
	 
	int niter = 0, iter_max = 100 * dim;
	constexpr REAL grad_err_tol = 1E-08;
	constexpr REAL rel_sol_change_tol = 1E-08;
	constexpr REAL small_number = std::is_same_v<REAL, float> ? 1e-5 : 1e-8;

	double sign = IsMinimize ? -1 : 1;

	Mat<double> x0(re0.unc_space, dim, 1, false, true);
	Mat<double> x1(re1.unc_space, dim, 1, false, true);
	SetVal(re0.unc_space, 1.0, dim);
	SetVal(re1.unc_space, 1.0, dim);
	rmat G0, &H0 = *(rmat*)NULL;
	rmat G1, &H1 = *(rmat*)NULL;

	if (Init) SetVal(re0.unc_space, Init, dim);
	
	double f0 = 0, f1 = 0;
	REAL grad_err0 = 0, grad_err1 = 0, rel_sol_change0 = 0, rel_sol_change1 = 0;

	f0 = sign * func(Param, re0, G0, H0);
	re1.neval++;

	if (norm(G0) < small_number) G0 = G0 + small_number;
	grad_err0 = norm(G0);
	rel_sol_change0 = norm(x0 / (abs(x0) + small_number), 1);
	x1 = x0 + sign * 0.1 * G0 / norm(G0);

	f1 = sign * func(Param, re1, G1, H1);
	re1.neval++;
	grad_err1 = norm(G1);
	rel_sol_change1 = norm((x1 - x0) / (abs(x0) + small_number), 1);

	while (niter < iter_max)
	{
		niter++;

		rmat dG = G1 - G0;
		double dG_norm = norm(dG);
		double r1 = abs(trace((x1 - x0).t() * dG)) / (dG_norm * dG_norm);

		x0 = x1; f0 = f1; G0 = G1; grad_err0 = grad_err1; rel_sol_change0 = rel_sol_change1;

		x1 = x1 + sign * r1 * G1;

		f1 = sign * func(Param, re1, G1, H1);
		re1.neval++;

		grad_err1 = norm(G1);
		rel_sol_change1 = norm((x1 - x0) / (abs(x0) + small_number), 1);

		if (grad_err1 < grad_err_tol && rel_sol_change1 < rel_sol_change_tol)
			break;
	}

	re1.lnL = f1;
	return re1;
}

/* Down-Hill Simplex algorithm */
template<typename REAL>
TARGET CPOINT CPOINT::DownHillSimplex(void* Param, double (*func)(void* Param, CPOINT&, rmat&, rmat&), int dim, double sep, int nrep, bool IsMinimize)
{
	CPOINT xx[20] = { 0 };
	rmat &G = *(rmat*)NULL, &H = *(rmat*)NULL;

	int neval = 0;
	double sign = IsMinimize ? -1 : 1;
	for (int i = 0; i <= dim; ++i)
	{
		xx[i].dim = dim;
		SetVal(xx[i].unc_space, 0.01, dim);
		if (i > 0) xx[i].unc_space[i - 1] = 0.1;
		xx[i].lnL = sign * func(Param, xx[i], G, H);
		neval++;
	}

	for (int kk = 0; kk < nrep; ++kk)
	{
		//Order
		for (int searchcount = 0; ; ++searchcount)
		{
			CPOINT::Order(xx, dim);
			if (searchcount >= MAX_ITER_DOWNHILL || 
				(CPOINT::IsBreak(xx, dim, LIKELIHOOD_TERM) && abs(xx[0].real_space[0] - xx[1].real_space[0]) < LIKELIHOOD_TERM * 1e-3))
				break;

			//Reflect
			CPOINT x0 = xx[0];
			for (int i = 1; i < dim; ++i)
				x0 += xx[i];
			x0 /= dim;
			
			CPOINT xr = x0 + (x0 - xx[dim]);
			xr.lnL = sign * func(Param, xr, G, H);
			neval++;

			//Expansion
			//best
			if (xr > xx[0])
			{
				CPOINT xe = x0 + (xr - x0) * 2; 
				xe.lnL = sign * func(Param, xe, G, H);
				neval++;
				SetVal(xx + 1, xx, dim);
				xx[0] = xe > xr ? xe : xr;
				continue;
			}

			//better than second worst
			if (xr > xx[dim - 1])
			{
				xx[dim] = xr;
				continue;
			}

			//worse than second worst
			//Contraction
			CPOINT xc = x0 + (xx[dim] - x0) * 0.5;
			xc.lnL = sign * func(Param, xc, G, H);
			neval++;
			if (xc > xx[1])
			{
				xx[dim] = xc;
				continue;
			}

			//Reduction
			for (int i = 1; i <= dim; ++i)
			{
				xx[i] = (xx[0] + xx[i]) * 0.5;
				xx[i].lnL = sign * func(Param, xx[i], G, H);
				neval++;
			}
			continue;
		}

		//step 2
		for (int i = 1; i <= dim; ++i)
		{
			xx[i] = xx[0];
			xx[i].unc_space[i - 1] += xx[i].unc_space[i - 1] * sep;
			xx[i].lnL = sign * func(Param, xx[i], G, H);
			neval++;
		}

		//Order
		sep /= 2;
	}

	CPOINT::Order(xx, dim);
	xx[0].neval = neval;
	return xx[0];
}

#endif

#ifndef _MEMORY
/* Initialize */
TARGET MEMORY::MEMORY()
{
	cblock = 0;
	nblocks = 2;
	blocks = new MEMBLOCK[nblocks];
	SetZero(blocks, nblocks);
	blocks[0].size = 65536;
	blocks[0].bucket = new bool[blocks[0].size];
	InitLock(lock);
}

/* Uninitialize */
TARGET MEMORY::~MEMORY()
{
	if (cblock != 0xFFFFFFFF)
		for (int i = 0; i <= cblock; ++i)
			DEL(blocks[i].bucket);
	DEL(blocks);
	cblock = 0xFFFFFFFF;
    if (nblocks) UnInitLock(lock);
	nblocks = 0;
}

// Clear entries
TARGET void MEMORY::ClearMemory()
{
	for (int i = 0; i <= cblock; ++i)
		blocks[i].used = 0;
	cblock = 0;
}

// Expand number of blocks
TARGET void MEMORY::Expand()
{
	MEMBLOCK* nblock = new MEMBLOCK[(uint64)nblocks << 1];
	SetVal(nblock, blocks, cblock);
	SetZero(nblock + cblock, (nblocks << 1) - cblock);
	DEL(blocks); blocks = nblock;
	nblocks <<= 1;
}

// Allocate a small piece of memory
TARGET byte* MEMORY::Alloc(int size)
{
	//if (islock) 
	Lock(lock);
	while (size + blocks[cblock].used > blocks[cblock].size)
	{
		if (++cblock == nblocks) Expand();
		if (!blocks[cblock].bucket)
		{
			blocks[cblock].size = std::max(std::min(blocks[cblock - 1].size << 1, BIG_FILE ? 256 * 1024 * 1024 : 8 * 1024 * 1024), size);
			blocks[cblock].bucket = new bool[blocks[cblock].size];
		}
	}
	byte* addr = (byte*)(blocks[cblock].bucket + blocks[cblock].used);
	blocks[cblock].used += size;
	//if (islock) 
	UnLock(lock);
	return addr;
}
#endif

/* Virtual memory management class */
#ifndef _VMEMORY

	/* Do nothing */
	TARGET VMEMORY::VMEMORY()
	{
		size = 0;
		blocksize = 0;
		base_addr = tail_addr = head_addr = NULL;
		InitLock(lock);
	}

	/* Assign value but do not alloc */
	TARGET VMEMORY::VMEMORY(int64 _size)
	{
		//estimate file size
		int64 decomp_size = _size == 0 ?
			TOTLEN_COMPRESS * (g_format_val < 0 ? 40 : 1) : _size;

		if (decomp_size < 0x10000000)		// < 256 MiB
			blocksize = 0x200000;           // 2 MiB block
		else
		{
			blocksize = 1ull << CeilLog2(decomp_size / 500);
			BIG_FILE = true;
		}

		ApplyClearSlot(true);

		size = 0;
		head_addr = tail_addr = base_addr;
	}

	/* Do nothing */
	TARGET VMEMORY::~VMEMORY()
	{
        UnInitLock(lock);
	}

	/* Alloc or Clear a Virtual Memory slot, each have 16 TiB */
	TARGET void VMEMORY::ApplyClearSlot(bool isapply)
	{
		static atomic<byte> slots[16] = { 0 };

		if (isapply)
		{
			for (int i = 0; i < sizeof(slots); ++i)
			{
				byte exp = 0;
				if (slots[i].compare_exchange_strong(exp, 1))
				{
					base_addr = (byte*)(0x200000000000ull + i * 0x40000000000ull); // 4 TiB
					break;
				}
			}
		}
		else
		{
			slots[((uint64)base_addr - 0x200000000000ull) / 0x40000000000ull] = 0; // 4 TiB
			size = 0;
			blocksize = 0;
			base_addr = tail_addr = head_addr = NULL;
		}
	}

	/* Alloc a blocksize memory at the end */
	TARGET void VMEMORY::Alloc(byte* ed)
	{
		if (ed == NULL) ed = tail_addr + 1;

		Lock(lock);
		while (ed > tail_addr)
		{
#ifdef _WIN64
			if (!VirtualAlloc(tail_addr, blocksize, MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE))
#else		
			if (mmap(tail_addr, blocksize, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANON | MAP_FIXED, -1, 0) == (void*)-1)
#endif	
				Exit("\nError: cannot alloc %lld Mib memory.", blocksize / 1024 / 1024);

			VMEM_SIZE += blocksize;
			VMEM_NBLOCK++;
			tail_addr += blocksize;
			size += blocksize;
		}
		UnLock(lock);
	}

	/* Alloc a range of addresses in virtual memory */
	TARGET void VMEMORY::AllocRange(byte* st, byte* ed)
	{
		int stid = (st - base_addr) / blocksize;
		int edid = (ed - 1 - base_addr) / blocksize;

		if (flag.size() <= edid)
		{
			Lock(lock);
			if (flag.size() <= edid)
				flag.insert(flag.end(), edid + 1 - flag.size(), false);
			UnLock(lock);
		}

		bool alloc = false;
		for (int i = stid; i <= edid; ++i)
			if (flag[i] == false)
				alloc = true;

		if (!alloc) return;

		Lock(lock);
		for (int i = stid; i <= edid; ++i)
		{
			if (flag[i]) continue;
			
			byte* taddr = base_addr + i * blocksize;
#ifdef _WIN64
			if (!VirtualAlloc(taddr, blocksize, MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE))
#else		
			if (mmap(taddr, blocksize, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANON | MAP_FIXED, -1, 0) == (void*)-1)
#endif	
				Exit("\nError: cannot alloc %lld Mib memory.", blocksize / 1024 / 1024);

			VMEM_SIZE += blocksize;
			VMEM_NBLOCK++;
			size += blocksize;
			flag[i] = true;
		}

		tail_addr = base_addr + flag.size() * blocksize;
		UnLock(lock);
	}

	/* Unalloc all memory */
	TARGET void VMEMORY::UnAllocAll()
	{
		Lock(lock);

		while (size)
		{
			tail_addr -= blocksize;
#ifdef _WIN64
			if (!VirtualFree(tail_addr, 0, MEM_RELEASE))
#else
			if (munmap(tail_addr, blocksize))
#endif
				Exit("\nError: cannot unalloc virtual memory at %llx.", tail_addr);

			VMEM_SIZE -= blocksize;
			VMEM_NBLOCK --;
			size -= blocksize;

			if (size == 0)
				ApplyClearSlot(false);
		}

		vector<byte>().swap(flag);
		UnLock(lock);
	}

	/* If the data size in a block is zero, call UnAllocAddr */
	TARGET void VMEMORY::SubUnAlloc(atomic<int>& _size, int subv, int blockid)
	{
		if (_size.fetch_sub(subv) - subv == 0)
			UnAllocAddr(base_addr + blockid * blocksize);
	}

	/* Unalloc a blocksize memory */
	TARGET void VMEMORY::UnAllocAddr(byte* addr)
	{
		Lock(lock);

		if (size)
		{
			//head_addr has been added before calling
#ifdef _WIN64
			if (!VirtualFree(addr, 0, MEM_RELEASE))
#else
			if (munmap(addr, blocksize))
#endif
				Exit("\nError: cannot unalloc virtual memory at %llx.", head_addr - blocksize);

			VMEM_SIZE -= blocksize;
			VMEM_NBLOCK--;
			size -= blocksize;

			if (size == 0)
				ApplyClearSlot(false);
		}

		UnLock(lock);
	}

#endif

#ifndef _BUCKET

	/* Do nothing */
	TARGET BUCKET::BUCKET()
	{
		//new(this) VMEMORY();
		coffset = 0;
		//expanding = false;
	}

	/* Release memory */
	TARGET BUCKET::~BUCKET()
	{
		if (base_addr)
			UnAllocAll();

		coffset = 0; 
		size = 0;
		blocksize = 0;
		base_addr = NULL;
		tail_addr = NULL;
		head_addr = NULL;
		//expanding = false;

		vector<byte>().swap(flag);
	}

	/* Create bucket for genotype or allelic depth bucket */
	TARGET void BUCKET::CreateBucket()
	{
		offset.~LIST();//ok
		new(&offset) LIST<OFFSET>();
		new(this) VMEMORY(0);
		coffset = 0;
		//expanding = false;
	}

	/* Create bucket for genotype bucket in haplotype extraction */
	TARGET void BUCKET::CreateBucketGT(LIST<HAPLO_DUMMY_LOCUS>& hlocus)
	{
		coffset = 0;
		//expanding = false;

		for (int64 l = 0; l < hlocus.size; ++l)
		{
			uint64 gtsize = CeilLog2(hlocus[l].gsize);
			offset.Push(OFFSET { coffset, gtsize });
			coffset += (gtsize * nind + 7) >> 3; //xoffset
		}

		new(this) VMEMORY(coffset + 7);
		Alloc(base_addr + coffset + 7);
	}

	/* Create genotype table for unphase genotypes */
	TARGET void BUCKET::CreateBucketGT(SLOCUS* loc)
	{
		coffset = 0;
		//expanding = false;

		for (int64 l = 0; l < nloc; ++l)
		{
			uint64 gtsize = CeilLog2((int)loc[l].ngeno);//in bits 
			offset.Push(OFFSET{ coffset, gtsize });
			coffset += (gtsize * nind + 7) >> 3; //xoffset
		}

		new(this) VMEMORY(coffset + 7);
		Alloc(base_addr + coffset + 7);
	}

	/* Create genotype table for non-vcf/bcf input files */
	TARGET void BUCKET::CreateBucketGT(LOCUS* loc)
	{
		coffset = 0;
		//expanding = false;

		for (int64 l = 0; l < nloc; ++l)
		{
			uint64 gtsize = CeilLog2((int)loc[l].ngeno);//in bits 
			offset.Push(OFFSET{ coffset, gtsize });
			coffset += (gtsize * nind + 7) >> 3; //xoffset
		}

		new(this) VMEMORY(coffset + 7);
		Alloc(base_addr + coffset + 7);
	}

	/* Filter individual for genotype bucket */
	template<typename REAL>
	TARGET void BUCKET::FilterIndividualGT(BUCKET* obucket, int nfilt, bool progress, int nthread)
	{
		coffset = 0;
		//expanding = false;

		for (int64 l = 0; l < nloc; ++l)
		{
			uint64 gtsize = geno_bucket.offset[l].size;
			offset.Push(OFFSET{ coffset, gtsize });
			coffset += (gtsize * nfilt + 7) >> 3; //xoffset 3
		}

		if (true)
		{
			//move genotype id and allele depth table
			new(this) VMEMORY(coffset + 7);

			//dump obucket to a temp file
			time_t time1;
			time(&time1);
			tm* time2 = localtime(&time1);
			char* tmpname = new char[PATH_LEN];
			FILE* tmpfile = FOpen(tmpname, "wb+", "%svcfpop_%04d_%02d_%02d_%02d_%02d_%02d%s%d%s",
				g_tmpdir_val.c_str(), time2->tm_year + 1900, time2->tm_mon + 1, time2->tm_mday, time2->tm_hour, time2->tm_min, time2->tm_sec, ".filter", 1, ".tmp");
			fwrite(obucket->base_addr, 1, obucket->coffset, tmpfile);
			fclose(tmpfile);
			obucket->UnAllocAll();

			//map the temp file to memory
			FileMapping fmap;
			obucket->base_addr = fmap.MapingReadOnlyFile(tmpname);
			LIST<OFFSET>& ooffset = obucket->offset;

#pragma omp parallel  for num_threads(nthread)  schedule(dynamic, 1)
			for (int64 l = 0; l < nloc; ++l)
			{
				threadid = omp_get_thread_num();

				//check read pos
				uint64 gtsize = offset[l].size;
				uint64 gtlinesize = (gtsize * nfilt + 7) >> 3; //xoffset 4
				byte* write_pos = base_addr + offset[l].offset, * write_end = write_pos + gtlinesize;
				AllocRange(write_pos, write_end);

				GENO_READER rg(0, l, obucket);
				GENO_WRITER wg(l, this);

				for (int i = 0; i < nind; ++i)
				{
					uint gid = rg.Read();
					if (ainds<REAL>[i]) wg.Write(gid);
				}
				wg.FinishWrite();

				if (progress) PROGRESS_VALUE++;
			}

			//unmap the temp file
			fmap.UnMapingReadOnlyFile();

			//remove the temp file
			remove(tmpname);
			DEL(tmpname);
		}
		else
		{
			//calculate nbytes in each memblock in obucket
			byte* obase_addr = obucket->base_addr;
			int64 oblocksize = obucket->blocksize;
			LIST<OFFSET>& ooffset = obucket->offset;

			int nblock = (obucket->tail_addr - obucket->head_addr) / oblocksize;
			atomic<int>* nbytes = new atomic<int>[nblock];
			for (int j = 0; j < nblock - 1; ++j)
				nbytes[j] = oblocksize;
			nbytes[nblock - 1] = obucket->coffset - oblocksize * (nblock - 1);

			new(this) VMEMORY(coffset + 7);

#pragma omp parallel  for num_threads(nthread)  schedule(dynamic, 1)
			for (int64 l = 0; l < nloc; ++l)
			{
				threadid = omp_get_thread_num();
				byte* readpos = obase_addr + ooffset[l].offset;

				//check read pos
				uint64 gtsize = offset[l].size;
				uint64 gtlinesize = (gtsize * nfilt + 7) >> 3; //xoffset 4
				byte* write_pos = base_addr + offset[l].offset, * write_end = write_pos + gtlinesize;
				AllocRange(write_pos, write_end);

				GENO_READER rg(0, l, obucket);
				GENO_WRITER wg(l, this);

				for (int i = 0; i < nind; ++i)
				{
					uint gid = rg.Read();
					if (ainds<REAL>[i]) wg.Write(gid);
				}
				wg.FinishWrite();

				//Release memory if finish read a block
				int block_st = (readpos - obase_addr) / oblocksize;
				int block_ed = (readpos + gtlinesize - 1 - obase_addr) / oblocksize;

				if (block_st == block_ed)
					obucket->SubUnAlloc(nbytes[block_st], gtlinesize, block_st);
				else
				{
					obucket->SubUnAlloc(nbytes[block_st], obase_addr + (block_st + 1) * oblocksize - readpos, block_st);
					for (int j = block_st + 1; j < block_ed; ++j)
						obucket->SubUnAlloc(nbytes[j], nbytes[j].load(), j);
					obucket->SubUnAlloc(nbytes[block_ed], readpos + gtlinesize - (obase_addr + block_ed * oblocksize), block_ed);
				}

				if (progress) PROGRESS_VALUE++;
			}

			for (int j = 0; j < nblock; ++j)
				if (nbytes[j] > 0)
					obucket->SubUnAlloc(nbytes[j], nbytes[j].load(), j);

			DEL(nbytes);
		}
	}

	/* Filter individual for allele depth bucket */
	template<typename REAL>
	TARGET void BUCKET::FilterIndividualAD(BUCKET* obucket, int nfilt, bool progress, int nthread)
	{
		coffset = 0;
		//expanding = false;

		for (int64 l = 0; l < nloc; ++l)
		{
			int k2 = GetLoc(l).k;
			uint64 adsize = obucket->offset[l].size;
			offset.Push(OFFSET{ coffset, adsize });
			coffset += (adsize * k2 * nfilt + 7) >> 3;// xoffset 8
		}

		if (true)
		{
			//move genotype id and allele depth table
			new(this) VMEMORY(coffset + 7);

			//dump obucket to a temp file
			time_t time1;
			time(&time1);
			tm* time2 = localtime(&time1);
			char* tmpname = new char[PATH_LEN];
			FILE* tmpfile = FOpen(tmpname, "wb+", "%svcfpop_%04d_%02d_%02d_%02d_%02d_%02d%s%d%s",
				g_tmpdir_val.c_str(), time2->tm_year + 1900, time2->tm_mon + 1, time2->tm_mday, time2->tm_hour, time2->tm_min, time2->tm_sec, ".filter", 1, ".tmp");
			fwrite(obucket->base_addr, 1, obucket->coffset, tmpfile);
			fclose(tmpfile);
			obucket->UnAllocAll();

			//map the temp file to memory
			FileMapping fmap;
			obucket->base_addr = fmap.MapingReadOnlyFile(tmpname);
			LIST<OFFSET>& ooffset = obucket->offset;

#pragma omp parallel  for num_threads(nthread)  schedule(dynamic, 1)
			for (int64 l = 0; l < nloc; ++l)
			{
				threadid = omp_get_thread_num();
				int k2 = GetLoc(l).k;

				//check read pos
				uint64 adsize = offset[l].size;
				uint64 adlinesize = (adsize * k2 * nfilt + 7) >> 3; //xoffset 12
				byte* write_pos = base_addr + offset[l].offset, * write_end = write_pos + adlinesize;
				AllocRange(write_pos, write_end);

				GENO_READER rd(0, l, obucket);
				GENO_WRITER wd(l, this);

				for (int i = 0; i < nind; ++i)
				{
					for (int k = 0; k < k2; ++k)
					{
						uint depth = rd.Read();
						if (ainds<REAL>[i]) wd.Write(depth);
					}
				}
				wd.FinishWrite();

				if (progress) PROGRESS_VALUE++;
			}

			//unmap the temp file
			fmap.UnMapingReadOnlyFile();

			//remove the temp file
			remove(tmpname);
			DEL(tmpname);
		}
		else
		{
			//calculate nbytes in each memblock in obucket
			byte* obase_addr = obucket->base_addr;
			int64 oblocksize = obucket->blocksize;
			LIST<OFFSET>& ooffset = obucket->offset;

			int nblock = (obucket->tail_addr - obucket->head_addr) / oblocksize;
			atomic<int>* nbytes = new atomic<int>[nblock];
			for (int j = 0; j < nblock - 1; ++j)
				nbytes[j] = oblocksize;
			nbytes[nblock - 1] = obucket->coffset - oblocksize * (nblock - 1);

			new(this) VMEMORY(coffset + 7);

#pragma omp parallel  for num_threads(nthread)  schedule(dynamic, 1)
			for (int64 l = 0; l < nloc; ++l)
			{
				threadid = omp_get_thread_num();
				int k2 = GetLoc(l).k;
				byte* readpos = obase_addr + ooffset[l].offset;

				//check read pos
				uint64 adsize = offset[l].size;
				uint64 adlinesize = (adsize * k2 * nfilt + 7) >> 3; //xoffset 12
				byte* write_pos = base_addr + offset[l].offset, * write_end = write_pos + adlinesize;
				AllocRange(write_pos, write_end);

				GENO_READER rd(0, l, obucket);
				GENO_WRITER wd(l, this);

				for (int i = 0; i < nind; ++i)
				{
					for (int k = 0; k < k2; ++k)
					{
						uint depth = rd.Read();
						if (ainds<REAL>[i]) wd.Write(depth);
					}
				}
				wd.FinishWrite();

				//Release memory if finish read a block
				int block_st = (readpos - obase_addr) / oblocksize;
				int block_ed = (readpos + adlinesize - 1 - obase_addr) / oblocksize;

				if (block_st == block_ed)
					obucket->SubUnAlloc(nbytes[block_st], adlinesize, block_st);
				else
				{
					obucket->SubUnAlloc(nbytes[block_st], obase_addr + (block_st + 1) * oblocksize - readpos, block_st);
					for (int j = block_st + 1; j < block_ed; ++j)
						obucket->SubUnAlloc(nbytes[j], nbytes[j].load(), j);
					obucket->SubUnAlloc(nbytes[block_ed], readpos + adlinesize - (obase_addr + block_ed * oblocksize), block_ed);
				}

				if (progress) PROGRESS_VALUE++;
			}

			for (int j = 0; j < nblock; ++j)
				if (nbytes[j] > 0)
					obucket->SubUnAlloc(nbytes[j], nbytes[j].load(), j);

			DEL(nbytes);
		}
	}

	/* Filter locus for genotype bucket */
	TARGET void BUCKET::FilterLocusGT(BUCKET* obucket, bool progress, int nthread, MEMORY* mem, SLOCUS* nslocus, uint64* nlocus_pos)
	{
		//genotype id offset
		coffset = 0;
		//expanding = false;

		uint* newid = new uint[nloc];
		memset(newid, 0xFF, nloc * sizeof(uint));
		
		//calculate new offsets
		for (int64 l = 0, nl = 0; l < nloc; ++l)
		{
			if (!GetLoc(l).flag_pass) continue;

			newid[l] = nl++;
			uint64 gtsize = obucket->offset[l].size;
			uint64 gtlinesize = (gtsize * nind + 7) >> 3; //xoffset
			offset.Push(OFFSET{ coffset, gtsize });
			coffset += gtlinesize;
		}
		
		if (true)
		{
			//move genotype id and allele depth table
			new(this) VMEMORY(coffset + 7);
            
			//dump obucket to a temp file
			time_t time1;
			time(&time1);
			tm* time2 = localtime(&time1);
			char* tmpname = new char[PATH_LEN];
			FILE* tmpfile = FOpen(tmpname, "wb+", "%svcfpop_%04d_%02d_%02d_%02d_%02d_%02d%s%d%s",
				g_tmpdir_val.c_str(), time2->tm_year + 1900, time2->tm_mon + 1, time2->tm_mday, time2->tm_hour, time2->tm_min, time2->tm_sec, ".filter", 1, ".tmp");
			fwrite(obucket->base_addr, 1, obucket->coffset, tmpfile);
			fclose(tmpfile);
			obucket->UnAllocAll();
			
			//map the temp file to memory
			FileMapping fmap;
			byte* obase_addr = fmap.MapingReadOnlyFile(tmpname);
			LIST<OFFSET>& ooffset = obucket->offset;

#pragma omp parallel  for num_threads(nthread)  schedule(dynamic, 1)
			for (int64 l = 0; l < nloc; ++l)
			{
				threadid = omp_get_thread_num();
				byte* readpos = obase_addr + ooffset[l].offset;
				uint64 gtsize = ooffset[l].size;
				uint64 gtlinesize = (gtsize * nind + 7) >> 3; //xoffset
				
				if (GetLoc(l).flag_pass)
				{
					int64 nl = newid[l];
					byte* write_pos = base_addr + offset[nl].offset, * write_end = write_pos + gtlinesize;
					
					AllocRange(write_pos, write_end);
					memcpy(base_addr + offset[nl].offset, readpos, gtlinesize);
					
					//deep copy locus
					new(&nslocus[nl]) SLOCUS(mem[threadid], slocus[l]);
					if (uselocpos) nlocus_pos[nl] = GetLocPos(l);
				}
				
				if (progress) PROGRESS_VALUE++;
			}
			
			//unmap the temp file
			fmap.UnMapingReadOnlyFile();
			
			//remove the temp file
			remove(tmpname);
			DEL(tmpname);
		}
		else 
		{
			//calculate nbytes in each memblock in obucket
			byte* obase_addr = obucket->base_addr;
			int64 oblocksize = obucket->blocksize;
			LIST<OFFSET>& ooffset = obucket->offset;

			int nblock = (obucket->tail_addr - obucket->head_addr) / oblocksize;
			atomic<int>* nbytes = new atomic<int>[nblock];
			for (int j = 0; j < nblock - 1; ++j)
				nbytes[j] = oblocksize;
			nbytes[nblock - 1] = obucket->coffset - oblocksize * (nblock - 1);

			//move genotype id and allele depth table
			new(this) VMEMORY(coffset + 7);

#pragma omp parallel  for num_threads(nthread)  schedule(dynamic, 1)
			for (int64 l = 0; l < nloc; ++l)
			{
				threadid = omp_get_thread_num();
				byte* readpos = obase_addr + ooffset[l].offset;
				uint64 gtsize = ooffset[l].size;
				uint64 gtlinesize = (gtsize * nind + 7) >> 3; //xoffset

				if (GetLoc(l).flag_pass)
				{
					int64 nl = newid[l];
					byte* write_pos = base_addr + offset[nl].offset, * write_end = write_pos + gtlinesize;

					AllocRange(write_pos, write_end);
					memmove(base_addr + offset[nl].offset, readpos, gtlinesize);

					//deep copy locus
					new(&nslocus[nl]) SLOCUS(mem[threadid], slocus[l]);
					if (uselocpos) nlocus_pos[nl] = GetLocPos(l);
				}

				//Release memory if finish read a block
				int block_st = (readpos - obase_addr) / oblocksize;
				int block_ed = (readpos + gtlinesize - 1 - obase_addr) / oblocksize;

				if (block_st == block_ed)
					obucket->SubUnAlloc(nbytes[block_st], gtlinesize, block_st);
				else
				{
					obucket->SubUnAlloc(nbytes[block_st], obase_addr + (block_st + 1) * oblocksize - readpos, block_st);
					for (int j = block_st + 1; j < block_ed; ++j)
						obucket->SubUnAlloc(nbytes[j], nbytes[j].load(), j);
					obucket->SubUnAlloc(nbytes[block_ed], readpos + gtlinesize - (obase_addr + block_ed * oblocksize), block_ed);
				}

				if (progress) PROGRESS_VALUE++;
			}

			for (int j = 0; j < nblock; ++j)
				if (nbytes[j] > 0)
					obucket->SubUnAlloc(nbytes[j], nbytes[j].load(), j);

			DEL(nbytes);
		}
		
		DEL(newid);
	}

	/* Filter locus for allele depth bucket */
	TARGET void BUCKET::FilterLocusAD(BUCKET* obucket, bool progress, int nthread)
	{
		//genotype id offset
		coffset = 0;
		//expanding = false;

		uint* newid = new uint[nloc];
		memset(newid, 0xFF, nloc * sizeof(uint));

		//calculate new offsets
		for (int64 l = 0, nl = 0; l < nloc; ++l)
		{
			if (!GetLoc(l).flag_pass) continue;

			newid[l] = nl++;
			uint64 adsize = obucket->offset[l].size;
			uint64 adlinesize = (adsize * GetLoc(l).k * nind + 7) >> 3;//xoffset

			offset.Push(OFFSET{ coffset, adsize });
			coffset += adlinesize;
		}

		if (true)
		{
			//move genotype id and allele depth table
			new(this) VMEMORY(coffset + 7);

			//dump obucket to a temp file
			time_t time1;
			time(&time1);
			tm* time2 = localtime(&time1);
			char* tmpname = new char[PATH_LEN];
			FILE* tmpfile = FOpen(tmpname, "wb+", "%svcfpop_%04d_%02d_%02d_%02d_%02d_%02d%s%d%s",
				g_tmpdir_val.c_str(), time2->tm_year + 1900, time2->tm_mon + 1, time2->tm_mday, time2->tm_hour, time2->tm_min, time2->tm_sec, ".filter", 1, ".tmp");
			fwrite(obucket->base_addr, 1, obucket->coffset, tmpfile);
			fclose(tmpfile);
			obucket->UnAllocAll();

			//map the temp file to memory
			FileMapping fmap;
			byte* obase_addr = fmap.MapingReadOnlyFile(tmpname);
			LIST<OFFSET>& ooffset = obucket->offset;

#pragma omp parallel  for num_threads(nthread)  schedule(dynamic, 1)
			for (int64 l = 0; l < nloc; ++l)
			{
				threadid = omp_get_thread_num();
				byte* readpos = obase_addr + ooffset[l].offset;
				uint64 adsize = obucket->offset[l].size;
				uint64 adlinesize = (adsize * GetLoc(l).k * nind + 7) >> 3; //xoffset

				if (GetLoc(l).flag_pass)
				{
					int64 nl = newid[l];
					byte* write_pos = base_addr + offset[nl].offset, * write_end = write_pos + adlinesize;

					AllocRange(write_pos, write_end);
					memcpy(base_addr + offset[nl].offset, readpos, adlinesize);
				}

				if (progress) PROGRESS_VALUE++;
			}

			//unmap the temp file
			fmap.UnMapingReadOnlyFile();

			//remove the temp file
			remove(tmpname);
			DEL(tmpname);
		}
		else 
		{
			//calculate nbytes in each memblock in obucket
			byte* obase_addr = obucket->base_addr;
			int64 oblocksize = obucket->blocksize;
			LIST<OFFSET>& ooffset = obucket->offset;

			int nblock = (obucket->tail_addr - obucket->head_addr) / oblocksize;
			atomic<int>* nbytes = new atomic<int>[nblock];
			for (int j = 0; j < nblock - 1; ++j)
				nbytes[j] = oblocksize;
			nbytes[nblock - 1] = obucket->coffset - oblocksize * (nblock - 1);

			//move genotype id and allele depth table
			new(this) VMEMORY(coffset + 7);

#pragma omp parallel  for num_threads(nthread)  schedule(dynamic, 1)
			for (int64 l = 0; l < nloc; ++l)
			{
				threadid = omp_get_thread_num();
				byte* readpos = obase_addr + ooffset[l].offset;
				uint64 adsize = obucket->offset[l].size;
				uint64 adlinesize = (adsize * GetLoc(l).k * nind + 7) >> 3; //xoffset

				if (GetLoc(l).flag_pass)
				{
					int64 nl = newid[l];
					byte* write_pos = base_addr + offset[nl].offset, * write_end = write_pos + adlinesize;

					AllocRange(write_pos, write_end);
					memmove(base_addr + offset[nl].offset, readpos, adlinesize);
				}

				//Release memory if finish read a block
				int block_st = (readpos - obase_addr) / oblocksize;
				int block_ed = (readpos + adlinesize - 1 - obase_addr) / oblocksize;

				if (block_st == block_ed)
					obucket->SubUnAlloc(nbytes[block_st], adlinesize, block_st);
				else
				{
					obucket->SubUnAlloc(nbytes[block_st], obase_addr + (block_st + 1) * oblocksize - readpos, block_st);
					for (int j = block_st + 1; j < block_ed; ++j)
						obucket->SubUnAlloc(nbytes[j], nbytes[j].load(), j);
					obucket->SubUnAlloc(nbytes[block_ed], readpos + adlinesize - (obase_addr + block_ed * oblocksize), block_ed);
				}

				if (progress) PROGRESS_VALUE++;
			}

			for (int j = 0; j < nblock; ++j)
				if (nbytes[j] > 0)
					obucket->SubUnAlloc(nbytes[j], nbytes[j].load(), j);

			DEL(nbytes);
		}

		DEL(newid);
	}

	/* Allocate offset for genotype bucket */
	TARGET OFFSET BUCKET::AddOffsetGT(uint64 gtsize, uint64 id)
	{
		OFFSET re = OFFSET{ 
			coffset.fetch_add((gtsize * nind + 7) >> 3),
			gtsize };

		offset.Place(re, id, rwlock);

		//insufficient size, expand
		if (coffset > size)
			Alloc(base_addr + coffset);

		return re;
	}

	/* Allocate offset for allele depth bucket */
	TARGET OFFSET BUCKET::AddOffsetAD(uint64 adsize, uint64 id, uint k)
	{
		OFFSET re = OFFSET{
			coffset.fetch_add((adsize * k * nind + 7) >> 3),
			adsize };

		offset.Place(re, id, rwlock);

		//insufficient size, expand
		if (coffset > size)
			Alloc(base_addr + coffset);

		return re;
	}

	/* Use another bucket to replace this */
	TARGET void BUCKET::Replace(BUCKET& ref)
	{
		if (base_addr)
			UnAllocAll();

		//take vmem of ref
		size = ref.size;
		blocksize = ref.blocksize;
		base_addr = ref.base_addr;
		tail_addr = ref.tail_addr;
		head_addr = ref.head_addr;

		coffset = ref.coffset.load();
		//expanding = false;
		offset = ref.offset;
		flag.swap(ref.flag);

		//unref vmem of ref
		ref.size = 0;
		ref.blocksize = 0;
		ref.base_addr = NULL;
		ref.tail_addr = NULL;
		ref.head_addr = NULL;
		ref.coffset = 0;
		//ref.expanding = false;
	}

#endif

#ifndef _INCBUFFER

	/* Initialize, alloc 64 Kib */
	TARGET INCBUFFER::INCBUFFER()
	{
		len = INCREASE_BUFFER;
		data = (char*)malloc(len);
		derived = data;
	}

	/* Uninitialize */
	TARGET INCBUFFER::~INCBUFFER()
	{
		if (data) free(data);
		data = NULL;
		len = 0;
		derived = data;
	}

	/* Expand and memmove */
	TARGET void INCBUFFER::Move(char* src, int nlen)
	{
		Expand(nlen + 1);
		memmove(data, src, nlen);
		data[nlen] = '\0';
	}

	/* Read a line from a file and place in data */
	TARGET int64 INCBUFFER::Gets(FILE* handle)
	{
		int64 nread = 0, pos = FTell(handle);

		for (;;)
		{
			data[0] = '\0';

			if (g_format_val < 0)
				gzgets((gzFile)handle, data, (int)len);
			else
				fgets(data, (int)len, handle);

			nread = FTell(handle) - pos;
			if (nread >= len - 1)
			{
				Expand(len * 2);
				FSeek(handle, pos, SEEK_SET);
			}
			else
				return nread;
		}
	}

	/* Read a line from a file and place in p2 */
	TARGET int64 INCBUFFER::Gets(FILE* handle, char*& p2)
	{
		int64 offset = p2 == NULL ? 0 : p2 - data;
		int64 pos = FTell(handle);

		for (;;)
		{
			int64 rlen = len - offset;
			p2[0] = '\0';

			if (g_format_val < 0)
				gzgets((gzFile)handle, p2, (int)rlen);
			else
				fgets(p2, (int)rlen, handle);

			int64 nread = FTell(handle) - pos;

			if (nread >= rlen - 1)
			{
				Expand(len * 2);
				FSeek(handle, pos, SEEK_SET);
			}
			else
				return nread;
		}
	}

	/* Read a line from a file and place in p2 */
	TARGET int64 INCBUFFER::Gets(FILE* handle, int maxlen, char*& p2)
	{
		int64 offset = p2 == NULL ? 0 : p2 - data;

		if (offset + maxlen > len)
			Expand(len * 2);

		int64 pos = FTell(handle);
		data[offset] = '\0';

		if (g_format_val < 0)
			gzgets((gzFile)handle, data + offset, (int)len);
		else
			fgets(data + offset, (int)len, handle);

		return FTell(handle) - pos;
	}

	/* Expand buffer if nlen is greater than len*/
	TARGET void INCBUFFER::Expand(int nlen)
	{
		if (nlen > len)
		{
			char* tdata = (char*)realloc(data, nlen);

			if (tdata == NULL)
			{
				tdata = (char*)malloc(nlen);
				memmove(tdata, data, len);
				free(data);
			}

			derived += (int64)(tdata - data);

			data = tdata;
			len = nlen;
		}
	}

#endif 

#ifndef _VCFBUFFER
	/*	
	int64 geno_len;									//number of bytes to read in genotype
	char* line_pos;									//pointer of current line
	char* geno_pos;									//pointer of genotype
	int64 line_len;									//length of current line = read_len
	int64 data_size;								//number of bytes readed in this buffer
	int64 jloc;
	*/

	/* Read data from file */
	TARGET int64 VCFBUFFER::Read1In(FILE* handle)
	{
		if (g_format_val < 0)
			data_size = gzread((gzFile)handle, data, (uint)(len - 1));
		else
			data_size = fread(data, 1, len - 1, handle);

		data[data_size] = '\0';
		return data_size;
	}

	/* Read data from file, expand buffer if necessary, and set progress */
	TARGET void VCFBUFFER::Read1(FILE* handle)
	{
		uint64 opos = FTell(handle);
		Read1In(handle);
		
		while (CountChar(data, '\n', data_size) < 1)
		{
			//insufficient buffer, expand buffer, read again
			Expand((len - 1) * 16 + 1);
			FSeek(handle, opos, SEEK_SET);
			Read1In(handle);
		}

		//set line pos
		line_pos = data;

		//save offset
		file_offset = FOffset(handle);

		//add progress
		PROGRESS_VALUE += file_offset;

		//fix 1st line_len
		line_len = 0;
		jloc = 0;
	}

	/* Read data from file, expand buffer if necessary, and set progress */
	TARGET void VCFBUFFER::Read2(FILE* handle)
	{
		int64 remain_size = data_size - (line_pos - data);
		memmove(data, line_pos, remain_size);
		data_size = remain_size + FRead(data + remain_size, 1, len - 1 - remain_size, handle);
		data[data_size] = '\0';

		//save offset and add progress
		uint64 coffset = FOffset(handle);
		PROGRESS_VALUE += coffset - file_offset;
		file_offset = coffset;

		//set line pos
		line_pos = data;

		//fix 1st line_len
		line_len = 0;
	}
#endif

/* Initialize double reduction rates */
TARGET void InitAlpha()
{
	memset(ALPHA, 0, sizeof(ALPHA));

	//rcs
	for (int i = 1; i <= N_DRE_MODELT; ++i)
		for (int j = 1; j <= N_MAX_PLOIDY; ++j)
			ALPHA[i][j][0] = 1;

	//prcs
	ALPHA[2][4][0] = 6.0 / 7;		ALPHA[2][4][1] = 1.0 / 7;
	ALPHA[2][6][0] = 8.0 / 11;		ALPHA[2][6][1] = 3.0 / 11;
	ALPHA[2][8][0] = 40.0 / 65;		ALPHA[2][8][1] = 24.0 / 65;		ALPHA[2][8][2] = 1.0 / 65;
	ALPHA[2][10][0] = 168.0 / 323;	ALPHA[2][10][1] = 140.0 / 323;	ALPHA[2][10][2] = 15.0 / 323;

	//ces
	ALPHA[3][4][0] = 5.0 / 6;		ALPHA[3][4][1] = 1.0 / 6;
	ALPHA[3][6][0] = 7.0 / 10;		ALPHA[3][6][1] = 3.0 / 10;
	ALPHA[3][8][0] = 83.0 / 140;	ALPHA[3][8][1] = 54.0 / 140;	ALPHA[3][8][2] = 3.0 / 140;
	ALPHA[3][10][0] = 127.0 / 252;	ALPHA[3][10][1] = 110.0 / 252;	ALPHA[3][10][2] = 15.0 / 252;

	for (int E = 0; E <= 100; ++E)
	{
		double e = (double)E / 100;
		for (int p = 4; p <= N_MAX_PLOIDY; p += 2)
		{
			int D = p / 4;
			double apie[4] = { 0 };

			for (int i = 0; i <= D; ++i)
				apie[i] = BINOMIAL[p / 2][i] * BINOMIAL[p / 2 - i][p / 2 - 2 * i] * pow(2.0, p / 2 - 2 * i) / BINOMIAL[p][p / 2];

			for (int dd = 0; dd <= D; ++dd)
			{
				ALPHA[4 + E][p][dd] = 0;
				for (int j = dd; j <= D; ++j)
					for (int i = j; i <= D; ++i)
						ALPHA[4 + E][p][dd] += apie[i] * BINOMIAL[i][j] * pow(2.0, (double)-i) * BINOMIAL[j][dd] * pow(e, dd) * pow((double)1 - e, (double)j - dd);
			}
		}
	}
}

/* Calculate task in multiple threads */
TARGET void RunThreads(void (*Func) (int), void (*GuardFunc1) (int), void (*GuardFunc2) (int),
	int64 ntot, int64 nct, const char* info, int nthreads, bool isfirst, int nchars)
{
	if (isfirst)
	{
		printf("%s", info);
		fflush(stdout);
		PROGRESS_VALUE = PROGRESS_NOUTPUTED = PROGRESS_NOUTPUTED2 = 0;
		PROGRESS_TOTAL = ntot;
	}
	progress1 = progress2 = 0;
	SetZero(state_lock, NBUF);
	for (int i = 0; i < NBUF; ++i)
		state_lock[i] = (int64)i << 2;

	PROGRESS_CEND = PROGRESS_VALUE + nct;
	PROGRESS_CSTART = PROGRESS_VALUE;
	VLA_NEW(threads, thread, nthreads + 2);
	VLA_NEW(flags, bool, nthreads + 2);

	for (int i = 0; i < nthreads + 2; ++i)
	{
		flags[i] = true;
		if (i == nthreads && GuardFunc1)
			threads[i] = thread(GuardFunc1, i);
		else if (i == nthreads + 1 && GuardFunc2)
			threads[i] = thread(GuardFunc2, i);
		else if (i < nthreads && Func)
			threads[i] = thread(Func, i);
		else
			flags[i] = false;
	}

	while (PROGRESS_VALUE != PROGRESS_CEND)
	{
		Sleep(30);
		while (PROGRESS_VALUE * (double)nchars / PROGRESS_TOTAL > PROGRESS_NOUTPUTED)
		{
			if ((PROGRESS_NOUTPUTED + PROGRESS_NOUTPUTED2) % 10 == 0)
				printf("%lld", ((PROGRESS_NOUTPUTED + PROGRESS_NOUTPUTED2) / 10) % 10);
			else if ((PROGRESS_NOUTPUTED + PROGRESS_NOUTPUTED2) % 5 == 0)
				printf("o");
			else
				printf(".");
			PROGRESS_NOUTPUTED++;
			fflush(stdout);
		}
	}

	if (PROGRESS_CEND == PROGRESS_TOTAL)
	{
		while (PROGRESS_NOUTPUTED < nchars)
		{
			if ((PROGRESS_NOUTPUTED + PROGRESS_NOUTPUTED2) % 10 == 0)
				printf("%lld", ((PROGRESS_NOUTPUTED + PROGRESS_NOUTPUTED2) / 10) % 10);
			else if ((PROGRESS_NOUTPUTED + PROGRESS_NOUTPUTED2) % 5 == 0)
				printf("o");
			else
				printf(".");
			PROGRESS_NOUTPUTED++;
			fflush(stdout);
		}
	}

	for (int i = 0; i < nthreads + 2; ++i)
		if (flags[i])
			threads[i].join();

	VLA_DELETE(threads);
	VLA_DELETE(flags);
}

#pragma pack(pop)
