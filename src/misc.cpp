/* Misc */

#include "vcfpop.h"
#pragma pack(push, 1)

/*
TARGET void ReadWriteLock::lock_shared()
{
	for (;;)
	{
		int old_state = state;
		if ((old_state >> 16 == 0) &&	//if all writer finish
			state.compare_exchange_weak(old_state, old_state + 1)) //add nreader
			break;
	}
}

TARGET void ReadWriteLock::unlock_shared()
{
	state--;
}

TARGET void ReadWriteLock::lock()
{
	//add nwriter, disable new reader enter
	state += 0x10000;

	for (;;)
	{
		int old_state = state;
		if ((old_state & 0xFFFF) == 0 &&  //if all reader finish
			state.compare_exchange_weak(old_state, old_state | 0xFFFF)) //set low-16 bits to 0xFFFF
			break;
	}
}

TARGET void ReadWriteLock::unlock()
{
	//set low-16 bits to 0x0000 and minus nwriter 
	state -= 0x1FFFF;
}

*/

/* Count number of alleles from frequency array*/
TARGET int CountK(double* fre, int k2)
{
	int re = 0;
	for (int i = 0; i < k2; ++i)
		if (fre[i] > MIN_FREQ) re++;
	return re;
}

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

/* Clear temp files */
TARGET void ClearTempFiles(string& path)
{
#ifdef _WIN64
	_finddata_t file;
	string pattern = path + "vcfpop_*.tmp";
	int64 handle = _findfirst(pattern.c_str(), &file);
	if (handle != -1) do
	{
		remove((path + file.name).c_str());
	} while (_findnext(handle, &file) != -1);
	_findclose(handle);
#else
	DIR* dir = opendir(path.c_str());
	if (!dir)
		Exit("\nError: cannot open temp dir %s. \n", path.c_str());
	while (true)
	{
		dirent* p = readdir(dir);
		if (!p) break;
		if (strstr(p->d_name, "vcfpop_*.tmp"))
			remove(p->d_name);
	}
	closedir(dir);
#endif
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
	fprintf(f1, "%s\t%0.3lf s\n", text, duration);
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

	fprintf(f1, "nind\t%d\n", nind);
	fprintf(f1, "npop\t%d\n", npop);
	fprintf(f1, "nloc\t%d\n", nloc);
	fprintf(f1, "raw\t%s\n", sizestr);
	fprintf(f1, "memory\t%s\n", memstr);

	fclose(f1);
}

/* Pause console */
TARGET void RunRscript(string script)
{
	if (!FileExists(g_rscript_val.c_str()))
		Exit("\nError: Rscript binary executable is not found.\n");

	if (!FileExists((EXEDIR + "rscripts" + PATH_DELIM + script).c_str()))
		Exit("\nError: R script %s is not found.\n", script.c_str());

	printf("\nPloting figure using %s ... ", script.c_str());
	string cmd = "\"\"" + g_rscript_val + "\" \"" +
		EXEDIR + script + "\" \"" +
		OUTFILE + "\"\"";

	std::system(cmd.c_str());
}

/* Exit program with a message */
TARGET void Exit(const char* fmt, ...)
{
	static bool first = true;
	if (first)
	{
		first = false;
		va_list ap;
		int n = 0;
		va_start(ap, fmt);
		n = vprintf(fmt, ap);
		va_end(ap);
		Pause();
		exit(0);
	}
	else
		Sleep(200000);
}

#ifndef _CPOINT
TARGET CPOINT::CPOINT(int de)
{
	diff = 0;
	dim = de;
	SetVal(image, 0.0, 19);
}

TARGET bool CPOINT::operator >(CPOINT& a)
{
	return li > a.li;
}

TARGET bool CPOINT::operator >=(CPOINT& a)
{
	return li >= a.li;
}

TARGET bool CPOINT::operator <(CPOINT& a)
{
	return li < a.li;
}

TARGET bool CPOINT::operator <=(CPOINT& a)
{
	return li <= a.li;
}

TARGET bool CPOINT::operator ==(CPOINT& a)
{
	return li == a.li;
}

TARGET bool CPOINT::operator !=(CPOINT& a)
{
	return li != a.li;
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
		re.image[i] += b.image[i];
	return re;
}

TARGET CPOINT CPOINT::operator -(const CPOINT& b)
{
	CPOINT re(*this);
	for (int i = 0; i < dim; ++i)
		re.image[i] -= b.image[i];
	return re;
}

TARGET CPOINT CPOINT::operator *(const double b)
{
	CPOINT re(*this);
	for (int i = 0; i < dim; ++i)
		re.image[i] *= b;
	return re;
}

TARGET CPOINT CPOINT::operator /(const double b)
{
	CPOINT re(*this);
	for (int i = 0; i < dim; ++i)
		re.image[i] /= b;
	return re;
}

TARGET CPOINT& CPOINT::operator +=(CPOINT& a)
{
	for (int i = 0; i < dim; ++i)
		image[i] += a.image[i];
	return *this;
}

TARGET CPOINT& CPOINT::operator -=(CPOINT& a)
{
	for (int i = 0; i < dim; ++i)
		image[i] -= a.image[i];
	return *this;
}

TARGET CPOINT& CPOINT::operator *=(double a)
{
	for (int i = 0; i < dim; ++i)
		image[i] *= a;
	return *this;
}

TARGET CPOINT& CPOINT::operator /=(double a)
{
	for (int i = 0; i < dim; ++i)
		image[i] /= a;
	return *this;
}

/* Distance in real space */
TARGET double CPOINT::DistanceReal(CPOINT& a)
{
	double s = 0;
	for (int i = 0; i < dim; ++i)
		s += (real[i] - a.real[i]) * (real[i] - a.real[i]);
	return MySqrt(s);
}

/* Distance in image space */
TARGET double CPOINT::DistanceImage(CPOINT& a)
{
	double s = 0;
	for (int i = 0; i < dim; ++i)
		s += (image[i] - a.image[i]) * (image[i] - a.image[i]);
	return MySqrt(s);
}

/* Calculate real points */
TARGET void CPOINT::Image2Real()
{
	if (confine && !diff && dim % 2 == 0)
	{
		if (dim == 1)
		{
			real[0] = 1 / (1 + exp(-image[0]));
			real[1] = 1 - real[0];
		}
		else if (dim == 2)
		{
			double p1 = 1 / (1 + exp(-image[0]));
			double q1 = 1 / (1 + exp(-image[1]));
			double p0 = 1 - p1;
			double q0 = 1 - q1;
			real[0] = p1 * q1;
			real[1] = p0 * q1 + p1 * q0;
			real[2] = p0 * q0;
		}
		else if (dim == 4)
		{
			double p2 = 1 / (1 + exp(-image[0]));
			double p1 = (1 - p2) / (1 + exp(-image[1]));
			double q2 = 1 / (1 + exp(-image[2]));
			double q1 = (1 - q2) / (1 + exp(-image[3]));
			double p0 = 1 - p2 - p1;
			double q0 = 1 - q2 - q1;
			real[0] = p2 * q2;
			real[1] = p2 * q1 + p1 * q2;
			real[2] = p2 * q0 + p0 * q2 + p1 * q1;
			real[3] = p0 * q1 + p1 * q0;
			real[4] = p0 * q0;
		}
		else if (dim == 6)
		{
			double p3 = 1 / (1 + exp(-image[0]));
			double p2 = (1 - p3) / (1 + exp(-image[1]));
			double p1 = (1 - p3 - p2) / (1 + exp(-image[2]));
			double q3 = 1 / (1 + exp(-image[3]));
			double q2 = (1 - q3) / (1 + exp(-image[4]));
			double q1 = (1 - q3 - q2) / (1 + exp(-image[5]));
			double p0 = 1 - p3 - p2 - p1;
			double q0 = 1 - q3 - q2 - q1;
			real[0] = p3 * q3;
			real[1] = p3 * q2 + p2 * q3;
			real[2] = p3 * q1 + p2 * q2 * p1 * q3;
			real[3] = p3 * q0 + p2 * q1 + p1 * q2 + p0 * q3;
			real[4] = p2 * q0 + p1 * q1 + p0 * q2;
			real[5] = p1 * q0 + p0 * q1;
			real[6] = p0 * q0;
		}
		else if (dim == 8)
		{
			double p4 = 1 / (1 + exp(-image[0]));
			double p3 = (1 - p4) / (1 + exp(-image[1]));
			double p2 = (1 - p4 - p3) / (1 + exp(-image[2]));
			double p1 = (1 - p4 - p3 - p2) / (1 + exp(-image[3]));
			double q4 = 1 / (1 + exp(-image[4]));
			double q3 = (1 - q4) / (1 + exp(-image[5]));
			double q2 = (1 - q4 - q3) / (1 + exp(-image[6]));
			double q1 = (1 - q4 - q3 - q2) / (1 + exp(-image[7]));
			double p0 = 1 - p4 - p3 - p2 - p1;
			double q0 = 1 - q4 - q3 - q2 - q1;
			real[0] = p4 * q4;
			real[1] = p4 * q3 + p3 * q4;
			real[2] = p4 * q2 + p3 * q3 * p2 * q4;
			real[3] = p4 * q1 + p3 * q2 + p2 * q3 + p1 * q4;
			real[4] = p4 * q0 + p3 * q1 + p2 * q2 + p1 * q3 + p0 * q4;
			real[5] = p3 * q0 + p2 * q1 + p1 * q2 + p0 * q3;
			real[6] = p2 * q0 + p1 * q1 + p0 * q2;
			real[7] = p1 * q0 + p0 * q1;
			real[8] = p0 * q0;
		}
	}
	else
	{
		real[dim] = 1;
		for (int i = 0; i < dim; ++i)
		{
			real[i] = real[dim] / (1 + exp(-image[i]));
			real[dim] -= real[i];
		}
	}
}

/* Calculate real points */
TARGET void CPOINT::Image2RealSelfing()
{
	real[0] = 1.0 / (1 + exp(image[0]));
	if (real[0] <     MIN_FREQ) real[0] =     MIN_FREQ;
	if (real[0] > 1 - MIN_FREQ) real[0] = 1 - MIN_FREQ;
}

/* Break iteration if distance between points are smaller than eps */
TARGET bool CPOINT::IsBreak(CPOINT xx[9], double eps)
{
	return xx[0].DistanceReal(xx[xx[0].dim]) < eps;
}

/* Sort points */
TARGET void CPOINT::Order(CPOINT xx[9])
{
	int dim = xx[0].dim;
	for (int i = 0; i <= dim; ++i)
		for (int j = i + 1; j <= dim; ++j)
			if (xx[i] < xx[j])
				Swap(xx[i], xx[j]);
}

/* Down-Hill Simplex algorithm */
TARGET CPOINT CPOINT::DownHillSimplex(int dim, int diff, bool confine, double sep, int nrep, double (*Likelihood)(CPOINT&, void**), void* *Param)
{
	CPOINT xx[20] = { 0 };

	for (int i = 0; i <= dim; ++i)
	{
		xx[i].dim = dim;
		xx[i].diff = diff;
		xx[i].confine = confine;
		if (i > 0) xx[i].image[i - 1] = 0.1;
		xx[i].li = Likelihood(xx[i], Param);
	}

	double likestop = LIKELIHOOD_TERM;
	for (int kk = 0; kk < nrep; ++kk)
	{
		//Order
		for (int searchcount = 0; ; ++searchcount)
		{
			CPOINT::Order(xx);
			if (searchcount >= MAX_ITER_DOWNHILL || 
				(CPOINT::IsBreak(xx, LIKELIHOOD_TERM) && abs(xx[0].real[0] - xx[1].real[0]) < LIKELIHOOD_TERM * 1e-3))
				break;

			//Reflect
			CPOINT x0 = xx[0];
			for (int i = 1; i < dim; ++i)
				x0 += xx[i];
			x0 /= dim;
			
			CPOINT xr = x0 + (x0 - xx[dim]);
			xr.li = Likelihood(xr, Param); 

			//Expansion
			//best
			if (xr > xx[0])
			{
				CPOINT xe = x0 + (xr - x0) * 2; 
				xe.li = Likelihood(xe, Param); 
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
			xc.li = Likelihood(xc, Param); 
			if (xc > xx[1])
			{
				xx[dim] = xc;
				continue;
			}

			//Reduction
			for (int i = 1; i <= dim; ++i)
			{
				xx[i] = (xx[0] + xx[i]) * 0.5;
				xx[i].li = Likelihood(xx[i], Param);
			}
			continue;
		}

		//step 2
		for (int i = 1; i <= dim; ++i)
		{
			xx[i] = xx[0];
			xx[i].image[i - 1] *= (1 - sep);
			xx[i].li = Likelihood(xx[i], Param);
		}

		//Order
		sep /= 2;
		likestop /= 2;
	}

	CPOINT::Order(xx);
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
			delete[] blocks[i].bucket;
	if (blocks) delete[] blocks;
	blocks = NULL;
	cblock = 0xFFFFFFFF;
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
	delete[] blocks; blocks = nblock;
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
			blocks[cblock].size = Max(Min(blocks[cblock - 1].size << 1, BIG_FILE ? 256 * 1024 * 1024 : 8 * 1024 * 1024), size);
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
	TARGET VMEMORY::VMEMORY(int64 size)
	{
		//estimate file size
		int64 decomp_size = size == 0 ?
			TOTLEN_COMPRESS * (g_format_val < 0 ? 40 : 1) : size;

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

		InitLock(lock);
	}

	/* Do nothing */
	TARGET VMEMORY::~VMEMORY()
	{

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
					base_addr = 0x100000000000ull + (byte*)(i * 0x40000000000ull); // 4 TiB
					break;
				}
			}
		}
		else
		{
			slots[((uint64)base_addr - 0x100000000000ull) / 0x40000000000ull] = 0; // 4 TiB
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
			if (mmap(tail_addr, blocksize, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANON, -1, 0) == (void*)-1)
#endif	
				Exit("\nError: cannot alloc %lld Mib memory.", blocksize / 1024 / 1024);

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
			if (mmap(taddr, blocksize, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANON, -1, 0) == (void*)-1)
#endif	
				Exit("\nError: cannot alloc %lld Mib memory.", blocksize / 1024 / 1024);

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

			size -= blocksize;

			if (size == 0)
				ApplyClearSlot(false);
		}

		vector<byte>().swap(flag);
		UnLock(lock);
	}

	/* If the data size in a block is zero, call UnAllocAddr */
	TARGET void VMEMORY::SubUnAlloc(atomic<int>& size, int subv, int blockid)
	{
		if (size.fetch_sub(subv) - subv == 0)
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
		new(this) VMEMORY();
		coffset = 0;
		expanding = false;
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
		expanding = false;

		vector<byte>().swap(flag);
	}

	/* Create bucket for genotype or allelic depth bucket */
	TARGET void BUCKET::CreateBucket()
	{
		offset.~LIST();//ok
		new(&offset) LIST<OFFSET>();
		new(this) VMEMORY(0);
		coffset = 0;
		expanding = false;
	}

	/* Create bucket for genotype bucket in haplotype extraction */
	TARGET void BUCKET::CreateBucketGT(LIST<HAPLO_DUMMY_LOCUS>& hlocus)
	{
		coffset = 0;
		expanding = false;

		for (int64 l = 0; l < hlocus.size; ++l)
		{
			uint64 gsize = CeilLog2(hlocus[l].gsize);
			offset.Push(OFFSET { coffset, gsize });
			coffset += (gsize * nind + 7) >> 3; // xoffset 1
		}

		new(this) VMEMORY(coffset + 7);
		Alloc(base_addr + coffset + 7);
	}

	/* Create genotype table for non-vcf/bcf input files */
	TARGET void BUCKET::CreateBucketGT(LOCUS* loc)
	{
		coffset = 0;
		expanding = false;

		for (int64 l = 0; l < nloc; ++l)
		{
			uint64 gsize = CeilLog2((int)loc[l].ngeno);//in bits 
			offset.Push(OFFSET{ coffset, gsize });
			coffset += (gsize * nind + 7) >> 3;// xoffset 2
		}

		new(this) VMEMORY(coffset + 7);
		Alloc(base_addr + coffset + 7);
	}

	/* Filter individual for genotype bucket */
	TARGET void BUCKET::FilterIndividualGT(BUCKET* obucket, int nfilt, bool progress, int nthread)
	{
		coffset = 0;
		expanding = false;

		for (int64 l = 0; l < nloc; ++l)
		{
			uint64 gtsize = geno_bucket.offset[l].size;
			offset.Push(OFFSET{ coffset, gtsize });
			coffset += (gtsize * nfilt + 7) >> 3;// xoffset 7
		}

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
			uint64 gtlinesize = (gtsize * nfilt + 7) >> 3; //xoffset 11
			byte* write_pos = base_addr + offset[l].offset, * write_end = write_pos + gtlinesize;
			AllocRange(write_pos, write_end);

			GENO_READER rg(0, l, obucket);
			GENO_WRITER wg(l, this);

			for (int i = 0; i < nind; ++i)
			{
				uint gid = rg.Read();
				if (ainds[i]) wg.Write(gid);
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

		delete[] nbytes;
	}

	/* Filter individual for allele depth bucket */
	TARGET void BUCKET::FilterIndividualAD(BUCKET* obucket, int nfilt, bool progress, int nthread)
	{
		coffset = 0;
		expanding = false;

		for (int64 l = 0; l < nloc; ++l)
		{
			int k2 = GetLoc(l).k;
			uint64 adsize = obucket->offset[l].size;
			offset.Push(OFFSET{ coffset, adsize });
			coffset += (adsize * k2 * nfilt + 7) >> 3;// xoffset 8
		}

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
			byte* write_pos = base_addr + offset[l].offset, *write_end = write_pos + adlinesize;
			AllocRange(write_pos, write_end);

			GENO_READER rd(0, l, obucket);
			GENO_WRITER wd(l, this);

			for (int i = 0; i < nind; ++i)
			{
				for (int k = 0; k < k2; ++k)
				{
					uint depth = rd.Read();
					if (ainds[i]) wd.Write(depth);
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

		delete[] nbytes;
	}

	/* Filter locus for genotype bucket */
	TARGET void BUCKET::FilterLocusGT(BUCKET* obucket, bool progress, int nthread, MEMORY* mem, SLOCUS* nslocus, uint64* nlocus_pos)
	{
		//genotype id offset
		coffset = 0;
		expanding = false;

		uint* newid = new uint[nloc];
		memset(newid, 0xFF, nloc * sizeof(uint));

		//calculate new offsets
		for (int64 l = 0, nl = 0; l < nloc; ++l)
		{
			if (!GetLoc(l).flag_pass) continue;

			newid[l] = nl++;
			uint64 gtsize = obucket->offset[l].size;
			uint64 gtlinesize = (gtsize * nind + 7) >> 3;//xoffset 9
			offset.Push(OFFSET{ coffset, gtsize });
			coffset += gtlinesize;
		}

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
			uint64 gtlinesize = (gtsize * nind + 7) >> 3; //xoffset 11

			if (GetLoc(l).flag_pass)
			{
				int64 nl = newid[l];
				byte* write_pos = base_addr + offset[nl].offset, *write_end = write_pos + gtlinesize;

				AllocRange(write_pos, write_end);
				memmove(base_addr + offset[nl].offset, readpos, gtlinesize);

				//deep copy locus
				new(&nslocus[nl]) SLOCUS(mem[threadid], slocus[l]);
				if (haplotype) nlocus_pos[nl] = GetLocPos(l);
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

		delete[] newid;
		delete[] nbytes;
	}

	/* Filter locus for allele depth bucket */
	TARGET void BUCKET::FilterLocusAD(BUCKET* obucket, bool progress, int nthread)
	{
		//genotype id offset
		coffset = 0;
		expanding = false;

		uint* newid = new uint[nloc];
		memset(newid, 0xFF, nloc * sizeof(uint));

		//calculate new offsets
		for (int64 l = 0, nl = 0; l < nloc; ++l)
		{
			if (!GetLoc(l).flag_pass) continue;

			newid[l] = nl++;
			uint64 adsize = obucket->offset[l].size;
			uint64 adlinesize = (adsize * GetLoc(l).k * nind + 7) >> 3;//xoffset 10
			offset.Push(OFFSET{ coffset, adsize });
			coffset += adlinesize;
		}

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
			uint64 adlinesize = (adsize * GetLoc(l).k * nind + 7) >> 3; //xoffset 12

			if (GetLoc(l).flag_pass)
			{
				int64 nl = newid[l];
				byte* write_pos = base_addr + offset[nl].offset, *write_end = write_pos + adlinesize;

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

		delete[] newid;
		delete[] nbytes;
	}

	/* Allocate offset for genotype bucket */
	TARGET OFFSET BUCKET::AddOffsetGT(uint64 gsize, uint64 id)
	{
		OFFSET re = OFFSET{ 
			coffset.fetch_add((gsize * nind + 7) >> 3), 
			gsize };
		offset.Place(re, id, lock, expanding);

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
		offset.Place(re, id, lock, expanding);

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
		expanding = false;;
		offset = ref.offset;
		flag.swap(ref.flag);

		//unref vmem of ref
		ref.size = 0;
		ref.blocksize = 0;
		ref.base_addr = NULL;
		ref.tail_addr = NULL;
		ref.head_addr = NULL;
		ref.coffset = 0;
		ref.expanding = false;
	}

#endif

#ifndef _INCBUFFER

	/* Initialize, alloc 64 Kib */
	TARGET INCBUFFER::INCBUFFER()
	{
		len = INCREASE_BUFFER;
		data = (char*)malloc(len);
	}

	/* Uninitialize */
	TARGET INCBUFFER::~INCBUFFER()
	{
		if (data) free(data);
		data = NULL;
		len = 0;
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
		int64 nread = 0, rlen = len - offset, pos = FTell(handle);

		for (;;)
		{
			p2[0] = '\0';

			if (g_format_val < 0)
				gzgets((gzFile)handle, p2, (int)rlen);
			else
				fgets(p2, (int)rlen, handle);

			nread = FTell(handle) - pos;

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
	VLA_NEW(hThread, thread, nthreads + 2);
	VLA_NEW(fThread, bool, nthreads + 2);

	for (int i = 0; i < nthreads + 2; ++i)
	{
		fThread[i] = true;
		if (i == nthreads && GuardFunc1)
			hThread[i] = thread(GuardFunc1, i);
		else if (i == nthreads + 1 && GuardFunc2)
			hThread[i] = thread(GuardFunc2, i);
		else if (i < nthreads && Func)
			hThread[i] = thread(Func, i);
		else
			fThread[i] = false;
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
		if (fThread[i])
			hThread[i].join();

	VLA_DELETE(hThread);
	VLA_DELETE(fThread);
}

#pragma pack(pop)
