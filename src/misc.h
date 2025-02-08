/* Misc */

#include "vcfpop.h"

#pragma pack(push, 1)

struct CPOINT;
struct MEMBLOCK;
struct MEMORY;
template <typename T> struct LIST;
template <typename T, typename T2> struct TABLE_ENTRY;
template <typename T, typename T2> struct TABLE;
struct VMEMORY;
struct BUCKET;
struct INCBUFFER;
struct VCFBUFFER;

/* ReadWrite Lock */
struct ReadWriteLock
{
private:
	atomic<uint> state;

public:
	TARGET ReadWriteLock()
	{
		state = 0;
	}

	TARGET void lock_shared()
	{
		for (;;)
		{
			uint old_state = state.load();
			uint new_state = old_state + 1;

			if ((old_state >> 16) == 0 && // is not writing
				state.compare_exchange_weak(old_state, new_state)) //if state is equal to old_state
				break;
		}
	}

	TARGET void unlock_shared()
	{
		state--;
	}

	TARGET void lock()
	{
		for (;;)
		{
			uint old_state = state.load();
			uint new_state = old_state + 0x10000;

			//competite for write lock
			if ((old_state >> 16) == 0 && // if is not writing
				state.compare_exchange_weak(old_state, new_state)) //if state is equal to old_state
			{
				//wait until preiouvs read complete
				while (state & 0xFFFF);
				break;
			}
		}
	}

	TARGET void unlock() 
	{
		state -= 0x10000;
	}
};

/* Initialize a lock */
TARGET inline void InitLock(LOCK& x)
{
#ifdef LOCKFREE

#else 
	#ifdef _WIN64
		InitializeCriticalSection(&x);
	#else
        pthread_mutex_init(&x, NULL);
	#endif
#endif
}

/* Uninitialize a lock */
TARGET inline void UnInitLock(LOCK& x)
{
#ifdef LOCKFREE

#else
    #ifdef _WIN64
        DeleteCriticalSection(&x);
    #else
        pthread_mutex_destroy(&x);
    #endif
#endif
}

/* Lock a lock */
TARGET inline void Lock(LOCK& x)
{
#ifdef LOCKFREE
	while (x.test_and_set());
#else 
	#ifdef _WIN64
		EnterCriticalSection(&x);
	#else
		pthread_mutex_lock(&x);
	#endif
#endif
}

/* UnLock a lock */
TARGET inline void UnLock(LOCK& x)
{
#ifdef LOCKFREE
	x.clear();
#else 
	#ifdef _WIN64
		LeaveCriticalSection(&x);
	#else
		pthread_mutex_unlock(&x);
	#endif
#endif
}

/* Count number of alleles from frequency array*/
template<typename REAL>
TARGET int CountK(REAL* fre, int k2)
{
	int re = 0;
	for (int i = 0; i < k2; ++i)
		if (fre[i] > MIN_FREQ) re++;
	return re;
}

/* Set a bit in a variable */
TARGET void SetBit(byte& b, int pos, bool val);

/* Get a bit in a variable */
TARGET bool GetBit(byte b, int pos);

/* Get number of unique elements */
TARGET int GetNalleles(ushort* alleles, int ploidy);

/* Get current directory */
TARGET string GetCurDir();

/* Set current directory */
TARGET void SetCurDir(string dir);

/* Get absoulte path */
TARGET string GetAbsPath(string file);

/* Get parent path */
TARGET string GetParentPath(string file);

/* Clear temp files */
TARGET void ClearTempFiles(string& dir);

/* Pause console */
TARGET void Pause(void);

/* Start a timer */
TARGET void tic();

/* Show elapse time */
TARGET void toc();

/* Get time */
TARGET timepoint GetNow();

/* Get elapse time since begin */
TARGET double GetElapse(timepoint& begin);

/* Record time when begin evaluation */
TARGET void EvaluationBegin();

/* Record time when evaluation is finished and append results to the evaluation file */
TARGET void EvaluationEnd(const char* text);

/* Return memory usage */
TARGET int64 GetMemoryUsage();

/* Return memory usage string */
TARGET void GetUsageString(char* str, int64 val);

/* Write number of individuals, populations, loci, raw data size and memory usage */
TARGET void EvaluationStat();

/* Run a Rscript */
TARGET void RunRscript(string script);

/* Exit program with a message */
TARGET void Exit(const char* fmt, ...);

/* Point in a simplex to optimize by Down-Hill Simplex algorithm */	
struct CPOINT
{
	double unc_space[8];			//Unconstrained optimized parameter
	double real_space[9];			//Constrained real coefficients
	double lnL;						//Logarithm of likelihood
	int dim;						//Dimensions
	int neval;						//Different ploidy levels

	TARGET CPOINT(int de = 4);

	TARGET bool operator >(CPOINT& a);

	TARGET bool operator >=(CPOINT& a);

	TARGET bool operator <(CPOINT& a);

	TARGET bool operator <=(CPOINT& a);

	TARGET bool operator ==(CPOINT& a);

	TARGET bool operator !=(CPOINT& a);

	TARGET CPOINT operator =(const CPOINT& a);

	TARGET CPOINT operator +(const CPOINT& b);

	TARGET CPOINT operator -(const CPOINT& b);

	TARGET CPOINT operator *(const double b);

	TARGET CPOINT operator /(const double b);

	TARGET CPOINT& operator +=(CPOINT& a);

	TARGET CPOINT& operator -=(CPOINT& a);

	TARGET CPOINT& operator *=(double a);

	TARGET CPOINT& operator /=(double a);

	/* Distance in real space */
	TARGET double DistanceReal(CPOINT& a);

	/* Distance in unconstrained space */
	TARGET double DistanceUnc(CPOINT& a);

	/* Break iteration if distance between points are smaller than eps */
	TARGET bool static IsBreak(CPOINT* xx, int dim, double eps = 1e-7);

	/* Sort points */
	TARGET void static Order(CPOINT* xx, int dim);

	/* Gradient descent algorithm, need G and H */
	template<typename REAL>
	TARGET CPOINT static GradientDescent(void* Param, double (*func)(void* Param, CPOINT&, rmat&, rmat&), int dim, bool IsMinimize = false, double* Init = NULL);

	/* Downhill Simplex algorithm, do not need G and H*/
	template<typename REAL>
	TARGET CPOINT static DownHillSimplex(void* Param, double (*func)(void* Param, CPOINT&, rmat&, rmat&), int dim, double sep = 0.0001, int nrep = 15, bool IsMinimize = false);
};

/* Block of local memory management class */
struct MEMBLOCK
{
	bool* bucket;							//Data
	int size;								//Bucket size in bytes
	int used;								//Number of used bytes
};

/* Local memory management class for locus and genotype */
struct MEMORY
{
	MEMBLOCK* blocks;						//Blocks array
	int nblocks;							//Blocks size
	int cblock;								//Current block index
	LOCK lock; 								//Lock during expansion and allocation

	/* Initialize */
	TARGET MEMORY();

	/* Uninitialize */
	TARGET ~MEMORY();

	// Clear entries
	TARGET void ClearMemory();

	// Expand number of blocks
	TARGET void Expand();

	// Allocate a small piece of memory
	template <typename T>
	TARGET T* Alloc(T*& addr, int size, bool islock = true)
	{
		if (islock) Lock(lock);
		size *= sizeof(T);
		while (size + blocks[cblock].used > blocks[cblock].size)
		{
			if (++cblock == nblocks) Expand();
			if (!blocks[cblock].bucket)
			{
				blocks[cblock].size = std::max(std::min(blocks[cblock - 1].size << 1, BIG_FILE ? 256 * 1024 * 1024 : 8 * 1024 * 1024), size);
				blocks[cblock].bucket = new bool[blocks[cblock].size];
			}
		}
		addr = (T*)(blocks[cblock].bucket + blocks[cblock].used);
		blocks[cblock].used += size;
		if (islock) UnLock(lock);
		return addr;
	}

	// Allocate a small piece of memory
	TARGET byte* Alloc(int size);

	// Free a small piece of memory
	template <typename T>
	TARGET void Free(T* addr, int size, bool islock = true)
	{
		if (islock) Lock(lock);
		if (addr + size == (T*)(blocks[cblock].bucket + blocks[cblock].used))
			blocks[cblock].used -= size * sizeof(T);
		if (islock) UnLock(lock);
	}

	// ReAllocate a small piece of memory
	template <typename T>
	TARGET void ReAlloc(T*& addr, int size, T* oaddr, int osize)
	{
		Lock(lock);
		Free(oaddr, osize, false);
		Alloc(addr, size, false);
		if (oaddr != addr)
			SetVal(addr, oaddr, osize);
		UnLock(lock);
	}

	// Move from new mem to odd mem, free new mem, and increase old mem
	template <typename T>
	TARGET void Move(T*& oaddr, int osize, T* naddr, int nsize)
	{
		Lock(lock);

		// oaddr and naddr are adjacent and naddr is the newest
		if (oaddr + osize == naddr &&
			naddr + nsize == (T*)(blocks[cblock].bucket + blocks[cblock].used))
		{
			SetVal(oaddr, naddr, nsize); //move mem
			blocks[cblock].used -= osize * sizeof(T);
		}
		else
			oaddr = naddr;

		UnLock(lock);
	}
};

/* Fast list class */
template <typename T>
struct LIST
{
	T* bucket;							//Bucket
	uint bucket_size;					//Size of bucket
	uint size;							//Number of entries
	MEMORY* memory;						//Memory management class

	/* Deep copy a list */
	TARGET LIST<T>& operator=(LIST<T>& ref)
	{
		if ((uint64)memory == 0xFFFFFFFFFFFFFFFF) memory = NULL;
		if ((uint64)bucket == 0xFFFFFFFFFFFFFFFF) bucket = NULL;

		if (memory && bucket)
			memory->Free(bucket, bucket_size);

		if (!memory && bucket)
			free(bucket);

		bucket = NULL;
		bucket_size = ref.bucket_size;
		size = ref.size;
		memory = ref.memory;

		if (memory)
			bucket = ref.bucket;
		else
			bucket = (T*)malloc(bucket_size * sizeof(T));

		SetVal(bucket, ref.bucket, bucket_size);
		return *this;
	}

	/* Set Zero */
	/*TARGET LIST()
	{
		SetZero(this, 1);
	}
	*/

	/* Initialize */
	TARGET LIST(MEMORY* _memory = NULL)
	{
		memory = _memory;
		bucket_size = 0x8;
		size = 0;
		bucket = memory ? (T*)memory->Alloc(bucket_size * sizeof(T)) : (T*)malloc(bucket_size * sizeof(T));
		SetFF(bucket, bucket_size);
	}

	/* Clear entries but do not shrink bucket */
	TARGET void Clear()
	{
		SetFF(bucket, bucket_size);
		size = 0;
	}

	/* Uninitialize and free memory */
	TARGET ~LIST()
	{
		if (memory && bucket)
			memory->Free(bucket, bucket_size);

		if (!memory && bucket) 
			free(bucket);

		bucket = NULL;
		memory = NULL;
		size = 0;
		bucket_size = 0;
	}

	/* Expand list */
	TARGET void Expand()
	{
		uint nsize = std::max(bucket_size * 2, 4u);
		T* nbucket = bucket;

		if (memory)
		{
			memory->ReAlloc(nbucket, nsize, bucket, bucket_size);
			SetFF(nbucket + bucket_size, nsize - bucket_size);
			bucket = nbucket;
			bucket_size = nsize;
		}
		else
		{
			nbucket = (T*)malloc(nsize * sizeof(T));
			memmove(nbucket, bucket, bucket_size * sizeof(T));
			SetFF(nbucket + bucket_size, nsize - bucket_size);
			T* obucket = bucket;
			bucket = nbucket;
			bucket_size = nsize;
			free(obucket);
		}
	}

	/* Access by index */
	TARGET T& operator[](uint64 i)
	{
		if (i >= size)
			Exit("\nError: LIST class access an unexist element.\n");
		return bucket[i];
	}

	/* Add an entry to the end */
	TARGET void Push(T& val)
	{
		if (size >= bucket_size) Expand();
		bucket[size++] = val;
	}

	/* Add an entry to the end */
	TARGET void Push(T&& val)
	{
		if (size >= bucket_size) Expand();
		bucket[size++] = val;
	}

	/* Add an entry at a place id */
	TARGET void Place(T& val, uint id, ReadWriteLock& rwlock)
	{
		if (id >= bucket_size)
		{
			rwlock.lock();
			if (id >= bucket_size)
				Expand();
			rwlock.unlock();
		}

		rwlock.lock_shared(); 
		AtomicMax(*(atomic<uint>*) & size, id + 1);
		bucket[id] = val;
		rwlock.unlock_shared();
	}

	/* Return and delete the entry at the end */
	TARGET T Pop()
	{
		return bucket[--size];
	}

	/* Remove an entry by index */
	TARGET void Erase(uint i)
	{
		if (i >= bucket_size) Exit("\nError: LIST class erase an unexist element.\n");
		SetVal(bucket + i, bucket + i + 1, size - i - 1);
		SetFF(bucket + size - 1, 1);
		size--;
	}

	/* Alloc memory */
	TARGET void SetSize(uint nsize)
	{
		nsize = std::max(bucket_size, nsize);
		T* nbucket = bucket;

		if (memory)
		{
			memory->ReAlloc(nbucket, nsize, bucket, bucket_size);
			SetFF(nbucket + bucket_size, nsize - bucket_size);
			bucket = nbucket;
			bucket_size = nsize;
		}
		else
		{
			nbucket = (T*)malloc(nsize * sizeof(T));
			memmove(nbucket, bucket, bucket_size * sizeof(T));
			SetFF(nbucket + bucket_size, nsize - bucket_size);
			T* obucket = bucket;
			bucket = nbucket;
			bucket_size = nsize;
			free(obucket);
		}
	}

	/* Set _nsize empty elements */
	TARGET void SetEmpty(uint _nsize)
	{
		SetSize(_nsize);
		SetZero(bucket, _nsize);
		size = _nsize;
	}
};

/* Hash table entry: key-val pair */
template <typename T, typename T2>
struct TABLE_ENTRY
{
	T key;								//Hash key
	T2 val;								//Table value
};

/* Fast hash table class */
template <typename T, typename T2>
struct TABLE
{
	TABLE_ENTRY<T,T2>* bucket;		//Key-val pair bucket
	uint* index;					//The list saves the position in bucket of ith entry
	uint mask;						//mask + 1 is the bucket size
	uint size;						//Number of entries
	MEMORY* memory;					//Memory management class

	/* Set Zero */
	TARGET TABLE()
	{
		SetZero(this, 1);
	}

	/* Initialize a table */
	TARGET TABLE(bool haslist, MEMORY* _memory, uint _size = 0)
	{
		memory = _memory;
		if (_size <= 4)
			mask = 3;
		else
			mask = (1 << CeilLog2((int)(_size + (_size >> 1)))) - 1;
		
		size = 0;
		byte* buf = memory ? 
			memory->Alloc((sizeof(TABLE_ENTRY<T, T2>) + (haslist ? sizeof(uint) : 0)) * (mask + 1)) : 
			(byte*)malloc((sizeof(TABLE_ENTRY<T, T2>) + (haslist ? sizeof(uint) : 0)) * (mask + 1));

		bucket = (TABLE_ENTRY<T, T2>*)buf;
		index = (uint*)(haslist ? bucket + (mask + 1) : NULL);
		Clear();
	}

	/* Clear entries but do not release memory */
	TARGET void Clear()
	{
		SetFF((byte*)bucket, (sizeof(TABLE_ENTRY<T, T2>) + (index ? sizeof(uint) : 0)) * (mask + 1));
		size = 0;
	}

	/* Free memory */
	TARGET ~TABLE()
	{
		if (!memory && bucket) 
			free(bucket);

		if (memory && bucket)
			memory->Free((byte*)bucket, (sizeof(TABLE_ENTRY<T, T2>) + (index ? sizeof(uint) : 0)) * (mask + 1));

		bucket = NULL;
		index = NULL;
		mask = 0;
		size = 0;
		memory = NULL;
	}

	/* Expand table */
	TARGET void Expand()
	{
		uint mask2 = ((mask + 1) << 1) - 1;
		uint olen = (sizeof(TABLE_ENTRY<T, T2>) + (index ? sizeof(uint) : 0)) * (mask + 1);
		uint nlen = (sizeof(TABLE_ENTRY<T, T2>) + (index ? sizeof(uint) : 0)) * (mask2 + 1);

		byte*obuf = (byte*)bucket;
		byte* nbuf = memory ? memory->Alloc(nlen) : (byte*)malloc(nlen);

		SetFF(nbuf, nlen);
		TABLE_ENTRY<T, T2>* nbucket = (TABLE_ENTRY<T, T2>*)nbuf;
		uint* nindex = (uint*)(index ? nbucket + (mask2 + 1) : 0);

		//move to new bucket
		if (index)
		{
			//get an entry from old table
			for (uint i = 0; i < size; ++i)
			{
				//find an empty emtry in the new table
				uint st = bucket[index[i]].key & mask2;
				for (uint j = 0; j <= mask2; ++j)
					if (nbucket[(j + st) & mask2].key == (T)-1)
					{
						//insert
						nbucket[nindex[i] = ((j + st) & mask2)] = bucket[index[i]];
						break;
					}
			}
		}
		else for (uint i = 0; i <= mask; ++i)
		{
			//get an entry from old table
			T key = bucket[i].key;
			if (key == (T)-1) continue;

			//find an empty slot in the new table
			uint st = key & mask2;
			for (uint j = 0; j <= mask2; ++j)
				if (nbucket[(j + st) & mask2].key == (T)-1)
				{
					//insert
					nbucket[(j + st) & mask2] = bucket[i];
					break;
				}
		}
		if (memory)
		{
			memory->Move(obuf, olen, nbuf, nlen);
			bucket = (TABLE_ENTRY<T, T2>*)obuf;
		}
		else
		{
			free(bucket);
			bucket = nbucket;
		}
		index = index ? (uint*)(bucket + (mask2 + 1)) : NULL;
		mask = mask2;
	}

	/* Get array [k] with k allele sizes for SMM distance model
	TARGET void GetLength(T2*& alen, MEMORY& locus_mem)
	{
		//For non-vcf input, allele identifier is the size
		//Return an array [size]
		//Table is the allele length table, key:value = allele:length

		locus_mem.Alloc(alen, size);
		for (uint i = 0; i < size; ++i)
			alen[i] = bucket[index[i]].key;
	} */

	/* Return value if entry exists, otherwise return null */
	TARGET T2 Try(T key)
	{
		uint st = key & mask;
		for (uint i = 0; i <= mask; ++i)
		{
			uint idx = (i + st) & mask;
			T k = bucket[idx].key;
			if (k == (T)-1) return (T2)NULL;
			if (k == key) return bucket[idx].val;
		}
		return (T2)NULL;
	}

	/* Use table as list where value is id, insert if entry doesn't exists, and return id */
	TARGET T2 PushIndex(T key)
	{
	restart:
		uint st = key & mask;
		for (uint i = 0; i <= mask; ++i)
		{
			uint idx = (i + st) & mask;
			T k = bucket[idx].key;
			if (k == (T)-1)
			{
				if (mask != 3 && (uint)(size + (size >> 1)) > mask)
				{
					Expand();
					goto restart;
				}
				bucket[idx].key = key;
				bucket[idx].val = (T2)size;
				if (index) index[size] = idx;
				return (T2)size++;
			}
			if (k == key) return bucket[idx].val;
		}

		if (mask == 3 && size == 4)
		{
			Expand();
			goto restart;
		}

		return (T2)-1;
	}

	/* Has an entry with the key */
	TARGET bool ContainsKey(T key)
	{
		uint st = key & mask;
		for (uint i = 0; i <= mask; ++i)
		{
			uint idx = (i + st) & mask;
			T k = bucket[idx].key;
			if (k == (T)-1) return false;
			if (k == key) return true;
		}
		return false;
	}

	/* Add one to the entry, if not exists then create this entry */
	TARGET void AddCount(T key)
	{
	restart:
		//search from hash % (mask + 1)
		uint st = key & mask;

		for (uint i = 0; i <= mask; i++)
		{
			uint id1 = (i + st + 0) & mask;
			T k1 = bucket[id1].key;

			if (k1 == (T)-1)
			{
				//insert a new entry
				if (mask != 3 && (uint)(size + (size >> 1)) > mask)
				{
					Expand();
					goto restart;
				}
				bucket[id1].key = key;
				bucket[id1].val = 1;
				if (index) index[size] = id1;
				size++;
				return;
			}
			if (k1 == key)
			{
				bucket[id1].val++;
				return;
			}
		}

		if (mask == 3 && size == 4)
		{
			Expand();
			goto restart;
		}

		Exit("\nError: key not found.\n");
	}

	/* Access entry by index */
	TARGET T2& operator()(uint64 i)
	{
		if (i >= size)
			Exit("\nError: index exceed limit in interal LIST class.\n");
		return bucket[index[i]].val;
	}

	/* Access entry by index */
	TARGET TABLE_ENTRY<T,T2>& GetEntry(uint64 i)
	{
		if (i >= size)
			Exit("\nError: index exceed limit in interal LIST class.\n");
		return bucket[index[i]];
	}

	/* Access entry by key, single thread, not need to lock */
	TARGET T2& operator[](T key)
	{
	restart:
		//search from hash % (mask + 1)
		uint st = key & mask;

		for (uint i = 0; i <= mask; i++)
		{
			uint id1 = (i + st + 0) & mask;
			T k1 = bucket[id1].key;

			if (k1 == (T)-1)
			{
				//insert a new entry
				if (mask != 3 && (uint)(size + (size >> 1)) > mask)
				{
					Expand();
					goto restart;
				}
				bucket[id1].key = key;
				if (index) index[size] = id1;
				size++;
				return bucket[id1].val;
			}
			if (k1 == key)
                return bucket[id1].val;
		}

		if (mask == 3 && size == 4)
		{
			Expand();
			goto restart;
		}

		Exit("\nError: key not found.\n");
		return bucket[0].val;
	}
};

#pragma pack(pop)

/* Virtual memory management class */
struct alignas(8) VMEMORY
{
	int64 size;								//number of bytes in this virtual memory
	int64 blocksize;						//number of bytes in each block estimated data size / 800
	byte* base_addr;						//first alloced address
	byte* tail_addr;						//next alloced address, base_addr + nallocated * blocksize
	byte* head_addr;						//next removed address, base_addr + nremoved * blocksize
	vector<byte> flag;						//for discrete allocation
	LOCK lock;
    
	/* Do nothing */
	TARGET VMEMORY();

	/* Assign value but do not alloc */
	TARGET VMEMORY(int64 size);

	/* Do nothing */
	TARGET ~VMEMORY();

	/* Alloc or Clear a Virtual Memory slot, each have 16 TiB */
	TARGET void ApplyClearSlot(bool isapply);

	/* Alloc a blocksize memory at the end */
	TARGET void Alloc(byte* ed);

	/* Alloc a blocksize memory at the end */
	TARGET void AllocRange(byte* st, byte* ed);

	/* Unalloc all memory */
	TARGET void UnAllocAll();

	/* If the data size in a block is zero, call UnAllocAddr */
	TARGET void SubUnAlloc(atomic<int>& _size, int subv, int blockid);

	/* Unalloc a blocksize memory */
	TARGET void UnAllocAddr(byte* addr);
};

/* Genotype and allelic depth bucket */
struct alignas(8) BUCKET : public VMEMORY
{
    atomic<uint64> coffset;							//Size of used bucket memory
	LIST<OFFSET> offset;							//Genotype index data at each locus
	ReadWriteLock rwlock;
	//bool expanding;

	/* Do nothing */
	TARGET BUCKET();

	/* Release memory */
	TARGET ~BUCKET();

	/* Create bucket for genotype or allelic depth bucket */
	TARGET void CreateBucket();

	/* Create bucket for genotype bucket in haplotype extraction */
	TARGET void CreateBucketGT(LIST<HAPLO_DUMMY_LOCUS>& hlocus);

	/* Create genotype table for unphase genotypes */
	TARGET void CreateBucketGT(SLOCUS* locus);

	/* Create genotype table for non-vcf/bcf input files */
	TARGET void CreateBucketGT(LOCUS* locus);

	/* Filter individual for genotype bucket */
	template<typename REAL>
	TARGET void FilterIndividualGT(BUCKET* obucket, int nfilter, bool progress, int nthread);

	/* Filter individual for allele depth bucket */
	template<typename REAL>
	TARGET void FilterIndividualAD(BUCKET* obucket, int nfilter, bool progress, int nthread);

	/* Filter locus for genotype bucket */
	TARGET void FilterLocusGT(BUCKET* obucket, bool progress, int nthread, MEMORY* mem, SLOCUS* nslocus, uint64* nlocus_pos);

	/* Filter locus for allele depth bucket */
	TARGET void FilterLocusAD(BUCKET* obucket, bool progress, int nthread);

	/* Allocate offset for genotype bucket */
	TARGET OFFSET AddOffsetGT(uint64 gsize, uint64 id);

	/* Allocate offset for allele depth bucket */
	TARGET OFFSET AddOffsetAD(uint64 adsize, uint64 id, uint k);

	/* Use another bucket to replace this */
	TARGET void Replace(BUCKET& ref);
};

/* Buffer with increaseing size */
struct INCBUFFER
{
	char* data;										//buffer data
	int64 len;										//buffer size
	char* derived;

	/* Initialize, alloc 64 Kib */
	TARGET INCBUFFER();

	/* Uninitialize */
	TARGET ~INCBUFFER();

	/* Expand and memmove */
	TARGET void Move(char* src, int nlen);

	/* Read a line from a file */
	TARGET int64 Gets(FILE* handle);

	/* Read a line from a file */
	TARGET int64 Gets(FILE* handle, int maxlen, char*& p2);

	/* Read a line from a file */
	TARGET int64 Gets(FILE* handle, char*& p2);

	/* Expand buffer if nlen is greater than len*/
	TARGET void Expand(int nlen);
};

struct VCFBUFFER : public INCBUFFER
{
	int64 geno_len;									//number of bytes to read in genotype
	char* line_pos;									//pointer of current line
	char* geno_pos;									//pointer of genotype
	int64 line_len;									//length of current line = read_len
	int64 data_size;								//number of bytes readed in this buffer
	int64 jloc;
	int64 file_offset;								//previous offset of file pointer, for progress bar

	/* Read data from file */
	TARGET int64 Read1In(FILE* handle);

	/* Read data from file, expand buffer if necessary, and set progress */
	TARGET void Read1(FILE* handle);

	/* Read data from file, expand buffer if necessary, and set progress */
	TARGET void Read2(FILE* handle);
};

/* Initialize double reduction rates */
TARGET void InitAlpha();

/* Get aligned 32 bytes address */
template<typename T>
TARGET T Align64(T addr)
{
	return (T)(((uint64)addr + 63) & 0xFFFFFFFFFFFFFFC0);
}

/* Get aligned 32 bytes address */
template<typename T>
TARGET T Align32(T addr)
{
	return (T)(((uint64)addr + 31) & 0xFFFFFFFFFFFFFFE0);
}

/* Get aligned 16 bytes address */
template<typename T>
TARGET T Align16(T addr)
{
	return (T)(((uint64)addr + 15) & 0xFFFFFFFFFFFFFFF0);
}

/* Get aligned 16 bytes address */
template<typename T>
TARGET T Align8(T addr)
{
	return (T)(((uint64)addr + 7) & 0xFFFFFFFFFFFFFFF8);
}

/* Get aligned 16 bytes address */
template<typename T>
TARGET T Align(T addr, int nalign)
{
	return (T)(((uint64)addr + nalign - 1) & (~(nalign - 1)));
}

/* Calculate task in multiple threads */
TARGET void RunThreads(void (*Func) (int), void (*GuardFunc1) (int), void (*GuardFunc2) (int),
	int64 ntot, int64 nct, const char* info, int nthreads, bool isfirst, int nprogress = g_progress_val);

