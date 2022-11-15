/* Hash Functions */

#include "vcfpop.h"

template TARGET void HashHaplotype<double>(IND<double>* ti, int64 st, int64 ed, HASH* hash, int& ploidy);
template TARGET void HashHaplotype<float >(IND<float >* ti, int64 st, int64 ed, HASH* hash, int& ploidy);

/* Get crc32 for 64 bits data*/
uint crc32_u64(uint prev, uint64 val)
{
#ifdef __aarch64__
	return __crc32cd(prev, val);
#else
	return (uint)_mm_crc32_u64((uint)prev, val);
#endif
}

#ifdef HASH64

/* 64-bit hash to reduce hash collision */

/* Get crc32 value for a genotype */
TARGET HASH HashGenotype(ushort* allele, int ploidy)
{
	uint crc1 = crc32_u64(0xB385FDA9, ploidy), crc2 = crc32_u64(0x7FED5EA6, ploidy);
	switch (ploidy)
	{
	case 1:
		crc1 = crc32_u64(crc1, (uint64)*(ushort*)allele);
		crc2 = crc32_u64(crc2, (uint64)*(ushort*)allele);
		break;
	case 2:
		crc1 = crc32_u64(crc1, (uint64)*(uint*  )allele);
		crc2 = crc32_u64(crc2, (uint64)*(uint*  )allele);
		break;
	case 3:
		crc1 = crc32_u64(crc1, (uint64)*(uint*  )allele | ((uint64)allele[2] << 32));
		crc2 = crc32_u64(crc2, (uint64)*(uint*  )allele | ((uint64)allele[2] << 32));
		break;
	case 4:
		crc1 = crc32_u64(crc1, (uint64)*(uint64*)allele);
		crc2 = crc32_u64(crc2, (uint64)*(uint64*)allele);
		break;
	case 5:
		crc1 = crc32_u64(crc32_u64(crc1, (uint64)*(uint64*)allele), (uint64)*(ushort*)(allele + 4));
		crc2 = crc32_u64(crc32_u64(crc2, (uint64)*(uint64*)allele), (uint64)*(ushort*)(allele + 4));
		break;
	case 6:
		crc1 = crc32_u64(crc32_u64(crc1, (uint64)*(uint64*)allele), (uint64)*(uint*  )(allele + 4));
		crc2 = crc32_u64(crc32_u64(crc2, (uint64)*(uint64*)allele), (uint64)*(uint*  )(allele + 4));
		break;
	case 7:
		crc1 = crc32_u64(crc32_u64(crc1, (uint64)*(uint64*)allele), (uint64)*(uint*  )(allele + 4) | ((uint64)allele[6] << 32));
		crc2 = crc32_u64(crc32_u64(crc2, (uint64)*(uint64*)allele), (uint64)*(uint*  )(allele + 4) | ((uint64)allele[6] << 32));
		break;
	case 8:
		crc1 = crc32_u64(crc32_u64(crc1, (uint64)*(uint64*)allele), (uint64)*(uint64*)(allele + 4));
		crc2 = crc32_u64(crc32_u64(crc2, (uint64)*(uint64*)allele), (uint64)*(uint64*)(allele + 4));
		break;
	case 9:
		crc1 = crc32_u64(crc32_u64(crc32_u64(crc1, (uint64)*(uint64*)allele), (uint64)*(uint64*)(allele + 4)), (uint64)*(ushort*)(allele + 8));
		crc2 = crc32_u64(crc32_u64(crc32_u64(crc2, (uint64)*(uint64*)allele), (uint64)*(uint64*)(allele + 4)), (uint64)*(ushort*)(allele + 8));
		break;
	case 10:
		crc1 = crc32_u64(crc32_u64(crc32_u64(crc1, (uint64)*(uint64*)allele), (uint64)*(uint64*)(allele + 4)), (uint64)*(uint*)(allele + 8));
		crc2 = crc32_u64(crc32_u64(crc32_u64(crc2, (uint64)*(uint64*)allele), (uint64)*(uint64*)(allele + 4)), (uint64)*(uint*)(allele + 8));
		break;
	default: break;
	}
	return ((uint64)crc1 << 32) | (uint64)crc2;
}

/* Get crc32 value for a string */
TARGET HASH HashString(char* str, int len)
{
	uint crc1 = 0xB385FDA9, crc2 = 0x7FED5EA6;
	if (len == -1) len = (int)strlen(str);
	for (; len >= 8; len -= 8, str += 8)
	{
		crc1 = crc32_u64(crc1, *(uint64*)str);
		crc2 = crc32_u64(crc2, *(uint64*)str);
	}

	switch (len)
	{
	case 1:
		crc1 = crc32_u64(crc1, (uint64)*(byte  *)str);
		crc2 = crc32_u64(crc2, (uint64)*(byte  *)str);
		break;
	case 2:
		crc1 = crc32_u64(crc1, (uint64)*(ushort*)str);
		crc2 = crc32_u64(crc2, (uint64)*(ushort*)str);
		break;
	case 3:
		crc1 = crc32_u64(crc1, (uint64)*(ushort*)str | ((uint64)*(byte*  )(str + 2) << 16));
		crc2 = crc32_u64(crc2, (uint64)*(ushort*)str | ((uint64)*(byte*  )(str + 2) << 16));
		break;
	case 4:
		crc1 = crc32_u64(crc1, (uint64)*(uint*  )str);
		crc2 = crc32_u64(crc2, (uint64)*(uint*  )str);
		break;
	case 5:
		crc1 = crc32_u64(crc1, (uint64)*(uint*  )str | ((uint64)*(byte*  )(str + 4) << 32));
		crc2 = crc32_u64(crc2, (uint64)*(uint*  )str | ((uint64)*(byte*  )(str + 4) << 32));
		break;
	case 6:
		crc1 = crc32_u64(crc1, (uint64)*(uint*  )str | ((uint64)*(ushort*)(str + 4) << 32));
		crc2 = crc32_u64(crc2, (uint64)*(uint*  )str | ((uint64)*(ushort*)(str + 4) << 32));
		break;
	case 7:
		crc1 = crc32_u64(crc1, (uint64)*(uint*  )str | (((uint64)*(uint* )(str + 3) >> 8) << 32));
		crc2 = crc32_u64(crc2, (uint64)*(uint*  )str | (((uint64)*(uint* )(str + 3) >> 8) << 32));
		break;
	default:
		break;
	}
	return ((uint64)crc1 << 32) | (uint64)crc2;
}

/* Get crc32 value for unsigned long integer */
TARGET HASH HashULong(uint64 val)
{
	return ((uint64)crc32_u64(0xB385FDA9, val) << 32) | (uint64)crc32_u64(0x7FED5EA6, val);
}

/* Get crc32 value for two unsigned integer */
TARGET HASH HashUInt(uint val, uint val2)
{
	return HashULong(((uint64)val << 32) | (uint64)val2);
}

#else

/* 32-bit hash */

/* Get crc32 value for a genotype */
TARGET HASH HashGenotype(ushort* allele, int ploidy)
{
	HASH crc1 = crc32_u64(0xB385FDA9, (uint64)ploidy);
	switch (ploidy)
	{
	case 1:
		crc1 = crc32_u64(crc1, (uint64)*(ushort*)allele);
		break;
	case 2:
		crc1 = crc32_u64(crc1, (uint64)*(uint*  )allele);
		break;
	case 3:
		crc1 = crc32_u64(crc1, (uint64)*(uint*  )allele | ((uint64)allele[2] << 32));
		break;
	case 4:
		crc1 = crc32_u64(crc1, (uint64)*(uint64*)allele);
		break;
	case 5:
		crc1 = crc32_u64(crc32_u64(crc1, (uint64)*(uint64*)allele), (uint64)*(ushort*)(allele + 4));
		break;
	case 6:
		crc1 = crc32_u64(crc32_u64(crc1, (uint64)*(uint64*)allele), (uint64)*(uint*  )(allele + 4));
		break;
	case 7:
		crc1 = crc32_u64(crc32_u64(crc1, (uint64)*(uint64*)allele), (uint64)*(uint*  )(allele + 4) | ((uint64)allele[6] << 32));
		break;
	case 8:
		crc1 = crc32_u64(crc32_u64(crc1, (uint64)*(uint64*)allele), (uint64)*(uint64*)(allele + 4));
		break;
	case 9:
		crc1 = crc32_u64(crc32_u64(crc32_u64(crc1, (uint64)*(uint64*)allele), (uint64)*(uint64*)(allele + 4)), (uint64)*(ushort*)(allele + 8));
		break;
	case 10:
		crc1 = crc32_u64(crc32_u64(crc32_u64(crc1, (uint64)*(uint64*)allele), (uint64)*(uint64*)(allele + 4)), (uint64)*(uint*  )(allele + 8));
		break;
	default: break;
	}
	return (uint)crc1;
}

/* Get crc32 value for a string */
TARGET HASH HashString(char* str, int len)
{
	HASH crc1 = 0xB385FDA9;
	if (len == -1) len = (int)strlen(str);

	for (; len >= 8; len -= 8, str += 8)
		crc1 = crc32_u64(crc1, *(uint64*)str);

	switch (len)
	{
	case 1:
		crc1 = crc32_u64(crc1, (uint64)*(byte  *)str);
		break;
	case 2:
		crc1 = crc32_u64(crc1, (uint64)*(ushort*)str);
		break;
	case 3:
		crc1 = crc32_u64(crc1, (uint64)*(ushort*)str | ((uint64)*(byte*  )(str + 2) << 16));
		break;
	case 4:
		crc1 = crc32_u64(crc1, (uint64)*(uint*  )str);
		break;
	case 5:
		crc1 = crc32_u64(crc1, (uint64)*(uint*  )str | ((uint64)*(byte*  )(str + 4) << 32));
		break;
	case 6:
		crc1 = crc32_u64(crc1, (uint64)*(uint*  )str | ((uint64)*(ushort*)(str + 4) << 32));
		break;
	case 7:
		crc1 = crc32_u64(crc1, (uint64)*(uint*  )str | (((uint64)*(uint* )(str + 3) >> 8) << 32));
		break;
	default:
		break;
	}
	return (uint)crc1;
}

/* Get crc32 value for unsigned long integer */
TARGET HASH HashULong(uint64 val)
{
	return crc32_u64(0xB385FDA9, val);
}

/* Get crc32 value for unsigned long integer */
TARGET uint64 Hash64ULong(uint64 val)
{
	return ((uint64)crc32_u64(0xB385FDA9, val) << 32) | (uint64)crc32_u64(0x159A55E5, val);
}

/* Get crc32 value for two unsigned integer */
TARGET HASH HashUInt(uint val, uint val2)
{
	return HashULong(((uint64)val << 32) | (uint64)val2);
}

#endif

/* Initialize crypt table for Huang 2015 maximum-likelihood polyploid relatedness estimator */
TARGET void InitCryptTable()
{
	if (cryptTable != NULL)  return;
	cryptTable = new uint[0x500];

	uint dwHih, dwLow, seed = 0x00100001, index1 = 0, index2 = 0, i;
	for (index1 = 0; index1 < 0x100; ++index1)
	{
		for (index2 = index1, i = 0; i < 5; ++i, index2 += 0x100)
		{
			seed = (seed * 125 + 3) % 0x2AAAAB;
			dwHih = (seed & 0xFFFF) << 16;
			seed = (seed * 125 + 3) % 0x2AAAAB;
			dwLow = (seed & 0xFFFF);
			cryptTable[index2] = dwHih | dwLow;
		}
	}
}

/* Get IBS model of a genotype */
TARGET int GetSingleIBS(int* x, int ploidy)
{
	int n = 0, t[8], tt = 0;
	for (int i = 0; i < ploidy; )
	{
		int j = i + 1;
		for (; j < ploidy && x[j] == x[i]; ++j);
		t[n++] = j - i;
		i = j;
	}

	for (int i = 0; i < n; ++i)
		for (int j = i + 1; j < n; ++j)
			if (t[j] < t[i])
				Swap(t[i], t[j]);

	for (int i = 0; i < n; ++i)
		tt = tt * 10 + t[i];

	return tt;
}

/* Get hash of IBS mode of a pair of genotypes */
TARGET uint HashString32(char* s1, char* s2, char* s3)
{
	uint dwSeed1 = 0x7FED7FED, dwSeed2 = 0xEEEEEEEE;
	uint* cryptTable2 = (uint*)cryptTable;
	byte b1;
	for (int i = 0; s1[i]; ++i)
	{
		b1 = (byte)s1[i];
		dwSeed1 = cryptTable2[b1 + 0x100] ^ (dwSeed1 + dwSeed2);
		dwSeed2 = b1 + dwSeed1 + dwSeed2 + (dwSeed2 << 5) + 3;
	}
	for (int i = 0; s2[i]; ++i)
	{
		b1 = (byte)s2[i];
		dwSeed1 = cryptTable2[b1 + 0x100] ^ (dwSeed1 + dwSeed2);
		dwSeed2 = b1 + dwSeed1 + dwSeed2 + (dwSeed2 << 5) + 3;
	}
	for (int i = 0; s3[i]; ++i)
	{
		b1 = (byte)s3[i];
		dwSeed1 = cryptTable2[b1 + 0x100] ^ (dwSeed1 + dwSeed2);
		dwSeed2 = b1 + dwSeed1 + dwSeed2 + (dwSeed2 << 5) + 3;
	}
	return dwSeed1;
}

/* Get IBS mode of a pair of genotypes */
TARGET uint GetHuang2015Hash(int* x, int* y, int p)
{
	//xibs + yibs + crosssame + crossfactors
	int crosssame = 0;
	int crfa[8], crfb[8], crfan = 0, crfbn = 0;
	for (int i = 0; i < p; ++i)
	{
		int c1 = 0, c2 = 0;
		for (int j = 0; j < p; ++j)
		{
			if (x[i] == y[j])
			{
				crosssame++;
				c1++;
			}
			if (x[j] == y[i])
				c2++;
		}
		if (c1) crfa[crfan++] = c1;
		if (c2) crfb[crfbn++] = c2;
	}

	for (int i = 0; i < crfan; ++i)
		for (int j = i + 1; j < crfan; ++j)
			if (crfa[i] > crfa[j])
				Swap(crfa[i], crfa[j]);

	for (int i = 0; i < crfbn; ++i)
		for (int j = i + 1; j < crfbn; ++j)
			if (crfb[i] > crfb[j])
				Swap(crfb[i], crfb[j]);

	char s1[20], s2[20], s3[20];
	s2[0] = '\0';
	s3[0] = '\0';

	char* ss1 = s1;
	AppendInt(ss1, GetSingleIBS(x, p)); *ss1++ = ',';
	AppendInt(ss1, GetSingleIBS(y, p)); *ss1++ = ',';
	AppendInt(ss1, crosssame); *ss1++ = '\0';

	for (int i = 0; i < crfan; ++i)
	{
		s2[i * 2] = (char)('0' + crfa[i]);
		s2[i * 2 + 1] = ',';
	}
	if (crfan) s2[crfan * 2 - 1] = 0;
	for (int i = 0; i < crfbn; ++i)
	{
		s3[i * 2] = (char)('0' + crfb[i]);
		s3[i * 2 + 1] = ',';
	}
	if (crfbn) s3[crfbn * 2 - 1] = 0;
	return HashString32(s1, s2, s3);
}

/* Get hash of a pair of genotype id */
TARGET HASH HashDyadGenotypeIndex(HASH ha)
{
	return crc32_u64(0xB385FDA9, ha);
}

/* Get hash of a haplotype */
template<typename REAL>
TARGET void HashHaplotype(IND<REAL>* ti, int64 st, int64 ed, HASH* hash, int& ploidy)
{
	//allele should be sorted
#ifdef HASH64
	HASH Seed1[N_MAX_PLOIDY] = { 0xB385FDA97FED5EA6ull, 0xB385FDA97FED5EA6ull, 0xB385FDA97FED5EA6ull, 0xB385FDA97FED5EA6ull, 0xB385FDA97FED5EA6ull, 0xB385FDA97FED5EA6ull, 0xB385FDA97FED5EA6ull, 0xB385FDA97FED5EA6ull, 0xB385FDA97FED5EA6ull, 0xB385FDA97FED5EA6ull };
	uint* Seed2 = (uint*)Seed1;
#else
	HASH Seed1[N_MAX_PLOIDY] = { 0xB385FDA9, 0xB385FDA9, 0xB385FDA9, 0xB385FDA9, 0xB385FDA9, 0xB385FDA9, 0xB385FDA9, 0xB385FDA9, 0xB385FDA9, 0xB385FDA9 };
#endif
	uint64 alleles[N_MAX_PLOIDY] = { 0 };

	for (int64 l = st, hc = 0; l <= ed; ++l)
	{
		//load four alleles and hash
		GENOTYPE& gt = ti->GetGenotype(GetLocId(l), GetLoc(l).GetGtab());
		if (hc == 0) ploidy = gt.Ploidy();

		if (gt.Nalleles() == 0)
		{
			SetFF(hash, ploidy);
			return;
		}

		ushort* als = gt.GetAlleleArray();
		for (int vi = 0; vi < ploidy; ++vi)
			alleles[vi] = alleles[vi] << 16 | als[vi];

		if (++hc % 4 == 0)
		{
			for (int vi = 0; vi < ploidy; ++vi)
			{
#ifdef HASH64
				Seed2[(vi << 1)]     = crc32_u64(Seed2[(vi << 1)    ], alleles[vi]);
				Seed2[(vi << 1) + 1] = crc32_u64(Seed2[(vi << 1 + 1)], alleles[vi]);
#else
				Seed1[vi] = (HASH)crc32_u64(Seed1[vi], alleles[vi]);
#endif
			}
			SetZero(Seed1, N_MAX_PLOIDY);
		}
	}

	for (int vi = 0; vi < ploidy; ++vi)
	{
#ifdef HASH64
		Seed2[(vi << 1)    ] = crc32_u64(Seed2[(vi << 1)    ], alleles[vi]);
		Seed2[(vi << 1) + 1] = crc32_u64(Seed2[(vi << 1 + 1)], alleles[vi]);
#else
		Seed1[vi] = crc32_u64(Seed1[vi], alleles[vi]);
#endif
	}

	SetVal(hash, Seed1, ploidy);
}

/* 32 bit Integer Hashing */
TARGET uint MurmurHash2(uint data, uint seed)
{
	uint m = 0x5bd1e995;
	uint s = seed ^ sizeof(uint);

	uint a = data;

	a *= m;
	a ^= a >> 24;
	a *= m;
	s *= m;
	a ^= s;

	a ^= a >> 13;
	a *= m;
	a ^= a >> 15;

	return a;
}

/* 64 bit Integer Hashing */
TARGET uint64 MurmurHash64(uint64 data, uint64 seed)
{
	uint* d = (uint*)&data;
	uint* s = (uint*)&seed;

	d[1] = (~d[0]) ^ d[1];
	s[1] = (~s[0]) ^ s[1];

	return ((uint64)MurmurHash2(d[1], s[1]) << 32) |
		((uint64)MurmurHash2(d[0], s[0]));
}

/* 64 bit Integer Hashing */
TARGET uint MurmurHash32(uint64 data, uint64 seed)
{
	uint a32 = (~(uint)data) ^ ((uint)(data >> 32));
	uint s32 = (~(uint)seed) ^ ((uint)(seed >> 32));

	return MurmurHash2(a32, s32);
}

/* Mix high and low 32 bits */
TARGET uint Mix(uint64 x)
{
	return (~(uint)x) ^ ((uint)(x >> 32));
}