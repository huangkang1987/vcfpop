/* Allelic Depth Analysis Functions, TEST, DO NOT USE */

#pragma once
#include "vcfpop.h"

#define extern 
#undef extern 

template struct IND<double>;
template struct IND<float >;

template TARGET void IND<double>::SetAlleleDepth(int64 l, uint* depth, int K, int indid);
template TARGET void IND<float >::SetAlleleDepth(int64 l, uint* depth, int K, int indid);
template TARGET void IND<double>::SetAlleleDepth(int64 l, uint* depth, int K);
template TARGET void IND<float >::SetAlleleDepth(int64 l, uint* depth, int K);
template TARGET void IND<double>::SetAlleleDepth(int64 l, uint* depth, int K, OFFSET* _offset, byte* bucket);
template TARGET void IND<float >::SetAlleleDepth(int64 l, uint* depth, int K, OFFSET* _offset, byte* bucket);
template TARGET void IND<double>::GetAlleleDepth(int64 l, uint* depth);
template TARGET void IND<float >::GetAlleleDepth(int64 l, uint* depth);

/* Set allele sequencing depth, for ad, TEST */
template<typename REAL>
TARGET void IND<REAL>::SetAlleleDepth(int64 l, uint* depth, int K, int indid)
{
	uint size = (uint)ad_bucket.offset[l].size, mask = ((1u << size) - 1u);
	uint64 offset = size * K * indid;
	uint* pos = (uint*)(ad_bucket.base_addr + ad_bucket.offset[l].offset + (offset >> 3));
	offset &= 7;
	for (int k = 0; k < K; ++k)
	{
		*pos = (*pos & (~(mask << offset))) | (depth[k] << offset);
		offset += size;
		if (offset > 7)
		{
			pos = (uint*)((byte*)pos + 1);
			offset -= 8;
		}
	}
}

/* Set allele sequencing depth, for ad, TEST */
template<typename REAL>
TARGET void IND<REAL>::SetAlleleDepth(int64 l, uint* depth, int K)
{
	uint size = (uint)ad_bucket.offset[l].size, mask = ((1u << size) - 1u);
	uint64 offset = size * K * indid;
	uint* pos = (uint*)(ad_bucket.base_addr + ad_bucket.offset[l].offset + (offset >> 3));
	offset &= 7;
	for (int k = 0; k < K; ++k)
	{
		*pos = (*pos & (~(mask << offset))) | (depth[k] << offset);
		offset += size;
		if (offset > 7)
		{
			pos = (uint*)((byte*)pos + 1);
			offset -= 8;
		}
	}
}

/* Set allele sequencing depth, for ad, TEST */
template<typename REAL>
TARGET void IND<REAL>::SetAlleleDepth(int64 l, uint* depth, int K, OFFSET* _offset, byte* bucket)
{
	uint size = (uint)_offset[l].size;
	uint mask = ((1u << size) - 1u);
	uint64 offset = size * K * indid;
	uint* pos = (uint*)(bucket + _offset[l].offset + (offset >> 3));
	offset &= 7;
	for (int k = 0; k < K; ++k)
	{
		*pos = (*pos & (~(mask << offset))) | (depth[k] << offset);
		offset += size;
		if (offset > 7)
		{
			pos = (uint*)((byte*)pos + 1);
			offset -= 8;
		}
	}
}

/* Set allele sequencing depth, for ad, TEST */
template<typename REAL>
TARGET void IND<REAL>::GetAlleleDepth(int64 l, uint* depth)
{
	uint K = GetLoc(l).k;
	uint size = (uint)ad_bucket.offset[l].size;
	uint mask = ((1u << size) - 1u);
	uint64 offset = size * K * indid;
	uint* pos = (uint*)(ad_bucket.base_addr + ad_bucket.offset[l].offset + (offset >> 3));
	offset &= 7;
	for (uint k = 0; k < K; ++k)
	{
		depth[k] = (*pos >> offset) & mask;
		offset += size;
		if (offset > 7)
		{
			pos = (uint*)((byte*)pos + 1);
			offset -= 8;
		}
	}
}