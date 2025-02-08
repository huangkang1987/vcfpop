/* Math  Functions */

#include "vcfpop.h"

template struct RNGSIMD<double>;
template struct RNGSIMD<float >;
template TARGETSIMD void RNGSIMD<double>::Integer<uint  >(uint  * re, int64 n, uint   minv, uint   maxv);
template TARGETSIMD void RNGSIMD<double>::Integer<uint64>(uint64* re, int64 n, uint64 minv, uint64 maxv);
template TARGETSIMD void RNGSIMD<float >::Integer<uint  >(uint  * re, int64 n, uint   minv, uint   maxv);
template TARGETSIMD void RNGSIMD<float >::Integer<uint64>(uint64* re, int64 n, uint64 minv, uint64 maxv);

#ifndef _RNGSIMD

/* Initialize rng */
template<typename REAL>
TARGETSIMD RNGSIMD<REAL>::RNGSIMD()
{

}

/* Initialize rng */
template<typename REAL>
TARGETSIMD RNGSIMD<REAL>::RNGSIMD(uint64 seed, uint64 salt)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: new (this) RNGNEO<REAL>(seed, salt); return;
#else
	case 4: new (this) RNG512<REAL>(seed, salt); return;
	case 3: new (this) RNGAVX<REAL>(seed, salt); return;
	case 2: new (this) RNGSSE<REAL>(seed, salt); return;
#endif
	}
	if constexpr (std::is_same_v<REAL, double>)
	{
		uint64* x = (uint64*)(data + 0);
		uint64* y = (uint64*)(data + 512);
		uint64 a[64]; uint* aa = (uint*)a;

		//REP(32) { a[kk] = vld1q_u64(((uint64[]) { seed, seed + 1 })); seed += 2; }
		REP(64) a[kk] = seed + kk;

		//s = vdupq_n_u64(salt);
		uint64 s = salt;
		uint* ss = (uint*)&s;

		//m = vdupq_n_u32(0x5bd1e995);
		uint m = 0x5bd1e995;

		//REP(32) a[kk] = veorq_u64(a[kk], vshlq_n_u64(vmvnq_u32(a[kk]), 32));
		REP(64) a[kk] = (a[kk] ^ ((~a[kk]) << 32));

		//s = veorq_u64(s, vshlq_n_u64(vmvnq_u32(s), 32));
		s = (s ^ ((~s) << 32));

		// uint s = s ^ 4;
		//s = veorq_u32(s, vdupq_n_u32(4));
		s = s ^ 0x0000010000000100ULL;

		// a *= m;
		//REP(32) a[kk] = vmulq_u32(a[kk], m);
		REP(128) aa[kk] = aa[kk] * m;

		// a ^= a >> 24;
		//REP(32) a[kk] = veorq_u64(a[kk], vshrq_n_u32(a[kk], 24));
		REP(64) a[kk] = a[kk] ^ ((a[kk] & 0xFFFFF000FFFFF000ULL) >> 24);

		// a *= m;
		//REP(32) a[kk] = vmulq_u32(a[kk], m);
		REP(128) aa[kk] = aa[kk] * m;

		// s *= m;
		//s = vmulq_u32(s, m);
		REP(2) ss[kk] = ss[kk] * m;

		// a ^= s;
		//REP(32) a[kk] = veorq_u64(a[kk], s);
		REP(64) a[kk] = a[kk] ^ s;

		// a ^= a >> 13;
		//REP(32) a[kk] = veorq_u64(a[kk], vshrq_n_u32(a[kk], 13));
		REP(64) a[kk] = a[kk] ^ (a[kk] & 0xFFFFE000FFFFE000ULL) >> 13;

		// a *= m;
		//REP(32) a[kk] = vmulq_u32(a[kk], m);
		REP(128) aa[kk] = aa[kk] * m;

		// a ^= a >> 15;
		//REP(32) a[kk] = veorq_u64(a[kk], vshrq_n_u32(a[kk], 15));
		REP(64) a[kk] = a[kk] ^ (a[kk] & 0xFFFF8000FFFF8000ULL) >> 15;


		// original
		//REP(32) x[kk] = veorq_u64(vdupq_n_u64(0x159A55E5075BCD15), a[kk]);
		REP(64) x[kk] = 0x159A55E5075BCD15ULL ^ a[kk];

		//REP(32) a[kk] = vshlq_n_u64(a[kk], 6);
		REP(64) a[kk] = a[kk] << 6;

		//REP(32) y[kk] = veorq_u64(vdupq_n_u64(0x054913331F123BB5), a[kk]);
		REP(64) y[kk] = 0x054913331F123BB5 ^ a[kk];
	}
	else
	{
		uint* x = (uint*)(data + 0);
		uint* y = (uint*)(data + 256);
		uint* z = (uint*)(data + 512);

		uint a[16], s = Mix(salt), m = 0x5bd1e995;

		// REP(16) { a[kk] = vld1q_u32(((uint[]) { Mix(seed + 0), Mix(seed + 1), Mix(seed + 2), Mix(seed + 3) })); seed += 4; }
		REP(64) a[kk] = Mix(seed + kk);

		// uint s = s ^ 4;
		//s = veorq_u32(s, vdupq_n_u32(4));
		s = s ^ 4;

		// a *= m;
		//REP(16) a[kk] = vmulq_u32(a[kk], m);
		REP(64) a[kk] = a[kk] * m;

		// a ^= a >> 24;
		//REP(16) a[kk] = veorq_u32(a[kk], vshrq_n_u32(a[kk], 24));
		REP(64) a[kk] = a[kk] ^ (a[kk] >> 24);

		// a *= m;
		//REP(16) a[kk] = vmulq_u32(a[kk], m);
		REP(64) a[kk] = a[kk] * m;

		// s *= m;
		//s = vmulq_u32(s, m);
		s = s * m;

		// a ^= s;
		//REP(16) a[kk] = veorq_u32(a[kk], s);
		REP(64) a[kk] = a[kk] ^ s;

		// a ^= a >> 13;
		//REP(16) a[kk] = veorq_u32(a[kk], vshrq_n_u32(a[kk], 13));
		REP(64) a[kk] = a[kk] ^ (a[kk] >> 13);

		// a *= m;
		//REP(16) a[kk] = vmulq_u32(a[kk], m);
		REP(64) a[kk] = a[kk] * m;

		// a ^= a >> 15;
		//REP(16) a[kk] = veorq_u32(a[kk], vshrq_n_u32(a[kk], 15));
		REP(64) a[kk] = a[kk] ^ (a[kk] >> 15);

		// original
		//REP(16) x[kk] = veorq_u32(vdupq_n_u32(0x075BCD15), a[kk]);
		REP(64) x[kk] = 0x075BCD15 ^ a[kk];

		//REP(16) a[kk] = vshlq_n_u32(a[kk], 3);
		REP(64) a[kk] = a[kk] << 3;

		//REP(16) y[kk] = veorq_u32(vdupq_n_u32(0x159A55E5), a[kk]);
		REP(64) y[kk] = 0x159A55E5 ^ a[kk];

		//REP(16) a[kk] = vshlq_n_u32(a[kk], 3);
		REP(64) a[kk] = a[kk] << 3;

		//REP(16) z[kk] = veorq_u32(vdupq_n_u32(0x1F123BB5), a[kk]);
		REP(64) z[kk] = 0x1F123BB5 ^ a[kk];
	}
}

/* Draw 64 64-bit integers in [0,n), 64*n frequencies are in arr */
template<typename REAL>
TARGETSIMD void RNGSIMD<REAL>::Poly(REAL* arr, int n, int64* re)
{
	if constexpr (std::is_same_v<REAL, double>)
	{
		switch (SIMD_TYPE)
		{
#ifdef __aarch64__
		case 2: ((RNGNEO<double>*)this)->Poly((float64x2_t*)arr, n, re); return;
#else
		case 4: ((RNG512<double>*)this)->Poly((__m512d*)arr, n, (__m512i*)re); return;
		case 3: ((RNGAVX<double>*)this)->Poly((__m256d*)arr, n, (__m256i*)re); return;
		case 2: ((RNGSSE<double>*)this)->Poly((__m128d*)arr, n, (__m128i*)re); return;
#endif
		}
	}
	else
	{
		switch (SIMD_TYPE)
		{
#ifdef __aarch64__
		case 2: ((RNGNEO<float>*)this)->Poly((float32x4_t*)arr, n, re); return;
#else
		case 4: ((RNG512<float>*)this)->Poly((__m512*)arr, n, (__m512i*)re); return;
		case 3: ((RNGAVX<float>*)this)->Poly((__m256*)arr, n, (__m256i*)re); return;
		case 2: ((RNGSSE<float>*)this)->Poly((__m128*)arr, n, (__m128i*)re); return;
#endif
		}
	}

	if constexpr (std::is_same_v<REAL, double>)
	{
		uint64* x = (uint64*)(data + 0);
		uint64* y = (uint64*)(data + 512);

		double t[64], s[64];
		double one = 1.0;
		uint64 mask1 = 0x000FFFFFFFFFFFFF;
		uint64 mask2 = 0x3FF0000000000000;
		uint64* r = (uint64*)t; uint64* re64 = (uint64*)re; 

		REP(64) s[kk] = 0;

		for (int i64 = 0; i64 < n * 64; i64 += 64)
			REP(64) s[kk] = s[kk] + arr[kk + i64];

		XorShift();

		REP(64) r[kk] = x[kk] + y[kk];

		REP(64) r[kk] = r[kk] & mask1;

		REP(64) r[kk] = r[kk] | mask2;

		REP(64) t[kk] = t[kk] - one;

		REP(64) t[kk] = t[kk] * s[kk];

		uint64 midx[64], nidx = 0, ninc = 1;
		uint64 f[64], b[64];
		REP(64) midx[kk] = n - 1;
		REP(64) f[kk] = 0;

		for (int i64 = 0; i64 < n * 64; i64 += 64)
		{
			REP(64) b[kk] = t[kk] < arr[kk + i64];

			REP(64) t[kk] = t[kk] - arr[kk + i64];

			REP(64) b[kk] = (~f[kk]) & b[kk];

			REP(64) f[kk] = f[kk] | b[kk];

			REP(64) midx[kk] = b[kk] ? nidx : midx[kk];

			nidx = nidx + ninc;
		}

		REP(64) re64[kk] = midx[kk];
	}
	else
	{
		uint* z = (uint*)(data + 512);

		float t[64], s[64]; 
		float one = 1.0f;
		uint mask1 = 0x007FFFFF;
		uint mask2 = 0x3F800000;
		uint* r = (uint*)t; uint64* re64 = (uint64*)re; 

		REP(64) s[kk] = 0;

		for (int i64 = 0; i64 < n * 64; i64 += 64)
			REP(64) s[kk] = s[kk] + arr[kk + i64];

		XorShift();

		REP(64) r[kk] = z[kk] & mask1;

		REP(64) r[kk] = r[kk] | mask2;

		REP(64) t[kk] = t[kk] - one;

		REP(64) t[kk] = t[kk] * s[kk];

		uint midx[64], nidx = 0, ninc = 1;
		uint f[64], b[64];
		REP(64) midx[kk] = n - 1;
		REP(64) f[kk] = 0;

		for (int i64 = 0; i64 < n * 64; i64 += 64)
		{
			REP(64) b[kk] = t[kk] < arr[kk + i64];

			REP(64) t[kk] = t[kk] - arr[kk + i64];

			REP(64) b[kk] = (~f[kk]) & b[kk];

			REP(64) f[kk] = f[kk] | b[kk];

			REP(64) midx[kk] = b[kk] ? nidx : midx[kk];

			nidx = nidx + ninc;
		}

		REP(64) re64[kk] = midx[kk];
	}
}

/* Draw uniform distriubted intergers */
template<typename REAL>
TARGETSIMD void RNGSIMD<REAL>::XorShift()
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: ((RNGNEO<REAL>*)this)->XorShift(); return;
#else
	case 4: ((RNG512<REAL>*)this)->XorShift(); return;
	case 3: ((RNGAVX<REAL>*)this)->XorShift(); return;
	case 2: ((RNGSSE<REAL>*)this)->XorShift(); return;
#endif
	}

	if constexpr (std::is_same_v<REAL, double>)
	{
		uint64* x = (uint64*)(data + 0);
		uint64* y = (uint64*)(data + 512);
		uint64 a[64], b[64];

		//REP(32) a[kk] = x[kk];
		REP(64) a[kk] = x[kk];

		//REP(32) b[kk] = y[kk];
		REP(64) b[kk] = y[kk];

		//REP(32) x[kk] = b[kk];
		REP(64) x[kk] = b[kk];

		//REP(32) a[kk] = veorq_u64(a[kk], vshlq_n_u64(a[kk], 23));
		REP(64) a[kk] = a[kk] ^ (a[kk] << 23);

		//REP(32) a[kk] = veorq_u64(a[kk], vshrq_n_u64(a[kk], 18));
		REP(64) a[kk] = a[kk] ^ (a[kk] >> 18);

		//REP(32) a[kk] = veorq_u64(a[kk], b[kk]);
		REP(64) a[kk] = a[kk] ^ b[kk];

		//REP(32) a[kk] = veorq_u64(a[kk], vshrq_n_u64(b[kk], 5));
		REP(64) a[kk] = a[kk] ^ (b[kk] >> 5);

		//REP(32) y[kk] = a[kk];
		REP(64) y[kk] = a[kk];

		// use x + y
	}
	else
	{
		uint* x = (uint*)(data + 0);
		uint* y = (uint*)(data + 256);
		uint* z = (uint*)(data + 512);
		uint u[64];

		//REP(16) u[kk] = vshlq_n_u32(x[kk], 16);
		REP(64) u[kk] = x[kk] << 16;

		//REP(16) x[kk] = veorq_u32(x[kk], u[kk]);
		REP(64) x[kk] = x[kk] ^ u[kk];

		//REP(16) u[kk] = vshrq_n_u32(x[kk], 5);
		REP(64) u[kk] = x[kk] >> 5;

		//REP(16) x[kk] = veorq_u32(x[kk], u[kk]);
		REP(64) x[kk] = x[kk] ^ u[kk];

		//REP(16) u[kk] = vshlq_n_u32(x[kk], 1);
		REP(64) u[kk] = x[kk] << 1;

		//REP(16) x[kk] = veorq_u32(x[kk], u[kk]);
		REP(64) x[kk] = x[kk] ^ u[kk];

		//REP(16) u[kk] = x[kk];
		REP(64) u[kk] = x[kk];

		//REP(16) x[kk] = y[kk];
		REP(64) x[kk] = y[kk];

		//REP(16) y[kk] = z[kk];
		REP(64) y[kk] = z[kk];

		//REP(16) z[kk] = veorq_u32(u[kk], x[kk]);
		REP(64) z[kk] = u[kk] ^ x[kk];

		//REP(16) z[kk] = veorq_u32(z[kk], y[kk]);
		REP(64) z[kk] = z[kk] ^ y[kk];

		// use z
	}
}

/* Draw uniform distriubted integers */
template<typename REAL>
template<typename INT>
TARGETSIMD void RNGSIMD<REAL>::Integer(INT* re, int64 n, INT minv, INT maxv)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: ((RNGNEO<REAL>*)this)->Integer(re, n, minv, maxv); return;
#else
	case 4: ((RNG512<REAL>*)this)->Integer(re, n, minv, maxv); return;
	case 3: ((RNGAVX<REAL>*)this)->Integer(re, n, minv, maxv); return;
	case 2: ((RNGSSE<REAL>*)this)->Integer(re, n, minv, maxv); return;
#endif
	}

	if constexpr (std::is_same_v<REAL, double>)
	{
		uint64* x = (uint64*)(data + 0);
		uint64* y = (uint64*)(data + 512);
		uint64* rei = (uint64*)re;
		int64 i = 0;
		INT modv = maxv - minv;

		for (; i <= n - 512 / sizeof(INT); i += 512 / sizeof(INT))
		{
			XorShift();
			REP(512 / sizeof(uint64)) rei[kk] = x[kk] + y[kk];
			rei += 512 / sizeof(uint64);
		}

		if (i != n)
		{
			uint64_t re2[512 / sizeof(uint64)];
			XorShift();
			REP(512 / sizeof(uint64)) re2[kk] = x[kk] + y[kk];
			SetVal((int*)rei, (int*)re2, n - i);
		}

		if (maxv != (INT)-1 || minv != 0)
		{
			for (i = 0; i < n; ++i)
				re[i] = re[i] % modv + minv;
		}
	}
	else
	{
		uint* z = (uint*)(data + 512);
		uint* rei = (uint*)re;
		int64 i = 0;
		INT modv = maxv - minv;

		for (; i <= n - 256 / sizeof(INT); i += 256 / sizeof(INT))
		{
			XorShift();
			REP(256 / sizeof(uint)) rei[kk] = z[kk];
			rei += 256 / sizeof(uint);
		}

		if (i != n)
		{
			XorShift();
			SetVal((int*)rei, (int*)z, n - i);
		}

		if (maxv != (INT)-1 || minv != 0)
		{
			for (i = 0; i < n; ++i)
				re[i] = re[i] % modv + minv;
		}
	}
}

/* Draw uniform distriubted real numbers */
template<typename REAL>
TARGETSIMD void RNGSIMD<REAL>::Uniform(REAL* re, int n, REAL minv, REAL maxv)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: ((RNGNEO<REAL>*)this)->Uniform(re, n, minv, maxv); return;
#else
	case 4: ((RNG512<REAL>*)this)->Uniform(re, n, minv, maxv); return;
	case 3: ((RNGAVX<REAL>*)this)->Uniform(re, n, minv, maxv); return;
	case 2: ((RNGSSE<REAL>*)this)->Uniform(re, n, minv, maxv); return;
#endif
	}

	if constexpr (std::is_same_v<REAL, double>)
	{
		uint64* x = (uint64*)(data + 0);
		uint64* y = (uint64*)(data + 512);
		uint64* rei = (uint64*)re;
		REAL*& ref = *(REAL**)&rei;

		int i = 0;
		double range = maxv - minv;

		uint64 mask1 = 0x000FFFFFFFFFFFFF;
		uint64 mask2 = 0x3FF0000000000000;
		double v1 = minv - range;
		double v2 = range;

		for (; i <= n - 64; i += 64)
		{
			XorShift();

			REP(64) rei[kk] = x[kk] + y[kk];

			REP(64) rei[kk] = rei[kk] & mask1;

			REP(64) rei[kk] = rei[kk] | mask2;

			REP(64) ref[kk] = ref[kk] * v2;

			REP(64) ref[kk] = ref[kk] + v1;

			rei += 64;
			//ref += 64;
		}

		if (i != n)
		{
			uint64 rei2[64];
			double* ref2 = (double*)rei2;

			XorShift();

			REP(64) rei2[kk] = x[kk] + y[kk];

			REP(64) rei2[kk] = rei2[kk] & mask1;

			REP(64) rei2[kk] = rei2[kk] | mask2;

			REP(64) ref2[kk] = ref2[kk] * v2;

			REP(64) ref2[kk] = ref2[kk] + v1;

			SetVal((double*)rei, (double*)rei2, n - i);
		}
	}
	else
	{
		uint* z = (uint*)(data + 512);
		uint* rei = (uint*)re;
		REAL*& ref = *(REAL**)&rei;

		int i = 0;
		float range = maxv - minv;

		uint mask1 = 0x007FFFFF;
		uint mask2 = 0x3F800000;
		float v1 = minv - range;
		float v2 = range;

		for (; i <= n - 64; i += 64)
		{
			XorShift();

			REP(64) rei[kk] = z[kk] & mask1;

			REP(64) rei[kk] = rei[kk] | mask2;

			REP(64) ref[kk] = ref[kk] * v2;

			REP(64) ref[kk] = ref[kk] + v1;

			rei += 64;
			//ref += 64;
		}

		if (i != n)
		{
			uint rei2[64];
			float* ref2 = (float*)rei2;

			XorShift();

			REP(64) rei2[kk] = z[kk] & mask1;

			REP(64) rei2[kk] = rei2[kk] | mask2;

			REP(64) ref2[kk] = ref2[kk] * v2;

			REP(64) ref2[kk] = ref2[kk] + v1;

			SetVal((float*)rei, (float*)rei2, n - i);
		}
	}
}

/* Draw uniform distriubted real numbers */
template<typename REAL>
TARGETSIMD void RNGSIMD<REAL>::Normal(REAL* re, int n, REAL mean, REAL sd)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: ((RNGNEO<REAL>*)this)->Normal(re, n, mean, sd); return;
#else
	case 4: ((RNG512<REAL>*)this)->Normal(re, n, mean, sd); return;
	case 3: ((RNGAVX<REAL>*)this)->Normal(re, n, mean, sd); return;
	case 2: ((RNGSSE<REAL>*)this)->Normal(re, n, mean, sd); return;
#endif
	}
	if constexpr (std::is_same_v<REAL, double>)
	{
		uint64* x = (uint64*)(data + 0);
		uint64* y = (uint64*)(data + 512);
		uint64* rei = (uint64*)re;
		REAL*& ref = *(REAL**)&rei;

		int i = 0;

		uint64 mask1 = 0x000FFFFFFFFFFFFF;
		uint64 mask2 = 0x3FF0000000000000;
		double v1 = -1;

		double min_freq = MIN_FREQ;
		double pi2 = 2.0 * M_PI;
		double mu = mean;
		double s = sd;

		for (; i <= n - 64; i += 64)
		{
			XorShift();

			REP(64) rei[kk] = x[kk] + y[kk];
			REP(64) rei[kk] = rei[kk] & mask1;
			REP(64) rei[kk] = rei[kk] | mask2;
			REP(64) ref[kk] = ref[kk] + v1;

			double u1, u2, u3, u4;
			for (int j = 0; j < 32; ++j)
			{
				u1 = std::max(ref[j], min_freq);
				u2 = ref[j + 32] * pi2;

				u1 = sqrt(-2.0 * log(u1));
				u3 = cos(u2);
				u4 = sin(u2);

				ref[j     ] = u1 * u3;
				ref[j + 32] = u1 * u4;
			}

			REP(64) ref[kk] = ref[kk] * s;
			REP(64) ref[kk] = ref[kk] + mu;

			rei += 64;
			//ref += 64;
		}

		if (i != n)
		{
			uint64 rei2[64];
			double* ref2 = (double*)rei2;

			XorShift();

			REP(64) rei2[kk] = x[kk] + y[kk];
			REP(64) rei2[kk] = rei2[kk] & mask1;
			REP(64) rei2[kk] = rei2[kk] | mask2;
			REP(64) ref2[kk] = ref2[kk] + v1;

			double u1, u2, u3, u4;
			for (int j = 0; j < 32; ++j)
			{
				u1 = std::max(ref2[j], min_freq);
				u2 = ref2[j + 32] * pi2;

				u1 = sqrt(-2.0 * log(u1));
				u3 = cos(u2);
				u4 = sin(u2);

				ref2[j] = u1 * u3;
				ref2[j + 32] = u1 * u4;
			}

			REP(64) ref2[kk] = ref2[kk] * s;
			REP(64) ref2[kk] = ref2[kk] + mu;

			SetVal((double*)rei, (double*)rei2, n - i);
		}
	}
	else
	{
		uint* z = (uint*)(data + 512);
		uint* rei = (uint*)re;
		REAL*& ref = *(REAL**)&rei;

		int i = 0;

		uint mask1 = 0x007FFFFF;
		uint mask2 = 0x3F800000;
		float v1 = -1;
		float min_freq = (float)MIN_FREQ;
		float pi2 = (float)(2.0 * M_PI);
		float mu = mean;
		float s = sd;

		for (; i <= n - 64; i += 64)
		{
			XorShift();

			REP(64) rei[kk] = z[kk] & mask1;
			REP(64) rei[kk] = rei[kk] | mask2;
			REP(64) ref[kk] = ref[kk] + v1;

			float u1, u2, u3, u4;
			for (int j = 0; j < 32; ++j)
			{
				u1 = std::max(ref[j], min_freq);
				u2 = ref[j + 32] * pi2;

				u1 = sqrt(-2.0 * log(u1));
				u3 = cos(u2);
				u4 = sin(u2);

				ref[j     ] = u1 * u3;
				ref[j + 32] = u1 * u4;
			}

			REP(64) ref[kk] = ref[kk] * s;
			REP(64) ref[kk] = ref[kk] + mu;

			rei += 64;
			//ref += 64;
		}

		if (i != n)
		{
			uint rei2[64];
			float* ref2 = (float*)rei2;
			
			XorShift();

			REP(64) rei2[kk] = z[kk] & mask1;
			REP(64) rei2[kk] = rei2[kk] | mask2;
			REP(64) ref2[kk] = ref2[kk] + v1;

			float u1, u2, u3, u4;
			for (int j = 0; j < 32; ++j)
			{
				u1 = std::max(ref2[j], min_freq);
				u2 = ref2[j + 32] * pi2;

				u1 = sqrt(-2.0 * log(u1));
				u3 = cos(u2);
				u4 = sin(u2);

				ref2[j     ] = u1 * u3;
				ref2[j + 32] = u1 * u4;
			}

			REP(64) ref2[kk] = ref2[kk] * s;
			REP(64) ref2[kk] = ref2[kk] + mu;

			SetVal((float*)rei, (float*)rei2, n - i);
		}
	}
}
#endif

/* Find the index of the mimumum element */
TARGET int64 GetMinIdx(double* A, int64 n, double& val)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return GetMinIdxNEO(A, n, val);
#else
	case 4: return GetMinIdx512(A, n, val);
	case 3: return GetMinIdxAVX(A, n, val);
	case 2: return GetMinIdxSSE(A, n, val);
#endif
	}

	val = DBL_MAX;
	int64 idx = -1;
	for (int64 i = 0; i < n; ++i)
	{
		if (A[i] > val) continue;
		val = A[i];
		idx = i;
	}
	return idx;
}

/* Find the index of the mimumum element */
TARGET int64 GetMinIdx(float* A, int64 n, float& val)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return GetMinIdxNEO(A, n, val);
#else
	case 4: return GetMinIdx512(A, n, val);
	case 3: return GetMinIdxAVX(A, n, val);
	case 2: return GetMinIdxSSE(A, n, val);
#endif
	}

	val = FLT_MAX;
	int64 idx = -1;
	for (int64 i = 0; i < n; ++i)
	{
		if (A[i] > val) continue;
		val = A[i];
		idx = i;
	}
	return idx;
}

/* Find maximum and minimum element of A */
TARGET void GetMinMaxVal(double* A, int64 n, double& minv, double& maxv)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return GetMinMaxValNEO(A, n, minv, maxv);
#else
	case 4: return GetMinMaxVal512(A, n, minv, maxv);
	case 3: return GetMinMaxValAVX(A, n, minv, maxv);
	case 2: return GetMinMaxValSSE(A, n, minv, maxv);
#endif
	}
	
	minv = DBL_MAX;
	maxv = -DBL_MAX;
	for (int64 i = 0; i < n; ++i)
	{
		if (A[i] < minv) minv = A[i];
		if (A[i] > maxv) maxv = A[i];
	}
}

/* Find maximum and minimum element of A */
TARGET void GetMinMaxVal(float* A, int64 n, float& minv, float& maxv)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return GetMinMaxValNEO(A, n, minv, maxv);
#else
	case 4: return GetMinMaxVal512(A, n, minv, maxv);
	case 3: return GetMinMaxValAVX(A, n, minv, maxv);
	case 2: return GetMinMaxValSSE(A, n, minv, maxv);
#endif
	}
	
	minv = FLT_MAX;
	maxv = -FLT_MAX;
	for (int64 i = 0; i < n; ++i)
	{
		if (A[i] < minv) minv = A[i];
		if (A[i] > maxv) maxv = A[i];
	}
}

/* Find minimum element of A */
TARGET double GetMinVal(double* A, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return GetMinValNEO(A, n);
#else
	case 4: return GetMinVal512(A, n);
	case 3: return GetMinValAVX(A, n);
	case 2: return GetMinValSSE(A, n);
#endif
	}
	
	double val = DBL_MAX;
	for (int64 i = 0; i < n; ++i)
	{
		if (A[i] > val) continue;
		val = A[i];
	}
	return val;
}

/* Find minimum element of A */
TARGET float GetMinVal(float* A, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return GetMinValNEO(A, n);
#else
	case 4: return GetMinVal512(A, n);
	case 3: return GetMinValAVX(A, n);
	case 2: return GetMinValSSE(A, n);
#endif
	}
	
	float val = FLT_MAX;
	for (int64 i = 0; i < n; ++i)
	{
		if (A[i] > val) continue;
		val = A[i];
	}
	return val;
}

/* Find minimum element of A */
TARGET int64 GetMinVal(int64* A, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return GetMinValNEO(A, n);
#else
	case 4: return GetMinVal512(A, n);
	case 3: return GetMinValAVX(A, n);
	case 2: return GetMinValSSE(A, n);
#endif
	}

	int64 val = 0x7FFFFFFFFFFFFFFF;
	for (int64 i = 0; i < n; ++i)
	{
		if (A[i] > val) continue;
		val = A[i];
	}
	return val;
}

/* Find Maximum element of A */
TARGET double GetMaxVal(double* A, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return GetMaxValNEO(A, n);
#else
	case 4: return GetMaxVal512(A, n);
	case 3: return GetMaxValAVX(A, n);
	case 2: return GetMaxValSSE(A, n);
#endif
	}
	
	double val = -DBL_MAX;
	for (int64 i = 0; i < n; ++i)
	{
		if (A[i] < val) continue;
		val = A[i];
	}
	return val;
}

/* Find Maximum element of A */
TARGET float GetMaxVal(float* A, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return GetMaxValNEO(A, n);
#else
	case 4: return GetMaxVal512(A, n);
	case 3: return GetMaxValAVX(A, n);
	case 2: return GetMaxValSSE(A, n);
#endif
	}
	
	float val = -FLT_MAX;
	for (int64 i = 0; i < n; ++i)
	{
		if (A[i] < val) continue;
		val = A[i];
	}
	return val;
}

/* Find Maximum element of A */
TARGET double GetMaxVal(double* A, int64 n, int64 sep)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return GetMaxValNEO(A, n, sep);
#else
	case 4: return GetMaxVal512(A, n, sep);
	case 3: return GetMaxValAVX(A, n, sep);
	case 2: return GetMaxValSSE(A, n, sep);
#endif
	}

	double val = -DBL_MAX;
	for (int64 i = 0; i < n; ++i, A += sep)
	{
		if (*A < val) continue;
		val = *A;
	}
	return val;
}

/* Find Maximum element of A */
TARGET float GetMaxVal(float* A, int64 n, int64 sep)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return GetMaxValNEO(A, n, sep);
#else
	case 4: return GetMaxVal512(A, n, sep);
	case 3: return GetMaxValAVX(A, n, sep);
	case 2: return GetMaxValSSE(A, n, sep);
#endif
	}

	float val = -FLT_MAX;
	for (int64 i = 0; i < n; ++i, A += sep)
	{
		if (*A < val) continue;
		val = *A;
	}
	return val;
}

/* A[i] = B[i] */
TARGET void SetVal(uint* A, ushort* B, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SetValNEO(A, B, n);
#else
	case 4: return SetVal512(A, B, n);
	case 3: return SetValAVX(A, B, n);
	case 2: return SetValSSE(A, B, n);
#endif
	}
	
	for (int64 i = 0; i < n; ++i)
		A[i] = B[i];
}

/* log(prod(A[i++])) */
TARGET double LogProd(double* A, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return LogProdNEO(A, n);
#else
	case 4: return LogProd512(A, n);
	case 3: return LogProdAVX(A, n);
	case 2: return LogProdSSE(A, n);
#endif
	}

	int64 slog = 0; double prod = 1;
	for (int64 i = 0; i < n; ++i)
		ChargeLog(slog, prod, A[i]);

	CloseLog(slog, prod);
	return prod;
}

/* log(prod(A[i++])) */
TARGET double LogProd(float* A, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return LogProdNEO(A, n);
#else
	case 4: return LogProd512(A, n);
	case 3: return LogProdAVX(A, n);
	case 2: return LogProdSSE(A, n);
#endif
	}

	int64 slog = 0; double prod = 1;
	for (int64 i = 0; i < n; ++i)
		ChargeLog(slog, prod, A[i]);

	CloseLog(slog, prod);
	return prod;
}

/* log(prod(A[i += sep])) */
TARGET double LogProd(double* A, int64 n, int64 sep)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return LogProdNEO(A, n, sep);
#else
	case 4: return LogProd512(A, n, sep);
	case 3: return LogProdAVX(A, n, sep);
	case 2: return LogProdSSE(A, n, sep);
#endif
	}

	int64 slog = 0; double prod = 1;
	for (int64 i = 0; i < n; ++i, A += sep)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);
	return prod;
}

/* log(prod(A[i += sep])) */
TARGET double LogProd(float* A, int64 n, int64 sep)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return LogProdNEO(A, n, sep);
#else
	case 4: return LogProd512(A, n, sep);
	case 3: return LogProdAVX(A, n, sep);
	case 2: return LogProdSSE(A, n, sep);
#endif
	}

	int64 slog = 0; double prod = 1;
	for (int64 i = 0; i < n; ++i, A += sep)
		ChargeLog(slog, prod, *A);

	CloseLog(slog, prod);
	return prod;
}

/* log(prod(A[i += sep] / B[i += sep])) */
TARGET double LogProdDiv(double* A, double* B, int64 n, int64 sep)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return LogProdDivNEO(A, B, n, sep);
#else
	case 4: return LogProdDiv512(A, B, n, sep);
	case 3: return LogProdDivAVX(A, B, n, sep);
	case 2: return LogProdDivSSE(A, B, n, sep);
#endif
	}

	int64 slog1 = 0; double prod1 = 1;
	int64 slog2 = 0; double prod2 = 1;
	for (int64 i = 0; i < n; ++i, A += sep, B += sep)
	{
		ChargeLog(slog1, prod1, *A);
		ChargeLog(slog2, prod2, *B);
	}

	CloseLog(slog1, prod1);
	CloseLog(slog2, prod2);
	return prod1 - prod2;
}

/* log(prod(A[i += sep] / B[i += sep])) */
TARGET double LogProdDiv(float* A, float* B, int64 n, int64 sep)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return LogProdDivNEO(A, B, n, sep);
#else
	case 4: return LogProdDiv512(A, B, n, sep);
	case 3: return LogProdDivAVX(A, B, n, sep);
	case 2: return LogProdDivSSE(A, B, n, sep);
#endif
	}

	int64 slog1 = 0; double prod1 = 1;
	int64 slog2 = 0; double prod2 = 1;
	for (int64 i = 0; i < n; ++i, A += sep, B += sep)
	{
		ChargeLog(slog1, prod1, *A);
		ChargeLog(slog2, prod2, *B);
	}

	CloseLog(slog1, prod1);
	CloseLog(slog2, prod2);

	return prod1 - prod2;
}

/* Count non-zero elements */
TARGET int64 CountNonZero(byte* A, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return CountNonZeroNEO(A, n);
#else
	case 4: return CountNonZero512(A, n);
	case 3: return CountNonZeroAVX(A, n);
	case 2: return CountNonZeroSSE(A, n);
#endif
	}

	int64 re = 0;
	for (int64 i = 0; i < n; ++i)
		if (A[i]) re++;
	return re;
}

/* Sum of A */
TARGET double Sum(double* A, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumNEO(A, n);
#else
	case 4: return Sum512(A, n);
	case 3: return SumAVX(A, n);
	case 2: return SumSSE(A, n);
#endif
	}

	double re = 0;
	for (int64 i = 0; i < n; ++i)
		re += A[i];
	return re;
}

/* Sum of A */
TARGET double Sum(float* A, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumNEO(A, n);
#else
	case 4: return Sum512(A, n);
	case 3: return SumAVX(A, n);
	case 2: return SumSSE(A, n);
#endif
	}

	double re = 0;
	for (int64 i = 0; i < n; ++i)
		re += A[i];
	return re;
}

/* Sum of A */
TARGET float Sumx(float* A, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumNEOx(A, n);
#else
	case 4: return Sum512x(A, n);
	case 3: return SumAVXx(A, n);
	case 2: return SumSSEx(A, n);
#endif
	}

	float re = 0;
	for (int64 i = 0; i < n; ++i)
		re += A[i];
	return re;
}

/* Sum of A */
TARGET int64 Sum(byte* A, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumNEO(A, n);
#else
	case 4: return Sum512(A, n);
	case 3: return SumAVX(A, n);
	case 2: return SumSSE(A, n);
#endif
	}

	uint64 re = 0;
	for (int64 i = 0; i < n; ++i)
		re += A[i];
	return re;
}

/* re += A[i += sep] */
TARGET double Sum(double* A, int64 n, int64 sep)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumNEO(A, n, sep);
#else
	case 4: return Sum512(A, n, sep);
	case 3: return SumAVX(A, n, sep);
	case 2: return SumSSE(A, n, sep);
#endif
	}

	double re = 0;
	for (int64 i = 0; i < n; ++i, A += sep)
		re += *A;
	return re;
}

/* re += A[i += sep] */
TARGET double Sum(float* A, int64 n, int64 sep)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumNEO(A, n, sep);
#else
	case 4: return Sum512(A, n, sep);
	case 3: return SumAVX(A, n, sep);
	case 2: return SumSSE(A, n, sep);
#endif
	}

	double re = 0;
	for (int64 i = 0; i < n; ++i, A += sep)
		re += *A;
	return re;
}

/* re += A[i += sep] */
TARGET float Sumx(float* A, int64 n, int64 sep)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumNEOx(A, n, sep);
#else
	case 4: return Sum512x(A, n, sep);
	case 3: return SumAVXx(A, n, sep);
	case 2: return SumSSEx(A, n, sep);
#endif
	}

	float re = 0;
	for (int64 i = 0; i < n; ++i, A += sep)
		re += *A;
	return re;
}

/* Product of A */
TARGET double Prod(double* A, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return ProdNEO(A, n);
#else
	case 4: return Prod512(A, n);
	case 3: return ProdAVX(A, n);
	case 2: return ProdSSE(A, n);
#endif
	}
	
	double re = 1;
	for (int64 i = 0; i < n; ++i)
		re *= A[i];
	return re;
}

/* Product of A */
TARGET double Prod(float* A, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return ProdNEO(A, n);
#else
	case 4: return Prod512(A, n);
	case 3: return ProdAVX(A, n);
	case 2: return ProdSSE(A, n);
#endif
	}
	
	double re = 1;
	for (int64 i = 0; i < n; ++i)
		re *= A[i];
	return re;
}

/* Product of A */
TARGET float Prodx(float* A, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return ProdNEOx(A, n);
#else
	case 4: return Prod512x(A, n);
	case 3: return ProdAVXx(A, n);
	case 2: return ProdSSEx(A, n);
#endif
	}
	
	float re = 1;
	for (int64 i = 0; i < n; ++i)
		re *= A[i];
	return re;
}

/* re *= A[i += sep] */
TARGET double Prod(double* A, int64 n, int64 sep)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return ProdNEO(A, n, sep);
#else
	case 4: return Prod512(A, n, sep);
	case 3: return ProdAVX(A, n, sep);
	case 2: return ProdSSE(A, n, sep);
#endif
	}
	
	double re = 1;
	for (int64 i = 0; i < n; ++i)
		re *= A[i * sep];
	return re;
}

/* re *= A[i += sep] */
TARGET double Prod(float* A, int64 n, int64 sep)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return ProdNEO(A, n, sep);
#else
	case 4: return Prod512(A, n, sep);
	case 3: return ProdAVX(A, n, sep);
	case 2: return ProdSSE(A, n, sep);
#endif
	}
	
	double re = 1;
	for (int64 i = 0; i < n; ++i)
		re *= A[i * sep];
	return re;
}

/* re *= A[i += sep] */
TARGET float Prodx(float* A, int64 n, int64 sep)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return ProdNEOx(A, n, sep);
#else
	case 4: return Prod512x(A, n, sep);
	case 3: return ProdAVXx(A, n, sep);
	case 2: return ProdSSEx(A, n, sep);
#endif
	}
	
	float re = 1;
	for (int64 i = 0; i < n; ++i)
		re *= A[i * sep];
	return re;
}

/* Sum of squared A */
TARGET double SumSquare(double* A, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumSquareNEO(A, n);
#else
	case 4: return SumSquare512(A, n);
	case 3: return SumSquareAVX(A, n);
	case 2: return SumSquareSSE(A, n);
#endif
	}
	
	double re = 0;
	for (int64 i = 0; i < n; ++i)
		re += A[i] * A[i];
	return re;
}

/* Sum of squared A */
TARGET double SumSquare(float* A, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumSquareNEO(A, n);
#else
	case 4: return SumSquare512(A, n);
	case 3: return SumSquareAVX(A, n);
	case 2: return SumSquareSSE(A, n);
#endif
	}
	
	double re = 0;
	for (int64 i = 0; i < n; ++i)
		re += (double)A[i] * (double)A[i];
	return re;
}

/* Sum of squared A */
TARGET int64 SumSquare(byte* A, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumSquareNEO(A, n);
#else
	case 4: return SumSquare512(A, n);
	case 3: return SumSquareAVX(A, n);
	case 2: return SumSquareSSE(A, n);
#endif
	}
	
	uint64 re = 0;
	for (int64 i = 0; i < n; ++i)
		re += A[i] * A[i];
	return re;
}

/* Sum of A and sum of squared A */
TARGET void SumSumSquare(double* A, int64 n, double& sum, double& sumsq)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumSumSquareNEO(A, n, sum, sumsq);
#else
	case 4: return SumSumSquare512(A, n, sum, sumsq);
	case 3: return SumSumSquareAVX(A, n, sum, sumsq);
	case 2: return SumSumSquareSSE(A, n, sum, sumsq);
#endif
	}
	
	sum = sumsq = 0;
	for (int64 i = 0; i < n; ++i)
	{
		sum += A[i];
		sumsq += A[i] * A[i];
	}
}

/* Sum of A and sum of squared A */
TARGET void SumSumSquare(float* A, int64 n, double& sum, double& sumsq)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumSumSquareNEO(A, n, sum, sumsq);
#else
	case 4: return SumSumSquare512(A, n, sum, sumsq);
	case 3: return SumSumSquareAVX(A, n, sum, sumsq);
	case 2: return SumSumSquareSSE(A, n, sum, sumsq);
#endif
	}
	
	sum = sumsq = 0;
	for (int64 i = 0; i < n; ++i)
	{
		sum += (double)A[i];
		sumsq += (double)A[i] * (double)A[i];
	}
}

/* re = sum(A1[i++] * B[j += sep]) / Sum(A2[i++] * B[j += sep]) */
TARGET double SumProdDiv(double* A1, double* A2, double* B, int64 sep, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumProdDivNEO(A1, A2, B, sep, n);
#else
	case 4: return SumProdDiv512(A1, A2, B, sep, n);
	case 3: return SumProdDivAVX(A1, A2, B, sep, n);
	case 2: return SumProdDivSSE(A1, A2, B, sep, n);
#endif
	}

	double re1 = 0, re2 = 0;
	for (int64 i = 0; i < n; ++i, A1++, A2++, B += sep)
	{
		volatile double v1 = (double)*A1 * (double)*B;
		volatile double v2 = (double)*A2 * (double)*B;
		re1 += v1;
		re2 += v2;
	}
	return re1 / re2;
}

/* re = sum(A1[i++] * B[j += sep]) / Sum(A2[i++] * B[j += sep]) */
TARGET double SumProdDiv(double* A1, float* A2, float* B, int64 sep, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumProdDivNEO(A1, A2, B, sep, n);
#else
	case 4: return SumProdDiv512(A1, A2, B, sep, n);
	case 3: return SumProdDivAVX(A1, A2, B, sep, n);
	case 2: return SumProdDivSSE(A1, A2, B, sep, n);
#endif
	}

	double re1 = 0, re2 = 0;
	for (int64 i = 0; i < n; ++i, A1++, A2++, B += sep)
	{
		volatile double v1 = (double)*A1 * (double)*B;
		volatile double v2 = (double)*A2 * (double)*B;
		re1 += v1;
		re2 += v2;
	}
	return re1 / re2;
}

/* re = sum(A1[i++] * B[j += sep]) / Sum(A2[i++] * B[j += sep]) */
TARGET double SumProdDiv(float* A1, float* A2, float* B, int64 sep, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumProdDivNEO(A1, A2, B, sep, n);
#else
	case 4: return SumProdDiv512(A1, A2, B, sep, n);
	case 3: return SumProdDivAVX(A1, A2, B, sep, n);
	case 2: return SumProdDivSSE(A1, A2, B, sep, n);
#endif
	}

	double re1 = 0, re2 = 0;
	for (int64 i = 0; i < n; ++i, A1++, A2++, B += sep)
	{
		volatile double v1 = (double)*A1 * (double)*B;
		volatile double v2 = (double)*A2 * (double)*B;
		re1 += v1;
		re2 += v2;
	}
	return re1 / re2;
}

/* re = sum(A1[i++] * B[j += sep]) / Sum(A2[i++] * B[j += sep]) */
TARGET float SumProdDivx(float* A1, float* A2, float* B, int64 sep, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumProdDivNEOx(A1, A2, B, sep, n);
#else
	case 4: return SumProdDiv512x(A1, A2, B, sep, n);
	case 3: return SumProdDivAVXx(A1, A2, B, sep, n);
	case 2: return SumProdDivSSEx(A1, A2, B, sep, n);
#endif
	}

	float re1 = 0, re2 = 0;
	for (int64 i = 0; i < n; ++i, A1++, A2++, B += sep)
	{
		volatile float v1 = (float)*A1 * (float)*B;
		volatile float v2 = (float)*A2 * (float)*B;
		re1 += v1;
		re2 += v2;
	}
	return re1 / re2;
}

/* re = sum(A[i++] * B[j += sep]) */
TARGET double SumProd(double* A, double* B, int64 sep, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumProdNEO(A, B, sep, n);
#else
	case 4: return SumProd512(A, B, sep, n);
	case 3: return SumProdAVX(A, B, sep, n);
	case 2: return SumProdSSE(A, B, sep, n);
#endif
	}

	double re = 0;
	for (int64 i = 0; i < n; ++i, A++, B += sep)
	{
		volatile double v1 = *A * *B;
		re += v1;
	}
	return re;
}

/* re = sum(A[i++] * B[j += sep]) */
TARGET double SumProd(float* A, float* B, int64 sep, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumProdNEO(A, B, sep, n);
#else
	case 4: return SumProd512(A, B, sep, n);
	case 3: return SumProdAVX(A, B, sep, n);
	case 2: return SumProdSSE(A, B, sep, n);
#endif
	}

	double re = 0;
	for (int64 i = 0; i < n; ++i, A++, B += sep)
		re += (double)*A * (double)*B;
	return re;
}

/* re = sum(A[i++] * B[j += sep]) */
TARGET float SumProdx(float* A, float* B, int64 sep, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumProdNEOx(A, B, sep, n);
#else
	case 4: return SumProd512x(A, B, sep, n);
	case 3: return SumProdAVXx(A, B, sep, n);
	case 2: return SumProdSSEx(A, B, sep, n);
#endif
	}

	float re = 0;
	for (int64 i = 0; i < n; ++i, A++, B += sep)
	{
		volatile float v1 = *A * *B;
		re += v1;
	}
	return re;
}

/* re = sum(A[i] * B[i]) */
TARGET double SumProd(double* A, double* B, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumProdNEO(A, B, n);
#else
	case 4: return SumProd512(A, B, n);
	case 3: return SumProdAVX(A, B, n);
	case 2: return SumProdSSE(A, B, n);
#endif
	}

	double re = 0;
	for (int64 i = 0; i < n; ++i)
		re += A[i] * B[i];
	return re;
}

/* re = sum(A[i] * B[i]) */
TARGET double SumProd(float* A, float* B, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumProdNEO(A, B, n);
#else
	case 4: return SumProd512(A, B, n);
	case 3: return SumProdAVX(A, B, n);
	case 2: return SumProdSSE(A, B, n);
#endif
	}

	double re = 0;
	for (int64 i = 0; i < n; ++i)
		re += (double)A[i] * (double)B[i];
	return re;
}

/* re = sum(A[i] * B[i]) */
TARGET float SumProdx(float* A, float* B, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumProdNEOx(A, B, n);
#else
	case 4: return SumProd512x(A, B, n);
	case 3: return SumProdAVXx(A, B, n);
	case 2: return SumProdSSEx(A, B, n);
#endif
	}
	
	float re = 0;
	for (int64 i = 0; i < n; ++i)
	{
		volatile float v1 = A[i] * B[i];
		re += v1;
	}
	return re;
}

/* re = sum(A[i] * B[i] * C[i]) */
TARGET double SumProd(double* A, double* B, double* C, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumProdNEO(A, B, C, n);
#else
	case 4: return SumProd512(A, B, C, n);
	case 3: return SumProdAVX(A, B, C, n);
	case 2: return SumProdSSE(A, B, C, n);
#endif
	}

	double re = 0;
	for (int64 i = 0; i < n; ++i)
		re += A[i] * B[i] * C[i];
	return re;
}

/* re = sum(A[i] * B[i] * C[i]) */
TARGET float SumProd(float* A, float* B, float* C, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumProdNEO(A, B, C, n);
#else
	case 4: return SumProd512(A, B, C, n);
	case 3: return SumProdAVX(A, B, C, n);
	case 2: return SumProdSSE(A, B, C, n);
#endif
	}

	float re = 0;
	for (int64 i = 0; i < n; ++i)
	{
		volatile float v1 = A[i] * B[i] * C[i];
		re += v1;
	}
	return re;
}

/* re = sum(A[i] * A[i] * B[i]) */
TARGET double SumSqProd(double* A, double* B, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumSqProdNEO(A, B, n);
#else
	case 4: return SumSqProd512(A, B, n);
	case 3: return SumSqProdAVX(A, B, n);
	case 2: return SumSqProdSSE(A, B, n);
#endif
	}

	double re = 0;
	for (int64 i = 0; i < n; ++i)
		re += A[i] * A[i] * B[i];
	return re;
}

/* re = sum(A[i] * A[i] * B[i]) */
TARGET float SumSqProd(float* A, float* B, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return SumSqProdNEO(A, B, n);
#else
	case 4: return SumSqProd512(A, B, n);
	case 3: return SumSqProdAVX(A, B, n);
	case 2: return SumSqProdSSE(A, B, n);
#endif
	}

	float re = 0;
	for (int64 i = 0; i < n; ++i)
	{
		volatile float v1 = A[i] * A[i] * B[i];
		re += v1;
	}
	return re;
}

/* Add B into A, A[i] += B[i] */
TARGET void Add(double* A, double* B, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return AddNEO(A, B, n);
#else
	case 4: return Add512(A, B, n);
	case 3: return AddAVX(A, B, n);
	case 2: return AddSSE(A, B, n);
#endif
	}
	
	for (int64 i = 0; i < n; ++i)
		A[i] += B[i];
}

/* Add B into A, A[i] += B[i] */
TARGET void Add(float* A, float* B, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return AddNEO(A, B, n);
#else
	case 4: return Add512(A, B, n);
	case 3: return AddAVX(A, B, n);
	case 2: return AddSSE(A, B, n);
#endif
	}

	for (int64 i = 0; i < n; ++i)
		A[i] += B[i];
}

/* Add B into A, A[i] += B[i] */
TARGET void Add(int64* A, int64* B, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return AddNEO(A, B, n);
#else
	case 4: return Add512(A, B, n);
	case 3: return AddAVX(A, B, n);
	case 2: return AddSSE(A, B, n);
#endif
	}

	for (int64 i = 0; i < n; ++i)
		A[i] += B[i];
}

/* Add B into A, A[i] += B[i] */
TARGET void Add(int* A, int* B, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return AddNEO(A, B, n);
#else
	case 4: return Add512(A, B, n);
	case 3: return AddAVX(A, B, n);
	case 2: return AddSSE(A, B, n);
#endif
	}

	for (int64 i = 0; i < n; ++i)
		A[i] += B[i];
}

/* Add B into A, A[i] += B */
TARGET void Add(int* A, int B, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return AddNEO(A, B, n);
#else
	case 4: return Add512(A, B, n);
	case 3: return AddAVX(A, B, n);
	case 2: return AddSSE(A, B, n);
#endif
	}

	for (int64 i = 0; i < n; ++i)
		A[i] += B;
}

/* Add B into A, A[i] += B */
TARGET void Add(double* A, double B, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return AddNEO(A, B, n);
#else
	case 4: return Add512(A, B, n);
	case 3: return AddAVX(A, B, n);
	case 2: return AddSSE(A, B, n);
#endif
	}

	for (int64 i = 0; i < n; ++i)
		A[i] += B;
}

/* Add B into A, A[i] += B */
TARGET void Add(float* A, float B, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return AddNEO(A, B, n);
#else
	case 4: return Add512(A, B, n);
	case 3: return AddAVX(A, B, n);
	case 2: return AddSSE(A, B, n);
#endif
	}

	for (int64 i = 0; i < n; ++i)
		A[i] += B;
}

/* A[i] = B[i] * C[i] */
TARGET void Mul(double* A, double* B, double* C, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return MulNEO(A, B, C, n);
#else
	case 4: return Mul512(A, B, C, n);
	case 3: return MulAVX(A, B, C, n);
	case 2: return MulSSE(A, B, C, n);
#endif
	}

	for (int64 i = 0; i < n; ++i)
		A[i] = B[i] * C[i];
}

/* A[i] = B[i] * C[i] */
TARGET void Mul(float* A, float* B, float* C, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return MulNEO(A, B, C, n);
#else
	case 4: return Mul512(A, B, C, n);
	case 3: return MulAVX(A, B, C, n);
	case 2: return MulSSE(A, B, C, n);
#endif
	}

	for (int64 i = 0; i < n; ++i)
		A[i] = B[i] * C[i];
}

/* A[i] = B[i] * C */
TARGET void Mul(double* A, double* B, double C, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return MulNEO(A, B, C, n);
#else
	case 4: return Mul512(A, B, C, n);
	case 3: return MulAVX(A, B, C, n);
	case 2: return MulSSE(A, B, C, n);
#endif
	}

	for (int64 i = 0; i < n; ++i)
		A[i] = B[i] * C;
}

/* A[i] = B[i] * C */
TARGET void Mul(float* A, float* B, float C, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return MulNEO(A, B, C, n);
#else
	case 4: return Mul512(A, B, C, n);
	case 3: return MulAVX(A, B, C, n);
	case 2: return MulSSE(A, B, C, n);
#endif
	}

	for (int64 i = 0; i < n; ++i)
		A[i] = B[i] * C;
}

/* A[i] *= B */
TARGET void Mul(double* A, double B, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return MulNEO(A, B, n);
#else
	case 4: return Mul512(A, B, n);
	case 3: return MulAVX(A, B, n);
	case 2: return MulSSE(A, B, n);
#endif
	}
	
	for (int64 i = 0; i < n; ++i)
		A[i] *= B;
}

/* A[i] *= B */
TARGET void Mul(float* A, float B, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return MulNEO(A, B, n);
#else
	case 4: return Mul512(A, B, n);
	case 3: return MulAVX(A, B, n);
	case 2: return MulSSE(A, B, n);
#endif
	}
	
	for (int64 i = 0; i < n; ++i)
		A[i] *= B;
}

/* A[i] = B / C[i] */
TARGET void Div(double* A, double B, double* C, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return DivNEO(A, B, C, n);
#else
	case 4: return Div512(A, B, C, n);
	case 3: return DivAVX(A, B, C, n);
	case 2: return DivSSE(A, B, C, n);
#endif
	}
	
	for (int64 i = 0; i < n; ++i)
		A[i] = B / C[i];
}

/* A[i] = B / C[i] */
TARGET void Div(float* A, float B, float* C, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return DivNEO(A, B, C, n);
#else
	case 4: return Div512(A, B, C, n);
	case 3: return DivAVX(A, B, C, n);
	case 2: return DivSSE(A, B, C, n);
#endif
	}

	for (int64 i = 0; i < n; ++i)
		A[i] = B / C[i];
}

/* A[i] = B[i] / C[i] */
TARGET void Div(double* A, double* B, double* C, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return DivNEO(A, B, C, n);
#else
	case 4: return Div512(A, B, C, n);
	case 3: return DivAVX(A, B, C, n);
	case 2: return DivSSE(A, B, C, n);
#endif
	}

	for (int64 i = 0; i < n; ++i)
		A[i] = B[i] / C[i];
}

/* A[i] = B[i] / C[i] */
TARGET void Div(float* A, float* B, float* C, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return DivNEO(A, B, C, n);
#else
	case 4: return Div512(A, B, C, n);
	case 3: return DivAVX(A, B, C, n);
	case 2: return DivSSE(A, B, C, n);
#endif
	}

	for (int64 i = 0; i < n; ++i)
		A[i] = B[i] / C[i];
}

/* A[i] += B[i] * C[i] */
TARGET void AddProd(double* A, double* B, double* C, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return AddProdNEO(A, B, C, n);
#else
	case 4: return AddProd512(A, B, C, n);
	case 3: return AddProdAVX(A, B, C, n);
	case 2: return AddProdSSE(A, B, C, n);
#endif
	}

	for (int64 i = 0; i < n; ++i)
		A[i] += B[i] * C[i];
}

/* A[i] += B[i] * C[i] */
TARGET void AddProd(float* A, float* B, float* C, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return AddProdNEO(A, B, C, n);
#else
	case 4: return AddProd512(A, B, C, n);
	case 3: return AddProdAVX(A, B, C, n);
	case 2: return AddProdSSE(A, B, C, n);
#endif
	}

	for (int64 i = 0; i < n; ++i)
		A[i] += B[i] * C[i];
}

/* A[i] += B[i] * C */
TARGET void AddProd(double* A, double* B, double C, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return AddProdNEO(A, B, C, n);
#else
	case 4: return AddProd512(A, B, C, n);
	case 3: return AddProdAVX(A, B, C, n);
	case 2: return AddProdSSE(A, B, C, n);
#endif
	}

	for (int64 i = 0; i < n; ++i)
		A[i] += B[i] * C;
}

/* A[i] += B[i] * C */
TARGET void AddProd(double* A, float* B, double C, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return AddProdNEO(A, B, C, n);
#else
	case 4: return AddProd512(A, B, C, n);
	case 3: return AddProdAVX(A, B, C, n);
	case 2: return AddProdSSE(A, B, C, n);
#endif
	}

	for (int64 i = 0; i < n; ++i)
		A[i] += B[i] * C;
}

/* A[i] += B[i] * C */
TARGET void AddProd(float* A, float* B, float C, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return AddProdNEO(A, B, C, n);
#else
	case 4: return AddProd512(A, B, C, n);
	case 3: return AddProdAVX(A, B, C, n);
	case 2: return AddProdSSE(A, B, C, n);
#endif
	}

	for (int64 i = 0; i < n; ++i)
		A[i] += B[i] * C;
}

/* Set the sum of A to one */
TARGET void Unify(double* A, int64 n)
{
	switch (SIMD_TYPE)
	//A[i] = A[i] * invs + MIN_FREQ * invs
	{
#ifdef __aarch64__
	case 2: return UnifyNEO(A, n);
#else
	case 4: return Unify512(A, n);
	case 3: return UnifyAVX(A, n);
	case 2: return UnifySSE(A, n);
#endif
	}

	double invsum = 1.0 / (Sum(A, n) + n * MIN_FREQ);
	for (int64 i = 0; i < n; ++i)
		A[i] = (A[i] + MIN_FREQ) * invsum;
}

/* Set the sum of A to one */
TARGET void Unify(float* A, int64 n)
{
	switch (SIMD_TYPE)
		//A[i] = A[i] * invs + MIN_FREQ * invs
	{
#ifdef __aarch64__
	case 2: return UnifyNEO(A, n);
#else
	case 4: return Unify512(A, n);
	case 3: return UnifyAVX(A, n);
	case 2: return UnifySSE(A, n);
#endif
	}

	double invsum = 1.0 / (Sum(A, n) + n * MIN_FREQ);
	for (int64 i = 0; i < n; ++i)
		A[i] = ((double)A[i] + MIN_FREQ) * invsum;
}

/* Find next position of val in string A*/
TARGET char* StrNextIdx(char* A, char val, int64 rep, int64 n)
{
	if (!n) return NULL;
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return StrNextIdxNEO(A, val, rep, n);
#else
	case 4: return StrNextIdx512(A, val, rep, n);
	case 3: return StrNextIdxAVX(A, val, rep, n);
	case 2: return StrNextIdxSSE(A, val, rep, n);
#endif
	}

	A++; n--;
	for (int64 i = 0; i < n; ++i, A++)
		if (*A == val && !--rep)
			return A;
	return NULL;
}

/* Count val in string A */
TARGET int64 CountChar(char* A, char val, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return CountCharNEO(A, val, n);
#else
	case 4: return CountChar512(A, val, n);
	case 3: return CountCharAVX(A, val, n);
	case 2: return CountCharSSE(A, val, n);
#endif
	}

	int64 re = 0;
	for (int64 i = 0; i < n; ++i)
		if (A[i] == val) 
			re++;
	return re;
}

/* Quadratic form A' D A with D being a diagonal matrix, A is m*n, D is n*n, ColMajor */
TARGET void DiagQuadForm(double* res, double* A, double* D, int64 m, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return DiagQuadFormNEO(res, A, D, m, n);
#else
	case 4: return DiagQuadForm512(res, A, D, m, n);
	case 3: return DiagQuadFormAVX(res, A, D, m, n);
	case 2: return DiagQuadFormSSE(res, A, D, m, n);
#endif
	}
	
	Mat<double> a(A, n, m, false, true);
	Col<double> d(D, n, false, true);
	Mat<double> r(res, m, m, false, true);
	
	r = a.t() * diagmat(d) * a;
}

/* Quadratic form A' D A with D being a diagonal matrix, A is m*n, D is n*n, ColMajor */
TARGET void DiagQuadForm(float* res, float* A, float* D, int64 m, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return DiagQuadFormNEO(res, A, D, m, n);
#else
	case 4: return DiagQuadForm512(res, A, D, m, n);
	case 3: return DiagQuadFormAVX(res, A, D, m, n);
	case 2: return DiagQuadFormSSE(res, A, D, m, n);
#endif
	}
	
	Mat<float> a(A, n, m, false, true);
	Col<float> d(D, n, false, true);
	Mat<float> r(res, m, m, false, true);
	
	r = a.t() * diagmat(d) * a;
}

/* Quadratic form A' D B with D being a diagonal matrix, A is m*n, D is n*n, B is n*1, ColMajor */
TARGET void DiagQuadForm(double* res, double* A, double* D, double* B, int64 m, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return DiagQuadFormNEO(res, A, D, B, m, n);
#else
	case 4: return DiagQuadForm512(res, A, D, B, m, n);
	case 3: return DiagQuadFormAVX(res, A, D, B, m, n);
	case 2: return DiagQuadFormSSE(res, A, D, B, m, n);
#endif
	}
	
	Mat<double> a(A, n, m, false, true);
	Col<double> d(D, n, false, true);
	Mat<double> b(B, n, 1, false, true);
	Mat<double> r(res, m, 1, false, true);
	
	r = a.t() * diagmat(d) * b;
}

/* Quadratic form A D B with D being a diagonal matrix, A is m*n, D is n*n, B is n*1, ColMajor */
TARGET void DiagQuadForm(float* res, float* A, float* D, float* B, int64 m, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return DiagQuadFormNEO(res, A, D, B, m, n);
#else
	case 4: return DiagQuadForm512(res, A, D, B, m, n);
	case 3: return DiagQuadFormAVX(res, A, D, B, m, n);
	case 2: return DiagQuadFormSSE(res, A, D, B, m, n);
#endif
	}
	
	Mat<float> a(A, n, m, false, true);
	Col<float> d(D, n, false, true);
	Mat<float> b(B, n, 1, false, true);
	Mat<float> r(res, m, 1, false, true);
	
	r = a.t() * diagmat(d) * b;
}

/* Quadratic form A D A' with D being a diagonal matrix, A is 1*n, D is n*n, ColMajor */
TARGET void DiagQuadForm(double* res, double* A, double* D, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return DiagQuadFormNEO(res, A, D, n);
#else
	case 4: return DiagQuadForm512(res, A, D, n);
	case 3: return DiagQuadFormAVX(res, A, D, n);
	case 2: return DiagQuadFormSSE(res, A, D, n);
#endif
	}
	
	Mat<double> a(A, n, 1, false, true);
	Col<double> d(D, n, false, true);
	Mat<double> r(res, 1, 1, false, true);
	
	r = a.t() * diagmat(d) * a;
}

/* Quadratic form A D A' with D being a diagonal matrix, A is 1*n, D is n*n, ColMajor */
TARGET void DiagQuadForm(float* res, float* A, float* D, int64 n)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return DiagQuadFormNEO(res, A, D, n);
#else
	case 4: return DiagQuadForm512(res, A, D, n);
	case 3: return DiagQuadFormAVX(res, A, D, n);
	case 2: return DiagQuadFormSSE(res, A, D, n);
#endif
	}
	
	Mat<float> a(A, n, 1, false, true);
	Col<float> d(D, n, false, true);
	Mat<float> r(res, 1, 1, false, true);
	
	r = a.t() * diagmat(d) * a;
}

/* Matrix Muplification, A is m*n, B is n*p, ColMajor */
TARGET void MatrixMul(double* res, double* A, double* B, int64 m, int64 n, int64 p)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return MatrixMulNEO(res, A, B, m, n, p);
#else
	case 4: return MatrixMul512(res, A, B, m, n, p);
	case 3: return MatrixMulAVX(res, A, B, m, n, p);
	case 2: return MatrixMulSSE(res, A, B, m, n, p);
#endif
	}
	
	Mat<double> a(A, n, m, false, true);
	Mat<double> b(B, n, p, false, true);
	Mat<double> r(res, m, p, false, true);

	r = a.t() * b;
}

/* Matrix Muplification, A is m*n, B is n*p, ColMajor */
TARGET void MatrixMul(float* res, float* A, float* B, int64 m, int64 n, int64 p)
{
	switch (SIMD_TYPE)
	{
#ifdef __aarch64__
	case 2: return MatrixMulNEO(res, A, B, m, n, p);
#else
	case 4: return MatrixMul512(res, A, B, m, n, p);
	case 3: return MatrixMulAVX(res, A, B, m, n, p);
	case 2: return MatrixMulSSE(res, A, B, m, n, p);
#endif
	}
	
	Mat<float> a(A, n, m, false, true);
	Mat<float> b(B, n, p, false, true);
	Mat<float> r(res, m, p, false, true);

	r = a.t() * b;
}
