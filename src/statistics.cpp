/* Statistics Functions */

#include "vcfpop.h"
#include "gcem.hpp"
#include "stats.hpp"
#include <limits>

using std::numeric_limits;

/*
template<typename REAL>
struct numeric_constants
{
	static REAL pi = static_cast<REAL>(3.1415926535897932384626433832795029L);
	static REAL pi_2 = static_cast<REAL>(1.5707963267948966192313216916397514L);
	static REAL pi_3 = static_cast<REAL>(1.0471975511965977461542144610931676L);
	static REAL pi_4 = static_cast<REAL>(0.7853981633974483096156608458198757L);
	static REAL 1_pi = static_cast<REAL>(0.3183098861837906715377675267450287L);
	static REAL 2_sqrtpi = static_cast<REAL>(1.1283791670955125738961589031215452L);
	static REAL sqrt2 = static_cast<REAL>(1.4142135623730950488016887242096981L);
	static REAL sqrt3 = static_cast<REAL>(1.7320508075688772935274463415058723L);
	static REAL sqrtpio2 = static_cast<REAL>(1.2533141373155002512078826424055226L);
	static REAL sqrt1_2 = static_cast<REAL>(0.7071067811865475244008443621048490L);
	static REAL lnpi = static_cast<REAL>(1.1447298858494001741434273513530587L);
	static REAL gamma_e = static_cast<REAL>(0.5772156649015328606065120900824024L);
	static REAL euler = static_cast<REAL>(2.7182818284590452353602874713526625L);
};
*/

template struct RNG<double>;
template struct RNG<float >;

template TARGET void RNG<double>::Permute<VESSEL<double>*>(VESSEL<double>** a, int n);
template TARGET void RNG<double>::Permute<VESSEL<float >*>(VESSEL<float >** a, int n);
template TARGET void RNG<double>::Permute<int64>(int64* a, int64 n);
template TARGET void RNG<double>::Permute<ushort>(ushort* a, int n);
template TARGET void RNG<double>::Permute<int>(int* a, int n);
template TARGET void RNG<double>::GetRandSeq<int64>(int64* a, int64 n);
template TARGET void RNG<double>::Dirichlet<double, double>(double* res, double* a, int n, double f);
template TARGET void RNG<double>::Dirichlet<double, float >(double* res, float * a, int n, double f);
template TARGET void RNG<double>::Dirichlet<float , double>(float * res, double* a, int n, double f);
template TARGET void RNG<double>::Dirichlet<float , float >(float * res, float * a, int n, double f);
template TARGET void RNG<double>::Dirichlet<double, double, int  >(double* res, double* a, int  * b, int n);
template TARGET void RNG<double>::Dirichlet<double, double, int64>(double* res, double* a, int64* b, int n);
template TARGET void RNG<double>::Dirichlet<float , float , int  >(float * res, float * a, int  * b, int n);
template TARGET void RNG<double>::Dirichlet<float , double, int64>(float * res, double* a, int64* b, int n);
template TARGET void RNG<double>::Uniform<double>(double* arr, int n, double min, double max);
template TARGET void RNG<double>::Uniform<float >(float * arr, int n, double min, double max);
template TARGET void RNG<double>::Normal<double>(double* arr, int n, double mean, double std);
template TARGET void RNG<double>::Normal<float >(float * arr, int n, double mean, double std);
template TARGET void RNG<double>::Integer<uint64>(uint64* arr, int n, uint64 mean, uint64 std);
template TARGET void RNG<double>::Integer<uint  >(uint  * arr, int n, uint64 mean, uint64 std);


#ifndef _RNG_FP64

TARGET RNG<double>::RNG()
{

}

TARGET RNG<double>::RNG(uint64 s, uint64 salt)
{
	s = MurmurHash64(s, salt);

	x = 0x159A55E5075BCD15 ^ (s);       //123456789, 362436069
	y = 0x054913331F123BB5 ^ (s << 6);  //521288629, 88675123

	*(uint*)&U1 = 0xFFFFFFFF;
}

/* Draw a uniform distriubted interger */
TARGET uint64 RNG<double>::XorShift()
{
	//XorShift
	uint64 a = x, b = y;

	x = b;
	a ^= a << 23;
	a ^= a >> 18;
	a ^= b;
	a ^= b >> 5;
	y = a;

	return a + b;
}

/* Draw a polynormial distriubted integer */
TARGET int RNG<double>::Poly(double* a, int n, int sep)
{
	//row unify
	if (sep == 1)
	{
		double s = Sum(a, n) + MIN_FREQ * n;
		double t = Uniform(0, s);
		for (int i = 0; i < n; ++i)
		{
			if (t < a[i] + MIN_FREQ) return i;
			t -= a[i] + MIN_FREQ;
		}
	}
	else
	{
		double s = Sum(a, n, sep) + MIN_FREQ * n;
		double t = Uniform(0, s);
		for (int i = 0; i < n; ++i)
		{
			if (t < a[i * sep] + MIN_FREQ) return i;
			t -= a[i * sep] + MIN_FREQ;
		}
	}
	return n - 1;
}

/* Draw a polynormial distriubted integer with propoirtions in natural logarithm */
TARGET int RNG<double>::PolyLog(double* a, int n, int sep)
{
	if (sep == 1)
	{
		//proportional polynomial distribution, will overwrite a
		double maxval = GetMaxVal(a, n), s = MIN_FREQ * n;
		for (int i = 0; i < n; ++i)
		{
			double diff = a[i] - maxval;
			a[i] = diff < -23 ? MIN_FREQ : exp(diff);
			s += a[i];
		}

		double r = Uniform(0, s);
		for (int i = 0; i < n; ++i)
		{
			if (r < a[i]) return i;
			r -= a[i];
		}
	}
	else
	{
		//proportional polynomial distribution, will overwrite a
		double maxval = GetMaxVal(a, n, sep), s = MIN_FREQ * n;
		for (int i = 0; i < n; ++i)
		{
			double diff = a[i * sep] - maxval;
			a[i * sep] = diff < -23 ? MIN_FREQ : exp(diff);
			s += a[i * sep];
		}

		double r = Uniform(0, s);
		for (int i = 0; i < n; ++i)
		{
			if (r < a[i * sep]) return i;
			r -= a[i * sep];
		}
	}
	return n - 1;
}

/* Get a random sequence from 0 ~ n-1 */
template<typename T>
TARGET void RNG<double>::GetRandSeq(T* td, T n)
{
	for (T i = 0; i < n; ++i)
		td[i] = i;

	Permute(td, n);
}

/* Draw a uniform distriubted real number */
TARGET double RNG<double>::Uniform(double min, double max)
{
	uint64 u = XorShift(), r = 0x3FF0000000000000;
	double& re = *(double*)&r;
	r |= u & 0x000FFFFFFFFFFFFF;
	return (re - 1.0) * (max - min) + min;;
}

/* Draw uniform distriubted real numbers */
template<typename T1>
TARGET void RNG<double>::Uniform(T1* arr, int n, double min, double max)
{
	for (int i = 0; i < n; ++i)
		arr[i] = (T1)Uniform(min, max);
}

/* Draw a normal distriubted real number */
TARGET double RNG<double>::Normal(double mean, double std)
{
	//normal distribution
	if (*(uint*)&U1 != 0xFFFFFFFF)
	{
		double re = U1 * MySin(U2);
		*(uint*)&U1 = 0xFFFFFFFF;
		return re * std + mean;
	}

	volatile double v1 = Uniform();
	volatile double v2 = Uniform();

	U1 = MySqrt(-2.0 * log(v1));
	U2 = 2.0 * M_PI * v2;
	return U1 * MyCos(U2) * std + mean;
}

/* Draw normal distriubted real numbers */
template<typename T1>
TARGET void RNG<double>::Normal(T1* arr, int n, double mean, double std)
{
	for (int i = 0; i < n; ++i)
		arr[i] = (T1)Normal(mean, std);
}

/* Draw a real number from gamma distribution */
TARGET double RNG<double>::Gamma(double alpha, double beta)
{
	//gamma distribution
	if (alpha < 1)
	{
		//bug fixed on 20220816 to keep code sequence
		volatile double v1 = Gamma(1.0 + alpha, beta);
		volatile double v2 = pow(Uniform(), 1.0 / alpha);
		return v1 * v2;
	}
	double t, v, u;
	double d = alpha - 0.333333333333333;
	double c = 1.0 / (sqrt(d) * 3.0);

	for (;;)
	{
		do
		{
			t = Normal();
			volatile double v1 = c * t;
			v = 1.0 + v1;
		} while (v <= 0);

		v = v * v * v;
		u = Uniform();

		if (u < 1.0 - 0.0331 * t * t * t * t) break;
		if (log(u) < 0.5 * t * t + d * (1.0 - v + log(v))) break;
	}
	return beta * d * v;
}

/* Draw a real number from beta distribution */
TARGET double RNG<double>::Beta(double a, double b)
{
	//Beta distribution
	volatile double v1 = Gamma(a);
	volatile double v2 = Gamma(b);
	return v1 / (v1 + v2);
}

/* Draw a vector from Dirichlet distribution D(a1 f, a2 f, ...) */
template<typename T1, typename T2>
TARGET void RNG<double>::Dirichlet(T1* res, T2* a, int n, double f)
{
	//Dirichlet distribution
	double s = 0;
	for (int i = 0; i < n; ++i)
	{
		double v = Gamma(a[i] * f);
		res[i] = v;
		s += v;
	}
	Mul(res, (T1)(1.0 / s), n);
}

/* Draw a vector from Dirichlet distribution D(a1 + b1, a2 + b2, ...) */
template<typename T1, typename T2, typename T3>
TARGET void RNG<double>::Dirichlet(T1* res, T2* a, T3* b, int n)
{
	//Dirichlet distribution
	double s = 0;
	for (int i = 0; i < n; ++i)
	{
		double v = Gamma((double)a[i] + (double)b[i]);
		res[i] = v;
		s += v;
	}
	Mul(res, (T1)(1.0 / s), n);
}

/* Shuffle an array */
template<typename T, typename T2>
TARGET void RNG<double>::Permute(T* val, T2 n)
{
	//https://lemire.me/blog/2016/06/30/fast-random-shuffling/
	for (T2 i = n; i > 1; --i)
		Swap(val[i - 1], val[Next(i)]);
}

/* Draw uniform distriubted intergers */
template<typename T1>
TARGET void RNG<double>::Integer(T1* arr, int n, uint64 min, uint64 max)
{
	for (int i = 0; i < n; ++i)
		arr[i] = (T1)Next(min, max);
}

/* Draw a uniform distriubted interger */
TARGET uint64 RNG<double>::Next(uint64 min, uint64 max)
{
	// will not equal to max
	return XorShift() % (max - min) + min;
}

/* Draw a uniform distriubted interger */
TARGET uint64 RNG<double>::Next(uint64 max)
{
	// will not equal to max
	return XorShift() % max;
}

/* Draw a uniform distriubted interger and avoid sample av */
TARGET uint64 RNG<double>::NextAvoid(uint64 max, uint64 av)
{
	uint64 a = Next(max - 1);
	if (a >= av) a++;
	return a;
}
#endif

#ifndef _RNG_FP32
	/* Initialize rng */
	TARGET RNG<float>::RNG()
	{
	}

	/* Initialize rng */
	TARGET RNG<float>::RNG(uint64 seed, uint64 salt)
	{
		uint s = MurmurHash32(seed, salt);

		x = 0x075BCD15 ^ (s);
		y = 0x159A55E5 ^ (s << 3);
		z = 0x1F123BB5 ^ (s << 6);

		//*(uint*)&U1 = 0xFFFFFFFF;
	}

	/* Draw a uniform distriubted interger */
	TARGET uint RNG<float>::XorShift()
	{
		//XorShift96
		uint t;
		x ^= x << 16;
		x ^= x >> 5;
		x ^= x << 1;
		t = x;
		x = y;
		y = z;
		z = t ^ x ^ y;
		return z;
	}

	 /* Draw a polynormial distriubted integer */
	TARGET int RNG<float>::Poly(float* a, int n, int sep)
	{
		if (sep == 1)
		{
			float s = 0;
			s = Sumx(a, n) + MIN_FREQ * n;

			float t = Uniform(0, s);
			for (int i = 0; i < n; ++i)
			{
				if (t < a[i] + MIN_FREQ) return i;
				t -= a[i] + MIN_FREQ;
			}
		}
		else
		{
			float s = 0;
			s = Sumx(a, n, sep) + MIN_FREQ * n;

			float t = Uniform(0, s);
			for (int i = 0; i < n; ++i)
			{
				if (t < a[i * sep] + MIN_FREQ) return i;
				t -= a[i * sep] + MIN_FREQ;
			}
		}
		return n - 1;
	}

	/* Draw a uniform distriubted real number */
	TARGET float RNG<float>::Uniform(float min, float max)
	{
		uint u = XorShift(), r = 0x3F800000;
		float& re = *(float*)&r;
		r |= u & 0x007FFFFF;
		return (re - 1.0f) * (max - min) + min;
	}

#ifdef _xxx
	/* Draw a polynormial distriubted integer with propoirtions in natural logarithm */
	TARGET int RNG<float>::PolyLog(float* a, int n, int sep)
	{
		if (sep == 1)
		{
			//proportional polynomial distribution, will overwrite a
			float maxval = GetMaxVal(a, n), s = (float)MIN_FREQ * n;
			for (int i = 0; i < n; ++i)
			{
				float diff = a[i] - maxval;
				a[i] = diff < -23 ? (REAL)MIN_FREQ : (REAL)exp(diff);
				s += a[i];
			}

			float r = Uniform(0, s);
			for (int i = 0; i < n; ++i)
			{
				if (r < a[i]) return i;
				r -= a[i];
			}
		}
		else
		{
			//proportional polynomial distribution, will overwrite a
			float maxval = GetMaxVal(a, n, sep), s = (float)MIN_FREQ * n;
			for (int i = 0; i < n; ++i)
			{
				float diff = a[i * sep] - maxval;
				a[i * sep] = diff < -23 ? (REAL)MIN_FREQ : (REAL)exp(diff);
				s += a[i * sep];
			}

			float r = Uniform(0, s);
			for (int i = 0; i < n; ++i)
			{
				if (r < a[i * sep]) return i;
				r -= a[i * sep];
			}
		}
		return n - 1;
	}

	/* Get a random sequence from 0 ~ n-1 */
	template<typename T>
	TARGET void RNG<float>::GetRandSeq(T* td, T n)
	{
		for (T i = 0; i < n; ++i)
			td[i] = i;

		Permute(td, n);
		/*
		for (int i = 0; i < n; ++i)
			td[i] = (int)((XorShift() << 16) | i);
		
		//QuickSort(td, 0, n - 1);
		std::sort(td, td + n);

		for (int i = 0; i < n; ++i)
			td[i] &= 0xFFFF;
		*/
	}

	/* Draw uniform distriubted real numbers */
	TARGET void RNG<float>::Uniform(float* arr, int n, float min, float max)
	{
		for (int i = 0; i < n; ++i)
			arr[i] = Uniform(min, max);
	}

	/* Draw a normal distriubted real number */
	TARGET float RNG<float>::Normal(float mean, float std)
	{
		//normal distribution
		if (*(uint*)&U1 != 0xFFFFFFFF)
		{
			float re = U1 * MySin(U2);
			*(uint*)&U1 = 0xFFFFFFFF;
			return re * std + mean;
		}

		volatile float v1 = Uniform();
		volatile float v2 = Uniform();
		U1 = MySqrt(-2.0 * log(std::max(MIN_FREQ, v1)));
		U2 = 2.0 * M_PI * v2;
		return U1 * MyCos(U2) * std + mean;
	}

	/* Draw normal distriubted real numbers */
	TARGET void RNG<float>::Normal(float* arr, int n, float mean, float std)
	{
		for (int i = 0; i < n; ++i)
			arr[i] = Normal(mean, std);
	}

	/* Draw a real number from gamma distribution */
	TARGET float RNG<float>::Gamma(float alpha, float beta)
	{
		//gamma distribution
		if (alpha < 1)
		{
			//bug fixed on 20220816 to keep code sequence
			volatile float v1 = Gamma(1.0 + alpha, beta);
			volatile float v2 = pow(Uniform(), 1.0 / alpha);
			return v1 * v2;
		}
		float t, v, u;
		float d = alpha - 0.333333333333333;
		float c = 1.0 / (sqrt(d) * 3.0);

		for (;;)
		{
			do
			{
				t = Normal();
				v = 1.0 + c * t;
			} while (v <= 0);

			v = v * v * v;
			u = Uniform();

			if (u < 1.0 - 0.0331 * t * t * t * t) break;
			if (log(u) < 0.5 * t * t + d * (1.0 - v + log(v))) break;
		}
		return beta * d * v;
	}

	/* Draw a real number from beta distribution */
	TARGET float RNG<float>::Beta(float a, float b)
	{
		//Beta distribution
		volatile double v1 = Gamma(a);
		volatile double v2 = Gamma(b);
		return v1 / (v1 + v2);
	}

	/* Draw a vector from Dirichlet distribution D(a1 f, a2 f, ...) */
	template<typename T1, typename T2>
	TARGET void RNG<float>::Dirichlet(T1 * res, T2 * a, int n, double f)
	{
		//Dirichlet distribution
		double s = 0;
		for (int i = 0; i < n; ++i)
		{
			double v = Gamma((double)a[i] * f);
			res[i] = (T1)v;
			s += v;
		}
		Mul(res, (T1)(1.0 / s), n);
	}

	/* Draw a vector from Dirichlet distribution D(a1, a2, ...) */
	template<typename T1, typename T2>
	TARGET void RNG<float>::Dirichlet(T1 * res, T2 * a, int n)
	{
		//Dirichlet distribution
		double s = 0;
		for (int i = 0; i < n; ++i)
		{
			double v = Gamma((double)a[i]);
			res[i] = (T1)v;
			s += v;
		}
		Mul(res, (T1)(1.0 / s), n);
	}

	/* Draw a vector from Dirichlet distribution D(a1 + b1, a2 + b2, ...) */
	template<typename T1, typename T2, typename T3>
	TARGET void RNG<float>::Dirichlet(T1 * res, T2 * a, T3 * b, int n)
	{
		//Dirichlet distribution
		double s = 0;
		for (T3 i = 0; i < n; ++i)
		{
			double v = Gamma((double)a[i] + (double)b[i]);
			res[i] = (T1)v;
			s += v;
		}
		Mul(res, (T1)(1.0 / s), n);
	}

	/* Shuffle an array */
	template<typename T, typename T2>
	TARGET void RNG<float>::Permute(T * val, T2 n)
	{
		//https://lemire.me/blog/2016/06/30/fast-random-shuffling/
		for (T2 i = n; i > 1; --i)
			Swap(val[i - 1], val[Next(i)]);
	}

	/* Draw a uniform distriubted interger */
	TARGET uint RNG<float>::Next(uint min, uint max)
	{
		// will not equal to max
		return XorShift() % (max - min) + min;
	}

	/* Draw a uniform distriubted interger */
	TARGET uint RNG<float>::Next(uint max)
	{
		// will not equal to max
		return XorShift() % max;
	}

	/* Draw a uniform distriubted interger and avoid sample av */
	TARGET uint RNG<float>::NextAvoid(uint max, uint av)
	{
		uint a = Next(max - 1);
		if (a >= av) a++;
		return a;
	}
#endif
#endif

namespace GSL
{
	struct ChebSeries {

		double* c;
		size_t order;
		double a;
		double b;
		size_t order_sp;
		double* f;
	};

	static double erfc_xlt1_data[20] = {
	  1.06073416421769980345174155056,
	 -0.42582445804381043569204735291,
	  0.04955262679620434040357683080,
	  0.00449293488768382749558001242,
	 -0.00129194104658496953494224761,
	 -0.00001836389292149396270416979,
	  0.00002211114704099526291538556,
	 -5.23337485234257134673693179020e-7,
	 -2.78184788833537885382530989578e-7,
	  1.41158092748813114560316684249e-8,
	  2.72571296330561699984539141865e-9,
	 -2.06343904872070629406401492476e-10,
	 -2.14273991996785367924201401812e-11,
	  2.22990255539358204580285098119e-12,
	  1.36250074650698280575807934155e-13,
	 -1.95144010922293091898995913038e-14,
	 -6.85627169231704599442806370690e-16,
	  1.44506492869699938239521607493e-16,
	  2.45935306460536488037576200030e-18,
	 -9.29599561220523396007359328540e-19
	};
	static ChebSeries erfc_xlt1_cs = {
	  erfc_xlt1_data,
	  19,
	  -1, 1,
	  12
	};

	static double erfc_x15_data[25] = {
	  0.44045832024338111077637466616,
	 -0.143958836762168335790826895326,
	  0.044786499817939267247056666937,
	 -0.013343124200271211203618353102,
	  0.003824682739750469767692372556,
	 -0.001058699227195126547306482530,
	  0.000283859419210073742736310108,
	 -0.000073906170662206760483959432,
	  0.000018725312521489179015872934,
	 -4.62530981164919445131297264430e-6,
	  1.11558657244432857487884006422e-6,
	 -2.63098662650834130067808832725e-7,
	  6.07462122724551777372119408710e-8,
	 -1.37460865539865444777251011793e-8,
	  3.05157051905475145520096717210e-9,
	 -6.65174789720310713757307724790e-10,
	  1.42483346273207784489792999706e-10,
	 -3.00141127395323902092018744545e-11,
	  6.22171792645348091472914001250e-12,
	 -1.26994639225668496876152836555e-12,
	  2.55385883033257575402681845385e-13,
	 -5.06258237507038698392265499770e-14,
	  9.89705409478327321641264227110e-15,
	 -1.90685978789192181051961024995e-15,
	  3.50826648032737849245113757340e-16
	};
	static ChebSeries erfc_x15_cs = {
	  erfc_x15_data,
	  24,
	  -1, 1,
	  16
	};

	static double erfc_x510_data[20] = {
	  1.11684990123545698684297865808,
	  0.003736240359381998520654927536,
	 -0.000916623948045470238763619870,
	  0.000199094325044940833965078819,
	 -0.000040276384918650072591781859,
	  7.76515264697061049477127605790e-6,
	 -1.44464794206689070402099225301e-6,
	  2.61311930343463958393485241947e-7,
	 -4.61833026634844152345304095560e-8,
	  8.00253111512943601598732144340e-9,
	 -1.36291114862793031395712122089e-9,
	  2.28570483090160869607683087722e-10,
	 -3.78022521563251805044056974560e-11,
	  6.17253683874528285729910462130e-12,
	 -9.96019290955316888445830597430e-13,
	  1.58953143706980770269506726000e-13,
	 -2.51045971047162509999527428316e-14,
	  3.92607828989125810013581287560e-15,
	 -6.07970619384160374392535453420e-16,
	  9.12600607264794717315507477670e-17
	};
	static ChebSeries erfc_x510_cs = {
	  erfc_x510_data,
	  19,
	  -1, 1,
	  12
	};

	static double gstar_a_data[30] = {
	  2.16786447866463034423060819465,
	 -0.05533249018745584258035832802,
	  0.01800392431460719960888319748,
	 -0.00580919269468937714480019814,
	  0.00186523689488400339978881560,
	 -0.00059746524113955531852595159,
	  0.00019125169907783353925426722,
	 -0.00006124996546944685735909697,
	  0.00001963889633130842586440945,
	 -6.3067741254637180272515795142e-06,
	  2.0288698405861392526872789863e-06,
	 -6.5384896660838465981983750582e-07,
	  2.1108698058908865476480734911e-07,
	 -6.8260714912274941677892994580e-08,
	  2.2108560875880560555583978510e-08,
	 -7.1710331930255456643627187187e-09,
	  2.3290892983985406754602564745e-09,
	 -7.5740371598505586754890405359e-10,
	  2.4658267222594334398525312084e-10,
	 -8.0362243171659883803428749516e-11,
	  2.6215616826341594653521346229e-11,
	 -8.5596155025948750540420068109e-12,
	  2.7970831499487963614315315444e-12,
	 -9.1471771211886202805502562414e-13,
	  2.9934720198063397094916415927e-13,
	 -9.8026575909753445931073620469e-14,
	  3.2116773667767153777571410671e-14,
	 -1.0518035333878147029650507254e-14,
	  3.4144405720185253938994854173e-15,
	 -1.0115153943081187052322643819e-15
	};
	static ChebSeries gstar_a_cs = {
	  gstar_a_data,
	  29,
	  -1, 1,
	  17
	};

	static double gstar_b_data[] = {
	  0.0057502277273114339831606096782,
	  0.0004496689534965685038254147807,
	 -0.0001672763153188717308905047405,
	  0.0000615137014913154794776670946,
	 -0.0000223726551711525016380862195,
	  8.0507405356647954540694800545e-06,
	 -2.8671077107583395569766746448e-06,
	  1.0106727053742747568362254106e-06,
	 -3.5265558477595061262310873482e-07,
	  1.2179216046419401193247254591e-07,
	 -4.1619640180795366971160162267e-08,
	  1.4066283500795206892487241294e-08,
	 -4.6982570380537099016106141654e-09,
	  1.5491248664620612686423108936e-09,
	 -5.0340936319394885789686867772e-10,
	  1.6084448673736032249959475006e-10,
	 -5.0349733196835456497619787559e-11,
	  1.5357154939762136997591808461e-11,
	 -4.5233809655775649997667176224e-12,
	  1.2664429179254447281068538964e-12,
	 -3.2648287937449326771785041692e-13,
	  7.1528272726086133795579071407e-14,
	 -9.4831735252566034505739531258e-15,
	 -2.3124001991413207293120906691e-15,
	  2.8406613277170391482590129474e-15,
	 -1.7245370321618816421281770927e-15,
	  8.6507923128671112154695006592e-16,
	 -3.9506563665427555895391869919e-16,
	  1.6779342132074761078792361165e-16,
	 -6.0483153034414765129837716260e-17
	};
	static ChebSeries gstar_b_cs = {
	  gstar_b_data,
	  29,
	  -1, 1,
	  18
	};

	static double lopxmx_data[20] = {
	 -1.12100231323744103373737274541,
	  0.19553462773379386241549597019,
	 -0.01467470453808083971825344956,
	  0.00166678250474365477643629067,
	 -0.00018543356147700369785746902,
	  0.00002280154021771635036301071,
	 -2.8031253116633521699214134172e-06,
	  3.5936568872522162983669541401e-07,
	 -4.6241857041062060284381167925e-08,
	  6.0822637459403991012451054971e-09,
	 -8.0339824424815790302621320732e-10,
	  1.0751718277499375044851551587e-10,
	 -1.4445310914224613448759230882e-11,
	  1.9573912180610336168921438426e-12,
	 -2.6614436796793061741564104510e-13,
	  3.6402634315269586532158344584e-14,
	 -4.9937495922755006545809120531e-15,
	  6.8802890218846809524646902703e-16,
	 -9.5034129794804273611403251480e-17,
	  1.3170135013050997157326965813e-17
	};
	static ChebSeries lopxmx_cs = {
	  lopxmx_data,
	  19,
	  -1, 1,
	  9
	};

	TARGET double ChebEval(ChebSeries* cs, double x)
	{
		int j;
		double y = (2.0 * x - cs->a - cs->b) / (cs->b - cs->a);
		double d = 0.0, dd = 0.0, y2 = 2.0 * y, e = 0.0;

		for (j = (int)cs->order; j >= 1; j--)
		{
			double temp = d;
			d = y2 * d - dd + cs->c[j];
			e += fabs(y2 * temp) + fabs(dd) + fabs(cs->c[j]);
			dd = temp;
		}

		{
			double temp = d;
			d = y * d - dd + 0.5 * cs->c[0];
			e += fabs(y * temp) + fabs(dd) + 0.5 * fabs(cs->c[0]);
		}

		return d;
	}

	/* Evaluate the continued fraction for exprel */
	TARGET double ExprelNCF(double N, double x)
	{
		const double RECUR_BIG = 1.3407807929942596e+154;// SQRT_DBL_MAX;
		const int maxiter = 5000;
		int n = 1;
		double Anm2 = 1.0, Bnm2 = 0.0, Anm1 = 0.0, Bnm1 = 1.0;
		double a1 = 1.0, b1 = 1.0, a2 = -x, b2 = N + 1;
		double an, bn;

		double fn;
		double An = b1 * Anm1 + a1 * Anm2;   /* A1 */
		double Bn = b1 * Bnm1 + a1 * Bnm2;   /* B1 */

		/* One explicit step, before we get to the main pattern. */
		n++;
		Anm2 = Anm1; Bnm2 = Bnm1; Anm1 = An; Bnm1 = Bn;
		An = b2 * Anm1 + a2 * Anm2;   /* A2 */
		Bn = b2 * Bnm1 + a2 * Bnm2;   /* B2 */

		fn = An / Bn;

		while (n < maxiter)
		{
			double old_fn;
			double del;
			n++;
			Anm2 = Anm1;
			Bnm2 = Bnm1;
			Anm1 = An;
			Bnm1 = Bn;
			an = ((n) & 1 ? ((n - 1) / 2) * x : -(N + (n / 2) - 1) * x);
			bn = N + n - 1;
			An = bn * Anm1 + an * Anm2;
			Bn = bn * Bnm1 + an * Bnm2;

			if (fabs(An) > RECUR_BIG || fabs(Bn) > RECUR_BIG) {
				An /= RECUR_BIG;
				Bn /= RECUR_BIG;
				Anm1 /= RECUR_BIG;
				Bnm1 /= RECUR_BIG;
				Anm2 /= RECUR_BIG;
				Bnm2 /= RECUR_BIG;
			}

			old_fn = fn;
			fn = An / Bn;
			del = old_fn / fn;

			if (fabs(del - 1.0) < 2.0 * DBL_EPSILON) break;
		}

		if (n == maxiter)
			return NAN;
		else
			return fn;
	}

	/* Log(1 + x) - x */
	TARGET double LogOnePlusXMinusX(double x)
	{
		if (x <= -1.0)
			return NAN;
		else if (fabs(x) < 7.4009597974140505e-04)//ROOT5_DBL_EPSILON) 
		{
			double c1 = -0.5, c2 = 1.0 / 3.0, c3 = -1.0 / 4.0;
			double c4 = 1.0 / 5.0, c5 = -1.0 / 6.0, c6 = 1.0 / 7.0;
			double c7 = -1.0 / 8.0, c8 = 1.0 / 9.0, c9 = -1.0 / 10.0;
			double t = c5 + x * (c6 + x * (c7 + x * (c8 + x * c9)));
			return x * x * (c1 + x * (c2 + x * (c3 + x * (c4 + x * t))));
		}
		else if (fabs(x) < 0.5)
			return x * x * ChebEval(&lopxmx_cs, 0.5 * (8.0 * x + 1.0) / (x + 2.0));
		else
			return log(1.0 + x) - x;
	}

	/* Log(1 + x) - x */
	TARGET double LogOnePlusX(double x)
	{
		if (x <= -1.0)
			return NAN;
		else if (fabs(x) < 2.4607833005759251e-03)//ROOT6_DBL_EPSILON)
		{
			double c1 = -0.5, c2 = 1.0 / 3.0, c3 = -1.0 / 4.0;
			double c4 = 1.0 / 5.0, c5 = -1.0 / 6.0, c6 = 1.0 / 7.0;
			double c7 = -1.0 / 8.0, c8 = 1.0 / 9.0, c9 = -1.0 / 10.0;
			double t = c5 + x * (c6 + x * (c7 + x * (c8 + x * c9)));
			return x * (1.0 + x * (c1 + x * (c2 + x * (c3 + x * (c4 + x * t)))));
		}
		else
			return log(1 + x);
	}

	/* Complementary error function */
	TARGET double ErfC(double x)
	{
		double ax = fabs(x), e_val;

		if (ax <= 1.0)
			e_val = ChebEval(&erfc_xlt1_cs, 2.0 * ax - 1.0);
		else if (ax <= 5.0)
			e_val = exp(-x * x) * ChebEval(&erfc_x15_cs, 0.5 * (ax - 3.0));
		else if (ax < 10.0)
			e_val = exp(-x * x) / ax * ChebEval(&erfc_x510_cs, (2.0 * ax - 15.0) / 5.0);
		else
		{
			static double P[] = {
				2.97886562639399288862, 7.409740605964741794425, 6.1602098531096305440906,
				5.019049726784267463450058, 1.275366644729965952479585264, 0.5641895835477550741253201704
			};
			static double Q[] = {
				3.3690752069827527677, 9.608965327192787870698, 17.08144074746600431571095,
				12.0489519278551290360340491, 9.396034016235054150430579648, 2.260528520767326969591866945, 1.0
			};
			double num = 0.0, den = 0.0;
			int i;

			num = P[5];
			for (i = 4; i >= 0; --i)
				num = ax * num + P[i];
			den = Q[6];
			for (i = 5; i >= 0; --i)
				den = ax * den + Q[i];

			return num / den * exp(-ax * ax);
		}

		if (x < 0.0)
			return 2.0 - e_val;
		else
			return e_val;
	}

	/* Natual logarithm of complementary error function*/
	TARGET double LogErfC(double x)
	{
		if (x * x < 10.0 * 2.4607833005759251e-03)// ROOT6_DBL_EPSILON) 
		{
			double y = x / M_SQRTPI;
			double c3 = (4.0 - M_PI) / 3.0, c4 = 2.0 * (1.0 - M_PI / 3.0), c5 = -0.001829764677455021;
			double c6 = 0.02629651521057465, c7 = -0.01621575378835404, c8 = 0.00125993961762116;
			double c9 = 0.00556964649138, c10 = -0.0045563339802, c11 = 0.0009461589032;
			double c12 = 0.0013200243174, c13 = -0.00142906, c14 = 0.00048204;
			double series = c8 + y * (c9 + y * (c10 + y * (c11 + y * (c12 + y * (c13 + c14 * y)))));
			series = y * (1.0 + y * (1.0 + y * (c3 + y * (c4 + y * (c5 + y * (c6 + y * (c7 + y * series)))))));
			return -2.0 * series;
		}
		else if (x > 8.0)
		{
			static double P[] = {
				2.97886562639399288862, 7.409740605964741794425, 6.1602098531096305440906,
				5.019049726784267463450058, 1.275366644729965952479585264, 0.5641895835477550741253201704
			};
			static double Q[] = {
				3.3690752069827527677, 9.608965327192787870698, 17.08144074746600431571095,
				12.0489519278551290360340491, 9.396034016235054150430579648, 2.260528520767326969591866945, 1.0
			};
			double num = 0.0, den = 0.0;
			int i;

			num = P[5];
			for (i = 4; i >= 0; --i)
				num = x * num + P[i];

			den = Q[6];
			for (i = 5; i >= 0; --i)
				den = x * den + Q[i];

			return log(num / den) - x * x;
		}
		else
			return log(ErfC(x));
	}

	/* Incomplete gamma function detail: The dominant part */
	TARGET double GammaIncD(double a, double x)
	{
		if (a < 10.0)
			return exp(a * log(x) - x - LogGamma(a + 1.0));
		else
		{
			double ln_term, gstar;
			if (x < 0.5 * a)
				ln_term = log(x / a) - x / a + 1.0;
			else
				ln_term = LogOnePlusX((x - a) / a);

			if (x <= 0.0)
				gstar = NAN;
			else if (x < 0.5)
				gstar = exp(LogGamma(x) - (x - 0.5) * log(x) + x - 0.5 * (M_LN2 + M_LNPI));
			else if (x < 2.0)
				gstar = ChebEval(&gstar_a_cs, 4.0 / 3.0 * (x - 0.5) - 1.0);
			else if (x < 10.0)
				gstar = ChebEval(&gstar_b_cs, 0.25 * (x - 2.0) - 1.0) / (x * x) + 1.0 + 1.0 / (12.0 * x);
			else if (x < 1.0 / 1.2207031250000000e-04)//ROOT4_DBL_EPSILON)
			{
				const double y = 1.0 / (x * x);
				const double c0 = 1.0 / 12.0;
				const double c1 = -1.0 / 360.0;
				const double c2 = 1.0 / 1260.0;
				const double c3 = -1.0 / 1680.0;
				const double c4 = 1.0 / 1188.0;
				const double c5 = -691.0 / 360360.0;
				const double c6 = 1.0 / 156.0;
				const double c7 = -3617.0 / 122400.0;
				const double ser = c0 + y * (c1 + y * (c2 + y * (c3 + y * (c4 + y * (c5 + y * (c6 + y * c7))))));
				gstar = exp(ser / x);
			}
			else if (x < 1.0 / DBL_EPSILON)
			{
				double xi = 1.0 / x;
				gstar = 1.0 + xi / 12.0 * (1.0 + xi / 24.0 * (1.0 - xi * (139.0 / 180.0 + 571.0 / 8640.0 * xi)));
			}
			else
				gstar = 1.0;

			return exp(a * ln_term) / sqrt(2.0 * M_PI * a) / gstar;
		}
	}

	/* Incomplete gamma function detail: The dominant part of incomplete gamma function */
	TARGET double LogGammaIncD(double a, double x)
	{
		if (a < 10.0)
			return a * log(x) - x - LogGamma(a + 1.0);
		else
		{
			double ln_term, gstar;
			if (x < 0.5 * a)
				ln_term = log(x / a) - x / a + 1.0;
			else
				ln_term = LogOnePlusX((x - a) / a);

			if (x <= 0.0)
				gstar = NAN;
			else if (x < 0.5)
				gstar = exp(LogGamma(x) - (x - 0.5) * log(x) + x - 0.5 * (M_LN2 + M_LNPI));
			else if (x < 2.0)
				gstar = ChebEval(&gstar_a_cs, 4.0 / 3.0 * (x - 0.5) - 1.0);
			else if (x < 10.0)
				gstar = ChebEval(&gstar_b_cs, 0.25 * (x - 2.0) - 1.0) / (x * x) + 1.0 + 1.0 / (12.0 * x);
			else if (x < 1.0 / 1.2207031250000000e-04)//ROOT4_DBL_EPSILON)
			{
				const double y = 1.0 / (x * x);
				const double c0 = 1.0 / 12.0;
				const double c1 = -1.0 / 360.0;
				const double c2 = 1.0 / 1260.0;
				const double c3 = -1.0 / 1680.0;
				const double c4 = 1.0 / 1188.0;
				const double c5 = -691.0 / 360360.0;
				const double c6 = 1.0 / 156.0;
				const double c7 = -3617.0 / 122400.0;
				const double ser = c0 + y * (c1 + y * (c2 + y * (c3 + y * (c4 + y * (c5 + y * (c6 + y * c7))))));
				gstar = exp(ser / x);
			}
			else if (x < 1.0 / DBL_EPSILON)
			{
				double xi = 1.0 / x;
				gstar = 1.0 + xi / 12.0 * (1.0 + xi / 24.0 * (1.0 - xi * (139.0 / 180.0 + 571.0 / 8640.0 * xi)));
			}
			else
				gstar = 1.0;

			return a * ln_term - 0.5 * log(2.0 * M_PI * a) - log(gstar);
		}
	}

	/* Incomplete gamma function detail: Q large x asymptotic */
	TARGET double GammaIncQLargeX(double a, double x)
	{
		int nmax = 5000, n;
		double D = GammaIncD(a, x);

		double sum = 1.0, term = 1.0, last = 1.0;
		for (n = 1; n < nmax; n++)
		{
			term *= (a - n) / x;
			if (fabs(term / last) > 1.0) break;
			if (fabs(term / sum) < DBL_EPSILON) break;
			sum += term;
			last = term;
		}

		if (n == nmax)
			return NAN;

		return D * (a / x) * sum;
	}

	/* Incomplete gamma function detail: Q large x asymptotic */
	TARGET double LogGammaIncQLargeX(double a, double x)
	{
		int nmax = 5000, n;
		double logD = LogGammaIncD(a, x);

		double sum = 1.0, term = 1.0, last = 1.0;
		for (n = 1; n < nmax; n++)
		{
			term *= (a - n) / x;
			if (fabs(term / last) > 1.0) break;
			if (fabs(term / sum) < DBL_EPSILON) break;
			sum += term;
			last = term;
		}

		if (n == nmax)
			return NAN;

		return logD + log(a) - log(x) + log(sum);
	}

	/* Incomplete gamma function detail: Uniform asymptotic for x near a, a and x large. */
	TARGET double GammaIncQAsymp(double a, double x)
	{
		double rta = sqrt(a), eps = (x - a) / a;
		double ln_term = LogOnePlusXMinusX(eps);
		double eta = (eps >= 0.0 ? 1 : -1) * sqrt(-2.0 * ln_term);

		double R, c0, c1;
		double erfc = 2 * (1 - NormalDistCDF(eta * rta));

		if (fabs(eps) < 7.4009597974140505e-04)//ROOT5_DBL_EPSILON) 
		{
			c0 = -1.0 / 3.0 + eps * (1.0 / 12.0 - eps * (23.0 / 540.0 - eps * (353.0 / 12960.0 - eps * 589.0 / 30240.0)));
			c1 = -1.0 / 540.0 - eps / 288.0;
		}
		else
		{
			const double rt_term = sqrt(-2.0 * ln_term / (eps * eps));
			const double lam = x / a;
			c0 = (1.0 - 1.0 / rt_term) / eps;
			c1 = -(eta * eta * eta * (lam * lam + 10.0 * lam + 1.0) - 12.0 * eps * eps * eps) / (12.0 * eta * eta * eta * eps * eps * eps);
		}

		R = exp(-0.5 * a * eta * eta) / (M_SQRT2 * M_SQRTPI * rta) * (c0 + c1 / a);

		return 0.5 * erfc + R;
	}

	/* Incomplete gamma function detail: Uniform asymptotic for x near a, a and x large. */
	TARGET double LogGammaIncQAsymp(double a, double x)
	{
		return log(GammaIncQAsymp(a, x));
	}

	/* Incomplete gamma function detail: regularized lower series expansion */
	TARGET double GammaIncPSeries(double a, double x)
	{
		int nmax = 10000;

		double D = GammaIncD(a, x);

		if (x > 0.995 * a && a > 1e5)
			return D * ExprelNCF(a, x);
		if (x > (a + nmax))
			return NAN;

		double sum = 1.0, term = 1.0, remainder;
		int n, nlow = (x > a) ? (x - a) : 0;

		for (n = 1; n < nlow; n++)
		{
			term *= x / (a + n);
			sum += term;
		}

		for (; n < nmax; n++)
		{
			term *= x / (a + n);
			sum += term;
			if (fabs(term / sum) < DBL_EPSILON)
				break;
		}

		remainder = (x / (a + n)) * term / (1.0 - x / (a + n + 1.0));

		if (n == nmax && fabs(remainder / sum) > 1.4901161193847656e-08)//SQRT_DBL_EPSILON
			return NAN;

		return D * sum;
	}

	/* Incomplete gamma function detail: regularized upper series expansion */
	TARGET double GammaIncQSeries(double a, double x)
	{
		double term1, sum = 1.0;

		const double pg21 = -2.404113806319188570799476;  /* PolyGamma[2,1] */
		const double lnx = log(x);
		const double el = M_EULER + lnx;
		const double c1 = -el;
		const double c2 = M_PI * M_PI / 12.0 - 0.5 * el * el;
		const double c3 = el * (M_PI * M_PI / 12.0 - el * el / 6.0) + pg21 / 6.0;
		const double c4 = -0.04166666666666666667 * (-1.758243446661483480 + lnx) * (-0.764428657272716373 + lnx) * (0.723980571623507657 + lnx) * (4.107554191916823640 + lnx);
		const double c5 = -0.0083333333333333333 * (-2.06563396085715900 + lnx) * (-1.28459889470864700 + lnx) * (-0.27583535756454143 + lnx) * (1.33677371336239618 + lnx) * (5.17537282427561550 + lnx);
		const double c6 = -0.0013888888888888889 * (-2.30814336454783200 + lnx) * (-1.65846557706987300 + lnx) * (-0.88768082560020400 + lnx) * (0.17043847751371778 + lnx) * (1.92135970115863890 + lnx) * (6.22578557795474900 + lnx);
		const double c7 = -0.00019841269841269841 * (-2.5078657901291800 + lnx) * (-1.9478900888958200 + lnx) * (-1.3194837322612730 + lnx) * (-0.5281322700249279 + lnx) * (0.5913834939078759 + lnx) * (2.4876819633378140 + lnx) * (7.2648160783762400 + lnx);
		const double c8 = -0.00002480158730158730 * (-2.677341544966400 + lnx) * (-2.182810448271700 + lnx) * (-1.649350342277400 + lnx) * (-1.014099048290790 + lnx) * (-0.191366955370652 + lnx) * (0.995403817918724 + lnx) * (3.041323283529310 + lnx) * (8.295966556941250 + lnx);
		const double c9 = -2.75573192239859e-6 * (-2.8243487670469080 + lnx) * (-2.3798494322701120 + lnx) * (-1.9143674728689960 + lnx) * (-1.3814529102920370 + lnx) * (-0.7294312810261694 + lnx) * (0.1299079285269565 + lnx) * (1.3873333251885240 + lnx) * (3.5857258865210760 + lnx) * (9.3214237073814600 + lnx);
		const double c10 = -2.75573192239859e-7 * (-2.9540329644556910 + lnx) * (-2.5491366926991850 + lnx) * (-2.1348279229279880 + lnx) * (-1.6741881076349450 + lnx) * (-1.1325949616098420 + lnx) * (-0.4590034650618494 + lnx) * (0.4399352987435699 + lnx) * (1.7702236517651670 + lnx) * (4.1231539047474080 + lnx) * (10.342627908148680 + lnx);

		term1 = a * (c1 + a * (c2 + a * (c3 + a * (c4 + a * (c5 + a * (c6 + a * (c7 + a * (c8 + a * (c9 + a * c10)))))))));

		int nmax = 5000, n;
		double t = 1.0;

		for (n = 1; n < nmax; n++)
		{
			t *= -x / (n + 1.0);
			sum += (a + 1.0) / (a + n + 1.0) * t;
			if (fabs(t / sum) < DBL_EPSILON)
				break;
		}

		if (n == nmax)
			return NAN;

		return term1 + (1.0 - term1) * a / (a + 1.0) * x * sum;
	}

	/* Incomplete gamma function detail: regularized upper series expansion */
	TARGET double LogGammaIncQSeries(double a, double x)
	{
		return log(GammaIncQSeries(a, x));
	}

	/* Incomplete gamma function detail: Continued fraction which occurs in evaluation of Q(a,x) or Gamma(a,x) */
	TARGET double GammaIncFCF(double a, double x)
	{
		int    nmax = 5000, n;
		double s = DBL_EPSILON * DBL_EPSILON * DBL_EPSILON;
		double hn = 1.0, Cn = 1.0 / s, Dn = 1.0;

		for (n = 2; n < nmax; n++)
		{
			double an;
			double delta;

			if (n & 1)
				an = 0.5 * (n - 1) / x;
			else
				an = (0.5 * n - a) / x;

			Dn = 1.0 + an * Dn;
			if (fabs(Dn) < s)
				Dn = s;
			Cn = 1.0 + an / Cn;
			if (fabs(Cn) < s)
				Cn = s;
			Dn = 1.0 / Dn;
			delta = Cn * Dn;
			hn *= delta;
			if (fabs(delta - 1.0) < DBL_EPSILON) break;
		}

		if (n == nmax)
			return NAN;
		else
			return hn;
	}

	/* Incomplete gamma function detail: Continued fraction which occurs in evaluation of Q(a,x) or Gamma(a,x) */
	TARGET double LogGammaIncFCF(double a, double x)
	{
		int    nmax = 5000, n;
		double s = DBL_EPSILON * DBL_EPSILON * DBL_EPSILON;
		double loghn = 1.0, Cn = 1.0 / s, Dn = 1.0;

		int64 slog = 0; double prod = 1;
		OpenLog(slog, prod);
		for (n = 2; n < nmax; n++)
		{
			double an;
			double delta;

			if (n & 1)
				an = 0.5 * (n - 1) / x;
			else
				an = (0.5 * n - a) / x;

			Dn = 1.0 + an * Dn;
			if (fabs(Dn) < s)
				Dn = s;
			Cn = 1.0 + an / Cn;
			if (fabs(Cn) < s)
				Cn = s;
			Dn = 1.0 / Dn;
			delta = Cn * Dn;
			//loghn += log(delta);
			ChargeLog(slog, prod, delta);

			if (fabs(delta - 1.0) < DBL_EPSILON) break;
		}
		CloseLog(slog, prod);
		loghn = prod;

		if (n == nmax)
			return NAN;
		else
			return loghn;
	}

	/* Regularized upper incomplete gamma function */
	TARGET double GammaIncQ(double a, double x)
	{
		if (a < 0.0 || x < 0.0)
			return NAN;
		else if (x == 0.0)
			return 1.0;
		else if (a == 0.0)
			return 0.0;
		else if (x <= 0.5 * a)
			return 1.0 - GammaIncPSeries(a, x);
		else if (a >= 1.0e+06 && (x - a) * (x - a) < a)
			return GammaIncQAsymp(a, x);
		else if (a < 0.2 && x < 5.0)
			return GammaIncQSeries(a, x);
		else if (a <= x)
		{
			if (x <= 1.0e+06)
				return GammaIncD(a, x) * (a / x) * GammaIncFCF(a, x);
			else
				return GammaIncQLargeX(a, x);
		}
		else
		{
			if (x > a - sqrt(a))
				return GammaIncD(a, x) * (a / x) * GammaIncFCF(a, x);
			else
				return 1.0 - GammaIncPSeries(a, x);
		}
	}

	/* Regularized upper incomplete gamma function */
	TARGET double LogGammaIncQ(double a, double x)
	{
		if (a < 0.0 || x < 0.0)
			return NAN;
		else if (x == 0.0)
			return 0;
		else if (a == 0.0)
			return NAN;
		else if (x <= 0.5 * a)
			return log(1.0 - GammaIncPSeries(a, x));
		else if (a >= 1.0e+06 && (x - a) * (x - a) < a)
			return LogGammaIncQAsymp(a, x);
		else if (a < 0.2 && x < 5.0)
			return LogGammaIncQSeries(a, x);
		else if (a <= x)
		{
			if (x <= 1.0e+06)
				return LogGammaIncD(a, x) + log(a) - log(x) + LogGammaIncFCF(a, x);
			else
				return LogGammaIncQLargeX(a, x);
		}
		else
		{
			if (x > a - sqrt(a))
				return LogGammaIncD(a, x) + log(a) - log(x) + LogGammaIncFCF(a, x);
			else
				return log(1.0 - GammaIncPSeries(a, x));
		}
	}
}

namespace TR1
{
	/* Sign of Gamma(x) */
	double lgamma_sign(double x)
	{
		if (x > double(0))
			return double(1);
		else
		{
			const double sin_fact = std::sin(M_PI * x);
			if (sin_fact > double(0))
				return (1);
			else if (sin_fact < double(0))
				return -double(1);
			else
				return double(0);
		}
	}

	/* Bernoulli number */
	double Bernoulli(unsigned int n)
	{
		static const double num[28] = {
		  double(1UL),                        -double(1UL) / double(2UL),
		  double(1UL) / double(6UL),             double(0UL),
		  -double(1UL) / double(30UL),           double(0UL),
		  double(1UL) / double(42UL),            double(0UL),
		  -double(1UL) / double(30UL),           double(0UL),
		  double(5UL) / double(66UL),            double(0UL),
		  -double(691UL) / double(2730UL),       double(0UL),
		  double(7UL) / double(6UL),             double(0UL),
		  -double(3617UL) / double(510UL),       double(0UL),
		  double(43867UL) / double(798UL),       double(0UL),
		  -double(174611) / double(330UL),       double(0UL),
		  double(854513UL) / double(138UL),      double(0UL),
		  -double(236364091UL) / double(2730UL), double(0UL),
		  double(8553103UL) / double(6UL),       double(0UL)
		};

		if (n == 0)
			return double(1);

		if (n == 1)
			return -double(1) / double(2);

		if (n % 2 == 1)
			return double(0);

		if (n < 28)
			return num[n];

		double fact = double(1);
		if ((n / 2) % 2 == 0)
			fact *= double(-1);
		for (unsigned int k = 1; k <= n; ++k)
			fact *= k / (double(2) * M_PI);
		fact *= double(2);

		double sum = double(0);
		for (unsigned int i = 1; i < 1000; ++i)
		{
			double term = pow(double(i), -double(n));
			if (term < numeric_limits<double>::epsilon())
				break;
			sum += term;
		}

		return fact * sum;
	}

	/* Hurwitz zeta function */
	double HurwitzZeta(double a, double s)
	{
		double zeta = double(0);
		const double eps = numeric_limits<double>::epsilon();
		const double max_bincoeff = numeric_limits<double>::max_exponent10 * log(double(10)) - double(1);

		const unsigned int maxit = 10000;
		for (unsigned int i = 0; i < maxit; ++i)
		{
			bool punt = false;
			double sgn = double(1);
			double term = double(0);
			for (unsigned int j = 0; j <= i; ++j)
			{
				double bincoeff = lgamma(double(1 + i)) - lgamma(double(1 + j)) - lgamma(double(1 + i - j));
				if (bincoeff > max_bincoeff)
				{
					punt = true;
					break;
				}
				bincoeff = exp(bincoeff);
				term += sgn * bincoeff * pow(double(a + j), -s);
				sgn *= double(-1);
			}
			if (punt)
				break;
			term /= double(i + 1);
			if (abs(term / zeta) < eps)
				break;
			zeta += term;
		}

		zeta /= s - double(1);

		return zeta;
	}

	/* Digamma function by series expansion */
	double Psi_Series(double x)
	{
		double sum = -M_EULER - double(1) / x;
		const unsigned int max_iter = 100000;
		for (unsigned int k = 1; k < max_iter; ++k)
		{
			const double term = x / (k * (k + x));
			sum += term;
			if (abs(term / sum) < numeric_limits<double>::epsilon())
				break;
		}
		return sum;
	}

	/* Digamma function for large argument */
	double Psi_Asymp(double x)
	{
		double sum = log(x) - double(0.5L) / x;
		const double xx = x * x;
		double xp = xx;
		const unsigned int max_iter = 100;
		for (unsigned int k = 1; k < max_iter; ++k)
		{
			const double term = Bernoulli(2 * k) / (2 * k * xp);
			sum -= term;
			if (abs(term / sum) < numeric_limits<double>::epsilon())
				break;
			xp *= xx;
		}
		return sum;
	}

	/* Digamma function */
	double Psi(double x)
	{
		const int n = static_cast<int>(x + 0.5L);
		const double eps = double(4) * numeric_limits<double>::epsilon();
		if (n <= 0 && abs(x - double(n)) < eps)
			return numeric_limits<double>::quiet_NaN();
		else if (x < double(0))
		{
			const double pi = M_PI;
			return Psi(double(1) - x) - pi * cos(pi * x) / sin(pi * x);
		}
		else if (x > double(100))
			return Psi_Asymp(x);
		else
			return Psi_Series(x);
	}

	/* Polygamma function */
	double Psi(unsigned int n, double x)
	{
		if (x <= double(0))
			return NAN;
		else if (n == 0)
			return Psi(x);
		else
		{
			const double hzeta = HurwitzZeta(double(n + 1), x);
			const double ln_nfact = lgamma(double(n + 1));
			double result = exp(ln_nfact) * hzeta;
			if (n % 2 == 1)
				result = -result;
			return result;
		}
	}

	/* Hpergeometric function 1F1(a,c;x) by series expansion. */
	double Hyperg1F1_Series(double a, double c, double x)
	{
		const double eps = numeric_limits<double>::epsilon();

		double term = double(1);
		double Fac = double(1);
		const unsigned int max_iter = 100000;
		unsigned int i;
		for (i = 0; i < max_iter; ++i)
		{
			term *= (a + double(i)) * x / ((c + double(i)) * double(1 + i));
			if (abs(term) < eps)
				break;
			Fac += term;
		}
		if (i == max_iter)
			return NAN;

		return Fac;
	}

	/* Hypergeometric function 1F1(a,b;c;x) by Luke */
	double Hyperg1F1_Luke(double a, double c, double xin)
	{
		const double big = pow(numeric_limits<double>::max(), double(0.16L));
		const int nmax = 20000;
		const double eps = numeric_limits<double>::epsilon();
		const double x = -xin;
		const double x3 = x * x * x;
		const double t0 = a / c;
		const double t1 = (a + double(1)) / (double(2) * c);
		const double t2 = (a + double(2)) / (double(2) * (c + double(1)));
		double F = double(1);
		double prec;

		double Bnm3 = double(1);
		double Bnm2 = double(1) + t1 * x;
		double Bnm1 = double(1) + t2 * x * (double(1) + t1 / double(3) * x);

		double Anm3 = double(1);
		double Anm2 = Bnm2 - t0 * x;
		double Anm1 = Bnm1 - t0 * (double(1) + t2 * x) * x + t0 * t1 * (c / (c + double(1))) * x * x;

		int n = 3;
		while (1)
		{
			double npam1 = double(n - 1) + a;
			double npcm1 = double(n - 1) + c;
			double npam2 = double(n - 2) + a;
			double npcm2 = double(n - 2) + c;
			double tnm1 = double(2 * n - 1);
			double tnm3 = double(2 * n - 3);
			double tnm5 = double(2 * n - 5);
			double F1 = (double(n - 2) - a) / (double(2) * tnm3 * npcm1);
			double F2 = (double(n) + a) * npam1 / (double(4) * tnm1 * tnm3 * npcm2 * npcm1);
			double F3 = -npam2 * npam1 * (double(n - 2) - a) / (double(8) * tnm3 * tnm3 * tnm5 * (double(n - 3) + c) * npcm2 * npcm1);
			double E = -npam1 * (double(n - 1) - c) / (double(2) * tnm3 * npcm2 * npcm1);

			double An = (double(1) + F1 * x) * Anm1 + (E + F2 * x) * x * Anm2 + F3 * x3 * Anm3;
			double Bn = (double(1) + F1 * x) * Bnm1 + (E + F2 * x) * x * Bnm2 + F3 * x3 * Bnm3;
			double r = An / Bn;

			prec = abs((F - r) / F);
			F = r;

			if (prec < eps || n > nmax)
				break;

			if (abs(An) > big || abs(Bn) > big)
			{
				An /= big;
				Bn /= big;
				Anm1 /= big;
				Bnm1 /= big;
				Anm2 /= big;
				Bnm2 /= big;
				Anm3 /= big;
				Bnm3 /= big;
			}
			else if (abs(An) < double(1) / big
				|| abs(Bn) < double(1) / big)
			{
				An *= big;
				Bn *= big;
				Anm1 *= big;
				Bnm1 *= big;
				Anm2 *= big;
				Bnm2 *= big;
				Anm3 *= big;
				Bnm3 *= big;
			}

			++n;
			Bnm3 = Bnm2;
			Bnm2 = Bnm1;
			Bnm1 = Bn;
			Anm3 = Anm2;
			Anm2 = Anm1;
			Anm1 = An;
		}

		if (n >= nmax)
			return NAN;

		return F;
	}

	/* Hypergeometric function 1F1(a;c;x) */
	double Hypergeometric1F1(double a, double c, double x)
	{
		const double c_nint = static_cast<int>(c + double(0.5L));
		if (isnan(a) || isnan(c) || isnan(x))
			return numeric_limits<double>::quiet_NaN();
		else if (c_nint == c && c_nint <= 0)
			return numeric_limits<double>::infinity();
		else if (a == double(0))
			return double(1);
		else if (c == a)
			return exp(x);
		else if (x < double(0))
			return Hyperg1F1_Luke(a, c, x);
		else
			return Hyperg1F1_Series(a, c, x);
	}

	/* Hypergeometric function 2F1 */
	double Hyperg2F1_Series(double a, double b, double c, double x)
	{
		const double eps = numeric_limits<double>::epsilon();

		double term = double(1);
		double Fabc = double(1);
		const unsigned int max_iter = 100000;
		unsigned int i;
		for (i = 0; i < max_iter; ++i)
		{
			term *= (a + double(i)) * (b + double(i)) * x / ((c + double(i)) * double(1 + i));
			if (abs(term) < eps)
				break;
			Fabc += term;
		}
		if (i == max_iter)
			return NAN;

		return Fabc;
	}

	/* Hypergeometric function 2F1(a,b;c;x) by Luke */
	double Hyperg2F1_Luke(double a, double b, double c, double xin)
	{
		const double big = pow(numeric_limits<double>::max(), double(0.16L));
		const int nmax = 20000;
		const double eps = numeric_limits<double>::epsilon();
		const double x = -xin;
		const double x3 = x * x * x;
		const double t0 = a * b / c;
		const double t1 = (a + double(1)) * (b + double(1)) / (double(2) * c);
		const double t2 = (a + double(2)) * (b + double(2)) / (double(2) * (c + double(1)));

		double F = double(1);

		double Bnm3 = double(1);
		double Bnm2 = double(1) + t1 * x;
		double Bnm1 = double(1) + t2 * x * (double(1) + t1 / double(3) * x);

		double Anm3 = double(1);
		double Anm2 = Bnm2 - t0 * x;
		double Anm1 = Bnm1 - t0 * (double(1) + t2 * x) * x + t0 * t1 * (c / (c + double(1))) * x * x;

		int n = 3;
		while (1)
		{
			const double npam1 = double(n - 1) + a;
			const double npbm1 = double(n - 1) + b;
			const double npcm1 = double(n - 1) + c;
			const double npam2 = double(n - 2) + a;
			const double npbm2 = double(n - 2) + b;
			const double npcm2 = double(n - 2) + c;
			const double tnm1 = double(2 * n - 1);
			const double tnm3 = double(2 * n - 3);
			const double tnm5 = double(2 * n - 5);
			const double n2 = n * n;
			const double F1 = (double(3) * n2 + (a + b - double(6)) * n + double(2) - a * b - double(2) * (a + b)) / (double(2) * tnm3 * npcm1);
			const double F2 = -(double(3) * n2 - (a + b + double(6)) * n + double(2) - a * b) * npam1 * npbm1 / (double(4) * tnm1 * tnm3 * npcm2 * npcm1);
			const double F3 = (npam2 * npam1 * npbm2 * npbm1 * (double(n - 2) - a) * (double(n - 2) - b)) / (double(8) * tnm3 * tnm3 * tnm5 * (double(n - 3) + c) * npcm2 * npcm1);
			const double E = -npam1 * npbm1 * (double(n - 1) - c) / (double(2) * tnm3 * npcm2 * npcm1);

			double An = (double(1) + F1 * x) * Anm1 + (E + F2 * x) * x * Anm2 + F3 * x3 * Anm3;
			double Bn = (double(1) + F1 * x) * Bnm1 + (E + F2 * x) * x * Bnm2 + F3 * x3 * Bnm3;
			const double r = An / Bn;

			const double prec = abs((F - r) / F);
			F = r;

			if (prec < eps || n > nmax)
				break;

			if (abs(An) > big || abs(Bn) > big)
			{
				An /= big;
				Bn /= big;
				Anm1 /= big;
				Bnm1 /= big;
				Anm2 /= big;
				Bnm2 /= big;
				Anm3 /= big;
				Bnm3 /= big;
			}
			else if (abs(An) < double(1) / big || abs(Bn) < double(1) / big)
			{
				An *= big;
				Bn *= big;
				Anm1 *= big;
				Bnm1 *= big;
				Anm2 *= big;
				Bnm2 *= big;
				Anm3 *= big;
				Bnm3 *= big;
			}

			++n;
			Bnm3 = Bnm2;
			Bnm2 = Bnm1;
			Bnm1 = Bn;
			Anm3 = Anm2;
			Anm2 = Anm1;
			Anm1 = An;
		}

		if (n >= nmax)
			return NAN;

		return F;
	}

	/* Hypergeometric function 2F1(a,b;c;x) by the Reflection formulae */
	double Hyperg2F1_Reflect(double a, double b, double c, double x)
	{
		const double d = c - a - b;
		const int intd = floor(d + double(0.5L));
		const double eps = numeric_limits<double>::epsilon();
		const double toler = double(1000) * eps;
		const double log_max = log(numeric_limits<double>::max());
		const bool d_integer = (abs(d - intd) < toler);

		if (d_integer)
		{
			const double ln_omx = log(double(1) - x);
			const double ad0 = abs(d);
			double F1, F2;

			double d1, d2;
			if (d >= double(0))
			{
				d1 = d;
				d2 = double(0);
			}
			else
			{
				d1 = double(0);
				d2 = d;
			}

			const double lng_c = lgamma(c);

			if (ad0 < eps)
				F1 = double(0);
			else
			{

				bool ok_d1 = true;
				double lng_ad, lng_ad1, lng_bd1;
				//try
				{
					lng_ad = lgamma(ad0);
					lng_ad1 = lgamma(a + d1);
					lng_bd1 = lgamma(b + d1);
					ok_d1 = IsNormal(lng_ad + lng_ad1 + lng_bd1);
				}
				//catch (...)
				{
					//    ok_d1 = false;
				}

				if (ok_d1)
				{
					double sum1 = double(1);
					double term = double(1);
					double ln_pre1 = lng_ad + lng_c + d2 * ln_omx - lng_ad1 - lng_bd1;

					for (int i = 1; i < ad0; ++i)
					{
						const int j = i - 1;
						term *= (a + d2 + j) * (b + d2 + j) / (double(1) + d2 + j) / i * (double(1) - x);
						sum1 += term;
					}

					if (ln_pre1 > log_max)
						return NAN;
					else
						F1 = exp(ln_pre1) * sum1;
				}
				else
					F1 = double(0);
			}

			bool ok_d2 = true;
			double lng_ad2, lng_bd2;
			//try
			{
				lng_ad2 = lgamma(a + d2);
				lng_bd2 = lgamma(b + d2);
				ok_d2 = IsNormal(lng_ad2 + lng_bd2);
			}
			//catch(...)
			{
				//ok_d2 = false;
			}

			if (ok_d2)
			{
				const int maxiter = 2000;
				const double psi_1 = -M_EULER;
				const double psi_1pd = Psi(double(1) + ad0);
				const double psi_apd1 = Psi(a + d1);
				const double psi_bpd1 = Psi(b + d1);

				double psi_term = psi_1 + psi_1pd - psi_apd1 - psi_bpd1 - ln_omx;
				double fact = double(1);
				double sum2 = psi_term;
				double ln_pre2 = lng_c + d1 * ln_omx - lng_ad2 - lng_bd2;

				int j;
				for (j = 1; j < maxiter; ++j)
				{
					const double term1 = double(1) / double(j) + double(1) / (ad0 + j);
					const double term2 = double(1) / (a + d1 + double(j - 1)) + double(1) / (b + d1 + double(j - 1));
					psi_term += term1 - term2;
					fact *= (a + d1 + double(j - 1)) * (b + d1 + double(j - 1)) / ((ad0 + j) * j) * (double(1) - x);
					const double delta = fact * psi_term;
					sum2 += delta;
					if (abs(delta) < eps * abs(sum2))
						break;
				}
				if (j == maxiter)
					return NAN;

				if (sum2 == double(0))
					F2 = double(0);
				else
					F2 = exp(ln_pre2) * sum2;
			}
			else
				F2 = double(0);

			const double sgn_2 = (intd % 2 == 1 ? -double(1) : double(1));
			const double F = F1 + sgn_2 * F2;

			return F;
		}
		else
		{
			bool ok1 = true;
			double sgn_g1ca = double(0), ln_g1ca = double(0);
			double sgn_g1cb = double(0), ln_g1cb = double(0);
			//try
			{
				sgn_g1ca = lgamma_sign(c - a);
				ln_g1ca = lgamma(c - a);
				sgn_g1cb = lgamma_sign(c - b);
				ln_g1cb = lgamma(c - b);
				ok1 = IsNormal(sgn_g1ca + ln_g1ca + sgn_g1cb + ln_g1cb);
			}
			//catch(...)
			{
				//ok1 = false;
			}

			bool ok2 = true;
			double sgn_g2a = double(0), ln_g2a = double(0);
			double sgn_g2b = double(0), ln_g2b = double(0);
			//try
			{
				sgn_g2a = lgamma_sign(a);
				ln_g2a = lgamma(a);
				sgn_g2b = lgamma_sign(b);
				ln_g2b = lgamma(b);
				ok2 = IsNormal(sgn_g2a + ln_g2a + sgn_g2b + ln_g2b);
			}
			//catch(...)
			{
				//ok2 = false;
			}

			const double sgn_gc = lgamma_sign(c);
			const double ln_gc = lgamma(c);
			const double sgn_gd = lgamma_sign(d);
			const double ln_gd = lgamma(d);
			const double sgn_gmd = lgamma_sign(-d);
			const double ln_gmd = lgamma(-d);

			const double sgn1 = sgn_gc * sgn_gd * sgn_g1ca * sgn_g1cb;
			const double sgn2 = sgn_gc * sgn_gmd * sgn_g2a * sgn_g2b;

			double pre1, pre2;
			if (ok1 && ok2)
			{
				double ln_pre1 = ln_gc + ln_gd - ln_g1ca - ln_g1cb;
				double ln_pre2 = ln_gc + ln_gmd - ln_g2a - ln_g2b + d * log(double(1) - x);
				if (ln_pre1 < log_max && ln_pre2 < log_max)
				{
					pre1 = exp(ln_pre1);
					pre2 = exp(ln_pre2);
					pre1 *= sgn1;
					pre2 *= sgn2;
				}
				else
					return NAN;
			}
			else if (ok1 && !ok2)
			{
				double ln_pre1 = ln_gc + ln_gd - ln_g1ca - ln_g1cb;
				if (ln_pre1 < log_max)
				{
					pre1 = exp(ln_pre1);
					pre1 *= sgn1;
					pre2 = double(0);
				}
				else
					return NAN;
			}
			else if (!ok1 && ok2)
			{
				double ln_pre2 = ln_gc + ln_gmd - ln_g2a - ln_g2b + d * log(double(1) - x);
				if (ln_pre2 < log_max)
				{
					pre1 = double(0);
					pre2 = exp(ln_pre2);
					pre2 *= sgn2;
				}
				else
					return NAN;
			}
			else
			{
				pre1 = double(0);
				pre2 = double(0);
				return NAN;
			}

			const double F1 = Hyperg2F1_Series(a, b, double(1) - d, double(1) - x);
			const double F2 = Hyperg2F1_Series(c - a, c - b, double(1) + d, double(1) - x);
			const double F = pre1 * F1 + pre2 * F2;

			return F;
		}
	}

	/* Hypergeometric function 2F1(a,b;c;x) */
	double Hypergeometric2F1(double a, double b, double c, double x)
	{
		const double a_nint = static_cast<int>(a + double(0.5L));
		const double b_nint = static_cast<int>(b + double(0.5L));
		const double c_nint = static_cast<int>(c + double(0.5L));
		const double toler = double(1000) * numeric_limits<double>::epsilon();
		if (abs(x) >= double(1))
			return NAN;
		else if (isnan(a) || isnan(b) || isnan(c) || isnan(x))
			return numeric_limits<double>::quiet_NaN();
		else if (c_nint == c && c_nint <= double(0))
			return numeric_limits<double>::infinity();
		else if (abs(c - b) < toler || abs(c - a) < toler)
			return pow(double(1) - x, c - a - b);
		else if (a >= double(0) && b >= double(0) && c >= double(0) && x >= double(0) && x < double(0.995L))
			return Hyperg2F1_Series(a, b, c, x);
		else if (abs(a) < double(10) && abs(b) < double(10))
		{
			if (a < double(0) && abs(a - a_nint) < toler)
				return Hyperg2F1_Series(a_nint, b, c, x);
			else if (b < double(0) && abs(b - b_nint) < toler)
				return Hyperg2F1_Series(a, b_nint, c, x);
			else if (x < -double(0.25L))
				return Hyperg2F1_Luke(a, b, c, x);
			else if (x < double(0.5L))
				return Hyperg2F1_Series(a, b, c, x);
			else
				if (abs(c) > double(10))
					return Hyperg2F1_Series(a, b, c, x);
				else
					return Hyperg2F1_Reflect(a, b, c, x);
		}
		else
			return Hyperg2F1_Luke(a, b, c, x);
	}
};

/* Gamma function */
TARGET double Gamma(double x)
{
	return tgamma(x);
}

/* Natural logarithm of Gamma function */
TARGET double LogGamma(double x)
{
	return lgamma(x);
}

/* Regularized incomplete Gamma function */
TARGET double GammaIncRegularized(double a, double x)
{
	double tv = gcem::incomplete_gamma(a, x);
	if (tv < 0.9)
		return 1 - tv;
	else
		return GSL::GammaIncQ(a, x);
}

/* Incomplete Gamma function */
TARGET double GammaInc(double a, double x)
{
	double tv = gcem::incomplete_gamma(a, x);
	double tg = tgamma(a);
	if (tv < 0.9)
		return tg - tg * tv;
	else
		return tg * GSL::GammaIncQ(a, x);
}

/* Natural logarithm of regularized incomplete Gamma function */
TARGET double LogGammaIncRegularized(double a, double x)
{
	return GSL::LogGammaIncQ(a, x);
}

/* Natural logarithm of incomplete Gamma function */
TARGET double LogGammaInc(double a, double x)
{
	return GSL::LogGammaIncQ(a, x) + lgamma(a);
}

/* Beta function */
TARGET double Beta(double a, double b)
{
	return tgamma(a) * tgamma(b) / tgamma(a + b);
}

/* Natural logarithm of Beta function */
TARGET double LogBeta(double a, double b)
{
	return lgamma(a) + lgamma(b) - lgamma(a + b);
}

/* Regularized incomplete Beta function */
TARGET double BetaIncRegularized(double a, double b, double z)
{
	return gcem::incomplete_beta(a, b, z);
}

/* Incomplete Beta function */
TARGET double BetaInc(double a, double b, double z)
{
	return gcem::incomplete_beta(a, b, z) * tgamma(a) * tgamma(b) / tgamma(a + b);
}

/* Natural logarithm of regularized incomplete Beta function */
TARGET double LogBetaIncRegularized(double a, double b, double z)
{
	// Hypergeometric2F1[a, 1 - b, a + 1, z]*z^a*/a
	return log(TR1::Hypergeometric2F1(a, 1 - b, a + 1, z)) + a * log(z) - log(a) - lgamma(a) - lgamma(b) + lgamma(a + b);
}

/* Natural logarithm of incomplete Beta function */
TARGET double LogBetaInc(double a, double b, double z)
{
	// Hypergeometric2F1[a, 1 - b, a + 1, z]*z^a*/a
	return log(TR1::Hypergeometric2F1(a, 1 - b, a + 1, z)) + a * log(z) - log(a);
}

/* Two tailled probabiliy of normal distribution */
TARGET double MinusLogPNormal(double x)
{
	return GSL::LogErfC(abs(x) / M_SQRT2);
	/*
	x = -abs(x);

	if (x < -6.9217225)
	{
		double i = 1.0 / x, i2 = i * i, i4 = i2 * i2, i6 = i2 * i4, i8 = i4 * i4;
		double i10 = i4 * i6, i12 = i6 * i6, i14 = i6 * i8, i16 = i8 * i8;
		double i18 = i8 * i10, i20 = i10 * i10;
		return 0.22579135264472743236
			- 6.1389139250000000000e8 * i20
			+ 3.2002080111111111111e7 * i18
			- 1.8595041250000000000e6 * i16
			+ 122028.14285714285714 * i14
			- 9200.8333333333333333 * i12
			+ 816.20000000000000000 * i10
			- 88.250000000000000000 * i8
			+ 12.333333333333333333 * i6
			- 2.5000000000000000000 * i4
			+ i2 + 0.5 * x * x - log(-i);
	}
	else if (x < -3)
	{
		double a = x + 4, a2 = a * a, a3 = a2 * a, a4 = a2 * a2;
		double a5 = a2 * a3, a6 = a3 * a3, a7 = a3 * a4, a8 = a4 * a4;
		double a9 = a4 * a5, a10 = a5 * a5, a11 = a5 * a6, a12 = a6 * a6;
		double a13 = a6 * a7, a14 = a7 * a7, a15 = a7 * a8, a16 = a8 * a8;
		double a17 = a8 * a9, a18 = a9 * a9, a19 = a9 * a10, a20 = a10 * a10;
		return 9.6669543059673455184
			- 4.2256071444894710728 * a
			+ 0.47666358080128868442 * a2
			- 0.0029760565512764043280 * a3
			- 0.00039610735131649454517 * a4
			- 0.000052106440317749825969 * a5
			- 6.5710966093077351106e-6 * a6
			- 7.7200912155870303342e-7 * a7
			- 8.0876155893300327799e-8 * a8
			- 6.8060680666457862119e-9 * a9
			- 2.7383593411643506794e-10 * a10
			+ 5.3399855404108932981e-11 * a11
			+ 1.8310629242590046417e-11 * a12
			+ 3.5947319257399716912e-12 * a13
			+ 5.5787836709990877572e-13 * a14
			+ 7.1757387846409384591e-14 * a15
			+ 7.2913006900839955228e-15 * a16
			+ 4.3053950051899499599e-16 * a17
			- 3.7605122218165314753e-17 * a18
			- 1.8710785773783494169e-17 * a19
			- 4.1111804893252978308e-18 * a20;
	}
	else if (x < -0.873)
	{
		double a = x + 2, a2 = a * a, a3 = a2 * a, a4 = a2 * a2;
		double a5 = a2 * a3, a6 = a3 * a3, a7 = a3 * a4, a8 = a4 * a4;
		double a9 = a4 * a5, a10 = a5 * 5, a11 = a5 * a6, a12 = a6 * a6;
		double a13 = a6 * a7, a14 = a7 * a7, a15 = a7 * a8, a16 = a8 * a8;
		double a17 = a8 * a9, a18 = a9 * a9, a19 = a9 * a10, a20 = a10 * a10;
		return 3.0900371531220866394
			- 2.3732155328228408673 * a
			+ 0.44286044979295937168 * a2
			- 0.0098926435485943021389 * a3
			- 0.0016425830777477838861 * a4
			- 0.00024249157212790694704 * a5
			- 0.000028597249249507913641 * a6
			- 1.7645404227607055909e-6 * a7
			+ 2.9807834350019200181e-7 * a8
			+ 1.3943185588600379280e-7 * a9
			+ 3.1678520808621568654e-8 * a10
			+ 4.8626172534760406799e-9 * a11
			+ 3.8674654641876361599e-10 * a12
			- 5.3330021778376005693e-11 * a13
			- 3.0560940174483102109e-11 * a14
			- 7.5992006467586071649e-12 * a15
			- 1.2523658586839622891e-12 * a16
			- 1.0973167322183581916e-13 * a17
			+ 1.2765176605609572190e-14 * a18
			+ 8.2691331868880010361e-15 * a19
			+ 2.1496979842761700291e-15 * a20;
	}
	else
	{
		double a = x, a2 = a * a, a3 = a2 * a, a4 = a2 * a2;
		double a5 = a2 * a3, a6 = a3 * a3, a7 = a3 * a4, a8 = a4 * a4;
		double a9 = a4 * a5, a10 = a5 * 5, a11 = a5 * a6, a12 = a6 * a6;
		double a13 = a6 * a7, a14 = a7 * a7, a15 = a7 * a8, a16 = a8 * a8;
		double a17 = a8 * a9, a18 = a9 * a9, a19 = a9 * a10, a20 = a10 * a10;
		return -0.79788456080286535588 * a
			+ 0.31830988618379067154 * a2
			- 0.036335602357498360115 * a3
			- 0.0047821117522591190687 * a4
			+ 0.000036980737188481841574 * a5
			+ 0.00021202574144674573625 * a6
			+ 0.000052160001839763418049 * a7
			+ 1.6168147195571145235e-6 * a8
			- 2.8513348301409302801e-6 * a9
			- 9.3056390789070481885e-7 * a10
			- 7.7091198959303376376e-8 * a11
			+ 4.2907451210093349687e-8 * a12
			+ 1.8531497952965519450e-8 * a13
			+ 2.4937657642194420731e-9 * a14
			- 6.1421743521933407055e-10 * a15
			- 3.7724994306416603112e-10 * a16
			- 7.0366311936526674419e-11 * a17
			+ 6.9617193221560459909e-12 * a18
			+ 7.5996918277492621596e-12 * a19
			+ 1.8429302830526136512e-12 * a20;
	}
	*/
}

/* Minus Log(P) for large Chi-square static */
TARGET double MinusLogPChi2(double x2, double df)
{
	return -LogGammaIncRegularized(df / 2, x2 / 2);
}

/* Minus Log(P) for large F-static */
TARGET double MinusLogPF(double f, double df1, double df2)
{
	return -LogBetaIncRegularized(df2 / 2, df1 / 2, df2 / (f * df1 + df2));
}

/* Minus Log(P) for large t-static */
TARGET double MinusLogPT(double t, double df)
{
	t = fabs(t);
	double p1 = TDistCDF(t, df), p2 = 1 - p1, p = 2 * std::min(p1, p2);

	if (p > 1e-100) return -log(p);

	return -LogBetaIncRegularized(0.5 * df, 0.5, df / (df + t * t));
}

/* Cumulative distribution function for Chi-square distribution */
TARGET double ChiSquareDistCDF(double x2, double df)
{
	return stats::pchisq(x2, df, false);
	/*
	if (!x2) return 1;
	if (x2 < 0.0) x2 = -x2;
	return GammaInc(df / 2.0, x2 / 2.0);
	*/
}

// Cumulative distribution function for standard normal distribution
TARGET double NormalDistCDF(double x)
{
	return stats::pnorm(x, false);
	/*
	double A1 = 0.31938153;
	double A2 = -0.356563782;
	double A3 = 1.781477937;
	double A4 = -1.821255978;
	double A5 = 1.330274429;
	double RSQRT2PI = 0.3989422804014327;
	double K = 1.0 / (1.0 + 0.2316419 * abs(x));
	double cnd = RSQRT2PI * exp(-0.5 * x * x) * (K * (A1 + K * (A2 + K * (A3 + K * (A4 + K * A5)))));
	if (x > 0) cnd = 1.0 - cnd;
	return cnd;
	*/
}

// Cumulative distribution function for T-distribution
TARGET double TDistCDF(double t, double df) 
{
	return stats::pt(t, df, false);
	/*
	double val = BetaInc(0.5 * df, 0.5, df / (df + t * t));

	if (t > 0)
		return 1 - 0.5 * val;
	else
		return 0.5 * val;
	*/
}

// Cumulative distribution function for F-distribution
TARGET double FDistCDF(double f, double df1, double df2)
{
	return stats::pf(f, df1, df2, false);
	//return BetaInc(df1 / 2, df2 / 2, f / (f + df2 / df1));
}

/* Input a vector, return proportion of grids with exp < threshold and two mininum indices */
TARGET double FindMinIndex(double* exp, int m, int& i1, int& i2, double threshold)
{
	int nle = 0;
	double minval1 = DBL_MAX, minval2 = DBL_MAX;
	i1 = i2 = 0xFFFFFFFF;
	for (int i = 0; i < m; ++i)
	{
		if (exp[i] < threshold)
			nle++;
		if (exp[i] < minval1)
		{
			minval2 = minval1;
			i2 = i1;

			minval1 = exp[i];
			i1 = i;
		}
		else if (exp[i] < minval2)
		{
			minval2 = exp[i];
			i2 = i;
		}
	}

	return nle / (double)m;
}

/* Input a column, return proportion of grids with exp < threshold and two mininum indices */
TARGET double FindMinIndex(double* exp, int m, int n, int& i1, int& j1, int& i2, int& j2, double threshold)
{
	int nle = 0;
	double minval1 = DBL_MAX, minval2 = DBL_MAX;
	i1 = j1 = i2 = j2 = 0xFFFFFFFF;
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			if (exp[i * n + j] < threshold)
				nle++;
			if (exp[i * n + j] < minval1)
			{
				minval2 = minval1;
				i2 = i1;
				j2 = j1;

				minval1 = exp[i * n + j];
				i1 = i;
				j1 = j;
			}
			else if (exp[i * n + j] < minval2)
			{
				minval2 = exp[i * n + j];
				i2 = i;
				j2 = j;
			}
		}
	}
	return nle / (double)(m * n);
}

/* strcmp two rows */
TARGET int CompareRow(double* obs, int n, int r1, int r2)
{
	for (int j = 0; j < n; ++j)
		if (obs[r1 * n + j] > obs[r2 * n + j]) return 1;
		else if (obs[r1 * n + j] < obs[r2 * n + j]) return -1;
	return 0;
}

/* strcmp two columns */
TARGET int CompareCol(double* obs, int m, int n, int c1, int c2)
{
	for (int i = 0; i < m; ++i)
		if (obs[i * n + c1] > obs[i * n + c2]) return 1;
		else if (obs[i * n + c1] < obs[i * n + c2]) return -1;
	return 0;
}

/* Sort columns */
TARGET void SortCol(double* obs, int m, int n)
{
	for (int i = 0; i < n; ++i)
		for (int j = i + 1; j < n; ++j)
			if (CompareCol(obs, m, n, i, j) > 0)
				for (int k = 0; k < m; ++k)
					Swap(obs[k * n + i], obs[k * n + j]);
}

/* Chi-square test, combine tables with any grid with an expectition < 5 */
TARGET void CombineTable(double* obs, int m, int n, double& g, int& df, double& p, bool test, double* obs2, double* exp, double* exp2, double* rowsum, double* colsum)
{
	{
		SortCol(obs, m, n);
		for (int i = 0; i < m; ++i)
			rowsum[i] = Sum(obs + i * n, n);
		for (int j = 0; j < n; ++j)
			colsum[j] = Sum(obs + j, m, n);

		double invtot = Sum(rowsum, m);

		if (invtot < 0.5) return;
		invtot = 1.0 / invtot;

		for (int i = 0; i < m; ++i)
			for (int j = 0; j < n; ++j)
				exp[i * n + j] = rowsum[i] * colsum[j] * invtot;
	}

	if (m < 2 || n < 2) df = 0;
	else for (;;)
	{
		for (int i = 0; i < m; ++i)
			rowsum[i] = Sum(obs + i * n, n);
		for (int j = 0; j < n; ++j)
			colsum[j] = Sum(obs + j, m, n);

		int i1 = 0, i2 = 0, j1 = 0, j2 = 0;
		double rate = FindMinIndex(exp, m, n, i1, j1, i2, j2, 5.0);

		if (exp[i1 * n + j1] > 1.0 && rate < 0.2)
		{
			df = (m - 1) * (n - 1);
			break;
		}

		if (m == 2 && n == 2)
		{
			df = 0;
			break;
		}

		FindMinIndex(rowsum, m, i1, i2, 5.0);
		FindMinIndex(colsum, n, j1, j2, 5.0);

		int m1 = m, n1 = n;
		if (m > 2 && n > 2)
		{
			if (m >= n)
			{
				m1--;
				j1 = j2 = 0x7FFFFFFF;
			}
			else
			{
				n1--;
				i1 = i2 = 0x7FFFFFFF;
			}
		}
		else if (m > 2)
		{
			m1--;
			j1 = j2 = 0x7FFFFFFF;
		}
		else if (n > 2)
		{
			n1--;
			i1 = i2 = 0x7FFFFFFF;
		}

		SetZero(obs2, m1 * n1);
		SetZero(exp2, m1 * n1);
		for (int i = 0; i < m; ++i)
		{
			int r = i == i1 ? i2 : i;
			r = r > i1 ? r - 1 : r;
			for (int j = 0; j < n; ++j)
			{
				int c = j == j1 ? j2 : j;
				c = c > j1 ? c - 1 : c;
				obs2[r * n1 + c] += obs[i * n + j];
				exp2[r * n1 + c] += exp[i * n + j];
			}
		}

		m = m1; n = n1;
		Swap(obs, obs2);
		Swap(exp, exp2);
	}

	if (df > 0)
	{
		g = 0;
		for (int i = 0; i < m; ++i)
			for (int j = 0; j < n; ++j)
				if (obs[i * n + j] > 0)
					g += 2 * obs[i * n + j] * log(obs[i * n + j] / exp[i * n + j]);

		if (test) 
			p = 1 - ChiSquareDistCDF(g, df);
	}
	else
	{
		g = 0;
		p = NAN;
	}
}

/* Factorial */
TARGET double Factorial(int n)
{
	static double f[] =
	{
		1.00000000000000E00, 1.00000000000000E00, 2.00000000000000E00, 6.00000000000000E00, 2.40000000000000E01, 1.20000000000000E02, 7.20000000000000E02, 5.04000000000000E03,
		4.03200000000000E04, 3.62880000000000E05, 3.62880000000000E06, 3.99168000000000E07, 4.79001600000000E08, 6.22702080000000E09, 8.71782912000000E10, 1.30767436800000E12,
		2.09227898880000E13, 3.55687428096000E14, 6.40237370572800E15, 1.21645100408832E17, 2.43290200817664E18, 5.10909421717094E19, 1.12400072777761E21, 2.58520167388850E22,
		6.20448401733239E23, 1.55112100433310E25, 4.03291461126606E26, 1.08888694504184E28, 3.04888344611714E29, 8.84176199373970E30, 2.65252859812191E32, 8.22283865417792E33,
		2.63130836933694E35, 8.68331761881189E36, 2.95232799039604E38, 1.03331479663861E40, 3.71993326789901E41, 1.37637530912263E43, 5.23022617466601E44, 2.03978820811974E46,
		8.15915283247898E47, 3.34525266131638E49, 1.40500611775288E51, 6.04152630633738E52, 2.65827157478845E54, 1.19622220865480E56, 5.50262215981209E57, 2.58623241511168E59,
		1.24139155925361E61, 6.08281864034268E62, 3.04140932017134E64, 1.55111875328738E66, 8.06581751709439E67, 4.27488328406002E69, 2.30843697339241E71, 1.26964033536583E73,
		7.10998587804863E74, 4.05269195048772E76, 2.35056133128288E78, 1.38683118545690E80, 8.32098711274139E81, 5.07580213877225E83, 3.14699732603879E85, 1.98260831540444E87,
		1.26886932185884E89, 8.24765059208247E90, 5.44344939077443E92, 3.64711109181887E94, 2.48003554243683E96, 1.71122452428141E98, 1.19785716699699E100, 8.50478588567862E101,
		6.12344583768861E103, 4.47011546151269E105, 3.30788544151939E107, 2.48091408113954E109, 1.88549470166605E111, 1.45183092028286E113, 1.13242811782063E115, 8.94618213078298E116,
		7.15694570462638E118, 5.79712602074737E120, 4.75364333701284E122, 3.94552396972066E124, 3.31424013456535E126, 2.81710411438055E128, 2.42270953836727E130, 2.10775729837953E132,
		1.85482642257398E134, 1.65079551609085E136, 1.48571596448176E138, 1.35200152767840E140, 1.24384140546413E142, 1.15677250708164E144, 1.08736615665674E146, 1.03299784882391E148,
		9.91677934870949E149, 9.61927596824822E151, 9.42689044888325E153, 9.33262154439442E155, 9.33262154439442E157, 9.42594775983836E159, 9.61446671503513E161, 9.90290071648618E163,
		1.02990167451456E166, 1.08139675824029E168, 1.14628056373471E170, 1.22652020319614E172, 1.32464181945183E174, 1.44385958320249E176, 1.58824554152274E178, 1.76295255109025E180,
		1.97450685722107E182, 2.23119274865981E184, 2.54355973347219E186, 2.92509369349301E188, 3.39310868445190E190, 3.96993716080872E192, 4.68452584975429E194, 5.57458576120761E196,
		6.68950291344913E198, 8.09429852527344E200, 9.87504420083360E202, 1.21463043670253E205, 1.50614174151114E207, 1.88267717688893E209, 2.37217324288005E211, 3.01266001845766E213,
		3.85620482362581E215, 4.97450422247729E217, 6.46685548922047E219, 8.47158069087881E221, 1.11824865119600E224, 1.48727070609068E226, 1.99294274616152E228, 2.69047270731805E230,
		3.65904288195255E232, 5.01288874827499E234, 6.91778647261948E236, 9.61572319694109E238, 1.34620124757175E241, 1.89814375907617E243, 2.69536413788816E245, 3.85437071718007E247,
		5.55029383273931E249, 8.04792605747199E251, 1.17499720439091E254, 1.72724589045464E256, 2.55632391787286E258, 3.80892263763057E260, 5.71338395644586E262, 8.62720977423324E264,
		1.31133588568345E267, 2.00634390509568E269, 3.08976961384735E271, 4.78914290146339E273, 7.47106292628289E275, 1.17295687942641E278, 1.85327186949373E280, 2.94670227249504E282,
		4.71472363599206E284, 7.59070505394722E286, 1.22969421873945E289, 2.00440157654530E291, 3.28721858553430E293, 5.42391066613158E295, 9.00369170577843E297, 1.50361651486500E300,
		2.52607574497320E302, 4.26906800900471E304, 7.25741561530800E306
	};
	if (n <= 170) return f[n];
	return NAN;
}

/* Binomial coefficient */
TARGET double Binomial(int n, int r)
{
	if (n <= 170)
		return Factorial(n) / Factorial(r) / Factorial(n - r);
	else
		return exp(LogGamma(n + 1) - LogGamma(r + 1) - LogGamma(n - r + 1));
}

/* Natural logarithm of binomial coefficient */
TARGET double LogBinomial(int n, int r)
{
	return LogGamma(n + 1) - LogGamma(r + 1) - LogGamma(n - r + 1);
}

/* Initialize BINOMIAL global variable */
TARGET void InitBinomial()
{
	for (int i = 0; i <= N_MAX_PLOIDY; ++i)
		for (int j = 0; j <= i; ++j)
			BINOMIAL[i][j] = (int)(Binomial(i, j) + 0.5);
}
