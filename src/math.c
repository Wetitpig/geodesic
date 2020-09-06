#include <math.h>
#include "geodesic.h"

long double sqr(long double operand)
{
	return operand * operand;
}

const long double fac_val[33] = {1.0l,
1.0l, 2.0l,
3.0l, 8.0l,
15.0l, 48.0l,
105.0l, 384.0l,
945.0l, 3840.0l,
10395.0l, 46080.0l,
135135.0l, 645120.0l,
2027025.0l, 10321920.0l,
34459425.0l, 185794560.0l,
654729075.0l, 3715891200.0l,
13749310575.0l, 81749606400.0l,
316234143225.0l, 1961990553600.0l,
7905853580625.0l, 51011754393600.0l,
213458046676875.0l, 1428329123020800.0l,
6190283353629375.0l, 42849873690624000.0l,
191898783962510625.0l, 1371195958099968000.0l
};

long double double_fac(int x)
{
	switch (x)
	{
	case -3:
		return -1;
		break;
	case -1:
		return 1;
		break;
	default:
		return fac_val[x];
		break;
	}
}

long double atan2_modified(long double y, long double x)
{
	if (x > 0)
		return atanl(y / x);
	else if (x < 0) {
		if (y < 0)
			return atanl(y / x) - M_PI_L;
		else
			return atanl(y / x) + M_PI_L;
	}
	else {
		if (y > 0)
			return M_PI_L / 2;
		else if (y < 0)
			return -M_PI_L / 2;
		else
			return NAN;
	}
}

long double normalise_a(long double x)
{
	return fmodl(x + 2 * M_PI_L, 2 * M_PI_L);
}

long double normalise_c(long double x)
{
	if (x > M_PI_L)
		return x - 2 * M_PI_L;
	else if (x < -M_PI_L)
		return x + 2 * M_PI_L;
	return x;
}
