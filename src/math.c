#include <math.h>
#include "geodesic.h"

long double sqr(long double operand)
{
	return operand * operand;
}

const long double double_fac[17] = {-1.0l/1.0l,
1.0l/2.0l,
1.0/8.0l,
1.0l/16.0l,
5.0l/128.0l,
7.0l/256.0l,
21.0l/1024.0l,
33.0l/2048.0l,
429.0l/32768.0l,
715.0l/65536.0l,
2431.0l/262144.0l,
4199.0l/524288.0l,
29393.0l/4194304.0l,
52003.0l/8388608.0l,
185725.0l/33554432.0l,
334305.0l/67108864.0l,
9694845.0l/2147483648.0l
};

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
