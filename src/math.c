#include <math.h>
#include "geodesic.h"

long double sqr(long double operand)
{
	return operand * operand;
}

long double double_fac(int x)
{
	if (x == -3)
		return -1;
	else {
		long double y = 1;
		while (x > 0) {
			y *= x;
			x -= 2;
		}
		return y;
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
