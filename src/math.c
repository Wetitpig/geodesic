#include <math.h>
#include "geodesic.h"

long double sqr(long double operand)
{
	return operand * operand;
}


long double atan2_modified(long double y, long double x)
{
	if (x > 0)
		return atanl(y / x);
	else if (x < 0) {
		if (y < 0)
			return atanl(y / x) - M_PI;
		else
			return atanl(y / x) + M_PI;
	}
	else {
		if (y > 0)
			return M_PI_2;
		else if (y < 0)
			return -M_PI_2;
		else
			return NAN;
	}
}

long double normalise_a(long double x)
{
	return fmodl(x + 2 * M_PI, 2 * M_PI);
}

long double normalise_c(long double x)
{
	if (x > M_PI)
		return x - 2 * M_PI;
	else if (x < -M_PI)
		return x + 2 * M_PI;
	return x;
}
