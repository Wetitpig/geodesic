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
