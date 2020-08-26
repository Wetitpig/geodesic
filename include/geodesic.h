#include <math.h>

#ifndef __HAVE_GEODESIC_MATH_H__

#define __HAVE_GEODESIC_MATH_H__

#define RADIUS 6371.0088
#define RAD M_PI / 180

#define RAD_MAJ 6378.137
#define RAD_MIN 6356.752314245

#define FLAT 1/298.257223563

#define NORMALISE_A(x) fmodl(x + 360, 360)
#define NORMALISE_C(x) x - 2 * M_PI * (x > M_PI)

struct Coordinates {
	long double lat;
	long double lon;
};

struct Vector {
	long double s;
	long double theta;
};

long double sqr(long double operand);
long double atan2_modified(long double y, long double x);

#endif
