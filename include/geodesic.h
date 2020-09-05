#include <math.h>

#ifndef __HAVE_GEODESIC_MATH_H__

#define __HAVE_GEODESIC_MATH_H__

#define M_PI_L 3.14159265358979323846264338327950l

#define RADIUS 6371.0088l
#define RAD (M_PI_L / 180l)

#define RAD_MAJ 6378.137l
#define FLAT (1.0l/298.257223563l)

#define RAD_MIN (RAD_MAJ * (1.0l - FLAT))
#define FLAT_3 ((RAD_MAJ - RAD_MIN) / (RAD_MAJ + RAD_MIN))
#define ECC (FLAT * (2.0l - FLAT))
#define ECC2 (ECC / (1.0l - ECC))

struct Coordinates {
	long double lat;
	long double lon;
};

struct Vector {
	long double s;
	long double theta;
};

long double sqr(long double operand);
long double double_fac(int x);
long double atan2_modified(long double y, long double x);
long double normalise_a(long double x);
long double normalise_c(long double x);

#endif
