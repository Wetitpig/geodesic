#include "geodesic.h"

#ifndef __HAVE_MPBLOCK_H__

#define __HAVE_MPBLOCK_H__

void mpblock_area(struct Coordinates *vertex, int i, long double s, long double a, long double *res);
void mpblock_inverse(struct Coordinates *location, struct Coordinates *location2, long double *res);

#endif
