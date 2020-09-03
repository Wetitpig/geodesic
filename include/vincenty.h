#include "geodesic.h"

#ifndef __HAVE_VINCENTY_H__

#define __HAVE_VINCENTY_H__

long double sin_reduced_latitude(long double lat);
long double cos_reduced_latitude(long double lat);

void vincenty_inverse(struct Coordinates *location, struct Coordinates *location2, long double *res, int count);
void vincenty_direct(struct Coordinates *point, struct Vector *add, long double *res);

#endif
