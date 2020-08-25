#include "geodesic.h"

#ifndef __HAVE_HAVERSINE_H__

#define __HAVE_HAVERSINE_H__

long double haversine_inverse_distance(struct Coordinates *location);
long double haversine_inverse_bearing(struct Coordinates *start, struct Coordinates *end);

struct Coordinates haversine_direct(struct Coordinates *point, struct Vector *add);

#endif
