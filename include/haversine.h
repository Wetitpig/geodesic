#include "constants.h"

#ifndef __HAVE_HAVERSINE_H__

#define __HAVE_HAVERSINE_H__

long double haversine_distance(struct Coordinates *location);
long double haversine_bearing(struct Coordinates *start, struct Coordinates *end);

#endif
