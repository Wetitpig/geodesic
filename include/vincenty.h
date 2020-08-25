#include "geodesic.h"

#ifndef __HAVE_VINCENTY_H__

#define __HAVE_VINCENTY_H__

struct vincenty_result {
	long double distance;
	long double start;
	long double end;
};

struct vincenty_result vincenty_inverse(struct Coordinates *location);
struct vincenty_result vincenty_direct(struct Coordinates *point, struct Vector *add);

#endif
