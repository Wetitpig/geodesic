#include "constants.h"

#ifndef __HAVE_VINCENTY_H__

#define __HAVE_VINCENTY_H__

struct vincenty_result {
	long double distance;
	long double start;
	long double end;
};

void vincenty(struct vincenty_result *result, struct Coordinates *location);

#endif
