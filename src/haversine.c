#include <math.h>
#include "constants.h"
#include "functions.h"
#include "haversine.h"

long double haversine_distance(struct Coordinates *location)
{
	long double londiff, latdiff;
	long double a, b;

	latdiff = (location + 1)->lat - location->lat;
	londiff = (location + 1)->lon - location->lon;

	a = sqr(sin(latdiff / 2));
	b = sinl(londiff / 2);
	a += cosl((location + 1)->lat) * cosl(location->lat) * sqr(b);

	return 2 * atan2_modified(sqrtl(a), sqrtl(1-a)) * RADIUS;
}

long double haversine_bearing(struct Coordinates *start, struct Coordinates *end)
{
	long double londiff, y, x;

	londiff = end->lon - start->lon;
	y = sinl(londiff) * cosl(end->lat);
	x = cosl(start->lat) * sinl(end->lat) - sinl(start->lat) * cosl(end->lat) * cosl(londiff);
	return atan2_modified(y, x) / (RAD);
}
