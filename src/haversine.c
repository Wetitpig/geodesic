#include <math.h>
#include "geodesic.h"
#include "haversine.h"

long double haversine_inverse_distance(struct Coordinates *location)
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

struct Coordinates haversine_direct(struct Coordinates *point, struct Vector *add)
{
	struct Coordinates result;
	long double delta = add->s / RADIUS;

	result.lat = asinl(sinl(point->lat) * cosl(delta) + cosl(point->lat) * sinl(delta) * cosl(add->theta));
	result.lon = point->lon + atan2_modified(sinl(add->theta) * sinl(delta) * cosl(point->lat), cosl(delta) - sinl(point->lat) * sinl(result.lat));
	result.lon = NORMALISE_C(result.lon);
	return result;
}
