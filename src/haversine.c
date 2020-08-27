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
	a = a + cosl((location + 1)->lat) * cosl(location->lat) * sqr(b);

	return 2 * atan2_modified(sqrtl(a), sqrtl(1-a)) * RADIUS;
}

long double haversine_bearing(struct Coordinates *start, struct Coordinates *end)
{
	long double londiff, y, x;

	londiff = end->lon - start->lon;
	y = sinl(londiff) * cosl(end->lat);
	x = cosl(start->lat) * sinl(end->lat) - sinl(start->lat) * cosl(end->lat) * cosl(londiff);
	return normalise_a(atan2_modified(y, x));
}

struct Coordinates haversine_direct(struct Coordinates *point, struct Vector *add)
{
	struct Coordinates result;
	long double delta = add->s / RADIUS, y, x;

	result.lat = sinl(point->lat) * cosl(delta) + cosl(point->lat) * sinl(delta) * cosl(add->theta);

	x = cosl(delta) * sinl(point->lat) * result.lat;
	result.lat = asinl(result.lat);
	y = sinl(add->theta) * sinl(delta) * cosl(point->lat);
	if (x == 0 && y == 0) {
		y = cosl(delta) - sinl(point->lat) * sinl(result.lat);
		x = cosl(point->lat) * cosl(result.lat);
		result.lon = acosl(y / x);
		if (result.lon < 0)
			result.lon = M_PI_2 - result.lon;
	}
	else
		result.lon = atan2_modified(y, x);

	result.lon = normalise_c(result.lon + point->lon);
	return result;
}
