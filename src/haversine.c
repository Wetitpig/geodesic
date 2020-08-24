#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "constants.h"

#include "functions.h"

long double distance(struct Coordinates *location)
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

long double bearing(struct Coordinates *start, struct Coordinates *end)
{
	long double londiff, y, x;

	londiff = end->lon - start->lon;
	y = sinl(londiff) * cosl(end->lat);
	x = cosl(start->lat) * sinl(end->lat) - sinl(start->lat) * cosl(end->lat) * cosl(londiff);
	return atan2_modified(y, x) / (RAD);
}

int main(int argc, char **argv)
{
	if (argc == 2) {
		printf("Usage:\n\t%s [coordinate 1] [coordinate 2] ...\n", argv[0]);
		return 1;
	}

	struct Coordinates *location = malloc(sizeof(struct Coordinates) * 2);

	long double c, total = 0, start, end;
	int i, j;

	scan(argc, argv[1], location);

	location->lat = location->lat * RAD;
	location->lon = location->lon * RAD;

	puts("{");

	for (i = 2; i < argc || argc == 1; ++i) {
		if (scan(argc, argv[i], (location + 1)) == -1)
			break;

		(location + 1)->lat = (location + 1)->lat * RAD;
		(location + 1)->lon = (location + 1)->lon * RAD;

		c = distance(location);
		start = NORMALISE(bearing(location, location + 1));
		end = NORMALISE(bearing(location + 1, location) - 180);

		total += c;
		memcpy(location, location + 1, sizeof(struct Coordinates));

		print(i - 2, c, start, end);
	}

	free(location);

	printf("  \"total_distance\": %Lf\n}\n", total);

	return 0;
}
