#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "constants.h"

#include "functions.h"

int main(int argc, char **argv)
{
	if (argc == 2) {
		printf("Usage:\n\t%s [coordinate 1] [coordinate 2] ...\n", argv[0]);
		return 1;
	}

	struct Coordinates *location = malloc(sizeof(struct Coordinates) * 2);

	double latdiff, londiff, a, b, c, total = 0;
	int i, j;

	scan(argc, argv[1], location);
	if (location->lat < -90 || location->lat > 90 || location->lon < -180 || location->lon > 180) {
		printf("Coordinate %lf,%lf out of range.\nAbort.\n", location->lat, location->lon);
		goto freeing;
	}

	location->lat *= RAD;
	location->lon *= RAD;

	puts("{");

	for (i = 2; i < argc || argc == 1; ++i) {
		if (scan(argc, argv[i], (location + 1)) == -1)
			break;
		if ((location + 1)->lat < -90 || (location + 1)->lat > 90 || (location + 1)->lon < -180 || (location + 1)->lon > 180) {
			printf("Coordinate %lf,%lf out of range.\nAbort.\n", (location + 1)->lat, (location + 1)->lon);
			goto freeing;
		}

		(location + 1)->lat *= RAD;
		(location + 1)->lon *= RAD;

		latdiff = (location + 1)->lat - location->lat;
		londiff = (location + 1)->lon - location->lon;

		a = sin(latdiff / 2);
		b = sin(londiff / 2);
		a *= a;
		a += cos((location + 1)->lat) * cos(location->lat) * b * b;
		c = 2 * atan2(sqrt(a), sqrt(1-a)) * RADIUS;

		total += c;
		memcpy(location, location + 1, sizeof(struct Coordinates));

		printf("  \"%d\": %lf,\n", i - 1, c);
	}

	free(location);

	printf("  \"total_distance\": %lf\n}\n", total);

	return 0;
freeing:
	free(location);
	return 1;
}
