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

	long double latdiff, londiff, a, b, c, total = 0;
	int i, j;

	scan(argc, argv[1], location);

	location->lat *= RAD;
	location->lon *= RAD;

	puts("{");

	for (i = 2; i < argc || argc == 1; ++i) {
		if (scan(argc, argv[i], (location + 1)) == -1)
			break;

		(location + 1)->lat *= RAD;
		(location + 1)->lon *= RAD;

		latdiff = (location + 1)->lat - location->lat;
		londiff = (location + 1)->lon - location->lon;

		a = sqr(sin(latdiff / 2));
		b = sinl(londiff / 2);
		a += cosl((location + 1)->lat) * cosl(location->lat) * sqr(b);
		c = 2 * atan2l(sqrtl(a), sqrtl(1-a)) * RADIUS;

		total += c;
		memcpy(location, location + 1, sizeof(struct Coordinates));

		printf("  \"%d\": %Lf,\n", i - 2, c);
	}

	free(location);

	printf("  \"total_distance\": %Lf\n}\n", total);

	return 0;
}
