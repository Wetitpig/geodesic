#include <stdio.h>
#include <stdlib.h>
#include "constants.h"
#include "functions.h"

int scan(int argc, char *argv, struct Coordinates *loc)
{
	int count = 0;
	if (argc == 1)
		count = scanf("%Lf,%Lf", &loc->lat, &loc->lon);
	else
		sscanf(argv, "%Lf,%Lf", &loc->lat, &loc->lon);

	if (loc->lat < -90 || loc->lat > 90 || loc->lon < -180 || loc->lon > 180) {
		printf("Coordinate %Lf,%Lf out of range.\nAbort.\n", loc->lat, loc->lon);
		exit(1);
        }

	return count;
}

void print(int order, long double s, long double start, long double end)
{
	printf("  \"%d\": {\n    \"distance\": %Lf,\n    \"start_azimuth\": %Lf,\n    \"end_azimuth\": %Lf\n  },\n", order, s, start, end);
	return;
}

