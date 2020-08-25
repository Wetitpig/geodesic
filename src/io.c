#include <stdio.h>
#include <stdlib.h>
#include "geodesic.h"
#include "mathio.h"

int scan_coordinates(int argc, char *argv, struct Coordinates *loc)
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

	loc->lat = loc->lat * RAD;
	loc->lon = loc->lon * RAD;

	return count;
}

int scan_vector(int argc, char *argv, struct Vector *vector)
{
	int count = 0;
	if (argc == 1)
		count = scanf("%Lf:%Lf", &vector->s, &vector->theta);
	else
		count = sscanf(argv, "%Lf:%Lf", &vector->s, &vector->theta);

	if (vector->theta < 0 || vector->theta > 360) {
		printf("Bearing %Lf out of range. Abort.\n", vector->theta);
		exit(1);
	}

	vector->theta = vector->theta * RAD;

	return count;
}

