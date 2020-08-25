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
		count = sscanf(argv, "%Lf,%Lf", &loc->lat, &loc->lon);
	if (count == 2) {
		if (loc->lat < -90 || loc->lat > 90 || loc->lon < -180 || loc->lon > 180) {
			printf("Coordinate out of range: %Lf, %Lf.\nAbort.\n", loc->lat, loc->lon);
			exit(1);
	        }
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

	if (count == 2) {
		if (vector->theta < 0 || vector->theta > 360) {
			printf("Bearing out of range: %Lf. Abort.\n", vector->theta);
			exit(1);
		}
	}

	vector->theta = vector->theta * RAD;

	return count;
}

void start_print(int i)
{
	if (i == 1)
		puts("[");
	else
		puts(",");

	puts("  {");
	return;
}
