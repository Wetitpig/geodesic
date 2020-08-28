#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include "geodesic.h"
#include "io.h"

int scan_coordinates(FILE *in, struct Coordinates *loc)
{
	int count = fscanf(in, "%LF,%LF", &loc->lat, &loc->lon);

	if (count == 2 && (loc->lat < -90 || loc->lat > 90 || loc->lon < -180 || loc->lon > 180)) {
		fprintf(stderr, "Coordinate out of range: %LF, %LF.\n", loc->lat, loc->lon);
		error(NULL);
	}

	loc->lat = loc->lat * RAD;
	loc->lon = loc->lon * RAD;

	return count;
}

int scan_vector(FILE *in, struct Vector *vector)
{
	int count = fscanf(in, "%LF:%LF", &vector->s, &vector->theta);

	if (count == 2 && (vector->theta < 0 || vector->theta > 360)) {
		fprintf(stderr, "Bearing out of range: %LF.", vector->theta);
		error(NULL);
	}

	vector->theta = vector->theta * RAD;

	return count;
}

void start_print(char *out, int i)
{
	if (i == 1)
		*out = '[';
	else
		*out = ',';

	strcat(out, "\n  {\n");

	if (i % 10000 == 0)
		fprintf(stderr, "%d values processed.\n", i);
	return;
}

void error(char *msg)
{
	if (msg != NULL)
		fputs(msg, stderr);
	fputs("\nAbort.\n", stderr);
	exit(1);
}
