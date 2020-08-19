#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define RADIUS 6371.0088
#define RAD M_PI / 180

struct Coordinates {
	double lat;
	double lon;
};

int scan(int argc, char *argv, struct Coordinates *loc)
{
	int count = 0;
	if (argc == 1)
		count = scanf("%lf,%lf", &loc->lat, &loc->lon);
	else
		sscanf(argv, "%lf,%lf", &loc->lat, &loc->lon);
	return count;
}

int main(int argc, char **argv)
{
	if (argc == 2) {
		printf("Usage:\n\t%s [coordinate 1] [coordinate 2] ...\n", argv[0]);
		return 1;
	}

	struct Coordinates *location0 = malloc(sizeof(struct Coordinates));
	struct Coordinates *location1 = malloc(sizeof(struct Coordinates));

	double latdiff, londiff, a, b, c, total = 0;
	int i, j;

	scan(argc, argv[1], location0);
	if (location0->lat < -90 || location0->lat > 90 || location0->lon < -180 || location0->lon > 180) {
		printf("Coordinate %lf,%lf out of range.\nAbort.\n", location0->lat, location0->lon);
		goto freeing;
	}

	location0->lat *= RAD;
	location0->lon *= RAD;

	puts("{");

	for (i = 2; (argc > 1 && i < argc) || argc == 1; ++i) {
		if (scan(argc, argv[i], location1) == -1)
			break;
		if (location1->lat < -90 || location1->lat > 90 || location1->lon < -180 || location1->lon > 180) {
			printf("Coordinate %lf,%lf out of range.\nAbort.\n", location1->lat, location1->lon);
			goto freeing;
		}

		location1->lat *= RAD;
		location1->lon *= RAD;

		latdiff = location1->lat - location0->lat;
		londiff = location1->lon - location0->lon;

		a = sin(latdiff / 2);
		b = sin(londiff / 2);
		a *= a;
		a += cos(location1->lat) * cos(location0->lat) * b * b;
		c = 2 * atan2(sqrt(a), sqrt(1-a));
		c *= RADIUS;

		total += c;

		location0->lat = location1->lat;
		location0->lon = location1->lon;

		printf("  \"%d\": %lf,\n", i - 1, c);
	}

	free(location0);
	free(location1);

	printf("  \"total_distance\": %lf\n}\n", total);

	return 0;
freeing:
	free(location0);
	free(location1);
	return 1;
}
