#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define RADIUS 6371.0088
#define RAD M_PI / 180

struct Coordinates {
	double lat;
	double lon;
	int fail;
};

static inline void outrange(struct Coordinates loc)
{
	if (loc.lat < -90 || loc.lat > 90 || loc.lon < -180 || loc.lon > 180) {
		printf("Coordinate %lf,%lf out of range.\nAbort.\n", loc.lat, loc.lon);
		exit(1);
	}
}

static inline struct Coordinates scan(int argc, char *argv)
{
	struct Coordinates loc;
	int count;
	if (argc == 1)
		loc.fail = scanf("%lf,%lf", &loc.lat, &loc.lon);
	else
		sscanf(argv, "%lf,%lf", &loc.lat, &loc.lon);
	return loc;
}

int main(int argc, char **argv)
{
	if (argc == 2) {
		printf("Usage:\n\t%s [coordinate 1] [coordinate 2] ...\n", argv[0]);
		return 1;
	}

	struct Coordinates location0, location1;
	double latdiff, londiff, a, b, c, total = 0;
	int i, j;

	location0 = scan(argc, argv[1]);
	outrange(location0);

	location0.lat *= RAD;
	location0.lon *= RAD;

	puts("{");

	for (i = 2; (argc > 1 && i < argc) || argc == 1; ++i) {
		location1 = scan(argc, argv[i]);
		outrange(location1);

		if (location1.fail == -1)
			break;

		location1.lat *= RAD;
		location1.lon *= RAD;

		latdiff = location1.lat - location0.lat;
		londiff = location1.lon - location0.lon;

		a = sin(latdiff / 2);
		b = sin(londiff / 2);
		a *= a;
		a += cos(location1.lat) * cos(location0.lat) * b * b;
		c = 2 * atan2(sqrt(a), sqrt(1-a));
		c *= RADIUS;

		total += c;

		location0.lat = location1.lat;
		location0.lon = location1.lon;

		printf("  \"%d\": %lf,\n", i - 1, c);
	}

	printf("  \"total_distance\": %lf\n}\n", total);

	return 0;
}
