#include <math.h>
#include <stdio.h>

#define RADIUS 6371.0088
#define RAD M_PI / 180

static inline int valid(double lat, double lon)
{
	return lat < -90 || lat > 90 || lon < -180 || lon > 180;
}


static inline void outrange(double lat, double lon)
{
	printf("Coordinate %lf,%lf out of range.\nAbort.\n", lat, lon);
	return;
}

int main(int argc, char **argv)
{
	if (argc < 3) {
		printf("Usage:\n\t%s [coordinate 1] [coordinate 2] ...\n", argv[0]);
		return 1;
	}

	double lat[2], lon[2];
	double latdiff, londiff, a, b, c, total = 0;
	int i;

	sscanf(argv[1], "%lf,%lf", &lat[0], &lon[0]);
	if (valid(lat[0], lon[0])) {
		outrange(lat[0], lon[0]);
		return 1;
	}
	lat[0] *= RAD;
	lon[0] *= RAD;

	puts("{");

	for (i = 2; i < argc; ++i) {
		sscanf(argv[i], "%lf,%lf", &lat[1], &lon[1]);
		if (valid(lat[1], lon[1])) {
			outrange(lat[1], lon[1]);
			return 1;
		}

		lat[1] *= RAD;
		lon[1] *= RAD;

		latdiff = lat[1] - lat[0];
		londiff = lon[1] - lon[0];

		a = sin(latdiff / 2);
		b = sin(londiff / 2);
		a *= a;
		a += cos(lat[1]) * cos(lat[0]) * b * b;
		c = 2 * atan2(sqrt(a), sqrt(1-a));
		c *= RADIUS;

		total += c;

		lat[0] = lat[1];
		lon[0] = lon[1];

		printf("  \"%d\": %lf,\n", i - 1, c);
	}

	printf("  \"total_distance\": %lf\n}\n", total);

	return 0;
}
