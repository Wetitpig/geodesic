#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"

#include "functions.h"

int main(int argc, char **argv)
{
	if (argc == 2) {
		printf("Usage:\n\t%s [coordinate 1] [coordinate 2] ...\n", argv[0]);
		return 1;
	}

	struct Coordinates *location = malloc(sizeof(struct Coordinates) * 2);
	double londiff, lambda, oldlambda, u1, u2;
	double ssig, csig, sig, salp, calp, cos2, C;
	double usq, a, b, dsig, s;
	double total = 0;

	scan(argc, argv[1], location);
	if (location->lat < -90 || location->lat > 90 || location->lon < -180 || location->lon > 180) {
		printf("Coordinate %lf,%lf out of range.\nAbort.\n", location->lat, location->lon);
		goto freeing;
	}

	location->lat *= RAD;
	location->lon *= RAD;

	puts("{");

	for (int i = 2; i < argc || argc == 1; ++i) {
		if (scan(argc, argv[i], (location + 1)) == -1)
			break;
		if ((location + 1)->lat < -90 || (location + 1)->lat > 90 || (location + 1)->lon < -180 || (location + 1)->lon > 180) {
			printf("Coordinate %lf,%lf out of range.\nAbort.\n", (location + 1)->lat, (location + 1)->lon);
			goto freeing;
		}

		(location + 1)->lat *= RAD;
		(location + 1)->lon *= RAD;

		londiff = (location + 1)->lon - location->lon;
		lambda = londiff;
		u1 = atan((1 - FLATTENING) * tan(location->lat));
		u2 = atan((1 - FLATTENING) * tan((location + 1)->lat));

		do {
			oldlambda = lambda;

			ssig = cos(u2) * sin(lambda);
			csig = cos(u1) * sin(u2) - sin(u1) * cos(u2) * cos(lambda);
			ssig *= ssig;
			csig *= csig;
			ssig = sqrt(ssig + csig);

			csig = sin(u1) * sin(u2) + cos(u1) * cos(u2) * cos(lambda);

			sig = atan2(ssig, csig);

			salp = cos(u1) * cos(u2) * sin(lambda) / ssig;
			calp = 1 - salp * salp;
			cos2 = csig - 2 * sin(u1) * sin(u2) / calp;

			C = FLATTENING / 16 * calp * (4 + FLATTENING * (4 - 3 * calp));
			lambda = londiff + (1 - C) * FLATTENING * salp * (sig + C * ssig * (cos2 + C * csig * (2 * cos2 * cos2 - 1)));
		} while (fabs(oldlambda - lambda) >= pow(10,-12));

		usq = calp * (RAD_MAJ * RAD_MAJ - RAD_MIN * RAD_MIN) / (RAD_MIN * RAD_MIN);
		a = 1 + usq / 16384 * (4096 + usq * (-768 + usq * (320 - 175 * usq)));
		b = usq / 1024 * (256 + usq * (-128 + usq * (74 - 47 * usq)));

		dsig = b * ssig * (cos2 + b / 4 * (csig * (2 * cos2 * cos2 - 1) - b / 6 * cos2 * (4 * ssig * ssig - 3) * (4 * cos2 * cos2 - 3)));
		s = RAD_MIN * a * (sig - dsig);

		total += s;
		memcpy(location, location + 1, sizeof(struct Coordinates));
		printf("  \"%d\": %lf,\n", i - 1, s);
	}

	free(location);
	printf("  \"total_distance\": %lf\n}\n", total);

	return 0;

freeing:
	free(location);
	return 1;
}
