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
	long double londiff, lambda, oldlambda, u1, u2;
	long double ssig, csig, sig, salp, calp, cos2, C;
	long double usq, k1, a, b, dsig, s;
	long double total = 0;

	scan(argc, argv[1], location);
	if (location->lat < -90 || location->lat > 90 || location->lon < -180 || location->lon > 180) {
		printf("Coordinate %Lf,%Lf out of range.\nAbort.\n", location->lat, location->lon);
		goto freeing;
	}

	location->lat *= RAD;
	location->lon *= RAD;

	puts("{");

	for (int i = 2; i < argc || argc == 1; ++i) {
		if (scan(argc, argv[i], (location + 1)) == -1)
			break;
		if ((location + 1)->lat < -90 || (location + 1)->lat > 90 || (location + 1)->lon < -180 || (location + 1)->lon > 180) {
			printf("Coordinate %Lf,%Lf out of range.\nAbort.\n", (location + 1)->lat, (location + 1)->lon);
			goto freeing;
		}

		(location + 1)->lat *= RAD;
		(location + 1)->lon *= RAD;

		londiff = (location + 1)->lon - location->lon;
		lambda = londiff;
		u1 = atanl((1 - FLAT) * tanl(location->lat));
		u2 = atanl((1 - FLAT) * tanl((location + 1)->lat));

		do {
			oldlambda = lambda;

			ssig = cosl(u2) * sinl(lambda);
			csig = cosl(u1) * sinl(u2) - sinl(u1) * cosl(u2) * cosl(lambda);
			ssig *= ssig;
			csig *= csig;
			ssig = sqrtl(ssig + csig);

			csig = sinl(u1) * sinl(u2) + cosl(u1) * cosl(u2) * cosl(lambda);

			sig = atan2l(ssig, csig);

			salp = cosl(u1) * cosl(u2) * sinl(lambda) / ssig;
			calp = 1 - salp * salp;
			cos2 = csig - 2 * sinl(u1) * sinl(u2) / calp;

			C = FLAT / 16 * calp * (4 + FLAT * (4 - 3 * calp));
			lambda = londiff + (1 - C) * FLAT * salp * (sig + C * ssig * (cos2 + C * csig * (2 * cos2 * cos2 - 1)));
		} while (fabsl(oldlambda - lambda) >= powl(10,-12));

		usq = calp * (RAD_MAJ * RAD_MAJ / (RAD_MIN * RAD_MIN) - 1);
		k1 = sqrt(1 + usq);
		k1 = (k1 - 1) / (k1 + 1);
		a = (1 + k1 * k1 / 4) / (1 - k1);
		b = k1 * (1 - 3 / 8 * k1 * k1);

		dsig = b * ssig * (cos2 + b / 4 * (csig * (2 * cos2 * cos2 - 1) - b / 6 * cos2 * (4 * ssig * ssig - 3) * (4 * cos2 * cos2 - 3)));
		s = RAD_MIN * a * (sig - dsig);

		total += s;
		memcpy(location, location + 1, sizeof(struct Coordinates));
		printf("  \"%d\": %Lf,\n", i - 2, s);
	}

	free(location);
	printf("  \"total_distance\": %Lf\n}\n", total);

	return 0;

freeing:
	free(location);
	return 1;
}
