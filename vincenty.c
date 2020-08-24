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
	long double londiff, lambda, oldvalue, u1, u2;
	long double ssig, csig, sig, salp, calp, cos2, C;
	long double usq, k1, a, b, dsig, s;

	long double d = 0;
	long double total = 0;

	scan(argc, argv[1], location);

	location->lat = location->lat * RAD;
	location->lon = location->lon * RAD;

	puts("{");

	for (int i = 2; i < argc || argc == 1; ++i) {
		if (scan(argc, argv[i], (location + 1)) == -1)
			break;

		(location + 1)->lat = (location + 1)->lat * RAD;
		(location + 1)->lon = (location + 1)->lon * RAD;
		londiff = (location + 1)->lon - location->lon;
		lambda = londiff;

		if (fabsl(location->lat) == M_PI_2)
			u1 = location->lat > 0 ? M_PI_2 : -M_PI_2;
		else
			u1 = atanl((1 - FLAT) * tanl(location->lat));

		if (fabsl((location + 1)->lat) == M_PI_2)
			u2 = (location + 1)->lat > 0 ? M_PI_2 : -M_PI_2;
		else
			u2 = atanl((1 - FLAT) * tanl((location + 1)->lat));

		do {
			oldvalue = lambda;

			ssig = sqrtl(sqr(cosl(u2) * sinl(lambda)) + sqr(cosl(u1) * sinl(u2) - sinl(u1) * cosl(u2) * cosl(lambda)));
			csig = sinl(u1) * sinl(u2) + cosl(u1) * cosl(u2) * cosl(lambda);

			sig = atan2l(ssig, csig);

			salp = cosl(u1) * cosl(u2) * sinl(lambda) / ssig;
			calp = 1 - salp * salp;
			cos2 = csig - 2 * sinl(u1) * sinl(u2) / calp;
			if (cos2 < -1 || isnan(cos2))
				cos2 = -1;

			C = FLAT / 16 * calp * (4 + FLAT * (4 - 3 * calp));

			lambda = londiff + (1 - C) * FLAT * salp * (sig + C * ssig * (cos2 + C * csig * (2 * cos2 * cos2 - 1)));

			if (fabsl(lambda) > M_PI) {
				d = 1;
				break;
			}
		} while (fabsl(oldvalue - lambda) >= powl(10,-12));

		if (d == 1) {
			londiff = (londiff > 0 ? M_PI : -M_PI) - londiff;
			lambda = 0;

			calp = 0.5;
			cos2 = 0;

			sig = M_PI - fabsl(u1 + u2);
			ssig = sinl(sig);
			csig = cosl(sig);

			do {
				oldvalue = salp;

				C = FLAT / 16 * calp * (4 + FLAT * (4 - 3 * calp));
				cos2 = csig - 2 * sinl(u1) * sinl(u2) / calp;

				d = (1 - C) * FLAT * (sig + C * ssig * (cos2 + C * csig * (2 * cos2 * cos2 - 1)));
				salp = (londiff - asinl(lambda)) / d;
				calp = 1 - salp * salp;

				lambda = salp * ssig / (cosl(u1) * cosl(u2));

				ssig = sqr(cosl(u2) * lambda) + sqr(cosl(u1) * sinl(u2) + sinl(u1) * cosl(u2) * sqrtl(1 - sqr(lambda)));
				csig = -sqrtl(1 - ssig);
				ssig = sqrtl(ssig);
				sig = M_PI - asinl(ssig);

			} while (fabsl(oldvalue - salp) >= powl(10,-12));
		}

		usq = calp * (RAD_MAJ * RAD_MAJ / (RAD_MIN * RAD_MIN) - 1);
		k1 = sqrtl(1 + usq);
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
}
