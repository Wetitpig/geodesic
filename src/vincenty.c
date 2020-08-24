#include <math.h>
#include "constants.h"
#include "functions.h"
#include "vincenty.h"

void vincenty(struct vincenty_result *result, struct Coordinates *location)
{
	long double londiff, lambda, oldvalue, u1, u2;
	long double ssig, csig, sig, salp, calp, cos2, C;
	long double usq, k1, a, b, dsig, s;
	long double start, end;

	long double d = 0;

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

		sig = atan2_modified(ssig, csig);

		salp = cosl(u1) * cosl(u2) * sinl(lambda) / ssig;
		calp = 1 - sqr(salp);
		cos2 = csig - 2 * sinl(u1) * sinl(u2) / calp;
		if (cos2 < -1 || isnan(cos2))
			cos2 = -1;

		C = FLAT / 16 * calp * (4 + FLAT * (4 - 3 * calp));

		lambda = londiff + (1 - C) * FLAT * salp * (sig + C * ssig * (cos2 + C * csig * (2 * sqr(cos2) - 1)));

		if (fabsl(lambda) > M_PI) {
			d = 1;
			break;
		}
	} while (fabsl(oldvalue - lambda) >= powl(10,-15));

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

			d = (1 - C) * FLAT * (sig + C * ssig * (cos2 + C * csig * (2 * sqr(cos2) - 1)));
			salp = (londiff - asinl(lambda)) / d;
			calp = 1 - sqr(salp);

			lambda = salp * ssig / (cosl(u1) * cosl(u2));

			ssig = sqr(cosl(u2) * lambda) + sqr(cosl(u1) * sinl(u2) + sinl(u1) * cosl(u2) * sqrtl(1 - sqr(lambda)));
			csig = -sqrtl(1 - ssig);
			ssig = sqrtl(ssig);
			sig = M_PI - asinl(ssig);

		} while (fabsl(oldvalue - salp) >= powl(10,-15));
	}

	usq = calp * (sqr(RAD_MAJ) / sqr(RAD_MIN) - 1);
	k1 = sqrtl(1 + usq);
	k1 = (k1 - 1) / (k1 + 1);
	a = (1 + sqr(k1) / 4) / (1 - k1);
	b = k1 * (1 - 3 / 8 * sqr(k1));

	dsig = b * ssig * (cos2 + b / 4 * (csig * (2 * sqr(cos2) - 1) - b / 6 * cos2 * (4 * sqr(ssig) - 3) * (4 * sqr(cos2) - 3)));
	s = RAD_MIN * a * (sig - dsig);

	if (fabsl(oldvalue - salp) > fabsl(oldvalue - lambda)) {
		a = cosl(u2) * sinl(lambda);
		b = cosl(u1) * sinl(u2) - sinl(u1) * cosl(u2) * cosl(lambda);
		start = NORMALISE(atan2_modified(a, b) / (RAD));

		a = cosl(u1) * sinl(lambda);
		b = cosl(u1) * sinl(u2) * cosl(lambda) - sinl(u1) * cosl(u2);
		end = NORMALISE(atan2_modified(a, b) / (RAD));
	}
	else {
		a = salp / cosl(u1);
		b = sqrtl(1 - sqr(a));
		if (cosl(u1) * sinl(u2) + sinl(u1) * cosl(u2) * cosl(lambda) < 0)
			b = b * -1;

		start = fmodl(atan2_modified(a, b) / (RAD) + 360, 360);
		end = fmodl(atan2_modified(salp, -sinl(u1) * ssig + cosl(u1) * csig * b) / (RAD) + 360, 360);
	}

	result->distance = s;
	result->start = start;
	result->end = end;
	return;
}
