#include <math.h>
#include "geodesic.h"
#include "vincenty.h"

long double reduced_latitude(long double lat)
{
	if (lat == M_PI_2)
		return M_PI_2;
	else if (lat == -M_PI_2)
		return -M_PI_2;
	else
		return atanl((1 - FLAT) * tanl(lat));
}

struct Coordinates helmert(long double calp)
{
	struct Coordinates ab;
	long double usq, k1;

	usq = calp * (sqr(RAD_MAJ) / sqr(RAD_MIN) - 1);
	k1 = sqrtl(1 + usq);
	k1 = (k1 - 1) / (k1 + 1);
	ab.lat = (1 + sqr(k1) / 4) / (1 - k1);
	ab.lon = k1 * (1 - 3 / 8 * sqr(k1));
	return ab;
}

static inline long double Cfromcalp(long double calp)
{
	return FLAT / 16 * calp * (4 + FLAT * (4 - 3 * calp));
}

struct vincenty_result vincenty_inverse(struct Coordinates *location, struct Coordinates *location2)
{
	long double londiff, lambda, u1, u2;
	long double ssig, csig, sig, salp, calp, cos2, c;
	long double a, b, dsig, s;
	int i = 0;

	long double oldvalue[2];
	oldvalue[0] = 0;

	struct Coordinates ab;

	long double d = 0;

	londiff = location2->lon - location->lon;
	lambda = londiff;

	u1 = reduced_latitude(location->lat);
	u2 = reduced_latitude(location2->lat);

	do {
		oldvalue[1] = oldvalue[0];
		oldvalue[0] = lambda;

		ssig = sqrtl(sqr(cosl(u2) * sinl(lambda)) + sqr(cosl(u1) * sinl(u2) - sinl(u1) * cosl(u2) * cosl(lambda)));
		csig = sinl(u1) * sinl(u2) + cosl(u1) * cosl(u2) * cosl(lambda);

		sig = atan2_modified(ssig, csig);

		salp = cosl(u1) * cosl(u2) * sinl(lambda) / ssig;
		calp = 1 - sqr(salp);
		cos2 = csig - 2 * sinl(u1) * sinl(u2) / calp;
		if (cos2 < -1 || isnan(cos2))
			cos2 = -1;

		c = Cfromcalp(calp);

		lambda = londiff + (1 - c) * FLAT * salp * (sig + c * ssig * (cos2 + c * csig * (2 * sqr(cos2) - 1)));

		if (i++ > 8192 || (oldvalue[1] == lambda && fabsl(lambda) > M_PI)) {
			d = 1;
			break;
		}
	} while (fabsl(oldvalue[0] - lambda) >= powl(10,-15));

	if (d == 1) {
		londiff = (londiff > 0 ? M_PI : -M_PI) - londiff;
		lambda = 0;

		i = 16384;
		calp = 0.5;
		cos2 = 0;

		sig = M_PI - fabsl(u1 + u2);
		ssig = sinl(sig);
		csig = cosl(sig);

		do {
			oldvalue[0] = salp;
			c = Cfromcalp(calp);
			cos2 = csig - 2 * sinl(u1) * sinl(u2) / calp;

			d = (1 - c) * FLAT * (sig + c * ssig * (cos2 + c * csig * (2 * sqr(cos2) - 1)));
			salp = (londiff - asinl(lambda)) / d;
			calp = 1 - sqr(salp);

			lambda = salp * ssig / (cosl(u1) * cosl(u2));

			ssig = sqr(cosl(u2) * lambda) + sqr(cosl(u1) * sinl(u2) + sinl(u1) * cosl(u2) * sqrtl(1 - sqr(lambda)));
			csig = -sqrtl(1 - ssig);
			ssig = sqrtl(ssig);
			sig = M_PI - asinl(ssig);

		} while (fabsl(oldvalue[0] - salp) >= powl(10,-15));
	}

	ab = helmert(calp);
	a = ab.lat;
	b = ab.lon;

	dsig = b * ssig * (cos2 + b / 4 * (csig * (2 * sqr(cos2) - 1) - b / 6 * cos2 * (4 * sqr(ssig) - 3) * (4 * sqr(cos2) - 3)));
	s = RAD_MIN * a * (sig - dsig);

	if (i != 16384) {
		a = cosl(u2) * sinl(lambda);
		b = cosl(u1) * sinl(u2) - sinl(u1) * cosl(u2) * cosl(lambda);

		c = cosl(u1) * sinl(lambda);
		d = cosl(u1) * sinl(u2) * cosl(lambda) - sinl(u1) * cosl(u2);
	}
	else {
		a = salp / cosl(u1);
		b = sqrtl(1 - sqr(a));
		if (cosl(u1) * sinl(u2) + sinl(u1) * cosl(u2) * cosl(lambda) < 0)
			b = b * -1;

		c = salp;
		d = -sinl(u1) * ssig + cosl(u1) * csig * b;
	}

	struct vincenty_result result;
	result.distance = s;
	result.start = normalise_a(atan2_modified(a, b));
	result.end = normalise_a(atan2_modified(c, d));
	return result;
}

struct vincenty_result vincenty_direct(struct Coordinates *point, struct Vector *add)
{
	long double u1, s1, salp, a, b;
	long double sigma, s2, dsig, oldvalue;
	long double lambda, c;

	struct Coordinates ab;

	u1 = reduced_latitude(point->lat);
	s1 = atan2_modified(tanl(u1), cosl(add->theta));
	salp = cosl(u1) * sinl(add->theta);

	ab = helmert(1 - sqr(salp));
	a = ab.lat;
	b = ab.lon;
	sigma = add->s / (RAD_MIN * a);

	do {
		oldvalue = sigma;

		s2 = 2 * s1 + sigma;
		dsig = b * sinl(sigma) * (cosl(s2) + b / 4 * (cosl(sigma) * (2 * sqr(cosl(s2)) - 1) - b / 6 * cosl(s2) * (4 * sqr(sinl(sigma)) - 3) * (4 * sqr(cosl(s2)) - 3)));
		sigma = add->s / (RAD_MIN * a) + dsig;
	} while (fabsl(oldvalue - sigma) >= powl(10, -15));

	struct vincenty_result result;
	result.distance = atan2_modified(sinl(u1) * cosl(sigma) + cosl(u1) * sinl(sigma) * cosl(add->theta), (1 - FLAT) * sqrtl(sqr(salp) + sqr(sinl(u1) * sinl(sigma) - cosl(u1) * cosl(sigma) * cosl(add->theta))));
	lambda = atan2_modified(sinl(sigma) * sinl(add->theta), cosl(u1) * cosl(sigma) - sinl(u1) * sinl(sigma) * cosl(add->theta));

	c = Cfromcalp(1 - sqr(salp));
	result.start = lambda - (1 - c) * FLAT * salp * (sigma + c * sinl(sigma) * (cosl(s2) + c * cosl(sigma) * (2 * sqr(cosl(s2)) - 1)));
	result.start = result.start + point->lon;
	result.start = normalise_c(result.start);

	result.end = normalise_a(atan2_modified(salp, -sinl(u1) * sinl(sigma) + cosl(u1) * cosl(sigma) * cosl(add->theta)));

	return result;
}

