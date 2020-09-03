#include <math.h>
#include "geodesic.h"
#include "vincenty.h"

const long double Alookup[5] = {25.0/16384, 1.0/256, 1.0/64, 1.0/4, 1};

const long double Blookup[8][4] = {
{-1.0/2,3.0/16,-1.0/32,19.0/2048},
{-1.0/16,1.0/32,-9.0/2048,7.0/4096},
{-1.0/48,3.0/256,-3.0/2048},
{-5.0/512,3.0/512,-11.0/16384},
{-7.0/1280,7.0/2048},
{-7.0/2048,9.0/4096},
{-33.0/14336},
{-429.0/212144.0}
};

const long double B1lookup[8][4] = {
{1.0/2,-9.0/32,205.0/1536,-4879.0/73728},
{5.0/16,-37.0/96,1335.0/4096,-86171.0/363840},
{29.0/96,-75.0/128,2901.0/4096},
{539.0/1536,-2391.0/2560,1082857.0/737280},
{3467.0/7680.0,-28223.0/18432},
{38081.0/61440,-733437.0/286720},
{459485.0/516096},
{109167851.0/82575360}
};

long double tan_reduced_latitude(long double lat)
{
	return (1 - FLAT) * tanl(lat);
}

long double reduced_latitude(long double lat)
{
	if (isnan(tanl(lat)))
		return copysignl(M_PI_2, lat);
	else
		return atanl(tan_reduced_latitude(lat));
}

long double sin_reduced_latitude(long double lat)
{
	if (isnan(tanl(lat)))
		return copysignl(1, lat);
	else
		return copysignl(1 / hypotl(1, 1 / tan_reduced_latitude(lat)), lat);
}

long double cos_reduced_latitude(long double lat)
{
	if (isnan(tanl(lat)))
		return 0;
	else
		return 1 / hypotl(1, tan_reduced_latitude(lat));
}

long double helmertA(long double calp)
{
	long double A = 0;
	long double usq, k1;
	int k;

	usq = calp * ECC2;
	k1 = sqrtl(1 + usq);
	k1 = (k1 - 1) / (k1 + 1);
	for (k = 0; k < 5; k++) {
		A *= sqr(k1);
		A += Alookup[k];
	}
	A /= (1 - k1);
	return A;
}

long double helmertB(long double calp, long double sig)
{
	long double B = 0;
	long double usq, k1, c;
	int k, j, i;

	usq = calp * ECC2;
	k1 = sqrtl(1 + usq);
	k1 = (k1 - 1) / (k1 + 1);
	for (k = 1; k < 9; k++) {
		c = 0;
		i = 0;
		for (j = k; j < 9; j += 2)
			c += Blookup[k - 1][i++] * powl(k1, j);
		B += c * sinl(2 * k * sig);
	}
	return B;
}

long double helmertB1(long double calp, long double sig)
{
	long double B = 0;
	long double usq, k1, c;
	int k, j, i;

	usq = calp * ECC2;
	k1 = sqrtl(1 + usq);
	k1 = (k1 - 1) / (k1 + 1);
	for (k = 1; k < 9; k++) {
		c = 0;
		i = 0;
		for (j = k; j < 9; j += 2)
			c += B1lookup[k - 1][i++] * powl(k1, j);
		B += c * sinl(2 * k * sig);
	}
	return B;
}

static inline long double Cfromcalp(long double calp)
{
	return FLAT / 16 * calp * (4 + FLAT * (4 - 3 * calp));
}

void vincenty_inverse(struct Coordinates *location, struct Coordinates *location2, long double *res, int count)
{
	long double londiff, lambda, u1, u2, cu1, cu2, su1, su2;
	long double ssig, csig, sig, salp, calp, cos2, c;
	long double a, b, dsig, s;
	int i = 0;

	long double oldvalue[2];
	oldvalue[0] = 0;

	long double d = 0;

	londiff = location2->lon - location->lon;
	lambda = londiff;

	u1 = reduced_latitude(location->lat);
	u2 = reduced_latitude(location2->lat);
	cu1 = cos_reduced_latitude(location->lat);
	cu2 = cos_reduced_latitude(location2->lat);
	su1 = sin_reduced_latitude(location->lat);
	su2 = sin_reduced_latitude(location2->lat);

	do {
		oldvalue[1] = oldvalue[0];
		oldvalue[0] = lambda;

		ssig = hypotl(cu2 * sinl(lambda), cu1 * su2 - su1 * cu2 * cosl(lambda));
		csig = su1 * su2 + cu1 * cu2 * cosl(lambda);

		sig = atan2_modified(ssig, csig);

		salp = cu1 * cu2 * sinl(lambda) / ssig;
		calp = 1 - sqr(salp);
		cos2 = csig - 2 * su1 * su2 / calp;
		if (cos2 < -1 || isnan(cos2))
			cos2 = -1;

		c = Cfromcalp(calp);

		lambda = londiff + (1 - c) * FLAT * salp * (sig + c * ssig * (cos2 + c * csig * (2 * sqr(cos2) - 1)));

		if (i++ > 8192 || (oldvalue[1] == lambda && fabsl(lambda) > M_PI)) {
			d = 1;
			break;
		}
	} while (fabsl(oldvalue[0] - lambda) >= powl(2,-48));

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
			cos2 = csig - 2 * su1 * su2 / calp;

			d = (1 - c) * FLAT * (sig + c * ssig * (cos2 + c * csig * (2 * sqr(cos2) - 1)));
			salp = (londiff - asinl(lambda)) / d;
			calp = 1 - sqr(salp);

			lambda = salp * ssig / (cu1 * cu2);

			ssig = sqr(cu2 * lambda) + sqr(cu1 * su2 + su1 * cu2 * sqrtl(1 - sqr(lambda)));
			csig = -sqrtl(1 - ssig);
			ssig = sqrtl(ssig);
			sig = M_PI - asinl(ssig);

		} while (fabsl(oldvalue[0] - salp) >= powl(2,-48));
	}

	a = helmertA(calp);

	ssig = hypotl(cu2 * sinl(lambda), cu1 * su2 - su1 * cu2 * cosl(lambda));
	csig = su1 * su2 + cu1 * cu2 * cosl(lambda);
	sig = atan2_modified(ssig, csig);
	b = helmertB(calp, sig);

	s = RAD_MIN * a * (sig + b);

	if (i != 16384) {
		a = cu2 * sinl(lambda);
		b = cu1 * su2 - su1 * cu2 * cosl(lambda);

		c = cu1 * sinl(lambda);
		d = cu1 * su2 * cosl(lambda) - su1 * cu2;
	}
	else {
		a = salp / cu1;
		b = sqrtl(1 - sqr(a));
		if ((cu1 * su2 + su1 * cu2 * cosl(lambda)) < 0)
			b *= -1;

		c = salp;
		d = -su1 * ssig + cu1 * csig * b;
	}

	if (count > 2) {
		*res = s;
		*(res + 1) = normalise_a(atan2_modified(a, b));
		*(res + 2) = normalise_a(atan2_modified(c, d));
	}
	if (count > 3)
		*(res + 3) = salp;
	return;
}

void vincenty_direct(struct Coordinates *point, struct Vector *add, long double *res)
{
	long double u1, su1, cu1, s1, s2, salp, a;
	long double sigma, oldvalue;
	long double lambda, c;

	u1 = reduced_latitude(point->lat);
	su1 = sin_reduced_latitude(point->lat);
	cu1 = cos_reduced_latitude(point->lat);
	s1 = atan2_modified(tan_reduced_latitude(point->lat), cosl(add->theta));
	salp = cu1 * sinl(add->theta);

	a = helmertA(1 - sqr(salp));
	sigma = add->s / (RAD_MIN * a);
	sigma += helmertB1(1 - sqr(salp), sigma);
	s2 = 2 * s1 + sigma;

	*res = atan2_modified(su1 * cosl(sigma) + cu1 * sinl(sigma) * cosl(add->theta), (1 - FLAT) * hypotl(salp, su1 * sinl(sigma) - cu1 * cosl(sigma) * cosl(add->theta)));
	lambda = atan2_modified(sinl(sigma) * sinl(add->theta), cu1 * cosl(sigma) - su1 * sinl(sigma) * cosl(add->theta));

	c = Cfromcalp(1 - sqr(salp));
	*(res + 1) = lambda - (1 - c) * FLAT * salp * (sigma + c * sinl(sigma) * (cosl(s2) + c * cosl(sigma) * (2 * sqr(cosl(s2)) - 1)));
	*(res + 1) = normalise_c(*(res + 1) + point->lon);

	*(res + 2) = normalise_a(atan2_modified(salp, -su1 * sinl(sigma) + cu1 * cosl(sigma) * cosl(add->theta)));

	return;
}

