#include <math.h>
#include "geodesic.h"
#include "vincenty.h"

long double tan_reduced_latitude(long double lat)
{
	return (1.0l - FLAT) * tanl(lat);
}

long double reduced_latitude(long double lat)
{
	if (isnan(tanl(lat)))
		return copysignl(M_PI_L / 2.0l, lat);
	else
		return atanl(tan_reduced_latitude(lat));
}

long double sin_reduced_latitude(long double lat)
{
	if (isnan(tanl(lat)))
		return copysignl(1.0l, lat);
	else
		return copysignl(1.0l / hypotl(1.0l, 1.0l / tan_reduced_latitude(lat)), lat);
}

long double cos_reduced_latitude(long double lat)
{
	if (isnan(tanl(lat)))
		return 0.0l;
	else
		return 1.0l / hypotl(1.0l, tan_reduced_latitude(lat));
}

long double helmertA(long double k1)
{
	long double A = 0;
	int k;

	for (k = 0; k < 5; k++)
		A += sqr(double_fac[k] * powl(k1, k));
	return A;
}

long double helmertB(long double k1, long double sig)
{
	long double B = 0;
	long double c, c1 = 0;
	int k, j;

	for (k = 1; k < 9; k++) {
		c = double_fac[k];
		for (j = 1; (2 * j + k) < 9 && c1 != c; j++) {
			c1 = c;
			c -= double_fac[j + k] * double_fac[j] * powl(k1, 2 * j);
		}
		B += c / k * sinl(2.0l * k * sig) * powl(k1, k);
	}
	return B;
}

static inline long double Cfromcalp(long double calp)
{
	return FLAT / 16.0l * calp * (4.0l + FLAT * (4.0l - 3.0l * calp));
}

void vincenty_inverse(struct Coordinates *location, struct Coordinates *location2, long double *res, int count)
{
	long double londiff, lambda, u1, u2, cu1, cu2, su1, su2;
	long double ssig, csig, sig, salp, calp, cos2, c;
	long double a, b, k1, s;
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
		cos2 = csig - 2.0l * su1 * su2 / calp;
		if (cos2 < -1 || isnan(cos2))
			cos2 = -1.0l;

		c = Cfromcalp(calp);

		lambda = londiff + (1.0l - c) * FLAT * salp * (sig + c * ssig * (cos2 + c * csig * (2.0l * sqr(cos2) - 1.0l)));

		if (i++ > 8192 || (oldvalue[1] == lambda && fabsl(lambda) > M_PI_L)) {
			d = 1;
			break;
		}
	} while (fabsl(oldvalue[0] - lambda) >= powl(2,-48));

	if (d == 1) {
		londiff = (londiff > 0 ? M_PI_L : -M_PI_L) - londiff;
		lambda = 0;

		i = 16384;
		calp = 0.5l;
		cos2 = 0;

		sig = M_PI_L - fabsl(u1 + u2);
		ssig = sinl(sig);
		csig = cosl(sig);

		do {
			oldvalue[0] = salp;
			c = Cfromcalp(calp);
			cos2 = csig - 2.0l * su1 * su2 / calp;

			d = (1.0l - c) * FLAT * (sig + c * ssig * (cos2 + c * csig * (2.0l * sqr(cos2) - 1.0l)));
			salp = (londiff - asinl(lambda)) / d;
			calp = 1.0l - sqr(salp);

			lambda = salp * ssig / (cu1 * cu2);

			ssig = sqr(cu2 * lambda) + sqr(cu1 * su2 + su1 * cu2 * sqrtl(1 - sqr(lambda)));
			csig = -sqrtl(1.0l - ssig);
			ssig = sqrtl(ssig);
			sig = M_PI_L - asinl(ssig);

		} while (fabsl(oldvalue[0] - salp) >= powl(2,-48));
	}

	k1 = calp * ECC2;
	k1 = sqrtl(1.0l + k1);
	k1 = (k1 - 1.0l) / (k1 + 1.0l);

	s = RAD_MIN * (helmertA(k1) * sig - helmertB(k1, sig)) / (1.0l - k1);

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
			b *= -1.0l;

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
	long double u1, su1, cu1, s1, s2, salp, a, b, k1, k2, usq;
	long double sigma, sig0, oldvalue;
	long double lambda, c;
	int i;

	u1 = reduced_latitude(point->lat);
	su1 = sin_reduced_latitude(point->lat);
	cu1 = cos_reduced_latitude(point->lat);
	s1 = atan2_modified(tan_reduced_latitude(point->lat), cosl(add->theta));
	salp = cu1 * sinl(add->theta);

	usq = (1.0l - sqr(salp)) * ECC2;
	k1 = sqrtl(1.0l + usq);
	k1 = (k1 - 1.0l) / (k1 + 1.0l);
	k2 = 1.0l - k1;

	a = helmertA(k1);
	sig0 = add->s / RAD_MIN;
	sigma = k2 * sig0 / a;
	do {
		oldvalue = sigma;
		sigma += (sig0 - (a * sigma - helmertB(k1, sigma)) / k2) / sqrt(1.0l + usq * sqr(sinl(sigma)));
	} while (fabsl(oldvalue - sigma) >= powl(2, -48));
	s2 = 2.0l * s1 + sigma;

	*res = atan2_modified(su1 * cosl(sigma) + cu1 * sinl(sigma) * cosl(add->theta), (1 - FLAT) * hypotl(salp, su1 * sinl(sigma) - cu1 * cosl(sigma) * cosl(add->theta)));
	lambda = atan2_modified(sinl(sigma) * sinl(add->theta), cu1 * cosl(sigma) - su1 * sinl(sigma) * cosl(add->theta));

	c = Cfromcalp(1.0l - sqr(salp));
	*(res + 1) = lambda - (1 - c) * FLAT * salp * (sigma + c * sinl(sigma) * (cosl(s2) + c * cosl(sigma) * (2.0l * sqr(cosl(s2)) - 1.0l)));
	*(res + 1) = normalise_c(*(res + 1) + point->lon);

	*(res + 2) = normalise_a(atan2_modified(salp, -su1 * sinl(sigma) + cu1 * cosl(sigma) * cosl(add->theta)));

	return;
}

