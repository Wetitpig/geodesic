#include <math.h>
#include "geodesic.h"
#include "vincenty.h"

const long double B1lookup[8][4] = {
{1.0l/2.0l,-9.0l/32.0l,205.0l/1536.0l,-4879.0l/73728.0l},
{5.0l/16.0l,-37.0l/96.0l,1335.0l/4096.0l,-86171.0l/363840.0l},
{29.0l/96.0l,-75.0l/128.0l,2901.0l/4096.0l},
{539.0l/1536.0l,-2391.0l/2560.0l,1082857.0l/737280.0l},
{3467.0l/7680.0l,-28223.0l/18432.0l},
{38081.0l/61440.0l,-733437.0l/286720.0l},
{459485.0l/516096.0l},
{109167851.0l/82575360.0l}
};

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

long double helmertA(long double calp)
{
	long double A = 0;
	long double usq, k1;
	int k;

	usq = calp * ECC2;
	k1 = sqrtl(1.0l + usq);
	k1 = (k1 - 1.0l) / (k1 + 1.0l);
	for (k = 0; k < 8; k++)
		A += sqr(double_fac(2 * k - 3) / double_fac(2 * k) * powl(k1, k));
	A /= (1 - k1);
	return A;
}

long double helmertB(long double calp, long double sig, long double A)
{
	long double B = 0;
	long double usq, k1, c;
	int k, j;

	usq = calp * ECC2;
	k1 = sqrtl(1.0l + usq);
	k1 = (k1 - 1.0l) / (k1 + 1.0l);
	A *= (1 - k1);
	for (k = 1; k < 17; k++) {
		c = double_fac(2 * k - 3) / double_fac(2 * k) * powl(k1, k);
		for (j = 1; (2 * j + k) < 17; j++)
			c -= double_fac(2 * (j + k) - 3) / double_fac(2 * (j + k)) * double_fac(2 * j - 3) / double_fac(2 * j) * powl(k1, 2 * j + k);
		B += c / (A * k) * sinl(2.0l * k * sig);
	}
	return B;
}

long double helmertB1(long double calp, long double sig)
{
	long double B = 0;
	long double usq, k1, c;
	int k, j, i;

	usq = calp * ECC2;
	k1 = sqrtl(1.0l + usq);
	k1 = (k1 - 1.0l) / (k1 + 1.0l);
	for (k = 1; k < 9; k++) {
		c = 0;
		i = 0;
		for (j = k; j < 9; j += 2)
			c += B1lookup[k - 1][i++] * powl(k1, j);
		B += c * sinl(2.0l * k * sig);
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

	a = helmertA(calp);
	b = helmertB(calp, sig, a);

	s = RAD_MIN * a * (sig - b);

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
	long double u1, su1, cu1, s1, s2, salp, a;
	long double sigma, oldvalue;
	long double lambda, c;

	u1 = reduced_latitude(point->lat);
	su1 = sin_reduced_latitude(point->lat);
	cu1 = cos_reduced_latitude(point->lat);
	s1 = atan2_modified(tan_reduced_latitude(point->lat), cosl(add->theta));
	salp = cu1 * sinl(add->theta);

	a = helmertA(1.0l - sqr(salp));
	sigma = add->s / (RAD_MIN * a);
	sigma += helmertB1(1 - sqr(salp), sigma);
	s2 = 2.0l * s1 + sigma;

	*res = atan2_modified(su1 * cosl(sigma) + cu1 * sinl(sigma) * cosl(add->theta), (1 - FLAT) * hypotl(salp, su1 * sinl(sigma) - cu1 * cosl(sigma) * cosl(add->theta)));
	lambda = atan2_modified(sinl(sigma) * sinl(add->theta), cu1 * cosl(sigma) - su1 * sinl(sigma) * cosl(add->theta));

	c = Cfromcalp(1.0l - sqr(salp));
	*(res + 1) = lambda - (1 - c) * FLAT * salp * (sigma + c * sinl(sigma) * (cosl(s2) + c * cosl(sigma) * (2.0l * sqr(cosl(s2)) - 1.0l)));
	*(res + 1) = normalise_c(*(res + 1) + point->lon);

	*(res + 2) = normalise_a(atan2_modified(salp, -su1 * sinl(sigma) + cu1 * cosl(sigma) * cosl(add->theta)));

	return;
}

