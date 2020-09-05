#include <math.h>
#include "math.h"
#include "io.h"
#include "mpblock.h"

int isblock(struct Coordinates *vertex)
{
	int block = 1;
	block &= (vertex->lat == (vertex + 1)->lat && (vertex + 2)->lat == (vertex + 3)->lat && (vertex + 1)->lon == (vertex + 2)->lon && vertex->lon == (vertex + 3)->lon);
	block |= (vertex->lon == (vertex + 1)->lon && (vertex + 2)->lon == (vertex + 3)->lon && (vertex + 1)->lat == (vertex + 2)->lat && vertex->lat == (vertex + 3)->lat);
	return block;
}

int ispolariso(struct Coordinates *vertex)
{
	int k;
	for (k = 0; k < 3; k++) {
		if (fabsl((vertex + k)->lat) / RAD == 90 && (vertex + (k + 1) % 3)->lat == (vertex + (k + 2) % 3)->lat)
			break;
	}
	return k;
}

long double ellipblock(long double lat0, long double lat1, long double lon0, long double lon1)
{
	long double londiff = fabsl(normalise_c(lon1 - lon0));
	long double latdiff = sinl(lat1) / (1.0l - ECC * sqr(sinl(lat1))) - sinl(lat0) / (1.0l - ECC * sqr(sinl(lat0)));
	latdiff += logl((1.0l + sqrtl(ECC) * sinl(lat1)) * (1.0l - sqrtl(ECC) * sinl(lat0)) / (1.0l - sqrtl(ECC) * sinl(lat1)) / (1.0l + sqrtl(ECC) * sinl(lat0))) / (2.0l * sqrtl(ECC));

	return sqr(RAD_MIN) * londiff * latdiff / 2.0l;
}

long double parallel_length(long double lon0, long double lon1, long double lat)
{
	return cosl(lat) * fabsl(normalise_c(lon1 - lon0)) / sqrtl(1.0l - ECC * sqr(sinl(lat))) * RAD_MAJ;
}

long double meridian_arc(long double lat0, long double lat1)
{
	int k, j;
	long double c = 0, m;

	for (j = 0; j < 11; j++)
		c += sqr(double_fac(2 * j - 3) / double_fac(2 * j) * powl(FLAT_3, j));
	m = c * (lat1 - lat0);

	for (k = 1; k < 9; k++) {
		c = 0;
		for (j = 0; j < 8; j++)
			c += double_fac(2 * j - 3) / double_fac(2 * j) * double_fac(2 * (j + k) - 3) / double_fac(2 * (j + k)) * powl(FLAT_3, k + 2 * j);
		c *= powl(-1.0l, k) * (1 - sqr(2 * k)) / k;
		m += c * (sin(2.0l * k * lat1) - sin(2.0l * k * lat0));
	}

	return (RAD_MAJ + RAD_MIN) / 2 * m;
}

void mpblock_area(struct Coordinates *vertex, int i, long double s, long double a, long double *res)
{
	int h, k;

	if (i == 4 && isblock(vertex)) {
		long double lat[2], lon[2];
		if (vertex->lat < (vertex + 2)->lat) {
			lat[0] = vertex->lat;
			lat[1] = (vertex + 2)->lat;
		}
		else {
			lat[0] = (vertex + 2)->lat;
			lat[1] = vertex->lat;
		}

		if (vertex->lon < (vertex + 2)->lon) {
			lon[0] = vertex->lon;
			lon[1] = (vertex + 2)->lon;
		}
		else {
			lon[0] = (vertex + 2)->lon;
			lon[1] = vertex->lon;
		}

		if (a == 1)
			*(res + 1) = ellipblock(lat[0], lat[1], lon[0], lon[1]);
		if (s == 1) {
			*res = parallel_length(lon[0], lon[1], lat[0]);
			*res += parallel_length(lon[0], lon[1], lat[1]);
			*res += meridian_arc(lat[0], lat[1]) * 2.0l;
		}
	}
	else if (i == 3 && (k = ispolariso(vertex)) < 3) {
		long double lon[2];

		if ((vertex + k + 1)->lon < (vertex + ((k + 2) % 3))->lon) {
			lon[0] = (vertex + k + 1)->lon;
			lon[1] = (vertex + ((k + 2) % 3))->lon;
		}
		else {
			lon[1] = (vertex + k + 1)->lon;
			lon[0] = (vertex + ((k + 2) % 3))->lon;
		}

		if (a == 1)
			*(res + 1) = ellipblock((vertex + k + 1)->lat, (vertex + k)->lat, lon[0], lon[1]);
		if (s == 1) {
			*res = parallel_length(lon[0], lon[1], (vertex + k + 1)->lat);
			*res += meridian_arc((vertex + k + 1)->lat, (vertex + k)->lat) * 2.0l;
		}
	}
	else
		error("Not a meridian-parallel block.");
	return;
}

void mpblock_inverse(struct Coordinates *location, struct Coordinates *location2, long double *res)
{
	if (location->lon == location2->lon || location->lat / RAD == 90 || location2->lat / RAD == 90) {
		long double lat[2];
		if (location->lat < location2->lat) {
			lat[0] = location->lat;
			lat[1] = location2->lat;
			*(res + 1) = 0;
		}
		else {
			lat[1] = location->lat;
			lat[0] = location2->lat;
			*(res + 1) = M_PI_L;
		}
		*res = meridian_arc(lat[0], lat[1]);
	}
	else if (location->lat == location2->lat) {
		long double conjugate = normalise_c(location->lon + M_PI_L);

		if (location->lon > 0) {
			if (location2->lon < location->lon && location2->lon > conjugate)
				*(res + 1) = 3 * M_PI_L / 2.0l;
			else if (location2->lon > location->lon || location2->lon < conjugate)
				*(res + 1) = M_PI_L / 2.0l;
			else
				*(res + 1) = M_PI_L * (location->lat < 0);
		}
		else if (location->lon < 0) {
			if (location2->lon > location->lon && location2->lon < conjugate)
				*(res + 1) = M_PI_L / 2.0l;
			else if (location2->lon < location->lon || location2->lon > conjugate)
				*(res + 1) = 3 * M_PI_L / 2.0l;
			else
				*(res + 1) = M_PI_L * (location->lat < 0);
		}
		else if (location2->lon / RAD != 180)
			*(res + 1) = (location->lon < 0) * M_PI_L + M_PI_L / 2;
		else
			*(res + 1) = M_PI_L * (location->lat < 0);

		*res = parallel_length(location->lon, location2->lon, location->lat);
	}
	else
		error("Not a meridian arc or a parallel arc");
	return;
}
