#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include "vincenty.h"
#include "karney.h"

const long double Clookup[6][6] = {
{
2.0/3 - ECC2/15 + 4*ECC2*ECC2/105 - 8*ECC2*ECC2*ECC2/315 + 64*ECC2*ECC2*ECC2*ECC2/3465 - 128*ECC2*ECC2*ECC2*ECC2*ECC2/9009,
(- 1.0/20 + ECC2/35 - 2*ECC2*ECC2/105 + 16*ECC2*ECC2*ECC2/1155 - 32*ECC2*ECC2*ECC2*ECC2/3003)*ECC2,
(1.0/42 - ECC2/63 + 8*ECC2*ECC2/693 - 80*ECC2*ECC2*ECC2/9009)*ECC2*ECC2,
(- 1.0/72 + ECC2/99 - 10*ECC2*ECC2/1287)*ECC2*ECC2*ECC2,
(1.0/110 - ECC2/143)*ECC2*ECC2*ECC2*ECC2,
- ECC2*ECC2*ECC2*ECC2*ECC2/156
},
{
0,
(1.0/180 - ECC2/315 + 2*ECC2*ECC2/945 - 16*ECC2*ECC2*ECC2/10395 + 32*ECC2*ECC2*ECC2*ECC2/27027)*ECC2,
(- 1.0/252 + ECC2/378 - 4*ECC2*ECC2/2079 + 40*ECC2*ECC2*ECC2/27027)*ECC2*ECC2,
(1.0/360 - ECC2/495 + 2*ECC2*ECC2/1287)*ECC2*ECC2*ECC2,
(- 1.0/495 + 2*ECC2/1287)*ECC2*ECC2*ECC2*ECC2,
5*ECC2*ECC2*ECC2*ECC2*ECC2/3276
},
{
0,
0,
(1.0/2100 - ECC2/3150 + 4*ECC2*ECC2/17325 - 8*ECC2*ECC2*ECC2/45045)*ECC2*ECC2,
(- 1.0/1800 + ECC2/2475 - 2*ECC2*ECC2/6435)*ECC2*ECC2*ECC2,
(1.0/1925 - 2*ECC2/5005)*ECC2*ECC2*ECC2*ECC2,
- ECC2*ECC2*ECC2*ECC2*ECC2/2184
},
{
0,
0,
0,
(1.0/17640 - ECC2/24255 + 2*ECC2*ECC2/63063)*ECC2*ECC2*ECC2,
(- 1.0/10780 + ECC2/14014)*ECC2*ECC2*ECC2*ECC2,
5*ECC2*ECC2*ECC2*ECC2*ECC2/45864
},
{
0,
0,
0,
0,
(1.0/124740 - ECC2/162162)*ECC2*ECC2*ECC2*ECC2,
- ECC2*ECC2*ECC2*ECC2*ECC2/58968
},
{
0,
0,
0,
0,
0,
ECC2*ECC2*ECC2*ECC2*ECC2/792792
}
};

long double sigeval(long double lat, long double azimuth)
{
	long double y, x;
	y = sin_reduced_latitude(lat);
	x = cos_reduced_latitude(lat) * cosl(azimuth);
	return atan2_modified(y, x);
}

long double I4(long double csig, long double salp)
{
	int k, j;
	long double csigp[12];
	csigp[0] = 1;
	csigp[1] = csig;

	long double i = 0, c;
	for (k = 0; k < 6; k++) {
		c = 0;
		for (j = k; j < 6; j++)
			c += powl(1 - sqr(salp), j) * Clookup[k][j];
		csigp[2 * k + 2] = 2 * csig * csigp[2 * k + 1] - csigp[2 * k];
		csigp[2 * k + 3] = 2 * csig * csigp[2 * k + 2] - csigp[2 * k + 1];
		i += c * csigp[2 * k + 1];
	}
	return i;
}

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

long double ellipblock(struct Coordinates *vertex)
{
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

	long double londiff = lon[1] - lon[0];
	long double latdiff = sinl(lat[1]) / (1 - ECC * sqr(sinl(lat[1]))) - sinl(lat[0]) / (1 - ECC * sqr(sinl(lat[0])));
	latdiff += logl((1 + sqrtl(ECC) * sinl(lat[1])) * (1 - sqrtl(ECC) * sinl(lat[0])) / (1 - sqrtl(ECC) * sinl(lat[1])) / (1 + sqrtl(ECC) * sinl(lat[0]))) / (2 * sqrtl(ECC));

	return sqr(RAD_MIN) * londiff * latdiff / 2;
}

void karney(struct Coordinates *vertex, int i, int s, int a, long double *res)
{
	long double *inter = malloc(sizeof(long double) * 8);
	long double prev, next, excess = 0;
	long double area, darea = 0, interarea;
	long double sig0, sig1;

	long double perimeter = 0;

	int h, k;

	if (a == 1) {
		if (i == 4 && isblock(vertex)) {
			area = ellipblock(vertex);
			a = 0;
		}
		else if (i == 3 && (k = ispolariso(vertex)) < 3) {
			struct Coordinates *triangle = malloc(sizeof(struct Coordinates) * 4);

			triangle->lat = (vertex + k)->lat;
			triangle->lon = (vertex + ((k + 2) % 3))->lon;
			memcpy(triangle + 1, vertex + ((k + 1) % 3), sizeof(struct Coordinates));
			for (h = 2; h < 4; h++) {
				(triangle + h)->lat = (vertex + ((k + h) % 3))->lat;
				(triangle + h)->lon = (triangle + (3 - h))->lon;
			}

			area = ellipblock(triangle);
			free(triangle);
			a = 0;
		}
	}

	struct Coordinates *coor0, *coor1;
	for (h = 0; h < i; h++) {
		coor0 = vertex + h;
		coor1 = vertex + h + 1;
		vincenty_inverse(coor0, coor1, inter, 4);

		if (s == 1)
			perimeter += *inter;

		if (a == 1) {
			next = *(inter + 1);
			if (h == 0)
				vincenty_inverse(vertex + i - 1, vertex, inter + 4, 4);
			prev = *(inter + 6);
			excess += normalise_a(next - prev + M_PI);

			if (coor0->lon != coor1->lon && !isnan(tanl(coor0->lat)) && !isnan(tanl(coor1->lat))) {
				interarea = sinl(2 * asinl(*(inter + 3))) / 2;

				sig1 = sigeval(coor1->lat, *(inter + 2));
				sig0 = sigeval(coor0->lat, *(inter + 1));
				interarea *= I4(cosl(sig0), *(inter + 3)) - I4(cosl(sig1), *(inter + 3));

				darea += interarea;
			}
			memcpy(inter + 4, inter, sizeof(long double) * 4);
		}
	}

	free(inter);
	if (a == 1) {
		excess -= (i - 2) * M_PI;
		area = (sqr(RAD_MAJ) + sqr(RAD_MIN) * atanhl(sqrtl(ECC)) / sqrtl(ECC)) * excess / 2;
		area += darea * ECC * sqr(RAD_MAJ);
	}

	*res = perimeter;
	*(res + 1) = area;
	return;
}
