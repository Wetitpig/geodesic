#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include "math.h"
#include "vincenty.h"
#include "karney.h"

const long double Clookup[6][6] = {
{
2.0l/3.0l - ECC2/15.0l + 4.0l*ECC2*ECC2/105.0l - 8.0l*ECC2*ECC2*ECC2/315.0l + 64.0l*ECC2*ECC2*ECC2*ECC2/3465.0l - 128.0l*ECC2*ECC2*ECC2*ECC2*ECC2/9009.0l,
(- 1.0l/20.0l + ECC2/35.0l - 2.0l*ECC2*ECC2/105.0l + 16.0l*ECC2*ECC2*ECC2/1155.0l - 32.0l*ECC2*ECC2*ECC2*ECC2/3003.0l)*ECC2,
(1.0l/42.0l - ECC2/63.0l + 8.0l*ECC2*ECC2/693.0l - 80.0l*ECC2*ECC2*ECC2/9009.0l)*ECC2*ECC2,
(- 1.0l/72.0l + ECC2/99.0l - 10.0l*ECC2*ECC2/1287.0l)*ECC2*ECC2*ECC2,
(1.0l/110.0l - ECC2/143.0l)*ECC2*ECC2*ECC2*ECC2,
- ECC2*ECC2*ECC2*ECC2*ECC2/156.0l
},
{
0,
(1.0l/180.0l - ECC2/315.0l + 2.0l*ECC2*ECC2/945.0l - 16.0l*ECC2*ECC2*ECC2/10395.0l + 32.0l*ECC2*ECC2*ECC2*ECC2/27027.0l)*ECC2,
(- 1.0l/252.0l + ECC2/378.0l - 4.0l*ECC2*ECC2/2079.0l + 40.0l*ECC2*ECC2*ECC2/27027.0l)*ECC2*ECC2,
(1.0l/360.0l - ECC2/495.0l + 2.0l*ECC2*ECC2/1287.0l)*ECC2*ECC2*ECC2,
(- 1.0l/495.0l + 2.0l*ECC2/1287.0l)*ECC2*ECC2*ECC2*ECC2,
5.0l*ECC2*ECC2*ECC2*ECC2*ECC2/3276.0l
},
{
0,
0,
(1.0/2100.0l - ECC2/3150.0l + 4.0l*ECC2*ECC2/17325.0l - 8.0l*ECC2*ECC2*ECC2/45045.0l)*ECC2*ECC2,
(- 1.0l/1800.0l + ECC2/2475.0l - 2.0l*ECC2*ECC2/6435.0l)*ECC2*ECC2*ECC2,
(1.0l/1925.0l - 2.0l*ECC2/5005.0l)*ECC2*ECC2*ECC2*ECC2,
- ECC2*ECC2*ECC2*ECC2*ECC2/2184.0l
},
{
0,
0,
0,
(1.0l/17640.0l - ECC2/24255.0l + 2.0l*ECC2*ECC2/63063.0l)*ECC2*ECC2*ECC2,
(- 1.0l/10780.0l + ECC2/14014.0l)*ECC2*ECC2*ECC2*ECC2,
5.0l*ECC2*ECC2*ECC2*ECC2*ECC2/45864.0l
},
{
0,
0,
0,
0,
(1.0l/124740.0l - ECC2/162162.0l)*ECC2*ECC2*ECC2*ECC2,
- ECC2*ECC2*ECC2*ECC2*ECC2/58968.0l
},
{
0,
0,
0,
0,
0,
ECC2*ECC2*ECC2*ECC2*ECC2/792792.0l
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
		csigp[2 * k + 2] = 2.0l * csig * csigp[2 * k + 1] - csigp[2 * k];
		csigp[2 * k + 3] = 2.0l * csig * csigp[2 * k + 2] - csigp[2 * k + 1];
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

	long double londiff = fabsl(normalise_c(lon[1] - lon[0]));
	long double latdiff = sinl(lat[1]) / (1.0l - ECC * sqr(sinl(lat[1]))) - sinl(lat[0]) / (1.0l - ECC * sqr(sinl(lat[0])));
	latdiff += logl((1.0l + sqrtl(ECC) * sinl(lat[1])) * (1.0l - sqrtl(ECC) * sinl(lat[0])) / (1.0l - sqrtl(ECC) * sinl(lat[1])) / (1.0l + sqrtl(ECC) * sinl(lat[0]))) / (2.0l * sqrtl(ECC));

	return sqr(RAD_MIN) * londiff * latdiff / 2.0l;
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
			excess += normalise_a(next - prev + M_PI_L);

			if (coor0->lon != coor1->lon && !isnan(tanl(coor0->lat)) && !isnan(tanl(coor1->lat))) {
				interarea = sinl(2.0l * asinl(*(inter + 3))) / 2.0l;

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
		excess -= (i - 2) * M_PI_L;
		area = (sqr(RAD_MAJ) + sqr(RAD_MIN) * atanhl(sqrtl(ECC)) / sqrtl(ECC)) * excess / 2.0l;
		area += darea * ECC * sqr(RAD_MAJ);
	}

	*res = perimeter;
	*(res + 1) = area;
	return;
}
