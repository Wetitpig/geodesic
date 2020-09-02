#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include "vincenty.h"
#include "karney.h"

const long double Clookup[6][6] =
{
{
2.0/3 - ECC2/15 + 4*ECC2*ECC2/105 - 8*ECC2*ECC2*ECC2/315 + 64*ECC2*ECC2*ECC2*ECC2/3465 - 128*ECC2*ECC2*ECC2*ECC2*ECC2/9009,
- 1.0/20 + ECC2/35 - 2*ECC2*ECC2/105 + 16*ECC2*ECC2*ECC2/1155 - 32*ECC2*ECC2*ECC2*ECC2/3003,
1.0/42 - ECC2/63 + 8*ECC2*ECC2/693 - 80*ECC2*ECC2*ECC2/9009,
- 1.0/72 + ECC2/99 - 10*ECC2*ECC2/1287,
1.0/110 - ECC2/143,
- 1.0/156
},
{
0,
1.0/180 - ECC2/315 + 2*ECC2*ECC2/945 - 16*ECC2*ECC2*ECC2/10395 + 32*ECC2*ECC2*ECC2*ECC2/27027,
- 1.0/252 + ECC2/378 - 4*ECC2*ECC2/2079 + 40*ECC2*ECC2*ECC2/27027,
1.0/360 - ECC2/495 + 2*ECC2*ECC2/1287,
- 1.0/495 + 2*ECC2/1287,
5.0/3276
},
{
0,
0,
1.0/2100 - 1*ECC2/3150 + 4*ECC2*ECC2/17325 - 8*ECC2*ECC2*ECC2/45045,
- 1.0/1800 + ECC2/2475 - 2*ECC2*ECC2/6435,
1.0/1925 - 2*ECC2/5005,
- 1.0/2184
},
{
0,
0,
0,
1.0/17640 - ECC2/24255 + 2*ECC2*ECC2/63063,
- 1.0/10780 + ECC2/14014,
5.0/45864
},
{
0,
0,
0,
0,
1.0/124740 - ECC2/162162,
- 1.0/58968
},
{
0,
0,
0,
0,
0,
1.0/792792
}
};

long double I4(long double sig, long double salp)
{
	long double i = 0, c;
	int k, j;
	for (k = 0; k < 6; k++) {
		c = 0;	
		for (j = k; j < 6; j++)
			c += powl(ECC2 * (1 - sqr(salp)), j) * Clookup[k][j];
		i += c * cosl((2 * k + 1) * sig);
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
	long double slat = sinl(lat[1]);
	long double latdiff = slat / (1 - ECC * sqr(slat)) + logl((1 + sqrtl(ECC) * slat) / (1 - sqrtl(ECC) * slat)) / (2 * sqrtl(ECC));
	slat = sinl(lat[0]);
	latdiff -= slat / (1 - ECC * sqr(slat)) + logl((1 + sqrtl(ECC) * slat) / (1 - sqrtl(ECC) * slat)) / (2 * sqrtl(ECC));

	return sqr(RAD_MIN) * londiff * latdiff / 2;
}

void karney(struct Coordinates *vertex, int i, int s, int a, long double *res)
{
	long double *inter = malloc(sizeof(long double) * 8);
	long double prev, next, excess = 0;
	long double area, darea = 0, interarea;
	long double csig0, csig1;

	long double perimeter = 0;

	int h, k;

/*	if (a == 1) {
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
*/
	for (h = 0; h < i; h++) {
		vincenty_inverse(vertex + h, vertex + h + 1, inter, 4);

		if (s == 1)
			perimeter += *inter;

		if (a == 1) {
			next = *(inter + 1);
			if (h == 0) {
				vincenty_inverse(vertex, vertex + i - 1, inter + 4, 4);
				prev = *(inter + 5);
			}
			else
				prev = normalise_a(*(inter + 6) - M_PI);
			excess += normalise_a(next - prev);

			interarea = *(inter + 3) * sqrtl(1 - sqr(*(inter + 3)));
//			if ((vertex + h)->lat < (vertex + h + 1)->lat) {
//				csig0 = 1 / (hypotl(1, tan_reduced_latitude((vertex + h + 1)->lat)) * hypotl(1, sqrtl(ECC2) * cosl((vertex + h + 1)->lat)));
//				csig1 = 1 / (hypotl(1, tan_reduced_latitude((vertex + h)->lat)) * hypotl(1, sqrtl(ECC2) * cosl((vertex + h)->lat)));
//			}
//			else {
				csig1 = 1 / (hypotl(1, tan_reduced_latitude((vertex + h + 1)->lat)) * hypotl(1, sqrtl(ECC2) * cosl((vertex + h + 1)->lat)));
				csig0 = 1 / (hypotl(1, tan_reduced_latitude((vertex + h)->lat)) * hypotl(1, sqrtl(ECC2) * cosl((vertex + h)->lat)));
//			}
			interarea *= I4(acosl(csig1), *(inter + 3)) - I4(acosl(csig0), *(inter + 3));

			darea += interarea;
			memcpy(inter + 4, inter, sizeof(long double) * 4);
		}
	}

	free(inter);
	if (a == 1) {
		excess -= (i - 2) * M_PI;
		area = (sqr(RAD_MAJ) / 2 + sqr(RAD_MIN) * atanhl(sqrtl(ECC)) / (2 * sqrtl(ECC))) * excess;
		area += darea * ECC * sqr(RAD_MAJ);
	}

	*res = perimeter;
	*(res + 1) = area;
	return;
}
