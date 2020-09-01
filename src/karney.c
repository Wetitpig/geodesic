#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include "vincenty.h"
#include "karney.h"
#include "karney_lookup.h"

#define E2C2 ECC * sqr(c)

long double I4(long double sig, long double salp)
{
	long double i = 0, c = 0;
	int k, j;
	for (k = 0; k < 6; k++) {
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
			for (h = 1; h < 4; h++)
				(triangle + h)->lat = (vertex + ((k + h) % 3))->lat;
			triangle->lon = (vertex + ((k + 2) % 3))->lon;
			(triangle + 1)->lon = (vertex + ((k + 1) % 3))->lon;
			(triangle + 2)->lon = (triangle + 1)->lon;
			(triangle + 3)->lon = triangle->lon;

			area = ellipblock(triangle);
			free(triangle);
			a = 0;
		}
	}

	for (h = 0; h < i; h++) {
		vincenty_inverse(vertex + h, vertex + ((h + 1) % i), inter, 4);

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
			if ((vertex + h)->lon > (vertex + ((h + 1) % i))->lon)
				interarea *= I4(atan2_modified(tan_reduced_latitude((vertex + h)->lat), cos(*(inter + 1))), *(inter + 3)) - I4(atan2_modified(tan_reduced_latitude((vertex + ((h + 1) % i))->lat), cos(*(inter + 2))), *(inter + 3));
			else
				interarea *= I4(atan2_modified(tan_reduced_latitude((vertex + ((h + 1) % i))->lat), cos(*(inter + 2))), *(inter + 3)) - I4(atan2_modified(tan_reduced_latitude((vertex + h)->lat), cos(*(inter + 1))), *(inter + 3));

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
