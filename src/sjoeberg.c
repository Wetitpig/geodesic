#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vincenty.h"
#include "sjoeberg.h"

struct Coordinates z(struct Coordinates *vertex, struct Coordinates *vertex2, int k, long double *c)
{
	struct Coordinates z;

	if (vertex->lat > vertex2->lat)
		*c = cosl(reduced_latitude(vertex->lat));
	else
		*c = cosl(reduced_latitude(vertex2->lat));

	z.lat = sqr(*c) * (1 - ECC) / sqr(cosl(vertex->lat));
	z.lon = sqr(*c) * (1 - ECC) / sqr(cosl(vertex2->lat));
	return z;
}

int combi(int n, int r)
{
	int nr;
	int fac[3];
	fac[0] = 1;
	fac[1] = 1;
	fac[2] = 1;

	if (n != r) {
		nr = n - r;
		do
			fac[2] = fac[2] * nr;
		while (--nr);
	}

	do
		fac[0] = fac[0] * n;
	while (--n);

	if (r != 0) {
		do
			fac[1] = fac[1] * r;
		while (--r);
	}

	return fac[0] / (fac[1] * fac[2]);
}

long double E(long double z, int k, long double *c)
{
	int j;
	long double hsum = 0, h, p, f, S, h0;
	long double feval[k + 1];

	p = ECC * sqr(*c);
	f = (1 - p) * p;
	p = 1 - 2 * p;
	S = f + p * z - sqr(z);
	if (S < 0)
		S = 0;
	h0 = sqr(*c) * (1 - ECC);

	for (j = 0; j <= k; j++) {
		switch(j)
		{
			case 0:
			feval[0] = (2 * z - p) / sqrtl(sqr(p) + 4 * f);
			if (feval[0] > 1)
				feval[0] = 1;
			if (feval[0] < -1)
				feval[0] = -1;
			feval[0] = asinl(feval[0]);
			break;

			case 1:
			feval[1] = logl((2 * f + p * z - sqrtl(f * S)) / z) / sqrtl(f);
			break;

			default:
			feval[j] = -sqrtl(S) / ((j - 1) * f * powl(z, (j - 1)));
			feval[j] = feval[j] - p * (2 * j - 3) / (2 * (j - 1) * f) * feval[j - 1];
			feval[j] = feval[j] + (j - 2) / ((j - 1) * f) * feval[j - 2];
			break;
		}

		h = powl(-h0, j) * combi(k, j) * feval[j];
		hsum = hsum + h;
	}

	return hsum / 2;
}

struct Vector sjoeberg(struct Coordinates *vertex, int i, int s, int a)
{
	struct vincenty_result inter[2];

	long double prev, next, excess = 0;
	long double area, darea = 0, interarea;
	struct Coordinates zres;
	long double z0, z1;
	long double *c = malloc(sizeof(long double));

	long double perimeter = 0;

	int h, k;
	for (h = 0; h < i; h++) {
		inter[0] = vincenty_inverse(vertex + h, vertex + ((h + 1) % i));

		if (s == 1)
			perimeter += inter[0].distance;

		if (a == 1) {
			next = inter[0].start;
			if (h == 0) {
				inter[1] = vincenty_inverse(vertex, vertex + i - 1);
				prev = inter[1].start;
			}
			else
				prev = normalise_a(inter[1].end - M_PI);

			excess += normalise_a(next - prev);

			for (k = 1; k < 6; k++) {
				interarea = powl(ECC, k) * (k + 1) / (2 * k + 1);
				zres = z(vertex + h, vertex + ((h + 1) % i), k, c);
				z0 = zres.lat;
				z1 = zres.lon;
				interarea = interarea * (E(z1, k, c) - E(z0, k, c));
				darea = interarea + darea;
			}

			inter[1] = inter[0];
		}
	}

	free(c);

	excess = excess - (i - 2) * M_PI;
	area = sqr(RAD_MIN) * (excess + darea);

	struct Vector res;
	res.s = perimeter;
	res.theta = area;
	return res;
}
