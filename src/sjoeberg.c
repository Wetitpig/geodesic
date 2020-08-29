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

long double combi(int n, int r)
{
	int nr = n - r;
	long fac[3];
	fac[0] = 1;
	fac[1] = 1;
	fac[2] = 1;

	do
		fac[0] = fac[0] * n;
	while (--n);

	do
		fac[1] = fac[1] * r;
	while (--r);

	do
		fac[2] = fac[2] * nr;
	while (--nr);

	return fac[0] / (fac[1] * fac[2]);
}

long double E(long double z, int k, long double *c)
{
	int j;
	long double hsum = 0, h, p, f, S, h0;
	long double feval[k + 1];

	p = sqr(ECC * *c);
	f = (1 - p) * p;
	p = 1 - 2 * p;
	S = f + p * z - sqr(z);
	h0 = sqr(*c) * (1 - ECC);

	for (j = 0; j <= k; j++) {
		switch(j)
		{
			case 0:
			feval[0] = asinl((2 * z - p) / sqrtl(sqr(p) + 4 * f));
			break;

			case 1:
			feval[1] = (logl(2 * f + p * z - sqrtl(f * S)) - logl(z)) / sqrtl(f);
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

long double sjoeberg_area(struct Coordinates *vertex, int i)
{
	long double prev, next, excess = 0;
	long double area, darea = 0, interarea;
	struct vincenty_result inter[2];
	struct Coordinates zres;
	long double z0, z1;
	long double *c = malloc(sizeof(long double));

	int h, k;
	for (h = 0; h < i; h++) {
		inter[0] = vincenty_inverse(vertex + h, vertex + ((h + 1) % i));
		next = inter[0].start;
		if (h == 0) {
			inter[1] = vincenty_inverse(vertex, vertex + i - 1);
			prev = inter[1].start;
		}
		else
			prev = normalise_a(inter[1].end - M_PI);

		excess += normalise_a(next - prev);

		for (k = 1; k < 11; k++) {
			interarea = powl(ECC, k) * (k + 1) / (2 * k + 1);
			zres = z(vertex + h, vertex + ((h + 1) % i), k, c);
			z0 = zres.lat;
			z1 = zres.lon;
			interarea = interarea * (E(z1, k, c) - E(z0, k, c));
			darea = interarea + darea;
		}

		inter[1] = inter[0];
	}

	free(c);

	excess = fabsl(excess) - (i - 2) * M_PI;
	return sqr(RAD_MIN) * (excess + darea);
}

long double sjoeberg_perimeter(struct Coordinates *vertex, int i)
{
	struct vincenty_result inter;
	int k;
	long double perimeter = 0;
	for (k = 0; k < i; k++) {
		inter = vincenty_inverse(vertex + k, vertex + ((k + 1) % i));
		perimeter += inter.distance;
	}
	return perimeter;
}
