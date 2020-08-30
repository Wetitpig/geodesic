#include <stdlib.h>
#include <math.h>
#include "vincenty.h"
#include "sjoeberg.h"

int binom(int n, int r)
{
	int nr;
	int fac[2];
	fac[0] = 1;
	fac[1] = 1;

	nr = n - r;
	while (nr)
		fac[1] = fac[1] * nr--;

	while (n > r)
		fac[0] = fac[0] * n--;

	return fac[0] / fac[1];
}

long double E(struct Coordinates *vertex, int k, long double c)
{
	int j;
	long double hsum = 0, z, p, f, S, h0;
	long double feval[k + 1];

	z = sqr(c) * (1 - ECC) / sqr(cosl(vertex->lat));
	p = 1 - 2 * ECC * sqr(c);
	f = (1 - ECC * sqr(c)) * ECC * sqr(c);
	S = f + p * z - sqr(z);
	if (S < 0)
		S = 0;
	h0 = sqr(c) * (1 - ECC);

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
			j--;
			feval[j + 1] = - sqrtl(S) / (j * f * powl(z, j)) - p * (2 * j - 1) / (2 * j * f) * feval[j] + (j - 1) / (j * f) * feval[j - 1];
			j++;
			break;
		}

		hsum = hsum + powl(-h0, j) * binom(k, j) * feval[j];
	}

	return hsum / 2;
}

struct Vector sjoeberg(struct Coordinates *vertex, int i, int s, int a)
{
	struct vincenty_result inter[2];

	long double prev, next, excess = 0;
	long double area, darea = 0, interarea;
	long double c, c2;

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

			if ((vertex + h)->lon != (vertex + ((h + 1) % i))->lon && fabsl((vertex + h)->lat) != M_PI_2 && fabsl((vertex + ((h + 1) % i))->lat) != M_PI_2) {
				for (k = 1; k < 6; k++) {
					interarea = powl(ECC, k) * (k + 1) / (2 * k + 1);
					c = cosl(reduced_latitude((vertex + h)->lat)) * sinl(inter[0].start);
					c = c + cosl(reduced_latitude((vertex + ((h + 1) % i))->lat)) * sinl(inter[0].end);
					c = c / 2;

					interarea = interarea * (E(vertex + h, k, c) - E(vertex + ((h + 1) % i), k, c));
					darea = interarea + darea;
				}
			}

			inter[1] = inter[0];
		}
	}

	excess = excess - (i - 2) * M_PI;
	area = sqr(RAD_MIN) * (excess + darea);

	struct Vector res;
	res.s = perimeter;
	res.theta = area;
	return res;
}
