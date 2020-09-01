#include <math.h>
#include "greatcircle.h"
#include "geodesic.h"
#include "haversine.h"

long double greatcircle_area(struct Coordinates *vertex, int i)
{
	long double next, prev, excess = 0;
	int k;
	for (k = 0; k < i; k++) {
		next = haversine_bearing(vertex + k, vertex + k + 1);
		prev = haversine_bearing(vertex + k, vertex + ((k + i - 1) % i));
		excess += normalise_a(next - prev);
	}

	excess = excess - (i - 2) * M_PI;

	return fabsl(sqr(RADIUS) * excess);
}

long double greatcircle_perimeter(struct Coordinates *vertex, int i)
{
	long double perimeter = 0;
	int k;
	for (k = 0; k < i; k++)
		perimeter += haversine_inverse_distance(vertex + k);

	return perimeter;
}
