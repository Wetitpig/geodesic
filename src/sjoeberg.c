#include <math.h>
#include "vincenty.h"
#include "sjoeberg.h"

long double sjoeberg_area(struct Coordinates *vertex, int i)
{
	return 0;
}

long double sjoeberg_perimeter(struct Coordinates *vertex, int i)
{
	struct vincenty_result inter;
	int k;
	long double perimeter = 0;
	for (k = 0; k < i; k++) {
		inter = vincenty_inverse(vertex + k);
		perimeter += inter.distance;
	}
	return perimeter;
}
