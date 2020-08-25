#include "geodesic.h"

#ifndef __HAVE_MATHIO_H__

#define __HAVE_MATHIO_H__

int scan_coordinates(int argc, char *argv, struct Coordinates *loc);
int scan_vector(int argc, char *argv, struct Vector *vector);

void start_print(int i);
void error();

long double sqr(long double operand);
long double atan2_modified(long double y, long double x);

#endif
