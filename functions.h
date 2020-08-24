int scan(int argc, char *argv, struct Coordinates *loc)
{
	int count = 0;
	if (argc == 1)
		count = scanf("%Lf,%Lf", &loc->lat, &loc->lon);
	else
		sscanf(argv, "%Lf,%Lf", &loc->lat, &loc->lon);

	if (loc->lat < -90 || loc->lat > 90 || loc->lon < -180 || loc->lon > 180) {
		printf("Coordinate %Lf,%Lf out of range.\nAbort.\n", loc->lat, loc->lon);
		exit(1);
        }

	return count;
}

long double sqr(long double operand)
{
	return operand * operand;
}


long double atan2_modified(long double y, long double x)
{
	if (x > 0)
		return atanl(y / x);
	else if (x < 0) {
		if (y < 0)
			return atanl(y / x) - M_PI;
		else
			return atanl(y / x) + M_PI;
	}
	else {
		if (y > 0)
			return M_PI_2;
		else if (y < 0)
			return -M_PI_2;
		else
			return NAN;
	}
}
