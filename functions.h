int scan(int argc, char *argv, struct Coordinates *loc)
{
	int count = 0;
	if (argc == 1)
		count = scanf("%lf,%lf", &loc->lat, &loc->lon);
	else
		sscanf(argv, "%lf,%lf", &loc->lat, &loc->lon);
	return count;
}
