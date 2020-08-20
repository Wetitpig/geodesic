int scan(int argc, char *argv, struct Coordinates *loc)
{
	int count = 0;
	if (argc == 1)
		count = scanf("%Lf,%Lf", &loc->lat, &loc->lon);
	else
		sscanf(argv, "%Lf,%Lf", &loc->lat, &loc->lon);
	return count;
}
