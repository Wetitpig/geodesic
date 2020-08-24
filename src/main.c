#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <strings.h>

#include "constants.h"
#include "functions.h"
#include "haversine.h"
#include "vincenty.h"

void help(char *name)
{
	puts("Usage:");
	printf("\n\t%s [options] [coordinate 1] [coordinate 2] ...\n", name);
	printf("\t-h Show usage.\n\t-f [haversine|vincenty] Set formula to haversine or vincenty.\n\t-s Compute distance.\n\t-o Compute azimuths.\n\n");
	puts("More info in README.md");
	return;
}


int main(int argc, char **argv)
{
	int i, j;
	int distance = 0, azimuth = 0;

	while ((i = getopt(argc, argv, "f:soh")) != -1) {
		switch (i)
		{
			case 'f':
			if (strcmp(optarg, "haversine") == 0)
				j = 1;
			else if (strcmp(optarg, "vincenty") == 0)
				j = 2;
			else {
				puts("Invalid formula. Abort.");
				exit(1);
			}
			break;
			case 's':
			distance = 1;
			break;
			case 'o':
			azimuth = 1;
			break;
			case 'h':
			help(argv[0]);
			exit(1);
			break;
			default:
			exit(1);
			break;
		}
	}

	if (distance == 0 && azimuth == 0) {
		puts("Nothing to be shown. Abort.");
		exit(1);
	}

	struct Coordinates *location = malloc(sizeof(struct Coordinates) * 2);

	long double c, total = 0, start, end;
	struct vincenty_result *res;

	scan(optind == argc, argv[optind], location);

	location->lat = location->lat * RAD;
	location->lon = location->lon * RAD;

	putchar('{');

	for (i = 1; optind == argc || (optind + i) < argc; ++i) {
		if (optind != argc)
			scan(argc, argv[optind + i], location + 1);
		else if (scan(1, NULL, location + 1) == -1)
			break;

		(location + 1)->lat = (location + 1)->lat * RAD;
		(location + 1)->lon = (location + 1)->lon * RAD;

		if (azimuth == 1 && i != 1)
			printf(",");

		printf("\n  \"%d\": {\n", i - 1);

		switch (j)
		{
			case 1:
			if (distance == 1)
				c = haversine_distance(location);
			if (azimuth == 1) {
				start = NORMALISE(haversine_bearing(location, location + 1));
				end = NORMALISE(haversine_bearing(location + 1, location) - 180);
			}
			break;
			case 2:
			res = malloc(sizeof(struct vincenty_result));
			vincenty(res, location);
			c = res->distance;
			start = res->start;
			end = res->end;
			free(res);
			break;
		}

		if (distance == 1) {
			total += c;
			printf("    \"distance\": %Lf", c);
			if (azimuth == 1)
				printf(",\n");
			else
				printf("\n  },");
		}

		if (azimuth == 1)
			printf("    \"start_azimuth\": %Lf,\n    \"end_azimuth\": %Lf\n  }", start, end);

		memcpy(location, location + 1, sizeof(struct Coordinates));
	}

	free(location);

	if (distance == 1) {
		if (azimuth == 1 && total != 0)
			printf(",");
		printf("\n  \"total_distance\": %Lf\n}\n", total);
	}
	else
		printf("\n}\n");

	exit(0);
}
