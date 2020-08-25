#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <strings.h>

#include "geodesic.h"
#include "mathio.h"
#include "haversine.h"
#include "vincenty.h"

void help(char *name)
{
	puts("Usage:");
	printf("\n\t%s [options] [coordinate 1] [coordinate 2] ...\n", name);
	printf("\t-h Show usage.\n\t-p [direct|inverse] Solve direct or inverse problem. \n\t-f [haversine|vincenty] Set formula to haversine or Vincenty's.\n\t-s Compute distances / coordinates.\n\t-o Compute azimuths.\n\n");
	puts("More info in README.md.");
	return;
}


int main(int argc, char **argv)
{
	int i, j = 0, p = 0, count;
	int distance = 0, azimuth = 0;

	while ((i = getopt(argc, argv, "p:f:soh")) != -1) {
		switch (i)
		{
			case 'p':
			if (strcmp(optarg, "direct") == 0)
				p = 1;
			else if (strcmp(optarg, "inverse") == 0)
				p = 2;
			else {
				puts("Invalid problem.");
				error();
			}
			break;
			case 'f':
			if (strcmp(optarg, "haversine") == 0)
				j = 1;
			else if (strcmp(optarg, "vincenty") == 0)
				j = 2;
			else {
				puts("Invalid formula.");
				error();
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
		puts("Nothing to be shown.");
		error();
	}

	if (j == 0) {
		puts("No formula defined.");
		error();
	}

	struct vincenty_result res;

	switch (p)
	{
		case 2:
		{
		long double c, total = 0, start, end;
		struct Coordinates *location = malloc(sizeof(struct Coordinates) * 2);

		if (scan_coordinates(optind == argc, argv[optind], location) != 2) {
			puts("Incorrect format of coordinates.");
			error();
		}

		for (i = 1; optind == argc || (optind + i) < argc; ++i) {
			if ((count = scan_coordinates(optind == argc, argv[optind + i], location + 1)) != 2) {
				if (optind == argc && count == -1)
					break;
				else {
					puts("Incorrect format of coordinates.");
					error();
				}
			}

			start_print(i);
			switch (j)
			{
				case 0:
				if (distance == 1)
					c = haversine_inverse_distance(location);
				if (azimuth == 1) {
					start = NORMALISE_A(haversine_inverse_bearing(location, location + 1));
					end = NORMALISE_A(haversine_inverse_bearing(location + 1, location) - 180);
				}
				break;
				case 1:
				res = vincenty_inverse(location);
				c = res.distance;
				start = NORMALISE_A(res.start);
				end = NORMALISE_A(res.end);
				break;
			}

			if (distance == 1) {
				total += c;
				printf("    \"distance\": %Lf", c);
				if (azimuth == 1)
					printf(",\n");
				else
					printf("\n  }");
			}

			if (azimuth == 1)
				printf("    \"start_azimuth\": %Lf,\n    \"end_azimuth\": %Lf\n  }", start, end);

			memcpy(location, location + 1, sizeof(struct Coordinates));
		}

		free(location);

		if (distance == 1 && i != 1)
			printf("\n  {\n    \"total_distance\": %Lf\n  }", total);
		break;
		}

		case 1:
		{
		long double end;

		struct Coordinates *point = malloc(sizeof(struct Coordinates) * 2);
		struct Vector *add = malloc(sizeof(struct Vector));

		if (scan_coordinates(optind == argc, argv[optind], point) != 2) {
			puts("Incorrect format of coordinates.");
			error();
		}


		for (i = 1; optind == argc || (optind + i) < argc; i++) {
			if ((count = scan_vector(argc == optind, argv[optind + i], add)) != 2) {
				if (optind == argc && count == -1)
					break;
				else {
					puts("Incorrect format of vector.");
					error();
				}
			}

			start_print(i);
			switch (j)
			{
				case 0:
				*(point + 1) = haversine_direct(point, add);
				if (azimuth == 1)
					end = NORMALISE_A(haversine_inverse_bearing(point + 1, point) - 180);
				break;
				case 1:
				res = vincenty_direct(point, add);
				(point + 1)->lat = res.distance;
				(point + 1)->lon = res.start;
				end = NORMALISE_A(res.end);
				break;
			}

			if (distance == 1) {
				printf("    \"coordinate\": [%Lf,%Lf]", (point + 1)->lat / (RAD), (point + 1)->lon / (RAD));
				if (azimuth == 1)
					printf(",\n");
			}

			if (azimuth == 1)
				printf("    \"azimuth\": %Lf", end);

			printf("\n  }");
			memcpy(point, point + 1, sizeof(struct Coordinates));
		}

		free(point);
		free(add);

		break;
		}
		default:
		puts("No problem defined.");
		error();
		break;
	}

	if (i == 1)
		putchar('[');
	printf("\n]\n");
	exit(0);
}
