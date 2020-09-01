#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <strings.h>

#include "geodesic.h"
#include "io.h"
#include "haversine.h"
#include "vincenty.h"
#include "greatcircle.h"
#include "karney.h"

void help(char *name)
{
	fputs("Usage:\n", stderr);
	fprintf(stderr, "\t%s [options] [coordinate 1] [coordinate 2] ...\n", name);
	fputs("\t-h Show usage.\n", stderr);

	fputs("Algorithm options:\n", stderr);
	fputs("\t-p [direct|inverse|polygon] Solve direct or inverse problem or polygon problem. \n", stderr);
	fputs("\t-f [sphere|ellipsoid] Set projection to sphere or ellipsoid.\n", stderr);

	fputs("Computation options:\n", stderr);
	fputs("\t-s Compute distances (a.k.a. perimeter) / coordinates.\n", stderr);
	fputs("\t-a Compute azimuths (for direct / inverse problem). Compute area (for polygon problem).\n", stderr);

	fputs("IO options:\n", stderr);
	fputs("\t-i [-|FILE] Input from stdin or FILE. stdin is assumed for - or missing argument.\n", stderr);
	fputs("\t-o [-|FILE] Output to stdout or FILE. stdout is assumed for - or missing argument.\n", stderr);
	fputs("\t-k [int] Precision of output values.\n", stderr);
	fputs("\nMore info in README.md.\n", stderr);
	return;
}


int main(int argc, char **argv)
{
	int i, j = 0, p = 0, count, precision = 6;
	int distance = 0, azimuth = 0;

	FILE *in, *out;
	in = stdin;
	out = stdout;

	while ((i = getopt(argc, argv, "i:o:k:p:f:sah")) != -1) {
		switch (i)
		{
			case 'i':
			if (strcmp(optarg, "-") != 0) {
				in = fopen(optarg, "r");
				if (in == NULL) {
					perror(optarg);
					error(NULL);
				}
			}
			break;

			case 'o':
			if (strcmp(optarg, "-") != 0)
				out = fopen(optarg, "w");
			break;

			case 'k':
			precision = atoi(optarg);
			if (precision < 0)
				error("Precision must be positive value");
			break;

			case 'p':
			if (strcmp(optarg, "direct") == 0)
				p = 1;
			else if (strcmp(optarg, "inverse") == 0)
				p = 2;
			else if (strcmp(optarg, "polygon") == 0)
				p = 3;
			else
				error("Invalid problem.");
			break;

			case 'f':
			if (strcmp(optarg, "sphere") == 0)
				j = 1;
			else if (strcmp(optarg, "ellipsoid") == 0)
				j = 2;
			else
				error("Invalid formula.");
			break;

			case 's':
			distance = 1;
			break;

			case 'a':
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

	if (distance == 0 && azimuth == 0)
		error("Nothing to be shown.");

	if (j == 0)
		error("No formula defined.");

	long double *res = malloc(sizeof(long double) * 3);
	char *writeout = calloc(1024, sizeof(char));

	switch (p)
	{
		case 2:
		{
		long double c, total = 0, start, end;
		struct Coordinates *location = malloc(sizeof(struct Coordinates) * 2);

		if (scan_coordinates(in, location) != 2)
			error("Incorrect format of coordinates.");

		for (i = 1; (count = scan_coordinates(in, location + 1)) == 2; i++) {
			start_print(writeout, i);

			switch (j)
			{
				case 1:
				if (distance == 1)
					c = haversine_inverse_distance(location);
				if (azimuth == 1) {
					start = haversine_bearing(location, location + 1);
					end = normalise_a(haversine_bearing(location + 1, location) - M_PI);
				}
				break;
				case 2:
				vincenty_inverse(location, location + 1, res, 3);
				c = *res;
				start = *(res + 1);
				end = *(res + 2);
				break;
			}

			if (isnan(c) && memcmp(location, location + 1, sizeof(struct Coordinates)) == 0)
				c = 0;

			if (distance == 1) {
				total += c;
				sprintf(writeout, "%s    \"distance\": %.*LF", writeout, precision, c);
				if (azimuth == 1)
					strcat(writeout, ",\n");
				else
					strcat(writeout, "\n  }");
			}

			if (azimuth == 1)
				sprintf(writeout, "%s    \"start_azimuth\": %.*LF,\n    \"end_azimuth\": %.*LF\n  }", writeout, precision, start / RAD, precision, end / RAD);

			memcpy(location, location + 1, sizeof(struct Coordinates));

			fputs(writeout, out);
			memset(writeout, 0, 1024);
		}

		free(location);
		if (count != -1)
			error("Incorrect format of coordinates.");
		if (distance == 1 && i != 1)
			fprintf(out, ",\n  {\n    \"total_distance\": %.*LF\n  }", precision, total);
		break;
		}

		case 1:
		{
		long double end;

		struct Coordinates *point = malloc(sizeof(struct Coordinates) * 2);
		struct Vector *add = malloc(sizeof(struct Vector));

		if (scan_coordinates(in, point) != 2)
			error("Incorrect format of coordinates.");

		for (i = 1; ((count = scan_vector(in, add)) == 2); i++) {
			start_print(writeout, i);
			switch (j)
			{
				case 1:
				*(point + 1) = haversine_direct(point, add);
				if (azimuth == 1)
					end = normalise_a(haversine_bearing(point + 1, point) - M_PI);
				break;
				case 2:
				vincenty_direct(point, add, res);
				(point + 1)->lat = *res;
				(point + 1)->lon = *(res + 1);
				end = *(res + 2);
				break;
			}

			if ((isnan((point + 1)->lat) || isnan((point + 1)->lon)) && add->s == 0)
				memcpy(point + 1, point, sizeof(struct Coordinates));

			if (distance == 1) {
				sprintf(writeout, "%s    \"coordinate\": [%.*LF,%.*LF]", writeout, precision, (point + 1)->lat / RAD, precision, (point + 1)->lon / RAD);
				if (azimuth == 1)
					strcat(writeout, ",\n");
			}

			if (azimuth == 1) {
				if (add->s == 0 && isnan(end))
					end = add->theta;
				sprintf(writeout, "%s    \"azimuth\": %.*LF", writeout, precision, end / RAD);
			}

			strcat(writeout, "\n  }");
			memcpy(point, point + 1, sizeof(struct Coordinates));

			fputs(writeout, out);
			memset(writeout, 0, 1024);
		}

		free(point);
		free(add);
		if (count != -1)
			error("Incorrect format of vector.");

		break;
		}

		case 3:
		{
		struct Coordinates *vertex = malloc(sizeof(struct Coordinates));
		long double area, perimeter, esterror;

		i = 0;

		while ((count = scan_coordinates(in, vertex + i)) == 2)
			vertex = realloc(vertex, sizeof(struct Coordinates) * (++i + 1));
		if (count != -1)
			error("Incorrect format of coordinates.");
		if (memcmp(vertex + --i, vertex, sizeof(struct Coordinates)) != 0)
			error("Not a closed polygon.");

		if (i > 1) {
			switch (j)
			{
				case 1:
				if (azimuth == 1)
					area = greatcircle_area(vertex, i);
				if (distance == 1)
					perimeter = greatcircle_perimeter(vertex, i);
				break;
				case 2:
				karney(vertex, i, distance, azimuth, res);
				perimeter = *res;
				area = *(res + 1);
				break;
			}
		}
		else {
			area = 0;
			perimeter = 0;
		}

		free(vertex);

		memcpy(writeout, "{\n", 2);
		if (azimuth == 1) {
			sprintf(writeout, "%s  \"area\": %.*LF", writeout, precision, area);
			if (distance == 1)
				strcat(writeout, ",\n");
			else
				strcat(writeout, "\n}\n");
		}
		if (distance == 1)
			sprintf(writeout, "%s  \"perimeter\": %.*LF\n}\n", writeout, precision, perimeter);
		fputs(writeout, out);
		break;
		}

		default:
		error("No problem defined.");
		break;
	}

	free(res);
	free(writeout);
	if (in != stdin)
		fclose(in);

	if (p < 3) {
		if (i == 1)
			fputc('[', out);
		fputs("\n]\n", out);
	}
	if (out != stdout)
		fclose(out);

	exit(0);
}
