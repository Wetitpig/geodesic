#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <strings.h>

#include "geodesic.h"
#include "io.h"
#include "haversine.h"
#include "vincenty.h"

void help(char *name)
{
	fputs("Usage:\n", stderr);
	fprintf(stderr, "\t%s [options] [coordinate 1] [coordinate 2] ...\n", name);
	fputs("\t-h Show usage.\n", stderr);

	fputs("Computation options:\n", stderr);
	fputs("\t-p [direct|inverse] Solve direct or inverse problem. \n", stderr);
	fputs("\t-f [haversine|vincenty] Set formula to haversine or Vincenty's.\n", stderr);
	fputs("\t-s Compute distances / coordinates.\n", stderr);
	fputs("\t-a Compute azimuths.\n", stderr);

	fputs("IO options:\n", stderr);
	fputs("\t-i [-|FILE] Input from stdin or FILE. stdin is assumed for - or missing argument.\n", stderr);
	fputs("\t-o [-|FILE] Output to stdout or FILE. stdout is assumed for - or missing argument.\n", stderr);
	fputs("\nMore info in README.md.\n", stderr);
	return;
}


int main(int argc, char **argv)
{
	int i, j = 0, p = 0, count;
	int distance = 0, azimuth = 0;

	FILE *in, *out;
	in = stdin;
	out = stdout;

	while ((i = getopt(argc, argv, "i:o:p:f:sah")) != -1) {
		switch (i)
		{
			case 'i':
			if (strcmp(optarg, "-") != 0) {
				in = fopen(optarg, "r");
				if (in == NULL) {
					perror(optarg);
					error();
				}
			}
			break;

			case 'o':
			if (strcmp(optarg, "-") != 0)
				out = fopen(optarg, "w");
			break;

			case 'p':
			if (strcmp(optarg, "direct") == 0)
				p = 1;
			else if (strcmp(optarg, "inverse") == 0)
				p = 2;
			else {
				fputs("Invalid problem.", stderr);
				error();
			}
			break;

			case 'f':
			if (strcmp(optarg, "haversine") == 0)
				j = 1;
			else if (strcmp(optarg, "vincenty") == 0)
				j = 2;
			else {
				fputs("Invalid formula.", stderr);
				error();
			}
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

	if (distance == 0 && azimuth == 0) {
		fputs("Nothing to be shown.", stderr);
		error();
	}

	if (j == 0) {
		fputs("No formula defined.", stderr);
		error();
	}

	struct vincenty_result res;
	char *writeout = calloc(1024, sizeof(char));

	switch (p)
	{
		case 2:
		{
		long double c, total = 0, start, end;
		struct Coordinates *location = malloc(sizeof(struct Coordinates) * 2);

		if (scan_coordinates(in, location) != 2) {
			fputs("Incorrect format of coordinates.", stderr);
			error();
		}

		for (i = 1; (count = scan_coordinates(in, location + 1)) == 2; i++) {
			start_print(writeout, i);

			switch (j)
			{
				case 1:
				if (distance == 1)
					c = haversine_inverse_distance(location);
				if (azimuth == 1) {
					start = haversine_bearing(location, location + 1);
					end = haversine_bearing(location + 1, location);
				}
				break;
				case 2:
				res = vincenty_inverse(location);
				c = res.distance;
				start = res.start;
				end = res.end;
				break;
			}

			if (isnan(c) && location->lat == (location + 1)->lat && location->lon == (location + 1)->lon)
				c = 0;

			if (distance == 1) {
				total += c;
				sprintf(writeout, "%s    \"distance\": %LF", writeout, c);
				if (azimuth == 1)
					strcat(writeout, ",\n");
				else
					strcat(writeout, "\n  }");
			}

			if (azimuth == 1)
				sprintf(writeout, "%s    \"start_azimuth\": %LF,\n    \"end_azimuth\": %LF\n  }", writeout, start / RAD, end / RAD);

			memcpy(location, location + 1, sizeof(struct Coordinates));

			fputs(writeout, out);
			memset(writeout, 0, 1024);
		}

		free(location);
		if (count != -1) {
			fputs("Incorrect format of coordinates.", stderr);
			error();
		}
		if (distance == 1 && i != 1)
			fprintf(out, ",\n  {\n    \"total_distance\": %LF\n  }", total);
		break;
		}

		case 1:
		{
		long double end;

		struct Coordinates *point = malloc(sizeof(struct Coordinates) * 2);
		struct Vector *add = malloc(sizeof(struct Vector));

		if (scan_coordinates(in, point) != 2) {
			fputs("Incorrect format of coordinates.", stderr);
			error();
		}

		for (i = 1; ((count = scan_vector(in, add)) == 2); i++) {
			start_print(writeout, i);
			switch (j)
			{
				case 1:
				*(point + 1) = haversine_direct(point, add);
				if (azimuth == 1)
						end = haversine_bearing(point + 1, point);
				break;
				case 2:
				res = vincenty_direct(point, add);
				(point + 1)->lat = res.distance;
				(point + 1)->lon = res.start;
				end = res.end;
				break;
			}

			if ((isnan((point + 1)->lat) || isnan((point + 1)->lon)) && add->s == 0)
				memcpy(point + 1, point, sizeof(struct Coordinates));

			if (distance == 1) {
				sprintf(writeout, "%s    \"coordinate\": [%LF,%LF]", writeout, (point + 1)->lat / RAD, (point + 1)->lon / RAD);
				if (azimuth == 1)
					strcat(writeout, ",\n");
			}

			if (azimuth == 1) {
				if (add->s == 0 && isnan(end))
					end = add->theta;
				sprintf(writeout, "%s    \"azimuth\": %LF", writeout, end / RAD);
			}

			strcat(writeout, "\n  }");
			memcpy(point, point + 1, sizeof(struct Coordinates));

			fputs(writeout, out);
			memset(writeout, 0, 1024);
		}

		free(point);
		free(add);
		if (count != -1) {
			fputs("Incorrect format of vector.", stderr);
			error();
		}

		break;
		}
		default:
		fputs("No problem defined.", stderr);
		error();
		break;
	}

	free(writeout);
	if (in != stdin)
		fclose(in);

	if (i == 1)
		fputc('[', out);
	fputs("\n]\n", out);
	if (out != stdout)
		fclose(out);

	exit(0);
}
