#ifndef __HAVE_CONSTANTS_H__

#define __HAVE_CONSTANTS_H__

#define RADIUS 6371.0088
#define RAD M_PI / 180

#define RAD_MAJ 6378.137
#define RAD_MIN 6356.752314245

#define FLAT 1/298.257223563

struct Coordinates {
	long double lat;
	long double lon;
};

#endif
