# Geodesic
A program for calculating distance between coordinates.

## Formulae
* [Haversine formula](https://en.wikipedia.org/wiki/Haversine_formula)
* [incenty formula](https://en.wikipedia.org/wiki/Vincenty%27s_formulae)

## Compilation
* All: `make`
* haversine/vincenty: `make [haversine|vincenty]`

All binaries will be located in `./bin` directory.

## Usage
Coordinates can be inputted from both command line and standard input (stdin).

### Format of coordinates
`[coordinates1] [coordinates2] ... [coordinatesN]`, where
each coordinate is in the form of `latitude,longitude` in decimals.

For latitudes, positive is assumed for north.
For longitudes, positive is assumed for east.

### Output
Output is in JSON format:
```
{
  "0": {
    "distance": xxx,
    "start_azimuth": xxx,
    "end_azimuth": xxx
  },
  "1": {
    "distance": xxx,
    "start_azimuth": xxx,
    "end_azimuth": xxx
  },
  ...
  "N": {
    "distance": xxx,
    "start_azimuth": xxx,
    "end_azimuth": xxx
  },
  "total_distance": xxx
}
```
Where
each numbered pair specifies the following:
* `distance` as the distance between the nth point and the (n+1)th point.
* `start_azimuth` as the bearing at the nth point, on the line from the nth point to the (n+1)th point.
* `end_azimuth` as the bearing at the (n+1)th point, on the line from the nth point to the (n+1)th point.
 
`"total_distance"` specifies the total distance of the line joining all points.

### Units
* All distances are provided in kilometres.
* All angles are provided in degrees.
