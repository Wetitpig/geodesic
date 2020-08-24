# Geodesic
Programs for calculating distance between coordinates.

## Formula
* [Haversine formula](https://en.wikipedia.org/wiki/Haversine_formula)
* [Vincenty's formula](https://en.wikipedia.org/wiki/Vincenty%27s_formulae)

## Compilation
```
make
make clean
```
A binary will be located in `./bin` directory (`bin/geodesic`).

## Usage
```
geodesic [-f formula] [-s] [-o] [-h] [coordinate 1] [coordinate 2]
```
### Options
`-f [haversine|vincenty]` Formula for evaluating the azimuths and distances.
Current available formula are haversine and Vincenty's formula. Haversine formula is more resource-efficient, while Vincenty's formula is more accurate. Resource consumption effects are not noticeable for just a few coordinates.

`-s` Show distance between coordinates.

`-o` Show azimuth from one coordinate to the next.

### Arguments
Coordinates can be inputted from both command line and standard input (stdin) in the form of:

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
* `distance` as the distance between the nth coordinate and the (n+1)th coordinate.
* `start_azimuth` as the bearing at the nth coordinate, on the line from the nth coordinate to the (n+1)th coordinate.
* `end_azimuth` as the bearing at the (n+1)th coordinate, on the line from the nth coordinate to the (n+1)th coordinate.
 
`"total_distance"` specifies the total distance of the line joining all coordinates.

### Units
* All distances are provided in kilometres.
* All angles are provided in degrees.
