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
geodesic [-p problem] [-f formula] [-s] [-a] [-h] [arguments]
```
### Options
`-p [direct|inverse]` Solve direct problem or inverse problem.
* Direct problem: Given a coordinate and a vector of distance and initial bearing, evaluate the destination coordinate and the final bearing.
* Inverse problem: Given 2 coordinates, evaluate the distance, the initial bearing and the final bearing.

`-f [haversine|vincenty]` Select formula for evaluating the azimuths and distances.
Current available formula are haversine and Vincenty's formula. Haversine formula is more resource-efficient, while Vincenty's formula is more accurate. Resource consumption effects are not noticeable for just a few coordinates.

`-s` Show distance between coordinates for inverse problems, or show the destination coordinate for direct problems.

`-a` Show azimuth from one coordinate to the next for inverse problems, or show the final bearing for direct problems.

### Arguments
All arguments can be given from the command line or standard input (stdin).

### Output
Output is in JSON format.

### Inverse problem

#### Arguments
Coordinates are given in the form of:

`[coordinates1] [coordinates2] ... [coordinatesN]`, where
each coordinate is in the form of `latitude,longitude` in decimals.

#### Output
```
[
  {
    "distance": xxx,
    "start_azimuth": xxx,
    "end_azimuth": xxx
  },
  {
    "distance": xxx,
    "start_azimuth": xxx,
    "end_azimuth": xxx
  },
  ...
  {
    "distance": xxx,
    "start_azimuth": xxx,
    "end_azimuth": xxx
  },
  {
    "total_distance": xxx
  }
]
```
Where
each set of pairs in array specifies the following:
* `distance` as the distance between the nth coordinate and the (n+1)th coordinate.
* `start_azimuth` as the bearing at the nth coordinate, on the line from the nth coordinate to the (n+1)th coordinate.
* `end_azimuth` as the bearing at the (n+1)th coordinate, on the line from the nth coordinate to the (n+1)th coordinate.
 
`"total_distance"` specifies the total distance of the line joining all coordinates.

### Direct problem

#### Arguments
Vectors are given in the form of:

`[start coordinate] [vector1] [vector2] ... [vectorN]`, where
each vector is in the form of `distance:bearing`, in decimals and degrees respectively.

#### Output
```
[
  {
    "coordinate": [xxx,xxx],
    "azimuth": xxx
  },
  {
    "coordinate": [xxx,xxx],
    "azimuth": xxx
  },
  ...
  {
    "coordinate": [xxx,xxx],
    "azimuth": xxx
  }
]
```
Where
each set of pairs in array specifies the following:
* `coordinate` as the destination coordinate at the distance and bearing from the (n-1)th coordinate.
* `azimuth` as the bearing at `coordinate` (nth coordinate), on the line from (n-1)th coordinate to nth coordinate.

### Units
* All distances are provided in kilometres.
* All angles are provided in degrees.
* For latitudes, positive is assumed for north.
* For longitudes, positive is assumed for east.
* For bearings, full bearing must be used.

## References
1. [Movable Type Scripts](https://www.movable-type.co.uk/scripts/latlong.html)
2. [Vincenty's formulae @ Wikipedia](https://en.wikipedia.org/wiki/Vincenty%27s_formulae)
3. [Geodetic Inverse Solution between Antinodal Points](https://geographiclib.sourceforge.io/geodesic-papers/vincenty75b.pdf)
