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
geodesic [-i file] [-o file] [-k precision] [-p problem] [-f projection] [-s] [-a] [-h]
```
### Computation Options
`-p [direct|inverse|polygon]` Solve direct problem or inverse problem.
* Direct problem: Given a coordinate and a vector of distance and initial bearing, evaluate the destination coordinate and the final bearing.
* Inverse problem: Given 2 coordinates, evaluate the distance, the initial bearing and the final bearing.
* Polygon problem: Given a set of coordinates, evaluate the perimeter and area.

`-f [sphere|ellipsoid]` Select projection for evaluating the azimuths and distances.
* Sphere projection employs the haversine algorithm. Error is higher (error ~0.5%).
* Ellipsoid projection employs the Vincenty's algorithm. It is more resource-intensive (2x sphere projection).

`-s` Different computations for different problems.
* Distance between coordinates for inverse problems.
* The destination coordinate for direct problems.
* The perimeter for polygon problems.

`-a` Different computations for different problems.
* Initial and final bearings from one coordinate to the next for inverse problems.
* The final bearing for direct problems.
* The area for polygon problems.

### IO Options
`-i [FILE|-]` Input file. `-` signifies standard in (stdin). Defaults to stdin without the option.

`-o [FILE|-]` Output file. `-` signifies standard out (stdout). Defaults to stdout without the option.

`-k [integer]` Set the number of decimal places in output. Default is 6.

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

### Polygon problem

### Arguments
Coordinates are given in the form of:

`[coordinates1] [coordinates2] ... [coordinatesN] [coordinates1]`, where
each coordinate is in the form of `latitude,longitude` in decimals.

Coordinates of a polygon should be specified in the **counter-clockwise** direction.

For polygons that include the poles, coordinates should be specified in the **eastward** direction for the North Pole, and the **westward** direction for the South Pole.

Note that `coordinate1` must be placed at the end to close the polygon.

#### Output
```
{
  "perimeter": xxx,
  "area": xxx
}
```
Where
* `perimeter` denotes the perimeter of the polygon.
* `area` denotes the area of the polygon.

### Units
* All distances are provided in kilometres.
* All areas are provided in square kilometres.
* All angles are provided in degrees.
* For latitudes, positive is assumed for north.
* For longitudes, positive is assumed for east.
* For bearings, full bearing must be used.

### Notes
* `allow_nan` should be set to `true` for parsing JSON, because azimuths may return a `NaN` value for zero vectors depending on the algorithm.

## References
1. [Movable Type Scripts](https://www.movable-type.co.uk/scripts/latlong.html)
2. [Vincenty's formulae @ Wikipedia](https://en.wikipedia.org/wiki/Vincenty%27s_formulae)
3. [Geodetic Inverse Solution between Antinodal Points](https://geographiclib.sourceforge.io/geodesic-papers/vincenty75b.pdf)
4. [Determination of areas on the plane, sphere and ellipsoid](https://www.tandfonline.com/doi/abs/10.1179/sre.2006.38.301.583?journalCode=ysre20)
