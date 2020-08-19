# Haversine
A program for calculating distance between coordinates.

# Compilation
`gcc -o haversine haversine.c -lm`

# Usage
Coordinates can be inputted from both command line and standard input (stdin).

## Format of coordinates
`[coordinates1] [coordinates2] ... [coordinatesN]`, where
each coordinate is in the form of `latitude,longitude`

For latitudes, positive is assumed for north.
For longitudes, positive is assumed for east.

## Output
Output is in JSON format:
```
{
  "1": xxx,
  "2": xxx,
  ...
  "N": xxx,
  "total": xxx
}
```
Where
each numbered pair specifies the distance between the nth point and (n+1)th point.
 
`"total"` specifies the total distance of the line joining all points.
 
All distances are provided in kilometres.
