# SNAC

Structural Numeric Analysis in C++

A structural analysis solver nobody needed in C++. Just a proof of concept for now, can only add point and distributed loads on elements, and only outputs the solved u vector.

## Compile with correct libraries

```
$ g++ -std=c++17 -O2 -larmadillo -llapack -lblas
```

## Example input file

```
POINTS
0 0
72 96
264 96
264 0

ELEMENTS
0 1 29000 400 10
1 2 29000 400 10
2 3 29000 400 10

BOUNDS
0 X 0
0 Y 0
0 R 0
3 Y 0

LOADING
ELEMENT 0 DISTRIBUTED -0.125 -0.125
ELEMENT 1 POINT -10 96

```

Output:

```
Solved u
        0
        0
        0
   0.2095
  -0.1587
  -0.0023
   0.2095
  -0.0020
   0.0034
   0.5350
        0
   0.0034
```
