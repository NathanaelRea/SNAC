# SNAC

Structural Numeric Analysis in C++

A very simple structural analysis solver nobody needed, and especially not in C++. Just a proof of concept for now.

Lmitations:
2D Linear Small Deformations
Only point and distributed loads on elements
Only outputs the solved u vector and local forces
Input file needs to be in correct order

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
Solved displacement vector
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

Local Forces for Element 0
   3.1116e+00
   1.7334e+01
   9.6683e+02
  -3.1116e+00
  -2.3337e+00
   2.1321e+02

Local Forces for Element 1
   5.6843e-14
   3.8895e+00
  -2.1321e+02
  -5.6843e-14
   6.1105e+00
  -1.2506e-12

Local Forces for Element 2
   6.1105e+00
   3.5527e-15
   6.8212e-13
  -6.1105e+00
  -3.5527e-15
            0
```
