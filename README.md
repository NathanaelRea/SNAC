# SNAC

Structural Numeric Analysis in C++

A very simple structural analysis solver nobody needed, and especially not in C++. Just a proof of concept for now.

Lmitations:
2D Linear Small Deformations
Only point and distributed loads on elements
Only outputs the solved u vector and local forces

## Compile with correct libraries

```
$ g++ -std=c++17 -O2 -larmadillo -llapack -lblas
```

## Input Syntax
```
POINTS
<X> <Y>

ELEMENTS
<P#1> <P#2> <E> <I> <A>

BOUNDS
<P#> <DIR> <DISP>

LOADING
ELEMENT <E#> DISTRIBUTED <W1> <W2>
ELEMENT <E#> POINT <LOAD> <DIST>
POINT <P#> <DIR> <LOAD>
```
Points and Elements are 0 indexed


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
Point Displacement (X,Y,R)
Point 0
(0.000000, 0.000000, 0.000000)
Point 1
(0.209509, -0.158741, -0.002346)
Point 2
(0.209509, -0.002023, 0.003391)
Point 3
(0.535009, 0.000000, 0.003391)

Local Element End Forces (A,V,M)
Element 0
(3.111614, 17.333710, 966.832579)
(-3.111614, -2.333710, 213.212670)
Element 1
(0.000000, 3.889517, -213.212670)
(-0.000000, 6.110483, -0.000000)
Element 2
(6.110483, 0.000000, 0.000000)
(-6.110483, -0.000000, 0.000000)

```
