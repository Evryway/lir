# lir
Largest Interior Rectangle implementation in C# for Unity.

## REQUIREMENTS:

Unity.Mathematics is required for int2 and float2x2 usage (2D integer vectors and float 2x2 matrices, respectively).

The simplest way to get this package is to open the "Package Manager" window in Unity, click the plus icon, select "add package from git URL ..." and then enter "com.unity.mathematics" (without quotes).

## USAGE:

## LargestInteriorRectangle.cs 

the LargestInteriorRectangle class contains various static methods for processing Vector2 arrays.
The primary methods are:

#### CalculateInteriorCells(Vector2[] vs, out float[] xs, out float[] ys, out int[,] cells)

This method takes a list of Vector2s representing a simple polygon (ordered counter-clockwise, concave, non-intersecting) and
generates three output arrays : the xs (ordered, filtered x coordinates for each vertex), the ys
(ordered, filtered y coordinates for each vertex) and the 2D cells array representing each cell of the
projection of the xs and ys arrays into a 2D rectangle array. Each cell is marked as exterior (0) or interior
(1) based on it's position relative to the polygon.

#### CalculateLargestInteriorRectangle(float[] xs, float ys[], int[,] cells, out Bound2D best)

this method takes the xs, ys and cells arrays created by CalculateInteriorCells and calculates
the axis-aligned Largest Interior Rectangle. The bound of this rectangle is output in best.

Along with these methods, there are a variety of other geometry processing methods:

#### CalculateConcavePolygon


## Bound2D.cs

This class represents a 2D bound. The bound is not required to be axis aligned. The bound is
constructed from a centre point, a 2D size and the primary axis (the first element of the 2D size).

## Extensions.cs

This class contains some Vector2 extension methods to simplify the code (Rotate for simple
2D rotations, and F3() to make logging more concise).

## LargestInteriorRectangleTests.cs

This class contains a selection of tests across the LargestInteriorRectangle methods.


please see https://www.evryway.com/largest-interior/ for background details.


