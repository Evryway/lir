using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;

using LIR = Evryway.LargestInteriorRectangle;

namespace Evryway
{

    public static class LargestInteriorRectangleTests
    {

        // this test calculates the largest interior rectangle for a defined set of vertices
        // which form a simple polygon.

        public static void TestLargestInteriorRectangle()
        {

            Vector2[] vs = new Vector2[]
            {
                new Vector2(9,-7),
                new Vector2(12,-6),
                new Vector2(8,3),
                new Vector2(10,6),
                new Vector2(12,7),
                new Vector2(1,9),
                new Vector2(-8,7),
                new Vector2(-6,6),
                new Vector2(-4,6),
                new Vector2(-6,2),
                new Vector2(-6,0),
                new Vector2(-7,-5),
                new Vector2(-2,-7),
                new Vector2(1,-3),
                new Vector2(5,-7),
                new Vector2(8,-4),

            };

            LIR.CalculateInteriorCells(vs, out var xs, out var ys, out int[,] cells);
            LIR.CalculateLargestInteriorRectangle(xs, ys, cells, out var best);

            Debug.Log($"{best.centre.F3()}, {best.size.F3()}, {best.axis_a.F3()}");
        }


        public static void TestCovarianceMatrix()
        {
            /*
            Vector2[] v2s = new Vector2[]
            {
                new Vector2(1,1),
                new Vector2(3,0),
                new Vector2(-1,-1),
            };
            */

            //https://www.wikihow.com/Calculate-Covariance
            Vector2[] v2s = new Vector2[]
            {
                new Vector2(1,8),
                new Vector2(3,6),
                new Vector2(2,9),
                new Vector2(5,4),
                new Vector2(8,3),
                new Vector2(7,3),
                new Vector2(12,2),
                new Vector2(2,7),
                new Vector2(4,7),
            };

            var m = LIR.CovarianceMatrix(v2s);

            Debug.Log($"{m.c0.x}, {m.c0.y}\n{m.c1.x}, {m.c1.y}");
        }


        public static void TestCalculateEigenValues()
        {
            //float2x2 A = new float2x2(5, -3, -6, 2);
            //float2x2 A = new float2x2(1, 2, 4, 3);
            float2x2 A = new float2x2(5, 4, 3, 2);
            TestCalculateEigenValues(A);
        }

        public static void TestCalculateEigenValues(float2x2 A)
        {
            var ok = LIR.CalculateEigenValues(A, out float v1, out float v2);
            Debug.Log($"{ok} {v1} {v2}");
        }



        public static void TestCalculateEigenVectors()
        {
            //float2x2 A = new float2x2(5, -3, -6, 2);
            float2x2 A = new float2x2(1, 2, 4, 3);
            TestCalculateEigenVectors(A);
        }

        public static void TestCalculateEigenVectors(float2x2 A)
        {

            var ok = LIR.CalculateEigenValues(A, out float v1, out float v2);
            if (!ok) return;

            var vec2a = LIR.CalculateEigenVector(A, v1);
            var vec2b = LIR.CalculateEigenVector(A, v2);

            Debug.Log($"{v1} : {vec2a.x} {vec2a.y}");
            Debug.Log($"{v2} : {vec2b.x} {vec2b.y}");

        }



        public static void TestPrimaryAxis()
        {
            Vector2[] v2s = new Vector2[]
            {
                new Vector2(0,0),
                new Vector2(4,2),
                new Vector2(7,3.5f),
                new Vector2(-10,-5),
                new Vector2(30,15),
            };
            TestPrimaryAxis(v2s);
        }



        public static Vector2 TestPrimaryAxis(Vector2[] vs)
        {
            var ok = LIR.CalculatePrimaryAxis(vs, out Vector2 axis, out float eigenvalue);

            //Debug.Log($"{ok} {axis.x} {axis.y}");
            return axis;

        }



        public static bool TestConvexHull()
        {
            bool ok = false;

            Vector2[] hull_vs;
            Vector2[] vs = new Vector2[0];
            ok = LIR.CalculateConvexHull(vs, out hull_vs);
            if (ok || hull_vs.Length > 0) Debug.LogWarning("convex hull should fail and have no verts.");

            vs = new Vector2[] { new Vector2(1, 1) };
            LIR.CalculateConvexHull(vs, out hull_vs);
            if (ok || hull_vs.Length != 1) Debug.LogWarning("convex hull should fail and have 1 vert.");


            // all duplicate points.
            vs = new Vector2[] {
                new Vector2(1,1),
                new Vector2(1,1),
                new Vector2(1,1),
                new Vector2(1,1),
            };
            ok = LIR.CalculateConvexHull(vs, out hull_vs);
            if (ok || hull_vs.Length != 1) Debug.LogWarning("convex hull should fail and have 1 vert.");

            // two points.
            vs = new Vector2[]
            {
                new Vector2(-1,0),
                new Vector2(1,0),
            };
            ok = LIR.CalculateConvexHull(vs, out hull_vs);
            if (ok || hull_vs.Length != 2) Debug.LogWarning("convex hull should fail and have 2 verts.");


            // three points, colinear.
            vs = new Vector2[]
            {
                new Vector2(-1,0),
                new Vector2(0.5f, 0),
                new Vector2(1,0),
            };
            ok = LIR.CalculateConvexHull(vs, out hull_vs);
            if (ok || hull_vs.Length != 2) Debug.LogWarning("convex hull should fail with 2 verts");

            // six points, colinear and duplicates.
            vs = new Vector2[]
            {
                new Vector2(-1,0),
                new Vector2(-1,0),
                new Vector2(0.5f, 0),
                new Vector2(-0.5f, 0),
                new Vector2(1,0),
                new Vector2(1,0),
            };
            ok = LIR.CalculateConvexHull(vs, out hull_vs);
            if (ok || hull_vs.Length != 2) Debug.LogWarning("convex hull should fail with 2 verts");


            // six points, colinear
            vs = new Vector2[]
            {
                new Vector2(-1,0),
                new Vector2(-0.6f,0),
                new Vector2(0.3f, 0),
                new Vector2(-0.8f, 0),
                new Vector2(1.2f,0),
                new Vector2(1.9f,0),
            };
            ok = LIR.CalculateConvexHull(vs, out hull_vs);
            if (ok || hull_vs.Length != 2) Debug.LogWarning("convex hull should fail with 2 verts");



            // three points, valid.
            vs = new Vector2[]
            {
                new Vector2(-1,-1),
                new Vector2(1,-1),
                new Vector2(0,1),
            };
            ok = LIR.CalculateConvexHull(vs, out hull_vs);
            if (!ok || hull_vs.Length != 3) Debug.LogWarning("convex hull should pass with 3 verts");

            //Debug.Log("tests complete.");

            return true;
        }


        public static bool TestSmallestEnclosingRectangle(Vector2[] vs, out Bound2D bound)
        {
            var ok = LIR.CalculateSmallestEnclosingRectangle(vs, out bound);
            if (!ok)
            {
                Debug.LogWarning("failed smallest enclosing rectangle test.");
                return false;
            }
            return true;
        }



    }
}
