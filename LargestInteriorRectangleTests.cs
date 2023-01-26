using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;

using LIR = Evryway.LargestInteriorRectangle;

namespace Evryway
{

    public static class LargestInteriorRectangleTests
    {
        public static bool TestPolygon(Vector2[] vs, bool strict = true)
        {
            var ok = true;
            if (strict)
            {
                // construct polygon info and check that we're convex and counter-clockwise.
                var info = PolygonInfo.Create(vs);
                ok = info.valid && !info.clockwise;
                if (!ok)
                {
                    if (!info.valid) Debug.LogError("invalid points.");
                    if (info.clockwise) Debug.LogError("winding order is not counter-clockwise.");
                    return false;
                }
            }


            ok = LIR.CalculateInteriorCells(vs, out var xs, out var ys, out int[,] cells);
            if (!ok)
            {
                Debug.LogError("points are not valid.");
                return false;
            }
            ok |= LIR.CalculateLargestInteriorRectangle(xs, ys, cells, out var best);
            if (!ok)
            {
                Debug.LogError("failed to calculate Largest Interior Rectangle.");
                return false;
            }

            Debug.Log("Test passed.");
            Debug.Log($"{best.centre.F3()}, {best.size.F3()}, {best.axis_a.F3()}");
            return true;
        }
        
        public static Vector2[] Reverse(IEnumerable<Vector2> src)
        {
            var rev = new List<Vector2>(src);
            rev.Reverse();
            return rev.ToArray();
        }

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

            TestPolygon(vs);

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


        // specific raised issues
        // https://github.com/Evryway/lir/issues/4
        public static bool TestIssue4a_IndexOutOfRange()
        {
            Vector2[] vs = new Vector2[] {
                new Vector2(-6.012309074f,  -2.394122601f ),
                new Vector2( 0.561122477f,  -2.394122601f  ),
                new Vector2( 0.56188339f,   -8.094342232f   ),
                new Vector2( 0.562643945f, -13.79456139f  ),
                new Vector2(-6.010787487f, -13.79456139f ),
                new Vector2(-9.057126045f, -13.79471302f ),
                new Vector2(-9.05653286f,   -9.563955307f  ),
                new Vector2(-6.010194302f,  -9.563802719f ),
                new Vector2(-6.010483265f,  -8.186531067f ),
                new Vector2(-10.22067451f,  -8.238635063f ),
                new Vector2(-10.22199726f,  -2.395780802f ),
                //new Vector2(-6.012309074f,  -2.394122601f ),
            };

            vs = Reverse(vs);
            return TestPolygon(vs);
        }

        public static bool TestIssue4b_IndexOutOfRange()
        {
            Vector2[] vs = new Vector2[] {
                new Vector2(7.525075341f, 18.03584103f),
                new Vector2(8.873593212f, 18.03584103f),
                new Vector2(14.71882603f, 18.04310901f),
                new Vector2(20.09109498f, 18.04807795f),
                new Vector2(20.09109579f, 11.30548939f),
                new Vector2(20.08379922f, 7.390946768f),
                new Vector2(20.07960213f, 4.904943037f),
                new Vector2(14.8398342f, 4.913811439f),
                new Vector2(14.81097304f, 4.913967396f),
                new Vector2(14.81516979f, 7.399971046f),
                new Vector2(14.77711302f, 11.31475831f),
                new Vector2(8.897990992f, 11.31980019f),
                new Vector2(7.517957187f, 11.32239571f),
                new Vector2(7.517823608f, 11.16280541f),
                new Vector2(7.512232756f, 6.440598359f),
                new Vector2(2.821927548f, 6.425158617f),
                new Vector2(2.781571511f, 9.5992437f),
                new Vector2(2.761771935f, 11.32380416f),
                new Vector2(-2.768077453f, 11.33181232f),
                new Vector2(-2.756213328f, 18.02582882f),
                new Vector2(2.708230927f, 18.03272605f),
                new Vector2(2.707296381f, 19.4777793f),
                new Vector2(5.773485434f, 19.47942581f),
                new Vector2(5.774419978f, 18.03437257f),
                //new Vector2(7.525075341f, 18.03584103f),
            };

            vs = Reverse(vs);
            return TestPolygon(vs);
        }

    }
}
