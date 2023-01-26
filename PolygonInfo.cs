using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;

namespace Evryway
{

    public class PolygonInfo 
    {
        public bool valid;
        public bool convex;
        public bool clockwise;
        public double area;

        public string info { get => $"valid:{valid} convex:{convex} clockwise:{clockwise} area:{area}"; }

        const float twopi = Mathf.PI * 2.0f;
        const float ang_delta_epsilon = 0.001f;
        static Vector2[] edges = new Vector2[65536];
        public static PolygonInfo invalid = new PolygonInfo { valid = false, convex = false, clockwise = false, area = 0 };

        // points are all on the XY plane.
        public static PolygonInfo Create(Vector2[] points)
        {
            if (points.Length < 3) return invalid;
            var area = PolygonArea(points);
            if (area == 0) return invalid;

            int pc = points.Length;
            if (pc > edges.Length) edges = new Vector2[pc * 2];
            for (int i = 0; i < pc; i++)
            {
                var ip1 = (i + 1) % pc;
                var pi = points[i];
                var pip1 = points[ip1];
                var edge = new Vector2(pip1.x - pi.x, pip1.y - pi.y);
                if (math.length(edge) <= 0.0f)
                {
                    return invalid;
                }
                edges[i] = edge;
            }
            return CreateFromEdges(edges, pc, area);
        }

        public static PolygonInfo CreateFromEdges(Vector2[] edges, int count, double area)
        {
            bool convex = true;

            float ang_dir_expected = 0.0f;
            double ang_delta_total = 0.0;
            for (int i = 0; i < count; i++)
            {
                var ip1 = (i + 1) % count;
                var e1 = edges[i];
                var e2 = edges[ip1];
                var ang1 = math.atan2(e1.y, e1.x);
                var ang2 = math.atan2(e2.y, e2.x);
                var ang_delta = ang2 - ang1;
                if (ang_delta < -Mathf.PI) ang_delta += (twopi);
                else if (ang_delta > Mathf.PI) ang_delta -= (twopi);

                // clockwise goes negative (see Atan2 deets above). anticlockwise goes positive.
                var ang_dir = ang_delta >= 0.0f ? 1.0f : -1.0f;
                if (i == 0)
                {
                    // first edge, set expected direction. -1 ? clockwise. 1 ? anticlockwise.
                    ang_dir_expected = ang_dir;
                }
                else
                {
                    // if our expected direction varies, we're not convex.
                    if (ang_dir_expected != ang_dir)
                    {
                        convex = false;
                    }
                }
                ang_delta_total += ang_delta;
            }

            // ang_delta_total should be either -2*PI or 2*PI (within error bounds)
            // if it's not, we've got a poly that's self-winding more than once (e.g. star, bowtie, that sort of thing)
            // it's probably self-intersecting in that case ...
            if (convex)
            {
                var c = (ang_delta_total * ang_dir_expected);
                if (math.abs(c - twopi) > ang_delta_epsilon)
                {
                    convex = false;
                }
            }
            // if ang_delta_total is negative, we've gone clockwise around.
            bool clockwise = ang_delta_total < 0.0f ? true : false;
            bool valid = count > 4 || ang_delta_total != 0;
            
            return new PolygonInfo { valid = valid, convex = convex, clockwise = clockwise, area = area };
        }

        public static double PolygonArea(Vector2[] points)
        {
            var pc = points.Length;
            double sum = 0;
            for (int i = 0 ; i < pc ; i++)
            {
                var p = points[i];
                var j = (i+1) % pc;
                var q = points[j];
                sum += (q.x-p.x) * (q.y+p.y);
            }
            sum /= 2;
            return sum;
        }


    }
}
