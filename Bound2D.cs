using System.Collections;
using System.Collections.Generic;
using UnityEngine;

// 2D bounds.
// NOT XY axis aligned.

namespace Evryway
{

    public class Bound2D
    {
        public Vector2 centre { get; private set; }
        public Vector2 axis_a { get; private set; }     // first axis (unit length)
        public Vector2 axis_b { get; private set; }     // second axis (unit length)
        public float length_a { get; private set; }     // first axis length
        public float length_b { get; private set; }     // second axis length (may be larger)
        public Vector2 size { get; private set; }       // (length_a,length_b)
        public Vector2 extents { get => size * 0.5f; }    // half-length, e.g. from centre.

        public Vector2[] corners { get; private set; }  // four corners. in order: BL, BR, TR, TL (assuming axis a is X axis and axis b is Y axis)
                                                        // (- axis_a * extents.x - axis_b * extents.y), 
                                                        // (+ axis_a * extents.x - axis_b * extents.y), 
                                                        // (+ axis_a * extents.x + axis_b * extents.y), 
                                                        // (- axis_a * extents.x + axis_b * extents.y), 
        public float area { get; private set; }
        public bool major_axis_is_a { get => length_a >= length_b; }
        public float angle { get; private set; }

        public Vector2 bl { get => corners[0]; }
        public Vector2 br { get => corners[1]; }
        public Vector2 tr { get => corners[2]; }
        public Vector2 tl { get => corners[3]; }

        public Bound2D(Vector2 centre, Vector2 axis, Vector2 size)
        {
            this.centre = centre;
            this.axis_a = axis.normalized;
            this.axis_b = new Vector2(-axis_a.y, axis_a.x);
            this.size = size;
            this.length_a = size.x;
            this.length_b = size.y;
            Cache();

        }

        // make the bound aligned such that the major axis (longest length) is axis_a, and
        // axis_a.x is positive.
        
        public void AlignMajor()
        {
            bool dirty = false;
            if (!major_axis_is_a)
            {
                // rotate the axes.
                var axis_t = axis_a;
                axis_a = axis_b;
                axis_b = -axis_t;
                // swap lengths.
                length_a = size.y;
                length_b = size.x;
                size = new Vector2(length_a, length_b);
                dirty = true;
            }

            if (axis_a.x <= 0)
            {
                axis_a = -axis_a;
                axis_b = -axis_b;
                dirty = true;
            }

            if (dirty)
            {
                Cache();
            }

        }

        void Cache()
        {
            area = length_a * length_b;
            angle = Vector2.SignedAngle(Vector2.right, axis_a);
            corners = new Vector2[]
            {
                centre - (axis_a * extents.x) - (axis_b * extents.y),
                centre + (axis_a * extents.x) - (axis_b * extents.y),
                centre + (axis_a * extents.x) + (axis_b * extents.y),
                centre - (axis_a * extents.x) + (axis_b * extents.y),
            };
        }

        public string info { get => $"area: {area}, axes: ({axis_a.x}, {axis_a.y}) ; ({axis_b.x}, {axis_b.y}), size: ({length_a}, {length_b}) angle: {angle}"; }

    }
}
