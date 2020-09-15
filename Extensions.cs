using UnityEngine;

namespace Evryway
{

    // utility functions, extending Vector2.

    public static class Vector2Extensions
    {
        // rotate a Vector2. positive is Counter-clockwise.
        public static Vector2 Rotate(this Vector2 v, float angle_in_degrees)
        {
            float rad = angle_in_degrees * Mathf.Deg2Rad;
            float cosa = Mathf.Cos(rad);
            float sina = Mathf.Sin(rad);
            var rx = (v.x * cosa) - (v.y * sina);
            var ry = (v.x * sina) + (v.y * cosa);
            return new Vector2(rx, ry);
        }

        // rotate a Vector2 array. positive is Counter-clockwise.
        public static void Rotate(this Vector2[] array, float angle_in_degrees)
        {
            float rad = angle_in_degrees * Mathf.Deg2Rad;
            float cosa = Mathf.Cos(rad);
            float sina = Mathf.Sin(rad);

            for (int i = 0; i < array.Length; i++)
            {
                var v = array[i];
                var rx = (v.x * cosa) - (v.y * sina);
                var ry = (v.x * sina) + (v.y * cosa);
                array[i] = new Vector2(rx, ry);
            }
        }



        public static string F(this Vector2 source, int p) => source.ToString($"F{p}");
        public static string F3(this Vector2 source) => source.ToString("F3");
    }
}
