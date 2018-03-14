using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace Assets.SpatialPartition.SmoothingKernel
{
    public class SmoothingKernel
    {
        //default kernel
        public static float Poly6(Vector3 ptlDistance, float krnRadius)
        {
            if (ptlDistance.magnitude > krnRadius)
            {
                //Debug.Log("Out of range");
                return 0;
            }
            else
            {
                return (
                        (315 / (64 * Mathf.PI * Mathf.Pow(krnRadius, 9.0F)))
                        * Mathf.Pow(
                            (krnRadius * krnRadius - ptlDistance.magnitude * ptlDistance.magnitude),
                            3.0F
                          )
                        );
            }
            
        }

        //gradient of poly6
        public static Vector3 GradientPoly6(Vector3 ptlDistance, float krnRadius)
        {
            if (ptlDistance.magnitude > krnRadius)
            {
                //Debug.Log("Out of range");
                return (new Vector3(0, 0, 0));
            }
            else
            {
                return (
                    ((-945) / (32 * Mathf.PI * Mathf.Pow(krnRadius, 9.0F)))
                    * ptlDistance
                    * Mathf.Pow(
                        (krnRadius * krnRadius - ptlDistance.magnitude * ptlDistance.magnitude),
                        2.0F
                        )
                   );
            }
        }

        //laplacian of poly6
        public static float LaplacianPloy6(Vector3 ptlDistance, float krnRadius)
        {
            if (ptlDistance.magnitude > krnRadius)
            {
                //Debug.Log("Out of range");
                return 0;
            }
            else
            {
                return (
                    ((-945) / (32 * Mathf.PI * Mathf.Pow(krnRadius, 9.0F)))
                    * (krnRadius * krnRadius - ptlDistance.magnitude * ptlDistance.magnitude)
                    * (3 * krnRadius * krnRadius - 7 * ptlDistance.magnitude * ptlDistance.magnitude)
                   );
            }
        }

        //smoothing kernel for pressure
        public static float Spiky(Vector3 ptlDistance, float krnRadius)
        {
            if (ptlDistance.magnitude > krnRadius)
            {
                //Debug.Log("Out of range");
                return 0;
            }
            else
            {
                return (
                        (15 / (Mathf.PI * Mathf.Pow(krnRadius, 6.0F)))
                        * Mathf.Pow(
                            (krnRadius - ptlDistance.magnitude),
                            3.0F
                          )
                        );
            }
        }

        public static Vector3 GradientSpiky(Vector3 ptlDistance, float krnRadius)
        {
            if (ptlDistance.magnitude > krnRadius)
            {
                //Debug.Log("Out of range");
                return (new Vector3(0,0,0));
            }
            else
            {
                return (
                    ((-45) / (Mathf.PI * Mathf.Pow(krnRadius, 6.0F)))
                    * (ptlDistance.normalized)
                    * Mathf.Pow(
                        (krnRadius - ptlDistance.magnitude),
                        2.0F
                        )
                    );
            }
        }

        public static float LaplacianSpiky(Vector3 ptlDistance, float krnRadius)
        {
            if (ptlDistance.magnitude > krnRadius)
            {
                //Debug.Log("Out of range");
                return 0;
            }
            else
            {
                return (
                    ((-90) / (Mathf.PI * Mathf.Pow(krnRadius, 6.0F)))
                    * (1 / ptlDistance.magnitude)
                    * ((krnRadius - ptlDistance.magnitude))
                    * ((krnRadius - 2 * ptlDistance.magnitude))
                    );
            }
        }

        //smoothing kernel for viscosity
        public static float Viscosity(Vector3 ptlDistance, float krnRadius)
        {
            if (ptlDistance.magnitude > krnRadius)
            {
                //Debug.Log("Out of range");
                return 0;
            }
            else
            {
                return (
                        (15 / (2 * Mathf.PI * Mathf.Pow(krnRadius, 3.0F)))
                        * (
                            ((-1 * Mathf.Pow(ptlDistance.magnitude, 3.0F)) / (2 * Mathf.Pow(krnRadius, 3.0F))) +
                            (Mathf.Pow(ptlDistance.magnitude, 2.0F) / Mathf.Pow(krnRadius, 2.0F)) +
                            (krnRadius / (2 * ptlDistance.magnitude)) - 1
                          )
                        );
            }
        }

        public static Vector3 GradientViscosity(Vector3 ptlDistance, float krnRadius)
        {
            if (ptlDistance.magnitude > krnRadius)
            {
                //Debug.Log("Out of range");
                return (new Vector3(0, 0, 0));
            }
            else
            {
                return (
                    (15 / (2 * Mathf.PI * Mathf.Pow(krnRadius, 3.0F)))
                    * ptlDistance
                    * (
                        ((-3 * ptlDistance.magnitude) / (2 * Mathf.Pow(krnRadius, 3.0F))) +
                        (2 / Mathf.Pow(krnRadius, 2.0F)) -
                        (krnRadius / (2 * Mathf.Pow(ptlDistance.magnitude, 3.0F)))
                        )
                    );
            }
        }

        public static float LaplacianViscosity(Vector3 ptlDistance, float krnRadius)
        {
            if (ptlDistance.magnitude > krnRadius)
            {
                //Debug.Log("Out of range");
                return 0;
            }
            else
            {
                return (
                    (45 / (Mathf.PI * Mathf.Pow(krnRadius, 6.0F)))
                    * (krnRadius - ptlDistance.magnitude)
                    );
            }
        }
    }
}

