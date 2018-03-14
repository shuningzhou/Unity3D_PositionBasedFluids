using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using UnityEngine;

namespace Assets.SpatialPartition.Utils
{
    namespace FMaterial
    {
        public class Water
        {
            //smoothing kernel radius
            private static float fSKh;

            public static float GetH()
            {
                return fSKh;
            }
            public static void SetH(float h)
            {
                fSKh = h;
            }

            private static float fPtlMass;

            public static float GetPMass()
            {
                return fPtlMass;
            }
            public static void SetPMass(float mass)
            {
                fPtlMass = mass;
            }

            public const float PARTICLE_REST_DENSITY = 998.29F;

            public const float GAS_STIFFNESS = 3.0F;
            public const float VISCOSITY_COEFFICIENT = 3.5f;//3.5f;

            public const float SURFACE_TENSION_THRESHOLD = 7.065f;
            public const float SURFACE_TENSION_COEFFICIENT = 0.0728f;
        }

    }

    public class Consts
    {
        //Used in hash function, large primes
        public const int iHashP1 = 73856093;
        public const int iHashP2 = 19349663;
        public const int iHashP3 = 83492791;

        public const int NUM_OF_PARTICLES = 100;
        public const float FLUID_VOLUME = 0.1f;

        public const int KERNEL_PARTICLE = 100;
    }

    public class MyParticle
    {
        public int id;

        public float   fMass;
        public Vector3 v3Position;

        public Vector3 v3Velocity;
        public Vector3 h_v3Velocity;

        public GameObject model;

        public float fPressure;
        public float fMassDensity;

        public Vector3 InternalForce;
        public Vector3 ExternalForce;

        public Vector3 SurfaceNormal;
        public float LaplacianColorField;

        public List<int> neighboursID;
    }

    public class Functions
    {
        public static int sgn(float val)
        {
            if (val == 0)
                return 0;
            if (val > 0)
                return 1;
            else
                return -1;
        }
    }
}
