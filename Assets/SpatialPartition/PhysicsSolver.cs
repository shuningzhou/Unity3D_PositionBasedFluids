using System.Collections;
using System.Collections.Generic;
using UnityEngine;

using Assets.SpatialPartition.SmoothingKernel;
using Assets.SpatialPartition.Utils;
using Assets.SpatialPartition.Utils.FMaterial;

public class PhysicsSolver
{
    public void SetMassDensity(MyParticle ptl, List<MyParticle> neighbours)
    {
        float result = 0;
        foreach (MyParticle p in neighbours)
        {

            float temp = (p.fMass * SmoothingKernel.Poly6(ptl.v3Position - p.v3Position, Water.GetH()));

            result = result + temp;

            /*
            if (neighbours.Count > 18)
            {
                Debug.Log("Count: " + neighbours.Count);
                Debug.Log("ID: " + ptl.id + "  result: " + result + " temp: " + temp + " dist:" + (ptl.v3Position - p.v3Position).magnitude);
            }
            */

            //result = result + (p.fMass * SmoothingKernel.Poly6(ptl.v3Position - p.v3Position, Water.GetH()));
            //Debug.Log("kernel result: " + SmoothingKernel.Poly6(ptl.v3Position - p.v3Position, Water.GetH()));

            //Debug.Log("kernel result: " + result);
        }

        ptl.fMassDensity = result;
    }

    public void SetPressure(MyParticle ptl, List<MyParticle> neighbours)
    {
        //ptl.fPressure = (Mathf.Pow(ptl.fMassDensity / Water.PARTICLE_REST_DENSITY, 7) - 1) * 1;

        ptl.fPressure =  Water.GAS_STIFFNESS * (ptl.fMassDensity - Water.PARTICLE_REST_DENSITY);
    }

    public Vector3 PressureForce(MyParticle ptl, List<MyParticle> neighbours)
    {
        Vector3 result = new Vector3(0,0,0);

        //Excluding target particle is required
        foreach (MyParticle p in neighbours)
        {
            //Debug.Log("Particle location: " + ptl.v3Position + "  distance: " + (ptl.v3Position - p.v3Position));

            /*
            Vector3 pressureF = (
                ((ptl.fPressure + p.fPressure) / 2) *
                (p.fMass / p.fMassDensity) *
                SmoothingKernel.GradientSpiky(ptl.v3Position - p.v3Position, Water.GetH()));
                */

            Vector3 pressureF = (
                ((ptl.fPressure + p.fPressure) / 2) *
                (p.fMass / p.fMassDensity) *
                SmoothingKernel.GradientSpiky((ptl.v3Position - p.v3Position), Water.GetH()));

            //result = result + pressureF;

            result = result - pressureF;

            /*
            //Debug.Log("force = " + pressureF + "   result = " + result);
            if (ptl.id == 0)
            {
                Debug.Log("part3: " + (-45) / (Mathf.PI * Mathf.Pow(Water.GetH(), 6.0F)));

                Debug.Log("part2: " + (ptl.v3Position - p.v3Position).normalized);

                Debug.Log("part: " + Mathf.Pow(
                        (Water.GetH() - (ptl.v3Position - p.v3Position).magnitude),
                        2.0F
                        ));

                Debug.Log("kernel result: " + SmoothingKernel.GradientSpiky(ptl.v3Position - p.v3Position, Water.GetH()) + " meanP: " + ((ptl.fPressure + p.fPressure) / 2) + " m/d: " + (p.fMass / p.fMassDensity));

                Debug.Log("result: " + result + "  pressure: " + pressureF);
            }
            */
        }

        /*
        Debug.Log("++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        Debug.Log("Final result: " + result);
        Debug.Log("++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        */

        //result = result * -1;

        //Debug.Log("Final result  =   " + result + "   count = " + neighbours.Count);

        return result;
    }

    public Vector3 ViscosityForce(MyParticle ptl, List<MyParticle> neighbours)
    {
        Vector3 result = new Vector3(0, 0, 0);

        //Excluding target particle is required
        foreach (MyParticle p in neighbours)
        {
            Vector3 temp = (
                (p.v3Velocity - ptl.v3Velocity) *
                (p.fMass / p.fMassDensity) *
                SmoothingKernel.LaplacianViscosity(ptl.v3Position - p.v3Position, Water.GetH()));

            result = result + temp;

            if (ptl.id == 10)
            {
                //Debug.Log("ID: " + ptl.id + "  result: " + result + " temp: " + temp + " dist:" + SmoothingKernel.LaplacianViscosity(ptl.v3Position - p.v3Position, Water.GetH()));
            }

            
        }

        //if (ptl.id == 10)
            //Debug.Log("+++++++++++++++++++++++++++++++++++++++++++++");

        result = result * Water.VISCOSITY_COEFFICIENT;
        return result;
    }

    public Vector3 InternalForce(MyParticle ptl, List<MyParticle> neighbours)
    {
        return (PressureForce(ptl, neighbours) + ViscosityForce(ptl, neighbours));
    }

    public Vector3 Gravity(MyParticle ptl)
    {
        return ((new Vector3(0, -9.82f, 0)) * ptl.fMassDensity);  //vector -> downward gravitional acceleration
    }

    public Vector3 SurfaceTension(MyParticle ptl, List<MyParticle> neighbours)
    {
        float fColorField = 0;
        Vector3 SurfaceNormal = new Vector3(0, 0, 0);
        float SurfaceNormalD = 0;

        foreach (MyParticle p in neighbours)
        {
            float temp1 = p.fMass / p.fMassDensity;
            Vector3 temp2 = ptl.v3Position - p.v3Position;

            fColorField = fColorField +
                (temp1 * SmoothingKernel.Poly6(temp2, Water.GetH()));

            SurfaceNormal = SurfaceNormal +
                (temp1 * SmoothingKernel.GradientPoly6(temp2, Water.GetH()));

            SurfaceNormalD = SurfaceNormalD +
                (temp1 * SmoothingKernel.LaplacianPloy6(temp2, Water.GetH()));

        }

        Vector3 SurfaceTension = new Vector3(0, 0, 0);

        if (SurfaceNormal.magnitude >= Water.SURFACE_TENSION_THRESHOLD)
        {
            Debug.Log("ID: " + ptl.id);

            SurfaceTension = (-1) * Water.SURFACE_TENSION_COEFFICIENT *
                 SurfaceNormalD *
                SurfaceNormal.normalized;
        }

        return SurfaceTension;
    }

    public Vector3 ExternalForce(MyParticle ptl, List<MyParticle> neighbours)
    {
        return (Gravity(ptl) + SurfaceTension(ptl, neighbours));
    }
}
