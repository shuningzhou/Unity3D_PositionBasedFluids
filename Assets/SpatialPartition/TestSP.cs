using System.Collections;
using System.Collections.Generic;
using UnityEngine;

using Assets.SpatialPartition.Utils;
using Assets.SpatialPartition;
using Assets.SpatialPartition.Utils.FMaterial;
using Assets.SpatialPartition.SmoothingKernel;

using System;

public class TestSP : MonoBehaviour
{
#if false
    SpatialHashTable ht;
    PhysicsSolver phySolver;

    float fTimeStep;

    //MyParticle p1;

    //MySphereCollider Container;

    MyCapsuleCollider Container;

    MyParticle[] ps;
    Vector3 offset;
    float step;
    int counterY;
    int counterZ;
    bool locker;

    Material newMat;

    int num;

    //contains result of nearest neighbour search
    List<MyParticle> Test_ParticleQueryReasult;

    float fptlRadius;
    float sphereScale;

    bool firstTime;

    void SetParticleMass()
    {
        float m = (Water.PARTICLE_REST_DENSITY * Consts.FLUID_VOLUME) / Consts.NUM_OF_PARTICLES; //0.02f;

        //(Water.PARTICLE_REST_DENSITY * Consts.FLUID_VOLUME) / Consts.NUM_OF_PARTICLES;
        Debug.Log("mass: " + m);

        Water.SetPMass(m);
    }

    void SetSKH()
    {
        float h = Mathf.Pow(((3 * Consts.FLUID_VOLUME * Consts.KERNEL_PARTICLE) / (4 * Mathf.PI * Consts.NUM_OF_PARTICLES)), 1.0f / 3.0f); //0.0457

        //Mathf.Pow(((3 * Consts.FLUID_VOLUME * Consts.KERNEL_PARTICLE) / (4 * Mathf.PI * Consts.NUM_OF_PARTICLES)), 1.0f / 3.0f);

        Debug.Log("kernel radius: " + h);

        Water.SetH(h);
    }

    void SetParticleSize()
    {
        float cbrtNum = (Water.GetPMass() * 4) / (3 * Mathf.PI * Water.PARTICLE_REST_DENSITY);
        fptlRadius = (float)Math.Pow(cbrtNum, (1.0 / 3.0));
        Debug.Log("particle radius:  " + fptlRadius);

        sphereScale = fptlRadius * 2;

        //Debug.Log("Kernel h:  " + Mathf.Pow(((3 * 0.1f * 20) / (4 * Mathf.PI * 64)), 1.0f / 3.0f));
    }

    void SetStep()
    {
        num = (int)Mathf.Pow(Consts.NUM_OF_PARTICLES, 1.0f / 3.0f);
        step = fptlRadius * 2;//(Mathf.Pow((Water.GetPMass() * Consts.NUM_OF_PARTICLES / Water.PARTICLE_REST_DENSITY), 1.0f / 3.0f) / num);

        Debug.Log("step: " + step + "  num: " + num);
    }

    // Use this for initialization
    void Start()
    {
        firstTime = true;


        SetParticleMass();
        SetSKH();

        /*
        Container = new MySphereCollider(new Vector3(10f, 10f, 10f), 1f);
        Container.Model = GameObject.CreatePrimitive(PrimitiveType.Sphere);
        Container.Model.transform.localScale = new Vector3(2f, 2f, 2f);
        Container.Model.transform.position = new Vector3(10f, 10f, 10f);
        newMat = Resources.Load("MyMat2", typeof(Material)) as Material;
        Container.Model.GetComponent<Renderer>().material = newMat;
        */

        Container = new MyCapsuleCollider(new Vector3(10f, 9.5f, 10f), new Vector3(10f, 10.5f, 10f), 1f);

        SetParticleSize();

        SetStep();

        fTimeStep = 0.02f;

        offset = new Vector3(0, 0, 0);

        //step = fptlRadius * 2.0f;
        locker = true;


        ps = new MyParticle[Consts.NUM_OF_PARTICLES];
        ht = new SpatialHashTable(Consts.NUM_OF_PARTICLES, Water.GetH());
        phySolver = new PhysicsSolver();

        for (int i = 0; i < ps.Length; i++)
        {
            ps[i] = new MyParticle();


            //initialization
            ps[i].fMass = Water.GetPMass();
            //Debug.Log("Initial mass = " + ps[i].fMass);
            ps[i].v3Velocity = new Vector3(0, 0, 0);
            ps[i].fMassDensity = Water.PARTICLE_REST_DENSITY;
            ps[i].fPressure = 0;
            ps[i].h_v3Velocity = ps[i].v3Velocity;

            ps[i].id = i;

            ps[i].InternalForce = new Vector3(0, 0, 0);
            ps[i].ExternalForce = new Vector3(0, 0, 0);

            /*
            ps[i].v3Position = this.transform.position + offset;
            offset.y += step;

            ps[i].model = GameObject.CreatePrimitive(PrimitiveType.Sphere);
            ps[i].model.transform.position = ps[i].v3Position;

            ps[i].model.transform.localScale = new Vector3(sphereScale, sphereScale, sphereScale);
            */

            
            ps[i].v3Position = this.transform.position + offset;
            offset.y += step;
            counterY++;


            //assigning model
            ps[i].model = GameObject.CreatePrimitive(PrimitiveType.Sphere);
            ps[i].model.transform.position = ps[i].v3Position;

            ps[i].model.transform.localScale = new Vector3(sphereScale, sphereScale, sphereScale);

            if (counterY == num)
            {
                offset.y = 0;
                offset.x += step;
                counterY = 0;

                counterZ++;
            }

            if (counterZ == num)
            {
                offset.x = 0;
                offset.y = 0;
                offset.z += step;

                counterZ = 0;
            }
            

            ht.AddParticle(ps[i]);
        }

        //Debug.Log("Kernel: " + SmoothingKernel.Poly6(new Vector3(0f, 0f, 0f), 1f));
        //Debug.Log("Magnitude: " + (new Vector3(0f, 0f, 0f)).magnitude);

        Debug.Log("Kernel : " + SmoothingKernel.Poly6(new Vector3(0f, 0f, 0f), Water.GetH()));
        Debug.Log("Kernel : " + SmoothingKernel.Poly6(new Vector3(0, 0.228f, 0), Water.GetH()));
        Debug.Log("Kernel : " + SmoothingKernel.Poly6(new Vector3(0.1f, 0.1f, 0.1f), Water.GetH()));
        //Debug.Log(SmoothingKernel.Spiky(new Vector3(0, 0, 0), 0.13f));
        //Debug.Log(SmoothingKernel.GradientSpiky(new Vector3(0.0001f, 0.0001f, 0.0001f), 0.13f).magnitude);
        //Debug.Log(SmoothingKernel.LaplacianSpiky(new Vector3(0, 0, 0), 0.13f));

        newMat = Resources.Load("MyMat", typeof(Material)) as Material;
        ps[0].model.GetComponent<Renderer>().material = newMat;

        /*
        List<MyParticle> neighbours = ht.NeighboursQuery(ps[43], false);

        foreach (MyParticle p in neighbours)
        {
            ps[p.id].model.transform.position = new Vector3(ps[p.id].model.transform.position.x, ps[p.id].model.transform.position.y, ps[p.id].model.transform.position.z - 1.5f);
        }
        Debug.Log("nei count: " + neighbours.Count);
        */
    }

#if true


    private void Update()
    {
        if (Input.GetKeyUp(KeyCode.Space))
        {
            locker = !locker;
            Debug.Log("Locker: " + locker);
        }

        if (Input.GetKeyUp(KeyCode.P))
        {
            List<MyParticle> neighbours = ht.NeighboursQuery(ps[0], false);
            foreach (MyParticle p in neighbours)
            {
                p.model.GetComponent<Renderer>().material = newMat;
            }
        }
    }



    private void FixedUpdate()
    {

        if (locker == false)
        {
            //Density and pressure
            foreach (MyParticle p in ps)
            {

                //Vector3 oldPos = p.v3Position; //store previous position

                //compute density and pressure
                List<MyParticle> neighbours = ht.NeighboursQuery(p, false);
                phySolver.SetMassDensity(p, neighbours);
                phySolver.SetPressure(p, neighbours);

                //Debug.Log("ID: " + p.id + "  Pressure: " + p.fPressure + "   " + "Density: " + p.fMassDensity + "  Num of n: " + neighbours.Count);

                if (p.id == 0)
                {
                    //Debug.Log("ID: " + p.id + "  Pressure: " + p.fPressure + "   " + "Density: " + p.fMassDensity + "  Num of n: " + neighbours.Count);
                }
            }

#if true



            //Internal force
            foreach (MyParticle p in ps)
            {
                List<MyParticle> neighboursF = ht.NeighboursQuery(p, true);
                Vector3 PressureForce = phySolver.PressureForce(p, neighboursF);
                Vector3 ViscosityForce = phySolver.ViscosityForce(p, neighboursF);

                p.InternalForce = PressureForce;//ViscosityForce;

                //PressureForce + ViscosityForce

                //Debug.Log("ID: " + p.id + " Pressure force: " + PressureForce + " Count: " + neighboursF.Count);

                if (p.id == 0)
                {
                    //Debug.Log("Pressure force: " + PressureForce);

                    //Debug.Log("Pressure Force: " + PressureForce + "Viscosity force: " + ViscosityForce);
                }
            }

            //External force
            foreach (MyParticle p in ps)
            {
                List<MyParticle> neighboursF = ht.NeighboursQuery(p, true);
                Vector3 Gravity = phySolver.Gravity(p);
                Vector3 SurfaceTension = phySolver.SurfaceTension(p, neighboursF);

                p.ExternalForce = Gravity; //Gravity + SurfaceTension;
            }

            //Translation & time integrator
            foreach (MyParticle p in ps)
            {
                Vector3 oldPos = p.v3Position; //store previous position
                Vector3 oldHVel = p.h_v3Velocity;

                Vector3 Force = p.InternalForce + p.ExternalForce;

                    //p.InternalForce + p.ExternalForce;

                //Vector3 acc = p.InternalForce / 1000 + (p.ExternalForce / p.fMassDensity);

                Vector3 acc = Force / p.fMassDensity;

                if(p.id == 0)
                //Debug.Log("Force: " + Force);

                if (firstTime)
                {
                    oldHVel = p.v3Velocity - (0.5f * fTimeStep * acc);
                }

                p.h_v3Velocity = oldHVel + fTimeStep * acc;

                p.v3Position = p.v3Position + fTimeStep * p.h_v3Velocity;
                p.model.transform.position = p.v3Position;

                //Debug.Log("Position: " + p.v3Position + "  AccI: " + p.InternalForce / p.fMassDensity + "  Force: " + Force.magnitude + "  Velocity: " + p.v3Velocity);

                /*
                p.v3Velocity = p.v3Velocity + fTimeStep * acc;

                p.v3Position = p.v3Position + fTimeStep * p.v3Velocity;
                p.model.transform.position = p.v3Position;
                */



                /*
                Vector3 repulsiveOffset = Vector3.zero;

                List<MyParticle> neighboursF = ht.NeighboursQuery(p, true);
                // push away particles that are too close
                foreach (MyParticle pt in neighboursF)
                {
                    Vector3 diff = pt.v3Position - p.v3Position;
                    float md2 = 0.0168f;// * 0.0168f;
                    //float mD2 = maxDist * maxDist;
                    float dot = diff.sqrMagnitude;
                    if (dot < md2)
                    {
                        repulsiveOffset += 100f * diff.normalized * (md2 - dot);
                    }

                }

                p.v3Position -= repulsiveOffset;
                p.model.transform.position = p.v3Position;
                */

                /*
                int result = Container.CheckParticle(p.v3Position);
                if (result != -1)
                {
                    p.v3Position = Container.v3ContactPoint;
                    p.model.transform.position = p.v3Position;

                    p.h_v3Velocity = p.h_v3Velocity - ((Vector3.Dot(p.v3Velocity, Container.v3SurfaceNormal)) * Container.v3SurfaceNormal);
                }
                */

                int result = Container.CollisionCheck(p.v3Position);
                if (result != -1)
                {
                    p.v3Position = Container.contactPoint;
                    p.model.transform.position = p.v3Position;

                    p.h_v3Velocity = p.h_v3Velocity - ((Vector3.Dot(p.v3Velocity, Container.surfaceNormal)) * Container.surfaceNormal);
                }

                p.v3Velocity = (p.h_v3Velocity + oldHVel) / 2;

                ht.UpdateParticle(p, oldPos);
            }

            if (firstTime == true)
            {
                firstTime = false;
            }
#endif
        }

    }
#endif
#endif
}

