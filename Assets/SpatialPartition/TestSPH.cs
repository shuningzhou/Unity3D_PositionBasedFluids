#define CONTAINER_ON
//#undef CONTAINER_ON

using System.Collections;
using System.Collections.Generic;
using UnityEngine;

using Assets.SpatialPartition.Utils;
using Assets.SpatialPartition;
using Assets.SpatialPartition.Utils.FMaterial;

using System;

public class TestSPH : MonoBehaviour
{

    MyParticle[] ps;
    SpatialHashTable ht;

    //MyCapsuleCollider Container;
    MyBoxCollider Container;

    bool bActionLock = true;

    //consts
    //-----------------------------
    private float REST_DENSITY = 998.29f;//998.29F;
    private float MATERIAL_DENSITY = 998.29f;

    private float GAS_STIFFNESS = 25f;//15f;//18.0f;
    private float VISCOSITY_COEFFICIENT = 20f;//15f; //3.5f;

    private float SURFACE_TENSION_COEFFICIENT = 0.1456f;//0.0728f;

    private int TOTAL_PARTICLES = 300;
    private int PARTICLES_PER_VOLUME = 300;   //how many particles are there in a unit of volume of water e.g. 100 particles to form 1m3 water
    private float VOLUME = 1f;
    private int KERNEL_PARTICLE = 10;
    
    private float TIME_STEP = 0.02f;

    private Vector3 G = new Vector3(0, -9.8f, 0);
    //-----------------------------

    private float SURFACE_TENSION_THRESHOLD;

    private float fParticleMass;
    private float fSupportRadius;   //for smoothing kernel

    private float POLY6_C;
    private float GL_POLY6_C;  //gradient + laplacian
    private float G_PRESSURE_C;
    private float L_VISCOSITY_C;

    private float KRad;
    private float KRad2; //k ^ 2

    Material newMat;

    // Use this for initialization
    void Start()
    {
        //step 1
        //Initialising material

        //step 2
        //Initialising particles
        //--------------------------------------------------------------------------------------------------------------
        ps = new MyParticle[TOTAL_PARTICLES];
        fParticleMass = MATERIAL_DENSITY * VOLUME / PARTICLES_PER_VOLUME;
        Debug.Log("Initial mass = " + fParticleMass);

        /*
        Vector3 offset = new Vector3(0, 0, 0);
        float step = Mathf.Pow((fParticleMass * 3) / (4 * Mathf.PI * REST_DENSITY), (1.0f / 3.0f));
        */
        //set physical field
        for (int i = 0; i < ps.Length; i++)
        {
            ps[i] = new MyParticle();
            ps[i].id = i;

            ps[i].fMass = fParticleMass;

            ps[i].v3Velocity = Vector3.zero;
            ps[i].h_v3Velocity = new Vector3(Mathf.Infinity, Mathf.Infinity, Mathf.Infinity);

            ps[i].fMassDensity = REST_DENSITY;
            ps[i].fPressure = 0;

            ps[i].InternalForce = ps[i].ExternalForce = Vector3.zero;
            ps[i].SurfaceNormal = Vector3.zero;

            ps[i].neighboursID = null;
            /*
            ps[i].v3Position = this.transform.position + offset;
            offset.y += step;

            ps[i].model = GameObject.CreatePrimitive(PrimitiveType.Sphere);
            ps[i].model.transform.position = ps[i].v3Position;

            ps[i].model.transform.localScale = new Vector3(step, step, step);
            */
        }

        //set surface tension threshold
        SURFACE_TENSION_THRESHOLD = 0.5f;//Mathf.Sqrt(REST_DENSITY / KERNEL_PARTICLE);
        Debug.Log("Surface threshold: " + SURFACE_TENSION_THRESHOLD);

        //set position and size
        /*
        float l = Mathf.Pow(VOLUME, 1.0f / 3.0f);
        float a = Mathf.Pow(NUM_OF_PARTICLES, 1.0f / 3.0f);
        float step = l / a;
        */


        float step = Mathf.Pow((fParticleMass * 3) / (4 * Mathf.PI * MATERIAL_DENSITY), (1.0f / 3.0f)) * 2;
        int counter = 0;

        int a = (int)Mathf.Pow(PARTICLES_PER_VOLUME, 1.0f / 3.0f);
        float l = step * a;

        Vector3 center = this.transform.position;
        Vector3 bMin = center - (new Vector3(l / 2.0f, l / 2.0f, l / 2.0f));
        Vector3 bMax = center + (new Vector3(l / 2.0f, l / 2.0f, l / 2.0f));

        newMat = Resources.Load("GPUInstancing", typeof(Material)) as Material;

        for (float i = bMin.x; i <= bMax.x; i += step)
        {
            for (float j = bMin.y; j <= bMax.y; j += step)
            {
                for (float k = bMin.z; k <= bMax.z; k += step)
                {
                    if (counter < ps.Length)
                    {
                        ps[counter].v3Position = new Vector3(i, j, k);
                       
                        ps[counter].model = GameObject.CreatePrimitive(PrimitiveType.Sphere);
                        ps[counter].model.transform.position = ps[counter].v3Position;

                        ps[counter].model.transform.localScale = new Vector3(step, step, step);

                        ps[counter].model.GetComponent<Renderer>().material = newMat;

                        counter++;
                    }
                }
            }
        }
        
        for (int i = counter; i < ps.Length; i++)
        {
            ps[i].v3Position = new Vector3(bMax.x + step * (ps.Length - i - 1), bMax.y, bMax.z + step);

            ps[i].model = GameObject.CreatePrimitive(PrimitiveType.Sphere);
            ps[i].model.transform.position = ps[i].v3Position;

            ps[i].model.transform.localScale = new Vector3(step, step, step);

            ps[counter].model.GetComponent<Renderer>().material = newMat;
        }
        //--------------------------------------------------------------------------------------------------------------

        //step 3
        //Initialising kernel
        //--------------------------------------------------------------------------------------------------------------
        fSupportRadius = Mathf.Pow(((3 * VOLUME * KERNEL_PARTICLE) / (4 * Mathf.PI * PARTICLES_PER_VOLUME)), 1.0f / 3.0f);

        float k6 = Mathf.Pow(fSupportRadius, 6.0f);  //k ^ 9
        float k9 = Mathf.Pow(fSupportRadius, 9.0f);

        KRad = fSupportRadius;
        KRad2 = fSupportRadius * fSupportRadius;

        POLY6_C = (315 / (64 * Mathf.PI * k9));
        GL_POLY6_C = (-945) / (32 * Mathf.PI * k9);

        G_PRESSURE_C = (-45) / (Mathf.PI * k6);
        L_VISCOSITY_C = (45) / (Mathf.PI * k6);

        Debug.Log("fSupportRadius = " + fSupportRadius + "  step: " + step);
        //--------------------------------------------------------------------------------------------------------------


        //step 4
        //Initialising hash table
        //--------------------------------------------------------------------------------------------------------------
        ht = new SpatialHashTable(TOTAL_PARTICLES, fSupportRadius);

        foreach (MyParticle p in ps)
        {
            ht.AddParticle(p);
        }
        //--------------------------------------------------------------------------------------------------------------


        //step 5
        //Creating collision objects
        //--------------------------------------------------------------------------------------------------------------
#if (CONTAINER_ON)
        //Container = new MyCapsuleCollider(new Vector3(10f, 9f, 10f), new Vector3(10f, 10f, 10f), 0.4f);
        Container = new MyBoxCollider(new Vector3(9.9f, 9.9f, 10f), new Vector3(0.5f, 1.0f, 0.5f));
#endif
        //--------------------------------------------------------------------------------------------------------------


        //step 6 
        //Initialising the leap frog integrator 
    }


    //------------------------smoothing-kernel----------------------------------------
    public float Poly6(Vector3 Rij)
    {
        // ||r|| <= h is originally held by query machenism 

        float h2r2_3 = Mathf.Pow(KRad2 - Rij.magnitude * Rij.magnitude, 3.0f);

        //Debug.Log("c: " + POLY6_C + " h2r2: " + h2r2_3);

        return POLY6_C * h2r2_3;
    }

    public Vector3 GradientPoly6(Vector3 Rij)
    {
        float h2r2_2 = Mathf.Pow(KRad2 - Rij.magnitude * Rij.magnitude, 2.0f);

        return GL_POLY6_C * Rij * h2r2_2;
    }

    public float LaplacianPoly6(Vector3 Rij)
    {
        float h2r2 = KRad2 - Rij.magnitude * Rij.magnitude;
        float h3p2r7p2 = 3 * KRad2 - 7 * Rij.magnitude * Rij.magnitude;

        return GL_POLY6_C * h2r2 * h3p2r7p2;
    }

    public Vector3 GradientPressure(Vector3 Rij)
    {
        float hr2 = Mathf.Pow(KRad - Rij.magnitude, 2.0f);

        return G_PRESSURE_C * Rij.normalized * hr2;
    }

    public float LaplacianViscosity(Vector3 Rij)
    {
        float hr = KRad - Rij.magnitude;

        return L_VISCOSITY_C * hr;
    }

    //--------------------------------------------------------------------------------

#if true

    private void FixedUpdate()
    {
        if (Input.GetKeyUp(KeyCode.Space))
        {
            bActionLock = !bActionLock;
            Debug.Log("Locker: " + bActionLock);
        }

        if (Input.GetKeyUp(KeyCode.Z))
        {
            float avgDensity = 0;
            foreach(MyParticle p in ps)
            {
                avgDensity += p.fMassDensity;
            }
            Debug.Log("Average density: " + (avgDensity / ps.Length));
        }

#if (CONTAINER_ON)
        if (Input.GetKey(KeyCode.D))
        {
            Container.UpdatePosition(new Vector3(0.05f, 0,0));
        }
        if (Input.GetKey(KeyCode.A))
        {
            Container.UpdatePosition(new Vector3(-0.05f, 0, 0));
        }
        if (Input.GetKey(KeyCode.W))
        {
            Container.UpdatePosition(new Vector3(0, 0.05f, 0));
        }
        if (Input.GetKey(KeyCode.S))
        {
            Container.UpdatePosition(new Vector3(0, -0.05f, 0));
        }
#endif

        //animation procedure
        if (bActionLock == false)
        {
            //i.computing density & pressure
            //--------------------------------------------------------------------------------
            foreach (MyParticle p in ps)
            {
                p.neighboursID = ht.NeighboursQuery(p, false);
                //density
                p.fMassDensity = 0;
                foreach(int id in p.neighboursID)
                {
                    p.fMassDensity += (ps[id].fMass * Poly6(p.v3Position - ps[id].v3Position));

                    /*
                    if (p.id == 0)
                    {
                        Debug.Log("Density: " + p.fMassDensity);
                    }*/
                }


                //pressure
                p.fPressure = GAS_STIFFNESS * (p.fMassDensity - REST_DENSITY);

               
                //if (p.id == 0)
                     //Debug.Log("ID: " + p.id + " Pressure: " + p.fPressure+ " Density: " + p.fMassDensity + " Counted: " + p.neighboursID.Count);
            }
            
#if true

            //ii.computing internal & external force
            foreach (MyParticle p in ps)
            {
                //pressure & viscosity
                Vector3 pressureForce = Vector3.zero;
                Vector3 viscosityForce = Vector3.zero;

                foreach (int id in p.neighboursID)
                {
                    if (id == p.id)
                        continue;

                    float pipjmj = (p.fPressure + ps[id].fPressure) * ps[id].fMass;
                    float dj2 = 2 * ps[id].fMassDensity;

                    Vector3 distance = p.v3Position - ps[id].v3Position;

                    pressureForce -= ((pipjmj / dj2) * GradientPressure(distance));

                    //-------------------------------------------------------------------------------------------------------------

                    Vector3 diffVel = ps[id].v3Velocity - p.v3Velocity;
                    float md = ps[id].fMass / ps[id].fMassDensity;

                    viscosityForce += (diffVel * md * LaplacianViscosity(distance));
                }

                // Debug.Log("AP: " + accumulatorP);

                //Debug.Log("pressure: " + pressureForce + "  Acc: " + (pressureForce / p.fMassDensity));

                //if(p.id == 0)
                //Debug.Log("viscosity: " + viscosityForce + "  Acc: " + (viscosityForce / p.fMassDensity));

                viscosityForce *= VISCOSITY_COEFFICIENT;

                p.InternalForce = pressureForce + viscosityForce;
                //pressureForce + viscosityForce;
            }

            foreach (MyParticle p in ps)
            {
                //------------------------------------------------------------------------------------------
                //iii.computing external force

                Vector3 gravityForce = G * p.fMassDensity;
                Vector3 surfaceTension = Vector3.zero;

                //initialise quantities to 0
                p.SurfaceNormal = Vector3.zero;
                p.LaplacianColorField = 0.0f;

                //Surface tension
                //------------------------------------------------------------------------------------
                foreach (int id in p.neighboursID)
                {
                    //-------------------------------------------------------------------------------------------------------------
                    //surface nromal
                    float md = ps[id].fMass / ps[id].fMassDensity;
                    Vector3 distance = p.v3Position - ps[id].v3Position;

                    p.SurfaceNormal += (md * GradientPoly6(distance));
                    p.LaplacianColorField += (md * LaplacianPoly6(distance));
                }

                //Debug.Log("SURFACE_TENSION_THRESHOLD" + SURFACE_TENSION_THRESHOLD + "   " + p.SurfaceNormal.magnitude);
                if (p.SurfaceNormal.magnitude >= SURFACE_TENSION_THRESHOLD)
                {

                    //Debug.Log("aaaaaaaaaaaaaaa");

                    //visualise surface tracking result
                    //Material newMat = Resources.Load("MyMat2", typeof(Material)) as Material;
                    //p.model.GetComponent<Renderer>().material = newMat;

                    surfaceTension = (-SURFACE_TENSION_COEFFICIENT) * p.LaplacianColorField * p.SurfaceNormal.normalized;
                    //Debug.Log("surfaceTension: " + surfaceTension);
                }


                //------------------------------------------------------------------------------------
                p.ExternalForce = gravityForce + surfaceTension;//gravityForce;
            }
            


            //iv.time integration & collision handling
            foreach (MyParticle p in ps)
            {
                //
                Vector3 oldPos = p.v3Position;
                //


                Vector3 force = p.InternalForce + p.ExternalForce;
                //p.InternalForce + p.ExternalForce;

                Vector3 acc = force / p.fMassDensity;

                //leap-frog
                //--------------------------------------------------------------------------------
                Vector3 oldHVel;

                if (p.h_v3Velocity.Equals(new Vector3(Mathf.Infinity, Mathf.Infinity, Mathf.Infinity)))
                {
                    oldHVel = p.v3Velocity - (0.5f * TIME_STEP * acc);
                }
                else
                {
                    oldHVel = p.h_v3Velocity;
                }

                p.h_v3Velocity = oldHVel + TIME_STEP * acc;

                p.v3Position = p.v3Position + TIME_STEP * p.h_v3Velocity;
                p.model.transform.position = p.v3Position;
                //--------------------------------------------------------------------------------

#if (CONTAINER_ON)
                int result = Container.CollisionCheck(p.v3Position);
                if (result != -1)
                {
                    p.v3Position = Container.contactPoint;
                    p.model.transform.position = p.v3Position;

                    p.h_v3Velocity = p.h_v3Velocity - ((Vector3.Dot(p.v3Velocity, Container.surfaceNormal)) * Container.surfaceNormal);
                }
#endif

                p.v3Velocity = (p.h_v3Velocity + oldHVel) / 2;

                ht.UpdateParticle(p, oldPos);
            }
#endif
            }
        }


#endif

    }



//legacy code
/*
               //pressure & viscosity
               Vector3 accumulatorP = Vector3.zero;
               Vector3 accumulatorV = Vector3.zero;

               foreach (MyParticle n in neighbours)
               {
                   float pipjmj = (p.fPressure + n.fPressure) * n.fMass;
                   float dj2 = 2 * n.fMassDensity;

                   accumulatorP = accumulatorP + ((pipjmj/dj2) * GradientPressure(p.v3Position - n.v3Position));

                   //-------------------------------------------------------------------------------------------------------------

                   Vector3 diffVel = n.v3Velocity - p.v3Velocity;
                   float md = n.fMass / n.fMassDensity;

                   accumulatorV = accumulatorV + (diffVel * md * LaplacianViscosity(p.v3Position - n.v3Position));
               }

               Vector3 pressureForce = (-1 * accumulatorP);
              // Debug.Log("AP: " + accumulatorP);

               //Debug.Log("pressure: " + pressureForce + "  Acc: " + (pressureForce / p.fMassDensity));

               Vector3 viscosityForce = (VISCOSITY_COEFFICIENT * accumulatorV);
               //Debug.Log("viscosity: " + viscosityForce + "  Acc: " + (viscosityForce / p.fMassDensity));

               p.InternalForce = pressureForce + viscosityForce;
               */
//(pressureForce + viscosityForce);