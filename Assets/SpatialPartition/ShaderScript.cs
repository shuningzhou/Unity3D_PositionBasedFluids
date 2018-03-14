using System.Collections;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using UnityEngine;

using Assets.SpatialPartition;


#if false
public struct CSParticle
{
    public uint id;
    public float fMass;
    public Vector3 v3Position;
    public Vector3 v3Velocity;
    public Vector3 v3PredictedPos;
    public Vector3 v3PtlPosCorrection;
    public Vector3 omega; //vorticity;
    public float lambda;  //constraint parameters
}



public class ShaderScript : MonoBehaviour
{
    //particle structure

    private int iPtlMemSize;

    CSSpatialHashTable hashTable;

    uint[] neighboursList;       //2D array represented in 1D to send to shader
    Vector2[] neighboutsListIdx; //(start index, the number of neighbours);

    CSParticle[] particles;

    GameObject[] goLists;

    //particle container - structured buffer
    ComputeBuffer sb_particles;

    ComputeBuffer sb_neighboursList;
    ComputeBuffer sb_neighboursListIndex;


    //update lock controlled by the user
    private bool bActionLock = true;

    //-------------------constant---------------------------
    private const int SOLVER_ITERATION = 1;

    private const float GRAVITY = -0.098f;

    private const int WARP_SIZE = 256;

    private float REST_DENSITY = 998.29f;//998.29F;
    private float MATERIAL_DENSITY = 998.29f;

    private float GAS_STIFFNESS = 25f;//15f;//18.0f;
    private float VISCOSITY_COEFFICIENT = 15f;//15f; //3.5f;

    private float SURFACE_TENSION_COEFFICIENT = 0.1456f;//0.0728f;

    //----------------------------------------------------------------------
    private static int sideLength = 5;  //side length of generation cube
    //----------------------------------------------------------------------

    private int TOTAL_PARTICLES = sideLength * sideLength * sideLength;
    private int PARTICLES_PER_VOLUME = sideLength * sideLength * sideLength;   //how many particles are there in a unit of volume of water e.g. 100 particles to form 1m3 water
    private float VOLUME = 1f;
    private int KERNEL_PARTICLE = 20;

    private float TIME_STEP = 0.02f;

    private float KRAD;

    private float RELXATION_PARAM = 0;
    
    private Vector3 G = new Vector3(0, -9.8f, 0);
    //------------------------------------------------------
    private int iWarpCount;

    //Id of the kernel used
    private int i_cs_InitializerID;

    private int i_cs_TestGravID;

    private int i_cs_LambdaID;

    private int i_cs_PosCorrectionID;

    private int i_cs_UpdatePredPosID;

    // Compute shader used to initialize the Particles.
    public ComputeShader cs_Initializer;

    public ComputeShader cs_TestGrav;

    public ComputeShader cs_Lambda;

    Material newMat;
    Material newMat2;

    void OnDestroy()
    {
        if (sb_particles != null)
            sb_particles.Release();
        if (sb_neighboursList != null)
            sb_neighboursList.Release();
        if (sb_neighboursListIndex != null)
            sb_neighboursListIndex.Release();
    }

    private void Start()
    {
        KRAD = Mathf.Pow(((3 * VOLUME * KERNEL_PARTICLE) / (4 * Mathf.PI * PARTICLES_PER_VOLUME)), 1.0f / 3.0f);
        float mass = REST_DENSITY * VOLUME / PARTICLES_PER_VOLUME;
        float step = Mathf.Pow((mass * 3) / (4 * Mathf.PI * REST_DENSITY), (1.0f / 3.0f));

        //initialize-hash-table-------------------------------------------
        hashTable = new CSSpatialHashTable(TOTAL_PARTICLES, KRAD); //radius is not correct
        //----------------------------------------------------------------

        // Calculate the number of warp needed to handle all the particles
        if (TOTAL_PARTICLES <= 0)
            TOTAL_PARTICLES = 1;
        iWarpCount = Mathf.CeilToInt((float)TOTAL_PARTICLES / WARP_SIZE);

        particles = new CSParticle[TOTAL_PARTICLES];
        goLists = new GameObject[TOTAL_PARTICLES];

        iPtlMemSize = Marshal.SizeOf(new CSParticle());

        // Create the ComputeBuffer holding the Particles
        sb_particles = new ComputeBuffer(TOTAL_PARTICLES, iPtlMemSize);
        sb_particles.SetData(particles);

        // Find the id of the kernel
        i_cs_InitializerID = cs_Initializer.FindKernel("Init");

        // Bind the ComputeBuffer to the shader and the compute shader
        cs_Initializer.SetBuffer(i_cs_InitializerID, "particleBuffer", sb_particles);
        //material.SetBuffer("particleBuffer", particleBuffer);

        //set float variable in compute shader

        cs_Initializer.SetFloat("mass", mass);
        cs_Initializer.SetVector("coordCenter", this.transform.position);
        cs_Initializer.SetInt("sideLength", sideLength);
        cs_Initializer.SetFloat("step", step);

        cs_Initializer.Dispatch(i_cs_InitializerID, iWarpCount, 1, 1);

        sb_particles.GetData(particles);

        newMat = Resources.Load("GPUInstancing", typeof(Material)) as Material;
        foreach (CSParticle p in particles)
        {
            //Debug.Log("Mass: " + p.fMass);
            Debug.Log("Pos: " + p.v3Position);
            Debug.Log("Mass: " + p.fMass);

            goLists[p.id] = GameObject.CreatePrimitive(PrimitiveType.Sphere);
            goLists[p.id].transform.position = p.v3Position;

            //GPU Instancing
            goLists[p.id].GetComponent<Renderer>().material = newMat;
            goLists[p.id].transform.localScale = new Vector3(step, step, step);
        }

        /*
        //--------------------------------------------------------------------------------
        hashTable.NeighboursQuery(particles, ref neighboursList, ref neighboutsListIdx, false);

        int z = 30;
        Vector3 pos;

        for (int i = (int)neighboutsListIdx[z].x; i < (int)neighboutsListIdx[z + 1].x; i++)
        {
            int idx = (int)neighboursList[i];
            pos = goLists[idx].transform.position;
            goLists[idx].transform.position = new Vector3(pos.x - 5f, pos.y, pos.z);
        }
        goLists[z].GetComponent<Renderer>().material = newMat2;

        sb_neighboursList = new ComputeBuffer(neighboursList.Length, sizeof(uint));
        sb_neighboursList.SetData(neighboursList);

        int v2Size = Marshal.SizeOf(new Vector2());
        sb_neighboursListIndex = new ComputeBuffer(neighboutsListIdx.Length, v2Size);
        sb_neighboursListIndex.SetData(neighboutsListIdx);

        //sb_particles.SetData(particles);

        i_cs_LambdaID = cs_Lambda.FindKernel("CalLambda");

        cs_Lambda.SetBuffer(i_cs_LambdaID, "particleBuffer", sb_particles);
        cs_Lambda.SetBuffer(i_cs_LambdaID, "neighboursList", sb_neighboursList);
        cs_Lambda.SetBuffer(i_cs_LambdaID, "neighboursListIndex", sb_neighboursListIndex);

        cs_Lambda.SetFloat("KRAD", KRAD);
        cs_Lambda.SetFloat("REST_DENSITY", REST_DENSITY);
        cs_Lambda.SetFloat("RELXATION_PARAM", RELXATION_PARAM);

        cs_Lambda.Dispatch(i_cs_LambdaID, iWarpCount, 1, 1);

        sb_particles.GetData(particles);

        foreach (CSParticle p in particles)
        {
            Debug.Log("My lambda: " + p.lambda);
        }

        i_cs_PosCorrectionID = cs_Lambda.FindKernel("CalPositionCorrection");

        cs_Lambda.SetBuffer(i_cs_PosCorrectionID, "particleBuffer", sb_particles);
        cs_Lambda.SetBuffer(i_cs_PosCorrectionID, "neighboursList", sb_neighboursList);
        cs_Lambda.SetBuffer(i_cs_PosCorrectionID, "neighboursListIndex", sb_neighboursListIndex);

        cs_Lambda.Dispatch(i_cs_PosCorrectionID, iWarpCount, 1, 1);

        sb_particles.GetData(particles);

        foreach (CSParticle p in particles)
        {
            Debug.Log("My PosCorr: " + p.v3PtlPosCorrection);
        }
        */
        //--------------------------------------------------------------------------------

        /*legacy hash table
        //hashing------------------------------------------------------------------
        hashTable.AddParticles(particles);

        //-------------------------------------------------------------------------
        
        neighbourIdx = hashTable.NeighboursQuery(particles, false);
        */

        /*
         * for neighbour index test only
        foreach (uint[] idxs in neighbourIdx)
        {
            Debug.Log("--------------------");
            foreach (uint idx in idxs)
            {
                Debug.Log("Neighbours id: " + idx);
            }
            Debug.Log("--------------------");
        }
        */
        //newMat2 = Resources.Load("MyMat2", typeof(Material)) as Material;

        i_cs_LambdaID = cs_Lambda.FindKernel("CalLambda");
        i_cs_PosCorrectionID = cs_Lambda.FindKernel("CalPositionCorrection");

        cs_Lambda.SetFloat("KRAD", KRAD);
        cs_Lambda.SetFloat("REST_DENSITY", REST_DENSITY);
        cs_Lambda.SetFloat("RELXATION_PARAM", RELXATION_PARAM);

        i_cs_TestGravID = cs_TestGrav.FindKernel("TestGravity");
        cs_TestGrav.SetFloat("timestep", TIME_STEP);
        cs_TestGrav.SetFloat("gravityConst", GRAVITY);

        i_cs_UpdatePredPosID = cs_TestGrav.FindKernel("UpdatePredictedPos");

		//---------------------------------------------------------------------------------
		hashTable.NeighboursQuery(particles, ref neighboursList, ref neighboutsListIdx, false);

		sb_neighboursList = new ComputeBuffer(neighboursList.Length, sizeof(uint));
		sb_neighboursList.SetData(neighboursList);

		sb_neighboursListIndex = new ComputeBuffer(neighboutsListIdx.Length, v2Size);
		sb_neighboursListIndex.SetData(neighboutsListIdx);

		cs_TestGrav.SetBuffer(i_cs_TestGravID, "particleBuffer", sb_particles);

		cs_Lambda.SetBuffer(i_cs_LambdaID, "particleBuffer", sb_particles);
		cs_Lambda.SetBuffer(i_cs_LambdaID, "neighboursList", sb_neighboursList);
		cs_Lambda.SetBuffer(i_cs_LambdaID, "neighboursListIndex", sb_neighboursListIndex);

		cs_Lambda.SetBuffer(i_cs_PosCorrectionID, "particleBuffer", sb_particles);
		cs_Lambda.SetBuffer(i_cs_PosCorrectionID, "neighboursList", sb_neighboursList);
		cs_Lambda.SetBuffer(i_cs_PosCorrectionID, "neighboursListIndex", sb_neighboursListIndex);

		cs_TestGrav.SetBuffer(i_cs_UpdatePredPosID, "particleBuffer", sb_particles);

		sb_particles.SetData(particles);
    }

    private void Update()
    {
        //neighbours search test
        
        /*
        if (Input.GetKey(KeyCode.W))
        {
            Vector3 pos;
           
        }
        */

    }

    int v2Size = Marshal.SizeOf(new Vector2());

    private void FixedUpdate()
    {
        if (Input.GetKeyUp(KeyCode.Space))
        {
            bActionLock = !bActionLock;
            Debug.Log("Locker: " + bActionLock);


        }



        /*simple test
        //animation procedure
        if (bActionLock == false)
        {
            //putting set buffer call to here would make the data in gpu adapte to any external changes on particles but slower which is terrible.
            //cs_TestGrav.SetBuffer(i_cs_TestGravID, "particleBuffer", sb_particles);
            cs_TestGrav.Dispatch(i_cs_TestGravID, iWarpCount, 1, 1);

            sb_particles.GetData(particles);

            foreach(CSParticle p in particles)
            {
                goLists[p.id].transform.position = p.v3Position;
            }
        }
        */

        if (bActionLock == false)
        {
            //putting set buffer call to here would make the data in gpu accomodate to any external changes on particles but slower which is terrible.
            //cs_TestGrav.SetBuffer(i_cs_TestGravID, "particleBuffer", sb_particles);
            //1-4
            cs_TestGrav.Dispatch(i_cs_TestGravID, iWarpCount, 1, 1);
            sb_particles.GetData(particles);

			CSParticle p;
			for(int i = 0; i < particles.Length; i++)
			{
				p = particles[i];
				goLists[p.id].transform.position = p.v3Position;

				if (p.id == 0) 
				{
					Debug.Log ("lambda " + p.lambda);
					//Debug.Log ("posCorr" + p.v3PtlPosCorrection);
					//Debug.Log ("neighbour count: " + neighboutsListIdx [p.id].y);
				}
			}

            //5-7
			neighboursList = null;
			neighboutsListIdx = null;
            hashTable.NeighboursQuery(particles, ref neighboursList, ref neighboutsListIdx, false);

			sb_neighboursList.Release();
			sb_neighboursListIndex.Release();

            sb_neighboursList = new ComputeBuffer(neighboursList.Length, sizeof(uint));
            sb_neighboursList.SetData(neighboursList);

			sb_neighboursListIndex = new ComputeBuffer(neighboutsListIdx.Length, v2Size);
            sb_neighboursListIndex.SetData(neighboutsListIdx);

			cs_Lambda.SetBuffer(i_cs_LambdaID, "neighboursList", sb_neighboursList);
			cs_Lambda.SetBuffer(i_cs_LambdaID, "neighboursListIndex", sb_neighboursListIndex);
            //8-19

            int iter = 0;
            while (iter < SOLVER_ITERATION)
            {
                //9 - 11
                
                cs_Lambda.Dispatch(i_cs_LambdaID, iWarpCount, 1, 1);

                //sb_particles.GetData(particles);

                //12 - 15  //lack of collision
                //sb_particles.SetData(particles);

                
                //cs_Lambda.Dispatch(i_cs_PosCorrectionID, iWarpCount, 1, 1);

                //sb_particles.GetData(particles);

                //16 - 18
                //sb_particles.SetData(particles);
                
                //cs_TestGrav.Dispatch(i_cs_UpdatePredPosID, iWarpCount, 1, 1);
				//sb_particles.GetData(particles);

                iter++;
            }
        }
    }
}
#endif