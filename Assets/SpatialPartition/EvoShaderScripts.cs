using System.Collections;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using UnityEngine;

/* particle data structure layout
struct Particle
{
	float4 position;
	float4 velocity;
	float4 predictedPos;
	uint id;
	uint cellId;
	float pad;
	float pad2;
};
*/
//Structured Resource View -> StructuredBuffer
//Unordered Access View -> RWStructuredBuffer

public struct CSParticle
{
    public Vector4 position;
    public Vector4 velocity;
    public Vector4 predictedPos;
    public uint id;
    public uint cellId;
    public float pad;
    public float pad2;
};

//used for testing, the same as the uint2 in hlsl
public struct uint2
{
    public uint x;
    public uint y;
}

public struct uint3
{
    public uint x;
    public uint y;
    public uint z;
}

//param set 1
//KRAD = 0.1;
//solver iter = 4;
//density 6578;
//viscosity 0.16;
//vort 0.00012;
//s_corr 0.00012;

public class EvoShaderScripts : MonoBehaviour {

    bool executionLock = true;

    //for debugging
  
    //-------------------animation parameters----------------
    private const int SOLVER_ITERATION = 2;
    private Vector3 G = new Vector3(0, -9.8f, 0);

    private float REST_DENSITY = 1000;//998.29F;
    private static int particleAmount = 65536;//79872;//65536;
    private static int sideLength = 48;
    private int TOTAL_PARTICLES = particleAmount;
    private float TIME_STEP = 0.016f;
    private float KRAD;
    private float RELXATION_PARAM = 500;

    private float Poly6Constant;
    private float GradSpikyConstant;

    //------------------------------------------------------
    //fixed to 128 for simplicity and align with shared memory size
    private const int THREADS_GROUP_SIZE = 128;
    private int iWarpCount;
    private int iCellWarpCount;

    //structured buffer for compute shader-----------------------------
    //ping-pong particle buffer----------------
    ComputeBuffer sb_particles;
    ComputeBuffer sb_sortedParticles;
    //-----------------NNS---------------------
    ComputeBuffer sb_particleCellId;
    ComputeBuffer sb_particleInsertIdx;
    ComputeBuffer sb_gridCounter;
    ComputeBuffer sb_gridPrefixSum;
    //-------------NNS(Prefix Sum)-------------
    ComputeBuffer sb_scanaux0;
    ComputeBuffer sb_scanaux1;
    //-------------Simulation------------------
    ComputeBuffer sb_poscorr; //position correction
    ComputeBuffer sb_omega;
    ComputeBuffer sb_lambda;
    ComputeBuffer sb_density;
    //-----------------------------------------
    //structured buffer end--------------------------------------------

    //Id of the kernel used--------------------------------------------
    //Initialization-------------------------
    private int i_cs_InitializerID;
    //Position&Velocity update---------------
    private int i_cs_UpdateExternelForceID;
    private int i_cs_OmegaID;
    private int i_cs_UpdateVelocityID;
    private int i_cs_UpdatePositionID;
    //Animation------------------------------
    private int i_cs_LambdaID;
    private int i_cs_PosCorrectionID;
    private int i_cs_UpdatePredPosID;
    //NNS------------------------------------
    private int i_cs_UpdateGridID;
    private int i_cs_SortPtlID;
    private int i_cs_ResetNNSBufferID;
    //NNS(Prefix sum)------------------------
    private int i_cs_ScanInBucketID;
    private int i_cs_ScanBucketResultID;
    private int i_cs_ScanAddBueckResultID;
    //Kernel ID End----------------------------------------------------


    // Compute shader used to initialize the Particles.
    public ComputeShader cs_Initializer;
    public ComputeShader cs_TestGrav;
    public ComputeShader cs_Lambda;
    public ComputeShader cs_NNS;

    Material newMat;
    Material newMat2;
    //-----------------------------------------------------
    //grid parameters - should at least 128 cell large
    private int cellNumX = 32;
    private int cellNumY = 16;
    private int cellNumZ = 16;

    [SerializeField]
    private Vector3 gridOrigin;
    [SerializeField]
    private Vector3 gridCenter;
    [SerializeField]
    private Vector3 gridSize;        //should be int3
    [SerializeField]
    private float cellSize;
    [SerializeField]
    private float particleSize;
    [SerializeField]
    private float particleMass;

    //collision box
    public Vector3 boundBoxCenter;
    public Vector3 boundBoxExtent;
    Vector3 tempBBCenter;
    Vector3 tempBBExtent;


    //rendering---------------------------------------
    private Material m_depthMaterial;
    private Material m_particleMaterial;
    private Material m_blurredDepthMaterial;
    private Material m_normalMaterial;
    private Material m_surfaceShadingMaterial;
    private Material m_thicknessMaterial;

    public RenderTexture tex_depthTexture;
    public RenderTexture tex_blurredDepthTexture;
    public RenderTexture tex_tempBlurredDepthTexture;
    public RenderTexture tex_normalTexture;
    public RenderTexture tex_thicknessTexture;

    public Cubemap m_cubemap;

    private string shader_points = "PBF/PointSprite";
    private string shader_depth = "PBF/SphereDepth";
    private string shader_blur = "PBF/SphereBlurredDepth";
    private string shader_norm = "PBF/SphereNormal";
    private string shader_surf = "PBF/SurfaceShading";
    private string shader_thic = "PBF/Thickness";
    //------------------------------------------------

   
    void OnDestroy()
    {
        if (sb_particles != null)
            sb_particles.Release();
        if (sb_sortedParticles != null)
            sb_sortedParticles.Release();

        if (sb_particleCellId != null)
            sb_particleCellId.Release();
        if (sb_particleInsertIdx != null)
            sb_particleInsertIdx.Release();
        if (sb_gridCounter != null)
            sb_gridCounter.Release();
        if (sb_gridPrefixSum != null)
            sb_gridPrefixSum.Release();

        if (sb_scanaux0 != null)
            sb_scanaux0.Release();
        if (sb_scanaux1 != null)
            sb_scanaux1.Release();

        if (sb_poscorr != null)
            sb_poscorr.Release();
        if (sb_omega != null)
            sb_omega.Release();
        if (sb_lambda != null)
            sb_lambda.Release();
        if (sb_density != null)
            sb_density.Release();
    }

    // Use this for initialization
    void Start () {
        KRAD = 0.2f;//Mathf.Pow(((3 * VOLUME * KERNEL_PARTICLE) / (4 * Mathf.PI * PARTICLES_PER_VOLUME)), 1.0f / 3.0f);
        particleMass = 1f;//REST_DENSITY * VOLUME / PARTICLES_PER_VOLUME;
        particleSize = (0.2f / 2.2f);//Mathf.Pow((particleMass * 3) / (4 * Mathf.PI * REST_DENSITY), (1.0f / 3.0f));

        Debug.Log("kernel radius" + KRAD);
        Debug.Log("diameter" + particleSize);

        //kenrel constant terms
        Poly6Constant = 315.0f / (64.0f * Mathf.PI * Mathf.Pow(KRAD, 9.0f));
        GradSpikyConstant = -45.0f / (Mathf.PI * Mathf.Pow(KRAD, 6.0f));

        // Calculate the number of warp needed to handle all the particles
        if (TOTAL_PARTICLES <= 0)
            TOTAL_PARTICLES = 1;
        iWarpCount = Mathf.CeilToInt((float)TOTAL_PARTICLES / THREADS_GROUP_SIZE);

        GridInit();
        ComputeBufferInit();
        KernelIDInit();
        KernelConstantInit();

        // Bind the ComputeBuffer to the shader and the compute shader
        cs_Initializer.SetBuffer(i_cs_InitializerID, "particleBuffer", sb_particles);

        //set float variable in compute shader
        cs_Initializer.Dispatch(i_cs_InitializerID, iWarpCount, 1, 1);
       

        //create material for rendering
        m_particleMaterial = new Material(Shader.Find(shader_points));
        m_particleMaterial.hideFlags = HideFlags.HideAndDontSave;

        m_depthMaterial = new Material(Shader.Find(shader_depth));
        m_depthMaterial.hideFlags = HideFlags.HideAndDontSave;

        m_blurredDepthMaterial = new Material(Shader.Find(shader_blur));
        m_blurredDepthMaterial.hideFlags = HideFlags.HideAndDontSave;

        m_normalMaterial = new Material(Shader.Find(shader_norm));
        m_normalMaterial.hideFlags = HideFlags.HideAndDontSave;

        m_surfaceShadingMaterial = new Material(Shader.Find(shader_surf));
        m_surfaceShadingMaterial.hideFlags = HideFlags.HideAndDontSave;

        m_thicknessMaterial = new Material(Shader.Find(shader_thic));
        m_thicknessMaterial.hideFlags = HideFlags.HideAndDontSave;

#if false
        newMat = Resources.Load("GPUInstancing", typeof(Material)) as Material;
        foreach (CSParticle p in particles)
        {
            //Debug.Log("Mass: " + p.fMass);
            //Debug.Log("Pos: " + p.v3Position);
            //Debug.Log("Mass: " + p.fMass);

            goLists[p.id] = GameObject.CreatePrimitive(PrimitiveType.Cube);

            goLists[p.id].transform.position = p.v3Position;
            goLists[p.id].hideFlags = HideFlags.HideAndDontSave;
            //GPU Instancing
            goLists[p.id].GetComponent<Renderer>().material = newMat;
            goLists[p.id].transform.localScale = new Vector3(particleSize, particleSize, particleSize);
        }
#endif
        //nns test
#if false

        //testing only
        uint[] counter = new uint[64];
        uint[] prefixSum = new uint[64];

        uint2[] gridIdx = new uint2[256];

        uint[] sortedparticle = new uint[TOTAL_PARTICLES];
        //-----------------

        
        cs_NNS.SetBuffer(i_cs_UpdateGridID, "_ParticleBuffer_R", sb_particles);
        cs_NNS.SetBuffer(i_cs_UpdateGridID, "_ParticleGridIdx_W", sb_particleGridIdx);
        cs_NNS.SetBuffer(i_cs_UpdateGridID, "_GridCounter_W", sb_gridCounter);

        cs_NNS.SetBuffer(i_cs_PrefixSumID, "_GridCounter_W", sb_gridCounter);
        cs_NNS.SetBuffer(i_cs_PrefixSumID, "_GridPrefixSum_W", sb_gridPrefixSum);

        cs_NNS.SetBuffer(i_cs_SortPtlID, "_ParticleGridIdx_W", sb_particleGridIdx);
        cs_NNS.SetBuffer(i_cs_SortPtlID, "_GridPrefixSum_W", sb_gridPrefixSum);
        cs_NNS.SetBuffer(i_cs_SortPtlID, "_ParticleBuffer_R", sb_particles);
        cs_NNS.SetBuffer(i_cs_SortPtlID, "_SortedParticleBuffer_W", sb_sortedParticles);

        cs_NNS.Dispatch(i_cs_UpdateGridID, iWarpCount, 1, 1);     //dispatched per particle
        cs_NNS.Dispatch(i_cs_PrefixSumID, iCellWarpCount, 1, 1);  //dispatched per cell
        cs_NNS.Dispatch(i_cs_SortPtlID, iWarpCount, 1, 1);        //dispatched per particle

        sb_gridCounter.GetData(counter);
        sb_gridPrefixSum.GetData(prefixSum);
        sb_sortedParticles.GetData(sortedparticle);
        sb_particleGridIdx.GetData(gridIdx);
          
        //uint[] ids = new uint[125];

        //InitUintArray(ids);

       //Debug.Log("output length: " + sortedparticle.Length);
       //foreach (CSParticle p in sortedparticle)
       //for (int i = 0; i < sortedparticle.Length; i++)
       //{
           //Debug.Log(p.id);
           //Debug.Log(sortedparticle[i].id);
           //ids[sortedparticle[i].id]++;
       //}

       /*
       for(int i = 0; i < ids.Length; i++)
       {
           if (ids[i] != 1)
           {
               Debug.Log("count: " + ids[i] + "  " + "id: " + i);
           }
       }
       */

       Debug.Log("aaaaaaaaaaaaaaaaaaaaaaaa");
       
       for (int i = 0; i < gridIdx.Length; i++)
       {
           Debug.Log("id: " + i + "  " + gridIdx[i].x + "  " + gridIdx[i].y);
       }

       Debug.Log("aaaaaaaaaaaaaaaaaaaaaaaa");

       foreach (uint i in counter)
       {
           Debug.Log(i);
       }
       
       Debug.Log("aaaaaaaaaaaaaaaaaaaaaaaa");
       for (int i = 0; i < prefixSum.Length; i++)
       {
           Debug.Log("cell Id: " + i + "  " + "sum: " + prefixSum[i]);
       }
       Debug.Log("aaaaaaaaaaaaaaaaaaaaaaaa");


        //NNS testing code
        /*
        foreach (CSParticle p in sortedparticle)
        {
            Vector3 pos = p.v3Position;
            goLists[p.id].transform.position = new Vector3(pos.x, pos.y, pos.z - 1 * gridIdx[p.id].x);
        }
        */
        
       uint testCell = gridIdx[255].x;
        Vector3 testPos = particles[255].v3Position;
        Vector3 h = new Vector3(KRAD, KRAD, KRAD);

        Vector3 BBmin = testPos - h;
        Vector3 BBmax = testPos + h;

        Debug.Log("query cell index: " + testCell);

        uint3 BBCellMin = testCellIdGen(BBmin);
        uint3 BBCellMax = testCellIdGen(BBmax);

        List<int> n = new List<int>();

        for (int i = (int)BBCellMin.x; i <= (int)BBCellMax.x; i++)
        {
            for (int j = (int)BBCellMin.y; j <= (int)BBCellMax.y; j++)
            {
                for (int k = (int)BBCellMin.z; k <= (int)BBCellMax.z; k++)
                {
                    if (i >= 0 && i < 4 && j >= 0 && j < 4 && k >= 0 && k < 4)
                    {
                       Debug.Log("neighbour cell index: " +
                       (new Vector3(i, j, k)) + "   " + 
                       (i + j * 4 + k * 4 * 4)
                       );

                        n.Add((i + j * 4 + k * 4 * 4));
                    }
                    else
                    {
                        Debug.Log("neighbour cell index: " +
                        (new Vector3(i, j, k)) + "   " +
                        (i + j * 4 + k * 4 * 4)
                        );

                    }
                }
            }
        }

        /*
        for (int i = 0; i < gridIdx.Length; i++)
        {
            int id = (int)gridIdx[i].x;
            if (id == 1 || id == 5 || id == 4 || id == 16 || id == 0 || id == 20 || id == 17 || id == 21)
            //if(id == 9)
            {
                Vector3 pos = particles[i].v3Position;
                goLists[i].transform.position = new Vector3(pos.x, pos.y, pos.z - 1);
            }
        }
        */

        /*
        foreach (int i in n)
        {
            Debug.Log("neighbour cell: " + i);

            foreach (CSParticle p in sortedparticle)
            {
                if (gridIdx[p.id].x == i)
                {
                    Vector3 pos = p.v3Position;
                    goLists[p.id].transform.position = new Vector3(pos.x, pos.y, pos.z - 1);
                }
            }

        }
        */

        
        newMat2 = Resources.Load("MyMat2", typeof(Material)) as Material;
        goLists[255].GetComponent<Renderer>().material = newMat2;
        foreach (int cellId in n)
        {
            Debug.Log("cellId: " + cellId);
            for (int i = (int)prefixSum[cellId]; i < (int)(prefixSum[cellId] + counter[cellId]); i++)
            {
                uint pid = sortedparticle[i];

                Vector3 pos = goLists[pid].transform.position;

                if ((pos - particles[255].v3Position).magnitude <= KRAD)
                {
                    goLists[pid].transform.position = new Vector3(pos.x - 2f, pos.y, pos.z);
                }
            }
        }
        
#endif

        //lambda test
#if false
        /*
        cs_TestGrav.SetBuffer(i_cs_TestGravID, "_ParticleBuffer_W", sb_particles);
        cs_TestGrav.SetBuffer(i_cs_UpdatePredPosID, "_ParticleBuffer_W", sb_particles);
        cs_TestGrav.SetBuffer(i_cs_UpdatePV, "_ParticleBuffer_W", sb_particles);
        */

        cs_NNS.SetBuffer(i_cs_UpdateGridID, "_ParticleBuffer_R", sb_particles);
        cs_NNS.SetBuffer(i_cs_UpdateGridID, "_ParticleGridIdx_W", sb_particleGridIdx);
        cs_NNS.SetBuffer(i_cs_UpdateGridID, "_GridCounter_W", sb_gridCounter);

        cs_NNS.SetBuffer(i_cs_PrefixSumID, "_GridCounter_W", sb_gridCounter);
        cs_NNS.SetBuffer(i_cs_PrefixSumID, "_GridPrefixSum_W", sb_gridPrefixSum);

        cs_NNS.SetBuffer(i_cs_SortPtlID, "_ParticleGridIdx_W", sb_particleGridIdx);
        cs_NNS.SetBuffer(i_cs_SortPtlID, "_GridPrefixSum_W", sb_gridPrefixSum);
        cs_NNS.SetBuffer(i_cs_SortPtlID, "_ParticleBuffer_R", sb_particles);
        cs_NNS.SetBuffer(i_cs_SortPtlID, "_SortedParticleBuffer_W", sb_sortedParticles);

        cs_NNS.Dispatch(i_cs_UpdateGridID, iWarpCount, 1, 1);     //dispatched per particle
        cs_NNS.Dispatch(i_cs_PrefixSumID, iCellWarpCount, 1, 1);  //dispatched per cell
        cs_NNS.Dispatch(i_cs_SortPtlID, iWarpCount, 1, 1);        //dispatched per particle

        cs_Lambda.SetBuffer(i_cs_LambdaID, "_ParticleBuffer_RW", sb_particles);
        cs_Lambda.SetBuffer(i_cs_LambdaID, "_ParticleGridIdx_R", sb_particleGridIdx);
        cs_Lambda.SetBuffer(i_cs_LambdaID, "_GridCounter_R", sb_gridCounter);
        cs_Lambda.SetBuffer(i_cs_LambdaID, "_GridPrefixSum_R", sb_gridPrefixSum);
        cs_Lambda.SetBuffer(i_cs_LambdaID, "_SortedParticleBuffer_R", sb_sortedParticles);

        cs_Lambda.SetBuffer(i_cs_PosCorrectionID, "_ParticleBuffer_RW", sb_particles);
        cs_Lambda.SetBuffer(i_cs_PosCorrectionID, "_ParticleGridIdx_R", sb_particleGridIdx);
        cs_Lambda.SetBuffer(i_cs_PosCorrectionID, "_GridCounter_R", sb_gridCounter);
        cs_Lambda.SetBuffer(i_cs_PosCorrectionID, "_GridPrefixSum_R", sb_gridPrefixSum);
        cs_Lambda.SetBuffer(i_cs_PosCorrectionID, "_SortedParticleBuffer_R", sb_sortedParticles);

        cs_Lambda.Dispatch(i_cs_LambdaID, iWarpCount, 1, 1);
        cs_Lambda.Dispatch(i_cs_PosCorrectionID, iWarpCount, 1, 1);



        sb_particles.GetData(particles);

        foreach (CSParticle p in particles)
        {
            Debug.Log("lambda: " + p.lambda);
            Debug.Log("position correction: " + p.v3PtlPosCorrection);
        }
#endif
    }

    //nns testing use only
    uint3 testCellIdGen(Vector3 pos)
    {
        uint3 result;
        result.x = (uint)Mathf.Floor((pos.x - gridOrigin.x) / cellSize);
        result.y = (uint)Mathf.Floor((pos.y - gridOrigin.y) / cellSize);
        result.z = (uint)Mathf.Floor((pos.z - gridOrigin.z) / cellSize);

        return result;
    }

    void GridInit()
    {
        cellSize = 2 * KRAD;
        gridSize = new Vector3(cellNumX, cellNumY, cellNumZ);        //should be int3, num be cells
        gridOrigin = this.transform.position -  (new Vector3(KRAD, KRAD, KRAD));
        gridCenter = gridOrigin + (gridSize / 2) * cellSize;

        float size = gridSize.x * gridSize.y * gridSize.z;
        iCellWarpCount = Mathf.CeilToInt(size / THREADS_GROUP_SIZE);

        tempBBCenter = boundBoxCenter = gridCenter;

        Vector3 offset = new Vector3(particleSize, particleSize, particleSize);
        tempBBExtent = boundBoxExtent = (gridSize / 2) * cellSize - offset;

        Debug.Log("boundBoxCenter: " + boundBoxCenter);
        Debug.Log("boundBoxExtent: " + boundBoxExtent);

    }

    void KernelConstantInit()
    {
        cs_Initializer.SetFloat("mass", particleMass);
        cs_Initializer.SetVector("coordCenter", this.transform.position);
        cs_Initializer.SetInt("sideLength", sideLength);
        cs_Initializer.SetFloat("step", particleSize);

        cs_NNS.SetInt("_numGridX", cellNumX);
        cs_NNS.SetInt("_numGridY", cellNumY);
        cs_NNS.SetInt("_numGridZ", cellNumZ);
        cs_NNS.SetVector("_gridOrigin", gridOrigin);
        cs_NNS.SetFloat("_cellSize", cellSize);
        cs_NNS.SetFloat("_cellCount", cellNumX * cellNumY * cellNumZ);

        cs_Lambda.SetFloat("_KRAD", KRAD);
        cs_Lambda.SetFloat("_REST_DENSITY", REST_DENSITY);
        cs_Lambda.SetFloat("_RELXATION_PARAM", RELXATION_PARAM);
        cs_Lambda.SetInt("_numGridX", cellNumX);
        cs_Lambda.SetInt("_numGridY", cellNumY);
        cs_Lambda.SetInt("_numGridZ", cellNumZ);
        cs_Lambda.SetVector("_gridOrigin", gridOrigin);
        cs_Lambda.SetFloat("_cellSize", cellSize);
        cs_Lambda.SetVector("_boundBoxCenter", boundBoxCenter);
        cs_Lambda.SetVector("_boundBoxExtent", boundBoxExtent);
        cs_Lambda.SetFloat("_Poly6Const", Poly6Constant);
        cs_Lambda.SetFloat("_GradSpikyConst", GradSpikyConstant);
        cs_Lambda.SetFloat("_ptlMass", particleMass);

        cs_TestGrav.SetFloat("_timestep", TIME_STEP);
        cs_TestGrav.SetVector("_gravity", G);
        cs_TestGrav.SetFloat("_KRAD", KRAD);
        cs_TestGrav.SetInt("_numGridX", cellNumX);
        cs_TestGrav.SetInt("_numGridY", cellNumY);
        cs_TestGrav.SetInt("_numGridZ", cellNumZ);
        cs_TestGrav.SetFloat("_Poly6Const", Poly6Constant);
        cs_TestGrav.SetFloat("_GradSpikyConst", GradSpikyConstant);
    }

    void ComputeBufferInit()
    {
        //stride is hardcoded to prevent any inconsistency of data size between c# and hlsl
        //----------------particles-------------------------------------------
        sb_particles = new ComputeBuffer(TOTAL_PARTICLES, 64); //64 = 3xfloat4 + 2xuint + 2xfloat
        sb_sortedParticles = new ComputeBuffer(TOTAL_PARTICLES, 64);
        //----------------NNS-------------------------------------------------
        int cellNum = (int)(gridSize.x * gridSize.y * gridSize.z);
        sb_particleCellId = new ComputeBuffer(TOTAL_PARTICLES, 4); //4bytes uint
        sb_particleInsertIdx = new ComputeBuffer(TOTAL_PARTICLES, 4); //4bytes uint
        sb_gridCounter = new ComputeBuffer(cellNum, 4); //4bytes uint
        sb_gridPrefixSum = new ComputeBuffer(cellNum, 4); //4bytes uint
        sb_scanaux0 = new ComputeBuffer(cellNum, 4); //4bytes uint
        sb_scanaux1 = new ComputeBuffer(cellNum, 4); //4bytes uint
        //----------------Animation-------------------------------------------
        sb_poscorr = new ComputeBuffer(TOTAL_PARTICLES, 16); //position correction float4 - 4 bytes float
        sb_omega = new ComputeBuffer(TOTAL_PARTICLES, 16); //float4 - 4 bytes float
        sb_lambda = new ComputeBuffer(TOTAL_PARTICLES, 4);
        sb_density = new ComputeBuffer(TOTAL_PARTICLES, 4);
    }

    void KernelIDInit()
    {
        i_cs_InitializerID = cs_Initializer.FindKernel("Init");
        //NNS--------------------------------------------------------
        i_cs_UpdateGridID = cs_NNS.FindKernel("UpdateGrid");
        i_cs_SortPtlID = cs_NNS.FindKernel("SortPtl");
        i_cs_ResetNNSBufferID = cs_NNS.FindKernel("ResetNNSBuffer");
        //NNS(Prefix Sum)--------------------------------------------
        i_cs_ScanInBucketID = cs_NNS.FindKernel("CSScanInBucket");
        i_cs_ScanBucketResultID = cs_NNS.FindKernel("CSScanBucketResult");
        i_cs_ScanAddBueckResultID = cs_NNS.FindKernel("CSScanAddBucketResult");
        //Position Update--------------------------------------------
        i_cs_UpdateExternelForceID = cs_TestGrav.FindKernel("UpdateExternelForce");
        i_cs_UpdateVelocityID = cs_TestGrav.FindKernel("UpdateVelocity");
        i_cs_UpdatePositionID = cs_TestGrav.FindKernel("UpdatePosition");
        i_cs_OmegaID = cs_TestGrav.FindKernel("CalOmega");
        //Animation--------------------------------------------------
        i_cs_LambdaID = cs_Lambda.FindKernel("CalLambda");
        i_cs_PosCorrectionID = cs_Lambda.FindKernel("CalPositionCorrection");
        i_cs_UpdatePredPosID = cs_Lambda.FindKernel("UpdatePredictedPos");
    }

    void Update()
    {
        if (Input.GetKeyUp(KeyCode.Space))
        {
            executionLock = !executionLock;
            Debug.Log("Locker: " + executionLock);
        }

        if (Input.GetKey(KeyCode.A))
        {
            if (tempBBExtent.x > 0)
            {
                tempBBCenter.x -= 0.02f;
                tempBBExtent.x -= 0.02f;
                cs_Lambda.SetVector("_boundBoxCenter", tempBBCenter);
                cs_Lambda.SetVector("_boundBoxExtent", tempBBExtent);
            }
        }

        if (Input.GetKey(KeyCode.D))
        {
            if (tempBBExtent.x < boundBoxExtent.x)
            {
                tempBBCenter.x += 0.02f;
                tempBBExtent.x += 0.02f;
                cs_Lambda.SetVector("_boundBoxCenter", tempBBCenter);
                cs_Lambda.SetVector("_boundBoxExtent", tempBBExtent);
            }
        }
    }

    void FixedUpdate()
    {
        if (executionLock == false)
        {
            //Update gravity force---------------------------------------------
            cs_TestGrav.SetBuffer(i_cs_UpdateExternelForceID, "_ParticleBuffer_RW", sb_particles);
            cs_TestGrav.Dispatch(i_cs_UpdateExternelForceID, iWarpCount, 1, 1);
            //-----------------------------------------------------------------



            //Nearest Neightbour Search--------------------------------------------------------------------
            //reset work space-------------------------------------------------
            cs_NNS.SetBuffer(i_cs_ResetNNSBufferID, "_GridPrefixSum_W", sb_gridPrefixSum);
            cs_NNS.SetBuffer(i_cs_ResetNNSBufferID, "_GridCounter_W", sb_gridCounter);
            //clean workspace(also initialize workspace in first loop)
            cs_NNS.Dispatch(i_cs_ResetNNSBufferID, iCellWarpCount, 1, 1);
            //-----------------------------------------------------------------

            //grid allocation--------------------------------------------------
            cs_NNS.SetBuffer(i_cs_UpdateGridID, "_ParticleBuffer_RW", sb_particles);
            cs_NNS.SetBuffer(i_cs_UpdateGridID, "_ParticleCellId_W", sb_particleCellId);
            cs_NNS.SetBuffer(i_cs_UpdateGridID, "_ParticleInsertIdx_W", sb_particleInsertIdx);
            cs_NNS.SetBuffer(i_cs_UpdateGridID, "_GridCounter_W", sb_gridCounter);
            
            cs_NNS.Dispatch(i_cs_UpdateGridID, iWarpCount, 1, 1);     //dispatched per particle
            //-----------------------------------------------------------------

            //parallel prefix sum----------------------------------------------
            cs_NNS.SetBuffer(i_cs_ScanInBucketID, "Input", sb_gridCounter);
            cs_NNS.SetBuffer(i_cs_ScanInBucketID, "Result", sb_scanaux0);

            cs_NNS.Dispatch(i_cs_ScanInBucketID, iCellWarpCount, 1, 1); //dispatched per grid cell

            cs_NNS.SetBuffer(i_cs_ScanBucketResultID, "Input", sb_scanaux0);
            cs_NNS.SetBuffer(i_cs_ScanBucketResultID, "Result", sb_scanaux1);

            cs_NNS.Dispatch(i_cs_ScanBucketResultID, 1, 1, 1);

            cs_NNS.SetBuffer(i_cs_ScanAddBueckResultID, "Input", sb_scanaux0);
            cs_NNS.SetBuffer(i_cs_ScanAddBueckResultID, "Input1", sb_scanaux1);
            cs_NNS.SetBuffer(i_cs_ScanAddBueckResultID, "Result", sb_gridPrefixSum);

            cs_NNS.Dispatch(i_cs_ScanAddBueckResultID, iCellWarpCount, 1, 1);
            //-----------------------------------------------------------------

            //sort
            cs_NNS.SetBuffer(i_cs_SortPtlID, "_ParticleCellId_W", sb_particleCellId);
            cs_NNS.SetBuffer(i_cs_SortPtlID, "_ParticleInsertIdx_W", sb_particleInsertIdx);
            cs_NNS.SetBuffer(i_cs_SortPtlID, "_GridPrefixSum_W", sb_gridPrefixSum);
            cs_NNS.SetBuffer(i_cs_SortPtlID, "_ParticleBuffer_RW", sb_particles);
            cs_NNS.SetBuffer(i_cs_SortPtlID, "_SortedParticleBuffer_W", sb_sortedParticles);

            cs_NNS.Dispatch(i_cs_SortPtlID, iWarpCount, 1, 1);        //dispatched per particle
            //-----------------------------------------------------------------

            //swap ping-pong buffer
            ComputeBuffer t;
            t = sb_particles;
            sb_particles = sb_sortedParticles;
            sb_sortedParticles = t;
            //Nearest Neightbour Search End----------------------------------------------------------------



            //Animation Start------------------------------------------------------------------------------
            int iter = 0;
            while (iter < SOLVER_ITERATION)
            {
                //calculate lambda -- position correction force--------------
                cs_Lambda.SetBuffer(i_cs_LambdaID, "_ParticleBuffer_RW", sb_particles);
                cs_Lambda.SetBuffer(i_cs_LambdaID, "_GridCounter_R", sb_gridCounter);
                cs_Lambda.SetBuffer(i_cs_LambdaID, "_GridPrefixSum_R", sb_gridPrefixSum);
                cs_Lambda.SetBuffer(i_cs_LambdaID, "_Lambda_RW", sb_lambda);
                cs_Lambda.SetBuffer(i_cs_LambdaID, "_Density_RW", sb_density);

                cs_Lambda.Dispatch(i_cs_LambdaID, iWarpCount, 1, 1);
                //-----------------------------------------------------------

                //calulate position correction-------------------------------
                cs_Lambda.SetBuffer(i_cs_PosCorrectionID, "_ParticleBuffer_RW", sb_particles);
                cs_Lambda.SetBuffer(i_cs_PosCorrectionID, "_GridCounter_R", sb_gridCounter);
                cs_Lambda.SetBuffer(i_cs_PosCorrectionID, "_GridPrefixSum_R", sb_gridPrefixSum);
                cs_Lambda.SetBuffer(i_cs_PosCorrectionID, "_Lambda_RW", sb_lambda);
                cs_Lambda.SetBuffer(i_cs_PosCorrectionID, "_PosCorr_RW", sb_poscorr);

                cs_Lambda.Dispatch(i_cs_PosCorrectionID, iWarpCount, 1, 1);
                //-----------------------------------------------------------

                //apply correction and collision-----------------------------
                cs_Lambda.SetBuffer(i_cs_UpdatePredPosID, "_ParticleBuffer_RW", sb_particles);
                cs_Lambda.SetBuffer(i_cs_UpdatePredPosID, "_PosCorr_RW", sb_poscorr);

                cs_Lambda.Dispatch(i_cs_UpdatePredPosID, iWarpCount, 1, 1);
                //-----------------------------------------------------------

                iter++;
            }

            //update velocity---------------------------------------------
            cs_TestGrav.SetBuffer(i_cs_UpdateVelocityID, "_ParticleBuffer_RW", sb_particles);

            cs_TestGrav.Dispatch(i_cs_UpdateVelocityID, iWarpCount, 1, 1);
            //------------------------------------------------------------

            //calculate omega---------------------------------------------
            cs_TestGrav.SetBuffer(i_cs_OmegaID, "_ParticleBuffer_RW", sb_particles);
            cs_TestGrav.SetBuffer(i_cs_OmegaID, "_GridCounter_R", sb_gridCounter);
            cs_TestGrav.SetBuffer(i_cs_OmegaID, "_GridPrefixSum_R", sb_gridPrefixSum);
            cs_TestGrav.SetBuffer(i_cs_OmegaID, "_Omega_RW", sb_omega);

            cs_TestGrav.Dispatch(i_cs_OmegaID, iWarpCount, 1, 1);
            //------------------------------------------------------------

            //update position---------------------------------------------
            cs_TestGrav.SetBuffer(i_cs_UpdatePositionID, "_ParticleBuffer_RW", sb_particles);
            cs_TestGrav.SetBuffer(i_cs_UpdatePositionID, "_GridCounter_R", sb_gridCounter);
            cs_TestGrav.SetBuffer(i_cs_UpdatePositionID, "_GridPrefixSum_R", sb_gridPrefixSum);
            cs_TestGrav.SetBuffer(i_cs_UpdatePositionID, "_Omega_RW", sb_omega);
            cs_TestGrav.SetBuffer(i_cs_UpdatePositionID, "_Density_R", sb_density);

            cs_TestGrav.Dispatch(i_cs_UpdatePositionID, iWarpCount, 1, 1);
            //------------------------------------------------------------






            //error probe
            //sb_particles.GetData(particles);
            //foreach (CSParticle p in particles) { Debug.Log("id: " + p.id + "  velocity: " + p.velocity); }
            //sb_gridPrefixSum.GetData(test);
            //for (uint i = 0; i < test.Length; i++) { Debug.Log("id: " + i + " sum: " + test[i]); }
            //sb_poscorr.GetData(testf4);
            //for (uint i = 0; i < testf4.Length; i++) { Debug.Log("poscorr " + testf4[i]); }
            //sb_density.GetData(testf);
            //for (uint i = 0; i < testf.Length; i++) { Debug.Log("density " + testf[i]); }
            //sb_lambda.GetData(testf);
            //for (uint i = 0; i < testf.Length; i++) { Debug.Log("lambda " + testf[i]); }



#if false
            sb_particles.GetData(particles);
            foreach (CSParticle p in particles)
            {
                //Debug.Log("lambda: " + p.lambda);
                //Debug.Log("position correction: " + p.v3PtlPosCorrection);
                Debug.Log("omega: " + p.omega);
                //goLists[p.id].transform.position = p.v3Position;
            }
#endif



#if false
            sb_particles.GetData(particles);
            foreach (CSParticle p in particles)
            {
                //Debug.Log("lambda: " + p.lambda);
                //Debug.Log("position correction: " + p.v3PtlPosCorrection);
                Debug.Log("omega: " + p.omega);
                //goLists[p.id].transform.position = p.v3Position;
            }
#endif
#if false
            sb_particles.GetData(particles);
            foreach (CSParticle p in particles)
            {
                Debug.Log("-----------");
                Debug.Log("ID: " + p.id);
                Debug.Log("velocity: " + p.v3Velocity);
                Debug.Log("Density: " + p.fdensity);
                Debug.Log("Corr: " + p.v3PtlPosCorrection);
                Debug.Log("Position: " + p.v3Position);
                //goLists[p.id].transform.position = p.v3Position;
            }
#endif
#if false
            uint[] counter = new uint[64];
            uint[] prefixSum = new uint[64];

            uint2[] gridIdx = new uint2[256];

            uint[] sortedparticle = new uint[TOTAL_PARTICLES];

            sb_gridCounter.GetData(counter);
            sb_gridPrefixSum.GetData(prefixSum);
            sb_sortedParticles.GetData(sortedparticle);
            sb_particleGridIdx.GetData(gridIdx);

            Debug.Log("aaaaaaaaaaaaaaaaaaaaaaaa");

            for (int i = 0; i < gridIdx.Length; i++)
            {
                Debug.Log("id: " + i + "  " + gridIdx[i].x + "  " + gridIdx[i].y);
            }

            Debug.Log("aaaaaaaaaaaaaaaaaaaaaaaa");

            for (int i = 0; i < counter.Length; i++)
            {
                Debug.Log("id: " + i + "  " + "counter: " + counter[i]);
            }

            Debug.Log("aaaaaaaaaaaaaaaaaaaaaaaa");
            for (int i = 0; i < prefixSum.Length; i++)
            {
                Debug.Log("cell Id: " + i + "  " + "sum: " + prefixSum[i]);
            }
#endif

            //ScreenPositionUpdate();

            //step lock
            //executionLock = true;
        }
    }

    //float pointScale = 800.0f / Mathf.Tan(45.0f * (0.5f * Mathf.PI / 180.0f));
    //rendering
    //point sprite
    public float m_particleRadius = 0.1f;
    //smoothing parameters
    public float m_blurScale = 3f;//0.001f;
    public int m_blurRadius = 11;
    public float m_blurDepthFalloff = 15;
    public int m_filterRadius = 7;
    //surface shading parameters
    public Vector4 m_surfaceColor = new Vector4(0.5f,0.5f,0.5f,1);
	public float m_refractionIndex = 1.5f;
    public int m_shiness = 1;  //specular reflection

    //rendering step control
    public bool r_isBlur = true;
    public bool r_isNormal = true;
    public bool r_isSurfaceShading = true;
    public bool r_isTransparent = true;
    public bool r_isRaw = false;
    void OnRenderObject()
    {
        if (r_isRaw)
        {
            PointSprite();
        }
        else
        {
            Depth();
            BlurredDepth(r_isBlur);
            Normal(r_isNormal);
            Thickness();
            SurfaceShading(r_isSurfaceShading);
        }
 
        Graphics.SetRenderTarget(Camera.main.targetTexture);
    }

    void PointSprite()
    {
        // point sprite
        //particleMaterial.SetFloat("_pointScale", pointScale);
        m_particleMaterial.SetBuffer("_ParticlesBuffer", sb_particles);
        m_particleMaterial.SetFloat("_ptlRadius", m_particleRadius);
        m_particleMaterial.SetPass(0);
        //Debug.Log("running");
        Graphics.DrawProcedural(MeshTopology.Points, particleAmount, 1);
    }

    void Depth()
    {
        //save depth as texture------------------------------------------
        m_depthMaterial.SetBuffer("_ParticlesBuffer", sb_particles);
        m_depthMaterial.SetFloat("_ptlRadius", m_particleRadius);
        m_depthMaterial.SetPass(0);

        Graphics.SetRenderTarget(tex_depthTexture);
        GL.Clear(true, true, Color.white);

        Graphics.DrawProcedural(MeshTopology.Points, particleAmount, 1);
        //---------------------------------------------------------------
    }

    void BlurredDepth(bool sSwitch)
    {
        //bilateral filtering (separated bilateral filter, not strictly correct but much faster)
        if (sSwitch)
        {
            m_blurredDepthMaterial.SetInt("radius", m_blurRadius);
            m_blurredDepthMaterial.SetFloat("blurDepthFalloff", m_blurDepthFalloff);
            m_blurredDepthMaterial.SetInt("filterRadius", m_filterRadius);

            //smoothing on x axis
            m_blurredDepthMaterial.SetTexture("_DepthTex", tex_depthTexture);
            m_blurredDepthMaterial.SetFloat("scaleX", 1.0f / 512 * m_blurScale);
            m_blurredDepthMaterial.SetFloat("scaleY", 0.0f);
            Graphics.Blit(tex_depthTexture, tex_tempBlurredDepthTexture, m_blurredDepthMaterial, -1);

            //smoothing on y axis
            m_blurredDepthMaterial.SetTexture("_DepthTex", tex_tempBlurredDepthTexture);
            m_blurredDepthMaterial.SetFloat("scaleX", 0.0f);
            m_blurredDepthMaterial.SetFloat("scaleY", 1.0f / 512 * m_blurScale);
            Graphics.Blit(tex_tempBlurredDepthTexture, tex_blurredDepthTexture, m_blurredDepthMaterial, -1);
            //debugging only
            //Graphics.Blit(tex_tempBlurredDepthTexture, Camera.main.targetTexture, m_blurredDepthMaterial, -1);
        }
        else
        {
            Graphics.Blit(tex_depthTexture, tex_blurredDepthTexture);
            //debugging only
            //Graphics.Blit(tex_depthTexture, Camera.main.targetTexture);
        }
        //---------------------------------------------------------------------------------------
    }

    void Normal(bool sSwitch)
    {
        //calculate normal from depth------------------------------------------------------------
        if (sSwitch)
        {
            //all texture used in here is 512x512
            m_normalMaterial.SetInt("textureWidth", tex_normalTexture.width);
            m_normalMaterial.SetInt("textureHeight", tex_normalTexture.height);

            //smoothed depth texture
            m_normalMaterial.SetTexture("_DepthTex", tex_blurredDepthTexture);

            Graphics.Blit(tex_blurredDepthTexture, tex_normalTexture, m_normalMaterial, -1);
            //debugging only 
            //Graphics.Blit(tex_blurredDepthTexture, Camera.main.targetTexture, m_normalMaterial, -1);
        }
        else
        {
            Graphics.Blit(tex_blurredDepthTexture, tex_normalTexture);
            //debugging only
            //Graphics.Blit(tex_blurredDepthTexture, Camera.main.targetTexture);
        }
        //---------------------------------------------------------------------------------------
    }

    void SurfaceShading(bool sSwitch)
    {
        //opaque surface shading-----------------------------------------------------------------
        if (sSwitch)
        {
            //m_surfaceShadingMaterial.SetTexture("_BGTex", Camera.main.targetTexture);
            m_surfaceShadingMaterial.SetTexture("_DepthTex", tex_depthTexture);
            m_surfaceShadingMaterial.SetTexture("_NormalTex", tex_normalTexture);
            m_surfaceShadingMaterial.SetTexture("_ThicknessTex", tex_thicknessTexture);
            m_surfaceShadingMaterial.SetTexture("_Cube", m_cubemap);

            m_surfaceShadingMaterial.SetVector("color", m_surfaceColor);
            m_surfaceShadingMaterial.SetFloat("refractionIndex", m_refractionIndex);
            m_surfaceShadingMaterial.SetInt("shininess", m_shiness);

            m_surfaceShadingMaterial.SetMatrix("viewMatrix", Camera.main.worldToCameraMatrix);

            if(r_isTransparent)
                m_surfaceShadingMaterial.SetInt("isTransparent", 1);
            else
                m_surfaceShadingMaterial.SetInt("isTransparent", 0);

            Graphics.Blit(tex_normalTexture, Camera.main.targetTexture, m_surfaceShadingMaterial, -1);
        }
        else
        {
            Graphics.Blit(tex_normalTexture, Camera.main.targetTexture);
        }
        //---------------------------------------------------------------------------------------
    }

    void Thickness()
    {
        //Thickness------------------------------------------------------------------------------
        m_thicknessMaterial.SetBuffer("_ParticlesBuffer", sb_particles);
        m_thicknessMaterial.SetFloat("_ptlRadius", 0.12f);
        m_thicknessMaterial.SetPass(0);

        Graphics.SetRenderTarget(tex_thicknessTexture);
        GL.Clear(true, true, Color.black);

        Graphics.DrawProcedural(MeshTopology.Points, particleAmount, 1);

        //---------------------------------------------------------------------------------------
        m_blurredDepthMaterial.SetInt("radius", 11);
        m_blurredDepthMaterial.SetFloat("blurDepthFalloff", 2);
        m_blurredDepthMaterial.SetInt("filterRadius", 5);

        //smoothing on x axis
        m_blurredDepthMaterial.SetTexture("_DepthTex", tex_thicknessTexture);
        m_blurredDepthMaterial.SetFloat("scaleX", 1.0f / tex_depthTexture.width * 2);
        m_blurredDepthMaterial.SetFloat("scaleY", 0.0f);
        Graphics.Blit(tex_thicknessTexture, tex_tempBlurredDepthTexture, m_blurredDepthMaterial, -1);

        //smoothing on y axis
        m_blurredDepthMaterial.SetTexture("_DepthTex", tex_tempBlurredDepthTexture);
        m_blurredDepthMaterial.SetFloat("scaleX", 0.0f);
        m_blurredDepthMaterial.SetFloat("scaleY", 1.0f / tex_depthTexture.height * 2);
        Graphics.Blit(tex_tempBlurredDepthTexture, tex_thicknessTexture, m_blurredDepthMaterial, -1);
        //debugging only
        //Graphics.Blit(tex_tempBlurredDepthTexture, Camera.main.targetTexture, m_blurredDepthMaterial, -1);
    }

    Vector3 CubeSize = new Vector3(0.05f, 0.05f, 0.05f);
    void OnDrawGizmos()
    {
       Gizmos.DrawWireCube(gridCenter, gridSize * cellSize);

        /*
        if (executionLock == false)
        {
            sb_particles.GetData(particles);

            Gizmos.color = Color.cyan;
            foreach (CSParticle p in particles)
            {
                Gizmos.DrawWireCube(p.v3Position, CubeSize);
            }
        }
        */
        
    }

    /*
    void ScreenPositionUpdate()
    {
        //send data back to main memory
        sb_particles.GetData(particles);

        foreach (CSParticle p in particles)
        {
            //Debug.Log("lambda: " + p.lambda);
            //Debug.Log("position correction: " + p.v3PtlPosCorrection);

            goLists[p.id].transform.position = p.v3Position;
        }
    }
    */

}
