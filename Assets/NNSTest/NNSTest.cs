using System.Collections;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using UnityEngine;

public struct NNSParticle
{
    public Vector4 v3Position;
    public Vector4 v3Velocity;
    public Vector4 v3PredictedPos;
    public uint id;
    public uint cellId;
    public float pad;
    public float pad2;
}

public class NNSTest : MonoBehaviour {

    NNSParticle[] particles;

    GameObject[] goLists;
    
    //-------------------animation parameters----------------

    private static int particleAmount = 256;
    private static int sideLength = 8;
    private int TOTAL_PARTICLES = particleAmount;
    private float KRAD;

    //------------------------------------------------------
    //fixed to 128 for simplicity and align with shared memory size
    private const int WARP_SIZE = 128;
    private int iWarpCount;
    private int iCellWarpCount;

    //structured buffer
    ComputeBuffer sb_particles;

    //NNS working space
    ComputeBuffer sb_particleCellId;
    ComputeBuffer sb_particleInsertIdx;
    ComputeBuffer sb_gridCounter;
    ComputeBuffer sb_gridPrefixSum;
    ComputeBuffer sb_sortedParticles;

    //ComputeBuffer sb_scanbuf0;
    //ComputeBuffer sb_scanbuf1;
    //ComputeBuffer sb_scanaux;

    ComputeBuffer sb_scanaux0;
    ComputeBuffer sb_scanaux1;

    //ComputeBuffer sb_gridInfo;
    //ComputeBuffer sb_sortedParticles;

    //intermediate variable dont need initialization - no get/set
    ComputeBuffer sb_positionCorrection;
    ComputeBuffer sb_omega;
    ComputeBuffer sb_lambda;
    ComputeBuffer sb_density;

    //Id of the kernel used
    private int i_cs_InitializerID;
    private int i_cs_UpdateGridID;
    private int i_cs_SortPtlID;
    private int i_cs_ResetNNSBufferID;

    private int i_cs_ScanInBucketID;
    private int i_cs_ScanBucketResultID;
    private int i_cs_ScanAddBueckResultID;

    // Compute shader used to initialize the Particles.
    public ComputeShader cs_NNS;
    public ComputeShader cs_Initializer;

    public Material newMat;
    public Material newMat2;
    //-----------------------------------------------------
    //grid parameters - must be pow of 2 for simplicity
    private int cellNumX = 8;
    private int cellNumY = 8;
    private int cellNumZ = 8;

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
   
    void OnDestroy()
    {
        if (sb_particles != null)
            sb_particles.Release();
        if (sb_particleCellId != null)
            sb_particleCellId.Release();
        if (sb_particleInsertIdx != null)
            sb_particleInsertIdx.Release();
        if (sb_gridCounter != null)
          sb_gridCounter.Release();
        if (sb_gridPrefixSum != null)
          sb_gridPrefixSum.Release();
        if (sb_sortedParticles != null)
          sb_sortedParticles.Release();

        /*
        if (sb_scanbuf0 != null)
            sb_scanbuf0.Release();
        if (sb_scanbuf1 != null)
            sb_scanbuf1.Release();
        if (sb_scanaux != null)
            sb_scanaux.Release();
            */
        if (sb_scanaux0 != null)
            sb_scanaux0.Release();
        if (sb_scanaux1 != null)
            sb_scanaux1.Release();
    }

    // Use this for initialization
    void Start () {
        GL.Clear(true, true, new Color(1f,1f,1f,1f));

        //GameObject aa = GameObject.CreatePrimitive(PrimitiveType.Cube);


#if true
        KRAD = 0.24f;//Mathf.Pow(((3 * VOLUME * KERNEL_PARTICLE) / (4 * Mathf.PI * PARTICLES_PER_VOLUME)), 1.0f / 3.0f);
        particleMass = 1f;//REST_DENSITY * VOLUME / PARTICLES_PER_VOLUME;
        particleSize = (0.24f / 2.2f);//Mathf.Pow((particleMass * 3) / (4 * Mathf.PI * REST_DENSITY), (1.0f / 3.0f));
        
        Debug.Log("kernel radius"  + KRAD);
        Debug.Log("diameter" + particleSize);

        // Calculate the number of warp needed to handle all the particles
        if (TOTAL_PARTICLES <= 0)
            TOTAL_PARTICLES = 1;
        iWarpCount = Mathf.CeilToInt((float)TOTAL_PARTICLES / WARP_SIZE);

        particles = new NNSParticle[TOTAL_PARTICLES];
        goLists = new GameObject[TOTAL_PARTICLES];

        GridInit();

        ComputeBufferInit();
        KernelIDInit();
        KernelConstantInit();

        // Bind the ComputeBuffer to the shader and the compute shader
        cs_Initializer.SetBuffer(i_cs_InitializerID, "particleBuffer", sb_particles);

        //set float variable in compute shader
        cs_Initializer.Dispatch(i_cs_InitializerID, iWarpCount, 1, 1);
        sb_particles.GetData(particles);

        //newMat = Resources.Load("/Materials/MyMat", typeof(Material)) as Material;
        foreach (NNSParticle p in particles)
        {
            uint id = p.id;
            //Debug.Log("Mass: " + p.fMass);
            //Debug.Log("Pos: " + p.v3Position);
            //Debug.Log("Mass: " + p.fMass);
            //Debug.Log("id: " + id);
            goLists[id] = GameObject.CreatePrimitive(PrimitiveType.Cube);

            goLists[id].transform.position = p.v3Position;
            goLists[id].hideFlags = HideFlags.HideInHierarchy;

            goLists[id].GetComponent<Renderer>().material = newMat;
            goLists[id].transform.localScale = new Vector3(particleSize, particleSize, particleSize);
        }

        //testing only
        uint[] counter = new uint[512];
        uint[] prefixSum = new uint[512];

        uint[] cellIds = new uint[256];
        uint[] insertIdx = new uint[256];

        NNSParticle[] sortedparticle = new NNSParticle[TOTAL_PARTICLES];
        //-----------------
        //only one go don't need to reset
        cs_NNS.SetBuffer(i_cs_ResetNNSBufferID, "_GridCounter_W", sb_gridCounter);
        cs_NNS.SetBuffer(i_cs_ResetNNSBufferID, "_GridPrefixSum_W", sb_gridPrefixSum);

        cs_NNS.Dispatch(i_cs_ResetNNSBufferID, iCellWarpCount, 1, 1);

        //float updateGridStart = Time.realtimeSinceStartup * 1000;

        cs_NNS.SetBuffer(i_cs_UpdateGridID, "_ParticleBuffer_RW", sb_particles);
        cs_NNS.SetBuffer(i_cs_UpdateGridID, "_ParticleCellId_W", sb_particleCellId);
        cs_NNS.SetBuffer(i_cs_UpdateGridID, "_ParticleInsertIdx_W", sb_particleInsertIdx);
        cs_NNS.SetBuffer(i_cs_UpdateGridID, "_GridCounter_W", sb_gridCounter);

        cs_NNS.Dispatch(i_cs_UpdateGridID, iWarpCount, 1, 1);     //dispatched per particle

        //sb_particleCellId.GetData(cellIds);
        //float updateGridEnd = Time.realtimeSinceStartup * 1000;
        //Debug.Log("UpdateGrid Time: " + (updateGridEnd - updateGridStart));

        //float myPrefixSumStart = Time.realtimeSinceStartup;

        //legacy prefix sum - extremely slow
        //cs_NNS.SetBuffer(i_cs_PrefixSumID, "_GridCounter_W", sb_gridCounter);
        //cs_NNS.SetBuffer(i_cs_PrefixSumID, "_GridPrefixSum_W", sb_gridPrefixSum);

        //cs_NNS.Dispatch(i_cs_PrefixSumID, iCellWarpCount, 1, 1);  //dispatched per cell

        //new prefix sum---------------------------------------------------
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

        //sb_gridPrefixSum.GetData(prefixSum);
        //float myPrefixSumEnd = Time.realtimeSinceStartup;

        cs_NNS.SetBuffer(i_cs_SortPtlID, "_ParticleCellId_W", sb_particleCellId);
        cs_NNS.SetBuffer(i_cs_SortPtlID, "_ParticleInsertIdx_W", sb_particleInsertIdx);
        cs_NNS.SetBuffer(i_cs_SortPtlID, "_GridPrefixSum_W", sb_gridPrefixSum);
        cs_NNS.SetBuffer(i_cs_SortPtlID, "_ParticleBuffer_RW", sb_particles);
        cs_NNS.SetBuffer(i_cs_SortPtlID, "_SortedParticleBuffer_W", sb_sortedParticles);

        cs_NNS.Dispatch(i_cs_SortPtlID, iWarpCount, 1, 1);        //dispatched per particle

        //swap ping-pong buffer
        SwapBuffer(sb_particles, sb_sortedParticles);

        sb_particles.GetData(particles);
        sb_particleCellId.GetData(cellIds);
        sb_particleInsertIdx.GetData(insertIdx);
        sb_gridCounter.GetData(counter);
        sb_gridPrefixSum.GetData(prefixSum);
        sb_sortedParticles.GetData(sortedparticle);

        /*
        //directx sdk prefixsum test--------------------------------------------
        int testSize = 512;
        sb_scanbuf0 = new ComputeBuffer(testSize, sizeof(uint));
        sb_scanbuf1 = new ComputeBuffer(testSize, sizeof(uint));
        sb_scanaux = new ComputeBuffer(testSize, sizeof(uint));
        uint[] prefixSumTest = new uint[testSize];
        uint[] psResult = new uint[testSize];
        for (int i = 0; i < testSize; i++)
        {
            prefixSumTest[i] = 1;
        }

        //directx sdk prefixsum test
        sb_scanbuf0.SetData(counter);

        int scanInBucket = cs_NNS.FindKernel("CSScanInBucket");
        int scanBucketResult = cs_NNS.FindKernel("CSScanBucketResult");
        int scanAddBucketResult = cs_NNS.FindKernel("CSScanAddBucketResult");

        //float dxPrefixSumStart = Time.realtimeSinceStartup; 

        cs_NNS.SetBuffer(scanInBucket, "Input", sb_scanbuf0);
        cs_NNS.SetBuffer(scanInBucket, "Result", sb_scanbuf1);

        cs_NNS.Dispatch(scanInBucket, testSize/128, 1, 1);

        cs_NNS.SetBuffer(scanBucketResult, "Input", sb_scanbuf1);
        cs_NNS.SetBuffer(scanBucketResult, "Result", sb_scanaux);

        cs_NNS.Dispatch(scanBucketResult, 1, 1, 1);

        cs_NNS.SetBuffer(scanAddBucketResult, "Input", sb_scanbuf1);
        cs_NNS.SetBuffer(scanAddBucketResult, "Input1", sb_scanaux);
        cs_NNS.SetBuffer(scanAddBucketResult, "Result", sb_scanbuf0);

        cs_NNS.Dispatch(scanAddBucketResult, testSize/128, 1, 1);

        //sb_scanbuf0.GetData(psResult);
        //float dxPrefixSumEnd = Time.realtimeSinceStartup;

        //Debug.Log("My Prefix Sum: " + (myPrefixSumEnd * 1000 - myPrefixSumStart*1000) + "DX Prefix Sum: " + (dxPrefixSumEnd*1000 - dxPrefixSumStart *1000));

        sb_scanbuf0.GetData(psResult);
        //------------------------------------------------------------------
        */


        //uint[] ids = new uint[125];

        //InitUintArray(ids);

        //Debug.Log("output length: " + sortedparticle.Length);
        //foreach (NNSParticle p in sortedparticle)
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

        for (int i = 0; i < sortedparticle.Length; i++)
        {
            //Debug.Log("id: " + i + "  particleId: " + sortedparticle[i].id + " " + "cellId: " + sortedparticle[i].cellId);
        }

        //Debug.Log("aaaaaaaaaaaaaaaaaaaaaaaa");

        for (int i = 0; i < cellIds.Length; i++)
        {
            //Debug.Log("id: " + i + "  " + cellIds[i] + "  " + insertIdx[i]);
        }

        //Debug.Log("aaaaaaaaaaaaaaaaaaaaaaaa");

        for (int i = 0; i < counter.Length; i++)
        {
            //Debug.Log("cellid: " + i + "  " + "elements count: " + counter[i]);
        }
        
        //Debug.Log("aaaaaaaaaaaaaaaaaaaaaaaa");
        for (int i = 0; i < prefixSum.Length; i++)
        {
            Debug.Log("cell Id: " + i + "  " + "sum: " + prefixSum[i]); // + " sum2: " + psResult[i]);
        }
        //Debug.Log("aaaaaaaaaaaaaaaaaaaaaaaa");

        uint testId = 255;
       uint testCell = cellIds[testId];
        Vector3 testPos = particles[testId].v3Position;
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
                    if (i >= 0 && i < 8 && j >= 0 && j < 8 && k >= 0 && k < 8)
                    {
                       Debug.Log("neighbour cell index: " +
                       (new Vector3(i, j, k)) + "   " + 
                       (i + j * 8 + k * 8 * 8)
                       );

                        n.Add((i + j * 8 + k * 8 * 8));
                    }
                    else
                    {
                        /*
                        Debug.Log("neighbour cell index: " +
                        (new Vector3(i, j, k)) + "   " +
                        (i + j * 8 + k * 8 * 8)
                        );
                        */
                    }
                }
            }
        }
        
        //newMat2 = Resources.Load("/Materials/MyMat2", typeof(Material)) as Material;
        goLists[testId].GetComponent<Renderer>().material = newMat2;
        foreach (int cellId in n)
        {
            Debug.Log("cellId: " + cellId);
            for (int i = (int)prefixSum[cellId]; i < (int)(prefixSum[cellId] + counter[cellId]); i++)
            {
                uint pid = sortedparticle[i].id;

                Vector3 pos = goLists[pid].transform.position;

                Vector3 PtlPos = new Vector3(particles[testId].v3Position.x, particles[testId].v3Position.y, particles[testId].v3Position.z);

                if ((pos - PtlPos).magnitude <= KRAD)
                {
                    goLists[pid].transform.position = new Vector3(pos.x - 5 * KRAD, pos.y, pos.z);
                }
            }
        }
#endif
    }

    //since compute buffer are just pointers, swapping overhead is minimal
    void SwapBuffer(ComputeBuffer a, ComputeBuffer b)
    {
        ComputeBuffer t;
        t = a;
        a = b;
        b = t;
        
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
        cellSize = KRAD;
        gridSize = new Vector3(cellNumX, cellNumY, cellNumZ);        //should be int3, num be cells
        gridOrigin = this.transform.position -  (new Vector3(KRAD, KRAD, KRAD));
        gridCenter = gridOrigin + (gridSize / 2) * cellSize;

        float size = gridSize.x * gridSize.y * gridSize.z;
        iCellWarpCount = Mathf.CeilToInt(size / WARP_SIZE);

        
        Vector3 offset = new Vector3(particleSize, particleSize, particleSize);
       
        Debug.Log("boundBoxCenter: " + boundBoxCenter);
        Debug.Log("boundBoxExtent: " + boundBoxExtent);

    }

    void KernelConstantInit()
    {
        cs_NNS.SetInt("_numGridX", cellNumX);
        cs_NNS.SetInt("_numGridY", cellNumY);
        cs_NNS.SetInt("_numGridZ", cellNumZ);
        cs_NNS.SetVector("_gridOrigin", gridOrigin);
        cs_NNS.SetFloat("_cellSize", cellSize);
        cs_NNS.SetFloat("_cellCount", cellNumX * cellNumY * cellNumZ);

        cs_Initializer.SetFloat("mass", particleMass);
        cs_Initializer.SetVector("coordCenter", this.transform.position);
        cs_Initializer.SetInt("sideLength", sideLength);
        cs_Initializer.SetFloat("step", particleSize);

    }

    void ComputeBufferInit()
    {
        //----------------particles-------------------------------------------
        //64 = 4 * float4
        sb_particles = new ComputeBuffer(TOTAL_PARTICLES, 64);
        sb_sortedParticles = new ComputeBuffer(TOTAL_PARTICLES, 64);
        //----------------NNS-------------------------------------------------
        int cellNum = (int)(gridSize.x * gridSize.y * gridSize.z);
        sb_particleCellId  = new ComputeBuffer(TOTAL_PARTICLES, sizeof(uint));
        sb_particleInsertIdx = new ComputeBuffer(TOTAL_PARTICLES, sizeof(uint));
        sb_gridCounter = new ComputeBuffer(cellNum, sizeof(uint));
        sb_gridPrefixSum = new ComputeBuffer(cellNum, sizeof(uint));
        sb_scanaux0 = new ComputeBuffer(cellNum, sizeof(uint));
        sb_scanaux1 = new ComputeBuffer(cellNum, sizeof(uint));
    }

    void KernelIDInit()
    {
        i_cs_UpdateGridID = cs_NNS.FindKernel("UpdateGrid");
        //i_cs_PrefixSumID = cs_NNS.FindKernel("PrefixSum");
        i_cs_SortPtlID = cs_NNS.FindKernel("SortPtl");
        i_cs_ResetNNSBufferID = cs_NNS.FindKernel("ResetNNSBuffer");

        i_cs_ScanInBucketID = cs_NNS.FindKernel("CSScanInBucket");
        i_cs_ScanBucketResultID = cs_NNS.FindKernel("CSScanBucketResult");
        i_cs_ScanAddBueckResultID = cs_NNS.FindKernel("CSScanAddBucketResult");

        i_cs_InitializerID = cs_Initializer.FindKernel("Init");
    }
    
    void OnDrawGizmos()
    {
       Gizmos.DrawWireCube(gridCenter, gridSize * cellSize);
    }

}
