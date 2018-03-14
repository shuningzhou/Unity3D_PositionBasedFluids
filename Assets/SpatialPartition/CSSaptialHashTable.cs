using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Assets.SpatialPartition;

#if false
public class CSSpatialHashTable
{
    private const int iHashP1 = 73856093;
    private const int iHashP2 = 19349663;
    private const int iHashP3 = 83492791;

    private int iHashTableSize;
    private List<CSParticle>[] l_PtlHashTable;
    private float kernelRadius;

    public CSSpatialHashTable(int numOfParticles, float krnRad)
    {
        kernelRadius = krnRad;
        iHashTableSize = Prime(numOfParticles * 2);
        l_PtlHashTable = new List<CSParticle>[iHashTableSize];
        Debug.Log("Hash table size:  " + iHashTableSize);

    }

    public int Prime(int n)
    {
        //meaningless, just for assertion
        //int origin = n;

        int t = 0, r = 0, c, d;

        while (true)
        {
            n++;
            t = n;

            /* Calculating reverse of number */
            while (t != 0)
            {
                r *= 10;  /* Compound assignment operator r*=10 => r=r*10 */
                r += t % 10;
                t /= 10;
            }

            /* if reverse equals original then it is palindrome */
            if (r == n)
            {
                d = (int)Mathf.Sqrt(n);

                /* Checking prime */
                for (c = 2; c <= d; c++)
                {
                    if (n % c == 0)
                        break;
                }
                if (c == d + 1)
                    break;
            }
            r = 0;
        }
        //Debug.Log("origin: " + origin);
        //Debug.Log("prime: " + n);
        return n;
    }

    private Vector3 DiscretizedPoint(Vector3 vOPosition)
    {
        return (new Vector3(Mathf.Floor(vOPosition.x / kernelRadius), Mathf.Floor(vOPosition.y / kernelRadius), Mathf.Floor(vOPosition.z / kernelRadius)));
    }

    private int Hash(Vector3 vPosition)
    {
        int temp = ((int)(vPosition.x * iHashP1))
                ^ ((int)(vPosition.y * iHashP2))
                ^ ((int)(vPosition.z * iHashP3));

        //Debug.Log("temp = " + temp);

        int index = temp
                    % iHashTableSize;

        return Mathf.Abs(index);//Mathf.Abs(temp);
    }

    public void AddP(CSParticle ptl)
    {
        Vector3 discretizedPos = DiscretizedPoint(ptl.v3PredictedPos);
        int idx = Hash(discretizedPos);

        //if the list hasn't be initialized, initializing first
        if (l_PtlHashTable[idx] == null)
        {
            l_PtlHashTable[idx] = new List<CSParticle>();
        }

        l_PtlHashTable[idx].Add(ptl);

#if false
        Debug.Log("Size: " + l_PtlHashTable[idx].Count);
        foreach (MyParticle p in l_PtlHashTable[idx])
        {
            Debug.Log("Pos: " + p.v3Position.x + " " + p.v3Position.y + " " + p.v3Position.z);
        }
#endif
    }

    private void AddParticles(CSParticle[] ps)
    {
        foreach (CSParticle p in ps)
        {
            AddP(p);
        }
    }

    private void Clear()
    {
        foreach (List<CSParticle> ps in l_PtlHashTable)
        {
            if (ps != null)
                ps.Clear();
        }
    }

    public void NeighboursQuery(CSParticle[] ps, ref uint[] nList, ref Vector2[] nIndex, bool bExclusionFlag)
    {
        if (l_PtlHashTable.Length != 0)
        {
            Clear();
        }

        AddParticles(ps);

        //uint[][] result = new uint[ps.Length][];

        List<uint> neighboursList = new List<uint>();
        Vector2[] neighboursIndex = new Vector2[ps.Length];

        GetNeighbours(ps, neighboursList, neighboursIndex, bExclusionFlag);

        nList = neighboursList.ToArray();
        nIndex = neighboursIndex;

        //Debug.Log("nList size: " + nList.Length);
        //Debug.Log("nIndex size: " + nIndex.Length);
    }

    //flag is used to exclude the original particle from query results
    public void GetNeighbours(CSParticle[] ps, List<uint> nList, Vector2[] nIndex, bool bExclusionFlag)
    {
        Vector3 h = new Vector3(kernelRadius, kernelRadius, kernelRadius);

        List<Vector3> l_ps;

        //min and max bounding box of sphere
        Vector3 BBmin;
        Vector3 BBmax;

        List<CSParticle> l_temp;

        HashSet<uint> hs_result;

        uint iIndexCounter = 0;

        foreach (CSParticle ptl in ps)
        {
            l_ps = new List<Vector3>();

            //min and max bounding box of sphere
            BBmin = DiscretizedPoint(ptl.v3Position - h);
            BBmax = DiscretizedPoint(ptl.v3Position + h);

            for (int i = (int)BBmin.x; i <= (int)BBmax.x; i++)
            {
                for (int j = (int)BBmin.y; j <= (int)BBmax.y; j++)
                {
                    for (int k = (int)BBmin.z; k <= (int)BBmax.z; k++)
                    {
                        l_ps.Add(new Vector3(i, j, k));
                    }
                }
            }

            hs_result = new HashSet<uint>();

            foreach (Vector3 v3 in l_ps)
            {
                int index = Hash(v3);

                if (index >= 0 && index < l_PtlHashTable.Length)
                {
                    //Debug.Log(index);

                    l_temp = l_PtlHashTable[index];
                    if (l_temp != null && l_temp.Count != 0)  //safe because of lazy evaluation
                    {
                        foreach (CSParticle p in l_temp)
                        {
                            //Debug.Log("distance:  " + (ptl.v3Position - p.v3Position).magnitude);
                            if ((ptl.v3Position - p.v3Position).magnitude <= kernelRadius)
                            {
                                if (bExclusionFlag && (ptl.id == p.id))
                                {
                                    continue;
                                }

                                //Debug.Log("Added");
                                hs_result.Add(p.id);

                            }
                        }
                    }
                }
            }

            nList.AddRange(hs_result);
            nIndex[ptl.id].x = iIndexCounter;
            nIndex[ptl.id].y = hs_result.Count;

            iIndexCounter += (uint)hs_result.Count;

            /*legacy test
            if (ptl.id == 0)
            {
                Debug.Log("------------------------------------------");

                Debug.Log("index:  " + iIndexCounter);

                foreach (uint idx in nList)
                {
                    Debug.Log("ID: " + idx);
                }

                Debug.Log("------------------------------------------");
            }
            */
        }
    }
}

#if false
public class CSSpatialHashTable
{
    private const int iHashP1 = 73856093;
    private const int iHashP2 = 19349663;
    private const int iHashP3 = 83492791;

    private int iHashTableSize;
    private List<CSParticle>[] l_PtlHashTable;
    private float kernelRadius;

    public CSSpatialHashTable(int numOfParticles, float krnRad)
    {
        kernelRadius = krnRad;
        iHashTableSize = Prime(numOfParticles * 2);
        l_PtlHashTable = new List<CSParticle>[iHashTableSize];
        Debug.Log("Hash table size:  " + iHashTableSize);

    }

    public int Prime(int n)
    {
        //meaningless, just for assertion
        //int origin = n;

        int t = 0, r = 0, c, d;

        while (true)
        {
            n++;
            t = n;

            /* Calculating reverse of number */
            while (t != 0)
            {
                r *= 10;  /* Compound assignment operator r*=10 => r=r*10 */
                r += t % 10;
                t /= 10;
            }

            /* if reverse equals original then it is palindrome */
            if (r == n)
            {
                d = (int)Mathf.Sqrt(n);

                /* Checking prime */
                for (c = 2; c <= d; c++)
                {
                    if (n % c == 0)
                        break;
                }
                if (c == d + 1)
                    break;
            }
            r = 0;
        }
        //Debug.Log("origin: " + origin);
        //Debug.Log("prime: " + n);
        return n;
    }

    private Vector3 DiscretizedPoint(Vector3 vOPosition)
    {
        Vector3 posWithoutOffset = new Vector3(vOPosition.x - 10, vOPosition.y - 10, vOPosition.z - 10);
        return (new Vector3(Mathf.Floor(posWithoutOffset.x / kernelRadius), Mathf.Floor(posWithoutOffset.y / kernelRadius), Mathf.Floor(posWithoutOffset.z / kernelRadius)));
    }

    private int Hash(Vector3 vPosition)
    {
        int temp = ((int)(vPosition.x * iHashP1))
                ^ ((int)(vPosition.y * iHashP2))
                ^ ((int)(vPosition.z * iHashP3));

        //Debug.Log("temp = " + temp);

        int index = temp
                    % iHashTableSize;

        return Mathf.Abs(index);//Mathf.Abs(temp);
    }

    public void AddP(CSParticle ptl)
    {
        Vector3 discretizedPos = DiscretizedPoint(ptl.v3PredictedPos);
        int idx = Hash(discretizedPos);

        //if the list hasn't be initialized, initializing first
        if (l_PtlHashTable[idx] == null)
        {
            l_PtlHashTable[idx] = new List<CSParticle>();
        }

        l_PtlHashTable[idx].Add(ptl);

#if false
        Debug.Log("Size: " + l_PtlHashTable[idx].Count);
        foreach (MyParticle p in l_PtlHashTable[idx])
        {
            Debug.Log("Pos: " + p.v3Position.x + " " + p.v3Position.y + " " + p.v3Position.z);
        }
#endif
    }

    public void AddParticles(CSParticle[] ps)
    {
        foreach (CSParticle p in ps)
        {
            AddP(p);
        }
    }

    public void Clear()
    {
        foreach (List<CSParticle> ps in l_PtlHashTable)
        {
            ps.Clear();
        }
    }

    public uint[][] NeighboursQuery(CSParticle[] ps, bool bExclusionFlag)
    {
        uint[][] result = new uint[ps.Length][];

        List<uint> temp;

        foreach (CSParticle p in ps)
        {
            temp = GetNeighbours(p, bExclusionFlag);

            result[p.id] = temp.ToArray();
        }

        return result;
    }

    //flag is used to exclude the original particle from query results
    public List<uint> GetNeighbours(CSParticle ptl, bool bExclusionFlag)
    {
        Vector3 h = new Vector3(kernelRadius, kernelRadius, kernelRadius);

        HashSet<uint> hs_result = new HashSet<uint>();
        List<Vector3> l_ps = new List<Vector3>();

        //min and max bounding box of sphere
        Vector3 BBmin = DiscretizedPoint(ptl.v3Position - h);
        Vector3 BBmax = DiscretizedPoint(ptl.v3Position + h);

        for (int i = (int)BBmin.x; i <= (int)BBmax.x; i++)
        {
            for (int j = (int)BBmin.y; j <= (int)BBmax.y; j++)
            {
                for (int k = (int)BBmin.z; k <= (int)BBmax.z; k++)
                {
                    l_ps.Add(new Vector3(i, j, k));
                }
            }
        }

        List<CSParticle> l_temp;

        foreach (Vector3 v3 in l_ps)
        {
            int index = Hash(v3);

            if (index >= 0 && index < l_PtlHashTable.Length)
            {
                //Debug.Log(index);

                l_temp = l_PtlHashTable[index];
                if (l_temp != null && l_temp.Count != 0)  //safe because of lazy evaluation
                {
                    foreach (CSParticle p in l_temp)
                    {
                        //Debug.Log("distance:  " + (ptl.v3Position - p.v3Position).magnitude);
                        if ((ptl.v3Position - p.v3Position).magnitude <= kernelRadius)
                        {
                            if (bExclusionFlag && (ptl.id == p.id))
                            {
                                continue;
                            }

                            //Debug.Log("Added");
                            hs_result.Add(p.id);
                        }
                    }
                }
            }
        }

        return (new List<uint>(hs_result));
    }
}
#endif
#endif