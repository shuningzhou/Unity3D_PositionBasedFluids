using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Assets.SpatialPartition.Utils;

public class SpatialHashTable
{
    private int iHashTableSize;
    private List<MyParticle>[] l_PtlHashTable;
    private float kernelRadius;

    public SpatialHashTable(int numOfParticles, float krnRad)
    {
        kernelRadius = krnRad;
        iHashTableSize = Prime(numOfParticles * 2);
        l_PtlHashTable = new List<MyParticle>[iHashTableSize];
        Debug.Log("Hash table size:  " + iHashTableSize);

#if false
        Prime(2);
        Prime(5);

        Vector3 v3 = new Vector3(90, 110, 30);

        Debug.Log("Hashed index: " + Hash(v3));
#endif
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
        int temp = ((int)(vPosition.x * Consts.iHashP1))
                ^ ((int)(vPosition.y * Consts.iHashP2))
                ^ ((int)(vPosition.z * Consts.iHashP3));

        //Debug.Log("temp = " + temp);

        int index = temp
                    % iHashTableSize;

        //Debug.Log("Discretized position: " + "(" + vPosition.x + "," + vPosition.y + "," + vPosition.z + ")" + "  index: " + index);
        ///Debug.Log("size: " + iHashTableSize);

        return Mathf.Abs(index);//Mathf.Abs(temp);
#if false

        int temp = Mathf.Abs((((int)vPosition.x * Consts.iHashP1))
                ^ (((int)vPosition.y * Consts.iHashP2))
                ^ (((int)vPosition.z * Consts.iHashP3)))
                    % iHashTableSize;

        //Debug.Log("Discretized position: " + "(" + vPosition.x + "," + vPosition.y + "," + vPosition.z + ")" + "  index: " + temp);
        //Debug.Log("size: " + iHashTableSize);

        return temp;//Mathf.Abs(temp);
#endif
    }

    //^ -> xor
    //return (long)(vPosition.x *  Consts.iHashP1) ^ (long)(vPosition.y * Consts.iHashP2) ^ (long)(vPosition.z * Consts.iHashP3) 

    public void AddParticle(MyParticle ptl)
    {
        Vector3 discretizedPos = DiscretizedPoint(ptl.v3Position);
        int idx = Hash(discretizedPos);

        //if the list hasn't be initialized, initializing first
        if (l_PtlHashTable[idx] == null)
        {
            l_PtlHashTable[idx] = new List<MyParticle>();
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

    public void UpdateParticle(MyParticle ptl, Vector3 oldPos)
    {
        Vector3 dOldPos = DiscretizedPoint(oldPos);
        Vector3 dNewPos = DiscretizedPoint(ptl.v3Position);

        //Debug.Log("Old pos: " + dOldPos + "   new pos: " + dNewPos);

        if (dOldPos != dNewPos)
        {
            int idx = Hash(dOldPos);
            l_PtlHashTable[idx].Remove(ptl);

            AddParticle(ptl);
        }
    }

    //flag is used to exclude the original particle from query results
    public List<int> NeighboursQuery(MyParticle ptl, bool bExclusionFlag)
    {
        Vector3 h = new Vector3(kernelRadius, kernelRadius, kernelRadius);

        HashSet<int> hs_result = new HashSet<int>();
        List<Vector3> l_ps = new List<Vector3>();

        //min and max bounding box of sphere
        Vector3 BBmin = DiscretizedPoint(ptl.v3Position - h);
        Vector3 BBmax = DiscretizedPoint(ptl.v3Position + h);

        //Debug.Log("BBmin:  " + "(" + BBmin.x + "," + BBmin.y + "," +BBmin.z + ")");
        //Debug.Log("BBmax:  " + "(" + BBmax.x + "," + BBmax.y + "," + BBmax.z + ")");

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

        List<MyParticle> l_temp;

        foreach (Vector3 v3 in l_ps)
        {
            int index = Hash(v3);
                
            if (index >= 0 && index < l_PtlHashTable.Length)
            {
                //Debug.Log(index);

                l_temp = l_PtlHashTable[index];
                if (l_temp != null && l_temp.Count != 0)  //safe because of lazy evaluation
                {
                    foreach (MyParticle p in l_temp)
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

        return (new List<int>(hs_result));
    }
}
