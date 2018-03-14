using System.Collections;
using System.Collections.Generic;
using UnityEngine;

using Assets.SpatialPartition.Utils;

public class MyBoxCollider
{
    private Vector3 v3Center;
    private Vector3 v3Extent;

    public Vector3 contactPoint;
    public Vector3 surfaceNormal;



    private GameObject model;

    public MyBoxCollider(Vector3 ctr, Vector3 sideLgh)
    {
        //for now its just axis aligned  -> general oriented later

        v3Center = ctr;
        v3Extent = sideLgh;

        model = GameObject.CreatePrimitive(PrimitiveType.Cube);
        model.transform.position = ctr;

        //only vertical capsule for now
        model.transform.localScale = (sideLgh * 2);

        Material newMat = Resources.Load("MyMat2", typeof(Material)) as Material;
        model.GetComponent<Renderer>().material = newMat;


    }

    public int CollisionCheck(Vector3 ptlPos)
    {
        contactPoint = Vector3.zero;
        surfaceNormal = Vector3.zero;

        Vector3 pLocal = ptlPos - v3Center;

        Vector3 cpLocal = VMin(v3Extent, VMax((-1 * v3Extent), pLocal));

        contactPoint = v3Center + cpLocal;

        surfaceNormal = VSgn(cpLocal - pLocal).normalized;

        int i = Functions.sgn(VGetMax(VAbs(pLocal) - v3Extent));

        return i;
    }

    public Vector3 VMax(Vector3 v1, Vector3 v2)
    {
        Vector3 max = Vector3.zero;

        max.x = Mathf.Max(v1.x, v2.x);
        max.y = Mathf.Max(v1.y, v2.y);
        max.z = Mathf.Max(v1.z, v2.z);

        return max;
    }

    public Vector3 VMin(Vector3 v1, Vector3 v2)
    {
        Vector3 min = Vector3.zero;

        min.x = Mathf.Min(v1.x, v2.x);
        min.y = Mathf.Min(v1.y, v2.y);
        min.z = Mathf.Min(v1.z, v2.z);

        return min;
    }

    public Vector3 VSgn(Vector3 v)
    {
        Vector3 result = Vector3.zero;

        result.x = Functions.sgn(v.x);
        result.y = Functions.sgn(v.y);
        result.z = Functions.sgn(v.z);

        return result;
    }

    public float VGetMax(Vector3 v)
    {
        return Mathf.Max(v.x, Mathf.Max(v.y, v.z));
    }

    public Vector3 VAbs(Vector3 v)
    {
        return (new Vector3(Mathf.Abs(v.x), Mathf.Abs(v.y), Mathf.Abs(v.z)));
    }

    public void UpdatePosition(Vector3 distance)
    {
        v3Center += distance;

        model.transform.position += distance;
    }
}
