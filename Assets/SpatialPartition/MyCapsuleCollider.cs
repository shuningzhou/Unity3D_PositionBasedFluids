using System.Collections;
using System.Collections.Generic;
using UnityEngine;

using Assets.SpatialPartition.Utils;

public class MyCapsuleCollider
{
    private Vector3 v3Center1;
    private Vector3 v3Center2;
    private float radius;

    public Vector3 contactPoint;
    public Vector3 surfaceNormal;

    private GameObject model;

    public MyCapsuleCollider(Vector3 ctr1, Vector3 ctr2, float rad)
    {
        //for now its just axis aligned  -> general oriented later

        v3Center1 = ctr1;
        v3Center2 = ctr2;
        radius = rad;

        model = GameObject.CreatePrimitive(PrimitiveType.Cube);
        model.transform.position = (ctr1 + ctr2) / 2;

        //only vertical capsule for now
        model.transform.localScale = new Vector3(rad * 2, rad + (ctr1 - ctr2).magnitude, rad * 2);

        Material newMat = Resources.Load("MyMat2", typeof(Material)) as Material;
        model.GetComponent<Renderer>().material = newMat;
    }

    public int CollisionCheck(Vector3 ptlPos)
    {
        contactPoint = Vector3.zero;
        surfaceNormal = Vector3.zero;

        float t = 0;
        Vector3 cloestPoint = Vector3.zero;

        t = -1 * (
                    Vector3.Dot((v3Center1 - ptlPos),(v3Center2 - v3Center1)) *
                    Mathf.Pow((v3Center2 - v3Center1).magnitude, 2.0f)
                 );

        t = Mathf.Min(1, Mathf.Max(0, t));

        cloestPoint = v3Center1 + t * (v3Center2 - v3Center1);

        contactPoint = cloestPoint + radius * ((ptlPos - cloestPoint).normalized);

        int i = Functions.sgn(Capsule(ptlPos, cloestPoint));

        surfaceNormal = i * ((cloestPoint - ptlPos).normalized);

        return i;
    }

    public float Capsule(Vector3 ptlPosition, Vector3 cloestPoint)
    {
        return ((cloestPoint - ptlPosition).magnitude - radius);
    }

    public void UpdatePosition(Vector3 distance)
    {
        v3Center1 += distance;
        v3Center2 += distance;

        model.transform.position += distance;
    }
}