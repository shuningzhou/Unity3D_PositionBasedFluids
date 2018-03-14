using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Assets.SpatialPartition.Utils;

public class MySphereCollider{
    public Vector3 v3ContactPoint;
    public Vector3 v3SurfaceNormal;
    private Vector3 v3PrimCenter; //the center of the primitive

    public GameObject Model;

    private float fRadius;
    private float fPntDepth;  //penetration depth

    public MySphereCollider(Vector3 ctr, float r)
    {
        v3PrimCenter = ctr;
        fRadius = r;
    }

    /*
    // Use this for initialization
    void Start()
    {
        //v3PrimCenter = this.transform.position;
        //fRadius = (this.GetComponent<Renderer>().bounds.size.x) / 2;
    }
    */

    public int CheckParticle(Vector3 position)
    {
        Vector3 diffXC = position - v3PrimCenter;
        Vector3 diffCX = -diffXC;
        v3ContactPoint = v3PrimCenter + (fRadius * diffXC.normalized);

        fPntDepth = Mathf.Abs((diffCX).magnitude - fRadius);

        int i = Functions.sgn(SphereFunc(position));

        v3SurfaceNormal = i * diffCX.normalized;
            
        return i;
    }

    float SphereFunc(Vector3 position)
    {
        return (Mathf.Pow((position - v3PrimCenter).magnitude, 2.0f) - (fRadius * fRadius));
    }

    private void Update()
    {
        
    }
}
