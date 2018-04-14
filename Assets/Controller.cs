using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Controller : MonoBehaviour {

    private float mSpeed = 0.2f;
    private float rSpeed = 3f;

    private Vector3 initPosition;
    private Vector3 initRotation;
    private float yaw = 0.0f;
    private float pitch = 0.0f;

    // Use this for initialization
    void Start () {
        //record current position
        initPosition = transform.position;
        initRotation = transform.eulerAngles;

    }
	
	// Update is called once per frame
	void LateUpdate ()
    {
        //reset position
        if (Input.GetKeyUp(KeyCode.R))
        {
            transform.position = initPosition;
            transform.eulerAngles = initRotation;
        }

        //move by key
        if (Input.GetKey(KeyCode.W))
        {
            transform.Translate(Vector3.forward * mSpeed, Space.Self);
        }
        if (Input.GetKey(KeyCode.A))
        {
            transform.Translate(Vector3.left * mSpeed, Space.Self);
        }
        if (Input.GetKey(KeyCode.S))
        {
            transform.Translate(Vector3.back * mSpeed, Space.Self);
        }
        if (Input.GetKey(KeyCode.D))
        {
            transform.Translate(Vector3.right * mSpeed, Space.Self);
        }
        if (Input.GetKey(KeyCode.LeftControl))
        {
            transform.Translate(Vector3.up * mSpeed, Space.Self);
        }
        if (Input.GetKey(KeyCode.LeftShift))
        {
            transform.Translate(Vector3.down * mSpeed, Space.Self);
        }

        if (Input.GetMouseButton(1))
        {
            //rotate by mouse
            yaw = (yaw + rSpeed * Input.GetAxis("Mouse X")) % 360;
            pitch = (pitch - rSpeed * Input.GetAxis("Mouse Y")) % 360;

            transform.eulerAngles = new Vector3(pitch, yaw, 0.0f);
        }

    }
}
