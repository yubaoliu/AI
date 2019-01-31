using UnityEngine;
using System.Collections;
using System.IO;
using System;
using System.Net.Sockets;
using System.Threading;
using UnityEngine.Networking;
using System.Net;
using System.Text;

public class cameraMove : MonoBehaviour
{

    int sequence = 0;
    int oldsequence = 0;
    static int TotalFrameNumber = 1501;
    Matrix4x4[] cameraPos = new Matrix4x4[1501];

    void Start()
    {
        readFrame(System.Environment.CurrentDirectory + "/Assets/data/cam_outdoor.txt");
    }
    int calculateCameraRotateTransformFromProjectionMatrix(Matrix4x4 P, ref Matrix4x4 cameraRotation, ref Vector4 cameraPosition)//P is the projection matrix
    {
        Vector4 T = new Vector4(P.m03, P.m13, P.m23, P.m33);
        Matrix4x4 R = Matrix4x4.zero;
        R.m33 = 1.0f;
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                R[i + j * 4] = P[i + j * 4];
            }
        }
        cameraRotation = R.inverse;
        cameraPosition = R.inverse * T * (-1.0f);
        return 0;
    }
    int computingEulerAnglesFromaRotationMatrix(Matrix4x4 R, ref Vector4 eulerangle)
    {
        float R11 = R.m00;
        float R12 = R.m01;
        float R13 = R.m02;
        float R21 = R.m10;
        float R22 = R.m11;
        float R23 = R.m12;
        float R31 = R.m20;
        float R32 = R.m21;
        float R33 = R.m22;
        float thet1 = 0f, thet2 = 0f, psi1 = 0f, psi2 = 0f, phi1 = 0f, phi2 = 0f;
        float thet = 0f, phi = 0f, psi = 0;

        if (!(R31 + 1.0f < 0.0000000001f) || !(R31 - 1 < 0.0000000001f))
        {
            thet1 = (-1.0f) * Mathf.Asin(R31);
            thet2 = Mathf.PI - thet1;

            psi1 = Mathf.Atan2(R32 / Mathf.Cos(thet1), R33 / Mathf.Cos(thet1));
            psi2 = Mathf.Atan2(R32 / Mathf.Cos(thet2), R33 / Mathf.Cos(thet2));

            phi1 = Mathf.Atan2(R21 / Mathf.Cos(thet1), R11 / Mathf.Cos(thet1));
            phi2 = Mathf.Atan2(R21 / Mathf.Cos(thet2), R11 / Mathf.Cos(thet2));
            //both 1 and 2 are valid, here use 1
            eulerangle.x = Mathf.Rad2Deg * psi1;
            eulerangle.y = Mathf.Rad2Deg * thet1;
            eulerangle.z = Mathf.Rad2Deg * phi1;
        }
        else
        {
            phi = 0;
            if (R31 + 1.0f < 0.0000000001f)
            {
                thet = Mathf.PI / 2;
                psi = phi + Mathf.Atan2(R12, R13);
            }
            else
            {
                thet = (-1.0f) * Mathf.PI / 2.0f;
                psi = (-1.0f) * psi + Mathf.Atan2((-1.0f) * R12, (-1.0f) * R13);
            }
            eulerangle.x = Mathf.Rad2Deg * psi;
            eulerangle.y = Mathf.Rad2Deg * thet;
            eulerangle.z = Mathf.Rad2Deg * phi;
        }

        return 0;
    }

    int frameIndex = 0;
    // Update is called once per frame
    void Update()
    {
        Matrix4x4 projectionMatrix = cameraPos[++frameIndex % TotalFrameNumber];
        Matrix4x4 cameraRotation = Matrix4x4.zero;
        Vector4 cameraPosition = Vector4.zero;
        calculateCameraRotateTransformFromProjectionMatrix(projectionMatrix, ref cameraRotation, ref cameraPosition);

        int scale = 1;
        int basement = 50;
        transform.position = new Vector3(scale * cameraPosition.x, basement+scale * cameraPosition.y, -scale * cameraPosition.z);

        Vector4 eulerangle = Vector4.zero;
        computingEulerAnglesFromaRotationMatrix(cameraRotation, ref eulerangle);

        transform.rotation = new Quaternion(-eulerangle.x, -eulerangle.y, eulerangle.z, 0);
        transform.Rotate(Vector3.right, 180);

    }
    int readFrame(string filePath)
    {
        StreamReader sr = new StreamReader(filePath, Encoding.Default);
        string line;
        Matrix4x4 RT = Matrix4x4.zero;
        while ((line = sr.ReadLine()) != null)
        {
            sequence = Convert.ToInt32(line);
            Debug.Log("sequence: " + sequence);
            for (int i = 0; i < 3; i++)
            {
                string[] data = new string[4];
                line = sr.ReadLine();
                data = line.Split(' ');
                for (int j = 0; j < 4; j++)
                {
                    RT[i + j * 4] = Convert.ToSingle(data[j]);
                }
            }
            RT.m33 = 1;
            cameraPos[sequence] = RT;
        }
        return 0;
    }
}
