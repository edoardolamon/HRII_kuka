/** 
 \author (Edoardo Lamon)
	\file
    \brief Planning and Implementation of a Cartesian Impedence Control

	Get this one running on your favorite system as proof of concept
NOTE: This sample, as the corresponding FRI (Fast Research inteface) is subject to radical change

 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef __WIN32
#include "stdafx.h"
#endif



#include <iostream>
#include <cstdlib>
#include <cmath>
#include <limits.h>
#include "friudp.h"
#include "friremote.h"

#include <vector>
#include <string>
#include <eigen3/Eigen/Dense>

#ifndef M_PI 
#define M_PI 3.14159
#endif

using namespace std;

float _degToRad(float degAngle)
{
    return degAngle*M_PI/180;
}

void _computeCubicTrajectory(vector<float> &a, float finalJntPos, float initialJntPos, float avgJntVel=0.174, float finalJntVel=0., float initialJntVel=0.)
{
    //float finalTime = abs((finalJntPos - initialJntPos) / avgJntVel);
    float finalTime = 5.;
    cout << "Planning time: " << finalTime << " seconds" << endl;

    a.push_back(initialJntPos); //a[0]
    a.push_back(initialJntVel); //a[1]

    float a2 = (-3*(initialJntPos - finalJntPos) - (2*initialJntVel + finalJntVel)) / (pow(finalTime, 2));
    a.push_back(a2);

    float a3 = (2*(initialJntPos - finalJntPos) + (initialJntVel + finalJntVel)) / (pow(finalTime, 3));
    a.push_back(a3);
}


float _nextStepCubic(vector<float> &cubicCoeffs, float t)
{
    return cubicCoeffs[3]*pow(t,3.) + cubicCoeffs[2]*pow(t,2.) + cubicCoeffs[1]*t + cubicCoeffs[0]; // static_cast<float>
}

void _nextStepCircular(Eigen::Vector3f &nextPos, float t, float period, float radius, vector<float> &center)
{
    float omega = (2*M_PI)/(period);
    float x = - radius*cos(omega*t) + center[0];
    float y = radius*sin(omega*t) + center[1];

    nextPos(0) = x;
    nextPos(1) = y;
}

void _calcKukaJacob(Eigen::MatrixXf &J_kuka, float jointPos[])
{
    if (!(J_kuka.rows() == FRI_CART_VEC && J_kuka.cols() == LBR_MNJ))  // (6x7)
    {
        cerr << "Jacobian wrong dimensions. Jacobian should be 6x7" << endl;
        J_kuka(FRI_CART_VEC,LBR_MNJ);
    }
//    if (!(jointPos.size() == LBR_MNJ))
//    {
//        cerr << "Joint vector wrong dimensions. It should be 6x1" << endl;
//    }

    float C1,S1,S2,C2,S3,C3,S4,C4,S5,C5,S6,C6,S7,C7;
    float l1=0.310, l2=0.4, l3=0.39, l4=0.075;

    C1=cos(jointPos[0]);
    S1=sin(jointPos[0]);
    C2=cos(jointPos[1]+M_PI/2);
    S2=sin(jointPos[1]+M_PI/2);
    C3=cos(jointPos[2]);
    S3=sin(jointPos[2]);
    C4=cos(jointPos[3]);
    S4=sin(jointPos[3]);
    C5=cos(jointPos[4]);
    S5=sin(jointPos[4]);
    C6=cos(jointPos[5]);
    S6=sin(jointPos[5]);
    C7=cos(jointPos[6]);
    S7=sin(jointPos[6]);

    J_kuka.row(0) << - l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3))) - l2*C2*S1, -C1*(l3*(C4*S2 - C2*C3*S4) + l2*S2 + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))) , C2*S1*(l3*(C4*S2 - C2*C3*S4) + l2*S2 + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))) - S2*(l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3))) + l2*C2*S1) , (C1*C3 - S1*S2*S3)*(l3*(C4*S2 - C2*C3*S4) + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))) - C2*S3*(l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3)))) , (S4*(C1*S3 + C3*S1*S2) + C2*C4*S1)*(l3*(C4*S2 - C2*C3*S4) + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))) - (l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3))))*(C4*S2 - C2*C3*S4) , l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))*(S5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) - C5*(C1*C3 - S1*S2*S3)) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3)))*(S5*(S2*S4 + C2*C3*C4) + C2*C5*S3) , 0;
    J_kuka.row(1) << l2*C1*C2 - l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3))) - l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4), -S1*(l3*(C4*S2 - C2*C3*S4) + l2*S2 + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))), - S2*(l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3))) - l2*C1*C2) - C1*C2*(l3*(C4*S2 - C2*C3*S4) + l2*S2 + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))), (C3*S1 + C1*S2*S3)*(l3*(C4*S2 - C2*C3*S4) + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))) - C2*S3*(l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3)))), (l3*(C4*S2 - C2*C3*S4) + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4)))*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - (C4*S2 - C2*C3*S4)*(l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3)))), l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))*(S5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) - C5*(C3*S1 + C1*S2*S3)) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3)))*(S5*(S2*S4 + C2*C3*C4) + C2*C5*S3), 0;
    J_kuka.row(2) << 0 , S1*(l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3))) + l2*C2*S1) - C1*(l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3))) - l2*C1*C2), C2*S1*(l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3))) - l2*C1*C2) + C1*C2*(l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3))) + l2*C2*S1), (C1*C3 - S1*S2*S3)*(l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3)))) - (l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3))))*(C3*S1 + C1*S2*S3), (S4*(C1*S3 + C3*S1*S2) + C2*C4*S1)*(l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3)))) - (l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3))))*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4), l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3)))*(S5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) - C5*(C1*C3 - S1*S2*S3)) - l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3)))*(S5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) - C5*(C3*S1 + C1*S2*S3)), 0;
    J_kuka.row(3) << 0 , S1 , C1*C2 , - C3*S1 - C1*S2*S3, C1*C2*C4 - S4*(S1*S3 - C1*C3*S2), C5*(C3*S1 + C1*S2*S3) - S5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4), S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3)) - C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4);
    J_kuka.row(4) << 0 , -C1 , C2*S1 , C1*C3 - S1*S2*S3, S4*(C1*S3 + C3*S1*S2) + C2*C4*S1, S5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) - C5*(C1*C3 - S1*S2*S3), C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3));
    J_kuka.row(5) << 1 , 0 , S2 , C2*S3 , C4*S2-C2*C3*S4 , - S5*(S2*S4+C2*C3*C4) - C2*C5*S3, S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4);
}

int
#ifdef __WIN32

_tmain
#else
#ifdef _WRS_KERNEL
friFirstApp
#else
main
#endif
#endif

(int argc, char *argv[])
{

	cout << "Opening FRI Version " 
		<< FRI_MAJOR_VERSION << "." << FRI_SUB_VERSION << "." <<FRI_DATAGRAM_ID_CMD << "." <<FRI_DATAGRAM_ID_MSR 
		<< " Interface for First Sample" << endl;
	{
		// do checks, whether the interface - and the host meets the requirements
		// Note:: This Check remains in friRempte.cpp -- should go to your code ...
		FRI_PREPARE_CHECK_BYTE_ORDER;
		if (!FRI_CHECK_BYTE_ORDER_OK) 
		{
			cerr << "Byte order on your system is not appropriate - expect deep trouble" <<endl;
		}
		if (!FRI_CHECK_SIZES_OK)
		{
			cout << "Sizes of datastructures not appropriate - expect even deeper trouble" << endl;

		}
	}
        // GENERAL VARIABLES
//        friRemote friInst(49948,"192.168.0.10"); // HRII Kuka
        friRemote friInst(49938,"192.168.0.20");  // ADVR Kuka
		FRI_QUALITY lastQuality = FRI_QUALITY_BAD;
        FRI_CTRL lastCtrlScheme = FRI_CTRL_OTHER;  // added

        bool firstLoop = true;
        double timeCounter = 0.;
        const float eps = 0.0001;

/////// JOINT IMPEDENCE VARIABLES
//        float myJntStiffness[LBR_MNJ]={100,100,100,100,100,100,100};
        float nullJntStiffness[LBR_MNJ]={0,0,0,0,0,0,0};
        float newJntVals[LBR_MNJ];
        float msrdJntVals[LBR_MNJ];
        float desJntVals[LBR_MNJ] = {-0.04, -0.58, 0.07, 2., 0.1, 0.99, 0.77};

        vector< vector<float> > cubicJntCoeffs(LBR_MNJ);
//        vector< vector<float> > cubicJntCoeffs;
//        cout << "Test cubicJntCoeffs dimension: " << cubicJntCoeffs[0].size() << endl;

//////// CARTESIAN IMPEDNCE VARIABLES
//        float myCartStiffness[FRI_CART_VEC]={1500.0,1500.0,0.0,150.0,150.0,150.0};
        float nullCartStiffness[FRI_CART_VEC]={0.0,0.0,0.0,0.0,0.0,0.0};
        float desCartVals[FRI_CART_FRM_DIM]; //= {-0.781831, 0.279979, 0.557092, 0.574816, 0.331044, 0.943566};
        float newCartVals[FRI_CART_FRM_DIM], msrdCartVals[FRI_CART_FRM_DIM];
        Eigen::VectorXf newCartVector(FRI_CART_FRM_DIM), msrdCartVector(FRI_CART_FRM_DIM);
        Eigen::MatrixXf jacobian(FRI_CART_VEC,LBR_MNJ);
        Eigen::Affine3f newCartAffine;
        Eigen::Matrix3f newCartRotation;
        Eigen::Vector3f newCartTranslation;
        Eigen::Vector3f desCartTranslation(0.1, 0.0, 0.0);

        vector< vector<float> > cubicCartCoeffs(3);
        const float radius = 0.05;
        const float period = 5.;
        vector<float> center(2);
        Eigen::VectorXf externalForce(FRI_CART_VEC);
                        externalForce << 0.0,0.0,5.0,0.0,0.0,0.0;
        float externalForceVector[FRI_CART_VEC] = {0.0,0.0,5.0,0.0,0.0,0.0};
//        float externalDeltaForceVector[FRI_CART_VEC] = {0.0,0.0,0.001,0.0,0.0,0.0};
        Eigen::VectorXf externalTorque(LBR_MNJ);

        printf("START CONTROL\n");
        friInst.doDataExchange();
        double deltaT = friInst.getSampleTime();
        printf("SAMPLE TIME: %f \n", deltaT);
        cout << endl;

//////* enter main loop - wait until we enter stable command mode */
        for(;;)
		{
			/// perform some arbitrary handshake to KRL -- possible in monitor mode already
			// send to krl int a value
			friInst.setToKRLInt(0,1);
			if ( friInst.getQuality() >= FRI_QUALITY_OK)
			{
				// send a second marker
				friInst.setToKRLInt(0,10);
			}

			//
            // just mirror the real value

			friInst.setToKRLReal(0,friInst.getFrmKRLReal(1));

            /// Connection quality check
            if ( friInst.getQuality() != lastQuality)
            {
                cout << endl;
                cout << "Quality change detected "<< friInst.getQuality()<< " \n";
                cout << friInst.getMsrBuf().intf;

                lastQuality = friInst.getQuality();
            }

            switch (friInst.getCurrentControlScheme())
            {

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// JOINT IMPEDANCE CONTROL
                case FRI_CTRL_JNT_IMP:

                // Read the current joint position
                for (int i = 0; i < LBR_MNJ; i++)
                {
                    newJntVals[i] = friInst.getMsrCmdJntPosition()[i];
                    msrdJntVals[i] = friInst.getMsrMsrJntPosition()[i];
                }

                // Read the current cartesian position
                for (int i = 0; i < FRI_CART_FRM_DIM; i++)
                {
                    newCartVals[i] = friInst.getMsrCmdCartPosition()[i];
                    msrdCartVals[i] = friInst.getMsrCartPosition()[i];

                    int row = i / 4;
                    int col = i % 4;
                    newCartAffine(row,col) = friInst.getMsrCmdCartPosition()[i];
                }

                /// First loop computation
                if (firstLoop)
                {

                    firstLoop = false;
                }

                if ( friInst.getState() == FRI_STATE_CMD)
                {
                    if ( friInst.isPowerOn() )
                    {
                        // Update time counter
                        timeCounter += deltaT;

                        // Compute jacobian
//                        _calcKukaJacob(jacobian,newJntVals);
//                        cout << "Jacobian = " << jacobian << endl;
                        externalTorque = jacobian.transpose() * externalForce;
//                        cout << "Input torque : \n" << externalTorque.transpose() << endl;

    //					for (int i = 0; i < LBR_MNJ; i++)
    //					{
    //                        msrdJntVals[i] = friInst.getMsrCmdJntPosition()[i];

    //                        if(abs(desJntVals[i] - msrdJntVals[i]) > eps)
    //                        if ( abs(desJntVals[i] - newJntVals[i]) > eps)
    //                        {
    //                            newJntVals[i] = _nextStepCubic(cubicJntCoeffs[i],timeCounter);  // compute trajectory
    //                        }

    //                        if (i=0)
    //                        {
    //                            cout << "Joint " << i <<" next desired position: "<< newJntVals[i] << endl;
    //                        }

                            // perform some sort of sine wave motion

    //						newJntVals[i]+=(float)sin( timeCounter * M_PI * 0.02) * (float)(10./180.*M_PI);
    //					}
//                        newJntVals[2] = _nextStepArmonic(desJntVals[2],msrdJntVals[2],timeCounter);
                    }
                    else
                    {
                        timeCounter = 0.;
                    }
                }
                else
                {
                    timeCounter = 0.;
                }

                friInst.doJntImpedanceControl(NULL, NULL, NULL, NULL);
//                friInst.doJntImpedanceControl(newJntVals, myJntStiffness, NULL, NULL);
//                  friInst.doJntImpedanceControl(newJntVals, myJntStiffness, NULL, externalTorque.data());
                break;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CARTESIAN IMPEDANCE CONTROL
                case FRI_CTRL_CART_IMP:

                for (int i = 0; i < LBR_MNJ; i++)
                {
                    newJntVals[i] = friInst.getMsrCmdJntPosition()[i];
                }

                /// First loop computation
                if (firstLoop)
                {
                    // Read the current cartesian position
                    for (int i = 0; i < FRI_CART_FRM_DIM; i++)
                    {
                        newCartVals[i] = friInst.getMsrCmdCartPosition()[i];
                        msrdCartVals[i] = friInst.getMsrCartPosition()[i];

                        int row = i / 4;
                        int col = i % 4;
                        newCartAffine(row,col) = friInst.getMsrCmdCartPosition()[i];
                    }

                    newCartRotation = newCartAffine.linear();
                    newCartTranslation = newCartAffine.translation();
                    desCartTranslation = desCartTranslation + newCartTranslation;

//                    cout << "Actual cartesian position: " << newCartVals[3] << " , " << newCartVals[7] << " , " << newCartVals[11] << endl;
//                    cout << "Measured cartesian position: " << msrdCartVals[3] << " , " << msrdCartVals[7] << " , " << msrdCartVals[11] << endl;

//                    cout << "Desired cartesian translation: " << " : " << desCartTranslation(0) << " , " << desCartTranslation(1) << " , " << desCartTranslation(2) << endl;

                    // Compute jacobian
                    _calcKukaJacob(jacobian,newJntVals);
                    cout << "Jacobian = " << jacobian << endl;
//                    cout << "Jacobian determinant = " << jacobian.determinant() << endl;

                    // Planning coeffs
//                    for (int i = 0; i < 3; i++)
//                    {
//                       _computeCubicTrajectory(cubicCartCoeffs[i], desCartTranslation(i), newCartTranslation(i));
//                        cout << "Trajectory coeffs of translation " << i << " : a0 = " << cubicCartCoeffs[0][i]
//                             << ", a1 = " << cubicCartCoeffs[i][1] << ", a2 = " << cubicCartCoeffs[i][2] << ", a3 = " << cubicCartCoeffs[i][3] << endl;
//                        cout << endl;
//                    }

//                    cout << newCartTranslation(0) << endl;
//                    cout << newCartTranslation(1) << endl;
                    center[0] = newCartTranslation(0) + radius;
                    center[1] = newCartTranslation(1);

                    firstLoop = false;
                    cout << "Plan computed" << endl;
                }

                if ( friInst.getState() == FRI_STATE_CMD)
                {
                    if ( friInst.isPowerOn() )
                    {
                        timeCounter += deltaT;
//
//                            vector<float> tmp_CartVals(newCartVals, newCartVals+6);     !!!!!!!!!
//                            _nextStepCircular(tmp_CartVals,timeCounter,period,radius);
//                            copy(tmp_CartVals.begin(),tmp_CartVals.end(),newCartVals);

                        // PLAN CARTESIAN POSITION TRAJECTORY
//                        for (int i = 0; i < 3; i++)
//                        {
//                            if ( abs(desCartTranslation(i) - newCartTranslation(i)) > eps)
//                            {
//                                newCartTranslation(i) = _nextStepCubic(cubicCartCoeffs[i],timeCounter);
//                            }
//                        }

                        // PLAN CARTESIAN CIRCLE TRAJECTORY
//                        _nextStepCircular(newCartTranslation,timeCounter,period,radius,center);

                        // CALCULATE JACOBIAN AND REBUILD THE INITIAL CARTESIAN VECTOR
                        _calcKukaJacob(jacobian,newJntVals);
//                        cout << "Jacobian = " << jacobian << endl;

                        newCartAffine.translation() = newCartTranslation;
                        for (int i = 0; i < FRI_CART_FRM_DIM; i++)
                        {
                            int row = i / 4;
                            int col = i % 4;
                            newCartVals[i] = newCartAffine(row,col);
                        }
//                        cout << "New cartesian position: " << newCartVals[3] << " , " << newCartVals[7] << " , " << newCartVals[11] << endl;
//                        cout << endl;

//                        externalTorque = jacobian.transpose() * externalForce;
//                        friInst.doCartesianImpedanceControl(NULL, NULL, NULL, externalTorque.data());
                    }
                    else
                    {
                        timeCounter = 0.;
                    }
                }
                else
                {
                    timeCounter = 0.;
                }

//                for (int i = 0; i < FRI_CART_VEC; i++)
//                {
//                    externalForceVector[i] += externalDeltaForceVector[i];
//                }

//                cout << "Input force : \n" << externalForce.transpose() << endl;


//                friInst.doCartesianImpedanceControl(NULL, NULL, NULL, NULL);
//                friInst.doCartesianImpedanceControl(newCartVals, myCartStiffness, NULL, NULL);
                friInst.doCartesianImpedanceControl(newCartVals, nullCartStiffness, NULL, externalForceVector);/*externalForce.data());*/

//                if (timeCounter < 10)
//                {
//                    friInst.doCartesianImpedanceControl(newCartVals, myCartStiffness, NULL, NULL);
//                    cout << "Waiting.. " << endl;
//                }
//                if (timeCounter < 10 && timeCounter < 15)
//                if (timeCounter < 10)
//                {
//                    friInst.doCartesianImpedanceControl(newCartVals, myCartStiffness, NULL, externalForceVector);
//                    cout << "Input force on z axis : \n" << externalForceVector[2] << endl;
//                }
//                else
//                    friInst.doCartesianImpedanceControl(newCartVals, myCartStiffness, NULL, NULL);
                break;


                default:
                friInst.doDataExchange();
                break;
            }

			// Call to data exchange - and the like 
//			friInst.doPositionControl(newJntVals);

			// have some debug information every n.th. step
//			int divider = (int)( (1./friInst.getSampleTime()) *2.0);
			
//			if ( friInst.getSequenceCount() % divider == 0)
//			{
//				cout << "krl interaction \n";
//				cout << friInst.getMsrBuf().krl;
//				cout << "intf stat interaction \n";
//				cout << friInst.getMsrBuf().intf.stat;
//				cout << "smpl " << friInst.getSampleTime();

//				cout << endl;
//			}

			// Stop request is issued from the other side
			if ( friInst.getFrmKRLInt(0) == -1) 
			{
				cout << "leaving \n";
				break;	  
			}
			//
			// Quality change leads to output of statistics
			// for informational reasons
			//
//			if ( friInst.getQuality() != lastQuality)
//			{
//				cout << "quality change detected "<< friInst.getQuality()<< " \n";
//				cout << friInst.getMsrBuf().intf;
//				cout << endl;
//				lastQuality=friInst.getQuality();
//			}
		}

		/* and leave it on */
    //}

	return EXIT_SUCCESS;
}
/* @} */
