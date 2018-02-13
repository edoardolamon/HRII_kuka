/** 
 \author (Edoardo Lamon)
    \file kuka_cartesianImpedenceController.cpp
    \brief Planning and control of a kuka LW4 through a Cartesian Impedence Controller
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

Eigen::MatrixXf _calcKukaJacob(float jointPos[]);
Eigen::VectorXf _getScrewFromAffine(Eigen::Affine3f &sourceAffine, Eigen::Affine3f &targetAffine);
Eigen::Affine3f _getAffineFromScrew(Eigen::VectorXf screw);
vector<float> _computeCubicTrajectory(float &relativePos, float avgJntVel=0.174);
vector<float> _computeQuinticTrajectory(float &relativePos);
float _nextStepCubic(vector<float> &cubicCoeffs, float t);
float _nextStepQuintic(vector<float> &quinticCoeffs, float t);


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
        cout << endl;

        bool firstLoopPlan = true;
        bool firstLoopCmd = true;
        double timeCounter = 0.;
        int k = 1;

        float newJntVals[LBR_MNJ],  msrdJntVals[LBR_MNJ];
        float newCartVals[FRI_CART_FRM_DIM],  msrdCartVals[FRI_CART_FRM_DIM];
        Eigen::Affine3f newCartAffine, msrdCartAffine, initialCartAffine, desCartAffine, unitAffine,
                        nextDesiredAffine, relativeNextDesiredAffine, errorAffine, relativeErrorAffine;
        unitAffine.linear().setIdentity(3,3);
        unitAffine.translation().setZero();
        Eigen::Matrix3f desRotation;
                        desRotation.setIdentity(3,3);
//                        desRotation << 1., 0., 0.,
//                                       0.,-1., 0.,
//                                       0., 0.,-1.;
        Eigen::Vector3f desTranslation;
                        desTranslation << 0., 0., 0.;
//        Eigen::Matrix3f msrdCartRotation;
//        Eigen::Vector3f msrdCartTranslation;
        Eigen::MatrixXf calcJacobian(FRI_CART_VEC,LBR_MNJ), msrdJacobian(FRI_CART_VEC,LBR_MNJ), measToWorldRotation(FRI_CART_VEC, FRI_CART_VEC);
        float * msrdJacobianVector;
        Eigen::VectorXf screwError(FRI_CART_VEC), screwErrorPrev(FRI_CART_VEC), relativeScrewError(FRI_CART_VEC),
                        desTwistError(FRI_CART_VEC), twistError(FRI_CART_VEC);
        measToWorldRotation.setZero(FRI_CART_VEC, FRI_CART_VEC);
        twistError.setZero(FRI_CART_VEC);
//        Eigen::MatrixXf measToWorldRotation = Eigen::MatrixXf::Zero(FRI_CART_VEC, FRI_CART_VEC);


        Eigen::VectorXf relativeDesiredScrew(FRI_CART_VEC), relativeNextDesiredScrew(FRI_CART_VEC);
        vector< vector<float> > cubicTrajCoeffs(FRI_CART_VEC);
        vector< vector<float> > quinticTrajCoeffs(FRI_CART_VEC);

        Eigen::MatrixXf myStiffness(FRI_CART_VEC,FRI_CART_VEC), myDamping(FRI_CART_VEC,FRI_CART_VEC);
        Eigen::VectorXf stiffnessVector(FRI_CART_VEC);
                        stiffnessVector << 1000.,1000.,1000.,100.,100.,100.;
        myStiffness = stiffnessVector.asDiagonal();
        cout << "-------------------Stiffness matrix------------------- \n" << myStiffness << endl;
        float dampingFactor = 0.7;
//        myDamping = 2*dampingFactor*myStiffness.array().sqrt();
//        myDamping = 2*dampingFactor*myStiffness.cwiseSqrt();  // should work in the same way as previous
//        myDamping.setIdentity(FRI_CART_VEC,FRI_CART_VEC);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> eigenSolver(myStiffness);
        Eigen::MatrixXf eigenvalDiag = eigenSolver.eigenvalues().asDiagonal();
        eigenvalDiag = eigenvalDiag.array().sqrt();
//        cout << "Eigenvalues: " << eigenSolver.eigenvalues() << endl;
//        cout << "Eigenvectors: \n" << eigenSolver.eigenvectors() << endl;
//        cout << "Eigenval : \n" << eigenvalDiag << endl;
        myDamping = 2 * dampingFactor * eigenSolver.eigenvectors() * eigenvalDiag * eigenSolver.eigenvectors().transpose();
        cout << "-------------------Damping matrix------------------- \n" << myDamping << endl;

        Eigen::VectorXf desiredForce(FRI_CART_VEC);
//                        desiredForce << 10.,0.,0.,0.,0.,0.;
        Eigen::VectorXf desiredTorque(LBR_MNJ);
        float commandedTorque[LBR_MNJ] = {0.,0.,0.,0.,0.,0.,0.};
        float smoothingFactor = 0.1;
        float nullJntStiffness[LBR_MNJ]={0,0,0,0,0,0,0};

        cout << endl;
        printf("START CONTROL\n");
        friInst.doDataExchange();
        double deltaT = friInst.getSampleTime();
        printf("Sample time: %f \n", deltaT);
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

            // just mirror the real value
			friInst.setToKRLReal(0,friInst.getFrmKRLReal(1));

            if ( lastCtrlScheme != friInst.getCurrentControlScheme())
            {
                cout << "Switching control scheme " << lastCtrlScheme;
                lastCtrlScheme = friInst.getCurrentControlScheme();
                cout << " to " << lastCtrlScheme << endl;
            }

            switch (friInst.getCurrentControlScheme())
            {
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
                    msrdCartAffine(row,col) = friInst.getMsrCartPosition()[i];
//                    initialCartAffine(row,col) = friInst.getMsrCartPosition()[i];
                }


//                msrdCartRotation = msrdCartAffine.linear();
//                msrdCartTranslation = msrdCartAffine.translation();

                /// First loop computation
                if (firstLoopPlan)
                {
                    // Desired pose initialization
                    initialCartAffine = msrdCartAffine;
                    cout << "Initial pose transformation: \n" << initialCartAffine.matrix() << endl;
                    cout << endl;
                    desCartAffine.linear() = initialCartAffine.linear() * desRotation;
                    desCartAffine.translation() = desTranslation + initialCartAffine.translation();
                    cout << "Desired pose transformation: \n" << desCartAffine.matrix() << endl;
                    cout << endl;
//                    calcJacobian = _calcKukaJacob(msrdJntVals);
//                    cout << "Jacobian = " << calcJacobian << endl;
//                    msrdJacobianVector = friInst.getJacobian();
//                    for (int i = 0; i < 42; i++)
//                    {
//                        int col = i % 7;
//                        int row = i / 7;
//                        msrdJacobian(row,col) = msrdJacobianVector[i];
//                    }
//                    cout << "Measured jacobian = " << msrdJacobian << endl;

                    // Planning
                    relativeDesiredScrew = _getScrewFromAffine(initialCartAffine, desCartAffine);
                    cout << "Screw computed: \n" << relativeDesiredScrew.transpose() << endl;

                    for (int i = 0; i < FRI_CART_VEC; i++)
                    {
                        cubicTrajCoeffs[i] = _computeCubicTrajectory(relativeDesiredScrew(i));
                        quinticTrajCoeffs[i] = _computeQuinticTrajectory(relativeDesiredScrew(i));

//                        cout << "Trajectory coeffs " << to_string(i) << " : a0 = " << cubicTrajCoeffs[i][0]
//                             << ", a1 = " << cubicTrajCoeffs[i][1] << ", a2 = " << cubicTrajCoeffs[i][2] << ", a3 = " << cubicTrajCoeffs[i][3] << endl;
//                        cout << "Trajectory coeffs " << to_string(i) << " : a0 = " << quinticTrajCoeffs[i][0]
//                             << ", a1 = " << quinticTrajCoeffs[i][1] << ", a2 = " << quinticTrajCoeffs[i][2] << ", a3 = " << quinticTrajCoeffs[i][3]
//                             << ", a4 = " << quinticTrajCoeffs[i][4] << ", a5 = " << quinticTrajCoeffs[i][5] << endl;
                    }

                    firstLoopPlan = false;
                }

                if ( friInst.getState() == FRI_STATE_CMD)
                {
                    if ( friInst.isPowerOn() )
                    {
                        // Update time counter
                        timeCounter += deltaT;

                        if (firstLoopCmd)
                        {
                            initialCartAffine = msrdCartAffine;
                        }

                        if (timeCounter == k*100*deltaT)
                        {
//                            cout << "---------------Measured cartesian pose--------------- \n" << msrdCartAffine.matrix() << endl;
                        }

                        // Compute next desired pose
                        for (int i = 0; i < FRI_CART_VEC; i++)
                        {
//                           relativeNextDesiredScrew(i) = _nextStepCubic(cubicTrajCoeffs[i], timeCounter);
                           relativeNextDesiredScrew(i) = _nextStepQuintic(quinticTrajCoeffs[i], timeCounter);
                        }
//                        if (timeCounter == k*100*deltaT)
//                        {
//                            cout << "Relative next desired screw planned: \n" << relativeNextDesiredScrew << endl;
//                        }
                        relativeNextDesiredAffine = _getAffineFromScrew(relativeNextDesiredScrew);
//                        cout << "Relative next desired screw converted into an affine matrix: \n" << relativeNextDesiredAffine.matrix() << endl;
                        nextDesiredAffine = initialCartAffine * relativeNextDesiredAffine;
//                        cout << "Next desired affine transformed wrt world frame: \n" << nextDesiredAffine.matrix() << endl;

                        // Compute position and orientation error
//                        relativeScrewError = _getScrewFromAffine(msrdCartAffine, nextDesiredAffine);
                        relativeScrewError = _getScrewFromAffine(msrdCartAffine, initialCartAffine);
                        measToWorldRotation.topLeftCorner(3,3) = msrdCartAffine.linear();
                        measToWorldRotation.bottomRightCorner(3,3) = msrdCartAffine.linear();
//                        if (timeCounter == k*100*deltaT)
//                        {
                            cout << "error:" << relativeScrewError.transpose() << endl;
//                        }
//                        relativeErrorAffine = _getAffineFromScrew(relativeScrewError);
//                        cout << "Relative next screw error converted into an affine matrix: \n" << relativeErrorAffine.matrix() << endl;
//                        errorAffine = msrdCartAffine * relativeErrorAffine;
//                        cout << "Next affine error transformed wrt world frame: \n" << errorAffine.matrix() << endl;
//                        screwError = _getScrewFromAffine(unitAffine, errorAffine);
                        screwError = measToWorldRotation * relativeScrewError;
//                        screwError(3) = 0.;
//                        screwError(4) = 0.;
//                        screwError(5) = 0.;
                        if (timeCounter == k*100*deltaT)
                        {
//                            cout << "Next screw error computed: \n" << screwError.transpose() << endl;
                        }
                        if (firstLoopCmd)                   
                        {
                            screwErrorPrev = screwError;
                        }
                        desTwistError = (screwError - screwErrorPrev) / deltaT;
                        cout << "v error:" << desTwistError.transpose() << endl;
                        twistError = smoothingFactor * desTwistError + (1 - smoothingFactor)*twistError;
                        screwErrorPrev = screwError;
                        cout << "v error filt:" << twistError.transpose() << endl;

//                        cout << "Control errors computed" << endl;

                        // Compute jacobian
                        calcJacobian = _calcKukaJacob(msrdJntVals);
//                        cout << "Jacobian computed: \n" << calcJacobian << endl;
                        msrdJacobianVector = friInst.getJacobian();
                        for (int i = 0; i < LBR_MNJ*FRI_CART_VEC; i++)
                        {
                            int col = i % LBR_MNJ;
                            int row = i / LBR_MNJ;
                            msrdJacobian(row,col) = msrdJacobianVector[i];
                        }
//                        cout << "KUKA Jacobian computed" << endl;

                        // Compute desired cartesian force through impedance behavior
                        desiredForce = myStiffness * screwError + myDamping * twistError;
                        if (timeCounter == k*100*deltaT)
                        {
//                           cout << "Desired force computed: \n" << desiredForce << endl;
                           k++;
                        }

                        // Compute desired joint torque
                        desiredTorque = calcJacobian.transpose() * desiredForce;
//                        desiredTorque.setZero();
//                        cout << "Desired torque computed" << endl;
//                        desiredTorque = msrdJacobian.transpose() * desiredForce;
//                        cout << "Input torque : \n" << desiredTorque.transpose() << endl;
                        for (int i = 0; i < LBR_MNJ; i++)
                        {
//                            commandedTorque[i] = smoothingFactor * desiredTorque(i) + (1 - smoothingFactor)*commandedTorque[i];
                            commandedTorque[i] = desiredTorque(i);
                        }
//                        cout << "Commanded torque computed" << endl;
//                        cout << endl;
                        if (firstLoopCmd)
                        {
                            firstLoopCmd = false;
                        }
//                        return 0;
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

//                friInst.doJntImpedanceControl(NULL, nullJntStiffness, NULL, NULL);
//                friInst.doJntImpedanceControl(NULL, NULL, NULL, NULL);
                friInst.doJntImpedanceControl(NULL, nullJntStiffness, NULL, commandedTorque);
                break;

                case FRI_CTRL_CART_IMP:

                break;

                default:
                friInst.doDataExchange();
                break;
            }


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
            if ( friInst.getQuality() != lastQuality)
            {
                cout << "quality change detected "<< friInst.getQuality()<< " \n";
                cout << friInst.getMsrBuf().intf;
                cout << endl;
                lastQuality=friInst.getQuality();
            }
		}

		/* and leave it on */
    //}

	return EXIT_SUCCESS;
}

// Methods implementation
Eigen::MatrixXf _calcKukaJacob(float jointPos[])
{
//    if (!(J_kuka.rows() == FRI_CART_VEC && J_kuka.cols() == LBR_MNJ))  // (6x7)
//    {
//        cerr << "Jacobian wrong dimensions. Jacobian should be 6x7" << endl;
//        J_kuka(FRI_CART_VEC,LBR_MNJ);
//    }

    float C1,S1,S2,C2,S3,C3,S4,C4,S5,C5,S6,C6,S7,C7;
    float l1=0.310, l2=0.4, l3=0.39, l4=0.075;
    Eigen::MatrixXf J_kuka(FRI_CART_VEC,LBR_MNJ);

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

//    J_kuka.row(0) << - l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3))) - l2*C2*S1, -C1*(l3*(C4*S2 - C2*C3*S4) + l2*S2 + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))) , C2*S1*(l3*(C4*S2 - C2*C3*S4) + l2*S2 + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))) - S2*(l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3))) + l2*C2*S1) , (C1*C3 - S1*S2*S3)*(l3*(C4*S2 - C2*C3*S4) + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))) - C2*S3*(l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3)))) , (S4*(C1*S3 + C3*S1*S2) + C2*C4*S1)*(l3*(C4*S2 - C2*C3*S4) + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))) - (l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3))))*(C4*S2 - C2*C3*S4) , l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))*(S5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) - C5*(C1*C3 - S1*S2*S3)) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3)))*(S5*(S2*S4 + C2*C3*C4) + C2*C5*S3) , 0;
//    J_kuka.row(1) << l2*C1*C2 - l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3))) - l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4), -S1*(l3*(C4*S2 - C2*C3*S4) + l2*S2 + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))), - S2*(l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3))) - l2*C1*C2) - C1*C2*(l3*(C4*S2 - C2*C3*S4) + l2*S2 + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))), (C3*S1 + C1*S2*S3)*(l3*(C4*S2 - C2*C3*S4) + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))) - C2*S3*(l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3)))), (l3*(C4*S2 - C2*C3*S4) + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4)))*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - (C4*S2 - C2*C3*S4)*(l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3)))), l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))*(S5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) - C5*(C3*S1 + C1*S2*S3)) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3)))*(S5*(S2*S4 + C2*C3*C4) + C2*C5*S3), 0;
//    J_kuka.row(2) << 0 , S1*(l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3))) + l2*C2*S1) - C1*(l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3))) - l2*C1*C2), C2*S1*(l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3))) - l2*C1*C2) + C1*C2*(l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3))) + l2*C2*S1), (C1*C3 - S1*S2*S3)*(l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3)))) - (l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3))))*(C3*S1 + C1*S2*S3), (S4*(C1*S3 + C3*S1*S2) + C2*C4*S1)*(l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3)))) - (l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3))))*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4), l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3)))*(S5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) - C5*(C1*C3 - S1*S2*S3)) - l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3)))*(S5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) - C5*(C3*S1 + C1*S2*S3)), 0;
//    J_kuka.row(3) << 0 , S1 , C1*C2 , - C3*S1 - C1*S2*S3, C1*C2*C4 - S4*(S1*S3 - C1*C3*S2), C5*(C3*S1 + C1*S2*S3) - S5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4), S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3)) - C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4);
//    J_kuka.row(4) << 0 , -C1 , C2*S1 , C1*C3 - S1*S2*S3, S4*(C1*S3 + C3*S1*S2) + C2*C4*S1, S5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) - C5*(C1*C3 - S1*S2*S3), C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3));
//    J_kuka.row(5) << 1 , 0 , S2 , C2*S3 , C4*S2-C2*C3*S4 , - S5*(S2*S4+C2*C3*C4) - C2*C5*S3, S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4);
    J_kuka <<- l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3))) - l2*C2*S1, -C1*(l3*(C4*S2 - C2*C3*S4) + l2*S2 + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))), C2*S1*(l3*(C4*S2 - C2*C3*S4) + l2*S2 + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))) - S2*(l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3))) + l2*C2*S1), (C1*C3 - S1*S2*S3)*(l3*(C4*S2 - C2*C3*S4) + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))) - C2*S3*(l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3)))), (S4*(C1*S3 + C3*S1*S2) + C2*C4*S1)*(l3*(C4*S2 - C2*C3*S4) + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))) - (l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3))))*(C4*S2 - C2*C3*S4), l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))*(S5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) - C5*(C1*C3 - S1*S2*S3)) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3)))*(S5*(S2*S4 + C2*C3*C4) + C2*C5*S3), 0,
           l2*C1*C2 - l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3))) - l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4), -S1*(l3*(C4*S2 - C2*C3*S4) + l2*S2 + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))), - S2*(l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3))) - l2*C1*C2) - C1*C2*(l3*(C4*S2 - C2*C3*S4) + l2*S2 + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))), (C3*S1 + C1*S2*S3)*(l3*(C4*S2 - C2*C3*S4) + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))) - C2*S3*(l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3)))), (l3*(C4*S2 - C2*C3*S4) + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4)))*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - (C4*S2 - C2*C3*S4)*(l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3)))), l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))*(S5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) - C5*(C3*S1 + C1*S2*S3)) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3)))*(S5*(S2*S4 + C2*C3*C4) + C2*C5*S3), 0,
           0, S1*(l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3))) + l2*C2*S1) - C1*(l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3))) - l2*C1*C2), C2*S1*(l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3))) - l2*C1*C2) + C1*C2*(l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3))) + l2*C2*S1), (C1*C3 - S1*S2*S3)*(l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3)))) - (l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3))))*(C3*S1 + C1*S2*S3), (S4*(C1*S3 + C3*S1*S2) + C2*C4*S1)*(l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3)))) - (l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3))))*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4), l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3)))*(S5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) - C5*(C1*C3 - S1*S2*S3)) - l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3)))*(S5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) - C5*(C3*S1 + C1*S2*S3)), 0,
           0, S1, C1*C2, - C3*S1 - C1*S2*S3, C1*C2*C4 - S4*(S1*S3 - C1*C3*S2), C5*(C3*S1 + C1*S2*S3) - S5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4), S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3)) - C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4),
           0, -C1, C2*S1, C1*C3 - S1*S2*S3, S4*(C1*S3 + C3*S1*S2) + C2*C4*S1, S5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) - C5*(C1*C3 - S1*S2*S3), C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3)),
           1, 0, S2, C2*S3, C4*S2 - C2*C3*S4, - S5*(S2*S4 + C2*C3*C4) - C2*C5*S3, S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4);

    return J_kuka;

}


Eigen::VectorXf _getScrewFromAffine(Eigen::Affine3f &sourceAffine, Eigen::Affine3f &targetAffine)
{
    Eigen::Affine3f targetFromSource = sourceAffine.inverse() * targetAffine;
    Eigen::AngleAxisf angleAxis(targetFromSource.linear());
    Eigen::VectorXf screw(FRI_CART_VEC);
//    screw << targetFromSource.translation();
//    screw << angleAxis.angle()*angleAxis.axis();
    Eigen::Vector3f w = angleAxis.angle()*angleAxis.axis();
//    screw(0) = targetFromSource.translation()(0);
//    screw(1) = targetFromSource.translation()(1);
//    screw(2) = targetFromSource.translation()(2);
//    screw(3) = w(0);
//    screw(4) = w(1);
//    screw(5) = w(2);
//    screw.head(3) = targetFromSource.translation();
//    screw.segment(3,3) = w;
      screw << targetFromSource.translation(), w;
      for(int i=0;i<3;i++)
      {
          if ( screw(3+i) > 0.03 )
              screw(3+i) = screw(3+i) - 0.03;
          else if ( screw(3+i) < -0.03 )
              screw(3+i) = screw(3+i) + 0.03;
          else
              screw(3+i) = 0;
      }

    return screw;
}

Eigen::Affine3f _getAffineFromScrew(Eigen::VectorXf screw)
{
//    cout << "-----------------------START _getAffineFromScrew-----------------------" << endl;
    float angle;
    Eigen::Vector3f axisAngle, axis;
    Eigen::Matrix3f rotation;
    Eigen::Affine3f screwAffine;
    Eigen::VectorXf myScrew(FRI_CART_VEC);

    myScrew = screw;
//    cout << "Input screw : " << screw.transpose() << endl;
    axisAngle(0) = myScrew(3);
    axisAngle(1) = myScrew(4);
    axisAngle(2) = myScrew(5);
//    axisAngle = myScrew.segment(3,3);
//    cout << "axisAngle: \n" << axisAngle.transpose() << endl;
    angle = axisAngle.norm();
    axis = axisAngle.normalized();
//    cout << "Angle is " << angle << " and axis is \n" << axis << endl;
    rotation = Eigen::AngleAxisf(angle,axis);
    screwAffine.linear() = rotation;
    screwAffine.translation() = screw.head(3);

//    cout << "-----------------------END _getAffineFromScrew-----------------------" << endl;
    return screwAffine;
}

vector<float> _computeCubicTrajectory(float &relativePos, float avgJntVel) /*float avgJntVel=0.174)*/
{
    //float finalTime = abs(relativePos / avgJntVel);
    float finalTime = 5.;
//    cout << "Planning time: " << finalTime << " seconds" << endl;

    vector<float> trajCoeffs(4);
    trajCoeffs[0] = 0; // initial position
    trajCoeffs[1] = 0; // initial velocity
    trajCoeffs[2] = (3*(relativePos) / (pow(finalTime, 2)));   /*(-3*(initialJntPos - finalJntPos) - (2*initialJntVel + finalJntVel)) / (pow(finalTime, 2));*/
    trajCoeffs[3] = (-2*(relativePos)) / (pow(finalTime, 3));  /*(2*(initialJntPos - finalJntPos) + (initialJntVel + finalJntVel)) / (pow(finalTime, 3));*/

    return trajCoeffs;
}


vector<float> _computeQuinticTrajectory(float &relativePos)
{
    float finalTime = 5.;
//    cout << "Planning time: " << finalTime << " seconds" << endl;

    vector<float> trajCoeffs(6);
    trajCoeffs[0] = 0; // initial position
    trajCoeffs[1] = 0; // initial velocity
    trajCoeffs[2] = 0; // intial acceleration/2
    trajCoeffs[3] = (20*(relativePos)) / 2*(pow(finalTime, 3));
    trajCoeffs[4] = (-30*(relativePos)) / 2*(pow(finalTime, 4));
    trajCoeffs[5] = (20*(relativePos)) / 2*(pow(finalTime, 5));

    return trajCoeffs;
}

float _nextStepCubic(vector<float> &cubicCoeffs, float t)
{
    return cubicCoeffs[3]*pow(t,3.) + cubicCoeffs[2]*pow(t,2.) + cubicCoeffs[1]*t + cubicCoeffs[0]; // static_cast<float>
}

float _nextStepQuintic(vector<float> &quinticCoeffs, float t)
{
    return quinticCoeffs[5]*pow(t,5.) + quinticCoeffs[4]*pow(t,4.) + quinticCoeffs[3]*pow(t,3.) + quinticCoeffs[2]*pow(t,2.) + quinticCoeffs[1]*t + quinticCoeffs[0];
}
