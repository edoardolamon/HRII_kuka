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


#include <vector>
#include <string>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <limits.h>
#include <signal.h>
#include "friudp.h"
#include "friremote.h"
#include "utils.h"
#include "XBotLogger/Logger.hpp"


#ifndef M_PI 
#define M_PI 3.14159
#endif

using namespace std;

volatile bool loop = false;

void handler(int sig)
{
    if (!loop)
        exit(1);
    loop = false;
}

int
#ifdef __WIN32

_tmain
#else
#ifdef _WRS_KERNEL
kuka_cartesianImpedenceController
#else
main
#endif
#endif

(int argc, char *argv[])
{
    signal(SIGINT, handler);
    XBot::MatLogger::Ptr logger = XBot::MatLogger::getLogger("/tmp/kuka_hrii_log");

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
        bool finalPosReached = true;
        double timeCounter = 0.;
        int k = 1;

        float newJntVals[LBR_MNJ],  msrdJntVals[LBR_MNJ];
//         float newCartVals[FRI_CART_FRM_DIM],  msrdCartVals[FRI_CART_FRM_DIM];
        Eigen::Affine3f newCartAffine, msrdCartAffine, initialCartAffine, desCartAffine, unitAffine,
                        nextDesiredAffine, relativeNextDesiredAffine, errorAffine, relativeErrorAffine;
        unitAffine.linear().setIdentity(3,3);
        unitAffine.translation().setZero();
        Eigen::Matrix3f desRotation;
                        desRotation.setIdentity(3,3);
//                         desRotation << 0.7071068,-0.7071068, 0.,
//                                        0.7071068, 0.7071068, 0.,
//                                                0.,       0., 1.;
        Eigen::Vector3f desTranslation;
        desTranslation.setZero(3);
//         desTranslation << 0., 0.2, 0.;
        Eigen::MatrixXf calcJacobian(FRI_CART_VEC,LBR_MNJ), msrdJacobian(FRI_CART_VEC,LBR_MNJ), measToWorldRotation(FRI_CART_VEC, FRI_CART_VEC), 
                        calcJacobianPinv(LBR_MNJ,FRI_CART_VEC), jntIdentity(LBR_MNJ, LBR_MNJ);
                        jntIdentity.setIdentity();
        float * msrdJacobianVector;
        Eigen::VectorXf screwError(FRI_CART_VEC), desScrewError(FRI_CART_VEC), screwErrorPrev(FRI_CART_VEC), relativeScrewError(FRI_CART_VEC),
                        desTwistError(FRI_CART_VEC), twistError(FRI_CART_VEC);
        Eigen::VectorXf msrdCartScrew(FRI_CART_VEC), desCartScrew(FRI_CART_VEC), nextDesCartScrew(FRI_CART_VEC);                
        measToWorldRotation.setZero(FRI_CART_VEC, FRI_CART_VEC);
        screwError.setZero(FRI_CART_VEC);
        twistError.setZero(FRI_CART_VEC);
        Eigen::VectorXf msrdJntPos(LBR_MNJ);

        Eigen::VectorXf relativeDesiredScrew(FRI_CART_VEC), relativeNextDesiredScrew(FRI_CART_VEC);
        vector< vector<float> > cubicTrajCoeffs(FRI_CART_VEC);
        vector< vector<float> > quinticTrajCoeffs(FRI_CART_VEC);
        
        // Cartesian Impedence
        Eigen::MatrixXf myStiffness(FRI_CART_VEC,FRI_CART_VEC), myDamping(FRI_CART_VEC,FRI_CART_VEC);
        Eigen::VectorXf stiffnessVector(FRI_CART_VEC);
                        stiffnessVector << 1500.,1500.,1500.,20.,20.,20.;
//                         stiffnessVector.setZero();
        myStiffness = stiffnessVector.asDiagonal();
        cout << "-------------------Stiffness matrix------------------- \n" << myStiffness << endl;
        float dampingFactor = 0.7;
//         myDamping = 2*dampingFactor*myStiffness.array().sqrt();
//         myDamping = 2*dampingFactor*myStiffness.cwiseSqrt();  // should work in the same way as previous
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> eigenSolver(myStiffness);
        Eigen::MatrixXf eigenvalDiag = eigenSolver.eigenvalues().asDiagonal();
        eigenvalDiag = eigenvalDiag.array().sqrt();
        myDamping = 2 * dampingFactor * eigenSolver.eigenvectors() * eigenvalDiag * eigenSolver.eigenvectors().transpose();
//         Eigen::VectorXf dampingVector(FRI_CART_VEC);
//         dampingVector << 80.,80.,80.,0.,0.,0.;
//         myDamping = dampingVector.asDiagonal();
        cout << "-------------------Damping matrix------------------- \n" << myDamping << endl;
        Eigen::VectorXf desiredForce(FRI_CART_VEC), addForce(FRI_CART_VEC), desAddForce(FRI_CART_VEC);
        desAddForce << 0.,0.,0.,0.,0.,0.;
        Eigen::VectorXf desiredTorque(LBR_MNJ), addTorque(LBR_MNJ), addTorqueFiltered(LBR_MNJ);
        addTorqueFiltered << 0.,0.,0.,0.,0.,0.,0.;
        float commandedTorque[LBR_MNJ] = {0.,0.,0.,0.,0.,0.,0.};
        float smoothingFactor = 0.1;
        float nullJntStiffness[LBR_MNJ]={0,0,0,0,0,0,0};
        
        // Circular trajectory
        const float radius = 0.05;
        const float period = 5.;
        Eigen::Vector3f center;
        
        // Second task: nullspace stiffness
        Eigen::VectorXf secondTaskTorque(LBR_MNJ), secondTaskJntEquilib(LBR_MNJ), secondTaskStiffnessVector(LBR_MNJ), secondTaskDampingVector(LBR_MNJ);
        Eigen::MatrixXf secondTaskStiffness(LBR_MNJ, LBR_MNJ), secondTaskDamping(LBR_MNJ, LBR_MNJ);
        secondTaskStiffnessVector << 20.,20.,20.,20.,20.,20.,20.;
//                         secondTaskStiffnessVector.setZero();
        secondTaskStiffness = secondTaskStiffnessVector.asDiagonal();
        secondTaskDampingVector << 0.5,0.5,0.5,0.5,0.5,0.5,0.5;
//                         secondTaskDampingVector.setZero();
        secondTaskDamping = secondTaskDampingVector.asDiagonal();
        
        // External forces estimated
        Eigen::VectorXf estExtJntTorque(LBR_MNJ), estExtCartForce(FRI_CART_VEC), wrenchError(FRI_CART_VEC), wrenchErrorPrev(FRI_CART_VEC);
        wrenchErrorPrev.setZero();
        Eigen::MatrixXf proportionalWrenchCoeff(FRI_CART_VEC, FRI_CART_VEC), integralWrenchCoeff(FRI_CART_VEC, FRI_CART_VEC);
        proportionalWrenchCoeff.setIdentity();
        integralWrenchCoeff.setIdentity();
        
        // Logger
        logger->log("StiffnessMatrix", myStiffness);
        logger->log("DampingMatrix", myDamping);
        logger->createVectorVariable("desScrewError", FRI_CART_VEC);
        logger->createVectorVariable("screwError", FRI_CART_VEC);
        logger->createVectorVariable("desTwistError", FRI_CART_VEC);
        logger->createVectorVariable("twistError", FRI_CART_VEC);
        logger->createMatrixVariable("nextDesiredAffine", 4, 4);
        logger->createMatrixVariable("msrdCartAffine", 4, 4);
        logger->createVectorVariable("desCartScrew", FRI_CART_VEC);
        logger->createVectorVariable("nextDesCartScrew", FRI_CART_VEC);
        logger->createVectorVariable("msrdCartScrew", FRI_CART_VEC);
        logger->createVectorVariable("desiredForce", FRI_CART_VEC);
        logger->createVectorVariable("estExtCartForce", FRI_CART_VEC);
        logger->createVectorVariable("desiredTorque", LBR_MNJ);
        logger->createVectorVariable("estExtJntTorque", LBR_MNJ);

        cout << endl;
        printf("START CONTROL\n");
        friInst.doDataExchange();
        double deltaT = friInst.getSampleTime();
        logger->add("deltaT", deltaT);
//        printf("Sample time: %f \n", deltaT);
        cout << endl;

//////* enter main loop - wait until we enter stable command mode */
//        for(;;)
        loop = true;
        while(loop)
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

                // Read the current joint position and external joint torques
                for (int i = 0; i < LBR_MNJ; i++)
                {
//                     newJntVals[i] = friInst.getMsrCmdJntPosition()[i];
                    msrdJntVals[i] = friInst.getMsrMsrJntPosition()[i];
                    msrdJntPos(i) = friInst.getMsrMsrJntPosition()[i];
                    estExtJntTorque(i) = friInst.getMsrEstExtJntTrq()[i];
                }

                // Read the current cartesian position
                for (int i = 0; i < FRI_CART_FRM_DIM; i++)
                {
//                     newCartVals[i] = friInst.getMsrCmdCartPosition()[i];
//                     msrdCartVals[i] = friInst.getMsrCartPosition()[i];

                    int row = i / 4;
                    int col = i % 4;
                    newCartAffine(row,col) = friInst.getMsrCmdCartPosition()[i];
                    msrdCartAffine(row,col) = friInst.getMsrCartPosition()[i];
                }
                
                // Read the external cartesian forces
                for (int i = 0; i < FRI_CART_VEC; i++)
                {
                    estExtCartForce(i) = friInst.getMsrEstExtTcpFT()[i];
                }

                /// First loop computation
                if (firstLoopPlan)
                {
                    // Desired pose initialization
                    initialCartAffine = msrdCartAffine;
                    cout << "Initial pose transformation: \n" << initialCartAffine.matrix() << endl;
                    cout << endl;
                    logger->add("nextDesiredAffine", initialCartAffine.matrix());
                    desCartAffine.linear() =  desRotation;
                    desCartAffine.translation() = desTranslation;
                    desCartAffine = initialCartAffine * desCartAffine;
// 		    
//                     cout << "Desired pose transformation: \n" << desCartAffine.matrix() << endl;
//                     cout << endl;
                    secondTaskJntEquilib = msrdJntPos;
//                     for (int i = 0; i < LBR_MNJ; i++)
//                     {
//                         secondTaskJntEquilib(i) = msrdJntVals[i];
//                     }
                     
                   // Planning a point to point trajectory
                   relativeDesiredScrew = utils::getScrewErrorFromAffine(initialCartAffine, desCartAffine);
//                    cout << "Screw computed: \n" << relativeDesiredScrew.transpose() << endl;

                   for (int i = 0; i < FRI_CART_VEC; i++)
                   {
//                        cubicTrajCoeffs[i] = utils::computeCubicTrajectory(relativeDesiredScrew(i));
                       quinticTrajCoeffs[i] = utils::computeQuinticTrajectory(relativeDesiredScrew(i));

//                        cout << "Trajectory coeffs " << to_string(i) << " : a0 = " << cubicTrajCoeffs[i][0]
//                             << ", a1 = " << cubicTrajCoeffs[i][1] << ", a2 = " << cubicTrajCoeffs[i][2] << ", a3 = " << cubicTrajCoeffs[i][3] << endl;
//                        cout << "Trajectory coeffs " << to_string(i) << " : a0 = " << quinticTrajCoeffs[i][0]
//                             << ", a1 = " << quinticTrajCoeffs[i][1] << ", a2 = " << quinticTrajCoeffs[i][2] << ", a3 = " << quinticTrajCoeffs[i][3]
//                             << ", a4 = " << quinticTrajCoeffs[i][4] << ", a5 = " << quinticTrajCoeffs[i][5] << endl;
                   }
                   
                   // Planning a circular trajectory (fix the center w.r.t. world frame)
                   center(0) = initialCartAffine.translation()(0) + radius;
                   center(1) = initialCartAffine.translation()(1);
                   center(2) = initialCartAffine.translation()(2);

                   firstLoopPlan = false;
                }

                if ( friInst.getState() == FRI_STATE_CMD)
                {
                    if ( friInst.isPowerOn() )
                    {
                        // Update time counter
                        timeCounter += deltaT;
                        logger->add("msrdCartAffine", msrdCartAffine.matrix());
                        msrdCartScrew = utils::getScrewFromAffine(msrdCartAffine);
                        logger->add("msrdCartScrew", msrdCartScrew);
                           
//                         desCartAffine = utils::nextStepCircular(initialCartAffine,timeCounter,period,radius,center);
//                         logger->add("desCartAffine", desCartAffine.matrix());
//                         relativeDesiredScrew = utils::getScrewErrorFromAffine(initialCartAffine, desCartAffine);
//                         for (int i = 0; i < FRI_CART_VEC; i++)
//                         {
//                              quinticTrajCoeffs[i] = utils::computeQuinticTrajectory(relativeDesiredScrew(i), 1);
//                         }
                        
                        // Compute next desired screw error
                        if (timeCounter <= 5.)
                        {
                             for (int i = 0; i < FRI_CART_VEC; i++)
                             {
                                 relativeNextDesiredScrew(i) = utils::nextStepQuintic(quinticTrajCoeffs[i], timeCounter);
                             }
                        }
                        else if (finalPosReached)
                        {
                             cout << "Desired position reached" << endl;
                             finalPosReached = false;
                        }

                       relativeNextDesiredAffine = utils::getAffineFromScrew(relativeNextDesiredScrew);
//                        cout << "Relative next desired screw converted into an affine matrix: \n" << relativeNextDesiredAffine.matrix() << endl;
                       nextDesiredAffine = initialCartAffine * relativeNextDesiredAffine;
                       // Circular trajectory
//                        nextDesiredAffine = utils::nextStepCircular(initialCartAffine,timeCounter,period,radius,center);
                       logger->add("nextDesiredAffine", nextDesiredAffine.matrix());
//                        cout << "Next desired affine transformed wrt world frame: \n" << nextDesiredAffine.matrix() << endl;
                       desCartScrew = utils::getScrewFromAffine(desCartAffine);
                       nextDesCartScrew = utils::getScrewFromAffine(nextDesiredAffine);
                       logger->add("desCartScrew", desCartScrew);
                       logger->add("nextDesCartScrew", nextDesCartScrew);

                        // Compute position and orientation error
                       relativeScrewError = utils::getScrewErrorFromAffine(msrdCartAffine, nextDesiredAffine);
//                         relativeScrewError = utils::getScrewErrorFromAffine(msrdCartAffine, initialCartAffine);
//                        logger->add("relativeScrewError", relativeScrewError);
                       measToWorldRotation.topLeftCorner(3,3) = msrdCartAffine.linear();
                       measToWorldRotation.bottomRightCorner(3,3) = msrdCartAffine.linear();
//                        cout << "Measured to world rotation: \n" << measToWorldRotation << endl;

                       desScrewError = measToWorldRotation * relativeScrewError;
                       logger->add("desScrewError", desScrewError);
                       screwError = smoothingFactor * desScrewError + (1 - smoothingFactor)*screwError;
                       logger->add("screwError", screwError);
                       if (firstLoopCmd)                   
                       {
                           screwErrorPrev = desScrewError;
                       }
                       desTwistError = (desScrewError - screwErrorPrev) / deltaT;
//                        cout << "Velocity error:" << desTwistError.transpose() << endl;
                       logger->add("desTwistError", desTwistError);
                       twistError = smoothingFactor * desTwistError + (1 - smoothingFactor)*twistError;
                       logger->add("twistError", twistError);
                       screwErrorPrev = desScrewError;
//                        cout << "Velocity error filtered:" << twistError.transpose() << endl;

                       // Compute desired cartesian force through impedance behavior
                       desiredForce = myStiffness * desScrewError + myDamping * twistError;
                       addForce = smoothingFactor * desAddForce + (1 - smoothingFactor) * addForce;
                       wrenchError = addForce - estExtCartForce; 
                       logger->add("desiredForce", desiredForce);
                       logger->add("estExtCartForce", estExtCartForce);

                       // Compute jacobian
                       calcJacobian = utils::calcKukaJacob(msrdJntVals);
//                        cout << "Jacobian computed: \n" << calcJacobian << endl;
//                        msrdJacobianVector = friInst.getJacobian();
//                        for (int i = 0; i < LBR_MNJ*FRI_CART_VEC; i++)
//                        {
//                            int col = i % LBR_MNJ;
//                            int row = i / LBR_MNJ;
//                            msrdJacobian(row,col) = msrdJacobianVector[i];
//                        }

                       // Compute desired joint torque
                       desiredTorque = calcJacobian.transpose() * desiredForce; //(desiredForce + addForce);
                       addTorque = calcJacobian.transpose() * 0.5 * (proportionalWrenchCoeff * wrenchError + integralWrenchCoeff * (wrenchError - wrenchErrorPrev) * deltaT);
                       calcJacobianPinv = calcJacobian.transpose() * (calcJacobian * calcJacobian.transpose()).inverse();
                       wrenchErrorPrev = wrenchError;
                       secondTaskTorque = (jntIdentity - calcJacobianPinv * calcJacobian) * secondTaskStiffness * (secondTaskJntEquilib - msrdJntPos);
//                        cout << "Desired torque computed" << endl;
//                        desiredTorque = msrdJacobian.transpose() * desiredForce;
//                        cout << "Input torque : \n" << desiredTorque.transpose() << endl;
//                        logger->add("desiredTorque", desiredTorque);
                       for (int i = 0; i < LBR_MNJ; i++)
                       {
//                            addTorqueFiltered(i) = smoothingFactor * addTorque(i) + (1 - smoothingFactor)*addTorqueFiltered(i);
                           addTorqueFiltered(i) = addTorque(i);
                           commandedTorque[i] = desiredTorque(i) + addTorqueFiltered(i);
//                            commandedTorque[i] = desiredTorque(i) + addTorqueFiltered(i) + secondTaskTorque(i);
                       }
                       logger->add("desiredTorque", desiredTorque);  
                       logger->add("estExtJntTorque", estExtJntTorque);
                       if (firstLoopCmd)
                       {
                           firstLoopCmd = false;
                       }
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

//                 friInst.doJntImpedanceControl(NULL, nullJntStiffness, NULL, NULL);
//                 friInst.doJntImpedanceControl(NULL, NULL, NULL, NULL);
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

    logger->flush();

	return EXIT_SUCCESS;
}
