/**
 \author Edoardo Lamon
        \file KukaLWRInterface_test.cpp
    \brief Test for KukaLWRInterface library
*/

#include "../include/KukaLWRInterface.h"

using namespace std;

int main(int argc, char** argv)
{   
    //Initialization
    KukaLWRInterface kukaInterface = KukaLWRInterface(49938,"192.168.0.20");
    
    double controlTime = 0.;
    double planTime = 5;
    
    //Planning variables
    bool planning = true;
    bool firstLoop = true;
    Eigen::Affine3f initialCartAffine, desCartAffine, msrdCartAffine, nextDesiredAffine, relativeNextDesiredAffine, errorAffine, relativeErrorAffine;
    vector< vector<float> > quinticTrajCoeffs(FRI_CART_VEC);
    Eigen::Matrix3f desRotation;
//     desRotation.setIdentity(3,3);
    desRotation << 0.7071068,-0.7071068, 0.,
                   0.7071068, 0.7071068, 0.,
                          0.,        0., 1.;
    Eigen::Vector3f desTranslation;
//     desTranslation.setZero(3);
    desTranslation << 0., 0.2, 0.;
    Eigen::VectorXf relativeDesiredScrew(FRI_CART_VEC), relativeNextDesiredScrew(FRI_CART_VEC);
    Eigen::VectorXf screwError(FRI_CART_VEC), screwErrorPrev(FRI_CART_VEC), relativeScrewError(FRI_CART_VEC),
                    desTwistError(FRI_CART_VEC), twistError(FRI_CART_VEC);
    Eigen::MatrixXf jacobian(FRI_CART_VEC,LBR_MNJ), measToWorldRotation(FRI_CART_VEC, FRI_CART_VEC);
    measToWorldRotation.setZero(FRI_CART_VEC, FRI_CART_VEC);
    twistError.setZero(FRI_CART_VEC);
    
    Eigen::VectorXf desiredForce(FRI_CART_VEC);
    Eigen::VectorXf desiredTorque(LBR_MNJ);
                    desiredTorque.setZero();
    float smoothingFactor = 0.1;
    
    
    Eigen::MatrixXf myStiffness(FRI_CART_VEC,FRI_CART_VEC), myDamping(FRI_CART_VEC,FRI_CART_VEC);
    Eigen::VectorXf stiffnessVector(FRI_CART_VEC);
                    stiffnessVector << 1500.,1500.,1500.,0.,0.,0.;
    myStiffness = stiffnessVector.asDiagonal();
    cout << "-------------------Stiffness matrix------------------- \n" << myStiffness << endl;
    float dampingFactor = 1.;
//     myDamping = 2*dampingFactor*myStiffness.array().sqrt();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> eigenSolver(myStiffness);
    Eigen::MatrixXf eigenvalDiag = eigenSolver.eigenvalues().asDiagonal();
    eigenvalDiag = eigenvalDiag.array().sqrt();
    myDamping = 2 * dampingFactor * eigenSolver.eigenvectors() * eigenvalDiag * eigenSolver.eigenvectors().transpose();
//     Eigen::VectorXf dampingVector(FRI_CART_VEC);
//     dampingVector << 80.,80.,80.,0.,0.,0.;
//     myDamping = dampingVector.asDiagonal();
    cout << "-------------------Damping matrix------------------- \n" << myDamping << endl;

    double sampleTime = kukaInterface.getSampleTime();
    for(;;)
    {
        kukaInterface.communicationQualityCheck();
        
        switch (kukaInterface.getCurrentControlScheme())
        {
            case FRI_CTRL_JNT_IMP:
                
                 if (planning)
                 {
                    //  Planning a point to point trajectory
                    initialCartAffine = kukaInterface.getMsrdCartPosition();
                    cout << "Initial pose transformation: \n" << initialCartAffine.matrix() << endl;
                    cout << endl;
                    desCartAffine.linear() =  desRotation;
                    desCartAffine.translation() = desTranslation;
                    desCartAffine = initialCartAffine * desCartAffine;
                    cout << "Desired pose transformation: \n" << desCartAffine.matrix() << endl;
                    cout << endl;
                    relativeDesiredScrew = utils::getScrewErrorFromAffine(initialCartAffine, desCartAffine);
                    cout << "Screw computed: \n" << relativeDesiredScrew.transpose() << endl;

                    for (int i = 0; i < FRI_CART_VEC; i++)
                    {
                        quinticTrajCoeffs[i] = utils::computeQuinticTrajectory(relativeDesiredScrew(i), planTime);
                    }
                
                    planning = false;
                }
                
                if (kukaInterface.controlState())
                {
                    controlTime += sampleTime;
                    
                    if (controlTime <= planTime)
                    {
                        for (int i = 0; i < FRI_CART_VEC; i++)
                        {
                            relativeNextDesiredScrew(i) = utils::nextStepQuintic(quinticTrajCoeffs[i], controlTime);
                        }
                        
                    }
                    
                    relativeNextDesiredAffine = utils::getAffineFromScrew(relativeNextDesiredScrew);
                    cout << "Relative next desired screw converted into an affine matrix: \n" << relativeNextDesiredAffine.matrix() << endl;
                    nextDesiredAffine = initialCartAffine * relativeNextDesiredAffine;
                    msrdCartAffine = kukaInterface.getMsrdCartPosition();
                    relativeScrewError = utils::getScrewErrorFromAffine(msrdCartAffine, nextDesiredAffine);
//                        logger->add("relativeScrewError", relativeScrewError);
                    measToWorldRotation.topLeftCorner(3,3) = msrdCartAffine.linear();
                    measToWorldRotation.bottomRightCorner(3,3) = msrdCartAffine.linear();
                    cout << "Measured to world rotation: \n" << measToWorldRotation << endl;
                    screwError = measToWorldRotation * relativeScrewError;
//                         logger->add(screwErrorName, screwError);
                    if (firstLoop)                   
                    {
                        screwErrorPrev = screwError;
                        firstLoop = false;
                    }
                    desTwistError = (screwError - screwErrorPrev) / kukaInterface.getSampleTime();
                    cout << "Velocity error:" << desTwistError.transpose() << endl;
//                         logger->add("desTwistError", desTwistError);
                    twistError = smoothingFactor * desTwistError + (1 - smoothingFactor)*twistError;
//                         logger->add("twistError", twistError);
                    screwErrorPrev = screwError;
                    cout << "Velocity error filtered:" << twistError.transpose() << endl;
                    desiredForce = utils::computeForceImpedenceControl(screwError, twistError, myStiffness, myDamping);
                    desiredTorque = kukaInterface.calcKukaJacob().transpose() * desiredForce;    
                }
                else
                    controlTime = 0.;           
                
                kukaInterface.doJntImpedenceControl(desiredTorque);
                break;

                default:
                kukaInterface.dataExchange();
                break;
        }
        
        // Stop request is issued from the other side
        if ( kukaInterface.onStop() == -1) 
        {
            cout << "Leaving \n";
            break;
        }
    }
    
    return EXIT_SUCCESS;
}