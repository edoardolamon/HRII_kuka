/**
 \author Edoardo Lamon
        \file utils_test.cpp
    \brief Test for utils library
*/

#include "utils.h"
#include "XBotLogger/Logger.hpp"

int main(int argc, char** argv)
{
    XBot::MatLogger::Ptr logger = XBot::MatLogger::getLogger("/tmp/utils_test_log");

//    cout << "------------------test relativeScrewError------------------" << endl;
//    Eigen::Affine3f testAffine, unitAffine;
//    Eigen::Matrix3f testRot;
//    Eigen::VectorXf relativeScrewError(FRI_CART_VEC);

//    unitAffine.translation().setZero();
//    unitAffine.linear().setIdentity(3,3);
//    testAffine.translation().setZero();
//    testRot << 1, 0, 0,
//               0, 1, 0,
//               0, 0, 1;
//    testAffine.linear() = testRot;
//    relativeScrewError = utils::getScrewErrorFromAffine(unitAffine, testAffine);
//    cout << "Input rotation: \n" << testAffine.linear() << endl;
//    cout << endl;
//    cout << "Output screw rotation: \n" << relativeScrewError.segment(3,3) << endl;

//    testRot <<-1, 0, 0,
//               0, 1, 0,
//               0, 0,-1;
//    testAffine.linear() = testRot;
//    relativeScrewError = utils::getScrewErrorFromAffine(unitAffine, testAffine);
//    cout << "Input rotation: \n" << testAffine.linear() << endl;
//    cout << endl;
//    cout << "Output screw rotation: \n" << relativeScrewError.segment(3,3) << endl;

//    testRot <<-1, 0, 0,
//               0,-1, 0,
//               0, 0, 1;
//    testAffine.linear() = testRot;
//    relativeScrewError = utils::getScrewErrorFromAffine(unitAffine, testAffine);
//    cout << "Input rotation: \n" << testAffine.linear() << endl;
//    cout << endl;
//    cout << "Output screw rotation: \n" << relativeScrewError.segment(3,3) << endl;
//    cout << "-----------------------------------------------------------" << endl;

//    cout << "------------------test getAffineFromScrew------------------" << endl;
//    Eigen::VectorXf testScrew(FRI_CART_VEC);
//    Eigen::Affine3f affineFromScrew;
//    testScrew << 0.,0.,0.,M_PI/4,0.,0.;
//    affineFromScrew = utils::getAffineFromScrew(testScrew);
//    cout << "Input screw rotation: \n" << testScrew.segment(3,3) << endl;
//    cout << endl;
//    cout << "Output rotation: \n" << affineFromScrew.linear() << endl;

//    testScrew << 0.,0.,0.,M_PI/sqrt(3),M_PI/sqrt(3),M_PI/sqrt(3);
//    affineFromScrew = utils::getAffineFromScrew(testScrew);
//    cout << "Input screw rotation: \n" << testScrew.segment(3,3) << endl;
//    cout << endl;
//    cout << "Output rotation: \n" << affineFromScrew.linear() << endl;

//    testScrew << 0.,0.,0.,2*M_PI/sqrt(3),2*M_PI/sqrt(3),2*M_PI/sqrt(3);
//    affineFromScrew = utils::getAffineFromScrew(testScrew);
//    cout << "Input screw rotation: \n" << testScrew.segment(3,3) << endl;
//    cout << endl;
//    cout << "Output rotation: \n" << affineFromScrew.linear() << endl;
//    cout << "-----------------------------------------------------------" << endl;

//     cout << "------------------test computeQuinticTrajectory------------------" << endl;
//     float initialPos = 0.;
//     float desPos = 1.;
//     vector<float> coeffs;
//     coeffs = utils::computeQuinticTrajectory(desPos);
//     logger->log("quinticCoeffs", coeffs);
//    cout << "Quintic coeffs :" << endl;
//    for (int i=0; i < 6; i++)
//    {
//        cout << "i: " << i;
//        cout << "; coeff: " << coeffs[i] << endl;
//    }
//    cout << endl;
//     cout << "-----------------------------------------------------------" << endl;

//     cout << "------------------test nextStepQuintic------------------" << endl;
    const int trajSize = 1000;
//     Eigen::VectorXf trajectory(trajSize);
    const float deltaT = 0.005;
// float nextPos = initialPos;
    float time = 0;
//     for (int i = 0; i < trajSize; i++)
//     {
//         logger->add("trajectory", nextPos);
//         trajectory(i) = nextPos;
//         time += deltaT;
//         if (abs(desPos  - trajectory(i)) > 0.0001)
//         {
//             nextPos = utils::nextStepQuintic(coeffs,time);
//         }
//     }
// 
//     logger->log("trajectory_batch", trajectory);

//     cout << "-----------------------------------------------------------" << endl;
    
    cout << "------------------test nextStepCircular------------------" << endl;
    logger->createMatrixVariable("trajectory",4,4);
    logger->createVectorVariable("trajectory_screw",6);
    const float radius = 0.05;
    const float period = 5.;
    Eigen::Vector3f center;
    center << radius,0,0.;
    Eigen::VectorXf nextCircScrew(6);
    Eigen::Affine3f nextCircAffine, initialCartAffine;
    initialCartAffine.linear().setIdentity(3,3);
    initialCartAffine.translation().setZero();
    initialCartAffine.translation() << 1., 2., 3.;
    nextCircScrew = utils::getScrewErrorFromAffine(initialCartAffine, initialCartAffine);
    
    logger->add("trajectory",initialCartAffine.matrix());
    logger->add("trajectory_screw", nextCircScrew);
    
    for (int i = 0; i < trajSize; i++)
    {
        time += deltaT;
        nextCircAffine = utils::nextStepCircular(initialCartAffine,time,period,radius,center);
	nextCircScrew = utils::getScrewErrorFromAffine(initialCartAffine, nextCircAffine);
	logger->add("trajectory", nextCircAffine.matrix());
	logger->add("trajectory_screw", nextCircScrew);
    }

    
    cout << "-----------------------------------------------------------" << endl;
    
    logger->flush();
    return 0;
}


