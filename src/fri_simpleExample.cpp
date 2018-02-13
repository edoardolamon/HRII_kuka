/** 
 \author (Edoardo Lamon)
	\file
    \brief Simple example for understanding how to give torque command to a KUKA LW4 using FRI

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
//#include <eigen3/unsupported/Eigen/MatrixFunctions>

#ifndef M_PI 
#define M_PI 3.14159
#endif

using namespace std;

void _calcKukaJacob(Eigen::MatrixXf &J_kuka, float jointPos[])
{
    if (!(J_kuka.rows() == FRI_CART_VEC && J_kuka.cols() == LBR_MNJ))  // (6x7)
    {
        cerr << "Jacobian wrong dimensions. Jacobian should be 6x7" << endl;
        J_kuka(FRI_CART_VEC,LBR_MNJ);
    }

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
        Eigen::MatrixXf calcJacobian(FRI_CART_VEC,LBR_MNJ);
        float * msrdJacobianVector;
        Eigen::MatrixXf msrdJacobian(FRI_CART_VEC,LBR_MNJ);
        float smoothingFactor = 0.1;

//        float myJntStiffness[LBR_MNJ]={100,100,100,100,100,100,100};
        float nullJntStiffness[LBR_MNJ]={0,0,0,0,0,0,0};
        float newJntVals[LBR_MNJ];
        float msrdJntVals[LBR_MNJ];

        Eigen::VectorXf externalForce(FRI_CART_VEC);
                        externalForce << 10.0,0.0,0.0,0.0,0.0,0.0;
        Eigen::VectorXf externalTorque(LBR_MNJ);
        float initialTorque[LBR_MNJ] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};

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

                /// First loop computation
                if (firstLoop)
                {
//                    _calcKukaJacob(calcJacobian,msrdJntVals);
//                    cout << "Jacobian = " << calcJacobian << endl;
//                    msrdJacobianVector = friInst.getJacobian();
//                    for (int i = 0; i < 42; i++)
//                    {
//                        int col = i % 7;
//                        int row = i / 7;
//                        msrdJacobian(row,col) = msrdJacobianVector[i];
//                    }
//                    cout << "Measured jacobian = " << msrdJacobian << endl;
                    firstLoop = false;
                }

                if ( friInst.getState() == FRI_STATE_CMD)
                {
                    if ( friInst.isPowerOn() )
                    {
                        // Update time counter
                        timeCounter += deltaT;

                        // Compute jacobian
                        _calcKukaJacob(calcJacobian,msrdJntVals);
                        msrdJacobianVector = friInst.getJacobian();
                        for (int i = 0; i < 42; i++)
                        {
                            int col = i % 7;
                            int row = i / 7;
                            msrdJacobian(row,col) = msrdJacobianVector[i];
                        }
                        externalTorque = calcJacobian.transpose() * externalForce;
//                        externalTorque = msrdJacobian.transpose() * externalForce;
//                        cout << "Input torque : \n" << externalTorque.transpose() << endl;
                        for (int i = 0; i < LBR_MNJ; i++)
                        {
                            initialTorque[i] = smoothingFactor * externalTorque(i) + (1 - smoothingFactor)*initialTorque[i];
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

//                friInst.doJntImpedanceControl(NULL, NULL, NULL, NULL);
//                friInst.doJntImpedanceControl(NULL, nullJntStiffness, NULL, NULL);
                friInst.doJntImpedanceControl(NULL, nullJntStiffness, NULL, initialTorque);
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
/* @} */
