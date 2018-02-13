#include "../include/KukaLWRInterface.h"

KukaLWRInterface::KukaLWRInterface(int kukaPort, char * kukaAddress):
_friInst(kukaPort, kukaAddress),
_msrdJntPos(LBR_MNJ), _cmdJntPos(LBR_MNJ)
{   
    cout << "Opening FRI Version " 
	 << FRI_MAJOR_VERSION << "." << FRI_SUB_VERSION << "." << FRI_DATAGRAM_ID_CMD << "." << FRI_DATAGRAM_ID_MSR 
	 << " Interface for First Sample" << endl;
	 
    positionInit();
    loggerInit();
    _sampleTime = 0.;
}

KukaLWRInterface::~KukaLWRInterface()
{

}

void KukaLWRInterface::loggerInit()
{
    _logger = XBot::MatLogger::getLogger("/tmp/kuka_hrii_log");
    _logger->createMatrixVariable("nextDesiredAffine", AFFINE_DIM, AFFINE_DIM);
    _logger->createMatrixVariable("msrdCartAffine", AFFINE_DIM, AFFINE_DIM);
    _logger->createVectorVariable("screwError", FRI_CART_VEC);
    _logger->createVectorVariable("relativeScrewError", FRI_CART_VEC);
    _logger->createVectorVariable("desForce", FRI_CART_VEC);
}

void KukaLWRInterface::positionInit()
{
    _msrdJntPos.setZero();
    _cmdJntPos.setZero();
    _msrdCartPos.linear().setIdentity();
    _msrdCartPos.translation().setZero();
    _cmdCartPos.linear().setIdentity();
    _cmdCartPos.translation().setZero();
}

double KukaLWRInterface::getSampleTime()
{
	_friInst.doDataExchange();
	_sampleTime = _friInst.getSampleTime();
// 	_logger->log("sampleTime", _sampleTime);
    
    return _sampleTime;
}

Eigen::VectorXf KukaLWRInterface::getMsrdJntPosition()
{
    for (int i = 0; i < LBR_MNJ; i++)
    {
        _msrdJntPos(i) = _friInst.getMsrMsrJntPosition()[i];
    }
}

Eigen::VectorXf KukaLWRInterface::getCmdJntPosition()
{
    for (int i = 0; i < LBR_MNJ; i++)
    {
        _cmdJntPos(i) = _friInst.getMsrCmdJntPosition()[i];
    }
}

Eigen::Affine3f KukaLWRInterface::getMsrdCartPosition()
{
    for (int i = 0; i < FRI_CART_FRM_DIM; i++)
    {
        int row = i / AFFINE_DIM;
        int col = i % AFFINE_DIM;
        _msrdCartPos(row,col) = _friInst.getMsrCartPosition()[i];
    }
}

Eigen::Affine3f KukaLWRInterface::getCmdCartPosition()
{
    for (int i = 0; i < FRI_CART_FRM_DIM; i++)
    {
        int row = i / AFFINE_DIM;
        int col = i % AFFINE_DIM;
        _cmdCartPos(row,col) = _friInst.getMsrCmdCartPosition()[i];
    }
}


// void KukaLWRInterface::updatePosition()
// {
//     // Read the current joint position
//     for (int i = 0; i < LBR_MNJ; i++)
//     {
//         _cmdJntPos(i) = _friInst.getMsrCmdJntPosition()[i];
//         _msrdJntPos(i) = _friInst.getMsrMsrJntPosition()[i];
//     }
// 
//     // Read the current cartesian position
//     for (int i = 0; i < FRI_CART_FRM_DIM; i++)
//     {
// 		int row = i / AFFINE_DIM;
// 		int col = i % AFFINE_DIM;
// 		_cmdCartPos(row,col) = _friInst.getMsrCmdCartPosition()[i];
// 		_msrdCartPos(row,col) = _friInst.getMsrCartPosition()[i];
//     }
// }

int KukaLWRInterface::getCurrentControlScheme()
{
	return _friInst.getCurrentControlScheme();
}

void KukaLWRInterface::dataExchange()
{
    _friInst.doDataExchange();
}

void KukaLWRInterface::communicationQualityCheck()
{
    /// perform some arbitrary handshake to KRL -- possible in monitor mode already
    // send to krl int a value
    _friInst.setToKRLInt(0,1);
    if ( _friInst.getQuality() >= FRI_QUALITY_OK)
    {
	// send a second marker
	_friInst.setToKRLInt(0,10);
    }
    
    // just mirror the real value
    _friInst.setToKRLReal(0,_friInst.getFrmKRLReal(1));

    if ( _lastCtrlScheme != _friInst.getCurrentControlScheme())
    {
	cout << "Switching control scheme " << _lastCtrlScheme;
	_lastCtrlScheme = _friInst.getCurrentControlScheme();
	cout << " to " << _lastCtrlScheme << endl;
    }
}

void KukaLWRInterface::onStart()
{

}

int KukaLWRInterface::onStop()
{
	//
	// Quality change leads to output of statistics
	// for informational reasons
	//
    if ( _friInst.getQuality() != _lastQuality)
    {
		cout << "quality change detected "<< _friInst.getQuality()<< " \n";
        cout << _friInst.getMsrBuf().intf;
        cout << endl;
        _lastQuality = _friInst.getQuality();
    }
    
    return _friInst.getFrmKRLInt(0);
}

bool KukaLWRInterface::controlState()
{
    return ( _friInst.isPowerOn() & _friInst.getState() == FRI_STATE_CMD);
}

Eigen::MatrixXf KukaLWRInterface::calcKukaJacob()
{
        float C1,S1,S2,C2,S3,C3,S4,C4,S5,C5,S6,C6,S7,C7;
        float l1=0.310, l2=0.4, l3=0.39, l4=0.075;
        Eigen::MatrixXf J_kuka(FRI_CART_VEC,LBR_MNJ);

        C1=cos(getMsrdJntPosition()(0));
        S1=sin(getMsrdJntPosition()(0));
        C2=cos(getMsrdJntPosition()(1)+M_PI/2);
        S2=sin(getMsrdJntPosition()(1)+M_PI/2);
        C3=cos(getMsrdJntPosition()(2));
        S3=sin(getMsrdJntPosition()(2));
        C4=cos(getMsrdJntPosition()(3));
        S4=sin(getMsrdJntPosition()(3));
        C5=cos(getMsrdJntPosition()(4));
        S5=sin(getMsrdJntPosition()(4));
        C6=cos(getMsrdJntPosition()(5));
        S6=sin(getMsrdJntPosition()(5));
        C7=cos(getMsrdJntPosition()(6));
        S7=sin(getMsrdJntPosition()(6));

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

void KukaLWRInterface::doJntImpedenceControl(Eigen::VectorXf desiredTorque)
{
    float nullJntStiffness[LBR_MNJ]={0,0,0,0,0,0,0};
    float commandedTorque[LBR_MNJ]; // = {0.,0.,0.,0.,0.,0.,0.};
    for (int i = 0; i < LBR_MNJ; i++)
    {
        commandedTorque[i] = desiredTorque(i);
    }
    _friInst.doJntImpedanceControl(NULL, nullJntStiffness, NULL, commandedTorque);
}







