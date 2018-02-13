#ifndef KUKALWRINTERFACE_H
#define KUKALWRINTERFACE_H

#include <eigen3/Eigen/Dense>

#include "friudp.h"
#include "friremote.h"
#include "utils.h"
#include "XBotLogger/Logger.hpp"

#ifndef M_PI 
#define M_PI 3.14159
#endif

#ifndef AFFINE_DIM 
#define AFFINE_DIM 4
#endif

class KukaLWRInterface
{

  public:
      
      KukaLWRInterface(int kukaPort, char * kukaAddress);
  
      ~KukaLWRInterface();
      
      void loggerInit();
      
      void positionInit();
      
      double getSampleTime();
      
      int getCurrentControlScheme();
      
      Eigen::VectorXf getMsrdJntPosition();
      
      Eigen::VectorXf getCmdJntPosition();
      
      Eigen::Affine3f getMsrdCartPosition();
      
      Eigen::Affine3f getCmdCartPosition();
    
//       void updatePosition();
      
      void dataExchange();
      
      void communicationQualityCheck();
      
      void onStart();
      
      int onStop();
      
      bool controlState();
      
      Eigen::MatrixXf calcKukaJacob();
      
      void doJntImpedenceControl(Eigen::VectorXf desiredTorque);
    
    
  protected:
      
      friRemote _friInst;
      
      XBot::MatLogger::Ptr _logger;
      
      FRI_QUALITY _lastQuality = FRI_QUALITY_BAD;
      FRI_CTRL _lastCtrlScheme = FRI_CTRL_OTHER;
      
      Eigen::VectorXf _msrdJntPos, _cmdJntPos;
      Eigen::Affine3f _msrdCartPos, _cmdCartPos, _nextDesCartPos;
      
      double _sampleTime;
      
      
  
};

#endif // KUKALWRINTERFACE_H
