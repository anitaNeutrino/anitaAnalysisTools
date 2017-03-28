/* -*- C++ -*-.*********************************************************************************************
 Author: Ben Strutt
 Email: strutt@physics.ucla.edu

 Description:
             A Cross Correlator to interact with the ROOTified ANITA-3 data and do some interferometry.
***********************************************************************************************************/

#ifndef INTERFEROMETRIC_MAP_H
#define INTERFEROMETRIC_MAP_H

// #define NUM_BINS_THETA 60
#define NUM_BINS_THETA 35
#define NUM_BINS_PHI 9
// #define THETA_RANGE 150
#define MIN_THETA -55
#define MAX_THETA 35
#define PHI_RANGE 22.5



#include "TH2D.h"



class InterferometricMap : public TH2D {

public:
  InterferometricMap();
  InterferometricMap(TString name, TString title, Int_t nBinsPhi, Double_t phiMin, Double_t phiMax, Int_t nBinsTheta, Double_t minTheta, Double_t maxTheta);
  InterferometricMap(TString name, TString title, Double_t phiMin); // constructor for coarse map, don't need extra info

  

protected:
  bool thetaAxisInSinTheta;
  void initializeInternals();
};




#endif // INTERFEROMETRIC_MAP_H
