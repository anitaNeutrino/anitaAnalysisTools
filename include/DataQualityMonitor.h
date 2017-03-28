/* -*- C++ -*-.*********************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description: 
             A class to take RawAnitaEvents and look for signs of bad data.
***********************************************************************************************************/

#ifndef DATAQUALITYMONITOR_H
#define DATAQUALITYMONITOR_H

// Ryan things
#include "UsefulAnitaEvent.h"
#include "AnitaEventCalibrator.h"
#include "AnitaGeomTool.h"
#include "UsefulAdu5Pat.h"

// standard c++ things
#include <iostream>

// ROOT things
#include "TChain.h"


// My things
#include "CrossCorrelator.h"


/**
 * @class DataQualityMonitor
 * @brief Looks for SURF saturation
 * 
 * A class to take RawAnitaEvents and look for signs of bad data.
*/
class DataQualityMonitor{

public:

  DataQualityMonitor();
  DataQualityMonitor(TChain* c);
  ~DataQualityMonitor();
  void setBranches(TChain* c);

  // tree variables
  TChain* dataQualityChain;
  Double_t maxAbsSecondDeriv[AnitaPol::kNotAPol][NUM_SEAVEYS];
  Double_t maxVolts[AnitaPol::kNotAPol][NUM_SEAVEYS];
  Int_t numPoints[AnitaPol::kNotAPol][NUM_SEAVEYS];
  UInt_t eventNumber;


  // quality variables
  Int_t phiAboveMaxVoltsThresh[NUM_PHI];
  Int_t numChannelsAboveSurfSaturation;
  Int_t numAboveVoltsBlastThresh;
  Int_t numPhiAboveMaxVoltsBlastThresh;


  Int_t processEntry(Long64_t entry, UInt_t eventNumberCheck = 0);  




  // cut variables, set in constructor, public for convenience... don't mess with them
  double maxVoltsBlastThresh;
  double saturationVolts;
  
private:


  
  
};
#endif

