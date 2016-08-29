/* -*- C++ -*-.*********************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description: 
             My ANITA-3 analysis cuts
***********************************************************************************************************/

#ifndef ANALYSISCUTS_H
#define ANALYSISCUTS_H

#include "FancyFFTs.h"
#include "AnitaConventions.h"
#include "RawAnitaHeader.h"



// #include "CrossCorrelator.h"

/**
 * @namespace Analysis Cuts 
 * @brief Put my analysis cuts in one place for multiple scripts.
 * 
 *
 */
namespace AnalysisCuts{

  
  typedef enum EStatus {
    kPass = 0,
    kFail = 1
  } Status_t; // to make the scripts easy to read with comparisons to pass/fail

  // CUT FLOW
  // Aiming for combined reduction of factor O(1e9) for thermal noise events
  // How many signal events will be left?
  // here we go!

  const double ratioCutHigh = 2.8;
  const double ratioCutLow = 1.14;  
  // Step 1: cut self triggered blasts (this is almost a data quality cut)
  Status_t applyBottomToTopRingPeakToPeakRatioCut(AnitaPol::AnitaPol_t pol, Double_t* peakToPeak, Double_t& maxRatio);
  
  

  const int maxAbsDeltaPhiSect = 2;
  Status_t L3TriggerDirectionCut(AnitaPol::AnitaPol_t pol, RawAnitaHeader* header, Double_t recoPhiDeg, Int_t& deltaPhiSect);

  
  const double deltaSolarPhiDegCut = 20;
  Status_t applySunPointingCut(Double_t deltaSolarPhiDeg);


  // these variables all come from the output of things in the defineThermalCut subfolder
  const int numFisherWeights = 3;
  const Double_t fisherWeights[numFisherWeights] = {-2.81448, 15.7929, 0.00783315};
  const Double_t fisherCutVal = -0.526251;
  Status_t applyThermalBackgroundCut(Double_t imagePeak, Double_t hilbertPeak, Double_t& fisher);
  

};



#endif //ANALYSISCUTS
