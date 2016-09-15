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
#include "CrossCorrelator.h"

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

  static const double maxVoltsLimit = 2000;
  static const double minVoltsLimit = -2000;
  static const double absMaxMinSumLimit = 500;
  Status_t applySurfSaturationCut(Double_t maxVolts[][NUM_SEAVEYS], Double_t minVolts[][NUM_SEAVEYS], Double_t& maxMaxVolts, Double_t& minMinVolts, Double_t& absSumMaxMin);
  Status_t applySurfSaturationCutBetter(Double_t theMaxVolts, Double_t theMinVolts, Double_t absSumMaxMin);


  static const double ratioCutHigh = 2.8;
  static const double ratioCutLow = 1.14;
  // Step 1: cut self triggered blasts (this is almost a data quality cut)
  Status_t applyBottomToTopRingPeakToPeakRatioCut(AnitaPol::AnitaPol_t pol, Double_t* peakToPeak, Double_t& maxRatio);
  Status_t applyBottomToTopRingPeakToPeakRatioCut(Double_t maxPeakToPeakRatio);


  static const int maxAbsDeltaPhiSect = 2;
  Status_t L3TriggerDirectionCut(AnitaPol::AnitaPol_t pol, RawAnitaHeader* header, Double_t recoPhiDeg, Int_t& deltaPhiSect);


  static const double deltaSolarPhiDegCut = 20;
  Status_t applySunPointingCut(Double_t deltaSolarPhiDeg);


  // these variables all come from the output of things in the defineThermalCut subfolder
  static const int numFisherWeights = 3;
  // 0.105275	0.000663689
  static const Double_t fisherWeights[numFisherWeights] = {-2.77328, 14.5089, 0.0104461};
  // static const Double_t fisherWeights[numFisherWeights] = {-2.80469, 15.9061, 0.00742922};
  // static const Double_t fisherWeights[numFisherWeights] = {-2.80469, 15.9061, 0.0074292};
  // static const Double_t fisherWeights[numFisherWeights] = {-2.81448, 15.7929, 0.00783315};
  static const Double_t fisherCutVal = 0.000950146; //-0.0101151; //-0.64206; //-0.526251;
  Status_t applyThermalBackgroundCut(Double_t imagePeak, Double_t hilbertPeak, Double_t& fisher);

  static const Double_t maxImagePeakRatio = 0.5;
  Status_t applyImagePeakRatioCut(Double_t p1, Double_t p2, Double_t& peakRatio);

  static const Double_t maxThetaDeg = 0;
  static const Double_t minThetaDeg = -30;
  Status_t applyThetaAngleCut(Double_t thetaDeg);


};



#endif //ANALYSISCUTS
