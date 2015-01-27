/* -*- C++ -*-.***************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description: 
             A Cross Correlator to interact with the ROOTified ANITA-3 data and do some interferometry. 
*************************************************************************************************************** */

#ifndef CROSSCORRELATOR_H
#define CROSSCORRELATOR_H


/* Ryan things */
#include "UsefulAnitaEvent.h"
#include "AnitaGeomTool.h"
#include "FFTtools.h"


/* ROOT things */
#include "TGraph.h"
#include "TH2D.h"
#include "Math/Interpolator.h"
#include "Math/InterpolationTypes.h"


/* standard c++ things */
#include <iostream>
#include <array>
#include <assert.h>

/* GPU definitions */
#define GLOBAL_COMBOS_PER_PHI_GPU 12
#define LOCAL_COMBOS_PER_PHI_GPU 27

/* Offline reconstruction definitions */
#define NUM_COMBOS 336

#define NUM_BINS_THETA 256
#define NUM_BINS_PHI 64
#define THETA_RANGE 150
#define PHI_RANGE 22.5
#define NUM_SAMPLES 256

#define NUM_POL 2
#define NUM_RING 3


class CrossCorrelator : public TObject{

public:
  CrossCorrelator();
  ~CrossCorrelator();



  /* Cross correlation */
  Double_t correlationWithOffset(TGraph* gr1, TGraph* gr2, Int_t offset);
  void getNormalizedInterpolatedTGraphs(UsefulAnitaEvent* realEvent);
  void correlateEventGPU(UsefulAnitaEvent* realEvent);
  TGraph* getCrossCorrelationGraph(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2);
  void correlateEvent(UsefulAnitaEvent* realEvent);
  Double_t* crossCorrelateFourier(TGraph* gr1, TGraph* gr2);


  /* Interformetry */
  void preCalculateDeltaTsAndCombos();
  Int_t getDeltaTExpected(Int_t ant1, Int_t ant2,Double_t phiWave, Double_t thetaWave);

  short* fillDeltaTLookupGPU();
  TH2D* makeImageGPU(AnitaPol::AnitaPol_t pol);
  TH2D* makeImageGPU(AnitaPol::AnitaPol_t pol, UInt_t phiSectorMask);

  TH2D* makeImage(AnitaPol::AnitaPol_t pol, Double_t& imagePeak, Double_t& peakPhiDeg, Double_t& peakThetaDeg);
  TH2D* makeImage(AnitaPol::AnitaPol_t pol);

  Double_t findImagePeak(TH2D* hist, Double_t& imagePeakTheta, Double_t& imagePeakPhi);
  void do5PhiSectorCombinatorics();

  /* Waveform manipulation */
  TGraph* interpolateWithStartTime(TGraph* grIn, Double_t startTime);
  TGraph* normalizeTGraph(TGraph* gr);

  
  Double_t correlationDeltaT;
  short* offsetIndGPU;
  
  /* For GPU style*/
  Int_t ant1Gpu[NUM_PHI][LOCAL_COMBOS_PER_PHI_GPU];
  Int_t ant2Gpu[NUM_PHI][LOCAL_COMBOS_PER_PHI_GPU];
  Double_t correlationsGPU[NUM_POL][NUM_PHI][GLOBAL_COMBOS_PER_PHI_GPU][NUM_SAMPLES];

  /* For wider reconstruction */
  std::array<std::vector<int>, NUM_SEAVEYS> ant2s;
  std::array<std::array<int, NUM_SEAVEYS>, NUM_SEAVEYS> comboIndices;
  std::array<std::array<std::array<Double_t, NUM_SAMPLES>, NUM_COMBOS>, NUM_POL> crossCorrelations;
  std::array<std::array<int, NUM_COMBOS>, NUM_POL> doneCrossCorrelations;

  TGraph* grs[NUM_POL][NUM_SEAVEYS];
  TGraph* grsInterp[NUM_POL][NUM_SEAVEYS];


  std::array<Double_t,NUM_SEAVEYS> rArrayFeed = {0.9675,0.7402,0.9675,0.7402,0.9675,0.7402,0.9675,0.7402,
  					       0.9675,0.7402,0.9675,0.7402,0.9675,0.7402,0.9675,0.7402,
  					       2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,
  					       2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,
  					       2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,
  					       2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447};

  // /* Flipped top ring high/low -> low/high */
  // std::array<Double_t,NUM_SEAVEYS> rArrayFeed = {0.7402,0.9675,0.7402,0.9675,0.7402,0.9675,0.7402,0.9675,
  // 					       0.7402,0.9675,0.7402,0.9675,0.7402,0.9675,0.7402,0.9675,
  // 					       2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,
  // 					       2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,
  // 					       2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,
  // 					       2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447};

  std::array<Double_t,NUM_SEAVEYS> phiArrayDegFeed = {0.0,22.5,45.0,67.5,90.0,112.5,135.0,157.5,
  						    180.0,202.5,225.0,247.5,270.0,292.5,315.0,337.5,
  						    0.0,22.5,45.0,67.5,90.0,112.5,135.0,157.5,
  						    180.0,202.5,225.0,247.5,270.0,292.5,315.0,337.5,
  						    0.0,22.5,45.0,67.5,90.0,112.5,135.0,157.5,180.0,
  						    202.5,225.0,247.5,270.0,292.5,315.0,337.5};

  std::array<Double_t,NUM_SEAVEYS> zArrayFeed = {-1.4407,-2.4135,-1.4407,-2.4135,-1.4407,-2.4135,-1.4407,-2.4135,
  					       -1.4407,-2.4135,-1.4407,-2.4135,-1.4407,-2.4135,-1.4407,-2.4135,
  					       -5.1090,-5.1090,-5.1090,-5.1090,-5.1090,-5.1090,-5.1090,-5.1090,
  					       -5.1090,-5.1090,-5.1090,-5.1090,-5.1090,-5.1090,-5.1090,-5.1090,
  					       -6.1951,-6.1951,-6.1951,-6.1951,-6.1951,-6.1951,-6.1951,-6.1951,
  					       -6.1951,-6.1951,-6.1951,-6.1951,-6.1951,-6.1951,-6.1951,-6.1951};

  // /* Flipped top ring high/low -> low/high */
  // std::array<Double_t,NUM_SEAVEYS> zArrayFeed = {-2.4135,-1.4407,-2.4135,-1.4407,-2.4135,-1.4407,-2.4135,-1.4407,
  // 					       -2.4135,-1.4407,-2.4135,-1.4407,-2.4135,-1.4407,-2.4135,-1.4407,
  // 					       -5.1090,-5.1090,-5.1090,-5.1090,-5.1090,-5.1090,-5.1090,-5.1090,
  // 					       -5.1090,-5.1090,-5.1090,-5.1090,-5.1090,-5.1090,-5.1090,-5.1090,
  // 					       -6.1951,-6.1951,-6.1951,-6.1951,-6.1951,-6.1951,-6.1951,-6.1951,
  // 					       -6.1951,-6.1951,-6.1951,-6.1951,-6.1951,-6.1951,-6.1951,-6.1951};


  std::array<Double_t,NUM_SEAVEYS> rArray={0.9675,0.7402,0.9675,0.7402,0.9675,0.7402,0.9675,0.7402,
					 0.9675,0.7402,0.9675,0.7402,0.9675,0.7402,0.9675,0.7402,
					 2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,
					 2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,
					 2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,
					 2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447};

  std::array<Double_t,NUM_SEAVEYS> phiArrayDeg = {0.0,22.5,45.0,67.5,90.0,112.5,135.0,157.5,
						180.0,202.5,225.0,247.5,270.0,292.5,315.0,337.5,
						0.0,22.5,45.0,67.5,90.0,112.5,135.0,157.5,
						180.0,202.5,225.0,247.5,270.0,292.5,315.0,337.5,
						0.0,22.5,45.0,67.5,90.0,112.5,135.0,157.5,180.0,
						202.5,225.0,247.5,270.0,292.5,315.0,337.5};

  std::array<Double_t,NUM_SEAVEYS> zArray = {-1.4407,-2.4135,-1.4407,-2.4135,-1.4407,-2.4135,-1.4407,-2.4135,
					   -1.4407,-2.4135,-1.4407,-2.4135,-1.4407,-2.4135,-1.4407,-2.4135,
					   -5.1090,-5.1090,-5.1090,-5.1090,-5.1090,-5.1090,-5.1090,-5.1090,
					   -5.1090,-5.1090,-5.1090,-5.1090,-5.1090,-5.1090,-5.1090,-5.1090,
					   -6.1951,-6.1951,-6.1951,-6.1951,-6.1951,-6.1951,-6.1951,-6.1951,
					   -6.1951,-6.1951,-6.1951,-6.1951,-6.1951,-6.1951,-6.1951,-6.1951};


  std::array<Double_t, NUM_BINS_THETA> tanThetaLookup;
  std::array<Double_t, NUM_BINS_THETA> cosThetaLookup;
  std::array<Double_t, NUM_BINS_PHI*NUM_PHI> sinPhiWaveLookup;
  std::array<Double_t, NUM_BINS_PHI*NUM_PHI> cosPhiWaveLookup;
  std::array<Double_t, NUM_SEAVEYS> cosPhiArrayLookup;
  std::array<Double_t, NUM_SEAVEYS> sinPhiArrayLookup;

  Int_t feedOffsetStep = 0;
  Double_t angleDeg = -10;
  Double_t dlPerStep = 0.1;
  static const Int_t numSteps = 1;
  Int_t getDeltaTExpected(Int_t ant1, Int_t ant2, Int_t phiBin, Int_t thetaBin); /* Slightly faster? */
  void fillDeltaTLookupWithOffsets();
  std::array<std::array<std::array<std::array<unsigned char, NUM_BINS_THETA>, NUM_PHI*NUM_BINS_PHI>, NUM_COMBOS>, numSteps>deltaTsVaryingPosition;

};
#endif
