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
#include "FancyFFTs.h"
#include "RootTools.h"

/* ROOT things */
#include "TGraph.h"
#include "TH2D.h"
#include "Math/Interpolator.h"
#include "Math/InterpolationTypes.h"

/* standard c++ things */
#include <iostream>
#include <assert.h>

/* GPU definitions */
#define GLOBAL_COMBOS_PER_PHI_GPU 12
#define LOCAL_COMBOS_PER_PHI_GPU 27

/* Offline reconstruction definitions */
#define NUM_COMBOS 336
//#define NUM_COMBOS 480
// #define NUM_COMBOS 672

/* Image definitions*/
#define NUM_BINS_THETA 256
#define NUM_BINS_PHI 64
#define THETA_RANGE 150
#define PHI_RANGE 22.5
#define NUM_SAMPLES 256

/* Anita Geometry definitions, shouldn't really be here */
#define NUM_POL 2
#define NUM_RING 3

#define SPEED_OF_LIGHT 2.99792458e8

class CrossCorrelator : public TObject{

public:
  /**********************************************************************************************************
  Constructor and destructor functions
  **********************************************************************************************************/
  CrossCorrelator();
  ~CrossCorrelator();


  /**********************************************************************************************************
  Waveform manipulation functions
  **********************************************************************************************************/
  void getNormalizedInterpolatedTGraphs(UsefulAnitaEvent* realEvent);
  TGraph* interpolateWithStartTime(TGraph* grIn, Double_t startTime);


  /**********************************************************************************************************
  All correlation functions
  **********************************************************************************************************/

  Double_t correlationWithOffset(TGraph* gr1, TGraph* gr2, Int_t offset);
  void correlateEvent(UsefulAnitaEvent* realEvent);
  void correlateEventGPU(UsefulAnitaEvent* realEvent);
  void doAllCrossCorrelations();
  Double_t* crossCorrelateFourier(TGraph* gr1, TGraph* gr2);


  /**********************************************************************************************************
  Calculate deltaT between two antennas (for a plane wave unless function name says otherwise)
  **********************************************************************************************************/
  Int_t getDeltaTExpected(Int_t ant1, Int_t ant2,Double_t phiWave, Double_t thetaWave);
  Int_t getDeltaTExpectedSpherical(Int_t ant1, Int_t ant2,Double_t phiWave, Double_t thetaWave, Double_t rWave);    


  Int_t getDeltaTExpected(Int_t ant1, Int_t ant2, Int_t phiBin, Int_t thetaBin); /* Slightly faster? */


  /**********************************************************************************************************
  Precalculate DeltaTs during initialization where appropriate
  **********************************************************************************************************/
  void do5PhiSectorCombinatorics();
  void fillDeltaTLookup();
  short* fillDeltaTLookupGPU();

  void writeDeltaTsFile(); /* Attempt to speed up initialization */
  Int_t readDeltaTsFile(); /* Attempt to speed up initialization */



  /**********************************************************************************************************
  Image generation functions.
  **********************************************************************************************************/

  TH2D* makeImage(AnitaPol::AnitaPol_t pol, UInt_t l3Trigger=0xffff);
  TH2D* makeImage(AnitaPol::AnitaPol_t pol, Double_t& imagePeak, 
		  Double_t& peakPhiDeg, Double_t& peakThetaDeg, UInt_t l3Trigger=0xffff);
  TH2D* makeImageSpherical(AnitaPol::AnitaPol_t pol, Double_t rWave);
  TH2D* makeImageSpherical(AnitaPol::AnitaPol_t pol, Double_t rWave, Double_t& imagePeak, 
			   Double_t& peakPhiDeg, Double_t& peakThetaDeg);
  TH2D* makeImageGPU(AnitaPol::AnitaPol_t pol);
  TH2D* makeImageGPU(AnitaPol::AnitaPol_t pol, UInt_t phiSectorMask);
  Double_t findImagePeak(TH2D* hist, Double_t& imagePeakTheta, Double_t& imagePeakPhi);



  /**********************************************************************************************************
  Functions to delete pointers to internal variables
  **********************************************************************************************************/
  void deleteCrossCorrelations();
  void deleteAllWaveforms();




  /**********************************************************************************************************
  Functions for debugging or testing
  **********************************************************************************************************/
  void correlateEventTest(Double_t phiDegSource, Double_t thetaDegSource, Double_t rSource);
  TGraph* getCrossCorrelationGraph(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2);


  /**********************************************************************************************************
  Variables
  **********************************************************************************************************/
  UInt_t lastEventNormalized;
  UInt_t eventNumber;
  Double_t correlationDeltaT;
  short* offsetIndGPU;
  
  Int_t ant1Gpu[NUM_PHI][LOCAL_COMBOS_PER_PHI_GPU];
  Int_t ant2Gpu[NUM_PHI][LOCAL_COMBOS_PER_PHI_GPU];
  Double_t correlationsGPU[NUM_POL][NUM_PHI][GLOBAL_COMBOS_PER_PHI_GPU][NUM_SAMPLES];
  std::vector<Int_t> ant2s[NUM_SEAVEYS];
  int comboIndices[NUM_SEAVEYS][NUM_SEAVEYS];
  Double_t* crossCorrelations[NUM_POL][NUM_COMBOS];
  int doneCrossCorrelations[NUM_POL][NUM_COMBOS];
  TGraph* grs[NUM_POL][NUM_SEAVEYS];
  TGraph* grsInterp[NUM_POL][NUM_SEAVEYS];
  Double_t interpRMS[NUM_POL][NUM_SEAVEYS];
  Double_t rArray[NUM_SEAVEYS];
  Double_t phiArrayDeg[NUM_SEAVEYS];
  Double_t zArray[NUM_SEAVEYS];
  Double_t tanThetaLookup[NUM_BINS_THETA];
  Double_t cosThetaLookup[NUM_BINS_THETA];
  Double_t sinPhiWaveLookup[NUM_BINS_PHI*NUM_PHI];
  Double_t cosPhiWaveLookup[NUM_BINS_PHI*NUM_PHI];
  Double_t cosPhiArrayLookup[NUM_SEAVEYS];
  Double_t sinPhiArrayLookup[NUM_SEAVEYS];
  unsigned char deltaTs[NUM_COMBOS][NUM_PHI*NUM_BINS_PHI][NUM_BINS_THETA];

  ClassDef(CrossCorrelator, 0);
};
#endif
