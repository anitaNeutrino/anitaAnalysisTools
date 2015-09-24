/* -*- C++ -*-.***************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description: 
             A Cross Correlator to interact with the ROOTified ANITA-3 data and do some interferometry. 
*************************************************************************************************************** */

#ifndef CROSSCORRELATOR_H
#define CROSSCORRELATOR_H

// Ryan things
#include "UsefulAnitaEvent.h"
#include "AnitaGeomTool.h"
// #include "FFTtools.h"
#include "FancyFFTs.h"

// My things
#include "RootTools.h"

// ROOT things
#include "TGraph.h"
#include "TH2D.h"
#include "Math/Interpolator.h"
#include "Math/InterpolationTypes.h"

// standard c++ things
#include <iostream>
#include <assert.h>


// Offline reconstruction definitions
#define NUM_COMBOS 336

// Image definitions
#define NUM_BINS_THETA 256
#define NUM_BINS_PHI 64
#define THETA_RANGE 150
#define PHI_RANGE 22.5
#define NUM_SAMPLES 256

// Anita Geometry definitions, shouldn't really be here
#define NUM_POL 2
#define NUM_RING 3
#define DELTA_PHI_SECT 2

#define SPEED_OF_LIGHT 2.99792458e8


/*! \class CrossCorrelator
\brief Pass CrossCorrelator an event and get interferometric maps with a single function.

Does all the heavy lifting of getting waveforms from a UsefulAnitaEvent, cross correlating them and producing interferometric maps. 
*/


class CrossCorrelator : public TObject{

public:
  /**********************************************************************************************************
  Constructor and destructor functions
  **********************************************************************************************************/
  CrossCorrelator(Int_t upsampleFactorTemp = 1);
  ~CrossCorrelator();
  void initializeVariables(Int_t upsampleFactorTemp);



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
  void doAllCrossCorrelations();
  Double_t* crossCorrelateFourier(TGraph* gr1, TGraph* gr2);
  std::vector<std::vector<Double_t> > getMaxCorrelationTimes();
  std::vector<std::vector<Double_t> > getMaxCorrelationValues();

  /**********************************************************************************************************
  Calculate deltaT between two antennas (for a plane wave unless function name says otherwise)
  **********************************************************************************************************/
  Int_t getDeltaTExpected(Int_t ant1, Int_t ant2,Double_t phiWave, Double_t thetaWave);
  Int_t getDeltaTExpectedSpherical(Int_t ant1, Int_t ant2,Double_t phiWave, Double_t thetaWave, Double_t rWave);

  /**********************************************************************************************************
  Precalculate DeltaTs during initialization where appropriate
  **********************************************************************************************************/
  void do5PhiSectorCombinatorics();
  void fillDeltaTLookup();

  void writeDeltaTsFile(); // Speed up initialization
  Int_t readDeltaTsFile(); // Speed up initialization



  /**********************************************************************************************************
  Image generation functions.
  **********************************************************************************************************/
  
  Double_t getPhi0();
  Bool_t useCombo(Int_t ant1, Int_t ant2, Int_t phiSector);
  TH2D* makeBlankImage(TString name, TString title);
  TH2D* makeImage(AnitaPol::AnitaPol_t pol, UInt_t l3Trigger=0xffff);
  TH2D* makeImage(AnitaPol::AnitaPol_t pol, Double_t& imagePeak, 
		  Double_t& peakPhiDeg, Double_t& peakThetaDeg, UInt_t l3Trigger=0xffff);
  TH2D* makeImageSpherical(AnitaPol::AnitaPol_t pol, Double_t rWave);
  TH2D* makeImageSpherical(AnitaPol::AnitaPol_t pol, Double_t rWave, Double_t& imagePeak, 
			   Double_t& peakPhiDeg, Double_t& peakThetaDeg);
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
  UInt_t eventNumber; ///< For tracking event number
  UInt_t lastEventNormalized; ///< Prevents cross-correlation of the same event twice
  Double_t nominalSamplingDeltaT; ///< ANITA-3 => 1./2.6 ns
  Double_t upsampleFactor; ///< Default = 1, 
  Double_t correlationDeltaT; ///< nominalSamplingDeltaT/upsampleFactor, deltaT of interpolation and cross correlation
  Int_t numSamplesUpsampled; ///< Number of samples in waveform after up sampling is applied
  Int_t numCombos;
  
  std::vector<Int_t> ant2s[NUM_SEAVEYS]; ///< Vector holding ant2 indices for each ant1 (for combinatorics)
  std::vector<Int_t> comboToAnt1s; ///< Vector mapping combination index to ant1
  std::vector<Int_t> comboToAnt2s; ///< Vector mapping combination index to ant1
  Int_t comboIndices[NUM_SEAVEYS][NUM_SEAVEYS]; ///< Array mapping ant1+ant2 to combo index
  Double_t* crossCorrelations[NUM_POL][NUM_COMBOS]; ///< Arrays for cross correlations
  TGraph* grs[NUM_POL][NUM_SEAVEYS]; ///< Raw waveforms obtained from a UsefulAnitaEvent
  TGraph* grsInterp[NUM_POL][NUM_SEAVEYS]; ///< Interpolated TGraphs 
  Double_t interpRMS[NUM_POL][NUM_SEAVEYS]; ///< RMS of interpolation
  std::vector<Double_t> rArray; ///< Vector of antenna radial positions
  std::vector<Double_t> phiArrayDeg; ///< Vector of antenna azimuth positions
  std::vector<Double_t> zArray; ///< Vector of antenna heights
  Char_t deltaTs[NUM_COMBOS][NUM_PHI*NUM_BINS_PHI][NUM_BINS_THETA]; ///< Lookup of deltaTs between antenna pairs for making an image (UChar_t to reduce size)

  ClassDef(CrossCorrelator, 0);
};
#endif
