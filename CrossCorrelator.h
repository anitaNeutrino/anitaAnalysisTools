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
#include "TThread.h"

// standard c++ things
#include <iostream>


// Offline reconstruction definitions
#define NUM_COMBOS 336
#define NUM_THREADS 4
 
// Typical number of samples in waveform
#define NUM_SAMPLES 256

// Image definitions
#define NUM_BINS_THETA 150
#define NUM_BINS_PHI 25
#define THETA_RANGE 150
#define PHI_RANGE 22.5

#define NUM_BINS_THETA_ZOOM 64
#define NUM_BINS_PHI_ZOOM 64
#define THETA_RANGE_ZOOM 6.4
#define PHI_RANGE_ZOOM 6.4

// Anita Geometry definitions, shouldn't really be here
#define NUM_POL 2
#define NUM_RING 3
#define DELTA_PHI_SECT 2

#define SPEED_OF_LIGHT 2.99792458e8

#define ALL_PHI_TRIGS 0xffff




/*! \class CrossCorrelator
\brief Pass CrossCorrelator an event and get interferometric maps with a single function.

Does all the heavy lifting of getting waveforms from a UsefulAnitaEvent, cross correlating them and producing interferometric maps. 
*/


class CrossCorrelator : public TObject{

public:

  /**********************************************************************************************************
  typdef enums: flags for making maps
  **********************************************************************************************************/
  typedef enum {
    kGlobal,
    kTriggered,
    kNumMapModes
  } mapMode_t;

  typedef enum{
    kZoomedOut,
    kZoomedIn,
    kNumZoomModes
  } zoomMode_t;
  TString mapModeNames[kNumMapModes];
  TString zoomModeNames[kNumZoomModes];

  
  /**********************************************************************************************************
  Constructor and destructor functions
  **********************************************************************************************************/
  CrossCorrelator();
  ~CrossCorrelator();
  void initializeVariables();  
  void printInfo();
  
  /**********************************************************************************************************
  Waveform manipulation functions
  **********************************************************************************************************/
  void getNormalizedInterpolatedTGraphs(UsefulAnitaEvent* realEvent, AnitaPol::AnitaPol_t pol);
  void doFFTs(AnitaPol::AnitaPol_t pol);
  TGraph* interpolateWithStartTime(TGraph* grIn, Double_t startTime);

  /**********************************************************************************************************
  All correlation functions
  **********************************************************************************************************/

  void correlateEvent(UsefulAnitaEvent* realEvent);
  void correlateEvent(UsefulAnitaEvent* realEvent, AnitaPol::AnitaPol_t pol);
  Double_t* crossCorrelateFourier(TGraph* gr1, TGraph* gr2);
  Double_t* crossCorrelateFourier(Int_t numSamplesTimeDomain,
				  std::complex<Double_t>* fft1, std::complex<Double_t>* fft2,
				  Int_t threadInd=0);

  // Double_t* crossCorrelateFourier(FFTWComplex* fft1, FFTWComplex* fft2);    
  std::vector<std::vector<Double_t> > getMaxCorrelationTimes();
  std::vector<std::vector<Double_t> > getMaxCorrelationValues();
  std::vector<Double_t> getMaxCorrelationTimes(AnitaPol::AnitaPol_t pol);
  std::vector<Double_t> getMaxCorrelationValues(AnitaPol::AnitaPol_t pol);

  void doAllCrossCorrelations(AnitaPol::AnitaPol_t pol);
  void doAllCrossCorrelationsThreaded(AnitaPol::AnitaPol_t pol); // Launches threads

  void doUpsampledCrossCorrelations(AnitaPol::AnitaPol_t pol, UInt_t l3TrigPattern);
  void doUpsampledCrossCorrelationsThreaded(AnitaPol::AnitaPol_t pol, UInt_t l3TrigPattern); // Launches threads

  // Actual functions executed by thread
  static void* doSomeCrossCorrelationsThreaded(void* voidPtrArgs);
  static void* doSomeUpsampledCrossCorrelationsThreaded(void* voidPtrArgs);  

  /**********************************************************************************************************
  Calculate deltaT between two antennas (for a plane wave unless function name says otherwise)
  **********************************************************************************************************/
  Double_t getDeltaTExpected(Int_t ant1, Int_t ant2,Double_t phiWave, Double_t thetaWave);
  Int_t getDeltaTExpectedSpherical(Int_t ant1, Int_t ant2,Double_t phiWave, Double_t thetaWave, Double_t rWave);

  /**********************************************************************************************************
  Precalculate DeltaTs during initialization where appropriate
  **********************************************************************************************************/
  Bool_t useCombo(Int_t ant1, Int_t ant2, Int_t phiSector);
  void fillCombosToUseIfNeeded(mapMode_t mapMode, UInt_t l3TrigPattern);
  void do5PhiSectorCombinatorics();
  void fillDeltaTLookup();
  void fillDeltaTLookupZoomed(Double_t zoomCenterPhiDeg, Double_t zoomCenterThetaDeg, UInt_t l3TrigPattern);
  Double_t getBin0PhiDeg();
  // void writeDeltaTsFile(); // Speed up initialization
  // Int_t readDeltaTsFile(); // Speed up initialization

  /**********************************************************************************************************
  Image generation functions.
  **********************************************************************************************************/

  // Create blank histogram with proper axis ranges and axis titles
  void createImageNameAndTitle(TString& name, TString& title, mapMode_t mapMode, zoomMode_t zoomMode,
			       Double_t rWave, AnitaPol::AnitaPol_t pol);
  
  TH2D* makeBlankImage(TString name, TString title);
  TH2D* makeBlankZoomedImage(TString name, TString title, Double_t zoomCenterPhiDeg, Double_t zoomCenterThetaDeg);
  
  // Workhorse function which is called by many aliases.
  // Don't call directly unless you want to pass a lot of redundant parameters
  TH2D* makeImage(AnitaPol::AnitaPol_t pol, Double_t rWave, Double_t& imagePeak, 
		  Double_t& peakPhiDeg, Double_t& peakThetaDeg, UInt_t l3TrigPattern,
		  mapMode_t mapMode, zoomMode_t zoomMode,
		  Double_t zoomCenterPhiDeg=0, Double_t zoomCenterThetaDeg=0);


  TH2D* makeGlobalImage(AnitaPol::AnitaPol_t pol, Double_t& imagePeak,
			Double_t& peakPhiDeg, Double_t& peakThetaDeg);

  TH2D* makeTriggeredImage(AnitaPol::AnitaPol_t pol, Double_t& imagePeak, Double_t& peakPhiDeg,
			   Double_t& peakThetaDeg, UInt_t l3TrigPattern);
  

  TH2D* makeGlobalSphericalImage(AnitaPol::AnitaPol_t pol, Double_t rWave);
  TH2D* makeGlobalSphericalImage(AnitaPol::AnitaPol_t pol, Double_t rWave,
				 Double_t& imagePeak, Double_t& peakPhiDeg, Double_t& peakThetaDeg);

  TH2D* makeTriggeredSphericalImage(AnitaPol::AnitaPol_t pol, Double_t rWave, UInt_t l3TrigPattern);
  TH2D* makeTriggeredSphericalImage(AnitaPol::AnitaPol_t pol, Double_t rWave, Double_t& imagePeak,
				    Double_t& peakPhiDeg, Double_t& peakThetaDeg, UInt_t l3TrigPattern);
  
  Double_t findImagePeak(TH2D* hist, Double_t& imagePeakTheta, Double_t& imagePeakPhi);

  

  // To be completed
  TH2D* makeZoomedImage(AnitaPol::AnitaPol_t pol, UInt_t l3TrigPattern,
			Double_t zoomCenterPhiDeg, Double_t zoomCenterThetaDeg);

  TH2D* makeZoomedImage(AnitaPol::AnitaPol_t pol, Double_t& imagePeak, Double_t& peakPhiDeg,
			Double_t& peakThetaDeg, UInt_t l3TrigPattern, Double_t zoomCenterPhiDeg,
			Double_t zoomCenterThetaDeg);

  static void* makeSomeOfImageThreaded(void* voidPtrArgs);
  

  /**********************************************************************************************************
  Functions to delete pointers to internal variables
  **********************************************************************************************************/
  void deleteCrossCorrelations(AnitaPol::AnitaPol_t pol);
  void deleteUpsampledCrossCorrelations(AnitaPol::AnitaPol_t pol);
  void deleteAllWaveforms(AnitaPol::AnitaPol_t pol);
  void deleteAllFFTs(AnitaPol::AnitaPol_t pol);
  void deleteAllPaddedFFTs(AnitaPol::AnitaPol_t pol);


  /**********************************************************************************************************
  Functions for debugging or testing
  **********************************************************************************************************/
  void correlateEventTest(Double_t phiDegSource, Double_t thetaDegSource, Double_t rSource);
  TGraph* getCrossCorrelationGraph(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2);
  TGraph* getUpsampledCrossCorrelationGraph(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2);
  TGraph* getCrossCorrelationGraphWorker(Int_t numSamps, AnitaPol::AnitaPol_t pol,
					 Int_t ant1, Int_t ant2);
  
  
  /**********************************************************************************************************
  Variables
  **********************************************************************************************************/
  UInt_t eventNumber[NUM_POL]; ///< For tracking event number
  UInt_t lastEventNormalized[NUM_POL]; ///< Prevents cross-correlation of the same event twice
  Double_t nominalSamplingDeltaT; ///< ANITA-3 => 1./2.6 ns
  Double_t upsampleFactor; ///< Default = 1, 
  Double_t correlationDeltaT; ///< nominalSamplingDeltaT/upsampleFactor, deltaT of interpolation and cross correlation
  Int_t numSamples; ///< Number of samples in waveform after padding
  Int_t numSamplesUpsampled; ///< Number of samples in waveform after padding and up sampling
  Int_t numCombos;
  
  std::vector<Int_t> ant2s[NUM_SEAVEYS]; ///< Vector holding ant2 indices for each ant1 (for combinatorics)
  std::vector<Int_t> comboToAnt1s; ///< Vector mapping combination index to ant1
  std::vector<Int_t> comboToAnt2s; ///< Vector mapping combination index to ant2
  std::vector<Int_t> combosToUseGlobal[NUM_PHI]; ///< Depends on L3 trigger for global image
  std::map<UInt_t, std::vector<Int_t> > combosToUseTriggered; ///< Depends on L3 trigger for triggered image  
  Int_t comboIndices[NUM_SEAVEYS][NUM_SEAVEYS]; ///< Array mapping ant1+ant2 to combo index
  Double_t* crossCorrelations[NUM_POL][NUM_COMBOS]; ///< Arrays for cross correlations
  Double_t* crossCorrelationsUpsampled[NUM_POL][NUM_COMBOS]; ///< Arrays for upsampled cross correlations
  TGraph* grs[NUM_POL][NUM_SEAVEYS]; ///< Raw waveforms obtained from a UsefulAnitaEvent
  TGraph* grsInterp[NUM_POL][NUM_SEAVEYS]; ///< Interpolated TGraphs
  std::complex<Double_t>* ffts[NUM_POL][NUM_SEAVEYS]; ///< FFTs of TGraphs
  std::complex<Double_t>* fftsPadded[NUM_POL][NUM_SEAVEYS]; ///< Padded with zeros.
  Double_t interpRMS[NUM_POL][NUM_SEAVEYS]; ///< RMS of interpolation
  std::vector<Double_t> rArray; ///< Vector of antenna radial positions
  std::vector<Double_t> phiArrayDeg; ///< Vector of antenna azimuth positions
  std::vector<Double_t> zArray; ///< Vector of antenna heights

  typedef Short_t dtIndex_t;
  dtIndex_t deltaTs[NUM_PHI*NUM_BINS_PHI][NUM_BINS_THETA][NUM_COMBOS]; ///< Lookup of deltaTs between antenna pairs for making an image.
  dtIndex_t deltaTsZoom[NUM_BINS_PHI_ZOOM][NUM_BINS_THETA_ZOOM][NUM_COMBOS]; ///< Lookup of deltaTs between antenna pairs for making a zoomed in image, must recalculated each event (probably)

  // pairs for making an image
  Int_t deltaTMax;
  Int_t deltaTMin;

  AnitaPol::AnitaPol_t threadPol;
  UInt_t threadL3TrigPattern;
  std::vector<TThread*> mapThreads;
  std::vector<TThread*> corrThreads;
  std::vector<TThread*> upsampledCorrThreads;

  struct threadArgs{
    Long_t threadInd;
    CrossCorrelator* ptr;
  };

private:
  // Messing with this will muck up the threading.
  std::vector<threadArgs> threadArgsVec;  
  
  
  ClassDef(CrossCorrelator, 0);
};
#endif