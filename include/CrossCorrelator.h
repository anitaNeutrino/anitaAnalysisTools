/* -*- C++ -*-.*********************************************************************************************
 Author: Ben Strutt
 Email: strutt@physics.ucla.edu

 Description:
             A Cross Correlator to interact with the ROOTified ANITA-3 data and do some interferometry.
***********************************************************************************************************/

#ifndef CROSSCORRELATOR_H
#define CROSSCORRELATOR_H

// Ryan things
// #include "UsefulAnitaEvent.h"
#include "FilteredAnitaEvent.h"
#include "AnalysisWaveform.h"
#include "AnitaEventCalibrator.h"
#include "AnitaGeomTool.h"
#include "UsefulAdu5Pat.h"
#include "FFTtools.h"

#include "AnitaEventSummary.h"

// My things
#include "RootTools.h"
#include "FancyFFTs.h"

// ROOT things
#include "TGraph.h"
#include "TH2D.h"
#include "TROOT.h" // for gDirectory pointer?

// standard c++ things
#include <iostream>



// Anita & Geometry definitions
#define NUM_POL AnitaPol::kNotAPol
#define NUM_RING AnitaRing::kNotARing
#define DEGREES_IN_CIRCLE 360

#define MAX_NUM_PEAKS 5
#define PEAK_PHI_DEG_RANGE 10
#define PEAK_THETA_DEG_RANGE 10
// #define PEAK_THETA_DEG_RANGE 180

#define SPEED_OF_LIGHT 2.99792458e8
#define SPEED_OF_LIGHT_NS 0.299792458


// Number of phi-sectors to cross correlate between
#define DELTA_PHI_SECT 2

// Offline reconstruction definitions
#define NUM_COMBOS 336

// Typical number of samples in waveform
#define NUM_SAMPLES 260
#define UPSAMPLE_FACTOR 6
#define NOMINAL_SAMPLING_DELTAT (1./2.6f)
#define PAD_FACTOR 2
#define GET_NUM_FREQS(n)((n)/2+1)



/**
 * @class CrossCorrelator
 * @brief A class to take in UsefulAnitaEvents or FiteredAnitaEvents and get interferometric maps with a single function.
 *
 */
class CrossCorrelator{

public:



  //--------------------------------------------------------------------------------------------------------
  // Public member functions
  //--------------------------------------------------------------------------------------------------------



  CrossCorrelator();
  ~CrossCorrelator();

  void initializeVariables();
  void printInfo();

  void getNormalizedInterpolatedTGraphs(FilteredAnitaEvent* realEvent, AnitaPol::AnitaPol_t pol);

  void renormalizeFourierDomain(AnitaPol::AnitaPol_t pol, Int_t ant);

  TGraph* interpolateWithStartTimeAndZeroMean(TGraph* grIn, Double_t startTime, Double_t dt, Int_t nSamp);
  void doFFTs(AnitaPol::AnitaPol_t pol);

  // template <class FilteredAnitaEvent> // needs eventNumber member and getGraph(int ant, AnitaPol::AnitaPol_t pol)
  void correlateEvent(FilteredAnitaEvent* realEvent);

  // template <class FilteredAnitaEvent> // needs eventNumber member and getGraph(int ant, AnitaPol::AnitaPol_t pol)
  void correlateEvent(FilteredAnitaEvent* realEvent, AnitaPol::AnitaPol_t pol);

  void getMaxCorrelationTimeValue(AnitaPol::AnitaPol_t pol, Int_t combo,Double_t& time, Double_t& value);
  void getMaxCorrelationTimeValue(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2, Double_t& time, Double_t& value);
  void getMaxUpsampledCorrelationTimeValue(AnitaPol::AnitaPol_t pol, Int_t combo, Double_t& time, Double_t& value);
  void getMaxUpsampledCorrelationTimeValue(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2, Double_t& time, Double_t& value);

  void doAllCrossCorrelations(AnitaPol::AnitaPol_t pol);
  void doUpsampledCrossCorrelations(AnitaPol::AnitaPol_t pol, Int_t phiSector);

  Bool_t useCombo(Int_t ant1, Int_t ant2, Int_t phiSector, Int_t deltaPhiSect);
  void fillCombosToUse();
  void do5PhiSectorCombinatorics();



  Double_t getInterpolatedUpsampledCorrelationValue(AnitaPol::AnitaPol_t pol, Int_t combo, Double_t deltaT);


  TGraph* getCrossCorrelationGraph(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2);
  TGraph* getUpsampledCrossCorrelationGraph(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2);
  TGraph* getCrossCorrelationGraphWorker(Int_t numSamps, AnitaPol::AnitaPol_t pol,
					 Int_t ant1, Int_t ant2);



  Double_t getCrossCorrelation(AnitaPol::AnitaPol_t pol, Int_t combo, Double_t deltaT);


  //--------------------------------------------------------------------------------------------------------
  // Public member variables
  //--------------------------------------------------------------------------------------------------------

  Double_t crossCorrelations[NUM_POL][NUM_COMBOS][NUM_SAMPLES*PAD_FACTOR]; //!< Cross correlations.

  Double_t crossCorrelationsUpsampled[NUM_POL][NUM_COMBOS][NUM_SAMPLES*PAD_FACTOR*UPSAMPLE_FACTOR*PAD_FACTOR]; //!< Upsampled cross correlations.

  std::complex<Double_t> ffts[NUM_POL][NUM_SEAVEYS][GET_NUM_FREQS(NUM_SAMPLES*PAD_FACTOR)]; //!< FFTs of evenly resampled waveforms.


  Double_t fVolts[NUM_POL][NUM_SEAVEYS][NUM_SAMPLES*PAD_FACTOR]; //!< Hold the filtered waveforms for padding...
  Double_t startTimes[NUM_POL][NUM_SEAVEYS];
  
  Double_t interpRMS[NUM_POL][NUM_SEAVEYS]; //!< RMS of interpolated TGraphs.
  Double_t interpRMS2[NUM_POL][NUM_SEAVEYS]; //!< RMS of interpolated TGraphs with extra zero padding.
  Int_t comboIndices[NUM_SEAVEYS][NUM_SEAVEYS]; //!< Array mapping ant1+ant2 to combo index

  UInt_t eventNumber[NUM_POL]; //!< For tracking event number
  Double_t nominalSamplingDeltaT; //!< ANITA-3 => 1./2.6 ns, deltaT for evenly resampling.
  Double_t correlationDeltaT; //!< nominalSamplingDeltaT/UPSAMPLE_FACTOR, deltaT of for interpolation.
  Int_t numSamples; //!< Number of samples in waveform after padding.
  Int_t numSamplesUpsampled; //!< Number of samples in waveform after padding and up sampling.
  Int_t numCombos; //!< Number of possible antenna pairs, counted during initialization. Should equal NUM_COMBOS.

  std::vector<Int_t> comboToAnt1s; //!< Vector mapping combo index to ant1.
  std::vector<Int_t> comboToAnt2s; //!< Vector mapping combo index to ant2.
  std::vector<Int_t> combosToUseGlobal[NUM_PHI]; //!< Depends on L3 trigger for global image


  Int_t multiplyTopRingByMinusOne; //!< For showing how I'm an idiot with respect to compiling the ANITA-3 prioritizer
  Int_t kOnlyThisCombo; //!< For debugging, only fill histograms with one particular antenna pair.
  Int_t kDeltaPhiSect; //!< Specifies how many phi-sectors around peak use in reconstruction.

private:

  //--------------------------------------------------------------------------------------------------------
  // Private member variables
  //--------------------------------------------------------------------------------------------------------


};
#endif
