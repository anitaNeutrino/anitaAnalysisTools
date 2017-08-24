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
// #define PEAK_THETA_DEG_RANGE 180

#define SPEED_OF_LIGHT 2.99792458e8
#define SPEED_OF_LIGHT_NS 0.299792458


// Number of phi-sectors to cross correlate between
#define DELTA_PHI_SECT 2

// Offline reconstruction definitions
#define NUM_COMBOS 336

// Typical number of samples in waveform
#define UPSAMPLE_FACTOR 6
#define NOMINAL_SAMPLING_DELTAT (1./2.6f)
#define PAD_FACTOR 2
#define GET_NUM_FREQS(n)((n)/2+1)

namespace Acclaim
{

  /**
   * @class CrossCorrelator
   * @brief A class to take in FiteredAnitaEvents and cross-correlate nearby channels
   *
   */
  class CrossCorrelator{

  public:

    //--------------------------------------------------------------------------------------------------------
    // Public member functions
    //--------------------------------------------------------------------------------------------------------

    CrossCorrelator();
    virtual ~CrossCorrelator();

    virtual void correlateEvent(const FilteredAnitaEvent* realEvent);
    virtual void correlateEvent(const FilteredAnitaEvent* realEvent, AnitaPol::AnitaPol_t pol);

    void doUpsampledCrossCorrelations(AnitaPol::AnitaPol_t pol, Int_t phiSector);  

    TGraph* getCrossCorrelationGraph(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2);
    TGraph* getUpsampledCrossCorrelationGraph(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2);
    TGraph* getCrossCorrelationGraphWorker(Int_t numSamps, AnitaPol::AnitaPol_t pol,
					   Int_t ant1, Int_t ant2);


    virtual Double_t getCrossCorrelation(AnitaPol::AnitaPol_t pol, Int_t combo, Double_t deltaT) const;

    void initializeVariables();
    void getNormalizedInterpolatedTGraphs(const FilteredAnitaEvent* realEvent, AnitaPol::AnitaPol_t pol, bool raw = false);
    void renormalizeFourierDomain(AnitaPol::AnitaPol_t pol, Int_t ant);
    void doFFTs(AnitaPol::AnitaPol_t pol);
    void doCrossCorrelations(AnitaPol::AnitaPol_t pol);

  
    Bool_t useCombo(Int_t ant1, Int_t ant2, Int_t phiSector, Int_t deltaPhiSect);
    void fillCombosToUse();
    void do5PhiSectorCombinatorics();

    Double_t getInterpolatedUpsampledCorrelationValue(AnitaPol::AnitaPol_t pol, Int_t combo, Double_t deltaT);




    //--------------------------------------------------------------------------------------------------------
    // Public member variables
    //--------------------------------------------------------------------------------------------------------

    std::complex<Double_t> ffts[AnitaPol::kNotAPol][NUM_SEAVEYS][GET_NUM_FREQS(NUM_SAMP*PAD_FACTOR)]; //!< FFTs of evenly resampled waveforms.
    Double_t crossCorrelations[AnitaPol::kNotAPol][NUM_COMBOS][NUM_SAMP*PAD_FACTOR]; //!< Cross correlations.
    Double_t crossCorrelationsUpsampled[AnitaPol::kNotAPol][NUM_COMBOS][NUM_SAMP*PAD_FACTOR*UPSAMPLE_FACTOR*PAD_FACTOR]; //!< Upsampled cross correlations.
  
    Double_t fVolts[AnitaPol::kNotAPol][NUM_SEAVEYS][NUM_SAMP*PAD_FACTOR]; //!< Hold the filtered waveforms for padding...
    Double_t startTimes[AnitaPol::kNotAPol][NUM_SEAVEYS];
  
    Double_t interpRMS[AnitaPol::kNotAPol][NUM_SEAVEYS]; //!< RMS of interpolated TGraphs.
    Double_t interpRMS2[AnitaPol::kNotAPol][NUM_SEAVEYS]; //!< RMS of interpolated TGraphs with extra zero padding.
    Int_t comboIndices[NUM_SEAVEYS][NUM_SEAVEYS]; //!< Array mapping ant1+ant2 to combo index

    UInt_t eventNumber[AnitaPol::kNotAPol]; //!< For tracking event number
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



  /**
   * @class TemplateCorrelator
   * @brief Cross-correlate all channels of an event with a template
   * 
   * This class overloads a few member functions while repurposing some member variable
   * so not really a case of straight forward inheritance...
   */
class TemplateCorrelator : public CrossCorrelator {

 public:
  TemplateCorrelator(Int_t run, UInt_t eventNumber);
  virtual ~TemplateCorrelator();
  void initTemplate(const FilteredAnitaEvent* fEv);
  void initTemplate(Int_t run, UInt_t eventNumber);
  virtual void correlateEvent(const FilteredAnitaEvent* fEv);
  virtual void correlateEvent(const FilteredAnitaEvent* fEv, AnitaPol::AnitaPol_t pol);

  TGraph* getCrossCorrelationGraph(AnitaPol::AnitaPol_t pol, Int_t ant) const;
  
  Double_t getPeakCorrelation(AnitaPol::AnitaPol_t pol, Double_t minOffset=-100, Double_t maxOffset=100, Double_t stepSize=NOMINAL_SAMPLING_DELTAT) const;
  
  virtual Double_t getCrossCorrelation(AnitaPol::AnitaPol_t pol, Int_t ant, Double_t deltaT) const;

 protected:
  // double templateAllChannelRMS;
  
};

}

#endif
