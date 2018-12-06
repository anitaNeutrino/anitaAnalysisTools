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
  class CorrelationSummary;

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

    TGraph* getCrossCorrelationGraph(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2) const;
    TGraph* getUpsampledCrossCorrelationGraph(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2) const;
    TGraph* getCrossCorrelationGraphWorker(Int_t numSamps, AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2) const;


    virtual Double_t getCrossCorrelation(AnitaPol::AnitaPol_t pol, Int_t combo, Double_t deltaT) const;
    void initializeVariables();
    void getNormalizedInterpolatedTGraphs(const FilteredAnitaEvent* realEvent, AnitaPol::AnitaPol_t pol, bool raw = false);
    void renormalizeFourierDomain(AnitaPol::AnitaPol_t pol, Int_t ant);
    void doFFTs(AnitaPol::AnitaPol_t pol);
    void doCrossCorrelations(AnitaPol::AnitaPol_t pol);

  
    Bool_t useCombo(Int_t ant1, Int_t ant2, Int_t phiSector, Int_t deltaPhiSect);
    void fillCombosToUse();
    void do5PhiSectorCombinatorics();

    std::shared_ptr<const Acclaim::CorrelationSummary> makeSummary(AnitaPol::AnitaPol_t pol, const FilteredAnitaEvent* event, double waisPhi, double waisTheta, Adu5Pat* pat);

    double correlationIndexToTime(bool upsampled, int corrIndex, AnitaPol::AnitaPol_t pol, int combo) const;
    double correlationIndexToTime(bool upsampled, int corrIndex, AnitaPol::AnitaPol_t pol, int ant1, int ant2) const;

    //--------------------------------------------------------------------------------------------------------
    // Public member variables
    //--------------------------------------------------------------------------------------------------------

    static constexpr Int_t numSamples = NUM_SAMP*PAD_FACTOR; //!< Number of samples in waveform after padding.
    static constexpr Int_t numSamplesUpsampled = numSamples*UPSAMPLE_FACTOR; //!< Number of samples in waveform after padding and up sampling.
    static constexpr Double_t nominalSamplingDeltaT = NOMINAL_SAMPLING_DELTAT; //!< ANITA-3 => 1./2.6 ns, deltaT for evenly resampling.
    static constexpr Double_t correlationDeltaT = nominalSamplingDeltaT/UPSAMPLE_FACTOR;; //!< nominalSamplingDeltaT/UPSAMPLE_FACTOR, deltaT of for interpolation.
    
    std::complex<Double_t> ffts[AnitaPol::kNotAPol][NUM_SEAVEYS][GET_NUM_FREQS(numSamples)]; //!< FFTs of evenly resampled waveforms.
    Double_t crossCorrelations[AnitaPol::kNotAPol][NUM_COMBOS][numSamples]; //!< Cross correlations.
    Double_t crossCorrelationsUpsampled[AnitaPol::kNotAPol][NUM_COMBOS][numSamplesUpsampled]; //!< Upsampled cross correlations.
  
    Double_t fVolts[AnitaPol::kNotAPol][NUM_SEAVEYS][numSamples]; //!< Hold the filtered waveforms for padding...
    Double_t startTimes[AnitaPol::kNotAPol][NUM_SEAVEYS];
  
    Double_t interpRMS[AnitaPol::kNotAPol][NUM_SEAVEYS]; //!< RMS of interpolated TGraphs.
    Double_t interpRMS2[AnitaPol::kNotAPol][NUM_SEAVEYS]; //!< RMS of interpolated TGraphs with extra zero padding.
    Int_t comboIndices[NUM_SEAVEYS][NUM_SEAVEYS]; //!< Array mapping ant1+ant2 to combo index

    UInt_t eventNumber[AnitaPol::kNotAPol]; //!< For tracking event number
    Int_t numCombos; //!< Number of possible antenna pairs, counted during initialization. Should equal NUM_COMBOS.

    std::vector<Int_t> comboToAnt1s; //!< Vector mapping combo index to ant1.
    std::vector<Int_t> comboToAnt2s; //!< Vector mapping combo index to ant2.
    std::array<std::vector<Int_t>, NUM_PHI>combosToUseGlobal; //!< Depends on L3 trigger for global image


    Int_t multiplyTopRingByMinusOne; //!< For showing how I'm an idiot with respect to compiling the ANITA-3 prioritizer
    Int_t kOnlyThisCombo; //!< For debugging, only fill histograms with one particular antenna pair.
    Int_t kDeltaPhiSect; //!< Specifies how many phi-sectors around peak use in reconstruction.

  private:

    //--------------------------------------------------------------------------------------------------------
    // Private member variables
    //--------------------------------------------------------------------------------------------------------


  };

}

#endif
