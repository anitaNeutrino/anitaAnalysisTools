/* -*- C++ -*-.*********************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description: 
             A Cross Correlator to interact with the ROOTified ANITA-3 data and do some interferometry. 
***********************************************************************************************************/

#ifndef CROSSCORRELATOR_H
#define CROSSCORRELATOR_H

// Ryan things
#include "UsefulAnitaEvent.h"
#include "AnitaEventCalibrator.h"
#include "AnitaGeomTool.h"
#include "UsefulAdu5Pat.h"

// #include "FFTtools.h"
#include "FancyFFTs.h"

// My things
#include "RootTools.h"

// ROOT things
#include "TGraph.h"
#include "TH2D.h"
#include "TThread.h"

// standard c++ things
#include <iostream>


// Number of phi-sectors to cross correlate between
#define DELTA_PHI_SECT 2

// Offline reconstruction definitions
#define NUM_COMBOS 336
#define NUM_THREADS 4

// Typical number of samples in waveform
#define NUM_SAMPLES 256
#define UPSAMPLE_FACTOR 40
#define NOMINAL_SAMPLING_DELTAT (1./2.6f)
#define PAD_FACTOR 2
#define GETNUMFREQS(n)((n)/2+1)

// Image definitions
#define NUM_BINS_THETA 100
#define NUM_BINS_PHI 15
#define THETA_RANGE 150
#define PHI_RANGE 22.5

#define NUM_BINS_THETA_ZOOM 100
#define NUM_BINS_PHI_ZOOM 100
// #define NUM_BINS_THETA_ZOOM 60
// #define NUM_BINS_PHI_ZOOM 60
#define ZOOM_BIN_SIZE_PHI 0.05
#define ZOOM_BIN_SIZE_THETA 0.05
#define THETA_RANGE_ZOOM (NUM_BINS_THETA_ZOOM*ZOOM_BIN_SIZE_THETA)
#define PHI_RANGE_ZOOM (NUM_BINS_PHI_ZOOM*ZOOM_BIN_SIZE_PHI)

// #define NUM_BINS_PHI_ZOOM_TOTAL 7200
// #define NUM_BINS_THETA_ZOOM_TOTAL 3000
#define NUM_BINS_PHI_ZOOM_TOTAL 7269
#define NUM_BINS_THETA_ZOOM_TOTAL 3069
// #define NUM_BINS_PHI_ZOOM_TOTAL 7229
// #define NUM_BINS_THETA_ZOOM_TOTAL 3029

// Anita Geometry definitions, shouldn't really be here
#define NUM_POL AnitaPol::kNotAPol
#define NUM_RING AnitaRing::kNotARing
#define DEGREES_IN_CIRCLE 360

#define MAX_NUM_PEAKS 2
#define PEAK_PHI_DEG_RANGE 10
#define PEAK_THETA_DEG_RANGE 10

#define SPEED_OF_LIGHT 2.99792458e8
#define SPEED_OF_LIGHT_NS 0.299792458

#define ALL_PHI_TRIGS 0xffff



/**
 * @class CrossCorrelator
 * @brief A class to take in UsefulAnitaEvents and get interferometric maps with a single function.
 * 
 * Does all the heavy lifting: gets waveforms from a UsefulAnitaEvent, cross correlates them, and produces interferometric maps. 
*/
class CrossCorrelator{

public:

    /**
   * @brief Flag to pass to CrossCorrelator when making a map telling it whether to use all phi-sectors or triggered phi-sectors.
   */  
  enum mapMode_t{
    kGlobal,
    kTriggered,
    kNumMapModes
  };

  /**
   * @brief Flag to pass to CrossCorrelator when making a map telling it whether to reconstruct all arrival directions or a finer binned close up of a particular region
   */
  enum zoomMode_t{
    kZoomedOut,
    kZoomedIn,
    kNumZoomModes
  };

  /**
   * @brief Container required to get threading to work inside a class, includes the thread index and the pointer to the class.
   *
   * All functions called by threads have to be static.
   * So we write any functions which we want to be threaded as static member functions, which take a pointer to the class as the first argument (requires that pointer).
   * We then use the thread index to figure out what portion of the work each thread should do.
   */
  struct threadArgs{
    Long_t threadInd; //!< The thread index
    CrossCorrelator* ptr; //!< Pointer to the CrossCorrelator
  };

  CrossCorrelator();
  ~CrossCorrelator();

  void initializeVariables();
  void printInfo();

  void getNormalizedInterpolatedTGraphs(UsefulAnitaEvent* realEvent, AnitaPol::AnitaPol_t pol);
  TGraph* interpolateWithStartTimeAndZeroMean(TGraph* grIn, Double_t startTime, Double_t dt, Int_t nSamp);
  void doFFTs(AnitaPol::AnitaPol_t pol);
  void correlateEvent(UsefulAnitaEvent* realEvent);
  void correlateEvent(UsefulAnitaEvent* realEvent, AnitaPol::AnitaPol_t pol);
  void reconstructEvent(UsefulAnitaEvent* usefulEvent, Int_t numFinePeaks=MAX_NUM_PEAKS, Int_t numCoarsePeaks=MAX_NUM_PEAKS);
  void findPeakValues(AnitaPol::AnitaPol_t pol, Int_t numPeaks, Double_t* peakValues,
		      Double_t* phiDegs, Double_t* thetaDegs);


  void getMaxCorrelationTimeValue(AnitaPol::AnitaPol_t pol, Int_t combo,Double_t& time, Double_t& value);
  void getMaxCorrelationTimeValue(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2, Double_t& time, Double_t& value);
  void getMaxUpsampledCorrelationTimeValue(AnitaPol::AnitaPol_t pol, Int_t combo, Double_t& time, Double_t& value);
  void getMaxUpsampledCorrelationTimeValue(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2, Double_t& time, Double_t& value);


  
  void doAllCrossCorrelationsThreaded(AnitaPol::AnitaPol_t pol);
  void doUpsampledCrossCorrelationsThreaded(AnitaPol::AnitaPol_t pol, Int_t phiSector);



  
  Double_t getDeltaTExpected(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2,
			     Double_t phiWave, Double_t thetaWave);
  // Double_t getDeltaTExpectedFast(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2,
  // 					Int_t phiIndex, Int_t thetaIndex);
  Double_t getOffAxisDelay(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2, Double_t phiWave);



  Bool_t useCombo(Int_t ant1, Int_t ant2, Int_t phiSector, Int_t deltaPhiSect);
  void fillCombosToUse();
  void do5PhiSectorCombinatorics();


  void fillDeltaTLookup();
  Double_t getBin0PhiDeg();








  void getCoarsePeakInfo(AnitaPol::AnitaPol_t pol, Int_t peakIndex, Double_t& value,
			 Double_t& phiDeg, Double_t& thetaDeg);
  void getFinePeakInfo(AnitaPol::AnitaPol_t pol, Int_t peakIndex, Double_t& value,
		       Double_t& phiDeg, Double_t& thetaDeg);  


  
  Double_t getInterpolatedUpsampledCorrelationValue(AnitaPol::AnitaPol_t pol, Int_t combo, Double_t deltaT);

  TH2D* getMap(AnitaPol::AnitaPol_t pol, Double_t& peakValue,
	       Double_t& peakPhiDeg, Double_t& peakThetaDeg,
	       UShort_t l3TrigPattern=ALL_PHI_TRIGS);
    

  TH2D* getZoomMap(AnitaPol::AnitaPol_t pol);
  
  void reconstruct(AnitaPol::AnitaPol_t pol, Double_t& imagePeak,
		   Double_t& peakPhiDeg, Double_t& peakThetaDeg);
  void reconstructZoom(AnitaPol::AnitaPol_t pol,
		       Double_t& imagePeak, Double_t& peakPhiDeg,
		       Double_t& peakThetaDeg, Double_t zoomCenterPhiDeg=0, Double_t zoomCenterThetaDeg=0);
  
  TH2D* makeGlobalImage(AnitaPol::AnitaPol_t pol, Double_t& imagePeak,
  			Double_t& peakPhiDeg, Double_t& peakThetaDeg);
  TH2D* makeGlobalImage(AnitaPol::AnitaPol_t pol);

  TH2D* makeTriggeredImage(AnitaPol::AnitaPol_t pol, Double_t& imagePeak, Double_t& peakPhiDeg,
  			   Double_t& peakThetaDeg, UShort_t l3TrigPattern);

  TH2D* makeZoomedImage(AnitaPol::AnitaPol_t pol, UShort_t l3TrigPattern,
  			Double_t zoomCenterPhiDeg, Double_t zoomCenterThetaDeg);
  TH2D* makeZoomedImage(AnitaPol::AnitaPol_t pol, Double_t& imagePeak, Double_t& peakPhiDeg,
  			Double_t& peakThetaDeg, UShort_t l3TrigPattern, Double_t zoomCenterPhiDeg,
  			Double_t zoomCenterThetaDeg);
  TH2D* makeZoomedImage(AnitaPol::AnitaPol_t pol, Double_t& imagePeak, Double_t& peakPhiDeg,
  			Double_t& peakThetaDeg, Double_t zoomCenterPhiDeg, Double_t zoomCenterThetaDeg);

  Int_t getPhiSectorOfAntennaClosestToPhiDeg(AnitaPol::AnitaPol_t pol, Double_t phiDeg);
  TGraph* makeCoherentlySummedWaveform(AnitaPol::AnitaPol_t pol, Double_t phiDeg,
				       Double_t thetaDeg, Int_t maxDeltaPhiSect, Double_t& snr);
  TGraph* makeUpsampledCoherentlySummedWaveform(AnitaPol::AnitaPol_t pol, Double_t phiDeg,
						Double_t thetaDeg, Int_t maxDeltaPhiSect, Double_t& snr);  
  TGraph* makeCoherentWorker(AnitaPol::AnitaPol_t pol, Double_t phiDeg, Double_t thetaDeg,
			     Int_t maxDeltaPhiSect, Double_t& snr, Int_t nSamp);
  static void* doSomeCrossCorrelationsThreaded(void* voidPtrArgs);
  static void* doSomeUpsampledCrossCorrelationsThreaded(void* voidPtrArgs);
  static void* makeSomeOfImageThreaded(void* voidPtrArgs);
  static void* makeSomeOfZoomImageThreaded(void* voidPtrArgs);
  
  void deleteAllWaveforms(AnitaPol::AnitaPol_t pol);

  TGraph* getCrossCorrelationGraph(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2);
  TGraph* getUpsampledCrossCorrelationGraph(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2);
  TGraph* getCrossCorrelationGraphWorker(Int_t numSamps, AnitaPol::AnitaPol_t pol,
					 Int_t ant1, Int_t ant2);
  static Int_t directlyInsertGeometry(TString pathToLindasFile, AnitaPol::AnitaPol_t pol);
  void insertPhotogrammetryGeometry();
  
  TString mapModeNames[kNumMapModes];//!< Maps text to the mapMode_t enum, used for histogram names/titles.
  TString zoomModeNames[kNumZoomModes];//!< Maps text to the zoomMode_t enum, used for histogram names/titles.
  
  // Double_t interpPreFactors[NUM_POL][NUM_COMBOS][NUM_BINS_THETA][NUM_PHI*NUM_BINS_PHI]; //!< The interpolation factor for neighbouring samples
  Double_t interpPreFactors[NUM_POL][NUM_COMBOS][NUM_PHI*NUM_BINS_PHI][NUM_BINS_THETA]; //!< The interpolation factor for neighbouring samples  
  // Int_t offsetLows[NUM_POL][NUM_COMBOS][NUM_BINS_THETA][NUM_PHI*NUM_BINS_PHI]; //!< The interpolation factor for neighbouring samples
  Int_t offsetLows[NUM_POL][NUM_COMBOS][NUM_PHI*NUM_BINS_PHI][NUM_BINS_THETA]; //!< The interpolation factor for neighbouring samples  

  Double_t crossCorrelations[NUM_POL][NUM_COMBOS][NUM_SAMPLES*PAD_FACTOR]; //!< Cross correlations.
  Double_t coarseMap[NUM_POL][NUM_BINS_PHI*NUM_PHI][NUM_BINS_THETA]; //!< Internal storage for the coarsely binned map
  
  Double_t partBAsZoom[NUM_POL][NUM_COMBOS][NUM_BINS_THETA_ZOOM_TOTAL]; //!< Yet more geometric caching
  Double_t part21sZoom[NUM_POL][NUM_COMBOS][NUM_BINS_PHI_ZOOM_TOTAL]; //!< Yet more geometric caching

  Double_t crossCorrelationsUpsampled[NUM_POL][NUM_COMBOS][NUM_SAMPLES*PAD_FACTOR*UPSAMPLE_FACTOR*PAD_FACTOR]; //!< Upsampled cross correlations.
  Double_t fineMap[NUM_POL][NUM_BINS_THETA_ZOOM][NUM_BINS_PHI_ZOOM]; //!< Internal storage for the finely binned map
  
  std::complex<Double_t> fftsPadded[NUM_POL][NUM_SEAVEYS][GETNUMFREQS(NUM_SAMPLES*PAD_FACTOR*UPSAMPLE_FACTOR)]; //!< FFTs of evenly resampled waveforms, padded with zeros so that the inverse fourier transform is interpolated.

  std::complex<Double_t> ffts[NUM_POL][NUM_SEAVEYS][GETNUMFREQS(NUM_SAMPLES*PAD_FACTOR)]; //!< FFTs of evenly resampled waveforms.
  TGraph* grs[NUM_POL][NUM_SEAVEYS]; //!< Raw waveforms obtained from the UsefulAnitaEvent.
  TGraph* grsResampled[NUM_POL][NUM_SEAVEYS]; //!< Evenly resampled TGraphs.
  Double_t interpRMS[NUM_POL][NUM_SEAVEYS]; //!< RMS of interpolated TGraphs.
  Double_t interpRMS2[NUM_POL][NUM_SEAVEYS]; //!< RMS of interpolated TGraphs with extra zero padding.  
  Int_t comboIndices[NUM_SEAVEYS][NUM_SEAVEYS]; //!< Array mapping ant1+ant2 to combo index

  UInt_t eventNumber[NUM_POL]; //!< For tracking event number
  UInt_t lastEventNormalized[NUM_POL]; //!< Prevents cross-correlation of the same event twice
  Double_t nominalSamplingDeltaT; //!< ANITA-3 => 1./2.6 ns, deltaT for evenly resampling.
  Double_t correlationDeltaT; //!< nominalSamplingDeltaT/UPSAMPLE_FACTOR, deltaT of for interpolation.
  Int_t numSamples; //!< Number of samples in waveform after padding.
  Int_t numSamplesUpsampled; //!< Number of samples in waveform after padding and up sampling.
  Int_t numCombos; //!< Number of possible antenna pairs, counted during initialization. Should equal NUM_COMBOS.
  
  std::vector<Int_t> comboToAnt1s; //!< Vector mapping combo index to ant1.
  std::vector<Int_t> comboToAnt2s; //!< Vector mapping combo index to ant2.
  std::vector<Int_t> combosToUseGlobal[NUM_PHI]; //!< Depends on L3 trigger for global image

  std::map<std::pair<UInt_t, Int_t>, std::vector<Int_t> > combosToUseTriggered; //!< Depends on L3 trigger for triggered image

  std::vector<Double_t> rArray[NUM_POL]; //!< Vector of antenna radial positions
  std::vector<Double_t> phiArrayDeg[NUM_POL]; //!< Vector of antenna azimuth positions
  std::vector<Double_t> zArray[NUM_POL]; //!< Vector of antenna heights

  Double_t coarseMapPeakValues[NUM_POL][MAX_NUM_PEAKS]; //!< Stores the peak of the interally stored map
  Double_t coarseMapPeakPhiDegs[NUM_POL][MAX_NUM_PEAKS]; //!< Stores the peak phi (degrees) of the interally stored map
  Double_t coarseMapPeakThetaDegs[NUM_POL][MAX_NUM_PEAKS]; //!< Stores the peak theta (degrees) of the interally stored map

  Double_t fineMapPeakValues[NUM_POL][MAX_NUM_PEAKS]; //!< Stores the peak of the interally stored map
  Double_t fineMapPeakPhiDegs[NUM_POL][MAX_NUM_PEAKS]; //!< Stores the peak phi (degrees) of the interally stored map
  Double_t fineMapPeakThetaDegs[NUM_POL][MAX_NUM_PEAKS]; //!< Stores the peak theta (degrees) of the interally stored map

  
  Double_t zoomedThetaWaves[NUM_BINS_THETA_ZOOM_TOTAL]; //!< Cached theta for zoomed image.
  Double_t zoomedTanThetaWaves[NUM_BINS_THETA_ZOOM_TOTAL]; //!< Cached tan(theta) for zoomed image.
  Double_t zoomedCosThetaWaves[NUM_BINS_THETA_ZOOM_TOTAL]; //!< Cached cos(theta) for zoomed image.
  Double_t zoomedPhiWaveLookup[NUM_BINS_PHI_ZOOM_TOTAL]; //!< Cached phi for zoomed image.  
  Double_t zoomedCosPartLookup[NUM_POL][NUM_SEAVEYS][NUM_BINS_PHI_ZOOM_TOTAL]; //!< Cached part of the deltaT calculation.
  Double_t offAxisDelays[NUM_POL][NUM_COMBOS][NUM_BINS_PHI_ZOOM_TOTAL]; //!< Off-axis delays for fine binned images.
  
  Double_t thetaWaves[NUM_BINS_THETA]; //!< Cached theta for image.
  Double_t phiWaveLookup[NUM_BINS_PHI*NUM_PHI]; //!< Cached phi for image.
  
  AnitaPol::AnitaPol_t threadPol; //!< Polarization to use in thread functions.
  UInt_t threadL3TrigPattern; //!< l3TrigPattern to use in thread functions.
  Int_t threadPhiSector; //!< phi-sector to use in thread functions.  
  Double_t threadImagePeak[NUM_THREADS]; //!< Store image peaks found by different threads.
  Double_t threadPeakPhiDeg[NUM_THREADS]; //!< Store phi of image peaks found by different threads.
  Double_t threadPeakThetaDeg[NUM_THREADS]; //!< Store theta of image peaks found by different threads.
  
  std::vector<TThread*> mapThreads; //!< TThreads for doing interferometric map making.
  std::vector<TThread*> corrThreads; //!< TThreads for doing cross correlations.
  std::vector<TThread*> upsampledCorrThreads; //!< TThreads for doing upsampled cross correlations.
   
  Int_t kOnlyThisCombo; //!< For debugging, only fill histograms with one particular antenna pair.
  Int_t kDeltaPhiSect; //!< Specifies how many phi-sectors around peak use in reconstruction.
  Int_t kUseOffAxisDelay; //!< Flag for whether or not to apply off axis delay to deltaT expected.
  Double_t maxDPhiDeg; //!< Variable for testing how wide an off axis angle is used in reconstruction

private:
  
  std::vector<threadArgs> threadArgsVec; //!< Vector of threadArgs, accessed by threaded functions so they can work out what portion of the work are supposed to be doing.

  Double_t aftForeOffset; //!< From AnitaGeomTool, defines the location of the antennas relative to the axis of the heading.
  Double_t minThetaDegZoom; //!< Emperically determined minimum possible zoomed theta (Degrees)
  Double_t minPhiDegZoom; //!< Emperically determined minimum possible zoomed phi (Degrees)
  Double_t zoomPhiMin; //!< For the current map
  Double_t zoomThetaMin; //!< For the current map

  
};
#endif
