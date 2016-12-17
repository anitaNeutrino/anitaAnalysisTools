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
#define MAX_THREADS 8

// Typical number of samples in waveform
#define NUM_SAMPLES 256
#define UPSAMPLE_FACTOR 6
// #define UPSAMPLE_FACTOR 10
// #define UPSAMPLE_FACTOR 40
#define NOMINAL_SAMPLING_DELTAT (1./2.6f)
#define PAD_FACTOR 2
#define GET_NUM_FREQS(n)((n)/2+1)

// Image definitions
// #define NUM_BINS_THETA 100
// #define NUM_BINS_PHI 15
// #define THETA_RANGE 150
// #define PHI_RANGE 22.5

// #define NUM_BINS_THETA 60
#define NUM_BINS_THETA 35
#define NUM_BINS_PHI 9
// #define THETA_RANGE 150
#define MIN_THETA -55
#define MAX_THETA 35
#define PHI_RANGE 22.5


// #define NUM_BINS_THETA_ZOOM 200
// #define NUM_BINS_PHI_ZOOM 200
// #define ZOOM_BINS_PER_DEGREE_PHI 20
// #define ZOOM_BINS_PER_DEGREE_THETA 20

// #define NUM_BINS_THETA_ZOOM 100
// #define NUM_BINS_PHI_ZOOM 200
// #define ZOOM_BINS_PER_DEGREE_PHI 20
// #define ZOOM_BINS_PER_DEGREE_THETA 20

#define NUM_BINS_THETA_ZOOM 40
#define NUM_BINS_PHI_ZOOM 40
#define ZOOM_BINS_PER_DEGREE_PHI 4
#define ZOOM_BINS_PER_DEGREE_THETA 4

// angular distance to traverse akima fits along peak
#define AKIMA_PEAK_INTERP_DELTA_THETA_DEG 0.05
#define AKIMA_PEAK_INTERP_DELTA_PHI_DEG 0.05

#define ZOOM_BIN_SIZE_PHI (1./ZOOM_BINS_PER_DEGREE_PHI)
#define ZOOM_BIN_SIZE_THETA (1./ZOOM_BINS_PER_DEGREE_THETA)
#define THETA_RANGE_ZOOM (NUM_BINS_THETA_ZOOM*ZOOM_BIN_SIZE_THETA)
#define PHI_RANGE_ZOOM (NUM_BINS_PHI_ZOOM*ZOOM_BIN_SIZE_PHI)

#define NUM_BINS_PHI_ZOOM_TOTAL (DEGREES_IN_CIRCLE*ZOOM_BINS_PER_DEGREE_PHI + NUM_BINS_PHI_ZOOM)
#define NUM_BINS_THETA_ZOOM_TOTAL ((MAX_THETA - MIN_THETA)*ZOOM_BINS_PER_DEGREE_THETA + NUM_BINS_THETA_ZOOM)

#define ALL_PHI_TRIGS 0xffff



/**
 * @class CrossCorrelator
 * @brief A class to take in UsefulAnitaEvents and get interferometric maps with a single function.
 *
 * Does all the heavy lifting: gets waveforms from a UsefulAnitaEvent, cross correlates them, and produces interferometric maps.
*/
class CrossCorrelator{

public:






  //--------------------------------------------------------------------------------------------------------
  // Classes declared inside this class
  //--------------------------------------------------------------------------------------------------------




  //--------------------------------------------------------------------------------------------------------
  /**
   * @class SimpleNotch
   * @ A class to hold two frequencies, a low notch edge and high notch edge. Should be ROOT read/writable.
   */
  class SimpleNotch : public TNamed{
  public:
    //------------------------------------------------------------------------------------------------------
    /**
     * @brief Default Constructor for ROOT IO
     *
     */
    SimpleNotch(){
      lowPassFreqMHz=0;
      highPassFreqMHz=0;
    }
    //------------------------------------------------------------------------------------------------------
    /**
     * @brief Proper Constructor
     *
     * @param name The name of the notch, for ROOT IO
     * @param title The title of the notch, for ROOT IO
     * @param theLowPassFreqMHz The low pass frequency in MHz, i.e. the low edge of the notch (frequencies less than this value ARE NOT filtered, and frequencies greater than and equal to this value ARE filtered)
     * @param theHighPassFreqMHz The high pass frequency in MHz, i.e. the high edge of the notch (frequencies less than this value ARE filtered, and frequencies greater than and equal to this value ARE NOT filtered)
     */
    SimpleNotch(TString name, TString title, Double_t theLowPassFreqMHz, Double_t theHighPassFreqMHz) : lowPassFreqMHz(theLowPassFreqMHz) , highPassFreqMHz(theHighPassFreqMHz){
      fName = name;
      fTitle = title;
      if(lowPassFreqMHz > highPassFreqMHz){
	std::cerr << "Warning in " << __FUNCTION__ << ", your highPassFreqMHz < lowPassFreqMHz!" << std::endl;
	printInfo(std::cerr);
	std::cerr << "This notch isn't going to do anything!" << std::endl;
      }
    }
    //------------------------------------------------------------------------------------------------------
    /**
     * @grief Get the notch edge values, retured by reference
     *
     * @param theLowPassFreqMHz gets the theLowPassFreqMHz
     * @param theHighPassFreqMHz get the highPassFreqMHz
     */

    void getNotchEdges(Double_t &theLowPassFreqMHz, Double_t& theHighPassFreqMHz) const{
      theLowPassFreqMHz =  lowPassFreqMHz;
      theHighPassFreqMHz = highPassFreqMHz;
    }
    void printInfo(std::ostream& output = std::cout){
      output << fName << "\t" << fTitle << ": Low Pass = " << lowPassFreqMHz
	     << " MHz, High Pass = " << highPassFreqMHz << " MHz" << std::endl;
    }

  private:
    Double_t lowPassFreqMHz;
    Double_t highPassFreqMHz;
    ClassDef(SimpleNotch, 1)
  };


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














  //--------------------------------------------------------------------------------------------------------
  // Public member functions
  //--------------------------------------------------------------------------------------------------------



  CrossCorrelator();
  ~CrossCorrelator();

  void initializeVariables();
  void printInfo();

  void getNormalizedInterpolatedTGraphs(UsefulAnitaEvent* realEvent, AnitaPol::AnitaPol_t pol);
  Double_t applyNotch(AnitaPol::AnitaPol_t pol, Int_t ant, const SimpleNotch& notch);
  Double_t applyFilterBits(AnitaPol::AnitaPol_t pol, Int_t ant, TBits* filterBits);

  void writeNotchesIfAble();
  void renormalizeFourierDomain(AnitaPol::AnitaPol_t pol, Int_t ant);

  TGraph* interpolateWithStartTimeAndZeroMean(TGraph* grIn, Double_t startTime, Double_t dt, Int_t nSamp);
  void doFFTs(AnitaPol::AnitaPol_t pol, TBits* filterBits=NULL);
  void correlateEvent(UsefulAnitaEvent* realEvent, TBits* filterBits=NULL);
  void correlateEvent(UsefulAnitaEvent* realEvent, AnitaPol::AnitaPol_t pol, TBits* filterBits=NULL);
  void reconstructEvent(UsefulAnitaEvent* usefulEvent, Int_t numFinePeaks=MAX_NUM_PEAKS, Int_t numCoarsePeaks=MAX_NUM_PEAKS, TBits* filterBits = NULL);
  AnitaPol::AnitaPol_t reconstructEventPeakPol(UsefulAnitaEvent* usefulEvent, Int_t numFinePeaks=MAX_NUM_PEAKS, Int_t numCoarsePeaks=MAX_NUM_PEAKS, TBits* filterBits = NULL);
  void findPeakValues(AnitaPol::AnitaPol_t pol, Int_t numPeaks, Double_t* peakValues,
		      Double_t* phiDegs, Double_t* thetaDegs);


  void getMaxCorrelationTimeValue(AnitaPol::AnitaPol_t pol, Int_t combo,Double_t& time, Double_t& value);
  void getMaxCorrelationTimeValue(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2, Double_t& time, Double_t& value);
  void getMaxUpsampledCorrelationTimeValue(AnitaPol::AnitaPol_t pol, Int_t combo, Double_t& time, Double_t& value);
  void getMaxUpsampledCorrelationTimeValue(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2, Double_t& time, Double_t& value);



  void doAllCrossCorrelationsThreaded(AnitaPol::AnitaPol_t pol);
  void akimaUpsampleCrossCorrelations(AnitaPol::AnitaPol_t pol, Int_t phiSector); // akima interp
  void doUpsampledCrossCorrelationsThreaded(AnitaPol::AnitaPol_t pol, Int_t phiSector);


  Double_t getDeltaTExpected(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2,
			     Double_t phiWave, Double_t thetaWave);
  // Double_t getDeltaTExpectedFast(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2,
  // 					Int_t phiIndex, Int_t thetaIndex);

  Double_t relativeOffAxisDelay(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2, Double_t phiDeg);
  Double_t singleAntennaOffAxisDelay(Double_t deltaPhiDeg);


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


  TH2D* getZoomMap(AnitaPol::AnitaPol_t pol, Int_t peakInd=0);

  void reconstruct(AnitaPol::AnitaPol_t pol, Double_t& imagePeak,
		   Double_t& peakPhiDeg, Double_t& peakThetaDeg);
  void reconstructZoom(AnitaPol::AnitaPol_t pol,
		       Double_t& imagePeak, Double_t& peakPhiDeg,
		       Double_t& peakThetaDeg, Double_t zoomCenterPhiDeg=0, Double_t zoomCenterThetaDeg=0,
		       Int_t peakIndex = 0);

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


  UInt_t addNotch(SimpleNotch simpleNotch);
  UInt_t setDefaultNotches();
  void printNotchInfo();











  //--------------------------------------------------------------------------------------------------------
  // Public member variables
  //--------------------------------------------------------------------------------------------------------


  TString mapModeNames[kNumMapModes];//!< Maps text to the mapMode_t enum, used for histogram names/titles.
  TString zoomModeNames[kNumZoomModes];//!< Maps text to the zoomMode_t enum, used for histogram names/titles.

  // Double_t interpPreFactors[NUM_POL][NUM_COMBOS][NUM_BINS_THETA][NUM_PHI*NUM_BINS_PHI]; //!< The interpolation factor for neighbouring samples
  Double_t interpPreFactors[NUM_POL][NUM_COMBOS][NUM_PHI*NUM_BINS_PHI][NUM_BINS_THETA]; //!< The interpolation factor for neighbouring samples
  // Double_t interpPreFactors[NUM_POL][NUM_PHI*NUM_BINS_PHI][NUM_COMBOS][NUM_BINS_THETA]; //!< The interpolation factor for neighbouring samples
  // Int_t offsetLows[NUM_POL][NUM_COMBOS][NUM_BINS_THETA][NUM_PHI*NUM_BINS_PHI]; //!< The interpolation factor for neighbouring samples
  Int_t offsetLows[NUM_POL][NUM_COMBOS][NUM_PHI*NUM_BINS_PHI][NUM_BINS_THETA]; //!< The interpolation factor for neighbouring samples
  // Int_t offsetLows[NUM_POL][NUM_PHI*NUM_BINS_PHI][NUM_COMBOS][NUM_BINS_THETA]; //!< The interpolation factor for neighbouring samples

  Double_t crossCorrelations[NUM_POL][NUM_COMBOS][NUM_SAMPLES*PAD_FACTOR]; //!< Cross correlations.
  Double_t coarseMap[NUM_POL][NUM_BINS_PHI*NUM_PHI][NUM_BINS_THETA]; //!< Internal storage for the coarsely binned map

  Double_t partBAsZoom[NUM_POL][NUM_COMBOS][NUM_BINS_THETA_ZOOM_TOTAL]; //!< Yet more geometric caching
  Double_t part21sZoom[NUM_POL][NUM_COMBOS][NUM_BINS_PHI_ZOOM_TOTAL]; //!< Yet more geometric caching

  Double_t crossCorrelationsUpsampled[NUM_POL][NUM_COMBOS][NUM_SAMPLES*PAD_FACTOR*UPSAMPLE_FACTOR*PAD_FACTOR]; //!< Upsampled cross correlations.
  Double_t fineMap[NUM_POL][MAX_NUM_PEAKS][NUM_BINS_THETA_ZOOM][NUM_BINS_PHI_ZOOM]; //!< Internal storage for the finely binned map

  std::complex<Double_t> ffts[NUM_POL][NUM_SEAVEYS][GET_NUM_FREQS(NUM_SAMPLES*PAD_FACTOR)]; //!< FFTs of evenly resampled waveforms.
  TGraph* grs[NUM_POL][NUM_SEAVEYS]; //!< Raw waveforms obtained from the UsefulAnitaEvent.
  TGraph* grsResampled[NUM_POL][NUM_SEAVEYS]; //!< Evenly resampled TGraphs.
  Double_t interpRMS[NUM_POL][NUM_SEAVEYS]; //!< RMS of interpolated TGraphs.
  Double_t interpRMS2[NUM_POL][NUM_SEAVEYS]; //!< RMS of interpolated TGraphs with extra zero padding.
  Int_t comboIndices[NUM_SEAVEYS][NUM_SEAVEYS]; //!< Array mapping ant1+ant2 to combo index

  UInt_t eventNumber[NUM_POL]; //!< For tracking event number
  UInt_t lastEventNormalized[NUM_POL]; //!< Prevents cross-correlation of the same event twice
  UInt_t lastEventUpsampleCorrelated[NUM_POL][NUM_COMBOS]; //!< Prevents upsampled cross-correlation of the same event twice
  Double_t nominalSamplingDeltaT; //!< ANITA-3 => 1./2.6 ns, deltaT for evenly resampling.
  Double_t correlationDeltaT; //!< nominalSamplingDeltaT/UPSAMPLE_FACTOR, deltaT of for interpolation.
  Int_t numSamples; //!< Number of samples in waveform after padding.
  Int_t numSamplesUpsampled; //!< Number of samples in waveform after padding and up sampling.
  Int_t numCombos; //!< Number of possible antenna pairs, counted during initialization. Should equal NUM_COMBOS.

  std::vector<Int_t> comboToAnt1s; //!< Vector mapping combo index to ant1.
  std::vector<Int_t> comboToAnt2s; //!< Vector mapping combo index to ant2.
  std::vector<Int_t> combosToUseGlobal[NUM_PHI]; //!< Depends on L3 trigger for global image

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
  Double_t dtFactors[NUM_BINS_THETA_ZOOM_TOTAL]; //!< Cached cos(theta)/c/dt for zoomed image.
  Double_t zoomedPhiWaveLookup[NUM_BINS_PHI_ZOOM_TOTAL]; //!< Cached phi for zoomed image.
  Double_t zoomedCosPartLookup[NUM_POL][NUM_SEAVEYS][NUM_BINS_PHI_ZOOM_TOTAL]; //!< Cached part of the deltaT calculation.
  Double_t offAxisDelays[NUM_POL][NUM_COMBOS][NUM_BINS_PHI_ZOOM_TOTAL]; //!< Off-axis delays for fine binned images.
  Double_t offAxisDelaysDivided[NUM_POL][NUM_COMBOS][NUM_BINS_PHI_ZOOM_TOTAL]; //!< Off-axis delays divided such to remove an operation from the inner loop of an image making function.

  Double_t thetaWaves[NUM_BINS_THETA]; //!< Cached theta for image.
  Double_t phiWaveLookup[NUM_BINS_PHI*NUM_PHI]; //!< Cached phi for image.


  Int_t multiplyTopRingByMinusOne; //!< For showing how I'm an idiot with respect to compiling the ANITA-3 prioritizer
  Int_t kOnlyThisCombo; //!< For debugging, only fill histograms with one particular antenna pair.
  Int_t kDeltaPhiSect; //!< Specifies how many phi-sectors around peak use in reconstruction.
  Int_t kUseOffAxisDelay; //!< Flag for whether or not to apply off axis delay to deltaT expected.
  Double_t maxDPhiDeg; //!< Variable for testing how wide an off axis angle is used in reconstruction


  Double_t aftForeOffset; //!< From AnitaGeomTool, defines the location of the antennas relative to the axis of the heading.
  Double_t minThetaDegZoom; //!< Minimum possible zoomed theta (Degrees)
  Double_t minPhiDegZoom; //!< Minimum possible zoomed phi (Degrees)
  Double_t zoomPhiMin[NUM_POL]; //!< For the current map
  Double_t zoomThetaMin[NUM_POL]; //!< For the current map

  Int_t getNumThreads();
  Int_t setNumThreads(Int_t numDesiredThreads);

private:

  //--------------------------------------------------------------------------------------------------------
  // Private member variables
  //--------------------------------------------------------------------------------------------------------

  Int_t numThreads;//!< Number of threads to use in threaded functions.

  std::vector<SimpleNotch> allChannelNotches; //!< Holds notches to be applied to all channels for all events (e.g. satellite filters).

  ROOT::Math::Interpolator* corrInterp;//(tVec,corVec,ROOT::Math::Interpolation::kAKIMA_PERIODIC);
  std::vector<double> ccTimes;
  double ccMaxTime;


  AnitaPol::AnitaPol_t threadPol; //!< Polarization to use in thread functions.
  Int_t threadPeakIndex; //!< Which fine peak are we on
  UInt_t threadL3TrigPattern; //!< l3TrigPattern to use in thread functions.
  Int_t threadPhiSector; //!< phi-sector to use in thread functions.
  std::vector<Double_t> threadImagePeak; //!< Store image peaks found by different threads.
  std::vector<Double_t> threadPeakPhiDeg; //!< Store phi of image peaks found by different threads.
  std::vector<Double_t> threadPeakThetaDeg; //!< Store theta of image peaks found by different threads.
  std::vector<Double_t> threadImagePeakZoom; //!< Store image peaks found by different threads.
  std::vector<Double_t> threadPeakPhiDegZoom; //!< Store phi of image peaks found by different threads.
  std::vector<Double_t> threadPeakThetaDegZoom; //!< Store theta of image peaks found by different threads.
  std::vector<Int_t> threadPeakPhiBinZoom; //!< Store phi bin of image peaks found by different threads.
  std::vector<Int_t> threadPeakThetaBinZoom; //!< Store theta bin of image peaks found by different threads.



};
#endif
