#ifndef INTERFEROMETRIC_MAP_MAKER_H
#define INTERFEROMETRIC_MAP_MAKER_H

#include "AnitaEventReconstructor.h"
#include "InterferometricMap.h"
#include "CrossCorrelator.h"

#define DEGREES_IN_CIRCLE 360

#define MAX_NUM_PEAKS 5
#define PEAK_PHI_DEG_RANGE 10
#define PEAK_THETA_DEG_RANGE 10

// Image definitions
// #define NUM_BINS_THETA 100
// #define NUM_BINS_PHI 15
// #define THETA_RANGE 150
// #define PHI_RANGE 22.5


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

#define ZOOM_BIN_SIZE_PHI (1./ZOOM_BINS_PER_DEGREE_PHI)
#define ZOOM_BIN_SIZE_THETA (1./ZOOM_BINS_PER_DEGREE_THETA)
#define THETA_RANGE_ZOOM (NUM_BINS_THETA_ZOOM*ZOOM_BIN_SIZE_THETA)
#define PHI_RANGE_ZOOM (NUM_BINS_PHI_ZOOM*ZOOM_BIN_SIZE_PHI)

#define NUM_BINS_PHI_ZOOM_TOTAL (DEGREES_IN_CIRCLE*ZOOM_BINS_PER_DEGREE_PHI + NUM_BINS_PHI_ZOOM)
#define NUM_BINS_THETA_ZOOM_TOTAL ((MAX_THETA - MIN_THETA)*ZOOM_BINS_PER_DEGREE_THETA + NUM_BINS_THETA_ZOOM)

#define ALL_PHI_TRIGS 0xffff




class DeltaTCache {

  Double_t* operator[](int phiBin){return &vals[NUM_BINS_THETA*phiBin];}
  
private:
  // Double_t vals[NUM_BINS_PHI*NUM_PHI][NUM_BINS_THETA];
  Double_t vals[NUM_BINS_PHI*NUM_PHI*NUM_BINS_THETA];    
  
  
};


/** 
 * Class to make interferometric maps from FilteredAnitaEvents...
 */
class InterferometricMapMaker : public AnitaEventReconstructor{


public:

    /**
   * @brief Flag to pass to InterferometricMapMaker when making a map telling it whether to use all phi-sectors or triggered phi-sectors.
   */
  enum mapMode_t{
    kGlobal,
    kTriggered,
    kNumMapModes
  };

  /**
   * @brief Flag to pass to InterferometricMapMaker when making a map telling it whether to reconstruct all arrival directions or a finer binned close up of a particular region
   */
  enum zoomMode_t{
    kZoomedOut,
    kZoomedIn,
    kNumZoomModes
  };

  
  virtual void process(const FilteredAnitaEvent * ev, UsefulAdu5Pat* pat ,AnitaEventSummary * summary) const;

  InterferometricMapMaker();
  virtual ~InterferometricMapMaker();

  void initializeInternals();


  void reconstructEvent(FilteredAnitaEvent* usefulEvent, Int_t numFinePeaks=MAX_NUM_PEAKS, Int_t numCoarsePeaks=MAX_NUM_PEAKS);

  AnitaPol::AnitaPol_t reconstructEventPeakPol(FilteredAnitaEvent* usefulEvent, Int_t numFinePeaks=MAX_NUM_PEAKS, Int_t numCoarsePeaks=MAX_NUM_PEAKS);

  void reconstructEvent(FilteredAnitaEvent* usefulEvent, UsefulAdu5Pat& usefulPat, AnitaEventSummary* eventSummary);

  void findPeakValues(AnitaPol::AnitaPol_t pol, Int_t numPeaks, Double_t* peakValues,
		      Double_t* phiDegs, Double_t* thetaDegs);

  Double_t getDeltaTExpected(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2,
			     Double_t phiWave, Double_t thetaWave);
  // Double_t getDeltaTExpectedFast(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2,
  // 					Int_t phiIndex, Int_t thetaIndex);

  Double_t relativeOffAxisDelay(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2, Double_t phiDeg);
  Double_t singleAntennaOffAxisDelay(Double_t deltaPhiDeg);


  void fillDeltaTLookup();
  Double_t getBin0PhiDeg();


  TGraph* makeCoherentlySummedWaveform(AnitaPol::AnitaPol_t pol, Double_t phiDeg,
				       Double_t thetaDeg, Int_t maxDeltaPhiSect, Double_t& snr);
  TGraph* makeUpsampledCoherentlySummedWaveform(AnitaPol::AnitaPol_t pol, Double_t phiDeg,
						Double_t thetaDeg, Int_t maxDeltaPhiSect, Double_t& snr);
  TGraph* makeCoherentWorker(AnitaPol::AnitaPol_t pol, Double_t phiDeg, Double_t thetaDeg,
			     Int_t maxDeltaPhiSect, Double_t& snr, Int_t nSamp);






  void getCoarsePeakInfo(AnitaPol::AnitaPol_t pol, Int_t peakIndex, Double_t& value,
			 Double_t& phiDeg, Double_t& thetaDeg);
  void getFinePeakInfo(AnitaPol::AnitaPol_t pol, Int_t peakIndex, Double_t& value,
		       Double_t& phiDeg, Double_t& thetaDeg);


  InterferometricMap* getInterferometricMap(AnitaPol::AnitaPol_t pol){
    InterferometricMap* h = coarseMaps[pol];
    coarseMaps[pol] = NULL;
    return h;
  }
  
  InterferometricMap* getMap(AnitaPol::AnitaPol_t pol, Double_t& peakValue,
			     Double_t& peakPhiDeg, Double_t& peakThetaDeg,
			     UShort_t l3TrigPattern=ALL_PHI_TRIGS);


  TH2D* getZoomMap(AnitaPol::AnitaPol_t pol, Int_t peakInd=0);

  void reconstruct(AnitaPol::AnitaPol_t pol, Double_t& imagePeak,
		   Double_t& peakPhiDeg, Double_t& peakThetaDeg);
  void reconstructZoom(AnitaPol::AnitaPol_t pol,
		       Double_t& imagePeak, Double_t& peakPhiDeg,
		       Double_t& peakThetaDeg, Double_t zoomCenterPhiDeg=0, Double_t zoomCenterThetaDeg=0,
		       Int_t peakIndex = 0);

  InterferometricMap* makeGlobalImage(AnitaPol::AnitaPol_t pol, Double_t& imagePeak,
				      Double_t& peakPhiDeg, Double_t& peakThetaDeg);
  InterferometricMap* makeGlobalImage(AnitaPol::AnitaPol_t pol);

  InterferometricMap* makeTriggeredImage(AnitaPol::AnitaPol_t pol, Double_t& imagePeak, Double_t& peakPhiDeg,
					 Double_t& peakThetaDeg, UShort_t l3TrigPattern);

  TH2D* makeZoomedImage(AnitaPol::AnitaPol_t pol, UShort_t l3TrigPattern,
  			Double_t zoomCenterPhiDeg, Double_t zoomCenterThetaDeg);
  TH2D* makeZoomedImage(AnitaPol::AnitaPol_t pol, Double_t& imagePeak, Double_t& peakPhiDeg,
  			Double_t& peakThetaDeg, UShort_t l3TrigPattern, Double_t zoomCenterPhiDeg,
  			Double_t zoomCenterThetaDeg);
  TH2D* makeZoomedImage(AnitaPol::AnitaPol_t pol, Double_t& imagePeak, Double_t& peakPhiDeg,
  			Double_t& peakThetaDeg, Double_t zoomCenterPhiDeg, Double_t zoomCenterThetaDeg);

  Int_t getPhiSectorOfAntennaClosestToPhiDeg(AnitaPol::AnitaPol_t pol, Double_t phiDeg);

  static Int_t directlyInsertGeometry(TString pathToLindasFile, AnitaPol::AnitaPol_t pol);
  void insertPhotogrammetryGeometry();

  UInt_t eventNumber[AnitaPol::kNotAPol];

  Double_t interpPreFactors[AnitaPol::kNotAPol][NUM_COMBOS][NUM_PHI*NUM_BINS_PHI][NUM_BINS_THETA]; //!< The interpolation factor for neighbouring samples
  Int_t offsetLows[AnitaPol::kNotAPol][NUM_COMBOS][NUM_PHI*NUM_BINS_PHI][NUM_BINS_THETA]; //!< The interpolation factor for neighbouring samples
  // std::vector<DeltaTCache> deltaTs[AnitaPol::kNotAPol];

  Double_t coarseMap[AnitaPol::kNotAPol][NUM_BINS_PHI*NUM_PHI][NUM_BINS_THETA]; //!< Internal storage for the coarsely binned map
  Double_t partBAsZoom[AnitaPol::kNotAPol][NUM_COMBOS][NUM_BINS_THETA_ZOOM_TOTAL]; //!< Yet more geometric caching
  Double_t part21sZoom[AnitaPol::kNotAPol][NUM_COMBOS][NUM_BINS_PHI_ZOOM_TOTAL]; //!< Yet more geometric caching

  Double_t fineMap[AnitaPol::kNotAPol][MAX_NUM_PEAKS][NUM_BINS_THETA_ZOOM][NUM_BINS_PHI_ZOOM]; //!< Internal storage for the finely binned map
  std::vector<Double_t> rArray[AnitaPol::kNotAPol]; //!< Vector of antenna radial positions
  std::vector<Double_t> phiArrayDeg[AnitaPol::kNotAPol]; //!< Vector of antenna azimuth positions
  std::vector<Double_t> zArray[AnitaPol::kNotAPol]; //!< Vector of antenna heights

  Double_t coarseMapPeakValues[AnitaPol::kNotAPol][MAX_NUM_PEAKS]; //!< Stores the peak of the interally stored map
  Double_t coarseMapPeakPhiDegs[AnitaPol::kNotAPol][MAX_NUM_PEAKS]; //!< Stores the peak phi (degrees) of the interally stored map
  Double_t coarseMapPeakThetaDegs[AnitaPol::kNotAPol][MAX_NUM_PEAKS]; //!< Stores the peak theta (degrees) of the interally stored map

  Double_t fineMapPeakValues[AnitaPol::kNotAPol][MAX_NUM_PEAKS]; //!< Stores the peak of the interally stored map
  Double_t fineMapPeakPhiDegs[AnitaPol::kNotAPol][MAX_NUM_PEAKS]; //!< Stores the peak phi (degrees) of the interally stored map
  Double_t fineMapPeakThetaDegs[AnitaPol::kNotAPol][MAX_NUM_PEAKS]; //!< Stores the peak theta (degrees) of the interally stored map

  Double_t zoomedThetaWaves[NUM_BINS_THETA_ZOOM_TOTAL]; //!< Cached theta for zoomed image.
  Double_t zoomedTanThetaWaves[NUM_BINS_THETA_ZOOM_TOTAL]; //!< Cached tan(theta) for zoomed image.
  Double_t zoomedCosThetaWaves[NUM_BINS_THETA_ZOOM_TOTAL]; //!< Cached cos(theta) for zoomed image.
  Double_t dtFactors[NUM_BINS_THETA_ZOOM_TOTAL]; //!< Cached cos(theta)/c/dt for zoomed image.
  Double_t zoomedPhiWaveLookup[NUM_BINS_PHI_ZOOM_TOTAL]; //!< Cached phi for zoomed image.
  Double_t zoomedCosPartLookup[AnitaPol::kNotAPol][NUM_SEAVEYS][NUM_BINS_PHI_ZOOM_TOTAL]; //!< Cached part of the deltaT calculation.
  Double_t offAxisDelays[AnitaPol::kNotAPol][NUM_COMBOS][NUM_BINS_PHI_ZOOM_TOTAL]; //!< Off-axis delays for fine binned images.
  Double_t offAxisDelaysDivided[AnitaPol::kNotAPol][NUM_COMBOS][NUM_BINS_PHI_ZOOM_TOTAL]; //!< Off-axis delays divided such to remove an operation from the inner loop of an image making function.

  Double_t thetaWaves[NUM_BINS_THETA]; //!< Cached theta for image.
  Double_t phiWaveLookup[NUM_BINS_PHI*NUM_PHI]; //!< Cached phi for image.

  Int_t kUseOffAxisDelay; //!< Flag for whether or not to apply off axis delay to deltaT expected.
  Double_t maxDPhiDeg; //!< Variable for testing how wide an off axis angle is used in reconstruction
  Int_t coherentDeltaPhi;


  Double_t minThetaDegZoom; //!< Minimum possible zoomed theta (Degrees)
  Double_t minPhiDegZoom; //!< Minimum possible zoomed phi (Degrees)
  Double_t zoomPhiMin[AnitaPol::kNotAPol]; //!< For the current map
  Double_t zoomThetaMin[AnitaPol::kNotAPol]; //!< For the current map

  TString mapModeNames[kNumMapModes];//!< Maps text to the mapMode_t enum, used for histogram names/titles.
  TString zoomModeNames[kNumZoomModes];//!< Maps text to the zoomMode_t enum, used for histogram names/titles.
  
  
  CrossCorrelator* cc;
  
private:
  void initializeVariables();

  // The axes of the interferometric maps store the uneven bins in theta (degrees)
  // these will replace the internal map storage...
  // we will imply some pointer ownership scheme with these.
  // I delete them if the event number changes and they've not been returned. You delete them otherwise...  

  InterferometricMap* coarseMaps[AnitaPol::kNotAPol]; // these guys do the whole 360 az, and defined elevation...
  
};




#endif
