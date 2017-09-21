#ifndef ANALYSIS_RECO_H
#define ANALYSIS_RECO_H

#include "AnitaEventReconstructor.h"
#include "InterferometricMap.h"
#include "InterferometryCache.h"
#include "CrossCorrelator.h"
#include "AnalysisSettings.h"
#include <list>

namespace AnitaResponse {
  class DeconvolutionMethod;
  class ResponseManager;
}
class NoiseMonitor;
class TruthAnitaEvent;

namespace Acclaim
{

  class InterferometryCache;


/** 
 * Class to make interferometric maps from FilteredAnitaEvents...
 */

class AnalysisReco : public TObject {

 public:

  // Constructor and destructor
  AnalysisReco();
  virtual ~AnalysisReco();  

  // Does the reconstruction.
  void process(const FilteredAnitaEvent * ev, AnitaEventSummary * summary, NoiseMonitor* noiseMonitor=NULL, TruthAnitaEvent* truth = NULL);

  // Master functions for delays between antennas.
  Double_t getDeltaTExpected(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2,
                             Double_t phiWave, Double_t thetaWave) const;
  Double_t relativeOffAxisDelay(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2, Double_t phiDeg) const;
  Double_t singleAntennaOffAxisDelay(Double_t deltaPhiDeg) const;

  // For MagicDisplay
  void drawSummary(TPad* pad, AnitaPol::AnitaPol_t pol);

  // Get the maps generated during the reconstruction 
  InterferometricMap* getMap(AnitaPol::AnitaPol_t pol);
  InterferometricMap* getZoomMap(AnitaPol::AnitaPol_t pol, UInt_t peakInd=0);
  TProfile2D* getHeatMap(AnitaPol::AnitaPol_t pol){return heatMaps[pol];}
  
  void reconstruct(AnitaPol::AnitaPol_t pol, const Adu5Pat* pat = NULL);
  void reconstructZoom(AnitaPol::AnitaPol_t pol, Int_t peakIndex, Double_t zoomCenterPhiDeg, Double_t zoomCenterThetaDeg, const Adu5Pat* pat = NULL);

  AnalysisWaveform* coherentlySum(const FilteredAnitaEvent* fEv, AnitaPol::AnitaPol_t pol, const std::vector<Int_t>& theAnts, Double_t peakPhiDeg, Double_t peakThetaDeg, Double_t* biggestPeakToPeak=NULL);
  AnalysisWaveform* coherentlySum(std::vector<const AnalysisWaveform*>& waves, std::vector<Double_t>& dts);


  // For calibration work
  static Int_t directlyInsertGeometry(TString pathToLindasFile, AnitaPol::AnitaPol_t pol);
  void insertPhotogrammetryGeometry();

  // Internal copies of the antenna positions from AnitaGeomTool
  std::vector<Double_t> rArray[AnitaPol::kNotAPol]; //!< Vector of antenna radial positions
  std::vector<Double_t> phiArrayDeg[AnitaPol::kNotAPol]; //!< Vector of antenna azimuth positions
  std::vector<Double_t> zArray[AnitaPol::kNotAPol]; //!< Vector of antenna heights
  
  AnalysisWaveform* getCoherentFiltered(AnitaPol::AnitaPol_t pol, Int_t peakInd=0, bool xPol=false);
  AnalysisWaveform* getCoherent(AnitaPol::AnitaPol_t pol, Int_t peakInd=0, bool xPol=false);
  AnalysisWaveform* getDeconvolved(AnitaPol::AnitaPol_t pol, Int_t peakInd=0, bool xPol=false);
  AnalysisWaveform* getDeconvolvedFiltered(AnitaPol::AnitaPol_t pol, Int_t peakInd=0, bool xPol=false);

  FilteredAnitaEvent* getEvMin();
  FilteredAnitaEvent* getEvMinDeco();
  FilteredAnitaEvent* getEvDeco();


  void directionAndAntennasToDeltaTs(const std::vector<Int_t>& theAnts, AnitaPol::AnitaPol_t pol,
                                     Double_t peakPhiDeg, Double_t peakThetaDeg, std::vector<double>& dts);

  static void setTriggerInfoFromPeakPhi(const RawAnitaHeader* header, AnitaPol::AnitaPol_t pol,
                                        Int_t peakPhiSector, AnitaEventSummary::PointingHypothesis& peak);
  static void fillChannelInfo(const FilteredAnitaEvent* fEv, AnitaEventSummary* sum);

  const std::vector<int>& phiSectorToCoherentAnts(int peakPhiSector){return phiSectorToAnts[peakPhiSector];}
  void wavesInCoherent(std::vector<const AnalysisWaveform*>& waves, std::vector<Double_t>& dts, std::vector<TGraphAligned*>& grs);

  CrossCorrelator* getCrossCorrelator(){return fCrossCorr;}
 protected:

  void initializeInternals();

  InterferometricMap* coarseMaps[AnitaPol::kNotAPol];
  InterferometricMap* fineMaps[AnitaPol::kNotAPol][AnitaEventSummary::maxDirectionsPerPol];
  TProfile2D* heatMaps[AnitaPol::kNotAPol];
  TProfile2D* makeHeatMap(const TString& name, const TString& title);

  size_t checkWavesAndDtsMatch(std::vector<const AnalysisWaveform*>& waves, std::vector<Double_t>& dts);
  void fillWaveformInfo(AnitaPol::AnitaPol_t pol, AnitaEventSummary::WaveformInfo& info, const FilteredAnitaEvent* fEv,
                        AnalysisWaveform** waveStore, InterferometricMap* h, NoiseMonitor* noiseMonitor);
  
  AnalysisWaveform* wfCoherentFiltered[AnitaPol::kNotAPol][AnitaEventSummary::maxDirectionsPerPol][AnitaPol::kNotAPol];
  AnalysisWaveform* wfCoherent[AnitaPol::kNotAPol][AnitaEventSummary::maxDirectionsPerPol][AnitaPol::kNotAPol];
  AnalysisWaveform* wfDeconvolved[AnitaPol::kNotAPol][AnitaEventSummary::maxDirectionsPerPol][AnitaPol::kNotAPol];
  AnalysisWaveform* wfDeconvolvedFiltered[AnitaPol::kNotAPol][AnitaEventSummary::maxDirectionsPerPol][AnitaPol::kNotAPol];
  std::vector<Int_t> phiSectorToAnts[NUM_PHI]; //!< Vector of antennas to use in coherently summed waveforms
  
  AnitaEventSummary summary;
  CrossCorrelator* fCrossCorr;
  bool spawnedCrossCorrelator;
  InterferometryCache dtCache;
  
  FilterStrategy* fMinFilter; //!< Would be no filters, but for ANITA-3 data we need an ALFA filter
  FilterStrategy* fMinDecoFilter; //!< Minimum filter + deconvolution filter

  FilteredAnitaEvent* fEvMin;
  FilteredAnitaEvent* fEvMinDeco;
  FilteredAnitaEvent* fEvDeco;

  void chooseAntennasForCoherentlySumming(int coherentDeltaPhi);
  void nicelyDeleteInternalFilteredEvents();
  
  ANALYSIS_SETTING(Int_t, CoherentDeltaPhi);
  Int_t fLastCoherentDeltaPhi; // for checking whether recalculation is needed
  ANALYSIS_SETTING(Int_t, WhichResponseDir);
  ANALYSIS_SETTING(Int_t, UseOffAxisDelay);
  ANALYSIS_SETTING(Int_t, ResponseNPad);
  ANALYSIS_SETTING(Int_t, NumPeaks);
  ANALYSIS_SETTING(Int_t, DoHeatMap);
  ANALYSIS_SETTING(Double_t, HeatMapHorizonKm);
  ANALYSIS_SETTING(Double_t, CoherentDtNs);
  ClassDef(AnalysisReco, 0)

};

}


#endif
