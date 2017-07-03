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
  void process(const FilteredAnitaEvent * ev, UsefulAdu5Pat* pat ,AnitaEventSummary * summary, NoiseMonitor* noiseMonitor=NULL);
  void process(const FilteredAnitaEvent * ev, Adu5Pat* pat, AnitaEventSummary * summary, NoiseMonitor* noiseMonitor=NULL);    

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
  
  void reconstruct(AnitaPol::AnitaPol_t pol);
  void reconstructZoom(AnitaPol::AnitaPol_t pol, Int_t peakIndex, Double_t zoomCenterPhiDeg, Double_t zoomCenterThetaDeg);

  AnalysisWaveform* coherentlySum(const FilteredAnitaEvent* fEv, AnitaPol::AnitaPol_t pol, const std::vector<Int_t>& theAnts, Double_t peakPhiDeg, Double_t peakThetaDeg);
  AnalysisWaveform* coherentlySum(std::vector<const AnalysisWaveform*>& waves, std::vector<Double_t>& dts);


  // For calibration work
  static Int_t directlyInsertGeometry(TString pathToLindasFile, AnitaPol::AnitaPol_t pol);
  void insertPhotogrammetryGeometry();

  // Internal copies of the antenna positions from AnitaGeomTool
  std::vector<Double_t> rArray[AnitaPol::kNotAPol]; //!< Vector of antenna radial positions
  std::vector<Double_t> phiArrayDeg[AnitaPol::kNotAPol]; //!< Vector of antenna azimuth positions
  std::vector<Double_t> zArray[AnitaPol::kNotAPol]; //!< Vector of antenna heights
  
  
 protected:

  void initializeInternals();
  
  // I delete maps left in memory next time process() is called
  // if getMap() or getZoomMap() is called, I set the internal pointer to NULL
  // and it is the user's responsbility to delete
  InterferometricMap* coarseMaps[AnitaPol::kNotAPol]; // these guys do the whole 360 az, and defined elevation...

  void fillWaveformInfo(AnitaPol::AnitaPol_t pol, AnitaEventSummary::WaveformInfo& info, const FilteredAnitaEvent* fEv,
                        AnalysisWaveform** waveStore, InterferometricMap* h, NoiseMonitor* noiseMonitor);
  
  InterferometricMap* fineMaps[AnitaPol::kNotAPol][AnitaEventSummary::maxDirectionsPerPol];
  AnalysisWaveform* wfCoherentFiltered[AnitaPol::kNotAPol][AnitaEventSummary::maxDirectionsPerPol][AnitaPol::kNotAPol];
  AnalysisWaveform* wfCoherent[AnitaPol::kNotAPol][AnitaEventSummary::maxDirectionsPerPol][AnitaPol::kNotAPol];
  AnalysisWaveform* wfDeconvolved[AnitaPol::kNotAPol][AnitaEventSummary::maxDirectionsPerPol][AnitaPol::kNotAPol];
  AnalysisWaveform* wfDeconvolvedFiltered[AnitaPol::kNotAPol][AnitaEventSummary::maxDirectionsPerPol][AnitaPol::kNotAPol];
  std::vector<Int_t> phiSectorToAnts[NUM_PHI]; //!< Vector of antennas to use in coherently summed waveforms
  
  AnitaEventSummary summary;
  // lazily generates CrossCorrelator when process() is called
  // (plan to add attachment function for external cross correlator)
  // if spawns one, sets the bool to true and deletes the cross correlator in destructor
  CrossCorrelator* fCrossCorr;

  bool spawnedCrossCorrelator;
  
  // in theory this could change if I end up making some settings dynamic, e.g. for MagicDisplay.
  InterferometryCache dtCache;
  
  FilterStrategy* fMinFilter; //!< Would be no filters, but for ANITA-3 data we need an ALFA filter
  FilterStrategy* fMinDecoFilter; //!< Minimum filter + deconvolution filter

  void chooseAntennasForCoherentlySumming(int coherentDeltaPhi);
  
  ANALYSIS_SETTING(Int_t, CoherentDeltaPhi);
  ANALYSIS_SETTING(Int_t, WhichResponseDir);
  ANALYSIS_SETTING(Int_t, UseOffAxisDelay);
  ANALYSIS_SETTING(Int_t, ResponseNPad);
  ANALYSIS_SETTING(Int_t, NumPeaks);  
  ClassDef(AnalysisReco, 0)

};

}


#endif
