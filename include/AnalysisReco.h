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

namespace Acclaim
{

  class InterferometryCache;


/** 
 * Class to make interferometric maps from FilteredAnitaEvents...
 */

class AnalysisReco : public TObject, public AnitaEventReconstructor {

 public:
  
  virtual void process(const FilteredAnitaEvent * ev, UsefulAdu5Pat* pat ,AnitaEventSummary * summary) const;
  virtual void process(const FilteredAnitaEvent * ev, Adu5Pat* pat, AnitaEventSummary * summary) const;    

  AnalysisReco();
  virtual ~AnalysisReco();

  void initializeInternals();


  Double_t getDeltaTExpected(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2,
                             Double_t phiWave, Double_t thetaWave) const;

  Double_t relativeOffAxisDelay(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2, Double_t phiDeg) const;
  Double_t singleAntennaOffAxisDelay(Double_t deltaPhiDeg) const;


  void drawSummary(TPad* pad, AnitaPol::AnitaPol_t pol);
  
  // TGraph* makeCoherentlySummedWaveform(AnitaPol::AnitaPol_t pol, Double_t phiDeg,
  //     				 Double_t thetaDeg, Int_t maxDeltaPhiSect, Double_t& snr);
  // TGraph* makeUpsampledCoherentlySummedWaveform(AnitaPol::AnitaPol_t pol, Double_t phiDeg,
  //     					  Double_t thetaDeg, Int_t maxDeltaPhiSect, Double_t&
  // snr);
  // TGraph* makeCoherentWorker(AnitaPol::AnitaPol_t pol, Double_t phiDeg, Double_t thetaDeg,
  //     		       Int_t maxDeltaPhiSect, Double_t& snr, Int_t nSamp);






  // void getCoarsePeakInfo(AnitaPol::AnitaPol_t pol, Int_t peakIndex, Double_t& value,
  //     		   Double_t& phiDeg, Double_t& thetaDeg);
  // void getFinePeakInfo(AnitaPol::AnitaPol_t pol, Int_t peakIndex, Double_t& value,
  //     		 Double_t& phiDeg, Double_t& thetaDeg);
  
  InterferometricMap* getMap(AnitaPol::AnitaPol_t pol);
  InterferometricMap* getZoomMap(AnitaPol::AnitaPol_t pol, UInt_t peakInd=0);

  void reconstruct(AnitaPol::AnitaPol_t pol) const;
  void reconstructZoom(AnitaPol::AnitaPol_t pol, Int_t peakIndex, Double_t zoomCenterPhiDeg, Double_t zoomCenterThetaDeg) const;

  AnalysisWaveform* coherentlySum(std::vector<const AnalysisWaveform*>& waves, std::vector<Double_t>& dts) const;
  AnalysisWaveform* coherentlySum(const FilteredAnitaEvent* fEv, const InterferometricMap* h, bool swapPol = false) const;
    
  static Int_t directlyInsertGeometry(TString pathToLindasFile, AnitaPol::AnitaPol_t pol);
  void insertPhotogrammetryGeometry();

  std::vector<Double_t> rArray[AnitaPol::kNotAPol]; //!< Vector of antenna radial positions
  std::vector<Double_t> phiArrayDeg[AnitaPol::kNotAPol]; //!< Vector of antenna azimuth positions
  std::vector<Double_t> zArray[AnitaPol::kNotAPol]; //!< Vector of antenna heights

  
 protected:
  
  // I delete maps left in memory next time process() is called
  // if getMap() or getZoomMap() is called, I set the internal pointer to NULL
  // and it is the user's responsbility to delete
  mutable InterferometricMap* coarseMaps[AnitaPol::kNotAPol]; // these guys do the whole 360 az, and defined elevation...

  void putInCoherentMap(std::map<Int_t, AnalysisWaveform*>& theMap, Int_t peakInd, AnalysisWaveform* coherentWave) const;
  void fillWaveformInfo(AnitaPol::AnitaPol_t pol,
                        Int_t peakInd,
                        AnitaEventSummary::WaveformInfo& info,
                        const FilteredAnitaEvent* fEv,
                        std::map<Int_t, AnalysisWaveform*>* waveStore,
                        InterferometricMap* h) const;
  
  mutable std::map<Int_t, InterferometricMap*> fineMaps[AnitaPol::kNotAPol]; // map of peak index to finely binned InterferometricMap
  mutable std::map<Int_t, AnalysisWaveform*> wfCoherentFiltered[AnitaPol::kNotAPol][AnitaPol::kNotAPol];
  mutable std::map<Int_t, AnalysisWaveform*> wfCoherent[AnitaPol::kNotAPol][AnitaPol::kNotAPol];
  mutable std::map<Int_t, AnalysisWaveform*> wfDeconvolved[AnitaPol::kNotAPol][AnitaPol::kNotAPol];
  mutable std::map<Int_t, AnalysisWaveform*> wfDeconvolvedFiltered[AnitaPol::kNotAPol][AnitaPol::kNotAPol];

  mutable AnitaEventSummary summary;
  // lazily generates CrossCorrelator when process() is called
  // (plan to add attachment function for external cross correlator)
  // if spawns one, sets the bool to true and deletes the cross correlator in destructor
  mutable CrossCorrelator* cc;

  mutable bool spawnedCrossCorrelator;
  
  // in theory this could change if I end up making some settings dynamic, e.g. for MagicDisplay.
  mutable InterferometryCache dtCache;
  
  FilterStrategy* fMinFilter; //!< Would be no filters, but for ANITA-3 data we need an ALFA filter
  FilterStrategy* fMinDecoFilter; //!< Minimum filter + deconvolution filter


  void swapPol(AnitaPol::AnitaPol_t& pol) const;
  
  ANALYSIS_SETTING(Int_t, CoherentDeltaPhi);
  ANALYSIS_SETTING(Int_t, WhichResponseDir);
  ANALYSIS_SETTING(Int_t, UseOffAxisDelay);
  ANALYSIS_SETTING(Int_t, ResponseNPad);  
  ClassDef(AnalysisReco, 0)
};

}


#endif
