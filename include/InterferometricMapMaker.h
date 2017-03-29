#ifndef INTERFEROMETRIC_MAP_MAKER_H
#define INTERFEROMETRIC_MAP_MAKER_H

#include "AnitaEventReconstructor.h"
#include "InterferometricMap.h"
#include "InterferometryCache.h"
#include "CrossCorrelator.h"




class InterferometryCache;


/** 
 * Class to make interferometric maps from FilteredAnitaEvents...
 */
class InterferometricMapMaker : public AnitaEventReconstructor{


public:
  
  virtual void process(const FilteredAnitaEvent * ev, UsefulAdu5Pat* pat ,AnitaEventSummary * summary) const;

  InterferometricMapMaker();
  virtual ~InterferometricMapMaker();

  void initializeInternals();


  Double_t getDeltaTExpected(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2,
			     Double_t phiWave, Double_t thetaWave) const;

  Double_t relativeOffAxisDelay(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2, Double_t phiDeg) const;
  Double_t singleAntennaOffAxisDelay(Double_t deltaPhiDeg) const;


  void drawSummary(TPad* pad, AnitaPol::AnitaPol_t pol);
  
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
  
  InterferometricMap* getMap(AnitaPol::AnitaPol_t pol);
  InterferometricMap* getZoomMap(AnitaPol::AnitaPol_t pol, UInt_t peakInd=0);

  void reconstruct(AnitaPol::AnitaPol_t pol) const;
  void reconstructZoom(AnitaPol::AnitaPol_t pol, Int_t peakIndex, Double_t zoomCenterPhiDeg, Double_t zoomCenterThetaDeg) const;

  AnalysisWaveform* coherentlySum(std::vector<const AnalysisWaveform*>& waves, std::vector<Double_t>& dts);

  static Int_t directlyInsertGeometry(TString pathToLindasFile, AnitaPol::AnitaPol_t pol);
  void insertPhotogrammetryGeometry();

  std::vector<Double_t> rArray[AnitaPol::kNotAPol]; //!< Vector of antenna radial positions
  std::vector<Double_t> phiArrayDeg[AnitaPol::kNotAPol]; //!< Vector of antenna azimuth positions
  std::vector<Double_t> zArray[AnitaPol::kNotAPol]; //!< Vector of antenna heights

  Int_t kUseOffAxisDelay; //!< Flag for whether or not to apply off axis delay to deltaT expected.
  Double_t maxDPhiDeg; //!< Variable for testing how wide an off axis angle is used in reconstruction
  Int_t coherentDeltaPhi;

  
  mutable InterferometricMap* coarseMaps[AnitaPol::kNotAPol]; // these guys do the whole 360 az, and defined elevation...
  mutable std::map<Int_t, InterferometricMap*> fineMaps[AnitaPol::kNotAPol]; // map of peak index to finely binned InterferometricMap
  
private:
  void initializeVariables();

  mutable CrossCorrelator* cc;  
  mutable bool spawnedCrossCorrelator;
  
  // in theory this could change if I end up making some settings dynamic, e.g. for MagicDisplay
  mutable InterferometryCache dtCache;
  // The axes of the interferometric maps store the uneven bins in theta (degrees)
  // these will replace the internal map storage...
  // we will imply some pointer ownership scheme with these.
  // I delete them if the event number changes and they've not been returned. You delete them otherwise...  

  
};




#endif
