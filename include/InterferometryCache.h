#ifndef INTERFEROMETRY_CACHE_H
#define INTERFEROMETRY_CACHE_H

#include "AnitaConventions.h"


#include <vector>


namespace Acclaim
{


  class AnalysisReco;
  class CrossCorrelator;
  
  /** 
   * Class to cache the deltaTs for fast map making. 
   * Really it's a glorified set of vectors, but I like to think of it as a lean, mean caching machine...
   */
  class InterferometryCache {

    friend class InterferometricMap;

  public:
    InterferometryCache();
    InterferometryCache(CrossCorrelator* cc, const AnalysisReco* reco);

    void init(CrossCorrelator* cc, const AnalysisReco* reco, bool forceCacheRecalculation = false);
  
    // pretty please inline
    inline int coarseIndex(int pol, int combo, int phiBin, int thetaBin){  
      return ((pol*numCombos + combo)*nCoarseBinsPhi + phiBin)*nCoarseBinsTheta + thetaBin;
    }

    inline double coarseDt(int pol, int combo, int phiBin, int thetaBin){
      return deltaTs[coarseIndex(pol, combo, phiBin, thetaBin)];
    }

    inline int partBAsIndex(int pol, int combo, int fineThetaBin){
      return (pol*numCombos + combo)*nFineBinsTheta + fineThetaBin;
    }

    inline int zoomedCosPartIndex(int pol, int ant, int finePhiBin){
      return (pol*NUM_SEAVEYS + ant)*nFineBinsPhi + finePhiBin;
    }

    // inline int offAxisDelayIndex(int pol, int combo, int finePhiBin){
    inline int part21sIndex(int pol, int combo, int finePhiBin){
      return (pol*numCombos + combo)*nFineBinsPhi + finePhiBin;
    }

    int kUseOffAxisDelay;
  private:  


    void populateCache(CrossCorrelator* cc, const AnalysisReco* reco);
    void populateFineCache(CrossCorrelator* cc, const AnalysisReco* reco);
    bool initialized;
  
    int nCoarseBinsPhi;
    int nCoarseBinsTheta;  
    int numCombos;
    std::vector<double> deltaTs; ///!< This is for the coarsely binned map

    int nFineBinsPhi;
    int nFineBinsTheta;

  
    // these represent the quite complicated caching of finely mapped dts so they scale well as the number of possible fine bins increases
    // I should give them accessor variables and make them private... but for now...

    std::vector<double> partBAsZoom; //[AnitaPol::kNotAPol][NUM_COMBOS][NUM_BINS_THETA_ZOOM_TOTAL]; //!< Yet more geometric caching
    std::vector<double> part21sZoom; //[AnitaPol::kNotAPol][NUM_COMBOS][NUM_BINS_PHI_ZOOM_TOTAL]; //!< Yet more geometric caching
    std::vector<double> zoomedThetaWaves; // NUM_BINS_PHI_ZOOM_TOTAL
    std::vector<double> zoomedTanThetaWaves; // NUM_BINS_THETA_ZOOM_TOTAL
    std::vector<double> zoomedCosThetaWaves; // NUM_BINS_THETA_ZOOM_TOTAL
    std::vector<double> dtFactors; //
    std::vector<double> zoomedPhiWaveLookup; //[NUM_BINS_PHI_ZOOM_TOTAL]; //!< Cached phi for zoomed image.
    std::vector<double> zoomedCosPartLookup; // [AnitaPol::kNotAPol][NUM_SEAVEYS][NUM_BINS_PHI_ZOOM_TOTAL]; //!< Cached part of the deltaT calculation.
    std::vector<double> offAxisDelays; //[AnitaPol::kNotAPol][NUM_COMBOS][NUM_BINS_PHI_ZOOM_TOTAL]; //!< Off-axis delays for fine binned images.
    std::vector<double> offAxisDelaysDivided; // [AnitaPol::kNotAPol][NUM_COMBOS][NUM_BINS_PHI_ZOOM_TOTAL]; //!< Off-axis delays divided such to remove an operation from the inner   
  };
}

#endif
