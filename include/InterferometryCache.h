#ifndef INTERFEROMETRY_CACHE_H
#define INTERFEROMETRY_CACHE_H

class InterferometricMapMaker;
class CrossCorrelator;

#include <vector>

/** 
 * Class to cache the deltaTs for fast map making. It's just a glorified array...
 * 
 */
class InterferometryCache {

public:
  InterferometryCache();
  InterferometryCache(CrossCorrelator* cc, InterferometricMapMaker* mm);
  void populateCache(CrossCorrelator* cc, InterferometricMapMaker* mm);
  

  // pretty please inline
  inline int coarseIndex(int pol, int combo, int phiBin, int thetaBin){  
    return ((pol*numCombos + combo)*nCoarseBinsPhi + phiBin)*nCoarseBinsTheta + thetaBin;
  }

  inline double dt(int pol, int combo, int phiBin, int thetaBin){
    return deltaTs[coarseIndex(pol, combo, phiBin, thetaBin)];
  }
  
private:


  int nCoarseBinsPhi;
  int nCoarseBinsTheta;  
  int numCombos;
  std::vector<double> deltaTs; ///!< This is for the coarsely binned map




  // std::vector<double> zoomedThetaWaves; // NUM_BINS_PHI_ZOOM_TOTAL
  // std::vector<double> zoomedTanThetaWaves; // NUM_BINS_THETA_ZOOM_TOTAL
  // std::vector<double> zoomedCosThetaWaves; // NUM_BINS_THETA_ZOOM_TOTAL
  // std::vector<double> dtFactors; //
  // Double_t zoomedTanThetaWaves[NUM_BINS_THETA_ZOOM_TOTAL]; //!< Cached tan(theta) for zoomed image.
  // Double_t zoomedCosThetaWaves[NUM_BINS_THETA_ZOOM_TOTAL]; //!< Cached cos(theta) for zoomed image.
  // Double_t dtFactors[NUM_BINS_THETA_ZOOM_TOTAL]; //!< Cached cos(theta)/c/dt for zoomed image.
  // Double_t zoomedPhiWaveLookup[NUM_BINS_PHI_ZOOM_TOTAL]; //!< Cached phi for zoomed image.
  // Double_t zoomedCosPartLookup[AnitaPol::kNotAPol][NUM_SEAVEYS][NUM_BINS_PHI_ZOOM_TOTAL]; //!< Cached part of the deltaT calculation.
  // Double_t offAxisDelays[AnitaPol::kNotAPol][NUM_COMBOS][NUM_BINS_PHI_ZOOM_TOTAL]; //!< Off-axis delays for fine binned images.
  // Double_t offAxisDelaysDivided[AnitaPol::kNotAPol][NUM_COMBOS][NUM_BINS_PHI_ZOOM_TOTAL]; //!< Off-axis delays divided such to remove an operation from the inner loop of an im

  
  // Double_t zoomedThetaWaves[NUM_BINS_THETA_ZOOM_TOTAL]; //!< Cached theta for zoomed image.
  // Double_t zoomedTanThetaWaves[NUM_BINS_THETA_ZOOM_TOTAL]; //!< Cached tan(theta) for zoomed image.
  // Double_t zoomedCosThetaWaves[NUM_BINS_THETA_ZOOM_TOTAL]; //!< Cached cos(theta) for zoomed image.
  // Double_t dtFactors[NUM_BINS_THETA_ZOOM_TOTAL]; //!< Cached cos(theta)/c/dt for zoomed image.
  // Double_t zoomedPhiWaveLookup[NUM_BINS_PHI_ZOOM_TOTAL]; //!< Cached phi for zoomed image.
  // Double_t zoomedCosPartLookup[AnitaPol::kNotAPol][NUM_SEAVEYS][NUM_BINS_PHI_ZOOM_TOTAL]; //!< Cached part of the deltaT calculation.
  // Double_t offAxisDelays[AnitaPol::kNotAPol][NUM_COMBOS][NUM_BINS_PHI_ZOOM_TOTAL]; //!< Off-axis delays for fine binned images.
  // Double_t offAxisDelaysDivided[AnitaPol::kNotAPol][NUM_COMBOS][NUM_BINS_PHI_ZOOM_TOTAL]; //!< Off-axis delays divided such to remove an operation from the inner loop of an image making function.
  
  
};


#endif
