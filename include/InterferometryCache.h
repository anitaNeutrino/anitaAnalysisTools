#ifndef INTERFEROMETRY_CACHE_H
#define INTERFEROMETRY_CACHE_H

#include "AnitaConventions.h"

class InterferometricMapMaker;
class CrossCorrelator;

#include <vector>

/** 
 * Class to cache the deltaTs for fast map making. 
 * Really it's a glorified set of vectors, but I like to think of it as a lean, mean caching machine...
 */
class InterferometryCache {

public:
  InterferometryCache();
  InterferometryCache(CrossCorrelator* cc, InterferometricMapMaker* mm);
  void populateCache(CrossCorrelator* cc, InterferometricMapMaker* mm);
  void populateFineCache(CrossCorrelator* cc, InterferometricMapMaker* mm);

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

  inline double fineDt(int pol, int combo, int phiBin, int thetaBin){
    return 0;
  }
  
  
private:


  int nCoarseBinsPhi;
  int nCoarseBinsTheta;  
  int numCombos;
  std::vector<double> deltaTs; ///!< This is for the coarsely binned map

  int nFineBinsPhi;
  int nFineBinsTheta;
public:  
  std::vector<double> partBAsZoom; //[AnitaPol::kNotAPol][NUM_COMBOS][NUM_BINS_THETA_ZOOM_TOTAL]; //!< Yet more geometric caching
  std::vector<double> part21sZoom; //[AnitaPol::kNotAPol][NUM_COMBOS][NUM_BINS_PHI_ZOOM_TOTAL]; //!< Yet more geometric caching
  // Double_t partBAsZoom[AnitaPol::kNotAPol][NUM_COMBOS][NUM_BINS_THETA_ZOOM_TOTAL]; //!< Yet more geometric caching
  // Double_t part21sZoom[AnitaPol::kNotAPol][NUM_COMBOS][NUM_BINS_PHI_ZOOM_TOTAL]; //!< Yet more geometric caching
  
  std::vector<double> zoomedThetaWaves; // NUM_BINS_PHI_ZOOM_TOTAL
  std::vector<double> zoomedTanThetaWaves; // NUM_BINS_THETA_ZOOM_TOTAL
  std::vector<double> zoomedCosThetaWaves; // NUM_BINS_THETA_ZOOM_TOTAL
  std::vector<double> dtFactors; //
  std::vector<double> zoomedPhiWaveLookup; //[NUM_BINS_PHI_ZOOM_TOTAL]; //!< Cached phi for zoomed image.
  std::vector<double> zoomedCosPartLookup; // [AnitaPol::kNotAPol][NUM_SEAVEYS][NUM_BINS_PHI_ZOOM_TOTAL]; //!< Cached part of the deltaT calculation.
  std::vector<double> offAxisDelays; //[AnitaPol::kNotAPol][NUM_COMBOS][NUM_BINS_PHI_ZOOM_TOTAL]; //!< Off-axis delays for fine binned images.
  std::vector<double> offAxisDelaysDivided; // [AnitaPol::kNotAPol][NUM_COMBOS][NUM_BINS_PHI_ZOOM_TOTAL]; //!< Off-axis delays divided such to remove an operation from the inner loop of an image making function.

  
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
