#ifndef INTERFEROMETRY_CACHE_H
#define INTERFEROMETRY_CACHE_H

#include "AnitaConventions.h"


#include <vector>


namespace Acclaim
{


  class AnalysisReco;
  class CrossCorrelator;
  
  /** 
   * @brief Class to cache the deltaTs or parts of their calculation for fast map making.
   * 
   * It's a glorified set of std::vectors with some indexing logic
   */
  class InterferometryCache {

    friend class InterferometricMap;

  public:
    /** 
     * Default constructor, not very useful
     */
    InterferometryCache();

    /** 
     * Useful constructor
     * 
     * @param cc is a pointer to a CrossCorrelator
     * @param reco is a pointer to an AnalysisReco
     */
    InterferometryCache(CrossCorrelator* cc, const AnalysisReco* reco);

    /** 
     * @brief Initialize the cache, only needs to be done once
     * 
     * @param cc is a pointer to the CrossCorrelator
     * @param reco is a pointer to the AnalysisReco
     * @param forceCacheRecalculation set true to calculate  the cached delays even if they've already been calculated.
     */    
    void init(CrossCorrelator* cc, const AnalysisReco* reco, bool forceCacheRecalculation = false);
  
    /** 
     * @brief  Get the internal array index for the #fDeltaTs vector (for coarsely binned InterferometericMap)
     * 
     * @param pol is the polarization
     * @param combo is an index for the antenna pair
     * @param phiBin is the azimuthal bin (counting from 0) of the coarsely binned InterferometricMap
     * @param thetaBin is the elevation bin (counting from 0) of the coarsely binned InterferometricMap
     * 
     * @return the index for #fDeltaTs
     */
    inline int coarseIndex(int pol, int combo, int phiBin, int thetaBin){  
      return ((pol*fNumCombos + combo)*fNCoarseBinsPhi + phiBin)*fNCoarseBinsTheta + thetaBin;
    }

    /** 
     * @brief Look up the delay for the antenna pair / polarization at a given angle
     * 
     * @param pol is the polarization
     * @param combo is an index for the antenna pair
     * @param phiBin is the azimuthal bin (counting from 0) of the coarsely binned InterferometricMap
     * @param thetaBin is the elevation bin (counting from 0) of the coarsely binned InterferometricMap
     * 
     * @return the delay between the antenna pair
     */
    inline double coarseDt(int pol, int combo, int phiBin, int thetaBin){
      return fDeltaTs[coarseIndex(pol, combo, phiBin, thetaBin)];
    }

    /** 
     * @brief Get the index for part of the deltaT calculation for the finely binned InterferometricMap 
     * 
     * At the current resolution caching the full set of dts for the finely binned InterferometricMap is too large to calculate over all azimuth/elevation
     * 
     * @param pol is the polarization
     * @param combo is an index for the antenna pair
     * @param fineThetaBin the theta bin of the antenna pair (considering the full range of theta at the finely binned resolution)
     * 
     * @return the index for the #fPartBAsZoom
     */
    inline int partBAsIndex(int pol, int combo, int fineThetaBin){
      return (pol*fNumCombos + combo)*fNFineBinsTheta + fineThetaBin;
    }

    /** 
     * @brief Get the index for part of the deltaT calculation for the finely binned InterferometricMap 
     * 
     * At the current resolution caching the full set of dts for the finely binned InterferometricMap is too large to calculate over all azimuth/elevation
     * 
     * @param pol is the polarization
     * @param ant is the antenna
     * @param finePhiBin the phi bin of the antenna pair (considering the full range of azimuth at the finely binned resolution)
     * 
     * @return the index for the #fZoomedCosPartLookup
     */    
    inline int zoomedCosPartIndex(int pol, int ant, int finePhiBin){
      return (pol*NUM_SEAVEYS + ant)*fNFineBinsPhi + finePhiBin;
    }

    /** 
     * Get the index for the phi part of the deltaT calculation for the finely binned InterferometricMap
     * 
     * @param pol is the polarization
     * @param combo is an index for the antenna pair
     * @param finePhiBin the phi bin of the antenna pair (considering the full range of azimuth at the finely binned resolution)
     * 
     * @return the index for the #fPart21sZoom
     */
    inline int part21sIndex(int pol, int combo, int finePhiBin){
      return (pol*fNumCombos + combo)*fNFineBinsPhi + finePhiBin;
    }

    /** 
     * Populate the #fDeltaTs cache for the coarsely binned InterferometricMap
     * 
     * @param cc is the CrossCorrelator
     * @param reco is the AnalysisReco
     */
    void populateCache(CrossCorrelator* cc, const AnalysisReco* reco);

    /** 
     * Populate the various caches for the finely binned InterferometricMap
     * 
     * @param cc is the CrossCorrelator
     * @param reco is the AnalysisReco
     */
    void populateFineCache(CrossCorrelator* cc, const AnalysisReco* reco);


    int fUseOffAxisDelay; ///< Should we use the off-axis delay? Set to 1 or 0 as it's used in a multiplication for speed    
    bool fInitialized; ///< Has the cache been initialized?
    int fNCoarseBinsPhi; ///< The number of phi bins in the coarsely binned interferometric map
    int fNCoarseBinsTheta; ///< The number of theta bins in the coarsely binned interferometric map
    int fNumCombos; ///< The number of possible antenna pairs to be correlated
    std::vector<double> fDeltaTs; ///< Delays between antenna pairs for the coarsely binned interferometric map

    int fNFineBinsPhi; ///< If a finely binned InterferometricMap extended over all 360 in phi (plus edge effects) it would have this many bins
    int fNFineBinsTheta; ///< If a finely binned InterferometricMap extended over all theta (plus edge effects) it would have this many bins

    /**
     * The following variables represent a partial caching of the finely binned InterferometricMap dts.
     * Calculating the full set of deltaTs would use too much memory.
     * These are partially cached such that the minimize the dt calculations in an inner loop in InterferometricMap.
     * I'm don't think the cached quantities are very physically meaningful and as such they have silly names.
     */
    
    //[AnitaPol::kNotAPol][NUM_COMBOS][NUM_BINS_THETA_ZOOM_TOTAL]
    std::vector<double> fPartBAsZoom; ///< Partial caching for the finely binned InterferometricMap 

    //[AnitaPol::kNotAPol][NUM_COMBOS][NUM_BINS_PHI_ZOOM_TOTAL]
    std::vector<double> fPart21sZoom; ///< Partial caching for the finely binned InterferometricMap 

    // NUM_BINS_PHI_ZOOM_TOTAL
    std::vector<double> fZoomedThetaWaves; ///< Values of theta (radians) for the finely binned InterferometricMap

    // NUM_BINS_THETA_ZOOM_TOTAL
    std::vector<double> fZoomedTanThetaWaves; ///< Values of tan(theta) for the finely binned InterferometricMap

    // NUM_BINS_THETA_ZOOM_TOTAL
    std::vector<double> fZoomedCosThetaWaves; ///< Values of cos(theta) for the finely binned InterferometricMap

    // NUM_BINS_THETA_ZOOM_TOTAL
    std::vector<double> fDtFactors; ///< Multiplier factor for the finely binned InterferometricMap

    // NUM_BINS_PHI_ZOOM_TOTAL
    std::vector<double> fZoomedPhiWaveLookup; ///< The value of phi in each bin in radians

    // NUM_BINS_PHI_ZOOM_TOTAL
    std::vector<double> fZoomedCosPartLookup; ///< Partial caching for the of the deltaT calculation

    // [AnitaPol::kNotAPol][NUM_COMBOS][NUM_BINS_PHI_ZOOM_TOTAL]
    std::vector<double> fOffAxisDelays; ///< Off-axis delays for use when making fine binned InterferometricMap

    // [AnitaPol::kNotAPol][NUM_COMBOS][NUM_BINS_PHI_ZOOM_TOTAL]    
    std::vector<double> fOffAxisDelaysDivided; ///< Same as #fOffAxisDelays, but delays scaled to remove an operation from an inner loop in InterfermetricMap
  };
}

#endif
