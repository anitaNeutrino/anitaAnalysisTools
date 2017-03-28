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
  std::vector<double> deltaTs;
  
  
};


#endif
