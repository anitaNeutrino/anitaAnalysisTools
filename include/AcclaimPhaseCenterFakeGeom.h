#ifndef ACCLAIM_PHASE_CENTER_FAKE_GEOM_H
#define ACCLAIM_PHASE_CENTER_FAKE_GEOM_H

#include "AnitaConventions.h"
#include "AcclaimPhaseCenterParameters.h"
#include <map>

class AnitaGeomTool;

namespace Acclaim {
  namespace PhaseCenter {


    /**
     * @class FakeGeomTool
     * @brief Caches the  "true" AnitaGeomTool positions for the minimizer
     */

    class FakeGeomTool {
    public:
      FakeGeomTool(const AnitaGeomTool* geom);
      inline Double_t getAntPhi(int ant, AnitaPol::AnitaPol_t pol=AnitaPol::kVertical) const {
	return get(fFittedPhi, fPhotoPhi, pol, ant);
      }
      inline Double_t getAntR(int ant, AnitaPol::AnitaPol_t pol=AnitaPol::kVertical) const {
	return get(fFittedR, fPhotoR, pol, ant);
      }
      inline Double_t getAntZ(int ant, AnitaPol::AnitaPol_t pol=AnitaPol::kVertical) const {
	return get(fFittedZ, fPhotoZ, pol, ant);
      }
      void restorePhotogrammetryPositionsInAnitaGeomTool(AnitaGeomTool* geom) const;
      void overwritePhotogrammetryPositionsInAnitaGeomTool(AnitaGeomTool* geom) const;

      void print() const;
    
    private:

      typedef std::map<std::pair<AnitaPol::AnitaPol_t, int>, double> PolAntMap;
      PolAntMap fPhotoR; /// pol/ant to radial (m)
      PolAntMap fPhotoZ; /// pol/ant to z position (m)
      PolAntMap fPhotoPhi; /// pol/ant to phi in radians

      typedef std::map<std::pair<AnitaPol::AnitaPol_t, PhysicalRing>, double> PolRingMap;
      PolRingMap fFittedR; /// pol/ant to radial (m)
      PolRingMap fFittedZ; /// pol/ant to z position (m)
      PolRingMap fFittedPhi; /// pol/ant to phi in radians

      inline double get(const PolRingMap& fittedVals, const PolAntMap& defaultVals,
			AnitaPol::AnitaPol_t pol, int ant) const {      
	auto fittedKey = std::make_pair(pol, antToPhysicalRing(ant));
	if(fittedVals.find(fittedKey)!=fittedVals.end()){
	  return fittedVals.at(fittedKey);
	}
	else{
	  auto defaultKey = std::make_pair(pol, ant);
	  return defaultVals.at(defaultKey);
	}
      }    
    };
    
  }
}

#endif // ACCLAIM_PHASE_CENTER_FAKE_GEOM_H
