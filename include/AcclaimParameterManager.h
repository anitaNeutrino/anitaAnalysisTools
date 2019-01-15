#ifndef ACCLAIM_FIT_PARAMS_H
#define ACCLAIM_FIT_PARAMS_H

#include "AnitaConventions.h"

class AnitaGeomTool;
class Adu5Pat;

namespace Acclaim {

  namespace PhaseCenterFit {

    enum class PhysicalRing {TopHigh, TopLow, Middle, Bottom};
    static PhysicalRing antToPhysicalRing(Int_t ant){
      int ring = 1 + ant/NUM_PHI;
      if(ant < NUM_PHI && (ant % 2)==0){
	ring -= 1;
      }
      return static_cast<PhysicalRing>(ring);
    }

    enum class ParameterSpace {None,
			       PitchRollHeading,
			       RingR,
			       RingPhi,
			       RingZ,
			       RingPhiRZ,
			       RingEllipse,
			       ExtraDeltaT
    };



    class EllipseParams {
    public:
      double x0 = 0; ///x-coordinate of ellipse center (m)
      double y0 = 0; ///y-coordinate of ellipse center (m)
      double alpha = 0; /// angle (radians) between x-axis in semi-major axis
      double Ra = 1; /// length of semi-major axis (m)
      double eccentricity = 0; /// = sqrt(1 - (Rb*Rb)/(Ra*Ra)) where Ra, Rb are the lengths of the major/minor axes
      double z = 0;

      EllipseParams(const double* params = nullptr);
      static int N(){return 6;}
      inline double Rb() const;
      void fill(const double* params);
      static const char* name(int p);
      void phiToEllipseXY(double phi, double& x, double& y) const;
      void phiToEllipseRZ(double phi, double& r, double& z) const;
    private:
      double tFromPhi(double phi) const;
    };



    /**
     * @class ParameterManager
     * @brief Class to handle the nitty gritty of how the fit params map to AnitaGeomTool/Adu5Pat
     */

    class ParameterManager {
    public:
      ParameterManager(ParameterSpace ps = ParameterSpace::None,  int nParams=0, const double* params=nullptr);
      virtual ~ParameterManager();
      void update(const double* params);
      void applyGeom(AnitaGeomTool* geom) const;
      void applyPat(Adu5Pat* pat) const ;
      const char* name(int param) const;
      int N() const {return fN;}

    private:
      const int fN;
      const double* fParams;
      const ParameterSpace fParamSpace;
      std::vector<EllipseParams> fEllipseParams;
  };


  }

}


std::ostream& operator<<(std::ostream& os, const Acclaim::PhaseCenterFit::PhysicalRing& r);


#endif // ACCLAIM_FIT_PARAMS_H
