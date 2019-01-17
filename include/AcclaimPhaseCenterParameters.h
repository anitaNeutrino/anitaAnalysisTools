#ifndef ACCLAIM_PHASE_CENTER_PARAMS_H
#define ACCLAIM_PHASE_CENTER_PARAMS_H

#include "AnitaConventions.h"
#include <array>

class AnitaGeomTool;
class Adu5Pat;

namespace ROOT {
  namespace Math {
    class Minimizer;
  }
}

namespace Acclaim {

  namespace PhaseCenter {

    enum class ParameterSpace {None,
			       PitchRoll,
			       PitchRollHeading,
			       RingR,
			       RingPhi,
			       RingZ,
			       RingPhiRZ,
			       RingEllipse,
			       ExtraDeltaT
    };
    const char* toCString(ParameterSpace ps);
    
    
    enum class PhysicalRing {TopHigh,
			     TopLow,
			     Middle,
			     Bottom};
    const char* toCString(PhysicalRing r);
    PhysicalRing antToPhysicalRing(int ant);
    constexpr std::array<PhysicalRing, 4> PhysicalRings {PhysicalRing::TopHigh,
							 PhysicalRing::TopLow,
							 PhysicalRing::Middle,
							 PhysicalRing::Bottom};


    
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
      ParameterManager(ParameterSpace ps = ParameterSpace::None);
      virtual ~ParameterManager();

      void update(const double* params);
      void setInputs(ROOT::Math::Minimizer* min, std::vector<double>& inputs) const;
      void applyDelay(double& dt, AnitaPol::AnitaPol_t pol, int ant1, int ant2) const;      
      void applyGeom(AnitaGeomTool* geom) const;
      void applyPat(Adu5Pat* pat) const ;

      const char* name(int param) const;
      int N() const {return fN;}
      ParameterSpace space() const {return fParamSpace;}
    private:
      int fN;
      const double* fParams;
      ParameterSpace fParamSpace;
      std::vector<EllipseParams> fEllipseParams;
  };


  }

}


std::ostream& operator<<(std::ostream& os, const Acclaim::PhaseCenter::PhysicalRing& r);


#endif // ACCLAIM_PHASE_CENTER_PARAMS_H
