#ifndef ACCLAIM_FIT_PARAMS_H
#define ACCLAIM_FIT_PARAMS_H

#include "TMath.h"
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
      EllipseParams(const double* params = nullptr){
	if(params){fill(params);}
      }
      double x0 = 0; ///x-coordinate of ellipse center (m)
      double y0 = 0; ///y-coordinate of ellipse center (m)
      double alpha = 0; /// angle (radians) between x-axis in semi-major axis
      double Ra = 1; /// length of semi-major axis (m)
      double eccentricity = 0; /// = sqrt(1 - (Rb*Rb)/(Ra*Ra)) where Ra, Rb are the lengths of the major/minor axes
      double z = 0;

      static int N() {
	return  6;
      }

      inline double Rb() const {
	return Ra* TMath::Sqrt(1 - eccentricity*eccentricity);
      }

      void fill(const double* params){
	x0           = params[0];
	y0           = params[1];
	alpha        = params[2];
	Ra           = params[3];
	eccentricity = params[4];
	z            = params[5];
      }

      static const char* name(int p){
	switch(p){
	case 0: return "x0";
	case 1: return "y0";
	case 2: return "alpha";
	case 3: return "Ra";
	case 4: return "eccentricity";
	case 5: return "z";
	default: return "Unknown!";
	}      
      }

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

      virtual void applyGeom(AnitaGeomTool* geom) const;
      virtual void applyPat(Adu5Pat* pat) const ;

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


#endif // ACCLAIM_FIT_PARAMS_H
