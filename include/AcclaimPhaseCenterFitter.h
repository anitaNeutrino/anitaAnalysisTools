#ifndef ACCLAIM_PHASE_CENTER_FITTER_H
#define ACCLAIM_PHASE_CENTER_FITTER_H

#include "AcclaimCorrelationSummary.h"
#include <map>
#include "Math/Functor.h"
#include "AnitaConventions.h"
#include "TMath.h"
#include "TVector2.h"

class AnitaGeomTool;
class TH2D;
class TChain;

namespace ROOT {
  namespace Math {
    class Minimizer;
  }
}

namespace Acclaim {

  /**
   * @class PhaseCenterFitter
   * @brief Well, it fits the phase centers!
   * 
   * Reads the output of make correlation summary
   */

  class PhaseCenterFitter {
  public:
    PhaseCenterFitter(const char* corrTreeFiles = nullptr);
    ~PhaseCenterFitter();

    enum class ParameterSpace {PitchRollHeading,
			       RingR,
			       RingPhi,
			       RingZ,
			       RingPhiRZ,
			       RingEllipse
    };

    enum class PhysicalRing {TopHigh, TopLow, Middle, Bottom};
    static PhysicalRing antToPhysicalRing(Int_t ant){
      int ring = 1 + ant/NUM_PHI;
      if(ant < NUM_PHI && (ant % 2)==0){
	ring -= 1;
      }
      return static_cast<PhysicalRing>(ring);
    }
    

    void setFitParameterSpace(ParameterSpace ps);
    const std::vector<double>& fit(AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, const char* outFileName = nullptr);
    void printResults() const;
    void setPrintOnEval(bool print){
      fPrintOnEval = print;
    }
    
    
  private:
    void readInSummaries();
    void makeFunctors();
    double eval(const double* params);
    TH2D* makeDdtHist(AnitaPol::AnitaPol_t pol, const TString& name, const char* title = nullptr);

    const double fCorrelationThreshold = 0.4;
    const double fDeltaDeltaTThreshold = 1.0; /// ns
    TChain* fChain = nullptr;
    CorrelationSummary* fSum = nullptr;
    std::vector<Acclaim::CorrelationSummary> fSummaries;
    ROOT::Math::Minimizer* fMin = nullptr;
    bool fPrintOnEval = false;
    AnitaPol::AnitaPol_t fFitPol = AnitaPol::kNotAPol;

    double fMinimum = 0;
    std::vector<double> fResults;

    double fInitial = 0;
    std::vector<double> fInputs;
    
    std::map<ParameterSpace, ROOT::Math::Functor> fFuncs;
    ParameterSpace fParamSpace;

    std::array<std::array<Int_t, NUM_SEAVEYS>, AnitaPol::kNotAPol> fNormalization;
    std::array<std::array<Double_t, NUM_SEAVEYS>, AnitaPol::kNotAPol> fDdts;

    bool fApplyParams = true;



  public:
    class EllipseParams {
    public:
      EllipseParams(const double* params = nullptr){
	if(params){fill(params);}
      }
      double x0 = 0; ///x-coordinate of ellipse center (m)
      double y0 = 0; ///y-coordinate of ellipse center (m)
      double alpha = 0; /// angle (radians) between x-axis in semi-major axis
      double Ra = 1; /// length of semi-major axis (m)
      double Rb = 1; /// length of semi-minor axis (m)

      void fill(const double* params){
	x0    = params[0];
	y0    = params[1];
	alpha = params[2];
	Ra    = params[3];
	Rb    = params[4];
      }

      static const char* name(int p){
	switch(p){
	case 0: return "x0";
	case 1: return "y0";
	case 2: return "alpha";
	case 3: return "Ra";
	case 4: return "Rb";
	default: return "Unknown!";
	}      
      }

      void phiToEllipseXY(double phi, double& x, double& y);
      
    private:
      double tFromPhi(double phi);
    };

    
      
  };


  
  
  class FakeGeomTool {
  public:
    FakeGeomTool(const AnitaGeomTool* geom);
    inline Double_t getAntPhi(Int_t ant, AnitaPol::AnitaPol_t pol=AnitaPol::kVertical) const {
      return get(fFittedPhi, fPhotoPhi, pol, ant);
    }
    inline Double_t getAntR(Int_t ant, AnitaPol::AnitaPol_t pol=AnitaPol::kVertical) const {
      return get(fFittedR, fPhotoR, pol, ant);
    }
    inline Double_t getAntZ(Int_t ant, AnitaPol::AnitaPol_t pol=AnitaPol::kVertical) const {
      return get(fFittedZ, fPhotoZ, pol, ant);
    }
    void restorePhotogrammetryPositionsInAnitaGeomTool(AnitaGeomTool* geom) const;
    void overwritePhotogrammetryPositionsInAnitaGeomTool(AnitaGeomTool* geom) const;

    void print() const;
    
  private:

    typedef std::map<std::pair<AnitaPol::AnitaPol_t, Int_t>, double> PolAntMap;
    PolAntMap fPhotoR; /// pol/ant to radial (m)
    PolAntMap fPhotoZ; /// pol/ant to z position (m)
    PolAntMap fPhotoPhi; /// pol/ant to phi in radians

    typedef std::map<std::pair<AnitaPol::AnitaPol_t, PhaseCenterFitter::PhysicalRing>, double> PolRingMap;
    PolRingMap fFittedR; /// pol/ant to radial (m)
    PolRingMap fFittedZ; /// pol/ant to z position (m)
    PolRingMap fFittedPhi; /// pol/ant to phi in radians

    inline double get(const PolRingMap& fittedVals, const PolAntMap& defaultVals,
				    AnitaPol::AnitaPol_t pol, Int_t ant) const {      
      auto fittedKey = std::make_pair(pol, PhaseCenterFitter::antToPhysicalRing(ant));
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

std::ostream& operator<<(std::ostream& os, const Acclaim::PhaseCenterFitter::PhysicalRing& r);

#endif 
