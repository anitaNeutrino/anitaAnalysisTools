#ifndef ACCLAIM_PHASE_CENTER_FITTER_H
#define ACCLAIM_PHASE_CENTER_FITTER_H

#include "AcclaimCorrelationSummary.h"
#include <map>
#include "Math/Functor.h"
#include "AnitaConventions.h"
#include "TMath.h"
#include "TVector2.h"
#include "AcclaimParameterManager.h"

class AnitaGeomTool;
class TH2D;
class TChain;

namespace ROOT {
  namespace Math {
    class Minimizer;
  }
}

namespace Acclaim {

  namespace PhaseCenterFit {

    

    class FakeGeomTool;
    
    /**
     * @class Minimizer
     * @brief Well, it fits the phase centers!
     * 
     * Reads the output of make correlation summary
     */

    class Minimizer {
    public:
      Minimizer(const char* corrTreeFiles = nullptr);
    
      void setParameterSpace(ParameterSpace ps);
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
      std::shared_ptr<TChain> fChain = nullptr;
      CorrelationSummary* fSum = nullptr;
      std::vector<Acclaim::CorrelationSummary> fSummaries;
      std::shared_ptr<ROOT::Math::Minimizer> fMin = nullptr;
      bool fPrintOnEval = false;
      AnitaPol::AnitaPol_t fFitPol = AnitaPol::kNotAPol;

      double fMinimum = 0;
      std::vector<double> fResults;

      double fInitial = 0;
      std::vector<double> fInputs;
    
      std::map<ParameterSpace, ROOT::Math::Functor> fFuncs;
      ParameterSpace   fParamSpace;
      ParameterManager fParamManager;

      std::array<std::array<Int_t, NUM_SEAVEYS>, AnitaPol::kNotAPol> fNormalization;
      std::array<std::array<Double_t, NUM_SEAVEYS>, AnitaPol::kNotAPol> fDdts;

      bool fApplyParams = true;
      bool fSaveResults = false;

      
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

      typedef std::map<std::pair<AnitaPol::AnitaPol_t, PhysicalRing>, double> PolRingMap;
      PolRingMap fFittedR; /// pol/ant to radial (m)
      PolRingMap fFittedZ; /// pol/ant to z position (m)
      PolRingMap fFittedPhi; /// pol/ant to phi in radians

      inline double get(const PolRingMap& fittedVals, const PolAntMap& defaultVals,
			AnitaPol::AnitaPol_t pol, Int_t ant) const {      
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

std::ostream& operator<<(std::ostream& os, const Acclaim::PhaseCenterFit::PhysicalRing& r);

#endif 
