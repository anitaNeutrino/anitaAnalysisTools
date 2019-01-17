#ifndef ACCLAIM_PHASE_CENTER_FITTER_H
#define ACCLAIM_PHASE_CENTER_FITTER_H

#include "AcclaimCorrelationSummary.h"
#include <map>
#include "Math/Functor.h"
#include "AnitaConventions.h"
#include "TMath.h"
#include "TVector2.h"
#include "AcclaimPhaseCenterParameters.h"

class AnitaGeomTool;
class TH2D;
class TChain;

namespace ROOT {
  namespace Math {
    class Minimizer;
  }
}

namespace Acclaim {

  namespace PhaseCenter {

    
    /**
     * @class Minimizer
     * @brief Fits the phases centers using ROOT's Minuit minimizer.
     * 
     * Reads the output of make correlation summary
     */

    class Minimizer {
    public:
      Minimizer(const char* corrTreeFiles = nullptr);
    
      void setParameterSpace(ParameterSpace ps);

      const std::vector<double>& fit(AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol);

      void printResults() const;

      void setPrintOnEval(bool print){
	fPrintOnEval = print;
      }

      enum class AllowedPairs {All,
			       SamePhiSector,
			       SamePhiSectorOrHorizontalNeigbour};
      static const char* toCString(AllowedPairs ap);

      void setAllowedPairs(AllowedPairs option) {
	fAllowedPairs = option;
      }

      static bool allowedPair(AllowedPairs ap, int ant1, int ant2);
    
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

      ParameterManager fParamManager;
      ROOT::Math::Functor fFunc;

      std::array<std::array<Int_t, NUM_SEAVEYS>, AnitaPol::kNotAPol> fNormalization;
      std::array<std::array<Double_t, NUM_SEAVEYS>, AnitaPol::kNotAPol> fDdts;

      bool fApplyParams = true;
      bool fSaveResults = false;
      AllowedPairs fAllowedPairs = AllowedPairs::All;

      bool allowedPair(int ant1, int ant2) const;
      template <class T> bool allowedPair(T t) const {return allowedPair(t.ant1, t.ant2);}
      
    };


  
  

  }  
}

#endif 
