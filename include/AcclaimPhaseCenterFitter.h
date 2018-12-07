#ifndef ACCLAIM_PHASE_CENTER_FITTER_H
#define ACCLAIM_PHASE_CENTER_FITTER_H

#include "AcclaimCorrelationSummary.h"
#include <map>
#include "Math/Functor.h"

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

    enum class ParameterSpace {PitchRoll};

    void SetFitParameterSpace(ParameterSpace ps);
    const std::vector<double>& fit();
    void printResults() const;
    
    
  private:
    void readInSummaries();
    void makeFunctors();
    double eval(const double* params);

    const double fCorrelationThreshold = 0.4;
    TChain* fChain = nullptr;
    CorrelationSummary* fSum = nullptr;
    std::vector<Acclaim::CorrelationSummary> fSummaries;
    ROOT::Math::Minimizer* fMin = nullptr;
    std::vector<double> fResults;
    std::map<ParameterSpace, ROOT::Math::Functor> fFuncs;
    ParameterSpace fParamSpace;
    
  };
}

#endif 
