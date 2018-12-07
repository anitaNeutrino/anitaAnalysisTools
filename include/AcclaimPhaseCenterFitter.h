#ifndef ACCLAIM_PHASE_CENTER_FITTER_H
#define ACCLAIM_PHASE_CENTER_FITTER_H

#include "AcclaimCorrelationSummary.h"

class TChain;

namespace ROOT {
  namespace Math {
    class Minimizer;
    class Functor;
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

    void fit();

  private:
    void readInSummaries();
    double evalPitchRoll(const double* params);

    const double fCorrelationThreshold = 0.4;
    TChain* fChain = nullptr;
    CorrelationSummary* fSum = nullptr;
    std::vector<Acclaim::CorrelationSummary> fSummaries;
    ROOT::Math::Minimizer* fMin = nullptr;
    ROOT::Math::Functor* fFunc = nullptr;
    
  };
}

#endif 
