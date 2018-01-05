#ifndef ACCLAIM_THERMAL_TREE_MAKER_H
#define ACCLAIM_THERMAL_TREE_MAKER_H

#include "TCut.h"

namespace Acclaim {

  /**
   * @class ThermalTreeMaker
   * @brief Reduces the bloated AnitaEventSummaries down to a lean, mean, TMVA-friendly tree.
   */

  class ThermalTreeMaker {
  public:    
    ThermalTreeMaker(const char* glob, const char* outfileName);
    void produce(const std::vector<const char*>& formulas, const std::vector<const TCut*> cuts);

  private:
    TString fSummaryGlob;
    TString fFileName;
    
  };



}







#endif
