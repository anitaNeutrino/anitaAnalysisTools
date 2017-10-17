#ifndef SUM_TREE_REDUCTION_SELECTOR_H
#define SUM_TREE_REDUCTION_SELECTOR_H

#include "SummarySelector.h"

class TProofOutputFile;

namespace Acclaim {

  /**
   * @class SumTreeReductionSelector
   * @brief An almost trivial proof of concept, for Tree merging with PROOF
   * 
   * May at some point be used to reduce to trees for clustering
   */  

  class SumTreeReductionSelector : public SummarySelector {
  public:
    AnitaEventSummary* fOutSum;
    TTree* fOutTree;
    TProofOutputFile* fProofOutFile;
    TFile* fOut;
    TNamed fOutFileName;
    TNamed fReducedSumTreeName;
    
    SumTreeReductionSelector(const char* outFileName="reduced", const char* reducedSumTreeName = "sumTree");

    // The minimum set required to be useful?
    virtual void   Begin(TTree *tree);    
    virtual void   SlaveBegin(TTree *tree);
    virtual Bool_t Process(Long64_t entry);
    virtual void   SlaveTerminate();
    virtual void   Terminate();

    ClassDef(SumTreeReductionSelector, 0);
  };


  
}




#endif
