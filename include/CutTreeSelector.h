#ifndef CUT_TREE_SELECTOR_H
#define CUT_TREE_SELECTOR_H

#include "SummarySelector.h"

class TList;
class TProofOutputFile;
class TTreeFormula;

namespace Acclaim {

  namespace AnalysisCuts{
    class AnalysisCut;
  }


  /**
   * @class CutTreeSelector
   * @brief A TSelector to parallelize the creation of trees for CutOptimizer
   */
  class CutTreeSelector : public SummarySelector {
  public:

    TTree* fOutTree;						/// TTree to create
    TProofOutputFile* fProofOutFile;				/// The proof output file
    TFile* fOut;						/// Created by the fProofOutFile
    TNamed fOutFileName;					/// The output file name (stored in TNamed.fTitle)
    TNamed fTreeName;						/// The output tree name (stored in TNamed.fTitle)

    std::vector<Int_t> fFormulaReturnTypes;			/// Return type of the trees
    std::vector<Float_t> fFloatVals;				/// Where the formula results are written in the case of a float-like variable
    std::vector<Int_t> fIntVals;				/// Where the formula results are written in the case of a int-like variable
    TList* fFormulaStrings;					/// Internal storage of formula strings, set these with setFormulaStrings()
    TList* fFormulas;						/// List of TTreeFormulas created from fFormulaStrings

    CutTreeSelector(const char* outFileName="CutTreeSelector.root", const char* reducedSumTreeName = "cutTree");
    void setFormulaStrings(const std::vector<const char*>& formulaStrings);

    virtual void   Begin(TTree *tree);
    virtual void   SlaveBegin(TTree *tree);
    virtual void   Init(TTree* tree);
    virtual Bool_t Process(Long64_t entry);
    virtual Bool_t Notify();
    virtual void   SlaveTerminate();
    virtual void   Terminate();    

    ClassDef(CutTreeSelector, 0);
  };
  
}


#endif
