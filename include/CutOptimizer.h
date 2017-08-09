/* -*- C++ -*-.*********************************************************************************************
 Author: Ben Strutt
 Email: strutt@physics.ucla.edu

 Description:
             Read in trees of ANITA event summaries and try to separate signal and background...
***********************************************************************************************************/

#ifndef CUT_OPTIMIZER_H
#define CUT_OPTIMIZER_H


#include "TString.h"
#include "TChain.h"
#include "TXMLEngine.h"
#include "TH1D.h"

namespace Acclaim{

class AnalysisCut;
class SummarySet;

  /** 
   * @class CutOptimizer
   * @brief A class to parse AnitaEventSummary trees, extract quantities of interest,
   *        feed them into a TMVA framework and separate them
   *
   */
class CutOptimizer{

 public:
  static void setDebug(bool db);
  
  CutOptimizer(const char* signalGlob, const char* backgroundGlob = NULL, bool save_trees = false);
  virtual ~CutOptimizer();
  void optimize(const std::vector<const Acclaim::AnalysisCut*>& signalSelection,
                const std::vector<const Acclaim::AnalysisCut*>& backgroundSelection,
                const std::vector<const char*>& formulaStrings,
                const char* outFileName = "");


 protected:

  enum BranchType{
    kUnassigned,
    kInt,
    kFloat
  };

  TFile* makeOutputFile(const char* outFileName);
  void generateSignalAndBackgroundTrees(const std::vector<const Acclaim::AnalysisCut*>& signalSelection,
                                        const std::vector<const Acclaim::AnalysisCut*>& backgroundSelection,
                                        const std::vector<const char*>& treeVars);
  BranchType setBranchFromFormula(TTree* t, const TTreeFormula* f, const char* formulaString, Int_t* intPtr, Float_t* floatPtr);

  TString fSignalGlob;
  TString fBackgroundGlob;
  TString fOutFileName;
  TTree* fSignalTree;
  TTree* fBackgroundTree;
  Bool_t fSaveTrees;
  std::vector<Float_t> fSignalFloatVals;
  std::vector<Float_t> fBackgroundFloatVals;
  std::vector<Int_t> fSignalIntVals;
  std::vector<Int_t> fBackgroundIntVals;




  /**
   * @class Class to get the results of the Fisher Discriminant into a useful form
   */
  class FisherResult {

    typedef std::map<int, double> WeightMap;
    typedef std::map<int, TString> ExpressionMap;
    
   public:
    FisherResult(const TString& fileName){
      getResultFromXML(fileName.Data());
    }
    TH1D* makeHist(int nBinsX, const TString& histName, const TString& histTitle, TTree* t, EColor col) const;
    
   protected:
    void getResultFromXML(const char* filename);
    void parseNode(TXMLEngine* xml, XMLNodePointer_t node, Int_t level);
    
    ExpressionMap fExpressions;
    WeightMap fWeights;
  };




  

  /** 
   * @class Dummy class to hold the TTreeFormula and trick TChain into notifying all the formulas
   */
  class FormulaHolder : public TObject {
   public:
    FormulaHolder(TChain* c);
    virtual ~FormulaHolder();
    virtual Bool_t Notify();
    virtual size_t add(const char* formulaString);
    TTreeFormula* at(UInt_t i) {return fForms.at(i);}
    const char* str(UInt_t i) const {return fFormStrs.at(i);}
    size_t N(){return fForms.size();}
   protected:
    TChain* fChain;
    std::vector<TTreeFormula*> fForms;
    std::vector<const char*> fFormStrs;
  };
  
};

}
#endif //CUT_OPTIMIZER
