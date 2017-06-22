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

  /** 
   * @class CutOptimizer
   * @brief A class to read in signal/background samples, feed them into a TMVA framework and separate
   *
   */
class CutOptimizer{

 public:
  static void setDebug(bool db);
  
  CutOptimizer(const TString& outFileName, const TString& signalTreeWildCards, const TString& backgroundTree);
  void optimize();

 protected:
  void makeOutputFile();
  void getSignalAndBackgroundTrees();

  
  TString fOutFileName;  
  TString fSignalName;
  TString fBackgroundName;
  TFile* fOutFile;

  TFile* fSignalFile;
  TTree* fSignalTree;

  TFile* fBackgroundFile;
  TTree* fBackgroundTree;


  class FisherResult {

    typedef std::map<int, double> WeightMap;
    typedef std::map<int, TString> ExpressionMap;
    
   public:
    FisherResult(const TString& fileName){
      getResultFromXML(fileName.Data());
    }
    TH1D* makeHist(int nBinsX, const TString& histName, const TString& histTitle, TTree* t, EColor col=kBlack) const;
    
   protected:
    void getResultFromXML(const char* filename);
    void parseNode(TXMLEngine* xml, XMLNodePointer_t node, Int_t level);
    
    ExpressionMap fExpressions;
    WeightMap fWeights;
  };

  
  
};

}
#endif //CUT_OPTIMIZER
