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

namespace Acclaim{

  /** 
   * @class CutOptimizer
   * @brief A class to read in signal/background samples, feed them into a TMVA framework and separate
   *
   */
class CutOptimizer{

 public:
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


  
  
};

}
#endif //CUT_OPTIMIZER
