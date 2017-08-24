#ifndef SummarySelector_h
#define SummarySelector_h

#include <TSelector.h>

#define SUMMARY_SELECTOR_MIN_ROOT_VERSION ROOT_VERSION(6,0,0)

#if ROOT_VERSION_CODE < SUMMARY_SELECTOR_MIN_ROOT_VERSION
#include <iostream>
namespace Acclaim{
class SummarySelector : public TSelector {
 public:
  SummarySelector(){
    std::cerr << "Use of summary selector requires ROOT version at least " << SUMMARY_SELECTOR_MIN_ROOT_VERSION << std::endl;
  }
};
}  

#else

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include "AnitaEventSummary.h"

namespace Acclaim {


class SummarySelector : public TSelector {
 public :
  TTreeReader     fReader;  //!the tree reader
  TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

  TTreeReaderValue<AnitaEventSummary> sum = {fReader, "sum"}; //!
  // TBranch* sumBranch; //!branch for AnitaEventSummary
  
  
  SummarySelector(TTree * /*tree*/ =0) { ;}
  virtual ~SummarySelector() { }
  virtual Int_t   Version() const { return 2; }
  virtual void    Begin(TTree *tree);
  virtual void    SlaveBegin(TTree *tree);
  virtual void    Init(TTree *tree);
  virtual Bool_t  Notify();
  virtual Bool_t  Process(Long64_t entry);
  virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
  virtual void    SetOption(const char *option) { fOption = option; }
  virtual void    SetObject(TObject *obj) { fObject = obj; }
  virtual void    SetInputList(TList *input) { fInput = input; }
  virtual TList  *GetOutputList() const { return fOutput; }
  virtual void    SlaveTerminate();
  virtual void    Terminate();

  ClassDef(SummarySelector,0);

};

}

#endif
#endif
