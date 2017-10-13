#ifndef SummarySelector_h
#define SummarySelector_h

#include <TSelector.h>

// For now, is approximate
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
#include "AnitaEventSummary.h"

class TH1D;

namespace Acclaim {

  /**
   * @class SummarySelector
   * @brief Template TSelector to be inherited from, designed for use with SummarySet.
   * 
   * The following methods *must* be defined in the cxx file:
   *    Begin():        Called every time a loop on the tree starts,
   *                    a convenient place to create your histograms.
   * 
   *    SlaveBegin():   Called after Begin(), when on PROOF called only on the
   *                    slave servers.
   * 
   *    Process():      Called for each event, in this function you decide what
   *                    to read and fill your histograms.
   * 
   *    SlaveTerminate: Called at the end of the loop on the tree, when on PROOF
   *                    called only on the slave servers.
   * 
   *    Terminate():    Called at the end of the loop on the tree,
   *                    a convenient place to draw/fit your histograms.
   * 
   * With the Acclaim::SummarySet, create a SummarySelector (or dereived) object.
   * Acclaim::SummarySet::Process(&SummarySelector).
   * Make sure to set SummarySet::SetUseProof(true) to get the full multithreaded PROOF goodness.
   */

  class SummarySelector : public TSelector {
  public :
    TTreeReader     fReader;					/// The tree reader
    TTree          *fChain = 0;					/// The analyzed TTree or TChain
    TTreeReaderValue<AnitaEventSummary> fSumReaderValue;	///The TTree reader value
    AnitaEventSummary* fSum;					/// Summary loaded with tree entry by GetEntry(entry)
    TH1D* fDemoHist;
    bool fDoDemoHist;
  
    SummarySelector(TTree * /*tree*/ =0, const char* sumBranchName = "sum")
      : fSumReaderValue({fReader, sumBranchName}), fDemoHist(NULL), fDoDemoHist(false) { ;}
    
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
    virtual TList*  GetOutputList() const { return fOutput; }
    virtual void    SlaveTerminate();
    virtual void    Terminate();

    ClassDef(SummarySelector,0);

  };

}

#endif
#endif
