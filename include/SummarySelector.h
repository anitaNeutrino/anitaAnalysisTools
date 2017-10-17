#ifndef SummarySelector_h
#define SummarySelector_h

#include <TSelector.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "AnitaEventSummary.h"
#include "AnalysisCuts.h"

// Use the TTree reader?
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
#define USE_TTREE_READER
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#endif

class TH1D;
class TBranch;

namespace Acclaim {

  // namespace AnalysisCuts {
  //   class AnalysisCut;
  // }

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

    /** Useful in derived classes */
    TTree*				 fChain;			/// The analyzed TTree or TChain
#if USE_TTREE_READER
    TTreeReader				 fReader;			/// The tree reader
    TTreeReaderValue<AnitaEventSummary>	 fSumReaderValue;		/// The TTree reader value
#else
    TString fSumBranchName;
    TBranch* fSumBranch; //!
#endif
    AnitaEventSummary*			 fSum;				/// AnitaEventSummary loaded with tree entry by GetEntry(entry)
    TList*				 fEventSelection;		/// A list of AnalysisCut::AnalysisCuts objects, for event selection (none means selecting all)

    /** For demonstration */
    TH1D*				 fSummarySelectorDemoHist;	/// A histogram of peak[1][0].value
    bool				 fDoSummarySelectorDemoHist;	/// Set this to true to generate and fill fSummarySelectorDemoHist

    SummarySelector(const char* sumBranchName = "sum");
    virtual ~SummarySelector();

    void addEventSelectionCut(const Acclaim::AnalysisCuts::AnalysisCut* analysisCut);
    
    virtual void   Begin(TTree *tree);
    virtual void   SlaveBegin(TTree *tree);
    virtual void   Init(TTree *tree);
    virtual Bool_t Notify();
    virtual Bool_t Process(Long64_t entry);
    virtual void   SlaveTerminate();
    virtual void   Terminate();
    
    virtual Int_t  Version() const					/// From ROOT
    {
      return 2;
    } 
    virtual Int_t  GetEntry(Long64_t entry, Int_t getall = 0)		/// From ROOT, Gets the local tree entry
    { 
      return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0;
    } 
    virtual void   SetOption(const char *option)			/// From ROOT, set the option
    {
      fOption = option;
    } 
    virtual void   SetObject(TObject *obj)				/// Set the current object, not sure what this does
    {
      fObject = obj;
    }
    virtual void   SetInputList(TList *input)				/// Set the input list, not sure what this does
    {
      fInput = input;
    }
    virtual TList* GetOutputList() const				/// Used to combine the objects
    {
      return fOutput;
    }

    ClassDef(SummarySelector,0);

  };

}

#endif
