
#define SumSelector_cxx
// The class definition in SumSelector.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("SumSelector.C")
// root> T->Process("SumSelector.C","some options")
// root> T->Process("SumSelector.C+")
//


#include "SummarySelector.h"

#if ROOT_VERSION_CODE >= SUMMARY_SELECTOR_MIN_ROOT_VERSION

#include <TH2.h>
#include <TStyle.h>
#include "AnalysisPlot.h"
// #include <TOutputListSelectorDataMap.h>

ClassImp(Acclaim::SummarySelector)

void Acclaim::SummarySelector::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the reader is initialized.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  fReader.SetTree(tree);

}

Bool_t Acclaim::SummarySelector::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}


void Acclaim::SummarySelector::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
  Info(__PRETTY_FUNCTION__, "here");
  TString option = GetOption();
}

void Acclaim::SummarySelector::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

  // const int nInputs = fInput->GetEntries();

  // // clone all the inputs...
  // for(int i=0; i < nInputs; i++){
  //   // For analysis plots
  //   TObject* thisInput= fInput->At(i);
  //   fOutput->Add(thisInput->Clone());
  // }
  // fOutput->Add(new TH1D("h", "h", 1024, 0, 1));
  auto p = new AnalysisPlot("h", "h", 1024, 0, 1);
  p->addCut(&AnalysisCuts::isAboveHorizontal);
  fOutput->Add(p);
  

}

Bool_t Acclaim::SummarySelector::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // When processing keyed objects with PROOF, the object is already loaded
  // and is available via the fObject pointer.
  //
  // This function should contain the \"body\" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.

  fReader.SetLocalEntry(entry);


  AnitaEventSummary* s = sum.Get();
  // AnitaEventSummary* s = sum->Get();

  for(int i=0; i < fOutput->GetEntries(); i++){
    // TH1D* plot = dynamic_cast<TH1D*>(fOutput->At(i));
    AnalysisPlot* plot = dynamic_cast<AnalysisPlot*>(fOutput->At(i));
    if(plot){
      plot->Fill(s, s->highestPeak().value);
    }
    // Acclaim::AnalysisPlot* plot = dynamic_cast<Acclaim::AnalysisPlot*>(fOutput->At(i));
    // if(plot){
    //   plot->Fill(s, s->higherPeak().value);
    // }
    
  }

  // std::cout << s->eventNumber << std::endl;

  return kTRUE;
}

void Acclaim::SummarySelector::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
  // std:: cerr << fOutput->GetEntries() << std::endl;
}

void Acclaim::SummarySelector::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  // std:: cerr << fOutput->GetEntries() << std::endl;

  for(int i=0; i < fOutput->GetEntries(); i++){
    TObject* thisOut = fOutput->At(i);
    // TOutputListSelectorDataMap* silly = dynamic_cast<TOutputListSelectorDataMap*>(thisOut);
    TNamed* named = dynamic_cast<TNamed*>(thisOut);
    if(named && gDirectory){
      named->SetBit(kCanDelete, false);
      gDirectory->Append(thisOut);
    }
  }
  
}

#endif
