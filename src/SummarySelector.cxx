#include "SummarySelector.h"

#if ROOT_VERSION_CODE >= SUMMARY_SELECTOR_MIN_ROOT_VERSION

#include <TH2.h>
#include <TStyle.h>
#include "TCanvas.h"
#include "AnalysisCuts.h"

ClassImp(Acclaim::SummarySelector);


/** 
 * Default constructor
 * 
 * @param tree Default parameter which is not used, for backward compatilbility reasons
 * @param sumBranchName the name of the branch of AnitaEventSummaries, default is "sum"
 */
Acclaim::SummarySelector::SummarySelector(const char* sumBranchName)
: fReader(), fChain(NULL), fSumReaderValue({fReader, sumBranchName}),
  fSum(NULL), fEventSelection(new TList), fSummarySelectorDemoHist(NULL),
  fDoSummarySelectorDemoHist(false) {

  fEventSelection->SetName("fEventSelection");
}



/** 
 * Destructor
 */
Acclaim::SummarySelector::~SummarySelector()
{

}



void  Acclaim::SummarySelector::setEventSelection(const std::vector<const AnalysisCuts::AnalysisCut*>& eventSelection)
{

  fEventSelection->Clear();


  // const TList?
  std::vector<const AnalysisCuts::AnalysisCut*>::const_iterator it = eventSelection.begin();
  for(; it != eventSelection.end(); ++it){
    AnalysisCuts::AnalysisCut* c = const_cast<AnalysisCuts::AnalysisCut*>(*it);
    fEventSelection->Add(c);
  }
}




/** 
 * The Init() function is called when the selector needs to initialize
 * a new tree or chain. Typically here the reader is initialized.
 * It is normally not necessary to make changes to the generated
 * code, but the routine can be extended by the user if needed.
 * Init() will be called many times when running on PROOF
 * (once per file to be processed).
 * 
 * @param tree the tree to initialize
 */
void Acclaim::SummarySelector::Init(TTree *tree)
{
  fReader.SetTree(tree);
}



Bool_t Acclaim::SummarySelector::Notify()
{
  
  return kTRUE;
}


/** 
 * The Begin() function is called at the start of the query.
 * When running with PROOF Begin() is only called on the client.
 * The tree argument is deprecated (on PROOF 0 is passed).
 */
void Acclaim::SummarySelector::Begin(TTree * /*tree*/)
{
  // TString option = GetOption();
  fInput->Add(fEventSelection);
}


/** 
 * The SlaveBegin() function is called after the Begin() function.
 * When running with PROOF SlaveBegin() is called on each slave server.
 * The tree argument is deprecated (on PROOF 0 is passed).
 */
void Acclaim::SummarySelector::SlaveBegin(TTree * /*tree*/)
{
  fEventSelection = dynamic_cast<TList*>(fInput->FindObject("fEventSelection")); 

  if(fDoSummarySelectorDemoHist){
    fSummarySelectorDemoHist = new TH1D("hDemo", "SummarySelector demo histogram (peak[1][0].value)", 1024, 0, 1);
    fOutput->Add(fSummarySelectorDemoHist);
  }
}


/** 
 * @brief Reads the AnitaEventSummary TTree entry and sets the fSum pointer.
 * Cycles through the fEventSelection, applying each AnalysisCut in turn.
 * 
 * The Process() function is called for each entry in the tree (or possibly
 * keyed object in the case of PROOF) to be processed. The entry argument
 * specifies which entry in the currently loaded tree is to be processed.
 * When processing keyed objects with PROOF, the object is already loaded
 * and is available via the fObject pointer.
 * 
 * This function should contain the \"body\" of the analysis. It can contain
 * simple or elaborate selection criteria, run algorithms on the data
 * of the event and typically fill histograms.
 * 
 * The processing can be stopped by calling Abort().
 * 
 * Use fStatus to set the return value of TTree::Process().
 * 
 * @param entry in the currently loaded tree
 * 
 * @return true if the entry passes all the fEventSelection values
 */
Bool_t Acclaim::SummarySelector::Process(Long64_t entry)
{

  fReader.SetLocalEntry(entry);
  fSum = fSumReaderValue.Get();

  if(fSummarySelectorDemoHist){
    fSummarySelectorDemoHist->Fill(fSum->peak[AnitaPol::kVertical][0].value);
  }

  Bool_t matchesSelection = true;
  TIter next(fEventSelection);
  while (TObject* obj = next()){
    const AnalysisCuts::AnalysisCut* eventSelection = dynamic_cast<const AnalysisCuts::AnalysisCut*>(obj);

    Int_t retVal = eventSelection->apply(fSum);
    matchesSelection = matchesSelection && retVal > 0;

    if(!matchesSelection){
      break;
    }    
  }
  
  return matchesSelection;
}


/** 
 * The SlaveTerminate() function is called after all entries or objects
 * have been processed. When running with PROOF SlaveTerminate() is called
 * on each slave server.
 */
void Acclaim::SummarySelector::SlaveTerminate()
{
  fEventSelection->Clear(); // don't delete globals
  delete fEventSelection;
  fEventSelection = NULL;
}



/** 
 * The Terminate() function is the last function to be called during
 * a query. It always runs on the client, it can be used to present
 * the results graphically or save the results to file.
 */
void Acclaim::SummarySelector::Terminate()
{

  if(fDoSummarySelectorDemoHist){
    fSummarySelectorDemoHist = dynamic_cast<TH1D *>(fOutput->FindObject("fSummarySelectorDemoHist"));  
    if (fSummarySelectorDemoHist) {
      TCanvas *c1 = new TCanvas();
      fSummarySelectorDemoHist->Draw("h");

      // Final update
      c1->cd();
      c1->Update();
    }
    else {
      Warning("Terminate", "histogram not found");
    }
  }  
}

#endif
