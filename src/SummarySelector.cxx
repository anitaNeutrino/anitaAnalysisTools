#include "SummarySelector.h"

#if ROOT_VERSION_CODE >= SUMMARY_SELECTOR_MIN_ROOT_VERSION

#include <TH2.h>
#include <TStyle.h>
#include "TCanvas.h"

ClassImp(Acclaim::SummarySelector);


/** 
 * Default constructor
 * 
 * @param tree Default parameter which is not used, for backward compatilbility reasons
 * @param sumBranchName the name of the branch of AnitaEventSummaries, default is "sum"
 */
Acclaim::SummarySelector::SummarySelector(TTree* tree, const char* sumBranchName)
: fReader(), fChain(NULL), fSumReaderValue({fReader, sumBranchName}),
  fSum(NULL), fSummarySelectorDemoHist(NULL), fDoSummarySelectorDemoHist(false) {
  (void) tree;
}



/** 
 * Destructor
 */
Acclaim::SummarySelector::~SummarySelector()
{
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
}


/** 
 * The SlaveBegin() function is called after the Begin() function.
 * When running with PROOF SlaveBegin() is called on each slave server.
 * The tree argument is deprecated (on PROOF 0 is passed).
 */
void Acclaim::SummarySelector::SlaveBegin(TTree * /*tree*/)
{
  // TString option = GetOption();

  if(fDoSummarySelectorDemoHist){
    fSummarySelectorDemoHist = new TH1D("hDemo", "SummarySelector demo histogram (peak[1][0].value)", 1024, 0, 1);
    fOutput->Add(fSummarySelectorDemoHist);
  }
}


/** 
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
 * @return The return value is currently not used.
 */
Bool_t Acclaim::SummarySelector::Process(Long64_t entry)
{

  fReader.SetLocalEntry(entry);
  fSum = fSumReaderValue.Get();

  if(fSummarySelectorDemoHist){
    fSummarySelectorDemoHist->Fill(fSum->peak[AnitaPol::kVertical][0].value);
  }
  return kTRUE;
}


/** 
 * The SlaveTerminate() function is called after all entries or objects
 * have been processed. When running with PROOF SlaveTerminate() is called
 * on each slave server.
 */
void Acclaim::SummarySelector::SlaveTerminate()
{
}



/** 
 * The Terminate() function is the last function to be called during
 * a query. It always runs on the client, it can be used to present
 * the results graphically or save the results to file.
 */
void Acclaim::SummarySelector::Terminate()
{

  if(fDoSummarySelectorDemoHist){
    TCanvas *c1 = new TCanvas("c1","Proof ProofEvent canvas",200,10,700,700);
    fSummarySelectorDemoHist = dynamic_cast<TH1D *>(fOutput->FindObject("fSummarySelectorDemoHist"));  
    if (fSummarySelectorDemoHist) {
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
