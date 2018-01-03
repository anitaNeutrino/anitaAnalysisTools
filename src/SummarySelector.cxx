#include "SummarySelector.h"

#include <TH2.h>
#include <TStyle.h>
#include "TCanvas.h"
#include "SummaryDraw.h"
#include "TTreeFormula.h"
#include "TMath.h"

ClassImp(Acclaim::SummarySelector);


/** 
 * Default constructor
 * 
 * @param tree Default parameter which is not used, for backward compatibility reasons
 * @param sumBranchName the name of the branch of AnitaEventSummaries, default is "sum"
 */
Acclaim::SummarySelector::SummarySelector(const char* sumBranchName)
: fChain(NULL),
#ifdef USE_TTREE_READER
  fReader(), fSumReaderValue(fReader, sumBranchName),
#else
  fSumBranchName(sumBranchName),
  fSumBranch(NULL),
#endif
  fSum(NULL), fEventSelection(new TList), fEventSelectionFormulas(new TList),
  fAnalysisCutTreeName("analysisCutTree"), fAnalysisCutTree(NULL),
  fAnalysisCutReturns(), fDoAnalysisCutTree(true),
  fSummarySelectorDemoHist(NULL), fDoSummarySelectorDemoHist(false),
  fEventNumber(0), fRun(0), fWeight(0)
{
  fEventSelection->SetName("fEventSelection");
}



/** 
 * Destructor
 */
Acclaim::SummarySelector::~SummarySelector()
{

}



/** 
 * Add an analysis cut to fEventSelection, takes care of const casting
 * 
 * @param analysisCut 
 */
void  Acclaim::SummarySelector::addEventSelectionCut(const TCut* analysisCut)
{
  fEventSelection->Add(const_cast<TCut*>(analysisCut));
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
#ifdef USE_TTREE_READER  
  fReader.SetTree(tree);
#else
  fChain = tree;
  fChain->SetBranchAddress(fSumBranchName, &fSum, &fSumBranch);
#endif


  fEventSelectionFormulas->SetOwner();
  fEventSelectionFormulas->Clear();
  fDirectionFormulaIndex = -1;

  const int nForm = fEventSelection->GetEntries();
  for(UInt_t i=0; i < nForm; i++){
    const TCut* eventSelection = dynamic_cast<const TCut*>(fEventSelection->At(i));

    TString formName = TString::Format("form_%s", eventSelection->GetName());
    TTreeFormula* f = new TTreeFormula(formName, eventSelection->GetTitle(), tree);
    fEventSelectionFormulas->Add(f);

    if(TString(eventSelection->GetName())=="highestPeak"){
      fDirectionFormulaIndex = i;
    }
  }

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

  fAnalysisCutReturns.resize(fEventSelection->GetEntries(), 0);

  if(fDoAnalysisCutTree){
    fAnalysisCutTree = new TTree(fAnalysisCutTreeName, fAnalysisCutTreeName);
    fAnalysisCutTree->Branch("eventNumber", &fEventNumber);
    fAnalysisCutTree->Branch("run", &fRun);
    fAnalysisCutTree->Branch("weight", &fWeight);

    TIter next(fEventSelection);
    int i=0;
    while (TObject* obj = next()){
      const TCut* eventSelection = dynamic_cast<const TCut*>(obj);
      TString name = eventSelection->GetName();
      fAnalysisCutTree->Branch(name, &fAnalysisCutReturns.at(i));
      i++;
    }
  }
}


/** 
 * @brief Reads the AnitaEventSummary TTree entry and sets the fSum pointer.
 * Cycles through the fEventSelection, applying each TCut in turn.
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
#ifdef USE_TTREE_READER
  fReader.SetLocalEntry(entry);
  fSum = fSumReaderValue.Get();
#else
  fChain->GetEntry(entry);
#endif  


  Int_t peakIteration = 0;
  TTreeFormula* peakDirectionFormula = dynamic_cast<TTreeFormula*>(fEventSelectionFormulas->At(fDirectionFormulaIndex));//FindObject("form_highestPeak"));
  for(int dirInstance=0; dirInstance < peakDirectionFormula->GetNdata(); dirInstance++){
    float val = peakDirectionFormula->EvalInstance(dirInstance);
    if(val > 0){
      peakIteration = dirInstance;
      // break;
    }
  }

  Bool_t matchesSelection = true;
  TIter next(fEventSelectionFormulas);
  int i=0;
  while (TObject* obj = next()){
    TTreeFormula* cutFormula = dynamic_cast<TTreeFormula*>(obj);

    Float_t cutVal = cutFormula->GetNdata() > 1 ? cutFormula->EvalInstance(peakIteration) : cutFormula->EvalInstance();
    // fAnalysisCutReturns.at(i) = fDirectionFormulaIndex; //cutVal;
    // fAnalysisCutReturns.at(i) = TMath::Nint(cutVal);//peakIteration + 1000*fDirectionFormulaIndex;
    fAnalysisCutReturns.at(i) = cutVal;
    i++;

    matchesSelection = matchesSelection && cutVal > 0;

    // we can break out of the loop early if we're not storing the results
    // of all the cuts, and the selection doesn't match all the cuts
    if(!fDoAnalysisCutTree && !matchesSelection){
      break;
    }
  }

  // if doing the demo hist, then fill event passes the cut
  if(matchesSelection && fSummarySelectorDemoHist){
    fSummarySelectorDemoHist->Fill(fSum->peak[AnitaPol::kVertical][0].value);
  }

  // If storing the analysis cut result, then update
  if(fDoAnalysisCutTree && fAnalysisCutTree){
    fEventNumber = fSum->eventNumber;
    fRun = fSum->run;
    fWeight = fSum->weight();
    fAnalysisCutTree->Fill();
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
  // fInput->Remove(fEventSelection);

  // fEventSelection->Clear(); // don't delete globals
  // delete fEventSelection;
  // fEventSelection = NULL;
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
