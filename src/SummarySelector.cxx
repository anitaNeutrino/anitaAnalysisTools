#include "SummarySelector.h"

#include <TH2.h>
#include <TStyle.h>
#include "TCanvas.h"
#include "SummaryDraw.h"
#include "TTreeFormula.h"
#include "TMath.h"
#include "TTreeFormulaManager.h"

ClassImp(Acclaim::SummarySelector);


/** 
 * Default constructor
 * 
 * @param tree Default parameter which is not used, for backward compatibility reasons
 * @param sumBranchName the name of the branch of AnitaEventSummaries, default is "sum"
 */
Acclaim::SummarySelector::SummarySelector(const char* sumBranchName)
  : TSelector(), fChain(NULL),
// #ifdef USE_TTREE_READER
//     fReader(), fSumReaderValue(fReader, sumBranchName),
// #else
//     fSumBranchName(sumBranchName),
//     fSumBranch(NULL),
// #endif
    // fSum(NULL),
    fCuts(new TList), fCutFormulas(NULL),
    fDemoForm(NULL),
    // fAnalysisCutTreeName("analysisCutTree"), fAnalysisCutTree(NULL),
    // fAnalysisCutReturns(), fDoAnalysisCutTree(true),
    fDemoHist(NULL), fDoDemoHist(false)
    // fEventNumber(0), fRun(0), fWeight(0)
{
  fInput = new TList;
  fCuts->SetName("fCuts");
}



/** 
 * Destructor
 */
Acclaim::SummarySelector::~SummarySelector()
{

}



/** 
 * Add an analysis cut to fCuts, takes care of const casting
 * 
 * @param analysisCut 
 */
void  Acclaim::SummarySelector::addCut(const TCut* analysisCut)
{
  fCuts->Add(const_cast<TCut*>(analysisCut));
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

  fTree = tree;

  if(!fCutFormulas){
    fCutFormulas = new TList;
  }
  if(fCutFormulas->GetEntries()==0){
    const int nForm = fCuts->GetEntries();
    fCutReturns.resize(nForm, std::vector<Int_t>());
    
    for(UInt_t i=0; i < nForm; i++){
      const TCut* eventSelection = dynamic_cast<const TCut*>(fCuts->At(i));

      TString formName = TString::Format("form_%s", eventSelection->GetName());
      TTreeFormula* f = new TTreeFormula(formName, eventSelection->GetTitle(), fTree);
      fCutFormulas->Add(f);
      // fManager->Add(f);
      fCutReturns.at(i).resize(f->GetNdata(), -1);
    }


    if(fDoDemoHist){
      // fDemoForm = new TTreeFormula("demo_form", "peak[1][0].value", fTree);
      fDemoForm = new TTreeFormula("demo_form", "run", fTree);      
    }
    
  }
  else{
    TIter next(fCutFormulas);
    while(TTreeFormula* form = dynamic_cast<TTreeFormula*>(next())){
      form->SetTree(tree);
      form->UpdateFormulaLeaves();
    }

    if(fDemoForm){
      fDemoForm->SetTree(tree);
      fDemoForm->UpdateFormulaLeaves();
    }
  }
}
  







Bool_t Acclaim::SummarySelector::Notify()
{
  TIter next(fCutFormulas);
  while(TTreeFormula* form = dynamic_cast<TTreeFormula*>(next())){
    form->UpdateFormulaLeaves();
  }
  
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
  
  fInput->Add(fCuts);
  
  if(fDoDemoHist){
    fInput->Add(new TNamed("fDoDemo", "hDemo"));
  }
  
  
}


/** 
 * The SlaveBegin() function is called after the Begin() function.
 * When running with PROOF SlaveBegin() is called on each slave server.
 * The tree argument is deprecated (on PROOF 0 is passed).
 */
void Acclaim::SummarySelector::SlaveBegin(TTree * /*tree*/)
{
  fCuts = dynamic_cast<TList*>(fInput->FindObject("fCuts"));

  TNamed* n = dynamic_cast<TNamed*>(fInput->FindObject("fDoDemo"));
  fDoDemoHist = n ? true : false;
  if(fDoDemoHist){
    fDemoHist = new TH1D("hDemo", "Demo histogram - run; run; number of events", 310, 130, 440);
    fOutput->Add(fDemoHist);
  }

  
}


/** 
 * @brief Reads the AnitaEventSummary TTree entry and sets the fSum pointer.
 * Cycles through the fCuts, applying each TCut in turn.
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
 * @return true if the entry passes all the fCuts values
 */
Bool_t Acclaim::SummarySelector::Process(Long64_t entry)
{

  fTree->LoadTree(entry);
  
  // #ifdef USE_TTREE_READER
//   fReader.SetLocalEntry(entry);
//   fSum = fSumReaderValue.Get();
// #else
//   fChain->GetEntry(entry);
// #endif  

  // Int_t peakIteration = 0;
  // TTreeFormula* peakDirectionFormula = dynamic_cast<TTreeFormula*>(fCutFormulas->At(fDirectionFormulaIndex));//FindObject("form_highestPeak"));
  // for(int dirInstance=0; dirInstance < peakDirectionFormula->GetNdata(); dirInstance++){
  //   float val = peakDirectionFormula->EvalInstance(dirInstance);
  //   if(val > 0){
  //     peakIteration = dirInstance;
  //     break;
  //   }
  // }

  Bool_t matchesSelection = true;
  TIter next(fCutFormulas);
  int i=0;
  while (TObject* obj = next()){
    TTreeFormula* form = dynamic_cast<TTreeFormula*>(obj);
    const int peakIteration = 0;
    Float_t cutVal = form->GetNdata() > 1 ? form->EvalInstance(peakIteration) : form->EvalInstance();
    fCutReturns.at(i).at(peakIteration) = cutVal;
    i++;

    matchesSelection = matchesSelection && cutVal > 0;

    // we can break out of the loop early if we're not storing the results
    // of all the cuts, and the selection doesn't match all the cuts
    // if(!fDoAnalysisCutTree && !matchesSelection){
    if(!matchesSelection){
      break;
    }
  }

  // if doing the demo hist, then fill event passes the cut
  if(matchesSelection && fDemoHist){
    fDemoHist->Fill(fDemoForm->EvalInstance());
  }

  // // If storing the analysis cut result, then update
  // if(fDoAnalysisCutTree && fAnalysisCutTree){
  //   fEventNumber = fSum->eventNumber;
  //   fRun = fSum->run;
  //   fWeight = fSum->weight();
  //   fAnalysisCutTree->Fill();
  // }

  return matchesSelection;
}


/** 
 * The SlaveTerminate() function is called after all entries or objects
 * have been processed. When running with PROOF SlaveTerminate() is called
 * on each slave server.
 */
void Acclaim::SummarySelector::SlaveTerminate()
{
  // fInput->Remove(fCuts);

  // fCuts->Clear(); // don't delete globals
  // delete fCuts;
  // fCuts = NULL;
}



/** 
 * The Terminate() function is the last function to be called during
 * a query. It always runs on the client, it can be used to present
 * the results graphically or save the results to file.
 */
void Acclaim::SummarySelector::Terminate()
{

  if(fDoDemoHist){
    fDemoHist = dynamic_cast<TH1D *>(fOutput->FindObject("hDemo"));  
    if (fDemoHist) {
      TCanvas *c1 = new TCanvas();
      fDemoHist->Draw("hist");

      // Final update
      c1->cd();
      c1->Update();
    }
    else {
      Warning("Terminate", "histogram not found");
    }
  }  
}
