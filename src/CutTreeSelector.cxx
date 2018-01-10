#include "CutTreeSelector.h"
#include "TProofOutputFile.h"
#include "TTreeFormula.h"
#include "DrawStrings.h"
#include "CutOptimizer.h"
#include "TList.h"


/**
 * Default constructor
 *

 * @param treeName is the name to give the output ttree
 */
Acclaim::CutTreeSelector::CutTreeSelector(const char* outFileName, const char* treeName)
  : fOutTree(NULL), fProofOutFile(NULL), fOutFileName("fOutFileName", outFileName), fTreeName("fTreeName", treeName),
    fFormulaStrings(new TList), fFormulas(NULL), fIterationFormula()
{

  fFormulaStrings->SetName("fFormulaStrings");
  fDoDemoHist = false;
}


/**
 * Not sure if this is necessary...
 * @return true
 */
Bool_t Acclaim::CutTreeSelector::Notify()  
{

  SummarySelector::Notify();
  TIter next(fFormulas);
  while(TTreeFormula* form = dynamic_cast<TTreeFormula*>(next())){
    form->Notify();
  }
  
  return kTRUE;
}


/**
 * @brief Set the formula strings to evaluate and put into the created output tree.
 *
 * Puts c-style strings into a TList of TObjStrings, which can be passed to slaves via fInput
 *
 * @param formulaStrings are the formula strings to evaluate inside Process()
 */
void  Acclaim::CutTreeSelector::setFormulaStrings(const std::vector<const char*>& formulaStrings)
{

  fFormulaStrings->Delete(); // clear and delete?

  std::vector<const char*>::const_iterator it = formulaStrings.begin();
  for(; it != formulaStrings.end(); ++it){
    fFormulaStrings->Add(new TObjString(*it));
  }
  fFormulaStrings->SetOwner();
}


/**
 * Add input here, the slave processes can only access things in fInput
 */
void Acclaim::CutTreeSelector::Begin(TTree* tree){
  SummarySelector::Begin(tree); // get the selection into memory also!
  fInput->Add(&fOutFileName);
  fInput->Add(&fTreeName);
  fInput->Add(fFormulaStrings);
}


/**
 * Called by all slave processes. Read fInput here, init fOutput here
 */
void Acclaim::CutTreeSelector::SlaveBegin(TTree* tree){

  fOutFileName = *(dynamic_cast<TNamed*>(fInput->FindObject("fOutFileName")));
  fTreeName = *(dynamic_cast<TNamed*>(fInput->FindObject("fTreeName")));
  fFormulaStrings = dynamic_cast<TList*>(fInput->FindObject("fFormulaStrings"));

  // https://root-forum.cern.ch/t/proof-tree-merging/8097/2
  // https://root.cern.ch/handling-large-outputs-root-files
  fProofOutFile = new TProofOutputFile(fOutFileName.GetTitle(), TProofOutputFile::kMerge);
  fOutput->Add(fProofOutFile);

  fOut = fProofOutFile->OpenFile("recreate");
  fOutTree = new TTree(fTreeName.GetTitle(), fTreeName.GetTitle());

  fFormulas = new TList();
  fFormulas->SetOwner(true);

  const size_t nForm = fFormulaStrings->GetEntries();
  fIntVals.resize(nForm, 0);
  fFloatVals.resize(nForm, 0);
  fFormulaReturnTypes.resize(nForm, 0);

  SummarySelector::SlaveBegin(tree); // Optional analysisCutTree needs to be booked after TProofOutputFile
}



/**
 * Allocate branches on first pass, reassign TTree formulas on each call
 *
 * @param tree is the ttree which is about to be processed
 */
void Acclaim::CutTreeSelector::Init(TTree* tree){

  SummarySelector::Init(tree);  

  fFormulas->SetOwner(true);
  fFormulas->Delete();
  
  if(fFormulas->GetEntries()==0){

    const size_t nForm = fFormulaStrings->GetEntries();
    TObjArray* fOutBranchList = fOutTree->GetListOfBranches();
    fIterationFormula.resize(nForm, false);

    for(UInt_t fInd=0; fInd < nForm; fInd++){
      const TString& formula = dynamic_cast<TObjString*>(fFormulaStrings->At(fInd))->String();

      TString formName = TString::Format("form%u", fInd);
      TTreeFormula* f = new TTreeFormula(formName, formula, tree);
      fFormulas->Add(f);

      // std::cout << fInd << "\t" << formula << "\t" << f->GetNdata() << std::endl;

      if(formula.Contains("Iteration$")){
	fIterationFormula[fInd] = true;
	// std::cout << "the iteration formula index is " << fIterationFormula << std::endl;
      }

      if(fOutBranchList->GetEntries() < (Int_t)nForm){
	TString bName = CutOptimizer::branchifyName(formula);

	fFormulaReturnTypes.at(fInd) = f->IsInteger();
	if(fFormulaReturnTypes.at(fInd) > 0){
	  fOutTree->Branch(bName, &fIntVals.at(fInd));
	}
	else{
	  fOutTree->Branch(bName, &fFloatVals.at(fInd));
	}
      }
    }

    // if(fIterationFormula==-1){
    //   std::cout << "the iteration formula index is " << fIterationFormula << std::endl;
    // }
  }
  // else{
  //   std::cout << tree->GetName() << std::endl;
  //   tree->Show(0);
  //   TIter next(fFormulas);
  //   while(TTreeFormula* form = dynamic_cast<TTreeFormula*>(next())){
  //     form->SetTree(tree);
  //     // form->UpdateFormulaLeaves();
  //   }
  // }
}




/**
 * Main event loop
 *
 * @param entry is the local tree entry
 *
 * @return always returns true
 */
Bool_t Acclaim::CutTreeSelector::Process(Long64_t entry){

  Bool_t matchesSelection = SummarySelector::Process(entry);
  if(matchesSelection){
    int fInd=0;
    TIter next(fFormulas);
    while(TObject* obj = next()){
      TTreeFormula* f = dynamic_cast<TTreeFormula*>(obj);

      fIntVals.at(fInd) = -9999;
      fFloatVals.at(fInd) = -9999;

      bool doneAnIteration = false;

      //***************************************************************
      // this necessary to force the TTree formula to load the data!!!!
      f->EvalInstance(); 
      //***************************************************************

      for(int i=0; i < fMaxNdata; i++){
	if(!doneAnIteration && fCumulativeCutReturns.at(i) > 0){

	  bool requireIteration = fIterationFormula[fInd] || (f->GetNdata() > 1 && i < f->GetNdata());
	  // std::cout << i << "\t" << fInd << "\t" << fCumulativeCutReturns.at(i) << "\t" << f->GetTitle() << "\t";

	  if(fFormulaReturnTypes.at(fInd) > 0){
	    fIntVals.at(fInd) = requireIteration ? f->EvalInstance(i) : f->EvalInstance();
	    // std::cout << "(int) = " << fIntVals.at(fInd) << "\t" << requireIteration << std::endl;
	  }
	  else{
	    fFloatVals.at(fInd) = requireIteration ? f->EvalInstance(i) : f->EvalInstance();
	    // std::cout << "(flt) = " << fFloatVals.at(fInd) << "\t" << requireIteration << std::endl;
	  }
	  doneAnIteration = true;
	}
      }
      // std::cout << fInd << "\t" << f->GetTitle() << "\t" << "the return vals are..." << std::endl;
      fInd++;
    }
    // std::cout << "will fill those values... " << std::endl;
    // std::cout << std::endl;
    
    fOutTree->Fill();
  }

  return kTRUE;
}



/**
 * Deallocate slave heap objects
 *
 */
void Acclaim::CutTreeSelector::SlaveTerminate(){  

  SummarySelector::SlaveTerminate();

  // fOutTree->Write();
  // delete fOutTree;
  // fOutTree = NULL;
  
  fOut->Write();
  fOut->Close();
  fOut = NULL;
  
  // fInput->Remove(fFormulas);
  // fFormulas->Delete();
  // delete fFormulas;
  // fFormulas = NULL;

  fFormulaReturnTypes.clear();
  fFloatVals.clear();
  fIntVals.clear();
}


/**
 * Deallocate master heap objects
 */
void Acclaim::CutTreeSelector::Terminate(){

  // fInput->Remove(fFormulaStrings);  
  // fFormulaStrings->Delete();
  // delete fFormulaStrings;
  // fFormulaStrings = NULL;  

  fProofOutFile = dynamic_cast<TProofOutputFile*>(fOutput->FindObject(fOutFileName.GetTitle()));
  // std::cout << fProofOutFile << std::endl;
  if(fProofOutFile){
    TFile* f = fProofOutFile->OpenFile("read");
    // std::cout << f << std::endl;    
    if(f){
      TTree* t = dynamic_cast<TTree*>(f->Get(fTreeName.GetTitle()));
      // std::cout << t << std::endl;          
      if(t){
	std::cout << "Created " << t->GetName() << " in file " << f->GetName() <<  " has "
		  << t->GetEntries() << " entries..." << std::endl;
	// t->Print();
      }
      f->Close();
    }
  }
}
