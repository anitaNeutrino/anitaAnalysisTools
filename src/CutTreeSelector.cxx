#include "CutTreeSelector.h"
#include "TProofOutputFile.h"
#include "TTreeFormula.h"
#include "AnalysisCuts.h"
#include "CutOptimizer.h"
#include "TList.h"


/**
 * Default constructor
 *
 * @param outFileName is the name to give the file containing the combined trees
 * @param treeName is the name to give the output ttree
 */
Acclaim::CutTreeSelector::CutTreeSelector(const char* outFileName, const char* treeName)
  : fOutTree(NULL), fProofOutFile(NULL), fOutFileName("fOutFileName", outFileName), fTreeName("fTreeName", treeName),
    fFormulaStrings(new TList), fFormulas(NULL)
{

  fFormulaStrings->SetName("fFormulaStrings");

}


/**
 * Not sure if this is necessary...
 * @return true
 */
Bool_t Acclaim::CutTreeSelector::Notify()
{
  fFormulas->Notify();
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
}


/**
 * Add input here, the slave processes can only access things in fInput
 */
void Acclaim::CutTreeSelector::Begin(TTree* ){
  fInput->Add(&fOutFileName);
  fInput->Add(&fTreeName);
  fInput->Add(fFormulaStrings);
}


/**
 * Called by all slave processes. Read fInput here, init fOutput here
 */
void Acclaim::CutTreeSelector::SlaveBegin(TTree* ){

  fOutFileName = *(dynamic_cast<TNamed*>(fInput->FindObject("fOutFileName")));
  fTreeName = *(dynamic_cast<TNamed*>(fInput->FindObject("fTreeName")));
  fFormulaStrings = dynamic_cast<TList*>(fInput->FindObject("fFormulaStrings"));

  // https://root-forum.cern.ch/t/proof-tree-merging/8097/2
  // https://root.cern.ch/handling-large-outputs-root-files
  fProofOutFile = new TProofOutputFile(fOutFileName.GetTitle(), TProofOutputFile::kMerge);
  fOutput->Add(fProofOutFile);

  fOut = fProofOutFile->OpenFile("recreate");
  fOutTree = new TTree(fTreeName.GetTitle(), fTreeName.GetTitle());

  if(fFormulas){
    fFormulas->Delete();
    delete fFormulas;
  }
  fFormulas = new TList();

  const size_t nForm = fFormulaStrings->GetEntries();
  fIntVals.resize(nForm, 0);
  fFloatVals.resize(nForm, 0);
  fFormulaReturnTypes.resize(nForm, 0);
}



/**
 * Allocate branches on first pass, reassign TTree formulas on each call
 *
 * @param tree is the ttree which is about to be processed
 */
void Acclaim::CutTreeSelector::Init(TTree* tree){

  SummarySelector::Init(tree);

  fFormulas->Delete();

  const size_t nForm = fFormulaStrings->GetEntries();
  TObjArray* fOutBranchList = fOutTree->GetListOfBranches();

  for(UInt_t i=0; i < nForm; i++){
    const TString& formula = dynamic_cast<TObjString*>(fFormulaStrings->At(i))->String();

    TString formName = TString::Format("form%u", i);
    TTreeFormula* f = new TTreeFormula(formName, formula, tree);
    fFormulas->Add(f);


    if(fOutBranchList->GetEntries() < (Int_t)nForm){
      TString bName = CutOptimizer::branchifyName(formula);

      fFormulaReturnTypes.at(i) = f->IsInteger();
      if(fFormulaReturnTypes.at(i) > 0){
	fOutTree->Branch(bName, &fIntVals.at(i));
      }
      else{
	fOutTree->Branch(bName, &fFloatVals.at(i));
      }
    }

  }
}




/**
 * Main event loop
 *
 * @param entry is the local tree entry
 *
 * @return always returns true
 */
Bool_t Acclaim::CutTreeSelector::Process(Long64_t entry){

  SummarySelector::Process(entry);

  Bool_t matchesSelection = true;
  for(UInt_t i=0; i < fEventSelection.size(); i++){
    Int_t retVal = fEventSelection[i]->apply(fSum);
    matchesSelection = matchesSelection && retVal > 0;
    if(!matchesSelection){
      break;
    }
  }

  if(matchesSelection){
    int i=0;
    TIter next(fFormulas);
    while(TObject* obj = next()){
      TTreeFormula* f = dynamic_cast<TTreeFormula*>(obj);
      if(fFormulaReturnTypes.at(i) > 0){
	fIntVals.at(i) = f->EvalInstance();
      }
      else{
	fFloatVals.at(i) = f->EvalInstance();
      }
      i++;
    }

    fOutTree->Fill();
  }

  return kTRUE;
}



/**
 * Deallocate slave heap objects
 *
 */
void Acclaim::CutTreeSelector::SlaveTerminate(){

  fOut->Write();
  fOut->Close();
  fOut = NULL;
  fOutTree = NULL; // should have been written to disk

  fFormulaStrings->Delete();
  delete fFormulaStrings;
  fFormulaStrings = NULL;

  fFormulas->Delete();
  delete fFormulas;
  fFormulas = NULL;

  fFormulaReturnTypes.clear();
  fFloatVals.clear();
  fIntVals.clear();
}


/**
 * Deallocate master heap objects
 */
void Acclaim::CutTreeSelector::Terminate(){

  TList* l = GetOutputList();
  fProofOutFile = dynamic_cast<TProofOutputFile*>(l->FindObject(fOutFileName.GetTitle()));
  if(fProofOutFile){
    TFile* f = fProofOutFile->OpenFile("read");
    TTree* t = dynamic_cast<TTree*>(f->Get(fTreeName.GetTitle()));
    if(t){
      std::cout << "Created " << t->GetName() << " in file " << f->GetName() <<  " has "
		<< t->GetEntries() << " entries..." << std::endl;
      t->Print();
    }
    f->Close();
  }
}
