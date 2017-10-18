#include "SumTreeReductionSelector.h"
#include "TProofOutputFile.h"

/** 
 * Default constructor
 * 
 * @param outFileName is the name to give the file containing the combined trees
 */
Acclaim::SumTreeReductionSelector::SumTreeReductionSelector(const char* outFileName, const char* reducedSumTreeName)
  : fOutSum(NULL), fOutTree(NULL), fProofOutFile(NULL), fOutFileName("fOutFileName", outFileName), fReducedSumTreeName("fReducedSumTreeName",  reducedSumTreeName)
{
  
}

void Acclaim::SumTreeReductionSelector::Begin(TTree* tree){
  SummarySelector::Begin(tree);

  // Does this mean the slaves can find this?
  fInput->Add(&fOutFileName);
  fInput->Add(&fReducedSumTreeName);
}


void Acclaim::SumTreeReductionSelector::SlaveBegin(TTree* tree){
  SummarySelector::SlaveBegin(tree);

  fOutFileName = *(dynamic_cast<TNamed*>(fInput->FindObject("fOutFileName")));
  fReducedSumTreeName = *(dynamic_cast<TNamed*>(fInput->FindObject("fReducedSumTreeName")));
  
  // https://root-forum.cern.ch/t/proof-tree-merging/8097/2
  // https://root.cern.ch/handling-large-outputs-root-files
  fProofOutFile = new TProofOutputFile(fOutFileName.GetTitle(), TProofOutputFile::kMerge);
  GetOutputList()->Add(fProofOutFile);

  fOut = fProofOutFile->OpenFile("recreate");
  fOutTree = new TTree("sumTree", "sumTree");
  fOutTree->Branch("sum", &fOutSum);
}


Bool_t Acclaim::SumTreeReductionSelector::Process(Long64_t entry){

  Bool_t matchSelection = SummarySelector::Process(entry);
  if(matchSelection){
    *fOutSum = *fSum;
    fOutTree->Fill();
  }
  return matchSelection;
}

void Acclaim::SumTreeReductionSelector::SlaveTerminate(){
  SummarySelector::SlaveTerminate();

  fOut->Write();
  fOut->Close();
  fOut = NULL;
  fOutTree = NULL; // should have been written to disk
}  

void Acclaim::SumTreeReductionSelector::Terminate(){
  SummarySelector::Terminate();

  TList* l = GetOutputList();
  fProofOutFile = dynamic_cast<TProofOutputFile*>(l->FindObject(fOutFileName.GetTitle()));
  if(fProofOutFile){
    TFile* f = fProofOutFile->OpenFile("read");
    TTree* t = dynamic_cast<TTree*>(f->Get(fReducedSumTreeName.GetTitle()));
    if(t){
      std::cout << "Created " << t->GetName() << " in file " << f->GetName() <<  " has "
		<< t->GetEntries() << " entries..." << std::endl;
      t->Print();
    }
    f->Close();
    // write to file, maybe...
  }
}
