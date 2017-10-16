#include "SumTreeReductionSelector.h"
#include "TProofOutputFile.h"

/** 
 * Default constructor
 * 
 * @param outFileName is the name to give the file containing the combined trees
 */
Acclaim::SumTreeReductionSelector::SumTreeReductionSelector(const char* outFileName, const char* reducedSumTreeName)
  : fOutSum(NULL), fOutTree(NULL), fOutFileName(outFileName), fReducedSumTreeName(reducedSumTreeName)
{
  
}


void Acclaim::SumTreeReductionSelector::SlaveBegin(TTree* ){

  // https://root-forum.cern.ch/t/proof-tree-merging/8097/2
  // https://root.cern.ch/handling-large-outputs-root-files
  fProofOutFile = new TProofOutputFile(fOutFileName, TProofOutputFile::kMerge);
  GetOutputList()->Add(fProofOutFile);

  fOut = fProofOutFile->OpenFile("recreate");
  fOutTree = new TTree("sumTree", "sumTree");
  fOutTree->Branch("sum", &fOutSum);

}


Bool_t Acclaim::SumTreeReductionSelector::Process(Long64_t entry){

  SummarySelector::Process(entry);
  
  *fOutSum = *fSum;
  fOutTree->Fill();

  return kTRUE;
}

void Acclaim::SumTreeReductionSelector::SlaveTerminate(){
  fOut->Write();
  fOut->Close();
  fOut = NULL;
  fOutTree = NULL; // should have been written to disk
}  

void Acclaim::SumTreeReductionSelector::Terminate(){

  TList* l = GetOutputList();
  fProofOutFile = dynamic_cast<TProofOutputFile*>(l->FindObject(fOutFileName));
  if(fProofOutFile){
    TFile* f = fProofOutFile->OpenFile("read");
    TTree* t = dynamic_cast<TTree*>(f->Get("sumTree"));
    if(t){
      std::cout << "The tree has " << t->GetEntries() << " entries" << std::endl;
      t->Scan("eventNumber:run:sum.trainingPeak().value");
    }
    // write to file, maybe...
  }
}




