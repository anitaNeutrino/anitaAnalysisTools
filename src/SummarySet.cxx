#include "SummarySet.h"
#include "TChain.h"
#include "AnitaEventSummary.h"
#include <iostream>

Acclaim::SummarySet::SummarySet(const char* pathToSummaryFiles, const char* treeName, const char* summaryBranchName)
    : fChain(NULL), fSum(NULL) {

  fChain = new TChain(treeName);
  fChain->Add(pathToSummaryFiles);
  fN = fChain->GetEntries();
  fChain->SetBranchAddress(summaryBranchName, &fSum);

  if(fN == 0){
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", no entries in TChain of AnitaEventSummary" << std::endl;
  }

}


Acclaim::SummarySet::~SummarySet(){
  delete fChain;
  fChain = NULL;
  fSum = NULL;
}


Long64_t Acclaim::SummarySet::getEntry(Long64_t entry){
  return fChain->GetEntry(entry);
}
