#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "AnitaEventSummary.h"

#include <iostream>
#include "ProgressBar.h"

int main(int argc, char* argv[]){

  if(argc != 3){
    std::cerr << argv[0] << "[Input file glob] [Output file name]";
    return 1;
  }
  TString inFileGlob = argv[1];
  TString outFileName = argv[2];

  TString treeName = "sumTree";
  std::cout << inFileGlob << "\t" << outFileName << std::endl;

  TFile* f = new TFile(outFileName, "create");
  if(f==NULL){
    return 2;
  }
  
  TTree* t = new TTree(treeName, treeName);
  AnitaEventSummary* sumOut = NULL;
  t->Branch("sum", &sumOut);
  
  
  TChain* c = new TChain(treeName);
  c->Add(inFileGlob);

  AnitaEventSummary* sumIn = NULL;
  c->SetBranchAddress("sum", &sumIn);
  

  Long64_t n = c->GetEntries();

  Acclaim::ProgressBar p(n);
  for(Long64_t entry=0; entry < n; entry++){
    c->GetEntry(entry);

    sumOut = sumIn;
    t->Fill();

    p.inc(entry, n);
  }

  t->BuildIndex("run", "eventNumber");
  
  f->Write();
  f->Close();

  return 0;
}
