#include "AnitaEventSummary.h"
#include "TTree.h"
#include "TFile.h"
#include "SummarySet.h"
#include "AnitaDataset.h"
#include "ProgressBar.h"
#include "AnalysisReco.h"
#include "FilterStrategy.h"
#include "InterferometricMap.h"

using namespace Acclaim;

/** 
 * Script to make a copy of summary trees, and hackily edit a couple of values if possible
 */
int main(int argc, char* argv[]){

  if(argc!=4){
    std::cerr << "Usage: " <<  argv[0] << " inputSumTree inputFlagTree run" << std::endl;
    return 1;
  }
 
  TString inputSumTree = argv[1];
  TString inputFlagTree = argv[2];
  int run = atoi(argv[3]);

  TString inGlob = inputSumTree;
  SummarySet ss(inGlob);
  ss.addFlagChain(inputFlagTree);

  const Long64_t N = ss.N();
  if(N > 0){    
    TString outFileName = TString::Format("doSineSub_all_%d.root", run);
    TFile* fOut = new TFile(outFileName, "recreate");
    TTree* sumTree = new TTree("sumTree", "sumTree");

    TNamed* settings = (TNamed*)(ss.getChain()->GetFile())->Get("AcclaimSettings.conf");
    settings->Write();
    delete settings;
    settings = NULL;
    
    AnitaEventSummary* sum = NULL;
    sumTree->Branch("sum", &sum);

    ProgressBar p(N);
    for(Long64_t entry=0; entry < N; entry++){
      ss.getEntry(entry);

      AnitaEventSummary* inSum = ss.summary();
      *sum = *inSum;
      
      sumTree->Fill();
      p.inc(entry, N);
    }
    fOut->Write();
    fOut->Close();
  }
  
  return 0;
}
