#include "AnitaEventSummary.h"
#include "TTree.h"
#include "TFile.h"
#include "SummarySet.h"
#include "AnitaDataset.h"
#include "ProgressBar.h"
#include "AnalysisReco.h"
#include "InterferometricMap.h"
#include "QualityCut.h"

using namespace Acclaim;


/** 
 * Script to make a copy of summary trees, and hackily edit a couple of values if possible
 */
int main(int argc, char* argv[]){

  if(argc!=2){
    std::cerr << "Usage: " <<  argv[0] << " inputFile" << std::endl;
    return 1;
  }

  const bool copyEverything = false;

  TString inputFile = argv[1];
  SummarySet ss(inputFile);

  const Long64_t N = ss.N();

  AnitaDataset* d = NULL;
  
  if(N > 0){    

    std::vector<TString> directories;
    RootTools::tokenize(directories, inputFile, "/");
    TString outFileName = directories.back();
    std::cout << outFileName << std::endl;
    TFile* fOut = new TFile(outFileName, "new");

    TTree* outTree = NULL;
    AnitaEventSummary::EventFlags* flags = NULL;
    AnitaEventSummary* sum = NULL;
    UInt_t eventNumber = 0;
    Int_t run = 0;
    
    if(copyEverything){
      outTree = new TTree("sumTree", "sumTree");
      outTree->Branch("sum", &sum);
      TList* keys = ss.getChain()->GetFile()->GetListOfKeys();
      for(int i=0; i < keys->GetEntries(); i++){
	TString name = keys->At(i)->GetName();
	if(name!="sumTree"){
	  TObject* thing = ss.getChain()->GetFile()->Get(name);
	  thing->Write();
	}
      }
    }
    else{
      outTree = new TTree("flagsTree", "flagsTree");
      outTree->Branch("flags", &flags);
      sum = new AnitaEventSummary();
      outTree->Branch("run", &run);
      outTree->Branch("eventNumber", &eventNumber);      
    }
    
    ProgressBar p(N);
    for(Long64_t entry=0; entry < N; entry++){
      ss.getEntry(entry);
      AnitaEventSummary* inSum = ss.summary();
      *sum = *inSum;

      run = inSum->run;
      if(d && d->getCurrRun()!=run){
	delete d;
	d = NULL;
      }
      if(!d){
	d = new AnitaDataset(run);
      }      


      d->getEvent(inSum->eventNumber);

      UsefulAnitaEvent* useful = d->useful();
      QualityCut::applyAll(useful, sum);

      if(!copyEverything){
	*flags = sum->flags;
	run = sum->run;
	eventNumber = sum->eventNumber;
      }
      
      outTree->Fill();
      p.inc(entry, N);
    }
    fOut->Write();
    fOut->Close();
  }
  
  return 0;
}
