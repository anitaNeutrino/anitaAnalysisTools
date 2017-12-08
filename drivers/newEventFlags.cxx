#include "AnitaEventSummary.h"
#include "TTree.h"
#include "TFile.h"
#include "SummarySet.h"
#include "AnitaDataset.h"
#include "ProgressBar.h"
#include "AnalysisReco.h"
#include "FilterStrategy.h"
#include "InterferometricMap.h"
#include "QualityCut.h"

using namespace Acclaim;


/** 
 * Script to make a copy of summary trees, and hackily edit a couple of values if possible
 */
int main(int argc, char* argv[]){

  if(argc!=5 && argc!=4 && argc!=3){
    std::cerr << "Usage: " <<  argv[0] << " inputDir inputGlob" << std::endl;
    std::cerr << "Usage: " <<  argv[0] << " inputDir inputGlob run" << std::endl;
    std::cerr << "Usage: " <<  argv[0] << " inputDir inputGlob run pulser" << std::endl;
    return 1;
  }

  std::vector<UInt_t> eventNumbers {26323378, 27569662, 27897245, 28254843, 28627679, 28776748, 28930237, 28934295, 29006913, 29008144, 32879955, 32953817, 34745811, 34757629, 35020101, 35116217, 35118197, 36710262, 37009416, 37009533, 37009810, 38273261, 38619569, 39464799, 40633613, 40646231, 40996918, 41492058, 41654702, 41799439, 42256877, 42699600, 44675867, 45985415, 47512906, 47831713, 48225288, 48242449, 49313992, 50827460, 51950538, 53275871, 55496464, 55504417, 55861013, 56463512, 56829451, 56894100, 57108465, 57471113, 57519585, 57523926, 57654756, 57839882, 57948519, 58002861, 58008286, 58110716, 58226572, 58310817, 58326495, 58332241, 58421469, 59967256, 60072506, 61104509, 61173194, 61215330, 61238684, 61497702, 61534104, 61784557, 61876418, 62235301, 62235933, 62420819, 62732869, 62849802, 63024583, 63031985, 63484004, 63714468, 63845258, 63845261};
  // std::vector<UInt_t> eventNumbers {60832108};
 
  TString inDir = argv[1];
  TString inFileBaseName = argv[2];
  int run = argc > 3 ? atoi(argv[3]) : -1;
  int type = argc > 4 ? atoi(argv[4]) : 0;

  // TString inGlob = inDir + "/" + inFileBaseName;// + TString::Format("_%d_*.root", run);
  TString inGlob = inDir + "/" + inFileBaseName;
  if(run > -1){
    inGlob += TString::Format("_%d.root", run);
  }
  else{
    inGlob += TString::Format("*.root", run);
  }
  SummarySet ss(inGlob);

  

  const Long64_t N = ss.N();
  if(N > 0){    
    // TString outFileName = inFileBaseName + TString::Format("_%d.root", run);
    // TString outFileName = TString::Format("editSumTree_%d.root", run);

    TString outFileName = inFileBaseName;
    if(run > -1){
      outFileName += TString::Format("_%d_%d.root", run, type);
    }
    else{
      outFileName += TString::Format(".root", run);
    }

    TFile* fOut = new TFile(outFileName, "recreate");
    TTree* flagTree = new TTree("flagTree", "flagTree");
    AnitaEventSummary::EventFlags* flags = NULL;
    UInt_t eventNumber;
    flagTree->Branch("eventNumber", &eventNumber);
    flagTree->Branch("flags", &flags);

    AnalysisSettings settings;
    Acclaim::AnalysisReco reco;
    FilterStrategy doNoFiltering;
    AnitaDataset* d = NULL;
    
    settings.apply(dynamic_cast<TObject*>(&reco));
    settings.write(fOut);

    const Long64_t nLoop = type ? N : eventNumbers.size();

    // ProgressBar p(N);
    // for(Long64_t entry=0; entry < N; entry++){
    //   ss.getEntry(entry);
    
    ProgressBar p(nLoop);

    for(UInt_t entryInd=0; entryInd < nLoop; entryInd++){

      if(type > 0){
	ss.getEntry(entryInd);
      }
      else{
	ss.getEvent(eventNumbers.at(entryInd));
      }

      AnitaEventSummary* inSum = ss.summary();
      
      if(type==0 || (type==1 && inSum->flags.pulser > 0) || (type==2 && inSum->flags.isPayloadBlast > 0)){

	if(d && d->getCurrRun() != inSum->run){
	  delete d;
	  d = NULL;
	}
      
	if(!d){
	  d = new AnitaDataset(inSum->run);
	}
	eventNumber = inSum->eventNumber;

	d->getEvent(inSum->eventNumber);
	FilteredAnitaEvent fEv(d->useful(), &doNoFiltering, d->gps(), d->header());

	Acclaim::QualityCut::applyAll(d->useful(), inSum);

	*flags = inSum->flags;
	reco.fillPowerFlags(&fEv, *flags);

	flagTree->Fill();
      }

      p.inc(entryInd, nLoop);
    }
    fOut->Write();
    fOut->Close();
  }
  
  return 0;
}
