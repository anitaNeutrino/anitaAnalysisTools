#include <getopt.h> // for getopt
#include <iostream>

#include <algorithm> //for std::sort 
#include "TString.h"
#include "TFile.h"
#include "TTree.h"

#include "AnitaDataset.h"
#include "RawAnitaHeader.h"
#include "CalibratedAnitaEvent.h"
#include "Adu5Pat.h"


/** 
 * Parse the command line input
 * 
 * @param argc from main
 * @param argv from main
 * @param eventNumbers filled with eventNumbers
 * 
 * @return the run passed on the command lnie
 */
int parseCommandLine(int argc, char* argv[],  std::vector<UInt_t>& eventNumbers);



/** 
 * A program to cut groups of event numbers 
 * from ANITA-3 runs and make a "fake run" out of just those.
 *
 * To make a "sub run" called run 1000 with events 83139414 60832108, do
 * 
 * ANITA_UTIL_INSTALL_DIR/bin/Acclaim/drivers/makeSubRun -r 1000 -e 83139414 60832108 
 */
int main(int argc, char* argv[]){

  Int_t outputRun = 0; ///< The run number given to the shortened run
  std::vector<UInt_t> eventNumbers; ///< The event numbers to trim
  outputRun = parseCommandLine(argc, argv, eventNumbers);  

  // minimize changing runs in case events aren't in order
  std::sort(eventNumbers.begin(), eventNumbers.end());

  TString outputDir = TString::Format("run%d/", outputRun);
  TString mkdirCommand("mkdir -p ");
  
  system(mkdirCommand + outputDir);

  // create only makes a file if one doesn't exit!
  const char* optionForTFile = "create";
  
  TFile* fNewHeadFile = new TFile(outputDir + TString::Format("headFile%d.root", outputRun), optionForTFile);

  if(!fNewHeadFile){
    return 1;
  }
  
  TTree* fHeadTree = new TTree("headTree", "headTree");
  RawAnitaHeader* header = NULL;
  fHeadTree->Branch("header", &header);


  TFile* fNewEventFile = new TFile(outputDir + TString::Format("calEventFile%d.root", outputRun), optionForTFile);
  if(!fNewEventFile){
    return 1;
  }
  
  TTree* fEventTree = new TTree("eventTree", "eventTree");
  CalibratedAnitaEvent* calEvent = NULL;
  fEventTree->Branch("event", &calEvent);

  TFile* fNewGpsFile = new TFile(outputDir + TString::Format("gpsEvent%d.root", outputRun), optionForTFile);
  if(!fNewGpsFile){
    return 1;
  }
  
  TTree* fGpsTree = new TTree("adu5PatTree", "adu5PatTree");
  Adu5Pat* pat = NULL;
  fGpsTree->Branch("pat", &pat);

  int runOfLastEvent = -1;
  AnitaDataset* d = NULL;
  AnitaVersion::set(3);
  for(UInt_t i=0; i < eventNumbers.size(); i++){
    UInt_t eventNumber = eventNumbers[i];
    Int_t run = AnitaDataset::getRunContainingEventNumber(eventNumber);

    if(run!=runOfLastEvent){
      if(d) {
	delete d;
	d = NULL;
      }
      d = new AnitaDataset(run);
    }

    d->getEvent(eventNumber);

    header = d->header();
    fHeadTree->Fill();

    calEvent = d->calibrated();
    fEventTree->Fill();

    pat = d->gps();
    fGpsTree->Fill();

    runOfLastEvent = run;

  }


  fNewHeadFile->Write();
  fNewHeadFile->Close();

  fNewEventFile->Write();
  fNewEventFile->Close();

  fNewGpsFile->Write();
  fNewGpsFile->Close();

  if(d){
    delete d;
    d = NULL;
  }  
  
  return 0;
}
















int parseCommandLine(int argc, char* argv[],  std::vector<UInt_t>& eventNumbers){
  Int_t outputRun = -1;
  bool expectingRun = false;
  bool expectingEventNumber = false;

  for(int i=1; i < argc; i++){
    TString opt(argv[i]);
    // std::cout << "argv[" << i << "] = " << opt << std::endl;
    if(opt=="-r"){
      expectingRun = true;
      expectingEventNumber = false;
      continue;
    }
    else if(opt=="-e"){
      expectingEventNumber = true;
      expectingRun = false;
      continue;
    }
    else{
      if(expectingRun && outputRun < 0){
	outputRun = atoi(argv[i]);
      }
      else if(expectingRun && outputRun > 0){
	std::cerr << "Warning already set run to " << outputRun << std::endl;
      }
      else if(expectingEventNumber){
	eventNumbers.push_back(atoi(argv[i]));
      }
      else{
	std::cerr << "Not sure what to do with argument " << opt << std::endl;
      }
    }
  }
  bool enoughArgs = true;
  if(outputRun < 0){
    std::cerr << "Didn't specify output run" << std::endl;
    enoughArgs = false;
  }
  if(eventNumbers.size()==0){
    std::cerr << "Didn't give any event numbers" << std::endl;
    enoughArgs = false;
  }
  
  if(!enoughArgs){
    std::cerr << "Aborting program due to insufficient arguments!" << std::endl;
    std::cerr << argv[0] << "-r run -e event1 event2 event3... " << std::endl;
    
  }
  
  return outputRun;
}
