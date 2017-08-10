#include "SummarySet.h"
#include "TChain.h"
#include "AnitaEventSummary.h"
#include <iostream>
#include "TH2D.h"
#include "AnalysisPlot.h"
#include "TFile.h"

Acclaim::SummarySet::SummarySet(const char* pathToSummaryFiles, const char* treeName, const char* summaryBranchName)
    : fPathToSummaryFiles(pathToSummaryFiles), fTreeName(treeName), fSummaryBranchName(summaryBranchName),
      fChain(NULL), fSum(NULL), fFirstTime(0), fFirstEventNumber(0), fLastTime(0), fLastEventNumber(0) {

  init();
}



Acclaim::SummarySet::~SummarySet(){
  delete fChain;
  fChain = NULL;

  if(fSum){
    delete fSum;
  }
  fSum = NULL;

  if(gFile && gFile->IsWritable()){
    TString name = "";
    TString title = "Acclaim::SummarySet ss(" + fPathToSummaryFiles + ", " + fTreeName + "," + fSummaryBranchName + ")";
    TNamed note(name, title);
    note.Write();
  }
}





void Acclaim::SummarySet::init(){


  fChain = new TChain(fTreeName);
  fChain->Add(fPathToSummaryFiles);
  fN = fChain->GetEntries();
  fChain->SetBranchAddress(fSummaryBranchName, &fSum);

  if(fN == 0){
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", no entries in TChain of AnitaEventSummary" << std::endl;
  }
  else{
    last();
    fLastTime = fSum->realTime;
    fLastEventNumber = fSum->eventNumber;

    first();
    fFirstTime = fSum->realTime;
    fFirstEventNumber = fSum->eventNumber;

    // std::cerr << "Just called " << __PRETTY_FUNCTION__ << ", time range is " << fFirstTime << " to " << fLastTime << "\t eventNumber range is " << fFirstEventNumber << " to " << fLastEventNumber << std::endl;
  }
}

Long64_t Acclaim::SummarySet::getEntry(Long64_t entry){
  return fChain->GetEntry(entry);
}



Acclaim::AnalysisProf* Acclaim::SummarySet::bookTimeAnalysisProf(const char* name, const char* title, int nx, int ny, double yMin, double yMax){
  AnalysisProf* h = new AnalysisProf(name, title, nx, getFirstTime(), getLastTime(), ny, yMin, yMax);
  return h;
}


Acclaim::AnalysisProf* Acclaim::SummarySet::bookEventNumberAnalysisProf(const char* name, const char* title, int nx, int ny, double yMin, double yMax){
  AnalysisProf* h = new AnalysisProf(name, title, nx, getFirstEventNumber(), getLastEventNumber(), ny, yMin, yMax);
  return h;
}



Acclaim::AnalysisPlot* Acclaim::SummarySet::bookTimeAnalysisPlot(const char* name, const char* title, int nx, int ny, double yMin, double yMax){
  AnalysisPlot* h = new AnalysisPlot(name, title, nx, getFirstTime(), getLastTime(), ny, yMin, yMax);
  return h;
}


Acclaim::AnalysisPlot* Acclaim::SummarySet::bookEventNumberAnalysisPlot(const char* name, const char* title, int nx, int ny, double yMin, double yMax){
  AnalysisPlot* h = new AnalysisPlot(name, title, nx, getFirstEventNumber(), getLastEventNumber(), ny, yMin, yMax);
  return h;
}



TH2D* Acclaim::SummarySet::bookTimeHistogram(const char* name, const char* title, int nx, int ny, double yMin, double yMax){
  TH2D* h = new TH2D(name, title, nx, getFirstTime(), getLastTime(), ny, yMin, yMax);
  return h;
}


TH2D* Acclaim::SummarySet::bookEventNumberHistogram(const char* name, const char* title, int nx, int ny, double yMin, double yMax){
  TH2D* h = new TH2D(name, title, nx, getFirstEventNumber(), getLastEventNumber(), ny, yMin, yMax);
  return h;
}


