#include "SummarySet.h"
#include "TChain.h"
#include "AnitaEventSummary.h"
#include <iostream>
#include "TH2D.h"
#include "AnalysisPlot.h"

Acclaim::SummarySet::SummarySet(const char* pathToSummaryFiles, const char* treeName, const char* summaryBranchName)
    : fChain(NULL), fSum(NULL), fFirstTime(0), fFirstEventNumber(0), fLastTime(0), fLastEventNumber(0) {

  init(pathToSummaryFiles, treeName, summaryBranchName);
}



Acclaim::SummarySet::~SummarySet(){
  delete fChain;
  fChain = NULL;

  if(fSum){
    delete fSum;
  }
  fSum = NULL;
}


void Acclaim::SummarySet::init(const char* pathToSummaryFiles, const char* treeName, const char* summaryBranchName){


  fChain = new TChain(treeName);
  fChain->Add(pathToSummaryFiles);
  fN = fChain->GetEntries();
  fChain->SetBranchAddress(summaryBranchName, &fSum);

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
  }
}

Long64_t Acclaim::SummarySet::getEntry(Long64_t entry){
  return fChain->GetEntry(entry);
}



Acclaim::AnalysisPlot* Acclaim::SummarySet::bookTimeAnalysisPlot(const char* name, const char* title, int nx, int ny, double yMin, double yMax){
  AnalysisPlot* h = new AnalysisPlot(name, title, nx, getFirstEventNumber(), getLastEventNumber(), ny, yMin, yMax);
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


