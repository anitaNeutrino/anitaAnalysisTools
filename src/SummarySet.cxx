#include "SummarySet.h"
#include "TChain.h"
#include "AnitaEventSummary.h"
#include <iostream>
#include "TH2D.h"
#include "AnalysisPlot.h"
#include "TFile.h"
#include "TProof.h"
#include <stdlib.h>
#include "SummarySelector.h"
#include "TCanvas.h"
#include "TSeqCollection.h"
#include "TROOT.h"

Acclaim::SummarySet::SummarySet(const char* pathToSummaryFiles, const char* treeName, const char* summaryBranchName, bool useProof)
    : fPathToSummaryFiles(pathToSummaryFiles), fTreeName(treeName), fSummaryBranchName(summaryBranchName),
      fChain(NULL), fSum(NULL), fFirstTime(0), fFirstEventNumber(0), fLastTime(0), fLastEventNumber(0), fUseProof(useProof) {

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
    TString title = "Acclaim::SummarySet ss(\"" + fPathToSummaryFiles + "\", \"" + fTreeName + "\", \"" + fSummaryBranchName + "\")";
    TNamed note(name, title);
    note.Write();
  }
}



/** 
 * Enable (disable) PROOF, if fUseProof is true (false)
 */
void Acclaim::SummarySet::initProof(){

  if(fUseProof && !gProof){
    TProof::Open("");
    const char* anitaUtilInstallDir = getenv("ANITA_UTIL_INSTALL_DIR");
    TString loadAnita = TString::Format("%s/share/Acclaim/loadAnita.C", anitaUtilInstallDir);
    gProof->Load(loadAnita);
  }
  
  if(fChain){
    fChain->SetProof(fUseProof);
  }
}


Long64_t Acclaim::SummarySet::Process(TList* things){
  (void) things;
  std::cout << fN << " entries" << std::endl;

  initProof();
  SummarySelector* sel = new SummarySelector();
  // AnalysisPlot* p = new AnalysisPlot("h", "asd", 100, 0, 1);
  TList* l = new TList();
  // l->Add(p);
  sel->SetInputList(l);
  fChain->Process(sel);

  TList* l2 = sel->GetOutputList();
  std::cerr << "there are " << l2->GetEntries() << " outputs" << std::endl;
  // const char* anitaUtilInstallDir = getenv("ANITA_UTIL_INSTALL_DIR");
  // TString loadSelector = TString::Format("%s/share/Acclaim/SumSelector.C", anitaUtilInstallDir);
  
  return 0; //fChain->Process(loadSelector);  
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




/** 
 * Little hack to make PROOF give me a histogram and canvas with different names
 * 
 * @param varexp the Draw expression
 */
void Acclaim::SummarySet::renameProofCanvas(const char* varexp){
  TCanvas* c = gPad->GetCanvas();
  TString canName = c->GetName();
  TString command = varexp;
  TString histName = "htemp";
  if(command.Contains(">>")){ // then we have a histogram name, let's go get it
    TObjArray* tkns = command.Tokenize(">>");
    TObjString* s1 = (TObjString*) tkns->At(1);
    TString s1Str = s1->String();
    TObjArray* tkns2 = s1Str.Tokenize("(");
    TObjString* s2 = (TObjString*) tkns2->At(0);
    histName = s2->String();
    delete tkns;
    delete tkns2;
  }

  // for some reason, proof does this...
  // and if I want the histogram on the command line, I need to rename the canvas
  if(canName==histName){
    TString defCanName = gROOT->GetDefCanvasName();
    TSeqCollection* cans = gROOT->GetListOfCanvases();
    bool alreadyHaveThisCanName = true;
    int i=1;
    TString newCanName;
    while(alreadyHaveThisCanName){
      newCanName = i == 1 ? defCanName : defCanName + TString::Format("_n%d", i);
      alreadyHaveThisCanName = cans->Contains(newCanName);
    }
    c->SetName(newCanName);
    c->SetTitle(newCanName);
  }      
}


Long64_t Acclaim::SummarySet::Draw(const char* varexp, const TCut &selection, Option_t *option, Long64_t nentries, Long64_t firstentry){
  initProof();
  Long64_t retVal = fChain->Draw(varexp, selection, option, nentries, firstentry);

  if(fUseProof){
    renameProofCanvas(varexp);
  }
  return retVal;
  
}

Long64_t Acclaim::SummarySet::Draw(const char* varexp, const char* selection, Option_t* option, Long64_t nentries, Long64_t firstentry){
  initProof();
  Long64_t retVal = fChain->Draw(varexp, selection, option, nentries, firstentry);
  if(fUseProof){
    renameProofCanvas(varexp);
  }
  return retVal;
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


