#include "TFile.h"
#include "TTree.h"
#include "UsefulAnitaEvent.h"
#include "CrossCorrelator.h"

void testCorrelatorInRoot(){


  // /* Tedious loading of libraries to make this ROOT macro work... */
  // gSystem->Load("libMathMore.so"); // prereq for Ryan's ROOT FFTW wrapper
  // gSystem->Load("libRootFftwWrapper.so"); // Ryan's ROOT FFTW wrapper, required by EventReaderRoot
  // gSystem->Load("libAnitaEvent.so"); // load event reader ROOT into CINT
  // //  gSystem->Load("libBensAnitaTools.so"); // load my correlator and other lovely things into CINT
  // gSystem->Load("libAnitaAnalysisTools.so"); // load my correlator and other lovely things into CINT




  //  TString dataPath = "~/UCL/ANITA/flight1415/root/";
  TString dataPath = "~/UCL/ANITA/newMonteCarlo/Energy_E21/";
  const Int_t run = 1; //352;//158;




  // Open ROOTified data files...
  //  TString fNameEvent = dataPath + TString::Format("run%d/calEventFile%d.root", run, run);
  TString fNameEvent = dataPath + TString::Format("run%d/SimulatedAnitaEventFile%d.root", run, run);
  TFile* fEvent = TFile::Open(fNameEvent);
  TTree* tEvent = (TTree*) fEvent->Get("eventTree");
  //  CalibratedAnitaEvent* calEvent = NULL;
  UsefulAnitaEvent* calEvent = NULL;
  tEvent->SetBranchAddress("event", &calEvent);

  //  TString fNameHead = dataPath + TString::Format("run%d/headFile%d.root", run, run);
  TString fNameHead = dataPath + TString::Format("run%d/SimulatedAnitaHeadFile%d.root", run, run);
  TFile* fHead = TFile::Open(fNameHead);
  TTree* tHead = (TTree*) fHead->Get("headTree");
  RawAnitaHeader* header = NULL;
  tHead->SetBranchAddress("header", &header);


  // UInt_t eventIWant = 10832108; //10106161;
  // UInt_t eventIWant = 10832108; //10106161;
  // Long64_t entryIWant = -1;
  // Long64_t numEntries = tHead->GetEntries();
  // for(Long64_t entry=0; entry<numEntries; entry++){
  //   tHead->GetEntry(entry);
  //   if(header->eventNumber == eventIWant){
  //     entryIWant = entry;
  //   }
  // }

  int entryIWant = 12;
  tHead->GetEntry(entryIWant);
  tEvent->GetEntry(entryIWant);

  // Make cross correlator
  CrossCorrelator* cc = new CrossCorrelator();

  //  UsefulAnitaEvent* usefulEvent = new UsefulAnitaEvent(calEvent);//, WaveCalType::kDefault, header);
  UsefulAnitaEvent* usefulEvent = calEvent;//, WaveCalType::kDefault, header);

  // cc->correlateEvent(usefulEvent); // Generates set of cross correlations
  cc->reconstructEvent(usefulEvent);

  // TH2D* hImage = cc->makeImage(AnitaPol::kVertical); // then generate and image

  TCanvas* c1 = new TCanvas();

  TH2D* hImage = cc->makeGlobalImage(AnitaPol::kVertical); // then generate and image
  // hImage->SetTitle("A pulser event reconstructed by Cross Correlator");
  // hImage->SetName(TString::Format("hImage%u", eventIWant));



  hImage->Draw("colz"); // ... and enjoy it, isn't it pretty?


  TCanvas* c2 = new TCanvas();
  TH2D* hZoom = cc->getZoomMap(AnitaPol::kVertical);
  hZoom->Draw("colz");

  TCanvas* c3 = new TCanvas();
  Double_t coarseValue, coarsePhiDeg, coarseThetaDeg;
  cc->getCoarsePeakInfo(AnitaPol::kVertical, 0, coarseValue, coarsePhiDeg, coarseThetaDeg);
  Double_t coarseSnr = 0;
  TGraph* gr = cc->makeCoherentlySummedWaveform(AnitaPol::kVertical, coarsePhiDeg, coarseThetaDeg, 0, coarseSnr);
  gr->Draw("al");

  TCanvas* c4 = new TCanvas();
  Double_t fineValue, finePhiDeg, fineThetaDeg;
  cc->getFinePeakInfo(AnitaPol::kVertical, 0, fineValue, finePhiDeg, fineThetaDeg);
  Double_t fineSnr = 0;
  TGraph* grFine = cc->makeUpsampledCoherentlySummedWaveform(AnitaPol::kVertical, finePhiDeg, fineThetaDeg, 0, fineSnr);
  grFine->Draw("al");


}
