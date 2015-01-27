// -*- C++ -*-.
/*****************************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description:
             Program to run over events near WAIS divide and look at Prioritizerd performance.
*************************************************************************************************************** */

#include "TFile.h"
#include "TChain.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "RawAnitaHeader.h"
#include "UsefulAdu5Pat.h"

#include "FancyTTreeInterpolator.h"

// using namespace std;

int main(){

  /* WAIS position things, taken from Steph's e-log */
  Double_t sourceLat = - (79 + (27.93728/60));
  Double_t sourceLon = -(112 + (6.74974/60));
  Double_t sourceAlt = 1813.42;

  Double_t cutTimeNs = 1200;
  // Double_t cutTimeLowNs = -1e3;
  // Double_t cutTimeHighNs = -8e2;
  //  Double_t cutTimeNs = 2e3;

  TChain* headChain = new TChain("headTree");

  TChain* gpsChain = new TChain("adu5PatTree");
  Int_t firstRun = 300; //1;
  Int_t lastRun = 380; //439;
  for(int run=firstRun; run<lastRun; run++){
  //  for(int run=340; run<341; run++){
    char fileName[1024];
    //sprintf(fileName, "root/run%d/eventHeadFile%d.root", run, run);
    sprintf(fileName, "root/run%d/headFile%d.root", run, run);
    headChain->Add(fileName);

    //sprintf(fileName, "root/run%d/gpsEvent%d.root", run, run);
    sprintf(fileName, "root/run%d/gpsFile%d.root", run, run);
    gpsChain->Add(fileName);
  }
  RawAnitaHeader* header = NULL;
  headChain->SetBranchAddress("header", &header);

  FancyTTreeInterpolator fancyGpsInterp(gpsChain, "pat->realTime");
  fancyGpsInterp.add("pat->heading", "pat->heading>-500"); // -1000ish is some kind of bad state/error condition
  fancyGpsInterp.add("pat->altitude");
  fancyGpsInterp.add("pat->latitude");
  fancyGpsInterp.add("pat->longitude");

  Long64_t numHeaderEntries = headChain->GetEntries();
  std::cout << numHeaderEntries << " " << gpsChain->GetEntries() << std::endl;

  TFile* outFile = new TFile("distanceFromWaisDividePlots.root", "recreate");    

  std::string title = "Distance from WAIS divide runs " + std::to_string(firstRun) + "-" + std::to_string(lastRun) +"; run; Distance (km)";
  TH2D* hDistanceFromWais = new TH2D("hDistanceFromWais", title.c_str(), lastRun+1-firstRun, firstRun, lastRun+1, 2048, 0, 2e3);
  for(Long64_t entry = 0; entry < numHeaderEntries; entry++){
    headChain->GetEntry(entry);

    if(entry==0 || entry == numHeaderEntries-1){
      if(entry==0) std::cout << "start time : ";
      else std::cout << "end time : ";
      std::cout << header->realTime << std::endl;
    }

    if(!(header->realTime >= fancyGpsInterp.fXmin && header->realTime <= fancyGpsInterp.fXmax)) continue;
    


    // Long64_t gpsEntry = gpsChain->GetEntryNumberWithIndex(header->realTime);
    // if(gpsEntry < 0 ) continue;

    // if(gpsEntry >= 0 ){
    // 	gpsChain->GetEntry(gpsEntry);
    // }

    // UsefulAdu5Pat usefulPat(pat);

    Adu5Pat pat2;
    pat2.heading = fancyGpsInterp.interp("pat->heading", header->realTime);
    pat2.altitude = fancyGpsInterp.interp("pat->altitude", header->realTime);
    pat2.latitude = fancyGpsInterp.interp("pat->latitude", header->realTime);
    pat2.longitude = fancyGpsInterp.interp("pat->longitude", header->realTime);
    UsefulAdu5Pat usefulPat2(&pat2);


    UInt_t triggerTimeNsExpected = usefulPat2.getTriggerTimeNsFromSource(sourceLat, sourceLon, sourceAlt);

    Double_t distKm = triggerTimeNsExpected*1e-9*C_LIGHT/1e3;
    hDistanceFromWais->Fill(header->run, distKm);
    // std::cout << distKm << std::endl;

  }


  outFile->Write();
  outFile->Close();

  return 0;

}
