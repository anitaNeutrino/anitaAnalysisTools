/*****************************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description: 
             Program to test the functionality of the CrossCorrelator class, nothing more complicated than that. Writes some output to /tmp.

*************************************************************************************************************** */


#include "AveragePowerSpectrum.h"
#include "FancyFFTs.h"
#include "RootTools.h"
#include "FFTtools.h"

#include "TGraph.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include <RawAnitaHeader.h>
#include <UsefulAnitaEvent.h>
#include <CalibratedAnitaEvent.h>

void testAveragePowerSpectrum();

int main(){

  testAveragePowerSpectrum();  
  return 0;

}

void testAveragePowerSpectrum(){

  char eventFileName[1024];

  Int_t run = 352;
  sprintf(eventFileName, "~/UCL/ANITA/calibratedFlight1415/run%d/calEventFile%d.root", run, run);

  TFile* eventFile = TFile::Open(eventFileName);
  TTree* eventTree = (TTree*) eventFile->Get("eventTree");

  char rawHeaderFileName[1024];
  sprintf(rawHeaderFileName, "~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);

  TFile* rawHeaderFile = TFile::Open(rawHeaderFileName);
  TTree* headTree = (TTree*) rawHeaderFile->Get("headTree");

  CalibratedAnitaEvent* event = NULL;  
  eventTree->SetBranchAddress("event", &event);

  RawAnitaHeader* header = NULL;
  headTree->SetBranchAddress("header", &header);

  TFile* outFile = new TFile("/tmp/testAveragePowerSpectrum.root","recreate");


  const Long64_t numEntries = 1000;
  // const Int_t numSamps = 256;
  const Int_t numSamps = 256;  
  const Double_t deltaT = 1./2.6;

  AveragePowerSpectrum aps("aps", "Antenna 16TH");
  AveragePowerSpectrum aps2("aps2", "Antenna 8MH");
  
  for(Long64_t entry=0; entry<numEntries; entry++){

    headTree->GetEntry(entry);
    eventTree->GetEntry(entry);
    UsefulAnitaEvent* realEvent = new UsefulAnitaEvent(event);

    TGraph* gr = realEvent->getGraph(16, AnitaPol::kHorizontal);
    TGraph* grInterp = RootTools::interpolateWithStartTime(gr, gr->GetX()[0], deltaT, numSamps);
    // RootTools::printTGraphInfo(grInterp);
    aps.add(grInterp);
    delete gr;
    delete grInterp;

    TGraph* gr2 = realEvent->getGraph(16+8, AnitaPol::kHorizontal);
    TGraph* grInterp2 = RootTools::interpolateWithStartTime(gr2, gr2->GetX()[0], deltaT, numSamps);
    aps2.add(grInterp2);
    delete gr2;
    delete grInterp2;

    delete realEvent;
    
  }
  
  TGraph* gr = aps.makeAvePowSpecTGraph();
  gr->Write();
  TGraph* gr2 = aps2.makeAvePowSpecTGraph_dB();
  gr2->Write();

  aps.fitAllRayleighHistograms();

  aps.Write();
  aps2.Write();  
  
  outFile->Write();
  outFile->Close();

}



