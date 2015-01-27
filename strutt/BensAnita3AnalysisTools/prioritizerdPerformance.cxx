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

  TChain* headChain = new TChain("headTree");
  TChain* gpsChain = new TChain("adu5PatTree");
  Int_t firstRun = 330;
  Int_t lastRun = 356;
  for(int run=firstRun; run<lastRun; run++){
    char fileName[1024];
    sprintf(fileName, "root/run%d/headFile%d.root", run, run);
    headChain->Add(fileName);

    sprintf(fileName, "root/run%d/gpsFile%d.root", run, run);
    gpsChain->Add(fileName);
  }
  RawAnitaHeader* header = NULL;
  headChain->SetBranchAddress("header", &header);

  /* Setting branch address for gpsChain no longer required because of FancyTTreeInterpolator */
  FancyTTreeInterpolator fancyGpsInterp(gpsChain, "pat->realTime");

  // -1000ish is some kind of bad state flag, also added wrapped value at 360 for interpolation
  fancyGpsInterp.add("pat->heading", "pat->heading>-500", 360.0); 
  fancyGpsInterp.add("pat->altitude");
  fancyGpsInterp.add("pat->latitude");
  fancyGpsInterp.add("pat->longitude");

  Long64_t numHeaderEntries = headChain->GetEntries();
  std::cout << numHeaderEntries << " " << gpsChain->GetEntries() << std::endl;

  TFile* outFile = new TFile("prioritizerdPerformancePlots.root", "recreate");

  TH1D* hGpsValidationHeading = new TH1D("hGpsValidationHeading", "heading from TTree - heading from FancyTTreeInterpolator; Heading difference (degrees); Fraction of events", 100, -10, 10);
  TH2D* hGpsValidationHeading2 = new TH2D("hGpsValidationHeading2", "heading from TTree vs heading from FancyTTreeInterpolator; ADU5 heading from TTree (degrees); ADU5 heading from TTree interpolator (degrees)", 1024, 0, numHeaderEntries, 1024, 0, 360);

  std::string h0Title = "Near WAIS divide: TriggerTimeNs for runs " 
    + std::to_string(firstRun) + "-" + std::to_string(lastRun) 
    + "; run; triggerTimeNs (ns)";

  TH2D* h0 = new TH2D("h0", h0Title.c_str(), lastRun+1-firstRun, firstRun, lastRun+1, 64, 0, 3e6);


  std::string h1Title = "triggerTimeNs vs TriggerTimeNsExpected WAIS pulser runs " + std::to_string(firstRun) + "-" + std::to_string(lastRun) + "; TriggerTimeNs (ns); Expected triggerTimeNs (ns)";
  TH2D* h1 = new TH2D("h1", h1Title.c_str(), 64, 1e6, 2e6, 64, 1e6, 2e6);
  
  std::string h2Title = "TriggerTimeNsExpected - triggerTimeNs WAIS pulser runs " + std::to_string(firstRun) + "-" + std::to_string(lastRun) + "; #deltat (ns); Fraction of events";
  TH1D* h2 = new TH1D("h2", h2Title.c_str(), 10000, -1e9, 1e9);

  
  
  std::string h3Title = "TriggerTimeNsExpected - triggerTimeNs WAIS pulser runs " + std::to_string(firstRun) + "-" + std::to_string(lastRun) + "; #deltat (ns); Fraction of events";
  TH1D* h3 = new TH1D("h3", h3Title.c_str(), 128, -2000, 2000);
  TH1D* h4 = new TH1D("h4", "Priority; Priority; Fraction of events", 10, 0, 10);
  TH2D* h15 = new TH2D("h15", "run vs Priority of WAIS pulses; Run; Priority", lastRun+1-firstRun, firstRun, lastRun+1, 10, 0, 10);

  TH1D* h5 = new TH1D("h5", "Priorities near WAIS divide, runs 331-354; Priority; Fraction of events", 10, 0, 10);
  TH2D* h6 = new TH2D("h6", "Phi-sector trigger HPOL vs phi-expected near WAIS divide, runs 331-354; Phi-Sector; Phi expected (degrees)", 16, 0, 16, 16, 0, 360);
  TH2D* h7 = new TH2D("h7", "GPU phi vs Expected phi near WAIS divide, runs 331-354; GPU peak phi (Degrees); Phi expected (degrees)", 32, 0, 360, 32, 0, 360);
  TH1D* h8 = new TH1D("h8", "GPU phi - Expected phi (+45) near WAIS divide, runs 331-354; GPU peak phi - Phi expected (degrees); Fraction of events", 360, -180, 180);
  TH2D* h9 = new TH2D("h9", "Heading - Expected phi near WAIS divide, runs 331-354; Entry; Heading - Phi expected (degrees)", 1024, 0, numHeaderEntries, 360, -180, 180);
  TH2D* h10 = new TH2D("h10", "L3 trigger HPOL vs GPU Peak PhiSector, runs 331-354; l3trigger; GPU Peak Phi", 16, 0, 16, 720, 0, 360);
  TH1D* h11 = new TH1D("h11", "GPU theta - Expected theta, near WAIS divide, runs 331-354; GPU Peak Theta - Theta expected (degrees); Fraction of events", 512, -150, 150);
  TH2D* h12 = new TH2D("h12", "GPU theta vs Expected theta, near WAIS divide, runs 331-354; GPU Peak Theta (degrees); Expected theta (degrees)", 150, -75, 75, 150, -100, 100);
  for(Long64_t entry = 0; entry < numHeaderEntries; entry++){
    headChain->GetEntry(entry);

    if(entry==0 || entry == numHeaderEntries-1){
      if(entry==0) std::cout << "start time : ";
      else std::cout << "end time : ";
      std::cout << header->realTime << std::endl;
    }

    if(!(header->realTime >= fancyGpsInterp.fXmin && header->realTime <= fancyGpsInterp.fXmax)) continue;

    if(header->trigType & 1){

      Adu5Pat pat2;
      pat2.heading = fancyGpsInterp.interp("pat->heading", header->realTime);
      pat2.altitude = fancyGpsInterp.interp("pat->altitude", header->realTime);
      pat2.latitude = fancyGpsInterp.interp("pat->latitude", header->realTime);
      pat2.longitude = fancyGpsInterp.interp("pat->longitude", header->realTime);
      UsefulAdu5Pat usefulPat2(&pat2);

      UInt_t triggerTimeNsExpected = usefulPat2.getTriggerTimeNsFromSource(sourceLat, sourceLon, sourceAlt);

      //      hGpsValidationHeading->Fill(usefulPat.heading - usefulPat2.heading);
      hGpsValidationHeading2->Fill(entry, usefulPat2.heading);


      h0->Fill(header->run, header->triggerTimeNs);
      h1->Fill(header->triggerTimeNs, triggerTimeNsExpected);
      h2->Fill(Double_t(triggerTimeNsExpected) - Double_t(header->triggerTimeNs));
      h3->Fill(Double_t(triggerTimeNsExpected) - Double_t(header->triggerTimeNs));

      if(TMath::Abs(Double_t(header->triggerTimeNs) - Double_t(triggerTimeNsExpected)) < cutTimeNs){

	h4->Fill(header->priority&0xf);
	h15->Fill(header->run, header->priority&0xf);
	// std::cout << header->l3TrigPatternH << ": ";
	Double_t thetaWave = -1;
	Double_t phiWave = -1;
	usefulPat2.getThetaAndPhiWaveAnita3(sourceLon, sourceLat, sourceAlt, thetaWave, phiWave);
	//usefulPat.getThetaAndPhiWave(pat->longitude, pat->latitude, 0, thetaWave, phiWave);
	Double_t thetaDeg = thetaWave*TMath::RadToDeg(); 
	Double_t phiDeg = phiWave*TMath::RadToDeg();
	phiDeg = phiDeg < 0 ? phiDeg + 360 : phiDeg;
	phiDeg = phiDeg >= 360 ? phiDeg - 360 : phiDeg;
	for(int phi=0; phi<16; phi++){
	  Int_t val = ((header->l3TrigPatternH >> phi) & 1);
	  if(val){
	    h6->Fill(phi, phiDeg);

	    //Int_t phiInd2=(header->prioritizerStuff&0x8fff)>>1;
	    Int_t phiInd2=(header->prioritizerStuff&0x7ff)>>1;
	    Int_t phiInd=(phiInd2&0x3f);
	    Int_t phiSector=(phiInd2>>6);
	    if(phiSector > 15){
	      std::cout << phiSector << " " << phiInd << " " << phiInd2 << std::endl;
	    }

	    Double_t gpuPhiDeg = 22.5*phiSector - 11.5 + 22.5*phiInd/64;
	    gpuPhiDeg = gpuPhiDeg < 0 ? gpuPhiDeg + 360 : gpuPhiDeg;

	    Double_t gpuPhiDeg2 = header->getPeakPhiDeg();
	    gpuPhiDeg2 = gpuPhiDeg2 < 0 ? gpuPhiDeg2 + 360 : gpuPhiDeg2;

	    h10->Fill(phi, gpuPhiDeg2);
	  }
	}
	Double_t gpuPhi = header->getPeakPhiDeg();
	if(gpuPhi < 0){
	  // std::cout << gpuPhi << std::endl;
	  gpuPhi += 360;
	}

	h7->Fill(header->getPeakPhiDeg(), phiDeg);

	Double_t deltaPhiGPUvsExpected = gpuPhi - phiDeg;
	deltaPhiGPUvsExpected = deltaPhiGPUvsExpected < -180 ? deltaPhiGPUvsExpected + 360 : deltaPhiGPUvsExpected;
	deltaPhiGPUvsExpected = deltaPhiGPUvsExpected > 180 ? deltaPhiGPUvsExpected - 360 : deltaPhiGPUvsExpected;
	h8->Fill(deltaPhiGPUvsExpected);
	h11->Fill(header->getPeakThetaDeg() - thetaDeg);
	h12->Fill(header->getPeakThetaDeg(), thetaDeg);

	Double_t deltaHeadingvsPhiExpected = usefulPat2.heading - phiDeg;
	deltaHeadingvsPhiExpected = deltaHeadingvsPhiExpected < -180 ? deltaHeadingvsPhiExpected + 360 : deltaHeadingvsPhiExpected;
	deltaHeadingvsPhiExpected = deltaHeadingvsPhiExpected > 180 ? deltaHeadingvsPhiExpected - 360 : deltaHeadingvsPhiExpected;
	h9->Fill(entry, deltaHeadingvsPhiExpected);
      }
      else{
	h5->Fill(header->priority&0xf);
      }
    }
  }

  const Int_t nx = 10;
  const char *pris[nx] = {"0","1","2","3","4","5","6","7","8","9"};
  for (int i=1;i<=nx;i++) {
    h5->GetXaxis()->SetBinLabel(i,pris[i-1]);
  }
  h5->GetXaxis()->SetLabelOffset(0.005);
  h5->GetXaxis()->SetLabelSize(0.04);




  lastRun = 356;
  TChain* headChain2 = new TChain("headTree");
  for(int run=firstRun; run<lastRun; run++){
    char fileName[1024];
    sprintf(fileName, "root/run%d/eventHeadFile%d.root", run, run);
    headChain2->Add(fileName);
  }
  headChain2->SetBranchAddress("header", &header);


  TH1D* h13 = new TH1D("h13", "Number of telemetered WAIS divide pulses; Run; Number of Events", lastRun-firstRun+1, firstRun, lastRun);
  h13->Sumw2();
  TH2D* h14 = new TH2D("h14", "Number of telemetered WAIS divide pulses vs Priority; Run; Priority", lastRun-firstRun+1, firstRun, lastRun, 10, 0, 10);
  for (int i=1;i<=nx;i++) {
    h14->GetYaxis()->SetBinLabel(i,pris[i-1]);
  }
  h14->GetYaxis()->SetLabelOffset(0.005);
  h14->GetYaxis()->SetLabelSize(0.04);

  TH1D* h16 = new TH1D("h16", "Number of telemetered events; Run; Number of Events", lastRun-firstRun+1, firstRun, lastRun);
  h15->Sumw2();


  numHeaderEntries = headChain2->GetEntries();
  std::cout << "numHeaderEntries = " << numHeaderEntries << std::endl;
  Double_t count1 = 0;
  Double_t count2 = 0;
  for(Long64_t entry = 0; entry < numHeaderEntries; entry++){
    headChain2->GetEntry(entry);

    if(entry==0 || entry == numHeaderEntries-1){
      if(entry==0) std::cout << "start time : ";
      else std::cout << "end time : ";
      std::cout << header->realTime << std::endl;
    }

    if(!(header->realTime >= fancyGpsInterp.fXmin && header->realTime <= fancyGpsInterp.fXmax)) continue;

    if(header->trigType & 1){
      
      h16->Fill(header->run);
      
      Adu5Pat pat2;
      pat2.heading = fancyGpsInterp.interp("pat->heading", header->realTime);
      pat2.altitude = fancyGpsInterp.interp("pat->altitude", header->realTime);
      pat2.latitude = fancyGpsInterp.interp("pat->latitude", header->realTime);
      pat2.longitude = fancyGpsInterp.interp("pat->longitude", header->realTime);
      UsefulAdu5Pat usefulPat2(&pat2);

      UInt_t triggerTimeNsExpected = usefulPat2.getTriggerTimeNsFromSource(sourceLat, sourceLon, sourceAlt);

      if(TMath::Abs(Double_t(header->triggerTimeNs) - Double_t(triggerTimeNsExpected)) < cutTimeNs){
	h13->Fill(header->run);
	h14->Fill(header->run, header->priority&0xf);
      }

      // std::cout << count1 << " " << count2 << std::endl;
    }

  }
  

  outFile->Write();
  outFile->Close();

  return 0;

}
