#include "CrossCorrelator.h"

#include <TMath.h>
#include <TApplication.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TROOT.h>

#include <array> 

int main(int argc, char *argv[]){

  TApplication* theApp = new TApplication("App", &argc, argv);
  gROOT->ProcessLine(".x ~/.rootlogon.C");
  
  /* Consider two antennas on middle ring */
  
  Int_t ant1 = 16;
  Int_t ant2 = 19;

  const int numRVals = 4;
  Double_t rVals[numRVals] = {10, 100, 1000, 1000};
  TH2D* hDts[numRVals] = {NULL};

  CrossCorrelator* cc = new CrossCorrelator();

  for(Int_t rInd = 0; rInd<numRVals; rInd++){
    Double_t rWave = rVals[rInd];    
    TString name = TString::Format("hDt%d", rInd);
    hDts[rInd] = new TH2D(name, name,
			  NUM_BINS_PHI*NUM_PHI, cc->phiArrayDeg[0]-PHI_RANGE/2, cc->phiArrayDeg[15]+PHI_RANGE/2,
			  NUM_BINS_THETA, -THETA_RANGE/2, THETA_RANGE/2);

    for(Int_t phiSector = 0; phiSector<NUM_PHI; phiSector++){
      for(Int_t phiInd = 0; phiInd < NUM_BINS_PHI; phiInd++){
	Int_t phiBin = phiSector*NUM_BINS_PHI + phiInd;
	Double_t phiDeg = -0.5*PHI_RANGE + phiBin*Double_t(PHI_RANGE)/NUM_BINS_PHI;
	Double_t phiWave = TMath::DegToRad()*phiDeg;
	for(Int_t thetaBin = 0; thetaBin < NUM_BINS_THETA; thetaBin++){
	  Double_t thetaDeg = THETA_RANGE*((Double_t)thetaBin/NUM_BINS_THETA - 0.5);
	  Double_t thetaWave = TMath::DegToRad()*thetaDeg;
	  Int_t offset = cc->getDeltaTExpectedSpherical(ant1, ant2, phiWave, thetaWave, rWave);
	  hDts[rInd]->SetBinContent(phiBin, thetaBin, offset);
	}
      }
    }
  }

  TCanvas* c1 = new TCanvas();
  c1->Divide(2, 2);
  for(Int_t rInd=0; rInd<numRVals; rInd++){
    c1->cd(rInd+1);
    hDts[rInd]->Draw("colz");
  }

  c1->Update();
  gSystem->ProcessEvents();
  std::cerr << "Press ctrl+c to quit." << std::endl;
  theApp->Run();
  
  return 0;
}
