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
  

  


  const int numRVals = 4;
  Double_t rVals[numRVals] = {5, 6, 7, 8};
  TH2D* hDts[numRVals] = {NULL};

  CrossCorrelator* cc = new CrossCorrelator();
  cc->correlateEventTest(30, -15, 7);

  for(Int_t rInd = 0; rInd<numRVals; rInd++){    
    hDts[rInd] = cc->makeImageSpherical(AnitaPol::kVertical, rVals[rInd]);
    hDts[rInd]->SetName(TString::Format("hSpherical%d", rInd));
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

  delete cc;
  
  return 0;
}
