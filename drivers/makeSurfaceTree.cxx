#include <iostream>
#include "ThermalChain.h"
#include "TMath.h"
#include "UsefulAdu5Pat.h"
#include "OutputConvention.h"
#include "ProgressBar.h"

using namespace Acclaim;

int main(int argc, char* argv[]){

  if(argc != 3){
    std::cerr << argv[0] << " 'glob' run " << std::endl;
    return 1;
  }
  const char* glob = argv[1];
  Int_t run = atoi(argv[2]);

  TCut runCut(TString::Format("run==%d", run));
  
  ThermalChain c(glob);
  c.setCut(runCut);
  const Long64_t n = c.N();

  const int nFakeArgv = 2;
  const char* fakeArgv[nFakeArgv] = {argv[0], argv[2]};
  OutputConvention oc(nFakeArgv, (char**)fakeArgv);
  TFile* fOut = oc.makeFile();
  TTree* surfaceTree = new TTree("surfaceTree", "surfaceTree");
  
  Double_t lon;
  Double_t lat;
  Double_t alt;
  Double_t adjust;
  Int_t onContinent;
  Int_t onIceShelf;
  Double_t iceThickness;

  surfaceTree->Branch("longitude", &lon);
  surfaceTree->Branch("latitude", &lat);
  surfaceTree->Branch("altitude", &alt);
  surfaceTree->Branch("thetaAdjustmentRequired", &adjust);
  surfaceTree->Branch("onContinent", &onContinent);
  surfaceTree->Branch("onIceShelf", &onIceShelf);
  surfaceTree->Branch("iceThickness", &iceThickness);
  surfaceTree->Branch("eventNumber3", &c.eventNumber);
  surfaceTree->Branch("run3", &c.run);

  ProgressBar p(n);
  for(Long64_t entry=0; entry < n; entry++){

    c.getEntry(entry);
    Adu5Pat pat = c.pat();
    UsefulAdu5Pat usefulPat(&pat);

    lon = -9999;
    lat = -9999;
    alt = -9999;
    adjust = -9999;
    onIceShelf = 0;
    onContinent = 0;
    iceThickness = 0;

    // usefulPat.setDebug(true);    
    usefulPat.traceBackToContinent3(c.peak_phi*TMath::DegToRad(), -c.peak_theta*TMath::DegToRad(), &lon, &lat, &alt, &adjust);//, maxThetaAdjust, 10);

    if(adjust != -9999){
      adjust *= TMath::RadToDeg();
      onContinent = RampdemReader::isOnContinent(lon, lat);
      onIceShelf = RampdemReader::isOnIceShelf(lon, lat);
      iceThickness = RampdemReader::iceThickness(lon, lat);
    }
    
    // std::cout << c.peak_theta << "\t" << lon << "\t" << lat << "\t" << alt << "\t" << adjust << std::endl;
    
    surfaceTree->Fill();
    p.inc(entry);
  }

  fOut->Write();
  fOut->Close();  
  
  return 0;
}
