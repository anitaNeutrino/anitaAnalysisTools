#include <iostream>
#include "ThermalChain.h"
#include "TMath.h"
#include "UsefulAdu5Pat.h"
#include "OutputConvention.h"

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
  TTree* hiCalTree = new TTree("hiCalTree", "hiCalTree");
  
  Int_t duringHiCal;
  Double_t hiCalPhi;
  Double_t hiCalTheta;

  hiCalTree->Branch("duringHiCal", &duringHiCal);
  hiCalTree->Branch("hiCalPhi", &hiCalPhi);
  hiCalTree->Branch("hiCalTheta", &hiCalTheta);
  hiCalTree->Branch("eventNumber2", &c.eventNumber);
  hiCalTree->Branch("run2", &c.run);

  for(Long64_t entry=0; entry < n; entry++){

    c.getEntry(entry);
    Adu5Pat pat = c.pat();
    UsefulAdu5Pat usefulPat(&pat);

    usefulPat.getThetaAndPhiWaveHiCal(hiCalTheta, hiCalPhi);

    if(hiCalTheta == -9999){
      duringHiCal = 0;
    }
    else{
      duringHiCal = 1;
      hiCalTheta*=TMath::RadToDeg();
      hiCalPhi*=TMath::RadToDeg();
    }
    
    hiCalTree->Fill();
    
  }

  fOut->Write();
  fOut->Close();  
  
  return 0;
}
