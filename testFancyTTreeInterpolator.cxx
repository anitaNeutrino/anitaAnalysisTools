#include "FancyTTreeInterpolator.h"

#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
//#include "UsefulAdu5Pat.h"
#include "TGraph.h"


int main(){

  // gSystem->Load("libAnitaEvent.so");

  TFile f("testFancyTTreeInterpolatorOutput.root", "recreate");

  TChain* gpsChain = new TChain("adu5PatTree");
  const Int_t startRun = 331;
  const Int_t endRun = 355;

  for(Int_t run=startRun; run<endRun; run++){
    char gpsFileName[1024];
    sprintf(gpsFileName, "root/run%d/gpsFile%d.root", run, run);
    gpsChain->Add(gpsFileName);
  }


  FancyTTreeInterpolator treeInterp(gpsChain, "pat->realTime");
  treeInterp.add("pat->latitude", "");
  std::cout << treeInterp.fStringToGraph.find("pat->latitude")->first << std::endl;
  int n = treeInterp.fStringToGraph.find("pat->latitude")->second->GetN();
  std::cout << n << std::endl;

  // for(int i=0; i<n; i++){
  //   std::cout << treeInterp.fStringToGraph.find("pat->latitude")->second->GetY()[i] << std::endl;
  // }

  try{// deliberate typo
    Double_t someNum = treeInterp.interp("pat->latitudea", 10000);
  }
  catch(...){
    std::cout << "Caught something!" << std::endl;
  }
  

  try{ // out of range time request
    Double_t someNum = treeInterp.interp("pat->latitude", 10000);
    
    // these lines never get executed...
    std::cout << someNum << std::endl;
    std::cout << treeInterp.fStringToGraph.find("pat->latitude")->second->GetX()[0] << std::endl;
    std::cout << treeInterp.fStringToGraph.find("pat->latitude")->second->GetX()[n-1] << std::endl;
  }
  catch(...){
    std::cout << "Caught something!" << std::endl;
  }


  std::shared_ptr<TGraph> gr = treeInterp.get("pat->latitude");
  gr->SetName("grLat");
  gr->Write();


  /* Test unwrapping functionality in heading, should be smooth */
  treeInterp.add("pat->heading", "pat->heading > -500", 360.0);
  std::shared_ptr<TGraph> gr2 = treeInterp.get("pat->heading");
  gr2->SetName("grHeadUnwrapped");
  gr2->Write();


  TGraph* gr3 = new TGraph();
  for(int i=0; i<gr2->GetN(); i++){
    Double_t xVal = gr2->GetX()[i];
    Double_t yVal = treeInterp.interp("pat->heading", xVal);
    gr3->SetPoint(i, xVal, yVal);
  }
  gr3->SetName("grHeadUnwrappedAndInterped");
  gr3->Write();


  /* Test unwrapping functionality in heading, check we get back to the intial graph */
  treeInterp.add("pat->heading", "pat->heading > -500");
  std::shared_ptr<TGraph> gr4 = treeInterp.get("pat->heading");
  gr4->SetName("grHeadNoWrap");
  gr4->Write();



  f.Close();

}


