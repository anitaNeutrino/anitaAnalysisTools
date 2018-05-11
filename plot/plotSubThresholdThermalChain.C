// #include "TGraphAntarctica.h"
// #include "TH2DAntarctica.h"
// #include "ProgressBar.h"
#include "DrawStrings.h"
#include "ThermalChain.h"

using namespace Acclaim;

void plotSubThresholdThermalChain(const char* glob, double deltaF = 0, double deltaF2 = 0){

  // // unknown-base singlets
  // std::vector<UInt_t> eventsH;
  // eventsH.push_back( 7613856);
  // eventsH.push_back( 9097075);
  // eventsH.push_back(11116669);
  // eventsH.push_back(11989349);
  // eventsH.push_back(15717147); // Peter's
  // eventsH.push_back(16952229);
  // eventsH.push_back(19459851);
  // eventsH.push_back(22345215);
  // eventsH.push_back(23695286);
  // eventsH.push_back(27142546);
  // eventsH.push_back(32907848);
  // eventsH.push_back(33484995);
  // eventsH.push_back(41529195);
  // eventsH.push_back(48837708);
  // eventsH.push_back(58592863);
  // eventsH.push_back(62273732);
  // eventsH.push_back(65187079);
  // eventsH.push_back(66313844);
  // eventsH.push_back(68298837);
  // eventsH.push_back(70013898);
  // eventsH.push_back(71171108);
  // eventsH.push_back(71766273);
  // eventsH.push_back(73726742);
  // eventsH.push_back(74592579);
  // eventsH.push_back(75277769);
  // eventsH.push_back(80840274);

  // std::vector<UInt_t> eventsV;
  // eventsV.push_back(83877990); // VPol

  // ThermalChain tc2(glob);
  // TString cutStr = "0";
  // for(h : eventsH){
  //   cutStr += TString::Format(" || eventNumber==%u", h);    
  // }
  // TCut hPolEvents = 

  const TCut subThresholdFisherCut("subThresholdFisherCut", ThermalTree::fisherDiscriminant + TString::Format(" > %lf", ThermalTree::fisherThreshold - deltaF));


  bool topCut = true;
  const TCut notAboveThresholdFisherCut("notAboveThresholdFisherCut", ThermalTree::fisherDiscriminant + TString::Format(" <= %lf", ThermalTree::fisherThreshold - deltaF2));
  
  ThermalChain tc(glob);  
  tc.addCut( !ThermalTree::isAboveHorizontal);
  tc.addCut(  ThermalTree::passAllQualityCuts);
  tc.addCut(  ThermalTree::isNotTaggedAsPulser);
  tc.addCut(  subThresholdFisherCut);
  if(TMath::Abs(deltaF) >  0 && topCut){
    tc.addCut(notAboveThresholdFisherCut);
  }
  tc.addCut( !ThermalTree::closeToHiCal);
  tc.addCut(  ThermalTree::closeToMC);

  TCut hackyCut("longitude < 90 && longitude > 0 && latitude < -73");
  tc.addCut(hackyCut);

  const Long64_t n = tc.N();

  std::cout << "There are " << n << " events to plot" << std::endl;
  
  TH2DAntarctica* h = new TH2DAntarctica("h", "h");
  ProgressBar p(n);
  for(Long64_t entry=0; entry < n; entry++){
    tc.getEntry(entry);
    // Acclaim::RootTools::vectorContainsValue(eventsH, )
    h->Fill(tc.longitude, tc.latitude);
    std::cout << tc.run << "\t" << tc.eventNumber << std::endl;
    p.inc(entry);
  }

  h->SetGrayScale(1);
  h->SetIcemask(1);
  auto c = new TCanvas();
  h->Draw("colz");
  h->ShowBackgroundColorAxis(false);
  c->SetLogz(1);

  // TString canName = TString::Format("projection_dF_%2.1lf_dF2_%2.1lf", deltaF, deltaF2);
  // if(!topCut){
  //   canName = TString::Format("projection_dF_%2.1lf", deltaF);
  // }
  // RootTools::saveCanvas(c, canName);
  
}
