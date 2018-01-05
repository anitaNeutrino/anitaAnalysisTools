#include "ThermalTreeMaker.h"
#include "SummaryDraw.h"
#include <iostream>
#include "TSystem.h" // require gSystem->Exit(0) to avoid segfault with PROOF

using namespace Acclaim;

int main(int argc, char* argv[]){

  if(argc != 2){
    std::cerr << argv[0] << " 'glob' " << std::endl;
    return 1;
  }
  const char* outFileName = argv[0];
  const char* glob = argv[1];

  ThermalTreeMaker ttm(glob, outFileName);

  std::vector<const char*> formulas;
  formulas.push_back("run");
  formulas.push_back("eventNumber");
  formulas.push_back("realTime");
  formulas.push_back("sum.weight()"); // Automatically filled with 1 if non-MC
  formulas.push_back("sum.mc.energy"); // The MC energy

  formulas.push_back("sum.highestPeak().value");  
  formulas.push_back("sum.highestPeak().phi");
  formulas.push_back("sum.highestPeak().theta");
  formulas.push_back("sum.highestPeak().dPhiTagged()");
  formulas.push_back("sum.highestPeak().dThetaTagged()");
  formulas.push_back("sum.highestPolAsInt()"); // debugging
  formulas.push_back("sum.highestPeakInd()"); // debugging

  formulas.push_back("TMath::Abs(sum.highestPeak().dPhiSun())"); // delta phi sun
  formulas.push_back("TMath::Abs(sum.highestPeak().minAbsHwAngle())");

  formulas.push_back("sum.highestDeconvolvedFiltered().fracPowerWindowGradient()");
  formulas.push_back("sum.highestCoherentFiltered().fracPowerWindowGradient()");

  formulas.push_back("sum.highestDeconvolvedFiltered().impulsivityMeasure");
  formulas.push_back("sum.highestCoherentFiltered().impulsivityMeasure");
  
  formulas.push_back("sum.highestDeconvolvedFiltered().peakHilbert");
  formulas.push_back("sum.highestCoherentFiltered().peakHilbert");


  std::vector<const TCut *> cuts;
  cuts.push_back(&Cuts::highestPeak); // Best peak?
  

  ttm.produce(formulas, cuts);

  gSystem->Exit(0);
  return 0;
}
