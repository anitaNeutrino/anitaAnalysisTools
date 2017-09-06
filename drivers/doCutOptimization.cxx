#include "CutOptimizer.h"
#include "AnalysisCuts.h"
#include <iostream>

using namespace Acclaim;

int main(int argc, char* argv[]){

  if(!(argc == 2 || argc==3)){
    std::cerr << argv[0] << " 'signalAndBackgroundGlob' " << std::endl;
    std::cerr << argv[0] << " 'signalGlob' 'backGroundGlob'" << std::endl;
    return 1;
  }
  const char* outFileName = argv[0];
  const char* signalGlob = argv[1];
  const char* backgroundGlob = argc >= 2 ? argv[2] : NULL;

  CutOptimizer co(signalGlob, backgroundGlob, true, true);




  std::vector<const Acclaim::AnalysisCut *> signalSelection;
  TString sg = signalGlob;
  if(sg.Contains("Wais")){
    // extra data quality cuts
    signalSelection.push_back(&AnalysisCuts::closeToWais); // Is this the right peak?
  }
  else{
    signalSelection.push_back(&AnalysisCuts::closeToMC); // Is this the right peak?
  }
  

  std::vector<const Acclaim::AnalysisCut *> backgroundSelection;
  backgroundSelection.push_back(&AnalysisCuts::isAboveHorizontal); // Upward pointing
  backgroundSelection.push_back(&AnalysisCuts::anita3QuietTime); // quiet
  backgroundSelection.push_back(&AnalysisCuts::isNotTaggedAsPulser); // not a pulser...

  const int nGen = 6;
  const AnalysisCut* genericDataQualityCuts[nGen] = {&AnalysisCuts::isGood, // not payload blast, SURF saturation, 
                                                     &AnalysisCuts::smallDeltaRough, // agreement between coarse/fine peak
                                                     &AnalysisCuts::goodGPS, // do we have GPS data?
                                                     &AnalysisCuts::nonZeroStokesI, // other badness
                                                     &AnalysisCuts::realSNR, // other badness (should be fixed with new noise monitor)
                                                     &AnalysisCuts::isRfTrigger
                                                     }; // for hardware angle related stuff
  
  for(unsigned i=0; i < nGen; i++){
    signalSelection.push_back(genericDataQualityCuts[i]);
    backgroundSelection.push_back(genericDataQualityCuts[i]);    
  }






  
  
  
  std::vector<CutOptimizer::FormulaString> treeFormulas;
  treeFormulas.push_back(CutOptimizer::FormulaString("run", false));
  treeFormulas.push_back(CutOptimizer::FormulaString("eventNumber", false)); // these get excluded from training by the CutOptimizer
  treeFormulas.push_back(CutOptimizer::FormulaString("realTime", false)); // these get excluded from training by the CutOptimizer
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.weight()", false)); // Automatically filled with 1 if non-MC
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.mc.energy", false)); // The MC energy
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.mcPeak().phi", false)); // debugging
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.mcPeak().theta", false)); // debugging
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.mcPeak().dPhiMC()", false)); // debugging
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.mcPeak().dThetaMC()", false)); // debugging

  treeFormulas.push_back(CutOptimizer::FormulaString("sum.mcPeak().absHwAngle()", true));  // definitionally can't do better than the hardware trigger, I think  

  // map info
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.mcPeak().value", true)); // map values
  treeFormulas.push_back(CutOptimizer::FormulaString("TMath::Abs(sum.mcPeak().dPhiSun())", true)); // delta phi sun
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.mcPeak().triggered", true));
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.mcPeak().triggered_xpol", true));

  // waveform info
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.mcCoherentFiltered().standardizedPeakMoment(1)", true));
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.mcDeconvolvedFiltered().standardizedPeakMoment(1)", true));  
  
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.mcDeconvolvedFiltered().snr", true));
  treeFormulas.push_back(CutOptimizer::FormulaString("TMath::Abs(sum.mcDeconvolved().impulsivityMeasure)", true));
  treeFormulas.push_back(CutOptimizer::FormulaString("TMath::SignBit(sum.mcDeconvolved().impulsivityMeasure)", true));
  treeFormulas.push_back(CutOptimizer::FormulaString("TMath::Abs(sum.mcDeconvolved().localMaxToMinTime)", true));
  treeFormulas.push_back(CutOptimizer::FormulaString("TMath::Abs(sum.mcDeconvolved().globalMaxToMinTime)", true));
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.mcDeconvolvedFiltered().peakHilbert", true));
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.mcDeconvolvedFiltered().linearPolFrac()", true));
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.mcDeconvolvedFiltered().circPolFrac()", true));

  treeFormulas.push_back(CutOptimizer::FormulaString("sum.mcDeconvolved().totalPower/sum.mcDeconvolvedFiltered().totalPower", true));
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.mcDeconvolved().peakHilbert/sum.mcDeconvolvedFiltered().peakHilbert", true));
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.mcCoherent().peakHilbert/sum.mcCoherentFiltered().peakHilbert", true));

  co.optimize(signalSelection, backgroundSelection, treeFormulas, outFileName);

  return 0;
}
