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
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.trainingPeak().phi", false)); // debugging
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.trainingPeak().theta", false)); // debugging
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.trainingPeak().dPhiTagged()", false)); // debugging
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.trainingPeak().dThetaTagged()", false)); // debugging

  // map info
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.trainingPeak().value", false)); // map values
  treeFormulas.push_back(CutOptimizer::FormulaString("TMath::Abs(sum.trainingPeak().dPhiSun())", false)); // delta phi sun
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.trainingPeak().triggered", false));
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.trainingPeak().triggered_xpol", false));
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.trainingPeak().masked", false));
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.trainingPeak().masked_xpol", false));
  treeFormulas.push_back(CutOptimizer::FormulaString("TMath::Abs(sum.trainingPeak().minAbsHwAngle())", false));
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.trainingPeak().absHwAngleLessThanAbsHwAngleXPol()", false));

  // waveform info
  // treeFormulas.push_back(CutOptimizer::FormulaString("sum.trainingCoherentFiltered().standardizedPeakMoment(1)", true));
  // treeFormulas.push_back(CutOptimizer::FormulaString("sum.trainingDeconvolvedFiltered().standardizedPeakMoment(1)", true));

  treeFormulas.push_back(CutOptimizer::FormulaString("sum.trainingDeconvolvedFiltered().impulsivityMeasure", true));
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.trainingDeconvolved().impulsivityMeasure", true));
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.trainingDeconvolvedFiltered().narrowestWidthsGradient()", true));
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.trainingDeconvolved().narrowestWidthsGradient()", true));

  treeFormulas.push_back(CutOptimizer::FormulaString("sum.trainingDeconvolvedFiltered().snr", true));
  // treeFormulas.push_back(CutOptimizer::FormulaString("TMath::Abs(sum.trainingDeconvolvedFiltered().localMaxToMinTime)", true));
  // treeFormulas.push_back(CutOptimizer::FormulaString("TMath::Abs(sum.trainingDeconvolvedFiltered().globalMaxToMinTime)", true));
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.trainingDeconvolvedFiltered().peakHilbert", true));
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.trainingDeconvolvedFiltered().linearPolFrac()", true));
  // treeFormulas.push_back(CutOptimizer::FormulaString("sum.trainingDeconvolvedFiltered().circPolFrac()", true));

  // treeFormulas.push_back(CutOptimizer::FormulaString("sum.trainingDeconvolved().totalPower/sum.trainingDeconvolvedFiltered().totalPower", true));
  // treeFormulas.push_back(CutOptimizer::FormulaString("sum.trainingDeconvolved().peakHilbert/sum.trainingDeconvolvedFiltered().peakHilbert", true));
  // treeFormulas.push_back(CutOptimizer::FormulaString("sum.trainingCoherent().peakHilbert/sum.trainingCoherentFiltered().peakHilbert", true));


  std::vector<const AnalysisCut*> waisCuts;
  waisCuts.push_back(&Acclaim::AnalysisCuts::isTaggedAsWaisPulser);
  co.addSpectatorTree("waisTree", backgroundGlob, waisCuts);

  std::vector<const AnalysisCut*> selectingBlastsCuts;
  selectingBlastsCuts.push_back(&Acclaim::AnalysisCuts::isTaggedAsPayloadBlast);
  co.addSpectatorTree("blastTree", backgroundGlob, selectingBlastsCuts);
  
  co.optimize(signalSelection, backgroundSelection, treeFormulas, outFileName);

  return 0;
}
