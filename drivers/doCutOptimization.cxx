#include "CutOptimizer.h"
#include "AnalysisCuts.h"
#include <iostream>
#include "TSystem.h" // require gSystem->Exit(0) to avoid segfault with PROOF

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

  AnalysisCuts::setMode(AnalysisCuts::kTraining);
  // AnalysisCuts::setMode(AnalysisCuts::kAcclaimAnalysis);  

  CutOptimizer co(signalGlob, backgroundGlob, true, true);

  std::vector<const AnalysisCuts::AnalysisCut *> signalSelection;
  TString sg = signalGlob;
  if(sg.Contains("Wais")){
    // extra data quality cuts
    signalSelection.push_back(&AnalysisCuts::closeToWais); // Is this the right peak?
  }
  else{
    signalSelection.push_back(&AnalysisCuts::closeToMC); // Is this the right peak?
  }

  std::vector<const AnalysisCuts::AnalysisCut *> backgroundSelection;
  backgroundSelection.push_back(&AnalysisCuts::isAboveHorizontal); // Upward pointing
  backgroundSelection.push_back(&AnalysisCuts::anita3QuietTime); // quiet

  const int nGen = 5;
  const AnalysisCuts::AnalysisCut* preThermalCuts[nGen] = {&AnalysisCuts::smallDeltaRough,
  							   &AnalysisCuts::goodGPS,
  							   &AnalysisCuts::realSNR,
							   &AnalysisCuts::nonZeroDenominators,
  							   &AnalysisCuts::isRfTrigger};
  // const AnalysisCuts::AnalysisCut* preThermalCuts[nGen] = {&AnalysisCuts::isGood,
  // 							   &AnalysisCuts::smallDeltaRough,
  // 							   &AnalysisCuts::goodGPS,
  // 							   &AnalysisCuts::realSNR,
  // 							   &AnalysisCuts::isRfTrigger};
  // 							   &AnalysisCuts::higherPeakHilbertAfterDedispersion,
  // 							   &AnalysisCuts::higherImpulsivityMeasureAfterDedispersion};

  for(unsigned i=0; i < nGen; i++){
    signalSelection.push_back(preThermalCuts[i]);
    backgroundSelection.push_back(preThermalCuts[i]);
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
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.trainingPolAsInt()", false)); // debugging
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.trainingPeakInd()", false)); // debugging
  
  // map info
  treeFormulas.push_back(CutOptimizer::FormulaString("TMath::Abs(sum.trainingPeak().dPhiSun())", true)); // delta phi sun
  treeFormulas.push_back(CutOptimizer::FormulaString("TMath::Abs(sum.trainingPeak().minAbsHwAngle())", true));

  treeFormulas.push_back(CutOptimizer::FormulaString("sum.trainingCoherentFiltered().fracPowerWindowGradient()/sum.trainingDeconvolvedFiltered().fracPowerWindowGradient()", true));
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.trainingDeconvolvedFiltered().fracPowerWindowGradient()", true));
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.trainingCoherentFiltered().fracPowerWindowGradient()", true));

  treeFormulas.push_back(CutOptimizer::FormulaString("sum.trainingDeconvolvedFiltered().impulsivityMeasure/sum.trainingCoherentFiltered().impulsivityMeasure()", true));
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.trainingDeconvolvedFiltered().impulsivityMeasure", true));
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.trainingCoherentFiltered().impulsivityMeasure", true));

// | HigherPeakHilbertAfterDedispersion            | 31381542 | 34802404 |  651135 |
// | HigherImpulsivityMeasureAfterDedispersion     | 31358120 | 78339170 |  646652 |
// | LowerFracPowerWindowGradientAfterDedispersion |  8298605 | 22180410 |  648464 |
// | FisherScoreAboveThreshold                     |   646643 |   954698 | 8298605 |

  

  std::vector<const AnalysisCuts::AnalysisCut*> waisCuts;
  waisCuts.push_back(&AnalysisCuts::isTaggedAsWaisPulser);
  co.addSpectatorTree("waisTree", backgroundGlob, waisCuts);

  std::vector<const AnalysisCuts::AnalysisCut*> selectingBlastsCuts;
  selectingBlastsCuts.push_back(&AnalysisCuts::isTaggedAsPayloadBlast);
  co.addSpectatorTree("blastTree", backgroundGlob, selectingBlastsCuts);
  
  co.optimize(signalSelection, backgroundSelection, treeFormulas, outFileName);

  gSystem->Exit(0);
  return 0;
}
