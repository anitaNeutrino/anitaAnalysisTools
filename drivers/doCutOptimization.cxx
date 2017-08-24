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
  // signalSelection.push_back(&AnalysisCuts::isTaggedAsWaisPulser); // match timing criteria
  // signalSelection.push_back(&AnalysisCuts::isGood); // Remove payload blasts and SURF saturation
  signalSelection.push_back(&AnalysisCuts::goodGPS); // // Apparently there are a few stray NaNs
  signalSelection.push_back(&AnalysisCuts::nonZeroStokesI); // Apparently there are a few stray NaNs
  // signalSelection.push_back(&AnalysisCuts::realSNR); // Apparently there are a few stray NaNs
  signalSelection.push_back(&AnalysisCuts::closeToMC); // Need to get near other
  
  std::vector<const Acclaim::AnalysisCut *> backgroundSelection;
  backgroundSelection.push_back(&AnalysisCuts::isAboveHorizontal); // Upward pointing
  backgroundSelection.push_back(&AnalysisCuts::isGood); // Remove payload blasts and SURF saturation
  backgroundSelection.push_back(&AnalysisCuts::goodGPS); // Apparently there are a few stray NaNs
  backgroundSelection.push_back(&AnalysisCuts::nonZeroStokesI); // Apparently there are a few stray NaNs
  // backgroundSelection.push_back(&AnalysisCuts::realSNR); // Apparently there are a few stray NaNs
  backgroundSelection.push_back(&AnalysisCuts::anita3QuietTime); // quiet time

  // meta-data for reviewing results
  std::vector<CutOptimizer::FormulaString> treeFormulas;
  treeFormulas.push_back(CutOptimizer::FormulaString("run", false));
  treeFormulas.push_back(CutOptimizer::FormulaString("eventNumber", false)); // these get excluded from training by the CutOptimizer
  treeFormulas.push_back(CutOptimizer::FormulaString("realTime", false)); // these get excluded from training by the CutOptimizer
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.weight()", false)); // Automatically filled with 1 if non-MC
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.mc.energy", false)); // Automatically filled with 1 if non-MC


  treeFormulas.push_back(CutOptimizer::FormulaString("sum.mcPeak().hwAngle", true));  // definitionally can't do better than the hardware trigger, I think

  // map info
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.mcPeak().value", true)); // map values
  treeFormulas.push_back(CutOptimizer::FormulaString("TMath::Abs(sum.mcPeak().dPhiSource(wais))", true)); // delta phi sun

  // waveform info
  treeFormulas.push_back(CutOptimizer::FormulaString("TMath::Abs(sum.mcDeconvolved().impulsivityMeasure)", true));
  treeFormulas.push_back(CutOptimizer::FormulaString("TMath::SignBit(sum.mcDeconvolved().impulsivityMeasure)", true));
  treeFormulas.push_back(CutOptimizer::FormulaString("TMath::Abs(sum.mcDeconvolved().localMaxToMinTime)", true));
  treeFormulas.push_back(CutOptimizer::FormulaString("TMath::Abs(sum.mcDeconvolved().globalMaxToMinTime)", true));
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.mcDeconvolved().peakHilbert", true));
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.mcDeconvolved().linearPolFrac()", true));
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.mcDeconvolved().circPolFrac()", true));


  treeFormulas.push_back(CutOptimizer::FormulaString("sum.mcDeconvolved().totalPower/sum.mcDeconvolvedFiltered().totalPower", true));
  treeFormulas.push_back(CutOptimizer::FormulaString("sum.mcDeconvolved().peakHilbert/sum.mcDeconvolvedFiltered().peakHilbert", true));

  co.optimize(signalSelection, backgroundSelection, treeFormulas, outFileName);

  return 0;
}
