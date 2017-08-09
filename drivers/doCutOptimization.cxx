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

  CutOptimizer co(signalGlob, backgroundGlob, true);

  std::vector<const Acclaim::AnalysisCut *> signalSelection;
  signalSelection.push_back(&AnalysisCuts::isTaggedAsWaisPulser); // match timing criteria
  signalSelection.push_back(&AnalysisCuts::isGood); // Remove payload blasts and SURF saturation
  signalSelection.push_back(&AnalysisCuts::goodGPS); // // Apparently there are a few stray NaNs

  std::vector<const Acclaim::AnalysisCut *> backgroundSelection;
  backgroundSelection.push_back(&AnalysisCuts::isAboveHorizontal); // Upward pointing
  backgroundSelection.push_back(&AnalysisCuts::isGood); // Remove payload blasts and SURF saturation
  backgroundSelection.push_back(&AnalysisCuts::goodGPS); // Apparently there are a few stray NaNs

  std::vector<const char*> treeFormulas;
  treeFormulas.push_back("run");         // these get excluded from training by the CutOptimizer
  treeFormulas.push_back("eventNumber"); // these get excluded from training by the CutOptimizer

  treeFormulas.push_back("sum.higherPeak().value"); // map values

  // Waveform info
  treeFormulas.push_back("sum.higherCoherent().peakHilbert");
  treeFormulas.push_back("sum.higherDeconvolved().peakHilbert");
  treeFormulas.push_back("sum.higherCoherentFiltered().peakHilbert");
  treeFormulas.push_back("sum.higherDeconvolvedFiltered().peakHilbert");
  
  co.optimize(signalSelection, backgroundSelection, treeFormulas, outFileName);

  return 0;
}
